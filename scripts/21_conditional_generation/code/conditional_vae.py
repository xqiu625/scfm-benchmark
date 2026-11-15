#!/usr/bin/env python3
"""
ACC 21: Conditional VAE for Severity Transition Generation

Train a conditional VAE on STATE embeddings to model severity transitions.
Key experiment: generate "recovered" cells from severe cell latent codes.

Architecture:
    Encoder: [2058 + 6 (condition)] -> 512 -> 256 -> 128 (mu, logvar)
    Decoder: [128 + 6 (condition)] -> 256 -> 512 -> 2058

Training: beta-VAE loss (reconstruction + 0.5 * KL divergence)

Usage:
    python conditional_vae.py
    python conditional_vae.py --epochs 200 --beta 0.5

Author: Xinru Qiu
Date: February 2026
"""

import argparse
import logging
from pathlib import Path
import json
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader, TensorDataset
import scanpy as sc
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, StandardScaler
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

BASE_DIR = Path("/bigdata/godziklab/shared/Xinru/302006")
EMB_DIR = BASE_DIR / "02_EMBEDDINGS"
RESULTS_DIR = BASE_DIR / "05_RESULTS/21_conditional_generation"
FIGURES_DIR = RESULTS_DIR / "figures"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Severity classes for conditioning
SEVERITY_CLASSES = ['healthy', 'mild', 'moderate', 'severe', 'fatal', 'recovered']
N_CONDITIONS = len(SEVERITY_CLASSES)


class ConditionalVAE(nn.Module):
    """
    Conditional Variational Autoencoder for cell embeddings.

    Conditions on severity class to enable controlled generation
    of cell states in embedding space.
    """

    def __init__(self, input_dim=2058, latent_dim=128, n_conditions=6,
                 hidden_dims=[512, 256]):
        super().__init__()

        self.input_dim = input_dim
        self.latent_dim = latent_dim
        self.n_conditions = n_conditions

        # Encoder: input + condition -> latent
        encoder_input = input_dim + n_conditions
        self.encoder = nn.Sequential(
            nn.Linear(encoder_input, hidden_dims[0]),
            nn.BatchNorm1d(hidden_dims[0]),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(hidden_dims[0], hidden_dims[1]),
            nn.BatchNorm1d(hidden_dims[1]),
            nn.ReLU(),
            nn.Dropout(0.1),
        )
        self.fc_mu = nn.Linear(hidden_dims[1], latent_dim)
        self.fc_logvar = nn.Linear(hidden_dims[1], latent_dim)

        # Decoder: latent + condition -> reconstruction
        decoder_input = latent_dim + n_conditions
        self.decoder = nn.Sequential(
            nn.Linear(decoder_input, hidden_dims[1]),
            nn.BatchNorm1d(hidden_dims[1]),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(hidden_dims[1], hidden_dims[0]),
            nn.BatchNorm1d(hidden_dims[0]),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(hidden_dims[0], input_dim),
        )

    def encode(self, x, condition):
        """Encode input + condition to latent distribution."""
        h = torch.cat([x, condition], dim=1)
        h = self.encoder(h)
        mu = self.fc_mu(h)
        logvar = self.fc_logvar(h)
        return mu, logvar

    def reparameterize(self, mu, logvar):
        """Reparameterization trick."""
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def decode(self, z, condition):
        """Decode latent + condition to reconstruction."""
        h = torch.cat([z, condition], dim=1)
        return self.decoder(h)

    def forward(self, x, condition):
        """Full forward pass: encode, sample, decode."""
        mu, logvar = self.encode(x, condition)
        z = self.reparameterize(mu, logvar)
        recon = self.decode(z, condition)
        return recon, mu, logvar, z

    def generate(self, z, condition):
        """Generate from latent code + condition."""
        return self.decode(z, condition)

    def transition(self, x, source_condition, target_condition):
        """
        Generate transition: encode with source condition,
        decode with target condition.
        """
        mu, logvar = self.encode(x, source_condition)
        z = mu  # Use mean for deterministic transition
        return self.decode(z, target_condition), z


def vae_loss(recon, x, mu, logvar, beta=0.5):
    """Beta-VAE loss: reconstruction + beta * KL divergence."""
    recon_loss = F.mse_loss(recon, x, reduction='mean')
    kl_loss = -0.5 * torch.mean(1 + logvar - mu.pow(2) - logvar.exp())
    return recon_loss + beta * kl_loss, recon_loss, kl_loss


def load_state_embeddings():
    """Load STATE embeddings with severity labels."""
    covid_path = EMB_DIR / "state/covid_state.h5ad"
    adata = sc.read_h5ad(covid_path)

    X = None
    for key in ['X_state', 'X_STATE', 'X_emb']:
        if key in adata.obsm:
            X = np.array(adata.obsm[key])
            break
    if X is None:
        X = adata.X.toarray() if hasattr(adata.X, 'toarray') else np.array(adata.X)

    severity = adata.obs['covid_severity'].astype(str).str.lower().values

    # Filter valid severities
    valid_mask = np.isin(severity, SEVERITY_CLASSES)
    X = X[valid_mask]
    severity = severity[valid_mask]

    logger.info(f"Loaded STATE embeddings: {X.shape}")
    for s in SEVERITY_CLASSES:
        logger.info(f"  {s}: {np.sum(severity == s)}")

    return X, severity


def severity_to_onehot(severity_labels, classes=SEVERITY_CLASSES):
    """Convert severity labels to one-hot vectors."""
    class_to_idx = {c: i for i, c in enumerate(classes)}
    indices = np.array([class_to_idx[s] for s in severity_labels])
    onehot = np.zeros((len(indices), len(classes)), dtype=np.float32)
    onehot[np.arange(len(indices)), indices] = 1.0
    return onehot, indices


def train_cvae(epochs=200, beta=0.5, latent_dim=128, batch_size=256, lr=1e-3):
    """Train conditional VAE."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    logger.info("=" * 60)
    logger.info("ACC 21: Conditional VAE Training")
    logger.info(f"Epochs: {epochs}, Beta: {beta}, Latent: {latent_dim}")
    logger.info("=" * 60)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logger.info(f"Device: {device}")

    # Load data
    X, severity = load_state_embeddings()
    conditions, class_indices = severity_to_onehot(severity)

    # Standardize embeddings
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X).astype(np.float32)

    # Train/val split
    X_train, X_val, c_train, c_val, sev_train, sev_val = train_test_split(
        X_scaled, conditions, severity, test_size=0.2, stratify=severity, random_state=42
    )

    train_dataset = TensorDataset(
        torch.FloatTensor(X_train),
        torch.FloatTensor(c_train),
    )
    val_dataset = TensorDataset(
        torch.FloatTensor(X_val),
        torch.FloatTensor(c_val),
    )

    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, drop_last=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)

    # Model
    input_dim = X_scaled.shape[1]
    model = ConditionalVAE(
        input_dim=input_dim,
        latent_dim=latent_dim,
        n_conditions=N_CONDITIONS,
    ).to(device)

    n_params = sum(p.numel() for p in model.parameters())
    logger.info(f"Model parameters: {n_params:,}")

    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs, eta_min=1e-5)

    # Training loop
    train_losses = []
    val_losses = []

    for epoch in range(epochs):
        model.train()
        epoch_loss = 0
        epoch_recon = 0
        epoch_kl = 0
        n_batches = 0

        for x_batch, c_batch in train_loader:
            x_batch = x_batch.to(device)
            c_batch = c_batch.to(device)

            optimizer.zero_grad()
            recon, mu, logvar, z = model(x_batch, c_batch)
            loss, recon_loss, kl_loss = vae_loss(recon, x_batch, mu, logvar, beta)

            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()

            epoch_loss += loss.item()
            epoch_recon += recon_loss.item()
            epoch_kl += kl_loss.item()
            n_batches += 1

        scheduler.step()

        avg_loss = epoch_loss / n_batches
        avg_recon = epoch_recon / n_batches
        avg_kl = epoch_kl / n_batches
        train_losses.append(avg_loss)

        # Validation
        model.eval()
        val_loss_total = 0
        val_n = 0
        with torch.no_grad():
            for x_batch, c_batch in val_loader:
                x_batch = x_batch.to(device)
                c_batch = c_batch.to(device)
                recon, mu, logvar, z = model(x_batch, c_batch)
                loss, _, _ = vae_loss(recon, x_batch, mu, logvar, beta)
                val_loss_total += loss.item()
                val_n += 1

        val_avg = val_loss_total / val_n
        val_losses.append(val_avg)

        if (epoch + 1) % 20 == 0:
            logger.info(f"Epoch {epoch+1}/{epochs} | Train: {avg_loss:.4f} "
                       f"(recon={avg_recon:.4f}, kl={avg_kl:.4f}) | Val: {val_avg:.4f}")

    # Save model
    torch.save({
        'model_state_dict': model.state_dict(),
        'scaler_mean': scaler.mean_,
        'scaler_scale': scaler.scale_,
        'input_dim': input_dim,
        'latent_dim': latent_dim,
        'train_losses': train_losses,
        'val_losses': val_losses,
    }, RESULTS_DIR / f"cvae_model_{timestamp}.pt")

    # Save training history
    with open(RESULTS_DIR / f"cvae_training_history_{timestamp}.json", 'w') as f:
        json.dump({
            'train_losses': train_losses,
            'val_losses': val_losses,
            'epochs': epochs,
            'beta': beta,
            'latent_dim': latent_dim,
            'n_params': n_params,
        }, f, indent=2)

    logger.info(f"Model saved to {RESULTS_DIR}")
    return model, scaler, X_scaled, severity, device


def run_severity_transition(model, scaler, X_scaled, severity, device):
    """
    Key experiment: generate "recovered" cells from severe cells.
    """
    logger.info("\n" + "=" * 60)
    logger.info("SEVERITY TRANSITION EXPERIMENT")
    logger.info("Severe -> Recovered transition")
    logger.info("=" * 60)

    model.eval()

    # Get severe cells
    severe_mask = severity == 'severe'
    recovered_mask = severity == 'recovered'

    X_severe = torch.FloatTensor(X_scaled[severe_mask]).to(device)
    X_recovered_real = X_scaled[recovered_mask]

    # Create conditions
    severe_cond = torch.zeros(X_severe.shape[0], N_CONDITIONS, device=device)
    severe_cond[:, SEVERITY_CLASSES.index('severe')] = 1.0

    recovered_cond = torch.zeros(X_severe.shape[0], N_CONDITIONS, device=device)
    recovered_cond[:, SEVERITY_CLASSES.index('recovered')] = 1.0

    # Generate transition
    with torch.no_grad():
        generated_recovered, latent_z = model.transition(X_severe, severe_cond, recovered_cond)
        generated_recovered = generated_recovered.cpu().numpy()
        latent_z = latent_z.cpu().numpy()

    # Evaluate quality
    from sklearn.neighbors import NearestNeighbors
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score

    # 1. L2 distance to nearest real recovered cell
    nn = NearestNeighbors(n_neighbors=1)
    nn.fit(X_recovered_real)
    distances, _ = nn.kneighbors(generated_recovered)
    mean_l2 = float(np.mean(distances))
    logger.info(f"  Mean L2 to nearest real recovered: {mean_l2:.4f}")

    # 2. Cosine similarity of transition vector to mean recovery direction
    severe_centroid = np.mean(X_scaled[severe_mask], axis=0)
    recovered_centroid = np.mean(X_recovered_real, axis=0)
    real_direction = recovered_centroid - severe_centroid
    real_direction = real_direction / (np.linalg.norm(real_direction) + 1e-8)

    generated_centroid = np.mean(generated_recovered, axis=0)
    gen_direction = generated_centroid - severe_centroid
    gen_direction = gen_direction / (np.linalg.norm(gen_direction) + 1e-8)

    cosine_sim = float(np.dot(real_direction, gen_direction))
    logger.info(f"  Cosine similarity (gen vs real recovery direction): {cosine_sim:.4f}")

    # 3. Classifier test: LogReg trained on real data classifies generated
    X_real_all = np.vstack([X_scaled[severe_mask], X_recovered_real])
    y_real = np.array([0] * np.sum(severe_mask) + [1] * np.sum(recovered_mask))

    clf = LogisticRegression(max_iter=2000, random_state=42)
    clf.fit(X_real_all, y_real)

    gen_probs = clf.predict_proba(generated_recovered)[:, 1]
    gen_preds = clf.predict(generated_recovered)
    frac_classified_recovered = float(np.mean(gen_preds == 1))
    logger.info(f"  Fraction classified as recovered: {frac_classified_recovered:.4f}")
    logger.info(f"  Mean probability(recovered): {float(np.mean(gen_probs)):.4f}")

    results = {
        'mean_l2_to_nearest_recovered': mean_l2,
        'cosine_sim_recovery_direction': cosine_sim,
        'frac_classified_recovered': frac_classified_recovered,
        'mean_prob_recovered': float(np.mean(gen_probs)),
        'n_severe_cells': int(np.sum(severe_mask)),
        'n_recovered_cells': int(np.sum(recovered_mask)),
    }

    with open(RESULTS_DIR / "transition_results.json", 'w') as f:
        json.dump(results, f, indent=2)

    # Save generated embeddings for UMAP visualization
    np.save(RESULTS_DIR / "generated_recovered_embeddings.npy", generated_recovered)
    np.save(RESULTS_DIR / "latent_z_severe.npy", latent_z)

    logger.info(f"  Results saved to {RESULTS_DIR}")
    return results


def main():
    parser = argparse.ArgumentParser(description='Conditional VAE')
    parser.add_argument('--epochs', type=int, default=200)
    parser.add_argument('--beta', type=float, default=0.5)
    parser.add_argument('--latent-dim', type=int, default=128)
    parser.add_argument('--batch-size', type=int, default=256)
    parser.add_argument('--lr', type=float, default=1e-3)
    args = parser.parse_args()

    torch.manual_seed(42)
    np.random.seed(42)

    model, scaler, X_scaled, severity, device = train_cvae(
        epochs=args.epochs, beta=args.beta,
        latent_dim=args.latent_dim, batch_size=args.batch_size, lr=args.lr
    )

    run_severity_transition(model, scaler, X_scaled, severity, device)


if __name__ == '__main__':
    main()
