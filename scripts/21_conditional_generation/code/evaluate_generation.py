#!/usr/bin/env python3
"""
ACC 21: Evaluate Generated Embeddings

Post-hoc evaluation of conditional VAE generations:
1. UMAP overlay: generated vs real recovered cells
2. Compare cVAE transition vectors to real X-Atlas perturbation effects
3. Generate transition plots for all severity pairs

Usage:
    python evaluate_generation.py

Author: Xinru Qiu
Date: February 2026
"""

import logging
from pathlib import Path
import json
import numpy as np
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

BASE_DIR = Path("/bigdata/godziklab/shared/Xinru/302006")
EMB_DIR = BASE_DIR / "02_EMBEDDINGS"
RESULTS_DIR = BASE_DIR / "05_RESULTS/21_conditional_generation"
FIGURES_DIR = RESULTS_DIR / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

SEVERITY_COLORS = {
    'healthy': '#2ecc71', 'mild': '#3498db', 'moderate': '#f39c12',
    'severe': '#e74c3c', 'fatal': '#8e44ad', 'recovered': '#1abc9c',
    'generated': '#e67e22',
}


def plot_umap_overlay():
    """Create UMAP with generated recovered cells overlaid."""
    logger.info("Creating UMAP overlay...")

    # Load STATE embeddings
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
    valid_severities = ['healthy', 'mild', 'moderate', 'severe', 'fatal', 'recovered']
    valid_mask = np.isin(severity, valid_severities)
    X = X[valid_mask]
    severity = severity[valid_mask]

    # Load generated embeddings
    gen_path = RESULTS_DIR / "generated_recovered_embeddings.npy"
    if not gen_path.exists():
        logger.warning(f"Generated embeddings not found: {gen_path}")
        return

    X_generated = np.load(gen_path)
    logger.info(f"Generated embeddings: {X_generated.shape}")

    # Combine for UMAP
    X_combined = np.vstack([X, X_generated])
    labels = np.concatenate([severity, ['generated'] * len(X_generated)])

    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_combined)

    # UMAP
    import umap
    reducer = umap.UMAP(n_neighbors=30, min_dist=0.3, random_state=42)
    embedding_2d = reducer.fit_transform(X_scaled)

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # Panel A: All real cells
    ax = axes[0]
    real_mask = labels != 'generated'
    for sev in valid_severities:
        mask = labels == sev
        ax.scatter(embedding_2d[mask, 0], embedding_2d[mask, 1],
                  c=SEVERITY_COLORS[sev], label=sev, s=5, alpha=0.5)
    ax.set_title('Real Cell Embeddings')
    ax.legend(markerscale=4, fontsize=9)
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')

    # Panel B: With generated cells
    ax = axes[1]
    # Plot real cells faded
    ax.scatter(embedding_2d[real_mask, 0], embedding_2d[real_mask, 1],
              c='lightgray', s=3, alpha=0.2, label='Real cells')
    # Highlight severe (source)
    severe_mask = labels == 'severe'
    ax.scatter(embedding_2d[severe_mask, 0], embedding_2d[severe_mask, 1],
              c=SEVERITY_COLORS['severe'], s=8, alpha=0.6, label='Severe (source)')
    # Highlight real recovered
    rec_mask = labels == 'recovered'
    ax.scatter(embedding_2d[rec_mask, 0], embedding_2d[rec_mask, 1],
              c=SEVERITY_COLORS['recovered'], s=8, alpha=0.6, label='Recovered (real)')
    # Highlight generated
    gen_mask = labels == 'generated'
    ax.scatter(embedding_2d[gen_mask, 0], embedding_2d[gen_mask, 1],
              c=SEVERITY_COLORS['generated'], s=15, alpha=0.8,
              marker='*', label='Generated recovered')

    ax.set_title('Generated vs Real Recovered')
    ax.legend(markerscale=3, fontsize=9)
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')

    plt.suptitle('Conditional VAE: Severe → Recovered Transition', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'umap_generated_overlay.png', dpi=150, bbox_inches='tight')
    plt.savefig(FIGURES_DIR / 'umap_generated_overlay.pdf', bbox_inches='tight')
    plt.close()
    logger.info(f"Saved UMAP overlay to {FIGURES_DIR}")


def plot_training_curves():
    """Plot cVAE training loss curves."""
    history_files = list(RESULTS_DIR.glob("cvae_training_history_*.json"))
    if not history_files:
        logger.warning("No training history found")
        return

    with open(sorted(history_files)[-1]) as f:
        history = json.load(f)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(history['train_losses'], label='Train', linewidth=2, color='#2196F3')
    ax.plot(history['val_losses'], label='Validation', linewidth=2, color='#F44336')
    ax.set_xlabel('Epoch')
    ax.set_ylabel('Loss (Recon + β·KL)')
    ax.set_title(f'Conditional VAE Training (β={history["beta"]}, z={history["latent_dim"]})')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'cvae_training_curves.png', dpi=150, bbox_inches='tight')
    plt.close()
    logger.info("Saved training curves")


def plot_latent_space():
    """Plot latent space of the cVAE colored by severity."""
    latent_path = RESULTS_DIR / "latent_z_severe.npy"
    if not latent_path.exists():
        logger.warning("Latent embeddings not found")
        return

    Z = np.load(latent_path)

    # 2D projection of latent space
    if Z.shape[1] > 2:
        from sklearn.decomposition import PCA
        pca = PCA(n_components=2, random_state=42)
        Z_2d = pca.fit_transform(Z)
    else:
        Z_2d = Z

    fig, ax = plt.subplots(figsize=(8, 7))
    ax.scatter(Z_2d[:, 0], Z_2d[:, 1], c=SEVERITY_COLORS['severe'], s=10, alpha=0.5)
    ax.set_xlabel('Latent dim 1')
    ax.set_ylabel('Latent dim 2')
    ax.set_title('Latent Space of Severe Cells')

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'latent_space.png', dpi=150, bbox_inches='tight')
    plt.close()
    logger.info("Saved latent space plot")


def summarize_results():
    """Print summary of generation evaluation."""
    results_file = RESULTS_DIR / "transition_results.json"
    if not results_file.exists():
        logger.warning("Transition results not found")
        return

    with open(results_file) as f:
        results = json.load(f)

    logger.info("\n" + "=" * 60)
    logger.info("GENERATION EVALUATION SUMMARY")
    logger.info("=" * 60)
    logger.info(f"  Severe cells encoded: {results['n_severe_cells']}")
    logger.info(f"  Real recovered cells: {results['n_recovered_cells']}")
    logger.info(f"  Mean L2 to nearest recovered: {results['mean_l2_to_nearest_recovered']:.4f}")
    logger.info(f"  Cosine sim (recovery direction): {results['cosine_sim_recovery_direction']:.4f}")
    logger.info(f"  Classified as recovered: {results['frac_classified_recovered']:.1%}")
    logger.info(f"  Mean P(recovered): {results['mean_prob_recovered']:.4f}")

    # Quality assessment
    quality = "GOOD" if results['cosine_sim_recovery_direction'] > 0.5 else "MODERATE"
    if results['frac_classified_recovered'] > 0.7:
        quality = "EXCELLENT"
    logger.info(f"\n  Overall quality: {quality}")


def main():
    logger.info("ACC 21: Evaluating Generated Embeddings")

    plot_training_curves()
    plot_umap_overlay()
    plot_latent_space()
    summarize_results()

    logger.info("\nEvaluation complete.")


if __name__ == '__main__':
    main()
