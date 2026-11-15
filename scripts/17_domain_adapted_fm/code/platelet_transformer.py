#!/usr/bin/env python3
"""
ACC 17B: Train Lightweight Transformer from Scratch on Platelet Data

Architecture: 6-layer transformer, 384 hidden dim, 6 heads (~20M params)
Pre-training: MGEM objective on platelet cells (48K) + X-Atlas perturbation data (500K subset)

This demonstrates the ability to design and train foundation models from scratch,
not just fine-tune existing ones.

Requirements:
    conda activate state_env

Usage:
    python platelet_transformer.py --epochs 100
    python platelet_transformer.py --epochs 100 --include-xatlas
    python platelet_transformer.py --epochs 100 --resume

Author: Xinru Qiu
Date: February 2026
"""

import os
import sys
import argparse
import logging
from pathlib import Path
import json
import math
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader, ConcatDataset
import scanpy as sc
from scipy.sparse import issparse

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


# =============================================================================
# Configuration
# =============================================================================

class Config:
    BASE_DIR = Path("/bigdata/godziklab/shared/Xinru/302006")
    DATA_DIR = BASE_DIR / "01_DATA"
    EMBEDDING_DIR = BASE_DIR / "02_EMBEDDINGS"
    RESULTS_DIR = BASE_DIR / "05_RESULTS/17_domain_adapted_fm/platelet_transformer"
    CHECKPOINT_DIR = BASE_DIR / "05_RESULTS/checkpoints/17_platelet_transformer"

    # X-Atlas data for pre-training augmentation
    XATLAS_DIR = DATA_DIR / "xatlas"

    # Model architecture
    N_GENES = 18903       # Total gene vocabulary size
    D_MODEL = 384         # Hidden dimension
    N_HEADS = 6           # Attention heads
    N_LAYERS = 6          # Transformer layers
    D_FF = 1536           # Feed-forward dimension (4x d_model)
    DROPOUT = 0.1
    MAX_SEQ_LEN = 2000    # Max genes per cell

    # Training
    BATCH_SIZE = 64
    LEARNING_RATE = 3e-4     # Reduced from 1e-3 to prevent NaN divergence
    EPOCHS = 100
    MASK_RATIO = 0.15
    WARMUP_EPOCHS = 10
    GRAD_ACCUM_STEPS = 4  # Effective batch size = 64 * 4 = 256

    SEED = 42
    INCLUDE_XATLAS = False
    XATLAS_MAX_CELLS = 500000  # Use 500K subset


# =============================================================================
# Model Architecture
# =============================================================================

class GeneEmbedding(nn.Module):
    """Gene embedding layer mapping gene indices to learned representations."""

    def __init__(self, n_genes, d_model, max_seq_len):
        super().__init__()
        self.gene_embedding = nn.Embedding(n_genes + 2, d_model, padding_idx=0)  # +2 for PAD and CLS
        self.expr_projection = nn.Linear(1, d_model)
        self.cls_token = nn.Parameter(torch.randn(1, 1, d_model) * 0.02)
        self.position_embedding = nn.Embedding(max_seq_len + 1, d_model)
        self.layer_norm = nn.LayerNorm(d_model)
        self.dropout = nn.Dropout(0.1)

    def forward(self, gene_ids, expr_values, lengths):
        """
        Args:
            gene_ids: [batch, seq_len] gene vocabulary indices
            expr_values: [batch, seq_len] expression values
            lengths: [batch] number of valid genes per cell
        Returns:
            embeddings: [batch, seq_len+1, d_model] (includes CLS token)
            padding_mask: [batch, seq_len+1] True for padding
        """
        batch_size, seq_len = gene_ids.shape

        # Gene embeddings + expression-weighted projection
        gene_emb = self.gene_embedding(gene_ids)  # [batch, seq, d_model]
        expr_proj = self.expr_projection(expr_values.unsqueeze(-1))  # [batch, seq, d_model]
        combined = gene_emb + expr_proj

        # Add CLS token
        cls_tokens = self.cls_token.expand(batch_size, -1, -1)
        combined = torch.cat([cls_tokens, combined], dim=1)  # [batch, seq+1, d_model]

        # Positional encoding
        positions = torch.arange(seq_len + 1, device=gene_ids.device).unsqueeze(0).expand(batch_size, -1)
        combined = combined + self.position_embedding(positions)

        combined = self.layer_norm(combined)
        combined = self.dropout(combined)

        # Padding mask (True = pad, for src_key_padding_mask)
        padding_mask = torch.arange(seq_len + 1, device=gene_ids.device).unsqueeze(0) > lengths.unsqueeze(1)

        return combined, padding_mask


class PlateletTransformer(nn.Module):
    """
    Lightweight transformer for single-cell gene expression.

    Architecture: 6-layer, 384 hidden, 6 heads (~20M params)
    Input: gene IDs + expression values
    Output: cell embedding (384 dims) from CLS token
    Pre-training: Masked gene expression modeling (MGEM)
    """

    def __init__(self, n_genes, d_model=384, n_heads=6, n_layers=6,
                 d_ff=1536, dropout=0.1, max_seq_len=2000):
        super().__init__()

        self.d_model = d_model
        self.embedding = GeneEmbedding(n_genes, d_model, max_seq_len)

        # Transformer encoder
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=n_heads,
            dim_feedforward=d_ff,
            dropout=dropout,
            activation='gelu',
            batch_first=True,
            norm_first=True  # Pre-norm for training stability
        )
        self.transformer = nn.TransformerEncoder(encoder_layer, num_layers=n_layers)

        # Pre-training head: predict masked expression values
        self.mgem_head = nn.Sequential(
            nn.Linear(d_model, d_model),
            nn.GELU(),
            nn.LayerNorm(d_model),
            nn.Linear(d_model, 1)
        )

        self._init_weights()

    def _init_weights(self):
        """Xavier uniform initialization."""
        for p in self.parameters():
            if p.dim() > 1:
                nn.init.xavier_uniform_(p)

    def forward(self, gene_ids, expr_values, lengths, return_embeddings=False):
        """
        Forward pass.

        Args:
            gene_ids: [batch, seq_len]
            expr_values: [batch, seq_len]
            lengths: [batch]
            return_embeddings: if True, return CLS embedding only

        Returns:
            If return_embeddings: cell_embedding [batch, d_model]
            Else: predicted_values [batch, seq_len]
        """
        # Embed genes
        x, padding_mask = self.embedding(gene_ids, expr_values, lengths)

        # Transformer
        x = self.transformer(x, src_key_padding_mask=padding_mask)

        if return_embeddings:
            # Return CLS token embedding
            return x[:, 0, :]  # [batch, d_model]

        # MGEM: predict expression values for all positions (loss applied on masked only)
        predictions = self.mgem_head(x[:, 1:, :]).squeeze(-1)  # [batch, seq_len]
        return predictions

    def get_cell_embeddings(self, gene_ids, expr_values, lengths):
        """Extract cell embeddings (CLS token)."""
        return self.forward(gene_ids, expr_values, lengths, return_embeddings=True)


# =============================================================================
# Dataset
# =============================================================================

class CellExpressionDataset(Dataset):
    """Dataset for cell expression data with MGEM masking."""

    def __init__(self, adata, gene_vocab, max_seq_len=2000, mask_ratio=0.15, is_pretrain=True):
        self.gene_vocab = gene_vocab
        self.max_seq_len = max_seq_len
        self.mask_ratio = mask_ratio
        self.is_pretrain = is_pretrain

        # Extract expression matrix
        if issparse(adata.X):
            self.X = adata.X.toarray().astype(np.float32)
        else:
            self.X = np.array(adata.X, dtype=np.float32)

        self.gene_names = [g.upper() for g in adata.var_names]
        self.n_cells = self.X.shape[0]

        # Pre-compute gene-to-vocab mapping
        self.gene_to_vocab = {}
        for i, gene in enumerate(self.gene_names):
            if gene in gene_vocab:
                self.gene_to_vocab[i] = gene_vocab[gene]

        logger.info(f"Dataset: {self.n_cells} cells, {len(self.gene_to_vocab)}/{len(self.gene_names)} genes in vocab")

    def __len__(self):
        return self.n_cells

    def __getitem__(self, idx):
        expr = self.X[idx]

        # Get expressed genes sorted by expression (descending)
        nonzero = np.where(expr > 0)[0]
        if len(nonzero) == 0:
            nonzero = np.arange(min(10, len(expr)))

        sorted_idx = nonzero[np.argsort(-expr[nonzero])]

        # Map to vocabulary and truncate
        gene_ids = []
        values = []
        for i in sorted_idx:
            if i in self.gene_to_vocab:
                gene_ids.append(self.gene_to_vocab[i])
                values.append(float(expr[i]))
                if len(gene_ids) >= self.max_seq_len:
                    break

        n_genes = len(gene_ids)
        if n_genes == 0:
            gene_ids = [1]  # CLS placeholder
            values = [0.0]
            n_genes = 1

        # Create MGEM mask
        mask_labels = [-100.0] * n_genes
        masked_values = list(values)

        if self.is_pretrain:
            n_mask = max(1, int(n_genes * self.mask_ratio))
            mask_pos = np.random.choice(n_genes, size=n_mask, replace=False)
            for pos in mask_pos:
                mask_labels[pos] = values[pos]
                r = np.random.random()
                if r < 0.8:
                    masked_values[pos] = 0.0
                elif r < 0.9:
                    masked_values[pos] = float(np.random.uniform(0, max(values) + 1))

        # Pad
        pad_len = self.max_seq_len - n_genes
        if pad_len > 0:
            gene_ids.extend([0] * pad_len)
            masked_values.extend([0.0] * pad_len)
            mask_labels.extend([-100.0] * pad_len)

        return {
            'gene_ids': torch.tensor(gene_ids[:self.max_seq_len], dtype=torch.long),
            'expr_values': torch.tensor(masked_values[:self.max_seq_len], dtype=torch.float32),
            'true_values': torch.tensor(values[:self.max_seq_len] + [0.0] * max(0, pad_len), dtype=torch.float32),
            'mask_labels': torch.tensor(mask_labels[:self.max_seq_len], dtype=torch.float32),
            'length': torch.tensor(n_genes, dtype=torch.long),
        }


def build_gene_vocab(adata_list):
    """Build gene vocabulary from one or more AnnData objects."""
    gene_set = set()
    for adata in adata_list:
        gene_set.update(g.upper() for g in adata.var_names)

    # Reserve 0=PAD, 1=CLS
    vocab = {gene: idx + 2 for idx, gene in enumerate(sorted(gene_set))}
    logger.info(f"Gene vocabulary: {len(vocab)} genes")
    return vocab


def load_xatlas_data(config: Config, max_cells=500000):
    """Load X-Atlas perturbation data for pre-training augmentation."""
    xatlas_dir = config.XATLAS_DIR
    if not xatlas_dir.exists():
        logger.warning(f"X-Atlas directory not found: {xatlas_dir}")
        return None

    chunks = sorted(xatlas_dir.glob("*.h5ad"))
    if not chunks:
        logger.warning("No X-Atlas chunks found")
        return None

    logger.info(f"Loading X-Atlas data ({len(chunks)} chunks, max {max_cells} cells)...")
    adatas = []
    total_cells = 0

    for chunk_file in chunks:
        if total_cells >= max_cells:
            break
        try:
            adata = sc.read_h5ad(chunk_file)
            remaining = max_cells - total_cells
            if adata.n_obs > remaining:
                adata = adata[:remaining].copy()
            adatas.append(adata)
            total_cells += adata.n_obs
            logger.info(f"  {chunk_file.name}: {adata.n_obs} cells (total: {total_cells})")
        except Exception as e:
            logger.warning(f"  Failed to load {chunk_file.name}: {e}")

    if adatas:
        combined = sc.concat(adatas)
        logger.info(f"X-Atlas combined: {combined.n_obs} cells, {combined.n_vars} genes")
        return combined

    return None


# =============================================================================
# Training
# =============================================================================

def run_train_from_scratch(config: Config, resume: bool = False):
    """Train PlateletTransformer from scratch."""

    logger.info("=" * 60)
    logger.info("ACC 17B: Training Platelet Transformer from Scratch")
    logger.info(f"Architecture: {config.N_LAYERS}L, {config.D_MODEL}d, {config.N_HEADS}h")
    logger.info("=" * 60)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logger.info(f"Device: {device}")
    if device.type == 'cuda':
        logger.info(f"GPU: {torch.cuda.get_device_name(0)}")

    config.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    config.CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)

    torch.manual_seed(config.SEED)
    np.random.seed(config.SEED)

    # Load platelet data
    possible_files = [
        config.EMBEDDING_DIR / "state/covid_state.h5ad",
        config.DATA_DIR / "covid_pbmc_adata.h5ad",
    ]
    platelet_adata = None
    for f in possible_files:
        if f.exists():
            platelet_adata = sc.read_h5ad(f)
            logger.info(f"Loaded platelet data: {platelet_adata.shape}")
            break
    if platelet_adata is None:
        raise FileNotFoundError("Could not find platelet data")

    # Optionally load X-Atlas data
    xatlas_adata = None
    if config.INCLUDE_XATLAS:
        xatlas_adata = load_xatlas_data(config, config.XATLAS_MAX_CELLS)

    # Build gene vocabulary
    adata_list = [platelet_adata]
    if xatlas_adata is not None:
        adata_list.append(xatlas_adata)
    gene_vocab = build_gene_vocab(adata_list)
    n_genes = len(gene_vocab) + 2  # +2 for PAD, CLS

    # Create datasets
    platelet_dataset = CellExpressionDataset(
        platelet_adata, gene_vocab, config.MAX_SEQ_LEN, config.MASK_RATIO
    )

    datasets = [platelet_dataset]
    if xatlas_adata is not None:
        xatlas_dataset = CellExpressionDataset(
            xatlas_adata, gene_vocab, config.MAX_SEQ_LEN, config.MASK_RATIO
        )
        datasets.append(xatlas_dataset)
        logger.info(f"Combined dataset: {len(platelet_dataset)} platelet + {len(xatlas_dataset)} X-Atlas cells")

    combined_dataset = ConcatDataset(datasets) if len(datasets) > 1 else datasets[0]

    dataloader = DataLoader(
        combined_dataset, batch_size=config.BATCH_SIZE,
        shuffle=True, num_workers=4, pin_memory=True, drop_last=True
    )

    # Initialize model
    model = PlateletTransformer(
        n_genes=n_genes,
        d_model=config.D_MODEL,
        n_heads=config.N_HEADS,
        n_layers=config.N_LAYERS,
        d_ff=config.D_FF,
        dropout=config.DROPOUT,
        max_seq_len=config.MAX_SEQ_LEN,
    ).to(device)

    n_params = sum(p.numel() for p in model.parameters())
    logger.info(f"Model parameters: {n_params:,} ({n_params/1e6:.1f}M)")

    # Optimizer
    optimizer = torch.optim.AdamW(model.parameters(), lr=config.LEARNING_RATE, weight_decay=0.01)

    # Cosine schedule with warmup
    total_steps = config.EPOCHS * len(dataloader)
    warmup_steps = config.WARMUP_EPOCHS * len(dataloader)

    def lr_lambda(step):
        if step < warmup_steps:
            return step / max(1, warmup_steps)
        progress = (step - warmup_steps) / max(1, total_steps - warmup_steps)
        return max(0.01, 0.5 * (1 + math.cos(math.pi * progress)))

    scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda)

    # Resume
    start_epoch = 0
    best_loss = float('inf')
    if resume:
        ckpt_path = config.CHECKPOINT_DIR / "platelet_tf_latest.pt"
        if ckpt_path.exists():
            ckpt = torch.load(ckpt_path, map_location=device)
            model.load_state_dict(ckpt['model_state_dict'])
            optimizer.load_state_dict(ckpt['optimizer_state_dict'])
            start_epoch = ckpt['epoch']
            best_loss = ckpt.get('best_loss', float('inf'))
            logger.info(f"Resumed from epoch {start_epoch}")

    # Training
    scaler = torch.amp.GradScaler('cuda') if device.type == 'cuda' else None
    loss_history = []

    for epoch in range(start_epoch, config.EPOCHS):
        model.train()
        epoch_loss = 0
        n_batches = 0

        optimizer.zero_grad()

        for step, batch in enumerate(dataloader):
            gene_ids = batch['gene_ids'].to(device)
            expr_values = batch['expr_values'].to(device)
            mask_labels = batch['mask_labels'].to(device)
            lengths = batch['length'].to(device)

            if scaler is not None:
                with torch.amp.autocast('cuda'):
                    predictions = model(gene_ids, expr_values, lengths)
                    # Loss on masked positions only
                    mask = mask_labels != -100.0
                    if mask.any():
                        loss = F.mse_loss(predictions[mask], mask_labels[mask])
                    else:
                        continue
                    loss = loss / config.GRAD_ACCUM_STEPS

                scaler.scale(loss).backward()

                if (step + 1) % config.GRAD_ACCUM_STEPS == 0:
                    scaler.unscale_(optimizer)
                    torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                    scaler.step(optimizer)
                    scaler.update()
                    optimizer.zero_grad()
                    scheduler.step()
            else:
                predictions = model(gene_ids, expr_values, lengths)
                mask = mask_labels != -100.0
                if mask.any():
                    loss = F.mse_loss(predictions[mask], mask_labels[mask])
                else:
                    continue
                loss = loss / config.GRAD_ACCUM_STEPS
                loss.backward()

                if (step + 1) % config.GRAD_ACCUM_STEPS == 0:
                    torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                    optimizer.step()
                    optimizer.zero_grad()
                    scheduler.step()

            epoch_loss += loss.item() * config.GRAD_ACCUM_STEPS
            n_batches += 1

        avg_loss = epoch_loss / max(n_batches, 1)
        loss_history.append(avg_loss)
        lr = optimizer.param_groups[0]['lr']

        logger.info(f"Epoch {epoch+1}/{config.EPOCHS} | Loss: {avg_loss:.4f} | LR: {lr:.2e}")

        # NaN detection: stop early to avoid wasting GPU hours
        if math.isnan(avg_loss) or math.isinf(avg_loss):
            logger.error(f"Training diverged at epoch {epoch+1} (loss={avg_loss}). "
                        f"Best loss before divergence: {best_loss:.4f}. Stopping.")
            break

        # Save best
        if avg_loss < best_loss:
            best_loss = avg_loss
            torch.save({
                'epoch': epoch + 1,
                'model_state_dict': model.state_dict(),
                'config': {
                    'n_genes': n_genes,
                    'd_model': config.D_MODEL,
                    'n_heads': config.N_HEADS,
                    'n_layers': config.N_LAYERS,
                    'd_ff': config.D_FF,
                    'max_seq_len': config.MAX_SEQ_LEN,
                },
                'gene_vocab': gene_vocab,
                'best_loss': best_loss,
            }, config.CHECKPOINT_DIR / "platelet_tf_best.pt")
            logger.info(f"  Saved best checkpoint (loss={best_loss:.4f})")

        # Save latest
        if (epoch + 1) % 10 == 0 or epoch == config.EPOCHS - 1:
            torch.save({
                'epoch': epoch + 1,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'config': {
                    'n_genes': n_genes,
                    'd_model': config.D_MODEL,
                    'n_heads': config.N_HEADS,
                    'n_layers': config.N_LAYERS,
                    'd_ff': config.D_FF,
                    'max_seq_len': config.MAX_SEQ_LEN,
                },
                'gene_vocab': gene_vocab,
                'best_loss': best_loss,
            }, config.CHECKPOINT_DIR / "platelet_tf_latest.pt")

    # Save history
    history = {
        'loss_history': loss_history,
        'best_loss': best_loss,
        'n_params': n_params,
        'n_cells_platelet': len(platelet_dataset),
        'n_cells_xatlas': len(xatlas_dataset) if xatlas_adata is not None else 0,
        'n_genes_vocab': len(gene_vocab),
    }
    with open(config.RESULTS_DIR / "train_history.json", 'w') as f:
        json.dump(history, f, indent=2)

    logger.info(f"\nTraining complete. Best loss: {best_loss:.4f}")
    logger.info(f"Parameters: {n_params/1e6:.1f}M")

    return model, gene_vocab


def main():
    parser = argparse.ArgumentParser(description='Train Platelet Transformer from Scratch')
    parser.add_argument('--epochs', type=int, default=100)
    parser.add_argument('--batch-size', type=int, default=64)
    parser.add_argument('--lr', type=float, default=1e-3)
    parser.add_argument('--d-model', type=int, default=384)
    parser.add_argument('--n-layers', type=int, default=6)
    parser.add_argument('--n-heads', type=int, default=6)
    parser.add_argument('--include-xatlas', action='store_true', default=False)
    parser.add_argument('--no-xatlas', action='store_true')
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()

    config = Config()
    config.EPOCHS = args.epochs
    config.BATCH_SIZE = args.batch_size
    config.LEARNING_RATE = args.lr
    config.D_MODEL = args.d_model
    config.N_LAYERS = args.n_layers
    config.N_HEADS = args.n_heads
    config.D_FF = args.d_model * 4
    config.INCLUDE_XATLAS = not args.no_xatlas

    run_train_from_scratch(config, resume=args.resume)


if __name__ == '__main__':
    main()
