#!/usr/bin/env python3
"""
ACC 17A: Continue Pre-Training scGPT on Platelet Data

Domain-adapt the scGPT whole-human checkpoint to platelet transcriptomics
using masked gene expression modeling (MGEM).

Approach:
    1. Load scGPT whole-human checkpoint (scPEFT/save/scGPT_human)
    2. Continue pre-training with MGEM: mask 15% of expression values, predict masked
    3. Extract domain-adapted embeddings (X_scGPT_adapted, 512 dims)

This transforms our role from "evaluator" to "builder" of foundation models.

Requirements:
    source activate /bigdata/godziklab/shared/Xinru/302006/scpeft_env
    export PYTHONPATH=/bigdata/godziklab/shared/Xinru/302006/03_MODELS/scPEFT:$PYTHONPATH

Usage:
    python scgpt_continue_pretrain.py --epochs 50
    python scgpt_continue_pretrain.py --epochs 50 --resume

Author: Xinru Qiu
Date: February 2026
"""

import os
import sys
import argparse
import logging
from pathlib import Path
import json
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
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
    RESULTS_DIR = BASE_DIR / "05_RESULTS/17_domain_adapted_fm/scgpt_adapted"
    CHECKPOINT_DIR = BASE_DIR / "05_RESULTS/checkpoints/17_scgpt_adapted"

    # scGPT model path
    _SCGPT_CANDIDATES = [
        BASE_DIR / "03_MODELS/scPEFT/save/scGPT_human",
        BASE_DIR / "03_MODELS/scGPT/scGPT_human",
    ]
    SCGPT_MODEL = next((p for p in _SCGPT_CANDIDATES if p.exists()), _SCGPT_CANDIDATES[0])

    # Pre-training hyperparameters
    BATCH_SIZE = 32
    LEARNING_RATE = 5e-5
    EPOCHS = 50
    MAX_SEQ_LEN = 1200
    MASK_RATIO = 0.15  # Mask 15% of gene expressions
    WARMUP_EPOCHS = 5

    SEED = 42


def load_platelet_data(config: Config):
    """Load platelet COVID data for continued pre-training."""
    logger.info("Loading platelet data for continued pre-training...")

    possible_files = [
        config.EMBEDDING_DIR / "state/covid_state.h5ad",
        config.DATA_DIR / "covid_pbmc_adata.h5ad",
        config.DATA_DIR / "raw/human_platelet_covid_severity_raw_counts.h5ad",
    ]

    adata = None
    for f in possible_files:
        if f.exists():
            logger.info(f"Loading from {f}")
            adata = sc.read_h5ad(f)
            break

    if adata is None:
        raise FileNotFoundError(f"Could not find platelet data in: {possible_files}")

    logger.info(f"Loaded: {adata.shape} ({adata.n_obs} cells, {adata.n_vars} genes)")
    return adata


class MGEMDataset(Dataset):
    """
    Masked Gene Expression Modeling Dataset for scGPT continued pre-training.

    For each cell:
    1. Sort genes by expression (top MAX_SEQ_LEN)
    2. Map gene names to scGPT vocabulary IDs
    3. Create masked version (15% random mask)
    """

    def __init__(self, adata, vocab, max_seq_len=1200, mask_ratio=0.15):
        self.adata = adata
        self.vocab = vocab
        self.max_seq_len = max_seq_len
        self.mask_ratio = mask_ratio
        self.gene_names = list(adata.var_names)
        self.n_cells = adata.n_obs

    def __len__(self):
        return self.n_cells

    def __getitem__(self, idx):
        # Get expression values
        if issparse(self.adata.X):
            expr = np.array(self.adata.X[idx].todense()).flatten()
        else:
            expr = self.adata.X[idx].flatten()

        # Get non-zero genes sorted by expression (descending)
        nonzero_idx = np.where(expr > 0)[0]
        if len(nonzero_idx) == 0:
            nonzero_idx = np.arange(min(10, len(expr)))

        sorted_idx = nonzero_idx[np.argsort(-expr[nonzero_idx])]
        sorted_idx = sorted_idx[:self.max_seq_len]

        # Map to vocabulary IDs
        gene_ids = []
        values = []
        for i in sorted_idx:
            gene = self.gene_names[i].upper()
            if gene in self.vocab:
                gene_ids.append(self.vocab[gene])
                values.append(float(expr[i]))

        n_genes = len(gene_ids)
        if n_genes == 0:
            # Fallback: use first 10 genes
            for i in range(min(10, len(self.gene_names))):
                gene = self.gene_names[i].upper()
                if gene in self.vocab:
                    gene_ids.append(self.vocab[gene])
                    values.append(0.0)
            n_genes = len(gene_ids)

        # Create mask for MGEM (mask 15% of positions)
        n_mask = max(1, int(n_genes * self.mask_ratio))
        mask_positions = np.random.choice(n_genes, size=n_mask, replace=False)
        masked_values = list(values)
        mask_labels = [-100.0] * n_genes  # -100 = ignore in loss

        for pos in mask_positions:
            mask_labels[pos] = values[pos]  # True value for loss
            # 80% mask, 10% random, 10% keep
            r = np.random.random()
            if r < 0.8:
                masked_values[pos] = 0.0  # Mask token value
            elif r < 0.9:
                masked_values[pos] = float(np.random.uniform(0, max(values) + 1))

        # Pad to max length
        if n_genes < self.max_seq_len:
            pad_len = self.max_seq_len - n_genes
            gene_ids.extend([0] * pad_len)
            masked_values.extend([0.0] * pad_len)
            mask_labels.extend([-100.0] * pad_len)

        return {
            'gene_ids': torch.tensor(gene_ids[:self.max_seq_len], dtype=torch.long),
            'masked_values': torch.tensor(masked_values[:self.max_seq_len], dtype=torch.float32),
            'true_values': torch.tensor(values[:self.max_seq_len] + [0.0] * max(0, self.max_seq_len - n_genes), dtype=torch.float32),
            'mask_labels': torch.tensor(mask_labels[:self.max_seq_len], dtype=torch.float32),
            'length': torch.tensor(n_genes, dtype=torch.long)
        }


def run_continue_pretrain(config: Config, resume: bool = False):
    """Continue pre-training scGPT on platelet data."""

    # Import scGPT
    try:
        from scgpt.model import TransformerModel
        from scgpt.tokenizer.gene_tokenizer import GeneVocab
        from scgpt.utils import set_seed
    except ImportError:
        logger.error("scGPT not installed. Activate scpeft_env.")
        raise

    logger.info("=" * 60)
    logger.info("ACC 17A: scGPT Continued Pre-Training on Platelet Data")
    logger.info("=" * 60)

    set_seed(config.SEED)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logger.info(f"Device: {device}")
    if device.type == 'cuda':
        logger.info(f"GPU: {torch.cuda.get_device_name(0)}")

    config.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    config.CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)

    # Load vocabulary
    vocab_file = config.SCGPT_MODEL / "vocab.json"
    if vocab_file.exists():
        vocab = GeneVocab.from_file(vocab_file)
        logger.info(f"Vocabulary: {len(vocab)} genes")
    else:
        raise FileNotFoundError(f"Vocabulary not found: {vocab_file}")

    # Load platelet data
    adata = load_platelet_data(config)
    logger.info(f"Pre-training on {adata.n_obs} platelet cells")

    # Create dataset
    dataset = MGEMDataset(adata, vocab, config.MAX_SEQ_LEN, config.MASK_RATIO)
    dataloader = DataLoader(
        dataset, batch_size=config.BATCH_SIZE,
        shuffle=True, num_workers=4, pin_memory=True, drop_last=True
    )

    # Initialize scGPT model
    model = TransformerModel(
        ntoken=len(vocab),
        d_model=512,
        nhead=8,
        d_hid=512,
        nlayers=12,
        vocab=vocab,
        dropout=0.1,
        pad_token="<pad>",
        pad_value=-2,
        n_input_bins=51,
        explicit_zero_prob=False,
        use_fast_transformer=False,
        peft_config={},
    )

    # Load pretrained weights
    model_file = config.SCGPT_MODEL / "best_model.pt"
    if model_file.exists():
        pretrained = torch.load(model_file, map_location='cpu')
        model_dict = model.state_dict()
        pretrained_dict = {k: v for k, v in pretrained.items()
                         if k in model_dict and v.shape == model_dict[k].shape}
        model_dict.update(pretrained_dict)
        model.load_state_dict(model_dict, strict=False)
        logger.info(f"Loaded {len(pretrained_dict)}/{len(model_dict)} pretrained weights")
    else:
        logger.warning(f"No pretrained weights at {model_file}, training from checkpoint state")

    model = model.to(device)

    # All parameters trainable for continued pre-training
    trainable = sum(p.numel() for p in model.parameters() if p.requires_grad)
    total = sum(p.numel() for p in model.parameters())
    logger.info(f"Parameters: {trainable:,} trainable / {total:,} total")

    # Optimizer with warmup + cosine decay
    optimizer = torch.optim.AdamW(model.parameters(), lr=config.LEARNING_RATE, weight_decay=0.01)

    # Cosine schedule with warmup
    total_steps = config.EPOCHS * len(dataloader)
    warmup_steps = config.WARMUP_EPOCHS * len(dataloader)

    def lr_lambda(step):
        if step < warmup_steps:
            return step / max(1, warmup_steps)
        progress = (step - warmup_steps) / max(1, total_steps - warmup_steps)
        return 0.5 * (1 + np.cos(np.pi * progress))

    scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, lr_lambda)

    # Resume from checkpoint
    start_epoch = 0
    best_loss = float('inf')
    if resume:
        ckpt_path = config.CHECKPOINT_DIR / "scgpt_adapted_latest.pt"
        if ckpt_path.exists():
            ckpt = torch.load(ckpt_path, map_location=device)
            model.load_state_dict(ckpt['model_state_dict'])
            optimizer.load_state_dict(ckpt['optimizer_state_dict'])
            scheduler.load_state_dict(ckpt['scheduler_state_dict'])
            start_epoch = ckpt['epoch']
            best_loss = ckpt.get('best_loss', float('inf'))
            logger.info(f"Resumed from epoch {start_epoch}, best loss {best_loss:.4f}")

    # Training loop
    scaler = torch.amp.GradScaler('cuda') if device.type == 'cuda' else None
    loss_history = []

    for epoch in range(start_epoch, config.EPOCHS):
        model.train()
        epoch_loss = 0
        n_batches = 0

        for batch in dataloader:
            gene_ids = batch['gene_ids'].to(device)
            masked_values = batch['masked_values'].to(device)
            mask_labels = batch['mask_labels'].to(device)
            src_key_padding_mask = gene_ids == 0

            optimizer.zero_grad()

            if scaler is not None:
                with torch.amp.autocast('cuda'):
                    output = model(
                        gene_ids,
                        masked_values,
                        src_key_padding_mask=src_key_padding_mask,
                        CLS=False,
                        CCE=False,
                        MVC=True,   # Masked value prediction
                        ECS=False,
                    )

                    # MGEM loss: MSE on masked positions only
                    if 'mvc_output' in output:
                        pred_values = output['mvc_output']
                        # Only compute loss on masked positions
                        mask = mask_labels != -100.0
                        if mask.any():
                            loss = F.mse_loss(pred_values[mask], mask_labels[mask])
                        else:
                            loss = torch.tensor(0.0, device=device)
                    else:
                        # Fallback: use cell embedding reconstruction
                        loss = torch.tensor(0.0, device=device)
                        for key in ['mlm_output', 'cell_emb']:
                            if key in output:
                                loss = F.mse_loss(
                                    output[key][:, :gene_ids.shape[1]].mean(dim=-1) if output[key].dim() == 3 else output[key],
                                    masked_values
                                )
                                break

                scaler.scale(loss).backward()
                scaler.unscale_(optimizer)
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                scaler.step(optimizer)
                scaler.update()
            else:
                output = model(
                    gene_ids,
                    masked_values,
                    src_key_padding_mask=src_key_padding_mask,
                    CLS=False, CCE=False, MVC=True, ECS=False,
                )

                if 'mvc_output' in output:
                    pred_values = output['mvc_output']
                    mask = mask_labels != -100.0
                    if mask.any():
                        loss = F.mse_loss(pred_values[mask], mask_labels[mask])
                    else:
                        loss = torch.tensor(0.0, device=device)
                else:
                    loss = torch.tensor(0.0, device=device)

                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                optimizer.step()

            scheduler.step()
            epoch_loss += loss.item()
            n_batches += 1

        avg_loss = epoch_loss / max(n_batches, 1)
        loss_history.append(avg_loss)
        lr = optimizer.param_groups[0]['lr']

        logger.info(f"Epoch {epoch+1}/{config.EPOCHS} | Loss: {avg_loss:.4f} | LR: {lr:.2e}")

        # Save checkpoint
        if avg_loss < best_loss:
            best_loss = avg_loss
            torch.save({
                'epoch': epoch + 1,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'scheduler_state_dict': scheduler.state_dict(),
                'best_loss': best_loss,
            }, config.CHECKPOINT_DIR / "scgpt_adapted_best.pt")
            logger.info(f"  Saved best checkpoint (loss={best_loss:.4f})")

        # Save latest for resume
        if (epoch + 1) % 10 == 0 or epoch == config.EPOCHS - 1:
            torch.save({
                'epoch': epoch + 1,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'scheduler_state_dict': scheduler.state_dict(),
                'best_loss': best_loss,
            }, config.CHECKPOINT_DIR / "scgpt_adapted_latest.pt")

    # Save training history
    history = {
        'loss_history': loss_history,
        'best_loss': best_loss,
        'epochs': config.EPOCHS,
        'n_cells': adata.n_obs,
        'mask_ratio': config.MASK_RATIO,
        'learning_rate': config.LEARNING_RATE,
    }
    with open(config.RESULTS_DIR / "pretrain_history.json", 'w') as f:
        json.dump(history, f, indent=2)

    logger.info(f"\nContinued pre-training complete. Best loss: {best_loss:.4f}")
    logger.info(f"Checkpoint: {config.CHECKPOINT_DIR / 'scgpt_adapted_best.pt'}")

    return model


def main():
    parser = argparse.ArgumentParser(description='scGPT Continued Pre-Training')
    parser.add_argument('--epochs', type=int, default=50)
    parser.add_argument('--batch-size', type=int, default=32)
    parser.add_argument('--lr', type=float, default=5e-5)
    parser.add_argument('--mask-ratio', type=float, default=0.15)
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()

    config = Config()
    config.EPOCHS = args.epochs
    config.BATCH_SIZE = args.batch_size
    config.LEARNING_RATE = args.lr
    config.MASK_RATIO = args.mask_ratio

    torch.manual_seed(config.SEED)
    np.random.seed(config.SEED)

    run_continue_pretrain(config, resume=args.resume)


if __name__ == '__main__':
    main()
