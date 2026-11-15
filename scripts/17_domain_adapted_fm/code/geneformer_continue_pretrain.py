#!/usr/bin/env python3
"""
ACC 17A: Continue Pre-Training Geneformer on Platelet Data

Domain-adapt the Geneformer-12L-30M checkpoint to platelet transcriptomics
using masked language modeling (MLM) on rank-value encoded gene tokens.

Approach:
    1. Load Geneformer-12L-30M (BertForMaskedLM, hidden_size=256)
    2. Tokenize platelet cells using rank-value encoding (Ensembl IDs)
    3. Continue pre-training with MLM: 15% masking, BERT 80/10/10 strategy
    4. Save adapted model in HuggingFace format for downstream embedding extraction

This transforms our role from "evaluator" to "builder" of foundation models.

Requirements:
    source activate /bigdata/godziklab/shared/Xinru/302006/scpeft_env
    export PYTHONPATH=/bigdata/godziklab/shared/Xinru/302006/03_MODELS/scPEFT:$PYTHONPATH

Usage:
    python geneformer_continue_pretrain.py --epochs 50
    python geneformer_continue_pretrain.py --epochs 50 --resume

Author: Xinru Qiu
Date: March 2026
"""

import os
import sys
import argparse
import logging
import pickle
from pathlib import Path
import json
import numpy as np
import torch
import torch.nn as nn
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
    RESULTS_DIR = BASE_DIR / "05_RESULTS/17_domain_adapted_fm/geneformer_adapted"
    CHECKPOINT_DIR = BASE_DIR / "05_RESULTS/checkpoints/17_geneformer_adapted"

    # Geneformer model path
    GENEFORMER_MODEL = BASE_DIR / "03_MODELS/geneformer"
    _PRETRAIN_CANDIDATES = [
        GENEFORMER_MODEL / "geneformer-V2-316M-native",
        GENEFORMER_MODEL / "Geneformer",
        BASE_DIR / "03_MODELS/Geneformer",
    ]
    PRETRAIN_CKPT = next((p for p in _PRETRAIN_CANDIDATES if p.exists()), _PRETRAIN_CANDIDATES[0])

    # Geneformer architecture (12L-30M)
    VOCAB_SIZE = 20275  # Matches geneformer-V2-316M-native config.json
    HIDDEN_SIZE = 256
    NUM_LAYERS = 12
    NUM_HEADS = 4
    INTERMEDIATE_SIZE = 1024
    MAX_POSITION_EMBEDDINGS = 4096  # Matches geneformer-V2-316M-native config.json

    # Pre-training hyperparameters
    BATCH_SIZE = 8
    LEARNING_RATE = 1e-5
    EPOCHS = 50
    MAX_SEQ_LEN = 512  # Reduced from 2048 for full pre-training memory (316M params)
    MASK_RATIO = 0.15  # Mask 15% of gene tokens
    WARMUP_EPOCHS = 5

    # Special tokens
    PAD_TOKEN_ID = 0
    MASK_TOKEN_ID = 1
    CLS_TOKEN_ID = 2
    SEP_TOKEN_ID = 3

    SEED = 42


def load_geneformer_token_dict(config: Config):
    """Load Geneformer gene-to-token dictionary."""
    candidates = [
        config.BASE_DIR / "03_MODELS/geneformer/geneformer-V2-316M-native/geneformer/token_dictionary_gc104M.pkl",
        config.BASE_DIR / "03_MODELS/geneformer/Geneformer/geneformer/token_dictionary_gc104M.pkl",
        config.BASE_DIR / "03_MODELS/scPEFT/geneformer_peft/geneformer/token_dictionary.pkl",
        config.GENEFORMER_MODEL / "geneformer" / "token_dictionary_gc104M.pkl",
        config.PRETRAIN_CKPT / "token_dictionary.pkl",
    ]
    for token_dict_file in candidates:
        if token_dict_file.exists():
            with open(token_dict_file, 'rb') as f:
                token_dict = pickle.load(f)
            logger.info(f"Loaded token dictionary: {len(token_dict)} genes from {token_dict_file}")
            return token_dict

    raise FileNotFoundError(
        f"Token dictionary not found at any expected path: {candidates}"
    )


def load_platelet_data(config: Config):
    """Load platelet COVID data for continued pre-training.

    Priority:
      1. Geneformer-ready file (var_names already Ensembl IDs)
      2. Raw file with ensembl_id column in .var -> remap var_names
      3. Generic preprocessed data (will warn about tokenization)
    """
    logger.info("Loading platelet data for continued pre-training...")

    # Try Geneformer-ready data first
    geneformer_ready = config.DATA_DIR / "geneformer_ready/covid_geneformer_ready.h5ad"
    if geneformer_ready.exists():
        adata = sc.read_h5ad(geneformer_ready)
        logger.info(f"Loaded Geneformer-ready data: {adata.shape} from {geneformer_ready}")
        return adata

    # Fall back to other data sources
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
        raise FileNotFoundError(f"Could not find platelet data in: {[geneformer_ready] + possible_files}")

    # Remap var_names to Ensembl IDs if needed
    first_var = str(adata.var_names[0])
    if not first_var.startswith("ENSG"):
        if 'ensembl_id' in adata.var.columns:
            import pandas as pd
            logger.info("Remapping var_names from gene symbols to Ensembl IDs via 'ensembl_id' column")
            if 'gene_symbol' not in adata.var.columns:
                adata.var['gene_symbol'] = adata.var_names.values
            adata.var_names = pd.Index(adata.var['ensembl_id'].values)
        else:
            logger.warning(
                "var_names are gene symbols and no 'ensembl_id' column found. "
                "Geneformer token dictionary uses Ensembl IDs - tokenization may have low match rate. "
                "Run prepare_data_for_geneformer.py first."
            )

    logger.info(f"Loaded: {adata.shape} ({adata.n_obs} cells, {adata.n_vars} genes)")
    return adata


# =============================================================================
# MLM Dataset
# =============================================================================

class GeneformerMLMDataset(Dataset):
    """
    Masked Language Modeling Dataset for Geneformer continued pre-training.

    For each cell:
    1. Sort genes by expression (rank-value encoding, Geneformer convention)
    2. Map Ensembl IDs to Geneformer token IDs
    3. Apply MLM masking: 15% with BERT 80/10/10 strategy
    """

    def __init__(self, adata, token_dict, max_len=2048, mask_ratio=0.15,
                 mask_token_id=1, cls_token_id=2, sep_token_id=3, pad_token_id=0,
                 vocab_size=25500):
        self.adata = adata
        self.token_dict = token_dict
        self.max_len = max_len
        self.mask_ratio = mask_ratio
        self.mask_token_id = mask_token_id
        self.cls_token_id = cls_token_id
        self.sep_token_id = sep_token_id
        self.pad_token_id = pad_token_id
        self.vocab_size = vocab_size
        self.gene_names = list(adata.var_names)
        self.n_cells = adata.n_obs

        # Pre-compute per-gene token IDs
        self.gene_token_ids = []
        matched = 0
        for name in self.gene_names:
            tid = token_dict.get(name)
            if tid is None:
                tid = token_dict.get(name.upper())
            self.gene_token_ids.append(tid)
            if tid is not None:
                matched += 1

        match_rate = matched / len(self.gene_names) if self.gene_names else 0
        logger.info(
            f"GeneformerMLMDataset: {matched}/{len(self.gene_names)} genes "
            f"matched token dictionary ({match_rate:.1%})"
        )
        if match_rate < 0.10:
            raise ValueError(
                f"Token match rate is only {match_rate:.1%} - almost no genes match. "
                f"var_names are likely gene symbols instead of Ensembl IDs. "
                f"Run prepare_data_for_geneformer.py first."
            )

    def __len__(self):
        return self.n_cells

    def __getitem__(self, idx):
        # Get expression values
        if issparse(self.adata.X):
            expr = np.array(self.adata.X[idx].todense()).flatten()
        else:
            expr = self.adata.X[idx].flatten()

        # Get non-zero genes sorted by expression (rank ordering)
        nonzero_idx = np.where(expr > 0)[0]
        if len(nonzero_idx) == 0:
            nonzero_idx = np.arange(min(10, len(expr)))
        sorted_idx = nonzero_idx[np.argsort(-expr[nonzero_idx])]

        # Convert to token IDs (skip genes not in dictionary)
        gene_tokens = []
        for i in sorted_idx:
            if len(gene_tokens) >= self.max_len - 2:  # Leave room for CLS + SEP
                break
            tid = self.gene_token_ids[i]
            if tid is not None:
                gene_tokens.append(tid)

        if len(gene_tokens) == 0:
            # Fallback: use first available gene tokens
            for tid in self.gene_token_ids:
                if tid is not None:
                    gene_tokens.append(tid)
                    if len(gene_tokens) >= 10:
                        break

        # Build input: [CLS] gene1 gene2 ... geneN [SEP]
        input_ids = [self.cls_token_id] + gene_tokens + [self.sep_token_id]
        seq_len = len(input_ids)

        # MLM masking: only mask gene tokens (not CLS/SEP)
        labels = [-100] * seq_len  # -100 = ignore in cross-entropy loss
        n_gene_tokens = len(gene_tokens)
        n_mask = max(1, int(n_gene_tokens * self.mask_ratio))
        mask_positions = np.random.choice(n_gene_tokens, size=n_mask, replace=False) + 1  # +1 for CLS offset

        for pos in mask_positions:
            labels[pos] = input_ids[pos]  # True token for loss computation
            r = np.random.random()
            if r < 0.8:
                input_ids[pos] = self.mask_token_id  # 80% -> [MASK]
            elif r < 0.9:
                input_ids[pos] = np.random.randint(4, self.vocab_size)  # 10% -> random token
            # else: 10% -> keep original

        # Attention mask
        attention_mask = [1] * seq_len

        # Pad to max_len
        if seq_len < self.max_len:
            pad_len = self.max_len - seq_len
            input_ids.extend([self.pad_token_id] * pad_len)
            attention_mask.extend([0] * pad_len)
            labels.extend([-100] * pad_len)

        return {
            'input_ids': torch.tensor(input_ids[:self.max_len], dtype=torch.long),
            'attention_mask': torch.tensor(attention_mask[:self.max_len], dtype=torch.long),
            'labels': torch.tensor(labels[:self.max_len], dtype=torch.long),
        }


# =============================================================================
# Training
# =============================================================================

def run_continue_pretrain(config: Config, resume: bool = False):
    """Continue pre-training Geneformer on platelet data with MLM."""

    from transformers import BertForMaskedLM, BertConfig

    logger.info("=" * 60)
    logger.info("ACC 17A: Geneformer Continued Pre-Training on Platelet Data")
    logger.info("=" * 60)

    torch.manual_seed(config.SEED)
    np.random.seed(config.SEED)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logger.info(f"Device: {device}")
    if device.type == 'cuda':
        logger.info(f"GPU: {torch.cuda.get_device_name(0)}")

    config.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    config.CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)

    # Load token dictionary
    token_dict = load_geneformer_token_dict(config)

    # Load platelet data
    adata = load_platelet_data(config)
    logger.info(f"Pre-training on {adata.n_obs} platelet cells")

    # Create MLM dataset
    dataset = GeneformerMLMDataset(
        adata, token_dict,
        max_len=config.MAX_SEQ_LEN,
        mask_ratio=config.MASK_RATIO,
        mask_token_id=config.MASK_TOKEN_ID,
        cls_token_id=config.CLS_TOKEN_ID,
        sep_token_id=config.SEP_TOKEN_ID,
        pad_token_id=config.PAD_TOKEN_ID,
        vocab_size=config.VOCAB_SIZE,
    )
    dataloader = DataLoader(
        dataset, batch_size=config.BATCH_SIZE,
        shuffle=True, num_workers=4, pin_memory=True, drop_last=True,
    )

    # Load pretrained Geneformer
    model_path = config.PRETRAIN_CKPT
    if model_path.exists():
        logger.info(f"Loading pretrained Geneformer from {model_path}")
        model = BertForMaskedLM.from_pretrained(
            model_path,
            output_attentions=False,
            output_hidden_states=False,
        )
    else:
        logger.warning(f"Pretrained checkpoint not found at {model_path}, creating new model")
        bert_config = BertConfig(
            vocab_size=config.VOCAB_SIZE,
            hidden_size=config.HIDDEN_SIZE,
            num_hidden_layers=config.NUM_LAYERS,
            num_attention_heads=config.NUM_HEADS,
            intermediate_size=config.INTERMEDIATE_SIZE,
            hidden_dropout_prob=0.02,
            attention_probs_dropout_prob=0.02,
            max_position_embeddings=config.MAX_POSITION_EMBEDDINGS,
        )
        model = BertForMaskedLM(bert_config)

    # All parameters trainable for continued pre-training (full fine-tuning, not LoRA)
    for param in model.parameters():
        param.requires_grad = True

    trainable = sum(p.numel() for p in model.parameters() if p.requires_grad)
    total = sum(p.numel() for p in model.parameters())
    logger.info(f"Parameters: {trainable:,} trainable / {total:,} total (full continued pre-training)")

    model = model.to(device)

    # Optimizer with warmup + cosine decay
    optimizer = torch.optim.AdamW(model.parameters(), lr=config.LEARNING_RATE, weight_decay=0.01)

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
    global_step = 0
    if resume:
        ckpt_path = config.CHECKPOINT_DIR / "geneformer_adapted_latest.pt"
        if ckpt_path.exists():
            ckpt = torch.load(ckpt_path, map_location=device)
            model.load_state_dict(ckpt['model_state_dict'])
            optimizer.load_state_dict(ckpt['optimizer_state_dict'])
            scheduler.load_state_dict(ckpt['scheduler_state_dict'])
            start_epoch = ckpt['epoch']
            best_loss = ckpt.get('best_loss', float('inf'))
            global_step = ckpt.get('global_step', start_epoch * len(dataloader))
            logger.info(f"Resumed from epoch {start_epoch}, best loss {best_loss:.4f}")
        else:
            logger.warning(f"No checkpoint found at {ckpt_path}, starting from scratch")

    # AMP scaler (compatible with PyTorch 2.1+)
    scaler = torch.amp.GradScaler('cuda') if device.type == 'cuda' else None
    loss_history = []

    logger.info(f"Starting training: {config.EPOCHS} epochs, {len(dataloader)} batches/epoch")
    logger.info(f"Warmup: {config.WARMUP_EPOCHS} epochs ({warmup_steps} steps)")

    for epoch in range(start_epoch, config.EPOCHS):
        model.train()
        epoch_loss = 0
        n_batches = 0

        for batch in dataloader:
            input_ids = batch['input_ids'].to(device)
            attention_mask = batch['attention_mask'].to(device)
            labels = batch['labels'].to(device)

            optimizer.zero_grad()

            if scaler is not None:
                with torch.amp.autocast('cuda'):
                    outputs = model(
                        input_ids=input_ids,
                        attention_mask=attention_mask,
                        labels=labels,
                    )
                    loss = outputs.loss

                scaler.scale(loss).backward()
                scaler.unscale_(optimizer)
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                scaler.step(optimizer)
                scaler.update()
            else:
                outputs = model(
                    input_ids=input_ids,
                    attention_mask=attention_mask,
                    labels=labels,
                )
                loss = outputs.loss

                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                optimizer.step()

            scheduler.step()
            epoch_loss += loss.item()
            n_batches += 1
            global_step += 1

        avg_loss = epoch_loss / max(n_batches, 1)
        loss_history.append(avg_loss)
        lr = optimizer.param_groups[0]['lr']

        logger.info(f"Epoch {epoch+1}/{config.EPOCHS} | MLM Loss: {avg_loss:.4f} | LR: {lr:.2e}")

        # Save best checkpoint (HuggingFace format + .pt for resume)
        if avg_loss < best_loss:
            best_loss = avg_loss
            # Save in HuggingFace format for easy loading with from_pretrained()
            best_hf_dir = config.CHECKPOINT_DIR / "best_model"
            model.save_pretrained(best_hf_dir)
            # Also save optimizer state for resume
            torch.save({
                'epoch': epoch + 1,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'scheduler_state_dict': scheduler.state_dict(),
                'best_loss': best_loss,
                'global_step': global_step,
            }, config.CHECKPOINT_DIR / "geneformer_adapted_best.pt")
            logger.info(f"  Saved best checkpoint (loss={best_loss:.4f})")

        # Save latest for resume
        if (epoch + 1) % 10 == 0 or epoch == config.EPOCHS - 1:
            torch.save({
                'epoch': epoch + 1,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'scheduler_state_dict': scheduler.state_dict(),
                'best_loss': best_loss,
                'global_step': global_step,
            }, config.CHECKPOINT_DIR / "geneformer_adapted_latest.pt")

    # Save final model in HuggingFace format
    final_hf_dir = config.RESULTS_DIR / "model"
    model.save_pretrained(final_hf_dir)
    logger.info(f"Saved final adapted model to {final_hf_dir}")

    # Save training history
    history = {
        'loss_history': loss_history,
        'best_loss': best_loss,
        'epochs': config.EPOCHS,
        'n_cells': adata.n_obs,
        'mask_ratio': config.MASK_RATIO,
        'learning_rate': config.LEARNING_RATE,
        'batch_size': config.BATCH_SIZE,
        'max_seq_len': config.MAX_SEQ_LEN,
        'model': 'Geneformer-12L-30M',
    }
    with open(config.RESULTS_DIR / "pretrain_history.json", 'w') as f:
        json.dump(history, f, indent=2)

    logger.info(f"\nContinued pre-training complete. Best MLM loss: {best_loss:.4f}")
    logger.info(f"Best checkpoint: {config.CHECKPOINT_DIR / 'best_model'}")
    logger.info(f"Final model: {final_hf_dir}")

    return model


def main():
    parser = argparse.ArgumentParser(description='Geneformer Continued Pre-Training')
    parser.add_argument('--epochs', type=int, default=50)
    parser.add_argument('--batch-size', type=int, default=16)
    parser.add_argument('--lr', type=float, default=1e-5)
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
