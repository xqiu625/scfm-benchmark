#!/usr/bin/env python3
"""
STATE SE-600M LoRA Fine-Tuning for COVID-19 Severity Classification

This script applies Parameter-Efficient Fine-Tuning (PEFT) using LoRA (Low-Rank Adaptation)
to the STATE SE-600M model. Instead of training all 603M parameters, we only train
~1-3% of parameters through low-rank adapter matrices.

Benefits over full fine-tuning:
    - 96-99% fewer trainable parameters (~5-20M vs 603M)
    - ~50% less GPU memory
    - Faster training
    - Reduces catastrophic forgetting of pre-trained knowledge
    - Can use more genes per cell (larger sequence length)

Architecture:
    [Raw Expression + Gene Names]
    → [ESM2 Gene Embeddings (5120 dim)]
    → [SE Encoder (frozen)]
    → [SE Transformer with LoRA adapters (16 layers)]
    → [Cell Embedding]
    → [Classification Head]
    → [Predictions]

Based on: scPEFT (Nature Machine Intelligence, 2025)
    - Paper: https://doi.org/10.1038/s42256-025-01170-z
    - Code: https://github.com/coffee19850519/scPEFT

Requirements:
    pip install peft arc-state torch scanpy

Usage:
    python state_lora_finetune.py --task binary --epochs 15
    python state_lora_finetune.py --task binary --epochs 15 --lora-r 8 --max-genes 512

Author: Xinru Qiu
Date: January 2026
"""

import sys
import os
import argparse
import logging
from pathlib import Path
import json
from typing import Dict, List, Optional

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from torch.optim import AdamW
from torch.optim.lr_scheduler import CosineAnnealingWarmRestarts
import scanpy as sc
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, accuracy_score, f1_score

# =============================================================================
# CRITICAL: Patch torch.load BEFORE any STATE imports
# =============================================================================
_original_torch_load = torch.load

def _patched_torch_load(f, *args, **kwargs):
    kwargs['weights_only'] = False
    return _original_torch_load(f, *args, **kwargs)

torch.load = _patched_torch_load

# Setup VCI aliases for backward compatibility
from state.emb import _setup_vci_aliases
_setup_vci_aliases()

import state.emb as emb
sys.modules['vci'] = emb
sys.modules['vci.nn'] = emb.nn
sys.modules['vci.nn.model'] = emb.nn.model

from state.emb.nn.model import StateEmbeddingModel

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


# =============================================================================
# Configuration
# =============================================================================

class Config:
    BASE_DIR = Path("/bigdata/godziklab/shared/Xinru/302006")
    MODEL_DIR = BASE_DIR / "03_MODELS/SE-600M"
    EMBEDDING_DIR = BASE_DIR / "02_EMBEDDINGS/state"
    RESULTS_DIR = BASE_DIR / "05_RESULTS/state_lora_finetune"
    CHECKPOINT_DIR = BASE_DIR / "05_RESULTS/checkpoints/state_lora_finetune"

    SE_CHECKPOINT = MODEL_DIR / "se600m_epoch16.ckpt"
    PROTEIN_EMBEDDINGS = MODEL_DIR / "protein_embeddings.pt"

    # Training - can use larger batch/genes due to LoRA efficiency
    BATCH_SIZE = 16
    LEARNING_RATE = 1e-4  # Higher LR ok for LoRA
    WEIGHT_DECAY = 0.01
    EPOCHS = 15
    MAX_GENES = 512  # Can use more genes with LoRA (less memory)

    # LoRA configuration
    LORA_R = 16          # Low-rank dimension (smaller = fewer params)
    LORA_ALPHA = 32      # Scaling factor (typically 2x r)
    LORA_DROPOUT = 0.1

    # Classifier
    HIDDEN_DIM = 512
    DROPOUT = 0.1

    N_FOLDS = 5
    SEED = 42

    TASK_CONFIGS = {
        'binary': {'n_classes': 2, 'label_map': {
            'healthy': 0, 'mild': 0, 'moderate': 0, 'recovered': 0,
            'severe': 1, 'fatal': 1
        }},
        '3class': {'n_classes': 3, 'label_map': {
            'healthy': 0, 'mild': 1, 'moderate': 1,
            'severe': 2, 'fatal': 2, 'recovered': 1
        }},
        '6class': {'n_classes': 6, 'label_map': {
            'healthy': 0, 'mild': 1, 'moderate': 2,
            'severe': 3, 'fatal': 4, 'recovered': 5
        }}
    }


# =============================================================================
# LoRA Implementation for STATE Transformer
# =============================================================================

class LoRALayer(nn.Module):
    """
    Low-Rank Adaptation layer.

    Implements: W' = W + BA where B is (d, r) and A is (r, d)
    Only A and B are trainable, W is frozen.
    """

    def __init__(self, original_layer: nn.Linear, r: int = 16, alpha: int = 32, dropout: float = 0.1):
        super().__init__()
        self.original_layer = original_layer
        self.r = r
        self.alpha = alpha
        self.scaling = alpha / r

        in_features = original_layer.in_features
        out_features = original_layer.out_features

        # Freeze original weights
        for param in self.original_layer.parameters():
            param.requires_grad = False

        # LoRA matrices
        self.lora_A = nn.Parameter(torch.zeros(r, in_features))
        self.lora_B = nn.Parameter(torch.zeros(out_features, r))
        self.dropout = nn.Dropout(dropout)

        # Initialize A with Kaiming, B with zeros (so initial output = original)
        nn.init.kaiming_uniform_(self.lora_A, a=np.sqrt(5))
        nn.init.zeros_(self.lora_B)

    def forward(self, x):
        # Original forward
        original_output = self.original_layer(x)

        # LoRA forward: x @ A^T @ B^T * scaling
        lora_output = self.dropout(x) @ self.lora_A.T @ self.lora_B.T * self.scaling

        return original_output + lora_output


def apply_lora_to_transformer(transformer_encoder, r: int = 16, alpha: int = 32, dropout: float = 0.1):
    """
    Apply LoRA adapters to transformer encoder layers.

    Targets the linear layers in the feed-forward network (linear1, linear2)
    and any attention projections if accessible.
    """
    lora_layers_added = 0

    for name, module in transformer_encoder.named_modules():
        # Target FFN layers (linear1, linear2 in each transformer layer)
        if isinstance(module, nn.Linear) and ('linear1' in name or 'linear2' in name):
            # Get parent module and attribute name
            parts = name.rsplit('.', 1)
            if len(parts) == 2:
                parent_name, attr_name = parts
                parent = transformer_encoder
                for part in parent_name.split('.'):
                    parent = getattr(parent, part)

                # Replace with LoRA layer
                original_layer = getattr(parent, attr_name)
                lora_layer = LoRALayer(original_layer, r=r, alpha=alpha, dropout=dropout)
                setattr(parent, attr_name, lora_layer)
                lora_layers_added += 1

    logger.info(f"Applied LoRA to {lora_layers_added} layers (r={r}, alpha={alpha})")
    return lora_layers_added


def count_parameters(model):
    """Count trainable and total parameters."""
    trainable = sum(p.numel() for p in model.parameters() if p.requires_grad)
    total = sum(p.numel() for p in model.parameters())
    return trainable, total


# =============================================================================
# Data Loading (same as v3)
# =============================================================================

def load_protein_embeddings(path: Path) -> Dict[str, torch.Tensor]:
    """Load ESM2 protein embeddings (on CPU for DataLoader compatibility)."""
    logger.info(f"Loading protein embeddings from {path}")
    embeddings = torch.load(str(path), map_location='cpu', weights_only=False)

    if isinstance(embeddings, dict):
        embeddings = {k.upper(): v.cpu() for k, v in embeddings.items()}
        logger.info(f"Loaded {len(embeddings)} gene embeddings")
        sample_key = next(iter(embeddings.keys()))
        logger.info(f"Embedding dimension: {embeddings[sample_key].shape}")

    return embeddings


def load_raw_covid_data(config: Config, task: str = 'binary'):
    """Load raw COVID expression data."""
    logger.info("Loading raw COVID data...")

    possible_files = [
        config.EMBEDDING_DIR / "covid_state.h5ad",
        config.BASE_DIR / "01_DATA" / "covid" / "covid_raw.h5ad",
    ]

    adata = None
    for f in possible_files:
        if f.exists():
            logger.info(f"Loading from {f}")
            adata = sc.read_h5ad(f)
            break

    if adata is None:
        raise FileNotFoundError(f"Could not find COVID data")

    logger.info(f"Loaded AnnData: {adata.shape}")

    X = adata.raw.X if hasattr(adata, 'raw') and adata.raw is not None else adata.X
    if hasattr(X, 'toarray'):
        X = X.toarray()
    X = np.array(X, dtype=np.float32)

    var_names = list(adata.raw.var_names if hasattr(adata, 'raw') and adata.raw is not None else adata.var_names)

    severity = adata.obs.get('covid_severity', adata.obs.get('severity'))
    label_map = config.TASK_CONFIGS[task]['label_map']

    labels, valid_indices = [], []
    for i, s in enumerate(severity):
        s_lower = str(s).lower().strip()
        if s_lower in label_map:
            labels.append(label_map[s_lower])
            valid_indices.append(i)

    labels = np.array(labels)
    X = X[valid_indices]

    logger.info(f"Valid cells: {len(labels)}, Expression shape: {X.shape}")
    return X, var_names, labels


# =============================================================================
# Dataset
# =============================================================================

class RawExpressionDataset(Dataset):
    """Dataset for SE model fine-tuning with expression data."""

    def __init__(self, expression_matrix: np.ndarray, gene_names: List[str],
                 labels: np.ndarray, protein_embeddings: Dict[str, torch.Tensor],
                 max_genes: int = 512):
        self.expression = torch.FloatTensor(expression_matrix)
        self.gene_names = [g.upper() for g in gene_names]
        self.labels = torch.LongTensor(labels)
        self.protein_embeddings = protein_embeddings
        self.max_genes = max_genes

        self.gene_to_idx = {g: i for i, g in enumerate(self.gene_names)}
        self.valid_genes = [g for g in self.gene_names if g in protein_embeddings]
        self.valid_gene_indices = [self.gene_to_idx[g] for g in self.valid_genes]

        logger.info(f"Genes with ESM2 embeddings: {len(self.valid_genes)}/{len(self.gene_names)}")

        if self.valid_genes:
            self.embed_dim = protein_embeddings[self.valid_genes[0]].shape[-1]
        else:
            raise ValueError("No genes found in protein embeddings!")

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        cell_expr = self.expression[idx]
        label = self.labels[idx]

        valid_expr = cell_expr[self.valid_gene_indices]

        if len(self.valid_genes) > self.max_genes:
            top_indices = torch.topk(valid_expr, self.max_genes).indices
        else:
            top_indices = torch.argsort(valid_expr, descending=True)

        gene_embeds, expr_weights = [], []
        for i in top_indices:
            gene_name = self.valid_genes[i]
            gene_embeds.append(self.protein_embeddings[gene_name])
            expr_weights.append(valid_expr[i])

        n_genes = len(gene_embeds)
        if n_genes < self.max_genes:
            padding = self.max_genes - n_genes
            zero_emb = torch.zeros(self.embed_dim)
            gene_embeds.extend([zero_emb] * padding)
            expr_weights.extend([0.0] * padding)

        gene_embeds = torch.stack(gene_embeds[:self.max_genes])
        expr_weights = torch.tensor(expr_weights[:self.max_genes], dtype=torch.float32)
        mask = expr_weights > 0

        return {'gene_embeds': gene_embeds, 'expr_weights': expr_weights, 'mask': mask, 'label': label}


# =============================================================================
# Model with LoRA
# =============================================================================

class SELoRAClassifier(nn.Module):
    """SE-600M with LoRA adapters and classification head."""

    def __init__(self, se_model: StateEmbeddingModel, n_classes: int,
                 hidden_dim: int = 512, dropout: float = 0.1,
                 lora_r: int = 16, lora_alpha: int = 32, lora_dropout: float = 0.1):
        super().__init__()

        self.se_model = se_model
        self.embed_dim = se_model.hparams.get('d_model', 2048)

        # Freeze entire SE model first
        for param in self.se_model.parameters():
            param.requires_grad = False
        logger.info("Froze SE-600M backbone (603M parameters)")

        # Apply LoRA to transformer encoder
        n_lora = apply_lora_to_transformer(
            self.se_model.transformer_encoder,
            r=lora_r, alpha=lora_alpha, dropout=lora_dropout
        )

        # Classification head (trainable)
        self.classifier = nn.Sequential(
            nn.Linear(self.embed_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.LayerNorm(hidden_dim // 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, n_classes)
        )

        trainable, total = count_parameters(self)
        logger.info(f"Trainable parameters: {trainable/1e6:.2f}M / {total/1e6:.1f}M ({100*trainable/total:.2f}%)")

    def forward(self, gene_embeds, expr_weights, mask):
        batch_size, max_genes, embed_dim = gene_embeds.shape

        # Encoder (frozen)
        gene_embeds_flat = gene_embeds.view(batch_size * max_genes, embed_dim)
        encoded = self.se_model.encoder(gene_embeds_flat)
        encoded = encoded.view(batch_size, max_genes, -1)

        # Expression weighting
        expr_weights_norm = expr_weights / (expr_weights.sum(dim=1, keepdim=True) + 1e-8)
        encoded_weighted = encoded * expr_weights_norm.unsqueeze(-1)

        # Transformer with LoRA (LoRA layers are trainable)
        encoded_t = encoded_weighted.transpose(0, 1)
        attn_mask = ~mask
        transformer_out = self.se_model.transformer_encoder(encoded_t, src_key_padding_mask=attn_mask)
        transformer_out = transformer_out.transpose(0, 1)

        # Masked mean pooling
        mask_expanded = mask.unsqueeze(-1).float()
        pooled = (transformer_out * mask_expanded).sum(dim=1)
        pooled = pooled / (mask_expanded.sum(dim=1) + 1e-8)

        return self.classifier(pooled)


def load_se_model(config: Config, device: torch.device) -> StateEmbeddingModel:
    """Load SE-600M model."""
    logger.info(f"Loading SE-600M from {config.SE_CHECKPOINT}")
    checkpoint = torch.load(str(config.SE_CHECKPOINT), map_location=device, weights_only=False)

    hparams = checkpoint['hyper_parameters']
    logger.info(f"Model config: d_model={hparams['d_model']}, nlayers={hparams['nlayers']}")

    model = StateEmbeddingModel(**hparams)
    model.load_state_dict(checkpoint['state_dict'], strict=False)
    model = model.to(device)

    return model


# =============================================================================
# Training
# =============================================================================

def train_epoch(model, loader, optimizer, device, scaler):
    """Train one epoch with mixed precision."""
    model.train()
    total_loss, correct, total = 0, 0, 0

    for batch in loader:
        gene_embeds = batch['gene_embeds'].to(device)
        expr_weights = batch['expr_weights'].to(device)
        mask = batch['mask'].to(device)
        labels = batch['label'].to(device)

        optimizer.zero_grad()

        with torch.amp.autocast('cuda'):
            logits = model(gene_embeds, expr_weights, mask)
            loss = F.cross_entropy(logits, labels)

        scaler.scale(loss).backward()
        scaler.unscale_(optimizer)
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        scaler.step(optimizer)
        scaler.update()

        total_loss += loss.item()
        correct += (logits.argmax(dim=1) == labels).sum().item()
        total += len(labels)

    return total_loss / len(loader), correct / total


def evaluate(model, loader, device, n_classes):
    """Evaluate model."""
    model.eval()
    all_probs, all_labels = [], []

    with torch.no_grad():
        for batch in loader:
            gene_embeds = batch['gene_embeds'].to(device)
            expr_weights = batch['expr_weights'].to(device)
            mask = batch['mask'].to(device)
            labels = batch['label'].to(device)

            with torch.amp.autocast('cuda'):
                logits = model(gene_embeds, expr_weights, mask)
                probs = F.softmax(logits, dim=1)

            all_probs.extend(probs.float().cpu().numpy())
            all_labels.extend(labels.cpu().numpy())

    all_probs = np.array(all_probs)
    all_labels = np.array(all_labels)
    preds = all_probs.argmax(axis=1)

    acc = accuracy_score(all_labels, preds)
    f1 = f1_score(all_labels, preds, average='weighted')
    auc = roc_auc_score(all_labels, all_probs[:, 1]) if n_classes == 2 else \
          roc_auc_score(all_labels, all_probs, multi_class='ovr', average='weighted')

    return {'accuracy': acc, 'f1': f1, 'auc': auc}


def run_lora_finetuning(config: Config, task: str = 'binary'):
    """Run LoRA fine-tuning."""
    logger.info("=" * 60)
    logger.info(f"STATE SE-600M LoRA Fine-Tuning")
    logger.info(f"Task: {task}")
    logger.info(f"LoRA config: r={config.LORA_R}, alpha={config.LORA_ALPHA}")
    logger.info("=" * 60)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logger.info(f"Device: {device}")
    if device.type == 'cuda':
        logger.info(f"GPU: {torch.cuda.get_device_name(0)}")

    config.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    config.CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)

    # Load data
    protein_embeddings = load_protein_embeddings(config.PROTEIN_EMBEDDINGS)
    expression, gene_names, labels = load_raw_covid_data(config, task)
    n_classes = config.TASK_CONFIGS[task]['n_classes']

    # Cross-validation
    skf = StratifiedKFold(n_splits=config.N_FOLDS, shuffle=True, random_state=config.SEED)
    fold_results = []

    for fold, (train_idx, val_idx) in enumerate(skf.split(expression, labels)):
        logger.info(f"\n{'='*40}\nFold {fold+1}/{config.N_FOLDS}\n{'='*40}")

        # Datasets
        train_data = RawExpressionDataset(expression[train_idx], gene_names, labels[train_idx],
                                          protein_embeddings, config.MAX_GENES)
        val_data = RawExpressionDataset(expression[val_idx], gene_names, labels[val_idx],
                                        protein_embeddings, config.MAX_GENES)

        train_loader = DataLoader(train_data, batch_size=config.BATCH_SIZE,
                                  shuffle=True, num_workers=4, pin_memory=True)
        val_loader = DataLoader(val_data, batch_size=config.BATCH_SIZE * 2,
                                shuffle=False, num_workers=4, pin_memory=True)

        # Model with LoRA
        se_model = load_se_model(config, device)
        model = SELoRAClassifier(
            se_model=se_model,
            n_classes=n_classes,
            hidden_dim=config.HIDDEN_DIM,
            dropout=config.DROPOUT,
            lora_r=config.LORA_R,
            lora_alpha=config.LORA_ALPHA,
            lora_dropout=config.LORA_DROPOUT
        ).to(device)

        # Only optimize trainable params (LoRA + classifier)
        trainable_params = [p for p in model.parameters() if p.requires_grad]
        optimizer = AdamW(trainable_params, lr=config.LEARNING_RATE, weight_decay=config.WEIGHT_DECAY)
        scheduler = CosineAnnealingWarmRestarts(optimizer, T_0=5, T_mult=2)
        scaler = torch.amp.GradScaler('cuda')

        best_auc = 0
        best_metrics = None
        patience, patience_counter = 5, 0

        for epoch in range(config.EPOCHS):
            train_loss, train_acc = train_epoch(model, train_loader, optimizer, device, scaler)
            scheduler.step()
            val_metrics = evaluate(model, val_loader, device, n_classes)

            logger.info(f"Epoch {epoch+1}/{config.EPOCHS} | "
                       f"Loss: {train_loss:.4f}, Train Acc: {train_acc:.4f} | "
                       f"Val AUC: {val_metrics['auc']:.4f}, Val Acc: {val_metrics['accuracy']:.4f}")

            if val_metrics['auc'] > best_auc:
                best_auc = val_metrics['auc']
                best_metrics = val_metrics.copy()
                patience_counter = 0
                torch.save({
                    'model_state_dict': model.state_dict(),
                    'lora_config': {'r': config.LORA_R, 'alpha': config.LORA_ALPHA},
                    'epoch': epoch,
                    'auc': best_auc
                }, config.CHECKPOINT_DIR / f"best_lora_{task}_fold{fold}.pt")
            else:
                patience_counter += 1
                if patience_counter >= patience:
                    logger.info(f"Early stopping at epoch {epoch+1}")
                    break

        fold_results.append({'fold': fold, 'auc': best_auc, 'accuracy': best_metrics['accuracy'], 'f1': best_metrics['f1']})
        logger.info(f"Fold {fold+1} Best AUC: {best_auc:.4f}")

        del model, se_model, optimizer, scaler
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

    # Summary
    aucs = [r['auc'] for r in fold_results]
    accs = [r['accuracy'] for r in fold_results]

    summary = {
        'task': task,
        'method': 'lora_finetuning',
        'lora_r': config.LORA_R,
        'lora_alpha': config.LORA_ALPHA,
        'max_genes': config.MAX_GENES,
        'mean_auc': float(np.mean(aucs)),
        'std_auc': float(np.std(aucs)),
        'mean_accuracy': float(np.mean(accs)),
        'std_accuracy': float(np.std(accs)),
        'fold_results': fold_results
    }

    logger.info(f"\n{'='*60}")
    logger.info(f"FINAL RESULTS (LoRA Fine-Tuning)")
    logger.info(f"AUC = {summary['mean_auc']:.4f} +/- {summary['std_auc']:.4f}")
    logger.info(f"Accuracy = {summary['mean_accuracy']:.4f} +/- {summary['std_accuracy']:.4f}")
    logger.info(f"{'='*60}")

    with open(config.RESULTS_DIR / f"results_lora_{task}.json", 'w') as f:
        json.dump(summary, f, indent=2)

    return summary


def main():
    parser = argparse.ArgumentParser(description='STATE SE-600M LoRA Fine-Tuning')
    parser.add_argument('--task', default='binary', choices=['binary', '3class', '6class'])
    parser.add_argument('--epochs', type=int, default=15)
    parser.add_argument('--batch-size', type=int, default=16)
    parser.add_argument('--lr', type=float, default=1e-4)
    parser.add_argument('--max-genes', type=int, default=512)
    parser.add_argument('--lora-r', type=int, default=16, help='LoRA rank (lower = fewer params)')
    parser.add_argument('--lora-alpha', type=int, default=32, help='LoRA scaling (typically 2x r)')
    parser.add_argument('--all-tasks', action='store_true')

    args = parser.parse_args()

    config = Config()
    config.EPOCHS = args.epochs
    config.BATCH_SIZE = args.batch_size
    config.LEARNING_RATE = args.lr
    config.MAX_GENES = args.max_genes
    config.LORA_R = args.lora_r
    config.LORA_ALPHA = args.lora_alpha

    torch.manual_seed(config.SEED)
    np.random.seed(config.SEED)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(config.SEED)

    if args.all_tasks:
        for task in ['binary', '3class', '6class']:
            run_lora_finetuning(config, task)
    else:
        run_lora_finetuning(config, args.task)


if __name__ == '__main__':
    main()
