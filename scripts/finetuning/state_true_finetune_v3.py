#!/usr/bin/env python3
"""
STATE SE-600M TRUE End-to-End Fine-Tuning for COVID-19 Severity Classification

This script performs TRUE fine-tuning where gradients flow through the SE-600M
transformer encoder (603M parameters), not just a classifier on frozen embeddings.

Architecture:
    [Raw Expression + Gene Names]
    → [Map to ESM2 Gene Embeddings (5120 dim)]
    → [SE-600M Encoder (5120→2048)]
    → [SE-600M Transformer (16 layers, 2048 dim)]
    → [Cell Embedding (2048 dim)]
    → [Classification Head]
    → [Predictions]

Key Difference from Embedding Classification:
    - Embedding Classification: Pre-computed frozen embeddings → Classifier (only classifier trained)
    - True Fine-tuning: Raw data → SE model (weights updated) → Classifier (everything trained)

Requirements:
    pip install arc-state torch scanpy omegaconf

Usage:
    python state_true_finetune_v3.py --task binary --epochs 15
    python state_true_finetune_v3.py --task binary --freeze-encoder  # Freeze SE, train classifier only

Author: Xinru Qiu
Date: January 2026
"""

import sys
import os
import argparse
import logging
from pathlib import Path
import json
from typing import Dict, List, Tuple, Optional

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
# Fixes PyTorch 2.6+ weights_only=True default issue
# =============================================================================
_original_torch_load = torch.load

def _patched_torch_load(f, *args, **kwargs):
    """Force weights_only=False for all torch.load calls."""
    kwargs['weights_only'] = False
    return _original_torch_load(f, *args, **kwargs)

torch.load = _patched_torch_load

# Setup VCI aliases for backward compatibility with SE checkpoints
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
    DATA_DIR = BASE_DIR / "01_DATA"
    RESULTS_DIR = BASE_DIR / "05_RESULTS/state_true_finetune_v3"
    CHECKPOINT_DIR = BASE_DIR / "05_RESULTS/checkpoints/state_true_finetune_v3"

    # SE-600M model files
    SE_CHECKPOINT = MODEL_DIR / "se600m_epoch16.ckpt"
    PROTEIN_EMBEDDINGS = MODEL_DIR / "protein_embeddings.pt"
    CONFIG_FILE = MODEL_DIR / "config.yaml"

    # Training hyperparameters
    BATCH_SIZE = 8  # Small for full model fine-tuning (memory intensive)
    LEARNING_RATE = 1e-5  # Smaller LR for fine-tuning pre-trained model
    WEIGHT_DECAY = 0.01
    EPOCHS = 15
    MAX_GENES = 256  # Top expressed genes per cell (reduced for memory - attention is O(n²))

    # Classifier head
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
# Data Loading
# =============================================================================

def load_protein_embeddings(path: Path, device: torch.device) -> Dict[str, torch.Tensor]:
    """
    Load ESM2 protein embeddings that map gene names to 5120-dim vectors.

    NOTE: Always load to CPU for DataLoader compatibility. Move to GPU during forward pass.

    Returns:
        Dictionary mapping gene names (uppercase) to embedding tensors (on CPU)
    """
    logger.info(f"Loading protein embeddings from {path}")
    # Always load to CPU - CUDA tensors can't be shared across DataLoader workers
    embeddings = torch.load(str(path), map_location='cpu', weights_only=False)

    # The protein_embeddings.pt is a dict mapping gene names to embeddings
    if isinstance(embeddings, dict):
        # Convert all keys to uppercase for consistent lookup
        # Keep tensors on CPU for DataLoader compatibility
        embeddings = {k.upper(): v.cpu() for k, v in embeddings.items()}
        logger.info(f"Loaded {len(embeddings)} gene embeddings (on CPU for DataLoader compatibility)")
        # Check embedding dimension
        sample_key = next(iter(embeddings.keys()))
        logger.info(f"Embedding dimension: {embeddings[sample_key].shape}")
    else:
        # If it's a tensor with index mapping, we need config.yaml
        logger.info(f"Embeddings type: {type(embeddings)}, shape: {embeddings.shape if hasattr(embeddings, 'shape') else 'N/A'}")

    return embeddings


def load_raw_covid_data(config: Config, task: str = 'binary'):
    """
    Load raw COVID expression data (counts) for fine-tuning.

    Returns:
        adata: AnnData with raw counts
        labels: numpy array of labels
        valid_indices: indices of cells with valid labels
    """
    logger.info("Loading raw COVID data...")

    # Try multiple possible locations for raw data
    possible_files = [
        config.EMBEDDING_DIR / "covid_state.h5ad",  # Has raw in .raw or X
        config.DATA_DIR / "covid" / "covid_raw.h5ad",
        config.DATA_DIR / "covid_severity.h5ad",
    ]

    adata = None
    for f in possible_files:
        if f.exists():
            logger.info(f"Loading from {f}")
            adata = sc.read_h5ad(f)
            break

    if adata is None:
        raise FileNotFoundError(f"Could not find raw COVID data in: {possible_files}")

    logger.info(f"Loaded AnnData: {adata.shape}")
    logger.info(f"obs columns: {list(adata.obs.columns)}")
    logger.info(f"obsm keys: {list(adata.obsm.keys())}")

    # Check if we have raw counts
    if hasattr(adata, 'raw') and adata.raw is not None:
        logger.info("Using adata.raw for raw counts")
        X = adata.raw.X
        var_names = list(adata.raw.var_names)
    else:
        logger.info("Using adata.X (assuming raw or normalized counts)")
        X = adata.X
        var_names = list(adata.var_names)

    # Convert to dense if sparse
    if hasattr(X, 'toarray'):
        logger.info("Converting sparse matrix to dense...")
        X = X.toarray()
    X = np.array(X, dtype=np.float32)

    # Get labels
    severity = adata.obs.get('covid_severity', adata.obs.get('severity'))
    if severity is None:
        raise KeyError("No severity column found in obs")

    label_map = config.TASK_CONFIGS[task]['label_map']

    labels = []
    valid_indices = []
    for i, s in enumerate(severity):
        s_lower = str(s).lower().strip()
        if s_lower in label_map:
            labels.append(label_map[s_lower])
            valid_indices.append(i)

    labels = np.array(labels)
    valid_indices = np.array(valid_indices)
    X = X[valid_indices]

    logger.info(f"Valid cells: {len(labels)}")
    logger.info(f"Expression shape: {X.shape}")
    logger.info(f"Gene names sample: {var_names[:5]}")

    return X, var_names, labels


# =============================================================================
# Dataset for True Fine-Tuning
# =============================================================================

class RawExpressionDataset(Dataset):
    """
    Dataset that provides raw expression data for SE model fine-tuning.

    For each cell:
    1. Select top expressed genes (or all non-zero)
    2. Look up ESM2 embeddings for those genes
    3. Return (gene_embeddings, expression_values, label)

    NOTE: All tensors are kept on CPU for DataLoader multiprocessing compatibility.
    Data is moved to GPU in the training loop.
    """

    def __init__(
        self,
        expression_matrix: np.ndarray,
        gene_names: List[str],
        labels: np.ndarray,
        protein_embeddings: Dict[str, torch.Tensor],
        max_genes: int = 2000
    ):
        self.expression = torch.FloatTensor(expression_matrix)
        self.gene_names = [g.upper() for g in gene_names]  # Uppercase for matching
        self.labels = torch.LongTensor(labels)
        self.protein_embeddings = protein_embeddings
        self.max_genes = max_genes

        # Pre-compute which genes have embeddings
        self.gene_to_idx = {g: i for i, g in enumerate(self.gene_names)}
        self.valid_genes = [g for g in self.gene_names if g in protein_embeddings]
        self.valid_gene_indices = [self.gene_to_idx[g] for g in self.valid_genes]

        logger.info(f"Genes with ESM2 embeddings: {len(self.valid_genes)}/{len(self.gene_names)}")

        # Get embedding dimension from first valid gene
        if self.valid_genes:
            sample_emb = protein_embeddings[self.valid_genes[0]]
            self.embed_dim = sample_emb.shape[-1]
            logger.info(f"ESM2 embedding dimension: {self.embed_dim}")
        else:
            raise ValueError("No genes found in protein embeddings!")

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        """
        Get gene embeddings weighted by expression for a single cell.

        Returns:
            gene_embeds: [max_genes, embed_dim] tensor of ESM2 embeddings
            expr_weights: [max_genes] expression values for weighting
            mask: [max_genes] boolean mask for valid positions
            label: class label
        """
        cell_expr = self.expression[idx]
        label = self.labels[idx]

        # Get expression values for valid genes only
        valid_expr = cell_expr[self.valid_gene_indices]

        # Select top expressed genes
        if len(self.valid_genes) > self.max_genes:
            top_indices = torch.topk(valid_expr, self.max_genes).indices
        else:
            # Use all valid genes, sorted by expression
            top_indices = torch.argsort(valid_expr, descending=True)

        # Collect embeddings and expression values
        gene_embeds = []
        expr_weights = []

        for i in top_indices:
            gene_name = self.valid_genes[i]
            emb = self.protein_embeddings[gene_name]
            gene_embeds.append(emb)
            expr_weights.append(valid_expr[i])

        # Pad if necessary
        n_genes = len(gene_embeds)
        if n_genes < self.max_genes:
            padding = self.max_genes - n_genes
            zero_emb = torch.zeros(self.embed_dim)
            gene_embeds.extend([zero_emb] * padding)
            expr_weights.extend([0.0] * padding)

        gene_embeds = torch.stack(gene_embeds[:self.max_genes])  # [max_genes, embed_dim]
        expr_weights = torch.tensor(expr_weights[:self.max_genes], dtype=torch.float32)  # [max_genes]
        mask = expr_weights > 0  # [max_genes] - True for real genes

        return {
            'gene_embeds': gene_embeds,
            'expr_weights': expr_weights,
            'mask': mask,
            'label': label
        }


# =============================================================================
# Model: SE-600M with Classification Head
# =============================================================================

class SEClassifier(nn.Module):
    """
    SE-600M model with classification head for end-to-end fine-tuning.

    Architecture:
        Gene embeddings (5120) → SE Encoder (2048) → SE Transformer → Pool → Classifier
    """

    def __init__(
        self,
        se_model: StateEmbeddingModel,
        n_classes: int,
        hidden_dim: int = 512,
        dropout: float = 0.1,
        freeze_encoder: bool = False
    ):
        super().__init__()

        self.se_model = se_model
        self.freeze_encoder = freeze_encoder

        # Get output dimension from SE model
        # SE model outputs 2048-dim embeddings
        self.embed_dim = se_model.hparams.get('d_model', 2048)

        if freeze_encoder:
            logger.info("Freezing SE encoder weights")
            for param in self.se_model.parameters():
                param.requires_grad = False
        else:
            # Enable gradient checkpointing to save memory
            if hasattr(self.se_model, 'transformer_encoder'):
                self.se_model.transformer_encoder.gradient_checkpointing = True
                logger.info("Enabled gradient checkpointing for transformer encoder")

        # Classification head
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

        logger.info(f"Classifier: {self.embed_dim} → {hidden_dim} → {hidden_dim//2} → {n_classes}")

    def forward(self, gene_embeds, expr_weights, mask):
        """
        Forward pass through SE model and classifier.

        Args:
            gene_embeds: [batch, max_genes, 5120] ESM2 gene embeddings
            expr_weights: [batch, max_genes] expression values
            mask: [batch, max_genes] boolean mask for valid positions

        Returns:
            logits: [batch, n_classes]
        """
        batch_size = gene_embeds.shape[0]

        # The SE model expects:
        # 1. Gene embeddings passed through encoder (5120 → 2048)
        # 2. Then through transformer encoder

        if self.freeze_encoder:
            with torch.no_grad():
                cell_embeds = self._encode_cells(gene_embeds, expr_weights, mask)
        else:
            cell_embeds = self._encode_cells(gene_embeds, expr_weights, mask)

        # Classification
        logits = self.classifier(cell_embeds)

        return logits

    def _encode_cells(self, gene_embeds, expr_weights, mask):
        """
        Encode cells using SE model.

        Args:
            gene_embeds: [batch, max_genes, 5120]
            expr_weights: [batch, max_genes]
            mask: [batch, max_genes]

        Returns:
            cell_embeds: [batch, d_model]
        """
        batch_size, max_genes, embed_dim = gene_embeds.shape

        # Project gene embeddings through SE encoder (5120 → 2048)
        # Reshape for batch processing
        gene_embeds_flat = gene_embeds.view(batch_size * max_genes, embed_dim)
        encoded = self.se_model.encoder(gene_embeds_flat)  # [batch*max_genes, 2048]
        encoded = encoded.view(batch_size, max_genes, -1)  # [batch, max_genes, 2048]

        # Weight by expression (optional - helps with expression-aware encoding)
        # Normalize expression weights
        expr_weights_norm = expr_weights / (expr_weights.sum(dim=1, keepdim=True) + 1e-8)
        expr_weights_norm = expr_weights_norm.unsqueeze(-1)  # [batch, max_genes, 1]

        # Apply expression weighting
        encoded_weighted = encoded * expr_weights_norm

        # Create attention mask for transformer
        # mask shape: [batch, max_genes] -> need [batch, max_genes] for src_key_padding_mask
        attn_mask = ~mask  # Invert: True where we want to mask (padding)

        # Pass through transformer encoder
        # SE transformer expects [seq_len, batch, d_model]
        encoded_t = encoded_weighted.transpose(0, 1)  # [max_genes, batch, 2048]

        # Use the transformer encoder
        transformer_out = self.se_model.transformer_encoder(
            encoded_t,
            src_key_padding_mask=attn_mask
        )  # [max_genes, batch, 2048]

        # Pool over genes (masked mean pooling)
        transformer_out = transformer_out.transpose(0, 1)  # [batch, max_genes, 2048]

        # Masked mean pooling
        mask_expanded = mask.unsqueeze(-1).float()  # [batch, max_genes, 1]
        pooled = (transformer_out * mask_expanded).sum(dim=1)  # [batch, 2048]
        pooled = pooled / (mask_expanded.sum(dim=1) + 1e-8)  # Normalize by valid count

        return pooled


def load_se_model(config: Config, device: torch.device) -> StateEmbeddingModel:
    """Load the SE-600M model from checkpoint."""
    logger.info("Loading SE-600M model...")

    checkpoint_path = config.SE_CHECKPOINT
    if not checkpoint_path.exists():
        raise FileNotFoundError(f"Checkpoint not found: {checkpoint_path}")

    logger.info(f"Loading from {checkpoint_path}")
    checkpoint = torch.load(str(checkpoint_path), map_location=device, weights_only=False)

    # Get hyperparameters
    hparams = checkpoint['hyper_parameters']
    logger.info(f"Model config: d_model={hparams['d_model']}, nlayers={hparams['nlayers']}")

    # Create model
    model = StateEmbeddingModel(**hparams)

    # Load weights
    result = model.load_state_dict(checkpoint['state_dict'], strict=False)
    logger.info(f"Loaded weights. Missing: {len(result.missing_keys)}, Unexpected: {len(result.unexpected_keys)}")

    model = model.to(device)
    n_params = sum(p.numel() for p in model.parameters()) / 1e6
    logger.info(f"SE-600M loaded: {n_params:.1f}M parameters")

    return model


# =============================================================================
# Training
# =============================================================================

def train_epoch(model, loader, optimizer, device, scaler, grad_clip=1.0):
    """Train for one epoch with mixed precision."""
    model.train()
    total_loss = 0
    correct = 0
    total = 0

    for batch in loader:
        gene_embeds = batch['gene_embeds'].to(device)
        expr_weights = batch['expr_weights'].to(device)
        mask = batch['mask'].to(device)
        labels = batch['label'].to(device)

        optimizer.zero_grad()

        # Mixed precision forward pass
        with torch.amp.autocast('cuda'):
            logits = model(gene_embeds, expr_weights, mask)
            loss = F.cross_entropy(logits, labels)

        # Scaled backward pass
        scaler.scale(loss).backward()

        # Gradient clipping (unscale first)
        scaler.unscale_(optimizer)
        torch.nn.utils.clip_grad_norm_(model.parameters(), grad_clip)

        # Optimizer step with scaler
        scaler.step(optimizer)
        scaler.update()

        total_loss += loss.item()
        preds = logits.argmax(dim=1)
        correct += (preds == labels).sum().item()
        total += len(labels)

    return total_loss / len(loader), correct / total


def evaluate(model, loader, device, n_classes):
    """Evaluate model on validation set with mixed precision."""
    model.eval()
    all_probs, all_labels = [], []

    with torch.no_grad():
        for batch in loader:
            gene_embeds = batch['gene_embeds'].to(device)
            expr_weights = batch['expr_weights'].to(device)
            mask = batch['mask'].to(device)
            labels = batch['label'].to(device)

            # Mixed precision inference
            with torch.amp.autocast('cuda'):
                logits = model(gene_embeds, expr_weights, mask)
                probs = F.softmax(logits, dim=1)

            all_probs.extend(probs.float().cpu().numpy())  # Convert to float32 for numpy
            all_labels.extend(labels.cpu().numpy())

    all_probs = np.array(all_probs)
    all_labels = np.array(all_labels)
    preds = all_probs.argmax(axis=1)

    acc = accuracy_score(all_labels, preds)
    f1 = f1_score(all_labels, preds, average='weighted')

    if n_classes == 2:
        auc = roc_auc_score(all_labels, all_probs[:, 1])
    else:
        auc = roc_auc_score(all_labels, all_probs, multi_class='ovr', average='weighted')

    return {'accuracy': acc, 'f1': f1, 'auc': auc}


def run_true_finetuning(config: Config, task: str = 'binary', freeze_encoder: bool = False):
    """
    Run TRUE end-to-end fine-tuning of SE-600M.

    Args:
        config: Configuration object
        task: Classification task ('binary', '3class', '6class')
        freeze_encoder: If True, freeze SE weights and only train classifier
    """
    logger.info("=" * 60)
    logger.info(f"STATE SE-600M TRUE Fine-Tuning")
    logger.info(f"Task: {task}")
    logger.info(f"Freeze encoder: {freeze_encoder}")
    logger.info("=" * 60)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logger.info(f"Device: {device}")

    if device.type == 'cuda':
        logger.info(f"GPU: {torch.cuda.get_device_name(0)}")
        logger.info(f"Memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")

    config.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    config.CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)

    # Load protein embeddings (gene name → ESM2 embedding)
    protein_embeddings = load_protein_embeddings(config.PROTEIN_EMBEDDINGS, device)

    # Load raw expression data
    expression, gene_names, labels = load_raw_covid_data(config, task)
    n_classes = config.TASK_CONFIGS[task]['n_classes']

    # Load SE-600M model
    se_model = load_se_model(config, device)

    # Cross-validation
    skf = StratifiedKFold(n_splits=config.N_FOLDS, shuffle=True, random_state=config.SEED)
    fold_results = []

    for fold, (train_idx, val_idx) in enumerate(skf.split(expression, labels)):
        logger.info(f"\n{'='*40}")
        logger.info(f"Fold {fold+1}/{config.N_FOLDS}")
        logger.info(f"{'='*40}")

        # Create datasets (data stays on CPU, moved to GPU in training loop)
        train_data = RawExpressionDataset(
            expression[train_idx], gene_names, labels[train_idx],
            protein_embeddings, config.MAX_GENES
        )
        val_data = RawExpressionDataset(
            expression[val_idx], gene_names, labels[val_idx],
            protein_embeddings, config.MAX_GENES
        )

        train_loader = DataLoader(
            train_data, batch_size=config.BATCH_SIZE,
            shuffle=True, num_workers=4, pin_memory=True
        )
        val_loader = DataLoader(
            val_data, batch_size=config.BATCH_SIZE * 2,
            shuffle=False, num_workers=4, pin_memory=True
        )

        # Create classifier model (wraps SE model)
        # Create fresh SE model for each fold to avoid weight leakage
        se_model_fold = load_se_model(config, device)

        model = SEClassifier(
            se_model=se_model_fold,
            n_classes=n_classes,
            hidden_dim=config.HIDDEN_DIM,
            dropout=config.DROPOUT,
            freeze_encoder=freeze_encoder
        ).to(device)

        # Count trainable parameters
        trainable_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
        total_params = sum(p.numel() for p in model.parameters())
        logger.info(f"Trainable parameters: {trainable_params/1e6:.1f}M / {total_params/1e6:.1f}M")

        # Optimizer - use different LR for encoder vs classifier
        if not freeze_encoder:
            optimizer = AdamW([
                {'params': model.se_model.parameters(), 'lr': config.LEARNING_RATE * 0.1},
                {'params': model.classifier.parameters(), 'lr': config.LEARNING_RATE}
            ], weight_decay=config.WEIGHT_DECAY)
        else:
            optimizer = AdamW(
                model.classifier.parameters(),
                lr=config.LEARNING_RATE,
                weight_decay=config.WEIGHT_DECAY
            )

        scheduler = CosineAnnealingWarmRestarts(optimizer, T_0=5, T_mult=2)

        # Mixed precision scaler
        scaler = torch.amp.GradScaler('cuda')

        best_auc = 0
        best_metrics = None
        patience = 5
        patience_counter = 0

        for epoch in range(config.EPOCHS):
            train_loss, train_acc = train_epoch(model, train_loader, optimizer, device, scaler)
            scheduler.step()

            val_metrics = evaluate(model, val_loader, device, n_classes)

            logger.info(
                f"Epoch {epoch+1}/{config.EPOCHS} | "
                f"Loss: {train_loss:.4f}, Train Acc: {train_acc:.4f} | "
                f"Val AUC: {val_metrics['auc']:.4f}, Val Acc: {val_metrics['accuracy']:.4f}"
            )

            if val_metrics['auc'] > best_auc:
                best_auc = val_metrics['auc']
                best_metrics = val_metrics.copy()
                patience_counter = 0

                # Save best model
                ckpt_path = config.CHECKPOINT_DIR / f"best_{task}_fold{fold}.pt"
                torch.save({
                    'model_state_dict': model.state_dict(),
                    'optimizer_state_dict': optimizer.state_dict(),
                    'epoch': epoch,
                    'auc': best_auc,
                    'freeze_encoder': freeze_encoder
                }, ckpt_path)
            else:
                patience_counter += 1
                if patience_counter >= patience:
                    logger.info(f"Early stopping at epoch {epoch+1}")
                    break

        fold_results.append({
            'fold': fold,
            'auc': best_auc,
            'accuracy': best_metrics['accuracy'],
            'f1': best_metrics['f1']
        })
        logger.info(f"Fold {fold+1} Best AUC: {best_auc:.4f}")

        # Clean up GPU memory
        del model, se_model_fold, optimizer, scaler
        torch.cuda.empty_cache() if device.type == 'cuda' else None

    # Summary
    aucs = [r['auc'] for r in fold_results]
    accs = [r['accuracy'] for r in fold_results]

    summary = {
        'task': task,
        'freeze_encoder': freeze_encoder,
        'method': 'embedding_classification' if freeze_encoder else 'true_finetuning',
        'mean_auc': float(np.mean(aucs)),
        'std_auc': float(np.std(aucs)),
        'mean_accuracy': float(np.mean(accs)),
        'std_accuracy': float(np.std(accs)),
        'fold_results': fold_results,
        'config': {
            'batch_size': config.BATCH_SIZE,
            'learning_rate': config.LEARNING_RATE,
            'epochs': config.EPOCHS,
            'max_genes': config.MAX_GENES,
            'hidden_dim': config.HIDDEN_DIM
        }
    }

    logger.info(f"\n{'='*60}")
    logger.info(f"FINAL RESULTS ({summary['method'].upper()})")
    logger.info(f"AUC = {summary['mean_auc']:.4f} +/- {summary['std_auc']:.4f}")
    logger.info(f"Accuracy = {summary['mean_accuracy']:.4f} +/- {summary['std_accuracy']:.4f}")
    logger.info(f"{'='*60}")

    # Save results
    suffix = 'frozen' if freeze_encoder else 'finetuned'
    with open(config.RESULTS_DIR / f"results_{task}_{suffix}.json", 'w') as f:
        json.dump(summary, f, indent=2)

    return summary


def main():
    parser = argparse.ArgumentParser(description='STATE SE-600M True Fine-Tuning')
    parser.add_argument('--task', default='binary', choices=['binary', '3class', '6class'],
                        help='Classification task')
    parser.add_argument('--epochs', type=int, default=15, help='Number of epochs')
    parser.add_argument('--batch-size', type=int, default=8, help='Batch size (default: 8 for memory)')
    parser.add_argument('--lr', type=float, default=1e-5, help='Learning rate')
    parser.add_argument('--max-genes', type=int, default=256,
                        help='Max genes per cell (default: 256, attention is O(n²))')
    parser.add_argument('--freeze-encoder', action='store_true',
                        help='Freeze SE encoder and only train classifier')
    parser.add_argument('--all-tasks', action='store_true',
                        help='Run all classification tasks')

    args = parser.parse_args()

    config = Config()
    config.EPOCHS = args.epochs
    config.BATCH_SIZE = args.batch_size
    config.LEARNING_RATE = args.lr
    config.MAX_GENES = args.max_genes

    torch.manual_seed(config.SEED)
    np.random.seed(config.SEED)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(config.SEED)

    if args.all_tasks:
        for task in ['binary', '3class', '6class']:
            run_true_finetuning(config, task, args.freeze_encoder)
    else:
        run_true_finetuning(config, args.task, args.freeze_encoder)


if __name__ == '__main__':
    main()
