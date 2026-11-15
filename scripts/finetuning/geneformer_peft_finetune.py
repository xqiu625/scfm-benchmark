#!/usr/bin/env python3
"""
Geneformer PEFT Fine-Tuning for COVID-19 Severity Classification

Uses Geneformer with LoRA adapters based on scPEFT methodology
(Nature Machine Intelligence, 2025).

The scPEFT approach for Geneformer:
    - Uses modified BertForSequenceClassification
    - Applies LoRA to attention layers via loralib
    - Freezes base model, trains only LoRA + classifier

Requirements:
    source activate /bigdata/godziklab/shared/Xinru/302006/scpeft_env
    export PYTHONPATH=/bigdata/godziklab/shared/Xinru/302006/03_MODELS/scPEFT:$PYTHONPATH

Usage:
    python geneformer_peft_finetune.py --task binary --epochs 25
    python geneformer_peft_finetune.py --task 3class --lr 0.0008
    python geneformer_peft_finetune.py --task binary --resume  # resume from checkpoint

Author: Xinru Qiu
Date: January 2026
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
from torch.utils.data import Dataset, DataLoader, WeightedRandomSampler
import scanpy as sc
from sklearn.model_selection import StratifiedKFold, StratifiedGroupKFold
from sklearn.metrics import roc_auc_score, accuracy_score, f1_score, precision_score, recall_score
from scipy.sparse import issparse

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


# =============================================================================
# Configuration
# =============================================================================

class Config:
    BASE_DIR = Path("/bigdata/godziklab/shared/Xinru/302006")
    DATA_DIR = BASE_DIR / "01_DATA"
    RESULTS_DIR = BASE_DIR / "05_RESULTS/geneformer_peft"
    CHECKPOINT_DIR = BASE_DIR / "05_RESULTS/checkpoints/geneformer_peft"

    # Geneformer model path -- check shared locations first
    _GENEFORMER_ROOT = BASE_DIR / "03_MODELS/geneformer/geneformer-V2-316M-native"
    _GENEFORMER_CANDIDATES = [
        _GENEFORMER_ROOT / "Geneformer-V1-10M",   # V1 model (10M cells, smaller - fits in GPU memory)
        _GENEFORMER_ROOT / "Geneformer-V2-104M",  # V2 model (104M cells, larger - may OOM)
        _GENEFORMER_ROOT,
        BASE_DIR / "geneformer_models/Geneformer",
    ]
    GENEFORMER_MODEL = _GENEFORMER_ROOT  # Root for token dictionary etc
    PRETRAIN_CKPT = next((p for p in _GENEFORMER_CANDIDATES if (p / "model.safetensors").exists() or (p / "pytorch_model.bin").exists()), _GENEFORMER_CANDIDATES[0])

    # Training
    BATCH_SIZE = 16  # V1 model is smaller, can use larger batch
    LEARNING_RATE = 8e-4  # From scPEFT Geneformer tutorial
    EPOCHS = 25
    MAX_INPUT_LEN = 2048  # Geneformer context length

    # LoRA config
    LORA_R = 16
    LORA_ALPHA = 32
    LORA_DROPOUT = 0.1

    N_FOLDS = 5
    SEED = 42
    RESUME = False  # Resume from checkpoint

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
        }},
        # In-hospital mortality (the project's primary binary task as of 2026-05-01).
        # fatal = 1; survived = 0 includes COVID-severe (ICU survivors) and
        # sepsis-severe (also ICU survivors), aligning the survival semantics
        # across diseases. healthy and unknown are EXCLUDED at load time —
        # healthy controls are not part of mortality prediction.
        'mortality': {'n_classes': 2, 'label_map': {
            'fatal': 1,
            'severe': 0, 'moderate': 0, 'mild': 0, 'recovered': 0,
            # 'healthy' and 'unknown' intentionally absent → cells excluded
        }},
        # Deprecated alias (kept for any pre-2026-05-01 result reproducibility).
        # Older definition excluded COVID-severe as ambiguous.
        'survived_vs_fatal': {'n_classes': 2, 'label_map': {
            'healthy': 0, 'mild': 0, 'moderate': 0, 'recovered': 0,
            'fatal': 1
        }},
    }


def load_covid_data(config: Config, task: str = 'binary'):
    """Load COVID expression data with Ensembl IDs for Geneformer tokenization.

    Priority:
      1. Geneformer-ready file (var_names already Ensembl IDs)
      2. Raw file with ensembl_id column in .var → remap var_names
      3. Raw file (gene symbols) — will warn about likely tokenization failure
    """
    logger.info("Loading COVID data...")

    # --- Try Geneformer-ready data first (prepared by prepare_data_for_geneformer.py) ---
    geneformer_ready = config.DATA_DIR / "geneformer_ready/covid_geneformer_ready.h5ad"
    if geneformer_ready.exists():
        adata = sc.read_h5ad(geneformer_ready)
        logger.info(f"Loaded Geneformer-ready data: {adata.shape} from {geneformer_ready}")
    else:
        # Fall back to generic preprocessed data
        data_file = config.DATA_DIR / "covid_pbmc_adata.h5ad"
        if not data_file.exists():
            data_file = config.BASE_DIR / "02_EMBEDDINGS/state/covid_state.h5ad"
        if not data_file.exists():
            raise FileNotFoundError(
                f"Data not found. Expected Geneformer-ready file at {geneformer_ready} "
                f"or raw data at {config.DATA_DIR / 'covid_pbmc_adata.h5ad'}"
            )
        adata = sc.read_h5ad(data_file)
        logger.info(f"Loaded raw data: {adata.shape} from {data_file}")

        # If var_names are gene symbols but ensembl_id column exists, remap
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
                    "Geneformer token dictionary uses Ensembl IDs — tokenization will likely fail. "
                    "Run prepare_data_for_geneformer.py first."
                )

    # Get severity labels
    severity = adata.obs.get('covid_severity', adata.obs.get('severity'))
    if severity is None:
        raise KeyError("No severity column found")

    label_map = config.TASK_CONFIGS[task]['label_map']

    labels = []
    valid_indices = []
    for i, s in enumerate(severity):
        s_lower = str(s).lower().strip()
        if s_lower in label_map:
            labels.append(label_map[s_lower])
            valid_indices.append(i)

    adata = adata[valid_indices].copy()
    labels = np.array(labels)

    logger.info(f"Valid cells: {len(labels)}")
    unique, counts = np.unique(labels, return_counts=True)
    for u, c in zip(unique, counts):
        logger.info(f"  Class {u}: {c} cells")

    if "Sample_ID" not in adata.obs.columns:
        raise KeyError(
            "Sample_ID missing from adata.obs — patient-level CV requires it. "
            "Re-tokenize from the canonical Option C input."
        )
    groups = adata.obs["Sample_ID"].astype(str).to_numpy()
    logger.info(f"Patient groups: {len(np.unique(groups))} unique Sample_IDs")
    return adata, labels, groups


def load_geneformer_token_dict(config: Config):
    """Load Geneformer gene-to-token dictionary."""
    candidates = [
        config.GENEFORMER_MODEL / "token_dictionary.pkl",
        config.GENEFORMER_MODEL / "geneformer" / "token_dictionary_gc104M.pkl",
        config.PRETRAIN_CKPT / "token_dictionary.pkl",
        config.BASE_DIR / "geneformer_models/Geneformer/geneformer/token_dictionary_gc104M.pkl",
    ]
    for token_dict_file in candidates:
        if token_dict_file.exists():
            import pickle
            with open(token_dict_file, 'rb') as f:
                token_dict = pickle.load(f)
            logger.info(f"Loaded token dictionary: {len(token_dict)} genes from {token_dict_file}")
            return token_dict

    logger.warning("Token dictionary not found at any expected path, using gene index fallback")
    return None


class GeneformerDataset(Dataset):
    """Dataset for Geneformer fine-tuning.

    Expects adata.var_names to be Ensembl IDs (ENSG...) that match the
    Geneformer token dictionary.  Genes not found in the dictionary are
    *skipped* (Geneformer convention — only matched genes contribute).
    """

    def __init__(self, adata, labels, token_dict=None, max_len=2048):
        self.adata = adata
        self.labels = labels
        self.token_dict = token_dict
        self.max_len = max_len
        self.gene_names = list(adata.var_names)

        # Pre-compute per-gene token IDs (None if not found)
        # This avoids repeated dict lookups in __getitem__
        self.gene_token_ids = []
        if token_dict is not None:
            matched = 0
            for name in self.gene_names:
                tid = token_dict.get(name)
                if tid is None:
                    # Also try upper-case (legacy gene-symbol fallback)
                    tid = token_dict.get(name.upper())
                self.gene_token_ids.append(tid)
                if tid is not None:
                    matched += 1

            match_rate = matched / len(self.gene_names) if self.gene_names else 0
            logger.info(
                f"GeneformerDataset: {matched}/{len(self.gene_names)} genes "
                f"matched token dictionary ({match_rate:.1%})"
            )
            if match_rate < 0.10:
                raise ValueError(
                    f"Token match rate is only {match_rate:.1%} — almost no genes match the "
                    f"Geneformer dictionary. This usually means var_names are gene symbols "
                    f"instead of Ensembl IDs. Run prepare_data_for_geneformer.py first."
                )

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        # Get expression values
        if issparse(self.adata.X):
            expr = np.array(self.adata.X[idx].todense()).flatten()
        else:
            expr = self.adata.X[idx].flatten()

        # Get non-zero genes sorted by expression (rank ordering like Geneformer)
        nonzero_idx = np.where(expr > 0)[0]
        if len(nonzero_idx) == 0:
            nonzero_idx = np.arange(min(10, len(expr)))

        sorted_idx = nonzero_idx[np.argsort(-expr[nonzero_idx])]

        # Convert to token IDs — skip genes not in dictionary
        input_ids = [2]  # CLS token
        for i in sorted_idx:
            if len(input_ids) >= self.max_len - 1:  # Leave room for SEP
                break
            if self.token_dict is not None:
                tid = self.gene_token_ids[i]
                if tid is not None:
                    input_ids.append(tid)
                # else: skip unknown gene (Geneformer convention)
            else:
                input_ids.append(i + 4)  # Offset by special tokens

        input_ids.append(3)  # SEP token

        # Padding
        attention_mask = [1] * len(input_ids)
        if len(input_ids) < self.max_len:
            pad_len = self.max_len - len(input_ids)
            input_ids.extend([0] * pad_len)
            attention_mask.extend([0] * pad_len)

        return {
            'input_ids': torch.tensor(input_ids, dtype=torch.long),
            'attention_mask': torch.tensor(attention_mask, dtype=torch.long),
            'labels': torch.tensor(self.labels[idx], dtype=torch.long)
        }


class LoRALinear(nn.Module):
    """LoRA adapter for linear layers."""

    def __init__(self, original_layer: nn.Linear, r: int = 16, alpha: int = 32, dropout: float = 0.1):
        super().__init__()
        self.original_layer = original_layer
        in_features = original_layer.in_features
        out_features = original_layer.out_features

        self.r = r
        self.alpha = alpha
        self.scaling = alpha / r

        # Freeze original weights
        for param in self.original_layer.parameters():
            param.requires_grad = False

        # LoRA matrices
        self.lora_A = nn.Parameter(torch.zeros(r, in_features))
        self.lora_B = nn.Parameter(torch.zeros(out_features, r))
        self.dropout = nn.Dropout(dropout)

        # Initialize
        nn.init.kaiming_uniform_(self.lora_A, a=np.sqrt(5))
        nn.init.zeros_(self.lora_B)

    def forward(self, x):
        original_output = self.original_layer(x)
        lora_output = self.dropout(x) @ self.lora_A.T @ self.lora_B.T * self.scaling
        return original_output + lora_output


def apply_lora_to_bert(model, r=16, alpha=32, dropout=0.1):
    """Apply LoRA to BERT attention layers."""
    lora_count = 0

    for name, module in model.named_modules():
        if isinstance(module, nn.Linear) and any(k in name for k in ['query', 'key', 'value', 'dense']):
            # Get parent module
            parent_name = '.'.join(name.split('.')[:-1])
            attr_name = name.split('.')[-1]

            parent = model
            for part in parent_name.split('.'):
                if part:
                    parent = getattr(parent, part)

            # Replace with LoRA
            lora_layer = LoRALinear(module, r=r, alpha=alpha, dropout=dropout)
            setattr(parent, attr_name, lora_layer)
            lora_count += 1

    logger.info(f"Applied LoRA to {lora_count} layers")
    return model


def run_geneformer_peft(config: Config, task: str = 'binary'):
    """Run Geneformer PEFT fine-tuning."""

    try:
        from transformers import BertForSequenceClassification, BertConfig
    except ImportError:
        logger.error("transformers not installed. Run: pip install transformers")
        raise

    logger.info("=" * 60)
    logger.info("Geneformer PEFT Fine-Tuning")
    logger.info(f"Task: {task}")
    logger.info(f"LoRA r={config.LORA_R}, alpha={config.LORA_ALPHA}")
    logger.info("=" * 60)

    torch.manual_seed(config.SEED)
    np.random.seed(config.SEED)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logger.info(f"Device: {device}")

    config.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    config.CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)

    # Load token dictionary
    token_dict = load_geneformer_token_dict(config)

    # Load data
    adata, labels, groups = load_covid_data(config, task)
    n_classes = config.TASK_CONFIGS[task]['n_classes']

    # Patient-level cross-validation (TODO §1.1 — same fix applied across the project)
    skf = StratifiedGroupKFold(n_splits=config.N_FOLDS, shuffle=True, random_state=config.SEED)
    fold_results = []

    # ---- Fix 5: collect results from previously completed folds on resume ----
    if config.RESUME:
        for prev_fold in range(config.N_FOLDS):
            done_file = config.CHECKPOINT_DIR / f"fold{prev_fold}_done_{task}.json"
            if done_file.exists():
                with open(done_file) as f:
                    fold_results.append(json.load(f))
                logger.info(f"Fold {prev_fold} already completed (loaded from {done_file.name})")

    for fold, (train_idx, val_idx) in enumerate(skf.split(np.zeros(len(labels)), labels, groups=groups)):
        # ---- Fix 5: skip completed folds on resume ----
        done_file = config.CHECKPOINT_DIR / f"fold{fold}_done_{task}.json"
        if config.RESUME and done_file.exists():
            logger.info(f"Skipping fold {fold} (already completed)")
            continue

        logger.info(f"\n{'='*60}")
        logger.info(f"Fold {fold+1}/{config.N_FOLDS}")
        logger.info(f"{'='*60}")

        train_adata = adata[train_idx].copy()
        val_adata = adata[val_idx].copy()
        train_labels = labels[train_idx]
        val_labels = labels[val_idx]

        # Create datasets
        train_dataset = GeneformerDataset(train_adata, train_labels, token_dict, config.MAX_INPUT_LEN)
        val_dataset = GeneformerDataset(val_adata, val_labels, token_dict, config.MAX_INPUT_LEN)

        # Weighted sampler for class imbalance
        class_counts = np.bincount(train_labels)
        class_weights = 1.0 / class_counts
        sample_weights = class_weights[train_labels]
        sampler = WeightedRandomSampler(sample_weights, len(sample_weights))

        train_loader = DataLoader(train_dataset, batch_size=config.BATCH_SIZE,
                                  sampler=sampler, num_workers=4, pin_memory=True)
        val_loader = DataLoader(val_dataset, batch_size=config.BATCH_SIZE,
                               shuffle=False, num_workers=4, pin_memory=True)

        # Load or create model
        model_path = config.PRETRAIN_CKPT
        if model_path.exists():
            logger.info(f"Loading pretrained model from {model_path}")
            model = BertForSequenceClassification.from_pretrained(
                model_path,
                num_labels=n_classes,
                output_attentions=False,
                output_hidden_states=False
            )
        else:
            logger.info("Creating new BERT model (no pretrained weights)")
            bert_config = BertConfig(
                vocab_size=25500,  # Geneformer vocab size
                hidden_size=256,
                num_hidden_layers=12,
                num_attention_heads=4,
                intermediate_size=1024,
                hidden_dropout_prob=0.02,
                attention_probs_dropout_prob=0.02,
                max_position_embeddings=2048,
                num_labels=n_classes
            )
            model = BertForSequenceClassification(bert_config)

        # Apply LoRA
        model = apply_lora_to_bert(model, r=config.LORA_R, alpha=config.LORA_ALPHA, dropout=config.LORA_DROPOUT)

        # Freeze non-LoRA parameters (except classifier)
        for name, param in model.named_parameters():
            if 'lora' not in name.lower() and 'classifier' not in name.lower():
                param.requires_grad = False

        # Count trainable parameters
        trainable = sum(p.numel() for p in model.parameters() if p.requires_grad)
        total = sum(p.numel() for p in model.parameters())
        logger.info(f"Trainable: {trainable:,} / {total:,} ({100*trainable/total:.2f}%)")

        model = model.to(device)

        # Optimizer
        optimizer = torch.optim.AdamW(
            filter(lambda p: p.requires_grad, model.parameters()),
            lr=config.LEARNING_RATE,
            weight_decay=0.01
        )

        # ---- Fix 4: Cosine annealing scheduler (matches scGPT PEFT) ----
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
            optimizer, T_max=config.EPOCHS, eta_min=1e-6
        )

        # Class weights for loss
        class_weights_tensor = torch.tensor(class_weights, dtype=torch.float32).to(device)
        criterion = nn.CrossEntropyLoss(weight=class_weights_tensor)

        # ---- Fix 3: Dual checkpoint paths (best + latest) ----
        ckpt_best = config.CHECKPOINT_DIR / f"ckpt_{task}_lora_fold{fold}_best.pt"
        ckpt_latest = config.CHECKPOINT_DIR / f"ckpt_{task}_lora_fold{fold}_latest.pt"

        # Training state
        best_auc = 0
        best_metrics = {}
        start_epoch = 0

        # Resume from *latest* checkpoint to avoid losing progress
        if config.RESUME and ckpt_latest.exists():
            ckpt = torch.load(ckpt_latest, map_location=device)
            model.load_state_dict(ckpt['model_state_dict'])
            optimizer.load_state_dict(ckpt['optimizer_state_dict'])
            scheduler.load_state_dict(ckpt['scheduler_state_dict'])
            start_epoch = ckpt['epoch']
            best_auc = ckpt['best_auc']
            best_metrics = ckpt.get('best_metrics', {})
            logger.info(f"Resumed from latest checkpoint: epoch {start_epoch}, best AUC {best_auc:.4f}")
        elif config.RESUME and ckpt_best.exists():
            ckpt = torch.load(ckpt_best, map_location=device)
            model.load_state_dict(ckpt['model_state_dict'])
            optimizer.load_state_dict(ckpt['optimizer_state_dict'])
            scheduler.load_state_dict(ckpt['scheduler_state_dict'])
            start_epoch = ckpt['epoch']
            best_auc = ckpt['best_auc']
            best_metrics = ckpt.get('best_metrics', {})
            logger.info(f"Resumed from best checkpoint: epoch {start_epoch}, best AUC {best_auc:.4f}")

        # ---- Fix 6: No early stopping — run all epochs ----
        for epoch in range(start_epoch, config.EPOCHS):
            # Train
            model.train()
            train_loss = 0

            for batch in train_loader:
                input_ids = batch['input_ids'].to(device)
                attention_mask = batch['attention_mask'].to(device)
                labels_batch = batch['labels'].to(device)

                optimizer.zero_grad()

                outputs = model(
                    input_ids=input_ids,
                    attention_mask=attention_mask,
                    labels=labels_batch
                )

                loss = criterion(outputs.logits, labels_batch)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                optimizer.step()

                train_loss += loss.item()

            scheduler.step()

            # Validate
            model.eval()
            val_preds = []
            val_probs = []
            val_labels_list = []

            with torch.no_grad():
                for batch in val_loader:
                    input_ids = batch['input_ids'].to(device)
                    attention_mask = batch['attention_mask'].to(device)
                    labels_batch = batch['labels']

                    outputs = model(
                        input_ids=input_ids,
                        attention_mask=attention_mask
                    )

                    probs = F.softmax(outputs.logits, dim=1)
                    preds = torch.argmax(outputs.logits, dim=1)

                    val_preds.extend(preds.cpu().numpy())
                    val_probs.extend(probs.cpu().numpy())
                    val_labels_list.extend(labels_batch.numpy())

            val_preds = np.array(val_preds)
            val_probs = np.array(val_probs)
            val_labels_arr = np.array(val_labels_list)

            # Metrics
            accuracy = accuracy_score(val_labels_arr, val_preds)
            f1 = f1_score(val_labels_arr, val_preds, average='macro')
            if n_classes == 2:
                auc = roc_auc_score(val_labels_arr, val_probs[:, 1])
            else:
                try:
                    auc = roc_auc_score(val_labels_arr, val_probs, multi_class='ovr', average='macro')
                except Exception:
                    auc = 0.5

            logger.info(f"Epoch {epoch+1}/{config.EPOCHS} - Loss: {train_loss/len(train_loader):.4f} "
                       f"- Val AUC: {auc:.4f} - Val Acc: {accuracy:.4f} - Val F1: {f1:.4f}")

            # Build checkpoint payload (shared between best and latest)
            ckpt_payload = {
                'epoch': epoch + 1,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'scheduler_state_dict': scheduler.state_dict(),
                'best_auc': best_auc,
                'best_metrics': best_metrics,
            }

            if auc > best_auc:
                best_auc = auc
                best_metrics = {
                    'auc': auc,
                    'accuracy': accuracy,
                    'f1': f1,
                    'epoch': epoch + 1
                }
                ckpt_payload['best_auc'] = best_auc
                ckpt_payload['best_metrics'] = best_metrics
                torch.save(ckpt_payload, ckpt_best)
                logger.info(f"  New best AUC — saved {ckpt_best.name}")

            # Always save latest checkpoint every epoch
            torch.save(ckpt_payload, ckpt_latest)

        # ---- Fix 5: Mark fold as done ----
        if not best_metrics:
            best_metrics = {'auc': 0.0, 'accuracy': 0.0, 'f1': 0.0}
        fold_result = {'fold': fold, **best_metrics}
        fold_results.append(fold_result)
        with open(done_file, 'w') as f:
            json.dump(fold_result, f, indent=2)
        logger.info(f"Fold {fold+1} Best AUC: {best_metrics.get('auc', 'N/A')}")

    # Summary
    aucs = [r['auc'] for r in fold_results]
    accs = [r['accuracy'] for r in fold_results]
    f1s = [r['f1'] for r in fold_results]

    summary = {
        'model': 'Geneformer',
        'peft_type': 'LoRA',
        'lora_r': config.LORA_R,
        'lora_alpha': config.LORA_ALPHA,
        'task': task,
        'mean_auc': float(np.mean(aucs)),
        'std_auc': float(np.std(aucs)),
        'mean_accuracy': float(np.mean(accs)),
        'std_accuracy': float(np.std(accs)),
        'mean_f1': float(np.mean(f1s)),
        'std_f1': float(np.std(f1s)),
        'fold_results': fold_results
    }

    logger.info(f"\n{'='*60}")
    logger.info(f"FINAL RESULTS (Geneformer LoRA)")
    logger.info(f"AUC = {summary['mean_auc']:.4f} +/- {summary['std_auc']:.4f}")
    logger.info(f"Accuracy = {summary['mean_accuracy']:.4f} +/- {summary['std_accuracy']:.4f}")
    logger.info(f"F1 = {summary['mean_f1']:.4f} +/- {summary['std_f1']:.4f}")
    logger.info(f"{'='*60}")

    # Save results
    with open(config.RESULTS_DIR / f"results_{task}_lora.json", 'w') as f:
        json.dump(summary, f, indent=2)

    return summary


def main():
    parser = argparse.ArgumentParser(description='Geneformer PEFT Fine-Tuning')
    parser.add_argument('--task', default='mortality',
                        choices=['mortality', 'binary', '3class', '6class', 'survived_vs_fatal'],
                        help='mortality (PRIMARY task as of 2026-05-01): fatal=1 vs survived=0 '
                             '(severe included as 0). Other tasks kept for reproducibility.')
    parser.add_argument('--epochs', type=int, default=25)
    parser.add_argument('--batch-size', type=int, default=20)
    parser.add_argument('--lr', type=float, default=8e-4)
    parser.add_argument('--lora-r', type=int, default=16)
    parser.add_argument('--lora-alpha', type=int, default=32)
    parser.add_argument('--max-len', type=int, default=2048)
    parser.add_argument('--resume', action='store_true', help='Resume from checkpoint')
    parser.add_argument('--all-tasks', action='store_true')

    args = parser.parse_args()

    config = Config()
    config.EPOCHS = args.epochs
    config.BATCH_SIZE = args.batch_size
    config.LEARNING_RATE = args.lr
    config.LORA_R = args.lora_r
    config.RESUME = args.resume
    config.LORA_ALPHA = args.lora_alpha
    config.MAX_INPUT_LEN = args.max_len

    if args.all_tasks:
        for task in ['binary', '3class', '6class']:
            run_geneformer_peft(config, task)
    else:
        run_geneformer_peft(config, args.task)


if __name__ == '__main__':
    main()
