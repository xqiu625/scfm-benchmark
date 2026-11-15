#!/usr/bin/env python3
"""
scGPT PEFT Fine-Tuning for COVID-19 Severity Classification

Uses scGPT's native PEFT support with adapter configurations from scPEFT
(Nature Machine Intelligence, 2025).

PEFT types:
    - LORA: Low-rank adaptation of attention layers
    - PREFIX: Learnable prefix tokens
    - ENCODER: Bottleneck adapters in transformer layers
    - TOKEN: Token-level adapters
    - HYBRID: Combination of multiple adapters

Requirements:
    source activate /bigdata/godziklab/shared/Xinru/302006/scpeft_env
    export PYTHONPATH=/bigdata/godziklab/shared/Xinru/302006/03_MODELS/scPEFT:$PYTHONPATH

Usage:
    python scgpt_peft_finetune.py --task binary --peft lora
    python scgpt_peft_finetune.py --task binary --peft prefix --epochs 25
    python scgpt_peft_finetune.py --task binary --peft lora --resume  # resume from checkpoint

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
from torch.utils.data import Dataset, DataLoader
import scanpy as sc
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, accuracy_score, f1_score
from scipy.sparse import issparse

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


# =============================================================================
# Configuration
# =============================================================================

class Config:
    BASE_DIR = Path("/bigdata/godziklab/shared/Xinru/302006")
    DATA_DIR = BASE_DIR / "01_DATA"
    RESULTS_DIR = BASE_DIR / "05_RESULTS/scgpt_peft"
    CHECKPOINT_DIR = BASE_DIR / "05_RESULTS/checkpoints/scgpt_peft"

    # scGPT model path -- check scPEFT location first, then standalone
    _SCGPT_CANDIDATES = [
        BASE_DIR / "03_MODELS/scPEFT/save/scGPT_human",
        BASE_DIR / "03_MODELS/scGPT/scGPT_human",
    ]
    SCGPT_MODEL = next((p for p in _SCGPT_CANDIDATES if p.exists()), _SCGPT_CANDIDATES[0])

    # Training
    BATCH_SIZE = 16
    LEARNING_RATE = 1e-4
    EPOCHS = 20
    MAX_SEQ_LEN = 1200  # Max genes per cell

    # PEFT config
    PEFT_TYPE = 'LORA'  # Options: False, 'LORA', 'PREFIX', 'ENCODER', 'TOKEN', 'HYBRID'

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
        }}
    }


def load_covid_data(config: Config, task: str = 'binary'):
    """Load COVID expression data for scGPT fine-tuning."""
    logger.info("Loading COVID data...")

    # Load preprocessed data
    data_file = config.DATA_DIR / "covid_pbmc_adata.h5ad"
    if not data_file.exists():
        # Try alternate location
        data_file = config.BASE_DIR / "02_EMBEDDINGS/state/covid_state.h5ad"

    if not data_file.exists():
        raise FileNotFoundError(f"Data file not found: {data_file}")

    adata = sc.read_h5ad(data_file)
    logger.info(f"Loaded: {adata.shape}")

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

    return adata, labels


class scGPTDataset(Dataset):
    """Dataset for scGPT fine-tuning."""

    def __init__(self, adata, labels, vocab, max_seq_len=1200):
        self.adata = adata
        self.labels = labels
        self.vocab = vocab
        self.max_seq_len = max_seq_len
        self.gene_names = list(adata.var_names)

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        # Get expression values
        if issparse(self.adata.X):
            expr = np.array(self.adata.X[idx].todense()).flatten()
        else:
            expr = self.adata.X[idx].flatten()

        # Get non-zero genes sorted by expression
        nonzero_idx = np.where(expr > 0)[0]
        if len(nonzero_idx) == 0:
            nonzero_idx = np.arange(min(10, len(expr)))

        sorted_idx = nonzero_idx[np.argsort(-expr[nonzero_idx])]
        sorted_idx = sorted_idx[:self.max_seq_len]

        # Convert to token IDs
        gene_ids = []
        values = []
        for i in sorted_idx:
            gene = self.gene_names[i].upper()
            if gene in self.vocab:
                gene_ids.append(self.vocab[gene])
                values.append(expr[i])

        # Pad to max length
        n_genes = len(gene_ids)
        if n_genes < self.max_seq_len:
            pad_len = self.max_seq_len - n_genes
            gene_ids.extend([0] * pad_len)  # PAD token
            values.extend([0.0] * pad_len)

        return {
            'gene_ids': torch.tensor(gene_ids, dtype=torch.long),
            'values': torch.tensor(values, dtype=torch.float32),
            'label': torch.tensor(self.labels[idx], dtype=torch.long),
            'length': torch.tensor(n_genes, dtype=torch.long)
        }


def run_scgpt_peft(config: Config, task: str = 'binary'):
    """Run scGPT PEFT fine-tuning."""

    # Import scGPT
    try:
        from scgpt.model import TransformerModel
        from scgpt.tokenizer.gene_tokenizer import GeneVocab
        from scgpt.utils import set_seed
        # Try to import PeftConfig from scPEFT or scgpt
        try:
            from scgpt.utils import PeftConfig, freeze_parameters, DownstreamTasks
            HAS_PEFT_CONFIG = True
        except ImportError:
            logger.warning("PeftConfig not found in scgpt. Using manual PEFT setup.")
            HAS_PEFT_CONFIG = False
    except ImportError:
        logger.error("scGPT not installed. Run: pip install scgpt")
        raise

    logger.info("=" * 60)
    logger.info("scGPT PEFT Fine-Tuning")
    logger.info(f"Task: {task}")
    logger.info(f"PEFT Type: {config.PEFT_TYPE}")
    logger.info("=" * 60)

    set_seed(config.SEED)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logger.info(f"Device: {device}")

    config.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    config.CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)

    # Load vocabulary
    vocab_file = config.SCGPT_MODEL / "vocab.json"
    if vocab_file.exists():
        vocab = GeneVocab.from_file(vocab_file)
        logger.info(f"Loaded vocabulary: {len(vocab)} genes")
    else:
        raise FileNotFoundError(f"Vocabulary not found: {vocab_file}")

    # Load data
    adata, labels = load_covid_data(config, task)
    n_classes = config.TASK_CONFIGS[task]['n_classes']

    # Cross-validation
    skf = StratifiedKFold(n_splits=config.N_FOLDS, shuffle=True, random_state=config.SEED)
    fold_results = []

    for fold, (train_idx, val_idx) in enumerate(skf.split(np.zeros(len(labels)), labels)):
        logger.info(f"\n{'='*60}")
        logger.info(f"Fold {fold+1}/{config.N_FOLDS}")
        logger.info(f"{'='*60}")

        train_adata = adata[train_idx].copy()
        val_adata = adata[val_idx].copy()
        train_labels = labels[train_idx]
        val_labels = labels[val_idx]

        # Create datasets
        train_dataset = scGPTDataset(train_adata, train_labels, vocab, config.MAX_SEQ_LEN)
        val_dataset = scGPTDataset(val_adata, val_labels, vocab, config.MAX_SEQ_LEN)

        train_loader = DataLoader(train_dataset, batch_size=config.BATCH_SIZE,
                                  shuffle=True, num_workers=4, pin_memory=True)
        val_loader = DataLoader(val_dataset, batch_size=config.BATCH_SIZE,
                               shuffle=False, num_workers=4, pin_memory=True)

        # Initialize model with PEFT config
        if HAS_PEFT_CONFIG:
            peft_config = PeftConfig(peft_type=config.PEFT_TYPE, use_default_settings=True).to_dict()
        else:
            # Manual PEFT config for LoRA
            peft_config = {
                'peft_type': config.PEFT_TYPE,
                'lora_r': 16,
                'lora_alpha': 32,
                'lora_dropout': 0.1,
            } if config.PEFT_TYPE else None

        # Model parameters from scGPT whole-human
        model = TransformerModel(
            ntoken=len(vocab),
            d_model=512,
            nhead=8,
            d_hid=512,
            nlayers=12,
            vocab=vocab,
            dropout=0.2,
            pad_token="<pad>",
            pad_value=-2,
            n_input_bins=51,
            explicit_zero_prob=False,
            use_fast_transformer=True,
            peft_config=peft_config,
            n_cls=n_classes,
        )

        # Load pretrained weights
        model_file = config.SCGPT_MODEL / "best_model.pt"
        if model_file.exists():
            pretrained = torch.load(model_file, map_location='cpu')
            # Load compatible weights
            model_dict = model.state_dict()
            pretrained_dict = {k: v for k, v in pretrained.items()
                             if k in model_dict and v.shape == model_dict[k].shape}
            model_dict.update(pretrained_dict)
            model.load_state_dict(model_dict, strict=False)
            logger.info(f"Loaded {len(pretrained_dict)}/{len(model_dict)} pretrained weights")

        # Freeze base model parameters for PEFT
        if config.PEFT_TYPE and HAS_PEFT_CONFIG:
            freeze_parameters(model, DownstreamTasks.Identification)
        elif config.PEFT_TYPE:
            # Manual freeze: only train PEFT parameters
            for name, param in model.named_parameters():
                if 'lora' not in name.lower() and 'adapter' not in name.lower() and 'cls' not in name.lower():
                    param.requires_grad = False

        # Count trainable parameters
        trainable = sum(p.numel() for p in model.parameters() if p.requires_grad)
        total = sum(p.numel() for p in model.parameters())
        logger.info(f"Trainable: {trainable:,} / {total:,} ({100*trainable/total:.2f}%)")

        model = model.to(device)

        # Optimizer and scheduler
        optimizer = torch.optim.AdamW(
            filter(lambda p: p.requires_grad, model.parameters()),
            lr=config.LEARNING_RATE,
            weight_decay=0.01
        )
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
            optimizer, T_max=config.EPOCHS, eta_min=1e-6
        )

        # Checkpoint path for this fold/task
        ckpt_path = config.CHECKPOINT_DIR / f"ckpt_{task}_{config.PEFT_TYPE}_fold{fold}.pt"

        # Training state
        best_auc = 0
        best_metrics = {}
        start_epoch = 0

        # Resume from checkpoint if requested
        if config.RESUME and ckpt_path.exists():
            ckpt = torch.load(ckpt_path, map_location=device)
            model.load_state_dict(ckpt['model_state_dict'])
            optimizer.load_state_dict(ckpt['optimizer_state_dict'])
            scheduler.load_state_dict(ckpt['scheduler_state_dict'])
            start_epoch = ckpt['epoch']
            best_auc = ckpt['best_auc']
            best_metrics = ckpt.get('best_metrics', {})
            logger.info(f"Resumed from checkpoint: epoch {start_epoch}, best AUC {best_auc:.4f}")

        for epoch in range(start_epoch, config.EPOCHS):
            # Train
            model.train()
            train_loss = 0

            for batch in train_loader:
                gene_ids = batch['gene_ids'].to(device)
                values = batch['values'].to(device)
                labels_batch = batch['label'].to(device)

                # Create padding mask
                src_key_padding_mask = gene_ids == 0

                optimizer.zero_grad()

                output = model(
                    gene_ids,
                    values,
                    src_key_padding_mask=src_key_padding_mask,
                    CLS=True,  # Classification task
                    CCE=False,
                    MVC=False,
                    ECS=False
                )

                logits = output['cls_output']
                loss = F.cross_entropy(logits, labels_batch)

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
                    gene_ids = batch['gene_ids'].to(device)
                    values = batch['values'].to(device)
                    labels_batch = batch['label']

                    src_key_padding_mask = gene_ids == 0

                    output = model(
                        gene_ids,
                        values,
                        src_key_padding_mask=src_key_padding_mask,
                        CLS=True,
                        CCE=False,
                        MVC=False,
                        ECS=False
                    )

                    logits = output['cls_output']
                    probs = F.softmax(logits, dim=1)
                    preds = torch.argmax(logits, dim=1)

                    val_preds.extend(preds.cpu().numpy())
                    val_probs.extend(probs.cpu().numpy())
                    val_labels_list.extend(labels_batch.numpy())

            val_preds = np.array(val_preds)
            val_probs = np.array(val_probs)
            val_labels_arr = np.array(val_labels_list)

            # Metrics
            accuracy = accuracy_score(val_labels_arr, val_preds)
            if n_classes == 2:
                auc = roc_auc_score(val_labels_arr, val_probs[:, 1])
            else:
                auc = roc_auc_score(val_labels_arr, val_probs, multi_class='ovr', average='macro')

            logger.info(f"Epoch {epoch+1}/{config.EPOCHS} - Loss: {train_loss/len(train_loader):.4f} "
                       f"- Val AUC: {auc:.4f} - Val Acc: {accuracy:.4f}")

            if auc > best_auc:
                best_auc = auc
                best_metrics = {
                    'auc': auc,
                    'accuracy': accuracy,
                    'epoch': epoch + 1
                }
                # Save checkpoint on best metric
                torch.save({
                    'epoch': epoch + 1,
                    'model_state_dict': model.state_dict(),
                    'optimizer_state_dict': optimizer.state_dict(),
                    'scheduler_state_dict': scheduler.state_dict(),
                    'best_auc': best_auc,
                    'best_metrics': best_metrics,
                }, ckpt_path)
                logger.info(f"Saved checkpoint: {ckpt_path.name}")

        fold_results.append({
            'fold': fold,
            **best_metrics
        })
        logger.info(f"Fold {fold+1} Best AUC: {best_metrics['auc']:.4f}")

    # Summary
    aucs = [r['auc'] for r in fold_results]
    accs = [r['accuracy'] for r in fold_results]

    summary = {
        'model': 'scGPT',
        'peft_type': config.PEFT_TYPE,
        'task': task,
        'mean_auc': float(np.mean(aucs)),
        'std_auc': float(np.std(aucs)),
        'mean_accuracy': float(np.mean(accs)),
        'std_accuracy': float(np.std(accs)),
        'fold_results': fold_results
    }

    logger.info(f"\n{'='*60}")
    logger.info(f"FINAL RESULTS (scGPT {config.PEFT_TYPE})")
    logger.info(f"AUC = {summary['mean_auc']:.4f} +/- {summary['std_auc']:.4f}")
    logger.info(f"Accuracy = {summary['mean_accuracy']:.4f} +/- {summary['std_accuracy']:.4f}")
    logger.info(f"{'='*60}")

    # Save results
    with open(config.RESULTS_DIR / f"results_{task}_{config.PEFT_TYPE}.json", 'w') as f:
        json.dump(summary, f, indent=2)

    return summary


def main():
    parser = argparse.ArgumentParser(description='scGPT PEFT Fine-Tuning')
    parser.add_argument('--task', default='binary', choices=['binary', '3class', '6class'])
    parser.add_argument('--peft', default='LORA', choices=['LORA', 'PREFIX', 'ENCODER', 'TOKEN', 'HYBRID', 'none'])
    parser.add_argument('--epochs', type=int, default=20)
    parser.add_argument('--batch-size', type=int, default=16)
    parser.add_argument('--lr', type=float, default=1e-4)
    parser.add_argument('--max-seq-len', type=int, default=1200)
    parser.add_argument('--resume', action='store_true', help='Resume from checkpoint')
    parser.add_argument('--all-tasks', action='store_true')

    args = parser.parse_args()

    config = Config()
    config.PEFT_TYPE = args.peft if args.peft != 'none' else False
    config.EPOCHS = args.epochs
    config.BATCH_SIZE = args.batch_size
    config.LEARNING_RATE = args.lr
    config.MAX_SEQ_LEN = args.max_seq_len
    config.RESUME = args.resume

    torch.manual_seed(config.SEED)
    np.random.seed(config.SEED)

    if args.all_tasks:
        for task in ['binary', '3class', '6class']:
            run_scgpt_peft(config, task)
    else:
        run_scgpt_peft(config, args.task)


if __name__ == '__main__':
    main()
