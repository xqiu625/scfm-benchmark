#!/usr/bin/env python3
"""
Re-run finetuned benchmark for scGPT only (with cleaned embeddings)

Uses covid_scgpt_clean.h5ad (platelet data only, Orion data removed)

Author: Project 302006
Date: January 2026
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, classification_report
from sklearn.pipeline import Pipeline
from sklearn.decomposition import PCA
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def load_scgpt_clean(base_dir):
    """Load cleaned scGPT embeddings"""
    filepath = Path(base_dir) / "02_EMBEDDINGS/scgpt/covid_scgpt_clean.h5ad"

    if not filepath.exists():
        raise FileNotFoundError(f"Cleaned scGPT file not found: {filepath}")

    logger.info(f"Loading cleaned scGPT from: {filepath}")
    adata = sc.read_h5ad(filepath)

    # Find embedding key
    emb_key = None
    for key in ['X_scGPT', 'X_scgpt', 'scgpt']:
        if key in adata.obsm:
            emb_key = key
            break

    if emb_key is None:
        raise ValueError(f"No scGPT embedding found. Available: {list(adata.obsm.keys())}")

    embeddings = adata.obsm[emb_key]
    logger.info(f"Embeddings shape: {embeddings.shape}")

    # Get labels
    severity_col = None
    for col in ['three_group', 'mild_vs_severe', 'covid_severity']:
        if col in adata.obs.columns:
            severity_col = col
            break

    if severity_col is None:
        raise ValueError("No severity column found")

    labels = adata.obs[severity_col].values
    logger.info(f"Using label column: {severity_col}")
    logger.info(f"Label distribution: {pd.Series(labels).value_counts().to_dict()}")

    # Extract patient-level groups for GroupKFold
    groups = None
    for col in ['Sample_ID', 'sample_id', 'Data_NO', 'data_no', 'batch', 'donor']:
        if col in adata.obs.columns:
            groups = adata.obs[col].values
            logger.info(f"Using group column: {col} ({len(np.unique(groups))} groups)")
            break
    if groups is None:
        logger.warning("No group column found; StratifiedGroupKFold will fall back to cell-level splitting")

    return embeddings, labels, adata.obs, groups


def run_benchmark(embeddings, labels, task='binary', n_folds=5, groups=None):
    """Run classification benchmark"""

    # Prepare labels
    if task == 'binary':
        # Map to binary: severe vs non-severe
        label_map = {
            'healthy': 'non_severe', 'mild': 'non_severe', 'moderate': 'non_severe',
            'severe': 'severe', 'fatal': 'severe',
            'non_severe_group': 'non_severe', 'severe_group': 'severe',
            'healthy_group': 'non_severe', 'mild_group': 'non_severe'
        }
        labels_mapped = np.array([label_map.get(str(l).lower().replace(' ', '_'), 'unknown') for l in labels])
        mask = labels_mapped != 'unknown'
        X = embeddings[mask]
        y = labels_mapped[mask]
        groups_valid = groups[mask] if groups is not None else None
    else:
        # 3-class: healthy, mild, severe
        label_map = {
            'healthy': 'healthy', 'mild': 'mild', 'severe': 'severe',
            'healthy_group': 'healthy', 'mild_group': 'mild', 'severe_group': 'severe',
            'non_severe_group': 'mild', 'moderate': 'mild', 'fatal': 'severe'
        }
        labels_mapped = np.array([label_map.get(str(l).lower().replace(' ', '_'), 'unknown') for l in labels])
        mask = labels_mapped != 'unknown'
        X = embeddings[mask]
        y = labels_mapped[mask]
        groups_valid = groups[mask] if groups is not None else None

    logger.info(f"Task: {task}, Samples: {len(X)}, Classes: {np.unique(y)}")

    # Encode labels
    le = LabelEncoder()
    y_encoded = le.fit_transform(y)

    # Cross-validation
    skf = StratifiedGroupKFold(n_splits=n_folds, shuffle=True, random_state=42)
    aucs = []

    for fold, (train_idx, test_idx) in enumerate(skf.split(X, y_encoded, groups=groups_valid)):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y_encoded[train_idx], y_encoded[test_idx]

        # Pipeline: StandardScaler + LogisticRegression
        pipeline = Pipeline([
            ('scaler', StandardScaler()),
            ('clf', LogisticRegression(max_iter=1000, class_weight='balanced', random_state=42))
        ])

        pipeline.fit(X_train, y_train)

        if len(np.unique(y_encoded)) == 2:
            y_prob = pipeline.predict_proba(X_test)[:, 1]
            auc = roc_auc_score(y_test, y_prob)
        else:
            y_prob = pipeline.predict_proba(X_test)
            auc = roc_auc_score(y_test, y_prob, multi_class='ovr', average='macro')

        aucs.append(auc)
        logger.info(f"  Fold {fold+1}: AUC = {auc:.4f}")

    mean_auc = np.mean(aucs)
    std_auc = np.std(aucs)
    logger.info(f"  Mean AUC: {mean_auc:.4f} ± {std_auc:.4f}")

    return mean_auc, std_auc, aucs


def run_with_pca(embeddings, labels, task='binary', n_components=256, n_folds=5, groups=None):
    """Run benchmark with PCA preprocessing"""
    logger.info(f"Running with PCA-{n_components}...")

    # Apply PCA
    pca = PCA(n_components=n_components, random_state=42)
    embeddings_pca = pca.fit_transform(embeddings)

    return run_benchmark(embeddings_pca, labels, task, n_folds, groups=groups)


def main():
    base_dir = Path("/bigdata/godziklab/shared/Xinru/302006")
    output_dir = base_dir / "05_RESULTS/scgpt_rerun"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load cleaned scGPT
    embeddings, labels, obs, groups = load_scgpt_clean(base_dir)

    results = []

    # Binary task
    logger.info("\n" + "="*60)
    logger.info("BINARY TASK (severe vs non-severe)")
    logger.info("="*60)

    auc, std, _ = run_benchmark(embeddings, labels, task='binary', groups=groups)
    results.append({'model': 'scGPT_clean', 'task': 'binary', 'method': 'raw', 'auc': auc, 'std': std})

    auc_pca, std_pca, _ = run_with_pca(embeddings, labels, task='binary', n_components=256, groups=groups)
    results.append({'model': 'scGPT_clean', 'task': 'binary', 'method': 'PCA-256', 'auc': auc_pca, 'std': std_pca})

    # 3-class task
    logger.info("\n" + "="*60)
    logger.info("3-CLASS TASK (healthy vs mild vs severe)")
    logger.info("="*60)

    auc3, std3, _ = run_benchmark(embeddings, labels, task='3class', groups=groups)
    results.append({'model': 'scGPT_clean', 'task': '3class', 'method': 'raw', 'auc': auc3, 'std': std3})

    auc3_pca, std3_pca, _ = run_with_pca(embeddings, labels, task='3class', n_components=256, groups=groups)
    results.append({'model': 'scGPT_clean', 'task': '3class', 'method': 'PCA-256', 'auc': auc3_pca, 'std': std3_pca})

    # Save results
    df = pd.DataFrame(results)
    output_file = output_dir / "scgpt_clean_finetuned_results.csv"
    df.to_csv(output_file, index=False)
    logger.info(f"\nResults saved to: {output_file}")

    # Print summary
    logger.info("\n" + "="*60)
    logger.info("SUMMARY: scGPT (Cleaned) Results")
    logger.info("="*60)
    print(df.to_string(index=False))

    # Compare with previous results
    logger.info("\n" + "="*60)
    logger.info("COMPARISON: Before vs After Cleaning")
    logger.info("="*60)
    logger.info("Previous scGPT (contaminated):")
    logger.info("  Binary: 0.735 AUC")
    logger.info("  Binary + PCA-256: 0.827 AUC (+7%)")
    logger.info(f"\nNew scGPT (cleaned):")
    logger.info(f"  Binary: {auc:.3f} AUC")
    logger.info(f"  Binary + PCA-256: {auc_pca:.3f} AUC")

    return 0


if __name__ == "__main__":
    sys.exit(main())
