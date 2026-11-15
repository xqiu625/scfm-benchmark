#!/usr/bin/env python3
"""
Re-run differential diagnosis for scGPT only (with cleaned embeddings)

5-class: COVID-19, Sepsis, SLE, Pneumonia, Healthy
Uses covid_scgpt_clean.h5ad (platelet data only)

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
from sklearn.metrics import roc_auc_score, confusion_matrix, classification_report
from sklearn.pipeline import Pipeline
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


# Disease mapping
DISEASE_MAP = {
    'covid-19': 'COVID-19', 'covid': 'COVID-19', 'covid19': 'COVID-19',
    'sepsis': 'Sepsis',
    'sle': 'SLE', 'lupus': 'SLE',
    'healthy': 'Healthy', 'hc': 'Healthy', 'control': 'Healthy',
    'pneumonia': 'Pneumonia', 'lung_infect': 'Pneumonia',
    'influenza': 'Influenza', 'flu': 'Influenza',
    'non_covid': 'Non_COVID'
}

TARGET_CLASSES = ['Healthy', 'COVID-19', 'Sepsis', 'SLE', 'Pneumonia']


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

    # Get disease labels
    disease_col = 'Disease'
    if disease_col not in adata.obs.columns:
        raise ValueError(f"Disease column not found")

    diseases = adata.obs[disease_col].values
    logger.info(f"Disease distribution: {pd.Series(diseases).value_counts().to_dict()}")

    # Extract patient-level groups for GroupKFold
    groups = None
    for col in ['Sample_ID', 'sample_id', 'Data_NO', 'data_no', 'batch', 'donor']:
        if col in adata.obs.columns:
            groups = adata.obs[col].values
            logger.info(f"Using group column: {col} ({len(np.unique(groups))} groups)")
            break
    if groups is None:
        logger.warning("No group column found; StratifiedGroupKFold will fall back to cell-level splitting")

    return embeddings, diseases, adata.obs, groups


def standardize_disease(diseases):
    """Standardize disease names"""
    standardized = []
    for d in diseases:
        d_lower = str(d).lower().strip()
        mapped = DISEASE_MAP.get(d_lower, None)
        if mapped is None:
            # Try partial match
            for key, val in DISEASE_MAP.items():
                if key in d_lower:
                    mapped = val
                    break
        standardized.append(mapped if mapped else 'Other')
    return np.array(standardized)


def run_differential_diagnosis(embeddings, diseases, n_folds=5, use_pca=False, n_components=256, groups=None):
    """Run 5-class differential diagnosis"""

    # Standardize disease names
    diseases_std = standardize_disease(diseases)

    # Filter to target classes
    mask = np.isin(diseases_std, TARGET_CLASSES)
    X = embeddings[mask]
    y = diseases_std[mask]
    groups_valid = groups[mask] if groups is not None else None

    logger.info(f"Filtered to {len(X)} samples")
    logger.info(f"Class distribution: {pd.Series(y).value_counts().to_dict()}")

    # Apply PCA if requested
    if use_pca:
        logger.info(f"Applying PCA-{n_components}")
        pca = PCA(n_components=n_components, random_state=42)
        X = pca.fit_transform(X)

    # Encode labels
    le = LabelEncoder()
    le.fit(TARGET_CLASSES)
    y_encoded = le.transform(y)

    # Cross-validation
    skf = StratifiedGroupKFold(n_splits=n_folds, shuffle=True, random_state=42)

    all_y_true = []
    all_y_pred = []
    all_y_prob = []
    fold_aucs = []

    for fold, (train_idx, test_idx) in enumerate(skf.split(X, y_encoded, groups=groups_valid)):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y_encoded[train_idx], y_encoded[test_idx]

        # Pipeline
        pipeline = Pipeline([
            ('scaler', StandardScaler()),
            ('clf', LogisticRegression(max_iter=1000, class_weight='balanced', random_state=42))
        ])

        pipeline.fit(X_train, y_train)

        y_pred = pipeline.predict(X_test)
        y_prob = pipeline.predict_proba(X_test)

        all_y_true.extend(y_test)
        all_y_pred.extend(y_pred)
        all_y_prob.extend(y_prob)

        # Fold AUC
        auc = roc_auc_score(y_test, y_prob, multi_class='ovr', average='macro')
        fold_aucs.append(auc)
        logger.info(f"  Fold {fold+1}: Macro AUC = {auc:.4f}")

    all_y_true = np.array(all_y_true)
    all_y_pred = np.array(all_y_pred)
    all_y_prob = np.array(all_y_prob)

    # Overall metrics
    macro_auc = roc_auc_score(all_y_true, all_y_prob, multi_class='ovr', average='macro')
    weighted_auc = roc_auc_score(all_y_true, all_y_prob, multi_class='ovr', average='weighted')

    # Per-class AUC
    per_class_auc = {}
    for i, cls in enumerate(le.classes_):
        y_true_binary = (all_y_true == i).astype(int)
        y_prob_class = all_y_prob[:, i]
        per_class_auc[cls] = roc_auc_score(y_true_binary, y_prob_class)

    logger.info(f"\nOverall Macro AUC: {macro_auc:.4f}")
    logger.info(f"Overall Weighted AUC: {weighted_auc:.4f}")
    logger.info(f"Per-class AUC: {per_class_auc}")

    return {
        'macro_auc': macro_auc,
        'weighted_auc': weighted_auc,
        'per_class_auc': per_class_auc,
        'fold_aucs': fold_aucs,
        'y_true': all_y_true,
        'y_pred': all_y_pred,
        'label_encoder': le
    }


def plot_confusion_matrix(results, output_path, title='scGPT (Cleaned)'):
    """Plot confusion matrix"""
    y_true = results['y_true']
    y_pred = results['y_pred']
    le = results['label_encoder']

    cm = confusion_matrix(y_true, y_pred)
    cm_norm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]

    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(cm_norm, annot=True, fmt='.2f', cmap='Blues',
                xticklabels=le.classes_, yticklabels=le.classes_, ax=ax)
    ax.set_xlabel('Predicted')
    ax.set_ylabel('True')
    ax.set_title(f'Confusion Matrix: {title}\nMacro AUC = {results["macro_auc"]:.3f}')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f"Saved: {output_path}")


def main():
    base_dir = Path("/bigdata/godziklab/shared/Xinru/302006")
    output_dir = base_dir / "05_RESULTS/scgpt_rerun"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load cleaned scGPT
    embeddings, diseases, obs, groups = load_scgpt_clean(base_dir)

    results_list = []

    # Raw embeddings
    logger.info("\n" + "="*60)
    logger.info("scGPT (Cleaned) - Raw Embeddings")
    logger.info("="*60)
    results_raw = run_differential_diagnosis(embeddings, diseases, use_pca=False, groups=groups)
    results_list.append({
        'model': 'scGPT_clean',
        'method': 'raw',
        'macro_auc': results_raw['macro_auc'],
        'weighted_auc': results_raw['weighted_auc'],
        **{f'auc_{k}': v for k, v in results_raw['per_class_auc'].items()}
    })

    plot_confusion_matrix(results_raw, output_dir / 'scgpt_clean_confusion_raw.png', 'scGPT Clean (Raw)')

    # PCA-256
    logger.info("\n" + "="*60)
    logger.info("scGPT (Cleaned) - PCA-256")
    logger.info("="*60)
    results_pca = run_differential_diagnosis(embeddings, diseases, use_pca=True, n_components=256, groups=groups)
    results_list.append({
        'model': 'scGPT_clean_PCA256',
        'method': 'PCA-256',
        'macro_auc': results_pca['macro_auc'],
        'weighted_auc': results_pca['weighted_auc'],
        **{f'auc_{k}': v for k, v in results_pca['per_class_auc'].items()}
    })

    plot_confusion_matrix(results_pca, output_dir / 'scgpt_clean_confusion_pca256.png', 'scGPT Clean (PCA-256)')

    # Save results
    df = pd.DataFrame(results_list)
    output_file = output_dir / "scgpt_clean_differential_diagnosis.csv"
    df.to_csv(output_file, index=False)
    logger.info(f"\nResults saved to: {output_file}")

    # Print summary
    logger.info("\n" + "="*60)
    logger.info("SUMMARY: scGPT (Cleaned) Differential Diagnosis")
    logger.info("="*60)
    print(df.to_string(index=False))

    # Compare with previous
    logger.info("\n" + "="*60)
    logger.info("COMPARISON: Before vs After Cleaning")
    logger.info("="*60)
    logger.info("Previous scGPT (contaminated):")
    logger.info("  Macro AUC: 0.973")
    logger.info("  Weighted AUC: 0.971")
    logger.info(f"\nNew scGPT (cleaned):")
    logger.info(f"  Macro AUC: {results_raw['macro_auc']:.3f}")
    logger.info(f"  Weighted AUC: {results_raw['weighted_auc']:.3f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
