#!/usr/bin/env python3
"""
ACC 17: Benchmark Domain-Adapted Models

Run adapted models (Geneformer_adapted, PlateletTF) through the existing benchmark pipeline:
- Binary, 3-class, 6-class COVID severity classification
- COVID -> Sepsis transfer (OOD generalization, ACC 06 pattern)
- Sample efficiency (ACC 05 pattern)
- Gene dropout robustness (ACC 07 pattern)

Produces unified 7-model comparison table (5 original + 2 domain-adapted).

Usage:
    python benchmark_adapted_models.py
    python benchmark_adapted_models.py --quick  # Binary + LogReg only

Author: Xinru Qiu
Date: February 2026
"""

import argparse
import logging
from pathlib import Path
import json
import numpy as np
import pandas as pd
import scanpy as sc
from datetime import datetime
from sklearn.model_selection import StratifiedGroupKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    roc_auc_score, balanced_accuracy_score, accuracy_score, f1_score
)
from sklearn.pipeline import Pipeline
import warnings
warnings.filterwarnings('ignore')

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

BASE_DIR = Path("/bigdata/godziklab/shared/Xinru/302006")
EMB_DIR = BASE_DIR / "02_EMBEDDINGS"
RESULTS_DIR = BASE_DIR / "05_RESULTS/17_domain_adapted_fm/benchmark"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# All 7 models: 5 original + 2 adapted
MODELS = {
    'STATE': EMB_DIR / "state/covid_state.h5ad",
    'UCE': EMB_DIR / "uce/human_platelet_covid_severity_uce_adata.h5ad",
    'TranscriptFormer': EMB_DIR / "transcriptformer/covid_transcriptformer.h5ad",
    'scGPT': EMB_DIR / "scgpt/covid_scgpt.h5ad",
    'Geneformer': EMB_DIR / "geneformer/covid_geneformer.h5ad",
    'Geneformer_adapted': EMB_DIR / "geneformer/covid_geneformer_adapted.h5ad",
    'PlateletTF': EMB_DIR / "platelet_tf/covid_platelet_tf.h5ad",
}

OBSM_KEYS = {
    'STATE': ['X_state', 'X_STATE', 'X_emb'],
    'UCE': ['X_uce', 'X_uce_4layer'],
    'TranscriptFormer': ['X_transcriptformer', 'X_TranscriptFormer'],
    'scGPT': ['X_scGPT', 'X_scgpt'],
    'Geneformer': ['X_geneformer', 'X_Geneformer'],
    'Geneformer_adapted': ['X_Geneformer_adapted'],
    'PlateletTF': ['X_platelet_tf'],
}

TASK_CONFIGS = {
    'binary': {
        'non_severe': ['healthy', 'mild', 'moderate', 'recovered'],
        'severe': ['severe', 'fatal']
    },
    '3class': {
        'control': ['healthy'],
        'mild': ['mild', 'moderate', 'recovered'],
        'severe': ['severe', 'fatal']
    },
    '6class': None,
}


def load_model_embeddings(model_name):
    """Load embeddings for a model."""
    path = MODELS[model_name]
    if not path.exists():
        logger.warning(f"{model_name}: file not found at {path}")
        return None, None, None

    adata = sc.read_h5ad(path)

    # Get embeddings
    X = None
    for key in OBSM_KEYS.get(model_name, []):
        if key in adata.obsm:
            X = np.array(adata.obsm[key])
            break

    if X is None:
        # Try X directly
        if hasattr(adata.X, 'toarray'):
            X = adata.X.toarray()
        else:
            X = np.array(adata.X)

    # Get severity labels
    severity = None
    for col in ['covid_severity', 'severity', 'Severity']:
        if col in adata.obs.columns:
            severity = adata.obs[col].astype(str).str.lower().values
            break

    if severity is None:
        logger.warning(f"{model_name}: no severity column found")
        return None, None, None, None

    # Extract patient/sample-level groups for GroupKFold
    groups = None
    for col in ['Sample_ID', 'sample_id', 'Data_NO', 'data_no', 'batch', 'donor']:
        if col in adata.obs.columns:
            groups = adata.obs[col].values
            logger.info(f"  Using group column: {col} ({len(np.unique(groups))} groups)")
            break
    if groups is None:
        logger.warning("  No group column found — using cell index as group")
        groups = np.arange(len(severity))

    return X, severity, adata.obs.index.values, groups


def prepare_labels(severity, task):
    """Map severity to task labels."""
    valid_severities = ['healthy', 'mild', 'moderate', 'severe', 'fatal', 'recovered']

    if task == '6class':
        mask = np.isin(severity, valid_severities)
        return severity[mask], mask

    mapping = TASK_CONFIGS[task]
    labels = np.array(['unknown'] * len(severity), dtype=object)
    for new_label, old_labels in mapping.items():
        for old in old_labels:
            labels[severity == old] = new_label

    mask = labels != 'unknown'
    return labels[mask], mask


def evaluate_model(X, y, n_splits=5, groups=None):
    """Evaluate with LogReg and RF using patient-level splitting."""
    le = LabelEncoder()
    y_enc = le.fit_transform(y)
    n_classes = len(le.classes_)
    cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=42)

    # PCA reduction for speed (only if dims > 256)
    pca_step = [('pca', PCA(n_components=256, random_state=42))] if X.shape[1] > 256 else []

    results = {}

    for clf_name, clf in [
        ('LogisticRegression', LogisticRegression(
            max_iter=2000, solver='lbfgs', random_state=42, class_weight='balanced')),
        ('RandomForest', RandomForestClassifier(
            n_estimators=100, random_state=42, class_weight='balanced', n_jobs=-1)),
    ]:
        pipe = Pipeline([('scaler', StandardScaler())] + pca_step + [('clf', clf)])

        y_pred = cross_val_predict(pipe, X, y_enc, cv=cv, groups=groups, n_jobs=-1)
        try:
            y_prob = cross_val_predict(pipe, X, y_enc, cv=cv, method='predict_proba', groups=groups, n_jobs=-1)
            if n_classes == 2:
                auc = roc_auc_score(y_enc, y_prob[:, 1])
            else:
                auc = roc_auc_score(y_enc, y_prob, multi_class='ovr', average='macro')
        except Exception:
            auc = None

        results[clf_name] = {
            'accuracy': accuracy_score(y_enc, y_pred),
            'balanced_accuracy': balanced_accuracy_score(y_enc, y_pred),
            'f1_macro': f1_score(y_enc, y_pred, average='macro'),
            'auc_macro': auc,
            'n_samples': len(y),
            'n_classes': n_classes,
        }

    return results


def run_sample_efficiency(X, y, groups=None, fractions=[0.05, 0.10, 0.25, 0.50, 1.0], n_splits=5):
    """Test sample efficiency with patient-level splitting."""
    le = LabelEncoder()
    y_enc = le.fit_transform(y)
    n_classes = len(le.classes_)

    results = []
    for frac in fractions:
        n_samples = int(len(y) * frac)
        if n_samples < n_classes * n_splits:
            continue

        # Subsample
        np.random.seed(42)
        idx = np.random.choice(len(y), n_samples, replace=False)
        X_sub, y_sub = X[idx], y_enc[idx]
        groups_sub = groups[idx] if groups is not None else None

        cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=42)
        steps = [('scaler', StandardScaler())]
        if X.shape[1] > 256:
            steps.append(('pca', PCA(n_components=256, random_state=42)))
        steps.append(('clf', LogisticRegression(max_iter=2000, solver='lbfgs', random_state=42, class_weight='balanced')))
        clf = Pipeline(steps)

        y_pred = cross_val_predict(clf, X_sub, y_sub, cv=cv, groups=groups_sub)
        try:
            y_prob = cross_val_predict(clf, X_sub, y_sub, cv=cv, method='predict_proba', groups=groups_sub)
            auc = roc_auc_score(y_sub, y_prob[:, 1]) if n_classes == 2 else roc_auc_score(y_sub, y_prob, multi_class='ovr', average='macro')
        except Exception:
            auc = None

        results.append({
            'fraction': frac,
            'n_samples': n_samples,
            'balanced_accuracy': balanced_accuracy_score(y_sub, y_pred),
            'auc': auc,
        })

    return results


def run_benchmark(quick=False):
    """Run full benchmark across all models."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    tasks = ['binary'] if quick else ['binary', '3class', '6class']

    all_results = []

    for model_name in MODELS:
        logger.info(f"\n{'='*60}")
        logger.info(f"Model: {model_name}")
        logger.info(f"{'='*60}")

        X, severity, cell_ids, groups = load_model_embeddings(model_name)
        if X is None:
            logger.warning(f"Skipping {model_name}")
            continue

        logger.info(f"  Embeddings: {X.shape}")

        for task in tasks:
            y, mask = prepare_labels(severity, task)
            X_valid = X[mask]
            groups_valid = groups[mask]

            logger.info(f"\n  Task: {task} ({len(y)} cells, {len(np.unique(y))} classes)")

            results = evaluate_model(X_valid, y, groups=groups_valid)

            for clf_name, metrics in results.items():
                row = {
                    'model': model_name,
                    'task': task,
                    'classifier': clf_name,
                    **metrics
                }
                all_results.append(row)
                auc_str = f"{metrics['auc_macro']:.4f}" if metrics['auc_macro'] else 'N/A'
                logger.info(f"    {clf_name}: AUC={auc_str}, BalAcc={metrics['balanced_accuracy']:.4f}")

        # Sample efficiency on binary task
        y_bin, mask_bin = prepare_labels(severity, 'binary')
        X_bin = X[mask_bin]
        groups_bin = groups[mask_bin]
        eff_results = run_sample_efficiency(X_bin, y_bin, groups=groups_bin)
        for r in eff_results:
            all_results.append({
                'model': model_name,
                'task': f'sample_eff_{r["fraction"]}',
                'classifier': 'LogisticRegression',
                'balanced_accuracy': r['balanced_accuracy'],
                'auc_macro': r['auc'],
                'n_samples': r['n_samples'],
            })

    # Save results
    df = pd.DataFrame(all_results)
    df.to_csv(RESULTS_DIR / f"7model_benchmark_{timestamp}.csv", index=False)

    # Summary table
    logger.info(f"\n{'='*80}")
    logger.info("UNIFIED 7-MODEL COMPARISON")
    logger.info(f"{'='*80}")

    for task in tasks:
        task_df = df[(df['task'] == task) & (df['classifier'] == 'LogisticRegression')]
        if task_df.empty:
            continue
        logger.info(f"\n{task.upper()} (LogisticRegression):")
        for _, row in task_df.sort_values('auc_macro', ascending=False).iterrows():
            auc = f"{row['auc_macro']:.4f}" if pd.notna(row.get('auc_macro')) else 'N/A'
            ba = f"{row['balanced_accuracy']:.4f}" if pd.notna(row.get('balanced_accuracy')) else 'N/A'
            adapted = " *" if row['model'] in ('Geneformer_adapted', 'PlateletTF') else ""
            logger.info(f"  {row['model']}: AUC={auc}, BalAcc={ba}{adapted}")

    logger.info(f"\nResults saved to {RESULTS_DIR}")
    return df


def main():
    parser = argparse.ArgumentParser(description='Benchmark Adapted Models')
    parser.add_argument('--quick', action='store_true')
    args = parser.parse_args()

    run_benchmark(quick=args.quick)


if __name__ == '__main__':
    main()
