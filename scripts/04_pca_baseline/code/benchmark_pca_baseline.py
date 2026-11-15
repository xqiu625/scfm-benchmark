#!/usr/bin/env python3
"""
PCA Baseline Benchmark for COVID Severity Prediction

Compares foundation model embeddings against traditional PCA baselines:
1. HVG → PCA (50, 100, 200, 500 dims) → LogisticRegression
2. HVG → Raw expression → LogisticRegression
3. HVG → Raw expression → RandomForest
4. HVG → Raw expression → XGBoost

This establishes whether foundation models provide value over simple baselines.

Usage:
    python benchmark_pca_baseline.py --tasks binary 3-class 6-class
    python benchmark_pca_baseline.py --tasks binary --quick

Author: Project 302006
Date: December 2025
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

from sklearn.model_selection import StratifiedGroupKFold, cross_val_predict
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score, f1_score,
    roc_auc_score, confusion_matrix
)
from sklearn.pipeline import Pipeline
from datetime import datetime
import json
import logging
import argparse

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_raw_expression_data(base_dir):
    """Load raw expression data for COVID severity prediction"""

    # Try multiple possible locations for COVID data with raw expression
    possible_files = [
        Path(base_dir) / "02_EMBEDDINGS/uce/human_platelet_covid_severity_uce_adata.h5ad",
        Path(base_dir) / "01_DATA/raw/human_platelet_covid_severity.h5ad",
        Path(base_dir) / "01_DATA/covid_severity.h5ad",
    ]

    file_path = None
    for fp in possible_files:
        if fp.exists():
            file_path = fp
            break

    if file_path is None:
        raise FileNotFoundError(f"No COVID data file found. Tried: {possible_files}")

    logger.info(f"Loading data from {file_path}")
    adata = sc.read_h5ad(file_path)

    # Get expression matrix
    X = adata.X
    if hasattr(X, 'toarray'):
        X = X.toarray()
    X = np.array(X, dtype=np.float32)

    # Get labels
    labels = None
    label_cols = ['covid_severity', 'severity', 'Category']
    for col in label_cols:
        if col in adata.obs.columns:
            raw_labels = adata.obs[col].values
            labels = np.array([str(x) if pd.notna(x) else 'unknown' for x in raw_labels])
            logger.info(f"  Using label column: {col}")
            break

    if labels is None:
        raise ValueError(f"No severity label column found. Available: {list(adata.obs.columns)}")

    # Extract patient/sample-level groups for GroupKFold
    groups = None
    for col in ['Sample_ID', 'sample_id', 'Data_NO', 'data_no', 'batch', 'donor']:
        if col in adata.obs.columns:
            groups = adata.obs[col].values
            logger.info(f"  Using group column: {col} ({len(np.unique(groups))} groups)")
            break
    if groups is None:
        logger.warning("  No group column found — using cell index as group (no patient-level splitting)")
        groups = np.arange(len(labels))

    logger.info(f"  Shape: {X.shape}")
    logger.info(f"  Label distribution: {pd.Series(labels).value_counts().to_dict()}")

    return X, labels, adata.var_names.tolist(), groups


def prepare_labels(task_name, labels):
    """Prepare labels for a specific classification task"""
    valid_severities = ['healthy', 'mild', 'moderate', 'severe', 'fatal', 'recovered']

    tasks = {
        '6-class': None,
        '3-class': {
            'control': ['healthy'],
            'mild': ['mild', 'moderate', 'recovered'],
            'severe': ['severe', 'fatal']
        },
        'binary': {
            'non_severe': ['healthy', 'mild', 'moderate', 'recovered'],
            'severe': ['severe', 'fatal']
        }
    }

    if task_name == '6-class':
        valid_mask = np.isin(labels, valid_severities)
        return labels, valid_mask

    task_mapping = tasks[task_name]
    new_labels = np.array(['unknown'] * len(labels), dtype=object)

    for new_label, old_labels in task_mapping.items():
        for old_label in old_labels:
            mask = labels == old_label
            new_labels[mask] = new_label

    valid_mask = new_labels != 'unknown'
    return new_labels, valid_mask


def evaluate_baseline(X, y, baseline_name, n_splits=5, groups=None):
    """Evaluate a baseline method using cross-validation with patient-level splitting"""
    le = LabelEncoder()
    y_encoded = le.fit_transform(y)
    n_classes = len(le.classes_)

    # Define baseline pipelines
    if baseline_name.startswith('PCA_'):
        n_components = int(baseline_name.split('_')[1])
        pipeline = Pipeline([
            ('scaler', StandardScaler()),
            ('pca', PCA(n_components=min(n_components, X.shape[1]), random_state=42)),
            ('classifier', LogisticRegression(
                max_iter=2000, random_state=42, n_jobs=-1,
                class_weight='balanced', solver='saga'
            ))
        ])
    elif baseline_name == 'Raw_LogReg':
        pipeline = Pipeline([
            ('scaler', StandardScaler()),
            ('classifier', LogisticRegression(
                max_iter=2000, random_state=42, n_jobs=-1,
                class_weight='balanced', solver='saga'
            ))
        ])
    elif baseline_name == 'Raw_RF':
        pipeline = Pipeline([
            ('scaler', StandardScaler()),
            ('classifier', RandomForestClassifier(
                n_estimators=100, random_state=42, n_jobs=-1,
                class_weight='balanced'
            ))
        ])
    elif baseline_name == 'Raw_XGBoost':
        try:
            from xgboost import XGBClassifier
            pipeline = Pipeline([
                ('scaler', StandardScaler()),
                ('classifier', XGBClassifier(
                    n_estimators=100, random_state=42, n_jobs=-1,
                    use_label_encoder=False, eval_metric='mlogloss'
                ))
            ])
        except ImportError:
            logger.warning("XGBoost not installed, using GradientBoosting instead")
            pipeline = Pipeline([
                ('scaler', StandardScaler()),
                ('classifier', GradientBoostingClassifier(
                    n_estimators=100, random_state=42
                ))
            ])
    else:
        raise ValueError(f"Unknown baseline: {baseline_name}")

    cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=42)

    # Get predictions
    y_pred = cross_val_predict(pipeline, X, y_encoded, cv=cv, groups=groups)

    # Get probabilities for AUC
    try:
        y_prob = cross_val_predict(pipeline, X, y_encoded, cv=cv, method='predict_proba', groups=groups)
        if n_classes == 2:
            auc_macro = roc_auc_score(y_encoded, y_prob[:, 1])
            auc_weighted = auc_macro
        else:
            auc_macro = roc_auc_score(y_encoded, y_prob, multi_class='ovr', average='macro')
            auc_weighted = roc_auc_score(y_encoded, y_prob, multi_class='ovr', average='weighted')
    except Exception as e:
        logger.warning(f"  Could not compute AUC: {e}")
        auc_macro = None
        auc_weighted = None

    metrics = {
        'accuracy': accuracy_score(y_encoded, y_pred),
        'balanced_accuracy': balanced_accuracy_score(y_encoded, y_pred),
        'f1_macro': f1_score(y_encoded, y_pred, average='macro'),
        'f1_weighted': f1_score(y_encoded, y_pred, average='weighted'),
        'auc_macro': auc_macro,
        'auc_weighted': auc_weighted,
        'n_samples': len(y),
        'n_classes': n_classes,
        'classes': le.classes_.tolist(),
        'confusion_matrix': confusion_matrix(y_encoded, y_pred).tolist()
    }

    return metrics


def run_baseline_benchmark(base_dir, tasks=None, quick=False, n_splits=5):
    """Run baseline benchmark for all specified tasks"""

    output_dir = Path(base_dir) / "05_RESULTS/covid_severity_benchmark"
    output_dir.mkdir(parents=True, exist_ok=True)

    if tasks is None:
        tasks = ['binary', '3-class', '6-class']

    # Load raw expression data
    X, labels, gene_names, groups = load_raw_expression_data(base_dir)

    # Define baselines
    if quick:
        baselines = ['PCA_50', 'PCA_100', 'Raw_LogReg']
    else:
        baselines = [
            'PCA_50', 'PCA_100', 'PCA_200', 'PCA_500',
            'Raw_LogReg', 'Raw_RF', 'Raw_XGBoost'
        ]

    logger.info("="*80)
    logger.info("PCA BASELINE BENCHMARK")
    logger.info(f"Quick mode: {quick}")
    logger.info(f"Tasks: {tasks}")
    logger.info(f"Baselines: {baselines}")
    logger.info(f"CV Folds: {n_splits}")
    logger.info(f"Expression matrix: {X.shape}")
    logger.info("="*80)

    all_results = {}

    for task_name in tasks:
        logger.info(f"\n{'='*60}")
        logger.info(f"Task: {task_name}")
        logger.info(f"{'='*60}")

        y, valid_mask = prepare_labels(task_name, labels)
        X_valid = X[valid_mask]
        y_valid = y[valid_mask]
        groups_valid = groups[valid_mask]

        logger.info(f"Samples: {len(y_valid)}, Classes: {np.unique(y_valid)}")

        all_results[task_name] = {}

        for baseline_name in baselines:
            logger.info(f"\n--- {baseline_name} ---")

            try:
                metrics = evaluate_baseline(X_valid, y_valid, baseline_name, n_splits=n_splits, groups=groups_valid)
                all_results[task_name][baseline_name] = metrics

                auc_str = f"{metrics['auc_macro']:.4f}" if metrics['auc_macro'] is not None else 'N/A'
                logger.info(f"  Accuracy: {metrics['accuracy']:.4f}")
                logger.info(f"  Balanced Acc: {metrics['balanced_accuracy']:.4f}")
                logger.info(f"  AUC Macro: {auc_str}")

            except Exception as e:
                logger.error(f"  Failed: {str(e)}")
                all_results[task_name][baseline_name] = {'error': str(e)}

    # Save results
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # JSON results
    json_file = output_dir / f'benchmark_PCA_baseline_{timestamp}.json'
    with open(json_file, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    logger.info(f"\nSaved: {json_file}")

    # CSV summary
    rows = []
    for task_name, task_results in all_results.items():
        for baseline_name, metrics in task_results.items():
            if 'error' not in metrics:
                rows.append({
                    'Task': task_name,
                    'Embedding': baseline_name,
                    'Classifier': 'Built-in',
                    'Accuracy': metrics['accuracy'],
                    'Balanced_Accuracy': metrics['balanced_accuracy'],
                    'F1_Macro': metrics['f1_macro'],
                    'F1_Weighted': metrics['f1_weighted'],
                    'AUC_Macro': metrics['auc_macro'],
                    'AUC_Weighted': metrics['auc_weighted'],
                    'N_Samples': metrics['n_samples'],
                    'N_Classes': metrics['n_classes']
                })

    df = pd.DataFrame(rows)
    csv_file = output_dir / f'benchmark_PCA_baseline_{timestamp}.csv'
    df.to_csv(csv_file, index=False)
    logger.info(f"Saved: {csv_file}")

    # Print summary
    print(f"\n{'='*80}")
    print("BASELINE BENCHMARK RESULTS SUMMARY")
    print(f"{'='*80}")

    for task in df['Task'].unique():
        task_df = df[df['Task'] == task].copy()
        task_df = task_df.sort_values('AUC_Macro', ascending=False)

        print(f"\n{task}:")
        print("-" * 60)
        for _, row in task_df.iterrows():
            auc_str = f"{row['AUC_Macro']:.4f}" if pd.notna(row['AUC_Macro']) else 'N/A'
            print(f"  {row['Embedding']:15s} | BalAcc: {row['Balanced_Accuracy']:.4f} | AUC: {auc_str}")

    print(f"\n{'='*80}")
    print("Compare with Foundation Models:")
    print("  STATE:            Binary AUC=0.894, 3-class AUC=0.893")
    print("  UCE:              Binary AUC=0.876, 3-class AUC=0.885, 6-class AUC=0.888")
    print("  TranscriptFormer: Binary AUC=0.838, 3-class AUC=0.851, 6-class AUC=0.849")
    print("  scGPT:            Binary AUC=0.833, 3-class AUC=0.747, 6-class AUC=0.810")
    print("  Geneformer:       Binary AUC=0.824, 3-class AUC=0.833")
    print(f"{'='*80}")

    return all_results


def main():
    parser = argparse.ArgumentParser(description='PCA Baseline Benchmark for COVID Severity')
    parser.add_argument('--base-dir', type=str,
                        default="/bigdata/godziklab/shared/Xinru/302006",
                        help='Base directory for the project')
    parser.add_argument('--tasks', nargs='+', default=['binary', '3-class', '6-class'],
                        choices=['binary', '3-class', '6-class'],
                        help='Classification tasks to run')
    parser.add_argument('--quick', action='store_true',
                        help='Quick mode: only PCA_50, PCA_100, Raw_LogReg')
    parser.add_argument('--cv-folds', type=int, default=5,
                        help='Number of cross-validation folds')
    args = parser.parse_args()

    run_baseline_benchmark(
        base_dir=args.base_dir,
        tasks=args.tasks,
        quick=args.quick,
        n_splits=args.cv_folds
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
