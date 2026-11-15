#!/usr/bin/env python3
"""
Benchmark scGPT Best Practices Embeddings

Compare scGPT embeddings generated with best practices vs original implementation.
Runs:
1. 6-class severity prediction
2. Binary severity prediction
3. 5-class differential diagnosis
4. Batch-shift generalization (COVID -> Sepsis)

Author: Project 302006
Date: January 2026
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
from datetime import datetime
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score, accuracy_score, f1_score
from sklearn.decomposition import PCA

import warnings
warnings.filterwarnings('ignore')

# Setup logging
import logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_embeddings(file_path: str, embed_key: str = "X_scgpt_bp") -> tuple:
    """Load embeddings and metadata from h5ad file."""
    logger.info(f"Loading {file_path}...")
    adata = sc.read_h5ad(file_path)

    if embed_key not in adata.obsm:
        available_keys = list(adata.obsm.keys())
        raise ValueError(f"Embedding key '{embed_key}' not found. Available: {available_keys}")

    embeddings = adata.obsm[embed_key]
    metadata = adata.obs.copy()

    logger.info(f"  Cells: {embeddings.shape[0]}")
    logger.info(f"  Embedding dim: {embeddings.shape[1]}")

    return embeddings, metadata


def run_classification(
    X: np.ndarray,
    y: np.ndarray,
    n_splits: int = 5,
    task_name: str = "Classification",
    groups: np.ndarray = None,
) -> dict:
    """Run cross-validated classification."""

    # Filter out classes with too few samples for stratified CV
    unique_classes, class_counts = np.unique(y, return_counts=True)
    min_count = n_splits  # need at least 1 sample per fold
    valid_classes = unique_classes[class_counts >= min_count]
    if len(valid_classes) < len(unique_classes):
        removed = set(unique_classes) - set(valid_classes)
        logger.warning(f"  Removing classes with <{min_count} samples: {removed}")
        valid_mask = np.isin(y, valid_classes)
        X = X[valid_mask]
        y = y[valid_mask]
        if groups is not None:
            groups = groups[valid_mask]

    # Encode labels
    le = LabelEncoder()
    y_encoded = le.fit_transform(y)
    n_classes = len(le.classes_)

    if n_classes < 2:
        logger.warning(f"  Skipping {task_name}: fewer than 2 classes after filtering")
        return {'task': task_name, 'n_classes': n_classes, 'auc_mean': float('nan'),
                'auc_std': float('nan'), 'acc_mean': float('nan'), 'f1_mean': float('nan'),
                'classes': list(le.classes_)}

    # Scale features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Cross-validation
    skf = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=42)

    aucs = []
    accs = []
    f1s = []

    for fold, (train_idx, test_idx) in enumerate(skf.split(X_scaled, y_encoded, groups=groups)):
        X_train, X_test = X_scaled[train_idx], X_scaled[test_idx]
        y_train, y_test = y_encoded[train_idx], y_encoded[test_idx]

        # Train model
        clf = LogisticRegression(
            max_iter=1000,
            multi_class='multinomial' if n_classes > 2 else 'auto',
            solver='lbfgs',
            random_state=42
        )
        clf.fit(X_train, y_train)

        # Predict
        y_pred = clf.predict(X_test)
        y_prob = clf.predict_proba(X_test)

        # Metrics
        acc = accuracy_score(y_test, y_pred)
        f1 = f1_score(y_test, y_pred, average='macro')

        try:
            if n_classes == 2:
                auc = roc_auc_score(y_test, y_prob[:, 1])
            else:
                auc = roc_auc_score(y_test, y_prob, multi_class='ovr', average='macro')
        except ValueError:
            # Some folds may still miss a class due to group-level splitting
            logger.warning(f"  Fold {fold}: AUC computation failed (missing class in test set), skipping")
            continue

        aucs.append(auc)
        accs.append(acc)
        f1s.append(f1)

    results = {
        'task': task_name,
        'n_classes': n_classes,
        'classes': list(le.classes_),
        'auc_mean': np.mean(aucs) if aucs else float('nan'),
        'auc_std': np.std(aucs) if aucs else float('nan'),
        'acc_mean': np.mean(accs) if accs else float('nan'),
        'f1_mean': np.mean(f1s) if f1s else float('nan'),
    }

    logger.info(f"  {task_name}: AUC={results['auc_mean']:.3f} (+/- {results['auc_std']:.3f})")

    return results


def run_batch_shift(
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_test: np.ndarray,
    y_test: np.ndarray,
    task_name: str = "Batch Shift"
) -> dict:
    """Run batch-shift generalization test."""

    # Encode labels
    le = LabelEncoder()
    le.fit(np.concatenate([y_train, y_test]))
    y_train_enc = le.transform(y_train)
    y_test_enc = le.transform(y_test)

    # Scale
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # Train on source, test on target
    clf = LogisticRegression(max_iter=1000, solver='lbfgs', random_state=42)
    clf.fit(X_train_scaled, y_train_enc)

    y_pred = clf.predict(X_test_scaled)
    y_prob = clf.predict_proba(X_test_scaled)

    auc = roc_auc_score(y_test_enc, y_prob[:, 1]) if len(le.classes_) == 2 else \
          roc_auc_score(y_test_enc, y_prob, multi_class='ovr', average='macro')

    results = {
        'task': task_name,
        'auc': auc,
        'train_size': len(y_train),
        'test_size': len(y_test)
    }

    logger.info(f"  {task_name}: AUC={auc:.3f}")

    return results


def main():
    """Main benchmark execution."""

    # Paths
    BASE_DIR = Path("/bigdata/godziklab/shared/Xinru/302006")
    EMBED_DIR = BASE_DIR / "02_EMBEDDINGS"
    OUTPUT_DIR = Path(__file__).parent.parent / "scGPT_results"
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Files to compare
    files = {
        'scGPT_original': EMBED_DIR / "scgpt" / "covid_scgpt_clean.h5ad",
        'scGPT_best_practices': EMBED_DIR / "scgpt_best_practices" / "covid_scgpt_best_practices.h5ad"
    }

    embed_keys = {
        'scGPT_original': 'X_scgpt',
        'scGPT_best_practices': 'X_scgpt_bp'
    }

    # Check which files exist
    available_files = {k: v for k, v in files.items() if v.exists()}

    if not available_files:
        logger.error("No embedding files found!")
        logger.error("Expected files:")
        for name, path in files.items():
            logger.error(f"  {name}: {path}")
        sys.exit(1)

    logger.info("=" * 60)
    logger.info("scGPT BEST PRACTICES BENCHMARK")
    logger.info("=" * 60)
    logger.info(f"Date: {datetime.now()}")
    logger.info(f"Available models: {list(available_files.keys())}")

    all_results = []

    for model_name, file_path in available_files.items():
        logger.info("\n" + "=" * 60)
        logger.info(f"MODEL: {model_name}")
        logger.info("=" * 60)

        # Load embeddings
        embed_key = embed_keys.get(model_name, 'X_scgpt')

        # Try multiple possible keys
        adata = sc.read_h5ad(file_path)
        if embed_key not in adata.obsm:
            # Try to find any embedding key
            possible_keys = [k for k in adata.obsm.keys() if 'scgpt' in k.lower() or k == 'X']
            if possible_keys:
                embed_key = possible_keys[0]
            elif adata.obsm:
                embed_key = list(adata.obsm.keys())[0]
            else:
                # Use X matrix directly if it's embeddings
                embeddings = adata.X
                if hasattr(embeddings, 'toarray'):
                    embeddings = embeddings.toarray()
        else:
            embeddings = adata.obsm[embed_key]

        if 'embeddings' not in dir():
            embeddings = adata.obsm[embed_key]

        metadata = adata.obs.copy()

        # Extract patient-level groups for GroupKFold
        groups = None
        for col in ['Sample_ID', 'sample_id', 'Data_NO', 'data_no', 'batch', 'donor']:
            if col in metadata.columns:
                groups = metadata[col].values
                logger.info(f"Using group column: {col} ({len(np.unique(groups))} groups)")
                break
        if groups is None:
            logger.warning("No group column found; StratifiedGroupKFold will fall back to cell-level splitting")

        logger.info(f"Embedding key: {embed_key}")
        logger.info(f"Shape: {embeddings.shape}")

        # Print available columns for debugging
        logger.info(f"Available columns: {list(metadata.columns)[:20]}...")

        # Find the correct column names (based on actual data)
        severity_6class_col = None
        for col in ['covid_severity', 'Severity_6class', 'Severity_Class', 'severity_6class', 'severity']:
            if col in metadata.columns:
                severity_6class_col = col
                break

        severity_binary_col = None
        for col in ['mild_vs_severe', 'Severity_binary', 'Severity_Binary', 'severity_binary', 'binary_severity']:
            if col in metadata.columns:
                severity_binary_col = col
                break

        disease_col = None
        for col in ['Disease', 'disease', 'Disease2', 'disease_state', 'condition']:
            if col in metadata.columns:
                disease_col = col
                break

        logger.info(f"Found columns - 6class: {severity_6class_col}, binary: {severity_binary_col}, disease: {disease_col}")

        # Task 1: 6-class severity
        logger.info("\nTask 1: 6-class Severity Prediction")
        if severity_6class_col:
            mask = metadata[severity_6class_col].notna()
            if mask.sum() > 0:
                result = run_classification(
                    embeddings[mask],
                    metadata.loc[mask, severity_6class_col].values,
                    task_name="6-class Severity",
                    groups=groups[mask] if groups is not None else None,
                )
                result['model'] = model_name
                all_results.append(result)
        else:
            logger.warning("  No 6-class severity column found")

        # Task 2: Binary severity (severe_group vs mild_group only)
        logger.info("\nTask 2: Binary Severity Prediction")
        if severity_binary_col:
            # Only use severe_group and mild_group (exclude intermediate and NaN)
            mask = metadata[severity_binary_col].isin(['severe_group', 'mild_group'])
            logger.info(f"  Binary severity values: {metadata[severity_binary_col].value_counts().to_dict()}")
            if mask.sum() > 0:
                result = run_classification(
                    embeddings[mask],
                    metadata.loc[mask, severity_binary_col].values,
                    task_name="Binary Severity",
                    groups=groups[mask] if groups is not None else None,
                )
                result['model'] = model_name
                all_results.append(result)
        else:
            logger.warning("  No binary severity column found")

        # Task 3: 5-class differential diagnosis
        logger.info("\nTask 3: 5-class Differential Diagnosis")
        if disease_col:
            # Define valid classes: COVID-19, Sepsis, SLE, HC (Healthy), Lung_infect (Pneumonia)
            # Map HC -> Healthy for consistency
            disease_mapping = {'HC': 'Healthy', 'Lung_infect': 'Pneumonia'}
            metadata['Disease_mapped'] = metadata[disease_col].replace(disease_mapping)

            valid_classes = ['COVID-19', 'Sepsis', 'SLE', 'Healthy', 'Pneumonia']
            mask = metadata['Disease_mapped'].isin(valid_classes)
            logger.info(f"  Disease values: {metadata[disease_col].unique().tolist()}")
            logger.info(f"  Samples per class: {metadata.loc[mask, 'Disease_mapped'].value_counts().to_dict()}")

            if mask.sum() > 0:
                result = run_classification(
                    embeddings[mask],
                    metadata.loc[mask, 'Disease_mapped'].values,
                    task_name="5-class Differential",
                    groups=groups[mask] if groups is not None else None,
                )
                result['model'] = model_name
                all_results.append(result)
            else:
                logger.warning(f"  No samples match valid classes")
        else:
            logger.warning("  No disease column found")

        # Task 4: Batch-shift generalization (COVID -> Sepsis)
        logger.info("\nTask 4: Batch-Shift Generalization (COVID -> Sepsis)")
        if disease_col and severity_binary_col:
            # Train on COVID with binary severity (exclude intermediate)
            covid_mask = (metadata[disease_col] == 'COVID-19') & \
                        (metadata[severity_binary_col].isin(['severe_group', 'mild_group']))
            # Test on Sepsis with binary severity
            sepsis_mask = (metadata[disease_col] == 'Sepsis') & \
                         (metadata[severity_binary_col].isin(['severe_group', 'mild_group']))

            logger.info(f"  COVID samples: {covid_mask.sum()}, Sepsis samples: {sepsis_mask.sum()}")

            # Check if both classes are present in test set
            sepsis_classes = metadata.loc[sepsis_mask, severity_binary_col].unique()
            logger.info(f"  Sepsis classes: {sepsis_classes.tolist()}")

            if covid_mask.sum() > 0 and sepsis_mask.sum() > 0 and len(sepsis_classes) >= 2:
                result = run_batch_shift(
                    embeddings[covid_mask],
                    metadata.loc[covid_mask, severity_binary_col].values,
                    embeddings[sepsis_mask],
                    metadata.loc[sepsis_mask, severity_binary_col].values,
                    task_name="COVID->Sepsis"
                )
                result['model'] = model_name
                all_results.append(result)
            else:
                logger.warning(f"  Skipping batch-shift: Sepsis only has {len(sepsis_classes)} class(es)")
        else:
            logger.warning("  Missing required columns for batch-shift test")

        # Clean up
        del adata, embeddings

    # Save results
    logger.info("\n" + "=" * 60)
    logger.info("SAVING RESULTS")
    logger.info("=" * 60)

    results_df = pd.DataFrame(all_results)
    output_file = OUTPUT_DIR / "scgpt_best_practices_benchmark_results.csv"
    results_df.to_csv(output_file, index=False)
    logger.info(f"Saved: {output_file}")

    # Print summary table
    logger.info("\n" + "=" * 60)
    logger.info("SUMMARY")
    logger.info("=" * 60)
    print(results_df.to_string(index=False))

    return 0


if __name__ == "__main__":
    sys.exit(main())
