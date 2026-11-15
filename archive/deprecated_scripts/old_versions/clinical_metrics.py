#!/usr/bin/env python3
"""
Clinical Metrics Analysis

Computes clinically-relevant metrics beyond AUC:
- Sensitivity (Recall) at various specificity thresholds
- Specificity at various sensitivity thresholds
- Positive Predictive Value (PPV)
- Negative Predictive Value (NPV)
- F1 Score, F2 Score
- AUPRC (Area Under Precision-Recall Curve)
- Calibration metrics (Brier score, Expected Calibration Error)
"""

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import (
    roc_auc_score, average_precision_score, brier_score_loss,
    precision_recall_curve, roc_curve, confusion_matrix,
    f1_score, precision_score, recall_score
)
from sklearn.calibration import calibration_curve
from pathlib import Path
from datetime import datetime
import logging
import warnings
warnings.filterwarnings('ignore')

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Paths
BASE_DIR = Path("/bigdata/godziklab/shared/Xinru/302006")
EMB_DIR = BASE_DIR / "02_EMBEDDINGS"
RESULTS_DIR = BASE_DIR / "05_RESULTS" / "clinical_metrics"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Models
MODELS = {
    'STATE': EMB_DIR / "state/covid_state.h5ad",
    'UCE': EMB_DIR / "uce/human_platelet_covid_severity_uce_adata.h5ad",
    'TranscriptFormer': EMB_DIR / "transcriptformer/covid_transcriptformer.h5ad",
    'scGPT': EMB_DIR / "scgpt/covid_scgpt.h5ad",
    'Geneformer': EMB_DIR / "geneformer/covid_geneformer.h5ad",
}

def load_embeddings(model_name):
    """Load embeddings and metadata."""
    path = MODELS[model_name]
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    logger.info(f"Loading {model_name} from {path}")
    adata = sc.read_h5ad(path)

    # Get embeddings
    if 'X_emb' in adata.obsm:
        X = adata.obsm['X_emb']
    elif 'X_uce' in adata.obsm:
        X = adata.obsm['X_uce']
    elif 'X_scGPT' in adata.obsm:
        X = adata.obsm['X_scGPT']
    else:
        X = adata.X
        if hasattr(X, 'toarray'):
            X = X.toarray()

    return X, adata

def get_severity_labels(adata, task='binary'):
    """Get severity labels."""
    if 'severity' in adata.obs.columns:
        severity = adata.obs['severity'].values
    elif 'Severity' in adata.obs.columns:
        severity = adata.obs['Severity'].values
    else:
        raise ValueError("No severity column found")

    if task == 'binary':
        severe_classes = ['severe', 'fatal', 'Severe', 'Fatal']
        y = np.array([1 if s in severe_classes else 0 for s in severity])
    else:
        le = LabelEncoder()
        y = le.fit_transform(severity)

    return y

def compute_clinical_metrics(y_true, y_prob, y_pred):
    """Compute comprehensive clinical metrics."""
    metrics = {}

    # Basic metrics
    metrics['auc'] = roc_auc_score(y_true, y_prob)
    metrics['auprc'] = average_precision_score(y_true, y_prob)
    metrics['brier_score'] = brier_score_loss(y_true, y_prob)

    # Classification metrics
    metrics['f1'] = f1_score(y_true, y_pred)
    metrics['precision'] = precision_score(y_true, y_pred)
    metrics['recall'] = recall_score(y_true, y_pred)

    # Confusion matrix derived
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    metrics['sensitivity'] = tp / (tp + fn) if (tp + fn) > 0 else 0
    metrics['specificity'] = tn / (tn + fp) if (tn + fp) > 0 else 0
    metrics['ppv'] = tp / (tp + fp) if (tp + fp) > 0 else 0
    metrics['npv'] = tn / (tn + fn) if (tn + fn) > 0 else 0

    # Sensitivity at fixed specificity
    fpr, tpr, thresholds = roc_curve(y_true, y_prob)
    for target_spec in [0.90, 0.95, 0.99]:
        target_fpr = 1 - target_spec
        idx = np.argmin(np.abs(fpr - target_fpr))
        metrics[f'sens_at_spec{int(target_spec*100)}'] = tpr[idx]

    # Specificity at fixed sensitivity
    for target_sens in [0.90, 0.95, 0.99]:
        idx = np.argmin(np.abs(tpr - target_sens))
        metrics[f'spec_at_sens{int(target_sens*100)}'] = 1 - fpr[idx]

    return metrics

def analyze_model(model_name, n_splits=5, n_repeats=3):
    """Analyze a single model with cross-validation."""
    try:
        X, adata = load_embeddings(model_name)
        y = get_severity_labels(adata, task='binary')

        logger.info(f"  Samples: {len(y)}, Positive: {y.sum()}, Negative: {(1-y).sum()}")

        all_metrics = []

        for repeat in range(n_repeats):
            cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42 + repeat)

            for fold, (train_idx, test_idx) in enumerate(cv.split(X, y)):
                X_train, X_test = X[train_idx], X[test_idx]
                y_train, y_test = y[train_idx], y[test_idx]

                # Scale
                scaler = StandardScaler()
                X_train_scaled = scaler.fit_transform(X_train)
                X_test_scaled = scaler.transform(X_test)

                # Train
                clf = LogisticRegression(max_iter=1000, random_state=42)
                clf.fit(X_train_scaled, y_train)

                # Predict
                y_prob = clf.predict_proba(X_test_scaled)[:, 1]
                y_pred = clf.predict(X_test_scaled)

                # Compute metrics
                metrics = compute_clinical_metrics(y_test, y_prob, y_pred)
                metrics['model'] = model_name
                metrics['repeat'] = repeat
                metrics['fold'] = fold
                all_metrics.append(metrics)

        return pd.DataFrame(all_metrics)

    except Exception as e:
        logger.error(f"Error processing {model_name}: {e}")
        return None

def main():
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    logger.info("=" * 80)
    logger.info("CLINICAL METRICS ANALYSIS")
    logger.info("=" * 80)

    all_results = []

    for model_name, path in MODELS.items():
        if not path.exists():
            logger.warning(f"Skipping {model_name}: file not found")
            continue

        logger.info(f"\n{'='*60}")
        logger.info(f"Model: {model_name}")
        logger.info(f"{'='*60}")

        results = analyze_model(model_name)
        if results is not None:
            all_results.append(results)

    if all_results:
        results_df = pd.concat(all_results, ignore_index=True)

        # Save detailed results
        output_file = RESULTS_DIR / f"clinical_metrics_detailed_{timestamp}.csv"
        results_df.to_csv(output_file, index=False)
        logger.info(f"\nDetailed results saved to {output_file}")

        # Summary statistics
        summary = results_df.groupby('model').agg({
            'auc': ['mean', 'std'],
            'auprc': ['mean', 'std'],
            'sensitivity': ['mean', 'std'],
            'specificity': ['mean', 'std'],
            'ppv': ['mean', 'std'],
            'npv': ['mean', 'std'],
            'sens_at_spec90': ['mean', 'std'],
            'sens_at_spec95': ['mean', 'std'],
            'brier_score': ['mean', 'std'],
        }).round(4)

        summary_file = RESULTS_DIR / f"clinical_metrics_summary_{timestamp}.csv"
        summary.to_csv(summary_file)
        logger.info(f"Summary saved to {summary_file}")

        # Print summary
        logger.info("\n" + "=" * 80)
        logger.info("SUMMARY")
        logger.info("=" * 80)

        for model in results_df['model'].unique():
            model_data = results_df[results_df['model'] == model]
            logger.info(f"\n{model}:")
            logger.info(f"  AUC: {model_data['auc'].mean():.4f} +/- {model_data['auc'].std():.4f}")
            logger.info(f"  AUPRC: {model_data['auprc'].mean():.4f} +/- {model_data['auprc'].std():.4f}")
            logger.info(f"  Sensitivity: {model_data['sensitivity'].mean():.4f}")
            logger.info(f"  Specificity: {model_data['specificity'].mean():.4f}")
            logger.info(f"  PPV: {model_data['ppv'].mean():.4f}")
            logger.info(f"  NPV: {model_data['npv'].mean():.4f}")
            logger.info(f"  Sens@Spec90: {model_data['sens_at_spec90'].mean():.4f}")
            logger.info(f"  Sens@Spec95: {model_data['sens_at_spec95'].mean():.4f}")
            logger.info(f"  Brier: {model_data['brier_score'].mean():.4f}")

    logger.info("\n" + "=" * 80)
    logger.info("Clinical Metrics Analysis Complete")
    logger.info("=" * 80)

if __name__ == "__main__":
    main()
