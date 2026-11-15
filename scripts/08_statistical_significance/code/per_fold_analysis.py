#!/usr/bin/env python3
"""
Per-Fold AUC Analysis

Get fold-level AUC values for:
1. Violin plots showing AUC distributions
2. Statistical tests (paired t-tests, Wilcoxon)
3. Confidence intervals
"""

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score, balanced_accuracy_score
from scipy import stats
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
RESULTS_DIR = BASE_DIR / "05_RESULTS" / "per_fold_analysis"
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
    """Load embeddings and labels."""
    path = MODELS[model_name]
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

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

    # Get severity
    if 'covid_severity' in adata.obs.columns:
        severity = adata.obs['covid_severity'].astype(str).str.lower().values
    elif 'severity' in adata.obs.columns:
        severity = adata.obs['severity'].values
    elif 'Severity' in adata.obs.columns:
        severity = adata.obs['Severity'].values
    else:
        raise ValueError("No severity column found")

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

    return X, severity, groups

def create_labels(severity, task='binary'):
    """Create labels for different tasks."""
    if task == 'binary':
        severe_classes = ['severe', 'fatal']
        y = np.array([1 if s in severe_classes else 0 for s in severity])
    elif task == '3-class':
        label_map = {
            'healthy': 0,
            'mild': 1, 'moderate': 1, 'recovered': 1,
            'severe': 2, 'fatal': 2
        }
        y = np.array([label_map.get(s, -1) for s in severity])
        mask = y >= 0
        return y[mask], mask
    elif task == '6-class':
        label_map = {
            'healthy': 0, 'mild': 1, 'moderate': 2,
            'severe': 3, 'fatal': 4, 'recovered': 5
        }
        y = np.array([label_map.get(s, -1) for s in severity])
        mask = y >= 0
        return y[mask], mask
    return y, np.ones(len(y), dtype=bool)

def get_per_fold_metrics(X, y, n_splits=10, n_repeats=5, groups=None):
    """Get per-fold AUC values with multiple repeats."""
    all_fold_aucs = []
    all_fold_bal_accs = []

    for repeat in range(n_repeats):
        cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=42 + repeat)

        for fold_idx, (train_idx, test_idx) in enumerate(cv.split(X, y, groups=groups)):
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
            y_prob = clf.predict_proba(X_test_scaled)
            y_pred = clf.predict(X_test_scaled)

            # Metrics
            if len(np.unique(y)) == 2:
                auc = roc_auc_score(y_test, y_prob[:, 1])
            else:
                auc = roc_auc_score(y_test, y_prob, multi_class='ovr', average='macro')

            bal_acc = balanced_accuracy_score(y_test, y_pred)

            all_fold_aucs.append(auc)
            all_fold_bal_accs.append(bal_acc)

    return np.array(all_fold_aucs), np.array(all_fold_bal_accs)

def pairwise_statistical_tests(fold_aucs_dict):
    """Perform pairwise statistical tests between models."""
    models = list(fold_aucs_dict.keys())
    results = []

    for i, model1 in enumerate(models):
        for model2 in models[i+1:]:
            aucs1 = fold_aucs_dict[model1]
            aucs2 = fold_aucs_dict[model2]

            # Use minimum length
            min_len = min(len(aucs1), len(aucs2))
            aucs1 = aucs1[:min_len]
            aucs2 = aucs2[:min_len]

            # Paired t-test
            t_stat, t_pval = stats.ttest_rel(aucs1, aucs2)

            # Wilcoxon signed-rank test
            try:
                w_stat, w_pval = stats.wilcoxon(aucs1, aucs2)
            except Exception:
                w_stat, w_pval = np.nan, np.nan

            # Effect size (Cohen's d)
            diff = aucs1 - aucs2
            cohens_d = np.mean(diff) / np.std(diff) if np.std(diff) > 0 else 0

            results.append({
                'model1': model1,
                'model2': model2,
                'mean_diff': np.mean(diff),
                'std_diff': np.std(diff),
                't_statistic': t_stat,
                't_pvalue': t_pval,
                'wilcoxon_statistic': w_stat,
                'wilcoxon_pvalue': w_pval,
                'cohens_d': cohens_d,
                'significant_t': t_pval < 0.05,
                'significant_w': w_pval < 0.05 if not np.isnan(w_pval) else False
            })

    return pd.DataFrame(results)

def main():
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    logger.info("=" * 80)
    logger.info("PER-FOLD AUC ANALYSIS")
    logger.info("=" * 80)

    tasks = ['binary', '3-class', '6-class']

    for task in tasks:
        logger.info(f"\n{'='*60}")
        logger.info(f"Task: {task}")
        logger.info(f"{'='*60}")

        fold_aucs_dict = {}
        all_results = []

        for model_name, path in MODELS.items():
            if not path.exists():
                logger.warning(f"Skipping {model_name}: file not found")
                continue

            try:
                logger.info(f"\n{model_name}:")
                X, severity, groups = load_embeddings(model_name)
                y, mask = create_labels(severity, task)

                if not mask.all():
                    X = X[mask]
                    groups = groups[mask]

                logger.info(f"  Samples: {len(y)}, Classes: {np.bincount(y)}")

                fold_aucs, fold_bal_accs = get_per_fold_metrics(X, y, n_splits=10, n_repeats=5, groups=groups)

                fold_aucs_dict[model_name] = fold_aucs

                # Store results
                for i, (auc, bal_acc) in enumerate(zip(fold_aucs, fold_bal_accs)):
                    all_results.append({
                        'model': model_name,
                        'task': task,
                        'fold': i,
                        'auc': auc,
                        'balanced_accuracy': bal_acc
                    })

                logger.info(f"  AUC: {np.mean(fold_aucs):.4f} ± {np.std(fold_aucs):.4f}")
                logger.info(f"  95% CI: [{np.percentile(fold_aucs, 2.5):.4f}, {np.percentile(fold_aucs, 97.5):.4f}]")

            except Exception as e:
                logger.error(f"Error processing {model_name}: {e}")
                continue

        # Save fold-level results
        if all_results:
            results_df = pd.DataFrame(all_results)
            results_df.to_csv(RESULTS_DIR / f"fold_aucs_{task}_{timestamp}.csv", index=False)

        # Statistical tests
        if len(fold_aucs_dict) >= 2:
            logger.info(f"\nStatistical Tests ({task}):")
            tests_df = pairwise_statistical_tests(fold_aucs_dict)
            tests_df.to_csv(RESULTS_DIR / f"statistical_tests_{task}_{timestamp}.csv", index=False)

            # Print significant comparisons
            sig_tests = tests_df[tests_df['significant_t']]
            for _, row in sig_tests.iterrows():
                winner = row['model1'] if row['mean_diff'] > 0 else row['model2']
                logger.info(f"  {row['model1']} vs {row['model2']}: "
                          f"p={row['t_pvalue']:.4f}, diff={row['mean_diff']:.4f} ({winner} wins)")

    logger.info(f"\nResults saved to {RESULTS_DIR}")

if __name__ == "__main__":
    main()
