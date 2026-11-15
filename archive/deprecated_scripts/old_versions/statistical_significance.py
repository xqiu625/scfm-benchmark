#!/usr/bin/env python3
"""
Statistical Significance Testing for Foundation Model Benchmark

Comprehensive statistical analysis including:
1. Bootstrap confidence intervals for all models
2. Pairwise model comparisons with multiple testing correction
3. Comparison against baselines (PCA, XGBoost)
4. Effect size calculations (Cohen's d)
5. Publication-ready summary tables
"""

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.decomposition import PCA
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, balanced_accuracy_score
from scipy import stats
from pathlib import Path
from datetime import datetime
import logging
import warnings
import json
warnings.filterwarnings('ignore')

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Paths
BASE_DIR = Path("/bigdata/godziklab/shared/Xinru/302006")
EMB_DIR = BASE_DIR / "02_EMBEDDINGS"
RESULTS_DIR = BASE_DIR / "05_RESULTS" / "statistical_significance"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Foundation Models
MODELS = {
    'STATE': EMB_DIR / "state/covid_state.h5ad",
    'UCE': EMB_DIR / "uce/human_platelet_covid_severity_uce_adata.h5ad",
    'TranscriptFormer': EMB_DIR / "transcriptformer/covid_transcriptformer.h5ad",
    'scGPT': EMB_DIR / "scgpt/platelet_scgpt_combined.h5ad",
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
    else:
        raise ValueError("No severity column found")

    return X, severity, adata


def load_expression_data():
    """Load raw expression data for baseline comparisons."""
    path = EMB_DIR / "uce/human_platelet_covid_severity_uce_adata.h5ad"
    adata = sc.read_h5ad(path)

    # Get expression matrix
    if hasattr(adata.X, 'toarray'):
        X = adata.X.toarray()
    else:
        X = np.array(adata.X)

    severity = adata.obs['covid_severity'].astype(str).str.lower().values
    return X, severity


def create_binary_labels(severity):
    """Create binary labels (severe/fatal vs others)."""
    severe_classes = ['severe', 'fatal']
    y = np.array([1 if s in severe_classes else 0 for s in severity])
    mask = np.array([s != 'unknown' for s in severity])
    return y[mask], mask


def create_3class_labels(severity):
    """Create 3-class labels."""
    label_map = {
        'healthy': 0,
        'mild': 1, 'moderate': 1, 'recovered': 1,
        'severe': 2, 'fatal': 2
    }
    y = np.array([label_map.get(s, -1) for s in severity])
    mask = y >= 0
    return y[mask], mask


def bootstrap_auc(X, y, n_bootstrap=1000, n_splits=5, random_state=42):
    """
    Bootstrap confidence intervals for AUC.
    Uses nested CV within each bootstrap sample.
    """
    np.random.seed(random_state)
    n_samples = len(y)
    bootstrap_aucs = []

    for b in range(n_bootstrap):
        # Bootstrap sample
        idx = np.random.choice(n_samples, size=n_samples, replace=True)
        X_boot, y_boot = X[idx], y[idx]

        # Skip if only one class
        if len(np.unique(y_boot)) < 2:
            continue

        # CV within bootstrap
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state + b)
        fold_aucs = []

        for train_idx, test_idx in cv.split(X_boot, y_boot):
            X_train, X_test = X_boot[train_idx], X_boot[test_idx]
            y_train, y_test = y_boot[train_idx], y_boot[test_idx]

            if len(np.unique(y_test)) < 2:
                continue

            scaler = StandardScaler()
            X_train_s = scaler.fit_transform(X_train)
            X_test_s = scaler.transform(X_test)

            clf = LogisticRegression(max_iter=1000, random_state=42)
            clf.fit(X_train_s, y_train)

            y_prob = clf.predict_proba(X_test_s)
            if len(np.unique(y)) == 2:
                auc = roc_auc_score(y_test, y_prob[:, 1])
            else:
                auc = roc_auc_score(y_test, y_prob, multi_class='ovr', average='macro')
            fold_aucs.append(auc)

        if fold_aucs:
            bootstrap_aucs.append(np.mean(fold_aucs))

    return np.array(bootstrap_aucs)


def cv_with_fold_scores(X, y, n_splits=10, n_repeats=5):
    """Get fold-level AUC scores for paired tests."""
    all_aucs = []

    for repeat in range(n_repeats):
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42 + repeat)

        for train_idx, test_idx in cv.split(X, y):
            X_train, X_test = X[train_idx], X[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]

            scaler = StandardScaler()
            X_train_s = scaler.fit_transform(X_train)
            X_test_s = scaler.transform(X_test)

            clf = LogisticRegression(max_iter=1000, random_state=42)
            clf.fit(X_train_s, y_train)

            y_prob = clf.predict_proba(X_test_s)
            if len(np.unique(y)) == 2:
                auc = roc_auc_score(y_test, y_prob[:, 1])
            else:
                auc = roc_auc_score(y_test, y_prob, multi_class='ovr', average='macro')

            all_aucs.append(auc)

    return np.array(all_aucs)


def get_pca_baseline_scores(X_expr, y, n_components=100, n_splits=10, n_repeats=5):
    """Get PCA baseline fold-level scores."""
    pca = PCA(n_components=min(n_components, X_expr.shape[1]))
    X_pca = pca.fit_transform(X_expr)
    return cv_with_fold_scores(X_pca, y, n_splits, n_repeats)


def get_xgboost_baseline_scores(X_expr, y, n_splits=10, n_repeats=3):
    """Get XGBoost baseline fold-level scores.
    Uses PCA to reduce dimensions first for speed.
    """
    all_aucs = []

    # Apply PCA first to reduce dimensions (makes XGBoost much faster)
    logger.info("    Applying PCA(500) for XGBoost baseline...")
    pca = PCA(n_components=min(500, X_expr.shape[1]))
    X_pca = pca.fit_transform(X_expr)

    for repeat in range(n_repeats):
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42 + repeat)

        for train_idx, test_idx in cv.split(X_pca, y):
            X_train, X_test = X_pca[train_idx], X_pca[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]

            # Use fewer estimators for speed
            clf = GradientBoostingClassifier(n_estimators=50, max_depth=4, random_state=42)
            clf.fit(X_train, y_train)

            y_prob = clf.predict_proba(X_test)
            if len(np.unique(y)) == 2:
                auc = roc_auc_score(y_test, y_prob[:, 1])
            else:
                auc = roc_auc_score(y_test, y_prob, multi_class='ovr', average='macro')

            all_aucs.append(auc)

    return np.array(all_aucs)


def pairwise_tests(scores_dict, alpha=0.05):
    """
    Perform pairwise statistical tests with Bonferroni correction.
    """
    models = list(scores_dict.keys())
    n_comparisons = len(models) * (len(models) - 1) // 2
    bonferroni_alpha = alpha / n_comparisons

    results = []

    for i, model1 in enumerate(models):
        for model2 in models[i+1:]:
            scores1 = scores_dict[model1]
            scores2 = scores_dict[model2]

            # Use minimum length
            min_len = min(len(scores1), len(scores2))
            s1, s2 = scores1[:min_len], scores2[:min_len]

            # Paired t-test
            t_stat, t_pval = stats.ttest_rel(s1, s2)

            # Wilcoxon signed-rank test
            try:
                w_stat, w_pval = stats.wilcoxon(s1, s2)
            except:
                w_stat, w_pval = np.nan, np.nan

            # Effect size (Cohen's d for paired samples)
            diff = s1 - s2
            cohens_d = np.mean(diff) / np.std(diff) if np.std(diff) > 0 else 0

            # Determine winner
            mean_diff = np.mean(s1) - np.mean(s2)
            winner = model1 if mean_diff > 0 else model2

            results.append({
                'model1': model1,
                'model2': model2,
                'mean1': np.mean(s1),
                'mean2': np.mean(s2),
                'mean_diff': mean_diff,
                'std_diff': np.std(diff),
                't_statistic': t_stat,
                't_pvalue': t_pval,
                'wilcoxon_pvalue': w_pval,
                'cohens_d': cohens_d,
                'effect_size': 'large' if abs(cohens_d) > 0.8 else ('medium' if abs(cohens_d) > 0.5 else 'small'),
                'significant_raw': t_pval < alpha,
                'significant_bonferroni': t_pval < bonferroni_alpha,
                'winner': winner if t_pval < alpha else 'tie',
                'bonferroni_alpha': bonferroni_alpha
            })

    return pd.DataFrame(results)


def main():
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    logger.info("=" * 80)
    logger.info("STATISTICAL SIGNIFICANCE TESTING")
    logger.info("=" * 80)

    # Load expression data for baselines
    logger.info("\nLoading expression data for baselines...")
    X_expr, severity_expr = load_expression_data()

    tasks = {
        'binary': create_binary_labels,
        '3class': create_3class_labels,
    }

    all_summary = []

    for task_name, label_func in tasks.items():
        logger.info(f"\n{'='*60}")
        logger.info(f"Task: {task_name}")
        logger.info(f"{'='*60}")

        scores_dict = {}
        bootstrap_ci = {}

        # Get baseline scores first
        y_expr, mask_expr = label_func(severity_expr)
        X_expr_masked = X_expr[mask_expr]

        logger.info(f"\nBaselines (n={len(y_expr)} samples):")

        # PCA baseline
        logger.info("  Computing PCA_100 baseline...")
        pca_scores = get_pca_baseline_scores(X_expr_masked, y_expr, n_components=100)
        scores_dict['PCA_100'] = pca_scores
        logger.info(f"    AUC: {np.mean(pca_scores):.4f} ± {np.std(pca_scores):.4f}")

        # XGBoost baseline
        logger.info("  Computing XGBoost baseline...")
        xgb_scores = get_xgboost_baseline_scores(X_expr_masked, y_expr)
        scores_dict['XGBoost'] = xgb_scores
        logger.info(f"    AUC: {np.mean(xgb_scores):.4f} ± {np.std(xgb_scores):.4f}")

        # Foundation models
        logger.info("\nFoundation Models:")
        for model_name in MODELS.keys():
            try:
                logger.info(f"\n  {model_name}:")
                X, severity, _ = load_embeddings(model_name)
                y, mask = label_func(severity)
                X_masked = X[mask]

                logger.info(f"    Samples: {len(y)}")

                # Fold-level scores for paired tests
                fold_scores = cv_with_fold_scores(X_masked, y)
                scores_dict[model_name] = fold_scores

                # Bootstrap CI (reduced to 500 for speed)
                logger.info(f"    Computing bootstrap CI (n=500)...")
                boot_aucs = bootstrap_auc(X_masked, y, n_bootstrap=500)
                ci_low, ci_high = np.percentile(boot_aucs, [2.5, 97.5])
                bootstrap_ci[model_name] = (ci_low, ci_high)

                mean_auc = np.mean(fold_scores)
                std_auc = np.std(fold_scores)

                logger.info(f"    AUC: {mean_auc:.4f} ± {std_auc:.4f}")
                logger.info(f"    95% CI: [{ci_low:.4f}, {ci_high:.4f}]")

                all_summary.append({
                    'task': task_name,
                    'model': model_name,
                    'mean_auc': mean_auc,
                    'std_auc': std_auc,
                    'ci_low': ci_low,
                    'ci_high': ci_high,
                    'n_samples': len(y)
                })

            except Exception as e:
                logger.error(f"    Error: {e}")
                continue

        # Pairwise statistical tests
        logger.info(f"\n--- Pairwise Statistical Tests ({task_name}) ---")
        tests_df = pairwise_tests(scores_dict)
        tests_df.to_csv(RESULTS_DIR / f"pairwise_tests_{task_name}_{timestamp}.csv", index=False)

        # Print significant results
        sig_tests = tests_df[tests_df['significant_bonferroni']].sort_values('mean_diff', ascending=False)
        logger.info(f"\nSignificant comparisons (Bonferroni α={tests_df['bonferroni_alpha'].iloc[0]:.4f}):")
        for _, row in sig_tests.head(10).iterrows():
            logger.info(f"  {row['model1']} vs {row['model2']}: "
                       f"Δ={row['mean_diff']:.4f}, p={row['t_pvalue']:.2e}, "
                       f"d={row['cohens_d']:.2f} ({row['effect_size']}), "
                       f"→ {row['winner']}")

        # Save fold-level scores
        scores_df = pd.DataFrame({model: pd.Series(scores) for model, scores in scores_dict.items()})
        scores_df.to_csv(RESULTS_DIR / f"fold_scores_{task_name}_{timestamp}.csv", index=False)

    # Save summary
    summary_df = pd.DataFrame(all_summary)
    summary_df.to_csv(RESULTS_DIR / f"summary_{timestamp}.csv", index=False)

    # Create publication-ready table
    logger.info("\n" + "=" * 60)
    logger.info("PUBLICATION SUMMARY TABLE")
    logger.info("=" * 60)

    pivot = summary_df.pivot(index='model', columns='task', values='mean_auc')
    pivot_std = summary_df.pivot(index='model', columns='task', values='std_auc')

    for task in pivot.columns:
        logger.info(f"\n{task.upper()} CLASSIFICATION:")
        task_data = summary_df[summary_df['task'] == task].sort_values('mean_auc', ascending=False)
        for _, row in task_data.iterrows():
            logger.info(f"  {row['model']:20s}: {row['mean_auc']:.4f} ± {row['std_auc']:.4f} "
                       f"[{row['ci_low']:.4f}, {row['ci_high']:.4f}]")

    logger.info(f"\n\nResults saved to: {RESULTS_DIR}")
    logger.info("=" * 80)
    logger.info("STATISTICAL SIGNIFICANCE TESTING COMPLETE")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
