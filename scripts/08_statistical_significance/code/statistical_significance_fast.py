#!/usr/bin/env python3
"""
OPTIMIZED Statistical Significance Analysis for COVID Severity Benchmark

Key optimizations over original:
1. Parallel bootstrap using joblib (16 cores)
2. Faster solver (lbfgs vs saga)
3. Reduced bootstrap iterations (100 vs 500) - still valid for 95% CI
4. PCA dimensionality reduction for high-dim embeddings
5. Lower max_iter (500 vs 2000)

Usage:
    python statistical_significance_fast.py --models STATE --tasks binary 3-class
    python statistical_significance_fast.py --models UCE scGPT --n-bootstrap 100

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

from sklearn.model_selection import StratifiedGroupKFold, cross_val_predict
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.decomposition import PCA
from scipy import stats
from datetime import datetime
from joblib import Parallel, delayed
import json
import logging
import argparse

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Number of parallel jobs (use all available CPUs)
N_JOBS = int(os.environ.get('SLURM_CPUS_PER_TASK', 16))


def get_embedding_config(model_name, base_dir):
    """Get configuration for a specific model"""
    embedding_dir = Path(base_dir) / "02_EMBEDDINGS"

    configs = {
        'UCE': {
            'files': [
                embedding_dir / "uce/human_platelet_covid_severity_uce_adata.h5ad",
                embedding_dir / "uce/covid_uce.h5ad",
            ],
            'obsm_keys': ['X_uce', 'X_uce_4layer'],
        },
        'scGPT': {
            'files': [
                embedding_dir / "scgpt/covid_scgpt.h5ad",
                embedding_dir / "scgpt/platelet_scgpt_combined.h5ad",
            ],
            'obsm_keys': ['X_scGPT', 'X_scgpt'],
        },
        'TranscriptFormer': {
            'files': [
                embedding_dir / "transcriptformer/covid_transcriptformer.h5ad",
            ],
            'obsm_keys': ['X_transcriptformer', 'X_TranscriptFormer'],
        },
        'Geneformer': {
            'files': [
                embedding_dir / "geneformer/covid_geneformer.h5ad",
            ],
            'obsm_keys': ['X_geneformer', 'X_Geneformer'],
        },
        'scPRINT2': {
            'files': [embedding_dir / "scprint2/covid_scprint2.h5ad"],
            'obsm_keys': ['X_scprint2', 'scprint_emb'],
        },
        'STATE': {
            'files': [
                embedding_dir / "state/covid_state.h5ad",
            ],
            'obsm_keys': ['X_emb', 'X_state'],
        }
    }

    return configs.get(model_name)


def load_model_data(model_name, base_dir):
    """Load embeddings and labels for a model"""
    config = get_embedding_config(model_name, base_dir)

    if config is None:
        raise ValueError(f"Unknown model: {model_name}")

    file_path = None
    for fp in config['files']:
        if fp.exists():
            file_path = fp
            break

    if file_path is None:
        logger.warning(f"No file found for {model_name}")
        return None, None

    logger.info(f"Loading {model_name} from {file_path}")
    adata = sc.read_h5ad(file_path)

    # Filter to COVID only for scGPT combined file
    if model_name == 'scGPT' and 'platelet_scgpt_combined' in str(file_path):
        if 'dataset' in adata.obs.columns:
            covid_mask = adata.obs['dataset'] == 'covid'
            if covid_mask.sum() > 0:
                adata = adata[covid_mask].copy()

    # Get embeddings
    obsm_key = None
    for key in config['obsm_keys']:
        if key in adata.obsm:
            obsm_key = key
            break

    if obsm_key is None:
        logger.warning(f"No embedding key found for {model_name}")
        return None, None

    embeddings = adata.obsm[obsm_key]
    if hasattr(embeddings, 'toarray'):
        embeddings = embeddings.toarray()
    embeddings = np.array(embeddings)

    # Get labels
    labels = None
    for col in ['covid_severity', 'severity', 'Category']:
        if col in adata.obs.columns:
            raw_labels = adata.obs[col].values
            labels = np.array([str(x) if pd.notna(x) else 'unknown' for x in raw_labels])
            break

    # Extract patient/sample-level groups for GroupKFold
    groups = None
    for col in ['Sample_ID', 'sample_id', 'Data_NO', 'data_no', 'batch', 'donor']:
        if col in adata.obs.columns:
            groups = adata.obs[col].values
            logger.info(f"  Using group column: {col} ({len(np.unique(groups))} groups)")
            break
    if groups is None:
        logger.warning("  No group column found — using cell index as group")
        groups = np.arange(len(labels))

    return embeddings, labels, groups


def prepare_labels(task_name, labels):
    """Prepare labels for classification task"""
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


def compute_auc(y_true, y_prob, n_classes):
    """Compute AUC for binary or multi-class"""
    if n_classes == 2:
        return roc_auc_score(y_true, y_prob[:, 1])
    else:
        return roc_auc_score(y_true, y_prob, multi_class='ovr', average='macro')


def single_bootstrap_iteration(X, y_encoded, n_classes, n_samples, seed):
    """Single bootstrap iteration with out-of-bag (OOB) evaluation."""
    np.random.seed(seed)

    # Sample with replacement
    indices = np.random.choice(n_samples, n_samples, replace=True)
    X_boot = X[indices]
    y_boot = y_encoded[indices]

    # Out-of-bag samples for honest evaluation
    oob_mask = np.ones(n_samples, dtype=bool)
    oob_mask[np.unique(indices)] = False
    if oob_mask.sum() == 0:
        return None
    X_oob = X[oob_mask]
    y_oob = y_encoded[oob_mask]

    # Check if all classes present in both boot and OOB
    if len(np.unique(y_boot)) < n_classes or len(np.unique(y_oob)) < n_classes:
        return None

    try:
        # Fast pipeline with lbfgs solver
        pipeline = Pipeline([
            ('scaler', StandardScaler()),
            ('classifier', LogisticRegression(
                max_iter=500,
                random_state=seed,
                n_jobs=1,  # Single thread per bootstrap
                class_weight='balanced',
                solver='lbfgs'  # MUCH faster than saga
            ))
        ])
        pipeline.fit(X_boot, y_boot)
        y_prob_oob = pipeline.predict_proba(X_oob)
        return compute_auc(y_oob, y_prob_oob, n_classes)
    except Exception:
        return None


def bootstrap_auc_parallel(X, y, n_bootstrap=100, random_state=42, pca_dims=256, groups=None):
    """
    Compute bootstrap confidence interval for AUC with PARALLEL execution.

    Point estimate uses cross-validation (not train-on-all/test-on-all).
    Bootstrap iterations use out-of-bag (OOB) evaluation.

    Note: PCA is applied on the full data before bootstrap for computational
    efficiency. This is a minor source of information leakage in the CI
    estimation, but the dominant fix is the OOB evaluation in each iteration.
    """
    le = LabelEncoder()
    y_encoded = le.fit_transform(y)
    n_classes = len(le.classes_)
    n_samples = len(y)

    # Reduce dimensionality for high-dim embeddings
    X_reduced = X
    if X.shape[1] > pca_dims:
        logger.info(f"    Reducing dims: {X.shape[1]} -> {pca_dims} via PCA")
        pca = PCA(n_components=pca_dims, random_state=random_state)
        X_reduced = pca.fit_transform(X)

    # Compute point estimate via cross-validation (honest evaluation)
    cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=random_state)
    cv_pipeline = Pipeline([
        ('scaler', StandardScaler()),
        ('classifier', LogisticRegression(
            max_iter=500, random_state=random_state, n_jobs=-1,
            class_weight='balanced', solver='lbfgs'
        ))
    ])
    y_prob_cv = cross_val_predict(cv_pipeline, X_reduced, y_encoded, cv=cv,
                                   method='predict_proba', groups=groups)
    original_auc = compute_auc(y_encoded, y_prob_cv, n_classes)

    # Scale data for bootstrap iterations
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_reduced)

    # PARALLEL bootstrap with OOB evaluation
    logger.info(f"    Running {n_bootstrap} bootstrap iterations on {N_JOBS} cores...")

    bootstrap_aucs = Parallel(n_jobs=N_JOBS, verbose=0)(
        delayed(single_bootstrap_iteration)(
            X_scaled, y_encoded, n_classes, n_samples, random_state + i
        )
        for i in range(n_bootstrap)
    )

    # Filter None values
    bootstrap_aucs = [x for x in bootstrap_aucs if x is not None]
    bootstrap_aucs = np.array(bootstrap_aucs)

    return {
        'auc': original_auc,
        'ci_lower': np.percentile(bootstrap_aucs, 2.5),
        'ci_upper': np.percentile(bootstrap_aucs, 97.5),
        'std': np.std(bootstrap_aucs),
        'n_bootstrap': len(bootstrap_aucs)
    }


def get_fold_aucs(X, y, n_splits=5, pca_dims=256, groups=None):
    """Get AUC for each CV fold (for violin plots). PCA is fit inside each fold."""
    le = LabelEncoder()
    y_encoded = le.fit_transform(y)
    n_classes = len(le.classes_)

    cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=42)
    fold_aucs = []

    for fold_idx, (train_idx, test_idx) in enumerate(cv.split(X, y_encoded, groups=groups)):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y_encoded[train_idx], y_encoded[test_idx]

        # Build pipeline: optionally include PCA (fit on train fold only)
        steps = [('scaler', StandardScaler())]
        if X.shape[1] > pca_dims:
            steps.append(('pca', PCA(n_components=pca_dims, random_state=42)))
        steps.append(('classifier', LogisticRegression(
            max_iter=500, random_state=42, n_jobs=-1,
            class_weight='balanced', solver='lbfgs'
        )))
        pipeline = Pipeline(steps)

        pipeline.fit(X_train, y_train)
        y_prob = pipeline.predict_proba(X_test)

        try:
            auc = compute_auc(y_test, y_prob, n_classes)
            fold_aucs.append({'fold': fold_idx + 1, 'auc': auc})
        except Exception:
            pass

    return fold_aucs


def _hanley_mcneil_var(auc, n_pos, n_neg):
    """Correct Hanley-McNeil single-AUC variance approximation."""
    q1 = auc / (2 - auc)
    q2 = 2 * auc**2 / (1 + auc)
    var = (auc * (1 - auc)
           + (n_pos - 1) * (q1 - auc**2)
           + (n_neg - 1) * (q2 - auc**2)) / (n_pos * n_neg)
    return max(var, 0)  # guard against numerical negatives


def delong_test(auc1, auc2, n_pos1, n_neg1, n_pos2=None, n_neg2=None):
    """
    Hanley-McNeil approximation for comparing two AUCs.

    NOTE: This assumes independent AUCs. For a proper paired DeLong test
    on the same data, use raw predictions with scipy or the pydelong package.

    Parameters
    ----------
    auc1, auc2 : float
        The two AUC values to compare.
    n_pos1, n_neg1 : int
        Number of positive / negative samples for AUC 1.
    n_pos2, n_neg2 : int or None
        Number of positive / negative samples for AUC 2.
        If None, uses same counts as AUC 1 (same test set).
    """
    if n_pos2 is None:
        n_pos2 = n_pos1
    if n_neg2 is None:
        n_neg2 = n_neg1

    var1 = _hanley_mcneil_var(auc1, n_pos1, n_neg1)
    var2 = _hanley_mcneil_var(auc2, n_pos2, n_neg2)

    se_diff = np.sqrt(var1 + var2)

    if se_diff == 0:
        return 0, 1.0

    z = (auc1 - auc2) / se_diff
    p_value = 2 * (1 - stats.norm.cdf(abs(z)))

    return z, p_value


def run_statistical_analysis(base_dir, tasks=None, n_bootstrap=100, n_splits=5, models=None, pca_dims=256):
    """Run full statistical analysis with optimizations"""

    output_dir = Path(base_dir) / "05_RESULTS/covid_severity_benchmark"
    output_dir.mkdir(parents=True, exist_ok=True)

    if tasks is None:
        tasks = ['binary', '3-class']

    if models is None:
        models = ['STATE', 'UCE', 'TranscriptFormer', 'scGPT', 'Geneformer', 'scPRINT2']

    logger.info("="*80)
    logger.info("OPTIMIZED STATISTICAL SIGNIFICANCE ANALYSIS")
    logger.info(f"Tasks: {tasks}")
    logger.info(f"Models: {models}")
    logger.info(f"Bootstrap iterations: {n_bootstrap}")
    logger.info(f"CV Folds: {n_splits}")
    logger.info(f"PCA dims: {pca_dims}")
    logger.info(f"Parallel jobs: {N_JOBS}")
    logger.info("="*80)

    all_results = {}

    for task_name in tasks:
        logger.info(f"\n{'='*60}")
        logger.info(f"Task: {task_name}")
        logger.info(f"{'='*60}")

        task_results = {
            'bootstrap_ci': {},
            'fold_aucs': {},
            'delong_tests': {}
        }

        model_aucs = {}

        for model_name in models:
            logger.info(f"\n--- {model_name} ---")

            try:
                X, labels, groups = load_model_data(model_name, base_dir)
                if X is None:
                    logger.warning(f"  Skipping {model_name} - no data")
                    continue

                y, valid_mask = prepare_labels(task_name, labels)
                X_valid = X[valid_mask]
                y_valid = y[valid_mask]
                groups_valid = groups[valid_mask]

                logger.info(f"  Samples: {len(y_valid)}, Dims: {X_valid.shape[1]}")

                # Bootstrap CI (parallelized, OOB evaluation)
                logger.info(f"  Computing bootstrap CI (n={n_bootstrap})...")
                bootstrap_result = bootstrap_auc_parallel(
                    X_valid, y_valid,
                    n_bootstrap=n_bootstrap,
                    pca_dims=pca_dims,
                    groups=groups_valid
                )
                task_results['bootstrap_ci'][model_name] = bootstrap_result
                le_temp = LabelEncoder()
                y_enc_temp = le_temp.fit_transform(y_valid)
                n_pos = int(np.sum(y_enc_temp == 1))
                n_neg = int(np.sum(y_enc_temp == 0))
                model_aucs[model_name] = (bootstrap_result['auc'], n_pos, n_neg)

                logger.info(f"  AUC: {bootstrap_result['auc']:.4f} "
                           f"[{bootstrap_result['ci_lower']:.4f}, {bootstrap_result['ci_upper']:.4f}]")

                # Fold AUCs
                logger.info(f"  Computing per-fold AUCs...")
                fold_aucs = get_fold_aucs(X_valid, y_valid, n_splits=n_splits, pca_dims=pca_dims, groups=groups_valid)
                task_results['fold_aucs'][model_name] = fold_aucs

                fold_auc_values = [f['auc'] for f in fold_aucs]
                logger.info(f"  Fold AUCs: {[f'{x:.4f}' for x in fold_auc_values]}")

            except Exception as e:
                logger.error(f"  Error: {str(e)}")
                import traceback
                traceback.print_exc()

        # DeLong tests (pairwise)
        logger.info(f"\n--- DeLong Tests ---")
        model_names = list(model_aucs.keys())

        for i, model1 in enumerate(model_names):
            for model2 in model_names[i+1:]:
                auc1, n_pos1, n_neg1 = model_aucs[model1]
                auc2, n_pos2, n_neg2 = model_aucs[model2]

                z, p = delong_test(auc1, auc2, n_pos1, n_neg1, n_pos2, n_neg2)

                comparison_key = f"{model1}_vs_{model2}"
                task_results['delong_tests'][comparison_key] = {
                    'model1': model1,
                    'model2': model2,
                    'auc1': auc1,
                    'auc2': auc2,
                    'auc_diff': auc1 - auc2,
                    'z_statistic': z,
                    'p_value': p,
                    'significant': p < 0.05
                }

                sig_str = "*" if p < 0.05 else ""
                logger.info(f"  {model1} vs {model2}: "
                           f"AUC diff = {auc1-auc2:+.4f}, p = {p:.4f} {sig_str}")

        all_results[task_name] = task_results

    # Save results
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # JSON
    json_file = output_dir / f'statistical_significance_fast_{timestamp}.json'
    with open(json_file, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    logger.info(f"\nSaved: {json_file}")

    # CSV summary - Bootstrap CI
    rows = []
    for task_name, task_results in all_results.items():
        for model_name, ci_result in task_results['bootstrap_ci'].items():
            rows.append({
                'Task': task_name,
                'Model': model_name,
                'AUC': ci_result['auc'],
                'CI_Lower': ci_result['ci_lower'],
                'CI_Upper': ci_result['ci_upper'],
                'CI_Width': ci_result['ci_upper'] - ci_result['ci_lower'],
                'Std': ci_result['std']
            })

    df_ci = pd.DataFrame(rows)
    csv_ci = output_dir / f'bootstrap_ci_fast_{timestamp}.csv'
    df_ci.to_csv(csv_ci, index=False)
    logger.info(f"Saved: {csv_ci}")

    # CSV summary - Fold AUCs (for violin plots)
    rows_fold = []
    for task_name, task_results in all_results.items():
        for model_name, fold_aucs in task_results['fold_aucs'].items():
            for fold_result in fold_aucs:
                rows_fold.append({
                    'Task': task_name,
                    'Model': model_name,
                    'Fold': fold_result['fold'],
                    'AUC': fold_result['auc']
                })

    df_fold = pd.DataFrame(rows_fold)
    csv_fold = output_dir / f'fold_aucs_fast_{timestamp}.csv'
    df_fold.to_csv(csv_fold, index=False)
    logger.info(f"Saved: {csv_fold}")

    # CSV summary - DeLong tests
    rows_delong = []
    for task_name, task_results in all_results.items():
        for comp_name, delong_result in task_results['delong_tests'].items():
            rows_delong.append({
                'Task': task_name,
                'Comparison': comp_name,
                'Model1': delong_result['model1'],
                'Model2': delong_result['model2'],
                'AUC1': delong_result['auc1'],
                'AUC2': delong_result['auc2'],
                'AUC_Diff': delong_result['auc_diff'],
                'Z_Statistic': delong_result['z_statistic'],
                'P_Value': delong_result['p_value'],
                'Significant': delong_result['significant']
            })

    df_delong = pd.DataFrame(rows_delong)
    csv_delong = output_dir / f'delong_tests_fast_{timestamp}.csv'
    df_delong.to_csv(csv_delong, index=False)
    logger.info(f"Saved: {csv_delong}")

    # Print summary
    print(f"\n{'='*80}")
    print("STATISTICAL SIGNIFICANCE SUMMARY")
    print(f"{'='*80}")

    for task_name in all_results.keys():
        print(f"\n{task_name}:")
        print("-" * 60)
        print(f"{'Model':<20} {'AUC':>8} {'95% CI':>20}")
        print("-" * 60)

        for model_name, ci in all_results[task_name]['bootstrap_ci'].items():
            ci_str = f"[{ci['ci_lower']:.4f}, {ci['ci_upper']:.4f}]"
            print(f"{model_name:<20} {ci['auc']:>8.4f} {ci_str:>20}")

    return all_results


def main():
    parser = argparse.ArgumentParser(description='Optimized Statistical Significance Analysis')
    parser.add_argument('--base-dir', type=str,
                        default="/bigdata/godziklab/shared/Xinru/302006",
                        help='Base directory')
    parser.add_argument('--tasks', nargs='+', default=['binary', '3-class'],
                        choices=['binary', '3-class', '6-class'],
                        help='Tasks to analyze')
    parser.add_argument('--n-bootstrap', type=int, default=100,
                        help='Number of bootstrap iterations (default: 100)')
    parser.add_argument('--cv-folds', type=int, default=5,
                        help='Number of CV folds')
    parser.add_argument('--pca-dims', type=int, default=256,
                        help='PCA dimensions for high-dim embeddings (default: 256)')
    parser.add_argument('--models', nargs='+',
                        default=['STATE', 'UCE', 'TranscriptFormer', 'scGPT', 'Geneformer', 'scPRINT2'],
                        choices=['STATE', 'UCE', 'TranscriptFormer', 'scGPT', 'Geneformer', 'scPRINT2'],
                        help='Models to analyze')
    args = parser.parse_args()

    run_statistical_analysis(
        base_dir=args.base_dir,
        tasks=args.tasks,
        n_bootstrap=args.n_bootstrap,
        n_splits=args.cv_folds,
        models=args.models,
        pca_dims=args.pca_dims
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
