#!/usr/bin/env python3
"""
Gene Programs Analysis using NMF on Foundation Model Embeddings

Learn latent gene programs from foundation model embeddings and analyze:
1. Disease severity as program composition shifts
2. Which programs distinguish severe from healthy
3. Perturbation effects in program space
"""

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.decomposition import NMF
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_predict, StratifiedGroupKFold
from sklearn.metrics import roc_auc_score, balanced_accuracy_score
import matplotlib.pyplot as plt
import seaborn as sns
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
RESULTS_DIR = BASE_DIR / "05_RESULTS" / "gene_programs"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Models and their embedding files
MODELS = {
    'STATE': EMB_DIR / "state/covid_state.h5ad",
    'UCE': EMB_DIR / "uce/human_platelet_covid_severity_uce_adata.h5ad",
}

# Perturbation data
PERT_FILES = {
    'STATE': {
        'HEK293T': EMB_DIR / "state/hek293t_state.h5ad",
        'HCT116': EMB_DIR / "state/hct116_state.h5ad",
    },
    'UCE': {
        'HEK293T': EMB_DIR / "uce/hek293t_uce.h5ad",
        'HCT116': EMB_DIR / "uce/hct116_uce.h5ad",
    }
}

def load_covid_embeddings(model_name):
    """Load COVID severity embeddings."""
    path = MODELS[model_name]
    logger.info(f"Loading {model_name} from {path}")
    adata = sc.read_h5ad(path)

    # Get embeddings
    if 'X_emb' in adata.obsm:
        X = adata.obsm['X_emb']
    elif 'X_uce' in adata.obsm:
        X = adata.obsm['X_uce']
    else:
        X = adata.X

    # Get severity labels
    if 'covid_severity' in adata.obs.columns:
        labels = adata.obs['covid_severity'].astype(str).str.lower().values
    elif 'severity' in adata.obs.columns:
        labels = adata.obs['severity'].values
    elif 'Severity' in adata.obs.columns:
        labels = adata.obs['Severity'].values
    else:
        raise ValueError("No severity column found")

    logger.info(f"  Shape: {X.shape}, Labels: {np.unique(labels)}")

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

    return X, labels, adata, groups

def learn_gene_programs(X, n_programs=50, random_state=42):
    """Learn gene programs using NMF."""
    logger.info(f"Learning {n_programs} gene programs with NMF...")

    # NMF requires non-negative input
    X_pos = np.abs(X)

    nmf = NMF(n_components=n_programs, init='nndsvda', random_state=random_state, max_iter=500)
    W = nmf.fit_transform(X_pos)  # cells x programs
    H = nmf.components_           # programs x embedding_dims

    reconstruction_error = nmf.reconstruction_err_
    logger.info(f"  Reconstruction error: {reconstruction_error:.4f}")
    logger.info(f"  Program loadings shape: {W.shape}")

    return W, H, nmf

def analyze_severity_in_program_space(W, labels):
    """Analyze how severity maps to program space."""
    logger.info("Analyzing severity in program space...")

    # Compute centroids
    unique_labels = np.unique(labels)
    centroids = {}
    for label in unique_labels:
        mask = labels == label
        centroids[label] = W[mask].mean(axis=0)

    # Define binary labels
    severe_labels = ['severe', 'fatal']
    non_severe_labels = ['healthy', 'mild', 'moderate', 'recovered']

    binary_labels = np.array(['severe' if l in severe_labels else 'non_severe' for l in labels])

    # Disease and recovery directions
    if 'severe' in centroids and 'healthy' in centroids:
        disease_direction = centroids['severe'] - centroids['healthy']
        disease_direction = disease_direction / np.linalg.norm(disease_direction)
    else:
        disease_direction = None

    if 'recovered' in centroids and 'severe' in centroids:
        recovery_direction = centroids['recovered'] - centroids['severe']
        recovery_direction = recovery_direction / np.linalg.norm(recovery_direction)
    else:
        recovery_direction = None

    return centroids, disease_direction, recovery_direction, binary_labels

def classify_with_programs(W, labels, n_splits=5, groups=None):
    """Classify severity using program loadings."""
    logger.info("Classifying severity using program loadings...")

    # Binary classification
    severe_labels = ['severe', 'fatal']
    binary_labels = np.array([1 if l in severe_labels else 0 for l in labels])

    cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=42)
    clf = LogisticRegression(max_iter=1000, random_state=42)

    # Scale features
    scaler = StandardScaler()
    W_scaled = scaler.fit_transform(W)

    # Cross-validation predictions
    y_prob = cross_val_predict(clf, W_scaled, binary_labels, cv=cv, groups=groups, method='predict_proba')
    y_pred = cross_val_predict(clf, W_scaled, binary_labels, cv=cv, groups=groups)

    auc = roc_auc_score(binary_labels, y_prob[:, 1])
    bal_acc = balanced_accuracy_score(binary_labels, y_pred)

    logger.info(f"  Program-based classification: AUC={auc:.4f}, BalAcc={bal_acc:.4f}")

    return auc, bal_acc

def identify_disease_programs(W, labels, n_top=10):
    """Identify which programs are most associated with disease severity."""
    logger.info("Identifying disease-associated programs...")

    severe_labels = ['severe', 'fatal']
    binary_labels = np.array([1 if l in severe_labels else 0 for l in labels])

    # Fit logistic regression
    scaler = StandardScaler()
    W_scaled = scaler.fit_transform(W)

    clf = LogisticRegression(max_iter=1000, random_state=42)
    clf.fit(W_scaled, binary_labels)

    # Get coefficients
    coefs = clf.coef_[0]

    # Rank programs by absolute coefficient
    program_importance = pd.DataFrame({
        'program': range(len(coefs)),
        'coefficient': coefs,
        'abs_coef': np.abs(coefs),
        'direction': ['SEVERE' if c > 0 else 'NON-SEVERE' for c in coefs]
    }).sort_values('abs_coef', ascending=False)

    logger.info(f"  Top {n_top} disease-associated programs:")
    for i, row in program_importance.head(n_top).iterrows():
        logger.info(f"    Program {row['program']}: coef={row['coefficient']:.4f} ({row['direction']})")

    return program_importance

def compare_n_programs(X, labels, program_counts=[10, 20, 30, 50, 100], groups=None):
    """Compare classification performance with different numbers of programs."""
    logger.info("Comparing different numbers of programs...")

    results = []
    for n_prog in program_counts:
        W, H, nmf = learn_gene_programs(X, n_programs=n_prog)
        auc, bal_acc = classify_with_programs(W, labels, groups=groups)
        results.append({
            'n_programs': n_prog,
            'auc': auc,
            'balanced_accuracy': bal_acc,
            'reconstruction_error': nmf.reconstruction_err_
        })
        logger.info(f"  n_programs={n_prog}: AUC={auc:.4f}")

    return pd.DataFrame(results)

def main():
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    logger.info("=" * 80)
    logger.info("GENE PROGRAMS ANALYSIS (NMF on Foundation Model Embeddings)")
    logger.info("=" * 80)

    all_results = []

    for model_name in MODELS:
        logger.info(f"\n{'='*60}")
        logger.info(f"Model: {model_name}")
        logger.info(f"{'='*60}")

        try:
            # Load embeddings
            X, labels, adata, groups = load_covid_embeddings(model_name)

            # Compare different numbers of programs
            comparison_df = compare_n_programs(X, labels, program_counts=[10, 20, 30, 50, 100], groups=groups)
            comparison_df['model'] = model_name
            all_results.append(comparison_df)

            # Use best number of programs (50 as default)
            W, H, nmf = learn_gene_programs(X, n_programs=50)

            # Analyze severity
            centroids, disease_dir, recovery_dir, binary_labels = analyze_severity_in_program_space(W, labels)

            # Identify disease programs
            program_importance = identify_disease_programs(W, labels)
            program_importance.to_csv(
                RESULTS_DIR / f"disease_programs_{model_name}_{timestamp}.csv",
                index=False
            )

            # Save program loadings
            np.save(RESULTS_DIR / f"program_loadings_{model_name}_{timestamp}.npy", W)
            np.save(RESULTS_DIR / f"program_components_{model_name}_{timestamp}.npy", H)

            # Save directions
            if disease_dir is not None:
                np.save(RESULTS_DIR / f"disease_direction_{model_name}_{timestamp}.npy", disease_dir)
            if recovery_dir is not None:
                np.save(RESULTS_DIR / f"recovery_direction_{model_name}_{timestamp}.npy", recovery_dir)

        except Exception as e:
            logger.error(f"Error processing {model_name}: {e}")
            continue

    # Save comparison results
    if all_results:
        results_df = pd.concat(all_results, ignore_index=True)
        results_df.to_csv(RESULTS_DIR / f"program_comparison_{timestamp}.csv", index=False)
        logger.info(f"\nResults saved to {RESULTS_DIR}")

    # Summary
    logger.info("\n" + "=" * 80)
    logger.info("SUMMARY")
    logger.info("=" * 80)

    if all_results:
        for model_name in MODELS:
            model_results = results_df[results_df['model'] == model_name]
            best_row = model_results.loc[model_results['auc'].idxmax()]
            logger.info(f"{model_name}: Best n_programs={best_row['n_programs']:.0f}, AUC={best_row['auc']:.4f}")

if __name__ == "__main__":
    main()
