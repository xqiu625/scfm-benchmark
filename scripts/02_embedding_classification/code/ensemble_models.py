#!/usr/bin/env python3
"""
Ensemble Foundation Models

Combine embeddings from multiple foundation models for improved classification.
Tests whether models capture complementary biological signals.
"""

import re
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, VotingClassifier, StackingClassifier
from sklearn.model_selection import cross_val_predict, StratifiedGroupKFold
from sklearn.metrics import roc_auc_score, balanced_accuracy_score
from sklearn.neural_network import MLPClassifier
from pathlib import Path
from datetime import datetime
from itertools import combinations
import logging
import warnings
warnings.filterwarnings('ignore')

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Paths
BASE_DIR = Path("/bigdata/godziklab/shared/Xinru/302006")
EMB_DIR = BASE_DIR / "02_EMBEDDINGS"
RESULTS_DIR = BASE_DIR / "05_RESULTS" / "ensemble_models"
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

    # Get embeddings — check all known FM keys before falling back
    X = None
    for key in ['X_emb', 'X_state', 'X_STATE', 'X_uce', 'X_UCE',
                'X_scGPT', 'X_scgpt', 'X_transcriptformer', 'X_TranscriptFormer',
                'X_geneformer', 'X_Geneformer', 'X_scprint2', 'scprint_emb']:
        if key in adata.obsm:
            X = adata.obsm[key]
            logger.info(f"  Using embedding key: {key} ({X.shape})")
            break
    if X is None:
        raise ValueError(f"No embedding key found for {model_name}. "
                         f"Available obsm keys: {list(adata.obsm.keys())}")

    # Get severity
    if 'covid_severity' in adata.obs.columns:
        severity = adata.obs['covid_severity'].astype(str).str.lower().values
    elif 'severity' in adata.obs.columns:
        severity = adata.obs['severity'].values
    elif 'Severity' in adata.obs.columns:
        severity = adata.obs['Severity'].values
    else:
        raise ValueError("No severity column found")

    # Get cell IDs for alignment — normalize scGPT's "-batch_N" suffix
    cell_ids = np.array([re.sub(r'-batch_\d+$', '', str(c)) for c in adata.obs.index])

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

    return X, severity, cell_ids, groups

def align_embeddings(embeddings_dict, labels_dict, cell_ids_dict, groups_dict):
    """Align embeddings across models using common cells.

    Cell IDs are normalized in load_embeddings() (e.g., scGPT's '-batch_N' suffix
    is stripped) so barcode intersection should work across all models.
    """
    # Find common cells
    all_cell_sets = [set(ids) for ids in cell_ids_dict.values()]
    common_cells = set.intersection(*all_cell_sets)
    common_cells = sorted(list(common_cells))

    logger.info(f"Common cells across all models: {len(common_cells)}")

    if len(common_cells) == 0:
        raise ValueError(
            "No common cells found across models after barcode normalization. "
            "Check cell index formats in each h5ad file."
        )

    aligned_X = {}
    aligned_y = None
    aligned_groups = None

    for model_name, X in embeddings_dict.items():
        cell_ids = cell_ids_dict[model_name]
        cell_to_idx = {c: i for i, c in enumerate(cell_ids)}

        indices = [cell_to_idx[c] for c in common_cells]
        aligned_X[model_name] = X[indices]

        if aligned_y is None:
            aligned_y = labels_dict[model_name][indices]
            aligned_groups = groups_dict[model_name][indices]

    return aligned_X, aligned_y, common_cells, aligned_groups

def create_binary_labels(severity):
    """Create binary severe/non-severe labels."""
    severe_classes = ['severe', 'fatal']
    binary = np.array([1 if s in severe_classes else 0 for s in severity])
    return binary

def evaluate_single_model(X, y, cv_folds=5, groups=None):
    """Evaluate a single model."""
    cv = StratifiedGroupKFold(n_splits=cv_folds, shuffle=True, random_state=42)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    clf = LogisticRegression(max_iter=1000, random_state=42)
    y_prob = cross_val_predict(clf, X_scaled, y, cv=cv, groups=groups, method='predict_proba')
    y_pred = cross_val_predict(clf, X_scaled, y, cv=cv, groups=groups)

    auc = roc_auc_score(y, y_prob[:, 1])
    bal_acc = balanced_accuracy_score(y, y_pred)

    return auc, bal_acc

def evaluate_concatenated(X_dict, y, model_names, cv_folds=5, groups=None):
    """Evaluate concatenated embeddings."""
    # Concatenate embeddings
    X_concat = np.hstack([X_dict[m] for m in model_names])
    logger.info(f"  Concatenated dims: {X_concat.shape[1]}")

    return evaluate_single_model(X_concat, y, cv_folds, groups=groups)

def evaluate_averaged(X_dict, y, model_names, cv_folds=5, groups=None):
    """Evaluate averaged embeddings (after PCA to common dimension)."""
    from sklearn.decomposition import PCA

    # Project each to common dimension (min dimension)
    min_dim = min(X_dict[m].shape[1] for m in model_names)
    common_dim = min(min_dim, 256)

    projected = []
    for model_name in model_names:
        X = X_dict[model_name]
        if X.shape[1] > common_dim:
            pca = PCA(n_components=common_dim, random_state=42)
            X_proj = pca.fit_transform(X)
        else:
            X_proj = X
        projected.append(X_proj)

    # Average
    X_avg = np.mean(projected, axis=0)

    return evaluate_single_model(X_avg, y, cv_folds, groups=groups)

def evaluate_voting_ensemble(X_dict, y, model_names, cv_folds=5, groups=None):
    """Evaluate voting ensemble (train separate classifiers, average predictions)."""
    cv = StratifiedGroupKFold(n_splits=cv_folds, shuffle=True, random_state=42)

    # Collect predictions from each model
    all_probs = []
    for model_name in model_names:
        X = X_dict[model_name]
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        clf = LogisticRegression(max_iter=1000, random_state=42)
        y_prob = cross_val_predict(clf, X_scaled, y, cv=cv, groups=groups, method='predict_proba')
        all_probs.append(y_prob[:, 1])

    # Average predictions
    avg_probs = np.mean(all_probs, axis=0)
    avg_preds = (avg_probs > 0.5).astype(int)

    auc = roc_auc_score(y, avg_probs)
    bal_acc = balanced_accuracy_score(y, avg_preds)

    return auc, bal_acc

def evaluate_stacking_ensemble(X_dict, y, model_names, cv_folds=5, groups=None):
    """Evaluate stacking ensemble (meta-learner on predictions)."""
    cv = StratifiedGroupKFold(n_splits=cv_folds, shuffle=True, random_state=42)

    # Get predictions from each model
    meta_features = []
    for model_name in model_names:
        X = X_dict[model_name]
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        clf = LogisticRegression(max_iter=1000, random_state=42)
        y_prob = cross_val_predict(clf, X_scaled, y, cv=cv, groups=groups, method='predict_proba')
        meta_features.append(y_prob[:, 1].reshape(-1, 1))

    # Combine predictions as meta-features
    X_meta = np.hstack(meta_features)

    # Train meta-learner
    meta_clf = LogisticRegression(max_iter=1000, random_state=42)
    y_prob_meta = cross_val_predict(meta_clf, X_meta, y, cv=cv, groups=groups, method='predict_proba')
    y_pred_meta = cross_val_predict(meta_clf, X_meta, y, cv=cv, groups=groups)

    auc = roc_auc_score(y, y_prob_meta[:, 1])
    bal_acc = balanced_accuracy_score(y, y_pred_meta)

    return auc, bal_acc

def main():
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    logger.info("=" * 80)
    logger.info("ENSEMBLE FOUNDATION MODELS ANALYSIS")
    logger.info("=" * 80)

    # Load all embeddings
    embeddings_dict = {}
    labels_dict = {}
    cell_ids_dict = {}
    groups_dict = {}
    available_models = []

    for model_name, path in MODELS.items():
        if not path.exists():
            logger.warning(f"Skipping {model_name}: file not found")
            continue
        try:
            X, severity, cell_ids, groups = load_embeddings(model_name)
            embeddings_dict[model_name] = X
            labels_dict[model_name] = severity
            cell_ids_dict[model_name] = cell_ids
            groups_dict[model_name] = groups
            available_models.append(model_name)
            logger.info(f"Loaded {model_name}: {X.shape}")
        except Exception as e:
            logger.warning(f"Error loading {model_name}: {e}")

    if len(available_models) < 2:
        logger.error("Need at least 2 models for ensemble")
        return

    # Align embeddings
    aligned_X, severity, common_cells, aligned_groups = align_embeddings(embeddings_dict, labels_dict, cell_ids_dict, groups_dict)
    y = create_binary_labels(severity)

    logger.info(f"\nUsing {len(common_cells)} common cells")
    logger.info(f"Labels: {np.bincount(y)}")

    all_results = []

    # Single model baselines
    logger.info("\n" + "=" * 60)
    logger.info("SINGLE MODEL BASELINES")
    logger.info("=" * 60)

    for model_name in available_models:
        auc, bal_acc = evaluate_single_model(aligned_X[model_name], y, groups=aligned_groups)
        all_results.append({
            'ensemble_type': 'single',
            'models': model_name,
            'n_models': 1,
            'auc': auc,
            'balanced_accuracy': bal_acc
        })
        logger.info(f"{model_name}: AUC={auc:.4f}, BalAcc={bal_acc:.4f}")

    # Pairwise ensembles
    logger.info("\n" + "=" * 60)
    logger.info("PAIRWISE ENSEMBLES")
    logger.info("=" * 60)

    for model_pair in combinations(available_models, 2):
        pair_name = "+".join(model_pair)
        logger.info(f"\n{pair_name}:")

        # Concatenated
        auc, bal_acc = evaluate_concatenated(aligned_X, y, model_pair, groups=aligned_groups)
        all_results.append({
            'ensemble_type': 'concatenate',
            'models': pair_name,
            'n_models': 2,
            'auc': auc,
            'balanced_accuracy': bal_acc
        })
        logger.info(f"  Concatenate: AUC={auc:.4f}")

        # Voting
        auc, bal_acc = evaluate_voting_ensemble(aligned_X, y, model_pair, groups=aligned_groups)
        all_results.append({
            'ensemble_type': 'voting',
            'models': pair_name,
            'n_models': 2,
            'auc': auc,
            'balanced_accuracy': bal_acc
        })
        logger.info(f"  Voting: AUC={auc:.4f}")

        # Stacking
        auc, bal_acc = evaluate_stacking_ensemble(aligned_X, y, model_pair, groups=aligned_groups)
        all_results.append({
            'ensemble_type': 'stacking',
            'models': pair_name,
            'n_models': 2,
            'auc': auc,
            'balanced_accuracy': bal_acc
        })
        logger.info(f"  Stacking: AUC={auc:.4f}")

    # Top 3 ensemble (STATE, UCE, TranscriptFormer if available)
    top_models = [m for m in ['STATE', 'UCE', 'TranscriptFormer'] if m in available_models]
    if len(top_models) >= 3:
        logger.info("\n" + "=" * 60)
        logger.info("TOP 3 ENSEMBLE (STATE + UCE + TranscriptFormer)")
        logger.info("=" * 60)

        trio_name = "+".join(top_models[:3])

        for ensemble_type, eval_func in [
            ('concatenate', evaluate_concatenated),
            ('voting', evaluate_voting_ensemble),
            ('stacking', evaluate_stacking_ensemble)
        ]:
            auc, bal_acc = eval_func(aligned_X, y, top_models[:3], groups=aligned_groups)
            all_results.append({
                'ensemble_type': ensemble_type,
                'models': trio_name,
                'n_models': 3,
                'auc': auc,
                'balanced_accuracy': bal_acc
            })
            logger.info(f"{ensemble_type}: AUC={auc:.4f}, BalAcc={bal_acc:.4f}")

    # All models ensemble
    if len(available_models) >= 4:
        logger.info("\n" + "=" * 60)
        logger.info(f"ALL MODELS ENSEMBLE ({len(available_models)} models)")
        logger.info("=" * 60)

        all_name = "+".join(available_models)

        for ensemble_type, eval_func in [
            ('voting', evaluate_voting_ensemble),
            ('stacking', evaluate_stacking_ensemble)
        ]:
            auc, bal_acc = eval_func(aligned_X, y, available_models, groups=aligned_groups)
            all_results.append({
                'ensemble_type': ensemble_type,
                'models': all_name,
                'n_models': len(available_models),
                'auc': auc,
                'balanced_accuracy': bal_acc
            })
            logger.info(f"{ensemble_type}: AUC={auc:.4f}, BalAcc={bal_acc:.4f}")

    # Save results
    results_df = pd.DataFrame(all_results)
    output_file = RESULTS_DIR / f"ensemble_results_{timestamp}.csv"
    results_df.to_csv(output_file, index=False)
    logger.info(f"\nResults saved to {output_file}")

    # Summary
    logger.info("\n" + "=" * 80)
    logger.info("SUMMARY - BEST ENSEMBLES")
    logger.info("=" * 80)

    # Best single model
    single = results_df[results_df['ensemble_type'] == 'single'].sort_values('auc', ascending=False).iloc[0]
    logger.info(f"Best single: {single['models']} (AUC={single['auc']:.4f})")

    # Best ensemble
    ensemble = results_df[results_df['n_models'] > 1].sort_values('auc', ascending=False).iloc[0]
    logger.info(f"Best ensemble: {ensemble['models']} [{ensemble['ensemble_type']}] (AUC={ensemble['auc']:.4f})")

    improvement = ensemble['auc'] - single['auc']
    logger.info(f"Improvement: +{improvement:.4f} ({improvement/single['auc']*100:.1f}%)")

if __name__ == "__main__":
    main()
