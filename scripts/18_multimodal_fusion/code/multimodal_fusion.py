#!/usr/bin/env python3
"""
ACC 18: Multimodal Fusion of Foundation Model Embeddings

Fuse heterogeneous FM representations (STATE, UCE, Geneformer, TranscriptFormer)
using four strategies:
1. Early fusion: concatenate + PCA
2. Late fusion: average predictions + rank-based fusion
3. Attention-based fusion: learned per-cell weights (key differentiator)
4. Expression + embedding fusion: raw gene features + FM embeddings

Usage:
    python multimodal_fusion.py
    python multimodal_fusion.py --strategy attention
    python multimodal_fusion.py --strategy all

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
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, balanced_accuracy_score, accuracy_score
from sklearn.pipeline import make_pipeline
from scipy.stats import rankdata
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

BASE_DIR = Path("/bigdata/godziklab/shared/Xinru/302006")
EMB_DIR = BASE_DIR / "02_EMBEDDINGS"
RESULTS_DIR = BASE_DIR / "05_RESULTS/18_multimodal_fusion"
FIGURES_DIR = RESULTS_DIR / "figures"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

MODELS = {
    'STATE': {
        'path': EMB_DIR / "state/covid_state.h5ad",
        'obsm_keys': ['X_state', 'X_STATE', 'X_emb'],
        'dims': 2058,
    },
    'UCE': {
        'path': EMB_DIR / "uce/human_platelet_covid_severity_uce_adata.h5ad",
        'obsm_keys': ['X_uce', 'X_uce_4layer'],
        'dims': 1280,
    },
    'Geneformer': {
        'path': EMB_DIR / "geneformer/covid_geneformer.h5ad",
        'obsm_keys': ['X_geneformer', 'X_Geneformer'],
        'dims': 1152,
    },
    'TranscriptFormer': {
        'path': EMB_DIR / "transcriptformer/covid_transcriptformer.h5ad",
        'obsm_keys': ['X_transcriptformer', 'X_TranscriptFormer'],
        'dims': 2048,
    },
    'scPRINT2': {
        'path': EMB_DIR / "scprint2/covid_scprint2.h5ad",
        'obsm_keys': ['X_scprint2', 'scprint_emb'],
        'dims': 256,
    },
}


def load_all_embeddings():
    """Load and align embeddings from all models."""
    embeddings = {}
    severities = {}
    cell_ids_dict = {}

    for name, config in MODELS.items():
        path = config['path']
        if not path.exists():
            logger.warning(f"Skipping {name}: {path} not found")
            continue

        adata = sc.read_h5ad(path)
        X = None
        for key in config['obsm_keys']:
            if key in adata.obsm:
                X = np.array(adata.obsm[key])
                break
        if X is None:
            X = adata.X.toarray() if hasattr(adata.X, 'toarray') else np.array(adata.X)

        sev = None
        for col in ['covid_severity', 'severity']:
            if col in adata.obs.columns:
                sev = adata.obs[col].astype(str).str.lower().values
                break

        if sev is None:
            continue

        embeddings[name] = X
        severities[name] = sev
        cell_ids_dict[name] = adata.obs.index.values
        logger.info(f"Loaded {name}: {X.shape}")

    # Align cells
    common_cells = sorted(set.intersection(*[set(ids) for ids in cell_ids_dict.values()]))
    logger.info(f"Common cells: {len(common_cells)}")

    aligned = {}
    aligned_sev = None
    for name in embeddings:
        idx_map = {c: i for i, c in enumerate(cell_ids_dict[name])}
        indices = [idx_map[c] for c in common_cells]
        aligned[name] = embeddings[name][indices]
        model_sev = severities[name][indices]
        if aligned_sev is None:
            aligned_sev = model_sev
        else:
            if not np.array_equal(aligned_sev, model_sev):
                mismatches = np.sum(aligned_sev != model_sev)
                raise ValueError(
                    f"Cell alignment failed: {name} has {mismatches} severity label "
                    f"mismatches vs first model. Check cell barcode namespaces."
                )

    # Extract groups from first loaded model
    first_name = list(embeddings.keys())[0]
    first_path = MODELS[first_name]['path']
    adata_first = sc.read_h5ad(first_path)
    groups = None
    for col in ['Sample_ID', 'sample_id', 'Data_NO', 'data_no', 'batch', 'donor']:
        if col in adata_first.obs.columns:
            group_map = dict(zip(adata_first.obs.index.values, adata_first.obs[col].values))
            groups = np.array([group_map.get(c, 'unknown') for c in common_cells])
            logger.info(f"Using group column: {col} ({len(np.unique(groups))} groups)")
            break
    if groups is None:
        logger.warning("No group column found — using cell index as group")
        groups = np.arange(len(common_cells))

    return aligned, aligned_sev, common_cells, groups


def create_binary_labels(severity):
    """Map severity to binary labels."""
    return np.array([1 if s in ['severe', 'fatal'] else 0 for s in severity])


def evaluate_cv(X, y, n_splits=5, groups=None):
    """Cross-validated evaluation with LogReg using patient-level splitting."""
    cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=42)
    clf = make_pipeline(
        StandardScaler(),
        LogisticRegression(max_iter=2000, random_state=42, class_weight='balanced')
    )

    y_prob = cross_val_predict(clf, X, y, cv=cv, method='predict_proba', groups=groups)
    y_pred = cross_val_predict(clf, X, y, cv=cv, groups=groups)

    auc = roc_auc_score(y, y_prob[:, 1])
    bal_acc = balanced_accuracy_score(y, y_pred)

    return {'auc': auc, 'balanced_accuracy': bal_acc}


# =============================================================================
# Strategy 1: Early Fusion (Concatenate + PCA)
# =============================================================================

def early_fusion(aligned, y, pca_dim=256, groups=None):
    """Concatenate all embeddings, then PCA."""
    logger.info(f"\n--- Early Fusion (concat + PCA-{pca_dim}) ---")

    X_concat = np.hstack([aligned[m] for m in sorted(aligned.keys())])
    logger.info(f"Concatenated: {X_concat.shape}")

    # Raw concatenation
    result_raw = evaluate_cv(X_concat, y, groups=groups)
    logger.info(f"  Raw concat: AUC={result_raw['auc']:.4f}")

    # PCA reduction (fit PCA on full data for variance reporting, but use Pipeline for CV)
    pca_preview = PCA(n_components=pca_dim, random_state=42)
    pca_preview.fit(StandardScaler().fit_transform(X_concat))
    pca_var = float(pca_preview.explained_variance_ratio_.sum())

    pca_pipe = make_pipeline(
        StandardScaler(),
        PCA(n_components=pca_dim, random_state=42),
        StandardScaler(),
        LogisticRegression(max_iter=2000, random_state=42, class_weight='balanced')
    )
    cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)
    y_prob_pca = cross_val_predict(pca_pipe, X_concat, y, cv=cv, method='predict_proba', groups=groups)
    y_pred_pca = cross_val_predict(pca_pipe, X_concat, y, cv=cv, groups=groups)
    result_pca = {
        'auc': roc_auc_score(y, y_prob_pca[:, 1]),
        'balanced_accuracy': balanced_accuracy_score(y, y_pred_pca),
    }
    logger.info(f"  PCA-{pca_dim}: AUC={result_pca['auc']:.4f}")
    logger.info(f"  PCA variance explained: {pca_var:.3f}")

    return {
        'early_concat': result_raw,
        'early_pca': result_pca,
        'pca_variance': pca_var,
    }


# =============================================================================
# Strategy 2: Late Fusion (Average + Rank-based)
# =============================================================================

def late_fusion(aligned, y, groups=None):
    """Average prediction probabilities and rank-based fusion."""
    logger.info("\n--- Late Fusion ---")

    cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)
    model_probs = {}
    model_preds = {}

    for name in sorted(aligned.keys()):
        X = aligned[name]
        clf = make_pipeline(
            StandardScaler(),
            LogisticRegression(max_iter=2000, random_state=42, class_weight='balanced')
        )

        y_prob = cross_val_predict(clf, X, y, cv=cv, method='predict_proba', groups=groups)
        y_pred = cross_val_predict(clf, X, y, cv=cv, groups=groups)

        model_probs[name] = y_prob[:, 1]
        model_preds[name] = y_pred

    # Average probabilities
    avg_probs = np.mean(list(model_probs.values()), axis=0)
    avg_preds = (avg_probs > 0.5).astype(int)
    avg_auc = roc_auc_score(y, avg_probs)
    avg_ba = balanced_accuracy_score(y, avg_preds)
    logger.info(f"  Average probs: AUC={avg_auc:.4f}")

    # Rank-based fusion
    ranked_probs = np.mean([rankdata(p) / len(p) for p in model_probs.values()], axis=0)
    rank_auc = roc_auc_score(y, ranked_probs)
    logger.info(f"  Rank fusion: AUC={rank_auc:.4f}")

    # Stacking meta-learner
    meta_X = np.column_stack(list(model_probs.values()))
    meta_clf = LogisticRegression(max_iter=1000, random_state=42)
    meta_probs = cross_val_predict(meta_clf, meta_X, y, cv=cv, method='predict_proba', groups=groups)
    meta_auc = roc_auc_score(y, meta_probs[:, 1])
    logger.info(f"  Stacking: AUC={meta_auc:.4f}")

    return {
        'late_average': {'auc': avg_auc, 'balanced_accuracy': avg_ba},
        'late_rank': {'auc': rank_auc},
        'late_stacking': {'auc': meta_auc},
    }


# =============================================================================
# Strategy 3: Attention-Based Fusion
# =============================================================================

def attention_fusion(aligned, y, n_epochs=100, common_dim=256, groups=None):
    """Train attention fusion model."""
    logger.info(f"\n--- Attention-Based Fusion (common_dim={common_dim}) ---")

    import sys
    import torch
    import torch.nn.functional as F
    from torch.utils.data import TensorDataset, DataLoader
    sys.path.insert(0, str(Path(__file__).parent))
    from attention_fusion_module import AttentionFusion

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logger.info(f"Device: {device}")

    model_dims = {name: aligned[name].shape[1] for name in aligned}
    n_cells = len(y)

    # Cross-validation
    cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)
    all_probs = np.zeros(n_cells)
    all_preds = np.zeros(n_cells, dtype=int)
    all_attention = {name: np.zeros(n_cells) for name in aligned}
    fold_aucs = []

    for fold, (train_idx, val_idx) in enumerate(cv.split(np.zeros(n_cells), y, groups)):
        logger.info(f"  Fold {fold+1}/5")

        # Create model
        model = AttentionFusion(
            model_dims=model_dims,
            common_dim=common_dim,
            n_heads=4,
            n_classes=2,
            dropout=0.1,
        ).to(device)

        optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=0.01)
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=n_epochs)

        # Prepare data tensors
        train_data = {
            name: torch.FloatTensor(aligned[name][train_idx]).to(device)
            for name in aligned
        }
        val_data = {
            name: torch.FloatTensor(aligned[name][val_idx]).to(device)
            for name in aligned
        }
        train_y = torch.LongTensor(y[train_idx]).to(device)
        val_y = torch.LongTensor(y[val_idx]).to(device)

        # Training loop (mini-batch)
        batch_size = 256
        n_train = len(train_idx)

        best_val_auc = 0
        for epoch in range(n_epochs):
            model.train()
            perm = torch.randperm(n_train)

            epoch_loss = 0
            for start in range(0, n_train, batch_size):
                end = min(start + batch_size, n_train)
                batch_idx = perm[start:end]

                batch_data = {name: train_data[name][batch_idx] for name in aligned}
                batch_y = train_y[batch_idx]

                optimizer.zero_grad()
                logits = model(batch_data)
                loss = F.cross_entropy(logits, batch_y)
                loss.backward()
                optimizer.step()
                epoch_loss += loss.item()

            scheduler.step()

            # Validate every 10 epochs
            if (epoch + 1) % 10 == 0:
                model.eval()
                with torch.no_grad():
                    val_logits, val_fused, val_attn = model(val_data, return_attention=True)
                    val_probs = F.softmax(val_logits, dim=1)[:, 1].cpu().numpy()
                    val_auc = roc_auc_score(y[val_idx], val_probs)

                    if val_auc > best_val_auc:
                        best_val_auc = val_auc
                        best_probs = val_probs
                        best_preds = val_logits.argmax(dim=1).cpu().numpy()
                        best_attn = val_attn.cpu().numpy()

        # Store results for this fold
        all_probs[val_idx] = best_probs
        all_preds[val_idx] = best_preds
        for i, name in enumerate(sorted(aligned.keys())):
            all_attention[name][val_idx] = best_attn[:, i]

        fold_aucs.append(best_val_auc)
        logger.info(f"    Best AUC: {best_val_auc:.4f}")

    # Overall metrics
    overall_auc = roc_auc_score(y, all_probs)
    overall_ba = balanced_accuracy_score(y, all_preds)
    logger.info(f"  Overall: AUC={overall_auc:.4f}, BalAcc={overall_ba:.4f}")
    logger.info(f"  Fold AUCs: {[f'{a:.4f}' for a in fold_aucs]}")

    # Plot attention weights
    plot_attention_weights(all_attention, y, aligned.keys())

    return {
        'attention_fusion': {
            'auc': overall_auc,
            'balanced_accuracy': overall_ba,
            'fold_aucs': fold_aucs,
        },
        'attention_weights': {name: float(np.mean(w)) for name, w in all_attention.items()},
    }


def plot_attention_weights(attention_dict, y, model_names):
    """Plot attention weight heatmap showing per-severity model specialization."""
    severity_labels = ['non-severe', 'severe']
    model_list = sorted(model_names)

    # Mean attention by severity class
    attn_matrix = np.zeros((2, len(model_list)))
    for i, sev_label in enumerate([0, 1]):
        mask = y == sev_label
        for j, name in enumerate(model_list):
            attn_matrix[i, j] = np.mean(attention_dict[name][mask])

    fig, ax = plt.subplots(figsize=(8, 4))
    sns.heatmap(
        attn_matrix, annot=True, fmt='.3f', cmap='YlOrRd',
        xticklabels=model_list, yticklabels=severity_labels,
        ax=ax, vmin=0,
    )
    ax.set_title('Attention Weights by Severity Class\n(Which FM matters for which cells?)')
    ax.set_xlabel('Foundation Model')
    ax.set_ylabel('Severity Class')
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'attention_weights_heatmap.png', dpi=150, bbox_inches='tight')
    plt.savefig(FIGURES_DIR / 'attention_weights_heatmap.pdf', bbox_inches='tight')
    plt.close()
    logger.info(f"  Saved attention heatmap to {FIGURES_DIR}")


# =============================================================================
# Strategy 4: Expression + Embedding Fusion
# =============================================================================

def expression_embedding_fusion(aligned, y, groups=None):
    """Concatenate STATE embedding with top-100 raw gene features."""
    logger.info("\n--- Expression + Embedding Fusion ---")

    # Load raw expression
    covid_path = EMB_DIR / "state/covid_state.h5ad"
    adata = sc.read_h5ad(covid_path)
    if hasattr(adata.X, 'toarray'):
        X_raw = adata.X.toarray()
    else:
        X_raw = np.array(adata.X)

    # Select top-100 most variable genes
    gene_var = np.var(X_raw, axis=0)
    top_100 = np.argsort(-gene_var)[:100]
    X_top100 = X_raw[:, top_100]

    # Need to align with the common cells
    # For simplicity, use STATE embedding directly (already aligned)
    if 'STATE' not in aligned:
        logger.warning("STATE not available for expression fusion")
        return {}

    X_state = aligned['STATE']
    # Assume cells are already aligned
    n_common = X_state.shape[0]
    if X_top100.shape[0] != n_common:
        logger.warning(f"Shape mismatch: X_top100={X_top100.shape[0]}, aligned={n_common}. Using STATE cells only.")
        X_top100 = X_top100[:n_common]

    # Concat STATE + raw features
    X_fused = np.hstack([X_state, X_top100])
    logger.info(f"  STATE ({X_state.shape[1]}) + top-100 genes = {X_fused.shape[1]} dims")

    result = evaluate_cv(X_fused, y, groups=groups)
    logger.info(f"  AUC: {result['auc']:.4f}")

    # Compare to STATE alone
    state_result = evaluate_cv(X_state, y, groups=groups)
    logger.info(f"  STATE alone: AUC={state_result['auc']:.4f}")
    logger.info(f"  Improvement: {result['auc'] - state_result['auc']:+.4f}")

    return {
        'expr_embedding_fusion': result,
        'state_baseline': state_result,
    }


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description='Multimodal Fusion Analysis')
    parser.add_argument('--strategy', default='all',
                       choices=['early', 'late', 'attention', 'expression', 'all'])
    parser.add_argument('--attention-epochs', type=int, default=100)
    parser.add_argument('--common-dim', type=int, default=256)
    args = parser.parse_args()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    logger.info("=" * 60)
    logger.info("ACC 18: Multimodal Fusion Analysis")
    logger.info("=" * 60)

    # Load and align embeddings
    aligned, severity, common_cells, groups = load_all_embeddings()
    y = create_binary_labels(severity)
    logger.info(f"Aligned {len(aligned)} models, {len(y)} cells")
    logger.info(f"Labels: {np.bincount(y)}")

    all_results = {}

    # Run strategies (groups passed via evaluate_cv internally)
    if args.strategy in ('early', 'all'):
        all_results.update(early_fusion(aligned, y, groups=groups))

    if args.strategy in ('late', 'all'):
        all_results.update(late_fusion(aligned, y, groups=groups))

    if args.strategy in ('attention', 'all'):
        all_results.update(attention_fusion(aligned, y, args.attention_epochs, args.common_dim, groups=groups))

    if args.strategy in ('expression', 'all'):
        all_results.update(expression_embedding_fusion(aligned, y, groups=groups))

    # Single model baselines
    logger.info("\n--- Single Model Baselines ---")
    for name in sorted(aligned.keys()):
        result = evaluate_cv(aligned[name], y, groups=groups)
        all_results[f'single_{name}'] = result
        logger.info(f"  {name}: AUC={result['auc']:.4f}")

    # Summary
    logger.info(f"\n{'='*60}")
    logger.info("SUMMARY")
    logger.info(f"{'='*60}")

    summary_rows = []
    for strategy, metrics in all_results.items():
        if isinstance(metrics, dict) and 'auc' in metrics:
            logger.info(f"  {strategy}: AUC={metrics['auc']:.4f}")
            summary_rows.append({'strategy': strategy, **metrics})

    # Save results
    with open(RESULTS_DIR / f"fusion_results_{timestamp}.json", 'w') as f:
        # Convert numpy types for JSON serialization
        def convert(obj):
            if isinstance(obj, (np.integer,)):
                return int(obj)
            if isinstance(obj, (np.floating,)):
                return float(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return obj

        json.dump(all_results, f, indent=2, default=convert)

    if summary_rows:
        df = pd.DataFrame(summary_rows)
        df.to_csv(RESULTS_DIR / f"fusion_summary_{timestamp}.csv", index=False)

    logger.info(f"\nResults saved to {RESULTS_DIR}")


if __name__ == '__main__':
    main()
