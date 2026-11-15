#!/usr/bin/env python3
"""
Step 5: 5-Modality Attention Fusion (STATE + UCE + Geneformer + TranscriptFormer + ESM_protein).

Extends ACC 18's 4-model attention fusion by adding per-cell ESM protein
embeddings as a 5th input modality. Runs all fusion strategies and compares
5-model vs 4-model baselines.

Usage:
    python step5_multimodal_fusion_5model.py
    python step5_multimodal_fusion_5model.py --strategy attention
    python step5_multimodal_fusion_5model.py --strategy all

Author: Xinru Qiu
Date: April 2026
"""

import argparse
import logging
from pathlib import Path
import json
import numpy as np
import pandas as pd
import scanpy as sc
from datetime import datetime
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, balanced_accuracy_score
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


# ============================================================================
# Configuration
# ============================================================================

FM_MODELS = {
    'STATE': {
        'obsm_keys': ['X_state', 'X_STATE', 'X_emb'],
        'dims': 2058,
    },
    'UCE': {
        'obsm_keys': ['X_uce', 'X_uce_4layer'],
        'dims': 1280,
    },
    'Geneformer': {
        'obsm_keys': ['X_geneformer', 'X_Geneformer'],
        'dims': 1152,
    },
    'TranscriptFormer': {
        'obsm_keys': ['X_transcriptformer', 'X_TranscriptFormer'],
        'dims': 2048,
    },
}


def load_all_embeddings(base_dir, esm_results_dir):
    """Load and align embeddings from all 5 modalities."""
    emb_dir = Path(base_dir) / "02_EMBEDDINGS"

    # Paths for FM embeddings
    fm_paths = {
        'STATE': emb_dir / "state/covid_state.h5ad",
        'UCE': emb_dir / "uce/human_platelet_covid_severity_uce_adata.h5ad",
        'Geneformer': emb_dir / "geneformer/covid_geneformer.h5ad",
        'TranscriptFormer': emb_dir / "transcriptformer/covid_transcriptformer.h5ad",
    }

    embeddings = {}
    severities = {}
    cell_ids_dict = {}

    # Load 4 FM embeddings
    for name, path in fm_paths.items():
        if not path.exists():
            logger.warning(f"Skipping {name}: {path} not found")
            continue

        adata = sc.read_h5ad(path)
        X = None
        for key in FM_MODELS[name]['obsm_keys']:
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
            logger.warning(f"Skipping {name}: no severity column")
            continue

        embeddings[name] = X
        severities[name] = sev
        cell_ids_dict[name] = adata.obs.index.values
        logger.info(f"Loaded {name}: {X.shape}")

    # Load ESM protein cell embeddings
    esm_h5ad = Path(esm_results_dir) / "cell_protein_embeddings.h5ad"
    if esm_h5ad.exists():
        adata_esm = sc.read_h5ad(esm_h5ad)
        if 'X_esm_protein' in adata_esm.obsm:
            X_esm = np.array(adata_esm.obsm['X_esm_protein'])
        else:
            X_esm = adata_esm.X.toarray() if hasattr(adata_esm.X, 'toarray') else np.array(adata_esm.X)

        embeddings['ESM_protein'] = X_esm
        cell_ids_dict['ESM_protein'] = adata_esm.obs.index.values
        logger.info(f"Loaded ESM_protein: {X_esm.shape}")

        # Severity from ESM adata obs if present, else will use from FM
        for col in ['covid_severity', 'severity']:
            if col in adata_esm.obs.columns:
                severities['ESM_protein'] = adata_esm.obs[col].astype(str).str.lower().values
                break
    else:
        logger.warning(f"ESM protein embeddings not found at {esm_h5ad}")
        logger.warning("Running with 4 models only (no ESM_protein)")

    # Align cells across all loaded modalities
    common_cells = sorted(set.intersection(*[set(ids) for ids in cell_ids_dict.values()]))
    logger.info(f"Common cells across {len(embeddings)} modalities: {len(common_cells)}")

    aligned = {}
    aligned_sev = None
    for name in embeddings:
        idx_map = {c: i for i, c in enumerate(cell_ids_dict[name])}
        indices = [idx_map[c] for c in common_cells]
        aligned[name] = embeddings[name][indices]
        if aligned_sev is None and name in severities:
            aligned_sev = severities[name][indices]

    # Extract patient-level groups for GroupKFold (from first loaded adata)
    groups = None
    first_path = next(iter(fm_paths.values()))
    if first_path.exists():
        adata_first = sc.read_h5ad(first_path)
        for col in ['Sample_ID', 'sample_id', 'Data_NO', 'data_no', 'batch', 'donor']:
            if col in adata_first.obs.columns:
                idx_map = {c: i for i, c in enumerate(adata_first.obs.index.values)}
                indices = [idx_map[c] for c in common_cells]
                groups = adata_first.obs[col].values[indices]
                logger.info(f"Using group column: {col} ({len(np.unique(groups))} groups)")
                break
        del adata_first
    if groups is None:
        logger.warning("No group column found; StratifiedGroupKFold will fall back to cell-level splitting")

    return aligned, aligned_sev, common_cells, groups


def create_binary_labels(severity):
    """Map severity to binary labels."""
    return np.array([1 if s in ['severe', 'fatal'] else 0 for s in severity])


def evaluate_cv(X, y, n_splits=5, groups=None):
    """Cross-validated evaluation with LogReg."""
    cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=42)
    clf = make_pipeline(
        StandardScaler(),
        PCA(n_components=min(256, X.shape[1]), random_state=42),
        LogisticRegression(max_iter=2000, random_state=42, class_weight='balanced', solver='lbfgs')
    )
    y_prob = cross_val_predict(clf, X, y, cv=cv, groups=groups, method='predict_proba')
    y_pred = cross_val_predict(clf, X, y, cv=cv, groups=groups)
    return {
        'auc': float(roc_auc_score(y, y_prob[:, 1])),
        'balanced_accuracy': float(balanced_accuracy_score(y, y_pred)),
    }


# ============================================================================
# Strategy 1: Early Fusion (Concatenate + PCA)
# ============================================================================

def early_fusion(aligned, y, pca_dim=256, groups=None):
    """Concatenate all embeddings, PCA, then classify."""
    logger.info(f"\n--- Early Fusion (concat + PCA-{pca_dim}) ---")

    X_concat = np.hstack([aligned[m] for m in sorted(aligned.keys())])
    logger.info(f"Concatenated: {X_concat.shape}")

    pca_pipe = make_pipeline(
        StandardScaler(),
        PCA(n_components=pca_dim, random_state=42),
        StandardScaler(),
        LogisticRegression(max_iter=2000, random_state=42, class_weight='balanced', solver='lbfgs')
    )
    cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)
    y_prob = cross_val_predict(pca_pipe, X_concat, y, cv=cv, groups=groups, method='predict_proba')
    y_pred = cross_val_predict(pca_pipe, X_concat, y, cv=cv, groups=groups)
    result = {
        'auc': float(roc_auc_score(y, y_prob[:, 1])),
        'balanced_accuracy': float(balanced_accuracy_score(y, y_pred)),
    }
    logger.info(f"  PCA-{pca_dim}: AUC={result['auc']:.4f}")

    return {'early_fusion_5model': result}


# ============================================================================
# Strategy 2: Late Fusion (Average + Rank + Stacking)
# ============================================================================

def late_fusion(aligned, y, groups=None):
    """Average prediction probabilities and rank-based fusion."""
    logger.info("\n--- Late Fusion ---")

    cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)
    model_probs = {}

    for name in sorted(aligned.keys()):
        X = aligned[name]
        clf = make_pipeline(
            StandardScaler(),
            PCA(n_components=min(256, X.shape[1]), random_state=42),
            LogisticRegression(max_iter=2000, random_state=42, class_weight='balanced', solver='lbfgs')
        )
        y_prob = cross_val_predict(clf, X, y, cv=cv, groups=groups, method='predict_proba')
        model_probs[name] = y_prob[:, 1]
        individual_auc = roc_auc_score(y, y_prob[:, 1])
        logger.info(f"  {name}: AUC={individual_auc:.4f}")

    # Average probabilities
    avg_probs = np.mean(list(model_probs.values()), axis=0)
    avg_preds = (avg_probs > 0.5).astype(int)
    avg_auc = float(roc_auc_score(y, avg_probs))
    avg_ba = float(balanced_accuracy_score(y, avg_preds))
    logger.info(f"  Average probs: AUC={avg_auc:.4f}")

    # Rank-based fusion
    ranked_probs = np.mean([rankdata(p) / len(p) for p in model_probs.values()], axis=0)
    rank_auc = float(roc_auc_score(y, ranked_probs))
    logger.info(f"  Rank fusion: AUC={rank_auc:.4f}")

    # Stacking meta-learner
    meta_X = np.column_stack(list(model_probs.values()))
    meta_clf = LogisticRegression(max_iter=1000, random_state=42)
    meta_probs = cross_val_predict(meta_clf, meta_X, y, cv=cv, groups=groups, method='predict_proba')
    meta_auc = float(roc_auc_score(y, meta_probs[:, 1]))
    logger.info(f"  Stacking: AUC={meta_auc:.4f}")

    return {
        'late_average_5model': {'auc': avg_auc, 'balanced_accuracy': avg_ba},
        'late_rank_5model': {'auc': rank_auc},
        'late_stacking_5model': {'auc': meta_auc},
    }


# ============================================================================
# Strategy 3: Attention-Based Fusion (main contribution)
# ============================================================================

def attention_fusion(aligned, y, n_epochs=100, common_dim=256, figures_dir=None, groups=None):
    """Train 5-model attention fusion and compare to 4-model baseline."""
    logger.info(f"\n--- Attention-Based Fusion (5-model, common_dim={common_dim}) ---")

    import sys
    import torch
    import torch.nn.functional as F

    # Import AttentionFusion from ACC 18
    acc18_code = Path(__file__).resolve().parent.parent.parent / "18_multimodal_fusion" / "code"
    sys.path.insert(0, str(acc18_code))
    from attention_fusion_module import AttentionFusion

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logger.info(f"Device: {device}")

    n_cells = len(y)

    def _train_attention(data_dict, label, tag=""):
        """Train attention fusion on given modalities and return results."""
        model_dims = {name: data_dict[name].shape[1] for name in data_dict}

        cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)
        all_probs = np.zeros(n_cells)
        all_preds = np.zeros(n_cells, dtype=int)
        all_attention = {name: np.zeros(n_cells) for name in data_dict}
        fold_aucs = []

        for fold, (train_idx, val_idx) in enumerate(cv.split(np.zeros(n_cells), label, groups=groups)):
            model = AttentionFusion(
                model_dims=model_dims,
                common_dim=common_dim,
                n_heads=4,
                n_classes=2,
                dropout=0.1,
            ).to(device)

            optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=0.01)
            scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=n_epochs)

            train_data = {
                name: torch.FloatTensor(data_dict[name][train_idx]).to(device)
                for name in data_dict
            }
            val_data = {
                name: torch.FloatTensor(data_dict[name][val_idx]).to(device)
                for name in data_dict
            }
            train_y = torch.LongTensor(label[train_idx]).to(device)

            batch_size = 256
            n_train = len(train_idx)
            best_val_auc = 0
            best_probs = None
            best_preds = None
            best_attn = None

            for epoch in range(n_epochs):
                model.train()
                perm = torch.randperm(n_train)
                for start in range(0, n_train, batch_size):
                    end = min(start + batch_size, n_train)
                    batch_idx = perm[start:end]
                    batch_data = {name: train_data[name][batch_idx] for name in data_dict}
                    batch_y = train_y[batch_idx]

                    optimizer.zero_grad()
                    logits = model(batch_data)
                    loss = F.cross_entropy(logits, batch_y)
                    loss.backward()
                    optimizer.step()

                scheduler.step()

                if (epoch + 1) % 10 == 0:
                    model.eval()
                    with torch.no_grad():
                        val_logits, val_fused, val_attn = model(val_data, return_attention=True)
                        val_probs = F.softmax(val_logits, dim=1)[:, 1].cpu().numpy()
                        val_auc = roc_auc_score(label[val_idx], val_probs)
                        if val_auc > best_val_auc:
                            best_val_auc = val_auc
                            best_probs = val_probs
                            best_preds = val_logits.argmax(dim=1).cpu().numpy()
                            best_attn = val_attn.cpu().numpy()

            all_probs[val_idx] = best_probs
            all_preds[val_idx] = best_preds
            for i, name in enumerate(sorted(data_dict.keys())):
                all_attention[name][val_idx] = best_attn[:, i]
            fold_aucs.append(best_val_auc)
            logger.info(f"  {tag} Fold {fold+1}/5: AUC={best_val_auc:.4f}")

        overall_auc = float(roc_auc_score(label, all_probs))
        overall_ba = float(balanced_accuracy_score(label, all_preds))
        logger.info(f"  {tag} Overall: AUC={overall_auc:.4f}, BalAcc={overall_ba:.4f}")

        return {
            'auc': overall_auc,
            'balanced_accuracy': overall_ba,
            'fold_aucs': [float(a) for a in fold_aucs],
        }, all_attention

    # --- 5-model fusion ---
    result_5, attn_5 = _train_attention(aligned, y, tag="5-model")

    # --- 4-model baseline (exclude ESM_protein) ---
    result_4 = None
    attn_4 = None
    if 'ESM_protein' in aligned:
        aligned_4 = {k: v for k, v in aligned.items() if k != 'ESM_protein'}
        result_4, attn_4 = _train_attention(aligned_4, y, tag="4-model")
        delta = result_5['auc'] - result_4['auc']
        logger.info(f"\n  Delta (5-model - 4-model): AUC {delta:+.4f}")

    # Plot attention weights
    if figures_dir is not None:
        _plot_attention_comparison(attn_5, attn_4, y, figures_dir)

    results = {'attention_fusion_5model': result_5}
    if result_4 is not None:
        results['attention_fusion_4model_baseline'] = result_4
    results['attention_weights_5model'] = {
        name: float(np.mean(w)) for name, w in attn_5.items()
    }

    return results


def _plot_attention_comparison(attn_5, attn_4, y, figures_dir):
    """Plot 5-model attention heatmap and comparison to 4-model."""
    figures_dir = Path(figures_dir)
    figures_dir.mkdir(parents=True, exist_ok=True)

    severity_labels = ['non-severe', 'severe']
    model_list = sorted(attn_5.keys())

    # 5-model heatmap
    attn_matrix = np.zeros((2, len(model_list)))
    for i, sev_label in enumerate([0, 1]):
        mask = y == sev_label
        for j, name in enumerate(model_list):
            attn_matrix[i, j] = np.mean(attn_5[name][mask])

    fig, ax = plt.subplots(figsize=(10, 4))
    sns.heatmap(
        attn_matrix, annot=True, fmt='.3f', cmap='YlOrRd',
        xticklabels=model_list, yticklabels=severity_labels,
        ax=ax, vmin=0,
    )
    ax.set_title('5-Modality Attention Weights by Severity Class\n(Including ESM Protein Embeddings)')
    ax.set_xlabel('Modality')
    ax.set_ylabel('Severity Class')
    plt.tight_layout()
    plt.savefig(figures_dir / 'attention_weights_5model_heatmap.png', dpi=150, bbox_inches='tight')
    plt.savefig(figures_dir / 'attention_weights_5model_heatmap.pdf', bbox_inches='tight')
    plt.close()

    # Bar chart comparing mean attention per modality
    fig, ax = plt.subplots(figsize=(8, 5))
    means = [float(np.mean(attn_5[name])) for name in model_list]
    colors = ['#4C72B0', '#DD8452', '#55A868', '#C44E52', '#8172B2']
    bars = ax.bar(model_list, means, color=colors[:len(model_list)])
    ax.set_ylabel('Mean Attention Weight')
    ax.set_title('Per-Modality Attention Weights (5-Model Fusion)')
    ax.set_ylim(0, max(means) * 1.3)
    for bar, val in zip(bars, means):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                f'{val:.3f}', ha='center', va='bottom', fontsize=10)
    plt.tight_layout()
    plt.savefig(figures_dir / 'attention_weights_5model_bar.png', dpi=150, bbox_inches='tight')
    plt.savefig(figures_dir / 'attention_weights_5model_bar.pdf', bbox_inches='tight')
    plt.close()

    logger.info(f"  Saved attention plots to {figures_dir}")


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description='5-Modality Multimodal Fusion')
    parser.add_argument('--base-dir', type=str,
                        default='/bigdata/godziklab/shared/Xinru/302006')
    parser.add_argument('--esm-results-dir', type=str, default=None,
                        help='Directory containing cell_protein_embeddings.h5ad')
    parser.add_argument('--output-dir', type=str, default=None)
    parser.add_argument('--strategy', default='all',
                        choices=['early', 'late', 'attention', 'all'])
    parser.add_argument('--attention-epochs', type=int, default=100)
    parser.add_argument('--common-dim', type=int, default=256)
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    results_dir = script_dir.parent / "results" if args.output_dir is None else Path(args.output_dir)
    esm_results_dir = results_dir if args.esm_results_dir is None else Path(args.esm_results_dir)
    figures_dir = results_dir / "figures"
    results_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    logger.info("=" * 60)
    logger.info("ACC 22: 5-Modality Multimodal Fusion (with ESM Protein)")
    logger.info("=" * 60)

    # Load and align all 5 modalities
    aligned, severity, common_cells, groups = load_all_embeddings(args.base_dir, esm_results_dir)
    y = create_binary_labels(severity)

    n_models = len(aligned)
    has_esm = 'ESM_protein' in aligned
    logger.info(f"Loaded {n_models} modalities, {len(y)} cells")
    logger.info(f"  ESM_protein included: {has_esm}")
    logger.info(f"  Labels: severe={y.sum()}, non-severe={(1-y).sum()}")
    for name in sorted(aligned.keys()):
        logger.info(f"  {name}: {aligned[name].shape}")

    all_results = {
        'n_models': n_models,
        'n_cells': len(y),
        'has_esm_protein': has_esm,
        'models': sorted(aligned.keys()),
    }

    # Run strategies
    if args.strategy in ('early', 'all'):
        all_results.update(early_fusion(aligned, y, groups=groups))

    if args.strategy in ('late', 'all'):
        all_results.update(late_fusion(aligned, y, groups=groups))

    if args.strategy in ('attention', 'all'):
        all_results.update(attention_fusion(
            aligned, y, n_epochs=args.attention_epochs,
            common_dim=args.common_dim, figures_dir=figures_dir, groups=groups,
        ))

    # Single-model baselines
    logger.info("\n--- Single Model Baselines ---")
    for name in sorted(aligned.keys()):
        result = evaluate_cv(aligned[name], y, groups=groups)
        all_results[f'single_{name}'] = result
        logger.info(f"  {name}: AUC={result['auc']:.4f}")

    # Summary table
    logger.info(f"\n{'='*60}")
    logger.info("SUMMARY")
    logger.info(f"{'='*60}")

    summary_rows = []
    for key, val in all_results.items():
        if isinstance(val, dict) and 'auc' in val:
            logger.info(f"  {key}: AUC={val['auc']:.4f}")
            summary_rows.append({'strategy': key, **val})

    # Save results
    output_path = results_dir / f"fusion_5model_results_{timestamp}.json"
    with open(output_path, 'w') as f:
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
        csv_path = results_dir / f"fusion_5model_summary_{timestamp}.csv"
        df.to_csv(csv_path, index=False)
        logger.info(f"Saved summary CSV: {csv_path}")

    logger.info(f"Saved results JSON: {output_path}")
    logger.info("\nDone. Compare 5-model vs 4-model attention fusion AUC.")


if __name__ == '__main__':
    main()
