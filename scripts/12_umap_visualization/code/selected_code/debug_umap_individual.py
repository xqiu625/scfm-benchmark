#!/usr/bin/env python3
"""
Debug UMAP: Generate individual UMAPs for each model
- Investigate 11 vs 12 batch issue
- Fix STATE and UCE skewed UMAPs
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.decomposition import PCA
import umap

import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300

SEVERITY_ORDER = ['healthy', 'mild', 'moderate', 'severe', 'fatal', 'recovered']

SEVERITY_COLORS = {
    'healthy': '#2ecc71',
    'mild': '#3498db',
    'moderate': '#f39c12',
    'severe': '#e74c3c',
    'fatal': '#8e44ad',
    'recovered': '#1abc9c',
}

# Batch ID to PMID mapping
BATCH_TO_PMID = {
    '302004data01': '32514174',  # Wilk et al.
    '302004data02': '32651212',  # Lee et al.
    '302004data04': '32810438',  # Schulte-Schrepping et al.
    '302004data05': '32788292',  # Arunachalam et al.
    '302004data15': '33657410',  # Ren et al.
    '302004data16': '33713619',  # Liu et al.
    '302004data17': '33879890',  # Stephenson et al.
    '302004data18': '32554932',  # Jiang et al.
    '302005data01': '33296687',  # Bernardes et al.
    '302005data02': '33494096',  # Combes et al.
    '302005data03': '34558746',  # Qiu et al.
}


def load_single_model(base_dir, model_name):
    """Load a single model's embeddings"""
    embedding_dir = Path(base_dir) / "02_EMBEDDINGS"

    configs = {
        'STATE': (embedding_dir / "state/covid_state.h5ad", ['X_emb', 'X_state']),
        'UCE': (embedding_dir / "uce/human_platelet_covid_severity_uce_adata.h5ad", ['X_uce', 'X_uce_4layer']),
        'scGPT': (embedding_dir / "scgpt/covid_scgpt_clean.h5ad", ['X_scGPT', 'X_scgpt']),
        'scGPT_BP': (embedding_dir / "scgpt_best_practices/covid_scgpt_best_practices.h5ad", ['X_scgpt_bp']),
        'TranscriptFormer': (embedding_dir / "transcriptformer/covid_transcriptformer.h5ad", ['X_transcriptformer']),
        'Geneformer': (embedding_dir / "geneformer/covid_geneformer.h5ad", ['X_geneformer']),
    }

    file_path, obsm_keys = configs[model_name]

    logger.info(f"Loading {model_name} from {file_path}")
    adata = sc.read_h5ad(file_path)

    # Find embedding key
    emb_key = None
    for key in obsm_keys:
        if key in adata.obsm:
            emb_key = key
            break
    if emb_key is None:
        for k in adata.obsm.keys():
            if 'emb' in k.lower() or model_name.lower().replace('_bp', '') in k.lower():
                emb_key = k
                break

    embeddings = np.array(adata.obsm[emb_key])

    # Get severity
    severity = None
    for col in ['covid_severity', 'severity', 'Severity', 'Category']:
        if col in adata.obs.columns:
            severity = np.array([str(x).lower() if pd.notna(x) else 'unknown' for x in adata.obs[col]])
            break

    # Get batch
    batch = None
    if 'Data_NO' in adata.obs.columns:
        batch = adata.obs['Data_NO'].values

    # Print detailed batch info
    logger.info(f"  Total cells: {len(adata)}")
    logger.info(f"  Embedding shape: {embeddings.shape}")

    if batch is not None:
        unique_batches = sorted(pd.Series(batch).dropna().unique())
        logger.info(f"  BATCH INFO:")
        logger.info(f"    Unique batches ({len(unique_batches)}): {unique_batches}")
        logger.info(f"    Batch value counts:")
        batch_counts = pd.Series(batch).value_counts().sort_index()
        for b, count in batch_counts.items():
            logger.info(f"      {b}: {count} cells")

    # Filter to valid severity
    valid_mask = np.isin(severity, SEVERITY_ORDER)
    logger.info(f"  Valid severity cells: {valid_mask.sum()}")

    return {
        'embeddings': embeddings[valid_mask],
        'severity': severity[valid_mask],
        'batch': batch[valid_mask] if batch is not None else None,
        'n_dims': embeddings.shape[1]
    }


def compute_umap_robust(embeddings, model_name, n_samples=10000, random_state=42):
    """Compute UMAP with model-specific robust preprocessing"""
    np.random.seed(random_state)

    logger.info(f"  Computing UMAP for {model_name}...")
    logger.info(f"    Input shape: {embeddings.shape}")

    # Check for NaN/Inf
    nan_count = np.isnan(embeddings).sum()
    inf_count = np.isinf(embeddings).sum()
    if nan_count > 0 or inf_count > 0:
        logger.warning(f"    Found {nan_count} NaN, {inf_count} Inf values!")
        embeddings = np.nan_to_num(embeddings, nan=0, posinf=0, neginf=0)

    # Step 1: PCA first for high-dimensional embeddings (STATE, UCE, TranscriptFormer)
    if model_name in ['STATE', 'UCE', 'TranscriptFormer'] and embeddings.shape[1] > 512:
        logger.info(f"    Applying PCA: {embeddings.shape[1]} -> 50 dims")
        pca = PCA(n_components=50, random_state=random_state)
        embeddings = pca.fit_transform(embeddings)
        logger.info(f"    PCA variance explained: {pca.explained_variance_ratio_.sum():.2%}")

    # Step 2: Robust outlier removal using Mahalanobis-like distance
    # Center and scale
    scaler = RobustScaler()
    X_scaled = scaler.fit_transform(embeddings)

    # Compute distances from median
    median = np.median(X_scaled, axis=0)
    distances = np.sqrt(np.sum((X_scaled - median)**2, axis=1))

    # Remove top 5% outliers for STATE/UCE, top 2% for others
    if model_name in ['STATE', 'UCE']:
        threshold = np.percentile(distances, 95)
    else:
        threshold = np.percentile(distances, 98)

    inlier_mask = distances <= threshold
    n_outliers = (~inlier_mask).sum()
    logger.info(f"    Removed {n_outliers} outliers ({n_outliers/len(embeddings)*100:.1f}%)")

    inlier_indices = np.where(inlier_mask)[0]

    # Step 3: Subsample
    if len(inlier_indices) > n_samples:
        idx = np.random.choice(inlier_indices, n_samples, replace=False)
    else:
        idx = inlier_indices

    X = X_scaled[idx]
    logger.info(f"    Final sample size: {len(X)}")

    # Step 4: UMAP with model-specific parameters
    if model_name in ['STATE', 'UCE']:
        # Spread out more, more global structure
        params = dict(n_neighbors=50, min_dist=0.5, spread=2.0, metric='euclidean')
    else:
        params = dict(n_neighbors=15, min_dist=0.1, spread=1.0, metric='euclidean')

    logger.info(f"    UMAP params: {params}")

    reducer = umap.UMAP(n_components=2, random_state=random_state, **params)
    coords = reducer.fit_transform(X)

    # Step 5: Center coordinates
    coords = coords - coords.mean(axis=0)

    return coords, idx


def plot_umap(coords, labels, label_type, colors_dict, model_name, output_dir, order=None, show_legend=True):
    """Plot single UMAP"""
    fig, ax = plt.subplots(figsize=(10, 8))

    if order is None:
        order = sorted(np.unique(labels))

    for label in order:
        if label in labels:
            mask = labels == label
            color = colors_dict.get(label, colors_dict.get(str(label), '#999999'))
            ax.scatter(coords[mask, 0], coords[mask, 1],
                      c=color, label=str(label), alpha=0.6, s=10, edgecolors='none')

    ax.set_title(f'{model_name} - {label_type}', fontsize=14, fontweight='bold')
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')

    if show_legend:
        ax.legend(loc='best', markerscale=2, fontsize=8)

    plt.tight_layout()

    fname = f'umap_{model_name.lower()}_{label_type.lower().replace(" ", "_")}'
    plt.savefig(output_dir / f'{fname}.png')
    plt.savefig(output_dir / f'{fname}.pdf')
    plt.close()
    logger.info(f"    Saved: {fname}.png/pdf")


def main():
    base_dir = Path("/bigdata/godziklab/shared/Xinru/302006")
    output_dir = base_dir / "302006_updategit/scripts/12_umap_visualization/figures/individual_umaps"
    output_dir.mkdir(parents=True, exist_ok=True)

    models = ['STATE', 'UCE', 'Geneformer', 'TranscriptFormer', 'scGPT', 'scGPT_BP']

    # Generate batch color palette for 12 batches
    batch_cmap = plt.cm.get_cmap('tab20', 20)

    all_batches = set()

    for model_name in models:
        logger.info(f"\n{'='*60}")
        logger.info(f"Processing {model_name}")
        logger.info(f"{'='*60}")

        try:
            data = load_single_model(base_dir, model_name)

            # Track all batches across models
            if data['batch'] is not None:
                unique_b = set(pd.Series(data['batch']).dropna().unique())
                all_batches.update(unique_b)

            # Compute UMAP
            coords, idx = compute_umap_robust(data['embeddings'], model_name)

            # Get labels for sampled cells
            severity_sampled = data['severity'][idx]

            # Plot severity UMAP (legend only for STATE)
            show_sev_legend = (model_name == 'STATE')
            plot_umap(coords, severity_sampled, '6class_severity',
                     SEVERITY_COLORS, model_name, output_dir, SEVERITY_ORDER, show_legend=show_sev_legend)

            # Plot batch UMAP if available
            if data['batch'] is not None:
                batch_sampled = data['batch'][idx]

                # Convert batch IDs to PMIDs
                batch_pmids = np.array([BATCH_TO_PMID.get(str(b), str(b)) for b in batch_sampled])
                unique_pmids = sorted(set(batch_pmids))
                batch_colors = {pmid: batch_cmap(i % 20) for i, pmid in enumerate(unique_pmids)}

                logger.info(f"    Plotting {len(unique_pmids)} batches (as PMIDs): {unique_pmids}")

                # Only show legend for STATE batch UMAP
                show_batch_legend = (model_name == 'STATE')
                plot_umap(coords, batch_pmids, 'batch',
                         batch_colors, model_name, output_dir, unique_pmids, show_legend=show_batch_legend)

        except Exception as e:
            logger.error(f"Error processing {model_name}: {e}")
            import traceback
            traceback.print_exc()

    # Summary
    logger.info(f"\n{'='*60}")
    logger.info("SUMMARY")
    logger.info(f"{'='*60}")
    logger.info(f"All unique batches across models: {sorted(all_batches)}")
    logger.info(f"Total unique batches: {len(all_batches)}")
    logger.info(f"\nFigures saved to: {output_dir}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
