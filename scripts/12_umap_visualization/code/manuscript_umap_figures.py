#!/usr/bin/env python3
"""
Manuscript UMAP Figures: 6 Models x 3 Coloring Schemes

Creates publication-quality UMAP visualizations for:
- 5 original models: STATE, UCE, Geneformer, TranscriptFormer, scGPT
- 1 additional: scGPT with best practices (batch effect removal)

Three coloring schemes:
1. 6-class COVID severity
2. Binary severity (severe vs non-severe)
3. 12 batches (Data_NO)

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

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from sklearn.preprocessing import StandardScaler
import umap

import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Publication-quality settings
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['pdf.fonttype'] = 42  # TrueType fonts for PDF

# Color palettes
SEVERITY_COLORS = {
    'healthy': '#2ecc71',      # Green
    'mild': '#3498db',         # Blue
    'moderate': '#f39c12',     # Orange
    'severe': '#e74c3c',       # Red
    'fatal': '#8e44ad',        # Purple
    'recovered': '#1abc9c',    # Teal
    'unknown': '#95a5a6'       # Gray
}

SEVERITY_ORDER = ['healthy', 'mild', 'moderate', 'severe', 'fatal', 'recovered']

BINARY_COLORS = {
    'non_severe': '#3498db',   # Blue
    'severe': '#e74c3c'        # Red
}

# 6 models for manuscript
MODEL_NAMES = ['STATE', 'UCE', 'Geneformer', 'TranscriptFormer', 'scGPT', 'scGPT_BP']
MODEL_DISPLAY_NAMES = {
    'STATE': 'STATE',
    'UCE': 'UCE',
    'Geneformer': 'Geneformer',
    'TranscriptFormer': 'TranscriptFormer',
    'scGPT': 'scGPT',
    'scGPT_BP': 'scGPT (Best Practices)'
}


def load_embeddings(base_dir):
    """Load all embedding files including scGPT best practices"""
    embedding_dir = Path(base_dir) / "02_EMBEDDINGS"

    configs = {
        'STATE': {
            'files': [embedding_dir / "state/covid_state.h5ad"],
            'obsm_keys': ['X_emb', 'X_state']
        },
        'UCE': {
            'files': [embedding_dir / "uce/human_platelet_covid_severity_uce_adata.h5ad"],
            'obsm_keys': ['X_uce', 'X_uce_4layer']
        },
        'scGPT': {
            'files': [
                embedding_dir / "scgpt/covid_scgpt_clean.h5ad",
                embedding_dir / "scgpt/covid_scgpt.h5ad"
            ],
            'obsm_keys': ['X_scGPT', 'X_scgpt']
        },
        'scGPT_BP': {
            'files': [embedding_dir / "scgpt_best_practices/covid_scgpt_best_practices.h5ad"],
            'obsm_keys': ['X_scgpt_bp', 'X_scGPT_bp']
        },
        'TranscriptFormer': {
            'files': [embedding_dir / "transcriptformer/covid_transcriptformer.h5ad"],
            'obsm_keys': ['X_transcriptformer', 'X_TranscriptFormer']
        },
        'Geneformer': {
            'files': [embedding_dir / "geneformer/covid_geneformer.h5ad"],
            'obsm_keys': ['X_geneformer', 'X_Geneformer']
        }
    }

    data = {}

    for model_name, config in configs.items():
        loaded = False
        for file_path in config['files']:
            if file_path.exists():
                try:
                    logger.info(f"Loading {model_name} from {file_path}")
                    adata = sc.read_h5ad(file_path)

                    # Find embedding key
                    emb_key = None
                    for key in config['obsm_keys']:
                        if key in adata.obsm:
                            emb_key = key
                            break

                    if emb_key is None:
                        # Try to find any embedding key
                        for k in adata.obsm.keys():
                            if 'emb' in k.lower() or model_name.lower().replace('_bp', '') in k.lower():
                                emb_key = k
                                break

                    if emb_key is None:
                        logger.warning(f"No embedding found for {model_name}")
                        continue

                    embeddings = adata.obsm[emb_key]

                    # Get severity labels
                    severity_col = None
                    for col in ['covid_severity', 'severity', 'Severity', 'Category']:
                        if col in adata.obs.columns:
                            severity_col = col
                            break

                    if severity_col:
                        severity = np.array([str(x).lower() if pd.notna(x) else 'unknown'
                                           for x in adata.obs[severity_col]])
                    else:
                        severity = np.array(['unknown'] * len(adata))

                    # Get batch labels (Data_NO - 12 batches)
                    batch_col = 'Data_NO'
                    if batch_col in adata.obs.columns:
                        batch = adata.obs[batch_col].values
                    else:
                        batch = np.array(['unknown'] * len(adata))

                    # Filter to valid severity labels
                    valid_mask = np.isin(severity, SEVERITY_ORDER)

                    data[model_name] = {
                        'embeddings': embeddings[valid_mask],
                        'severity': severity[valid_mask],
                        'batch': batch[valid_mask] if len(batch) == len(severity) else batch[valid_mask],
                        'n_dims': embeddings.shape[1]
                    }

                    # Log batch information
                    unique_batches = np.unique(batch[valid_mask])
                    logger.info(f"  {model_name}: {valid_mask.sum()} cells, {embeddings.shape[1]} dims")
                    logger.info(f"    Batches ({len(unique_batches)}): {sorted(unique_batches)}")
                    loaded = True
                    break

                except Exception as e:
                    logger.error(f"Error loading {model_name}: {e}")

        if not loaded:
            logger.warning(f"Could not load {model_name}")

    return data


def compute_umap(embeddings, n_samples=10000, random_state=42, model_name=None):
    """Compute UMAP coordinates with subsampling and robust preprocessing"""
    np.random.seed(random_state)

    # Step 1: Winsorize - clip extreme values per dimension
    # Use tighter clipping for STATE and UCE which have more outliers
    if model_name in ['STATE', 'UCE']:
        lower_pct, upper_pct = 2, 98  # More aggressive clipping
    else:
        lower_pct, upper_pct = 1, 99

    embeddings_clipped = embeddings.copy()
    for dim in range(embeddings.shape[1]):
        lower = np.percentile(embeddings[:, dim], lower_pct)
        upper = np.percentile(embeddings[:, dim], upper_pct)
        embeddings_clipped[:, dim] = np.clip(embeddings[:, dim], lower, upper)

    # Step 2: Remove outliers based on L2 norm
    # More aggressive for STATE and UCE
    if model_name in ['STATE', 'UCE']:
        outlier_pct = 97  # Remove top 3%
    else:
        outlier_pct = 99  # Remove top 1%

    norms = np.linalg.norm(embeddings_clipped, axis=1)
    norm_threshold = np.percentile(norms, outlier_pct)
    inlier_mask = norms <= norm_threshold

    n_outliers = (~inlier_mask).sum()
    if n_outliers > 0:
        logger.info(f"    Removed {n_outliers} outlier cells ({n_outliers/len(embeddings)*100:.1f}%)")

    # Get indices of inliers
    inlier_indices = np.where(inlier_mask)[0]

    # Step 3: Subsample from inliers
    if len(inlier_indices) > n_samples:
        sampled_inlier_idx = np.random.choice(len(inlier_indices), n_samples, replace=False)
        idx = inlier_indices[sampled_inlier_idx]
    else:
        idx = inlier_indices

    X = embeddings_clipped[idx]

    # Step 4: Robust scaling using median and IQR
    from sklearn.preprocessing import RobustScaler
    X_scaled = RobustScaler().fit_transform(X)

    # Step 5: Compute UMAP with tuned parameters
    # Use different parameters for STATE/UCE to spread out more
    if model_name in ['STATE', 'UCE']:
        n_neighbors = 30  # More neighbors for smoother layout
        min_dist = 0.3    # More spread out
    else:
        n_neighbors = 15
        min_dist = 0.1

    reducer = umap.UMAP(
        n_components=2,
        random_state=random_state,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric='euclidean',
        spread=1.5  # Increase spread for better visualization
    )
    coords = reducer.fit_transform(X_scaled)

    # Step 6: Center and clip UMAP coordinates to remove visual outliers
    # Clip to 1st-99th percentile to keep main mass centered
    for dim in range(2):
        lower = np.percentile(coords[:, dim], 1)
        upper = np.percentile(coords[:, dim], 99)
        coords[:, dim] = np.clip(coords[:, dim], lower, upper)

    # Center the coordinates
    coords = coords - coords.mean(axis=0)

    return coords, idx


def get_binary_labels(labels):
    """Convert 6-class to binary labels"""
    binary = np.array(['non_severe'] * len(labels), dtype=object)
    binary[np.isin(labels, ['severe', 'fatal'])] = 'severe'
    return binary


def create_manuscript_figure_6class(data, output_dir, n_samples=10000):
    """
    Create 2x3 UMAP figure with 6-class severity coloring
    """
    logger.info("Creating 6-class severity UMAP figure...")

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for idx, model_name in enumerate(MODEL_NAMES):
        ax = axes[idx]

        if model_name not in data:
            ax.text(0.5, 0.5, f'{MODEL_DISPLAY_NAMES.get(model_name, model_name)}\nNot Available',
                   ha='center', va='center', fontsize=12)
            ax.set_axis_off()
            continue

        embeddings = data[model_name]['embeddings']
        severity = data[model_name]['severity']
        n_dims = data[model_name]['n_dims']

        coords, sample_idx = compute_umap(embeddings, n_samples, model_name=model_name)
        severity_sampled = severity[sample_idx]

        # Plot each severity level
        for sev in SEVERITY_ORDER:
            mask = severity_sampled == sev
            if mask.sum() > 0:
                ax.scatter(coords[mask, 0], coords[mask, 1],
                          c=SEVERITY_COLORS[sev], label=sev,
                          alpha=0.6, s=8, edgecolors='none')

        ax.set_title(f'{MODEL_DISPLAY_NAMES.get(model_name, model_name)}\n({n_dims} dims)',
                    fontsize=11, fontweight='bold')
        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')

        # Only show legend on last panel
        if idx == 5:
            ax.legend(loc='upper right', fontsize=8, markerscale=2, framealpha=0.9)

    plt.suptitle('UMAP Visualization: COVID-19 Severity (6-class)\nFoundation Model Comparison',
                fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    output_path = Path(output_dir)
    plt.savefig(output_path / 'manuscript_umap_6models_6class.png')
    plt.savefig(output_path / 'manuscript_umap_6models_6class.pdf')
    plt.close()
    logger.info("Saved: manuscript_umap_6models_6class.png/pdf")


def create_manuscript_figure_binary(data, output_dir, n_samples=10000):
    """
    Create 2x3 UMAP figure with binary severity coloring
    """
    logger.info("Creating binary severity UMAP figure...")

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for idx, model_name in enumerate(MODEL_NAMES):
        ax = axes[idx]

        if model_name not in data:
            ax.text(0.5, 0.5, f'{MODEL_DISPLAY_NAMES.get(model_name, model_name)}\nNot Available',
                   ha='center', va='center', fontsize=12)
            ax.set_axis_off()
            continue

        embeddings = data[model_name]['embeddings']
        severity = data[model_name]['severity']
        n_dims = data[model_name]['n_dims']

        coords, sample_idx = compute_umap(embeddings, n_samples, model_name=model_name)
        binary_labels = get_binary_labels(severity[sample_idx])

        # Plot binary classes
        for label in ['non_severe', 'severe']:
            mask = binary_labels == label
            if mask.sum() > 0:
                ax.scatter(coords[mask, 0], coords[mask, 1],
                          c=BINARY_COLORS[label], label=label.replace('_', '-'),
                          alpha=0.6, s=8, edgecolors='none')

        ax.set_title(f'{MODEL_DISPLAY_NAMES.get(model_name, model_name)}\n({n_dims} dims)',
                    fontsize=11, fontweight='bold')
        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')

        # Only show legend on last panel
        if idx == 5:
            ax.legend(loc='upper right', fontsize=9, markerscale=2, framealpha=0.9)

    plt.suptitle('UMAP Visualization: COVID-19 Severity (Binary: Severe vs Non-Severe)\nFoundation Model Comparison',
                fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    output_path = Path(output_dir)
    plt.savefig(output_path / 'manuscript_umap_6models_binary.png')
    plt.savefig(output_path / 'manuscript_umap_6models_binary.pdf')
    plt.close()
    logger.info("Saved: manuscript_umap_6models_binary.png/pdf")


def create_manuscript_figure_batch(data, output_dir, n_samples=10000):
    """
    Create 2x3 UMAP figure with 12-batch coloring
    """
    logger.info("Creating batch-colored UMAP figure...")

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    # Use a colormap for 12 batches
    batch_cmap = plt.cm.get_cmap('tab20', 12)

    for idx, model_name in enumerate(MODEL_NAMES):
        ax = axes[idx]

        if model_name not in data:
            ax.text(0.5, 0.5, f'{MODEL_DISPLAY_NAMES.get(model_name, model_name)}\nNot Available',
                   ha='center', va='center', fontsize=12)
            ax.set_axis_off()
            continue

        embeddings = data[model_name]['embeddings']
        batch = data[model_name]['batch']
        n_dims = data[model_name]['n_dims']

        coords, sample_idx = compute_umap(embeddings, n_samples, model_name=model_name)
        batch_sampled = batch[sample_idx]

        # Get unique batches
        unique_batches = sorted([b for b in np.unique(batch_sampled) if str(b) != 'unknown'])

        # Plot each batch
        for i, batch_name in enumerate(unique_batches):
            mask = batch_sampled == batch_name
            if mask.sum() > 0:
                ax.scatter(coords[mask, 0], coords[mask, 1],
                          c=[batch_cmap(i % 12)], label=str(batch_name)[:10],
                          alpha=0.6, s=8, edgecolors='none')

        ax.set_title(f'{MODEL_DISPLAY_NAMES.get(model_name, model_name)}\n({n_dims} dims)',
                    fontsize=11, fontweight='bold')
        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')

        # Add batch count annotation
        ax.text(0.02, 0.98, f'{len(unique_batches)} batches',
               transform=ax.transAxes, fontsize=9, va='top',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.suptitle('UMAP Visualization: Colored by Batch (Data_NO)\nFoundation Model Comparison',
                fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    output_path = Path(output_dir)
    plt.savefig(output_path / 'manuscript_umap_6models_batch.png')
    plt.savefig(output_path / 'manuscript_umap_6models_batch.pdf')
    plt.close()
    logger.info("Saved: manuscript_umap_6models_batch.png/pdf")


def create_combined_figure(data, output_dir, n_samples=10000):
    """
    Create a combined figure with all 3 coloring schemes (18 panels total)
    """
    logger.info("Creating combined manuscript figure...")

    fig = plt.figure(figsize=(18, 18))
    gs = gridspec.GridSpec(3, 6, figure=fig, hspace=0.3, wspace=0.25)

    batch_cmap = plt.cm.get_cmap('tab20', 12)

    # Row labels
    row_titles = ['6-class Severity', 'Binary Severity', 'Batch (Data_NO)']

    for row in range(3):
        for col, model_name in enumerate(MODEL_NAMES):
            ax = fig.add_subplot(gs[row, col])

            if model_name not in data:
                ax.text(0.5, 0.5, 'N/A', ha='center', va='center', fontsize=10)
                ax.set_axis_off()
                continue

            embeddings = data[model_name]['embeddings']
            severity = data[model_name]['severity']
            batch = data[model_name]['batch']
            n_dims = data[model_name]['n_dims']

            coords, sample_idx = compute_umap(embeddings, n_samples, model_name=model_name)

            if row == 0:  # 6-class severity
                severity_sampled = severity[sample_idx]
                for sev in SEVERITY_ORDER:
                    mask = severity_sampled == sev
                    if mask.sum() > 0:
                        ax.scatter(coords[mask, 0], coords[mask, 1],
                                  c=SEVERITY_COLORS[sev], label=sev,
                                  alpha=0.6, s=5, edgecolors='none')

            elif row == 1:  # Binary severity
                binary_labels = get_binary_labels(severity[sample_idx])
                for label in ['non_severe', 'severe']:
                    mask = binary_labels == label
                    if mask.sum() > 0:
                        ax.scatter(coords[mask, 0], coords[mask, 1],
                                  c=BINARY_COLORS[label], label=label.replace('_', '-'),
                                  alpha=0.6, s=5, edgecolors='none')

            else:  # Batch
                batch_sampled = batch[sample_idx]
                unique_batches = sorted([b for b in np.unique(batch_sampled) if str(b) != 'unknown'])
                for i, batch_name in enumerate(unique_batches):
                    mask = batch_sampled == batch_name
                    if mask.sum() > 0:
                        ax.scatter(coords[mask, 0], coords[mask, 1],
                                  c=[batch_cmap(i % 12)],
                                  alpha=0.6, s=5, edgecolors='none')

            # Title only on first row
            if row == 0:
                ax.set_title(f'{MODEL_DISPLAY_NAMES.get(model_name, model_name)}\n({n_dims}d)',
                            fontsize=10, fontweight='bold')

            ax.set_xticks([])
            ax.set_yticks([])

            # Row label on first column
            if col == 0:
                ax.set_ylabel(row_titles[row], fontsize=11, fontweight='bold')

    plt.suptitle('Foundation Model Embedding Comparison\nCOVID-19 Platelet Data',
                fontsize=16, fontweight='bold', y=0.98)

    output_path = Path(output_dir)
    plt.savefig(output_path / 'manuscript_umap_combined_6models_3schemes.png')
    plt.savefig(output_path / 'manuscript_umap_combined_6models_3schemes.pdf')
    plt.close()
    logger.info("Saved: manuscript_umap_combined_6models_3schemes.png/pdf")


def main():
    # Paths
    base_dir = Path("/bigdata/godziklab/shared/Xinru/302006")
    output_dir = base_dir / "302006_updategit/scripts/12_umap_visualization/figures/manuscript"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Also save to main figures directory
    main_fig_dir = base_dir / "302006_updategit/figures"
    main_fig_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("MANUSCRIPT UMAP FIGURES")
    logger.info("=" * 60)

    # Load data
    logger.info("\nLoading embeddings...")
    data = load_embeddings(base_dir)

    if not data:
        logger.error("No data loaded!")
        return 1

    logger.info(f"\nLoaded {len(data)} models: {list(data.keys())}")

    # Generate figures
    logger.info("\n" + "=" * 60)
    logger.info("Generating manuscript figures...")
    logger.info("=" * 60)

    # Individual figures
    create_manuscript_figure_6class(data, output_dir)
    create_manuscript_figure_binary(data, output_dir)
    create_manuscript_figure_batch(data, output_dir)

    # Combined figure
    create_combined_figure(data, output_dir)

    # Copy to main figures directory
    import shutil
    for fig_file in output_dir.glob('manuscript_*.pdf'):
        shutil.copy(fig_file, main_fig_dir / fig_file.name)
    for fig_file in output_dir.glob('manuscript_*.png'):
        shutil.copy(fig_file, main_fig_dir / fig_file.name)

    logger.info("\n" + "=" * 60)
    logger.info(f"All figures saved to:")
    logger.info(f"  {output_dir}")
    logger.info(f"  {main_fig_dir}")
    logger.info("=" * 60)

    return 0


if __name__ == "__main__":
    sys.exit(main())
