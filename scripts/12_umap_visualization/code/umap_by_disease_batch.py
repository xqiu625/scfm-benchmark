#!/usr/bin/env python3
"""
UMAP Visualization by Disease and Data Source (Batch)

Investigates whether scGPT's crescent patterns correspond to
technical batches rather than biological signal.

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
import seaborn as sns
from sklearn.preprocessing import StandardScaler
import umap

import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 11
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['savefig.bbox'] = 'tight'

# Color palettes
DISEASE_COLORS = {
    'COVID-19': '#e74c3c',      # Red
    'Sepsis': '#3498db',        # Blue
    'SLE': '#9b59b6',           # Purple
    'Healthy': '#2ecc71',       # Green
    'Pneumonia': '#f39c12',     # Orange
    'Influenza': '#1abc9c',     # Teal
    'Non_COVID': '#95a5a6',     # Gray
    'unknown': '#7f8c8d'        # Dark gray
}

MODEL_NAMES = ['STATE', 'UCE', 'Geneformer', 'TranscriptFormer', 'scGPT']


def load_embeddings(base_dir):
    """Load all embedding files with metadata"""
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
                embedding_dir / "scgpt/covid_scgpt_clean.h5ad",  # Cleaned: platelet only
                embedding_dir / "scgpt/covid_scgpt.h5ad",
                embedding_dir / "scgpt/platelet_scgpt_combined.h5ad"
            ],
            'obsm_keys': ['X_scGPT', 'X_scgpt']
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
                        logger.warning(f"No embedding found for {model_name}")
                        continue

                    embeddings = adata.obsm[emb_key]

                    # Get disease labels
                    disease_col = None
                    for col in ['disease', 'Disease', 'condition', 'Condition', 'disease_state']:
                        if col in adata.obs.columns:
                            disease_col = col
                            break

                    if disease_col:
                        disease = adata.obs[disease_col].values
                    else:
                        disease = np.array(['unknown'] * len(adata))

                    # Get batch/source labels - use Data_NO (12 datasets)
                    batch_col = 'Data_NO'
                    if batch_col in adata.obs.columns:
                        batch = adata.obs[batch_col].values
                        logger.info(f"  Using batch column: {batch_col} ({len(np.unique(batch))} unique)")
                    else:
                        # Fallback to Sample_ID
                        if 'Sample_ID' in adata.obs.columns:
                            batch = adata.obs['Sample_ID'].values
                            logger.info(f"  Using batch column: Sample_ID")
                        else:
                            batch = np.array(['unknown'] * len(adata))

                    # Get severity labels
                    severity_col = None
                    for col in ['severity', 'Severity', 'severity_state', 'covid_severity']:
                        if col in adata.obs.columns:
                            severity_col = col
                            break

                    if severity_col:
                        severity = adata.obs[severity_col].values
                    else:
                        severity = np.array(['unknown'] * len(adata))

                    data[model_name] = {
                        'embeddings': embeddings,
                        'disease': disease,
                        'batch': batch,
                        'severity': severity,
                        'obs': adata.obs
                    }

                    logger.info(f"  {model_name}: {embeddings.shape[0]} cells, {embeddings.shape[1]} dims")
                    logger.info(f"  Disease labels: {np.unique(disease)[:5]}...")
                    logger.info(f"  Batch labels: {np.unique(batch)[:5]}...")

                    loaded = True
                    break

                except Exception as e:
                    logger.error(f"Error loading {model_name}: {e}")

        if not loaded:
            logger.warning(f"Could not load {model_name}")

    return data


def compute_umap(embeddings, n_samples=10000, random_state=42):
    """Compute UMAP coordinates with optional subsampling"""
    np.random.seed(random_state)

    if len(embeddings) > n_samples:
        idx = np.random.choice(len(embeddings), n_samples, replace=False)
    else:
        idx = np.arange(len(embeddings))

    X = embeddings[idx]
    X_scaled = StandardScaler().fit_transform(X)

    reducer = umap.UMAP(n_components=2, random_state=random_state, n_neighbors=30, min_dist=0.3)
    coords = reducer.fit_transform(X_scaled)

    return coords, idx


def fig_umap_by_disease(data, output_dir, n_samples=10000):
    """UMAP colored by disease for all 5 models"""
    logger.info("Creating UMAP by disease figure...")

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()

    for idx, model_name in enumerate(MODEL_NAMES):
        if model_name not in data:
            axes[idx].text(0.5, 0.5, f'{model_name}\nNot Available',
                          ha='center', va='center', fontsize=14)
            axes[idx].set_axis_off()
            continue

        embeddings = data[model_name]['embeddings']
        disease = data[model_name]['disease']

        coords, sample_idx = compute_umap(embeddings, n_samples)
        disease_sampled = disease[sample_idx]

        # Standardize disease names
        disease_std = []
        for d in disease_sampled:
            d_lower = str(d).lower()
            if 'covid' in d_lower:
                disease_std.append('COVID-19')
            elif 'sepsis' in d_lower:
                disease_std.append('Sepsis')
            elif 'sle' in d_lower or 'lupus' in d_lower:
                disease_std.append('SLE')
            elif 'healthy' in d_lower or 'control' in d_lower or d_lower == 'hc':
                disease_std.append('Healthy')
            elif 'pneumonia' in d_lower or 'lung_infect' in d_lower:
                disease_std.append('Pneumonia')
            elif 'flu' in d_lower or 'influenza' in d_lower:
                disease_std.append('Influenza')
            else:
                disease_std.append(str(d))
        disease_std = np.array(disease_std)

        # Plot
        ax = axes[idx]
        unique_diseases = np.unique(disease_std)

        for disease_name in unique_diseases:
            mask = disease_std == disease_name
            color = DISEASE_COLORS.get(disease_name, '#7f8c8d')
            ax.scatter(coords[mask, 0], coords[mask, 1],
                      c=color, label=disease_name, alpha=0.6, s=8)

        ax.set_title(f'{model_name}\n({embeddings.shape[1]} dims)', fontsize=14, fontweight='bold')
        ax.set_xlabel('Dimension 1')
        ax.set_ylabel('Dimension 2')
        ax.legend(loc='upper right', fontsize=8, markerscale=2)

    # Hide empty subplot
    axes[5].set_axis_off()

    plt.suptitle('UMAP Visualization: Colored by Disease\n(All 5 Foundation Models)',
                 fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_dir / 'umap_5models_by_disease.png')
    plt.savefig(output_dir / 'umap_5models_by_disease.pdf')
    plt.close()
    logger.info("Saved: umap_5models_by_disease.png/pdf")


def fig_umap_by_batch(data, output_dir, n_samples=10000):
    """UMAP colored by batch/source for all 5 models"""
    logger.info("Creating UMAP by batch figure...")

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()

    for idx, model_name in enumerate(MODEL_NAMES):
        if model_name not in data:
            axes[idx].text(0.5, 0.5, f'{model_name}\nNot Available',
                          ha='center', va='center', fontsize=14)
            axes[idx].set_axis_off()
            continue

        embeddings = data[model_name]['embeddings']
        batch = data[model_name]['batch']

        coords, sample_idx = compute_umap(embeddings, n_samples)
        batch_sampled = batch[sample_idx]

        # Plot
        ax = axes[idx]
        unique_batches = np.unique(batch_sampled)

        # Use a colormap for many batches
        if len(unique_batches) > 10:
            cmap = plt.cm.get_cmap('tab20', len(unique_batches))
            colors = [cmap(i) for i in range(len(unique_batches))]
        else:
            colors = plt.cm.tab10(np.linspace(0, 1, len(unique_batches)))

        for i, batch_name in enumerate(unique_batches):
            mask = batch_sampled == batch_name
            ax.scatter(coords[mask, 0], coords[mask, 1],
                      c=[colors[i]], label=str(batch_name)[:15], alpha=0.6, s=8)

        ax.set_title(f'{model_name}\n({embeddings.shape[1]} dims)', fontsize=14, fontweight='bold')
        ax.set_xlabel('Dimension 1')
        ax.set_ylabel('Dimension 2')

        # Only show legend if not too many batches
        if len(unique_batches) <= 10:
            ax.legend(loc='upper right', fontsize=7, markerscale=2)
        else:
            ax.text(0.02, 0.98, f'{len(unique_batches)} batches',
                   transform=ax.transAxes, fontsize=10, va='top')

    # Hide empty subplot
    axes[5].set_axis_off()

    plt.suptitle('UMAP Visualization: Colored by Batch/Source\n(All 5 Foundation Models)',
                 fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()

    output_dir = Path(output_dir)
    plt.savefig(output_dir / 'umap_5models_by_batch.png')
    plt.savefig(output_dir / 'umap_5models_by_batch.pdf')
    plt.close()
    logger.info("Saved: umap_5models_by_batch.png/pdf")


def fig_scgpt_detailed(data, output_dir, n_samples=10000):
    """Detailed scGPT analysis: side-by-side disease vs batch vs source_dataset"""
    if 'scGPT' not in data:
        logger.warning("scGPT data not available")
        return

    logger.info("Creating detailed scGPT figure...")

    embeddings = data['scGPT']['embeddings']
    obs = data['scGPT']['obs']

    coords, sample_idx = compute_umap(embeddings, n_samples)

    # Get relevant columns
    disease = obs['Disease'].values[sample_idx] if 'Disease' in obs.columns else np.array(['unknown'] * len(sample_idx))
    sample_id = obs['Sample_ID'].values[sample_idx] if 'Sample_ID' in obs.columns else np.array(['unknown'] * len(sample_idx))
    data_no = obs['Data_NO'].values[sample_idx] if 'Data_NO' in obs.columns else np.array(['unknown'] * len(sample_idx))
    source_dataset = obs['source_dataset'].values[sample_idx] if 'source_dataset' in obs.columns else np.array(['unknown'] * len(sample_idx))

    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    axes = axes.flatten()

    # Panel 1: By Disease
    ax = axes[0]
    for disease_name in np.unique(disease):
        if pd.isna(disease_name):
            continue
        mask = disease == disease_name
        color = DISEASE_COLORS.get(str(disease_name), '#7f8c8d')
        ax.scatter(coords[mask, 0], coords[mask, 1], c=color, label=str(disease_name), alpha=0.6, s=8)
    ax.set_title('Colored by Disease', fontsize=14, fontweight='bold')
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.legend(loc='upper right', fontsize=8, markerscale=2)

    # Panel 2: By Data_NO (12 datasets)
    ax = axes[1]
    unique_data_no = [d for d in np.unique(data_no) if not pd.isna(d)]
    cmap = plt.cm.get_cmap('tab20', len(unique_data_no))
    for i, dn in enumerate(unique_data_no):
        mask = data_no == dn
        ax.scatter(coords[mask, 0], coords[mask, 1], c=[cmap(i)], label=str(dn)[:15], alpha=0.6, s=8)
    ax.set_title('Colored by Data_NO (12 datasets)', fontsize=14, fontweight='bold')
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.legend(loc='upper right', fontsize=7, markerscale=2, ncol=2)

    # Panel 3: By source_dataset (HCT116, HEK293T, nan=platelet)
    ax = axes[2]
    source_colors = {'HCT116': '#e74c3c', 'HEK293T': '#3498db', 'platelet': '#2ecc71'}
    # Replace nan with 'platelet'
    source_clean = np.array(['platelet' if pd.isna(s) else str(s) for s in source_dataset])
    for src in np.unique(source_clean):
        mask = source_clean == src
        color = source_colors.get(src, '#7f8c8d')
        ax.scatter(coords[mask, 0], coords[mask, 1], c=color, label=src, alpha=0.6, s=8)
    ax.set_title('Colored by source_dataset\n(HCT116/HEK293T = cell lines, NOT platelets!)',
                 fontsize=14, fontweight='bold', color='red')
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.legend(loc='upper right', fontsize=10, markerscale=2)

    # Panel 4: Count summary
    ax = axes[3]
    ax.axis('off')

    # Create summary text
    total_cells = len(embeddings)
    platelet_cells = np.sum(source_clean == 'platelet') if 'source_dataset' in obs.columns else 'N/A'
    hct116_cells = np.sum(source_clean == 'HCT116') if 'source_dataset' in obs.columns else 'N/A'
    hek293t_cells = np.sum(source_clean == 'HEK293T') if 'source_dataset' in obs.columns else 'N/A'

    summary_text = f"""
    scGPT Embedding File Analysis
    ============================

    Total cells: {total_cells:,}

    By source_dataset:
    - Platelet (COVID study): {platelet_cells:,}
    - HCT116 (cell line): {hct116_cells:,}
    - HEK293T (cell line): {hek293t_cells:,}

    ⚠️  WARNING: scGPT file contains mixed data!
    The crescent patterns may be caused by
    non-platelet cell line data (HCT116, HEK293T)
    being included in the embedding file.

    Other models have ~47,000 cells (platelets only)
    scGPT has {total_cells:,} cells (mixed sources)
    """
    ax.text(0.1, 0.9, summary_text, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    plt.suptitle('scGPT Embedding Analysis: Why Do Crescents Appear?',
                 fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()

    output_dir = Path(output_dir)
    plt.savefig(output_dir / 'scgpt_detailed_analysis.png')
    plt.savefig(output_dir / 'scgpt_detailed_analysis.pdf')
    plt.close()
    logger.info("Saved: scgpt_detailed_analysis.png/pdf")


def main():
    base_dir = Path("/bigdata/godziklab/shared/Xinru/302006")
    output_dir = base_dir / "302006_updategit/scripts/12_umap_visualization/figures"
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Loading embeddings...")
    data = load_embeddings(base_dir)

    if not data:
        logger.error("No data loaded!")
        return 1

    logger.info(f"\nLoaded {len(data)} models: {list(data.keys())}")

    # Generate figures
    logger.info("\n" + "="*60)
    logger.info("Generating UMAP visualizations...")
    logger.info("="*60)

    fig_umap_by_disease(data, output_dir)
    fig_umap_by_batch(data, output_dir)
    fig_scgpt_detailed(data, output_dir)

    logger.info("\n" + "="*60)
    logger.info(f"All figures saved to: {output_dir}")
    logger.info("="*60)

    return 0


if __name__ == "__main__":
    sys.exit(main())
