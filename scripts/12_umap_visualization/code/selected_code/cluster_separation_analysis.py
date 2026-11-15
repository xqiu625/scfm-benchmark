#!/usr/bin/env python3
"""
Cluster Separation Analysis: Visualizing Why Metrics Differ

This script visualizes why scGPT can have high ARI/NMI but low Bio_ASW:
- Shows cluster overlap vs separation for each model
- Computes inter-cluster vs intra-cluster distances
- Creates visualizations to illustrate the concept

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
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import silhouette_samples
from sklearn.decomposition import PCA
from scipy.spatial.distance import cdist
import umap

import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42

MODEL_NAMES = ['STATE', 'UCE', 'Geneformer', 'TranscriptFormer', 'scGPT', 'scGPT_BP']
SEVERITY_ORDER = ['healthy', 'mild', 'moderate', 'severe', 'fatal', 'recovered']
SEVERITY_COLORS = {
    'healthy': '#2ecc71',
    'mild': '#3498db',
    'moderate': '#f39c12',
    'severe': '#e74c3c',
    'fatal': '#8e44ad',
    'recovered': '#1abc9c'
}


def load_embeddings(base_dir):
    """Load all embedding files"""
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
        for file_path in config['files']:
            if file_path.exists():
                try:
                    adata = sc.read_h5ad(file_path)

                    emb_key = None
                    for key in config['obsm_keys']:
                        if key in adata.obsm:
                            emb_key = key
                            break

                    if emb_key is None:
                        continue

                    embeddings = np.array(adata.obsm[emb_key])

                    severity_col = None
                    for col in ['covid_severity', 'severity', 'Severity']:
                        if col in adata.obs.columns:
                            severity_col = col
                            break

                    if severity_col:
                        severity = np.array([str(x).lower() if pd.notna(x) else 'unknown'
                                           for x in adata.obs[severity_col]])
                    else:
                        continue

                    valid_mask = np.isin(severity, SEVERITY_ORDER)

                    data[model_name] = {
                        'embeddings': embeddings[valid_mask],
                        'severity': severity[valid_mask]
                    }

                    logger.info(f"Loaded {model_name}: {valid_mask.sum()} cells")
                    break

                except Exception as e:
                    logger.error(f"Error loading {model_name}: {e}")

    return data


def compute_cluster_distances(X, labels):
    """Compute inter-cluster and intra-cluster distances"""
    unique_labels = np.unique(labels)
    le = LabelEncoder()
    labels_encoded = le.fit_transform(labels)

    # Compute centroids
    centroids = {}
    for label in unique_labels:
        mask = labels == label
        centroids[label] = X[mask].mean(axis=0)

    # Intra-cluster distances (average distance to centroid within cluster)
    intra_distances = {}
    for label in unique_labels:
        mask = labels == label
        if mask.sum() > 0:
            dists = np.linalg.norm(X[mask] - centroids[label], axis=1)
            intra_distances[label] = np.mean(dists)

    # Inter-cluster distances (distances between centroids)
    centroid_matrix = np.array([centroids[l] for l in unique_labels])
    inter_distances = cdist(centroid_matrix, centroid_matrix)

    # Average inter-cluster distance (excluding diagonal)
    mask = ~np.eye(len(unique_labels), dtype=bool)
    avg_inter = inter_distances[mask].mean()

    # Average intra-cluster distance
    avg_intra = np.mean(list(intra_distances.values()))

    # Separation ratio: higher = better separated clusters
    separation_ratio = avg_inter / (avg_intra + 1e-10)

    return {
        'avg_inter': avg_inter,
        'avg_intra': avg_intra,
        'separation_ratio': separation_ratio,
        'intra_by_class': intra_distances
    }


def compute_silhouette_by_class(X, labels):
    """Compute silhouette scores for each class"""
    le = LabelEncoder()
    labels_encoded = le.fit_transform(labels)

    if len(np.unique(labels_encoded)) < 2:
        return {}

    sil_samples = silhouette_samples(X, labels_encoded)

    result = {}
    for label in np.unique(labels):
        mask = labels == label
        result[label] = sil_samples[mask].mean()

    return result


def create_umap_comparison(data, output_dir, n_samples=3000):
    """Create UMAP visualization showing cluster separation"""
    logger.info("Creating UMAP comparison...")

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()

    for idx, model_name in enumerate(MODEL_NAMES):
        if model_name not in data:
            continue

        ax = axes[idx]

        X = data[model_name]['embeddings']
        severity = data[model_name]['severity']

        # Subsample
        np.random.seed(42)
        if len(X) > n_samples:
            sample_idx = np.random.choice(len(X), n_samples, replace=False)
            X_sub = X[sample_idx]
            sev_sub = severity[sample_idx]
        else:
            X_sub = X
            sev_sub = severity

        # Scale and compute UMAP
        X_scaled = StandardScaler().fit_transform(X_sub)

        reducer = umap.UMAP(n_neighbors=30, min_dist=0.3, random_state=42)
        X_umap = reducer.fit_transform(X_scaled)

        # Plot each severity class
        for sev in SEVERITY_ORDER:
            mask = sev_sub == sev
            if mask.sum() > 0:
                ax.scatter(X_umap[mask, 0], X_umap[mask, 1],
                          c=SEVERITY_COLORS[sev], label=sev,
                          alpha=0.5, s=10, edgecolors='none')

        # Compute separation metrics
        dist_metrics = compute_cluster_distances(X_scaled, sev_sub)

        ax.set_title(f'{model_name}\nSeparation Ratio: {dist_metrics["separation_ratio"]:.2f}',
                    fontsize=12, fontweight='bold')
        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')
        ax.set_xticks([])
        ax.set_yticks([])

    # Add legend to last plot
    handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(handles, labels, loc='center right', bbox_to_anchor=(1.02, 0.5),
              fontsize=10, title='Severity')

    plt.suptitle('Cluster Separation Comparison Across Models\n(Higher separation ratio = better separated clusters)',
                fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 0.92, 0.95])

    output_path = Path(output_dir)
    plt.savefig(output_path / 'cluster_separation_umap.png', bbox_inches='tight')
    plt.savefig(output_path / 'cluster_separation_umap.pdf', bbox_inches='tight')
    plt.close()
    logger.info("  Saved: cluster_separation_umap.png/pdf")


def create_separation_metrics_plot(data, output_dir, n_samples=5000):
    """Create bar plot showing separation metrics"""
    logger.info("Computing separation metrics...")

    results = []

    for model_name in MODEL_NAMES:
        if model_name not in data:
            continue

        X = data[model_name]['embeddings']
        severity = data[model_name]['severity']

        # Subsample
        np.random.seed(42)
        if len(X) > n_samples:
            sample_idx = np.random.choice(len(X), n_samples, replace=False)
            X = X[sample_idx]
            severity = severity[sample_idx]

        X_scaled = StandardScaler().fit_transform(X)

        dist_metrics = compute_cluster_distances(X_scaled, severity)

        results.append({
            'Model': model_name,
            'Inter-cluster Distance': dist_metrics['avg_inter'],
            'Intra-cluster Distance': dist_metrics['avg_intra'],
            'Separation Ratio': dist_metrics['separation_ratio']
        })

        logger.info(f"  {model_name}: Inter={dist_metrics['avg_inter']:.3f}, "
                   f"Intra={dist_metrics['avg_intra']:.3f}, "
                   f"Ratio={dist_metrics['separation_ratio']:.3f}")

    df = pd.DataFrame(results)

    # Create figure with 3 subplots
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    colors = sns.color_palette("husl", len(df))

    # Plot 1: Inter-cluster distance
    ax = axes[0]
    bars = ax.bar(range(len(df)), df['Inter-cluster Distance'], color=colors)
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df['Model'], rotation=45, ha='right')
    ax.set_ylabel('Distance')
    ax.set_title('Inter-cluster Distance\n(Higher = clusters further apart)', fontweight='bold')
    for bar, val in zip(bars, df['Inter-cluster Distance']):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
               f'{val:.2f}', ha='center', va='bottom', fontsize=9)

    # Plot 2: Intra-cluster distance
    ax = axes[1]
    bars = ax.bar(range(len(df)), df['Intra-cluster Distance'], color=colors)
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df['Model'], rotation=45, ha='right')
    ax.set_ylabel('Distance')
    ax.set_title('Intra-cluster Distance\n(Lower = tighter clusters)', fontweight='bold')
    for bar, val in zip(bars, df['Intra-cluster Distance']):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
               f'{val:.2f}', ha='center', va='bottom', fontsize=9)

    # Plot 3: Separation ratio
    ax = axes[2]
    bars = ax.bar(range(len(df)), df['Separation Ratio'], color=colors)
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df['Model'], rotation=45, ha='right')
    ax.set_ylabel('Ratio')
    ax.set_title('Separation Ratio (Inter/Intra)\n(Higher = better separated)', fontweight='bold')
    for bar, val in zip(bars, df['Separation Ratio']):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
               f'{val:.2f}', ha='center', va='bottom', fontsize=9)

    plt.suptitle('Cluster Separation Analysis: Why Bio_ASW Differs from ARI/NMI',
                fontsize=14, fontweight='bold')
    plt.tight_layout()

    output_path = Path(output_dir)
    plt.savefig(output_path / 'cluster_separation_metrics.png', bbox_inches='tight')
    plt.savefig(output_path / 'cluster_separation_metrics.pdf', bbox_inches='tight')
    plt.close()
    logger.info("  Saved: cluster_separation_metrics.png/pdf")

    return df


def create_silhouette_distribution(data, output_dir, n_samples=3000):
    """Create violin plot showing silhouette score distribution by class"""
    logger.info("Creating silhouette distribution plot...")

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()

    for idx, model_name in enumerate(MODEL_NAMES):
        if model_name not in data:
            continue

        ax = axes[idx]

        X = data[model_name]['embeddings']
        severity = data[model_name]['severity']

        # Subsample
        np.random.seed(42)
        if len(X) > n_samples:
            sample_idx = np.random.choice(len(X), n_samples, replace=False)
            X = X[sample_idx]
            severity = severity[sample_idx]

        X_scaled = StandardScaler().fit_transform(X)

        le = LabelEncoder()
        labels_encoded = le.fit_transform(severity)

        sil_samples = silhouette_samples(X_scaled, labels_encoded)

        # Create dataframe for plotting
        plot_df = pd.DataFrame({
            'Severity': severity,
            'Silhouette': sil_samples
        })

        # Order by severity
        plot_df['Severity'] = pd.Categorical(plot_df['Severity'],
                                              categories=SEVERITY_ORDER,
                                              ordered=True)

        # Violin plot
        sns.violinplot(data=plot_df, x='Severity', y='Silhouette', ax=ax,
                      palette=SEVERITY_COLORS, order=SEVERITY_ORDER)

        ax.axhline(y=0, color='red', linestyle='--', alpha=0.5)
        ax.set_title(f'{model_name}\nMean Silhouette: {sil_samples.mean():.3f}',
                    fontsize=11, fontweight='bold')
        ax.set_xlabel('')
        ax.set_ylabel('Silhouette Score')
        ax.tick_params(axis='x', rotation=45)

    plt.suptitle('Silhouette Score Distribution by Severity Class\n'
                '(Negative values = cell is closer to other clusters than its own)',
                fontsize=14, fontweight='bold')
    plt.tight_layout()

    output_path = Path(output_dir)
    plt.savefig(output_path / 'silhouette_distribution.png', bbox_inches='tight')
    plt.savefig(output_path / 'silhouette_distribution.pdf', bbox_inches='tight')
    plt.close()
    logger.info("  Saved: silhouette_distribution.png/pdf")


def create_explanation_figure(output_dir):
    """Create a conceptual figure explaining the metric differences"""
    logger.info("Creating explanation figure...")

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    np.random.seed(42)

    # Scenario 1: Well-separated clusters (high ASW, high ARI)
    ax = axes[0]
    for i, (cx, cy, color) in enumerate([(-2, 0, '#e74c3c'), (2, 0, '#3498db'), (0, 3, '#2ecc71')]):
        x = np.random.normal(cx, 0.4, 100)
        y = np.random.normal(cy, 0.4, 100)
        ax.scatter(x, y, c=color, alpha=0.6, s=30, edgecolors='none')
    ax.set_title('Well-Separated Clusters\n(High ASW, High ARI/NMI)', fontweight='bold')
    ax.set_xlim(-5, 5)
    ax.set_ylim(-3, 6)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.text(0, -2.5, 'STATE, UCE pattern', ha='center', fontsize=10, style='italic')

    # Scenario 2: Overlapping but distinguishable (low ASW, high ARI)
    ax = axes[1]
    for i, (cx, cy, color) in enumerate([(-0.8, 0, '#e74c3c'), (0.8, 0, '#3498db'), (0, 1.5, '#2ecc71')]):
        x = np.random.normal(cx, 0.8, 100)
        y = np.random.normal(cy, 0.8, 100)
        ax.scatter(x, y, c=color, alpha=0.6, s=30, edgecolors='none')
    ax.set_title('Overlapping but Distinguishable\n(Low ASW, High ARI/NMI)', fontweight='bold')
    ax.set_xlim(-5, 5)
    ax.set_ylim(-3, 6)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.text(0, -2.5, 'scGPT pattern', ha='center', fontsize=10, style='italic')

    # Scenario 3: Mixed/confused clusters (low ASW, low ARI)
    ax = axes[2]
    colors = ['#e74c3c', '#3498db', '#2ecc71']
    for i in range(300):
        x = np.random.normal(0, 1.5)
        y = np.random.normal(0, 1.5)
        ax.scatter(x, y, c=colors[i % 3], alpha=0.6, s=30, edgecolors='none')
    ax.set_title('Mixed/Confused Clusters\n(Low ASW, Low ARI/NMI)', fontweight='bold')
    ax.set_xlim(-5, 5)
    ax.set_ylim(-3, 6)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.text(0, -2.5, 'Poor embedding pattern', ha='center', fontsize=10, style='italic')

    plt.suptitle('Why Metrics Can Disagree: Cluster Geometry Matters\n'
                'ASW measures separation, ARI/NMI measure assignment correctness',
                fontsize=14, fontweight='bold')
    plt.tight_layout()

    output_path = Path(output_dir)
    plt.savefig(output_path / 'metric_explanation.png', bbox_inches='tight')
    plt.savefig(output_path / 'metric_explanation.pdf', bbox_inches='tight')
    plt.close()
    logger.info("  Saved: metric_explanation.png/pdf")


def main():
    base_dir = Path("/bigdata/godziklab/shared/Xinru/302006")
    output_dir = base_dir / "302006_updategit/scripts/12_umap_visualization/figures/selected_figures"
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("CLUSTER SEPARATION ANALYSIS")
    logger.info("=" * 60)

    # Load data
    logger.info("\nLoading embeddings...")
    data = load_embeddings(base_dir)

    if not data:
        logger.error("No data loaded!")
        return 1

    # Create visualizations
    create_explanation_figure(output_dir)
    create_umap_comparison(data, output_dir)
    df_metrics = create_separation_metrics_plot(data, output_dir)
    create_silhouette_distribution(data, output_dir)

    # Save metrics
    csv_path = output_dir / 'cluster_separation_metrics.csv'
    df_metrics.to_csv(csv_path, index=False)

    logger.info("\n" + "=" * 60)
    logger.info("COMPLETE")
    logger.info(f"Results saved to: {output_dir}")
    logger.info("=" * 60)

    return 0


if __name__ == "__main__":
    sys.exit(main())
