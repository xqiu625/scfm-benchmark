#!/usr/bin/env python3
"""
Batch Effect Analysis: 6 Models Comparison

Evaluates batch effect removal using multiple complementary metrics:
- Clustering quality: ARI, NMI, ASW (bio-conservation)
- Batch mixing: kBET, iLISI, batch ASW (batch-correction)

Models: STATE, UCE, Geneformer, TranscriptFormer, scGPT, scGPT_BP

Author: Project 302006
Date: January 2026

RETIRED from the active batch-effect workflow, 2026-10-09 -- not deleted,
not rerun, outputs not touched. `compute_batch_asw()` below computes one
GLOBAL `silhouette_score(X, batch)` with no biological stratification and
negates the raw value directly (not even `1 - abs(...)`) -- the same category
of defect fixed in `scripts/24_option_c_canonical_input/code/
embedding_quality_option_c.py::scib_batch_asw()` on 2026-10-08/09 (see
`docs/decisions.md`, 2026-10-08/09 entry, and `tests/test_batch_asw.py`).
Do not read any `Batch ASW` value this module produced as scIB-correct, and
do not call `compute_batch_asw()` from new code -- use
`embedding_quality_option_c.scib_batch_asw()` instead. This module is kept
for historical reference to whatever COVID-Option-C-era figures/tables it
already produced; it is not wired into the current 8-cohort rebuild.
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
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    silhouette_score,
    silhouette_samples
)
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import KMeans
import scipy.stats as stats

import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Publication settings
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['pdf.fonttype'] = 42

# Models
MODEL_NAMES = ['STATE', 'UCE', 'Geneformer', 'TranscriptFormer', 'scGPT', 'scGPT_BP']
MODEL_DISPLAY_NAMES = {
    'STATE': 'STATE',
    'UCE': 'UCE',
    'Geneformer': 'Geneformer',
    'TranscriptFormer': 'TranscriptFormer',
    'scGPT': 'scGPT',
    'scGPT_BP': 'scGPT (BP)'
}

SEVERITY_ORDER = ['healthy', 'mild', 'moderate', 'severe', 'fatal', 'recovered']


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
                        for k in adata.obsm.keys():
                            if 'emb' in k.lower() or model_name.lower().replace('_bp', '') in k.lower():
                                emb_key = k
                                break

                    if emb_key is None:
                        logger.warning(f"No embedding found for {model_name}")
                        continue

                    embeddings = np.array(adata.obsm[emb_key])

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

                    # Get batch labels (Data_NO)
                    batch_col = 'Data_NO'
                    if batch_col in adata.obs.columns:
                        batch = np.array([str(x) for x in adata.obs[batch_col].values])
                    else:
                        batch = np.array(['unknown'] * len(adata))

                    # Filter to valid severity labels
                    valid_mask = np.isin(severity, SEVERITY_ORDER)

                    data[model_name] = {
                        'embeddings': embeddings[valid_mask],
                        'severity': severity[valid_mask],
                        'batch': batch[valid_mask],
                        'n_dims': embeddings.shape[1]
                    }

                    logger.info(f"  {model_name}: {valid_mask.sum()} cells, {embeddings.shape[1]} dims")
                    logger.info(f"    Batches: {len(np.unique(batch[valid_mask]))}")
                    loaded = True
                    break

                except Exception as e:
                    logger.error(f"Error loading {model_name}: {e}")

        if not loaded:
            logger.warning(f"Could not load {model_name}")

    return data


def subsample_data(embeddings, severity, batch, n_samples=10000, random_state=42):
    """Subsample data for faster computation"""
    np.random.seed(random_state)
    if len(embeddings) > n_samples:
        idx = np.random.choice(len(embeddings), n_samples, replace=False)
        return embeddings[idx], severity[idx], batch[idx]
    return embeddings, severity, batch


def compute_silhouette_scores(X, labels, label_name="labels"):
    """Compute silhouette score for given labels"""
    le = LabelEncoder()
    labels_encoded = le.fit_transform(labels)

    # Need at least 2 unique labels
    if len(np.unique(labels_encoded)) < 2:
        return np.nan

    return silhouette_score(X, labels_encoded)


def compute_ari_nmi(X, true_labels, n_clusters=None):
    """
    Compute ARI and NMI by clustering embeddings and comparing to true labels.
    """
    le = LabelEncoder()
    true_encoded = le.fit_transform(true_labels)

    if n_clusters is None:
        n_clusters = len(np.unique(true_encoded))

    # Cluster the embeddings
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    pred_labels = kmeans.fit_predict(X)

    ari = adjusted_rand_score(true_encoded, pred_labels)
    nmi = normalized_mutual_info_score(true_encoded, pred_labels)

    return ari, nmi


def compute_kbet(X, batch_labels, k=50, n_samples=1000, random_state=42):
    """
    Compute kBET (k-nearest neighbor Batch Effect Test).

    kBET tests whether the local neighborhood composition matches
    the global batch composition. Lower rejection rate = better mixing.

    Returns: acceptance rate (1 - rejection_rate), higher is better
    """
    np.random.seed(random_state)

    le = LabelEncoder()
    batch_encoded = le.fit_transform(batch_labels)
    n_batches = len(np.unique(batch_encoded))

    # Global batch frequencies
    batch_counts = np.bincount(batch_encoded)
    global_freq = batch_counts / len(batch_encoded)

    # Subsample for speed
    if len(X) > n_samples:
        idx = np.random.choice(len(X), n_samples, replace=False)
        X_sub = X[idx]
        batch_sub = batch_encoded[idx]
    else:
        X_sub = X
        batch_sub = batch_encoded

    # Find k nearest neighbors
    nn = NearestNeighbors(n_neighbors=k+1, metric='euclidean')
    nn.fit(X)
    distances, indices = nn.kneighbors(X_sub)

    # Test each cell's neighborhood
    rejection_count = 0

    for i in range(len(X_sub)):
        # Get batch labels of neighbors (excluding self)
        neighbor_batches = batch_encoded[indices[i, 1:]]

        # Observed frequencies in neighborhood
        observed = np.bincount(neighbor_batches, minlength=n_batches) / k

        # Expected frequencies (global)
        expected = global_freq

        # Chi-squared test
        # Avoid division by zero
        expected_counts = expected * k
        observed_counts = observed * k

        # Only test if expected counts are sufficient
        if np.all(expected_counts > 0):
            chi2 = np.sum((observed_counts - expected_counts)**2 / expected_counts)
            p_value = 1 - stats.chi2.cdf(chi2, df=n_batches-1)

            if p_value < 0.05:
                rejection_count += 1

    rejection_rate = rejection_count / len(X_sub)
    acceptance_rate = 1 - rejection_rate

    return acceptance_rate


def compute_lisi(X, labels, perplexity=30, n_samples=1000, random_state=42):
    """
    Compute LISI (Local Inverse Simpson's Index).

    For batch labels: iLISI (integration LISI) - higher = better batch mixing
    For bio labels: cLISI (cell-type LISI) - lower = better bio conservation

    LISI measures the effective number of categories in local neighborhoods.
    """
    np.random.seed(random_state)

    le = LabelEncoder()
    labels_encoded = le.fit_transform(labels)
    n_labels = len(np.unique(labels_encoded))

    # Subsample for speed
    if len(X) > n_samples:
        idx = np.random.choice(len(X), n_samples, replace=False)
        X_sub = X[idx]
        labels_sub = labels_encoded[idx]
    else:
        X_sub = X
        labels_sub = labels_encoded

    # Number of neighbors based on perplexity
    k = min(3 * perplexity, len(X) - 1)

    # Find k nearest neighbors
    nn = NearestNeighbors(n_neighbors=k+1, metric='euclidean')
    nn.fit(X)
    distances, indices = nn.kneighbors(X_sub)

    # Compute LISI for each cell
    lisi_values = []

    for i in range(len(X_sub)):
        # Get distances and labels of neighbors
        neighbor_distances = distances[i, 1:]  # exclude self
        neighbor_labels = labels_encoded[indices[i, 1:]]

        # Convert distances to probabilities using Gaussian kernel
        sigma = np.median(neighbor_distances)
        if sigma == 0:
            sigma = 1.0

        probs = np.exp(-neighbor_distances**2 / (2 * sigma**2))
        probs = probs / probs.sum()

        # Compute label frequencies weighted by probabilities
        label_probs = np.zeros(n_labels)
        for j, label in enumerate(neighbor_labels):
            label_probs[label] += probs[j]

        # Simpson's index
        simpson = np.sum(label_probs**2)

        # Inverse Simpson's index (effective number of labels)
        if simpson > 0:
            lisi = 1 / simpson
        else:
            lisi = 1

        lisi_values.append(lisi)

    return np.mean(lisi_values)


def compute_batch_asw(X, batch_labels):
    """
    RETIRED 2026-10-09. Raises unconditionally -- this must not silently
    return the old, not-scIB-correct global/unstratified metric if some
    caller still imports it. Not scIB-correct: this computed one global
    silhouette_score(X, batch) with no biological stratification, so real
    cohort-level biology could be misread as a batch effect. Use
    `scripts/00_shared/code/batch_asw.py::scib_batch_asw(X, batch, bio)`
    instead (docs/decisions.md, 2026-10-08/09 entries; a 2026-10-09 review
    specifically required that the retired entry point "cannot silently
    emit the old misleading metric if invoked").
    """
    raise NotImplementedError(
        "compute_batch_asw() is retired and not scIB-correct (global, "
        "unstratified silhouette). Use scib_batch_asw(X, batch_labels, "
        "bio_labels) from scripts/00_shared/code/batch_asw.py instead -- it "
        "requires a biological stratification label, which this function "
        "never had, so the call sites above (analyze_all_models) also need "
        "a bio label wired in, not just a drop-in replacement."
    )


def analyze_all_models(data, n_samples=5000):
    """Compute all metrics for all models"""
    results = []

    for model_name in MODEL_NAMES:
        if model_name not in data:
            logger.warning(f"{model_name} not in data, skipping")
            continue

        logger.info(f"\n{'='*60}")
        logger.info(f"Analyzing {model_name}...")
        logger.info(f"{'='*60}")

        embeddings = data[model_name]['embeddings']
        severity = data[model_name]['severity']
        batch = data[model_name]['batch']

        # Subsample for computational efficiency
        X, sev, bat = subsample_data(embeddings, severity, batch, n_samples)

        # Scale embeddings
        X_scaled = StandardScaler().fit_transform(X)

        logger.info(f"  Computing metrics on {len(X)} cells...")
        logger.info(f"  Batch column: Data_NO with {len(np.unique(bat))} unique batches")
        logger.info(f"  Batches: {np.unique(bat)}")

        # 1. Bio-conservation metrics (clustering quality)
        logger.info("  Computing bio-conservation metrics...")

        # ASW for severity (biological signal)
        bio_asw = compute_silhouette_scores(X_scaled, sev)
        logger.info(f"    Bio ASW (severity): {bio_asw:.4f}")

        # ARI and NMI for severity clustering
        ari, nmi = compute_ari_nmi(X_scaled, sev)
        logger.info(f"    ARI (severity): {ari:.4f}")
        logger.info(f"    NMI (severity): {nmi:.4f}")

        # cLISI for severity (lower = better bio conservation)
        clisi = compute_lisi(X_scaled, sev, n_samples=min(1000, len(X)))
        logger.info(f"    cLISI (severity): {clisi:.4f}")

        # 2. Batch mixing metrics
        logger.info("  Computing batch mixing metrics...")

        # Batch ASW (negative, so higher = better mixing)
        batch_asw = compute_batch_asw(X_scaled, bat)
        logger.info(f"    Batch ASW (neg): {batch_asw:.4f}")

        # kBET acceptance rate (higher = better mixing)
        kbet = compute_kbet(X_scaled, bat, k=50, n_samples=min(1000, len(X)))
        logger.info(f"    kBET acceptance: {kbet:.4f}")

        # iLISI for batch (higher = better mixing)
        ilisi = compute_lisi(X_scaled, bat, n_samples=min(1000, len(X)))
        logger.info(f"    iLISI (batch): {ilisi:.4f}")

        results.append({
            'Model': model_name,
            'Bio_ASW': bio_asw,
            'ARI': ari,
            'NMI': nmi,
            'cLISI': clisi,
            'Batch_ASW_neg': batch_asw,
            'kBET': kbet,
            'iLISI': ilisi,
            'n_cells': len(X),
            'n_batches': len(np.unique(bat))
        })

    return pd.DataFrame(results)


def create_summary_figure(df, output_dir):
    """Create summary visualization of all metrics"""
    logger.info("\nCreating summary figures...")

    # Prepare data for plotting
    bio_metrics = ['Bio_ASW', 'ARI', 'NMI']
    batch_metrics = ['Batch_ASW_neg', 'kBET', 'iLISI']

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # Color palette
    colors = sns.color_palette("husl", len(df))

    # Plot bio-conservation metrics
    for i, metric in enumerate(bio_metrics):
        ax = axes[0, i]
        bars = ax.bar(range(len(df)), df[metric], color=colors)
        ax.set_xticks(range(len(df)))
        ax.set_xticklabels([MODEL_DISPLAY_NAMES.get(m, m) for m in df['Model']],
                          rotation=45, ha='right')
        ax.set_ylabel(metric)
        ax.set_title(f'{metric}\n(Bio-conservation, higher=better)', fontweight='bold')
        ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)

        # Add value labels
        for bar, val in zip(bars, df[metric]):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                   f'{val:.3f}', ha='center', va='bottom', fontsize=8)

    # Plot batch mixing metrics
    for i, metric in enumerate(batch_metrics):
        ax = axes[1, i]
        bars = ax.bar(range(len(df)), df[metric], color=colors)
        ax.set_xticks(range(len(df)))
        ax.set_xticklabels([MODEL_DISPLAY_NAMES.get(m, m) for m in df['Model']],
                          rotation=45, ha='right')
        ax.set_ylabel(metric)
        ax.set_title(f'{metric}\n(Batch mixing, higher=better)', fontweight='bold')
        ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)

        # Add value labels
        for bar, val in zip(bars, df[metric]):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                   f'{val:.3f}', ha='center', va='bottom', fontsize=8)

    plt.suptitle('Batch Effect Analysis: Bio-Conservation vs Batch Mixing\n6 Foundation Models Comparison',
                fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    output_path = Path(output_dir)
    plt.savefig(output_path / 'batch_effect_metrics_6models.png')
    plt.savefig(output_path / 'batch_effect_metrics_6models.pdf')
    plt.close()
    logger.info(f"  Saved: batch_effect_metrics_6models.png/pdf")


def create_single_radar(df, metrics, title, filename, output_dir):
    """Create a single radar chart for given metrics"""
    df_norm = df.copy()

    for metric in metrics:
        min_val = df[metric].min()
        max_val = df[metric].max()
        if max_val > min_val:
            df_norm[metric] = (df[metric] - min_val) / (max_val - min_val)
        else:
            df_norm[metric] = 0.5

    # Radar chart
    angles = np.linspace(0, 2 * np.pi, len(metrics), endpoint=False).tolist()
    angles += angles[:1]  # Complete the loop

    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(polar=True))

    colors = sns.color_palette("husl", len(df))

    for i, (idx, row) in enumerate(df_norm.iterrows()):
        values = [row[m] for m in metrics]
        values += values[:1]  # Complete the loop

        ax.plot(angles, values, 'o-', linewidth=2, label=MODEL_DISPLAY_NAMES.get(row['Model'], row['Model']),
               color=colors[i])
        ax.fill(angles, values, alpha=0.1, color=colors[i])

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(metrics, fontsize=10)
    ax.set_ylim(0, 1)
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.0))

    plt.title(title, fontsize=14, fontweight='bold', y=1.08)

    output_path = Path(output_dir)
    plt.savefig(output_path / f'{filename}.png', bbox_inches='tight')
    plt.savefig(output_path / f'{filename}.pdf', bbox_inches='tight')
    plt.close()
    logger.info(f"  Saved: {filename}.png/pdf")


def create_radar_chart(df, output_dir):
    """Create separate radar charts for bio-conservation and batch mixing metrics"""
    logger.info("Creating radar charts...")

    # 1. Bio-conservation radar (severity-based metrics)
    # Note: cLISI measures local neighborhood purity, LOWER = better bio conservation
    # So we negate it to make it consistent with other metrics (higher = better)
    df_bio = df.copy()
    df_bio['cLISI_neg'] = -df_bio['cLISI']  # Negate so higher = better

    bio_metrics = ['Bio_ASW', 'ARI', 'NMI', 'cLISI_neg']
    create_single_radar(
        df_bio, bio_metrics,
        'Bio-Conservation Metrics (Severity-Based)\n(Higher = Better preservation of biological signal)',
        'radar_bio_conservation_6models',
        output_dir
    )

    # 2. Batch mixing radar (11 data sources-based metrics)
    batch_metrics = ['Batch_ASW_neg', 'kBET', 'iLISI']
    create_single_radar(
        df, batch_metrics,
        'Batch Mixing Metrics (11 Data Sources)\n(Higher = Better batch integration)',
        'radar_batch_mixing_6models',
        output_dir
    )

    # 3. Also create combined radar for reference
    # Use cLISI_neg for consistency (higher = better)
    all_metrics = ['Bio_ASW', 'ARI', 'NMI', 'cLISI_neg', 'Batch_ASW_neg', 'kBET', 'iLISI']
    create_single_radar(
        df_bio, all_metrics,
        'Model Comparison: All Metrics\n(Higher = Better for all metrics)',
        'batch_effect_radar_6models',
        output_dir
    )


def create_heatmap(df, output_dir):
    """Create heatmap of all metrics"""
    logger.info("Creating heatmap...")

    metrics = ['Bio_ASW', 'ARI', 'NMI', 'cLISI', 'Batch_ASW_neg', 'kBET', 'iLISI']

    # Create matrix
    matrix = df[metrics].values

    # Normalize each metric to 0-1
    matrix_norm = np.zeros_like(matrix)
    for i in range(matrix.shape[1]):
        min_val = matrix[:, i].min()
        max_val = matrix[:, i].max()
        if max_val > min_val:
            matrix_norm[:, i] = (matrix[:, i] - min_val) / (max_val - min_val)
        else:
            matrix_norm[:, i] = 0.5

    # For cLISI, lower is better for bio conservation, so invert
    clisi_idx = metrics.index('cLISI')
    matrix_norm[:, clisi_idx] = 1 - matrix_norm[:, clisi_idx]

    fig, ax = plt.subplots(figsize=(12, 8))

    # Create heatmap
    im = ax.imshow(matrix_norm, cmap='RdYlGn', aspect='auto', vmin=0, vmax=1)

    # Labels
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels(metrics, rotation=45, ha='right')
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels([MODEL_DISPLAY_NAMES.get(m, m) for m in df['Model']])

    # Add values
    for i in range(len(df)):
        for j in range(len(metrics)):
            text = ax.text(j, i, f'{matrix[i, j]:.3f}',
                          ha='center', va='center', color='black', fontsize=9)

    # Colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Normalized Score (higher = better)', fontsize=10)

    # Divider between bio and batch metrics
    ax.axvline(x=3.5, color='white', linewidth=3)

    # Add annotations
    ax.text(1.5, -1.2, 'Bio-Conservation', ha='center', fontsize=11, fontweight='bold')
    ax.text(5, -1.2, 'Batch Mixing', ha='center', fontsize=11, fontweight='bold')

    plt.title('Batch Effect Analysis: All Metrics Comparison\n(Green = Better Performance)',
             fontsize=14, fontweight='bold', pad=20)
    plt.tight_layout()

    output_path = Path(output_dir)
    plt.savefig(output_path / 'batch_effect_heatmap_6models.png')
    plt.savefig(output_path / 'batch_effect_heatmap_6models.pdf')
    plt.close()
    logger.info(f"  Saved: batch_effect_heatmap_6models.png/pdf")


def main():
    # Paths
    base_dir = Path("/bigdata/godziklab/shared/Xinru/302006")
    output_dir = base_dir / "302006_updategit/scripts/12_umap_visualization/figures/manuscript"
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("BATCH EFFECT ANALYSIS: 6 MODELS")
    logger.info("=" * 60)
    logger.info("Metrics:")
    logger.info("  Bio-conservation: ASW, ARI, NMI, cLISI")
    logger.info("  Batch mixing: Batch ASW (neg), kBET, iLISI")
    logger.info("=" * 60)

    # Load data
    logger.info("\nLoading embeddings...")
    data = load_embeddings(base_dir)

    if not data:
        logger.error("No data loaded!")
        return 1

    logger.info(f"\nLoaded {len(data)} models: {list(data.keys())}")

    # Compute all metrics
    df_results = analyze_all_models(data, n_samples=5000)

    # Save results
    logger.info("\n" + "=" * 60)
    logger.info("RESULTS SUMMARY")
    logger.info("=" * 60)
    print(df_results.to_string(index=False))

    # Save CSV
    csv_path = output_dir / 'batch_effect_metrics_6models.csv'
    df_results.to_csv(csv_path, index=False)
    logger.info(f"\nSaved: {csv_path}")

    # Create visualizations
    create_summary_figure(df_results, output_dir)
    create_radar_chart(df_results, output_dir)
    create_heatmap(df_results, output_dir)

    # Copy to main figures directory
    import shutil
    main_fig_dir = base_dir / "302006_updategit/figures"
    for fig_file in output_dir.glob('batch_effect_*.png'):
        shutil.copy(fig_file, main_fig_dir / fig_file.name)
    for fig_file in output_dir.glob('batch_effect_*.pdf'):
        shutil.copy(fig_file, main_fig_dir / fig_file.name)
    shutil.copy(csv_path, main_fig_dir / csv_path.name)

    logger.info("\n" + "=" * 60)
    logger.info("COMPLETE")
    logger.info(f"Results saved to: {output_dir}")
    logger.info("=" * 60)

    return 0


if __name__ == "__main__":
    sys.exit(main())
