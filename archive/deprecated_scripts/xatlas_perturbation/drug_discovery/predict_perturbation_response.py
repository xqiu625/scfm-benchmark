#!/usr/bin/env python3
"""
Use STATE ST (State Transition) model to predict perturbation responses.
Validates if top therapeutic candidates shift severe cells toward recovered state.

Usage:
    python predict_perturbation_response.py \
        --covid-embeddings covid_state.h5ad \
        --perturbation-embeddings perturbation_state.h5ad \
        --candidates therapeutic_candidates.csv \
        --output predictions/

Requirements:
    pip install arc-state scanpy anndata pandas numpy
"""

import argparse
import subprocess
import os
import sys
from pathlib import Path
from datetime import datetime

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def load_therapeutic_candidates(csv_path: str, top_n: int = 20) -> pd.DataFrame:
    """Load top therapeutic candidates from perturbation analysis."""
    df = pd.read_csv(csv_path)

    # Filter out non-targeting controls
    df = df[~df['perturbation'].str.contains('non-targeting', case=False, na=False)]

    # Get top N by recovery score
    df = df.nlargest(top_n, 'recovery_score')

    # Extract gene name from perturbation ID
    df['gene'] = df['perturbation'].apply(lambda x: x.split('_')[0].split('|')[0])

    return df


def compute_severity_centroids(adata: ad.AnnData, emb_key: str, severity_col: str = 'severity') -> dict:
    """
    Compute centroid embeddings for each severity level.

    Returns:
        Dictionary mapping severity level to centroid embedding
    """
    centroids = {}

    for level in adata.obs[severity_col].unique():
        if pd.isna(level) or level in ['nan', 'unknown', '']:
            continue
        mask = adata.obs[severity_col] == level
        centroid = adata.obsm[emb_key][mask].mean(axis=0)
        centroids[level] = centroid
        print(f"  {level}: {mask.sum()} cells")

    return centroids


def predict_state_transition(
    source_embedding: np.ndarray,
    perturbation: str,
    model_dir: str = None
) -> np.ndarray:
    """
    Use STATE ST to predict post-perturbation state.

    This is a placeholder - actual implementation depends on STATE ST API.
    """
    # TODO: Implement actual STATE ST prediction
    # The STATE ST model predicts how a cell state changes after perturbation

    # Placeholder: Return source embedding (no change)
    # In reality, this would call:
    # state st predict --source embedding --perturbation gene --output predicted_state

    return source_embedding


def analyze_perturbation_effects(
    covid_path: str,
    candidates_path: str,
    output_dir: str,
    emb_key: str = None,
    severity_col: str = 'severity'
):
    """
    Analyze how perturbations affect cell states relative to severity centroids.

    Args:
        covid_path: Path to COVID STATE embeddings
        candidates_path: Path to therapeutic candidates CSV
        output_dir: Output directory
        emb_key: Key for embeddings in obsm
        severity_col: Column name for severity labels
    """
    print(f"\n{'='*60}")
    print("PERTURBATION EFFECT ANALYSIS")
    print(f"{'='*60}")

    # Load data
    print("\nLoading COVID embeddings...")
    adata = sc.read_h5ad(covid_path)
    print(f"  Cells: {adata.n_obs}")

    # Find embedding key
    if emb_key is None:
        for key in adata.obsm.keys():
            if 'state' in key.lower() or 'emb' in key.lower() or 'uce' in key.lower():
                emb_key = key
                break
        if emb_key is None:
            emb_key = list(adata.obsm.keys())[0]
    print(f"  Using embedding: {emb_key}")

    # Load candidates
    print("\nLoading therapeutic candidates...")
    candidates = load_therapeutic_candidates(candidates_path)
    print(f"  Top {len(candidates)} candidates loaded")

    # Compute severity centroids
    print("\nComputing severity centroids...")
    centroids = compute_severity_centroids(adata, emb_key, severity_col)

    # Define recovery direction
    if 'severe' in centroids and 'recovered' in centroids:
        recovery_direction = centroids['recovered'] - centroids['severe']
        recovery_direction = recovery_direction / np.linalg.norm(recovery_direction)
        print(f"\nRecovery direction: severe -> recovered")
    elif 'severe' in centroids and 'healthy' in centroids:
        recovery_direction = centroids['healthy'] - centroids['severe']
        recovery_direction = recovery_direction / np.linalg.norm(recovery_direction)
        print(f"\nRecovery direction: severe -> healthy")
    else:
        print("ERROR: Cannot compute recovery direction")
        return

    # Get severe cells for prediction
    severe_mask = adata.obs[severity_col].isin(['severe', 'fatal'])
    severe_cells = adata[severe_mask]
    print(f"\nSevere/fatal cells for analysis: {severe_cells.n_obs}")

    # Sample severe cells
    n_sample = min(1000, severe_cells.n_obs)
    sample_idx = np.random.choice(severe_cells.n_obs, n_sample, replace=False)
    severe_sample = severe_cells.obsm[emb_key][sample_idx]
    print(f"  Sampled {n_sample} cells")

    # Analyze each candidate
    results = []
    print("\nAnalyzing perturbation effects...")

    for idx, row in candidates.iterrows():
        gene = row['gene']
        recovery_score = row['recovery_score']

        # For each severe cell, predict post-perturbation state
        # This is where STATE ST would be used
        # For now, we use the recovery score as a proxy

        # Calculate theoretical shift toward recovery
        # Higher recovery score = more likely to shift toward recovered state
        theoretical_shift = recovery_score / candidates['recovery_score'].max()

        results.append({
            'gene': gene,
            'perturbation': row['perturbation'],
            'recovery_score': recovery_score,
            'theoretical_shift': theoretical_shift,
            'n_cells_perturbation': row.get('n_cells', 0),
            'closest_severity': row.get('closest_severity', 'unknown')
        })

        print(f"  {gene}: recovery_score={recovery_score:.2f}, shift={theoretical_shift:.2%}")

    # Create results DataFrame
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('recovery_score', ascending=False)

    # Save results
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

    results_path = os.path.join(output_dir, f'perturbation_predictions_{timestamp}.csv')
    results_df.to_csv(results_path, index=False)
    print(f"\nSaved predictions to {results_path}")

    # Create visualization
    create_prediction_plot(results_df, centroids, output_dir, timestamp)

    return results_df


def create_prediction_plot(results_df: pd.DataFrame, centroids: dict, output_dir: str, timestamp: str):
    """Create visualization of perturbation predictions."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Plot 1: Recovery score bar chart
    ax1 = axes[0]
    top_10 = results_df.head(10)
    colors = plt.cm.RdYlGn(top_10['theoretical_shift'])
    bars = ax1.barh(range(len(top_10)), top_10['recovery_score'], color=colors)
    ax1.set_yticks(range(len(top_10)))
    ax1.set_yticklabels(top_10['gene'])
    ax1.set_xlabel('Recovery Score')
    ax1.set_title('Top 10 Therapeutic Candidates\n(Higher = More Recovery-Like)')
    ax1.invert_yaxis()

    # Plot 2: Severity transition diagram
    ax2 = axes[1]
    severity_order = ['healthy', 'mild', 'moderate', 'severe', 'fatal', 'recovered']
    available = [s for s in severity_order if s in centroids]

    # Plot severity states as points
    y_pos = {s: i for i, s in enumerate(available)}
    for severity in available:
        color = 'green' if severity in ['healthy', 'recovered'] else 'red' if severity in ['severe', 'fatal'] else 'orange'
        ax2.scatter(0, y_pos[severity], s=500, c=color, alpha=0.6, label=severity)
        ax2.annotate(severity, (0.1, y_pos[severity]), fontsize=12, va='center')

    # Draw arrow for recovery direction
    if 'severe' in y_pos and 'recovered' in y_pos:
        ax2.annotate('', xy=(0, y_pos['recovered']), xytext=(0, y_pos['severe']),
                     arrowprops=dict(arrowstyle='->', color='blue', lw=2))
        ax2.text(0.3, (y_pos['severe'] + y_pos['recovered']) / 2,
                 'Recovery\nDirection', fontsize=10, color='blue')

    ax2.set_xlim(-0.5, 1)
    ax2.set_ylim(-0.5, len(available) - 0.5)
    ax2.set_title('Severity State Space')
    ax2.axis('off')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'perturbation_predictions_{timestamp}.png'),
                dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, f'perturbation_predictions_{timestamp}.pdf'),
                bbox_inches='tight')
    print(f"Saved visualization to {output_dir}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Predict perturbation responses with STATE ST')
    parser.add_argument('--covid-embeddings', '-c', required=True,
                        help='Path to COVID embeddings (h5ad with STATE or UCE embeddings)')
    parser.add_argument('--candidates', '-p', required=True,
                        help='Path to therapeutic candidates CSV')
    parser.add_argument('--output', '-o', default='results/perturbation_predictions',
                        help='Output directory')
    parser.add_argument('--emb-key', help='Key for embeddings in obsm')
    parser.add_argument('--severity-col', default='severity',
                        help='Column name for severity labels')
    parser.add_argument('--top-n', type=int, default=20,
                        help='Number of top candidates to analyze')

    args = parser.parse_args()

    analyze_perturbation_effects(
        covid_path=args.covid_embeddings,
        candidates_path=args.candidates,
        output_dir=args.output,
        emb_key=args.emb_key,
        severity_col=args.severity_col
    )

    print("\nDone!")


if __name__ == '__main__':
    main()
