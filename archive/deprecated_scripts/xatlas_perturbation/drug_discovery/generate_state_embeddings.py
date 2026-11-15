#!/usr/bin/env python3
"""
Generate STATE embeddings for COVID severity data and compare with UCE.
STATE SE (State Embedding) model from Arc Institute.

Usage:
    python generate_state_embeddings.py --input covid_data.h5ad --output state_embeddings.h5ad

Requirements:
    pip install arc-state scanpy anndata
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


def check_state_installation():
    """Check if STATE is installed."""
    try:
        result = subprocess.run(['state', '--version'], capture_output=True, text=True)
        print(f"STATE version: {result.stdout.strip()}")
        return True
    except FileNotFoundError:
        print("ERROR: STATE not installed. Run: pip install arc-state")
        return False


def generate_state_embeddings(input_path: str, output_path: str, model_dir: str = None):
    """
    Generate STATE SE embeddings using CLI.

    Args:
        input_path: Path to input h5ad file
        output_path: Path to output h5ad file
        model_dir: Path to STATE model checkpoint (optional)
    """
    print(f"\n{'='*60}")
    print("STATE EMBEDDING GENERATION")
    print(f"{'='*60}")
    print(f"Input: {input_path}")
    print(f"Output: {output_path}")

    # Load data to check format
    print("\nLoading input data...")
    adata = sc.read_h5ad(input_path)
    print(f"  Cells: {adata.n_obs}")
    print(f"  Genes: {adata.n_vars}")

    # Check for required columns
    if 'gene_name' not in adata.var.columns and 'gene_symbol' not in adata.var.columns:
        # Use var_names as gene names
        adata.var['gene_name'] = adata.var_names
        print("  Added gene_name column from var_names")

    gene_col = 'gene_name' if 'gene_name' in adata.var.columns else 'gene_symbol'

    # Save temporary file if needed
    temp_input = input_path
    if 'gene_name' not in adata.var.columns:
        temp_input = input_path.replace('.h5ad', '_state_ready.h5ad')
        adata.write_h5ad(temp_input)
        print(f"  Saved STATE-ready file: {temp_input}")

    # Build STATE command
    cmd = [
        'state', 'emb', 'transform',
        '--input', temp_input,
        '--output', output_path,
        '--gene-column', gene_col,
    ]

    if model_dir:
        cmd.extend(['--model-folder', model_dir])

    print(f"\nRunning STATE command:")
    print(f"  {' '.join(cmd)}")

    # Run STATE
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print("\nSTATE output:")
        print(result.stdout)
        if result.stderr:
            print("Warnings:", result.stderr)
    except subprocess.CalledProcessError as e:
        print(f"\nERROR running STATE:")
        print(e.stderr)
        raise

    # Verify output
    if os.path.exists(output_path):
        adata_out = sc.read_h5ad(output_path)
        print(f"\n{'='*60}")
        print("OUTPUT VERIFICATION")
        print(f"{'='*60}")
        print(f"Cells: {adata_out.n_obs}")
        print(f"Embedding keys: {list(adata_out.obsm.keys())}")

        # Find STATE embedding key
        state_key = None
        for key in adata_out.obsm.keys():
            if 'state' in key.lower() or 'emb' in key.lower():
                state_key = key
                break

        if state_key:
            emb_shape = adata_out.obsm[state_key].shape
            print(f"STATE embedding shape: {emb_shape}")

        return adata_out
    else:
        print(f"ERROR: Output file not created: {output_path}")
        return None


def compare_embeddings(uce_path: str, state_path: str, output_dir: str):
    """
    Compare UCE and STATE embeddings.

    Args:
        uce_path: Path to UCE embeddings
        state_path: Path to STATE embeddings
        output_dir: Output directory for comparison results
    """
    import matplotlib.pyplot as plt
    from sklearn.metrics import silhouette_score
    from sklearn.preprocessing import StandardScaler

    print(f"\n{'='*60}")
    print("EMBEDDING COMPARISON: UCE vs STATE")
    print(f"{'='*60}")

    # Load embeddings
    adata_uce = sc.read_h5ad(uce_path)
    adata_state = sc.read_h5ad(state_path)

    # Find embedding keys
    uce_key = 'X_uce' if 'X_uce' in adata_uce.obsm else list(adata_uce.obsm.keys())[0]
    state_key = None
    for key in adata_state.obsm.keys():
        if 'state' in key.lower() or 'emb' in key.lower():
            state_key = key
            break
    if state_key is None:
        state_key = list(adata_state.obsm.keys())[0]

    print(f"UCE key: {uce_key}, shape: {adata_uce.obsm[uce_key].shape}")
    print(f"STATE key: {state_key}, shape: {adata_state.obsm[state_key].shape}")

    # Get severity labels
    severity_col = None
    for col in ['severity', 'disease_severity', 'covid_severity']:
        if col in adata_uce.obs.columns:
            severity_col = col
            break

    if severity_col is None:
        print("WARNING: No severity column found")
        return

    # Calculate silhouette scores
    labels = adata_uce.obs[severity_col].astype(str)
    valid_mask = ~labels.isin(['nan', 'unknown', ''])

    uce_emb = adata_uce.obsm[uce_key][valid_mask]
    labels_clean = labels[valid_mask]

    # Scale embeddings
    scaler = StandardScaler()
    uce_scaled = scaler.fit_transform(uce_emb)

    # Sample for faster computation
    n_sample = min(10000, len(uce_scaled))
    idx = np.random.choice(len(uce_scaled), n_sample, replace=False)

    uce_sil = silhouette_score(uce_scaled[idx], labels_clean.iloc[idx])
    print(f"\nSilhouette Scores (severity separation):")
    print(f"  UCE: {uce_sil:.4f}")

    # If STATE has same cells, compute its score too
    if len(adata_state) == len(adata_uce):
        state_emb = adata_state.obsm[state_key][valid_mask]
        state_scaled = scaler.fit_transform(state_emb)
        state_sil = silhouette_score(state_scaled[idx], labels_clean.iloc[idx])
        print(f"  STATE: {state_sil:.4f}")

    # Create UMAP comparison
    os.makedirs(output_dir, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # UCE UMAP
    sc.pp.neighbors(adata_uce, use_rep=uce_key)
    sc.tl.umap(adata_uce)
    sc.pl.umap(adata_uce, color=severity_col, ax=axes[0], show=False, title='UCE Embeddings')

    # STATE UMAP
    sc.pp.neighbors(adata_state, use_rep=state_key)
    sc.tl.umap(adata_state)
    sc.pl.umap(adata_state, color=severity_col, ax=axes[1], show=False, title='STATE Embeddings')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'uce_vs_state_umap.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'uce_vs_state_umap.pdf'), bbox_inches='tight')
    print(f"\nSaved comparison plot to {output_dir}")

    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Generate STATE embeddings')
    parser.add_argument('--input', '-i', required=True, help='Input h5ad file')
    parser.add_argument('--output', '-o', required=True, help='Output h5ad file')
    parser.add_argument('--model-dir', '-m', help='Path to STATE model checkpoint')
    parser.add_argument('--compare-uce', help='Path to UCE embeddings for comparison')
    parser.add_argument('--output-dir', default='results/state_comparison', help='Output directory')

    args = parser.parse_args()

    # Check installation
    if not check_state_installation():
        sys.exit(1)

    # Generate embeddings
    adata = generate_state_embeddings(args.input, args.output, args.model_dir)

    # Compare if UCE path provided
    if args.compare_uce and adata is not None:
        compare_embeddings(args.compare_uce, args.output, args.output_dir)

    print("\nDone!")


if __name__ == '__main__':
    main()
