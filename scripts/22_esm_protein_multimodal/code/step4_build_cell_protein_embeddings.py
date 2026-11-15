#!/usr/bin/env python3
"""
Step 4: Build per-cell protein embeddings from ESM + gene expression.

For each cell, aggregates protein ESM embeddings weighted by expression level:
    cell_protein_embedding = sum(expression_i * esm_embedding_i) / sum(expression_i)

This creates a per-cell representation that captures:
- WHICH proteins are expressed (from single-cell expression)
- WHAT those proteins look like structurally/functionally (from ESM-2)

Output: H5AD file with per-cell protein embeddings in obsm['X_esm_protein']

Usage:
    python step4_build_cell_protein_embeddings.py
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
import logging
import json

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def load_esm_embeddings(results_dir):
    """Load ESM protein embeddings."""
    emb_path = results_dir / "esm_protein_embeddings.npz"
    data = np.load(emb_path)
    embeddings = {gene: data[gene] for gene in data.files}
    logger.info(f"Loaded ESM embeddings for {len(embeddings)} proteins")
    emb_dim = next(iter(embeddings.values())).shape[0]
    logger.info(f"  Embedding dim: {emb_dim}")
    return embeddings, emb_dim


def load_expression_data(base_dir):
    """Load expression data (same h5ad as step 1)."""
    state_path = Path(base_dir) / "02_EMBEDDINGS" / "state" / "covid_state.h5ad"
    logger.info(f"Loading expression data from {state_path}")
    adata = sc.read_h5ad(state_path)
    return adata


def build_expression_weighted_embeddings(adata, esm_embeddings, esm_dim):
    """
    Build per-cell protein embeddings via expression-weighted aggregation.

    For each cell:
        e_cell = sum_g(expr_g * esm_g) / sum_g(expr_g)

    where g iterates over genes with ESM embeddings.
    """
    # Find overlap between adata genes and ESM genes
    adata_genes = set(adata.var_names)
    esm_genes = set(esm_embeddings.keys())
    overlap = sorted(adata_genes & esm_genes)
    logger.info(f"Gene overlap: {len(overlap)} / {len(esm_genes)} ESM genes found in expression data")

    if len(overlap) < 10:
        raise ValueError(f"Too few overlapping genes ({len(overlap)}). Check gene name format.")

    # Build ESM embedding matrix for overlapping genes: (n_genes, esm_dim)
    esm_matrix = np.array([esm_embeddings[g] for g in overlap])  # (n_overlap, 1280)

    # Extract expression for overlapping genes: (n_cells, n_genes)
    gene_idx = [list(adata.var_names).index(g) for g in overlap]

    # Get expression matrix (handle sparse)
    if hasattr(adata.X, 'toarray'):
        expr_matrix = adata.X[:, gene_idx].toarray()
    else:
        expr_matrix = adata.X[:, gene_idx]

    # Normalize expression per cell (so weights sum to 1)
    # Use raw counts or log-normalized — either works since we normalize
    expr_matrix = np.maximum(expr_matrix, 0)  # ensure non-negative
    row_sums = expr_matrix.sum(axis=1, keepdims=True)
    row_sums = np.maximum(row_sums, 1e-8)  # avoid division by zero
    expr_weights = expr_matrix / row_sums  # (n_cells, n_genes), rows sum to 1

    # Weighted average: (n_cells, n_genes) @ (n_genes, esm_dim) = (n_cells, esm_dim)
    logger.info(f"Computing expression-weighted protein embeddings...")
    logger.info(f"  Expression matrix: {expr_weights.shape}")
    logger.info(f"  ESM matrix: {esm_matrix.shape}")

    cell_protein_embeddings = expr_weights @ esm_matrix  # (n_cells, 1280)

    logger.info(f"  Output: {cell_protein_embeddings.shape}")
    logger.info(f"  Mean norm: {np.linalg.norm(cell_protein_embeddings, axis=1).mean():.4f}")
    logger.info(f"  Zero rows: {(np.linalg.norm(cell_protein_embeddings, axis=1) < 1e-6).sum()}")

    return cell_protein_embeddings, overlap


def main():
    parser = argparse.ArgumentParser(description="Build per-cell protein embeddings")
    parser.add_argument('--base-dir', type=str,
                        default='/bigdata/godziklab/shared/Xinru/302006')
    parser.add_argument('--output-dir', type=str, default=None)
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    results_dir = script_dir.parent / "results" if args.output_dir is None else Path(args.output_dir)

    # Load ESM embeddings
    esm_embeddings, esm_dim = load_esm_embeddings(results_dir)

    # Load expression data
    adata = load_expression_data(args.base_dir)

    # Build per-cell protein embeddings
    cell_protein_embs, overlap_genes = build_expression_weighted_embeddings(
        adata, esm_embeddings, esm_dim
    )

    # Save as H5AD with protein embeddings in obsm
    output_path = results_dir / "cell_protein_embeddings.h5ad"

    # Create minimal adata with just the embeddings
    import anndata
    adata_out = anndata.AnnData(
        obs=adata.obs.copy(),
        obsm={'X_esm_protein': cell_protein_embs}
    )
    adata_out.write_h5ad(output_path)
    logger.info(f"Saved cell protein embeddings to {output_path}")
    logger.info(f"  Shape: {cell_protein_embs.shape}")
    logger.info(f"  Genes used: {len(overlap_genes)}")

    # Also save as numpy for easy loading
    np.save(results_dir / "cell_protein_embeddings.npy", cell_protein_embs)

    # Save metadata
    meta = {
        'n_cells': cell_protein_embs.shape[0],
        'esm_dim': cell_protein_embs.shape[1],
        'n_genes_used': len(overlap_genes),
        'genes_used': overlap_genes
    }
    with open(results_dir / "cell_protein_embedding_meta.json", 'w') as f:
        json.dump(meta, f, indent=2)

    logger.info("\nDone. Ready for step 5 (multimodal fusion).")


if __name__ == '__main__':
    main()
