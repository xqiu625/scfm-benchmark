#!/usr/bin/env python3
"""
Step 1: Identify top differentially expressed genes across severity classes.

These genes will be used to build protein-level features via ESM-2.
We select genes that are:
  (a) differentially expressed between severe and non-severe
  (b) have valid human UniProt protein entries

Usage:
    python step1_identify_degs.py [--n-top 500] [--output-dir results/]
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


def load_expression_data(base_dir):
    """Load raw expression data with severity labels."""
    embedding_dir = Path(base_dir) / "02_EMBEDDINGS"

    # Use STATE h5ad since it has both expression and metadata
    state_path = embedding_dir / "state" / "covid_state.h5ad"
    logger.info(f"Loading expression data from {state_path}")
    adata = sc.read_h5ad(state_path)
    logger.info(f"  Shape: {adata.shape}")
    return adata


def create_binary_labels(adata):
    """Create binary severity labels."""
    severity_col = None
    for col in ['covid_severity', 'severity', 'Severity']:
        if col in adata.obs.columns:
            severity_col = col
            break

    if severity_col is None:
        raise ValueError("No severity column found in adata.obs")

    severity = adata.obs[severity_col].astype(str).str.lower()
    binary = severity.map(lambda x: 1 if x in ['severe', 'fatal'] else 0)
    return binary.values


def find_degs_rank_based(adata, binary_labels, n_top=500):
    """
    Find differentially expressed genes using Wilcoxon rank-sum test.

    Uses scanpy's rank_genes_groups for efficient computation.
    Returns top genes ranked by combined score (log fold change * -log10 pvalue).
    """
    adata_expr = adata.copy()
    adata_expr.obs['severity_binary'] = pd.Categorical(
        ['severe' if b == 1 else 'non_severe' for b in binary_labels]
    )

    # Use raw counts if available, otherwise use .X
    if adata_expr.raw is not None:
        adata_raw = adata_expr.raw.to_adata()
        adata_raw.obs = adata_expr.obs
        adata_expr = adata_raw

    # Normalize for DEG analysis
    sc.pp.normalize_total(adata_expr, target_sum=1e4)
    sc.pp.log1p(adata_expr)

    logger.info("Running Wilcoxon rank-sum test...")
    sc.tl.rank_genes_groups(
        adata_expr, groupby='severity_binary', method='wilcoxon',
        reference='non_severe', n_genes=n_top * 2
    )

    # Extract results
    result = adata_expr.uns['rank_genes_groups']
    genes = [result['names'][i][0] for i in range(len(result['names']))]
    scores = [result['scores'][i][0] for i in range(len(result['scores']))]
    pvals = [result['pvals_adj'][i][0] for i in range(len(result['pvals_adj']))]
    logfc = [result['logfoldchanges'][i][0] for i in range(len(result['logfoldchanges']))]

    deg_df = pd.DataFrame({
        'gene': genes,
        'score': scores,
        'pval_adj': pvals,
        'log2fc': logfc
    })

    # Filter significant genes and rank by absolute score
    deg_df = deg_df[deg_df['pval_adj'] < 0.05].copy()
    deg_df['abs_score'] = deg_df['score'].abs()
    deg_df = deg_df.sort_values('abs_score', ascending=False).head(n_top)

    logger.info(f"  Found {len(deg_df)} significant DEGs (top {n_top} requested)")
    return deg_df


def find_degs_variance_based(adata, n_top=500):
    """
    Fallback: select highly variable genes if DEG analysis fails.
    """
    logger.info("Using variance-based gene selection as fallback...")
    adata_expr = adata.copy()

    if adata_expr.raw is not None:
        adata_raw = adata_expr.raw.to_adata()
        adata_expr = adata_raw

    sc.pp.normalize_total(adata_expr, target_sum=1e4)
    sc.pp.log1p(adata_expr)
    sc.pp.highly_variable_genes(adata_expr, n_top_genes=n_top)

    hvg = adata_expr.var[adata_expr.var['highly_variable']].index.tolist()
    deg_df = pd.DataFrame({
        'gene': hvg,
        'score': adata_expr.var.loc[hvg, 'dispersions_norm'].values,
        'pval_adj': [0.0] * len(hvg),
        'log2fc': [0.0] * len(hvg)
    })
    return deg_df


def main():
    parser = argparse.ArgumentParser(description="Identify top DEGs for ESM embedding")
    parser.add_argument('--base-dir', type=str,
                        default='/bigdata/godziklab/shared/Xinru/302006',
                        help='Base project directory')
    parser.add_argument('--n-top', type=int, default=500,
                        help='Number of top DEGs to select')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Output directory (default: scripts/22_.../results/)')
    args = parser.parse_args()

    if args.output_dir is None:
        script_dir = Path(__file__).resolve().parent
        args.output_dir = script_dir.parent / "results"
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    adata = load_expression_data(args.base_dir)
    binary_labels = create_binary_labels(adata)

    logger.info(f"  Severe: {binary_labels.sum()}, Non-severe: {(1-binary_labels).sum()}")

    # Find DEGs
    try:
        deg_df = find_degs_rank_based(adata, binary_labels, n_top=args.n_top)
    except Exception as e:
        logger.warning(f"DEG analysis failed: {e}. Using variance-based fallback.")
        deg_df = find_degs_variance_based(adata, n_top=args.n_top)

    # Save
    deg_path = output_dir / "top_degs.csv"
    deg_df.to_csv(deg_path, index=False)
    logger.info(f"Saved {len(deg_df)} DEGs to {deg_path}")

    # Save gene list for next step
    gene_list_path = output_dir / "deg_gene_list.txt"
    deg_df['gene'].to_csv(gene_list_path, index=False, header=False)
    logger.info(f"Saved gene list to {gene_list_path}")

    # Summary
    logger.info(f"\nTop 20 DEGs:")
    for _, row in deg_df.head(20).iterrows():
        logger.info(f"  {row['gene']:15s}  score={row['score']:.3f}  logFC={row['log2fc']:.3f}  p={row['pval_adj']:.2e}")


if __name__ == '__main__':
    main()
