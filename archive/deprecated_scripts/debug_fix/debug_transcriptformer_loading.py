#!/usr/bin/env python3
"""
Debug what TranscriptFormer actually sees when loading the data.
This mimics the exact logic in transcriptformer/data/dataloader.py
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path

def main():
    base_dir = Path("/bigdata/godziklab/shared/Xinru/302006")
    input_file = base_dir / "01_DATA/transcriptformer_ready/hek293t_transcriptformer_ready.h5ad"

    print("=" * 80)
    print("DEBUGGING TRANSCRIPTFORMER DATA LOADING")
    print("=" * 80)

    # Load data exactly as TranscriptFormer does
    print(f"\nLoading: {input_file}")
    adata = sc.read_h5ad(input_file)
    print(f"Shape: {adata.shape}")

    # TranscriptFormer's load_gene_features function:
    # 1. Gets gene_col_name from args (we pass "ensembl_id")
    # 2. If gene_col_name is in adata.var.columns, uses that
    # 3. Otherwise uses adata.var_names

    gene_col_name = "ensembl_id"

    print(f"\n--- Checking gene_col_name='{gene_col_name}' ---")
    print(f"Is '{gene_col_name}' in adata.var.columns? {gene_col_name in adata.var.columns}")

    if gene_col_name in adata.var.columns:
        gene_names = adata.var[gene_col_name].values
        print(f"Using adata.var['{gene_col_name}']")
    else:
        gene_names = adata.var_names.values
        print(f"Using adata.var_names")

    print(f"\nGene names (first 10): {gene_names[:10].tolist()}")
    print(f"Gene names dtype: {gene_names.dtype}")
    print(f"Total genes: {len(gene_names)}")

    # TranscriptFormer strips version numbers
    print("\n--- Stripping version numbers (TranscriptFormer logic) ---")
    gene_names_stripped = np.array([str(g).split('.')[0] for g in gene_names])
    print(f"After stripping (first 10): {gene_names_stripped[:10].tolist()}")

    # Check for duplicates
    unique_genes, counts = np.unique(gene_names_stripped, return_counts=True)
    duplicates = unique_genes[counts > 1]

    print(f"\nUnique genes: {len(unique_genes)}")
    print(f"Duplicate IDs: {len(duplicates)}")

    if len(duplicates) > 0:
        print(f"\n!!! FOUND DUPLICATES !!!")
        print(f"Showing first 20 duplicates:")
        for dup in duplicates[:20]:
            indices = np.where(gene_names_stripped == dup)[0]
            original = [gene_names[i] for i in indices]
            print(f"  '{dup}': count={len(indices)}, original={original}")
    else:
        print("\n✓ No duplicates found - this should work!")

    # Also check var_names directly
    print("\n--- Also checking adata.var_names directly ---")
    var_names = adata.var_names.values
    var_names_stripped = np.array([str(v).split('.')[0] for v in var_names])
    unique_var, var_counts = np.unique(var_names_stripped, return_counts=True)
    var_dups = unique_var[var_counts > 1]
    print(f"var_names (first 10): {var_names[:10].tolist()}")
    print(f"var_names unique: {len(unique_var)}")
    print(f"var_names duplicates: {len(var_dups)}")

    if len(var_dups) > 0:
        print(f"var_names duplicate examples:")
        for dup in var_dups[:10]:
            indices = np.where(var_names_stripped == dup)[0]
            print(f"  '{dup}': count={len(indices)}")

    # Check the raw layer too if exists
    print("\n--- Checking adata.raw ---")
    if adata.raw is not None:
        print(f"adata.raw exists")
        print(f"adata.raw.var_names (first 10): {adata.raw.var_names[:10].tolist()}")
        raw_var_stripped = np.array([str(v).split('.')[0] for v in adata.raw.var_names])
        unique_raw, raw_counts = np.unique(raw_var_stripped, return_counts=True)
        raw_dups = unique_raw[raw_counts > 1]
        print(f"adata.raw.var_names unique: {len(unique_raw)}")
        print(f"adata.raw.var_names duplicates: {len(raw_dups)}")
    else:
        print("adata.raw is None")


if __name__ == "__main__":
    main()
