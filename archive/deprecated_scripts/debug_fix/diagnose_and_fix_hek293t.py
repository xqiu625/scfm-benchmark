#!/usr/bin/env python3
"""
Diagnose and fix duplicate gene issues in HEK293T data for TranscriptFormer.

This script:
1. Loads the prepared HEK293T data
2. Checks for duplicates in the ensembl_id column
3. Shows exactly what TranscriptFormer sees
4. Creates a fixed version with guaranteed unique Ensembl IDs
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path

def main():
    base_dir = Path("/bigdata/godziklab/shared/Xinru/302006")
    input_file = base_dir / "01_DATA/transcriptformer_ready/hek293t_transcriptformer_ready.h5ad"
    output_file = base_dir / "01_DATA/transcriptformer_ready/hek293t_transcriptformer_ready_fixed.h5ad"

    print("=" * 80)
    print("DIAGNOSING HEK293T DATA FOR TRANSCRIPTFORMER")
    print("=" * 80)

    # Load data
    print(f"\nLoading: {input_file}")
    adata = sc.read_h5ad(input_file)
    print(f"Shape: {adata.shape}")

    # Check var structure
    print(f"\nvar columns: {list(adata.var.columns)}")
    print(f"var_names (first 10): {adata.var_names[:10].tolist()}")

    # Check ensembl_id column
    if 'ensembl_id' in adata.var.columns:
        ensembl_ids = adata.var['ensembl_id'].values
        print(f"\nensembl_id column (first 10): {ensembl_ids[:10].tolist()}")

        # Strip version numbers (this is what TranscriptFormer does)
        ensembl_ids_stripped = [str(eid).split('.')[0] for eid in ensembl_ids]

        # Check for duplicates
        unique_ids, counts = np.unique(ensembl_ids_stripped, return_counts=True)
        duplicates = unique_ids[counts > 1]

        print(f"\nTotal genes: {len(ensembl_ids)}")
        print(f"Unique Ensembl IDs (after version strip): {len(unique_ids)}")
        print(f"Number of duplicate IDs: {len(duplicates)}")

        if len(duplicates) > 0:
            print(f"\nDuplicate Ensembl IDs (showing first 20):")
            for dup_id in duplicates[:20]:
                indices = [i for i, eid in enumerate(ensembl_ids_stripped) if eid == dup_id]
                original_ids = [ensembl_ids[i] for i in indices]
                gene_symbols = [adata.var_names[i] for i in indices]
                print(f"  {dup_id}: indices={indices}, original={original_ids}, symbols={gene_symbols}")

    # Also check var_names for duplicates (TranscriptFormer might use these)
    print("\n" + "-" * 40)
    print("Checking var_names for duplicates:")
    var_names = adata.var_names.tolist()
    var_names_stripped = [str(v).split('.')[0] for v in var_names]
    unique_var, var_counts = np.unique(var_names_stripped, return_counts=True)
    var_duplicates = unique_var[var_counts > 1]

    print(f"Unique var_names (after version strip): {len(unique_var)}")
    print(f"Number of duplicate var_names: {len(var_duplicates)}")

    if len(var_duplicates) > 0:
        print(f"\nDuplicate var_names (showing first 20):")
        for dup_id in var_duplicates[:20]:
            indices = [i for i, v in enumerate(var_names_stripped) if v == dup_id]
            original_names = [var_names[i] for i in indices]
            print(f"  {dup_id}: count={len(indices)}, original={original_names[:5]}")

    # CREATE FIXED VERSION
    print("\n" + "=" * 80)
    print("CREATING FIXED VERSION")
    print("=" * 80)

    # Get ensembl IDs and strip versions
    ensembl_ids = adata.var['ensembl_id'].values
    ensembl_ids_stripped = np.array([str(eid).split('.')[0] for eid in ensembl_ids])

    # Find unique indices (keep first occurrence)
    _, unique_indices = np.unique(ensembl_ids_stripped, return_index=True)
    unique_indices = np.sort(unique_indices)  # Preserve original order

    print(f"\nKeeping {len(unique_indices)} unique genes (removing {len(ensembl_ids) - len(unique_indices)} duplicates)")

    # Subset to unique genes
    adata_fixed = adata[:, unique_indices].copy()

    # Update ensembl_id to stripped version
    adata_fixed.var['ensembl_id'] = ensembl_ids_stripped[unique_indices]

    # CRITICAL: Set var_names to the stripped Ensembl IDs
    # TranscriptFormer reads from var_names when checking for duplicates
    if 'gene_symbol' not in adata_fixed.var.columns:
        adata_fixed.var['gene_symbol'] = adata_fixed.var_names.values
    adata_fixed.var_names = pd.Index(adata_fixed.var['ensembl_id'].values)

    # Verify no duplicates
    final_var_names = adata_fixed.var_names.tolist()
    final_stripped = [str(v).split('.')[0] for v in final_var_names]
    final_unique = len(set(final_stripped))

    print(f"\nFinal check:")
    print(f"  Total genes: {len(final_var_names)}")
    print(f"  Unique after version strip: {final_unique}")
    print(f"  Duplicates: {len(final_var_names) - final_unique}")

    if len(final_var_names) == final_unique:
        print("\n✓ No duplicates - data is ready for TranscriptFormer!")
    else:
        print("\n✗ Still have duplicates - something is wrong")
        return

    # Save fixed version
    print(f"\nSaving fixed data to: {output_file}")
    adata_fixed.write_h5ad(output_file)

    print(f"\n✓ Fixed file saved!")
    print(f"\nNext step: Run embedding generation with the fixed file:")
    print(f"  python generate_transcriptformer_embeddings.py --dataset hek293t")
    print(f"\n  (Update the script to use the _fixed.h5ad file, or rename it)")


if __name__ == "__main__":
    main()
