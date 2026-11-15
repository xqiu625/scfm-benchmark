#!/usr/bin/env python3
"""
Properly fix HCT116 data for TranscriptFormer.

The issue: var_names and ensembl_id column are misaligned.
TranscriptFormer expects the gene column to match the matrix columns.

Solution: Set var_names directly to the ensembl_id values (in place).
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path

def main():
    base_dir = Path("/bigdata/godziklab/shared/Xinru/302006")
    input_file = base_dir / "01_DATA/transcriptformer_ready/hct116_transcriptformer_ready.h5ad"
    output_file = base_dir / "01_DATA/transcriptformer_ready/hct116_transcriptformer_fixed.h5ad"

    print("=" * 80)
    print("FIXING HCT116 DATA FOR TRANSCRIPTFORMER")
    print("=" * 80)

    # Load data
    print(f"\nLoading: {input_file}")
    adata = sc.read_h5ad(input_file)
    print(f"Shape: {adata.shape}")

    print(f"\nBefore fix:")
    print(f"  var_names[:5]: {adata.var_names[:5].tolist()}")
    print(f"  ensembl_id[:5]: {adata.var['ensembl_id'][:5].tolist()}")

    # The ensembl_id column IS the correct mapping for each gene
    # We need to set var_names to these Ensembl IDs

    # First, save original gene symbols/names
    adata.var['original_var_name'] = adata.var_names.values

    # Set var_names to the ensembl_id values
    ensembl_ids = adata.var['ensembl_id'].values.astype(str)

    # Strip version numbers
    ensembl_ids_stripped = np.array([eid.split('.')[0] for eid in ensembl_ids])

    # Check for duplicates BEFORE setting
    unique_ids, counts = np.unique(ensembl_ids_stripped, return_counts=True)
    n_dups = (counts > 1).sum()

    print(f"\nChecking for duplicates in ensembl_id column:")
    print(f"  Total IDs: {len(ensembl_ids_stripped)}")
    print(f"  Unique IDs: {len(unique_ids)}")
    print(f"  Duplicate IDs: {n_dups}")

    if n_dups > 0:
        print("\nRemoving duplicates (keeping first occurrence)...")
        _, unique_indices = np.unique(ensembl_ids_stripped, return_index=True)
        unique_indices = np.sort(unique_indices)
        adata = adata[:, unique_indices].copy()
        ensembl_ids_stripped = ensembl_ids_stripped[unique_indices]
        print(f"  Kept {len(unique_indices)} genes")

    # Set var_names to stripped Ensembl IDs
    adata.var_names = pd.Index(ensembl_ids_stripped)

    # Update ensembl_id column to match
    adata.var['ensembl_id'] = ensembl_ids_stripped

    print(f"\nAfter fix:")
    print(f"  var_names[:5]: {adata.var_names[:5].tolist()}")
    print(f"  ensembl_id[:5]: {adata.var['ensembl_id'][:5].tolist()}")

    # Final verification
    assert (adata.var_names == adata.var['ensembl_id']).all(), "var_names and ensembl_id should match!"

    # Check no duplicates in var_names
    assert len(adata.var_names) == len(set(adata.var_names)), "var_names should be unique!"

    print(f"\n✓ var_names and ensembl_id are now aligned and unique!")

    # Save
    print(f"\nSaving to: {output_file}")
    adata.write_h5ad(output_file)

    # Also overwrite the original file
    original_backup = input_file.with_suffix('.h5ad.bak')
    print(f"Backing up original to: {original_backup}")
    import shutil
    shutil.copy(input_file, original_backup)

    print(f"Overwriting original file: {input_file}")
    adata.write_h5ad(input_file)

    print(f"\n✓ Done! Both files saved.")
    print(f"\nNow run:")
    print(f"  python generate_transcriptformer_embeddings.py --dataset hct116")


if __name__ == "__main__":
    main()
