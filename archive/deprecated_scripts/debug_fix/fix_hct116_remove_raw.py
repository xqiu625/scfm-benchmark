#!/usr/bin/env python3
"""
Fix HCT116 data by removing adata.raw which has mismatched gene names.
Same fix as HEK293T.
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path

def main():
    base_dir = Path("/bigdata/godziklab/shared/Xinru/302006")
    input_file = base_dir / "01_DATA/transcriptformer_ready/hct116_transcriptformer_ready.h5ad"

    print("=" * 80)
    print("FIXING HCT116 - REMOVING adata.raw AND ALIGNING GENE IDs")
    print("=" * 80)

    # Load data
    print(f"\nLoading: {input_file}")
    adata = sc.read_h5ad(input_file)
    print(f"Shape: {adata.shape}")

    print(f"\nBefore fix:")
    print(f"  adata.var_names[:5]: {adata.var_names[:5].tolist()}")
    print(f"  adata.raw exists: {adata.raw is not None}")
    if adata.raw is not None:
        print(f"  adata.raw.shape: {adata.raw.shape}")
        print(f"  adata.raw.var_names[:5]: {adata.raw.var_names[:5].tolist()}")

    # Check if ensembl_id column exists and matches var_names
    if 'ensembl_id' in adata.var.columns:
        # Check if var_names are gene symbols (not Ensembl IDs)
        first_var = str(adata.var_names[0])
        if not first_var.startswith('ENSG'):
            print("\nvar_names are gene symbols, need to set to Ensembl IDs...")

            # Save original gene symbols
            adata.var['gene_symbol'] = adata.var_names.values

            # Get ensembl IDs and strip versions
            ensembl_ids = adata.var['ensembl_id'].values.astype(str)
            ensembl_ids_stripped = np.array([eid.split('.')[0] for eid in ensembl_ids])

            # Check for duplicates
            unique_ids, counts = np.unique(ensembl_ids_stripped, return_counts=True)
            n_dups = (counts > 1).sum()

            if n_dups > 0:
                print(f"  Found {n_dups} duplicate Ensembl IDs, removing...")
                _, unique_indices = np.unique(ensembl_ids_stripped, return_index=True)
                unique_indices = np.sort(unique_indices)
                adata = adata[:, unique_indices].copy()
                ensembl_ids_stripped = ensembl_ids_stripped[unique_indices]
                print(f"  Kept {len(unique_indices)} genes")

            # Set var_names to Ensembl IDs
            adata.var_names = pd.Index(ensembl_ids_stripped)
            adata.var['ensembl_id'] = ensembl_ids_stripped

    # Remove adata.raw
    print("\nRemoving adata.raw...")
    adata.raw = None

    print(f"\nAfter fix:")
    print(f"  adata.var_names[:5]: {adata.var_names[:5].tolist()}")
    print(f"  adata.raw exists: {adata.raw is not None}")

    # Verify ensembl_id column matches var_names
    if 'ensembl_id' in adata.var.columns:
        match = (adata.var_names == adata.var['ensembl_id']).all()
        print(f"  var_names == ensembl_id: {match}")

    # Save
    print(f"\nSaving to: {input_file}")
    adata.write_h5ad(input_file)

    print(f"\n✓ Done! HCT116 data fixed.")
    print(f"\nNow run:")
    print(f"  cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit")
    print(f"  python scripts/foundation_models/TranscriptFormer/generate_transcriptformer_embeddings.py --dataset hct116")


if __name__ == "__main__":
    main()
