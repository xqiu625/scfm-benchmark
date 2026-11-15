#!/usr/bin/env python3
"""
Fix HEK293T data by removing adata.raw which has mismatched gene names.

TranscriptFormer appears to be reading from adata.raw which has:
- Gene symbols instead of Ensembl IDs
- Different number of genes (24143 vs 24140)

Solution: Remove adata.raw entirely.
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path

def main():
    base_dir = Path("/bigdata/godziklab/shared/Xinru/302006")
    input_file = base_dir / "01_DATA/transcriptformer_ready/hek293t_transcriptformer_ready.h5ad"

    print("=" * 80)
    print("FIXING HEK293T - REMOVING adata.raw")
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

    # Remove adata.raw
    print("\nRemoving adata.raw...")
    adata.raw = None

    print(f"\nAfter fix:")
    print(f"  adata.raw exists: {adata.raw is not None}")

    # Verify ensembl_id column matches var_names
    if 'ensembl_id' in adata.var.columns:
        match = (adata.var_names == adata.var['ensembl_id']).all()
        print(f"  var_names == ensembl_id: {match}")

    # Save
    print(f"\nSaving to: {input_file}")
    adata.write_h5ad(input_file)

    print(f"\n✓ Done! adata.raw has been removed.")
    print(f"\nNow run:")
    print(f"  cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit")
    print(f"  python scripts/foundation_models/TranscriptFormer/generate_transcriptformer_embeddings.py --dataset hek293t")


if __name__ == "__main__":
    main()
