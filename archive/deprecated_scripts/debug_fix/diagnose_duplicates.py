#!/usr/bin/env python3
"""
Diagnose duplicate Ensembl IDs in prepared datasets
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
from collections import Counter

def diagnose_dataset(dataset_name):
    """Check for duplicates in a dataset"""
    print(f"\n{'='*80}")
    print(f"Diagnosing: {dataset_name}")
    print(f"{'='*80}")

    # Load data
    base_dir = Path("/bigdata/godziklab/shared/Xinru/302006/01_DATA/transcriptformer_ready")
    file_path = base_dir / f"{dataset_name}_transcriptformer_ready.h5ad"

    if not file_path.exists():
        print(f"File not found: {file_path}")
        return

    adata = sc.read_h5ad(file_path)
    print(f"Shape: {adata.shape}")

    # Get Ensembl IDs
    ensembl_ids = adata.var['ensembl_id'].values
    print(f"\nTotal genes: {len(ensembl_ids)}")

    # Check for None/NaN
    n_missing = sum(pd.isna(ensembl_ids))
    if n_missing > 0:
        print(f"WARNING: {n_missing} genes with missing Ensembl IDs")

    # Strip version numbers (what TranscriptFormer does)
    print("\nStripping version numbers (like TranscriptFormer)...")
    ensembl_ids_no_version = [str(eid).split('.')[0] if pd.notna(eid) else eid for eid in ensembl_ids]

    # Count occurrences
    counter = Counter(ensembl_ids_no_version)
    duplicates = {eid: count for eid, count in counter.items() if count > 1}

    if duplicates:
        print(f"\n❌ Found {len(duplicates)} duplicate Ensembl IDs:")
        print(f"   Total duplicate occurrences: {sum(duplicates.values())}")
        print("\nDuplicate examples (first 10):")
        for i, (eid, count) in enumerate(list(duplicates.items())[:10]):
            print(f"   {eid}: appears {count} times")

            # Show which genes have this ID
            indices = [i for i, x in enumerate(ensembl_ids_no_version) if x == eid]
            print(f"      Indices: {indices}")
            if hasattr(adata.var, 'gene_names'):
                gene_names = adata.var.iloc[indices]['gene_names'].tolist()
                print(f"      Gene names: {gene_names}")
    else:
        print("\n✓ No duplicate Ensembl IDs found")

    # Show sample IDs
    print(f"\nSample Ensembl IDs (first 10):")
    for i in range(min(10, len(ensembl_ids))):
        original = ensembl_ids[i]
        stripped = ensembl_ids_no_version[i]
        match = "✓" if original == stripped else "➜"
        print(f"   [{i}] {original} {match} {stripped}")

# Check all datasets
for dataset in ['covid', 'hek293t', 'hct116']:
    diagnose_dataset(dataset)

print(f"\n{'='*80}")
print("DIAGNOSIS COMPLETE")
print(f"{'='*80}")
