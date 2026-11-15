#!/usr/bin/env python3
"""
Inspect metadata columns in embedding files
Find the correct column name for batch/source
"""

import scanpy as sc
from pathlib import Path
import pandas as pd

base_dir = Path("/bigdata/godziklab/shared/Xinru/302006")
embedding_dir = base_dir / "02_EMBEDDINGS"

files = {
    'STATE': embedding_dir / "state/covid_state.h5ad",
    'UCE': embedding_dir / "uce/human_platelet_covid_severity_uce_adata.h5ad",
    'scGPT': embedding_dir / "scgpt/covid_scgpt.h5ad",
    'TranscriptFormer': embedding_dir / "transcriptformer/covid_transcriptformer.h5ad",
    'Geneformer': embedding_dir / "geneformer/covid_geneformer.h5ad",
}

for model, filepath in files.items():
    print(f"\n{'='*60}")
    print(f"MODEL: {model}")
    print(f"File: {filepath}")
    print('='*60)

    if not filepath.exists():
        print("  FILE NOT FOUND")
        continue

    try:
        adata = sc.read_h5ad(filepath)
        print(f"\nShape: {adata.shape[0]} cells x {adata.shape[1]} genes")

        print(f"\n--- All .obs columns ({len(adata.obs.columns)}) ---")
        for col in adata.obs.columns:
            unique_vals = adata.obs[col].nunique()
            print(f"  {col}: {unique_vals} unique values")

        print(f"\n--- Sample values for potential batch/source columns ---")
        potential_cols = ['batch', 'Batch', 'source', 'Source', 'dataset', 'Dataset',
                         'sample', 'Sample', 'donor', 'Donor', 'patient', 'Patient',
                         'study', 'Study', 'orig.ident', 'library', 'Library']

        for col in adata.obs.columns:
            if any(p.lower() in col.lower() for p in ['batch', 'source', 'sample', 'donor',
                                                        'patient', 'study', 'library', 'orig', 'ident']):
                print(f"\n  Column: '{col}'")
                print(f"    Unique values: {adata.obs[col].unique()[:10].tolist()}")

        print(f"\n--- Sample of first 5 rows ---")
        print(adata.obs.head())

    except Exception as e:
        print(f"  ERROR: {e}")
