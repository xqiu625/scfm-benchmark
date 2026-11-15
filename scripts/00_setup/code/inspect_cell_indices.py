#!/usr/bin/env python3
"""
Inspect cell barcode indices across all 5 FM h5ad files.
Diagnoses why ensemble alignment finds 0 common cells.

Usage:
    python scripts/inspect_cell_indices.py
"""

import scanpy as sc
import numpy as np
from pathlib import Path

BASE_DIR = Path("/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS")

FILES = {
    'STATE':            BASE_DIR / "state/covid_state.h5ad",
    'UCE':              BASE_DIR / "uce/human_platelet_covid_severity_uce_adata.h5ad",
    'scGPT':            BASE_DIR / "scgpt/covid_scgpt.h5ad",
    'TranscriptFormer': BASE_DIR / "transcriptformer/covid_transcriptformer.h5ad",
    'Geneformer':       BASE_DIR / "geneformer/covid_geneformer.h5ad",
}

print("=" * 80)
print("CELL INDEX INSPECTION ACROSS 5 FOUNDATION MODELS")
print("=" * 80)

indices = {}

for name, path in FILES.items():
    print(f"\n--- {name} ---")
    if not path.exists():
        print(f"  FILE NOT FOUND: {path}")
        continue

    adata = sc.read_h5ad(path, backed='r')
    idx = adata.obs.index
    indices[name] = set(idx)

    print(f"  File: {path.name}")
    print(f"  Cells: {len(idx)}")
    print(f"  Index name: {idx.name}")
    print(f"  Index dtype: {idx.dtype}")
    print(f"  First 5: {list(idx[:5])}")
    print(f"  Last 5:  {list(idx[-5:])}")
    print(f"  Example format: '{idx[0]}'")

    # Check for common patterns
    sample = str(idx[0])
    has_dash = '-' in sample
    has_underscore = '_' in sample
    has_colon = ':' in sample
    print(f"  Contains dash: {has_dash}, underscore: {has_underscore}, colon: {has_colon}")

    # Check obs columns
    print(f"  obs columns: {list(adata.obs.columns)[:15]}")

    adata.file.close()

# Pairwise overlap
print("\n" + "=" * 80)
print("PAIRWISE BARCODE OVERLAP")
print("=" * 80)

models = list(indices.keys())
for i in range(len(models)):
    for j in range(i + 1, len(models)):
        m1, m2 = models[i], models[j]
        overlap = len(indices[m1] & indices[m2])
        print(f"  {m1} ∩ {m2}: {overlap} common cells")

# Full intersection
if len(indices) >= 2:
    full_intersection = set.intersection(*indices.values())
    print(f"\n  ALL {len(indices)} models intersection: {len(full_intersection)} common cells")

    if len(full_intersection) == 0:
        # Try substring matching (e.g., strip suffixes like -1, -2)
        print("\n--- Attempting fuzzy matching (strip trailing -N suffix) ---")
        stripped = {}
        for name, idx_set in indices.items():
            stripped[name] = set()
            for barcode in idx_set:
                base = str(barcode).rsplit('-', 1)[0]  # strip -1, -2, etc.
                stripped[name].add(base)

        for i in range(len(models)):
            for j in range(i + 1, len(models)):
                m1, m2 = models[i], models[j]
                overlap = len(stripped[m1] & stripped[m2])
                print(f"  {m1} ∩ {m2} (stripped): {overlap}")

        full_stripped = set.intersection(*stripped.values())
        print(f"  ALL models (stripped): {len(full_stripped)}")
