#!/usr/bin/env python3
"""Normalize scPRINT-2 outputs to standard convention:
  - obs metadata copied from canonical input
  - embedding stored in obsm['X_scprint2']

scPRINT-2 output: embedding in obsm['scprint_emb']; obs partially preserved
                  but missing some columns; subset to ~47114 cells (dropped 863).
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import scanpy as sc

CANONICAL = Path("/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/human_platelet_covid_severity_raw_counts.h5ad")
OPT_C = Path("/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c")

print(f"Loading canonical: {CANONICAL}")
canon = sc.read_h5ad(CANONICAL)
print(f"  canonical: {canon.n_obs} cells, obs cols: {list(canon.obs.columns)[:5]}...")

# --- scPRINT-2 ---
scp_path = OPT_C / "scprint2_covid_canonical.h5ad"
print(f"\n--- scPRINT-2: {scp_path} ---")
scp = sc.read_h5ad(scp_path)
print(f"  scPRINT-2 pre: cells={scp.n_obs}, obsm keys: {list(scp.obsm.keys())}")

if "X_scprint2" not in scp.obsm:
    if "scprint_emb" not in scp.obsm:
        print(f"  scPRINT-2 ERROR: no scprint_emb found")
    else:
        # Align scPRINT-2 cells to canonical via index intersection
        common = scp.obs.index.intersection(canon.obs.index)
        print(f"  Common cells: {len(common)} (scprint2={scp.n_obs}, canonical={canon.n_obs})")
        if len(common) < 0.9 * scp.n_obs:
            print("  WARNING: less than 90% cell overlap - using scprint2's own obs")
            new = scp.copy()
            new.obsm["X_scprint2"] = np.asarray(scp.obsm["scprint_emb"], dtype=np.float32)
        else:
            scp_sub = scp[common]
            new = canon[common].copy()
            new.obsm["X_scprint2"] = np.asarray(scp_sub.obsm["scprint_emb"], dtype=np.float32)
        new.write_h5ad(scp_path, compression="gzip")
        print(f"  scPRINT-2 post: cells={new.n_obs}, X_scprint2.shape={new.obsm['X_scprint2'].shape}")

print("\nDone.")
