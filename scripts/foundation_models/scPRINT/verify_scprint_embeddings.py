#!/usr/bin/env python3
"""Sanity-check a generated scPRINT embedding h5ad.

    python verify_scprint_embeddings.py \
        /bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/scprint/covid_scprint.h5ad
"""
from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np


def check(path: str) -> None:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)

    a = ad.read_h5ad(p)
    print(f"file:        {p}")
    print(f"n_obs:       {a.n_obs}")
    print(f"n_vars:      {a.n_vars}")
    print(f"obsm keys:   {list(a.obsm.keys())}")
    print(f"obs columns: {list(a.obs.columns)}")

    if "X_scprint" not in a.obsm:
        raise SystemExit("FAIL: X_scprint key missing from .obsm")

    emb = np.asarray(a.obsm["X_scprint"])
    print(f"X_scprint:   shape={emb.shape}, dtype={emb.dtype}")
    print(f"             min={emb.min():.4f}, max={emb.max():.4f}, mean={emb.mean():.4f}")

    if np.isnan(emb).any():
        raise SystemExit("FAIL: NaN in X_scprint")
    if np.isinf(emb).any():
        raise SystemExit("FAIL: Inf in X_scprint")
    if emb.shape[0] != a.n_obs:
        raise SystemExit(f"FAIL: row mismatch {emb.shape[0]} vs n_obs {a.n_obs}")

    for col in ("Sample_ID", "covid_severity"):
        if col in a.obs.columns:
            print(f"obs.{col} dist: {a.obs[col].value_counts().to_dict()}")

    print("OK")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: verify_scprint_embeddings.py <path.h5ad>")
    check(sys.argv[1])
