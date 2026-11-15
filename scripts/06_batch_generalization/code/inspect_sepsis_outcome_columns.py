#!/usr/bin/env python3
"""Inspect the sepsis cells' obs columns to find the fatal-vs-survived outcome label.

The current `covid_severity` column for sepsis patients is all 'severe' (data
limitation), which makes binary outcome classification impossible. This script
scans every obs column for sepsis cells and surfaces likely candidates for a
fatal/survived outcome label.

Usage on HPCC:
    source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
    python scripts/06_batch_generalization/code/inspect_sepsis_outcome_columns.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import scanpy as sc

EMB = Path("/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/state/covid_state.h5ad")

KEYWORDS = (
    "surv", "fatal", "mort", "outcome", "death", "alive", "deceased",
    "status", "disposition", "30day", "90day", "discharge", "icu",
)


def main() -> int:
    print(f"Reading {EMB}")
    adata = sc.read_h5ad(EMB)

    if "Disease" in adata.obs.columns:
        disease_col = "Disease"
    elif "disease" in adata.obs.columns:
        disease_col = "disease"
    else:
        print("ERROR: No Disease/disease column on adata.obs")
        return 2

    sep_mask = adata.obs[disease_col].astype(str).str.contains(
        "sepsis", case=False, na=False
    )
    sepsis = adata[sep_mask]
    cov_mask = adata.obs[disease_col].astype(str).str.contains(
        "covid", case=False, na=False
    )
    covid = adata[cov_mask]

    print(f"COVID cells:  {covid.n_obs}")
    print(f"Sepsis cells: {sepsis.n_obs}")
    print(f"All obs columns ({len(adata.obs.columns)}): {list(adata.obs.columns)}")

    print("\n=== covid_severity distribution ===")
    if "covid_severity" in adata.obs.columns:
        print("COVID:  ", covid.obs["covid_severity"].value_counts().to_dict())
        print("Sepsis: ", sepsis.obs["covid_severity"].value_counts().to_dict())

    print("\n=== Outcome-keyword candidate columns (sepsis) ===")
    candidates = []
    for col in adata.obs.columns:
        col_lower = col.lower()
        if any(kw in col_lower for kw in KEYWORDS):
            vals = sepsis.obs[col].astype(str).value_counts().to_dict()
            n_unique = sepsis.obs[col].nunique(dropna=False)
            candidates.append(col)
            print(f"\n  [{col}]  unique={n_unique}")
            for k, v in list(vals.items())[:10]:
                print(f"    {k!r}: {v}")

    if not candidates:
        print("  (no columns matched keywords; printing top low-cardinality cols)")
        for col in adata.obs.columns:
            n_unique = sepsis.obs[col].nunique(dropna=False)
            if 2 <= n_unique <= 8:
                vals = sepsis.obs[col].astype(str).value_counts().to_dict()
                print(f"\n  [{col}]  unique={n_unique}")
                for k, v in list(vals.items())[:10]:
                    print(f"    {k!r}: {v}")

    print("\n=== Suggested next step ===")
    print("Pick the column above whose values look like fatal / survived, then run:")
    print("  bash scripts/06_batch_generalization/code/run_cross_disease_fatal.sh \\")
    print("    --sepsis-label-col <THAT_COLUMN> \\")
    print("    --sepsis-fatal-tokens 'died,fatal,deceased,non-survivor'")
    return 0


if __name__ == "__main__":
    sys.exit(main())
