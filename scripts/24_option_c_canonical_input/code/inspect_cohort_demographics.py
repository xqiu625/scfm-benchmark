#!/usr/bin/env python3
"""Cohort + demographic audit on the canonical raw AnnData.

Reports what's already in adata.obs (Layer 1). Layers 2 (GEO sample metadata)
and 3 (source-paper demographics) require manual entry — see the
`cohort_metadata_template.csv` companion file.

Addresses CZI Virtual Cells Workshop §10.1 demographic-bias audit
(Fahsbender et al. 2025, arXiv 2507.10502v2).

Usage:
    python inspect_cohort_demographics.py
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("cohort_audit")

INPUT = Path(
    "/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/"
    "human_platelet_covid_severity_raw_counts.h5ad"
)
OUTPUT_DIR = Path(
    "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/cohort_audit"
)

DEMOGRAPHIC_KEYWORDS = (
    "sex", "gender", "age", "ancestry", "ethnic", "race", "donor",
)
COHORT_KEYWORDS = (
    "data_no", "study", "cohort", "batch", "dataset", "source", "gse", "ega", "pmid", "orig.ident",
)


def report_columns(adata, keywords, label):
    print(f"\n=== {label} columns in obs ===")
    found = False
    for col in adata.obs.columns:
        if any(kw in col.lower() for kw in keywords):
            found = True
            vc = adata.obs[col].value_counts(dropna=False)
            print(f"\n[{col}] ({len(vc)} unique values, {adata.obs[col].isna().sum()} NaN)")
            for k, v in list(vc.items())[:20]:
                print(f"  {k!r}: {v}")
            if len(vc) > 20:
                print(f"  ... and {len(vc) - 20} more")
    if not found:
        print(f"  (no {label.lower()} columns found by keyword match)")


def cohort_x_severity_crosstab(adata, cohort_col: str, severity_col: str) -> pd.DataFrame:
    if cohort_col not in adata.obs.columns or severity_col not in adata.obs.columns:
        return pd.DataFrame()
    df = pd.crosstab(adata.obs[cohort_col], adata.obs[severity_col], dropna=False)
    df["total_cells"] = df.sum(axis=1)
    return df.sort_values("total_cells", ascending=False)


def cohort_x_patient_count(adata, cohort_col: str, patient_col: str = "Sample_ID") -> pd.DataFrame:
    if cohort_col not in adata.obs.columns or patient_col not in adata.obs.columns:
        return pd.DataFrame()
    counts = (
        adata.obs.groupby(cohort_col)[patient_col]
        .nunique()
        .reset_index()
        .rename(columns={patient_col: "n_patients"})
    )
    cells = adata.obs[cohort_col].value_counts().reset_index()
    cells.columns = [cohort_col, "n_cells"]
    return counts.merge(cells, on=cohort_col).sort_values("n_cells", ascending=False)


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if not INPUT.exists():
        log.error("Canonical raw input not found: %s", INPUT)
        return 1

    log.info("Reading %s", INPUT)
    adata = sc.read_h5ad(INPUT)
    log.info("Shape: %d cells × %d genes", adata.n_obs, adata.n_vars)

    print(f"\nAll obs columns ({len(adata.obs.columns)}):")
    for col in adata.obs.columns:
        print(f"  {col}")

    report_columns(adata, COHORT_KEYWORDS, "Cohort/study")
    report_columns(adata, DEMOGRAPHIC_KEYWORDS, "Demographic")

    # If Data_NO present, build cohort × severity crosstab
    if "Data_NO" in adata.obs.columns and "covid_severity" in adata.obs.columns:
        cross = cohort_x_severity_crosstab(adata, "Data_NO", "covid_severity")
        cross.to_csv(OUTPUT_DIR / "cohort_x_severity.csv")
        print("\n=== Cohort × covid_severity ===")
        print(cross.to_string())
        log.info("Saved: %s", OUTPUT_DIR / "cohort_x_severity.csv")

    # Cohort × patient counts
    for cohort_col in ("Data_NO", "Disease", "Disease2"):
        if cohort_col in adata.obs.columns:
            tbl = cohort_x_patient_count(adata, cohort_col)
            if not tbl.empty:
                tbl.to_csv(OUTPUT_DIR / f"{cohort_col.lower()}_patient_counts.csv", index=False)
                print(f"\n=== {cohort_col} → patients & cells ===")
                print(tbl.to_string(index=False))

    # Generate the metadata template for manual completion
    if "Data_NO" in adata.obs.columns:
        cohorts = sorted(adata.obs["Data_NO"].dropna().unique())
        template = pd.DataFrame({
            "cohort_id": cohorts,
            "geo_or_ega_accession": ["TBD"] * len(cohorts),
            "source_paper_pmid": ["TBD"] * len(cohorts),
            "n_patients_in_obs": [
                adata.obs[adata.obs["Data_NO"] == c]["Sample_ID"].nunique()
                if "Sample_ID" in adata.obs.columns else "?"
                for c in cohorts
            ],
            "n_cells_in_obs": [
                int((adata.obs["Data_NO"] == c).sum()) for c in cohorts
            ],
            "sex_distribution_M_F_unknown": ["TBD,TBD,TBD"] * len(cohorts),
            "median_age_IQR": ["TBD"] * len(cohorts),
            "ancestry_distribution": ["TBD"] * len(cohorts),
            "country_or_site": ["TBD"] * len(cohorts),
            "notes": [""] * len(cohorts),
        })
        out = OUTPUT_DIR / "cohort_metadata_template.csv"
        template.to_csv(out, index=False)
        log.info("Generated metadata template: %s", out)
        print(f"\nManual entry template generated at:\n  {out}\n"
              f"Fill in TBD fields from each source paper / GEO record.\n"
              f"This becomes Supplementary Table 1.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
