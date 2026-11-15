#!/usr/bin/env python3
"""
Export the canonical within-COVID cohort (same construction used throughout
this thread -- 240 patients/25 fatal, stratified 60 cells/patient, RNG_SEED
20260910) as plain matrix-market + CSV files for DecontX (R/celda), avoiding
any Python<->R anndata bridge dependency (zellkonverter etc.) that the
decontx_env may not have.

Outputs (OUT_DIR):
  counts.mtx        genes x cells, raw counts (matrix-market, celda/SCE
                    convention: rows=genes, cols=cells)
  genes.csv         var_names, one per row, same order as counts.mtx rows
  cells.csv         cell_id, patient_id, Sample_ID, y -- same order as
                    counts.mtx columns
"""
import os

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io as sio
import scipy.sparse as sp

RAW_H5AD = "/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/human_platelet_covid_severity_raw_counts.h5ad"
PATIENT_LIST_CSV = ("/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/"
                     "31_donor_level_revision/results_patient_v2/qc_baseline/patient_qc_distributions.csv")
SAMPLE_MAP_CSV = ("/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/"
                   "31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv")
OUT_DIR = "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/bib_revision_tier2/contamination_audit/decontx_io"
os.makedirs(OUT_DIR, exist_ok=True)

CELLS_PER_PATIENT = 60
RNG_SEED = 20260910


def load_canonical_cohort():
    patient_list = pd.read_csv(PATIENT_LIST_CSV, dtype={"patient_id": str})
    patient_list = patient_list[patient_list["disease"] == "COVID-19"]
    assert len(patient_list) == 240 and int(patient_list["y"].sum()) == 25
    y_by_patient = patient_list.set_index("patient_id")["y"].to_dict()
    covid_patient_ids = set(patient_list["patient_id"])

    sample_map = pd.read_csv(SAMPLE_MAP_CSV, dtype={"sample_id": str, "patient_id": str})
    sample_map = sample_map[sample_map["patient_id"].isin(covid_patient_ids)]
    sample_to_patient = sample_map.set_index("sample_id")["patient_id"].to_dict()

    print("Loading raw counts...", flush=True)
    a = sc.read_h5ad(RAW_H5AD)
    a.obs["patient_id"] = a.obs["Sample_ID"].astype(str).map(sample_to_patient)
    a = a[a.obs["patient_id"].notna()].copy()
    a.obs["y"] = a.obs["patient_id"].map(y_by_patient).astype(int)
    print(f"  {a.n_obs} cells, {a.obs['patient_id'].nunique()} patients", flush=True)

    rng = np.random.default_rng(RNG_SEED)
    keep_idx = []
    for pid, grp in a.obs.groupby("patient_id"):
        take = min(len(grp), CELLS_PER_PATIENT)
        keep_idx.extend(rng.choice(grp.index.to_numpy(), size=take, replace=False).tolist())
    a = a[keep_idx].copy()
    print(f"  after stratified subsample: {a.n_obs} cells", flush=True)
    a.obs_names_make_unique()
    return a


a = load_canonical_cohort()
X = a.X
if not sp.issparse(X):
    X = sp.csr_matrix(X)
X = X.tocsc().astype(np.float64)

# celda/SCE convention: genes (rows) x cells (columns)
sio.mmwrite(os.path.join(OUT_DIR, "counts.mtx"), X.T.tocsr())
pd.Series(a.var_names.astype(str)).to_csv(os.path.join(OUT_DIR, "genes.csv"), index=False, header=["gene"])
a.obs[["Sample_ID", "patient_id", "y"]].to_csv(os.path.join(OUT_DIR, "cells.csv"), index_label="cell_id")

print(f"\nSAVED to {OUT_DIR}:", flush=True)
print(f"  counts.mtx: {X.shape[1]} genes x {X.shape[0]} cells", flush=True)
print("  genes.csv, cells.csv", flush=True)
