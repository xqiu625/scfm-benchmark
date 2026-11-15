#!/usr/bin/env python3
"""
Contamination audit of this project's own consensus attribution gene set
(results_gene_attribution_gsea/cross_model/consensus_gene_set.csv, k=5
unanimous-sign core, 531 genes) -- prompted by the same concern
subpopulation_biomarker_design.md raised for Qiu et al. 2024's biomarker
list, applied here to OUR data only (not Qiu 2024's).

13 of the 531 consensus-core genes are canonical non-platelet lineage
markers (B-cell: CD79A/CD79B; NK: PRF1; monocyte/myeloid: S100A8/S100A12/
FCER1G/TYROBP; erythroid: SLC25A37/HBA1; T-cell: CD3E/GZMA/CCL5; one
platelet-reference marker, GP1BA, included as a positive control). All the
POSITIVE-sign flagged genes are B/NK/myeloid/erythroid; all the
NEGATIVE-sign ones are T-cell -- a pattern consistent with EITHER (a) a
technical ambient-RNA artifact scaling with library depth, or (b) genuine
immune dysregulation (COVID lymphopenia + myeloid activation is real,
well-documented severity biology) showing up as differential leukocyte
contamination in the platelet gate. Distinguishing these requires checking
whether each gene's per-patient expression tracks QC covariates (library
depth, %mito) the way a pure technical artifact would.

Method: per-patient mean raw expression of each flagged gene (canonical
within-COVID cohort, same cells as the attribution pipeline), correlated
(Spearman) against the same per-patient QC metrics already established in
this project (patient_qc_distributions.csv: median_counts, median_genes,
median_pct_mito) and against the mortality label itself, so the "tracks QC"
and "tracks outcome" correlations can be compared directly for each gene.
"""
import os

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

RAW_H5AD = "/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/human_platelet_covid_severity_raw_counts.h5ad"
PATIENT_LIST_CSV = ("/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/"
                     "31_donor_level_revision/results_patient_v2/qc_baseline/patient_qc_distributions.csv")
SAMPLE_MAP_CSV = ("/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/"
                   "31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv")
OUT_DIR = "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/bib_revision_tier2/contamination_audit"
os.makedirs(OUT_DIR, exist_ok=True)

FLAGGED_GENES = {
    "CD79B": "B_cell", "CD79A": "B_cell",
    "PRF1": "NK_cell",
    "S100A8": "Monocyte_Myeloid", "S100A12": "Monocyte_Myeloid",
    "FCER1G": "Monocyte_Myeloid", "TYROBP": "Monocyte_Myeloid",
    "SLC25A37": "Erythroid", "HBA1": "Erythroid",
    "CD3E": "T_cell", "GZMA": "T_cell", "CCL5": "T_cell",
    "GP1BA": "Platelet_reference_positive_control",
}

CELLS_PER_PATIENT = 60
RNG_SEED = 20260910


def load_canonical_cohort():
    patient_list = pd.read_csv(PATIENT_LIST_CSV, dtype={"patient_id": str})
    patient_list = patient_list[patient_list["disease"] == "COVID-19"]
    assert len(patient_list) == 240 and int(patient_list["y"].sum()) == 25
    y_by_patient = patient_list.set_index("patient_id")["y"].to_dict()
    qc_by_patient = patient_list.set_index("patient_id")[
        ["median_counts", "median_genes", "median_pct_mito"]]
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
    return a, y_by_patient, qc_by_patient


a, y_by_patient, qc_by_patient = load_canonical_cohort()
var_names = a.var_names.astype(str).to_numpy()
gene_to_col = {g: i for i, g in enumerate(var_names)}

found = {g: gene_to_col[g] for g in FLAGGED_GENES if g in gene_to_col}
missing = set(FLAGGED_GENES) - set(found)
print(f"\n{len(found)}/{len(FLAGGED_GENES)} flagged genes found in var_names; missing: {missing}",
      flush=True)

X = a.X
patient_ids = a.obs["patient_id"].to_numpy()
unique_patients = sorted(set(patient_ids))

rows = []
for gene, col in found.items():
    col_data = X[:, col].toarray().ravel() if hasattr(X[:, col], "toarray") else np.asarray(X[:, col]).ravel()
    per_patient_mean = pd.Series(col_data, index=patient_ids).groupby(level=0).mean()
    per_patient_mean = per_patient_mean.reindex(unique_patients)

    y = pd.Series({p: y_by_patient[p] for p in unique_patients})
    qc = qc_by_patient.reindex(unique_patients)

    rho_y, p_y = spearmanr(per_patient_mean, y)
    rho_counts, p_counts = spearmanr(per_patient_mean, qc["median_counts"])
    rho_genes, p_genes = spearmanr(per_patient_mean, qc["median_genes"])
    rho_mito, p_mito = spearmanr(per_patient_mean, qc["median_pct_mito"])

    rows.append({
        "gene": gene, "lineage": FLAGGED_GENES[gene],
        "rho_vs_mortality": rho_y, "p_vs_mortality": p_y,
        "rho_vs_median_counts": rho_counts, "p_vs_median_counts": p_counts,
        "rho_vs_median_genes": rho_genes, "p_vs_median_genes": p_genes,
        "rho_vs_pct_mito": rho_mito, "p_vs_pct_mito": p_mito,
    })
    print(f"{gene:10s} ({FLAGGED_GENES[gene]:12s}): "
          f"rho_vs_mortality={rho_y:+.3f} (p={p_y:.3f})  "
          f"rho_vs_depth={rho_counts:+.3f} (p={p_counts:.3f})  "
          f"rho_vs_pct_mito={rho_mito:+.3f} (p={p_mito:.3f})", flush=True)

df = pd.DataFrame(rows).sort_values("rho_vs_mortality", ascending=False)
out_path = os.path.join(OUT_DIR, "consensus_gene_contamination_audit.csv")
df.to_csv(out_path, index=False)
print(f"\nSAVED {out_path}", flush=True)
