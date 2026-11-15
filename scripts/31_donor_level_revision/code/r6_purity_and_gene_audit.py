#!/usr/bin/env python3
"""
Response to contamination_response_plan.md items 1, 2, and 6 -- all "cheap, on data
you already have," reusing the canonical 240-patient/10,737-cell platelet export
(decontx_io/counts.mtx, genes.csv, cells.csv) rather than re-loading the 5GB raw h5ad.

Item 1: per-cell platelet purity -- fraction of each cell's total counts coming from
a platelet-identity panel (PPBP, PF4, ITGA2B, ITGB3, GP9, GP1BA, TUBB1, NRGN, CLU,
SPARC). Reported as a distribution per cohort/Sample_ID group.

Item 2: cross-lineage co-expression flag -- fraction of each cell's total counts
coming from a non-platelet lineage panel (PTPRC, CD3E, CD3D, LYZ, CD14, FCGR3A,
MS4A1, CD79A, NKG7, HBB, HBA1). Fraction of cells above a stated threshold, per
cohort.

Item 6: depth/mortality/%mito audit for every gene actually named in the R6 section
of main.tex (not just the 13-gene consensus-audit subset) -- same per-patient
Spearman-correlation methodology as consensus_gene_contamination_audit.py, extended
to Geneformer's (RPLP2, RPL24, RPL18A, TIMP1), scGPT's (ALDOA, G6PD, PKM, SHMT2),
UCE's (AXL, C1QA, PLCG2, IL1R2), STATE's (TIMP1, RSU1, TNFSF13B, PFN1), the
consensus-agreement examples (AMFR, ARID5B, ODC1), and the sign-discordant genes
(PPBP, TMSB4X, CCL5). GP1BA included again as the positive control.
"""
import os

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
from scipy.stats import spearmanr

IO_DIR = "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/bib_revision_tier2/contamination_audit/decontx_io"
OUT_DIR = "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/bib_revision_tier2/contamination_audit"
PATIENT_LIST_CSV = ("/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/"
                     "31_donor_level_revision/results_patient_v2/qc_baseline/patient_qc_distributions.csv")

PLATELET_PANEL = ["PPBP", "PF4", "ITGA2B", "ITGB3", "GP9", "GP1BA", "TUBB1", "NRGN", "CLU", "SPARC"]
CROSS_LINEAGE_PANEL = ["PTPRC", "CD3E", "CD3D", "LYZ", "CD14", "FCGR3A", "MS4A1", "CD79A", "NKG7",
                        "HBB", "HBA1"]
CROSS_LINEAGE_THRESHOLD = 0.05  # fraction of total counts

R6_GENES = {
    "RPLP2": "Geneformer_ribosomal", "RPL24": "Geneformer_ribosomal", "RPL18A": "Geneformer_ribosomal",
    "TIMP1": "Geneformer_and_STATE",
    "ALDOA": "scGPT_glycolytic", "G6PD": "scGPT_glycolytic", "PKM": "scGPT_glycolytic",
    "SHMT2": "scGPT_glycolytic",
    "AXL": "UCE_myeloid", "C1QA": "UCE_myeloid", "PLCG2": "UCE_myeloid", "IL1R2": "UCE_myeloid",
    "RSU1": "STATE", "TNFSF13B": "STATE", "PFN1": "STATE",
    "AMFR": "consensus_positive", "ARID5B": "consensus_positive", "ODC1": "near_CCL5_consensus",
    "PPBP": "sign_discordant_platelet", "TMSB4X": "sign_discordant",
    "CCL5": "sign_consistent_negative",
    "GP1BA": "positive_control",
}

print("Loading canonical cohort export...", flush=True)
genes = pd.read_csv(os.path.join(IO_DIR, "genes.csv"))["gene"].tolist()
cells = pd.read_csv(os.path.join(IO_DIR, "cells.csv"))
gene_to_idx = {g: i for i, g in enumerate(genes)}

raw = sio.mmread(os.path.join(IO_DIR, "counts.mtx"))
raw_gc = (raw.tocsc() if raw.shape[0] == len(genes) else raw.T.tocsc())  # genes x cells
print(f"  {raw_gc.shape[0]} genes x {raw_gc.shape[1]} cells, {cells['patient_id'].nunique()} patients",
      flush=True)

total_counts_per_cell = np.asarray(raw_gc.sum(axis=0)).ravel()
total_counts_per_cell[total_counts_per_cell == 0] = 1  # avoid div-by-zero for empty cells

# ---------------------------------------------------------------------
# Item 1: per-cell platelet purity
# ---------------------------------------------------------------------
plt_idx = [gene_to_idx[g] for g in PLATELET_PANEL]
platelet_counts = np.asarray(raw_gc[plt_idx, :].sum(axis=0)).ravel()
platelet_purity = platelet_counts / total_counts_per_cell
cells["platelet_purity"] = platelet_purity

print("\n=== Item 1: per-cell platelet purity (fraction of counts from platelet panel) ===", flush=True)
print(f"Overall: mean={platelet_purity.mean():.3f}, median={np.median(platelet_purity):.3f}, "
      f"p10={np.percentile(platelet_purity,10):.3f}, p90={np.percentile(platelet_purity,90):.3f}",
      flush=True)
purity_by_cohort = cells.groupby("Sample_ID")["platelet_purity"].agg(["mean", "median", "count"])
purity_by_cohort.to_csv(os.path.join(OUT_DIR, "platelet_purity_by_cohort.csv"))
print(f"Per-cohort range: min mean={purity_by_cohort['mean'].min():.3f}, "
      f"max mean={purity_by_cohort['mean'].max():.3f} (n={len(purity_by_cohort)} cohorts)", flush=True)

# ---------------------------------------------------------------------
# Item 2: cross-lineage co-expression flag
# ---------------------------------------------------------------------
xl_idx = [gene_to_idx[g] for g in CROSS_LINEAGE_PANEL]
xl_counts = np.asarray(raw_gc[xl_idx, :].sum(axis=0)).ravel()
xl_fraction = xl_counts / total_counts_per_cell
cells["cross_lineage_fraction"] = xl_fraction
cells["cross_lineage_flag"] = xl_fraction > CROSS_LINEAGE_THRESHOLD

print(f"\n=== Item 2: cross-lineage co-expression (threshold={CROSS_LINEAGE_THRESHOLD}) ===", flush=True)
print(f"Overall: {cells['cross_lineage_flag'].mean()*100:.2f}% of cells flagged "
      f"({cells['cross_lineage_flag'].sum()} / {len(cells)})", flush=True)
flag_by_cohort = cells.groupby("Sample_ID")["cross_lineage_flag"].agg(["mean", "sum", "count"])
flag_by_cohort.to_csv(os.path.join(OUT_DIR, "cross_lineage_flag_by_cohort.csv"))
print(f"Per-cohort flagged-fraction range: min={flag_by_cohort['mean'].min()*100:.2f}%, "
      f"max={flag_by_cohort['mean'].max()*100:.2f}%", flush=True)

cells[["cell_id", "Sample_ID", "patient_id", "y", "platelet_purity", "cross_lineage_fraction",
       "cross_lineage_flag"]].to_csv(os.path.join(OUT_DIR, "per_cell_purity_and_crosslineage.csv"),
                                       index=False)

# ---------------------------------------------------------------------
# Item 6: depth/mortality/%mito audit for every R6-named gene
# ---------------------------------------------------------------------
patient_list = pd.read_csv(PATIENT_LIST_CSV, dtype={"patient_id": str})
patient_list = patient_list[patient_list["disease"] == "COVID-19"]
y_by_patient = patient_list.set_index("patient_id")["y"].to_dict()
qc_by_patient = patient_list.set_index("patient_id")[["median_counts", "median_pct_mito"]]

patient_ids = cells["patient_id"].astype(str).to_numpy()
rows = []
for gene, tag in R6_GENES.items():
    gi = gene_to_idx[gene]
    vals = np.asarray(raw_gc[gi, :].todense()).ravel()
    s = pd.Series(vals, index=patient_ids)
    per_patient = s.groupby(level=0).mean()
    y = pd.Series({p: y_by_patient[p] for p in per_patient.index if p in y_by_patient})
    common = per_patient.index.intersection(y.index)
    rho_y, p_y = spearmanr(per_patient.loc[common], y.loc[common])
    qc = qc_by_patient.reindex(common)
    rho_depth, p_depth = spearmanr(per_patient.loc[common], qc["median_counts"])
    rho_mito, p_mito = spearmanr(per_patient.loc[common], qc["median_pct_mito"])
    rows.append({"gene": gene, "r6_source": tag, "n_patients": len(common),
                 "rho_vs_mortality": rho_y, "p_vs_mortality": p_y,
                 "rho_vs_depth": rho_depth, "p_vs_depth": p_depth,
                 "rho_vs_pct_mito": rho_mito, "p_vs_pct_mito": p_mito})

df = pd.DataFrame(rows).sort_values("r6_source")
df.to_csv(os.path.join(OUT_DIR, "r6_full_gene_contamination_audit.csv"), index=False)
print(f"\n=== Item 6: depth/mortality/%mito audit, all {len(R6_GENES)} R6-named genes "
      f"(n={df['n_patients'].iloc[0]} patients) ===", flush=True)
for _, r in df.iterrows():
    flag = "  <-- depth-driven, null mortality" if (r.p_vs_depth < 0.05 and r.p_vs_mortality > 0.1) else \
           "  <-- direct mortality, no depth confound" if (r.p_vs_mortality < 0.05 and r.p_vs_depth > 0.1) else ""
    print(f"{r.gene:10s} ({r.r6_source:26s}): mortality rho={r.rho_vs_mortality:+.3f} (p={r.p_vs_mortality:.3f})  "
          f"depth rho={r.rho_vs_depth:+.3f} (p={r.p_vs_depth:.3f})  "
          f"pct_mito rho={r.rho_vs_pct_mito:+.3f} (p={r.p_vs_pct_mito:.3f}){flag}", flush=True)

print(f"\nSAVED platelet_purity_by_cohort.csv, cross_lineage_flag_by_cohort.csv, "
      f"per_cell_purity_and_crosslineage.csv, r6_full_gene_contamination_audit.csv", flush=True)
