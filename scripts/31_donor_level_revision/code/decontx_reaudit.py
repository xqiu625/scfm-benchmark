#!/usr/bin/env python3
"""
Consume run_decontx.R's output and answer the two questions this check was
built for:
  1. Do the flagged genes' contamination signature (strong depth
     correlation, weak/null mortality correlation) attenuate once DecontX-
     corrected expression is used instead of raw counts?
  2. Does contamination correction change which genes land in the
     masking-attribution pipeline's top-256-by-expression candidate set
     per cell (the actual gene universe the R6 consensus attribution was
     computed on) -- i.e., were some of the flagged non-platelet genes
     candidates in cells' masking sweeps ONLY because contamination
     inflated their raw expression above the top-256 cutoff?

Uses the exact same top_k_genes_per_cell selection logic (top-256 by
expression, stable sort) as masking_attribution_common.py, applied to both
the raw and DecontX-corrected matrices for a direct before/after
comparison.
"""
import os
import sys

import numpy as np
import pandas as pd
import scipy.io as sio
from scipy.stats import spearmanr

IO_DIR = "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/bib_revision_tier2/contamination_audit/decontx_io"
OUT_DIR = "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/bib_revision_tier2/contamination_audit"
PATIENT_LIST_CSV = ("/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/"
                     "31_donor_level_revision/results_patient_v2/qc_baseline/patient_qc_distributions.csv")

TOP_K = 256

FLAGGED_GENES = {
    "CD79B": "B_cell", "CD79A": "B_cell",
    "PRF1": "NK_cell",
    "S100A8": "Monocyte_Myeloid", "S100A12": "Monocyte_Myeloid",
    "FCER1G": "Monocyte_Myeloid", "TYROBP": "Monocyte_Myeloid",
    "SLC25A37": "Erythroid", "HBA1": "Erythroid",
    "CD3E": "T_cell", "GZMA": "T_cell", "CCL5": "T_cell",
    "GP1BA": "Platelet_reference_positive_control",
}


def top_k_genes_per_cell(X_csc_genes_x_cells, k=TOP_K):
    """X: genes x cells (as exported). Returns list of length n_cells, each
    a set of gene row-indices in that cell's top-k by expression (matches
    masking_attribution_common.py's stable-sort, descending-expression
    selection, restricted to nonzero entries)."""
    Xc = X_csc_genes_x_cells.tocsc()
    n_cells = Xc.shape[1]
    out = []
    for i in range(n_cells):
        col = Xc[:, i]
        row_idx = col.indices
        vals = col.data
        if len(row_idx) == 0:
            out.append(set())
            continue
        order = row_idx[np.argsort(-vals, kind="stable")][:k]
        out.append(set(order.tolist()))
    return out


genes = pd.read_csv(os.path.join(IO_DIR, "genes.csv"))["gene"].tolist()
cells = pd.read_csv(os.path.join(IO_DIR, "cells.csv"))
gene_to_idx = {g: i for i, g in enumerate(genes)}

print("Loading raw and decontaminated matrices...", flush=True)
raw = sio.mmread(os.path.join(IO_DIR, "counts.mtx")).tocsc()   # cells x genes as written by mmwrite(X.T...)
decont = sio.mmread(os.path.join(IO_DIR, "decontx_counts.mtx"))

# export_for_decontx.py wrote sio.mmwrite(..., X.T.tocsr()) where X was
# cells x genes, so counts.mtx is genes x cells (rows=genes). R's writeMM
# on decont_counts (also genes x cells, celda convention) matches this
# orientation directly -- both should already be genes x cells.
if raw.shape[0] == len(genes):
    raw_genes_x_cells = raw.tocsc()
else:
    raw_genes_x_cells = raw.T.tocsc()
if decont.shape[0] == len(genes):
    decont_genes_x_cells = decont.tocsc()
else:
    decont_genes_x_cells = decont.T.tocsc()
print(f"  raw: {raw_genes_x_cells.shape}, decont: {decont_genes_x_cells.shape}, "
      f"n_genes={len(genes)}, n_cells={len(cells)}", flush=True)
assert raw_genes_x_cells.shape == decont_genes_x_cells.shape, "shape mismatch between raw and decontaminated matrices"

contam = pd.read_csv(os.path.join(IO_DIR, "decontx_contamination_scores.csv"))
cells = cells.merge(contam, on="cell_id", how="left")
print(f"Contamination fraction: mean={cells['decontX_contamination'].mean():.4f}, "
      f"median={cells['decontX_contamination'].median():.4f}", flush=True)

# ---------------------------------------------------------------------
# Q1: flagged-gene correlations, raw vs decontaminated per-patient mean
# ---------------------------------------------------------------------
patient_list = pd.read_csv(PATIENT_LIST_CSV, dtype={"patient_id": str})
patient_list = patient_list[patient_list["disease"] == "COVID-19"]
qc_by_patient = patient_list.set_index("patient_id")[["median_counts", "median_genes", "median_pct_mito"]]
y_by_patient = patient_list.set_index("patient_id")["y"].to_dict()

rows = []
for gene in FLAGGED_GENES:
    if gene not in gene_to_idx:
        continue
    gi = gene_to_idx[gene]
    raw_vals = np.asarray(raw_genes_x_cells[gi, :].todense()).ravel()
    decont_vals = np.asarray(decont_genes_x_cells[gi, :].todense()).ravel()

    for label, vals in [("raw", raw_vals), ("decontaminated", decont_vals)]:
        s = pd.Series(vals, index=cells["patient_id"].to_numpy())
        per_patient = s.groupby(level=0).mean()
        unique_patients = sorted(per_patient.index)
        y = pd.Series({p: y_by_patient[p] for p in unique_patients if p in y_by_patient})
        common = per_patient.index.intersection(y.index)
        rho_y, p_y = spearmanr(per_patient.loc[common], y.loc[common])
        qc = qc_by_patient.reindex(common)
        rho_depth, p_depth = spearmanr(per_patient.loc[common], qc["median_counts"])
        rows.append({"gene": gene, "lineage": FLAGGED_GENES[gene], "matrix": label,
                      "rho_vs_mortality": rho_y, "p_vs_mortality": p_y,
                      "rho_vs_depth": rho_depth, "p_vs_depth": p_depth})

df1 = pd.DataFrame(rows)
df1.to_csv(os.path.join(OUT_DIR, "decontx_flagged_gene_correlations.csv"), index=False)
print("\n=== Q1: flagged-gene correlations, raw vs DecontX-corrected ===", flush=True)
for gene in FLAGGED_GENES:
    sub = df1[df1["gene"] == gene]
    if len(sub) == 0:
        continue
    r = sub[sub.matrix == "raw"].iloc[0]
    d = sub[sub.matrix == "decontaminated"].iloc[0]
    print(f"{gene:10s}: mortality rho raw={r.rho_vs_mortality:+.3f} (p={r.p_vs_mortality:.3f}) "
          f"-> decont={d.rho_vs_mortality:+.3f} (p={d.p_vs_mortality:.3f})   |   "
          f"depth rho raw={r.rho_vs_depth:+.3f} -> decont={d.rho_vs_depth:+.3f}", flush=True)

# ---------------------------------------------------------------------
# Q2: does the top-256 candidate gene set change after decontamination?
# ---------------------------------------------------------------------
print("\n=== Q2: top-256 candidate-gene-set membership, raw vs decontaminated ===", flush=True)
raw_topk = top_k_genes_per_cell(raw_genes_x_cells)
decont_topk = top_k_genes_per_cell(decont_genes_x_cells)

flagged_idx = {g: gene_to_idx[g] for g in FLAGGED_GENES if g in gene_to_idx}
rows2 = []
for gene, gi in flagged_idx.items():
    n_in_raw = sum(1 for s in raw_topk if gi in s)
    n_in_decont = sum(1 for s in decont_topk if gi in s)
    n_dropped = sum(1 for r, d in zip(raw_topk, decont_topk) if gi in r and gi not in d)
    n_gained = sum(1 for r, d in zip(raw_topk, decont_topk) if gi not in r and gi in d)
    rows2.append({"gene": gene, "lineage": FLAGGED_GENES[gene],
                  "n_cells_candidate_raw": n_in_raw, "n_cells_candidate_decont": n_in_decont,
                  "n_dropped_after_decont": n_dropped, "n_gained_after_decont": n_gained,
                  "n_cells_total": len(raw_topk)})
    print(f"{gene:10s}: candidate in {n_in_raw}/{len(raw_topk)} cells (raw) -> "
          f"{n_in_decont}/{len(raw_topk)} (decont); dropped={n_dropped}, gained={n_gained}", flush=True)

df2 = pd.DataFrame(rows2)
df2.to_csv(os.path.join(OUT_DIR, "decontx_candidate_gene_impact.csv"), index=False)
print(f"\nSAVED {os.path.join(OUT_DIR, 'decontx_flagged_gene_correlations.csv')} and "
      f"decontx_candidate_gene_impact.csv", flush=True)
