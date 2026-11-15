#!/usr/bin/env python3
"""
Follow-up to decontx_reaudit.py, addressing external review points A (genome-
wide |delta rho| reference distribution) and H (depth-outcome relationship),
and re-verifying D (level of analysis) explicitly with printed n.

Review: ~/Downloads/decontx_audit_review.md (2026-09-13). Point A: a Delta
rho <= 0.03 for the 13 flagged genes is only evidence of robustness if it is
small relative to how much OTHER genes moved under the same correction --
report the full 13,066-gene |Delta rho| distribution (vs depth AND vs
mortality) and where the 13 flagged genes fall in it. Point H: "tracks
depth, not mortality" needs the depth-outcome relationship itself stated,
since depth is not outcome-neutral in this cohort.
"""
import os

import numpy as np
import pandas as pd
import scipy.io as sio
from scipy.stats import spearmanr

IO_DIR = "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/bib_revision_tier2/contamination_audit/decontx_io"
OUT_DIR = "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/bib_revision_tier2/contamination_audit"
PATIENT_LIST_CSV = ("/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/"
                     "31_donor_level_revision/results_patient_v2/qc_baseline/patient_qc_distributions.csv")

FLAGGED_GENES = ["CD79B", "CD79A", "PRF1", "S100A8", "S100A12", "FCER1G", "TYROBP",
                  "SLC25A37", "HBA1", "CD3E", "GZMA", "CCL5", "GP1BA"]

genes = pd.read_csv(os.path.join(IO_DIR, "genes.csv"))["gene"].tolist()
cells = pd.read_csv(os.path.join(IO_DIR, "cells.csv"))
gene_to_idx = {g: i for i, g in enumerate(genes)}

print("Loading raw and decontaminated matrices...", flush=True)
raw = sio.mmread(os.path.join(IO_DIR, "counts.mtx"))
decont = sio.mmread(os.path.join(IO_DIR, "decontx_counts.mtx"))
raw_gc = (raw.tocsc() if raw.shape[0] == len(genes) else raw.T.tocsc())
decont_gc = (decont.tocsc() if decont.shape[0] == len(genes) else decont.T.tocsc())
print(f"  raw: {raw_gc.shape}, decont: {decont_gc.shape}", flush=True)

patient_list = pd.read_csv(PATIENT_LIST_CSV, dtype={"patient_id": str})
patient_list = patient_list[patient_list["disease"] == "COVID-19"]
qc_by_patient = patient_list.set_index("patient_id")[["median_counts"]]
y_by_patient = patient_list.set_index("patient_id")["y"].to_dict()

patient_ids = cells["patient_id"].to_numpy()
unique_patients = sorted(set(patient_ids))
y_s = pd.Series({p: y_by_patient[p] for p in unique_patients if p in y_by_patient})
depth_s = qc_by_patient.reindex(unique_patients)["median_counts"]
n_patients_common_y = len(y_s)
print(f"n_cells={len(patient_ids)}, n_unique_patients={len(unique_patients)}, "
      f"n_patients_with_mortality_label={n_patients_common_y} -- "
      f"ALL Q1-style correlations below are computed at this patient level, NOT per-cell.",
      flush=True)

# Filter out all-zero genes to avoid degenerate correlations
n_genes = len(genes)
rows = []
for gi in range(n_genes):
    raw_col = np.asarray(raw_gc[gi, :].todense()).ravel()
    dec_col = np.asarray(decont_gc[gi, :].todense()).ravel()
    if raw_col.sum() == 0 and dec_col.sum() == 0:
        continue
    raw_pp = pd.Series(raw_col, index=patient_ids).groupby(level=0).mean()
    dec_pp = pd.Series(dec_col, index=patient_ids).groupby(level=0).mean()

    common_y = raw_pp.index.intersection(y_s.index)
    if raw_pp.loc[common_y].std() == 0 or dec_pp.loc[common_y].std() == 0:
        rho_y_raw, rho_y_dec = np.nan, np.nan
    else:
        rho_y_raw, _ = spearmanr(raw_pp.loc[common_y], y_s.loc[common_y])
        rho_y_dec, _ = spearmanr(dec_pp.loc[common_y], y_s.loc[common_y])

    common_d = raw_pp.index.intersection(depth_s.dropna().index)
    if raw_pp.loc[common_d].std() == 0 or dec_pp.loc[common_d].std() == 0:
        rho_d_raw, rho_d_dec = np.nan, np.nan
    else:
        rho_d_raw, _ = spearmanr(raw_pp.loc[common_d], depth_s.loc[common_d])
        rho_d_dec, _ = spearmanr(dec_pp.loc[common_d], depth_s.loc[common_d])

    rows.append({"gene": genes[gi], "rho_mortality_raw": rho_y_raw, "rho_mortality_decont": rho_y_dec,
                 "rho_depth_raw": rho_d_raw, "rho_depth_decont": rho_d_dec})
    if gi % 2000 == 0:
        print(f"  ... {gi}/{n_genes} genes processed", flush=True)

df = pd.DataFrame(rows)
df["abs_delta_rho_mortality"] = (df["rho_mortality_decont"] - df["rho_mortality_raw"]).abs()
df["abs_delta_rho_depth"] = (df["rho_depth_decont"] - df["rho_depth_raw"]).abs()
df = df.dropna(subset=["abs_delta_rho_mortality", "abs_delta_rho_depth"])
df.to_csv(os.path.join(OUT_DIR, "decontx_genomewide_delta_rho.csv"), index=False)

print(f"\n=== Genome-wide |Delta rho| reference distribution (n={len(df)} genes with nonzero counts) ===",
      flush=True)
for col in ["abs_delta_rho_mortality", "abs_delta_rho_depth"]:
    q = df[col].quantile([0.5, 0.75, 0.90, 0.95, 0.99])
    print(f"{col}: median={q[0.5]:.4f}  p75={q[0.75]:.4f}  p90={q[0.90]:.4f}  "
          f"p95={q[0.95]:.4f}  p99={q[0.99]:.4f}  max={df[col].max():.4f}", flush=True)

print("\n=== Where the 13 flagged genes fall in the genome-wide |Delta rho| distribution ===", flush=True)
for g in FLAGGED_GENES:
    row = df[df["gene"] == g]
    if len(row) == 0:
        print(f"{g:10s}: not found / all-zero, skipped")
        continue
    row = row.iloc[0]
    pct_mort = (df["abs_delta_rho_mortality"] < row["abs_delta_rho_mortality"]).mean() * 100
    pct_depth = (df["abs_delta_rho_depth"] < row["abs_delta_rho_depth"]).mean() * 100
    print(f"{g:10s}: |d_rho_mortality|={row['abs_delta_rho_mortality']:.4f} (genome pctile={pct_mort:5.1f}%)   "
          f"|d_rho_depth|={row['abs_delta_rho_depth']:.4f} (genome pctile={pct_depth:5.1f}%)", flush=True)

# Point H: is depth itself outcome-associated in this cohort?
rho_depth_outcome, p_depth_outcome = spearmanr(depth_s.loc[y_s.index], y_s)
print(f"\n=== Point H: is median_counts (depth) itself associated with mortality in this "
      f"{len(y_s)}-patient cohort? ===", flush=True)
print(f"Spearman rho(median_counts, mortality) = {rho_depth_outcome:+.3f} (p={p_depth_outcome:.3f})", flush=True)
print("SAVED decontx_genomewide_delta_rho.csv", flush=True)
