#!/usr/bin/env python3
"""
Multi-cell-type DecontX rerun, step 4: subset the decontaminated matrix back
to platelet cells only, and redo the flagged-gene per-patient correlation
check (raw vs decontaminated) for the 148-patient subset -- this time with
DecontX run on a genuinely heterogeneous object (real T/B/NK/Mono/DC cells
from Ren/Stephenson/Qiu-sepsis), not a platelet-only matrix. Compares
against the original platelet-only DecontX run's result for the same genes.
"""
import os

import numpy as np
import pandas as pd
import scipy.io as sio
from scipy.stats import spearmanr

IO_DIR = "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/bib_revision_tier2/contamination_audit/decontx_multicelltype"
OUT_DIR = "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/bib_revision_tier2/contamination_audit"
PATIENT_LIST_CSV = ("/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/"
                     "31_donor_level_revision/results_patient_v2/qc_baseline/patient_qc_distributions.csv")

FLAGGED_GENES = {
    "CD79B": "B_cell", "CD79A": "B_cell",
    "PRF1": "NK_cell",
    "S100A8": "Monocyte_Myeloid", "S100A12": "Monocyte_Myeloid",
    "FCER1G": "Monocyte_Myeloid", "TYROBP": "Monocyte_Myeloid",
    "SLC25A37": "Erythroid", "HBA1": "Erythroid",
    "CD3E": "T_cell", "GZMA": "T_cell", "CCL5": "T_cell",
    "GP1BA": "Platelet_reference_positive_control",
}

genes = pd.read_csv(os.path.join(IO_DIR, "multicelltype_genes.csv"))["gene"].tolist()
cells = pd.read_csv(os.path.join(IO_DIR, "multicelltype_cells.csv"))
gene_to_idx = {g: i for i, g in enumerate(genes)}

print("Loading raw and decontaminated multi-cell-type matrices...", flush=True)
raw = sio.mmread(os.path.join(IO_DIR, "multicelltype_counts.mtx"))
decont = sio.mmread(os.path.join(IO_DIR, "multicelltype_decontx_counts.mtx"))
raw_gc = (raw.tocsc() if raw.shape[0] == len(genes) else raw.T.tocsc())
decont_gc = (decont.tocsc() if decont.shape[0] == len(genes) else decont.T.tocsc())
print(f"  raw: {raw_gc.shape}, decont: {decont_gc.shape}, n_cells={len(cells)}", flush=True)

platelet_mask = cells["is_platelet"].astype(bool).to_numpy()
print(f"Platelet cells: {platelet_mask.sum()} / {len(cells)} "
      f"across {cells.loc[platelet_mask, 'patient_id'].nunique()} patients", flush=True)

patient_list = pd.read_csv(PATIENT_LIST_CSV, dtype={"patient_id": str})
patient_list = patient_list[patient_list["disease"] == "COVID-19"]
y_by_patient = patient_list.set_index("patient_id")["y"].to_dict()
qc_by_patient = patient_list.set_index("patient_id")[["median_counts"]]

patient_ids_platelet = cells.loc[platelet_mask, "patient_id"].astype(str).to_numpy()

rows = []
for gene in FLAGGED_GENES:
    if gene not in gene_to_idx:
        print(f"{gene}: not in intersected gene panel, skipped", flush=True)
        continue
    gi = gene_to_idx[gene]
    raw_vals = np.asarray(raw_gc[gi, :].todense()).ravel()[platelet_mask]
    dec_vals = np.asarray(decont_gc[gi, :].todense()).ravel()[platelet_mask]

    for label, vals in [("raw", raw_vals), ("decontaminated_multicelltype", dec_vals)]:
        s = pd.Series(vals, index=patient_ids_platelet)
        per_patient = s.groupby(level=0).mean()
        y = pd.Series({p: y_by_patient[p] for p in per_patient.index if p in y_by_patient})
        common = per_patient.index.intersection(y.index)
        rho_y, p_y = spearmanr(per_patient.loc[common], y.loc[common])
        qc = qc_by_patient.reindex(common)
        rho_depth, p_depth = spearmanr(per_patient.loc[common], qc["median_counts"])
        rows.append({"gene": gene, "lineage": FLAGGED_GENES[gene], "matrix": label,
                      "n_patients": len(common),
                      "rho_vs_mortality": rho_y, "p_vs_mortality": p_y,
                      "rho_vs_depth": rho_depth, "p_vs_depth": p_depth})

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT_DIR, "decontx_multicelltype_flagged_gene_correlations.csv"), index=False)
print("\n=== Multi-cell-type DecontX: flagged-gene correlations, raw vs decontaminated "
      "(148-patient subset, platelet cells only) ===", flush=True)
for gene in FLAGGED_GENES:
    sub = df[df["gene"] == gene]
    if len(sub) == 0:
        continue
    r = sub[sub.matrix == "raw"].iloc[0]
    d = sub[sub.matrix == "decontaminated_multicelltype"].iloc[0]
    print(f"{gene:10s} (n={int(r.n_patients)} pts): mortality rho raw={r.rho_vs_mortality:+.3f} "
          f"(p={r.p_vs_mortality:.3f}) -> decont={d.rho_vs_mortality:+.3f} (p={d.p_vs_mortality:.3f})   |   "
          f"depth rho raw={r.rho_vs_depth:+.3f} -> decont={d.rho_vs_depth:+.3f}", flush=True)

print(f"\nSAVED {os.path.join(OUT_DIR, 'decontx_multicelltype_flagged_gene_correlations.csv')}", flush=True)
