#!/usr/bin/env python3
"""Supplementary figure: embedding-dimension dropout robustness.

Patient-level AUC as a function of randomly zeroed embedding dimensions
(0-70%), per model, from the post-audit mortality-endpoint rerun
(scripts/07_gene_dropout/results/gene_dropout_robustness_20260806_162334.csv;
legacy analysis pool, 5 patient-grouped folds, mean +/- 1 SD).

UCE's flat curve (delta +0.003 at 70% zeroing) is the direct visualization of
the main text's "UCE is the most dropout-resistant" statement
(Section r5_calibration). Stdlib + numpy + matplotlib only.
"""
import csv
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import MODEL_COLOR, save_fig

CSV = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir, os.pardir,
    "07_gene_dropout", "results",
    "gene_dropout_robustness_20260806_162334.csv"))

DISPLAY = {"STATE": "STATE", "UCE": "UCE", "Geneformer": "Geneformer",
           "scPRINT2": "scPRINT-2", "scGPT": "scGPT"}
ORDER = ["STATE", "UCE", "Geneformer", "scPRINT2", "scGPT"]


def main():
    rows = list(csv.DictReader(open(CSV)))
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    for m in ORDER:
        pts = sorted((int(r["Dropout_Pct"]), float(r["AUC_Mean"]), float(r["AUC_Std"]))
                     for r in rows if r["Model"] == m)
        x = [p[0] for p in pts]
        y = np.array([p[1] for p in pts])
        sd = np.array([p[2] for p in pts])
        name = DISPLAY[m]
        color = MODEL_COLOR[name]
        ax.plot(x, y, "-o", ms=4, lw=1.4, color=color, label=name)
        ax.fill_between(x, y - sd, y + sd, color=color, alpha=0.12)
    ax.set_xlabel("Embedding dimensions randomly zeroed (%)", fontsize=9)
    ax.set_ylabel("Patient-level AUC", fontsize=9)
    ax.set_title("Embedding-dimension dropout robustness\n"
                 "(mortality endpoint, legacy analysis pool; mean ± 1 SD over 5 patient-grouped folds)",
                 fontsize=9, loc="left")
    ax.legend(fontsize=8, frameon=False)
    ax.tick_params(labelsize=8)
    fig.tight_layout()
    save_fig(fig, "FigureS_dropout")


if __name__ == "__main__":
    main()
