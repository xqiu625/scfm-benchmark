"""Figure S-LOBO (new, v3 spec) -- leave-one-cohort-out generalization.

Left: one row per FM, 5 cohort-level held-out AUCs (point size prop. to
n_fatal_held_out) plus a vertical line at the pooled (concatenated held-out
predictions, 206 patients/25 fatal) AUC.
Right: small per-cohort n/fatal table.
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from style import MODEL_COLOR, MODEL_ORDER, check_close, save_fig

DATA_LOBO = "scripts/31_donor_level_revision/results_patient_v2/lobo_v2"
DATA_PI = "scripts/31_donor_level_revision/results_patient_v2/paired_inference_full"

POOLED_ANCHOR = {"STATE": 0.7525, "scPRINT-2": 0.7514, "Geneformer": 0.7481,
                  "UCE": 0.6813, "scGPT": 0.6376}


def main():
    per_cohort = pd.read_csv(os.path.join(DATA_LOBO, "lobo_per_cohort_auc.csv"))
    pooled = pd.read_csv(os.path.join(DATA_PI, "lobo_pooled.csv")).set_index("fm")
    cohort_n = pd.read_csv(os.path.join(DATA_PI, "lobo_pooled_per_cohort_n.csv"))

    for m, v in POOLED_ANCHOR.items():
        check_close("FigureS-LOBO", f"{m} pooled AUC", v, pooled.loc[m, "pooled_auc"],
                    os.path.join(DATA_PI, "lobo_pooled.csv"))

    order = sorted(MODEL_ORDER, key=lambda m: pooled.loc[m, "pooled_auc"])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"width_ratios": [2.2, 1]})

    y = np.arange(len(order))
    for i, m in enumerate(order):
        color = MODEL_COLOR[m]
        sub = per_cohort[per_cohort["fm"] == m]
        sizes = 18 + sub["n_fatal_held_out"].to_numpy() * 9
        jitter = np.random.default_rng(0).normal(0, 0.08, size=len(sub))
        ax1.scatter(sub["held_out_auc"], i + jitter, s=sizes, color=color, alpha=0.55,
                   edgecolors="none", zorder=2)
        p = pooled.loc[m, "pooled_auc"]
        ax1.plot([p, p], [i - 0.32, i + 0.32], color=color, lw=2.4, zorder=3)
    ax1.axvline(0.5, color="gray", ls=":", lw=0.8, zorder=0)
    ax1.set_yticks(y)
    ax1.set_yticklabels(order)
    ax1.set_xlabel("Held-out cohort AUC (dots, size scales with n fatal) / pooled AUC (vertical bar)")
    ax1.set_title("A", loc="left", fontweight="bold")

    cohorts = cohort_n[cohort_n["fm"] == "STATE"][["cohort_data_no", "n_patients_held_out",
                                                     "n_fatal_held_out"]].reset_index(drop=True)
    ax2.axis("off")
    tbl = ax2.table(
        cellText=cohorts.values.tolist(),
        colLabels=["cohort", "n patients", "n fatal"],
        loc="center", cellLoc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8)
    tbl.scale(1, 1.6)
    ax2.set_title("B  per-cohort n (held out)", loc="left", fontweight="bold", fontsize=10)

    fig.tight_layout()
    fig.text(0.5, -0.04,
              "LOBO is outcome-selected: folds require ≥1 fatal and ≥1 survivor in the held-out "
              "cohort; survivor-only cohorts (34 patients) train only. Pooled AUC over concatenated "
              "held-out predictions (206 patients, 25 fatal).",
              ha="center", fontsize=7.3, color="dimgray")
    save_fig(fig, "FigureS_LOBO")


if __name__ == "__main__":
    main()
