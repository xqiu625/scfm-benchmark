"""Figure 5 — cross-disease transfer. Two forest panels, TWO DIFFERENT SOURCE FILES
per figure_generation_spec.md 修订附录 v2.1 F5:

A: COVID->Sepsis (91 pts/17 fatal, GSE279452 test set) -- from the Chongqing-pool
   file, UNCHANGED by Job J (Jiang was never in this test set).
B: Sepsis->COVID (240 pts/25 fatal) -- from the Job J (Jiang-excluded) file's
   reverse rows, which is the LEGACY-only training arm (12->5 patients after
   excluding Jiang) tested on all 240 COVID patients. This is a different
   training regime from panel A by design (per spec), not an inconsistency.
LR = circle marker, RF = triangle marker, same panel. primary_* columns ONLY --
legacy 4-patient/now-2-patient test-set numbers (TF 0.897/STACK 0.858, and the
Job J file's own COVID->Sepsis rows tested on the tiny 5-patient legacy sepsis
set, which show degenerate 1.000-type AUCs) must not appear anywhere in this
figure -- supplementary-only if used at all.
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from style import MODEL_COLOR, MODEL_ORDER, check_close, save_fig

FWD_DATA = "scripts/31_donor_level_revision/results_patient_v2/cross_disease_chongqing"
FWD_CSV = os.path.join(FWD_DATA, "cross_disease_chongqing_20260904_233237.csv")
REV_DATA = "scripts/31_donor_level_revision/results_patient_v2/cross_disease"
REV_CSV = os.path.join(REV_DATA, "cross_disease_patient_jobJ_20260907_154957.csv")

# spec forbids legacy small-test-set numbers in this figure
FORBIDDEN_VALUES = {0.897, 0.858, 1.000}


def panel(ax, df, direction, title):
    sub = df[df["direction"] == direction]
    order = MODEL_ORDER
    y0 = np.arange(len(order))
    for i, m in enumerate(order):
        color = MODEL_COLOR[m]
        for clf, marker, dy in [("LogReg", "o", 0.12), ("RandomForest", "^", -0.12)]:
            row = sub[(sub["model"] == m) & (sub["classifier"] == clf)]
            if row.empty:
                continue
            row = row.iloc[0]
            y = i + dy
            ax.plot([row["primary_boot_ci_lo"], row["primary_boot_ci_hi"]], [y, y],
                    color=color, lw=1.3, alpha=0.8, zorder=1)
            ax.scatter([row["primary_auc"]], [y], s=70, color=color, marker=marker,
                      zorder=2, edgecolors="black", linewidths=0.5,
                      label=("LR" if clf == "LogReg" else "RF") if i == 0 else None)
    ax.axvline(0.5, color="gray", ls="--", lw=1, zorder=0)
    ax.axvline(0.65, color="gray", ls=":", lw=0.8, zorder=0)
    ax.set_yticks(y0)
    ax.set_yticklabels(order)
    ax.set_xlabel("Patient-level AUC")
    ax.set_title(title, loc="left", fontweight="bold")
    ax.legend(fontsize=7, frameon=False, loc="upper left")
    n_test_pt = sub["primary_n_test_patients"].iloc[0]
    n_fatal_pt = sub["primary_n_fatal_test_patients"].iloc[0]
    ax.text(0.98, 0.02, f"n={n_test_pt} test patients ({n_fatal_pt} fatal)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7, color="dimgray")


def main():
    fwd_df = pd.read_csv(FWD_CSV)
    rev_df = pd.read_csv(REV_CSV)

    # hard guard: forbidden legacy/degenerate test-set numbers must not appear
    # anywhere in what's actually plotted (fwd panel uses fwd_df only, rev panel
    # uses rev_df's Sepsis->COVID rows only -- Job J's own COVID->Sepsis rows,
    # tested on the tiny 5-patient legacy set, are never loaded for plotting).
    rev_plot_rows = rev_df[rev_df["direction"] == "Sepsis->COVID"]
    plotted_vals = (set(fwd_df["primary_auc"].round(3)) |
                     set(fwd_df["primary_boot_auc_mean"].round(3)) |
                     set(rev_plot_rows["primary_auc"].round(3)) |
                     set(rev_plot_rows["primary_boot_auc_mean"].round(3)))
    for bad in FORBIDDEN_VALUES:
        if any(abs(v - bad) < 0.002 for v in plotted_vals):
            raise SystemExit(f"FORBIDDEN legacy/degenerate value {bad} present in "
                              f"primary_* columns -- do not plot Figure 5, investigate.")

    # spec anchor cross-check: forward 0.644-0.684 all overlap (unchanged, Chongqing file)
    fwd_lr = fwd_df[(fwd_df["direction"] == "COVID->Sepsis") & (fwd_df["classifier"] == "LogReg")]
    check_close("Figure5", "forward AUC range low", 0.644, fwd_lr["primary_auc"].min(), FWD_CSV)
    check_close("Figure5", "forward AUC range high", 0.684, fwd_lr["primary_auc"].max(), FWD_CSV)

    # reverse anchors: Job J (Jiang-excluded) file -- new v2.1 anchors
    rev_lr = rev_plot_rows[rev_plot_rows["classifier"] == "LogReg"].set_index("model")
    rev_anchor_lr = {"STATE": 0.799, "scPRINT-2": 0.695, "Geneformer": 0.668,
                      "scGPT": 0.661, "UCE": 0.653}
    for m, v in rev_anchor_lr.items():
        check_close("Figure5", f"reverse {m}-LR AUC (post-Jiang)", v,
                    rev_lr.loc[m, "primary_auc"], REV_CSV)
    state_rf_rev = rev_plot_rows[(rev_plot_rows["model"] == "STATE") &
                                  (rev_plot_rows["classifier"] == "RandomForest")].iloc[0]
    check_close("Figure5", "STATE-RF reverse AUC (post-Jiang)", 0.653,
                state_rf_rev["primary_auc"], REV_CSV)
    check_close("Figure5", "STATE-LR reverse CI lo (post-Jiang)", 0.718,
                rev_lr.loc["STATE", "primary_boot_ci_lo"], REV_CSV)
    check_close("Figure5", "STATE-LR reverse CI hi (post-Jiang)", 0.876,
                rev_lr.loc["STATE", "primary_boot_ci_hi"], REV_CSV)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    panel(axes[0], fwd_df, "COVID->Sepsis", "A  COVID -> Sepsis")
    panel(axes[1], rev_plot_rows, "Sepsis->COVID", "B  Sepsis -> COVID (post-Jiang exclusion)")
    fig.tight_layout()
    save_fig(fig, "Figure5_cross_disease")


if __name__ == "__main__":
    main()
