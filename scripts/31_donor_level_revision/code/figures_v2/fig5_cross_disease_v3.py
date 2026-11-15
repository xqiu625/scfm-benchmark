"""Figure 5 — cross-disease transfer, v3 spec, 2 panels.

A: bidirectional dumbbell. Forward (COVID->Sepsis, 91-patient GSE279452 test
   pool, UNCHANGED source -- cross_disease_chongqing file, Jiang was never in
   this test set) vs reverse (Sepsis->COVID, 96-patient training pool
   [91 GSE279452 + 5 GSE167363], 240-patient test -- UPDATED source, the new
   96pool file) per model, LR only (matches the spec's explicit numbers).
B: adaptation curve (source unchanged), with a new gray unresolvable band
   [0.15, 0.85] at the 32-patient-test-set endpoint.
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from style import MODEL_COLOR, MODEL_ORDER, check_close, save_fig

FWD_CSV = ("scripts/31_donor_level_revision/results_patient_v2/cross_disease_chongqing/"
           "cross_disease_chongqing_20260904_233237.csv")
REV_CSV = ("scripts/31_donor_level_revision/results_patient_v2/cross_disease/"
           "cross_disease_patient_96pool_20260908_142954.csv")
ADAPT_DIR = "scripts/31_donor_level_revision/results_patient_v2/adaptation_curve_chongqing"

FWD_ANCHOR = {"UCE": 0.684, "STATE": 0.653, "scGPT": 0.650, "scPRINT-2": 0.647, "Geneformer": 0.644}
REV_ANCHOR = {"STATE": 0.7209, "Geneformer": 0.7118, "UCE": 0.6928, "scGPT": 0.6417, "scPRINT-2": 0.6186}


def panel_a(ax, fwd_df, rev_df):
    fwd_lr = fwd_df[(fwd_df["direction"] == "COVID->Sepsis") &
                     (fwd_df["classifier"] == "LogReg")].set_index("model")
    rev_lr = rev_df[(rev_df["direction"] == "Sepsis->COVID") &
                     (rev_df["classifier"] == "LogReg")].set_index("model")
    order = sorted(MODEL_ORDER, key=lambda m: rev_lr.loc[m, "primary_auc"])
    y = np.arange(len(order))
    for i, m in enumerate(order):
        color = MODEL_COLOR[m]
        f = fwd_lr.loc[m, "primary_auc"]
        r = rev_lr.loc[m, "primary_auc"]
        ax.plot([f, r], [i, i], color=color, lw=1.5, alpha=0.6, zorder=1)
        ax.scatter([f], [i], s=70, color="white", edgecolors=color, linewidths=2, zorder=2,
                  marker="o")
        ax.scatter([r], [i], s=70, color=color, zorder=3, marker="o")
    ax.set_yticks(y)
    ax.set_yticklabels(order)
    ax.axvline(0.5, color="gray", ls=":", lw=0.8, zorder=0)
    ax.set_xlabel("Patient-level AUC (LR)")
    ax.set_title("A", loc="left", fontweight="bold")
    from matplotlib.lines import Line2D
    legend_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white",
               markeredgecolor="dimgray", markeredgewidth=1.6, markersize=8,
               label="forward (COVID→Sepsis)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="dimgray",
               markeredgecolor="dimgray", markersize=8, label="reverse (Sepsis→COVID)"),
    ]
    ax.legend(handles=legend_handles, fontsize=7.5, frameon=False, loc="upper left")
    ax.text(0.98, 0.03,
            "forward: uninformative (mortality not learnable in the\n"
            "sepsis pool even in-domain, 91 pts/17 fatal)\n"
            "reverse: rank-indeterminate (top-2 reversed under\n"
            "severe-excluded; UCE shifts +0.086 with 5 added\n"
            "training patients; 240 pts/25 fatal)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6.4, color="dimgray")


def panel_b(ax, exp):
    for m in MODEL_ORDER + ["PCA-200+LR"]:
        sub = exp[exp["method"] == m]
        agg = sub.groupby("n_sepsis_req")["auc_point"].agg(["mean", "std"]).reset_index()
        agg = agg.sort_values("n_sepsis_req")
        x = agg["n_sepsis_req"].replace(-1, sub["n_sepsis"].max())
        x_plot = x.replace(0, 1)
        color = MODEL_COLOR.get(m, MODEL_COLOR.get(m.replace("+LR", ""), "black"))
        ax.plot(x_plot, agg["mean"], color=color, lw=1.8, marker="o", markersize=4, label=m)
        ax.fill_between(x_plot, agg["mean"] - agg["std"], agg["mean"] + agg["std"],
                        color=color, alpha=0.12)
    ax.axhspan(0.15, 0.85, color="gray", alpha=0.12, zorder=0)
    ax.text(0.98, 0.93, "unresolvable at this test-set size\n"
            "(≈32-patient test set, ~3 fatal events,\n"
            "95% interval spans ≈0.70)",
            transform=ax.transAxes, ha="right", va="top", fontsize=6.8, color="dimgray")
    ax.set_xscale("log")
    ax.axhline(0.5, color="gray", ls=":", lw=0.8, zorder=0)
    ax.set_xlabel("Sepsis training cells added (log scale; leftmost point = 0)")
    ax.set_ylabel("AUC (32-patient test set, 96-patient sepsis pool)")
    ax.set_title("B", loc="left", fontweight="bold")
    ax.legend(fontsize=6.3, frameon=False, ncol=2, loc="lower right")


def main():
    fwd_df = pd.read_csv(FWD_CSV)
    rev_df = pd.read_csv(REV_CSV)
    adapt_all = pd.read_csv(os.path.join(ADAPT_DIR, "all_rows.csv"))
    exp = adapt_all[adapt_all["test_set"] == "expanded_96pt"].copy()

    fwd_lr = fwd_df[(fwd_df["direction"] == "COVID->Sepsis") & (fwd_df["classifier"] == "LogReg")].set_index("model")
    rev_lr = rev_df[(rev_df["direction"] == "Sepsis->COVID") & (rev_df["classifier"] == "LogReg")].set_index("model")
    for m, v in FWD_ANCHOR.items():
        check_close("Figure5v3", f"forward {m}", v, fwd_lr.loc[m, "primary_auc"], FWD_CSV, tol=0.001)
    for m, v in REV_ANCHOR.items():
        check_close("Figure5v3", f"reverse {m}", v, rev_lr.loc[m, "primary_auc"], REV_CSV, tol=0.001)

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    panel_a(axes[0], fwd_df, rev_df)
    panel_b(axes[1], exp)
    fig.tight_layout()
    save_fig(fig, "Figure5_cross_disease")


if __name__ == "__main__":
    main()
