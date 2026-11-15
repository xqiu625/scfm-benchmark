"""Figure 4 — in-distribution (monocyte) positive control, v3 spec.

Per figure_table_spec_v3.md Sec.1 "Figure 4": relabel every "OOD" occurrence
to "underrepresented (platelets)"; "in-distribution (monocytes)" is kept.
Headline panel A uses the 148-patient/19-fatal cut (ALL_incl_Qiu_sepsis,
matches P0 Sec.3c main table); panel B's significance test is the
143-patient/17-fatal COVID-only subset (pre-existing, already reconciled --
the Qiu-sepsis n=5 addition is noise-only, excluded from the formal test,
see panel C). Both n's stated explicitly in-figure so the difference is
never silently blended.
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from style import MODEL_COLOR, MODEL_ORDER, check_close, save_fig

DATA = "scripts/31_donor_level_revision/results_patient_v2/monocyte"
PI = "scripts/31_donor_level_revision/results_patient_v2/paired_inference_full"

SPEC_ANCHOR = {"STATE": 0.696, "scGPT": 0.677, "Geneformer": 0.674,
               "scPRINT-2": 0.672, "UCE": 0.648, "PCA-200": 0.654}


def panel_a(ax, df):
    cut = df[df["cut"] == "ALL_incl_Qiu_sepsis"].set_index("fm")
    order = cut.sort_values("auc_platelet_mean", ascending=True).index.tolist()
    for i, m in enumerate(order):
        row = cut.loc[m]
        color = MODEL_COLOR[m]
        ax.plot([row["auc_monocyte_mean"], row["auc_platelet_mean"]], [i, i],
                color=color, lw=1.5, zorder=1, alpha=0.7)
        ax.scatter([row["auc_monocyte_mean"]], [i], s=70, color="white",
                  edgecolors=color, linewidths=2, zorder=2,
                  label="in-distribution (monocytes)" if i == 0 else None)
        ax.scatter([row["auc_platelet_mean"]], [i], s=70, color=color, zorder=3,
                  label="underrepresented (platelets)" if i == 0 else None)
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels(order)
    ax.set_xlabel("Patient-level repeat-mean AUC")
    ax.set_title("A", loc="left", fontweight="bold")
    ax.axvline(0.5, color="gray", ls=":", lw=0.8, zorder=0)
    ax.legend(loc="lower right", frameon=False, fontsize=7.5)
    ax.text(0.02, 0.02, "in-distribution (monocytes): 148 pts/19 fatal\nunderrepresented (platelets): 240 pts/25 fatal",
            transform=ax.transAxes, fontsize=7, color="dimgray", va="bottom")


def panel_b(ax, pw):
    order = pw.sort_values("delta_auc").reset_index(drop=True)
    y = np.arange(len(order))
    for i, row in order.iterrows():
        m = row["pair"].split(" - ")[0]
        color = MODEL_COLOR.get(m, "black")
        ax.plot([row["ci_lo"], row["ci_hi"]], [i, i], color=color, lw=1.4,
               alpha=0.5, zorder=1)
        ax.scatter([row["delta_auc"]], [i], s=80, color=color,
                  marker="D", zorder=2, alpha=0.9, edgecolors="none")
        ax.text(row["ci_hi"] + 0.01, i, f"p={row['p_boot']:.2f}", va="center",
               fontsize=6.5, color="dimgray")
    ax.axvline(0, color="gray", ls="--", lw=1, zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels([r["pair"].split(" - ")[0] for _, r in order.iterrows()])
    ax.set_xlabel("$\\Delta$ AUC (FM − PCA-200), in-distribution (monocytes, 143 pts/17 fatal)")
    ax.set_title("B  (no pair significant, seed-ensemble bootstrap p > 0.23)", loc="left",
                fontweight="bold", fontsize=9)


def panel_c(ax, strat):
    cuts = ["Ren_only", "Stephenson_only", "Qiu_sepsis_only"]
    labels = ["Ren", "Stephenson", "Qiu-sepsis"]
    for m in MODEL_ORDER + ["PCA-200"]:
        vals = [strat[(strat["fm"] == m) & (strat["cut"] == c)]["auc_repeat_mean"].iloc[0]
                for c in cuts]
        ax.plot(range(len(cuts)), vals, marker="o", color=MODEL_COLOR[m], label=m, lw=1.3,
               markersize=6)
    ax.set_xticks(range(len(cuts)))
    ax.set_xticklabels(labels)
    ax.axhline(0.5, color="gray", ls=":", lw=0.8)
    ax.text(2, 0.28, "n=5, noise", fontsize=7, color="dimgray", ha="center")
    ax.set_ylabel("Repeat-mean AUC")
    ax.set_title("C", loc="left", fontweight="bold")
    ax.legend(fontsize=6.5, frameon=False, ncol=2, loc="upper left")


def main():
    df = pd.read_csv(os.path.join(DATA, "indistribution_vs_platelet.csv"))
    pw = pd.read_csv(os.path.join(PI, "monocyte_vs_pca_5pairs.csv"))
    strat = pd.read_csv(os.path.join(DATA, "monocyte_source_stratified_auc.csv"))
    headline = pd.read_csv(os.path.join(DATA, "donor_level_auc_repeated.csv"))
    headline = headline[(headline["aggregation"] == "mean") & (headline["classifier"] == "LR")].set_index("fm")

    for m, v in SPEC_ANCHOR.items():
        check_close("Figure4v3", f"{m} monocyte AUC", v, headline.loc[m, "auc_repeat_mean"],
                    os.path.join(DATA, "donor_level_auc_repeated.csv"), tol=0.001)

    fig = plt.figure(figsize=(15, 4.5))
    ax1 = fig.add_subplot(1, 3, 1)
    ax2 = fig.add_subplot(1, 3, 2)
    ax3 = fig.add_subplot(1, 3, 3)
    panel_a(ax1, df)
    panel_b(ax2, pw)
    panel_c(ax3, strat)
    fig.tight_layout()
    fig.text(0.5, -0.05,
              "No foundation model significantly beats the from-scratch baseline in-distribution "
              "(seed-ensemble bootstrap, all p > 0.23); every method, including the baseline, scores "
              "lower on monocytes than on platelets (task property, not a foundation-model property).",
              ha="center", fontsize=7.5, color="dimgray")
    save_fig(fig, "Figure4_indist_vs_ood")


if __name__ == "__main__":
    main()
