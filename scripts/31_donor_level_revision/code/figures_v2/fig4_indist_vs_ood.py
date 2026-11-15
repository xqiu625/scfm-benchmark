"""Figure 4 — in-distribution (monocyte) vs OOD (platelet). 3 panels.

A: dumbbell/slope (monocyte COVID-only vs platelet, per model, sorted by platelet value)
B: forest plot (FigureYa90, direct reuse) of FM-vs-PCA-200 monocyte deltas
C: source-stratified point plot (Ren / Stephenson / Qiu-sepsis)
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from style import MODEL_COLOR, MODEL_ORDER, check_close, save_fig

DATA = "scripts/31_donor_level_revision/results_patient_v2/monocyte"


def panel_a(ax, df):
    covid_only = df[df["cut"] == "COVID_only_Ren_Stephenson"].set_index("fm")
    order = covid_only.sort_values("auc_platelet_mean", ascending=True).index.tolist()
    for i, m in enumerate(order):
        row = covid_only.loc[m]
        color = MODEL_COLOR[m]
        ax.plot([row["auc_monocyte_mean"], row["auc_platelet_mean"]], [i, i],
                color=color, lw=1.5, zorder=1, alpha=0.7)
        ax.scatter([row["auc_monocyte_mean"]], [i], s=70, color="white",
                  edgecolors=color, linewidths=2, zorder=2, label="Monocyte" if i == 0 else None)
        ax.scatter([row["auc_platelet_mean"]], [i], s=70, color=color, zorder=3,
                  label="Platelet" if i == 0 else None)
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels(order)
    ax.set_xlabel("Patient-level repeat-mean AUC")
    ax.set_title("A", loc="left", fontweight="bold")
    ax.axvline(0.5, color="gray", ls=":", lw=0.8, zorder=0)
    ax.legend(loc="lower right", frameon=False, fontsize=8)
    ax.text(0.02, 0.02, "monocyte: 143 pts/17 fatal (Ren+Stephenson)\nplatelet: 240 pts/25 fatal",
            transform=ax.transAxes, fontsize=7, color="dimgray", va="bottom")


def panel_b(ax, pw):
    # v2.1 (updated post-Job-Q, 2026-09-08): significance now from
    # paired_inference_full/monocyte_vs_pca_5pairs.csv (seed-ensemble paired
    # bootstrap, B=10,000, n=143 COVID-only). The per-seed p-values previously
    # shown here (STATE 0.002, scGPT 0.037 -- both "significant") are RETIRED and
    # do NOT survive honest inference: under the seed-ensemble bootstrap, NO
    # FM-vs-PCA-200 comparison reaches significance in-distribution (all p>=0.23).
    # This reverses the earlier framing -- report plainly, not smoothed over.
    order = pw.sort_values("delta_auc").reset_index(drop=True)
    y = np.arange(len(order))
    for i, row in order.iterrows():
        m = row["pair"].split(" - ")[0]
        color = MODEL_COLOR.get(m, "black")
        sig = row["p_boot"] < 0.05
        ax.plot([row["ci_lo"], row["ci_hi"]], [i, i], color=color, lw=1.4,
               alpha=0.9 if sig else 0.5, zorder=1)
        ax.scatter([row["delta_auc"]], [i], s=80, color=color,
                  marker="D", zorder=2, alpha=0.9 if sig else 0.5, edgecolors="none")
        ax.text(row["ci_hi"] + 0.01, i, f"p={row['p_boot']:.2f}", va="center",
               fontsize=6.5, color="dimgray")
    ax.axvline(0, color="gray", ls="--", lw=1, zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels([r["pair"].split(" - ")[0] for _, r in order.iterrows()])
    ax.set_xlabel("Δ AUC (FM − PCA-200), monocyte COVID-only")
    ax.set_title("B  (no pair significant, seed-ensemble bootstrap)", loc="left",
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
    pw = pd.read_csv(
        "scripts/31_donor_level_revision/results_patient_v2/paired_inference_full/monocyte_vs_pca_5pairs.csv")
    strat = pd.read_csv(os.path.join(DATA, "monocyte_source_stratified_auc.csv"))
    assert (pw["n_patients"] == 143).all(), "Job Q monocyte pairs should be the 143-patient COVID-only cut"
    assert (pw["p_boot"] >= 0.05).all(), (
        "Figure4 panel B text assumes no significant pair; Job Q data changed, update panel_b wording")

    spec_anchor = {"STATE": -0.0914, "Geneformer": -0.1191, "scPRINT-2": -0.0850,
                   "scGPT": -0.0554, "UCE": -0.0047, "PCA-200": -0.0714}
    covid_only = df[df["cut"] == "COVID_only_Ren_Stephenson"].set_index("fm")
    for m, spec_v in spec_anchor.items():
        check_close("Figure4", f"{m} delta (monocyte-platelet)", spec_v,
                    covid_only.loc[m, "delta_monocyte_minus_platelet"],
                    os.path.join(DATA, "indistribution_vs_platelet.csv"))

    fig = plt.figure(figsize=(15, 4.5))
    ax1 = fig.add_subplot(1, 3, 1)
    ax2 = fig.add_subplot(1, 3, 2)
    ax3 = fig.add_subplot(1, 3, 3)
    panel_a(ax1, df)
    panel_b(ax2, pw)
    panel_c(ax3, strat)
    fig.tight_layout()
    fig.text(0.5, -0.04,
              "UCE favors monocytes across all four configurations (seed × "
              "disease-pool); other models revert to the platelet-favoring "
              "direction (or zero) under the 10-seed ensemble.",
              ha="center", fontsize=7.5, color="dimgray")
    save_fig(fig, "Figure4_indist_vs_ood")


if __name__ == "__main__":
    main()
