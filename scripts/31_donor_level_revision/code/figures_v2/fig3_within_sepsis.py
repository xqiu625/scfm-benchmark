"""Figure 3 — within-sepsis null (91 patients / 17 fatal). Single panel violin+jitter+mean."""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from style import MODEL_COLOR, MODEL_ORDER, check_close, save_fig

DATA = "scripts/31_donor_level_revision/results_patient_v2/within_sepsis"
ORDER = MODEL_ORDER + ["PCA-200"]

SPEC_ANCHOR = {"UCE": 0.548, "scPRINT-2": 0.523, "Geneformer": 0.515,
               "STATE": 0.511, "PCA-200": 0.428}


def main():
    per_seed = pd.read_csv(os.path.join(DATA, "per_seed_auc_repeated.csv"))
    headline = pd.read_csv(os.path.join(DATA, "donor_level_auc_repeated.csv"))
    lr_mean = headline[(headline["aggregation"] == "mean") & (headline["classifier"] == "LR")].set_index("fm")

    for m in SPEC_ANCHOR:
        check_close("Figure3", f"{m} AUC", SPEC_ANCHOR[m], lr_mean.loc[m, "auc_repeat_mean"],
                    os.path.join(DATA, "donor_level_auc_repeated.csv"))

    # NOTE (v2.1, updated post-Job-Q): the per-seed Wilcoxon q-values previously
    # annotated here ("q=0.0037 vs 4/5 FMs") are RETIRED. Job Q's seed-ensemble
    # paired bootstrap (B=10,000, all 15 pairs) landed 2026-09-08:
    # paired_inference_full/within_sepsis_15pairs.csv. Result: NO pair reaches
    # significance after accounting for 15 comparisons; the sole nominal (uncorrected)
    # exception is UCE-scGPT (p=0.047). In particular PCA-200 is NOT significantly
    # different from any FM under honest inference -- the earlier "PCA sig. worse
    # than 4/5 FMs" claim does not survive the seed-ensemble bootstrap and is
    # retired without replacement, not just re-labeled.
    pairs_q = pd.read_csv(
        "scripts/31_donor_level_revision/results_patient_v2/paired_inference_full/within_sepsis_15pairs.csv")
    n_sig = int((pairs_q["p_boot"] < 0.05).sum())
    assert n_sig == 1 and pairs_q.loc[pairs_q["p_boot"].idxmin(), "pair"] == "UCE - scGPT", (
        "Figure3 note text assumes exactly 1 nominal pair (UCE-scGPT); Job Q data changed, update note.")

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    data = []
    for m in ORDER:
        seeds = per_seed[(per_seed["fm"] == m) & (per_seed["aggregation"] == "mean") &
                          (per_seed["classifier"] == "LR")]["auc"].to_numpy()
        data.append(seeds)
    parts = ax.violinplot(data, positions=np.arange(len(ORDER)), showextrema=False, widths=0.7)
    for i, (pc, m) in enumerate(zip(parts["bodies"], ORDER)):
        pc.set_facecolor(MODEL_COLOR[m])
        pc.set_alpha(0.35)
        pc.set_edgecolor(MODEL_COLOR[m])
    rng = np.random.default_rng(2)
    for i, (m, seeds) in enumerate(zip(ORDER, data)):
        jitter = rng.normal(0, 0.05, size=len(seeds))
        ax.scatter(i + jitter, seeds, s=14, color=MODEL_COLOR[m], alpha=0.7, zorder=3,
                   edgecolors="none")
        ax.scatter([i], [seeds.mean()], s=90, color=MODEL_COLOR[m], edgecolors="black",
                  linewidths=1, zorder=4, marker="D")
    ax.axhline(0.5, color="gray", ls="--", lw=1, zorder=0)
    ax.set_xticks(np.arange(len(ORDER)))
    ax.set_xticklabels(ORDER)
    ax.set_ylabel("Per-seed patient-level AUC (LR)")
    ax.set_ylim(0.25, 0.75)
    ax.text(0.02, 0.98, "28-day mortality; draws median 16 d post admission\n"
            "No pairwise comparison significant (seed-ensemble bootstrap, 15 pairs, B=10,000);\n"
            "UCE−scGPT nominal only (p=0.047, uncorrected).",
            transform=ax.transAxes, ha="left", va="top", fontsize=7.5, style="italic",
            color="dimgray")
    fig.tight_layout()
    save_fig(fig, "Figure3_within_sepsis")


if __name__ == "__main__":
    main()
