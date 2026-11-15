"""Figure S2 -- within-COVID full pairwise (15 pairs) delta-AUC, lower-triangle
heatmap. Magnitude from donor_level_pairwise_repeated.csv (unchanged, per
figure_generation_spec.md 修订附录 v2.1). Significance layer (updated
2026-09-08, Job Q landed): stars from paired_inference_full/within_covid_15pairs.csv
(seed-ensemble paired bootstrap, B=10,000) -- replaces the retired per-seed
Wilcoxon q-values.
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from style import MODEL_ORDER, save_fig

DATA = "scripts/31_donor_level_revision/results_patient_v2/within_covid"


def sig_lookup(pw_q):
    d = {}
    for _, row in pw_q.iterrows():
        a, b = row["pair"].split(" - ")
        d[frozenset([a, b])] = row["p_boot"]
    return d


def main():
    pw = pd.read_csv(os.path.join(DATA, "donor_level_pairwise_repeated.csv"))
    pw_q = pd.read_csv(
        "scripts/31_donor_level_revision/results_patient_v2/paired_inference_full/within_covid_15pairs.csv")
    assert len(pw_q) == 15, "Job Q within-COVID pairs should be exactly 15"
    p_by_pair = sig_lookup(pw_q)

    methods = MODEL_ORDER + ["PCA-200"]
    n = len(methods)
    mat = np.full((n, n), np.nan)
    for _, row in pw.iterrows():
        a, b = row["pair"].split(" - ")
        if a not in methods or b not in methods:
            continue
        ia, ib = methods.index(a), methods.index(b)
        i, j = max(ia, ib), min(ia, ib)  # lower triangle: row idx > col idx
        val = row["mean_delta"] if ia > ib else -row["mean_delta"]
        mat[i, j] = val

    fig, ax = plt.subplots(figsize=(6, 5.5))
    vmax = np.nanmax(np.abs(mat))
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    for i in range(n):
        for j in range(n):
            if np.isnan(mat[i, j]):
                continue
            p = p_by_pair.get(frozenset([methods[i], methods[j]]))
            star = "*" if (p is not None and p < 0.05) else ""
            ax.text(j, i, f"{mat[i, j]:+.3f}{star}", ha="center", va="center", fontsize=8,
                   color="white" if abs(mat[i, j]) > vmax * 0.5 else "black")
    ax.set_xticks(range(n))
    ax.set_xticklabels(methods, rotation=45, ha="right")
    ax.set_yticks(range(n))
    ax.set_yticklabels(methods)
    ax.set_title("Within-COVID pairwise ΔAUC (row − column), repeated-CV\n"
                  "* = p<0.05, seed-ensemble paired bootstrap (Job Q, B=10,000)", fontsize=10)
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("Δ AUC (mean, 10 seeds)")
    fig.tight_layout()
    save_fig(fig, "FigureS2_pairwise_heatmap")


if __name__ == "__main__":
    main()
