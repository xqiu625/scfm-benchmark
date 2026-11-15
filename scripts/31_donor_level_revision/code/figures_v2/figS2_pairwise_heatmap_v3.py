"""Figure S2 -- within-COVID full pairwise (15 pairs) delta-AUC, lower-triangle
heatmap, v3 spec: complete 6x6 matrix (5 FMs + PCA-200; the prior version was
missing STATE and PCA-200), exact DeLong throughout (paired_inference_full/
delong_exact.csv), raw p<0.05 shown with a bold cell border, BH (m=15) shown
with a star -- only UCE-Geneformer and STATE-UCE should survive BH (q=0.030,
per P0_onepager.md), verified in-script rather than hardcoded.
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from style import MODEL_ORDER, save_fig


def benjamini_hochberg(pvals, alpha=0.05):
    """Returns boolean array, same order as pvals, True = rejected (survives BH)."""
    pvals = np.asarray(pvals)
    m = len(pvals)
    order = np.argsort(pvals)
    sorted_p = pvals[order]
    thresh = (np.arange(1, m + 1) / m) * alpha
    passed = sorted_p <= thresh
    if not passed.any():
        k = 0
    else:
        k = np.max(np.where(passed)[0]) + 1
    reject_sorted = np.zeros(m, dtype=bool)
    reject_sorted[:k] = True
    reject = np.zeros(m, dtype=bool)
    reject[order] = reject_sorted
    return reject


def main():
    df = pd.read_csv(
        "scripts/31_donor_level_revision/results_patient_v2/paired_inference_full/delong_exact.csv")
    assert len(df) == 15, "expected exactly 15 pairs (6 choose 2)"

    df["bh_pass"] = benjamini_hochberg(df["delong_exact_p"].to_numpy())
    n_bh = df["bh_pass"].sum()
    bh_pairs = df.loc[df["bh_pass"], "pair"].tolist()
    print(f"BH (m=15) survivors: {n_bh} -> {bh_pairs}")
    expected_bh = {"STATE - UCE", "UCE - Geneformer"}
    if set(bh_pairs) != expected_bh:
        from style import flag_mismatch
        flag_mismatch("FigureS2v3", "BH survivors", str(expected_bh), str(set(bh_pairs)),
                      "paired_inference_full/delong_exact.csv")

    methods = MODEL_ORDER + ["PCA-200"]
    n = len(methods)
    mat = np.full((n, n), np.nan)
    sig_raw = np.zeros((n, n), dtype=bool)
    sig_bh = np.zeros((n, n), dtype=bool)
    for _, row in df.iterrows():
        a, b = row["pair"].split(" - ")
        ia, ib = methods.index(a), methods.index(b)
        i, j = max(ia, ib), min(ia, ib)
        val = row["delta_auc"] if ia > ib else -row["delta_auc"]
        mat[i, j] = val
        if row["delong_exact_p"] < 0.05:
            sig_raw[i, j] = True
        if row["bh_pass"]:
            sig_bh[i, j] = True

    fig, ax = plt.subplots(figsize=(6.3, 5.8))
    vmax = np.nanmax(np.abs(mat))
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    for i in range(n):
        for j in range(n):
            if np.isnan(mat[i, j]):
                continue
            star = "*" if sig_bh[i, j] else ""
            ax.text(j, i, f"{mat[i, j]:+.3f}{star}", ha="center", va="center", fontsize=8,
                   color="white" if abs(mat[i, j]) > vmax * 0.5 else "black")
            if sig_raw[i, j]:
                rect = plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                     edgecolor="black", linewidth=1.6, zorder=3)
                ax.add_patch(rect)
    ax.set_xticks(range(n))
    ax.set_xticklabels(methods, rotation=45, ha="right")
    ax.set_yticks(range(n))
    ax.set_yticklabels(methods)
    ax.set_title("Within-COVID pairwise $\\Delta$AUC (row − column), exact DeLong\n"
                  f"bold border = p<0.05; * = survives BH (m=15, {n_bh}/15 pairs)", fontsize=10)
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label(r"$\Delta$ AUC (ensemble-OOF)")
    fig.tight_layout()
    save_fig(fig, "FigureS2_pairwise_heatmap")


if __name__ == "__main__":
    main()
