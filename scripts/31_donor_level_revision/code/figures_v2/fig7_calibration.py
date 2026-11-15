"""Figure 7 — calibration, within-COVID, patient-level. Two-panel reliability diagram.

A = LR (donor_predictions_<FM>.csv 'mean' column), B = RF ('mean_rf').
10 equal-width bins; point area proportional to bin population; ECE = population-
weighted mean |predicted - observed| over the 10 bins (matches Eq. 3 in Methods).
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from style import MODEL_COLOR, MODEL_ORDER, save_fig

DATA = "scripts/31_donor_level_revision/results_patient_v2/within_covid"
N_BINS = 10


def reliability(y, p, n_bins=N_BINS):
    bins = np.linspace(0, 1, n_bins + 1)
    idx = np.clip(np.digitize(p, bins) - 1, 0, n_bins - 1)
    bin_mean_pred, bin_frac_pos, bin_n = [], [], []
    for b in range(n_bins):
        mask = idx == b
        n = mask.sum()
        if n == 0:
            continue
        bin_mean_pred.append(p[mask].mean())
        bin_frac_pos.append(y[mask].mean())
        bin_n.append(n)
    bin_n = np.array(bin_n)
    ece = np.sum(bin_n * np.abs(np.array(bin_mean_pred) - np.array(bin_frac_pos))) / bin_n.sum()
    return np.array(bin_mean_pred), np.array(bin_frac_pos), bin_n, ece


def panel(ax, score_col, title):
    for m in MODEL_ORDER + ["PCA-200"]:
        df = pd.read_csv(os.path.join(DATA, f"donor_predictions_{m}.csv"))
        pred, frac, n, ece = reliability(df["y"].to_numpy(), df[score_col].to_numpy())
        sizes = 15 + 60 * (n / n.max())
        ax.plot(pred, frac, color=MODEL_COLOR[m], lw=1.2, alpha=0.8, zorder=1)
        ax.scatter(pred, frac, s=sizes, color=MODEL_COLOR[m], zorder=2,
                  label=f"{m} (ECE={ece:.3f})", edgecolors="none")
    ax.plot([0, 1], [0, 1], color="gray", ls="--", lw=0.8, zorder=0)
    ax.set_xlabel("Mean predicted probability (bin)")
    ax.set_ylabel("Observed fatal fraction (bin)")
    ax.set_title(title, loc="left", fontweight="bold")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.legend(fontsize=6.5, frameon=False, loc="upper left")


def main():
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    panel(axes[0], "mean", "A  LR")
    panel(axes[1], "mean_rf", "B  RF")
    fig.tight_layout()
    save_fig(fig, "Figure7_calibration")


if __name__ == "__main__":
    main()
