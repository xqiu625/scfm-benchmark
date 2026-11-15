"""Figure 2 — within-COVID (240 patients / 25 fatal), 3 panels.

A: lollipop (borrows FigureYa19Lollipop's point-at-stem-end encoding, not its
   mutation-track layout) — repeat-mean AUC per model with per-seed jitter + sd cap.
B: single-panel multi-curve patient-level ROC (borrows FigureYa102multipanelROC's
   sensitivity-vs-100-specificity + diagonal + AUC-in-legend grammar, collapsed to
   one panel per the spec's explicit "单面板多曲线" instruction).
C: box + density (FigureYa227boxdensity, directly reused) on the 10 per-seed AUCs,
   the empirical case for why repeated-CV was adopted.
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from sklearn.metrics import roc_curve, roc_auc_score

sys.path.insert(0, os.path.dirname(__file__))
from style import DPI, MODEL_COLOR, MODEL_ORDER, check_close, save_fig, sig_stars

DATA = "scripts/31_donor_level_revision/results_patient_v2/within_covid"
METHODS_ALL = MODEL_ORDER + ["PCA-200", "null"]

# spec anchors (P0_onepager.md §3a / figure_generation_spec.md Figure 2 anchors)
SPEC_ANCHOR = {
    "Geneformer": 0.786, "STATE": 0.778, "scPRINT-2": 0.740,
    "PCA-200": 0.725, "scGPT": 0.720, "UCE": 0.661,
    "null": 0.570,  # STALE — see mismatch flag below; CSV says 0.5544
}


def load_headline():
    df = pd.read_csv(os.path.join(DATA, "donor_level_auc_repeated.csv"))
    lr_mean = df[(df["aggregation"] == "mean") & (df["classifier"] == "LR")].set_index("fm")
    null_df = pd.read_csv(os.path.join(DATA, "null_model_auc_repeated.csv"))
    null_row = null_df[(null_df["feature_set"] == "log_cell_count+disease") &
                        (null_df["variant"] == "cell-expanded")].iloc[0]
    return lr_mean, null_row


def load_per_seed():
    return pd.read_csv(os.path.join(DATA, "per_seed_auc_repeated.csv"))


def panel_a(ax, lr_mean, null_row, per_seed):
    order = sorted(METHODS_ALL, key=lambda m: (
        null_row["auc_repeat_mean"] if m == "null" else lr_mean.loc[m, "auc_repeat_mean"]))
    y = np.arange(len(order))
    for i, m in enumerate(order):
        mean = null_row["auc_repeat_mean"] if m == "null" else lr_mean.loc[m, "auc_repeat_mean"]
        sd = null_row["auc_repeat_sd"] if m == "null" else lr_mean.loc[m, "auc_repeat_sd"]
        color = MODEL_COLOR[m]
        ax.plot([0, mean], [i, i], color=color, lw=1.5, zorder=1)
        if m != "null":
            seeds = per_seed[(per_seed["fm"] == m) & (per_seed["aggregation"] == "mean") &
                              (per_seed["classifier"] == "LR")]["auc"].to_numpy()
            jitter = np.random.default_rng(0).normal(0, 0.08, size=len(seeds))
            ax.scatter(seeds, i + jitter, s=10, color=color, alpha=0.45, zorder=2,
                       edgecolors="none")
        ax.errorbar([mean], [i], xerr=[sd], fmt="o", color=color, ecolor=color,
                    elinewidth=1.2, capsize=3, markersize=8, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(order)
    ax.axvline(0.5, color="gray", ls=":", lw=0.8, zorder=0)
    ax.set_xlabel("Patient-level repeat-mean AUC (LR)")
    ax.set_xlim(0.35, 0.9)
    ax.set_title("A", loc="left", fontweight="bold")
    # significance annotations: Job L seed-ensemble paired bootstrap (B=10,000),
    # NOT the retired per-seed Wilcoxon q from donor_level_pairwise_repeated.csv.
    pb = pd.read_csv(os.path.join(
        "scripts/31_donor_level_revision/results_patient_v2/paired_inference",
        "seed_ensemble_paired_bootstrap.csv")).set_index("pair")

    def fmt_pair(pair, label, marginal=False):
        row = pb.loc[pair]
        p = row["p_boot"]
        star = "" if marginal else (" *" if p < 0.05 else " n.s.")
        tag = " (marginal)" if marginal else star
        return f"{label}: Δ={row['delta_auc']:+.3f} [{row['ci_lo']:.3f},{row['ci_hi']:.3f}] p={p:.3f}{tag}"

    lines = [
        fmt_pair("Geneformer - PCA-200", "Geneformer-PCA"),
        fmt_pair("STATE - PCA-200", "STATE-PCA"),
        fmt_pair("STATE - Geneformer", "STATE-Geneformer (could not be distinguished)"),
        fmt_pair("scPRINT-2 - PCA-200", "scPRINT-2-PCA"),
        fmt_pair("scGPT - PCA-200", "scGPT-PCA"),
        fmt_pair("UCE - PCA-200", "UCE-PCA", marginal=True),
    ]
    ax.text(0.98, 0.02, "\n".join(lines),
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6.3, color="dimgray")


def panel_b(ax, lr_mean):
    for m in MODEL_ORDER + ["PCA-200"]:
        pred = pd.read_csv(os.path.join(DATA, f"donor_predictions_{m}.csv"))
        fpr, tpr, _ = roc_curve(pred["y"], pred["mean"])
        auc_recomp = roc_auc_score(pred["y"], pred["mean"])
        anchor = lr_mean.loc[m, "auc_repeat_mean"]
        # NB: this is a single-seed-42 point AUC (donor_predictions is the seed-42
        # aggregation), so it will not exactly equal the 10-seed repeat-mean; that is
        # expected and not a mismatch -- the legend reports the repeat-mean per spec.
        ax.plot(fpr, tpr, color=MODEL_COLOR[m], lw=1.6,
                label=f"{m} (repeat-mean AUC={anchor:.3f})")
    ax.plot([0, 1], [0, 1], color="gray", ls="--", lw=0.8)
    ax.set_xlabel("1 - Specificity")
    ax.set_ylabel("Sensitivity")
    ax.set_title("B", loc="left", fontweight="bold")
    ax.legend(fontsize=6.5, loc="lower right", frameon=False)


def panel_c(ax, per_seed):
    order = MODEL_ORDER + ["PCA-200"]
    positions = np.arange(len(order))
    box_data = []
    for i, m in enumerate(order):
        seeds = per_seed[(per_seed["fm"] == m) & (per_seed["aggregation"] == "mean") &
                          (per_seed["classifier"] == "LR")]["auc"].to_numpy()
        box_data.append(seeds)
        color = MODEL_COLOR[m]
        # density (top half), scaled+placed above the box row
        if seeds.std() > 1e-8:
            kde = gaussian_kde(seeds, bw_method=0.5)
            xs = np.linspace(seeds.min() - 0.02, seeds.max() + 0.02, 200)
            dens = kde(xs)
            dens = dens / dens.max() * 0.35
            ax.fill_between(xs, i + 0.15, i + 0.15 + dens, color=color, alpha=0.35, lw=0)
        jitter = np.random.default_rng(1).normal(0, 0.05, size=len(seeds))
        ax.scatter(seeds, i - 0.15 + jitter, s=10, color=color, alpha=0.6, zorder=3,
                   edgecolors="none")
    bp = ax.boxplot(box_data, positions=positions - 0.15, vert=False, widths=0.18,
                     patch_artist=True, showfliers=False, zorder=2)
    for patch, m in zip(bp["boxes"], order):
        patch.set_facecolor(MODEL_COLOR[m])
        patch.set_alpha(0.5)
    ax.set_yticks(positions)
    ax.set_yticklabels(order)
    ax.set_xlabel("Per-seed AUC (LR, 10 seeds)")
    ax.set_title("C  (fold instability across seeds — the motivation for repeated-CV)",
                 loc="left", fontweight="bold", fontsize=9)


def main():
    lr_mean, null_row = load_headline()
    per_seed = load_per_seed()

    for m in SPEC_ANCHOR:
        csv_val = null_row["auc_repeat_mean"] if m == "null" else lr_mean.loc[m, "auc_repeat_mean"]
        check_close("Figure2", f"{m} AUC", SPEC_ANCHOR[m], csv_val,
                    os.path.join(DATA, "donor_level_auc_repeated.csv" if m != "null"
                                 else "null_model_auc_repeated.csv"),
                    note=("spec/P0 null anchor (0.570) matches the PRE-Task-4-fix, "
                          "pre-repeated-CV number; current repeated-CV value is 0.5544 "
                          "-- using the CSV value in the figure, flagging the stale anchor."
                          if m == "null" else ""))

    fig = plt.figure(figsize=(14, 4.2))
    ax1 = fig.add_subplot(1, 3, 1)
    ax2 = fig.add_subplot(1, 3, 2)
    ax3 = fig.add_subplot(1, 3, 3)
    panel_a(ax1, lr_mean, null_row, per_seed)
    panel_b(ax2, lr_mean)
    panel_c(ax3, per_seed)
    fig.tight_layout()
    fig.text(0.5, -0.04,
              "Contrasts use patient-level paired bootstrap (B=10,000) on the "
              "seed-ensemble out-of-fold scores.",
              ha="center", fontsize=7.5, color="dimgray")
    save_fig(fig, "Figure2_within_covid")


if __name__ == "__main__":
    main()
