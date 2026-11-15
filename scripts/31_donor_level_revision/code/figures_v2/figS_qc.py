"""Figure S-QC -- QC-metric distributions by outcome/disease + the QC-only
baseline vs. within-COVID model AUCs (evidence for external-review point A1:
QC metrics alone do not explain the mortality signal).

Panels A/B/C: box+density (FigureYa227boxdensity grammar, as used in Figure 2C)
for median_counts (nUMI, log y), median_genes, median_pct_mito -- each split by
outcome (survived/fatal) x disease (COVID-19/Sepsis), 4 groups per panel.
Panel D: lollipop, QC-only AUC vs. within-COVID model AUCs, 0.5 reference line.

COHORT-SCOPE NOTE (resolved, not assumed): qc_baseline_auc.csv's repeated-CV was
fit on the POOLED 252-patient cohort (COVID 240 + legacy sepsis 12, pre-Chongqing/
pre-Job-J) -- confirmed via patient_qc_distributions.csv's own disease column
(COVID-19: 240 rows, Sepsis: 12 rows, sums to 252). No per-patient OOF/prediction
file was saved for the QC-only model, so a true within-COVID-only re-stratified
AUC CANNOT be computed post-hoc without retraining (out of scope for a figure
script). Per the spec's own explicit fallback, Panel D's QC-only AUC is captioned
"pooled analysis cohort (252 patients, 27 fatal)" -- NOT mislabeled as within-COVID.
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

sys.path.insert(0, os.path.dirname(__file__))
from style import DISEASE_COLOR, MODEL_COLOR, MODEL_ORDER, check_close, save_fig

DATA = "scripts/31_donor_level_revision/results_patient_v2/qc_baseline"
WITHIN_COVID = "scripts/31_donor_level_revision/results_patient_v2/within_covid"

GROUPS = [("COVID-19", 0, "COVID\nsurvived"), ("COVID-19", 1, "COVID\nfatal"),
          ("Sepsis", 0, "Sepsis\nsurvived"), ("Sepsis", 1, "Sepsis\nfatal")]


def box_density_panel(ax, df, col, ylabel, log=False, title=""):
    positions = np.arange(len(GROUPS))
    box_data = []
    colors = []
    for disease, y, _ in GROUPS:
        vals = df[(df["disease"] == disease) & (df["y"] == y)][col].to_numpy()
        vals = vals[~np.isnan(vals)]
        box_data.append(vals)
        colors.append(DISEASE_COLOR[disease])
    for i, (vals, color) in enumerate(zip(box_data, colors)):
        if len(vals) > 3 and np.std(vals) > 1e-8:
            plot_vals = np.log10(vals) if log else vals
            kde = gaussian_kde(plot_vals, bw_method=0.5)
            lo, hi = plot_vals.min(), plot_vals.max()
            pad = (hi - lo) * 0.1 if hi > lo else 1
            xs = np.linspace(lo - pad, hi + pad, 150)
            dens = kde(xs)
            dens = dens / dens.max() * 0.32
            xs_disp = 10 ** xs if log else xs
            ax.fill_betweenx(xs_disp, i + 0.12, i + 0.12 + dens, color=color, alpha=0.35, lw=0)
        rng = np.random.default_rng(i)
        jitter = rng.normal(0, 0.06, size=len(vals))
        ax.scatter(i - 0.15 + jitter, vals, s=9, color=color, alpha=0.55, zorder=3,
                  edgecolors="none")
    bp = ax.boxplot(box_data, positions=positions - 0.15, widths=0.16,
                     patch_artist=True, showfliers=False, zorder=2)
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.45)
    if log:
        ax.set_yscale("log")
    ax.set_xticks(positions)
    ax.set_xticklabels([g[2] for g in GROUPS], fontsize=8)
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", fontweight="bold")


def panel_d(ax, qc_auc_mean, qc_auc_sd, cohort_note):
    lr_mean = pd.read_csv(os.path.join(WITHIN_COVID, "donor_level_auc_repeated.csv"))
    lr_mean = lr_mean[(lr_mean["aggregation"] == "mean") & (lr_mean["classifier"] == "LR")].set_index("fm")
    methods = ["QC-only"] + MODEL_ORDER + ["PCA-200"]
    means, sds, colors = [], [], []
    for m in methods:
        if m == "QC-only":
            means.append(qc_auc_mean); sds.append(qc_auc_sd); colors.append("#95A5A6")
        else:
            means.append(lr_mean.loc[m, "auc_repeat_mean"])
            sds.append(lr_mean.loc[m, "auc_repeat_sd"])
            colors.append(MODEL_COLOR[m])
    order = np.argsort(means)
    y = np.arange(len(methods))
    for rank, idx in enumerate(order):
        m, mean, sd, color = methods[idx], means[idx], sds[idx], colors[idx]
        ax.plot([0, mean], [rank, rank], color=color, lw=1.5, zorder=1)
        ax.errorbar([mean], [rank], xerr=[sd], fmt="o", color=color, ecolor=color,
                    elinewidth=1.2, capsize=3, markersize=9 if m == "QC-only" else 7, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([methods[i] for i in order])
    ax.axvline(0.5, color="gray", ls=":", lw=0.8, zorder=0)
    ax.set_xlabel("Patient-level AUC")
    ax.set_xlim(0.35, 0.9)
    ax.set_title("D", loc="left", fontweight="bold")
    ax.text(0.98, 0.02,
            f"QC-only: pooled analysis cohort ({cohort_note})\n"
            "Fatal patients show higher nUMI and %mito,\n"
            "yet QC alone does not predict outcome (AUC 0.53).",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7, color="dimgray")


def main():
    dist = pd.read_csv(os.path.join(DATA, "patient_qc_distributions.csv"))
    qc_auc = pd.read_csv(os.path.join(DATA, "qc_baseline_auc.csv")).iloc[0]

    check_close("FigureS_QC", "QC-only AUC", 0.528, qc_auc["auc_repeat_mean"],
                os.path.join(DATA, "qc_baseline_auc.csv"))

    # cohort-scope resolution (documented in module docstring): confirm 240+12=252
    disease_counts = dist["disease"].value_counts()
    n_covid, n_sepsis = disease_counts.get("COVID-19", 0), disease_counts.get("Sepsis", 0)
    if n_covid + n_sepsis != int(qc_auc["n_patients"]):
        from style import flag_mismatch
        flag_mismatch("FigureS_QC", "disease-count sum vs n_patients",
                      qc_auc["n_patients"], n_covid + n_sepsis,
                      os.path.join(DATA, "patient_qc_distributions.csv"))
    cohort_note = f"{int(qc_auc['n_patients'])} patients, {int(qc_auc['n_fatal_patients'])} fatal"

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    box_density_panel(axes[0, 0], dist, "median_counts", "Median nUMI (log)", log=True, title="A")
    box_density_panel(axes[0, 1], dist, "median_genes", "Median n genes", title="B")
    box_density_panel(axes[1, 0], dist, "median_pct_mito", "Median % mito", title="C")
    panel_d(axes[1, 1], qc_auc["auc_repeat_mean"], qc_auc["auc_repeat_sd"], cohort_note)
    fig.suptitle("Fatal patients show higher nUMI and %mito, yet QC alone does not "
                  "predict outcome (AUC 0.53).", fontsize=10, y=1.01)
    fig.tight_layout()
    save_fig(fig, "FigureS_QC_distributions")


if __name__ == "__main__":
    main()
