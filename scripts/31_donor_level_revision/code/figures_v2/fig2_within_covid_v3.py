"""Figure 2 — within-COVID (240 patients / 25 fatal), v3 spec, 2 panels.

Per reports/manuscript/revision/figure_table_spec_v3.md Sec.1:
A: six methods + three baselines, patient-level repeat-mean +/- repeat-sd,
   dot-and-whisker, sorted by value. Baselines (null/cohort-only/QC+count) as
   gray squares, five FMs as colored circles, PCA-200 as an open circle.
   Annotated: pre-specified primary (Geneformer vs PCA-200, winner-corrected
   permutation P=0.181) and the confirmatory leaders-vs-UCE range-statistic
   (P=0.011).
B: incremental-value forest plot (baseline = QC+count+cohort, AUC=0.627),
   points = DeLong ΔAUC, whiskers = 1.96*DeLong SE, text columns for DeLong p
   and nested-LR P; UCE row grayed (weakest under both readings).
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from style import MODEL_COLOR, MODEL_ORDER, check_close, save_fig

DATA_WC = "scripts/31_donor_level_revision/results_patient_v2/within_covid"
DATA_PI = "scripts/31_donor_level_revision/results_patient_v2/paired_inference_full"

# P0_onepager.md (round-17 final) Sec.3a anchors, verified against CSVs below
SPEC_ANCHOR_MAIN = {
    "Geneformer": 0.786, "STATE": 0.778, "scPRINT-2": 0.740,
    "PCA-200": 0.725, "scGPT": 0.720, "UCE": 0.661,
}
SPEC_ANCHOR_NULL = 0.554
SPEC_ANCHOR_COHORT = 0.486
SPEC_ANCHOR_QC = 0.618

INCR_ANCHOR = {
    # fm: (delta, delong_p, nested_lr_note)
    "STATE": (0.065, 0.00015, "P<0.0001"),
    "Geneformer": (0.053, 0.00014, "P<0.0001"),
    "PCA-200": (0.050, 0.0034, "P<0.0001 (pre-specified primary)"),
    "scPRINT-2": (0.046, 0.0054, "P<0.0001"),
    "scGPT": (0.026, 0.0063, "P<0.0001"),
    "UCE": (0.016, 0.072, "P=0.0022"),
}


def load_panel_a():
    df = pd.read_csv(os.path.join(DATA_WC, "donor_level_auc_repeated.csv"))
    lr_mean = df[(df["aggregation"] == "mean") & (df["classifier"] == "LR")].set_index("fm")
    null_df = pd.read_csv(os.path.join(DATA_WC, "null_model_auc_repeated.csv"))
    null_row = null_df[(null_df["feature_set"] == "log_cell_count+disease") &
                        (null_df["variant"] == "cell-expanded")].iloc[0]
    cohort_row = pd.read_csv(os.path.join(DATA_PI, "cohort_only_baseline.csv")).iloc[0]
    return lr_mean, null_row, cohort_row


def panel_a(ax, lr_mean, null_row, cohort_row):
    rows = {}
    for m in MODEL_ORDER + ["PCA-200"]:
        rows[m] = (lr_mean.loc[m, "auc_repeat_mean"], lr_mean.loc[m, "auc_repeat_sd"], "fm")
    rows["null (count+disease)"] = (null_row["auc_repeat_mean"], null_row["auc_repeat_sd"], "baseline")
    rows["cohort-only"] = (cohort_row["repeat_mean"], cohort_row["repeat_sd"], "baseline")
    rows["QC+cell count (as-fit)"] = (SPEC_ANCHOR_QC, 0.0, "baseline")

    order = sorted(rows, key=lambda k: rows[k][0])
    y = np.arange(len(order))
    for i, name in enumerate(order):
        mean, sd, kind = rows[name]
        if kind == "baseline":
            color = "#7F8C8D"
            marker = "s"
            face = color
        elif name == "PCA-200":
            color = MODEL_COLOR["PCA-200"]
            marker = "o"
            face = "white"
        else:
            color = MODEL_COLOR[name]
            marker = "o"
            face = color
        ax.plot([0, mean], [i, i], color=color, lw=1.2, zorder=1, alpha=0.6)
        if sd > 0:
            ax.errorbar([mean], [i], xerr=[sd], fmt="none", ecolor=color,
                        elinewidth=1.2, capsize=3, zorder=2)
        ax.scatter([mean], [i], s=70, color=face, edgecolors=color, linewidths=1.4,
                   marker=marker, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(order, fontsize=8.5)
    ax.axvline(0.5, color="gray", ls=":", lw=0.8, zorder=0)
    ax.set_xlabel("Patient-level repeat-mean AUC (ensemble-OOF, LR)")
    ax.set_xlim(0.35, 0.92)
    ax.set_title("A", loc="left", fontweight="bold")

    geneformer_y = order.index("Geneformer")
    uce_y = order.index("UCE")
    ax.text(0.905, geneformer_y, "pre-specified: Geneformer vs PCA-200\nwinner-corrected permutation P = 0.181 (n.s.)",
            ha="left", va="center", fontsize=6.6, color="dimgray")
    ax.text(0.905, uce_y, "confirmatory: leaders vs UCE\nrange-statistic permutation P = 0.011",
            ha="left", va="center", fontsize=6.6, color="dimgray")


def panel_b(ax, incr_df):
    order = incr_df.sort_values("incremental_delta_auc").index.tolist()
    y = np.arange(len(order))
    ax.set_xlim(-0.01, 0.10)
    for i, fm in enumerate(order):
        row = incr_df.loc[fm]
        delta = row["incremental_delta_auc"]
        se = row["delong_SE"]
        color = "#7F8C8D" if fm == "UCE" else MODEL_COLOR.get(fm, "#333333")
        marker = "o" if fm != "PCA-200" else "D"
        ax.errorbar([delta], [i], xerr=[1.96 * se], fmt="none", ecolor=color,
                    elinewidth=1.3, capsize=3, zorder=2)
        ax.scatter([delta], [i], s=70, color=color, marker=marker, zorder=3,
                   edgecolors="black" if fm == "UCE" else "none", linewidths=0.6)
        p_delong = row["delong_p"]
        nested_note = INCR_ANCHOR[fm][2]
        ax.text(1.02, i, f"DeLong p={p_delong:.4f}   nested-LR {nested_note}",
                transform=ax.get_yaxis_transform(), va="center", ha="left",
                fontsize=6.6, color="dimgray")
    ax.axvline(0, color="gray", ls=":", lw=0.8, zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels(order, fontsize=8.5)
    ax.set_xlabel(r"$\Delta$AUC over QC+count+cohort baseline (AUC=0.627)")
    ax.set_title("B", loc="left", fontweight="bold")


def main():
    lr_mean, null_row, cohort_row = load_panel_a()
    incr_df = pd.read_csv(os.path.join(DATA_PI, "incremental_value_test.csv")).set_index("fm")

    for m, v in SPEC_ANCHOR_MAIN.items():
        check_close("Figure2v3", f"{m} AUC", v, lr_mean.loc[m, "auc_repeat_mean"],
                    os.path.join(DATA_WC, "donor_level_auc_repeated.csv"))
    check_close("Figure2v3", "null AUC", SPEC_ANCHOR_NULL, null_row["auc_repeat_mean"],
                os.path.join(DATA_WC, "null_model_auc_repeated.csv"), tol=0.001)
    check_close("Figure2v3", "cohort-only AUC", SPEC_ANCHOR_COHORT, cohort_row["repeat_mean"],
                os.path.join(DATA_PI, "cohort_only_baseline.csv"), tol=0.001)
    for fm, (delta, p, _) in INCR_ANCHOR.items():
        check_close("Figure2v3", f"{fm} incremental delta", delta,
                    incr_df.loc[fm, "incremental_delta_auc"],
                    os.path.join(DATA_PI, "incremental_value_test.csv"), tol=0.001)
        check_close("Figure2v3", f"{fm} incremental delong_p", p,
                    incr_df.loc[fm, "delong_p"],
                    os.path.join(DATA_PI, "incremental_value_test.csv"), tol=0.0005)

    fig = plt.figure(figsize=(13, 5.2))
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)
    panel_a(ax1, lr_mean, null_row, cohort_row)
    panel_b(ax2, incr_df)
    fig.tight_layout()
    fig.text(0.5, -0.04,
              "Panel B baseline = QC metrics + cell count + cohort (AUC 0.627); increments over it "
              "by exact DeLong; nested likelihood-ratio P against an empirical permutation null; "
              "the pre-specified primary for this question is the PCA-200 row (nested-LR P < 0.0001).",
              ha="center", fontsize=7.3, color="dimgray", wrap=True)
    save_fig(fig, "Figure2_within_covid")


if __name__ == "__main__":
    main()
