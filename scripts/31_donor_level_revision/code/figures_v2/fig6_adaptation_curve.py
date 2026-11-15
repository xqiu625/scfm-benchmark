"""Figure 6 — adaptation curve, post-Jiang-exclusion (Job J) 96-patient sepsis
pool, 34-patient test set. Line + ribbon.

v2.1: test_set is now 'expanded_96pt' (was 'expanded_103pt' pre-Job-J -- Jiang's
7 legacy sepsis donors excluded from the training pool). STACK and
TranscriptFormer rows are gone entirely from this file (Job P deleted them --
they only ever had legacy-test-set rows, never a Chongqing embedding, so they
were never comparable on this axis; see the runbook's comparability-gap fix).
Endpoint (n=all) should land <=0.62 for every remaining method (verified below).
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from style import MODEL_COLOR, MODEL_ORDER, check_close, save_fig

DATA = "scripts/31_donor_level_revision/results_patient_v2/adaptation_curve_chongqing"
FORBIDDEN_VALUES = {0.897, 0.858}  # legacy TF/STACK endpoint numbers (any vintage)


def main():
    df = pd.read_csv(os.path.join(DATA, "all_rows.csv"))
    if "test_set" not in df.columns:
        raise SystemExit("test_set column missing -- wrong data source, aborting.")
    present_sets = set(df["test_set"].unique())
    if present_sets != {"expanded_96pt"}:
        raise SystemExit(f"Expected only 'expanded_96pt' rows post-Job-P, found "
                          f"{present_sets} -- data source changed, investigate before plotting.")
    exp = df[df["test_set"] == "expanded_96pt"].copy()

    present_methods = set(exp["method"].unique())
    for banned in ["STACK", "TranscriptFormer"]:
        if banned in present_methods:
            raise SystemExit(f"{banned} present -- should have been purged by Job P.")

    plotted_vals = set(exp["auc_point"].round(3))
    for bad in FORBIDDEN_VALUES:
        if any(abs(v - bad) < 0.002 for v in plotted_vals):
            raise SystemExit(f"FORBIDDEN legacy value {bad} present in expanded_96pt rows.")

    grid = sorted(exp["n_sepsis"].unique())
    # spec anchor v2.1: full-pool endpoint now <=0.62 for every method (was 0.41-0.60 pre-Job-J)
    full = exp[exp["n_sepsis_req"] == -1].groupby("method")["auc_point"].mean()
    if full.max() > 0.62:
        from style import flag_mismatch
        flag_mismatch("Figure6", "full-pool AUC max should be <=0.62", "<=0.62",
                      f"{full.max():.4f} ({full.idxmax()})",
                      os.path.join(DATA, "all_rows.csv"))

    fig, ax = plt.subplots(figsize=(7, 5))
    methods = MODEL_ORDER + ["PCA-200+LR"]
    for m in methods:
        sub = exp[exp["method"] == m]
        agg = sub.groupby("n_sepsis_req")["auc_point"].agg(["mean", "std"]).reset_index()
        agg = agg.sort_values("n_sepsis_req")
        # x-axis: log scale on n_sepsis, but n_sepsis_req=-1 means "all" -> use actual n
        x = agg["n_sepsis_req"].replace(-1, sub["n_sepsis"].max())
        x_plot = x.replace(0, 1)  # log scale can't show 0; place at 1 with a break note
        color = MODEL_COLOR.get(m, MODEL_COLOR.get(m.replace("+LR", ""), "black"))
        ax.plot(x_plot, agg["mean"], color=color, lw=1.8, marker="o", markersize=4, label=m)
        ax.fill_between(x_plot, agg["mean"] - agg["std"], agg["mean"] + agg["std"],
                        color=color, alpha=0.15)
    ax.set_xscale("log")
    ax.axhline(0.5, color="gray", ls="--", lw=1, zorder=0)
    ax.set_xlabel("Sepsis training cells (log scale; leftmost point = 0)")
    ax.set_ylabel("AUC (32-patient test set, expanded 96-patient sepsis pool)")
    ax.legend(fontsize=7.5, frameon=False, ncol=2)
    fig.tight_layout()
    save_fig(fig, "Figure6_adaptation_curve")


if __name__ == "__main__":
    main()
