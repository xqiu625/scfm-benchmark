#!/usr/bin/env python3
# =============================================================================
# DEPRECATED (2026-08-06): this script is superseded by
#   scripts/26_figures/code/regenerate_indexed_figures.py
# which is the canonical figure generator (its output filenames line up with
# the manuscript's \includegraphics numbering). This file is kept only because
# docs still reference it; do not use it for new figure runs.
# =============================================================================
"""
Regenerate the 5 data-driven main-text figures for the 302006 manuscript
under the 5-FM primary benchmark (STATE, UCE, Geneformer, scPRINT-2, scGPT).

Style follows the nature-figure skill: Arial, font.size=7, editable SVG/PDF,
no top/right spines, axis linewidth 0.8, NMI-pastel-style unified palette
with green/red reserved for directional cues only.

STACK and TranscriptFormer are dropped from every figure per the 2026-05-26
strategic restructuring. scGPT is shown on Fig 3a only (with its legacy
anti-correlated cross-disease value 0.352) because the Option C rerun on
the cleaner cohort is still pending for the other endpoints.
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Nature-figure style
# ---------------------------------------------------------------------------
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "legend.frameon": False,
    "legend.fontsize": 6.2,
    "axes.labelsize": 7.5,
    "axes.titlesize": 8,
})

ROOT = Path("/Users/xinruqiu/Desktop/302006/302006_updategit")
DATA = ROOT / "scripts/26_figures/data"
OUT = ROOT / "reports/manuscript/figures"
OUT.mkdir(parents=True, exist_ok=True)

# 5-FM primary set, in plot order
FMS = ["STATE", "UCE", "Geneformer", "scPRINT-2", "scGPT"]

# Drop-list (removed from every CSV before plotting)
DROP = {"STACK", "STACK ortholog", "TranscriptFormer", "TF", "TF-Metazoa"}

# Unified low-saturation palette (NMI-pastel family).
# Green/red reserved for directional cues (Fig 3a bidirectional, Fig 4 heatmap).
PALETTE = {
    "STATE":      "#3B6BB0",   # restrained blue
    "UCE":        "#5C9D7F",   # muted teal-green (allowed: identity, not directional)
    "Geneformer": "#7C5C9E",   # muted purple
    "scPRINT-2":  "#D4944B",   # warm ochre
    "scGPT":      "#A84A4A",   # restrained red (reserved for the anti-correlated reference)
    "PCA-200+LR": "#444444",   # neutral grey for baseline
}

# CSV-name harmonization (some CSVs use scprint2 / scPRINT-2 / scPRINT2 etc.)
def canon(name: str) -> str:
    s = str(name).strip().lower().replace("_", "").replace("-", "")
    if s in {"scprint2", "scprint", "scprint2native"}:
        return "scPRINT-2"
    if s == "state":
        return "STATE"
    if s == "uce":
        return "UCE"
    if s == "geneformer":
        return "Geneformer"
    if s == "scgpt":
        return "scGPT"
    if s.startswith("stack"):
        return "STACK"
    if s.startswith("transcriptformer") or s.startswith("tf"):
        return "TranscriptFormer"
    return name


def save_pub(fig, stem):
    fig.savefig(OUT / f"{stem}.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Fig 2 — Per-FM mortality AUC bars (LR + RF), PCA baseline overlay
# ---------------------------------------------------------------------------
def make_fig2_performance():
    df = pd.read_csv(DATA / "mortality_lr_rf.csv")
    df["fm"] = df["fm"].map(canon)
    df = df[~df["fm"].isin(DROP)]
    df = df[df["fm"].isin(FMS)]

    pca = pd.read_csv(DATA / "pca_baseline.csv")
    pca_lr200 = pca.loc[pca["method"] == "PCA-200+LR", "auc_mean"].iloc[0]

    # Pivot to FM × classifier
    fms_present = [f for f in FMS if f in df["fm"].unique()]
    classifiers = ["LR", "RF"]
    width = 0.36
    x = np.arange(len(FMS))

    fig, ax = plt.subplots(figsize=(3.6, 2.4))

    for i, clf in enumerate(classifiers):
        means, stds, colors = [], [], []
        for fm in FMS:
            row = df[(df["fm"] == fm) & (df["classifier"] == clf)]
            if len(row) == 0:
                means.append(np.nan)
                stds.append(0)
            else:
                means.append(row["auc_mean"].iloc[0])
                stds.append(row["auc_std"].iloc[0])
            colors.append(PALETTE.get(fm, "#888"))
        offset = (i - 0.5) * width
        ax.bar(
            x + offset,
            means,
            width,
            yerr=stds,
            color=[PALETTE.get(fm, "#888") for fm in FMS],
            edgecolor="black",
            linewidth=0.4,
            error_kw=dict(elinewidth=0.6, capsize=1.5, capthick=0.6, ecolor="#444"),
            alpha=0.55 if clf == "RF" else 1.0,
            label=clf,
        )

    # PCA-200+LR horizontal reference (label inside plot, left side, well below the title)
    ax.axhline(pca_lr200, ls="--", lw=0.8, color=PALETTE["PCA-200+LR"], zorder=0)
    ax.text(-0.45, pca_lr200, f"PCA-200+LR = {pca_lr200:.3f}",
            ha="left", va="bottom", fontsize=6, color=PALETTE["PCA-200+LR"])

    # Annotate scGPT slot with "data pending"
    if "scGPT" in FMS and df[df["fm"] == "scGPT"].empty:
        idx = FMS.index("scGPT")
        ax.text(idx, 0.05, "data\npending", ha="center", va="bottom",
                fontsize=5.8, color="#888", style="italic")

    ax.set_xticks(x)
    ax.set_xticklabels(FMS, rotation=0)
    ax.set_ylim(0, 0.95)
    ax.set_ylabel("Mortality AUC (5-fold patient-level CV)")
    ax.set_title("In-hospital mortality on platelet transcriptomes\n(COVID + Sepsis subset, n=42,095 cells / 322 donors)",
                 loc="left", fontsize=7.5, pad=4)

    # Custom legend: solid bar = LR, faded bar = RF
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor="#666", edgecolor="black", linewidth=0.4, label="LR"),
        Patch(facecolor="#666", edgecolor="black", linewidth=0.4, alpha=0.55, label="RF"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", fontsize=6.2,
              handlelength=1.0, handletextpad=0.4, borderaxespad=0.2)

    fig.tight_layout()
    save_pub(fig, "fig_performance")
    print("  ✓ fig_performance.png/pdf/svg")


# ---------------------------------------------------------------------------
# Fig 3a — Cross-disease bidirectional heatmap (severe-excluded, Job 23761950)
# ---------------------------------------------------------------------------
def make_fig3a_cross_disease():
    df = pd.read_csv(DATA / "cross_disease_severe_excluded.csv")
    df["model"] = df["model"].map(canon)
    df = df[~df["model"].isin(DROP)]
    df = df[df["classifier"] == "LogReg"]

    # FMs we have data for in this CSV: STATE, scPRINT-2, UCE.
    # Add scGPT manually from legacy manuscript value (0.352 forward).
    available = [f for f in FMS if f in df["model"].unique() or f == "scGPT"]
    fms_show = [f for f in available if f != "Geneformer"]  # Geneformer missing
    if "Geneformer" in FMS:
        fms_show = ["STATE", "UCE", "scPRINT-2", "scGPT"]  # Geneformer omitted with note

    # Matrix [n_fms × 2 directions]
    forward = []
    reverse = []
    fwd_lo, fwd_hi = [], []
    rev_lo, rev_hi = [], []
    for fm in fms_show:
        if fm == "scGPT":
            forward.append(0.352)
            reverse.append(0.555)
            fwd_lo.append(0.217); fwd_hi.append(0.503)
            rev_lo.append(0.395); rev_hi.append(0.667)
            continue
        f = df[(df["model"] == fm) & (df["direction"] == "COVID->Sepsis")]
        r = df[(df["model"] == fm) & (df["direction"] == "Sepsis->COVID")]
        forward.append(f["auc"].iloc[0] if len(f) else np.nan)
        reverse.append(r["auc"].iloc[0] if len(r) else np.nan)
        fwd_lo.append(f["boot_ci_lo"].iloc[0] if len(f) else np.nan)
        fwd_hi.append(f["boot_ci_hi"].iloc[0] if len(f) else np.nan)
        rev_lo.append(r["boot_ci_lo"].iloc[0] if len(r) else np.nan)
        rev_hi.append(r["boot_ci_hi"].iloc[0] if len(r) else np.nan)

    fig, ax = plt.subplots(figsize=(3.6, 2.5))
    y = np.arange(len(fms_show))
    height = 0.36

    fwd_err = np.array([
        [f - lo if not np.isnan(lo) else 0 for f, lo in zip(forward, fwd_lo)],
        [hi - f if not np.isnan(hi) else 0 for f, hi in zip(forward, fwd_hi)],
    ])
    rev_err = np.array([
        [r - lo if not np.isnan(lo) else 0 for r, lo in zip(reverse, rev_lo)],
        [hi - r if not np.isnan(hi) else 0 for r, hi in zip(reverse, rev_hi)],
    ])

    ax.barh(y - height / 2, forward, height, xerr=fwd_err,
            color=[PALETTE.get(fm, "#888") for fm in fms_show],
            edgecolor="black", linewidth=0.4,
            error_kw=dict(elinewidth=0.6, capsize=1.5, capthick=0.6, ecolor="#444"),
            label="COVID → Sepsis")
    ax.barh(y + height / 2, reverse, height, xerr=rev_err,
            color=[PALETTE.get(fm, "#888") for fm in fms_show],
            edgecolor="black", linewidth=0.4, alpha=0.55,
            error_kw=dict(elinewidth=0.6, capsize=1.5, capthick=0.6, ecolor="#444"),
            label="Sepsis → COVID")

    # Vertical line at chance (0.5)
    ax.axvline(0.5, ls=":", lw=0.7, color="#666", zorder=0)
    ax.text(0.5, len(fms_show) - 0.3, "chance", ha="center", va="bottom",
            fontsize=5.8, color="#666")

    ax.set_yticks(y)
    ax.set_yticklabels(fms_show)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("Cross-disease transfer AUC (200-iter patient bootstrap)")
    ax.set_title("Bidirectional cross-disease mortality transfer\n(severe-excluded label; Job 23761950)",
                 loc="left", fontsize=7.5, pad=4)

    # Annotate STATE as the only bidirectional generalizer (place inside the STATE row,
    # right-aligned, no arrow — avoids overlapping the title)
    state_idx = fms_show.index("STATE")
    ax.text(0.98, state_idx + 0.08, "only bidirectional\ngeneralizer",
            ha="right", va="center", fontsize=5.6, color="#222",
            fontstyle="italic")

    # Annotate scGPT as anti-correlated (place to the right of the bars, no arrow)
    sg_idx = fms_show.index("scGPT")
    ax.text(0.62, sg_idx, "anti-correlated\n(legacy reference)",
            ha="left", va="center", fontsize=5.6, color=PALETTE["scGPT"],
            fontstyle="italic")

    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor="#666", edgecolor="black", linewidth=0.4, label="COVID → Sepsis"),
        Patch(facecolor="#666", edgecolor="black", linewidth=0.4, alpha=0.55, label="Sepsis → COVID"),
    ]
    ax.legend(handles=legend_handles, loc="lower left", fontsize=6.2,
              handlelength=1.0, handletextpad=0.4, borderaxespad=0.4,
              bbox_to_anchor=(0.0, -0.02))

    fig.tight_layout()
    save_pub(fig, "fig_cross_disease_heatmap")
    print("  ✓ fig_cross_disease_heatmap.png/pdf/svg")


# ---------------------------------------------------------------------------
# Fig 3b — Per-FM activation-axis Spearman ρ
# ---------------------------------------------------------------------------
def make_fig3b_activation_axis():
    df = pd.read_csv(DATA / "activation_summary_per_fm.csv")
    df["fm"] = df["model"].map(canon)
    df = df[~df["fm"].isin(DROP)]
    df = df[df["fm"].isin(FMS)]
    # Split the chained set_index/reindex so the comprehension sees the post-set_index df
    df = df.set_index("fm")
    df = df.reindex([f for f in FMS if f in df.index])

    fig, ax = plt.subplots(figsize=(3.4, 2.4))
    x = np.arange(len(df))
    rhos = df["top_dim_spearman"].values
    bar_colors = [PALETTE[fm] for fm in df.index]

    ax.bar(x, rhos, width=0.65,
           color=bar_colors, edgecolor="black", linewidth=0.4)
    ax.axhline(0, color="black", lw=0.7)

    # Reference threshold |ρ| > 0.45
    ax.axhline(0.45, ls=":", lw=0.6, color="#888")
    ax.axhline(-0.45, ls=":", lw=0.6, color="#888")
    ax.text(len(df) - 0.45, 0.46, "|ρ|=0.45", ha="right", va="bottom",
            fontsize=5.5, color="#888")

    # Per-bar value labels
    for xi, (fm, r) in enumerate(zip(df.index, rhos)):
        va = "bottom" if r >= 0 else "top"
        ax.text(xi, r + (0.025 if r >= 0 else -0.025),
                f"{r:+.2f}", ha="center", va=va, fontsize=6)

    ax.set_xticks(x)
    ax.set_xticklabels(df.index, rotation=0, fontsize=6.5)
    ax.set_ylim(-0.7, 0.7)
    ax.set_ylabel("Spearman ρ\n(top dim vs. activation module score)")
    ax.set_title("Platelet activation axis per FM\n(+ : cytoskeleton/adhesion;  − : exhaustion)",
                 loc="left", fontsize=7.5, pad=4)

    fig.tight_layout()
    save_pub(fig, "fig_activation_axis")
    print("  ✓ fig_activation_axis.png/pdf/svg")


# ---------------------------------------------------------------------------
# Fig 4 — LOBO heatmap (FMs × held-out cohorts)
# ---------------------------------------------------------------------------
def make_fig4_lobo():
    df = pd.read_csv(DATA / "lobo.csv")
    df["fm"] = df["fm"].map(canon)
    df = df[~df["fm"].isin(DROP)]
    df = df[df["fm"].isin(FMS)]

    # Pivot to [FM x cohort] matrix
    M = df.pivot_table(index="fm", columns="held_cohort", values="auc")
    fms_order = [f for f in FMS if f in M.index]
    M = M.loc[fms_order]

    # Sort cohorts by mean AUC (descending) so easy cohorts on left, hard on right
    cohort_order = M.mean(axis=0).sort_values(ascending=False).index.tolist()
    M = M[cohort_order]

    fig, ax = plt.subplots(figsize=(4.5, 2.4))
    # Green-Red diverging for AUC around 0.5 chance; clamp 0.3..0.9
    cmap = plt.cm.RdYlGn
    im = ax.imshow(M.values, aspect="auto", cmap=cmap, vmin=0.3, vmax=0.9)

    # Cell text
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M.values[i, j]
            color = "white" if v < 0.45 or v > 0.78 else "black"
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    fontsize=5.5, color=color)

    ax.set_yticks(range(len(M.index)))
    ax.set_yticklabels(M.index)
    ax.set_xticks(range(len(M.columns)))
    ax.set_xticklabels(M.columns, rotation=45, ha="right", fontsize=5.8)
    ax.set_xlabel("Held-out cohort")
    ax.set_title("Leave-one-cohort-out (LOBO) generalization\n(mortality LR AUC; bold = best per cohort)",
                 loc="left", fontsize=7.5, pad=4)

    # Mark best per cohort with a bold border
    for j in range(M.shape[1]):
        best_i = np.nanargmax(M.values[:, j])
        rect = plt.Rectangle((j - 0.5, best_i - 0.5), 1, 1,
                             fill=False, edgecolor="black", lw=1.0)
        ax.add_patch(rect)

    # Colorbar
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("AUC", fontsize=6.5)
    cb.ax.tick_params(labelsize=5.5)
    cb.outline.set_linewidth(0.5)

    fig.tight_layout()
    save_pub(fig, "fig_lobo_heatmap")
    print("  ✓ fig_lobo_heatmap.png/pdf/svg")


# ---------------------------------------------------------------------------
# Fig 5 — Adaptation curve (AUC vs n added sepsis cells)
# ---------------------------------------------------------------------------
def make_fig5_adaptation_curve():
    df = pd.read_csv(DATA / "adaptation_curve_summary.csv")
    df["method"] = df["method"].map(canon)
    # Drop STACK and TF; keep STATE/UCE/Geneformer/scPRINT-2 + PCA-200+LR
    df = df[~df["method"].isin(DROP)]
    df["method"] = df["method"].replace({"PCA-200+LR": "PCA-200+LR"})
    keep = ["STATE", "UCE", "Geneformer", "scPRINT-2"]
    df = df[df["method"].isin(keep + ["PCA-200+LR"])]
    # Drop the "n_sepsis_req = -1" cap-on-train rows that signal asymptote-only entry
    df = df[df["n_sepsis"] >= 0]
    # Sort
    df = df.sort_values(["method", "n_sepsis"])

    fig, ax = plt.subplots(figsize=(3.8, 2.6))

    for method in keep + ["PCA-200+LR"]:
        sub = df[df["method"] == method]
        if sub.empty:
            continue
        color = PALETTE.get(method, "#888")
        ls = "--" if method == "PCA-200+LR" else "-"
        lw = 1.2 if method == "PCA-200+LR" else 1.0
        ax.plot(sub["n_sepsis"], sub["boot_auc_mean"],
                ls=ls, lw=lw, color=color, label=method, zorder=3 if method == "PCA-200+LR" else 2)
        # Shaded CI
        ax.fill_between(sub["n_sepsis"], sub["boot_ci_lo_mean"], sub["boot_ci_hi_mean"],
                        color=color, alpha=0.12, lw=0)

    ax.set_xscale("symlog", linthresh=10)
    ax.set_xlabel("# sepsis cells added to COVID training set (symlog)")
    ax.set_ylabel("Cross-disease AUC on sepsis test (bootstrap mean)")
    ax.set_title("Adaptation curve — PCA-200+LR baseline competitive at all scales",
                 loc="left", fontsize=7.5, pad=4)
    ax.set_ylim(0.55, 1.0)
    ax.axhline(0.5, ls=":", lw=0.6, color="#888", zorder=0)

    ax.legend(loc="lower right", fontsize=6.0, handlelength=1.5,
              ncol=1, columnspacing=0.7)

    fig.tight_layout()
    save_pub(fig, "fig_adaptation_curve")
    print("  ✓ fig_adaptation_curve.png/pdf/svg")


# ---------------------------------------------------------------------------
def main():
    print("Regenerating 5-FM main-text figures for 302006 manuscript...")
    print(f"  Data dir:    {DATA}")
    print(f"  Output dir:  {OUT}")
    print()
    make_fig2_performance()
    make_fig3a_cross_disease()
    make_fig3b_activation_axis()
    make_fig4_lobo()
    make_fig5_adaptation_curve()
    print()
    print("All 5 figures saved as PNG + PDF + SVG.")


if __name__ == "__main__":
    main()
