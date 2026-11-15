#!/usr/bin/env python3
"""
Indexed-figure regenerator for the 302006 manuscript.

Produces FigureN_<slug>.{png,pdf,svg} files in reports/manuscript/figures/ so that
\\includegraphics paths line up with the figure numbering used in the
manuscript prose.

Main figures (5-FM primary set: STATE, UCE, Geneformer, scPRINT-2, scGPT):
  Figure1_overview               schematic, copied from fig1/fig1_overview
  Figure2_performance            per-FM mortality LR + RF AUC bars
  Figure3_cross_disease          combined: (a) bidirectional transfer
                                          (b) activation-axis Spearman rho
  Figure4_lobo                   leave-one-cohort-out heatmap
  Figure5_adaptation             cross-disease adaptation curve

Supplementary:
  FigureS1_lora_frozen_workflow         (copied from fig_lora_workflow/)
  FigureS2_deepmlp_workflow             (copied from fig_deepmlp_workflow/)
  FigureS3_geneformer_lora_workflow     (copied from fig_geneformer_true_lora_workflow/)
  FigureS4_per_cohort_oof               per-cohort 5-fold OOF AUC heatmap
  FigureS5_delong_significance          DeLong pairwise dAUC + significance

STACK and TranscriptFormer are dropped from every figure per the 2026-05-26
strategic restructuring. scGPT is fully data-traced: its Option C canonical
embedding landed 2026-06-08 and all of its figure values come from the
cleaner-cohort rerun (the earlier legacy-value placeholder is resolved).
"""

from __future__ import annotations

import shutil
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

# ---------------------------------------------------------------------------
# Nature-figure style (mirrors the existing regenerator scripts)
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

# Resolve repo root relative to this script so paths survive the repo relocation
# (was hardcoded to the old ~/Desktop/302006 path). scripts/26_figures/code/<this> → parents[3].
ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "scripts/26_figures/data"
OUT = ROOT / "reports/manuscript/figures"
OUT.mkdir(parents=True, exist_ok=True)  # figures/ is gitignored; ensure it exists for regeneration

FMS = ["STATE", "UCE", "Geneformer", "scPRINT-2", "scGPT"]
DROP = {"STACK", "STACK ortholog", "TranscriptFormer", "TF", "TF-Metazoa"}

PALETTE = {
    "STATE":      "#3B6BB0",
    "UCE":        "#5C9D7F",
    "Geneformer": "#7C5C9E",
    "scPRINT-2":  "#D4944B",
    "scGPT":      "#A84A4A",
    "PCA-200+LR": "#444444",
}

# Cohort IDs are presented as source-study surnames in figures (2026-08-25 rename;
# mapping per docs/cohort-sources.md). Bernardes and Combes are COVID-19
# cohorts in the 302005 series; only Qiu is a sepsis cohort there. Jiang
# (302004data18) is also Sepsis-labeled but has zero fatal cells, so it never
# appears in AUC-based figures (verified 2026-08-25; record archived at
# archive/superseded_2026-09-21/docs_legacy/cohort-annotation-verification.md).
COHORT_LABELS = {
    "302004data01": "Wilk",
    "302004data02": "Lee",
    "302004data04": "Schulte-Schrepping",
    "302004data05": "Arunachalam",
    "302004data15": "Ren",
    "302004data16": "Liu",
    "302004data17": "Stephenson",
    "302004data18": "Jiang",
    "302005data01": "Bernardes",
    "302005data02": "Combes",
    "302005data03": "Qiu",
    "302005data04": "Mistry",
}


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


def save_pub(fig, stem: str):
    fig.savefig(OUT / f"{stem}.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Figure 1 — Study overview (copy existing schematic)
# ---------------------------------------------------------------------------
def make_figure1_overview():
    src_dir = OUT / "fig1"
    src_stem = "fig1_overview"
    dst_stem = "Figure1_overview"
    n_copied = 0
    for ext in ("png", "pdf", "svg"):
        src = src_dir / f"{src_stem}.{ext}"
        if src.exists():
            shutil.copy2(src, OUT / f"{dst_stem}.{ext}")
            n_copied += 1
    if n_copied == 0:
        print("  !! Figure1_overview source files missing under fig1/; run generate_fig1_overview_5fm.py first")
        return
    print(f"  ✓ Figure1_overview (copied {n_copied} formats from fig1/)")


# ---------------------------------------------------------------------------
# Figure 2 — Per-FM mortality AUC bars
# ---------------------------------------------------------------------------
def make_figure2_performance():
    df = pd.read_csv(DATA / "mortality_lr_rf.csv")
    df["fm"] = df["fm"].map(canon)
    df = df[~df["fm"].isin(DROP)]
    df = df[df["fm"].isin(FMS)]

    pca = pd.read_csv(DATA / "pca_baseline.csv")
    pca_lr200 = pca.loc[pca["method"] == "PCA-200+LR", "auc_mean"].iloc[0]

    classifiers = ["LR", "RF"]
    width = 0.36

    # Order bars by per-fold mean LR AUC (descending), matching tab:performance
    def _lr_mean(fm):
        row = df[(df["fm"] == fm) & (df["classifier"] == "LR")]
        return row["auc_mean"].iloc[0] if len(row) else np.nan
    fms_sorted = sorted(FMS, key=lambda f: -_lr_mean(f))
    x = np.arange(len(fms_sorted))

    fig, ax = plt.subplots(figsize=(3.6, 2.4))

    for i, clf in enumerate(classifiers):
        means, stds = [], []
        for fm in fms_sorted:
            row = df[(df["fm"] == fm) & (df["classifier"] == clf)]
            if len(row) == 0:
                means.append(np.nan)
                stds.append(0)
            else:
                means.append(row["auc_mean"].iloc[0])
                stds.append(row["auc_std"].iloc[0])
        offset = (i - 0.5) * width
        ax.bar(
            x + offset,
            means,
            width,
            yerr=stds,
            color=[PALETTE.get(fm, "#888") for fm in fms_sorted],
            edgecolor="black",
            linewidth=0.4,
            error_kw=dict(elinewidth=0.6, capsize=1.5, capthick=0.6, ecolor="#444"),
            alpha=0.55 if clf == "RF" else 1.0,
            label=clf,
        )

    # PCA-200+LR horizontal reference; label carried in the legend instead of
    # on-plot text (the previous left-side text overlapped the scPRINT-2 bar).
    ax.axhline(pca_lr200, ls="--", lw=0.8, color=PALETTE["PCA-200+LR"], zorder=0)

    if "scGPT" in fms_sorted and df[df["fm"] == "scGPT"].empty:
        idx = fms_sorted.index("scGPT")
        ax.text(idx, 0.05, "data\npending", ha="center", va="bottom",
                fontsize=5.8, color="#888", style="italic")

    ax.set_xticks(x)
    ax.set_xticklabels(fms_sorted, rotation=0)
    ax.set_ylim(0, 0.95)
    ax.set_ylabel("Mortality AUC (5-fold patient-level CV)")
    ax.set_title("In-hospital mortality on platelet transcriptomes\n(COVID + Sepsis subset, n=42,095 cells / 322 donors)",
                 loc="left", fontsize=7.5, pad=4)

    legend_handles = [
        Patch(facecolor="#666", edgecolor="black", linewidth=0.4, label="LR"),
        Patch(facecolor="#666", edgecolor="black", linewidth=0.4, alpha=0.55, label="RF"),
        mpl.lines.Line2D([0], [0], ls="--", lw=0.8, color=PALETTE["PCA-200+LR"],
                         label=f"PCA-200+LR = {pca_lr200:.3f}"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", fontsize=6.2,
              handlelength=1.0, handletextpad=0.4, borderaxespad=0.2)

    fig.tight_layout()
    save_pub(fig, "Figure2_performance")
    print("  ✓ Figure2_performance.{png,pdf,svg}")


# ---------------------------------------------------------------------------
# Figure 3 — Combined cross-disease + activation axis (2-panel composite)
# ---------------------------------------------------------------------------
def _draw_panel_a_cross_disease(ax):
    df = pd.read_csv(DATA / "cross_disease_severe_excluded.csv")
    df["model"] = df["model"].map(canon)
    df = df[~df["model"].isin(DROP)]
    df = df[df["classifier"] == "LogReg"]

    fms_show = ["STATE", "UCE", "scPRINT-2", "Geneformer", "scGPT"]
    forward, reverse = [], []
    fwd_lo, fwd_hi, rev_lo, rev_hi = [], [], [], []
    for fm in fms_show:
        # All shown FMs (incl. scGPT) are read from the snapshot. scGPT was previously
        # hardcoded (0.352/0.555); since 2026-06-08 the snapshot carries its C002 row
        # (build_cross_disease_snapshot.py), so the value is now fully data-traced.
        f = df[(df["model"] == fm) & (df["direction"] == "COVID->Sepsis")]
        r = df[(df["model"] == fm) & (df["direction"] == "Sepsis->COVID")]
        forward.append(f["auc"].iloc[0] if len(f) else np.nan)
        reverse.append(r["auc"].iloc[0] if len(r) else np.nan)
        fwd_lo.append(f["boot_ci_lo"].iloc[0] if len(f) else np.nan)
        fwd_hi.append(f["boot_ci_hi"].iloc[0] if len(f) else np.nan)
        rev_lo.append(r["boot_ci_lo"].iloc[0] if len(r) else np.nan)
        rev_hi.append(r["boot_ci_hi"].iloc[0] if len(r) else np.nan)

    y = np.arange(len(fms_show))
    height = 0.36

    def _err(vals, lo, hi):
        return np.array([
            [v - l if not np.isnan(l) else 0 for v, l in zip(vals, lo)],
            [h - v if not np.isnan(h) else 0 for v, h in zip(vals, hi)],
        ])
    fwd_err = _err(forward, fwd_lo, fwd_hi)
    rev_err = _err(reverse, rev_lo, rev_hi)

    ax.barh(y - height / 2, forward, height, xerr=fwd_err,
            color=[PALETTE.get(fm, "#888") for fm in fms_show],
            edgecolor="black", linewidth=0.4,
            error_kw=dict(elinewidth=0.6, capsize=1.5, capthick=0.6, ecolor="#444"))
    ax.barh(y + height / 2, reverse, height, xerr=rev_err,
            color=[PALETTE.get(fm, "#888") for fm in fms_show],
            edgecolor="black", linewidth=0.4, alpha=0.55,
            error_kw=dict(elinewidth=0.6, capsize=1.5, capthick=0.6, ecolor="#444"))

    ax.axvline(0.5, ls=":", lw=0.7, color="#666", zorder=0)

    ax.set_yticks(y)
    ax.set_yticklabels(fms_show)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.42)
    ax.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xlabel("Cross-disease transfer AUC")
    ax.set_title("(a) Bidirectional cross-disease mortality transfer",
                 loc="left", fontsize=7.5, pad=4)

    # Annotations placed cleanly outside bars (x in [1.05, 1.32])
    state_idx = fms_show.index("STATE")
    ax.text(1.06, state_idx, "strongest\nbidirectional",
            ha="left", va="center", fontsize=5.4, color="#222",
            fontstyle="italic")
    sg_idx = fms_show.index("scGPT")
    ax.text(1.06, sg_idx, "anti-correlated",
            ha="left", va="center", fontsize=5.4, color=PALETTE["scGPT"],
            fontstyle="italic")

    # Legend inside the panel, center right (clear of bars and annotations)
    legend_handles = [
        Patch(facecolor="#666", edgecolor="black", linewidth=0.4, label="COVID $\\rightarrow$ Sepsis"),
        Patch(facecolor="#666", edgecolor="black", linewidth=0.4, alpha=0.55, label="Sepsis $\\rightarrow$ COVID"),
    ]
    ax.legend(handles=legend_handles, fontsize=5.8,
              handlelength=1.0, handletextpad=0.4, borderaxespad=0.2,
              ncol=1, loc="center right")


def _draw_panel_b_activation(ax):
    df = pd.read_csv(DATA / "activation_summary_per_fm.csv")
    df["fm"] = df["model"].map(canon)
    df = df[~df["fm"].isin(DROP)]
    df = df[df["fm"].isin(FMS)]
    df = df.set_index("fm")
    df = df.reindex([f for f in FMS if f in df.index])

    x = np.arange(len(df))
    rhos = df["top_dim_spearman"].values
    bar_colors = [PALETTE[fm] for fm in df.index]

    ax.bar(x, rhos, width=0.65,
           color=bar_colors, edgecolor="black", linewidth=0.4)
    ax.axhline(0, color="black", lw=0.7)
    ax.axhline(0.45, ls=":", lw=0.6, color="#888")
    ax.axhline(-0.45, ls=":", lw=0.6, color="#888")
    ax.text(len(df) - 0.45, 0.46, "|$\\rho$|=0.45", ha="right", va="bottom",
            fontsize=5.5, color="#888")

    for xi, (_, r) in enumerate(zip(df.index, rhos)):
        va = "bottom" if r >= 0 else "top"
        ax.text(xi, r + (0.025 if r >= 0 else -0.025),
                f"{r:+.2f}", ha="center", va=va, fontsize=6)

    ax.set_xticks(x)
    ax.set_xticklabels(df.index, rotation=20, ha="right", fontsize=6.5)
    ax.set_ylim(-0.7, 0.7)
    ax.set_ylabel("Spearman $\\rho$\n(top dim vs. activation module score)")
    ax.set_title("(b) Platelet activation axis per FM",
                 loc="left", fontsize=7.5, pad=4)


def make_figure3_cross_disease():
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.7), gridspec_kw=dict(wspace=0.35))
    _draw_panel_a_cross_disease(axes[0])
    _draw_panel_b_activation(axes[1])
    fig.tight_layout()
    save_pub(fig, "Figure3_cross_disease")
    print("  ✓ Figure3_cross_disease.{png,pdf,svg}  (combined 2-panel)")


# ---------------------------------------------------------------------------
# Figure 4 — LOBO heatmap
# ---------------------------------------------------------------------------
def make_figure4_lobo():
    df = pd.read_csv(DATA / "lobo.csv")
    df["fm"] = df["fm"].map(canon)
    df = df[~df["fm"].isin(DROP)]
    df = df[df["fm"].isin(FMS)]

    M = df.pivot_table(index="fm", columns="held_cohort", values="auc")
    fms_order = [f for f in FMS if f in M.index]
    M = M.loc[fms_order]
    cohort_order = M.mean(axis=0).sort_values(ascending=False).index.tolist()
    M = M[cohort_order]

    fig, ax = plt.subplots(figsize=(4.5, 2.4))
    im = ax.imshow(M.values, aspect="auto", cmap=plt.cm.RdYlGn, vmin=0.3, vmax=0.9)

    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M.values[i, j]
            color = "white" if v < 0.45 or v > 0.78 else "black"
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    fontsize=5.5, color=color)

    ax.set_yticks(range(len(M.index)))
    ax.set_yticklabels(M.index)
    ax.set_xticks(range(len(M.columns)))
    ax.set_xticklabels([COHORT_LABELS.get(c, c) for c in M.columns], rotation=45, ha="right", fontsize=5.8)
    ax.set_xlabel("Held-out cohort")
    ax.set_title("Leave-one-cohort-out generalization\n(mortality LR AUC; bold border = best per cohort)",
                 loc="left", fontsize=7.5, pad=4)

    for j in range(M.shape[1]):
        best_i = np.nanargmax(M.values[:, j])
        ax.add_patch(plt.Rectangle((j - 0.5, best_i - 0.5), 1, 1,
                                   fill=False, edgecolor="black", lw=1.0))

    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("AUC", fontsize=6.5)
    cb.ax.tick_params(labelsize=5.5)
    cb.outline.set_linewidth(0.5)

    fig.tight_layout()
    save_pub(fig, "Figure4_lobo")
    print("  ✓ Figure4_lobo.{png,pdf,svg}")


# ---------------------------------------------------------------------------
# Figure 5 — Adaptation curve
# ---------------------------------------------------------------------------
def make_figure5_adaptation():
    df = pd.read_csv(DATA / "adaptation_curve_summary.csv")
    df["method"] = df["method"].map(canon)
    df = df[~df["method"].isin(DROP)]
    keep = ["STATE", "UCE", "Geneformer", "scPRINT-2", "scGPT"]
    df = df[df["method"].isin(keep + ["PCA-200+LR"])].copy()
    df["is_all"] = df["n_sepsis_req"] == -1
    df = df.sort_values(["method", "n_sepsis"])

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.6, 2.7),
                                   gridspec_kw={"width_ratios": [1.15, 1.0], "wspace": 0.32})

    # ---------------- Panel A: adaptation curves ----------------
    for method in keep + ["PCA-200+LR"]:
        sub = df[(df["method"] == method) & (~df["is_all"])].sort_values("n_sepsis")
        if sub.empty:
            continue
        color = PALETTE.get(method, "#888")
        ls = "--" if method == "PCA-200+LR" else "-"
        lw = 1.2 if method == "PCA-200+LR" else 1.0
        axA.plot(sub["n_sepsis"], sub["boot_auc_mean"],
                 ls=ls, lw=lw, color=color, label=method,
                 zorder=3 if method == "PCA-200+LR" else 2)
        axA.fill_between(sub["n_sepsis"], sub["boot_ci_lo_mean"], sub["boot_ci_hi_mean"],
                         color=color, alpha=0.12, lw=0)
        asym = df[(df["method"] == method) & (df["is_all"])]
        if not asym.empty:
            axA.axhline(asym["point_auc_mean"].iloc[0], color=color, ls=":", lw=0.7,
                        alpha=0.6, zorder=1)

    axA.set_xscale("symlog", linthresh=10)
    # Left-clip the meaningless negative region (n<0 cells); n=0 sits at the
    # left edge with a hair of linear-region padding, no bogus -10^1 tick.
    axA.set_xlim(-1, 1500)
    axA.set_xlabel("# sepsis cells added to COVID training pool (symlog)")
    axA.set_ylabel("Cross-disease AUC on sepsis test\n(bootstrap mean)")
    axA.set_title("A. Adaptation curve",
                  loc="left", fontsize=7.5, pad=4)
    axA.set_ylim(0.55, 1.0)
    axA.axhline(0.5, ls=":", lw=0.6, color="#888", zorder=0)
    axA.legend(loc="lower right", fontsize=5.6, handlelength=1.5)

    # ---------------- Panel B: three quantities vs PCA reference ----------------
    metrics = [("Y-int", 0), ("n=100", 100), ("asymp.", -1)]
    pca = df[df["method"] == "PCA-200+LR"]
    pca_refs = [pca.loc[pca["n_sepsis_req"] == req, "point_auc_mean"].iloc[0]
                for _, req in metrics]

    width = 0.26
    x = np.arange(len(keep))
    for i, (label, req) in enumerate(metrics):
        vals, errs = [], []
        for fm in keep:
            row = df[(df["method"] == fm) & (df["n_sepsis_req"] == req)]
            vals.append(row["point_auc_mean"].iloc[0])
            std = row["point_auc_std"].iloc[0]
            errs.append(0.0 if np.isnan(std) else std)
        offset = (i - 1) * width
        axB.bar(x + offset, vals, width, yerr=errs,
                color=[PALETTE[fm] for fm in keep], alpha=[0.45, 0.7, 1.0][i],
                edgecolor="black", linewidth=0.4,
                error_kw=dict(elinewidth=0.5, capsize=1.2, capthick=0.5, ecolor="#444"),
                label=label)
    for y, ls in zip(pca_refs, ["--", "-.", ":"]):
        axB.axhline(y, color="#444444", ls=ls, lw=0.8, alpha=0.8)
    _lab_kw = dict(fontsize=5.2, color="#444",
                   bbox=dict(facecolor="white", alpha=0.75, edgecolor="none", pad=0.4))
    axB.text(-0.45, pca_refs[0], f"PCA Y-int {pca_refs[0]:.3f}",
             ha="left", va="bottom", **_lab_kw)
    axB.text(len(keep) - 0.45, pca_refs[1], f"PCA n=100 {pca_refs[1]:.3f}",
             ha="right", va="top", **_lab_kw)
    axB.text(len(keep) - 0.45, pca_refs[2], f"PCA asymp. {pca_refs[2]:.3f}",
             ha="right", va="top", **_lab_kw)
    axB.set_xticks(x)
    axB.set_xticklabels(keep, rotation=20, ha="right", fontsize=6)
    axB.set_ylim(0.55, 1.0)
    axB.set_ylabel("AUC (point mean)")
    axB.set_title("B. Three quantities vs. PCA-200+LR reference",
                  loc="left", fontsize=7.5, pad=4)
    axB.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), fontsize=5.4, ncol=3,
               columnspacing=0.8)

    fig.tight_layout()
    save_pub(fig, "Figure5_adaptation")
    print("  ✓ Figure5_adaptation.{png,pdf,svg}  (2-panel)")


# ---------------------------------------------------------------------------
# Figure 6 — Calibration reliability diagram (LR vs RF panels, 5 FMs)
# ---------------------------------------------------------------------------
def make_figure6_calibration():
    # Per-cell out-of-fold predicted probabilities (headline severe-included
    # protocol; verified to reproduce tab:calibration Brier/ECE exactly).
    OOF = ROOT / "scripts/30_review_response_stats/data"
    bins = np.linspace(0, 1, 11)

    def _reliability(p, y):
        idx = np.digitize(p, bins[1:-1])
        xs, ys, ns, ece = [], [], [], 0.0
        for b in range(10):
            m = idx == b
            if m.sum() > 0:
                xs.append(p[m].mean())
                ys.append(y[m].mean())
                ns.append(m.sum())
                ece += m.mean() * abs(p[m].mean() - y[m].mean())
        return np.array(xs), np.array(ys), np.array(ns), ece

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(7.2, 3.1), sharex=True, sharey=True,
                                   gridspec_kw={"wspace": 0.18})
    for ax, col, title in [(axL, "p_lr", "A. Logistic regression"),
                           (axR, "p_rf", "B. Random forest")]:
        ax.plot([0, 1], [0, 1], ls=":", lw=0.7, color="#666", zorder=0,
                label="perfect calibration")
        for fm in FMS:
            df = pd.read_csv(OOF / f"oof_{fm}.csv")
            xs, ys, ns, ece = _reliability(df[col].values, df["y"].values)
            ax.plot(xs, ys, lw=0.7, color=PALETTE[fm], zorder=1, alpha=0.6)
            ax.scatter(xs, ys, s=3 + 30 * ns / ns.max(), color=PALETTE[fm],
                       edgecolor="black", linewidth=0.3, zorder=2,
                       label=f"{fm} ({ece:.3f})")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel("Mean predicted probability (10 equal-width bins)")
        ax.set_title(title, loc="left", fontsize=7.5, pad=4)
        # FM entries first, diagonal ("perfect calibration") last
        handles, labels_ = ax.get_legend_handles_labels()
        handles = handles[1:] + handles[:1]
        labels_ = labels_[1:] + labels_[:1]
        ax.legend(handles, labels_, loc="upper left", fontsize=5.8,
                  title="FM (ECE)", title_fontsize=5.8, borderaxespad=0.2)
    axL.set_ylabel("Observed mortality rate")

    fig.tight_layout()
    save_pub(fig, "Figure6_calibration")
    print("  ✓ Figure6_calibration.{png,pdf,svg}  (2-panel)")


# ---------------------------------------------------------------------------
# Supplementary figures — schematic copies (S1, S2, S3) + data plots (S4, S5)
# ---------------------------------------------------------------------------
def _copy_schematic(src_dir_name: str, src_stem: str, dst_stem: str):
    src_dir = OUT / src_dir_name
    n_copied = 0
    for ext in ("png", "pdf", "svg"):
        src = src_dir / f"{src_stem}.{ext}"
        if src.exists():
            shutil.copy2(src, OUT / f"{dst_stem}.{ext}")
            n_copied += 1
    if n_copied == 0:
        print(f"  !! {dst_stem}: source files not found under {src_dir_name}/{src_stem}.*")
    else:
        print(f"  ✓ {dst_stem} (copied {n_copied} formats from {src_dir_name}/)")


def make_supp_schematics():
    _copy_schematic("fig_lora_workflow",
                    "fig_lora_frozen_embedding_workflow",
                    "FigureS1_lora_frozen_workflow")
    _copy_schematic("fig_deepmlp_workflow",
                    "fig_deepmlp_frozen_embedding_workflow",
                    "FigureS2_deepmlp_workflow")
    _copy_schematic("fig_geneformer_true_lora_workflow",
                    "fig_geneformer_true_lora_workflow",
                    "FigureS3_geneformer_lora_workflow")


def make_figureS4_per_cohort_oof():
    df = pd.read_csv(DATA / "per_cohort_oof.csv")
    df["fm"] = df["fm"].map(canon)
    df = df[~df["fm"].isin(DROP)]
    df = df[df["fm"].isin(FMS)]

    M = df.pivot_table(index="fm", columns="cohort", values="auc")
    fms_order = [f for f in FMS if f in M.index]
    M = M.loc[fms_order]
    cohort_order = M.mean(axis=0).sort_values(ascending=False).index.tolist()
    M = M[cohort_order]

    fig, ax = plt.subplots(figsize=(4.6, 2.4))
    im = ax.imshow(M.values, aspect="auto", cmap=plt.cm.RdYlGn, vmin=0.3, vmax=0.9)

    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M.values[i, j]
            color = "white" if (np.isnan(v) or v < 0.45 or v > 0.78) else "black"
            txt = f"{v:.2f}" if not np.isnan(v) else "--"
            ax.text(j, i, txt, ha="center", va="center", fontsize=5.5, color=color)

    ax.set_yticks(range(len(M.index)))
    ax.set_yticklabels(M.index)
    ax.set_xticks(range(len(M.columns)))
    ax.set_xticklabels([COHORT_LABELS.get(c, c) for c in M.columns], rotation=45, ha="right", fontsize=5.8)
    ax.set_xlabel("Cohort (source study)")
    ax.set_title("Per-cohort standard 5-fold OOF mortality AUC\n(within-cohort, not LOBO)",
                 loc="left", fontsize=7.5, pad=4)

    for j in range(M.shape[1]):
        col = M.values[:, j]
        if np.isnan(col).all():
            continue
        best_i = int(np.nanargmax(col))
        ax.add_patch(plt.Rectangle((j - 0.5, best_i - 0.5), 1, 1,
                                   fill=False, edgecolor="black", lw=1.0))

    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("AUC", fontsize=6.5)
    cb.ax.tick_params(labelsize=5.5)
    cb.outline.set_linewidth(0.5)

    fig.tight_layout()
    save_pub(fig, "FigureS4_per_cohort_oof")
    print("  ✓ FigureS4_per_cohort_oof.{png,pdf,svg}")


def make_figureS5_delong():
    df = pd.read_csv(DATA / "delong_pairwise.csv")
    df["model_a"] = df["model_a"].map(canon)
    df["model_b"] = df["model_b"].map(canon)
    df = df[~df["model_a"].isin(DROP) & ~df["model_b"].isin(DROP)]
    df = df[df["model_a"].isin(FMS) & df["model_b"].isin(FMS)]

    models_present = sorted(set(df["model_a"]).union(df["model_b"]),
                            key=lambda x: FMS.index(x) if x in FMS else 99)

    delta_matrix = pd.DataFrame(np.nan, index=models_present, columns=models_present)
    p_matrix = pd.DataFrame(np.nan, index=models_present, columns=models_present)
    for _, row in df.iterrows():
        a, b = row["model_a"], row["model_b"]
        d, p = row["delta_auc"], row["p_value"]
        delta_matrix.loc[a, b] = d
        delta_matrix.loc[b, a] = -d
        p_matrix.loc[a, b] = p
        p_matrix.loc[b, a] = p

    fig, ax = plt.subplots(figsize=(4.0, 3.4))
    im = ax.imshow(delta_matrix.values, cmap=plt.cm.RdBu_r, vmin=-0.15, vmax=0.15)

    for i, _ in enumerate(models_present):
        for j, _ in enumerate(models_present):
            if i == j:
                ax.text(j, i, "—", ha="center", va="center", fontsize=6, color="black")
                continue
            d = delta_matrix.values[i, j]
            p = p_matrix.values[i, j]
            if np.isnan(d):
                continue
            star = ""
            if not np.isnan(p):
                if p < 1e-4:
                    star = "***"
                elif p < 1e-3:
                    star = "**"
                elif p < 0.05:
                    star = "*"
            txt = f"{d:+.3f}\n{star}" if star else f"{d:+.3f}"
            color = "white" if abs(d) > 0.08 else "black"
            ax.text(j, i, txt, ha="center", va="center", fontsize=5.5, color=color)

    ax.set_xticks(range(len(models_present)))
    ax.set_yticks(range(len(models_present)))
    ax.set_xticklabels(models_present, rotation=45, ha="right", fontsize=6.5)
    ax.set_yticklabels(models_present, fontsize=6.5)

    # Bold cell border = comparison survives the donor-clustered bootstrap with
    # Benjamini-Hochberg correction (tab:delong): only STATE's lead over each of
    # the other four FMs survives; all other pairs are ties under donor clustering.
    state_i = models_present.index("STATE")
    for j in range(len(models_present)):
        if j == state_i:
            continue
        for ii, jj in [(state_i, j), (j, state_i)]:
            ax.add_patch(plt.Rectangle((jj - 0.5, ii - 0.5), 1, 1,
                                       fill=False, edgecolor="black", lw=1.2))
    ax.set_xlabel("Model $B$")
    ax.set_ylabel("Model $A$")
    ax.set_title("DeLong pairwise $\\Delta$AUC = $\\mathrm{AUC}_A - \\mathrm{AUC}_B$\n(* $p<0.05$,  ** $p<10^{-3}$,  *** $p<10^{-4}$)",
                 loc="left", fontsize=7.5, pad=4)

    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label("$\\Delta$AUC", fontsize=6.5)
    cb.ax.tick_params(labelsize=5.5)
    cb.outline.set_linewidth(0.5)

    fig.tight_layout()
    save_pub(fig, "FigureS5_delong_significance")
    print("  ✓ FigureS5_delong_significance.{png,pdf,svg}")


# ---------------------------------------------------------------------------
def main():
    print(f"Output dir: {OUT}")
    print()
    print("Main figures:")
    make_figure1_overview()
    make_figure2_performance()
    make_figure3_cross_disease()
    make_figure4_lobo()
    make_figure5_adaptation()
    print()
    print("Supplementary figures:")
    make_supp_schematics()
    make_figureS4_per_cohort_oof()
    make_figureS5_delong()
    print()
    print("All indexed figures regenerated.")


if __name__ == "__main__":
    main()
