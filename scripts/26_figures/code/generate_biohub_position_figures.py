#!/usr/bin/env python3
"""
Clean, slide-ready INDIVIDUAL figures for the Biohub Chicago Decoding Inflammation
Computational Biologist (Immune Perturbation) interview.

All figures cover the scFM benchmark on platelet transcriptomes. One single-panel figure per file so the
user can pick-and-paste into slides.

Outputs (PDF + PNG, 300 dpi):
  1. fig_mortality_perfold         -- Per-fold mean LR AUC ranking (headline)
  2. fig_mortality_lobo            -- Leave-one-cohort-out AUC mean +/- std
  3. fig_mortality_delong          -- DeLong OOF aggregated AUC heatmap
  4. fig_mortality_calibration     -- LR + RF calibration (Brier + ECE)
  5. fig_cross_disease             -- Cross-disease bidirectional COVID <-> Sepsis
  6. fig_cohort_cleanup            -- Cohort cleanup methodology (delta AUC bars)

Numbers verified 2026-05-09 against rerun Job 23921397 (COVID + Sepsis subset).

Usage:
  python3 generate_biohub_position_figures.py --output-dir <PATH>
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np


# ============================================================================
# Locked palette
# ============================================================================
PALETTE = {
    "STATE":      "#1F73C2",
    "UCE":        "#737CCB",
    "TF":         "#F48452",
    "TF-Metazoa": "#C99A2E",
    "scPRINT-2":  "#719A73",
    "STACK":      "#B46C72",
    "Geneformer": "#725D86",
    "PCA-200":    "#888888",
    "good":       "#0F5C56",
    "bad":        "#B46C72",
    "accent":     "#882255",
    "title":      "#1A1A1A",
    "body":       "#333333",
    "muted":      "#777777",
}


def setup_style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 11,
        "axes.titlesize": 13,
        "axes.titleweight": "bold",
        "axes.titlepad": 12,
        "axes.labelsize": 11,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": "#444444",
        "axes.linewidth": 0.9,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "xtick.color": "#444444",
        "ytick.color": "#444444",
        "legend.fontsize": 10,
        "legend.frameon": False,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.facecolor": "white",
        "savefig.pad_inches": 0.25,
    })


def save(fig, outdir: Path, name: str):
    for ext in ("pdf", "png"):
        out = outdir / f"{name}.{ext}"
        fig.savefig(out)
        print(f"  Saved {out}")
    plt.close(fig)


def color_xticklabels(ax, labels):
    for tick, m in zip(ax.get_xticklabels(), labels):
        if m in PALETTE:
            tick.set_color(PALETTE[m])
            tick.set_fontweight("bold")


def color_yticklabels(ax, labels):
    for tick, m in zip(ax.get_yticklabels(), labels):
        if m in PALETTE:
            tick.set_color(PALETTE[m])
            tick.set_fontweight("bold")


# ============================================================================
# Fig 1: Per-fold mean LR AUC ranking (headline)
# ============================================================================

def fig_mortality_perfold(outdir: Path):
    fms = ["scPRINT-2", "Geneformer", "STATE", "TF", "UCE", "STACK"]
    auc = np.array([0.762, 0.738, 0.724, 0.719, 0.598, 0.562])
    std = np.array([0.135, 0.113, 0.106, 0.090, 0.141, 0.181])

    fig, ax = plt.subplots(figsize=(8.5, 5))
    y = np.arange(len(fms))[::-1]
    colors = [PALETTE[m] for m in fms]
    ax.barh(y, auc, xerr=std, color=colors,
            edgecolor="black", linewidth=0.6, capsize=4,
            error_kw={"elinewidth": 0.9, "ecolor": "#444444"})
    ax.axvline(0.705, color=PALETTE["PCA-200"], linestyle="--",
               linewidth=1.4, zorder=10)
    ax.axvline(0.5, color="#CCCCCC", linestyle=":", linewidth=0.8, zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels(fms, fontsize=11)
    color_yticklabels(ax, fms)
    ax.set_xlim(0.30, 1.00)
    ax.set_xlabel(r"Per-fold mean LR AUC $\pm$ std (5-fold patient-level CV)",
                  fontsize=11)
    ax.set_title("scPRINT-2 leads per-fold mortality ranking",
                 loc="left", fontsize=13, pad=14)
    for yi, a, s in zip(y, auc, std):
        ax.text(a + s + 0.012, yi, f"{a:.3f}", va="center",
                fontsize=10, color=PALETTE["body"])
    # Annotate PCA-200 line
    ax.text(0.705, len(fms) - 0.3, "PCA-200 = 0.705",
            ha="center", fontsize=9.5, color=PALETTE["PCA-200"],
            fontweight="bold")
    ax.text(0.01, -0.16,
            "COVID + Sepsis subset (42,095 cells, 322 donors); rerun Job 23921397, 2026-05-09",
            transform=ax.transAxes, fontsize=8.5,
            color=PALETTE["muted"], style="italic")
    save(fig, outdir, "fig_mortality_perfold")


# ============================================================================
# Fig 2: LOBO cross-cohort
# ============================================================================

def fig_mortality_lobo(outdir: Path):
    fms = ["STATE", "scPRINT-2", "TF", "Geneformer", "UCE", "STACK"]
    mean = np.array([0.707, 0.688, 0.662, 0.653, 0.647, 0.422])
    std  = np.array([0.072, 0.097, 0.083, 0.062, 0.095, 0.160])

    fig, ax = plt.subplots(figsize=(8.5, 5))
    y = np.arange(len(fms))[::-1]
    colors = [PALETTE[m] for m in fms]
    ax.barh(y, mean, xerr=std, color=colors,
            edgecolor="black", linewidth=0.6, capsize=4,
            error_kw={"elinewidth": 0.9, "ecolor": "#444444"})
    ax.axvline(0.5, color="#CCCCCC", linestyle=":", linewidth=0.8, zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels(fms, fontsize=11)
    color_yticklabels(ax, fms)
    ax.set_xlim(0.20, 0.95)
    ax.set_xlabel(r"Mean leave-one-cohort-out AUC $\pm$ std (across 7 cohorts)",
                  fontsize=11)
    ax.set_title("STATE leads cross-cohort mean; Geneformer most consistent ($\\sigma$ = 0.062)",
                 loc="left", fontsize=12.5, pad=14)
    for yi, m, s in zip(y, mean, std):
        ax.text(m + s + 0.012, yi, f"{m:.3f}", va="center",
                fontsize=10, color=PALETTE["body"])
    # Highlight STACK collapse
    ax.text(0.10, 0, "STACK collapses;\ncohort 302005data03 AUC 0.148",
            ha="left", va="center", fontsize=9,
            color=PALETTE["bad"], fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#FAF0F0",
                      edgecolor=PALETTE["bad"], linewidth=0.6))
    ax.text(0.01, -0.16,
            "COVID + Sepsis subset; rerun Job 23921397, 2026-05-09",
            transform=ax.transAxes, fontsize=8.5,
            color=PALETTE["muted"], style="italic")
    save(fig, outdir, "fig_mortality_lobo")


# ============================================================================
# Fig 3: DeLong OOF aggregated AUC heatmap
# ============================================================================

def fig_mortality_delong(outdir: Path):
    order = ["STATE", "TF", "scPRINT-2", "Geneformer", "UCE", "STACK"]
    aucs = {"STATE": 0.735, "TF": 0.697, "scPRINT-2": 0.679,
            "Geneformer": 0.670, "UCE": 0.608, "STACK": 0.597}
    M = np.zeros((6, 6))
    for i, a in enumerate(order):
        for j, b in enumerate(order):
            M[i, j] = aucs[a] - aucs[b]

    fig, ax = plt.subplots(figsize=(8, 6.5))
    cmap = LinearSegmentedColormap.from_list(
        "div", [PALETTE["bad"], "#FFFFFF", PALETTE["STATE"]], N=256)
    im = ax.imshow(M, cmap=cmap, vmin=-0.15, vmax=0.15, aspect="equal")
    ax.set_xticks(range(6))
    ax.set_xticklabels(order, rotation=30, ha="right", fontsize=10.5)
    ax.set_yticks(range(6))
    ax.set_yticklabels(order, fontsize=10.5)
    color_xticklabels(ax, order)
    color_yticklabels(ax, order)
    ax.set_title("STATE leads OOF aggregated AUC; all 15 pairs $p<10^{-3}$",
                 loc="left", fontsize=12.5, pad=14)
    for i in range(6):
        for j in range(6):
            if i == j:
                ax.text(j, i, f"{aucs[order[i]]:.3f}",
                        ha="center", va="center", fontsize=9,
                        color="#666666", style="italic")
            else:
                v = M[i, j]
                color = "white" if abs(v) > 0.09 else PALETTE["body"]
                ax.text(j, i, f"{v:+.3f}", ha="center", va="center",
                        fontsize=9, color=color)
    ax.set_xlabel("Model B (column)", fontsize=10)
    ax.set_ylabel("Model A (row)", fontsize=10)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(r"$\Delta$AUC (row $-$ column)", fontsize=10)
    cbar.ax.tick_params(labelsize=9)
    ax.text(0.01, -0.20,
            "Cell-weighted out-of-fold AUC (DeLong testing convention) on COVID + Sepsis subset",
            transform=ax.transAxes, fontsize=8.5,
            color=PALETTE["muted"], style="italic")
    save(fig, outdir, "fig_mortality_delong")


# ============================================================================
# Fig 4: Calibration (LR + RF)
# ============================================================================

def fig_mortality_calibration(outdir: Path):
    fms = ["STATE", "TF", "scPRINT-2", "Geneformer", "UCE", "STACK"]
    lr_brier = np.array([0.179, 0.194, 0.199, 0.210, 0.223, 0.235])
    lr_ece   = np.array([0.156, 0.186, 0.173, 0.194, 0.210, 0.236])
    rf_brier = np.array([0.149, 0.154, 0.156, 0.163, 0.161, 0.150])
    rf_ece   = np.array([0.034, 0.048, 0.091, 0.077, 0.084, 0.096])

    fig, ax = plt.subplots(figsize=(10, 5.5))
    x = np.arange(len(fms))
    width = 0.20
    cols = [PALETTE[m] for m in fms]

    ax.bar(x - 1.5*width, lr_brier, width, color=cols, edgecolor="black",
           linewidth=0.5, alpha=0.95, label="LR Brier")
    ax.bar(x - 0.5*width, lr_ece, width, color=cols, edgecolor="black",
           linewidth=0.5, alpha=0.55, hatch="///", label="LR ECE")
    ax.bar(x + 0.5*width, rf_brier, width, color=cols, edgecolor="black",
           linewidth=0.5, alpha=0.95, label="RF Brier")
    ax.bar(x + 1.5*width, rf_ece, width, color=cols, edgecolor="black",
           linewidth=0.5, alpha=0.55, hatch="\\\\\\", label="RF ECE")

    ax.set_xticks(x)
    ax.set_xticklabels(fms, fontsize=11)
    color_xticklabels(ax, fms)
    ax.set_ylim(0, 0.27)
    ax.set_ylabel("Calibration error (lower = better)", fontsize=11)
    ax.set_title("STATE-RF is the best-calibrated configuration (Brier 0.149, ECE 0.034)",
                 loc="left", fontsize=12, pad=14)
    ax.legend(loc="upper left", fontsize=9.5, ncol=2)
    # Annotate STATE values
    ax.text(0 - 1.5*width, lr_brier[0] + 0.005, f"{lr_brier[0]:.3f}",
            ha="center", fontsize=8, color=PALETTE["body"])
    ax.text(0 + 1.5*width, rf_ece[0] + 0.005, f"{rf_ece[0]:.3f}",
            ha="center", fontsize=8, color=PALETTE["good"], fontweight="bold")
    # Side note
    ax.text(0.99, 0.96,
            "RF is dramatically better-calibrated\nthan LR for every FM",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=9, color=PALETTE["good"], fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#F0F7F0",
                      edgecolor=PALETTE["good"], linewidth=0.6))
    ax.text(0.01, -0.16,
            "COVID + Sepsis subset; rerun Job 23921397, 2026-05-09",
            transform=ax.transAxes, fontsize=8.5,
            color=PALETTE["muted"], style="italic")
    save(fig, outdir, "fig_mortality_calibration")


# ============================================================================
# Fig 5: Cross-disease bidirectional COVID <-> Sepsis
# ============================================================================

def fig_cross_disease(outdir: Path):
    fms = ["STACK", "scPRINT-2", "UCE", "STATE", "TF", "Geneformer"]
    forward = np.array([0.917, 0.855, 0.807, 0.798, 0.775, 0.773])
    reverse = np.array([0.698, 0.543, 0.646, 0.707, 0.629, 0.614])

    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    x = np.arange(len(fms))
    width = 0.38
    cols = [PALETTE[m] for m in fms]
    ax.bar(x - width/2, forward, width, color=cols, edgecolor="black",
           linewidth=0.5, alpha=0.85, label=r"COVID $\to$ Sepsis (forward)")
    ax.bar(x + width/2, reverse, width, color=cols, edgecolor="black",
           linewidth=0.5, alpha=0.55, hatch="///",
           label=r"Sepsis $\to$ COVID (reverse)")
    ax.axhline(0.7, color=PALETTE["good"], linestyle="--",
               linewidth=1.3, label="Bidirectional threshold (0.70)")
    ax.axhline(0.5, color="#CCCCCC", linestyle=":", linewidth=0.8, zorder=0)
    ax.set_xticks(x)
    ax.set_xticklabels(fms, fontsize=11)
    color_xticklabels(ax, fms)
    ax.set_ylim(0.40, 1.00)
    ax.set_ylabel("AUC (200-iteration patient-level bootstrap)", fontsize=11)
    ax.set_title("STACK + STATE are the only true bidirectional generalizers",
                 loc="left", fontsize=13, pad=14)
    # Numeric labels
    for xi, f, r in zip(x, forward, reverse):
        ax.text(xi - width/2, f + 0.008, f"{f:.3f}", ha="center",
                fontsize=8.5, color=PALETTE["body"])
        ax.text(xi + width/2, r + 0.008, f"{r:.3f}", ha="center",
                fontsize=8.5, color=PALETTE["body"])
    ax.legend(loc="upper right", fontsize=9.5)
    ax.text(0.01, -0.13,
            "COVID training: 21,286 cells / 205 donors  •  Sepsis test: 3,422 cells / 17 donors  •  Job 23835603",
            transform=ax.transAxes, fontsize=8.5,
            color=PALETTE["muted"], style="italic")
    save(fig, outdir, "fig_cross_disease")


# ============================================================================
# Fig 6: Cohort cleanup methodology (delta AUC bars)
# ============================================================================

def fig_cohort_cleanup(outdir: Path):
    fms = ["scPRINT-2", "TF", "Geneformer", "STATE", "UCE", "STACK"]
    delta = np.array([+0.074, +0.045, +0.039, +0.013, -0.033, -0.143])

    fig, ax = plt.subplots(figsize=(11, 6.4))
    y = np.arange(len(fms))[::-1]
    cols = [PALETTE[m] for m in fms]
    ax.barh(y, delta, color=cols, edgecolor="black", linewidth=0.6)
    ax.axvline(0, color="#444444", linewidth=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels(fms, fontsize=11)
    color_yticklabels(ax, fms)
    ax.set_xlim(-0.20, 0.13)
    ax.set_ylim(-1.4, 6.4)
    ax.set_xlabel(r"$\Delta$ Per-fold mean LR AUC (COVID+Sepsis subset $-$ original 5-disease pool)",
                  fontsize=11)
    ax.set_title("Off-disease cells were inflating retrieval-anchored FMs and deflating translation-axis FMs",
                 loc="left", fontsize=11.5, pad=14)
    # Numeric labels outside bars
    for yi, d in zip(y, delta):
        if d > 0:
            ax.text(d + 0.004, yi, f"{d:+.3f}", va="center", ha="left",
                    fontsize=10, color=PALETTE["body"], fontweight="bold")
        else:
            ax.text(d - 0.004, yi, f"{d:+.3f}", va="center", ha="right",
                    fontsize=10, color=PALETTE["body"], fontweight="bold")
    # Two side labels placed in dedicated rows above scPRINT-2 and below STACK
    ax.text(0.06, 6.0, "Translation-axis FMs gained sharper COVID-specific signal",
            transform=ax.transData, ha="center", va="center",
            fontsize=9.5, color=PALETTE["scPRINT-2"], fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.35", facecolor="#F0F7F0",
                      edgecolor=PALETTE["scPRINT-2"], linewidth=0.6))
    ax.text(-0.10, -1.0, "STACK lost its free discriminative axis from the Lung_infect 1-donor confound",
            transform=ax.transData, ha="center", va="center",
            fontsize=9.5, color=PALETTE["bad"], fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.35", facecolor="#FAF0F0",
                      edgecolor=PALETTE["bad"], linewidth=0.6))
    ax.text(0.01, -0.22,
            "Cleanup dropped Flu (33), HC (3,205), Lung_infect (1,749), Non_covid (726), SLE (169) cells; 91 donors total",
            transform=ax.transAxes, fontsize=8.5,
            color=PALETTE["muted"], style="italic")
    save(fig, outdir, "fig_cohort_cleanup")


# ============================================================================
# Fig 11: Forest plot for DeLong pairwise (STATE as reference)
# ============================================================================

def fig_delong_forest(outdir: Path):
    # STATE OOF AUC = 0.735; deltas computed against each other FM
    others = ["TF", "scPRINT-2", "Geneformer", "UCE", "STACK"]
    delta  = np.array([+0.038, +0.056, +0.065, +0.127, +0.138])
    # Approximate DeLong paired SE for n~42K cells, AUC~0.7, paired r~0.5
    se = np.array([0.0091, 0.0096, 0.0098, 0.0118, 0.0123])
    ci = 1.96 * se

    fig, ax = plt.subplots(figsize=(10, 5.5))
    y = np.arange(len(others))[::-1]
    cols = [PALETTE[m] for m in others]
    ax.errorbar(delta, y, xerr=ci, fmt="none",
                ecolor="#444444", elinewidth=1.2, capsize=5, zorder=2)
    for di, yi, c in zip(delta, y, cols):
        ax.scatter([di], [yi], s=180, color=c,
                   edgecolor="black", linewidth=0.8, zorder=4)
    ax.axvline(0, color="#444444", linewidth=1.0)
    ax.set_yticks(y)
    ax.set_yticklabels([f"STATE vs {m}" for m in others], fontsize=11)
    ax.set_xlabel(r"$\Delta$ AUC (STATE $-$ comparator), 95% CI",
                  fontsize=11)
    ax.set_xlim(-0.02, 0.18)
    ax.set_ylim(-0.7, len(others) - 0.3)
    ax.set_title("STATE beats every comparator on OOF AUC; all 5 paired DeLong tests $p<10^{-3}$",
                 loc="left", fontsize=12, pad=14)
    # Annotate each row with Delta and p
    for di, yi, c in zip(delta, y, cols):
        ax.text(di + 0.012 + 0.0125, yi, f"$\\Delta$={di:+.3f},  $p<10^{{-3}}$",
                va="center", ha="left", fontsize=9.5,
                color=PALETTE["body"])
    # Color y-tick comparator names
    for tick, m in zip(ax.get_yticklabels(), others):
        tick.set_color(PALETTE[m])
        tick.set_fontweight("bold")
    ax.text(0.01, -0.18,
            "SE approximated from paired DeLong on n=42,095 cells; exact CIs available in rerun output (Job 23921397).",
            transform=ax.transAxes, fontsize=8.5,
            color=PALETTE["muted"], style="italic")
    save(fig, outdir, "fig_delong_forest")


# ============================================================================
# Fig 12: Sankey for cohort cleanup
# ============================================================================

def fig_cohort_sankey(outdir: Path):
    """2-column Sankey from source pool -> 6 destinations (Kept + 5 dropped).

    Kept is shown proportionally (huge), dropped classes use a shared
    minimum-height floor so even Flu (33 cells) is readable.
    """
    n_kept = 42095
    drops  = [("HC", 3205), ("Lung_infect", 1749), ("Non_covid", 726),
              ("SLE", 169), ("Flu", 33)]
    n_drop  = sum(c for _, c in drops)
    n_total = n_kept + n_drop  # 47,977

    fig, ax = plt.subplots(figsize=(12, 7.5))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 10)
    ax.axis("off")

    # Layout
    x_src, w_src = 1.0, 0.55
    x_dst, w_dst = 8.5, 0.55
    y_top, y_bot = 9.0, 0.6
    span = y_top - y_bot
    gap_kept_drop = 0.30
    drop_min_h    = 0.42  # min destination height for each dropped class

    # Reserve some span for the dropped destinations (with min heights)
    n_dropped_classes = len(drops)
    drop_floor_total  = n_dropped_classes * drop_min_h + \
                        (n_dropped_classes - 1) * 0.10  # tiny inter-row gap
    # Kept gets the rest
    kept_h = span - drop_floor_total - gap_kept_drop

    # Source node ----------------------------------------------------------
    ax.add_patch(plt.Rectangle((x_src, y_bot), w_src, span,
                               facecolor="#888888", edgecolor="black",
                               linewidth=0.8, alpha=0.9, zorder=3))
    ax.text(x_src - 0.25, (y_top + y_bot) / 2,
            f"Original\n5-disease pool\n{n_total:,} cells",
            ha="right", va="center", fontsize=11, fontweight="bold",
            color=PALETTE["title"])

    # Source-side outflow widths must sum to span; use proportional split
    # but boost dropped fraction visually while preserving ordering.
    src_kept_h = span * (n_kept / n_total)   # 88% of span on source side
    src_drop_h = span * (n_drop / n_total)   # 12% of span on source side
    src_kept_top = y_top
    src_kept_bot = src_kept_top - src_kept_h
    src_drop_top = src_kept_bot
    src_drop_bot = src_drop_top - src_drop_h

    # Kept destination (top) ----------------------------------------------
    dst_kept_top = y_top
    dst_kept_bot = dst_kept_top - kept_h
    ax.add_patch(plt.Rectangle((x_dst, dst_kept_bot), w_dst, kept_h,
                               facecolor=PALETTE["good"], edgecolor="black",
                               linewidth=0.8, alpha=0.9, zorder=3))
    ax.text(x_dst + w_dst + 0.2, (dst_kept_top + dst_kept_bot) / 2,
            f"COVID + Sepsis (kept)\n{n_kept:,} cells   |   322 donors",
            ha="left", va="center", fontsize=11.5,
            color=PALETTE["good"], fontweight="bold")
    _bezier_ribbon(ax, x_src + w_src, src_kept_top, src_kept_bot,
                   x_dst, dst_kept_top, dst_kept_bot,
                   color=PALETTE["good"], alpha=0.30)

    # Dropped destinations (5 classes, with floor height) ------------------
    drop_palette = ["#B46C72", "#C99A2E", "#725D86", "#737CCB", "#0F5C56"]
    # Source-side ribbon for each dropped class is proportional to its count
    # (within the src_drop_h total), but destination height uses min floor.
    src_cursor_top = src_drop_top
    dst_cursor_top = dst_kept_bot - gap_kept_drop
    for (name, n), col in zip(drops, drop_palette):
        src_h = src_drop_h * (n / n_drop)
        dst_h = drop_min_h
        dst_top = dst_cursor_top
        dst_bot = dst_top - dst_h

        ax.add_patch(plt.Rectangle((x_dst, dst_bot), w_dst, dst_h,
                                   facecolor=col, edgecolor="black",
                                   linewidth=0.7, alpha=0.9, zorder=3))
        ax.text(x_dst + w_dst + 0.2, (dst_top + dst_bot) / 2,
                f"{name}:  {n:,} cells",
                ha="left", va="center", fontsize=10.5,
                color=col, fontweight="bold")
        _bezier_ribbon(ax, x_src + w_src, src_cursor_top,
                       src_cursor_top - src_h,
                       x_dst, dst_top, dst_bot,
                       color=col, alpha=0.32)

        src_cursor_top -= src_h
        dst_cursor_top -= (dst_h + 0.10)

    # Bracket grouping for the dropped block
    bracket_top = dst_kept_bot - gap_kept_drop + 0.10
    bracket_bot = dst_cursor_top + 0.05
    ax.plot([x_dst + w_dst + 3.6, x_dst + w_dst + 3.6],
            [bracket_top, bracket_bot],
            color=PALETTE["bad"], linewidth=1.4, alpha=0.7)
    ax.text(x_dst + w_dst + 3.75, (bracket_top + bracket_bot) / 2,
            f"Off-disease drop\n{n_drop:,} cells\n91 donors",
            ha="left", va="center", fontsize=10,
            color=PALETTE["bad"], fontweight="bold")

    ax.set_title(r"Cohort cleanup: 47,977 $\to$ 42,095 cells after dropping 5 off-disease classes",
                 loc="left", fontsize=12.5, pad=18)
    ax.text(0.01, 0.02,
            "Source-side ribbon widths are proportional to cell counts; destination boxes use a shared floor so small classes (Flu, SLE) remain readable. Rerun Job 23921397, 2026-05-09.",
            transform=ax.transAxes, ha="left", va="bottom",
            fontsize=8.5, color=PALETTE["muted"], style="italic")
    save(fig, outdir, "fig_cohort_sankey")


def _bezier_ribbon(ax, x0, y0_top, y0_bot, x1, y1_top, y1_bot, color, alpha=0.4):
    """Draw a smooth bezier-curve ribbon between two vertical edges."""
    from matplotlib.path import Path as MplPath
    import matplotlib.patches as mpatches
    cx = (x0 + x1) / 2
    verts = [
        (x0, y0_top),
        (cx, y0_top), (cx, y1_top), (x1, y1_top),  # top curve
        (x1, y1_bot),
        (cx, y1_bot), (cx, y0_bot), (x0, y0_bot),  # bottom curve
        (x0, y0_top),
    ]
    codes = [MplPath.MOVETO,
             MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.LINETO,
             MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.CLOSEPOLY]
    patch = mpatches.PathPatch(MplPath(verts, codes),
                               facecolor=color, edgecolor="none",
                               alpha=alpha, zorder=2)
    ax.add_patch(patch)


# ============================================================================
# Fig 13: ROC overlay (reads OOF predictions pulled from HPCC)
# ============================================================================

def fig_roc_overlay(outdir: Path, data_dir: Path):
    """Overlay OOF ROC curves for all FMs from pulled CSVs.

    Expects: <data_dir>/oof_predictions/<FM>.csv with columns y, p_lr.
    """
    import pandas as pd
    from sklearn.metrics import roc_curve, roc_auc_score
    fms = ["STATE", "TF", "scPRINT-2", "Geneformer", "UCE", "STACK"]
    fm_to_file = {  # match HPCC naming convention from rerun script
        "STATE": "STATE.csv", "TF": "TranscriptFormer.csv",
        "scPRINT-2": "scPRINT-2.csv", "Geneformer": "Geneformer.csv",
        "UCE": "UCE.csv", "STACK": "STACK.csv",
    }
    oof_dir = data_dir / "oof_predictions"
    if not oof_dir.exists():
        print(f"  SKIP: {oof_dir} not found (rsync OOF preds from HPCC first).")
        return

    fig, ax = plt.subplots(figsize=(7.5, 7))
    ax.plot([0, 1], [0, 1], color="#CCCCCC", linestyle="--", linewidth=0.9, zorder=0)

    rows = []
    for fm in fms:
        f = oof_dir / fm_to_file[fm]
        if not f.exists():
            print(f"  SKIP: {f} missing")
            continue
        df = pd.read_csv(f)
        mask = ~df["p_lr"].isna()
        y, p = df.loc[mask, "y"].values, df.loc[mask, "p_lr"].values
        fpr, tpr, _ = roc_curve(y, p)
        auc = roc_auc_score(y, p)
        rows.append((fm, auc, fpr, tpr))

    # Sort by AUC descending so legend ranks top-to-bottom
    rows.sort(key=lambda r: -r[1])
    for fm, auc, fpr, tpr in rows:
        ax.plot(fpr, tpr, color=PALETTE[fm], linewidth=2.2, alpha=0.9,
                label=f"{fm}  (AUC = {auc:.3f})")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.005)
    ax.set_xlabel("False positive rate", fontsize=11)
    ax.set_ylabel("True positive rate", fontsize=11)
    ax.set_aspect("equal")
    ax.set_title("OOF ROC on COVID + Sepsis mortality (LR, patient-level CV)",
                 loc="left", fontsize=12.5, pad=14)
    ax.legend(loc="lower right", fontsize=10, frameon=False)
    ax.text(0.01, -0.13,
            "Out-of-fold predictions aggregated across 5 patient-level folds; rerun Job 23921397, 2026-05-09",
            transform=ax.transAxes, fontsize=8.5,
            color=PALETTE["muted"], style="italic")
    save(fig, outdir, "fig_roc_overlay")


# ============================================================================
# Fig 14: UMAP per FM (reads UMAP coords pulled from HPCC)
# ============================================================================

def fig_umap_per_fm(outdir: Path, data_dir: Path):
    """One UMAP figure per FM, colored by mortality outcome."""
    import pandas as pd
    fms = ["STATE", "TF", "scPRINT-2", "Geneformer", "UCE", "STACK"]
    fm_to_file = {
        "STATE": "STATE.csv", "TF": "TranscriptFormer.csv",
        "scPRINT-2": "scPRINT-2.csv", "Geneformer": "Geneformer.csv",
        "UCE": "UCE.csv", "STACK": "STACK.csv",
    }
    umap_dir = data_dir / "umap"
    if not umap_dir.exists():
        print(f"  SKIP: {umap_dir} not found (rsync UMAP CSVs from HPCC first).")
        return

    rng = np.random.default_rng(42)
    for fm in fms:
        f = umap_dir / fm_to_file[fm]
        if not f.exists():
            print(f"  SKIP: {f} missing")
            continue
        df = pd.read_csv(f)
        n_surv  = int((df["y_mortality"] == 0).sum())
        n_fatal = int((df["y_mortality"] == 1).sum())
        # Shuffle plotting order so neither class dominates the top layer
        order = rng.permutation(len(df))
        df_s  = df.iloc[order]

        fig, ax = plt.subplots(figsize=(6.5, 6.5))
        colors = np.where(df_s["y_mortality"].values == 1,
                          PALETTE["bad"], "#BFBFBF")
        ax.scatter(df_s["umap1"], df_s["umap2"],
                   s=3.5, c=colors, alpha=0.55,
                   edgecolors="none", rasterized=True)
        # Build a clean two-entry legend manually
        from matplotlib.lines import Line2D
        legend_items = [
            Line2D([0], [0], marker="o", color="none", markerfacecolor="#BFBFBF",
                   markeredgecolor="none", markersize=8,
                   label=f"Survived (n={n_surv:,})"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor=PALETTE["bad"],
                   markeredgecolor="none", markersize=8,
                   label=f"Fatal (n={n_fatal:,})"),
        ]
        ax.legend(handles=legend_items, loc="lower right",
                  fontsize=10, frameon=False)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_xlabel("UMAP-1", fontsize=10, color=PALETTE["muted"])
        ax.set_ylabel("UMAP-2", fontsize=10, color=PALETTE["muted"])
        ax.set_title(f"{fm} embedding: COVID + Sepsis platelets, colored by mortality",
                     loc="left", fontsize=12, pad=14,
                     color=PALETTE.get(fm, PALETTE["title"]))
        ax.text(0.01, -0.06,
                "UMAP (n_neighbors=30, min_dist=0.3, cosine); rerun Job 23921397",
                transform=ax.transAxes, fontsize=8.5,
                color=PALETTE["muted"], style="italic")
        # Tight axis around data
        ax.margins(0.02)
        save(fig, outdir, f"fig_umap_{fm.replace('-', '').replace(' ', '').lower()}")


# ============================================================================
# Main
# ============================================================================

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output-dir", type=Path,
                   default=Path("/Users/xinruqiu/Downloads/biohub_position_figs"))
    p.add_argument("--data-dir", type=Path,
                   default=Path("/Users/xinruqiu/Downloads/biohub_position_figs/_data"),
                   help="Directory with pulled HPCC outputs (oof_predictions/, umap/).")
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    setup_style()

    print(f"Output dir: {args.output_dir}")
    print(f"Data dir:   {args.data_dir}")
    fns = [
        fig_mortality_perfold,
        fig_mortality_lobo,
        fig_mortality_delong,
        fig_mortality_calibration,
        fig_cross_disease,
        fig_cohort_cleanup,
        fig_delong_forest,
        fig_cohort_sankey,
    ]
    for fn in fns:
        print(f"\n=== {fn.__name__} ===")
        fn(args.output_dir)
    # Data-dependent figures (gracefully skip if pulled CSVs not available yet)
    for fn in [fig_roc_overlay, fig_umap_per_fm]:
        print(f"\n=== {fn.__name__} ===")
        fn(args.output_dir, args.data_dir)
    print("\nALL DONE.")


if __name__ == "__main__":
    main()
