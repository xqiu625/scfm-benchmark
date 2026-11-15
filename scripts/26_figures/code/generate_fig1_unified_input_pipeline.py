#!/usr/bin/env python3
"""
Figure 1: Study design and unified-input standardized benchmarking pipeline.

Creates a clean, flat layout-reference schematic for BioRender redraw.
Outputs:
  reports/manuscript/figures/fig1/fig1_unified_input_pipeline.png
  reports/manuscript/figures/fig1/fig1_unified_input_pipeline.svg
"""

from __future__ import annotations

import textwrap
from pathlib import Path
import os

Path("/tmp/matplotlib").mkdir(parents=True, exist_ok=True)
Path("/tmp/fontconfig").mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/fontconfig")
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch, Rectangle


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "reports" / "manuscript" / "figures" / "fig1"
OUT.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica", "sans-serif"],
        "svg.fonttype": "none",
        "font.size": 8.0,
        "axes.linewidth": 0.0,
    }
)

COLORS = {
    "STATE": "#4472C4",
    "UCE": "#70AD47",
    "Geneformer": "#7B5EA7",
    "scPRINT-2": "#ED7D31",
    "scGPT": "#C00000",
    "baseline": "#D9D9D9",
    "gray_fill": "#EDEDED",
    "gray_edge": "#808080",
    "arrow": "#6F6F6F",
    "text": "#1F1F1F",
}


def wrap(text: str, width: int) -> str:
    return "\n".join(
        textwrap.wrap(
            text,
            width=width,
            break_long_words=False,
            break_on_hyphens=False,
        )
    )


def rounded_box(
    ax,
    x: float,
    y: float,
    w: float,
    h: float,
    face: str,
    edge: str,
    lw: float = 1.0,
    radius: float = 0.020,
):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0.004,rounding_size={radius}",
        linewidth=lw,
        edgecolor=edge,
        facecolor=face,
    )
    ax.add_patch(patch)
    return patch


def box_text(
    ax,
    x: float,
    y: float,
    w: float,
    h: float,
    title: str,
    body: str | None = None,
    title_size: float = 8.7,
    body_size: float = 7.1,
    title_color: str = COLORS["text"],
    body_color: str = COLORS["text"],
    title_weight: str = "bold",
    body_weight: str = "normal",
    line_spacing: float = 1.18,
):
    if body:
        ax.text(
            x + w / 2,
            y + h * 0.69,
            title,
            ha="center",
            va="center",
            fontsize=title_size,
            fontweight=title_weight,
            color=title_color,
            linespacing=line_spacing,
        )
        ax.text(
            x + w / 2,
            y + h * 0.34,
            body,
            ha="center",
            va="center",
            fontsize=body_size,
            fontweight=body_weight,
            color=body_color,
            linespacing=line_spacing,
        )
    else:
        ax.text(
            x + w / 2,
            y + h / 2,
            title,
            ha="center",
            va="center",
            fontsize=title_size,
            fontweight=title_weight,
            color=title_color,
            linespacing=line_spacing,
        )


def arrow(ax, start, end, lw: float = 1.5, mutation_scale: float = 12.0):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            color=COLORS["arrow"],
            linewidth=lw,
            mutation_scale=mutation_scale,
            shrinkA=3,
            shrinkB=3,
            connectionstyle="arc3,rad=0",
        )
    )


def double_arrow(ax, start, end, lw: float = 1.5, mutation_scale: float = 11.0):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="<->",
            color=COLORS["arrow"],
            linewidth=lw,
            mutation_scale=mutation_scale,
            shrinkA=0,
            shrinkB=0,
        )
    )


def tier_label(ax, text: str, y: float):
    ax.text(
        1.0,
        y,
        text,
        ha="left",
        va="center",
        fontsize=8.4,
        fontweight="bold",
        color="#4A4A4A",
        linespacing=1.08,
    )


def draw_platelet(ax, cx: float, cy: float):
    ax.add_patch(
        Ellipse(
            (cx, cy),
            width=6.0,
            height=2.6,
            angle=-10,
            facecolor="#F4F4F4",
            edgecolor=COLORS["gray_edge"],
            linewidth=1.2,
        )
    )
    ax.add_patch(
        Ellipse(
            (cx - 0.6, cy + 0.15),
            width=2.2,
            height=0.8,
            angle=-10,
            facecolor="white",
            edgecolor="#B0B0B0",
            linewidth=0.8,
        )
    )


def draw_donor_folds(ax, x: float, y: float):
    fold_w = 2.6
    gap = 0.45
    for idx in range(5):
        fx = x + idx * (fold_w + gap)
        ax.add_patch(Rectangle((fx, y - 0.85), fold_w, 1.7, facecolor="white", edgecolor=COLORS["gray_edge"], linewidth=0.8))
        for j in range(3):
            px = fx + 0.55 + j * 0.65
            ax.add_patch(Circle((px, y + 0.18), 0.16, facecolor="#595959", edgecolor="none"))
            ax.add_patch(Rectangle((px - 0.10, y - 0.45), 0.20, 0.45, facecolor="#595959", edgecolor="none"))
        ax.text(fx + fold_w / 2, y - 1.45, str(idx + 1), ha="center", va="center", fontsize=7.0, color="#4A4A4A")


def draw():
    fig, ax = plt.subplots(figsize=(16, 10), dpi=300)
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    left = 16.0
    right = 98.0
    width = right - left

    tier_label(ax, "1. INPUT", 84.5)
    tier_label(ax, "2. UNIFIED-INPUT\nSTANDARDIZATION +\nMODELS", 66.8)
    tier_label(ax, "3. EVALUATION ARMS", 36.0)

    # Tier 1
    draw_platelet(ax, left + 5.2, 86.6)
    ax.text(
        left + 5.2,
        79.4,
        "platelets: anucleate, absent from FM pretraining atlases",
        ha="center",
        va="center",
        fontsize=7.2,
        color="#404040",
    )
    t1 = (left + 14.0, 81.0, width - 14.0, 9.6)
    rounded_box(ax, *t1, face="white", edge=COLORS["gray_edge"], lw=1.0)
    ax.text(
        t1[0] + t1[2] / 2,
        86.6,
        "Human platelet scRNA-seq compendium | 47,977 cells · 413 donors · 11 cohorts",
        ha="center",
        va="center",
        fontsize=10.4,
        fontweight="bold",
        color=COLORS["text"],
    )
    ax.text(
        t1[0] + t1[2] / 2,
        83.8,
        "Primary task: COVID-19 + Sepsis subset | 42,095 cells · 322 donors · in-hospital mortality (fatal vs. survived)",
        ha="center",
        va="center",
        fontsize=8.4,
        color="#404040",
    )

    # Tier 2
    banner = (left, 71.2, width, 5.8)
    rounded_box(ax, *banner, face=COLORS["gray_fill"], edge=COLORS["gray_edge"], lw=1.0)
    box_text(
        ax,
        *banner,
        title="Unified-input standardization: identical raw AnnData input; each FM applies its own native preprocessing",
        title_size=9.8,
        title_weight="bold",
    )

    t1_center_x = t1[0] + t1[2] / 2
    arrow(ax, (t1_center_x, 81.0), (t1_center_x, 77.0), lw=2.6, mutation_scale=15)

    models = [
        ("STATE", "state-space, perturbation-aware"),
        ("UCE", "cross-species universal embedding"),
        ("Geneformer", "rank-value MLM (+ LoRA arm)"),
        ("scPRINT-2", "multi-task, ontology-conditioned"),
        ("scGPT", "generative MLM"),
    ]
    gap = 1.0
    baseline_w = 13.5
    model_w = (width - 5 * gap - baseline_w) / 5
    model_y = 60.0
    model_h = 8.6
    track_centers = []
    for idx, (name, desc) in enumerate(models):
        x = left + idx * (model_w + gap)
        track_centers.append((x + model_w / 2, model_y + model_h / 2))
        rounded_box(ax, x, model_y, model_w, model_h, face=COLORS[name], edge="#4D4D4D", lw=0.8)
        box_text(
            ax,
            x,
            model_y,
            model_w,
            model_h,
            title=name,
            body=wrap(desc, 28),
            title_size=8.9,
            body_size=7.0,
            title_color="white",
            body_color="white",
        )
        arrow(ax, (left + width / 2, banner[1]), (x + model_w / 2, model_y + model_h), lw=1.5)

    base_x = left + 5 * (model_w + gap)
    baseline = (base_x, model_y, baseline_w, model_h)
    track_centers.append((base_x + baseline_w / 2, model_y + model_h / 2))
    rounded_box(ax, *baseline, face=COLORS["baseline"], edge=COLORS["gray_edge"], lw=0.9)
    box_text(
        ax,
        *baseline,
        title="PCA-200 + LR baseline",
        body="raw counts, same protocol",
        title_size=8.2,
        body_size=7.0,
        title_weight="bold",
    )
    arrow(ax, (left + width / 2, banner[1]), (base_x + baseline_w / 2, model_y + model_h), lw=1.5)

    # Tier 3
    eval_banner = (left, 48.5, width, 5.8)
    rounded_box(ax, *eval_banner, face=COLORS["gray_fill"], edge=COLORS["gray_edge"], lw=1.0)
    draw_donor_folds(ax, left + 3.0, 51.6)
    ax.text(
        left + 22.0,
        51.4,
        "patient-level StratifiedGroupKFold, K = 5",
        ha="left",
        va="center",
        fontsize=9.6,
        fontweight="bold",
        color=COLORS["text"],
    )
    for idx, center in enumerate(track_centers):
        target_x = left + 6.0 + idx * ((width - 12.0) / 5)
        arrow(ax, (center[0], model_y), (target_x, eval_banner[1] + eval_banner[3]), lw=1.5)

    grid_top = 44.5
    box_w = (width - 2 * 1.4) / 3
    box_h = 11.2
    row_gap = 1.45
    eval_items = [
        (
            "Within-disease mortality",
            "per-fold & OOF AUC · DeLong + donor bootstrap",
            False,
        ),
        ("Calibration", "Brier · ECE · LR vs. RF", False),
        ("LOBO cross-cohort", "7 held-out cohorts · mean AUC + sigma", False),
        (
            "Cross-disease transfer",
            "COVID<->Sepsis bidirectional · bootstrap CIs",
            True,
        ),
        (
            "Adaptive",
            "sepsis labels added: 0 -> all · vs. PCA baseline",
            False,
        ),
        (
            "Sample efficiency + dimension dropout",
            "5-100% training · 0-70% dim zeroing",
            False,
        ),
        (
            "Embedding ablations",
            "top-PC removal · normalization · random dims",
            False,
        ),
        (
            "LoRA fine-tune (Geneformer)",
            "r=16, alpha=32 · patient-level CV",
            False,
        ),
        (
            "Interpretability",
            "activation axis · SHAP · GSEA · Jaccard",
            False,
        ),
    ]
    for idx, (title, body, has_double_arrow) in enumerate(eval_items):
        row, col = divmod(idx, 3)
        x = left + col * (box_w + 1.4)
        y = grid_top - row * (box_h + row_gap) - box_h
        rounded_box(ax, x, y, box_w, box_h, face="white", edge=COLORS["gray_edge"], lw=1.0)
        title_y = y + box_h * 0.66
        body_y = y + box_h * 0.34
        box_text(
            ax,
            x,
            y,
            box_w,
            box_h,
            title=wrap(title, 34),
            body=wrap(body, 45),
            title_size=8.2,
            body_size=7.0,
            title_color=COLORS["text"],
            body_color="#303030",
        )
        if has_double_arrow:
            double_arrow(ax, (x + box_w * 0.41, title_y - 2.2), (x + box_w * 0.59, title_y - 2.2), lw=1.4)

    png = OUT / "fig1_unified_input_pipeline.png"
    svg = OUT / "fig1_unified_input_pipeline.svg"
    fig.savefig(png, dpi=300, facecolor="white")
    fig.savefig(svg, facecolor="white")
    plt.close(fig)
    print(png)
    print(svg)


if __name__ == "__main__":
    draw()
