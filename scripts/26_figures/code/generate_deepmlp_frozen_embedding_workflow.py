#!/usr/bin/env python3
"""Workflow schematic for the frozen-embedding DeepMLP classifier head.

This figure illustrates the classifier-architecture comparison used in the
Option C mortality benchmark:

    pre-computed FM embeddings -> frozen design matrix -> LR or DeepMLP head

Only the downstream classifier head is trained. Foundation model weights and
pre-computed embeddings remain frozen.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle


TEXT = "#202020"
GREY = "#6B6B6B"
LINE = "#303030"
DATA = "#0F5C56"
EMB = "#725D86"
MLP = "#1F73C2"
DROP = "#B46C72"
BN = "#719A73"
LR = "#4A4A4A"
NOTE = "#C99A2E"


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "savefig.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def box(ax, xy, wh, title, subtitle="", fc="white", ec=LINE, title_color=TEXT,
        subtitle_color=GREY, title_size=10, subtitle_size=8, lw=1.0):
    x, y = xy
    w, h = wh
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.015,rounding_size=0.012",
        linewidth=lw, edgecolor=ec, facecolor=fc,
    )
    ax.add_patch(patch)
    if subtitle:
        ax.text(x + w / 2, y + h * 0.66, title, ha="center", va="center",
                fontsize=title_size, fontweight="bold", color=title_color,
                linespacing=1.0)
        ax.text(x + w / 2, y + h * 0.28, subtitle, ha="center", va="center",
                fontsize=subtitle_size, color=subtitle_color, linespacing=1.12)
    else:
        ax.text(x + w / 2, y + h / 2, title, ha="center", va="center",
                fontsize=title_size, fontweight="bold", color=title_color,
                linespacing=1.0)
    return patch


def arrow(ax, start, end, color=LINE, lw=1.25):
    ax.add_patch(
        FancyArrowPatch(
            start, end, arrowstyle="-|>", mutation_scale=13,
            linewidth=lw, color=color, shrinkA=3, shrinkB=3,
        )
    )


def step_label(ax, x, y, text):
    ax.text(x, y, text, fontsize=10, fontweight="bold", color=GREY, ha="left", va="center")


def draw_matrix(ax, x, y, w, h, rows=7, cols=9, color=EMB):
    ax.add_patch(Rectangle((x, y), w, h, facecolor="white", edgecolor=LINE, linewidth=0.9))
    for i in range(rows):
        for j in range(cols):
            alpha = 0.15 + 0.65 * ((i * 4 + j * 7) % 11) / 10
            ax.add_patch(
                Rectangle(
                    (x + j * w / cols, y + i * h / rows),
                    w / cols, h / rows,
                    facecolor=color, edgecolor="white", linewidth=0.25, alpha=alpha,
                )
            )
    ax.text(x + w * 0.60, y - 0.018, "X_fm  (n cells x embedding dim)", ha="center", fontsize=8, color=GREY)


def layer(ax, x, y, w, h, title, color, subtitle=""):
    return box(
        ax, (x, y), (w, h), title, subtitle,
        fc="#FFFFFF", ec=color, title_color=color,
        title_size=8.6, subtitle_size=7.2, lw=1.2,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path("reports/manuscript/figures/fig_deepmlp_workflow"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    setup_style()
    fig, ax = plt.subplots(figsize=(13.33, 7.5))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    fig.suptitle(
        "Frozen-embedding DeepMLP head: testing whether nonlinear classifiers add signal",
        fontsize=13.5,
        fontweight="bold",
        y=0.965,
    )
    ax.text(
        0.5, 0.905,
        "Foundation-model embeddings are fixed; only the downstream classifier head is trained inside each patient-held-out fold.",
        ha="center", fontsize=9.3, color=GREY,
    )

    # A. frozen feature extraction
    step_label(ax, 0.055, 0.82, "A. Feature extraction")
    box(
        ax, (0.055, 0.66), (0.19, 0.12),
        "Unified input\nAnnData",
        "cells + mortality labels\n+ Sample_ID",
        fc="#EEF8F6", ec=DATA, title_color=DATA,
    )
    box(
        ax, (0.31, 0.66), (0.20, 0.12),
        "Pre-computed\nFM embedding",
        ".obsm['X_state']\n.obsm['X_geneformer']\n...",
        fc="#F4F0F7", ec=EMB, title_color=EMB,
        title_size=9.2, subtitle_size=7.2,
    )
    arrow(ax, (0.245, 0.72), (0.31, 0.72), EMB)
    draw_matrix(ax, 0.095, 0.455, 0.33, 0.14, color=EMB)
    arrow(ax, (0.41, 0.66), (0.36, 0.60), EMB)
    box(
        ax, (0.095, 0.335), (0.33, 0.075),
        "Frozen design matrix",
        "no gradient into FM or embedding file",
        fc="#FAFAFA", ec=LR, title_color=LR,
    )
    arrow(ax, (0.12, 0.455), (0.12, 0.412), LR)

    # B. two heads compared
    step_label(ax, 0.49, 0.82, "B. Classifier heads compared on the same frozen X")
    box(
        ax, (0.54, 0.655), (0.16, 0.105),
        "LR baseline",
        "StandardScaler\n+ L2 logistic regression",
        fc="#F7F7F7", ec=LR, title_color=LR,
    )
    box(
        ax, (0.78, 0.655), (0.16, 0.105),
        "DeepMLP head",
        "BatchNorm + ReLU\n+ Dropout blocks",
        fc="#EEF5FF", ec=MLP, title_color=MLP,
    )
    arrow(ax, (0.425, 0.357), (0.54, 0.705), EMB)
    arrow(ax, (0.425, 0.357), (0.78, 0.705), EMB)

    # DeepMLP architecture detail
    ax.text(0.70, 0.585, "DeepMLP architecture", fontsize=10.5, fontweight="bold", color=MLP, ha="center")
    xs = [0.485, 0.595, 0.705, 0.815]
    ys = 0.475
    layer(ax, xs[0], ys, 0.085, 0.075, "Linear", MLP, "512")
    layer(ax, xs[1], ys, 0.085, 0.075, "BatchNorm", BN)
    layer(ax, xs[2], ys, 0.085, 0.075, "ReLU", MLP)
    layer(ax, xs[3], ys, 0.085, 0.075, "Dropout", DROP, "p=0.3")
    for i in range(3):
        arrow(ax, (xs[i] + 0.085, ys + 0.037), (xs[i + 1], ys + 0.037), GREY, lw=1.0)

    xs2 = [0.540, 0.650, 0.760]
    ys2 = 0.345
    layer(ax, xs2[0], ys2, 0.085, 0.075, "Linear", MLP, "256")
    layer(ax, xs2[1], ys2, 0.085, 0.075, "BatchNorm", BN)
    layer(ax, xs2[2], ys2, 0.085, 0.075, "Dropout", DROP, "p=0.3")
    for i in range(2):
        arrow(ax, (xs2[i] + 0.085, ys2 + 0.037), (xs2[i + 1], ys2 + 0.037), GREY, lw=1.0)
    arrow(ax, (0.86, 0.472), (0.82, 0.423), GREY, lw=1.0)

    box(
        ax, (0.61, 0.225), (0.19, 0.075),
        "Linear -> 2 logits",
        "fatal vs survived",
        fc="#FFFFFF", ec=MLP, title_color=MLP,
        title_size=9, subtitle_size=7.4,
    )
    arrow(ax, (0.80, 0.343), (0.72, 0.302), GREY, lw=1.0)

    # C. CV / conclusion
    step_label(ax, 0.055, 0.170, "C. Patient-held-out training and result")
    box(
        ax, (0.055, 0.035), (0.25, 0.095),
        "5-fold StratifiedGroupKFold",
        "groups = Sample_ID\nscaler fit inside each fold",
        fc="#F2F6FA", ec=MLP, title_color=MLP,
    )
    arrow(ax, (0.305, 0.083), (0.38, 0.083), MLP)
    box(
        ax, (0.38, 0.035), (0.25, 0.095),
        "Train only classifier head",
        "FM embeddings remain frozen\nno donor leakage",
        fc="#FFFFFF", ec=LINE, title_color=TEXT,
    )
    arrow(ax, (0.63, 0.083), (0.70, 0.083), MLP)
    box(
        ax, (0.70, 0.035), (0.25, 0.095),
        "Mortality result",
        "DeepMLP does not beat LR\nfor any FM (seeded rerun)",
        fc="#FFF8E8", ec=NOTE, title_color=NOTE,
    )

    # status legend
    ax.text(0.055, 0.295, "Parameter status", fontsize=9.8, fontweight="bold", color=TEXT)
    for i, (color, text) in enumerate([
        (LR, "Frozen: FM weights and pre-computed embedding matrix X"),
        (MLP, "Trainable: LR or DeepMLP classifier head"),
        (DROP, "Regularization: BatchNorm + Dropout control MLP variance"),
    ]):
        y = 0.270 - i * 0.035
        ax.add_patch(Rectangle((0.058, y - 0.012), 0.016, 0.016, facecolor=color, edgecolor=LINE, linewidth=0.35))
        ax.text(0.083, y - 0.004, text, fontsize=8.1, color=TEXT, ha="left", va="center")

    for ext in ("pdf", "png", "svg"):
        out = args.output_dir / f"fig_deepmlp_frozen_embedding_workflow.{ext}"
        fig.savefig(out, bbox_inches="tight")
        print(f"Saved {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
