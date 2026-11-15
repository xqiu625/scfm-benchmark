#!/usr/bin/env python3
"""Workflow schematic for the frozen-embedding LoRA-style classifier.

This figure intentionally distinguishes the January 2026 LoRA-classifier arm
from true transformer PEFT.  It visualizes the actual frozen-feature workflow:

    AnnData -> pre-computed FM embeddings -> frozen X -> low-rank classifier

Only the low-rank adapter matrices and classifier output are trainable.
Foundation-model transformer weights and pre-computed embeddings are not
updated.
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
LIGHT_GREY = "#F4F5F7"
LINE = "#303030"
GENEFORMER = "#725D86"
LORA = "#B46C72"
FROZEN = "#4A4A4A"
TRAIN = "#1F73C2"
DATA = "#0F5C56"
WARNING = "#C99A2E"


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


def box(
    ax,
    xy,
    wh,
    title,
    subtitle="",
    fc="white",
    ec=LINE,
    lw=1.0,
    title_size=10,
    subtitle_size=8,
    title_color=TEXT,
    subtitle_color=GREY,
    weight="bold",
):
    x, y = xy
    w, h = wh
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.015,rounding_size=0.012",
        linewidth=lw,
        edgecolor=ec,
        facecolor=fc,
    )
    ax.add_patch(patch)
    if subtitle:
        ax.text(
            x + w / 2,
            y + h * 0.62,
            title,
            ha="center",
            va="center",
            fontsize=title_size,
            fontweight=weight,
            color=title_color,
            linespacing=1.05,
        )
        ax.text(
            x + w / 2,
            y + h * 0.30,
            subtitle,
            ha="center",
            va="center",
            fontsize=subtitle_size,
            color=subtitle_color,
            linespacing=1.15,
        )
    else:
        ax.text(
            x + w / 2,
            y + h / 2,
            title,
            ha="center",
            va="center",
            fontsize=title_size,
            fontweight=weight,
            color=title_color,
            linespacing=1.05,
        )
    return patch


def arrow(ax, start, end, color=LINE, lw=1.2, style="-|>"):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle=style,
            mutation_scale=13,
            linewidth=lw,
            color=color,
            shrinkA=3,
            shrinkB=3,
        )
    )


def step_label(ax, x, y, label, color=GREY):
    ax.text(x, y, label, fontsize=8.5, color=color, fontweight="bold", ha="left", va="center")


def draw_matrix(ax, x, y, w, h, rows=7, cols=9, color=GENEFORMER):
    ax.add_patch(Rectangle((x, y), w, h, facecolor="white", edgecolor=LINE, linewidth=0.9))
    for i in range(rows):
        for j in range(cols):
            alpha = 0.15 + 0.65 * ((i * 3 + j * 5) % 11) / 10
            ax.add_patch(
                Rectangle(
                    (x + j * w / cols, y + i * h / rows),
                    w / cols,
                    h / rows,
                    facecolor=color,
                    edgecolor="white",
                    linewidth=0.25,
                    alpha=alpha,
                )
            )
    ax.text(x + w / 2, y - 0.032, "X_geneformer  (n cells x 1152)", ha="center", fontsize=8, color=GREY)


def draw_lora_equation(ax):
    y0 = 0.43
    ax.text(0.575, y0 + 0.210, "scPEFT-inspired\nLoRA head", ha="center", va="center", fontsize=10, fontweight="bold", color=LORA, linespacing=0.95)
    ax.text(0.575, y0 + 0.165, "low-rank adapter branch\nadded to frozen linear path", ha="center", va="center", fontsize=7.2, color=GREY, linespacing=1.05)

    # Main frozen path Wx
    box(ax, (0.430, y0 + 0.052), (0.09, 0.075), "W", "frozen", fc="#FFFFFF", ec=FROZEN, title_color=FROZEN)
    box(ax, (0.585, y0 + 0.052), (0.09, 0.075), "x", "embedding", fc="#FFFFFF", ec=GENEFORMER, title_color=GENEFORMER)
    arrow(ax, (0.536, y0 + 0.0895), (0.569, y0 + 0.0895), color=FROZEN)

    # Adapter path BAx
    box(ax, (0.410, y0 - 0.070), (0.060, 0.068), "LoRA A", "train", fc="#FFF5F7", ec=LORA, title_color=LORA, title_size=8.2)
    box(ax, (0.5175, y0 - 0.070), (0.060, 0.068), "LoRA B", "train", fc="#FFF5F7", ec=LORA, title_color=LORA, title_size=8.2)
    box(ax, (0.625, y0 - 0.070), (0.055, 0.068), "alpha/r", "scale", fc="#FFF5F7", ec=LORA, title_color=LORA, title_size=9)
    arrow(ax, (0.486, y0 - 0.036), (0.5015, y0 - 0.036), color=LORA)
    arrow(ax, (0.5935, y0 - 0.036), (0.609, y0 - 0.036), color=LORA)

    ax.text(
        0.55,
        y0 - 0.145,
        "logits = W x  +  (B A x) * alpha/r",
        ha="center",
        va="center",
        fontsize=11,
        fontweight="bold",
        color=TEXT,
    )
    ax.text(
        0.55,
        y0 - 0.205,
        "Backprop updates A, B, and output head only; no FM weights are touched.",
        ha="center",
        va="center",
        fontsize=7.8,
        color=GREY,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path("reports/manuscript/figures/fig_lora_workflow"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    setup_style()
    fig, ax = plt.subplots(figsize=(13.33, 7.5))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    fig.suptitle(
        "Frozen-embedding scPEFT-inspired LoRA classifier head",
        fontsize=13.5,
        fontweight="bold",
        y=0.965,
    )
    ax.text(
        0.5,
        0.905,
        "LoRA math from scPEFT is applied after frozen FM embeddings: W x + (B A x) * alpha/r.",
        ha="center",
        fontsize=9.5,
        color=GREY,
    )

    # Left: data and frozen embedding extraction.
    step_label(ax, 0.055, 0.82, "A. Feature extraction")
    box(
        ax,
        (0.055, 0.635),
        (0.18, 0.145),
        "Platelet AnnData\n(COVID + Sepsis)",
        "cells + labels + Sample_ID",
        fc="#EEF8F6",
        ec=DATA,
        title_color=DATA,
    )
    box(
        ax,
        (0.30, 0.635),
        (0.20, 0.145),
        "Geneformer\nembedding file",
        ".obsm['X_geneformer']\npre-computed once",
        fc="#F4F0F7",
        ec=GENEFORMER,
        title_color=GENEFORMER,
    )
    arrow(ax, (0.251, 0.7075), (0.284, 0.7075), color=GENEFORMER)
    draw_matrix(ax, 0.080, 0.45, 0.29, 0.135, color=GENEFORMER)
    arrow(ax, (0.40, 0.618), (0.355, 0.596), color=GENEFORMER)
    box(
        ax,
        (0.080, 0.315),
        (0.29, 0.075),
        "Frozen design matrix",
        "no gradient into embedding generator",
        fc="#FAFAFA",
        ec=FROZEN,
        title_color=FROZEN,
    )
    arrow(ax, (0.33, 0.448), (0.33, 0.407), color=FROZEN)

    # Middle: LoRA classifier math.
    step_label(ax, 0.45, 0.82, "B. scPEFT-inspired LoRA adapter head")
    draw_lora_equation(ax)
    arrow(ax, (0.385, 0.340), (0.409, 0.372), color=GENEFORMER)

    # Right: evaluation / CV.
    step_label(ax, 0.72, 0.82, "C. Patient-held-out evaluation")
    box(
        ax,
        (0.72, 0.655),
        (0.22, 0.13),
        "5-fold\nStratifiedGroupKFold",
        "groups = Sample_ID\nprevents donor leakage",
        fc="#F2F6FA",
        ec=TRAIN,
        title_color=TRAIN,
        title_size=9.2,
        subtitle_size=7.6,
    )
    box(
        ax,
        (0.72, 0.495),
        (0.22, 0.105),
        "Train",
        "A, B + output head\ninside each fold",
        fc="#EEF5FF",
        ec=TRAIN,
        title_color=TRAIN,
    )
    box(
        ax,
        (0.72, 0.345),
        (0.22, 0.105),
        "Evaluate",
        "held-out patients\nAUC by task",
        fc="#FFFFFF",
        ec=LINE,
        title_color=TEXT,
    )
    arrow(ax, (0.83, 0.639), (0.83, 0.616), color=TRAIN)
    arrow(ax, (0.83, 0.479), (0.83, 0.466), color=TRAIN)

    # Frozen vs trainable legend.
    ax.text(0.075, 0.225, "Parameter status", fontsize=10, fontweight="bold", color=TEXT, ha="left")
    legend_items = [
        (FROZEN, "Frozen: pre-computed embedding X and base linear W"),
        (LORA, "Trainable: low-rank A/B adapters"),
        (TRAIN, "Trainable: final classifier output"),
    ]
    for i, (color, text) in enumerate(legend_items):
        y = 0.185 - i * 0.045
        ax.add_patch(Rectangle((0.078, y - 0.012), 0.018, 0.018, facecolor=color, edgecolor=LINE, linewidth=0.4))
        ax.text(0.105, y - 0.003, text, fontsize=8.6, color=TEXT, ha="left", va="center")

    # Explicit distinction callout.
    box(
        ax,
        (0.51, 0.065),
        (0.40, 0.105),
        "scPEFT connection",
        "Same low-rank A/B adapter formula;\napplied to the classifier head, not transformer layers.",
        fc="#FFF8E8",
        ec=WARNING,
        title_color=WARNING,
        subtitle_color=TEXT,
        title_size=9.2,
        subtitle_size=8,
    )

    for ext in ("pdf", "png", "svg"):
        out = args.output_dir / f"fig_lora_frozen_embedding_workflow.{ext}"
        fig.savefig(out, bbox_inches="tight")
        print(f"Saved {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
