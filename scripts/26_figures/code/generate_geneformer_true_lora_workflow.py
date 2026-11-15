#!/usr/bin/env python3
"""Workflow schematic for true Geneformer scPEFT-style LoRA fine-tuning.

This figure is intentionally distinct from the frozen-embedding LoRA-style
classifier figure.  Here raw expression is re-tokenized into Geneformer input
IDs, a BertForSequenceClassification backbone is loaded, LoRA adapters are
inserted inside transformer Linear modules, and only LoRA + classifier-head
parameters are trainable.
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
TOKEN = "#725D86"
MODEL = "#4A4A4A"
LORA = "#B46C72"
TRAIN = "#1F73C2"
NOTE = "#C99A2E"
GREEN = "#719A73"


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
        ax.text(x + w / 2, y + h * 0.65, title, ha="center", va="center",
                fontsize=title_size, fontweight="bold", color=title_color,
                linespacing=1.0)
        ax.text(x + w / 2, y + h * 0.26, subtitle, ha="center", va="center",
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


def token_strip(ax, x, y, w, h):
    ax.add_patch(Rectangle((x, y), w, h, facecolor="white", edgecolor=TOKEN, linewidth=1.0))
    labels = ["CLS", "ENSG", "ENSG", "ENSG", "...", "SEP", "PAD"]
    widths = [0.10, 0.14, 0.14, 0.14, 0.12, 0.10, 0.12]
    xpos = x + 0.01
    for i, (lab, frac) in enumerate(zip(labels, widths)):
        ww = w * frac
        color = TOKEN if i not in (4, 6) else "#D7D0E0"
        ax.add_patch(Rectangle((xpos, y + 0.01), ww, h - 0.02, facecolor=color, edgecolor="white", linewidth=0.5, alpha=0.88))
        ax.text(xpos + ww / 2, y + h / 2, lab, ha="center", va="center", fontsize=7.3, color="white" if i not in (4, 6) else TEXT, fontweight="bold")
        xpos += ww + 0.004


def transformer_block(ax, x, y, w, h, label, lora=True):
    box(ax, (x, y), (w, h), "", fc="#F7F7F7", ec=MODEL,
        title_color=MODEL, title_size=8.5, subtitle_size=7.1)
    ax.text(x + w / 2, y + h * 0.76, label, ha="center", va="center",
            fontsize=8.5, fontweight="bold", color=MODEL)
    if lora:
        # Small red LoRA tags placed on target Linear modules.
        tags = [("Q", 0.05), ("K", 0.26), ("V", 0.47), ("Dense", 0.68)]
        for txt, frac in tags:
            tx = x + w * frac
            ty = y + h * 0.16
            ax.add_patch(FancyBboxPatch(
                (tx, ty), w * 0.14, h * 0.30,
                boxstyle="round,pad=0.005,rounding_size=0.006",
                facecolor="#FFF5F7", edgecolor=LORA, linewidth=0.9,
            ))
            ax.text(tx + w * 0.070, ty + h * 0.15, txt, ha="center", va="center",
                    fontsize=6.5, color=LORA, fontweight="bold")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path("reports/manuscript/figures/fig_geneformer_true_lora_workflow"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    setup_style()
    fig, ax = plt.subplots(figsize=(13.33, 7.5))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    fig.suptitle(
        "True Geneformer scPEFT-style LoRA fine-tuning",
        fontsize=14,
        fontweight="bold",
        y=0.965,
    )
    ax.text(
        0.5, 0.905,
        "Unlike frozen-embedding heads, this path re-tokenizes expression and trains LoRA adapters inside the Geneformer transformer.",
        ha="center", fontsize=9.3, color=GREY,
    )

    # A. data -> tokenization
    step_label(ax, 0.055, 0.82, "A. Re-tokenize raw expression for Geneformer")
    box(
        ax, (0.055, 0.67), (0.19, 0.11),
        "Geneformer-ready\nAnnData",
        "var_names = Ensembl IDs\nlabels + Sample_ID",
        fc="#EEF8F6", ec=DATA, title_color=DATA,
    )
    box(
        ax, (0.31, 0.67), (0.21, 0.11),
        "Rank-value\ntokenization",
        "sort nonzero genes by expression\nskip genes absent from dictionary",
        fc="#F4F0F7", ec=TOKEN, title_color=TOKEN,
        title_size=9.5, subtitle_size=7.4,
    )
    arrow(ax, (0.245, 0.725), (0.31, 0.725), TOKEN)
    token_strip(ax, 0.085, 0.515, 0.39, 0.07)
    arrow(ax, (0.415, 0.67), (0.34, 0.59), TOKEN)
    box(
        ax, (0.085, 0.395), (0.39, 0.08),
        "input_ids + attention_mask",
        "CLS + ranked gene tokens + SEP + PAD; max length 2048",
        fc="#FFFFFF", ec=TOKEN, title_color=TOKEN,
    )
    arrow(ax, (0.28, 0.515), (0.28, 0.475), TOKEN)

    # B. model with LoRA adapters
    step_label(ax, 0.57, 0.82, "B. Insert LoRA into Geneformer transformer")
    box(
        ax, (0.57, 0.685), (0.34, 0.095),
        "BertForSequenceClassification",
        "Geneformer-V1-10M backbone; binary classifier head",
        fc="#F7F7F7", ec=MODEL, title_color=MODEL,
        title_size=9.7, subtitle_size=7.7,
    )
    transformer_block(ax, 0.565, 0.540, 0.16, 0.105, "Layer 1")
    transformer_block(ax, 0.745, 0.540, 0.16, 0.105, "Layer 2")
    transformer_block(ax, 0.565, 0.400, 0.16, 0.105, "Layer ...")
    transformer_block(ax, 0.745, 0.400, 0.16, 0.105, "Layer 6")
    arrow(ax, (0.74, 0.685), (0.65, 0.645), MODEL)
    arrow(ax, (0.73, 0.540), (0.745, 0.592), MODEL, lw=0.8)
    arrow(ax, (0.65, 0.540), (0.65, 0.505), MODEL, lw=0.8)
    arrow(ax, (0.73, 0.400), (0.745, 0.452), MODEL, lw=0.8)

    box(
        ax, (0.60, 0.302), (0.27, 0.058),
        "Trainable output head",
        "fatal vs survived logits",
        fc="#EEF5FF", ec=TRAIN, title_color=TRAIN,
        title_size=9.0, subtitle_size=7.3,
    )
    arrow(ax, (0.735, 0.400), (0.735, 0.360), TRAIN)

    # C. training and CV
    step_label(ax, 0.055, 0.315, "C. Patient-held-out LoRA fine-tuning")
    box(
        ax, (0.055, 0.175), (0.24, 0.090),
        "5-fold StratifiedGroupKFold",
        "groups = Sample_ID\nno donor overlap",
        fc="#F2F6FA", ec=TRAIN, title_color=TRAIN,
    )
    box(
        ax, (0.335, 0.175), (0.24, 0.090),
        "Freeze backbone weights",
        "requires_grad = False\nunless name contains LoRA/classifier",
        fc="#FAFAFA", ec=MODEL, title_color=MODEL,
        title_size=9.2, subtitle_size=7.2,
    )
    box(
        ax, (0.615, 0.175), (0.24, 0.090),
        "Train only LoRA + head",
        "r=16, alpha=32, dropout=0.1\nAdamW, 25 epochs",
        fc="#FFF5F7", ec=LORA, title_color=LORA,
        title_size=9.2, subtitle_size=7.2,
    )
    arrow(ax, (0.295, 0.220), (0.335, 0.220), TRAIN)
    arrow(ax, (0.575, 0.220), (0.615, 0.220), TRAIN)

    # Parameter status / result
    ax.text(0.055, 0.135, "Parameter status", fontsize=9.8, fontweight="bold", color=TEXT)
    legend = [
        (MODEL, "Frozen: Geneformer backbone weights"),
        (LORA, "Trainable: LoRA A/B adapters inside Q/K/V/dense Linear modules"),
        (TRAIN, "Trainable: binary classification head"),
        (TOKEN, "Input: raw expression re-tokenized; no precomputed embedding matrix"),
    ]
    for i, (color, text) in enumerate(legend):
        y = 0.110 - i * 0.030
        ax.add_patch(Rectangle((0.058, y - 0.010), 0.015, 0.015,
                               facecolor=color, edgecolor=LINE, linewidth=0.3))
        ax.text(0.083, y - 0.002, text, fontsize=7.7, color=TEXT, ha="left", va="center")

    box(
        ax, (0.56, 0.045), (0.37, 0.095),
        "Domain-adaptation result",
        "V1-10M LoRA: 0.705 +/- 0.135 AUC\nsame-backbone zero-shot: 0.592 +/- 0.094 (+0.113)",
        fc="#FFF8E8", ec=NOTE, title_color=NOTE,
        title_size=9.2, subtitle_size=7.8,
    )

    for ext in ("pdf", "png", "svg"):
        out = args.output_dir / f"fig_geneformer_true_lora_workflow.{ext}"
        fig.savefig(out, bbox_inches="tight")
        print(f"Saved {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
