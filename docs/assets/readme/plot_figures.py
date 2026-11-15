#!/usr/bin/env python3
"""Render README graphics from a frozen coordinate table; no UMAP/model fitting."""

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
SOURCE = Path(
    "scripts/24_option_c_canonical_input/results/"
    "option_c_bio_interpretation/umap_scprint2.csv"
)
INK = "#172B44"
MUTED = "#526477"
TEAL = "#087E8B"
SEED = 42


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def workflow(out):
    fig, ax = plt.subplots(figsize=(12.8, 5.35))
    fig.subplots_adjust(0, 0, 1, 1)
    fig.set_facecolor("#F4F7FB")
    ax.set(xlim=(0, 1280), ylim=(535, 0))
    ax.axis("off")
    ax.text(40, 33, "SCFM-BENCHMARK  /  PORTABLE WORKFLOW", color=TEAL,
            fontsize=10, weight="bold", va="top")
    ax.text(40, 70, "From frozen embeddings to donor-level evaluation", color=INK,
            fontsize=22, weight="bold", va="top")
    ax.text(40, 118, "One set of cells. Shared folds. Inspectable predictions.",
            color=MUTED, fontsize=12, va="top")

    steps = [
        ("Supply inputs", "Cell metadata", "cell_id / donor_id",
         "cv_group / binary label", "Frozen embedding CSVs"),
        ("Align & validate", "Join rows by cell_id", "Check identities and values",
         "Reject inconsistent labels"),
        ("Split by group", "Shared folds for all models", "Keep cv_group intact",
         "Repeat across seeds"),
        ("Fit on training cells", "Scale / optional PCA", "Fit logistic regression",
         "Predict held-out cells"),
        ("Score donors", "Average cell probabilities", "within each held-out donor",
         "AUROC, AP and Brier", "Training-prior baseline"),
    ]
    for i, (title, *lines) in enumerate(steps):
        x = 40 + i * 246
        ax.add_patch(FancyBboxPatch((x, 177), 216, 201,
                     boxstyle="round,pad=0,rounding_size=12", facecolor="white",
                     edgecolor="#DCE4ED", linewidth=1))
        ax.text(x + 17, 196, f"0{i + 1}", color=TEAL, fontsize=11,
                weight="bold", va="top")
        ax.text(x + 17, 235, title, color=INK, fontsize=12.5,
                weight="bold", va="top")
        for j, line in enumerate(lines):
            ax.text(x + 17, 281 + j * 26, line, color=MUTED,
                    fontsize=10.2, va="top")
        if i < len(steps) - 1:
            ax.add_patch(FancyArrowPatch((x + 221, 275), (x + 241, 275),
                         arrowstyle="-|>", mutation_scale=13,
                         linewidth=1.4, color=TEAL))

    ax.add_patch(FancyBboxPatch((40, 411), 1200, 69,
                 boxstyle="round,pad=0,rounding_size=12", facecolor=INK,
                 edgecolor="none"))
    ax.text(61, 432, "Inspect every run", color="white", fontsize=13,
            weight="bold", va="top")
    ax.text(280, 434, "Predictions  /  Folds  /  Metrics  /  HTML report  /  Input and output hashes",
            color="#D9E6F1", fontsize=12, va="top")
    ax.text(40, 504, "Embeddings are prepared upstream. Training weights cells equally; metrics weight donors equally.",
            color=MUTED, fontsize=10, va="top")
    fig.savefig(out / "workflow.svg", metadata={"Date": None})
    fig.savefig(out / "workflow.png", dpi=160)
    plt.close(fig)


def historical_umap(df, out):
    cohorts = sorted(df.cohort.unique())
    cohort_colors = ["#087E8B", "#E69F00", "#0072B2", "#CC6677",
                     "#6654A4", "#009E73", "#882255", "#8C6D31", "#777777",
                     "#56B4E9", "#AA4499"]
    if len(cohorts) > len(cohort_colors):
        raise ValueError("Add distinct cohort colors before plotting more cohorts")
    cohort_palette = dict(zip(cohorts, cohort_colors))
    severity_palette = {"mild": "#4AA6A8", "moderate": "#D8A323",
                        "recovered": "#4477AA", "severe": "#8C6BB1",
                        "fatal": "#C44E52"}
    if set(df.severity) - severity_palette.keys():
        raise ValueError("Unexpected severity labels; preserve and label them explicitly")

    order = np.random.default_rng(SEED).permutation(len(df))
    points = df.iloc[order]
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 6.7), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.055, right=0.98, top=0.75, bottom=0.26, wspace=0.15)
    fig.set_facecolor("white")
    fig.text(0.055, 0.95, "ARCHIVED RESEARCH  /  SCIENTIFIC VISUALIZATION",
             fontsize=10, weight="bold", color=TEAL)
    fig.text(0.055, 0.89, "One embedding, two views", fontsize=23, weight="bold", color=INK)
    fig.text(0.055, 0.84,
             f"scPRINT-2  |  {len(df):,} stored cell coordinates  |  All cells shown in both panels",
             fontsize=12, color=MUTED)
    for ax, column, palette, title in zip(
            axes, ["cohort", "severity"], [cohort_palette, severity_palette],
            ["A   Source cohort", "B   Recorded severity"]):
        colors = [palette[value] for value in points[column]]
        ax.scatter(points.UMAP1, points.UMAP2, c=colors, s=2.8, alpha=0.55,
                   linewidths=0, rasterized=True)
        ax.set_title(title, loc="left", pad=16, fontsize=13, weight="bold", color=INK)
        ax.set_xlabel("UMAP 1", fontsize=10, color=MUTED)
        ax.set_ylabel("UMAP 2", fontsize=10, color=MUTED)
        ax.set_aspect("equal", adjustable="box")
        ax.margins(0.05)
        ax.tick_params(length=0, labelsize=8, colors=MUTED)
        for edge in ["top", "right"]:
            ax.spines[edge].set_visible(False)
        for edge in ["left", "bottom"]:
            ax.spines[edge].set_color("#D3DDE7")
        labels = [label for label in palette if label in set(df[column])]
        handles = [Line2D([], [], linestyle="", marker="o", markersize=5,
                          markeredgewidth=0, color=palette[label]) for label in labels]
        ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.16),
                  ncol=4 if column == "cohort" else 3, frameon=False,
                  fontsize=8.5, columnspacing=1.1,
                  handlelength=0.9, handletextpad=0.4)
    fig.text(0.055, 0.072, "Same coordinates, different annotations. Geometry is reused; no projection or model was refitted.",
             color=MUTED, fontsize=10)
    fig.text(0.055, 0.037, "Historical input predates the QC rebuild. Visual separation is not a performance score or proof of biological validity.",
             color=MUTED, fontsize=9)
    fig.savefig(out / "historical_umap.png", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True, help="New output directory")
    args = parser.parse_args()
    source = ROOT / SOURCE
    source_hash = sha256(source)
    df = pd.read_csv(source, usecols=["Cell_ID", "UMAP1", "UMAP2", "cohort", "severity"],
                     dtype={"Cell_ID": str, "cohort": str, "severity": str})
    if df.empty or df.isna().any().any() or df.Cell_ID.duplicated().any():
        raise ValueError("Coordinates require unique cell IDs and complete annotations")
    if not np.isfinite(df[["UMAP1", "UMAP2"]].to_numpy()).all():
        raise ValueError("Non-finite coordinates")
    args.out.mkdir(parents=True, exist_ok=False)
    plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none",
                         "svg.hashsalt": "scfm-readme", "savefig.facecolor": "auto"})
    workflow(args.out)
    historical_umap(df, args.out)
    if sha256(source) != source_hash:
        raise RuntimeError("Source coordinate table changed during rendering")
    manifest = {
        "source_file": str(SOURCE), "source_sha256": source_hash,
        "renderer_sha256": sha256(Path(__file__)), "n_cells": len(df),
        "annotation_counts": {key: df[key].value_counts().sort_index().to_dict()
                              for key in ["cohort", "severity"]},
        "plotting": {"all_rows_used": True, "same_coordinates_in_both_panels": True,
                     "draw_order_seed": SEED, "marker_size": 2.8, "alpha": 0.55,
                     "umap_refitted": False, "model_refitted": False},
        "versions": {"matplotlib": matplotlib.__version__, "numpy": np.__version__,
                     "pandas": pd.__version__},
        "outputs": {name: sha256(args.out / name)
                    for name in ["workflow.svg", "workflow.png", "historical_umap.png"]},
    }
    (args.out / "provenance.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
