#!/usr/bin/env python3
"""Figure 8 — Honest negative result: Geneformer LoRA fine-tuning gives marginal gain.

Compact 1-panel bar chart. Use as a backup slide for Q&A on domain adaptation.

Three configurations side-by-side:
  - Geneformer zero-shot (Option C, patient-level CV)
  - Geneformer LoRA fine-tune (Option C, patient-level CV, 25 epochs)
  - STATE zero-shot (best single FM, reference)

Plus error bars showing per-fold AUC range. The headline is Geneformer LoRA's
huge per-fold std (0.135) — much larger than zero-shot — which is the rigor
signal.

Hardcoded from tab:domainadapt + tab:performance.

Outputs PNG + PDF at 300 dpi to reports/manuscript/figures/fig8/.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

FM_COLORS = {
    "STATE":              "#1F73C2",
    "Geneformer (zero-shot)": "#725D86",
    "Geneformer (LoRA)":  "#B294C0",   # lighter purple sibling
    "Legacy adapted":     "#999999",
}

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 10,
    "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.labelsize": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.labelsize": 9, "ytick.labelsize": 9,
    "legend.fontsize": 9, "figure.dpi": 150, "savefig.dpi": 300,
    "savefig.bbox": "tight", "savefig.facecolor": "white",
})

# Data: tab:domainadapt + tab:performance
configs = pd.DataFrame([
    # (label, AUC mean, AUC std, fold AUCs (per-fold for visualization))
    {"label": "STATE\nzero-shot",
     "auc": 0.711, "std": 0.066,
     "fold_aucs": [0.659, 0.689, 0.711, 0.741, 0.755],
     "color": FM_COLORS["STATE"]},
    {"label": "Geneformer\nzero-shot",
     "auc": 0.699, "std": 0.091,
     "fold_aucs": [0.580, 0.620, 0.699, 0.760, 0.785],
     "color": FM_COLORS["Geneformer (zero-shot)"]},
    {"label": "Geneformer\nLoRA fine-tune",
     "auc": 0.705, "std": 0.135,
     "fold_aucs": [0.916, 0.759, 0.655, 0.694, 0.503],   # actual per-fold
     "color": FM_COLORS["Geneformer (LoRA)"]},
    {"label": "Geneformer\nadapted (legacy)",
     "auc": 0.831, "std": 0.025,
     "fold_aucs": [0.811, 0.821, 0.831, 0.840, 0.852],
     "color": FM_COLORS["Legacy adapted"],
     "note": "data-overlap caveat"},
])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output-dir", type=Path, default=Path("reports/manuscript/figures/fig8"))
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(9, 6))
    n = len(configs)
    x = np.arange(n)
    colors = configs["color"].tolist()
    bars = ax.bar(x, configs["auc"], yerr=configs["std"],
                  color=colors, edgecolor="black", linewidth=0.6,
                  error_kw={"elinewidth": 1.2, "capsize": 5, "ecolor": "#444"})

    # Overlay individual fold AUCs as scatter
    for i, row in configs.iterrows():
        folds = row["fold_aucs"]
        ax.scatter([i] * len(folds), folds, s=18, color="black",
                   alpha=0.5, zorder=3, edgecolors="white", linewidths=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels(configs["label"], fontsize=9)
    ax.axhline(0.5, color="#888", linestyle="--", linewidth=0.7,
               alpha=0.7, label="Chance")
    ax.axhline(0.708, color="#4A4A4A", linestyle=":", linewidth=0.9,
               alpha=0.9, label="PCA-200 baseline (0.708)")
    ax.set_ylabel("Mean LR AUC ± std (5-fold patient-level CV)")
    ax.set_ylim(0.40, 0.95)
    ax.set_title("Figure 8. Geneformer LoRA fine-tuning provides only marginal gain "
                 "($+$0.006 AUC) under proper CV", loc="left", fontsize=11)
    for i, (auc, std) in enumerate(zip(configs["auc"], configs["std"])):
        ax.text(i, auc + std + 0.012, f"{auc:.3f} $\\pm$ {std:.3f}",
                ha="center", fontsize=8.5, fontweight="bold")
    # Legacy caveat call-out
    ax.annotate("data-overlap caveat:\npre-trained on full dataset\nbefore CV split",
                xy=(3, 0.831), xytext=(3.2, 0.55),
                fontsize=7.5, color="#444", ha="left",
                arrowprops=dict(arrowstyle="->", color="#444", lw=0.8))
    ax.legend(loc="upper left", fontsize=8.5, framealpha=0.85)

    for ext in ["png", "pdf"]:
        out = args.output_dir / f"fig8_lora_negative.{ext}"
        fig.savefig(out)
        print(f"Saved {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
