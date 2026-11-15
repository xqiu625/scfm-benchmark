#!/usr/bin/env python3
"""
Fig 1 — Study-design schematic for the 302006 platelet-OOD FM benchmark.

5-FM primary set: STATE, UCE, Geneformer, scPRINT-2, scGPT.
No STACK, no TranscriptFormer, no mouse arm (mouse moved to 302035).

Four horizontal bands:
  1. INPUT (human platelets only)
  2. STANDARDIZED INPUT + 5 FMs (scGPT now primary, not legacy)
  3. EVALUATION ARMS (patient-level CV across mortality / calibration / LOBO /
     cross-disease / adaptation curve / ablations / LoRA / interpretability)
  4. HEADLINE FINDINGS (STATE bidirectional; PCA-200+LR baseline; scGPT 0.352)

Output: fig1/fig1_overview.png + .pdf + .svg (Nature-figure style: editable
text in SVG/PDF, Arial-first font stack, no top/right spines).
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 8,
})

ROOT = Path("/Users/xinruqiu/Desktop/302006/302006_updategit")
OUT = ROOT / "reports/manuscript/figures/fig1"
OUT.mkdir(parents=True, exist_ok=True)

# 5-FM palette (consistent with the data-driven figures)
FM_COLORS = {
    "STATE":      "#3B6BB0",
    "UCE":        "#5C9D7F",
    "Geneformer": "#7C5C9E",
    "scPRINT-2":  "#D4944B",
    "scGPT":      "#A84A4A",
}
TEXT_DARK = "#222222"
DIM_GREY = "#BBBBBB"
BAND_GREY = "#F4F4F4"

FIGSIZE = (12, 8.2)


def rounded_box(ax, x, y, w, h, label="", sublabel="", face="white", edge="black",
                lw=0.9, fontsize=9, sub_fontsize=7.5, title_weight="bold",
                title_color=None, sub_color="#444444", radius=0.025, alpha=1.0,
                ha="center"):
    box = FancyBboxPatch((x, y), w, h,
                         boxstyle=f"round,pad=0.003,rounding_size={radius}",
                         linewidth=lw, edgecolor=edge, facecolor=face, alpha=alpha)
    ax.add_patch(box)
    tcolor = title_color if title_color is not None else TEXT_DARK
    if sublabel:
        ax.text(x + w / 2, y + h * 0.66, label, ha=ha, va="center",
                fontsize=fontsize, color=tcolor, fontweight=title_weight)
        ax.text(x + w / 2, y + h * 0.32, sublabel, ha=ha, va="center",
                fontsize=sub_fontsize, color=sub_color, linespacing=1.15)
    else:
        ax.text(x + w / 2, y + h / 2, label, ha=ha, va="center",
                fontsize=fontsize, color=tcolor, fontweight=title_weight)


def arrow_down(ax, x, y_from, y_to):
    ax.add_patch(FancyArrowPatch(
        (x, y_from), (x, y_to),
        arrowstyle="-|>,head_width=4,head_length=5.6",
        color="#777777", linewidth=1.0, mutation_scale=8))


def band_label(ax, text, x, y):
    ax.text(x, y, text, fontsize=10, color="#555555", fontweight="bold",
            ha="left", va="center", style="italic")


def main():
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    # ----------------- BAND 1: INPUT (human-only) -----------------
    band_label(ax, "1. INPUT", 1, 95)
    rounded_box(ax, 4, 84, 93, 9, face=BAND_GREY, edge="#555555",
                label="Human platelet scRNA-seq  $\\bullet$  47,977 cells $\\cdot$ 413 donors $\\cdot$ 11 cohorts",
                sublabel="Primary task: COVID + Sepsis subset $\\bullet$ 42,095 cells $\\cdot$ 322 donors $\\bullet$ in-hospital mortality (fatal vs.\\ survived; severe-excluded)",
                fontsize=10.5, sub_fontsize=8.5)

    for xi in [20, 50, 82]:
        arrow_down(ax, xi, 83.5, 79.5)

    # ----------------- BAND 2: STANDARDIZATION + 5 FMs -----------------
    band_label(ax, "2. STANDARDIZED INPUT + 5 FOUNDATION MODELS", 1, 78)
    rounded_box(ax, 4, 70.5, 93, 5.5, face="#EAF3FC", edge="#0072B2",
                label="Unified-input standardization  $\\bullet$  identical raw AnnData $\\rightarrow$ each FM applies its native preprocessing",
                fontsize=10.5)

    fms = ["STATE", "UCE", "Geneformer", "scPRINT-2", "scGPT"]
    tags = ["state-space\nperturbation-aware",
            "cross-species\nuniversal embedding",
            "rank-value MLM\n+LoRA arm (V1-10M)",
            "multi-task\nontology-conditioned",
            "generative MLM\nmasked-gene prediction"]
    card_w, gap = 17.4, 1.0
    x0 = 4
    for i, (fm, tag) in enumerate(zip(fms, tags)):
        rounded_box(ax, x0 + i * (card_w + gap), 60, card_w, 8.5,
                    face=FM_COLORS[fm], edge="black", lw=0.6,
                    label=fm, sublabel=tag,
                    fontsize=10.5, sub_fontsize=7.5,
                    title_color="white", sub_color="white", alpha=0.92)

    # PCA baseline call-out below FMs
    rounded_box(ax, 4, 55.0, 93, 3.2, face="#FFFAEC", edge="#888",
                label="+ PCA-200 + LR baseline on raw counts  $\\bullet$  headline finding: beats UCE and scGPT;"
                      " ties or beats every FM except scPRINT-2 on cross-disease adaptation curve",
                fontsize=9, title_color="#444", title_weight="normal")

    for xi in [20, 50, 82]:
        arrow_down(ax, xi, 54.5, 51.0)

    # ----------------- BAND 3: EVALUATION ARMS -----------------
    band_label(ax, "3. EVALUATION ARMS (patient-level StratifiedGroupKFold, $K{=}5$)", 1, 50.5)
    eval_boxes = [
        ("Within-disease mortality",
         "Per-fold mean LR / RF AUC\nOOF cell-weighted AUC\nDeLong pairwise tests", "#D6EAF8"),
        ("Calibration",
         "Brier score  $\\cdot$  ECE\nLR vs RF comparison\nRF systematically better", "#D6EAF8"),
        ("LOBO cross-cohort",
         "7 held-out cohorts\nMean LOBO AUC + $\\sigma$\nPer-cohort breakdown", "#D6EAF8"),

        ("Cross-disease bidirectional",
         "COVID $\\rightarrow$ Sepsis  +  Sepsis $\\rightarrow$ COVID\n200-iter patient bootstrap CI\nsevere-excluded label", "#FCE4D6"),
        ("Adaptation curve",
         "$n_{\\mathrm{sepsis}} \\in \\{0,10,50,100,500,1000,\\mathrm{all}\\}$\nY-intercept / slope / asymptote\n+ PCA-200+LR baseline", "#FCE4D6"),
        ("Sample efficiency + dropout",
         "5 / 10 / 20 / 50 / 100\\% training\n0--70\\% embedding-dim dropout\nrepresentational redundancy", "#FCE4D6"),

        ("Embedding ablations",
         "PCA dim reduction\nTop-PC removal  $\\cdot$  Norm.\nRandom dim selection", "#E2F0CB"),
        ("LoRA fine-tune (Geneformer)",
         "V1-10M  $\\cdot$  $r{=}16,\\,\\alpha{=}32$\n37 LoRA Linear modules\nProper patient-level CV", "#E2F0CB"),
        ("Interpretability",
         "GSEA on top dim\n+ SHAP attribution\n+ Jaccard overlap", "#E2F0CB"),
    ]
    ebw, ebh = 30.0, 9.5
    egap_x, egap_y = 1.2, 0.7
    ex0, ey0 = 4, 38.0
    for idx, (title, sub, face) in enumerate(eval_boxes):
        r, c = divmod(idx, 3)
        x = ex0 + c * (ebw + egap_x)
        y = ey0 - r * (ebh + egap_y)
        rounded_box(ax, x, y, ebw, ebh, face=face, edge="#333333", lw=0.6,
                    label=title, sublabel=sub,
                    fontsize=9.5, sub_fontsize=7.8, sub_color="#222222")

    for xi in [20, 50, 82]:
        arrow_down(ax, xi, 19.5, 16.0)

    # ----------------- BAND 4: HEADLINE FINDINGS -----------------
    band_label(ax, "4. HEADLINE FINDINGS", 1, 15.5)
    findings = [
        ("STATE",
         "Strongest bidirectional\ncross-disease transfer\nforward 0.798 / reverse 0.707",
         FM_COLORS["STATE"]),
        ("scPRINT-2",
         "Best per-fold LR AUC 0.762\nBest zero-shot Y-intercept\non adaptation curve",
         FM_COLORS["scPRINT-2"]),
        ("PCA-200+LR baseline",
         "Within-disease 0.705\nAdaptation asymptote 0.964\npretrained FM superiority\nis conditional",
         "#555555"),
        ("scGPT — negative result",
         "Anti-correlated\ncross-disease 0.352\npublished as a finding",
         FM_COLORS["scGPT"]),
        ("Two-axis biology",
         "Cytoskeleton/adhesion\n(STATE, UCE)  vs.\ntranslation/coronavirus\n(scPRINT-2, Geneformer)",
         "#444444"),
    ]
    fw, fh = 18.0, 11.5
    fgap = 0.9
    fx0 = 4
    fy0 = 2.5
    for i, (title, sub, color) in enumerate(findings):
        x = fx0 + i * (fw + fgap)
        rounded_box(ax, x, fy0, fw, fh, face=color, edge="black", lw=0.6,
                    label=title, sublabel=sub,
                    fontsize=9.5, sub_fontsize=7.6,
                    title_color="white", sub_color="white", alpha=0.92)

    # ----------------- TITLE -----------------
    fig.suptitle("Out-of-distribution robustness benchmark of single-cell foundation models  $\\cdot$  in-hospital mortality on platelet transcriptomes",
                 fontsize=12, y=0.997, x=0.5, fontweight="bold")
    fig.text(0.5, 0.955,
             "5-FM primary set (STATE / UCE / Geneformer / scPRINT-2 / scGPT) + PCA-200+LR baseline; patient-level $K{=}5$ StratifiedGroupKFold throughout",
             fontsize=9, ha="center", color="#555555", style="italic")

    for ext in ("png", "pdf", "svg"):
        out = OUT / f"fig1_overview.{ext}"
        fig.savefig(out, bbox_inches="tight", dpi=600 if ext == "png" else None)
        print(f"  wrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
