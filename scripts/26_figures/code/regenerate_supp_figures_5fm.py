#!/usr/bin/env python3
"""
Generate new supplementary figures for the 302006 manuscript under the 5-FM
primary set (STATE, UCE, Geneformer, scPRINT-2, scGPT; STACK/TF dropped).

Adds two data-driven supplementary panels:
  - Supp Fig: per-cohort standard 5-fold OOF AUC heatmap (per_cohort_oof.csv)
  - Supp Fig: DeLong pairwise AUC-difference + significance heatmap (delong_pairwise.csv)

Style follows the nature-figure skill convention used in the main figures:
Arial, font.size=7, editable SVG/PDF, no top/right spines, lw 0.8.
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.8,
    "legend.frameon": False,
})

ROOT = Path("/Users/xinruqiu/Desktop/302006/302006_updategit")
DATA = ROOT / "scripts/26_figures/data"
OUT = ROOT / "reports/manuscript/figures"
OUT.mkdir(parents=True, exist_ok=True)

FMS = ["STATE", "UCE", "Geneformer", "scPRINT-2", "scGPT"]
DROP = {"STACK", "STACK ortholog", "TranscriptFormer", "TF", "TF-Metazoa"}

PALETTE = {
    "STATE":      "#3B6BB0",
    "UCE":        "#5C9D7F",
    "Geneformer": "#7C5C9E",
    "scPRINT-2":  "#D4944B",
    "scGPT":      "#A84A4A",
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


def save_pub(fig, stem):
    fig.savefig(OUT / f"{stem}.png", dpi=600, bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Supp Fig 4: per-cohort OOF AUC heatmap (standard 5-fold, not LOBO)
# ---------------------------------------------------------------------------
def make_supp_per_cohort_oof():
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
    ax.set_xticklabels(M.columns, rotation=45, ha="right", fontsize=5.8)
    ax.set_xlabel("Cohort (Data\\_NO)")
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
    save_pub(fig, "fig_supp_per_cohort_oof")
    print("  ✓ fig_supp_per_cohort_oof.png/pdf/svg")


# ---------------------------------------------------------------------------
# Supp Fig 5: DeLong pairwise AUC-difference heatmap with significance overlay
# ---------------------------------------------------------------------------
def make_supp_delong():
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
        # Place delta with sign: row=a, col=b means a-b (positive => a better)
        delta_matrix.loc[a, b] = d
        delta_matrix.loc[b, a] = -d
        p_matrix.loc[a, b] = p
        p_matrix.loc[b, a] = p

    fig, ax = plt.subplots(figsize=(4.0, 3.4))
    im = ax.imshow(delta_matrix.values, cmap=plt.cm.RdBu_r, vmin=-0.15, vmax=0.15)

    for i, ra in enumerate(models_present):
        for j, cb in enumerate(models_present):
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
    ax.set_xlabel("Model $B$")
    ax.set_ylabel("Model $A$")
    ax.set_title("DeLong pairwise $\\Delta$AUC = $\\mathrm{AUC}_A - \\mathrm{AUC}_B$\n(* $p<0.05$,  ** $p<10^{-3}$,  *** $p<10^{-4}$)",
                 loc="left", fontsize=7.5, pad=4)

    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label("$\\Delta$AUC", fontsize=6.5)
    cb.ax.tick_params(labelsize=5.5)
    cb.outline.set_linewidth(0.5)

    fig.tight_layout()
    save_pub(fig, "fig_supp_delong_significance")
    print("  ✓ fig_supp_delong_significance.png/pdf/svg")


# ---------------------------------------------------------------------------
def main():
    print("Generating supplementary figures (5-FM benchmark)...")
    make_supp_per_cohort_oof()
    make_supp_delong()
    print("Done.")


if __name__ == "__main__":
    main()
