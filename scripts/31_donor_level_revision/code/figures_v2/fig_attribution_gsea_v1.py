#!/usr/bin/env python3
"""Gene-attribution + GSEA cross-model figure (masking thread), v1.

Panel A — top mortality-attributed genes per model: union of each model's
top-8 genes by mean_attribution; cell color = that gene's rank percentile
within the model's own ranked list (darker = more important), +/- = sign of
the attribution. RANK-ONLY display: attribution magnitudes are not on a
common scale across models (classifier logit scales differ, across-gene SD
0.012-0.061), so no raw delta values are shown.

Panel B — GSEA pathway x model NES matrix: union-selected terms significant
(FDR q < 0.25) in >=2 models, ordered by sharing count, plus each model's
strongest model-specific term; dot = significant. scPRINT-2's empty column
is the direct visualization of its complete GSEA null (see R6 text: the
joint-occlusion dose control attributes this to occlusion insensitivity of
the architecture, with redundant encoding ruled out and biological absence
vs. method limitation open).

Inputs: results_gene_attribution_gsea/{*_masking/*_ranked.csv,
cross_model/gsea_nes_matrix_*.csv} (produced by
compare_masking_attribution_5models.py). Stdlib + numpy + matplotlib only.
"""
import csv
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import MODEL_COLOR, save_fig

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     os.pardir, os.pardir))
GSEA_DIR = os.path.join(BASE, "results_gene_attribution_gsea")
XMOD_DIR = os.path.join(GSEA_DIR, "cross_model")

MODELS = ["geneformer", "state", "scgpt", "uce", "scprint2"]
DISPLAY = {"geneformer": "Geneformer", "state": "STATE", "scgpt": "scGPT",
           "uce": "UCE", "scprint2": "scPRINT-2"}
GENESETS = ["hallmark", "kegg", "reactome"]
TOP_PER_MODEL = 7
SHARED_QUOTA = {"hallmark": 5, "kegg": 3, "reactome": 3}  # 11 shared rows
# + one strongest model-specific term per model = 16 rows total

# spec/P0-style anchors, verified against the CSVs at build time (printed, not asserted)
ANCHOR_SIG_COUNTS = {  # from the 2026-09-11 milestone summary, cross-checked below
    "geneformer": (6, 23, 2), "state": (5, 7, 2), "scgpt": (28, 63, 257),
    "uce": (23, 70, 82), "scprint2": (0, 0, 0),
}


def load_ranked(model):
    path = os.path.join(GSEA_DIR, f"{model}_masking",
                        f"{model}_masking_gene_attribution_ranked.csv")
    with open(path) as f:
        return {r["gene"]: float(r["mean_attribution"]) for r in csv.DictReader(f)}


def load_nes_matrix(gs):
    path = os.path.join(XMOD_DIR, f"gsea_nes_matrix_{gs}.csv")
    with open(path) as f:
        return list(csv.DictReader(f))


def panel_a(ax, ranked):
    union, seen = [], set()
    for m in MODELS:
        for g in sorted(ranked[m], key=ranked[m].get, reverse=True)[:TOP_PER_MODEL]:
            if g not in seen:
                seen.add(g)
                union.append(g)
    # rank of each gene in each model (1 = top); genes absent from a model's list = None
    rank_of = {m: {g: i + 1 for i, g in enumerate(
        sorted(ranked[m], key=ranked[m].get, reverse=True))} for m in MODELS}
    n_lists = {m: len(ranked[m]) for m in MODELS}

    def best_rank(g):
        return min(rank_of[m].get(g, 10 ** 9) for m in MODELS)

    union.sort(key=best_rank)
    G, M = len(union), len(MODELS)
    mat = np.full((G, M), np.nan)  # rank percentile (0 = top)
    signs = np.full((G, M), "", dtype=object)
    for j, m in enumerate(MODELS):
        for i, g in enumerate(union):
            if g in rank_of[m]:
                mat[i, j] = (rank_of[m][g] - 1) / max(n_lists[m] - 1, 1)
                signs[i, j] = "+" if ranked[m][g] > 0 else "−"

    cmap = plt.get_cmap("viridis_r").copy()
    cmap.set_bad("#E8E8E8")
    im = ax.imshow(np.ma.masked_invalid(mat), aspect="auto", cmap=cmap, vmin=0, vmax=0.5)
    for i in range(G):
        for j in range(M):
            if signs[i, j]:
                ax.text(j, i, signs[i, j], ha="center", va="center",
                        fontsize=7, color="white" if mat[i, j] < 0.25 else "dimgray")
            else:
                ax.text(j, i, "·", ha="center", va="center", fontsize=7, color="#BBBBBB")
    ax.set_xticks(range(M))
    ax.set_xticklabels([DISPLAY[m] for m in MODELS], fontsize=8, rotation=20, ha="right")
    for j, m in enumerate(MODELS):
        ax.get_xticklabels()[j].set_color(MODEL_COLOR[DISPLAY[m]])
    ax.set_yticks(range(G))
    ax.set_yticklabels(union, fontsize=7)
    ax.set_title("A  Top mortality-attributed genes per model\n"
                 "(rank within model; +/− = attribution sign)", fontsize=9, loc="left")
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    return im


def panel_b(ax):
    rows_by_gs = {gs: load_nes_matrix(gs) for gs in GENESETS}
    cols = [DISPLAY[m] for m in MODELS]

    # verify sig-count anchors against the matrices
    for m, col in zip(MODELS, cols):
        counts = tuple(sum(int(r[f"sig_{col}"]) for r in rows_by_gs[gs]) for gs in GENESETS)
        status = "OK" if counts == ANCHOR_SIG_COUNTS[m] else f"MISMATCH vs anchor {ANCHOR_SIG_COUNTS[m]}"
        print(f"sig counts {col}: {counts} ({status})", flush=True)

    selected, seen = [], set()
    for gs in GENESETS:  # shared (>=2 models) rows, quota per geneset
        rows = sorted(rows_by_gs[gs],
                      key=lambda r: (-int(r["n_models_sig"]),
                                     -max(abs(float(r[f"NES_{c}"])) if r[f"NES_{c}"] else 0.0
                                          for c in cols)))
        taken = 0
        for r in rows:
            if taken >= SHARED_QUOTA[gs]:
                break
            if int(r["n_models_sig"]) >= 2 and r["Term"] not in seen:
                seen.add(r["Term"])
                selected.append((gs, r))
                taken += 1
    # one strongest model-specific term per model (for the model-specific edge)
    specific = []
    for m, col in zip(MODELS, cols):
        best, best_key = None, None
        for gs in GENESETS:
            for r in rows_by_gs[gs]:
                if int(r["n_models_sig"]) == 1 and r[f"sig_{col}"] == "1" \
                        and r["Term"] not in seen:
                    key = abs(float(r[f"NES_{col}"]))
                    if best_key is None or key > best_key:
                        best, best_key = (gs, r), key
        if best is not None:
            seen.add(best[1]["Term"])
            specific.append(best)
    selected += specific

    P, M = len(selected), len(MODELS)
    nes = np.full((P, M), np.nan)
    sig = np.zeros((P, M), dtype=bool)
    labels = []
    for i, (gs, r) in enumerate(selected):
        labels.append(r["Term"] if len(r["Term"]) <= 44 else r["Term"][:43] + "…")
        for j, col in enumerate(cols):
            if r[f"NES_{col}"]:
                nes[i, j] = float(r[f"NES_{col}"])
                sig[i, j] = r[f"sig_{col}"] == "1"

    vmax = np.nanmax(np.abs(nes))
    im = ax.imshow(np.ma.masked_invalid(nes), aspect="auto", cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax)
    n_shared = P - len(specific)
    if 0 < n_shared < P:
        ax.axhline(n_shared - 0.5, color="black", lw=0.8, alpha=0.5)
    for i in range(P):
        for j in range(M):
            if sig[i, j]:
                ax.text(j, i, "●", ha="center", va="center", fontsize=6.5,
                        color="black")
    ax.set_xticks(range(M))
    ax.set_xticklabels(cols, fontsize=8, rotation=20, ha="right")
    for j, col in enumerate(cols):
        ax.get_xticklabels()[j].set_color(MODEL_COLOR[col])
    ax.set_yticks(range(P))
    ax.set_yticklabels(labels, fontsize=7)
    gs_marks = {}
    for i, (gs, _) in enumerate(selected):
        gs_marks.setdefault(gs, []).append(i)
    ylabel = "  ".join(f"{k}({len(v)})" for k, v in gs_marks.items())
    ax.set_title("B  GSEA pathway emphasis (NES; ● = FDR q < 0.25)\n"
                 f"rows: {ylabel}; below line: top model-specific term",
                 fontsize=9, loc="left")
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    return im


def main():
    ranked = {m: load_ranked(m) for m in MODELS}
    fig = plt.figure(figsize=(13, 7.0))
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)
    im1 = panel_a(ax1, ranked)
    im2 = panel_b(ax2)
    fig.tight_layout()
    cb1 = fig.colorbar(im1, ax=ax1, fraction=0.03, pad=0.01, ticks=[0, 0.25, 0.5])
    cb1.ax.set_yticklabels(["top", "25%", "50%+"])
    cb1.set_label("rank in model", fontsize=7)
    cb1.ax.tick_params(labelsize=7, length=0)
    cb2 = fig.colorbar(im2, ax=ax2, fraction=0.03, pad=0.01)
    cb2.set_label("NES", fontsize=7)
    cb2.ax.tick_params(labelsize=7, length=0)
    fig.text(0.5, -0.03,
             "Masking-based attribution on the canonical within-COVID cohort (240 patients, 25 fatal; "
             "10,737 cells); delta = full minus gene-occluded mortality-classifier logit. "
             "Cross-model comparison is rank-based only (classifier logit scales differ). "
             "scPRINT-2 column shows its complete GSEA null; the joint-occlusion dose control "
             "(N=10/30/128 flat) attributes this to occlusion insensitivity of the architecture, "
             "not to redundant gene encoding.",
             ha="center", fontsize=7.3, color="dimgray", wrap=True)
    save_fig(fig, "Figure_gene_attribution_gsea")


if __name__ == "__main__":
    main()
