#!/usr/bin/env python3
"""Cross-model comparison of the five masking-based gene-attribution ranked
lists and their GSEA results (Stage 5/6 of the unified masking design).

Reads (all under scripts/31_donor_level_revision/results_gene_attribution_gsea/):
  <model>_masking/<model>_masking_gene_attribution_ranked.csv
  <model>_masking/<model>_masking_gsea_{hallmark,kegg,reactome}.csv

Writes to results_gene_attribution_gsea/cross_model/:
  common_genes.csv                genes ranked by all five models (the
                                  post-hoc 5-way comparison set)
  cross_model_spearman.csv        5x5 Spearman rank correlation of
                                  mean_attribution on common genes
  top50_overlap.csv               5x5 Jaccard of top-50 genes
  top100_overlap.csv              5x5 Jaccard of top-100 genes
  consensus_ranking.csv           per-gene ranks in each model + mean rank
  sign_flip_genes.csv             genes top-100 in >=1 model whose sign
                                  differs across the models that rank it
  gsea_nes_matrix_<gs>.csv        term x model NES (+ FDR q per model)
  gsea_sig_jaccard_<gs>.csv       5x5 Jaccard of FDR<0.25 term sets
  gsea_consensus_terms.csv        FDR<0.25 in >=4 of 5 models
  gsea_model_specific_terms.csv   FDR<0.25 in exactly 1 model

Comparisons are RANK-ONLY: mean_attribution magnitudes are not on a common
scale across models (each model's mortality-direction classifier has its own
logit scale; across-gene SD spans 0.012-0.061), so no cross-model delta
magnitudes are computed or compared here.

Stdlib + numpy only (no pandas) so it runs with the local python3.12 that
has matplotlib but not pandas; numpy is used for arrays only, not I/O.
"""
import csv
import math
import os

import numpy as np

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir,
                    "results_gene_attribution_gsea")
BASE = os.path.normpath(BASE)
OUT_DIR = os.path.join(BASE, "cross_model")
os.makedirs(OUT_DIR, exist_ok=True)

MODELS = ["geneformer", "state", "scgpt", "uce", "scprint2"]
DISPLAY = {"geneformer": "Geneformer", "state": "STATE", "scgpt": "scGPT",
           "uce": "UCE", "scprint2": "scPRINT-2"}
GENESETS = ["hallmark", "kegg", "reactome"]
FDR_CUT = 0.25  # same threshold convention as the per-model GSEA runs


def load_ranked(model):
    path = os.path.join(BASE, f"{model}_masking",
                        f"{model}_masking_gene_attribution_ranked.csv")
    with open(path) as f:
        return {r["gene"]: float(r["mean_attribution"]) for r in csv.DictReader(f)}


def load_gsea(model, gs):
    path = os.path.join(BASE, f"{model}_masking",
                        f"{model}_masking_gsea_{gs}.csv")
    with open(path) as f:
        return {r["Term"]: {"NES": float(r["NES"]), "FDR": float(r["FDR q-val"])}
                for r in csv.DictReader(f)}


def avg_ranks(values):
    """Average-tie ranks, ascending (smallest value -> rank 1)."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        r = (i + j + 2) / 2.0
        for k in range(i, j + 1):
            ranks[order[k]] = r
        i = j + 1
    return ranks


def spearman(x, y):
    rx, ry = avg_ranks(x), avg_ranks(y)
    n = len(rx)
    mx, my = sum(rx) / n, sum(ry) / n
    cov = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    vx = sum((a - mx) ** 2 for a in rx)
    vy = sum((b - my) ** 2 for b in ry)
    return cov / math.sqrt(vx * vy) if vx > 0 and vy > 0 else float("nan")


def jaccard(a, b):
    u = a | b
    return len(a & b) / len(u) if u else float("nan")


def write_matrix(path, labels, mat, fmt="{:.4f}"):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["model"] + labels)
        for i, lab in enumerate(labels):
            w.writerow([lab] + [fmt.format(v) for v in mat[i]])


def main():
    ranked = {m: load_ranked(m) for m in MODELS}
    for m in MODELS:
        print(f"{DISPLAY[m]:12s} ranked genes: {len(ranked[m])}", flush=True)

    common = sorted(set.intersection(*[set(ranked[m]) for m in MODELS]))
    print(f"Common genes across all five ranked lists: {len(common)}", flush=True)
    with open(os.path.join(OUT_DIR, "common_genes.csv"), "w", newline="") as f:
        csv.writer(f).writerows([["gene"]] + [[g] for g in common])

    labels = [DISPLAY[m] for m in MODELS]

    # 1. pairwise Spearman on common genes (rank-only comparison)
    sp = np.zeros((5, 5))
    for i, a in enumerate(MODELS):
        for j, b in enumerate(MODELS):
            sp[i, j] = spearman([ranked[a][g] for g in common],
                                [ranked[b][g] for g in common])
    write_matrix(os.path.join(OUT_DIR, "cross_model_spearman.csv"), labels, sp, "{:.3f}")
    off = [sp[i, j] for i in range(5) for j in range(5) if i < j]
    print(f"Spearman range across pairs: {min(off):.3f} .. {max(off):.3f}", flush=True)

    # 2. top-N Jaccard overlaps (by each model's own mean_attribution order)
    for topn in (50, 100):
        tops = {m: set(sorted(ranked[m], key=ranked[m].get, reverse=True)[:topn])
                for m in MODELS}
        jac = np.array([[jaccard(tops[a], tops[b]) for b in MODELS] for a in MODELS])
        write_matrix(os.path.join(OUT_DIR, f"top{topn}_overlap.csv"), labels, jac, "{:.3f}")

    # 3. consensus ranking: per-model rank (1 = highest attribution) + mean
    rank_of = {m: {g: r for g, r in zip(common, avg_ranks([-ranked[m][g] for g in common]))}
               for m in MODELS}
    rows = []
    for g in common:
        rs = [rank_of[m][g] for m in MODELS]
        signs = "".join("+" if ranked[m][g] > 0 else "-" for m in MODELS)
        rows.append([g] + [f"{r:.0f}" for r in rs] + [f"{sum(rs)/len(rs):.1f}", signs])
    rows.sort(key=lambda r: float(r[6]))
    with open(os.path.join(OUT_DIR, "consensus_ranking.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene"] + [f"rank_{DISPLAY[m]}" for m in MODELS] + ["mean_rank", "signs"])
        w.writerows(rows)
    print("Top-10 consensus genes:", ", ".join(r[0] for r in rows[:10]), flush=True)

    # 4. sign flips: top-100 in >=1 model, ranked by >=3 models, signs differ
    flips = []
    for g in common:
        top100 = [m for m in MODELS if rank_of[m][g] <= 100]
        if not top100:
            continue
        signs = {m: (ranked[m][g] > 0) for m in MODELS}
        if len(set(signs.values())) > 1:
            flips.append([g, ",".join(top100),
                          "".join("+" if ranked[m][g] > 0 else "-" for m in MODELS)]
                         + [ranked[m][g] for m in MODELS])
    with open(os.path.join(OUT_DIR, "sign_flip_genes.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "top100_in", "signs"] + [f"attr_{DISPLAY[m]}" for m in MODELS])
        w.writerows(flips)
    print(f"Sign-flip genes (top-100 somewhere, disagreeing signs): {len(flips)}", flush=True)

    # 5. GSEA matrices + consensus/model-specific terms
    gsea = {(m, gs): load_gsea(m, gs) for m in MODELS for gs in GENESETS}
    all_consensus, all_specific = [], []
    for gs in GENESETS:
        terms = sorted(set.union(*[set(gsea[(m, gs)]) for m in MODELS]))
        sig = {t: [m for m in MODELS if t in gsea[(m, gs)] and gsea[(m, gs)][t]["FDR"] < FDR_CUT]
               for t in terms}
        with open(os.path.join(OUT_DIR, f"gsea_nes_matrix_{gs}.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["Term"] + [f"NES_{DISPLAY[m]}" for m in MODELS]
                       + [f"sig_{DISPLAY[m]}" for m in MODELS] + ["n_models_sig"])
            for t in terms:
                w.writerow([t]
                           + [f"{gsea[(m, gs)][t]['NES']:.3f}" if t in gsea[(m, gs)] else "" for m in MODELS]
                           + [int(t in gsea[(m, gs)] and gsea[(m, gs)][t]["FDR"] < FDR_CUT) for m in MODELS]
                           + [len(sig[t])])
        sets = {m: {t for t in terms if m in sig[t]} for m in MODELS}
        jac = np.array([[jaccard(sets[a], sets[b]) for b in MODELS] for a in MODELS])
        write_matrix(os.path.join(OUT_DIR, f"gsea_sig_jaccard_{gs}.csv"), labels, jac, "{:.3f}")
        n_sig = {m: len(sets[m]) for m in MODELS}
        print(f"{gs:9s} sig terms: " + ", ".join(f"{DISPLAY[m]}={n_sig[m]}" for m in MODELS),
              flush=True)
        for t in terms:
            rec = [gs, t, len(sig[t]), ",".join(DISPLAY[m] for m in sig[t]),
                   ",".join(f"{gsea[(m, gs)][t]['NES']:.3f}" for m in sig[t])]
            if len(sig[t]) >= 4:
                all_consensus.append(rec)
            elif len(sig[t]) == 1:
                all_specific.append(rec)
    for path, recs in [("gsea_consensus_terms.csv", all_consensus),
                       ("gsea_model_specific_terms.csv", all_specific)]:
        with open(os.path.join(OUT_DIR, path), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["geneset", "Term", "n_models_sig", "models", "NES_in_sig_models"])
            w.writerows(recs)
    print(f"Consensus terms (>=4/5 models): {len(all_consensus)}; "
          f"model-specific (exactly 1): {len(all_specific)}", flush=True)
    print(f"Outputs written to {OUT_DIR}", flush=True)


if __name__ == "__main__":
    main()
