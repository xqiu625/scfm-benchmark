#!/usr/bin/env python3
"""
scPRINT-2 multi-gene JOINT occlusion control, added 2026-09-11 in direct
response to a methodological concern about scprint2_masking_attribution.py's
null GSEA result: scPRINT-2's single-gene attribution dynamic range
(across-gene SD of mean_attribution = 0.0121) is ~1/5 of Geneformer's
(0.0609) and roughly half of scGPT/STATE/UCE's (0.019-0.022), while its
mean per-gene SD across cells (within-gene cell-to-cell noise = 0.1145) is
~9.5x its own across-gene signal -- by far the worst noise/signal ratio of
any of the 5 models (others range 0.4-1.7). This raises a real, testable
ambiguity: is scPRINT-2's null GSEA genuine "biological diffusion" (mortality
information is redundantly encoded across many genes, so removing any ONE
barely moves the score even though the information is real), or is it a
generic INSENSITIVITY of this model's forward pass / classifier direction to
single-token occlusion (in which case the null GSEA is a measurement-method
artifact, not a biological finding)?

Design: for each cell, using scPRINT-2's OWN already-computed single-gene
ranked list (scprint2_masking_gene_attribution_ranked.csv) as the importance
ranking, jointly mask (in ONE forward pass) the top-N most-POSITIVELY- and
top-N most-NEGATIVELY-attributed candidate genes present in that cell, at
two tiers N in {10, 30} for a dose-response check -- plus a random-N control
(3 draws, averaged) matched in size but NOT selected by attribution, drawn
from the same per-cell candidate set. Three comparisons per (N, polarity):
  joint_delta        = s_full - s_masked(topN)
  sum_individual      = sum of each gene's own single-gene mean_attribution
                         (from the existing ranked CSV) -- the naive
                         "expected if additive" prediction
  random_delta        = s_full - s_masked(random N), averaged over 3 draws

Interpretation (decided AFTER running, not assumed):
  |joint_delta| >> |sum_individual| AND |joint_delta| >> |random_delta|
      -> real redundant/compensatory encoding: individual genes carry real
         signal that's masked by compensation from co-expressed genes, so
         single-gene masking underestimates true importance. The null GSEA
         would then reflect a real methodological limitation of single-gene
         occlusion for THIS architecture, not an absence of biology.
  |joint_delta| ~ |sum_individual| and both stay small even at N=30, AND
  |random_delta| is comparably small
      -> generic insensitivity: even removing 30 genes at once barely moves
         the score. The classifier's decision direction is just not resolved
         by this occlusion granularity for this model -- caveat the null
         GSEA heavily as inconclusive rather than "no biology."
  |joint_delta| ~ |random_delta| (topN no more impactful than matched
  random-N, even jointly)
      -> the single-gene ranking itself carries no reliable signal at any
         granularity tested -- reinforces (via an independent test) that the
         null GSEA result is likely a genuine flat/uninformative attribution
         profile for this model, not just an artifact of testing one gene at
         a time.

Reuses the exact same cohort loading, gene universe, and candidate-selection
machinery as scprint2_masking_attribution.py (masking_attribution_common.py)
so the candidate genes available for masking here are IDENTICAL, cell for
cell, to the ones used in the original single-gene sweep.
"""
import os
import re
import sys
import time

import numpy as np
import pandas as pd
import torch
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from masking_attribution_common import (
    load_canonical_cohort, shared_gene_universe_5way, top_k_genes_per_cell,
    TOP_K_GENES,
)

from scprint2 import scPRINT2

CKPT_PATH = "/bigdata/godziklab/shared/Xinru/.huggingface/scprint2_ckpts/small-v2.ckpt"
OUT_DIR = "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/gene_attribution_gsea/scprint2_masking"
RANKED_CSV = f"{OUT_DIR}/scprint2_masking_gene_attribution_ranked.csv"
CTRL_OUT_DIR = "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/gene_attribution_gsea/scprint2_multigene_occlusion"
os.makedirs(CTRL_OUT_DIR, exist_ok=True)

HUMAN_GENE_RE = re.compile(r"^ENSG\d+$")
N_TIERS = [10, 30, 128]
N_RANDOM_DRAWS = 3
RNG_SEED_BASE = 20260911

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Device: {device}", flush=True)

_test_max_cells = os.environ.get("MASKING_TEST_MAX_CELLS")
_test_max_cells = int(_test_max_cells) if _test_max_cells else None
if _test_max_cells:
    print(f"[TESTING] MASKING_TEST_MAX_CELLS={_test_max_cells}", flush=True)
    CTRL_OUT_DIR = CTRL_OUT_DIR + "_TEST"
    os.makedirs(CTRL_OUT_DIR, exist_ok=True)

a, y_by_patient = load_canonical_cohort(max_cells=_test_max_cells)
var_names = a.var_names.astype(str).to_numpy()

print(f"Loading scPRINT-2 checkpoint: {CKPT_PATH} (expect ~7-10 min)...", flush=True)
t_load = time.time()
model = scPRINT2.load_from_checkpoint(CKPT_PATH, precpt_gene_emb=None, gene_pos_file=None)
print(f"Loaded in {(time.time() - t_load) / 60:.1f} min.", flush=True)

genes_list = list(model.genes)
human_ensembl_ids = {g for g in genes_list if HUMAN_GENE_RE.match(g)}
ensembl_to_global_idx = {g: i for i, g in enumerate(genes_list) if HUMAN_GENE_RE.match(g)}

valid_mask, sym2ens, _token_dict, _protein_embeds, _scgpt_vocab, _uce_pe, _chrom_df = \
    shared_gene_universe_5way(var_names, human_ensembl_ids)

if not torch.cuda.is_available():
    model = model.to(torch.float32)
model = model.to(device)
model.eval()

print(f"Selecting top-{TOP_K_GENES} candidate genes per cell...", flush=True)
top_k_idx = top_k_genes_per_cell(a.X, valid_mask, k=TOP_K_GENES)

col_to_global_idx = {}
for i, gene in enumerate(var_names):
    if valid_mask[i]:
        ens = sym2ens.get(gene)
        gi = ensembl_to_global_idx.get(ens)
        if gi is not None:
            col_to_global_idx[i] = gi

print(f"Loading existing single-gene ranked list: {RANKED_CSV}", flush=True)
ranked = pd.read_csv(RANKED_CSV)
gene_to_attr = dict(zip(ranked["gene"], ranked["mean_attribution"]))
gene_symbol_by_col = {i: g for i, g in enumerate(var_names)}
print(f"  {len(gene_to_attr)} genes have a single-gene attribution value "
      f"(>= 30 occurrences in the original sweep)", flush=True)

X = a.X


def cell_total_counts(cell_idx):
    row = X[cell_idx].toarray().ravel() if hasattr(X[cell_idx], "toarray") else np.asarray(X[cell_idx]).ravel()
    return float(np.clip(row, 0, None).sum())


def build_cell_tensors(cell_idx, gene_cols):
    row = X[cell_idx].toarray().ravel() if hasattr(X[cell_idx], "toarray") else np.asarray(X[cell_idx]).ravel()
    gene_pos, expr, cols = [], [], []
    for col in gene_cols:
        gi = col_to_global_idx.get(col)
        if gi is None:
            continue
        gene_pos.append(gi)
        expr.append(float(row[col]))
        cols.append(col)
    return (np.array(gene_pos, dtype=np.int64), np.array(expr, dtype=np.float32), cols)


n_cells = a.n_obs
y_cell = a.obs["y"].to_numpy()
patient_cell = a.obs["patient_id"].to_numpy()

# ---------------------------------------------------------------------
# Pass 1: refit the SAME classifier as the single-gene run (same cohort,
# same procedure) -- needed to get w, b for scoring. Not re-saved; the
# single-gene run's classifier_diagnostics.csv remains authoritative.
# ---------------------------------------------------------------------
BATCH_CELLS = 16
pooled_list = []
t0 = time.time()
with torch.no_grad():
    for start in range(0, n_cells, BATCH_CELLS):
        end = min(start + BATCH_CELLS, n_cells)
        gp_list, ex_list, depth_list = [], [], []
        for ci in range(start, end):
            gp, ex, _ = build_cell_tensors(ci, top_k_idx[ci])
            gp_list.append(gp)
            ex_list.append(ex)
            depth_list.append(cell_total_counts(ci))
        max_len = max(len(g) for g in gp_list)
        gp_2d = np.zeros((end - start, max_len), dtype=np.int64)
        ex_2d = np.zeros((end - start, max_len), dtype=np.float32)
        mz_2d = np.zeros((end - start, max_len), dtype=bool)
        for j, (gp, ex) in enumerate(zip(gp_list, ex_list)):
            gp_2d[j, :len(gp)] = gp
            ex_2d[j, :len(ex)] = ex
            mz_2d[j, :len(gp)] = True

        gp_t = torch.tensor(gp_2d, device=device)
        ex_t = torch.tensor(ex_2d, device=device)
        mz_t = torch.tensor(mz_2d, device=device)
        depth_t = torch.tensor(depth_list, dtype=torch.float32, device=device)

        out = model.forward(
            gene_pos=gp_t, expression=ex_t, mask_zeros=mz_t,
            req_depth=depth_t, do_class=False, do_mvc=False,
        )
        emb = out["output_cell_emb"]
        pooled_list.append(emb.detach().cpu().float().numpy())
        if (start // BATCH_CELLS) % 30 == 0:
            print(f"  pass1 embedding batch {end}/{n_cells} "
                  f"({(time.time() - t0) / 60:.1f} min)", flush=True)

pooled_all = np.concatenate(pooled_list, axis=0)
clf = LogisticRegression(penalty="l2", C=0.3, max_iter=2000)
clf.fit(pooled_all, y_cell)
cell_auc = roc_auc_score(y_cell, clf.predict_proba(pooled_all)[:, 1])
print(f"Refit classifier cell-level AUC = {cell_auc:.4f} "
      f"(sanity check -- should match the single-gene run's 0.857)", flush=True)

w = clf.coef_[0].astype(np.float32)
b = float(clf.intercept_[0])
w_t = torch.tensor(w, device=device)

# ---------------------------------------------------------------------
# Pass 2: joint-occlusion sweep.
# ---------------------------------------------------------------------
records = []
t0 = time.time()
with torch.no_grad():
    for i in range(n_cells):
        gp, ex, cols = build_cell_tensors(i, top_k_idx[i])
        n = len(cols)
        if n < max(N_TIERS):
            continue  # too few candidates even for the largest tier alone
        depth = cell_total_counts(i)
        cols_arr = np.array(cols)
        rng = np.random.default_rng(RNG_SEED_BASE + i)

        # rank this cell's candidates by existing single-gene attribution
        attrs = np.array([gene_to_attr.get(gene_symbol_by_col[c], np.nan) for c in cols_arr])
        has_attr = ~np.isnan(attrs)
        cand_with_attr = np.where(has_attr)[0]  # positions within cols_arr
        order_desc = cand_with_attr[np.argsort(-attrs[cand_with_attr], kind="stable")]
        order_asc = cand_with_attr[np.argsort(attrs[cand_with_attr], kind="stable")]

        row_idx = 1  # row 0 reserved for s_full
        row_meta = []  # (n_tier, arm, draw_i)
        mask_rows = []  # list of (positions-within-n to mask)

        cell_record = {"cell_idx": i}
        skip_cell = False
        for N in N_TIERS:
            if len(order_desc) < N or len(order_asc) < N:
                skip_cell = True
                break
            pos_positions = order_desc[:N]
            neg_positions = order_asc[:N]
            excluded = set(pos_positions.tolist()) | set(neg_positions.tolist())
            remaining = np.array([p for p in range(n) if p not in excluded])
            # At large N (e.g. 128 out of ~226 candidates), a disjoint
            # remainder of size >= N often doesn't exist (2*N can exceed n).
            # Fall back to sampling from the FULL candidate set (allowing
            # overlap with top_pos/top_neg) rather than skipping the cell --
            # the random-N control's purpose is "what does removing ANY N
            # genes do to the score", which doesn't require disjointness;
            # disjointness was only a (nice-to-have, now infeasible at this
            # N) way to avoid the random draw accidentally reusing the same
            # genes at the smaller tiers, where it was always achievable.
            if len(remaining) >= N:
                random_pool = remaining
                pool_is_disjoint = True
            else:
                random_pool = np.arange(n)
                pool_is_disjoint = False
            if len(random_pool) < N:
                skip_cell = True
                break

            mask_rows.append(pos_positions)
            row_meta.append((N, "top_pos", None))
            mask_rows.append(neg_positions)
            row_meta.append((N, "top_neg", None))
            for d in range(N_RANDOM_DRAWS):
                rand_positions = rng.choice(random_pool, size=N, replace=False)
                mask_rows.append(rand_positions)
                row_meta.append((N, "random", d))

            cell_record[f"sum_individual_topN_pos_N{N}"] = float(attrs[pos_positions].sum())
            cell_record[f"sum_individual_topN_neg_N{N}"] = float(attrs[neg_positions].sum())
            cell_record[f"random_pool_disjoint_N{N}"] = pool_is_disjoint

        if skip_cell:
            continue

        Bsz = 1 + len(mask_rows)
        gp_batch = np.tile(gp, (Bsz, 1))
        ex_batch = np.tile(ex, (Bsz, 1))
        mz_batch = np.ones((Bsz, n), dtype=bool)
        for r, positions in enumerate(mask_rows):
            mz_batch[r + 1, positions] = False

        gp_t = torch.tensor(gp_batch, device=device)
        ex_t = torch.tensor(ex_batch, device=device)
        mz_t = torch.tensor(mz_batch, device=device)
        depth_t = torch.full((Bsz,), depth, dtype=torch.float32, device=device)

        out = model.forward(
            gene_pos=gp_t, expression=ex_t, mask_zeros=mz_t,
            req_depth=depth_t, do_class=False, do_mvc=False,
        )
        emb = out["output_cell_emb"]
        logits = (emb @ w_t + b).cpu().numpy()
        s_full = logits[0]
        s_masked = logits[1:]
        deltas = s_full - s_masked

        rand_accum = {}
        for (N, arm, draw), d in zip(row_meta, deltas):
            if arm == "random":
                rand_accum.setdefault(N, []).append(d)
            else:
                cell_record[f"joint_delta_{arm}_N{N}"] = float(d)
        for N, vals in rand_accum.items():
            cell_record[f"joint_delta_random_N{N}"] = float(np.mean(vals))

        records.append(cell_record)

        if (i + 1) % 500 == 0:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed
            eta_min = (n_cells - i - 1) / rate / 60 if rate > 0 else float("nan")
            print(f"  occlusion sweep {i + 1}/{n_cells} "
                  f"({rate:.2f} cells/s, ETA {eta_min:.1f} min, "
                  f"{len(records)} cells with enough candidates)", flush=True)

print(f"Joint-occlusion sweep done in {(time.time() - t0) / 60:.1f} min "
      f"over {len(records)} / {n_cells} cells (rest skipped: too few candidates)", flush=True)

df = pd.DataFrame(records)
per_cell_path = f"{CTRL_OUT_DIR}/scprint2_multigene_occlusion_per_cell.csv"
df.to_csv(per_cell_path, index=False)
print(f"Saved per-cell results: {per_cell_path}", flush=True)

print("\n=== SUMMARY (mean +/- SD across cells) ===", flush=True)
summary_rows = []
for N in N_TIERS:
    for polarity, joint_col, sum_col in [
        ("positive", f"joint_delta_top_pos_N{N}", f"sum_individual_topN_pos_N{N}"),
        ("negative", f"joint_delta_top_neg_N{N}", f"sum_individual_topN_neg_N{N}"),
    ]:
        joint = df[joint_col]
        summed = df[sum_col]
        rand = df[f"joint_delta_random_N{N}"]
        ratio_vs_sum = joint.mean() / summed.mean() if summed.mean() != 0 else float("nan")
        row = {
            "N": N, "polarity": polarity,
            "joint_delta_mean": joint.mean(), "joint_delta_sd": joint.std(),
            "sum_individual_mean": summed.mean(),
            "random_delta_mean": rand.mean(), "random_delta_sd": rand.std(),
            "ratio_joint_over_sum": ratio_vs_sum,
            "ratio_joint_over_random_abs": abs(joint.mean()) / abs(rand.mean()) if rand.mean() != 0 else float("nan"),
            "n_cells": len(df),
        }
        summary_rows.append(row)
        print(f"N={N:2d} {polarity:8s}: joint={joint.mean():+.4f} (sd={joint.std():.4f})  "
              f"sum_individual={summed.mean():+.4f}  random={rand.mean():+.4f} (sd={rand.std():.4f})  "
              f"joint/sum={ratio_vs_sum:+.2f}x  |joint|/|random|={row['ratio_joint_over_random_abs']:.2f}x",
              flush=True)

summary_df = pd.DataFrame(summary_rows)
summary_path = f"{CTRL_OUT_DIR}/scprint2_multigene_occlusion_summary.csv"
summary_df.to_csv(summary_path, index=False)
print(f"\nSaved summary: {summary_path}", flush=True)
