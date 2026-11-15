#!/usr/bin/env python3
"""
scPRINT-2 masking-based ("leave-one-gene-out") gene attribution of the
mortality-associated direction (canonical 240-patient/25-fatal within-COVID
cohort). Fifth and final model in the unified masking design, after
Geneformer, STATE, scGPT, and UCE -- see masking_attribution_common.py's
docstring for the shared rationale.

Architecture facts this script depends on, all verified directly against
the ACTUAL environment this project uses in production
(scripts/foundation_models/scPRINT/run_scprint2_human_covid.sh):
`scprint2_env`, package literally named `scprint2` (checkpoint
small-v2.ckpt) -- NOT the older, now-archived `scprint_env` (package
`scprint==2.3.8`) whose file paths superficially matched an earlier
external code review's citations and caused a brief, caught-before-any-
code-was-written misidentification (see project memory for the full
course-correction).

- `type(model.transformer)` is `simpler_flash.flashformer.FlashTransformer`
  (confirmed by loading the actual checkpoint), NOT `Performer` -- the
  "FlashAttention is not installed" warning at import time refers only to
  an optional acceleration backend WITHIN simpler_flash, not to
  FlashTransformer's own availability.
- Unlike Geneformer/UCE, masking here needs NO bypass of forward() at all:
  `scPRINT2.forward(..., mask_zeros=...)` genuinely threads `mask_zeros`
  through to `self.transformer(..., mask_zeros=mask_zeros)`
  (model.py:~1378), and scprint2_env's own simpler_flash/flashformer.py
  (line 224) passes it straight to `src_key_padding_mask` in the
  underlying attention call. Verified polarity directly in THIS package's
  mha.py (not assumed from the older package): "True means to keep, False
  means to mask out" (mha.py:810).
- `model.cell_transformer` is present for this checkpoint (a separate
  6-layer FlashTransformer, d_model=256, that cross-attends cell tokens
  against the main transformer's output AFTER it runs) -- confirmed via
  `forward()`'s own `num = (...) if not self.cell_transformer else 0`
  logic, which evaluates to num=0 here. This means the main transformer
  processes ONLY gene positions (no prepended cell-token block), so
  mask_zeros only needs shape (batch, seq_len) matching gene_pos exactly
  -- simpler than STATE/UCE/scGPT's "+special tokens" mask shapes.
- `model.expr_encoder` is `EasyExprGNN` (not `ExprBasedFT`) for this
  checkpoint, and `model.normalization == "log"`: `_encoder()` applies a
  simple `log2(1+x)` transform to raw expression, independently per gene
  -- unlike "sum" normalization (used by the OLD scprint package and by
  scGPT/scPRINT's own weaker per-gene binning concerns), a log transform
  has NO cross-gene coupling at all, so there is no renormalization
  confound to guard against here even in principle -- not that it matters
  since this script uses flag-based masking (input tensors never change
  across masking variants), not deletion.
- `model.genes` is a 335,182-entry FLAT multi-species Ensembl-gene-ID
  vocabulary (confirmed by loading the checkpoint). The human
  (`^ENSG\\d+$` -- NOT a naive `startswith("ENSG")`, which also matches
  e.g. chicken's `ENSGALG...` prefix, a real bug caught during
  verification) subset is exactly 20,004 genes, occupying one CONTIGUOUS
  block (global indices 246612-266615) -- confirmed empirically, not
  assumed.
- `model.pos_encoder` is an active `PositionalEncoding()` module (unlike
  STATE, which has none) -- but it is applied as `enc.add_(self.
  pos_encoder(gene_pos))`, a function of gene IDENTITY, not of sequence
  POSITION within our candidate list -- so candidate gene order within
  the built sentence is immaterial to this term, and masking one gene
  does not disturb another gene's positional contribution.
- Checkpoint loading: `scPRINT2.load_from_checkpoint(ckpt_path,
  precpt_gene_emb=None, gene_pos_file=None)` (the exact call scprint2's
  own CLI uses, `scprint2/cli.py`). Loading takes ~7-10 minutes on this
  cluster (confirmed genuinely CPU-busy via `sstat`, not hung) -- budget
  for this when sizing SLURM walltime.
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
    aggregate_and_save, save_checkpoint_atomic, TOP_K_GENES, MIN_GENE_OCCURRENCE,
)

from scprint2 import scPRINT2

CKPT_PATH = "/bigdata/godziklab/shared/Xinru/.huggingface/scprint2_ckpts/small-v2.ckpt"
OUT_DIR = "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/gene_attribution_gsea/scprint2_masking"
os.makedirs(OUT_DIR, exist_ok=True)

HUMAN_GENE_RE = re.compile(r"^ENSG\d+$")

CHECKPOINT_PATH = f"{OUT_DIR}/masking_sweep_checkpoint.npz"
TIME_BUDGET_SECONDS = float(os.environ.get("TIME_BUDGET_SECONDS", 5400))
_script_start_time = time.time()

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Device: {device}", flush=True)

_test_max_cells = os.environ.get("MASKING_TEST_MAX_CELLS")
_test_max_cells = int(_test_max_cells) if _test_max_cells else None
if _test_max_cells:
    print(f"[TESTING] MASKING_TEST_MAX_CELLS={_test_max_cells}", flush=True)
    OUT_DIR = OUT_DIR + "_TEST"
    os.makedirs(OUT_DIR, exist_ok=True)

a, y_by_patient = load_canonical_cohort(max_cells=_test_max_cells)
var_names = a.var_names.astype(str).to_numpy()

# ---------------------------------------------------------------------
# Load model FIRST (the expensive ~7-10 min step): the shared gene
# universe needs model.genes, which is only obtainable this way for
# scPRINT-2 (no separate lightweight vocab file ships with it).
# ---------------------------------------------------------------------
print(f"Loading scPRINT-2 checkpoint: {CKPT_PATH} (expect ~7-10 min)...", flush=True)
t_load = time.time()
model = scPRINT2.load_from_checkpoint(CKPT_PATH, precpt_gene_emb=None, gene_pos_file=None)
print(f"Loaded in {(time.time() - t_load) / 60:.1f} min.", flush=True)

genes_list = list(model.genes)
human_ensembl_ids = {g for g in genes_list if HUMAN_GENE_RE.match(g)}
print(f"scPRINT-2 human gene vocabulary: {len(human_ensembl_ids)} genes", flush=True)
ensembl_to_global_idx = {g: i for i, g in enumerate(genes_list) if HUMAN_GENE_RE.match(g)}

valid_mask, sym2ens, _token_dict, _protein_embeds, _scgpt_vocab, _uce_pe, _chrom_df = \
    shared_gene_universe_5way(var_names, human_ensembl_ids)

if not torch.cuda.is_available():
    model = model.to(torch.float32)
model = model.to(device)
model.eval()

print(f"Selecting top-{TOP_K_GENES} candidate genes per cell...", flush=True)
top_k_idx = top_k_genes_per_cell(a.X, valid_mask, k=TOP_K_GENES)
n_per_cell = np.array([len(x) for x in top_k_idx])
print(f"Candidates per cell: min={n_per_cell.min()}, max={n_per_cell.max()}, "
      f"mean={n_per_cell.mean():.1f} (target {TOP_K_GENES})", flush=True)

col_to_global_idx = {}
for i, gene in enumerate(var_names):
    if valid_mask[i]:
        ens = sym2ens.get(gene)
        gi = ensembl_to_global_idx.get(ens)
        if gi is not None:
            col_to_global_idx[i] = gi
print(f"col_to_global_idx built for {len(col_to_global_idx)} / {int(valid_mask.sum())} "
      f"shared-universe genes", flush=True)

X = a.X


def cell_total_counts(cell_idx):
    """Full-transcriptome raw count total, used as the req_depth proxy
    (matches the "sequencing depth" role documented for this parameter)."""
    row = X[cell_idx].toarray().ravel() if hasattr(X[cell_idx], "toarray") else np.asarray(X[cell_idx]).ravel()
    return float(np.clip(row, 0, None).sum())


def build_cell_tensors(cell_idx, gene_cols):
    """Returns (gene_pos int64 array, expression float32 array, cols list)
    for the candidate genes present (nonzero, by construction of
    top_k_genes_per_cell) in this cell -- length n, no padding, no
    prepended special tokens (see module docstring: cell_transformer
    handles cell-level pooling separately)."""
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
# Pass 1: s_full embeddings (no grad) -> fit mortality classifier
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
        mz_2d = np.zeros((end - start, max_len), dtype=bool)  # True=keep
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
print(f"Pooled embeddings: {pooled_all.shape}", flush=True)

clf = LogisticRegression(penalty="l2", C=0.3, max_iter=2000)
clf.fit(pooled_all, y_cell)
cell_auc = roc_auc_score(y_cell, clf.predict_proba(pooled_all)[:, 1])

patient_score = pd.DataFrame({"patient_id": patient_cell, "p": clf.predict_proba(pooled_all)[:, 1]}) \
    .groupby("patient_id")["p"].mean()
patient_y = pd.Series(y_by_patient)
patient_auc = roc_auc_score(patient_y.loc[patient_score.index], patient_score)
print(f"Cell-level in-sample AUC = {cell_auc:.4f}; "
      f"patient-level (mean-pooled) in-sample AUC = {patient_auc:.4f}", flush=True)

w = clf.coef_[0].astype(np.float32)
b = float(clf.intercept_[0])
w_t = torch.tensor(w, device=device)

pd.DataFrame([{"model": "scprint2_masking", "cell_auc": cell_auc, "patient_auc": patient_auc,
               "n_cells": n_cells, "n_patients": a.obs["patient_id"].nunique(),
               "top_k_genes": TOP_K_GENES, "embedding_dim": pooled_all.shape[1]}]) \
    .to_csv(f"{OUT_DIR}/classifier_diagnostics.csv", index=False)

# ---------------------------------------------------------------------
# Pass 2: masking sweep -- flip mask_zeros False at one gene position per
# variant (all variants share the same sequence length, like scGPT --
# input tensors never change, only mask_zeros, so no renormalization
# concern even in principle given log-transform normalization).
# ---------------------------------------------------------------------
n_genes_total = len(var_names)
start_cell = 0
if os.path.exists(CHECKPOINT_PATH):
    _ckpt = np.load(CHECKPOINT_PATH)
    gene_sum = _ckpt["gene_sum"]
    gene_sumsq = _ckpt["gene_sumsq"]
    gene_count = _ckpt["gene_count"]
    start_cell = int(_ckpt["next_cell"])
    print(f"Resuming masking sweep from checkpoint: {start_cell}/{n_cells} cells "
          f"already done", flush=True)
else:
    gene_sum = np.zeros(n_genes_total, dtype=np.float64)
    gene_sumsq = np.zeros(n_genes_total, dtype=np.float64)
    gene_count = np.zeros(n_genes_total, dtype=np.int64)

t0 = time.time()
with torch.no_grad():
    for i in range(start_cell, n_cells):
        gp, ex, cols = build_cell_tensors(i, top_k_idx[i])
        n = len(cols)
        if n == 0:
            continue
        depth = cell_total_counts(i)

        Bsz = n + 1
        gp_batch = np.tile(gp, (Bsz, 1))
        ex_batch = np.tile(ex, (Bsz, 1))
        mz_batch = np.ones((Bsz, n), dtype=bool)  # True = keep
        for k in range(n):
            mz_batch[k + 1, k] = False  # mask out candidate k

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
        delta = s_full - s_masked
        cols_arr = np.array(cols)

        np.add.at(gene_sum, cols_arr, delta)
        np.add.at(gene_sumsq, cols_arr, delta ** 2)
        np.add.at(gene_count, cols_arr, 1)

        if (i + 1) % 200 == 0:
            elapsed = time.time() - t0
            rate = (i + 1 - start_cell) / elapsed
            eta_min = (n_cells - i - 1) / rate / 60 if rate > 0 else float("nan")
            print(f"  masking sweep {i + 1}/{n_cells} "
                  f"({rate:.2f} cells/s, ETA {eta_min:.1f} min)", flush=True)

            if time.time() - _script_start_time > TIME_BUDGET_SECONDS:
                save_checkpoint_atomic(CHECKPOINT_PATH, gene_sum=gene_sum,
                                        gene_sumsq=gene_sumsq, gene_count=gene_count,
                                        next_cell=i + 1)
                print(f"Time budget ({TIME_BUDGET_SECONDS:.0f}s) reached at cell "
                      f"{i + 1}/{n_cells}. Checkpoint saved to {CHECKPOINT_PATH}. "
                      f"Exiting cleanly for the sbatch wrapper to resubmit and resume.",
                      flush=True)
                sys.exit(0)

print(f"Masking sweep done in {(time.time() - t0) / 60:.1f} min", flush=True)

out_path = f"{OUT_DIR}/scprint2_masking_gene_attribution_ranked.csv"
aggregate_and_save(gene_sum, gene_sumsq, gene_count, var_names, out_path, MIN_GENE_OCCURRENCE)
if os.path.exists(CHECKPOINT_PATH):
    os.remove(CHECKPOINT_PATH)
    print(f"Removed checkpoint file (run complete).", flush=True)
