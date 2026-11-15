#!/usr/bin/env python3
"""
STATE (Arc Institute SE-600M) masking-based ("leave-one-gene-out") gene
attribution of the mortality-associated direction (canonical
240-patient/25-fatal within-COVID cohort). Companion to
geneformer_masking_attribution.py -- both import
masking_attribution_common.py so the candidate gene SET selected per cell
is byte-for-byte identical between the two models (see that module's
docstring for the full rationale).

Architecture facts this script depends on (all verified by reading
state_env's arc-state==0.9.27 source directly, not assumed):
- STATE's native sentence construction (VCIDatasetSentenceCollator.
  sample_cell_sentences) is, for a cell with more expressed genes than
  pad_length, effectively DETERMINISTIC top-(pad_length-1) by expression
  (random only for tie-breaking and for padding with unexpressed genes
  when a cell has fewer expressed genes than pad_length -- not the case
  for platelets). This script reimplements that same "rank by expression,
  take top-K" logic directly, restricted to the K=TOP_K_GENES candidate
  set shared with Geneformer, rather than STATE's native pad_length=2048.
- Crucially, unlike Geneformer's BertModel, STATE's transformer_encoder is
  called with `src_key_padding_mask=None` (state.emb.nn.model.
  StateEmbeddingModel.forward) -- the `mask` tensor is NOT used to exclude
  positions from attention, and pooling is CLS-token pooling (not
  mean-pooling), so a mask flag alone cannot "remove" a gene the way it
  does for Geneformer. Masking a gene here means literally DELETING its
  position from the sequence (built as a per-cell variable-length
  sentence, no padding needed since each cell's masking sweep is
  processed as its own batch of same-length variants: n_valid-length
  after any single deletion).
- pe_embedding does NOT need to be unfrozen for this method (no gradients
  used at all -- forward passes only).
- `dataset_nums` is accepted by StateEmbeddingModel.forward() but never
  referenced inside its body (verified via grep) -- only self.dataset_token
  (a single learned parameter, not per-dataset) affects the `embedding`
  output this script uses, so passing dataset_nums=None is safe.
- Checkpoint selection mirrors the project's own production invocation
  (glob.glob(model_folder + "/*.ckpt")[-1] -> se600m_epoch4.ckpt),
  identical to state_gene_attribution.py (the now-superseded gradient
  script) and to how state_covid_canonical.h5ad was originally built.
"""
import glob
import os
import sys
import time

import numpy as np
import pandas as pd
import torch

from state.emb.inference import Inference

# PyTorch >=2.6 defaults torch.load(weights_only=True), which rejects this
# checkpoint's embedded omegaconf.DictConfig during PyTorch Lightning's own
# load_from_checkpoint -> torch.load call chain (not something Inference's
# own API exposes a way to override). The checkpoint is Arc Institute's own
# first-party release (trusted source), so restore the pre-2.6 default for
# the duration of this script, per the error message's own recommendation.
_orig_torch_load = torch.load


def _patched_torch_load(*args, **kwargs):
    # force-override, not setdefault: Lightning's own pl_load passes
    # weights_only explicitly (even if None), so setdefault is a no-op --
    # confirmed by a first attempt that still hit the same UnpicklingError
    # with setdefault in place.
    kwargs["weights_only"] = False
    return _orig_torch_load(*args, **kwargs)


torch.load = _patched_torch_load

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from masking_attribution_common import (
    load_canonical_cohort, shared_gene_universe, top_k_genes_per_cell,
    aggregate_and_save, save_checkpoint_atomic, TOP_K_GENES, MIN_GENE_OCCURRENCE,
)

MODEL_FOLDER = "/bigdata/godziklab/shared/Xinru/302006/03_MODELS/SE-600M"
OUT_DIR = "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/gene_attribution_gsea/state_masking"
os.makedirs(OUT_DIR, exist_ok=True)
CLS_TOKEN_IDX_PLACEHOLDER = 3  # config.yaml's cls_token_idx; value is irrelevant since the
                               # embedding at this position gets overwritten by the learned
                               # CLS parameter regardless of what integer sits here

# Checkpoint/resume support, for chaining through short_gpu's 2h walltime cap
# instead of waiting in the (often heavily congested) gpu partition's queue
# for one uninterrupted multi-hour allocation. Pass1 (embeddings + classifier
# fit) is cheap (a few minutes) and deterministic given the same cohort/model
# -- so only pass2's masking-sweep accumulators are checkpointed; pass1 is
# simply redone on every restart rather than also being saved/restored.
# TIME_BUDGET_SECONDS should be set (via env var, by the sbatch wrapper) to
# comfortably less than the SLURM --time limit, leaving margin for pass1 and
# final I/O. Default here (5400s = 1h30m) assumes a 2h short_gpu allocation.
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
valid_mask, _sym2ens, _token_dict, protein_embeds = shared_gene_universe(var_names)

gene_symbols_by_global_idx = list(protein_embeds.keys())
global_idx_by_symbol = {s: i for i, s in enumerate(gene_symbols_by_global_idx)}
col_to_global = {}
for i, gene in enumerate(var_names):
    if valid_mask[i]:
        gi = global_idx_by_symbol.get(gene)
        if gi is not None:
            col_to_global[i] = gi
print(f"col_to_global built for {len(col_to_global)} / {int(valid_mask.sum())} "
      f"shared-universe genes", flush=True)

print(f"Selecting top-{TOP_K_GENES} candidate genes per cell...", flush=True)
top_k_idx = top_k_genes_per_cell(a.X, valid_mask, k=TOP_K_GENES)
n_per_cell = np.array([len(x) for x in top_k_idx])
print(f"Candidates per cell: min={n_per_cell.min()}, max={n_per_cell.max()}, "
      f"mean={n_per_cell.mean():.1f} (target {TOP_K_GENES}); "
      f"cells below target: {int((n_per_cell < TOP_K_GENES).sum())} / {len(n_per_cell)}", flush=True)

ckpt_files = glob.glob(os.path.join(MODEL_FOLDER, "*.ckpt"))
assert ckpt_files, f"no checkpoint found in {MODEL_FOLDER}"
checkpoint = ckpt_files[-1]
print(f"Using checkpoint: {checkpoint}", flush=True)

inferer = Inference(cfg=None, protein_embeds=protein_embeds)
inferer.load_model(checkpoint)
model = inferer.model
# Inference.load_model() casts the model to a mixed-precision dtype (bf16 on
# CUDA, via get_precision_config) for its own autocast-based encode() path.
# This script builds its own fp32 ids/counts tensors directly (no autocast),
# so force the model back to plain fp32 to avoid a dtype mismatch inside
# count_encoder's Linear layer (confirmed via a real crash: "mat1 and mat2
# must have the same dtype, but got Float and BFloat16").
model = model.float()
model.eval()
model.to(device)

X = a.X


def cell_expression_weights(row_idx):
    """log1p-normalized expression distribution over the FULL transcriptome
    for one cell -- matches VCIDatasetSentenceCollator.sample_cell_sentences's
    own log1p-then-normalize convention for the count channel (verified by
    reading its source: raw integer counts -> log1p -> divide by per-cell
    total)."""
    row = X[row_idx].toarray().ravel() if hasattr(X[row_idx], "toarray") else np.asarray(X[row_idx]).ravel()
    row = np.clip(row, 0, None)
    log1p = np.log1p(row)
    total = log1p.sum()
    if total <= 0:
        return np.zeros_like(log1p)
    return log1p / total


def build_full_sentence(cell_idx, gene_cols, ew):
    n = len(gene_cols)
    ids = np.empty(n + 1, dtype=np.int64)
    counts = np.empty(n + 1, dtype=np.float32)
    ids[0] = CLS_TOKEN_IDX_PLACEHOLDER
    counts[0] = 0.0
    for j, col in enumerate(gene_cols):
        ids[j + 1] = col_to_global[col]
        counts[j + 1] = 100.0 * ew[col]
    return ids, counts


def _dummy_batch(ids_2d, counts_2d):
    """Build the 9-tuple state._compute_embedding_for_batch expects.
    Only positions 0 (sentences), 5 (mask), 7 (counts) matter for the
    `embedding` output we use; batch[1] (task/query sentence 'X') must be a
    real tensor (unconditionally .to(device)'d) but its own output is
    never consumed by us, so a minimal valid dummy is used."""
    B, L = ids_2d.shape
    sentences = torch.tensor(ids_2d, dtype=torch.int64, device=device)
    mask = torch.ones((B, L), dtype=torch.bool, device=device)
    counts_t = torch.tensor(counts_2d, dtype=torch.float32, device=device)
    dummy_X = torch.zeros((B, 1), dtype=torch.int64, device=device)
    dummy_Y = torch.zeros((B, 1), dtype=torch.float32, device=device)
    dummy_weights = torch.zeros((B, L), dtype=torch.float32, device=device)
    return (sentences, dummy_X, dummy_Y, None, dummy_weights, mask, None, counts_t, None)


n_cells = a.n_obs
y_cell = a.obs["y"].to_numpy()
patient_cell = a.obs["patient_id"].to_numpy()

# ---------------------------------------------------------------------
# Pass 1: s_full embeddings (no grad) -> fit mortality classifier
# ---------------------------------------------------------------------
embedding_dim = None
pooled_all = None
BATCH_CELLS = 16  # cells sharing the same sequence length can't generally be
                   # batched together (different n per cell), so pass 1 pads
                   # to the per-batch max length instead
t0 = time.time()
pooled_list = []
with torch.no_grad():
    for start in range(0, n_cells, BATCH_CELLS):
        end = min(start + BATCH_CELLS, n_cells)
        cell_idxs = list(range(start, end))
        sentences_list, counts_list = [], []
        for ci in cell_idxs:
            ew = cell_expression_weights(ci)
            ids, counts = build_full_sentence(ci, top_k_idx[ci], ew)
            sentences_list.append(ids)
            counts_list.append(counts)
        max_len = max(len(s) for s in sentences_list)
        ids_2d = np.zeros((len(cell_idxs), max_len), dtype=np.int64)
        counts_2d = np.zeros((len(cell_idxs), max_len), dtype=np.float32)
        for j, (ids, counts) in enumerate(zip(sentences_list, counts_list)):
            ids_2d[j, :len(ids)] = ids
            counts_2d[j, :len(counts)] = counts
        batch = _dummy_batch(ids_2d, counts_2d)
        _, _, _, embedding, _ = model._compute_embedding_for_batch(batch)
        pooled_list.append(embedding.detach().cpu().float().numpy())
        if (start // BATCH_CELLS) % 30 == 0:
            print(f"  pass1 embedding batch {end}/{n_cells} "
                  f"({(time.time()-t0)/60:.1f} min elapsed)", flush=True)

pooled_all = np.concatenate(pooled_list, axis=0)
embedding_dim = pooled_all.shape[1]
print(f"Pooled embeddings: {pooled_all.shape}", flush=True)

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

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

pd.DataFrame([{"model": "state_masking", "cell_auc": cell_auc, "patient_auc": patient_auc,
               "n_cells": n_cells, "n_patients": a.obs["patient_id"].nunique(),
               "top_k_genes": TOP_K_GENES, "embedding_dim": embedding_dim}]) \
    .to_csv(f"{OUT_DIR}/classifier_diagnostics.csv", index=False)

# ---------------------------------------------------------------------
# Pass 2: masking sweep -- per cell, one full-sentence forward pass plus
# one batched forward pass over all single-gene-deletion variants (all
# sharing the same shorter length n, so they batch cleanly together;
# separate from the length-(n+1) full-sentence call).
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
        gene_cols = top_k_idx[i]
        n = len(gene_cols)
        if n == 0:
            continue
        ew = cell_expression_weights(i)
        ids_full, counts_full = build_full_sentence(i, gene_cols, ew)

        # s_full
        batch_full = _dummy_batch(ids_full[None, :], counts_full[None, :])
        _, _, _, emb_full, _ = model._compute_embedding_for_batch(batch_full)
        s_full = (emb_full @ w_t + b).item()

        # n deletion variants, each length n (one shorter than full)
        ids_masked = np.zeros((n, n), dtype=np.int64)
        counts_masked = np.zeros((n, n), dtype=np.float32)
        full_positions = np.arange(n + 1)
        for k in range(n):
            keep = np.delete(full_positions, k + 1)
            ids_masked[k] = ids_full[keep]
            counts_masked[k] = counts_full[keep]

        batch_masked = _dummy_batch(ids_masked, counts_masked)
        _, _, _, emb_masked, _ = model._compute_embedding_for_batch(batch_masked)
        s_masked = (emb_masked @ w_t + b).cpu().numpy()  # (n,)

        delta = s_full - s_masked
        cols = np.array(gene_cols)
        np.add.at(gene_sum, cols, delta)
        np.add.at(gene_sumsq, cols, delta ** 2)
        np.add.at(gene_count, cols, 1)

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

out_path = f"{OUT_DIR}/state_masking_gene_attribution_ranked.csv"
aggregate_and_save(gene_sum, gene_sumsq, gene_count, var_names, out_path, MIN_GENE_OCCURRENCE)
if os.path.exists(CHECKPOINT_PATH):
    os.remove(CHECKPOINT_PATH)
    print(f"Removed checkpoint file (run complete).", flush=True)
