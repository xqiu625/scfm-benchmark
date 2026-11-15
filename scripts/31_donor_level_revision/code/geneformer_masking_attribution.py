#!/usr/bin/env python3
"""
Geneformer masking-based ("leave-one-gene-out") gene attribution of the
mortality-associated direction (canonical 240-patient/25-fatal within-COVID
cohort). Replaces geneformer_gene_attribution.py's gradient x input method
with a model-agnostic perturbation method shared identically (by
construction, via masking_attribution_common.py) with
state_masking_attribution.py -- see that module's docstring for the full
design rationale.

Pipeline:
1. Load cohort (same restriction/subsample as the gradient pipeline).
2. Restrict to the gene universe usable by BOTH Geneformer and STATE, and
   select each cell's top TOP_K_GENES by raw expression from that universe
   (deterministic; identical candidate set STATE will use for the same
   cell).
3. Tokenize each cell as EXACTLY that K-gene set (not Geneformer's native
   rank-512 full tokenization), mean-pool -> pooled embedding.
4. Fit ONE in-sample L2 (C=0.3) logistic regression on these pooled
   embeddings -> mortality label (same convention as before, but embeddings
   are now K-gene-restricted, so this classifier is specific to this
   pipeline, not reused from the gradient-based run).
5. For each cell, batch (n_valid_genes + 1) forward passes: the original
   plus one variant per candidate gene with that gene's position zeroed out
   of the attention mask (token id left as-is; attention_mask=0 removes it
   from mean-pooling exactly as if it were never tokenized). Delta score
   per gene = s_full - s_masked_g, s = pooled_embedding . w + b (raw
   logit).
6. Aggregate to one ranked gene list (mean delta per gene across cell
   occurrences, min-occurrence filter) via the shared aggregation function,
   in the same schema geneformer_run_gsea.py already consumes.
"""
import os
import sys
import time

import numpy as np
import torch
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from transformers import BertModel

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from masking_attribution_common import (
    load_canonical_cohort, shared_gene_universe, top_k_genes_per_cell,
    aggregate_and_save, save_checkpoint_atomic, TOP_K_GENES, MIN_GENE_OCCURRENCE,
)

MODEL_PATH = "/bigdata/godziklab/shared/Xinru/302006/03_MODELS/geneformer/geneformer-V2-316M-native"
OUT_DIR = "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/gene_attribution_gsea/geneformer_masking"
os.makedirs(OUT_DIR, exist_ok=True)

# Checkpoint/resume support (see state_masking_attribution.py for the first
# use of this pattern and its rationale) -- lets a masking sweep be chained
# through short_gpu's 2h cap instead of waiting for one uninterrupted
# multi-hour gpu-partition allocation. Not used by the currently-running
# full-scale job (started before this was added); available for future
# reruns.
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
valid_mask, sym2ens, token_dict, _ = shared_gene_universe(var_names)

idx_to_token = {}
for i, gene in enumerate(var_names):
    if valid_mask[i]:
        ens = sym2ens.get(gene)
        tok = token_dict.get(ens)
        if tok is not None:
            idx_to_token[i] = tok
print(f"idx_to_token built for {len(idx_to_token)} / {int(valid_mask.sum())} "
      f"shared-universe genes", flush=True)

print(f"Selecting top-{TOP_K_GENES} candidate genes per cell...", flush=True)
top_k_idx = top_k_genes_per_cell(a.X, valid_mask, k=TOP_K_GENES)

n_cells = a.n_obs
all_ids = np.zeros((n_cells, TOP_K_GENES), dtype=np.int64)
all_mask = np.zeros((n_cells, TOP_K_GENES), dtype=np.int64)
all_gene_col = np.full((n_cells, TOP_K_GENES), -1, dtype=np.int64)
for i in range(n_cells):
    toks, cols = [], []
    for col in top_k_idx[i]:
        tok = idx_to_token.get(int(col))
        if tok is not None:
            toks.append(tok)
            cols.append(int(col))
    n = len(toks)
    all_ids[i, :n] = toks
    all_mask[i, :n] = 1
    all_gene_col[i, :n] = cols
    if (i + 1) % 2000 == 0:
        print(f"  tokenized {i + 1}/{n_cells}", flush=True)

n_per_cell = all_mask.sum(axis=1)
print(f"Candidate genes per cell after token mapping: min={n_per_cell.min()}, "
      f"max={n_per_cell.max()}, mean={n_per_cell.mean():.1f} (target {TOP_K_GENES})", flush=True)

model = BertModel.from_pretrained(MODEL_PATH).to(device).eval()
hidden_size = model.config.hidden_size
print(f"Model loaded. hidden_size={hidden_size}", flush=True)

y_cell = a.obs["y"].to_numpy()
patient_cell = a.obs["patient_id"].to_numpy()

# ---------------------------------------------------------------------
# Pass 1: s_full pooled embeddings (no grad) -> fit mortality classifier
# ---------------------------------------------------------------------
BATCH = 64
pooled_all = np.zeros((n_cells, hidden_size), dtype=np.float32)
with torch.no_grad():
    for start in range(0, n_cells, BATCH):
        end = min(start + BATCH, n_cells)
        ids = torch.tensor(all_ids[start:end], device=device)
        mask = torch.tensor(all_mask[start:end], device=device)
        out = model(input_ids=ids, attention_mask=mask).last_hidden_state
        m = mask.unsqueeze(-1).float()
        pooled = (out * m).sum(dim=1) / m.sum(dim=1).clamp(min=1)
        pooled_all[start:end] = pooled.cpu().numpy()
        if (start // BATCH) % 30 == 0:
            print(f"  pass1 embedding batch {end}/{n_cells}", flush=True)

clf = LogisticRegression(penalty="l2", C=0.3, max_iter=2000)
clf.fit(pooled_all, y_cell)
cell_auc = roc_auc_score(y_cell, clf.predict_proba(pooled_all)[:, 1])

import pandas as pd
patient_score = pd.DataFrame({"patient_id": patient_cell, "p": clf.predict_proba(pooled_all)[:, 1]}) \
    .groupby("patient_id")["p"].mean()
patient_y = pd.Series(y_by_patient)
patient_auc = roc_auc_score(patient_y.loc[patient_score.index], patient_score)
print(f"Cell-level in-sample AUC = {cell_auc:.4f}; "
      f"patient-level (mean-pooled) in-sample AUC = {patient_auc:.4f}", flush=True)

w = clf.coef_[0].astype(np.float32)
b = float(clf.intercept_[0])
w_t = torch.tensor(w, device=device)

pd.DataFrame([{"model": "geneformer_masking", "cell_auc": cell_auc, "patient_auc": patient_auc,
               "n_cells": n_cells, "n_patients": a.obs["patient_id"].nunique(),
               "top_k_genes": TOP_K_GENES}]) \
    .to_csv(f"{OUT_DIR}/classifier_diagnostics.csv", index=False)

# ---------------------------------------------------------------------
# Pass 2: masking sweep -- one batched forward pass per cell
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
        mask_row = all_mask[i]
        valid_positions = np.where(mask_row == 1)[0]
        n_valid = len(valid_positions)
        if n_valid == 0:
            continue
        Bsz = n_valid + 1
        ids_batch = np.tile(all_ids[i], (Bsz, 1))
        mask_batch = np.tile(mask_row, (Bsz, 1))
        for k, pos in enumerate(valid_positions):
            mask_batch[k + 1, pos] = 0

        ids_t = torch.tensor(ids_batch, device=device)
        mask_t = torch.tensor(mask_batch, device=device)
        out = model(input_ids=ids_t, attention_mask=mask_t).last_hidden_state
        m = mask_t.unsqueeze(-1).float()
        pooled = (out * m).sum(dim=1) / m.sum(dim=1).clamp(min=1)
        logits = (pooled @ w_t + b).cpu().numpy()

        s_full = logits[0]
        s_masked = logits[1:]
        delta = s_full - s_masked  # (n_valid,)
        cols = all_gene_col[i, valid_positions]

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

out_path = f"{OUT_DIR}/geneformer_masking_gene_attribution_ranked.csv"
aggregate_and_save(gene_sum, gene_sumsq, gene_count, var_names, out_path, MIN_GENE_OCCURRENCE)
if os.path.exists(CHECKPOINT_PATH):
    os.remove(CHECKPOINT_PATH)
    print(f"Removed checkpoint file (run complete).", flush=True)
