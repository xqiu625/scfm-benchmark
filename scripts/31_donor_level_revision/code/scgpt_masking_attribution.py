#!/usr/bin/env python3
"""
scGPT masking-based ("leave-one-gene-out") gene attribution of the
mortality-associated direction (canonical 240-patient/25-fatal within-COVID
cohort). Third model in the unified masking design, after Geneformer and
STATE -- see masking_attribution_common.py's docstring for the shared
rationale, and ~/Desktop/scFM_masking_attribution_illustration.md for the
full cross-model architecture writeup.

Architecture facts this script depends on, verified directly against our
vendored copy of bowang-lab/scGPT (03_MODELS/scgpt/scGPT), not assumed from
an external source alone:
- TransformerModel._encode sums the gene-identity token embedding and the
  (binned) expression-value embedding (total_embs = src + values), then
  calls `self.transformer_encoder(total_embs, src_key_padding_mask=...)` --
  the mask genuinely reaches attention (model.py's _encode method).
- Checked the one real risk this project has now hit twice for other
  models (a fast/flash attention path silently forbidding or mis-handling
  masks): this checkpoint's own args.json says fast_transformer=true, and
  TransformerModel.__init__ would route to a FastTransformerEncoderWrapper
  (linear backend, which explicitly REJECTS non-trailing padding: "Found
  padding tokens in the middle of the sequence") if flash-attn were
  installed -- but scpeft_env does NOT have flash_attn installed, so
  TransformerModel.__init__'s own fallback forces use_fast_transformer=False
  regardless of the checkpoint's training-time flag, using the standard
  PyTorch TransformerEncoderLayer/TransformerEncoder instead, which
  correctly supports an arbitrary (non-trailing) src_key_padding_mask
  pattern -- exactly what leave-one-gene-out masking needs.
- Pooling is CLS: cell_emb_style="cls" (this model's default), pooling is
  literally layer_output[:, 0, :] (the correct config knob, unlike
  scPRINT-2's dead cell_emb_style setting).
- Mask convention: PyTorch src_key_padding_mask, True = ignore. Masking a
  gene here means flipping ONE position's flag to True -- NOT replacing its
  value token (that only zeroes the expression signal; gene identity still
  attends, a strictly weaker intervention) and NOT physically dropping the
  (gene, value) pair (unnecessary sequence-length change; unlike STATE,
  scGPT's mask flag genuinely has teeth so no sequence surgery is needed).
- Binning trap: DataCollator bins expression values via percentile
  quantiles computed over WHATEVER gene set is already selected
  (scgpt.preprocess.binning, quantiles over the nonzero subset). If this
  script fed only the K candidate genes into binning, it would compute
  different bin edges than scGPT would see natively (a cell's FULL nonzero
  profile, which for platelets can be thousands of genes) -- the same
  "renormalisation couples every gene together" trap flagged for
  STATE/scPRINT-2. Fixed by binning the cell's FULL nonzero row ONCE
  (matching native usage), then selecting just the K candidates' resulting
  bin values for the short sequence actually fed to the model.
- Candidate gene set: drawn from the 3-way (Geneformer & STATE & scGPT)
  shared universe (masking_attribution_common.shared_gene_universe_3way),
  not the 2-way universe Geneformer/STATE's already-completed runs used --
  see that function's docstring for why those two runs were not redone
  (7-gene, 0.06% difference, judged negligible against hours of GPU time).
- torch.load(model_file, weights_only=False) set proactively, having
  already been bitten once by PyTorch's >=2.6 weights_only=True default
  while building STATE's script.
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import torch
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from masking_attribution_common import (
    load_canonical_cohort, shared_gene_universe_3way, top_k_genes_per_cell,
    aggregate_and_save, save_checkpoint_atomic, TOP_K_GENES, MIN_GENE_OCCURRENCE,
)

SCGPT_REPO = "/bigdata/godziklab/shared/Xinru/302006/03_MODELS/scgpt/scGPT"
sys.path.insert(0, SCGPT_REPO)
from scgpt.model import TransformerModel
from scgpt.tokenizer import GeneVocab
from scgpt.utils import load_pretrained
from scgpt.preprocess import binning

MODEL_DIR = "/bigdata/godziklab/shared/Xinru/302006/03_MODELS/scgpt/scGPT_human/scGPT_human"
OUT_DIR = "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/gene_attribution_gsea/scgpt_masking"
os.makedirs(OUT_DIR, exist_ok=True)

# Checkpoint/resume support (see state_masking_attribution.py for the first
# use of this pattern). scGPT's own full-scale run (job 28277291) already
# completed in 42m38s, comfortably inside short_gpu's 2h cap, so this was
# not needed for that run -- added here as future-proofing in case a rerun
# with a larger K or cohort takes longer and needs to chain through
# multiple short_gpu allocations.
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
valid_mask, _sym2ens, _token_dict, _protein_embeds, _scgpt_vocab_json = \
    shared_gene_universe_3way(var_names)

print(f"Selecting top-{TOP_K_GENES} candidate genes per cell...", flush=True)
top_k_idx = top_k_genes_per_cell(a.X, valid_mask, k=TOP_K_GENES)
n_per_cell = np.array([len(x) for x in top_k_idx])
print(f"Candidates per cell: min={n_per_cell.min()}, max={n_per_cell.max()}, "
      f"mean={n_per_cell.mean():.1f} (target {TOP_K_GENES})", flush=True)

# ---------------------------------------------------------------------
# Load model, mirroring scgpt.tasks.cell_emb.embed_data's official recipe
# ---------------------------------------------------------------------
vocab_file = f"{MODEL_DIR}/vocab.json"
model_config_file = f"{MODEL_DIR}/args.json"
model_file = f"{MODEL_DIR}/best_model.pt"
pad_token = "<pad>"
special_tokens = [pad_token, "<cls>", "<eoc>"]

vocab = GeneVocab.from_file(vocab_file)
for s in special_tokens:
    if s not in vocab:
        vocab.append_token(s)
vocab.set_default_index(vocab["<pad>"])

with open(model_config_file) as f:
    model_configs = json.load(f)

model = TransformerModel(
    ntoken=len(vocab),
    d_model=model_configs["embsize"],
    nhead=model_configs["nheads"],
    d_hid=model_configs["d_hid"],
    nlayers=model_configs["nlayers"],
    nlayers_cls=model_configs["n_layers_cls"],
    n_cls=1,
    vocab=vocab,
    dropout=model_configs["dropout"],
    pad_token=model_configs["pad_token"],
    pad_value=model_configs["pad_value"],
    do_mvc=True,
    do_dab=False,
    use_batch_labels=False,
    domain_spec_batchnorm=False,
    explicit_zero_prob=False,
    use_fast_transformer=True,  # auto-falls back to False: flash-attn not installed here
    fast_transformer_backend="flash",
    pre_norm=False,
)
ckpt = torch.load(model_file, map_location=device, weights_only=False)
load_pretrained(model, ckpt, verbose=False)
model.to(device)
model.eval()
print(f"Model loaded. use_fast_transformer={model.use_fast_transformer} "
      f"(expect False -- flash-attn not installed, safe fallback to standard "
      f"attention which supports arbitrary mask patterns)", flush=True)

cls_id = vocab["<cls>"]
pad_id = vocab[pad_token]
pad_value = model_configs["pad_value"]

col_to_scgpt_id = {}
n_silently_defaulted = 0
for i, gene in enumerate(var_names):
    if valid_mask[i]:
        tok = vocab[gene]
        # vocab.set_default_index(vocab["<pad>"]) above means an OOV gene
        # would silently resolve to <pad> here instead of raising -- valid_mask
        # was built from a separate raw-JSON membership check (shared_gene_
        # universe_3way), not this GeneVocab object directly, so confirm the
        # two agree rather than let a mismatch silently zero out a gene's
        # attribution (delta stuck at 0, no error).
        if tok == pad_id:
            n_silently_defaulted += 1
        col_to_scgpt_id[i] = tok
assert n_silently_defaulted == 0, (
    f"{n_silently_defaulted} shared-universe genes silently resolved to <pad> "
    f"in scGPT's vocab -- the raw-JSON membership check and the GeneVocab "
    f"object disagree, do not proceed")

X = a.X


def full_row_binned(cell_idx):
    """Bin the cell's FULL nonzero expression profile (matching native
    scGPT quantile computation over all nonzero genes, not just the K
    candidates) -- returns {var_names col -> bin value} for nonzero genes."""
    row = X[cell_idx].toarray().ravel() if hasattr(X[cell_idx], "toarray") else np.asarray(X[cell_idx]).ravel()
    nz = np.nonzero(row)[0]
    if len(nz) == 0:
        return {}
    binned = binning(row[nz].astype(np.float32), n_bins=51)
    return dict(zip(nz.tolist(), np.asarray(binned).tolist()))


def build_full_sentence(cell_idx, gene_cols):
    """[<cls>] + candidate genes present in gene_cols with nonzero expression.
    Returns (genes int64 array, values float32 array, cols list) -- all
    length-matched, cols excludes the CLS slot."""
    binned_map = full_row_binned(cell_idx)
    genes = [cls_id]
    values = [float(pad_value)]
    cols = []
    for col in gene_cols:
        if col in binned_map:
            genes.append(col_to_scgpt_id[col])
            values.append(float(binned_map[col]))
            cols.append(col)
    return np.array(genes, dtype=np.int64), np.array(values, dtype=np.float32), cols


n_cells = a.n_obs
y_cell = a.obs["y"].to_numpy()
patient_cell = a.obs["patient_id"].to_numpy()

# ---------------------------------------------------------------------
# Pass 1: s_full embeddings (no grad) -> fit mortality classifier
# ---------------------------------------------------------------------
BATCH_CELLS = 32
pooled_list = []
t0 = time.time()
with torch.no_grad():
    for start in range(0, n_cells, BATCH_CELLS):
        end = min(start + BATCH_CELLS, n_cells)
        cell_idxs = list(range(start, end))
        genes_list, values_list = [], []
        for ci in cell_idxs:
            genes, values, _ = build_full_sentence(ci, top_k_idx[ci])
            genes_list.append(genes)
            values_list.append(values)
        max_len = max(len(g) for g in genes_list)
        genes_2d = np.full((len(cell_idxs), max_len), pad_id, dtype=np.int64)
        values_2d = np.full((len(cell_idxs), max_len), float(pad_value), dtype=np.float32)
        for j, (g, v) in enumerate(zip(genes_list, values_list)):
            genes_2d[j, :len(g)] = g
            values_2d[j, :len(v)] = v
        genes_t = torch.tensor(genes_2d, device=device)
        values_t = torch.tensor(values_2d, device=device)
        mask_t = genes_t.eq(pad_id)
        out = model._encode(genes_t, values_t, src_key_padding_mask=mask_t)
        emb = out[:, 0, :]
        pooled_list.append(emb.cpu().numpy())
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

pd.DataFrame([{"model": "scgpt_masking", "cell_auc": cell_auc, "patient_auc": patient_auc,
               "n_cells": n_cells, "n_patients": a.obs["patient_id"].nunique(),
               "top_k_genes": TOP_K_GENES, "embedding_dim": pooled_all.shape[1]}]) \
    .to_csv(f"{OUT_DIR}/classifier_diagnostics.csv", index=False)

# ---------------------------------------------------------------------
# Pass 2: masking sweep -- flip src_key_padding_mask, one batched forward
# pass per cell (all variants share the same sequence length, unlike STATE)
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
        genes, values, cols = build_full_sentence(i, top_k_idx[i])
        n = len(cols)
        if n == 0:
            continue
        L = len(genes)  # = n + 1 (CLS + n candidates)
        Bsz = n + 1

        genes_batch = np.tile(genes, (Bsz, 1))
        values_batch = np.tile(values, (Bsz, 1))
        mask_batch = np.zeros((Bsz, L), dtype=bool)
        for k in range(n):
            mask_batch[k + 1, k + 1] = True  # position k+1 = candidate k (0 = CLS)

        genes_t = torch.tensor(genes_batch, device=device)
        values_t = torch.tensor(values_batch, device=device)
        mask_t = torch.tensor(mask_batch, device=device)
        out = model._encode(genes_t, values_t, src_key_padding_mask=mask_t)
        emb = out[:, 0, :]
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

out_path = f"{OUT_DIR}/scgpt_masking_gene_attribution_ranked.csv"
aggregate_and_save(gene_sum, gene_sumsq, gene_count, var_names, out_path, MIN_GENE_OCCURRENCE)
if os.path.exists(CHECKPOINT_PATH):
    os.remove(CHECKPOINT_PATH)
    print(f"Removed checkpoint file (run complete).", flush=True)
