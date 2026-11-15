#!/usr/bin/env python3
"""
UCE masking-based ("leave-one-gene-out") gene attribution of the
mortality-associated direction (canonical 240-patient/25-fatal within-COVID
cohort). Fourth model in the unified masking design, after Geneformer,
STATE, and scGPT -- see masking_attribution_common.py's docstring for the
shared rationale, and ~/Desktop/scFM_masking_attribution_illustration.md
for the full cross-model architecture writeup.

Architecture facts this script depends on, all verified directly against
our own vendored copy of snap-stanford/UCE (03_MODELS/uce/UCE), not
assumed from the earlier design-doc review alone:

- UCE's TransformerModel.forward(src, mask) does
  `src_key_padding_mask=(1-mask)` with NO dtype cast -- confirmed a real
  bug in the RELEASED, OFFICIAL inference path (evaluate.py's run_eval,
  the standard entry point via eval_single_anndata.py), not a stray/
  unused code branch: sample_cell_sentences (eval_data.py) builds `mask`
  via `torch.zeros(...)` (float32 by default) and
  `torch.concat((torch.ones(i), torch.zeros(remainder)))` -- so
  `(1-mask)` stays float32, and PyTorch treats a FLOAT src_key_padding_mask
  as an ADDITIVE BIAS on attention logits, not exclusion. Padding
  positions get +1.0 bias (mild up-weighting) instead of being ignored.
  THIS SCRIPT USES THE FIXED VERSION (`.bool()` cast) for both s_full and
  every s_masked variant -- a deliberate, necessary deviation from UCE's
  as-released behavior, not an optional preference: without the fix,
  flag-based masking would have almost no ablation effect (a masked gene
  would still contribute its actual embedding to attention, only losing a
  uniform +1.0 bias shared by every padding position), making the whole
  method non-functional. Because model.forward()'s buggy (1-mask) step
  is computed INSIDE that method with no way to inject a dtype fix from
  outside, this script does NOT call model.forward() -- it replicates its
  body directly (model.encoder / model.pos_encoder / model.transformer_encoder
  / model.decoder, called in the same order) with the one-line fix applied.
- Pooling: CLS, but read AFTER the decoder MLP -- embedding =
  gene_output[0, :, :] (position 0 along the sequence-first axis),
  L2-normalized. No override of the CLS slot's embedding afterward (unlike
  STATE's learned nn.Parameter overwrite) -- UCE's CLS is just index 3 in
  the SAME frozen ESM2-derived embedding table as every gene, going
  through the identical encoder/pos_encoder/transformer path.
- pos_encoder (positional encoding) IS active for UCE, unlike STATE (no
  positional encoding at all) -- irrelevant to this script's masking
  design though, since masking flips an attention flag rather than
  physically shortening the sequence (unlike STATE's forced deletion), so
  every position keeps its original index across all masking variants.
- Tokenization is NOT expression-rank ordered like the other 3 models.
  UCE's native sample_cell_sentences (a) samples genes STOCHASTICALLY
  WITH REPLACEMENT weighted by expression (so the same gene can appear
  multiple times in one native cell sentence), then (b) orders the
  sampled genes by CHROMOSOME (chromosomes themselves randomly shuffled
  per cell) then by GENOMIC START POSITION within each chromosome, with
  a chromosome-identity "open" token (chrom_code + CHROM_TOKEN_OFFSET)
  and a fixed "close" token (chrom_token_right_idx=2) bracketing each
  chromosome's gene block. This script builds a DETERMINISTIC sentence
  instead: the same shared top-K candidate genes as the other 3 models
  (expression-ranked SELECTION, via top_k_genes_per_cell), each used
  EXACTLY ONCE (no with-replacement duplication), re-ordered into UCE's
  native chromosome-then-position token structure (chromosomes visited in
  a fixed sorted order, not UCE's native per-cell random shuffle) --
  documented deviations from native stochastic behavior, structurally
  analogous to STATE's deterministic-sentence bypass of its own native
  sampler, necessary for a controlled, reproducible masking design.
- Checkpoint loading mirrors evaluate.py's run_eval() exactly: a
  placeholder-shaped pe_embedding is set BEFORE model.load_state_dict()
  (the checkpoint's own pe_embedding.weight entry is a dummy of the right
  shape, not real trained values -- UCE's protein embeddings are frozen
  and supplied externally at inference time), then REPLACED with the real
  frozen ESM2 + chromosome-identity table AFTER loading. The combined
  embedding table is base ESM2 tokens (all_tokens.torch, 143,574 rows)
  plus 1,895 deterministically-seeded (torch.manual_seed(23)) chromosome-
  identity vectors appended after -- replicated exactly from
  evaluate.py's get_ESM2_embeddings().
- Gene-to-token mapping requires TWO separate UCE-native lookups to
  succeed (verified via data_proc/data_utils.py's
  adata_path_to_prot_chrom_starts, which raises if either fails): (a) the
  gene's uppercased symbol in UCE's own human ESM2 protein-embedding
  dict, offset by species_offsets.pkl's human value (13,466, verified
  directly) to get its row in the shared multi-species embedding table;
  (b) the gene's uppercased symbol in species_chrom.csv (species=="human")
  for its chromosome + genomic start position. The chromosome CODE itself
  must come from a categorical fit over the FULL (all-species)
  species_chrom.csv column, THEN subsetted to human -- fitting the
  categorical on the human-only subset would reassign different codes,
  breaking the chrom-identity token's alignment with the checkpoint's
  trained chromosome-embedding rows. Replicated exactly (fit-then-subset
  order preserved) rather than reimplemented from a human-only view.
"""
import math
import os
import sys
import time

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from masking_attribution_common import (
    load_canonical_cohort, shared_gene_universe_4way, top_k_genes_per_cell,
    aggregate_and_save, save_checkpoint_atomic, TOP_K_GENES, MIN_GENE_OCCURRENCE,
)

UCE_REPO = "/bigdata/godziklab/shared/Xinru/302006/03_MODELS/uce/UCE"
sys.path.insert(0, UCE_REPO)
from model import TransformerModel

MODEL_DIR = f"{UCE_REPO}/model_files"
CKPT_PATH = f"{MODEL_DIR}/4layer_model.torch"
ALL_TOKENS_PATH = f"{MODEL_DIR}/all_tokens.torch"
OUT_DIR = "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/gene_attribution_gsea/uce_masking"
os.makedirs(OUT_DIR, exist_ok=True)

# Special token indices, verified from eval_single_anndata.py's own argparse
# defaults (this project's checkpoint uses these exact values, not
# re-derived from the checkpoint itself since they're architecture
# constants, not learned).
PAD_TOKEN_IDX = 0
CHROM_TOKEN_RIGHT_IDX = 2
CLS_TOKEN_IDX = 3
CHROM_TOKEN_OFFSET = 143574
HUMAN_OFFSET = 13466  # verified directly from species_offsets.pkl
TOKEN_DIM = 5120
D_MODEL = 1280
NHEAD = 20
D_HID = 5120
NLAYERS = 4
OUTPUT_DIM = 1280
DROPOUT = 0.05
TOTAL_VOCAB = 145469  # 143574 base tokens + 1895 chromosome-identity tokens

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
valid_mask, _sym2ens, _token_dict, _protein_embeds, _scgpt_vocab, uce_pe_raw, chrom_df = \
    shared_gene_universe_4way(var_names)

print(f"Selecting top-{TOP_K_GENES} candidate genes per cell...", flush=True)
top_k_idx = top_k_genes_per_cell(a.X, valid_mask, k=TOP_K_GENES)
n_per_cell = np.array([len(x) for x in top_k_idx])
print(f"Candidates per cell: min={n_per_cell.min()}, max={n_per_cell.max()}, "
      f"mean={n_per_cell.mean():.1f} (target {TOP_K_GENES})", flush=True)

# ---------------------------------------------------------------------
# Gene -> (global embedding row, chromosome code, genomic start) lookup,
# replicating data_proc.data_utils.adata_path_to_prot_chrom_starts exactly.
# ---------------------------------------------------------------------
pe_keys_upper = [k.upper() for k in uce_pe_raw.keys()]  # preserves original file order
pe_key_to_local_idx = {k: i for i, k in enumerate(pe_keys_upper)}

# Fit the chromosome categorical on the FULL (all-species) column first,
# THEN subset to human -- fitting on a human-only view would reassign
# different codes, breaking alignment with the checkpoint's trained
# chromosome-identity embedding rows.
chrom_df = chrom_df.copy()
chrom_df["spec_chrom"] = pd.Categorical(
    chrom_df["species"].astype(str) + "_" + chrom_df["chromosome"].astype(str)
)
chrom_df["chrom_code"] = chrom_df["spec_chrom"].cat.codes
human_chrom = chrom_df[chrom_df["species"] == "human"].copy()
human_chrom["gene_upper"] = human_chrom["gene_symbol"].str.upper()
human_chrom = human_chrom.drop_duplicates(subset="gene_upper").set_index("gene_upper")

col_to_global_idx = {}
col_to_chrom_code = {}
col_to_start = {}
for i, gene in enumerate(var_names):
    if not valid_mask[i]:
        continue
    gu = gene.upper()
    local_idx = pe_key_to_local_idx.get(gu)
    if local_idx is None or gu not in human_chrom.index:
        continue
    col_to_global_idx[i] = local_idx + HUMAN_OFFSET
    row = human_chrom.loc[gu]
    col_to_chrom_code[i] = int(row["chrom_code"])
    col_to_start[i] = int(row["start"])
print(f"col_to_global_idx built for {len(col_to_global_idx)} / {int(valid_mask.sum())} "
      f"shared-universe genes", flush=True)

# ---------------------------------------------------------------------
# Load model, mirroring evaluate.py's run_eval() checkpoint-loading order
# exactly: placeholder pe_embedding -> load_state_dict -> real pe_embedding.
# ---------------------------------------------------------------------
model = TransformerModel(token_dim=TOKEN_DIM, d_model=D_MODEL, nhead=NHEAD,
                          d_hid=D_HID, nlayers=NLAYERS, dropout=DROPOUT,
                          output_dim=OUTPUT_DIM)
empty_pe = torch.zeros(TOTAL_VOCAB, TOKEN_DIM)
model.pe_embedding = torch.nn.Embedding.from_pretrained(empty_pe)
state_dict = torch.load(CKPT_PATH, map_location="cpu", weights_only=False)
model.load_state_dict(state_dict, strict=True)

all_tokens = torch.load(ALL_TOKENS_PATH, map_location="cpu", weights_only=False)
if all_tokens.shape[0] == 143574:
    torch.manual_seed(23)
    chrom_tensors = torch.normal(mean=0, std=1, size=(1895, TOKEN_DIM))
    all_pe = torch.vstack((all_tokens, chrom_tensors))
else:
    all_pe = all_tokens
assert all_pe.shape[0] == TOTAL_VOCAB, f"unexpected combined vocab size {all_pe.shape[0]}"
all_pe.requires_grad = False
model.pe_embedding = torch.nn.Embedding.from_pretrained(all_pe)

model.to(device)
model.eval()
print(f"Model loaded. d_model={D_MODEL}, nlayers={NLAYERS}", flush=True)


def uce_forward(src, mask):
    """Replicates TransformerModel.forward()'s body exactly, with the
    float-vs-bool src_key_padding_mask dtype bug fixed (see module
    docstring). src: (L, B, token_dim), already pe_embedding'd +
    L2-normalized. mask: (B, L) float, 1.0=attend, 0.0=exclude (same
    convention as UCE's own sample_cell_sentences)."""
    src = model.encoder(src) * math.sqrt(model.d_model)
    src = model.pos_encoder(src)
    src_key_padding_mask = (1 - mask).bool()
    output = model.transformer_encoder(src, src_key_padding_mask=src_key_padding_mask)
    gene_output = model.decoder(output)
    embedding = gene_output[0, :, :]
    embedding = F.normalize(embedding, dim=1)
    return gene_output, embedding


def build_sentence(gene_cols):
    """[CLS] + candidate genes grouped by chromosome (fixed sorted order,
    not UCE's native per-cell random shuffle) then by genomic start
    position within each chromosome, with chrom-open/close marker tokens
    bracketing each block -- matches UCE's native token STRUCTURE, using
    a deterministic candidate SET (see module docstring for why).
    Returns (ids int64 array, gene_col_at_position dict: sequence
    position -> var_names column, for real gene positions only)."""
    chrom_groups = {}
    for col in gene_cols:
        code = col_to_chrom_code.get(col)
        if code is None:
            continue
        chrom_groups.setdefault(code, []).append(col)

    ids = [CLS_TOKEN_IDX]
    gene_col_at_position = {}
    for code in sorted(chrom_groups.keys()):
        genes_sorted = sorted(chrom_groups[code], key=lambda c: col_to_start[c])
        ids.append(int(code) + CHROM_TOKEN_OFFSET)
        for col in genes_sorted:
            gene_col_at_position[len(ids)] = col
            ids.append(col_to_global_idx[col])
        ids.append(CHROM_TOKEN_RIGHT_IDX)
    return np.array(ids, dtype=np.int64), gene_col_at_position


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
        ids_list, gcap_list = [], []
        for ci in range(start, end):
            ids, gcap = build_sentence(top_k_idx[ci])
            ids_list.append(ids)
            gcap_list.append(gcap)
        max_len = max(len(x) for x in ids_list)
        ids_2d = np.full((end - start, max_len), PAD_TOKEN_IDX, dtype=np.int64)
        mask_2d = np.zeros((end - start, max_len), dtype=np.float32)
        for j, ids in enumerate(ids_list):
            ids_2d[j, :len(ids)] = ids
            mask_2d[j, :len(ids)] = 1.0

        ids_t = torch.tensor(ids_2d, device=device)
        mask_t = torch.tensor(mask_2d, device=device)
        src = model.pe_embedding(ids_t)
        src = F.normalize(src, dim=2)
        src = src.permute(1, 0, 2)  # (B, L, token_dim) -> (L, B, token_dim)
        _, embedding = uce_forward(src, mask_t)
        pooled_list.append(embedding.cpu().numpy())
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

pd.DataFrame([{"model": "uce_masking", "cell_auc": cell_auc, "patient_auc": patient_auc,
               "n_cells": n_cells, "n_patients": a.obs["patient_id"].nunique(),
               "top_k_genes": TOP_K_GENES, "embedding_dim": pooled_all.shape[1]}]) \
    .to_csv(f"{OUT_DIR}/classifier_diagnostics.csv", index=False)

# ---------------------------------------------------------------------
# Pass 2: masking sweep -- flip the mask flag at one gene position per
# variant (all variants share the same sequence length, like scGPT --
# the fixed bool mask genuinely excludes the flagged position from
# attention, no physical deletion needed).
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
        ids, gcap = build_sentence(top_k_idx[i])
        positions = sorted(gcap.keys())
        n = len(positions)
        if n == 0:
            continue
        L = len(ids)
        Bsz = n + 1

        ids_batch = np.tile(ids, (Bsz, 1))
        mask_batch = np.ones((Bsz, L), dtype=np.float32)
        for k, pos in enumerate(positions):
            mask_batch[k + 1, pos] = 0.0  # exclude this gene's position

        ids_t = torch.tensor(ids_batch, device=device)
        mask_t = torch.tensor(mask_batch, device=device)
        src = model.pe_embedding(ids_t)
        src = F.normalize(src, dim=2)
        src = src.permute(1, 0, 2)
        _, embedding = uce_forward(src, mask_t)
        logits = (embedding @ w_t + b).cpu().numpy()

        s_full = logits[0]
        s_masked = logits[1:]
        delta = s_full - s_masked
        cols_arr = np.array([gcap[p] for p in positions])

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

out_path = f"{OUT_DIR}/uce_masking_gene_attribution_ranked.csv"
aggregate_and_save(gene_sum, gene_sumsq, gene_count, var_names, out_path, MIN_GENE_OCCURRENCE)
if os.path.exists(CHECKPOINT_PATH):
    os.remove(CHECKPOINT_PATH)
    print(f"Removed checkpoint file (run complete).", flush=True)
