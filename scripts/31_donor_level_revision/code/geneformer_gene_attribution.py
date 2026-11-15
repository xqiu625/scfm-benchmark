#!/usr/bin/env python3
"""
Geneformer gene-level attribution of the mortality-associated direction
(within-COVID cohort, 240 patients / 25 fatal), for cross-model GSEA
comparison (302006 project, gene-attribution follow-on analysis).

Pipeline:
1. Load raw within-COVID platelet counts, restrict to the canonical 240-patient
   cohort (patient_qc_distributions.csv, disease==COVID-19), stratified
   per-patient cell subsample for tractability.
2. Tokenize identically to the existing frozen-embedding pipeline
   (jobS_embed/geneformer_chongqing.sh): rank-value gene tokens via
   token_dictionary.pkl, max_len=512, native BertModel (scpeft_env),
   mean-pool last_hidden_state over the attention mask.
3. Forward pass (no grad) over all sampled cells -> pooled embeddings.
4. Fit ONE in-sample L2 (C=0.3) logistic regression on cell-level embeddings
   -> mortality label. Its weight vector w IS "the mortality-associated
   direction" in Geneformer's embedding space (project convention: cell-level
   fit, patient-level score = mean of cell probabilities).
5. Gradient attribution: for each sampled cell, forward again WITH gradients
   at the token-embedding-lookup layer (inputs_embeds path), scalar target =
   pooled_embedding . w (raw logit projection, monotonic in probability, no
   need for the sigmoid nonlinearity). Backprop to get d(target)/d(embeds).
   Attribution per token = grad . embeds ("gradient x input" -- a first-order
   Taylor / DeepLIFT-style attribution; chosen over raw gradient magnitude
   because it accounts for how large the actual embedding activation is, not
   just local sensitivity, and is the standard choice for this exact
   discrete-token/continuous-embedding architecture pattern).
6. Map token -> gene (reverse of the same token_dictionary.pkl used for
   tokenization), aggregate to one ranked gene list (mean attribution per
   gene across all cell-occurrences, min_occurrence filter).
7. Save the ranked gene list (for GSEA, run as a separate quick step) plus
   diagnostics (classifier AUC, attribution variance, per-gene occurrence
   counts) so a collapsed/degenerate signal would be caught here, not later.
"""
import os
import pickle
import sys

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import torch
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from transformers import BertModel

RAW_H5AD = "/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/human_platelet_covid_severity_raw_counts.h5ad"
PATIENT_LIST_CSV = ("/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/"
                     "31_donor_level_revision/results_patient_v2/qc_baseline/patient_qc_distributions.csv")
SAMPLE_MAP_CSV = ("/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/"
                   "31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv")
MODEL_PATH = "/bigdata/godziklab/shared/Xinru/302006/03_MODELS/geneformer/geneformer-V2-316M-native"
TOKEN_DICT = ("/bigdata/godziklab/shared/Xinru/302006/03_MODELS/geneformer/"
              "geneformer-V2-316M-native/geneformer/token_dictionary_gc104M.pkl")
OUT_DIR = "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/gene_attribution_gsea/geneformer"
os.makedirs(OUT_DIR, exist_ok=True)

MAX_LEN = 512
CELLS_PER_PATIENT = 60  # stratified subsample cap, per patient, for tractability
RNG_SEED = 20260910
MIN_GENE_OCCURRENCE = 30  # minimum cell-occurrences for a gene to be included in the ranked list

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Device: {device}", flush=True)

# ---------------------------------------------------------------------
# 1. Data loading and cohort restriction
# ---------------------------------------------------------------------
patient_list = pd.read_csv(PATIENT_LIST_CSV)
patient_list = patient_list[patient_list["disease"] == "COVID-19"]
assert len(patient_list) == 240 and int(patient_list["y"].sum()) == 25, \
    f"cohort mismatch: n={len(patient_list)}, fatal={int(patient_list['y'].sum())}"
y_by_patient = patient_list.set_index("patient_id")["y"].to_dict()
covid_patient_ids = set(patient_list["patient_id"].astype(str))
print(f"Canonical cohort: {len(covid_patient_ids)} patients, "
      f"{sum(y_by_patient.values())} fatal", flush=True)

sample_map = pd.read_csv(SAMPLE_MAP_CSV)
sample_map["patient_id"] = sample_map["patient_id"].astype(str)
sample_map = sample_map[sample_map["patient_id"].isin(covid_patient_ids)]
sample_to_patient = sample_map.set_index("sample_id")["patient_id"].to_dict()

print(f"Loading {RAW_H5AD} ...", flush=True)
a = sc.read_h5ad(RAW_H5AD)
a.obs["patient_id"] = a.obs["Sample_ID"].astype(str).map(sample_to_patient)
a = a[a.obs["patient_id"].notna()].copy()
a.obs["y"] = a.obs["patient_id"].map(y_by_patient).astype(int)
n_patients_present = a.obs["patient_id"].nunique()
print(f"Cells mapped to canonical cohort: {a.n_obs} cells, "
      f"{n_patients_present} / 240 patients present in this h5ad", flush=True)

# stratified per-patient cell subsample (cap, not proportional oversample --
# every patient contributes at most CELLS_PER_PATIENT cells regardless of
# class, so the fatal/survived cell-count ratio tracks the 25/215 PATIENT
# ratio rather than being skewed by cells-per-patient variation)
rng = np.random.default_rng(RNG_SEED)
keep_idx = []
for pid, grp in a.obs.groupby("patient_id"):
    n = len(grp)
    take = min(n, CELLS_PER_PATIENT)
    chosen = rng.choice(grp.index.to_numpy(), size=take, replace=False)
    keep_idx.extend(chosen.tolist())
a = a[keep_idx].copy()
print(f"After stratified per-patient subsample (cap={CELLS_PER_PATIENT}/patient): "
      f"{a.n_obs} cells, {a.obs['y'].sum()} from fatal patients' cells "
      f"(fatal patients contribute {a.obs.loc[a.obs.y==1,'patient_id'].nunique()} / "
      f"{a.obs['patient_id'].nunique()} patients)", flush=True)

# ---------------------------------------------------------------------
# 2. Tokenization. NOTE: this h5ad's var_names are gene SYMBOLS (verified
# via h5py peek: "LINC00115", "NOC2L", ...), but token_dictionary.pkl keys
# are Ensembl IDs ("ENSG..."). jobS_embed/geneformer_chongqing.sh's naive
# `if gene in token_dict` lookup only worked there because it embeds a
# DIFFERENT h5ad (chongqing_sepsis_platelet_raw_counts.h5ad) whose var_names
# happen to already be Ensembl IDs -- confirmed by first attempt here
# failing at 0/13066 mapped. Bridge via the pre-existing symbol->Ensembl
# table already used for an earlier (monocyte) pipeline in this project.
# ---------------------------------------------------------------------
with open(TOKEN_DICT, "rb") as f:
    token_dict = pickle.load(f)
SYMBOL_TO_ENSEMBL = ("/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/"
                     "option_c_monocyte/gene_symbol_to_ensembl.csv")
sym2ens = pd.read_csv(SYMBOL_TO_ENSEMBL).set_index("gene_symbol")["ensembl_id"].to_dict()
token_to_gene = {}  # reverse map for attribution (token -> gene SYMBOL, for readability)
idx_to_token = {}
n_mapped = 0
n_symbol_resolved = 0
for i, gene in enumerate(a.var_names.astype(str)):
    ens = sym2ens.get(gene)
    if ens is not None:
        n_symbol_resolved += 1
        if ens in token_dict:
            tok = token_dict[ens]
            idx_to_token[i] = tok
            token_to_gene[tok] = gene
            n_mapped += 1
print(f"Symbols resolved to Ensembl IDs: {n_symbol_resolved} / {a.n_vars}", flush=True)
print(f"Genes mapped to tokens: {n_mapped} / {a.n_vars}", flush=True)
assert n_mapped > 1000, "too few genes mapped -- tokenization broken, do not proceed"

X = a.X
def tokenize_cell(row_idx):
    row = X[row_idx].toarray().ravel() if sp.issparse(X) else np.asarray(X[row_idx]).ravel()
    nz = np.where(row > 0)[0]
    if len(nz) == 0:
        return [0] * MAX_LEN, [0] * MAX_LEN
    order = nz[np.argsort(-row[nz])]
    toks = []
    for i in order:
        t = idx_to_token.get(int(i))
        if t is not None:
            toks.append(t)
            if len(toks) >= MAX_LEN:
                break
    n = len(toks)
    if n < MAX_LEN:
        toks = toks + [0] * (MAX_LEN - n)
        mask = [1] * n + [0] * (MAX_LEN - n)
    else:
        mask = [1] * MAX_LEN
    return toks[:MAX_LEN], mask

# ---------------------------------------------------------------------
# 3. Model + forward pass (no grad) to get pooled embeddings for classifier fit
# ---------------------------------------------------------------------
model = BertModel.from_pretrained(MODEL_PATH).to(device).eval()
embed_layer = model.embeddings.word_embeddings
hidden_size = model.config.hidden_size
print(f"Model loaded. hidden_size={hidden_size}", flush=True)

n_cells = a.n_obs
all_ids = np.zeros((n_cells, MAX_LEN), dtype=np.int64)
all_mask = np.zeros((n_cells, MAX_LEN), dtype=np.int64)
for i in range(n_cells):
    toks, mask = tokenize_cell(i)
    all_ids[i] = toks
    all_mask[i] = mask
    if (i + 1) % 2000 == 0:
        print(f"  tokenized {i+1}/{n_cells}", flush=True)

y_cell = a.obs["y"].to_numpy()
patient_cell = a.obs["patient_id"].to_numpy()

BATCH = 32
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
        if (start // BATCH) % 50 == 0:
            print(f"  forward pass {end}/{n_cells}", flush=True)

print(f"Pooled embeddings: {pooled_all.shape}", flush=True)

# ---------------------------------------------------------------------
# 4. Fit the cell-level mortality classifier -> mortality direction w
# ---------------------------------------------------------------------
clf = LogisticRegression(C=0.3, solver="lbfgs", max_iter=2000)
clf.fit(pooled_all, y_cell)
p_cell = clf.predict_proba(pooled_all)[:, 1]
cell_auc = roc_auc_score(y_cell, p_cell)
patient_score = pd.Series(p_cell, index=range(n_cells)).groupby(
    pd.Series(patient_cell)).mean()
patient_y = pd.Series(y_cell, index=range(n_cells)).groupby(
    pd.Series(patient_cell)).max()
patient_auc = roc_auc_score(patient_y.reindex(patient_score.index),
                             patient_score) if patient_y.nunique() > 1 else float("nan")
print(f"Cell-level in-sample AUC = {cell_auc:.4f}; "
      f"patient-level (mean-pooled) in-sample AUC = {patient_auc:.4f}", flush=True)
w = torch.tensor(clf.coef_[0], dtype=torch.float32, device=device)  # (hidden_size,)

pd.DataFrame({"metric": ["cell_auc_in_sample", "patient_auc_in_sample", "n_cells", "n_patients"],
              "value": [cell_auc, patient_auc, n_cells, a.obs["patient_id"].nunique()]}
             ).to_csv(os.path.join(OUT_DIR, "classifier_diagnostics.csv"), index=False)

# ---------------------------------------------------------------------
# 5. Gradient (gradient x input) attribution at the token-embedding layer
# ---------------------------------------------------------------------
gene_attr_sum = {}
gene_attr_sumsq = {}
gene_attr_count = {}

model.eval()
for p in model.parameters():
    p.requires_grad_(False)

for start in range(0, n_cells, BATCH):
    end = min(start + BATCH, n_cells)
    ids = torch.tensor(all_ids[start:end], device=device)
    mask = torch.tensor(all_mask[start:end], device=device)
    embeds = embed_layer(ids).detach().clone().requires_grad_(True)
    out = model(inputs_embeds=embeds, attention_mask=mask).last_hidden_state
    m = mask.unsqueeze(-1).float()
    pooled = (out * m).sum(dim=1) / m.sum(dim=1).clamp(min=1)  # (batch, hidden)
    target = (pooled * w.unsqueeze(0)).sum(dim=1).sum()  # scalar, sum over batch is fine (independent per-row grads)
    target.backward()
    grad_x_input = (embeds.grad * embeds).sum(dim=-1).detach().cpu().numpy()  # (batch, MAX_LEN)

    ids_np = ids.cpu().numpy()
    mask_np = mask.cpu().numpy()
    for bi in range(end - start):
        for pos in range(MAX_LEN):
            if mask_np[bi, pos] == 0:
                continue
            tok = int(ids_np[bi, pos])
            gene = token_to_gene.get(tok)
            if gene is None:
                continue
            val = float(grad_x_input[bi, pos])
            gene_attr_sum[gene] = gene_attr_sum.get(gene, 0.0) + val
            gene_attr_sumsq[gene] = gene_attr_sumsq.get(gene, 0.0) + val * val
            gene_attr_count[gene] = gene_attr_count.get(gene, 0) + 1
    if (start // BATCH) % 50 == 0:
        print(f"  attribution pass {end}/{n_cells}", flush=True)

rows = []
for gene, s in gene_attr_sum.items():
    n = gene_attr_count[gene]
    if n < MIN_GENE_OCCURRENCE:
        continue
    mean_attr = s / n
    var_attr = gene_attr_sumsq[gene] / n - mean_attr ** 2
    rows.append({"gene": gene, "mean_attribution": mean_attr,
                 "sd_attribution": np.sqrt(max(var_attr, 0.0)), "n_occurrences": n})
attr_df = pd.DataFrame(rows).sort_values("mean_attribution", ascending=False).reset_index(drop=True)
print(f"Genes in ranked list (>= {MIN_GENE_OCCURRENCE} occurrences): {len(attr_df)}", flush=True)
print(f"Attribution value range: [{attr_df['mean_attribution'].min():.4f}, "
      f"{attr_df['mean_attribution'].max():.4f}], overall sd across genes = "
      f"{attr_df['mean_attribution'].std():.4f}", flush=True)

attr_df.to_csv(os.path.join(OUT_DIR, "geneformer_gene_attribution_ranked.csv"), index=False)
print(f"Saved {os.path.join(OUT_DIR, 'geneformer_gene_attribution_ranked.csv')}", flush=True)
print("\nTop 20 genes by mean attribution:")
print(attr_df.head(20).to_string(index=False))
print("\nBottom 10 genes by mean attribution:")
print(attr_df.tail(10).to_string(index=False))
