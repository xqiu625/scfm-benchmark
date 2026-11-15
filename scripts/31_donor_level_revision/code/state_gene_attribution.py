#!/usr/bin/env python3
"""
STATE (Arc Institute SE-600M) gene-level attribution of the mortality-associated
direction (within-COVID cohort, 240 patients / 25 fatal), for cross-model GSEA
comparison (302006 project, gene-attribution follow-on analysis). Companion to
geneformer_gene_attribution.py -- same cohort/subsample/classifier convention,
different architecture-specific attribution mechanics.

Architecture recon (verified by reading state_env's arc-state==0.9.27 source,
not assumed):
- STATE embeds a cell as a stochastically-sampled "cell sentence" of up to
  pad_length=2048 gene-identity indices into a GLOBAL 19790-gene index space
  (state.emb.data.loader.VCIDatasetSentenceCollator.sample_cell_sentences).
  Each index is looked up in a FROZEN nn.Embedding (model.pe_embedding, built
  from ESM2 protein embeddings, from_pretrained() so requires_grad=False by
  default -- must be explicitly re-enabled or no gradient reaches it).
- Position 0 is overwritten post-lookup with a learned CLS token (in-place),
  so it is excluded from attribution (its pre-overwrite value never actually
  reaches the transformer).
- `mask` (not the raw index value) is the authoritative real-vs-padding flag
  for a sentence position -- padding also uses index 0, which can coincide
  with a legitimately-sampled real gene (global index 0 = 'A1BG'), so masking
  by index value alone would be wrong.
- Checkpoint selection mirrors the project's own production invocation
  (`state emb transform --model-folder .../03_MODELS/SE-600M`, used to build
  02_EMBEDDINGS/option_c/state_covid_canonical.h5ad): the CLI does
  `glob.glob(model_folder + "/*.ckpt")[-1]`, which resolves to
  se600m_epoch4.ckpt on this filesystem (verified directly, not assumed) --
  replicated here identically for consistency with existing STATE results.
- Gene identity: protein_embeddings.pt keys are gene SYMBOLS directly (no
  Ensembl bridging needed, unlike Geneformer); var_names of the raw COVID
  h5ad overlap 12600/13066 (96.4%) with these keys.

Attribution: gradient x input at the pe_embedding lookup output (same method,
same layer role, as Geneformer's word_embeddings lookup), target = pooled
cell embedding projected onto one in-sample L2 logistic-regression mortality
direction (identical convention to Geneformer). Aggregated to one ranked gene
list per project convention (mean attribution per gene across occurrences,
min-occurrence filter).
"""
import glob
import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc
import torch
import torch.nn.functional as F
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

from state.emb.inference import Inference
from state.emb.data import create_dataloader

RAW_H5AD = "/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/human_platelet_covid_severity_raw_counts.h5ad"
PATIENT_LIST_CSV = ("/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/"
                     "31_donor_level_revision/results_patient_v2/qc_baseline/patient_qc_distributions.csv")
SAMPLE_MAP_CSV = ("/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/"
                   "31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv")
MODEL_FOLDER = "/bigdata/godziklab/shared/Xinru/302006/03_MODELS/SE-600M"
OUT_DIR = "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/gene_attribution_gsea/state"
os.makedirs(OUT_DIR, exist_ok=True)

CELLS_PER_PATIENT = 60  # same convention as geneformer_gene_attribution.py
RNG_SEED = 20260910
MIN_GENE_OCCURRENCE = 30
BATCH = 32  # cell-sentence batches are big (2048 x 5120 float32 per cell); keep modest

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Device: {device}", flush=True)

# ---------------------------------------------------------------------
# 1. Data loading and cohort restriction (identical to geneformer_gene_attribution.py)
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
a.obs_names_make_unique()

# ---------------------------------------------------------------------
# 2. Load STATE (SE-600M), matching the project's own production checkpoint
#    resolution exactly (glob.glob(...)[-1], not sorted -- same as
#    state._cli._emb._transform.run_emb_transform).
# ---------------------------------------------------------------------
ckpt_files = glob.glob(os.path.join(MODEL_FOLDER, "*.ckpt"))
assert ckpt_files, f"no checkpoint found in {MODEL_FOLDER}"
checkpoint = ckpt_files[-1]
print(f"Using checkpoint: {checkpoint} "
      f"(matches production state_covid_canonical.h5ad generation)", flush=True)

protein_embeds = torch.load(os.path.join(MODEL_FOLDER, "protein_embeddings.pt"),
                             weights_only=False, map_location="cpu")
print(f"Loaded protein embeddings: {len(protein_embeds)} genes", flush=True)
gene_symbols_by_global_idx = list(protein_embeds.keys())  # pe_embedding row order

inferer = Inference(cfg=None, protein_embeds=protein_embeds)
inferer.load_model(checkpoint)
model = inferer.model
model.eval()
model.to(device)

# pe_embedding is nn.Embedding.from_pretrained(...) -> frozen by default.
# Must explicitly re-enable requires_grad or no gradient reaches gene identity.
model.pe_embedding.weight.requires_grad_(True)
print(f"pe_embedding requires_grad after fix: {model.pe_embedding.weight.requires_grad}", flush=True)

cfg = inferer._vci_conf
dataloader = create_dataloader(
    cfg, adata=a, adata_name="covid_mortality_cohort", protein_embeds=protein_embeds,
    shuffle=False, gene_column=None, precision=None, workers=0,
)
n_batches = len(dataloader)
print(f"Dataloader built: {n_batches} batches of up to batch_size={cfg.model.batch_size}", flush=True)

# ---------------------------------------------------------------------
# 3. Pass 1 (no grad): pooled cell embeddings -> fit mortality classifier
# ---------------------------------------------------------------------
y_cell = a.obs["y"].to_numpy()
patient_cell = a.obs["patient_id"].to_numpy()

all_embeddings = []
with torch.no_grad():
    for i, batch in enumerate(dataloader):
        batch = [b.to(device) if isinstance(b, torch.Tensor) else b for b in batch]
        _, _, _, embedding, _ = model._compute_embedding_for_batch(batch)
        all_embeddings.append(embedding.detach().cpu().float().numpy())
        if (i + 1) % 20 == 0:
            print(f"  pass1 embedding batch {i+1}/{n_batches}", flush=True)
pooled_all = np.concatenate(all_embeddings, axis=0)
assert pooled_all.shape[0] == a.n_obs, f"embedding count mismatch: {pooled_all.shape[0]} vs {a.n_obs}"
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

w = torch.tensor(clf.coef_[0], dtype=torch.float32, device=device)
b = float(clf.intercept_[0])

pd.DataFrame([{"model": "state", "cell_auc": cell_auc, "patient_auc": patient_auc,
               "n_cells": a.n_obs, "n_patients": a.obs["patient_id"].nunique()}]) \
    .to_csv(f"{OUT_DIR}/classifier_diagnostics.csv", index=False)

# ---------------------------------------------------------------------
# 4. Pass 2 (grad): gradient x input attribution at pe_embedding output
# ---------------------------------------------------------------------
captured = {}


def _capture_first_call(module, inp, out):
    if "emb" not in captured:
        captured["emb"] = out


handle = model.pe_embedding.register_forward_hook(_capture_first_call)

gene_sum = np.zeros(len(gene_symbols_by_global_idx), dtype=np.float64)
gene_sumsq = np.zeros(len(gene_symbols_by_global_idx), dtype=np.float64)
gene_count = np.zeros(len(gene_symbols_by_global_idx), dtype=np.int64)

n_pos0_checked = 0
pos0_grad_absmax = 0.0

for i, batch in enumerate(dataloader):
    captured.clear()
    batch_dev = [b.to(device) if isinstance(b, torch.Tensor) else b for b in batch]
    batch_sentences_idx = batch_dev[0]  # (B, L) int, global gene index
    mask = batch_dev[5].to(torch.bool)  # (B, L) True = real gene position

    _, _, _, embedding, _ = model._compute_embedding_for_batch(batch_dev)
    logit = embedding @ w + b  # (B,)
    grad, = torch.autograd.grad(logit.sum(), captured["emb"], retain_graph=False)
    attr = (grad * captured["emb"].detach()).sum(dim=2)  # (B, L) gradient x input

    if i == 0:
        # sanity check (see module docstring): CLS overwrite should sever the
        # graph at position 0, so its gradient should be ~0 regardless of
        # which gene index happened to land there.
        pos0_grad_absmax = grad[:, 0, :].abs().max().item()
        n_pos0_checked = grad.shape[0]

    idx_np = batch_sentences_idx.cpu().numpy()
    mask_np = mask.cpu().numpy()
    attr_np = attr.detach().cpu().numpy()

    mask_np[:, 0] = False  # exclude CLS slot explicitly regardless of grad
    valid_b, valid_l = np.where(mask_np)
    valid_genes = idx_np[valid_b, valid_l]
    valid_attr = attr_np[valid_b, valid_l]

    np.add.at(gene_sum, valid_genes, valid_attr)
    np.add.at(gene_sumsq, valid_genes, valid_attr ** 2)
    np.add.at(gene_count, valid_genes, 1)

    if (i + 1) % 20 == 0:
        print(f"  pass2 attribution batch {i+1}/{n_batches}", flush=True)

handle.remove()
print(f"Position-0 (CLS slot) sanity check: max |grad| over {n_pos0_checked} cells "
      f"in first batch = {pos0_grad_absmax:.6g} (expect ~0 -- confirms CLS overwrite "
      f"severs the graph there, as reasoned in the module docstring)", flush=True)

# ---------------------------------------------------------------------
# 5. Aggregate to ranked gene list
# ---------------------------------------------------------------------
keep = gene_count >= MIN_GENE_OCCURRENCE
mean_attr = np.divide(gene_sum, gene_count, out=np.zeros_like(gene_sum), where=gene_count > 0)
var_attr = np.divide(gene_sumsq, gene_count, out=np.zeros_like(gene_sum), where=gene_count > 0) - mean_attr ** 2
sd_attr = np.sqrt(np.clip(var_attr, 0, None))

ranked = pd.DataFrame({
    "gene": [gene_symbols_by_global_idx[i] for i in np.where(keep)[0]],
    "mean_attribution": mean_attr[keep],
    "sd_attribution": sd_attr[keep],
    "n_occurrences": gene_count[keep],
}).sort_values("mean_attribution", ascending=False)

print(f"Genes in ranked list (>= {MIN_GENE_OCCURRENCE} occurrences): {len(ranked)}", flush=True)
print(f"Attribution value range: [{ranked.mean_attribution.min():.4f}, "
      f"{ranked.mean_attribution.max():.4f}], overall sd across genes = "
      f"{ranked.mean_attribution.std():.4f}", flush=True)

out_path = f"{OUT_DIR}/state_gene_attribution_ranked.csv"
ranked.to_csv(out_path, index=False)
print(f"Saved {out_path}", flush=True)

print("\nTop 20 genes by mean attribution:", flush=True)
print(ranked.head(20).to_string(index=False), flush=True)
print("\nBottom 10 genes by mean attribution:", flush=True)
print(ranked.tail(10).to_string(index=False), flush=True)
