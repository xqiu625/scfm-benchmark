#!/usr/bin/env python3
"""
Cheap standalone validation of state_gene_attribution.py's assumptions on a
tiny (30-cell) subset, before spending a full GPU job. Checks:
1. mask.sum(dim=1) is a plausible number of expressed genes per platelet cell
   (mask is claimed to be the authoritative real-vs-padding flag, not the
   raw index==0 value, since global index 0 = a real gene 'A1BG').
2. Position-0 (CLS slot) gradient really is ~0, confirming the in-place CLS
   overwrite severs the autograd graph there as reasoned in the pipeline
   script's docstring -- if this is wrong, position-0 attributions would
   need separate handling.
3. Gradients at real (masked-in) positions are non-zero and finite (catches
   the "frozen pe_embedding => zero gradient" failure mode before it wastes
   a GPU allocation).
4. Decoded gene symbols at a few real positions match highly-expressed genes
   in that same cell's raw counts (catches an index/gene-symbol misalignment
   bug, same class of bug that broke Geneformer's token dictionary earlier).
"""
import glob
import os

import numpy as np
import pandas as pd
import scanpy as sc
import torch

from state.emb.inference import Inference
from state.emb.data import create_dataloader

RAW_H5AD = "/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/human_platelet_covid_severity_raw_counts.h5ad"
MODEL_FOLDER = "/bigdata/godziklab/shared/Xinru/302006/03_MODELS/SE-600M"

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Device: {device}", flush=True)

print("Loading a 30-cell slice of the raw COVID h5ad...", flush=True)
a = sc.read_h5ad(RAW_H5AD)
a = a[:30].copy()
a.obs_names_make_unique()
var_names = a.var_names.to_numpy()

ckpt_files = glob.glob(os.path.join(MODEL_FOLDER, "*.ckpt"))
checkpoint = ckpt_files[-1]
print(f"Using checkpoint: {checkpoint}", flush=True)

protein_embeds = torch.load(os.path.join(MODEL_FOLDER, "protein_embeddings.pt"),
                             weights_only=False, map_location="cpu")
gene_symbols_by_global_idx = list(protein_embeds.keys())

inferer = Inference(cfg=None, protein_embeds=protein_embeds)
inferer.load_model(checkpoint)
model = inferer.model
model.eval()
model.to(device)
model.pe_embedding.weight.requires_grad_(True)
print(f"pe_embedding requires_grad: {model.pe_embedding.weight.requires_grad}", flush=True)

cfg = inferer._vci_conf
dataloader = create_dataloader(
    cfg, adata=a, adata_name="verify_slice", protein_embeds=protein_embeds,
    shuffle=False, gene_column=None, precision=None, workers=0,
)
print(f"Dataloader built, batch_size={cfg.model.batch_size}", flush=True)

batch = next(iter(dataloader))
batch_dev = [b.to(device) if isinstance(b, torch.Tensor) else b for b in batch]
batch_sentences_idx = batch_dev[0]
mask = batch_dev[5].to(torch.bool)

print(f"\n=== Check 1: mask plausibility ===", flush=True)
n_real = mask.sum(dim=1).cpu().numpy()
print(f"mask.sum(dim=1) per cell (n=~{len(n_real)}): min={n_real.min()}, "
      f"max={n_real.max()}, mean={n_real.mean():.1f}, pad_length={mask.shape[1]}", flush=True)
n_expr_raw = np.asarray((a.X[:len(n_real)] > 0).sum(axis=1)).ravel()
print(f"Raw counts: genes with nonzero expression per cell: min={n_expr_raw.min()}, "
      f"max={n_expr_raw.max()}, mean={n_expr_raw.mean():.1f}", flush=True)
print("(expect mask-count to be <= min(n_expressed, pad_length) and in a similar "
      "ballpark -- confirms mask correctly flags real vs padded positions)", flush=True)

captured = {}


def _capture_first_call(module, inp, out):
    if "emb" not in captured:
        captured["emb"] = out


handle = model.pe_embedding.register_forward_hook(_capture_first_call)
_, _, _, embedding, _ = model._compute_embedding_for_batch(batch_dev)
w = torch.randn(embedding.shape[1], device=device)  # arbitrary direction, just for grad-flow check
logit = embedding @ w
grad, = torch.autograd.grad(logit.sum(), captured["emb"], retain_graph=False)
handle.remove()

print(f"\n=== Check 2: position-0 (CLS slot) gradient ===", flush=True)
pos0_absmax = grad[:, 0, :].abs().max().item()
other_absmax = grad[:, 1:, :].abs().max().item()
print(f"max |grad| at position 0: {pos0_absmax:.6g}", flush=True)
print(f"max |grad| at positions 1..L: {other_absmax:.6g}", flush=True)
print("(expect position-0 grad to be ~0 relative to the rest -- confirms CLS "
      "overwrite severs the graph there)", flush=True)

print(f"\n=== Check 3: gradients at real positions are non-zero and finite ===", flush=True)
mask_np = mask.cpu().numpy()
grad_np = grad.detach().cpu().numpy()
grad_at_real = grad_np[mask_np]
print(f"n real (masked-in) positions in this batch: {grad_at_real.shape[0]}", flush=True)
print(f"grad at real positions: any_nan={np.isnan(grad_at_real).any()}, "
      f"any_inf={np.isinf(grad_at_real).any()}, "
      f"frac_exactly_zero={(grad_at_real == 0).all(axis=1).mean():.4f}, "
      f"abs_mean={np.abs(grad_at_real).mean():.6g}", flush=True)

print(f"\n=== Check 4: decoded gene symbols vs actual top-expressed genes (cell 0) ===", flush=True)
idx_np = batch_sentences_idx.cpu().numpy()
cell0_mask = mask_np[0].copy()
cell0_mask[0] = False  # exclude CLS
cell0_idx = idx_np[0][cell0_mask]
decoded_genes = set(gene_symbols_by_global_idx[i] for i in cell0_idx)
print(f"cell 0: {len(decoded_genes)} unique decoded gene symbols from sentence positions", flush=True)

row0 = np.asarray(a.X[0].todense()).ravel() if hasattr(a.X[0], "todense") else np.asarray(a.X[0]).ravel()
top_expressed_order = np.argsort(-row0)[:30]
top_expressed_genes = set(var_names[top_expressed_order]) & set(gene_symbols_by_global_idx)
overlap = decoded_genes & top_expressed_genes
print(f"top-30-expressed genes (that exist in protein_embeds): {sorted(top_expressed_genes)[:15]} ...", flush=True)
print(f"overlap between decoded cell-sentence genes and top-30-expressed genes: "
      f"{len(overlap)} -- {sorted(overlap)}", flush=True)
print("(expect substantial overlap -- STATE's own sampler should preferentially "
      "include highly-expressed genes in the cell sentence; near-zero overlap "
      "would indicate a global-index decode bug)", flush=True)

print("\nDone.", flush=True)
