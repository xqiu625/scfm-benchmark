#!/usr/bin/env python3
"""Job F Phase 3 (Geneformer): tokenize + extract embeddings for the monocyte
canonical object using the geneformer.sif container's Geneformer install
(older API than the platelet-cohort step2c script assumed: EmbExtractor needs
a pre-tokenized dataset file, not a live AnnData, and tokenization requires
adata.var['ensembl_id']).

Must run INSIDE the geneformer.sif container (sys.path needs /opt/Geneformer).
"""
import sys
sys.path.insert(0, "/opt/Geneformer")

import os
import shutil
import numpy as np
import pandas as pd
import scanpy as sc
from geneformer import TranscriptomeTokenizer, EmbExtractor

BASE = "/bigdata/godziklab/shared/Xinru/302006"
EMB_DIR = f"{BASE}/02_EMBEDDINGS/option_c_monocyte"
INPUT = f"{EMB_DIR}/monocyte_canonical_raw.h5ad"
MAPPING_CSV = f"{EMB_DIR}/gene_symbol_to_ensembl.csv"
WORK = f"{BASE}/302006_updategit/scripts/31_donor_level_revision/tmp_geneformer_monocyte_work"
OUTPUT = f"{EMB_DIR}/geneformer_monocyte_canonical.h5ad"

os.makedirs(WORK, exist_ok=True)
ready_dir = f"{WORK}/ready"
os.makedirs(ready_dir, exist_ok=True)

print("Loading monocyte canonical input...", flush=True)
adata = sc.read_h5ad(INPUT)
print(f"Loaded: {adata.shape}", flush=True)

gmap = pd.read_csv(MAPPING_CSV, dtype=str).set_index("gene_symbol")["ensembl_id"]
adata.var["ensembl_id"] = gmap.reindex(adata.var_names).values
n_before = adata.n_vars
adata = adata[:, adata.var["ensembl_id"].notna() & (adata.var["ensembl_id"] != "")].copy()
print(f"Genes with Ensembl IDs: {adata.n_vars}/{n_before}", flush=True)

# Container's bundled anndata/scipy combo hits a `.astype(copy=False)` TypeError on
# backed sparse fancy-indexing during tokenization (same issue
# generate_geneformer_embeddings_chunked.py works around) -- densify to avoid it.
if hasattr(adata.X, "toarray"):
    print("Densifying X to avoid backed-sparse tokenizer bug...", flush=True)
    adata.X = np.asarray(adata.X.todense())

adata.obs["orig_index"] = np.arange(adata.n_obs)
if "n_counts" not in adata.obs.columns:
    x = adata.X
    adata.obs["n_counts"] = np.asarray(x.sum(axis=1)).ravel()
adata.obs["filter_pass"] = True
for col in list(adata.obs.columns):
    if adata.obs[col].dtype.name == "category":
        adata.obs[col] = adata.obs[col].astype(str)
for col in list(adata.var.columns):
    if adata.var[col].dtype.name == "category":
        adata.var[col] = adata.var[col].astype(str)
adata.raw = None

# Write in shards -- the container's tokenizer hits a "two fancy-index arrays"
# h5py error on the full 248k-cell file (likely backed-mode indexing kicking in
# at this scale); matches why generate_geneformer_embeddings_chunked.py processes
# in ~5k-cell chunks. tokenize_data() iterates every h5ad in data_directory and
# concatenates results, so sharding here is a drop-in fix.
SHARD_SIZE = 5000
n_shards = int(np.ceil(adata.n_obs / SHARD_SIZE))
for s in range(n_shards):
    lo, hi = s * SHARD_SIZE, min((s + 1) * SHARD_SIZE, adata.n_obs)
    shard = adata[lo:hi].copy()
    shard.write_h5ad(f"{ready_dir}/monocyte_canonical_ready_shard{s:03d}.h5ad")
print(f"Wrote {n_shards} geneformer-ready shards ({SHARD_SIZE} cells each) to {ready_dir}", flush=True)

print("Tokenizing...", flush=True)
tk = TranscriptomeTokenizer(custom_attr_name_dict={"orig_index": "orig_index"}, nproc=8)
tk.tokenize_data(
    data_directory=ready_dir,
    output_directory=WORK,
    output_prefix="monocyte_canonical",
    file_format="h5ad",
)
tokenized = f"{WORK}/monocyte_canonical.dataset"
print(f"Tokenized dataset: {tokenized}", flush=True)

print("Extracting embeddings...", flush=True)
embex = EmbExtractor(
    model_type="Pretrained",
    emb_mode="cell",
    max_ncells=None,
    emb_layer=-1,
    forward_batch_size=32,
    nproc=8,
    emb_label=["orig_index"],
)
emb_df = embex.extract_embs(
    model_directory="ctheodoris/Geneformer",
    input_data_file=tokenized,
    output_directory=WORK,
    output_prefix="monocyte_canonical_emb",
)
print(f"Embedding table: {emb_df.shape}, columns include orig_index: "
      f"{'orig_index' in emb_df.columns}", flush=True)

label_cols = [c for c in ["orig_index"] if c in emb_df.columns]
emb_cols = [c for c in emb_df.columns if c not in label_cols]
emb_df = emb_df.set_index("orig_index")
emb_matrix = np.full((adata.n_obs, len(emb_cols)), np.nan, dtype=np.float32)
idx = emb_df.index.astype(int).values
emb_matrix[idx, :] = emb_df[emb_cols].values
n_missing = int(np.isnan(emb_matrix).any(axis=1).sum())
print(f"Cells missing an embedding after tokenizer filtering: {n_missing}/{adata.n_obs}", flush=True)

adata.obsm["X_geneformer"] = emb_matrix
adata.write_h5ad(OUTPUT)
print(f"Saved: {OUTPUT}", flush=True)

shutil.rmtree(WORK, ignore_errors=True)
print("DONE", flush=True)
