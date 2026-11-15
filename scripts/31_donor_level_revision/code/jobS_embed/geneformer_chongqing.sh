#!/bin/bash
#SBATCH --job-name=jobS_geneformer
#SBATCH --partition=short_gpu
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu04,gpu13,gpu14
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=2:00:00
#SBATCH --requeue
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/jobS_geneformer_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/jobS_geneformer_%j.err

set -euo pipefail
echo "=== Geneformer Chongqing (native scpeft_env, scPEFT tokenization) ==="
date

ENV_DIR=/bigdata/godziklab/shared/Xinru/302006/scpeft_env
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH=/bigdata/godziklab/shared/Xinru/302006/03_MODELS/scPEFT:${PYTHONPATH:-}
export PYTHONNOUSERSITE=1

INPUT=/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/chongqing_sepsis_platelet_raw_counts.h5ad
OUTPUT=/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c_chongqing/geneformer_chongqing_canonical.h5ad
MODEL_PATH=/bigdata/godziklab/shared/Xinru/302006/03_MODELS/geneformer/geneformer-V2-316M-native
TOKEN_DICT=/bigdata/godziklab/shared/Xinru/302006/03_MODELS/scPEFT/geneformer_peft/geneformer/token_dictionary.pkl
mkdir -p "$(dirname "$OUTPUT")"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

python - <<PY
import os, sys, pickle
import numpy as np, pandas as pd, scanpy as sc
import scipy.sparse as sp
import torch
from torch.utils.data import Dataset, DataLoader
from transformers import BertModel

INPUT = "$INPUT"
OUTPUT = "$OUTPUT"
MODEL_PATH = "$MODEL_PATH"
TOKEN_DICT = "$TOKEN_DICT"

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Device: {device}")

print(f"Loading: {INPUT}")
a = sc.read_h5ad(INPUT)
print(f"  cells={a.n_obs}, genes={a.n_vars}, var sample: {list(a.var_names[:5])}")

with open(TOKEN_DICT, "rb") as f:
    token_dict = pickle.load(f)
print(f"  Token dict: {len(token_dict)} entries")

model = BertModel.from_pretrained(MODEL_PATH).to(device).eval()
hidden_size = model.config.hidden_size
vocab_size = model.config.vocab_size
print(f"  Hidden size: {hidden_size}, Vocab size: {vocab_size}")

idx_to_token = {}
n_mapped = 0
for i, gene in enumerate(a.var_names.astype(str)):
    if gene in token_dict:
        tok = token_dict[gene]
        if tok < vocab_size:
            idx_to_token[i] = tok
            n_mapped += 1
print(f"Genes mapped to tokens: {n_mapped} / {a.n_vars}")
if n_mapped < 100:
    print("ERROR: too few genes mapped -- check var format")
    sys.exit(1)

class CellDataset(Dataset):
    def __init__(self, X, max_len=512):
        self.X = X
        self.max_len = max_len
    def __len__(self):
        return self.X.shape[0]
    def __getitem__(self, idx):
        row = self.X[idx].toarray().ravel() if sp.issparse(self.X) else np.asarray(self.X[idx]).ravel()
        nz = np.where(row > 0)[0]
        if len(nz) == 0:
            return {"input_ids": torch.zeros(self.max_len, dtype=torch.long),
                    "attention_mask": torch.zeros(self.max_len, dtype=torch.long)}
        order = nz[np.argsort(-row[nz])]
        toks = []
        for i in order:
            t = idx_to_token.get(int(i))
            if t is not None:
                toks.append(t)
                if len(toks) >= self.max_len:
                    break
        n = len(toks)
        if n < self.max_len:
            toks = toks + [0] * (self.max_len - n)
            mask = [1]*n + [0]*(self.max_len - n)
        else:
            mask = [1] * self.max_len
        return {"input_ids": torch.tensor(toks[:self.max_len], dtype=torch.long),
                "attention_mask": torch.tensor(mask[:self.max_len], dtype=torch.long)}

ds = CellDataset(a.X, max_len=512)
loader = DataLoader(ds, batch_size=32, shuffle=False, num_workers=4, pin_memory=True)
print(f"Extracting embeddings on {len(ds)} cells...")

all_embs = []
with torch.no_grad():
    for i, batch in enumerate(loader):
        ids = batch["input_ids"].to(device)
        mask = batch["attention_mask"].to(device)
        out = model(input_ids=ids, attention_mask=mask).last_hidden_state
        m = mask.unsqueeze(-1).float()
        pooled = (out * m).sum(dim=1) / m.sum(dim=1).clamp(min=1)
        all_embs.append(pooled.cpu().numpy())
        if (i+1) % 100 == 0:
            print(f"  Processed {(i+1) * loader.batch_size} cells")

emb = np.concatenate(all_embs, axis=0).astype(np.float32)
print(f"Embeddings shape: {emb.shape}")

new = a.copy()
new.obsm["X_geneformer"] = emb
new.write_h5ad(OUTPUT, compression="gzip")
print(f"Wrote: {OUTPUT}")
PY

ls -lh "$OUTPUT"
date
