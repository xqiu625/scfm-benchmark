"""
Priority (b): Geneformer embedding extraction on the shared 8-cohort
subsample (subsample_platelet_monocyte_cells.py output), per Cindy's
2026-09-19 subsampling decision. Pure embedding extraction -- no LOCO/CV
here (that is priority (c), run later against these saved embeddings so
every evaluation reuses the same extraction rather than re-embedding).

Reuses the confirmed-working native-BertModel setup validated in Phase 3
(geneformer_loco_before_after_decontx.py, reference_302006_envs memory):
scpeft_env, geneformer-V2-316M-native checkpoint, token_dictionary_gc104M.pkl
(Ensembl-keyed -- requires the symbol->Ensembl mapping first), rank-value
tokenization (max_len=512), mean-pooled last_hidden_state per cell.

Usage:
    python geneformer_extract_embeddings.py --which raw|decontx
"""
import argparse
import os
import pickle
import sys

import numpy as np
import pandas as pd
import torch
from scipy.io import mmread
from torch.utils.data import Dataset, DataLoader
from transformers import BertModel

SUBSAMPLE_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/fm_embeddings/subsampled_input'
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/fm_embeddings/geneformer'
MODEL_PATH = '/bigdata/godziklab/shared/Xinru/302006/03_MODELS/geneformer/geneformer-V2-316M-native'
TOKEN_DICT = (MODEL_PATH + '/geneformer/token_dictionary_gc104M.pkl')
SYMBOL_TO_ENSEMBL_CSV = '/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c_monocyte/gene_symbol_to_ensembl.csv'
MAX_LEN = 512
BATCH_SIZE = 32


class CellDataset(Dataset):
    def __init__(self, X, idx_to_token, max_len=MAX_LEN):
        self.X = X  # cells x genes, csr
        self.idx_to_token = idx_to_token
        self.max_len = max_len

    def __len__(self):
        return self.X.shape[0]

    def __getitem__(self, i):
        row = self.X[i].toarray().ravel()
        nz = np.where(row > 0)[0]
        if len(nz) == 0:
            return {'input_ids': torch.zeros(self.max_len, dtype=torch.long),
                    'attention_mask': torch.zeros(self.max_len, dtype=torch.long)}
        order = nz[np.argsort(-row[nz])]
        toks = []
        for g in order:
            t = self.idx_to_token.get(int(g))
            if t is not None:
                toks.append(t)
                if len(toks) >= self.max_len:
                    break
        n = len(toks)
        if n < self.max_len:
            toks = toks + [0] * (self.max_len - n)
            mask = [1] * n + [0] * (self.max_len - n)
        else:
            mask = [1] * self.max_len
        return {'input_ids': torch.tensor(toks[:self.max_len], dtype=torch.long),
                'attention_mask': torch.tensor(mask[:self.max_len], dtype=torch.long)}


def get_embeddings(mat, idx_to_token, device):
    ds = CellDataset(mat, idx_to_token)
    loader = DataLoader(ds, batch_size=BATCH_SIZE, shuffle=False, num_workers=4, pin_memory=True)
    model = BertModel.from_pretrained(MODEL_PATH).to(device).eval()
    all_embs = []
    with torch.no_grad():
        for i, batch in enumerate(loader):
            ids = batch['input_ids'].to(device)
            mask = batch['attention_mask'].to(device)
            out = model(input_ids=ids, attention_mask=mask).last_hidden_state
            m = mask.unsqueeze(-1).float()
            pooled = (out * m).sum(dim=1) / m.sum(dim=1).clamp(min=1)
            all_embs.append(pooled.cpu().numpy())
            if (i + 1) % 50 == 0:
                print(f'    {(i+1) * BATCH_SIZE} cells embedded', flush=True)
    return np.concatenate(all_embs, axis=0).astype(np.float32)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--which', required=True, choices=['raw', 'decontx'])
    ap.add_argument('--shard-dir', default=None,
                     help='if set, read shard_{shard:03d}_{counts.mtx,cells.csv} from here '
                          '(full-scale mode) instead of the subsampled input')
    ap.add_argument('--shard', type=int, default=None)
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f'Device: {device}', flush=True)

    if args.shard_dir is not None:
        assert args.shard is not None
        in_dir = f'{args.shard_dir}/{args.which}'
        mtx_path = f'{in_dir}/shard_{args.shard:03d}_counts.mtx'
        cells_path = f'{in_dir}/shard_{args.shard:03d}_cells.csv'
        genes_path = f'{args.shard_dir}/genes.csv'
        out_emb = f'{OUT_DIR}/geneformer_emb_{args.which}_shard{args.shard:03d}.npy'
        out_cells = f'{OUT_DIR}/geneformer_cells_{args.which}_shard{args.shard:03d}.csv'
        if os.path.exists(out_emb):
            print(f'shard {args.shard} already extracted, skipping', flush=True)
            return
    else:
        mtx_path = f'{SUBSAMPLE_DIR}/{args.which}_counts.mtx'
        cells_path = f'{SUBSAMPLE_DIR}/{args.which}_cells.csv'
        genes_path = f'{SUBSAMPLE_DIR}/genes.csv'
        out_emb = f'{OUT_DIR}/geneformer_emb_{args.which}.npy'
        out_cells = f'{OUT_DIR}/geneformer_cells_{args.which}.csv'

    print(f'Loading {args.which} counts.mtx + cells.csv + genes.csv...', flush=True)
    mat = mmread(mtx_path).tocsc()  # genes x cells
    cells = pd.read_csv(cells_path, low_memory=False)
    genes = pd.read_csv(genes_path)['gene'].tolist()
    assert mat.shape[1] == len(cells)
    print(f'  {mat.shape[1]} cells x {mat.shape[0]} genes, {cells.donor_id.nunique()} donors', flush=True)
    mat = mat.T.tocsr()  # cells x genes, for the Dataset

    print('Mapping gene symbols -> Ensembl -> Geneformer tokens...', flush=True)
    sym2ens = pd.read_csv(SYMBOL_TO_ENSEMBL_CSV).set_index('gene_symbol')['ensembl_id'].to_dict()
    with open(TOKEN_DICT, 'rb') as f:
        token_dict = pickle.load(f)
    idx_to_token = {}
    for i, g in enumerate(genes):
        ens = sym2ens.get(g)
        if ens is not None and ens in token_dict:
            idx_to_token[i] = token_dict[ens]
    print(f'  {len(idx_to_token)} / {len(genes)} genes mapped to Geneformer tokens', flush=True)
    assert len(idx_to_token) > 5000, 'too few genes mapped -- symbol/Ensembl mapping likely broken'

    print('Extracting Geneformer embeddings...', flush=True)
    emb = get_embeddings(mat, idx_to_token, device)
    print(f'  embeddings shape: {emb.shape}', flush=True)

    np.save(out_emb, emb)
    cells.to_csv(out_cells, index=False)
    print(f'\nSAVED -> {out_emb} ({out_cells} carries matching row order + labels)')


if __name__ == '__main__':
    main()
