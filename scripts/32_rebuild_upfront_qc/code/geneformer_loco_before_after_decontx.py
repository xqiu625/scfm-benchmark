"""
Before-vs-after-decontamination robustness rerun, step 2 of 2 (step 1 was
PCA-200, see pca200_loco_before_after_decontx.py): Geneformer embeddings on
the same within-COVID validation set, run through the same LOCO harness, on
BOTH the raw and decontaminated matrices.

Uses the confirmed-working native-BertModel setup for Geneformer inference
(reference_302006_envs memory, corrected 2026-09-10; NOT the older
jobS_embed/geneformer_chongqing.sh path, which uses a since-superseded token
dictionary that CUDA-crashes on this checkpoint):
  - scpeft_env (native transformers.BertModel, not Geneformer's own wrapper lib)
  - MODEL_PATH = 03_MODELS/geneformer/geneformer-V2-316M-native
  - token_dictionary_gc104M.pkl (Ensembl-keyed, NOT gene-symbol-keyed -- our
    gene lists are symbols, so a symbol->Ensembl mapping is required first;
    same gene_symbol_to_ensembl.csv used by masking_attribution_common.py)

Stratified per-donor subsample (cap 60 cells/donor, matching the established
convention in geneformer_gene_attribution.py) for GPU tractability -- running
Geneformer on the full 2.4M cells is not what the PCA-200 validation step
needed either; this is a per-cell-embedding-then-patient-aggregation model,
same as PCA-200's design, just with a capped per-donor cell budget instead of
using every cell.

Usage:
    python geneformer_loco_before_after_decontx.py --which raw
    python geneformer_loco_before_after_decontx.py --which decontx
"""
import argparse
import os
import pickle
import sys

import numpy as np
import pandas as pd
import torch
import scipy.sparse as sp
from scipy.io import mmread
from sklearn.linear_model import LogisticRegression
from torch.utils.data import Dataset, DataLoader
from transformers import BertModel

sys.path.insert(0, '/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/code')
from repeated_cv import run_loco, summarize_loco

VALIDATION_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/within_covid_validation'
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/geneformer_loco_validation'
MODEL_PATH = '/bigdata/godziklab/shared/Xinru/302006/03_MODELS/geneformer/geneformer-V2-316M-native'
TOKEN_DICT = (MODEL_PATH + '/geneformer/token_dictionary_gc104M.pkl')
SYMBOL_TO_ENSEMBL_CSV = '/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c_monocyte/gene_symbol_to_ensembl.csv'
CELLS_PER_DONOR = 60
RNG_SEED = 20260917
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


def stratified_subsample(cells_df, cap, seed):
    rng = np.random.default_rng(seed)
    keep = []
    for donor, grp in cells_df.groupby('donor_id'):
        n = len(grp)
        take = min(n, cap)
        chosen = rng.choice(grp.index.to_numpy(), size=take, replace=False)
        keep.extend(chosen.tolist())
    return np.sort(np.array(keep))


def get_embeddings(mat, genes, idx_to_token, device):
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


def fit_predict_factory():
    def fit_predict(X_tr, y_tr, X_te):
        clf = LogisticRegression(C=1.0, solver='lbfgs', max_iter=2000,
                                  class_weight='balanced', n_jobs=4)
        clf.fit(X_tr, y_tr)
        p = clf.predict_proba(X_te)[:, 1]
        return {'p_lr': p}
    return fit_predict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--which', required=True, choices=['raw', 'decontx'])
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f'Device: {device}', flush=True)

    print(f'Loading {args.which} counts.mtx + cells.csv + genes.csv...', flush=True)
    mat_full = mmread(f'{VALIDATION_DIR}/{args.which}/counts.mtx').tocsc()  # genes x cells
    cells = pd.read_csv(f'{VALIDATION_DIR}/{args.which}/cells.csv')
    genes = pd.read_csv(f'{VALIDATION_DIR}/genes.csv')['gene'].tolist()
    assert mat_full.shape[1] == len(cells)
    print(f'  {mat_full.shape[1]} cells x {mat_full.shape[0]} genes, {cells.donor_id.nunique()} donors', flush=True)

    sub_idx = stratified_subsample(cells, CELLS_PER_DONOR, RNG_SEED)
    cells = cells.iloc[sub_idx].reset_index(drop=True)
    mat = mat_full[:, sub_idx].T.tocsr()  # cells x genes, for the Dataset
    print(f'  after stratified subsample (cap={CELLS_PER_DONOR}/donor): {mat.shape[0]} cells', flush=True)

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
    emb = get_embeddings(mat, genes, idx_to_token, device)
    print(f'  embeddings shape: {emb.shape}', flush=True)
    np.save(f'{OUT_DIR}/geneformer_emb_{args.which}.npy', emb)
    cells.to_csv(f'{OUT_DIR}/geneformer_cells_{args.which}.csv', index=False)

    y = cells['event'].to_numpy().astype(int)
    patient_id = cells['donor_id'].to_numpy()
    cohort = cells['cohort'].to_numpy()

    print('Running LOCO (Geneformer embeddings + LR, fit on train folds only)...', flush=True)
    oof, held_out = run_loco(emb, y, cohort, fit_predict_factory())
    print(f'  held out cohorts in order: {held_out}', flush=True)

    report = summarize_loco(oof, patient_id, y, cohort, score_key='p_lr')
    report.to_csv(f'{OUT_DIR}/loco_report_{args.which}.csv', index=False)
    print(f'\n=== {args.which} Geneformer LOCO report ===')
    print(report.to_string(index=False))
    print(f'\nSAVED -> {OUT_DIR}/loco_report_{args.which}.csv')


if __name__ == '__main__':
    main()
