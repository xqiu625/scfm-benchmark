"""
Priority (b): scGPT embedding extraction on the shared 8-cohort subsample
(second of 5 FMs, after Geneformer).

Model loading mirrors scgpt_masking_attribution.py's proven recipe (job
28277291, confirmed working, 42m38s for a full-scale run) rather than
scGPT's own high-level `embed_data` helper: `embed_data` internally calls
`torch.load(model_file, map_location=device)` WITHOUT `weights_only=False`,
which the masking script's docstring flags as a real risk on PyTorch >=2.6's
weights_only=True default -- that script worked around it by loading the
model manually instead of trusting embed_data. Reused here for the same
reason, not re-risked.

Architecture facts (verified against the vendored bowang-lab/scGPT repo,
03_MODELS/scgpt/scGPT, same verification as scgpt_masking_attribution.py):
  - Pooling is CLS: out[:, 0, :] from model._encode(...).
  - use_fast_transformer=True in the config falls back to standard
    (non-flash) attention automatically since flash-attn is not installed
    in scpeft_env -- confirmed safe, not a masking-specific concern.
  - Binning (scgpt.preprocess.binning) is applied to each cell's FULL
    nonzero expression profile BEFORE gene selection, matching native
    scGPT usage (computing bin edges over a restricted gene subset would
    give different quantiles than the model saw in training/its own
    embed_data path).
  - max_length=1200 (scGPT's own default in get_batch_cell_embeddings) --
    cells with more nonzero genes than this are truncated to the
    top-(max_length-1) by expression (leaving room for <cls>), matching
    Geneformer's rank-value truncation convention used elsewhere in this
    project for consistency across FMs, though scGPT's own reference
    implementation does not specify a canonical truncation order for
    over-length cells.

Usage:
    python scgpt_extract_embeddings.py --which raw|decontx
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import torch
from scipy.io import mmread

SCGPT_REPO = '/bigdata/godziklab/shared/Xinru/302006/03_MODELS/scgpt/scGPT'
sys.path.insert(0, SCGPT_REPO)
from scgpt.model import TransformerModel
from scgpt.tokenizer import GeneVocab
from scgpt.utils import load_pretrained
from scgpt.preprocess import binning

SUBSAMPLE_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/fm_embeddings/subsampled_input'
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/fm_embeddings/scgpt'
MODEL_DIR = '/bigdata/godziklab/shared/Xinru/302006/03_MODELS/scgpt/scGPT_human/scGPT_human'
MAX_LENGTH = 1200
BATCH_CELLS = 32
N_BINS = 51


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
        out_emb = f'{OUT_DIR}/scgpt_emb_{args.which}_shard{args.shard:03d}.npy'
        out_cells = f'{OUT_DIR}/scgpt_cells_{args.which}_shard{args.shard:03d}.csv'
        if os.path.exists(out_emb):
            print(f'shard {args.shard} already extracted, skipping', flush=True)
            return
    else:
        mtx_path = f'{SUBSAMPLE_DIR}/{args.which}_counts.mtx'
        cells_path = f'{SUBSAMPLE_DIR}/{args.which}_cells.csv'
        genes_path = f'{SUBSAMPLE_DIR}/genes.csv'
        out_emb = f'{OUT_DIR}/scgpt_emb_{args.which}.npy'
        out_cells = f'{OUT_DIR}/scgpt_cells_{args.which}.csv'

    print(f'Loading {args.which} counts.mtx + cells.csv + genes.csv...', flush=True)
    mat = mmread(mtx_path).tocsr()  # genes x cells -> csr for row slicing after T
    cells = pd.read_csv(cells_path, low_memory=False)
    genes = pd.read_csv(genes_path)['gene'].tolist()
    assert mat.shape[1] == len(cells)
    print(f'  {mat.shape[1]} cells x {mat.shape[0]} genes, {cells.donor_id.nunique()} donors', flush=True)
    X = mat.T.tocsr()  # cells x genes

    print('Loading scGPT model + vocab...', flush=True)
    vocab_file = f'{MODEL_DIR}/vocab.json'
    model_config_file = f'{MODEL_DIR}/args.json'
    model_file = f'{MODEL_DIR}/best_model.pt'
    pad_token = '<pad>'
    special_tokens = [pad_token, '<cls>', '<eoc>']

    vocab = GeneVocab.from_file(vocab_file)
    for s in special_tokens:
        if s not in vocab:
            vocab.append_token(s)
    vocab.set_default_index(vocab['<pad>'])

    with open(model_config_file) as f:
        model_configs = json.load(f)

    model = TransformerModel(
        ntoken=len(vocab), d_model=model_configs['embsize'], nhead=model_configs['nheads'],
        d_hid=model_configs['d_hid'], nlayers=model_configs['nlayers'],
        nlayers_cls=model_configs['n_layers_cls'], n_cls=1, vocab=vocab,
        dropout=model_configs['dropout'], pad_token=model_configs['pad_token'],
        pad_value=model_configs['pad_value'], do_mvc=True, do_dab=False,
        use_batch_labels=False, domain_spec_batchnorm=False, explicit_zero_prob=False,
        use_fast_transformer=True, fast_transformer_backend='flash', pre_norm=False,
    )
    ckpt = torch.load(model_file, map_location=device, weights_only=False)
    load_pretrained(model, ckpt, verbose=False)
    model.to(device).eval()
    print(f'  Model loaded. use_fast_transformer={model.use_fast_transformer} '
          f'(expect False -- flash-attn not installed, safe fallback)', flush=True)

    cls_id = vocab['<cls>']
    pad_id = vocab[pad_token]
    pad_value = model_configs['pad_value']

    col_to_scgpt_id = {}
    n_matched = 0
    for i, g in enumerate(genes):
        tok = vocab[g]
        if tok != pad_id:
            col_to_scgpt_id[i] = tok
            n_matched += 1
    print(f'  {n_matched} / {len(genes)} genes matched scGPT vocab', flush=True)
    assert n_matched > 5000, 'too few genes matched -- vocab lookup likely broken'

    def build_sentence(cell_idx):
        row = X[cell_idx].toarray().ravel()
        nz = np.nonzero(row)[0]
        nz = np.array([c for c in nz if c in col_to_scgpt_id])
        if len(nz) == 0:
            return np.array([cls_id], dtype=np.int64), np.array([pad_value], dtype=np.float32)
        binned = binning(row[nz].astype(np.float32), n_bins=N_BINS)
        order = np.argsort(-row[nz])  # highest expression first, for truncation
        nz, binned = nz[order], np.asarray(binned)[order]
        if len(nz) > MAX_LENGTH - 1:
            nz, binned = nz[:MAX_LENGTH - 1], binned[:MAX_LENGTH - 1]
        gene_ids = [cls_id] + [col_to_scgpt_id[c] for c in nz]
        values = [float(pad_value)] + binned.tolist()
        return np.array(gene_ids, dtype=np.int64), np.array(values, dtype=np.float32)

    print('Extracting scGPT embeddings (CLS pooling)...', flush=True)
    n_cells = X.shape[0]
    pooled = []
    with torch.no_grad():
        for start in range(0, n_cells, BATCH_CELLS):
            end = min(start + BATCH_CELLS, n_cells)
            genes_list, values_list = [], []
            for ci in range(start, end):
                g, v = build_sentence(ci)
                genes_list.append(g)
                values_list.append(v)
            max_len = max(len(g) for g in genes_list)
            genes_2d = np.full((end - start, max_len), pad_id, dtype=np.int64)
            values_2d = np.full((end - start, max_len), float(pad_value), dtype=np.float32)
            for j, (g, v) in enumerate(zip(genes_list, values_list)):
                genes_2d[j, :len(g)] = g
                values_2d[j, :len(v)] = v
            genes_t = torch.tensor(genes_2d, device=device)
            values_t = torch.tensor(values_2d, device=device)
            mask_t = genes_t.eq(pad_id)
            out = model._encode(genes_t, values_t, src_key_padding_mask=mask_t)
            pooled.append(out[:, 0, :].cpu().numpy())
            if (start // BATCH_CELLS) % 50 == 0:
                print(f'    {end}/{n_cells} cells embedded', flush=True)

    emb = np.concatenate(pooled, axis=0).astype(np.float32)
    print(f'  embeddings shape: {emb.shape}', flush=True)

    np.save(out_emb, emb)
    cells.to_csv(out_cells, index=False)
    print(f'\nSAVED -> {out_emb} ({out_cells} carries matching row order + labels)')


if __name__ == '__main__':
    main()
