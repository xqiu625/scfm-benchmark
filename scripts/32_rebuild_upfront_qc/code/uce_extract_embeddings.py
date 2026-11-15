"""
Priority (b): UCE embedding extraction (third of 5 FMs, after Geneformer
and scGPT). Supports both the (superseded) subsampled input and, per the
2026-09-20 user override to run full-scale, sharded full-dataset input via
--shard-dir/--shard.

Unlike Geneformer/scGPT (hand-rolled tokenization+forward, reusing this
project's proven manual model-loading recipes), UCE uses its OFFICIAL,
UNMODIFIED CLI entry point (eval_single_anndata.py) via subprocess, run
from the UCE repo directory (relative default paths like
"./model_files/species_chrom.csv" require this cwd). This is deliberate,
not a shortcut: uce_masking_attribution.py's docstring documents a real bug
in UCE's released forward() (a float src_key_padding_mask is treated as an
additive attention bias instead of exclusion) that made GENE-MASKING
non-functional without a manual fix -- but for PLAIN embedding extraction
(no masking), the released, unmodified pipeline is the faithful choice: it
is what UCE's own published embeddings/benchmarks used, bugs and all, so
reimplementing around that quirk here would produce embeddings LESS
representative of the released model, not more correct ones.

Config matches uce_masking_attribution.py's verified setup (4-layer
checkpoint, d_model=1280) -- these also happen to be eval_single_anndata.py's
own defaults (nlayers=4, output_dim=1280, d_hid=5120, token_dim=5120),
confirmed via --help, so passed explicitly here only to avoid relying on
undocumented default resolution for --model_loc.

Real bug caught and fixed 2026-09-20: UCE's native sample_cell_sentences
does `np.random.choice(..., p=weights)` where weights are expression-based
sampling probabilities -- a cell whose total count RESTRICTED TO UCE's own
recognized gene set is exactly 0 (possible for cells with very low total
counts, e.g. some DecontX-decontaminated cells have fractional residual
sums as low as 0.03) makes every weight 0/0 = NaN, crashing the whole batch
with "probabilities contain NaN" (not a per-cell skip -- np.random.choice
raises immediately). Fixed by pre-filtering to cells with total raw count
sum >= MIN_TOTAL_COUNT before building UCE's input h5ad -- a UCE-specific
exclusion (Geneformer/scGPT embed these cells fine), so this script's
output can have FEWER cells than its input shard; the excluded cell_ids are
saved separately so this is traceable, not silently dropped.

Usage:
    python uce_extract_embeddings.py --which raw|decontx
    python uce_extract_embeddings.py --which raw|decontx --shard-dir <dir> --shard N
"""
import argparse
import os
import subprocess
import sys

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.io import mmread

SUBSAMPLE_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/fm_embeddings/subsampled_input'
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/fm_embeddings/uce'
UCE_REPO = '/bigdata/godziklab/shared/Xinru/302006/03_MODELS/uce/UCE'
MODEL_LOC = f'{UCE_REPO}/model_files/4layer_model.torch'
BATCH_SIZE = 25
MIN_TOTAL_COUNT = 1.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--which', required=True, choices=['raw', 'decontx'])
    ap.add_argument('--shard-dir', default=None)
    ap.add_argument('--shard', type=int, default=None)
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)

    if args.shard_dir is not None:
        assert args.shard is not None
        in_dir = f'{args.shard_dir}/{args.which}'
        mtx_path = f'{in_dir}/shard_{args.shard:03d}_counts.mtx'
        cells_path = f'{in_dir}/shard_{args.shard:03d}_cells.csv'
        genes_path = f'{args.shard_dir}/genes.csv'
        tag = f'{args.which}_shard{args.shard:03d}'
        out_emb = f'{OUT_DIR}/uce_emb_{tag}.npy'
        out_cells = f'{OUT_DIR}/uce_cells_{tag}.csv'
        out_excluded = f'{OUT_DIR}/uce_excluded_{tag}.csv'
        if os.path.exists(out_emb):
            print(f'shard {args.shard} already extracted, skipping', flush=True)
            return
    else:
        mtx_path = f'{SUBSAMPLE_DIR}/{args.which}_counts.mtx'
        cells_path = f'{SUBSAMPLE_DIR}/{args.which}_cells.csv'
        genes_path = f'{SUBSAMPLE_DIR}/genes.csv'
        tag = args.which
        out_emb = f'{OUT_DIR}/uce_emb_{tag}.npy'
        out_cells = f'{OUT_DIR}/uce_cells_{tag}.csv'
        out_excluded = f'{OUT_DIR}/uce_excluded_{tag}.csv'
    work_dir = f'{OUT_DIR}/{tag}_work/'
    os.makedirs(work_dir, exist_ok=True)

    print(f'Loading {args.which} counts.mtx + cells.csv + genes.csv...', flush=True)
    mat = mmread(mtx_path).tocsr()  # genes x cells
    cells = pd.read_csv(cells_path, low_memory=False)
    genes = pd.read_csv(genes_path)['gene'].tolist()
    assert mat.shape[1] == len(cells)
    print(f'  {mat.shape[1]} cells x {mat.shape[0]} genes, {cells.donor_id.nunique()} donors', flush=True)

    totals = np.asarray(mat.sum(axis=0)).ravel()
    keep = totals >= MIN_TOTAL_COUNT
    n_excluded = int((~keep).sum())
    if n_excluded:
        print(f'  excluding {n_excluded} cells with total count < {MIN_TOTAL_COUNT} '
              f'(UCE-specific: crashes np.random.choice with NaN weights otherwise)', flush=True)
        cells.loc[~keep, ['cell_id', 'donor_id', 'cohort']].to_csv(out_excluded, index=False)
    mat = mat[:, keep]
    cells = cells[keep].reset_index(drop=True)

    print('Writing input h5ad for UCE\'s native CLI...', flush=True)
    adata = AnnData(X=mat.T.tocsr(), obs=pd.DataFrame(index=cells['cell_id'].astype(str)),
                     var=pd.DataFrame(index=genes))
    input_path = f'{work_dir}input_{tag}.h5ad'
    adata.write(input_path)

    cmd = [
        'python3', 'eval_single_anndata.py',
        '--adata_path', input_path,
        '--dir', work_dir,
        '--species', 'human',
        '--batch_size', str(BATCH_SIZE),
        '--nlayers', '4',
        '--output_dim', '1280',
        '--d_hid', '5120',
        '--token_dim', '5120',
        '--model_loc', MODEL_LOC,
    ]
    print(f'Running UCE native CLI: {" ".join(cmd)}', flush=True)
    subprocess.run(cmd, cwd=UCE_REPO, check=True)

    out_h5ad = f'{work_dir}input_{tag}_uce_adata.h5ad'
    print(f'Loading UCE output from {out_h5ad}...', flush=True)
    import scanpy as sc
    result = sc.read_h5ad(out_h5ad)
    emb = np.asarray(result.obsm['X_uce'], dtype=np.float32)
    assert emb.shape[0] == len(cells), f'UCE returned {emb.shape[0]} cells, expected {len(cells)}'
    print(f'  embeddings shape: {emb.shape}', flush=True)

    np.save(out_emb, emb)
    cells.to_csv(out_cells, index=False)
    print(f'\nSAVED -> {out_emb} ({out_cells} carries matching row order + labels; '
          f'{n_excluded} cells excluded, see {out_excluded if n_excluded else "n/a"})')


if __name__ == '__main__':
    main()
