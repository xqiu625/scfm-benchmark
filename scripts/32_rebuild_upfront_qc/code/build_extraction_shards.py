"""
Priority (b) REVISED (user override 2026-09-20 -- explicit instruction to
NOT subsample, feed the FULL Platelet+Monocyte cell set to all 5 FMs;
supersedes Cindy's 2026-09-19 subsampling decision for this benchmark run,
see qc_caveats.md for the full record). Splits the full defined-event
Platelet+Monocyte dataset (1,186,293 cells, 612 donors -- ~96% of the
combined dataset's cells, since only 50,588 have an undefined outcome) into
N contiguous shards per FM extraction script to run as checkpointed
short_gpu array jobs -- a single job cannot hold ~1.19M cells within
short_gpu's 2h cap (the 61,441-cell subsample alone took Geneformer ~20 min,
so full scale is roughly 19x that per model).

Single-pass sharding, NOT N calls to read_mtx_columns: since ~96% of the
file's columns are kept here (unlike the Platelet/Monocyte cell-type
filtering step, which was a large reduction), calling the column-filtering
reader once per shard would re-read the entire 22-50G file N times.
Instead this streams the source counts.mtx exactly ONCE, bucketing each
chunk's rows into per-shard body-only text files (no header -- nnz isn't
known until the pass completes) as it goes, then prepends the correct
3-line MatrixMarket header to each shard's body afterward. Memory is
bounded by one chunk (20M rows), never the whole file or a whole shard.

Usage:
    python build_extraction_shards.py --which raw|decontx --n-shards 12
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

DATA_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/platelet_monocyte_benchmark_input_v2'
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/fm_embeddings/full_shards'
CHUNK_SIZE = 20_000_000


def stream_shard_mtx(src_path, col_to_shard, col_rank_in_shard, n_shards, shard_dir, n_row):
    body_paths = [f'{shard_dir}/.shard_{i:03d}.body' for i in range(n_shards)]
    handles = [open(p, 'w') for p in body_paths]
    nnz_counts = np.zeros(n_shards, dtype=np.int64)

    with open(src_path) as f:
        line = f.readline()
        assert line.startswith('%%MatrixMarket')
        n_header = 1
        line = f.readline()
        n_header += 1
        while line.startswith('%'):
            line = f.readline()
            n_header += 1
        src_n_row, src_n_col, nnz_total = (int(x) for x in line.split())
    assert src_n_row == n_row

    reader = pd.read_csv(src_path, sep=r'\s+', skiprows=n_header, header=None,
                          names=['i', 'j', 'v'], chunksize=CHUNK_SIZE,
                          dtype={'i': np.int64, 'j': np.int64, 'v': np.float64})
    n_seen = 0
    for chunk in reader:
        j0 = chunk['j'].to_numpy() - 1
        shard_id = col_to_shard[j0]
        keep = shard_id >= 0
        if keep.any():
            sub = chunk[keep].copy()
            sub['j'] = col_rank_in_shard[j0[keep]] + 1  # 1-indexed for MatrixMarket
            sub['shard'] = shard_id[keep]
            # vectorized per-shard write (pandas' C-level to_csv, not a Python
            # per-row loop -- a manual loop here would take hours at this scale)
            for s, grp in sub.groupby('shard'):
                grp[['i', 'j', 'v']].to_csv(handles[s], sep=' ', header=False, index=False)
                nnz_counts[s] += len(grp)
        n_seen += len(chunk)
        print(f'    {n_seen}/{nnz_total} entries streamed', flush=True)

    for h in handles:
        h.close()
    return nnz_counts, body_paths


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--which', required=True, choices=['raw', 'decontx'])
    ap.add_argument('--n-shards', type=int, required=True)
    args = ap.parse_args()
    shard_dir = f'{OUT_DIR}/{args.which}'
    os.makedirs(shard_dir, exist_ok=True)

    if os.path.exists(f'{shard_dir}/ALL_DONE'):
        print('All shards already built, skipping.', flush=True)
        return

    print(f'Loading {args.which} cells.csv...', flush=True)
    cells_full = pd.read_csv(f'{DATA_DIR}/{args.which}/cells.csv', low_memory=False)
    n_total = len(cells_full)
    keep = cells_full['event'].notna()
    idx_all = np.where(keep.to_numpy())[0]
    print(f'  {len(idx_all)} / {n_total} cells have a defined event (Platelet+Monocyte_Myeloid)',
          flush=True)

    shard_bounds = np.array_split(idx_all, args.n_shards)
    n_col_total = len(cells_full)
    col_to_shard = np.full(n_col_total, -1, dtype=np.int64)
    col_rank_in_shard = np.zeros(n_col_total, dtype=np.int64)
    for s, shard_idx in enumerate(shard_bounds):
        col_to_shard[shard_idx] = s
        col_rank_in_shard[shard_idx] = np.arange(len(shard_idx))
        cells_full.iloc[shard_idx].reset_index(drop=True).to_csv(
            f'{shard_dir}/shard_{s:03d}_cells.csv', index=False)

    genes = pd.read_csv(f'{DATA_DIR}/genes.csv')['gene']
    n_genes = len(genes)
    genes.to_csv(f'{OUT_DIR}/genes.csv', index=False)

    print(f'Streaming {args.which} counts.mtx ONCE, bucketing into {args.n_shards} shards...',
          flush=True)
    nnz_counts, body_paths = stream_shard_mtx(
        f'{DATA_DIR}/{args.which}/counts.mtx', col_to_shard, col_rank_in_shard,
        args.n_shards, shard_dir, n_genes)

    for s, shard_idx in enumerate(shard_bounds):
        header_path = f'{shard_dir}/.shard_{s:03d}.header'
        with open(header_path, 'w') as f:
            f.write('%%MatrixMarket matrix coordinate real general\n%\n')
            f.write(f'{n_genes} {len(shard_idx)} {nnz_counts[s]}\n')
        final_path = f'{shard_dir}/shard_{s:03d}_counts.mtx'
        os.system(f'cat {header_path} {body_paths[s]} > {final_path}.tmp && mv {final_path}.tmp {final_path}')
        os.remove(header_path)
        os.remove(body_paths[s])
        with open(f'{shard_dir}/shard_{s:03d}_DONE', 'w') as f:
            f.write('done\n')
        print(f'  shard {s}: {len(shard_idx)} cells, {nnz_counts[s]} nnz -> {final_path}', flush=True)

    with open(f'{shard_dir}/ALL_DONE', 'w') as f:
        f.write('done\n')
    print(f'\nDone. {args.n_shards} shards in {shard_dir}')


if __name__ == '__main__':
    main()
