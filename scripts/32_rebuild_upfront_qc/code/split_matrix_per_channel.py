"""
Split an already-built cohort counts.mtx/genes.csv/cells.csv (the
run_decontx_per_channel.R IO_DIR contract) into one small mtx file per
channel. Needed for cohorts whose combined matrix exceeds R/Matrix's 32-bit
non-zero-count ceiling (2,147,483,647) -- hit for real on Ren/GSE158055
(2,319,108,599 non-zeros): celda's own internal per-batch subsetting
(`x[, cellsInBatch]`) inside decontX(..., batch=) fails on a TsparseMatrix
whose 'i'/'j' slots exceed 2^31-1 in length, regardless of representation
(CSC's cumulative "p" slot fails even earlier, at construction). The only
robust fix is to never build one unified R Matrix object at that scale --
split per-channel here in Python (scipy has no such limitation) and run
DecontX per-channel directly in R (run_decontx_per_channel_looped.R),
bypassing celda's internal batch-splitting for oversized cohorts.

Usage:
    python split_matrix_per_channel.py --io-dir <dir with counts.mtx, genes.csv, cells.csv> \
        --out-dir <output dir for per-channel subdirectories>

Output: <out_dir>/<channel_id>/counts.mtx, cells.csv (that channel's rows of
the original cells.csv, same columns) per channel. genes.csv is not
duplicated -- callers read it once from --io-dir.
"""
import argparse
import os

import pandas as pd
from scipy.io import mmread, mmwrite


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--io-dir', required=True)
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    print('Loading cells.csv...', flush=True)
    cells = pd.read_csv(os.path.join(args.io_dir, 'cells.csv'))
    print(f'  {len(cells)} cells, {cells["channel_id"].nunique()} channels', flush=True)

    print('Loading counts.mtx (genes x cells)...', flush=True)
    mat = mmread(os.path.join(args.io_dir, 'counts.mtx')).tocsc()
    assert mat.shape[1] == len(cells), f'matrix has {mat.shape[1]} columns, cells.csv has {len(cells)} rows'

    channel_ids = cells['channel_id'].to_numpy()
    unique_channels = sorted(set(channel_ids))
    print(f'Splitting into {len(unique_channels)} per-channel files...', flush=True)

    manifest_rows = []
    for i, ch in enumerate(unique_channels):
        idx = (channel_ids == ch).nonzero()[0]
        sub = mat[:, idx]
        ch_dir = os.path.join(args.out_dir, ch)
        os.makedirs(ch_dir, exist_ok=True)
        mmwrite(os.path.join(ch_dir, 'counts.mtx'), sub)
        cells.iloc[idx].to_csv(os.path.join(ch_dir, 'cells.csv'), index=False)
        manifest_rows.append({'channel_id': ch, 'n_cells': len(idx), 'nnz': sub.nnz})
        if (i + 1) % 40 == 0 or i == len(unique_channels) - 1:
            print(f'  {i+1}/{len(unique_channels)} channels written', flush=True)

    pd.DataFrame(manifest_rows).to_csv(os.path.join(args.out_dir, 'split_manifest.csv'), index=False)
    print(f'\nDone. {len(unique_channels)} channels written to {args.out_dir}')
    print(f'Max per-channel nnz: {max(r["nnz"] for r in manifest_rows):,} '
          f'(all must stay well under 2,147,483,647 for the R side to work)')


if __name__ == '__main__':
    main()
