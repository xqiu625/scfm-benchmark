"""
Assemble Liu/GSE161918's counts.mtx / genes.csv / cells.csv for
run_decontx_per_channel.R, combining three things per lane:
  1. filter_liu_demuxlet_singlets.py's per-lane SNG barcode + donor_id
     assignment (<lane>_singlets.csv)
  2. Liu's per-lane RNA filtered_feature_bc_matrix.h5 (cellranger's own
     called-cell set for that GEM well)
  3. harmonized_lineage_labels.py's marker-score z-label, computed
     independently per lane (never pooled across lanes, per the
     project_302006_rebuild_architecture no-cross-channel-fitting rule)

A cell enters the final matrix only if it is BOTH a demuxlet singlet AND one
of cellranger's own called cells for that lane (in practice, for this dataset
specifically, this intersection turned out not to narrow anything further --
verified 2026-09-15 on the real data: every lane's RNA filtered_feature_bc_matrix.h5
has exactly 30,000 called barcodes, the SAME 30,000-barcode whitelist demuxlet
classified, so "filtered" here means Liu et al.'s own demuxlet-compatible
candidate pool, not an independent cellranger cell-call; the code still does
the intersection rather than trusting that this holds for every lane, and it
does not -- one lane, B3_10xlane15, had 10 demuxlet-SNG barcodes absent from
its own h5, a negligible, fully-traced 0.001% edge case, not a bug).

Real-data run (2026-09-15, all 55 lanes): 791,267 cells x 32,738 genes.
Peak memory for the final sparse.hstack + write step exceeded 32G (OOM-killed
on the batch partition after all 55 lanes' per-lane processing completed
successfully) -- needed the highmem partition at 220G. Harmonized lineage
composition skewed toward B_cell/Monocyte_Myeloid/NK_cell over T_cell and
carried a notably high HSPC share (~10%) -- plausible for a severe/critical
COVID-19 cohort (emergency myelopoiesis, immature low-density granulocytes
co-purifying with the PBMC fraction and scoring MPO/ELANE-positive) but not
independently confirmed; worth a biological sanity check before relying on
these labels downstream, not just accepting them because the code ran clean.

Usage:
    python build_liu_decontx_input.py \
        --singlets-dir <dir of *_singlets.csv from filter_liu_demuxlet_singlets.py> \
        --h5-dir <dir of *_RNA_filtered_feature_bc_matrix.h5> \
        --out-dir <output dir>

Output (matches run_decontx_per_channel.R's expected IO_DIR contract):
    <out_dir>/counts.mtx   -- genes (rows) x cells (cols), MatrixMarket
    <out_dir>/genes.csv    -- one column: gene
    <out_dir>/cells.csv    -- cell_id, channel_id, donor_id, harmonized_lineage
    <out_dir>/build_summary.csv -- per-lane cell counts at each filter stage,
        for a sanity check against liu_gse161918_channel_inventory.csv
"""
import argparse
import glob
import os
import re

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.io import mmwrite

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from harmonized_lineage_labels import label_channel


def lane_from_h5_filename(fname):
    # "GSM4929081_B1_10xlane1_RNA_filtered_feature_bc_matrix.h5" -> "B1_10xlane1"
    m = re.match(r'GSM\d+_(.+)_RNA_filtered_feature_bc_matrix\.h5$', os.path.basename(fname))
    assert m, f'unexpected h5 filename: {fname}'
    return m.group(1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--singlets-dir', required=True)
    ap.add_argument('--h5-dir', required=True)
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    h5_files = sorted(glob.glob(os.path.join(args.h5_dir, '*_RNA_filtered_feature_bc_matrix.h5')))
    if not h5_files:
        raise SystemExit(f'no RNA h5 files found under {args.h5_dir}')

    ref_genes = None
    X_blocks = []          # list of (genes x cells) sparse blocks, per lane
    cells_rows = []         # list of per-cell dicts across all lanes
    summary_rows = []

    for h5_path in h5_files:
        lane = lane_from_h5_filename(h5_path)
        singlets_path = os.path.join(args.singlets_dir, f'{lane}_singlets.csv')
        if not os.path.exists(singlets_path):
            print(f'{lane}: no singlets CSV found, skipping', flush=True)
            continue
        singlets = pd.read_csv(singlets_path)

        adata = sc.read_10x_h5(h5_path)
        adata.var_names_make_unique()
        n_called_cells = adata.n_obs

        if ref_genes is None:
            ref_genes = adata.var_names.astype(str).tolist()
        else:
            assert list(adata.var_names.astype(str)) == ref_genes, \
                f'{lane}: gene list differs from the reference lane -- cellranger reference mismatch'

        # intersect: demuxlet singlet barcodes that are ALSO cellranger-called cells
        called_barcodes = set(adata.obs_names.astype(str))
        singlets_in_called = singlets[singlets['barcode'].isin(called_barcodes)].copy()
        n_intersect = len(singlets_in_called)

        if n_intersect == 0:
            print(f'{lane}: 0 cells survive intersection (called={n_called_cells}, '
                  f'demuxlet_sng={len(singlets)}) -- skipping', flush=True)
            summary_rows.append({'lane': lane, 'n_called_cells': n_called_cells,
                                  'n_demuxlet_sng': len(singlets), 'n_intersect': 0})
            continue

        sub = adata[singlets_in_called['barcode'].to_numpy()].copy()

        labels = label_channel(sub)  # per-lane only, no cross-lane fitting

        block = sub.X.T.tocsr() if sparse.issparse(sub.X) else sparse.csr_matrix(sub.X.T)
        X_blocks.append(block)

        for barcode, cell_id, channel_id, donor_id, lineage in zip(
                singlets_in_called['barcode'], singlets_in_called['cell_id'],
                singlets_in_called['channel_id'], singlets_in_called['donor_id'],
                labels.to_numpy()):
            cells_rows.append({'cell_id': cell_id, 'channel_id': channel_id,
                                'donor_id': donor_id, 'harmonized_lineage': lineage})

        summary_rows.append({'lane': lane, 'n_called_cells': n_called_cells,
                              'n_demuxlet_sng': len(singlets), 'n_intersect': n_intersect})
        print(f'{lane}: cellranger-called={n_called_cells}, demuxlet SNG={len(singlets)}, '
              f'intersection kept={n_intersect}', flush=True)

    if not X_blocks:
        raise SystemExit('no lanes produced any cells -- nothing to write')

    combined = sparse.hstack(X_blocks, format='csr')
    cells_df = pd.DataFrame(cells_rows)
    assert combined.shape == (len(ref_genes), len(cells_df)), \
        f'shape mismatch: matrix {combined.shape} vs genes={len(ref_genes)} cells={len(cells_df)}'

    mmwrite(os.path.join(args.out_dir, 'counts.mtx'), combined)
    pd.DataFrame({'gene': ref_genes}).to_csv(os.path.join(args.out_dir, 'genes.csv'), index=False)
    cells_df.to_csv(os.path.join(args.out_dir, 'cells.csv'), index=False)
    pd.DataFrame(summary_rows).to_csv(os.path.join(args.out_dir, 'build_summary.csv'), index=False)

    print(f'\nDone. {combined.shape[1]} cells x {combined.shape[0]} genes across '
          f'{len(summary_rows)} lanes.')
    print('Harmonized lineage composition:')
    print(cells_df['harmonized_lineage'].value_counts())


if __name__ == '__main__':
    main()
