"""
Assemble Ren/GSE158055's counts.mtx / genes.csv / cells.csv for
run_decontx_per_channel.R.

Ren is structurally simpler than Liu: no multiplexing, one GEM well = one GSM
= one sampleID, already present as a single pooled matrix + per-cell
annotation (GSE158055_cell_annotation.csv.gz: cellName, sampleID, celltype,
majorType). Two things this script deliberately does NOT do, both load-bearing:
  1. Does not use Ren's own `celltype`/`majorType` as the z-label passed to
     DecontX -- per project_302006_rebuild_architecture, z must come from one
     portable, fixed marker panel applied identically to every channel/cohort
     (harmonized_lineage_labels.py), never a cohort's own published
     annotation. Ren's own labels are not used anywhere in this script.
  2. Does not filter cells beyond what's already in cell_annotation.csv --
     every barcode there already has a celltype/majorType, i.e. these are
     already Ren's own called real cells, not raw/empty droplets (unlike
     Liu, there is no separate raw-vs-filtered distinction to reconcile
     here).

Harmonized lineage labeling is still done per-channel (per sampleID, looping
over each of the 284 unique sampleIDs and calling label_channel()
independently), matching Liu's build script and harmonized_lineage_labels.py's
documented contract, even though all 284 channels' cells arrive in one
already-loaded matrix rather than 284 separate files.

Usage:
    python build_ren_decontx_input.py \
        --mtx-dir <dir with GSE158055_covid19_counts.mtx.gz, _barcodes.tsv.gz, _features.tsv.gz> \
        --annotation <GSE158055_cell_annotation.csv.gz path> \
        --out-dir <output dir>

Output (matches run_decontx_per_channel.R's expected IO_DIR contract):
    <out_dir>/counts.mtx   -- genes (rows) x cells (cols), MatrixMarket
    <out_dir>/genes.csv    -- one column: gene
    <out_dir>/cells.csv    -- cell_id, channel_id, harmonized_lineage
    <out_dir>/build_summary.csv -- per-channel (sampleID) cell counts and
        harmonized-lineage composition, for a sanity check
"""
import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy.io import mmread
from scipy.io import mmwrite
from scipy import sparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from harmonized_lineage_labels import label_channel


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mtx-dir', required=True)
    ap.add_argument('--annotation', required=True)
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    print('Loading annotation...', flush=True)
    ann = pd.read_csv(args.annotation)
    print(f'  {len(ann)} annotated cells, {ann["sampleID"].nunique()} unique sampleIDs', flush=True)

    print('Loading barcodes/features...', flush=True)
    with gzip.open(os.path.join(args.mtx_dir, 'GSE158055_covid19_barcodes.tsv.gz'), 'rt') as f:
        barcodes = [l.strip() for l in f]
    with gzip.open(os.path.join(args.mtx_dir, 'GSE158055_covid19_features.tsv.gz'), 'rt') as f:
        genes = [l.strip().split('\t')[0] for l in f]
    print(f'  {len(barcodes)} barcodes, {len(genes)} genes', flush=True)

    assert len(barcodes) == len(ann), \
        f'barcode count ({len(barcodes)}) != annotation row count ({len(ann)}) -- order must match exactly'
    assert list(barcodes) == list(ann['cellName']), \
        'barcodes.tsv order does not match cell_annotation.csv cellName order -- cannot assume alignment'

    print('Loading count matrix (large, sparse, gzipped -- this will take a while)...', flush=True)
    with gzip.open(os.path.join(args.mtx_dir, 'GSE158055_covid19_counts.mtx.gz'), 'rb') as f:
        mat = mmread(f).tocsc()  # genes x cells, as shipped
    print(f'  matrix shape: {mat.shape}', flush=True)
    assert mat.shape == (len(genes), len(barcodes)), \
        f'matrix shape {mat.shape} does not match genes={len(genes)} x barcodes={len(barcodes)}'

    sample_ids = ann['sampleID'].to_numpy()
    unique_samples = sorted(set(sample_ids))
    print(f'Labeling {len(unique_samples)} channels independently...', flush=True)

    cells_rows = []
    summary_rows = []
    for i, sid in enumerate(unique_samples):
        idx = np.where(sample_ids == sid)[0]
        sub_mat = mat[:, idx].T.tocsr()  # cells x genes for this channel, for AnnData/scanpy
        sub_adata = AnnData(X=sub_mat, var=pd.DataFrame(index=genes))
        sub_adata.var_names_make_unique()
        labels = label_channel(sub_adata)

        for j, cell_idx in enumerate(idx):
            cells_rows.append({'cell_id': barcodes[cell_idx], 'channel_id': sid,
                                'harmonized_lineage': labels.iloc[j]})
        lineage_counts = labels.value_counts().to_dict()
        summary_rows.append({'channel_id': sid, 'n_cells': len(idx), **lineage_counts})
        if (i + 1) % 20 == 0 or i == len(unique_samples) - 1:
            print(f'  {i+1}/{len(unique_samples)} channels labeled', flush=True)

    cells_df = pd.DataFrame(cells_rows)
    assert len(cells_df) == mat.shape[1], f'{len(cells_df)} labeled cells != {mat.shape[1]} matrix columns'

    print('Writing outputs...', flush=True)
    mmwrite(os.path.join(args.out_dir, 'counts.mtx'), mat)
    pd.DataFrame({'gene': genes}).to_csv(os.path.join(args.out_dir, 'genes.csv'), index=False)
    cells_df.to_csv(os.path.join(args.out_dir, 'cells.csv'), index=False)
    pd.DataFrame(summary_rows).to_csv(os.path.join(args.out_dir, 'build_summary.csv'), index=False)

    print(f'\nDone. {mat.shape[1]} cells x {mat.shape[0]} genes across {len(unique_samples)} channels.')
    print('Harmonized lineage composition:')
    print(cells_df['harmonized_lineage'].value_counts())


if __name__ == '__main__':
    main()
