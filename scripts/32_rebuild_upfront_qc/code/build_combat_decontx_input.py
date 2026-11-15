"""
Assemble COMBAT's per-channel counts.mtx/cells.csv files directly (the
split_matrix_per_channel.py output format that run_decontx_looped_per_channel.R
consumes), from COMBAT-CITESeq-DATA.h5ad.

Unlike Liu, COMBAT needs no demuxlet-style singlet filter: the deposit's own
obs already carries a clean, single, per-cell donor call (COMBAT_ID) -- there
is no separate SNG/DBL/AMB category to resolve, so no evidence the depositors'
donor assignment is unreliable the way Liu's untouched .best files were.
Unlike Ren, COMBAT's per-channel split isn't a hard requirement (verified
2026-09-16: `layers['raw']` has 1,255,272,732 non-zeros before even dropping
the antibody features, comfortably under R's 2,147,483,647 ceiling) -- writing
directly to per-channel files anyway, both for consistency with the now-proven
fast/parallel pipeline and to avoid materializing one unnecessary large
intermediate combined file the way Ren's build originally did before the
ceiling was discovered.

Two things specific to this h5ad's structure (verified directly, not assumed):
  1. `X` is normalized (non-integer values); `layers['raw']` is the actual
     raw-count layer DecontX needs.
  2. `var` mixes Gene Expression (20,615 features, real ENSG ids,
     feature_types==1) and Antibody Capture (192 features, "AB_..." names,
     feature_types==0) in one matrix -- subset to feature_types==1 before
     writing anything out; DecontX/the harmonized-lineage marker panel are
     GEX-only.

z (harmonized_lineage) is computed per-channel via harmonized_lineage_labels.py,
never COMBAT's own Annotation_* columns (per project_302006_rebuild_architecture
- z must come from one portable fixed marker panel, not a cohort's published
annotation). channel_id = Channel_ID (the GEM well, not Pool_ID, which is
coarser and not the physical ambient/doublet unit).

Usage:
    python build_combat_decontx_input.py --h5ad <path> --out-dir <split_per_channel output dir>
"""
import argparse
import os

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy import sparse
from scipy.io import mmwrite

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from harmonized_lineage_labels import label_channel


def read_obs_categorical(f, name):
    """Handles both anndata HDF5 encodings, verified directly against this
    file (2026-09-16): newer anndata stores a Group with sibling
    'codes'/'categories' datasets; THIS file uses the older encoding, where
    obs/<name> is a plain int8 codes dataset (with an attrs['categories']
    HDF5 object reference, not itself string data) and the actual string
    categories live separately at obs/__categories/<name>."""
    grp = f['obs'][name]
    if isinstance(grp, h5py.Group):
        codes = grp['codes'][:]
        cats = grp['categories'][:]
        cats = np.array([c.decode() if isinstance(c, bytes) else c for c in cats])
        return cats[codes]
    if '__categories' in f['obs'] and name in f['obs']['__categories']:
        codes = grp[:]
        cats = f['obs']['__categories'][name][:]
        cats = np.array([c.decode() if isinstance(c, bytes) else c for c in cats])
        return cats[codes]
    arr = grp[:]
    return np.array([c.decode() if isinstance(c, bytes) else c for c in arr])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--h5ad', required=True)
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    f = h5py.File(args.h5ad, 'r')

    print('Reading var (features)...', flush=True)
    var_names = np.array([v.decode() if isinstance(v, bytes) else v for v in f['var']['_index'][:]])
    feature_types = f['var']['feature_types'][:]  # int8: 1 = Gene Expression, 0 = Antibody Capture
    gex_mask = feature_types == 1
    genes = var_names[gex_mask]
    print(f'  {len(var_names)} total features, {gex_mask.sum()} Gene Expression, '
          f'{(~gex_mask).sum()} Antibody Capture -- keeping GEX only', flush=True)

    print('Reading obs (Channel_ID, COMBAT_ID, cell barcodes)...', flush=True)
    channel_id = read_obs_categorical(f, 'Channel_ID')
    donor_id = read_obs_categorical(f, 'COMBAT_ID')
    cell_barcodes = read_obs_categorical(f, '_index') if '_index' in f['obs'] else np.array(
        [f'cell{i}' for i in range(len(channel_id))])
    n_cells = len(channel_id)
    print(f'  {n_cells} cells, {len(set(channel_id))} channels', flush=True)

    print('Reading raw counts layer (sparse CSR, cells x all_features)...', flush=True)
    raw = f['layers']['raw']
    full = sparse.csr_matrix((raw['data'][:], raw['indices'][:], raw['indptr'][:]),
                              shape=(n_cells, len(var_names)))
    print(f'  full matrix: {full.shape}, nnz={full.nnz}', flush=True)
    gex = full[:, gex_mask].tocsc()  # cells x genes, GEX only
    del full
    print(f'  GEX-only matrix: {gex.shape}, nnz={gex.nnz}', flush=True)

    unique_channels = sorted(set(channel_id))
    manifest_rows = []
    for i, ch in enumerate(unique_channels):
        idx = np.where(channel_id == ch)[0]
        sub = gex[idx, :]  # cells x genes for this channel
        sub_adata = AnnData(X=sub.tocsr(), var=pd.DataFrame(index=genes))
        sub_adata.var_names_make_unique()
        labels = label_channel(sub_adata)

        ch_dir = os.path.join(args.out_dir, ch)
        os.makedirs(ch_dir, exist_ok=True)
        mmwrite(os.path.join(ch_dir, 'counts.mtx'), sub.T.tocsc())  # genes x cells, matches the R-side contract
        cells_df = pd.DataFrame({'cell_id': cell_barcodes[idx], 'channel_id': ch,
                                  'donor_id': donor_id[idx],
                                  'harmonized_lineage': labels.to_numpy()})
        cells_df.to_csv(os.path.join(ch_dir, 'cells.csv'), index=False)

        manifest_rows.append({'channel_id': ch, 'n_cells': len(idx), 'nnz': sub.nnz,
                               'n_donors': cells_df['donor_id'].nunique()})
        if (i + 1) % 10 == 0 or i == len(unique_channels) - 1:
            print(f'  {i+1}/{len(unique_channels)} channels written', flush=True)

    pd.DataFrame({'gene': genes}).to_csv(os.path.join(args.out_dir, '..', 'genes.csv'), index=False)
    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(os.path.join(args.out_dir, 'split_manifest.csv'), index=False)
    print(f'\nDone. {len(unique_channels)} channels written to {args.out_dir}')
    print(f'Max per-channel nnz: {manifest.nnz.max():,}')
    print(f'Total cells: {manifest.n_cells.sum()}, donors/channel range: '
          f'{manifest.n_donors.min()}-{manifest.n_donors.max()}')


if __name__ == '__main__':
    main()
