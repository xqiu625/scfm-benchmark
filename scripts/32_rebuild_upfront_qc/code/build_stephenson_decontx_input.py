"""
Assemble Stephenson/E-MTAB-10026's counts.mtx / genes.csv / cells.csv for
run_decontx_per_channel.R (single-combined-matrix + batch= driver -- GEX-only
nnz is 965,744,336, safely under R/Matrix's 2^31-1 ceiling, no per-channel
split needed).

Source: covid_portal_210320_with_raw.h5ad (ArrayExpress E-MTAB-10026,
publicly downloadable, "isPublic": true via the BioStudies API -- NOT
controlled access, unlike Schulte-Schrepping/EGA). `X` is normalized;
`layers['raw']` is the actual raw-count layer (verified directly, same
COMBAT precedent). `var['feature_types']` mixes Gene Expression and Antibody
Capture (CITE-seq) features -- subset to Gene Expression only, same as
COMBAT. Uses the OLDER anndata HDF5 categorical encoding throughout (obs/var
codes as plain int arrays, categories in obs/__categories/<name> --
verified directly, same read_obs_categorical() pattern as
build_combat_decontx_input.py), but obs's own row index is stored under the
non-standard column name "covid_index" (obs.attrs['_index'] ==
"covid_index"), not "_index" itself -- a real structural difference from
every other cohort's h5ad, caught by inspection before writing this script.

channel_id = sample_id (the GEM-well/library unit, 143 channels across 130
patients -- most patients have exactly one sample; a few have more, per
n_samples in T3_patient_sample_mapping/stephenson_patient_outcomes.csv).
donor_id = patient_id, matches canonical_donor_labels.csv's Stephenson rows
directly (e.g. "AP1"). draw_time = Days_from_onset, a per-cell categorical
field whose categories are numeric-day strings (empty string / not
resolvable for some cells) -- cast to float, NaN where blank.

Usage:
    python build_stephenson_decontx_input.py --h5ad <path> --out-dir <dir>

Output (matches run_decontx_per_channel.R's expected IO_DIR contract):
    <out_dir>/counts.mtx, genes.csv,
    cells.csv (cell_id, channel_id, donor_id, harmonized_lineage, draw_time)
"""
import argparse
import os

import h5py
import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse
from scipy.io import mmwrite

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from harmonized_lineage_labels import label_channel


def read_obs_categorical(f, name, group='obs'):
    grp = f[group][name]
    if isinstance(grp, h5py.Group):
        codes = grp['codes'][:]
        cats = grp['categories'][:]
        cats = np.array([c.decode() if isinstance(c, bytes) else c for c in cats])
        return cats[codes]
    if '__categories' in f[group] and name in f[group]['__categories']:
        codes = grp[:]
        cats = f[group]['__categories'][name][:]
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
    feature_types = read_obs_categorical(f, 'feature_types', group='var')
    gex_mask = feature_types == 'Gene Expression'
    genes = var_names[gex_mask]
    print(f'  {len(var_names)} total features, {gex_mask.sum()} Gene Expression', flush=True)

    print('Reading obs (sample_id, patient_id, Days_from_onset, cell barcodes)...', flush=True)
    channel_id = read_obs_categorical(f, 'sample_id')
    donor_id = read_obs_categorical(f, 'patient_id')
    days_raw = read_obs_categorical(f, 'Days_from_onset')
    draw_time = pd.to_numeric(pd.Series(days_raw).replace('', np.nan), errors='coerce').to_numpy()
    cell_barcodes = read_obs_categorical(f, 'covid_index')
    n_cells = len(channel_id)
    print(f'  {n_cells} cells, {len(set(channel_id))} channels, {len(set(donor_id))} donors', flush=True)

    print('Reading raw counts layer (sparse CSR, cells x all_features)...', flush=True)
    raw = f['layers']['raw']
    full = sparse.csr_matrix((raw['data'][:], raw['indices'][:], raw['indptr'][:]),
                              shape=(n_cells, len(var_names)))
    print(f'  full matrix: {full.shape}, nnz={full.nnz}', flush=True)
    gex = full[:, gex_mask].tocsc()  # cells x genes, GEX only
    del full
    print(f'  GEX-only matrix: {gex.shape}, nnz={gex.nnz}', flush=True)

    unique_channels = sorted(set(channel_id))
    X_blocks = []
    cells_rows = []
    for i, ch in enumerate(unique_channels):
        idx = np.where(channel_id == ch)[0]
        sub = gex[idx, :]  # cells x genes for this channel
        sub_adata = AnnData(X=sub.tocsr(), var=pd.DataFrame(index=genes))
        sub_adata.var_names_make_unique()
        labels = label_channel(sub_adata)

        X_blocks.append(sub.T.tocsc())  # genes x cells
        for bc, donor, dt, lineage in zip(cell_barcodes[idx], donor_id[idx], draw_time[idx],
                                           labels.to_numpy()):
            cells_rows.append({'cell_id': bc, 'channel_id': ch, 'donor_id': donor,
                                'harmonized_lineage': lineage, 'draw_time': dt})
        if (i + 1) % 20 == 0 or i == len(unique_channels) - 1:
            print(f'  {i+1}/{len(unique_channels)} channels processed', flush=True)

    combined = sparse.hstack(X_blocks, format='csc')
    cells_df = pd.DataFrame(cells_rows)
    assert combined.shape[1] == len(cells_df)

    mmwrite(os.path.join(args.out_dir, 'counts.mtx'), combined)
    pd.DataFrame({'gene': genes}).to_csv(os.path.join(args.out_dir, 'genes.csv'), index=False)
    cells_df.to_csv(os.path.join(args.out_dir, 'cells.csv'), index=False)

    print(f'\nDone. {combined.shape[1]} cells x {combined.shape[0]} genes across '
          f'{len(unique_channels)} channels, {cells_df.donor_id.nunique()} donors.')
    print('Harmonized lineage composition:')
    print(cells_df['harmonized_lineage'].value_counts())


if __name__ == '__main__':
    main()
