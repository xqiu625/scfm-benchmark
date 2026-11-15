#!/usr/bin/env python3
"""
Compute UMAP coordinates per FM on the COVID+Sepsis-only subset (matches the
rerun in scripts/27_covid_sepsis_only_rerun).

Saves one CSV per FM at:
  05_RESULTS/option_c_primary_benchmark/covid_sepsis_only/umap/<FM>.csv

Columns: cell_id, sample_id, disease, severity, y_mortality, umap1, umap2

Memory-safe: h5py-only access to FM obsm; never loads X.
"""
import os
import h5py
import numpy as np
import pandas as pd
import umap

ROOT    = '/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c'
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/option_c_primary_benchmark/covid_sepsis_only/umap'
os.makedirs(OUT_DIR, exist_ok=True)

EMBEDDINGS = [
    ('STATE',            'state_covid_canonical.h5ad',            'X_state'),
    ('UCE',              'uce_covid_canonical.h5ad',              'X_uce'),
    ('TranscriptFormer', 'transcriptformer_covid_canonical.h5ad', 'X_transcriptformer'),
    ('scPRINT-2',        'scprint2_covid_canonical.h5ad',         'X_scprint2'),
    ('STACK',            'stack_covid_canonical.h5ad',            'X_stack'),
    ('Geneformer',       'geneformer_covid_canonical.h5ad',       'X_geneformer'),
]

KEEP_DISEASE  = ['COVID-19', 'Sepsis']
KEEP_SEVERITY = ['fatal', 'mild', 'moderate', 'recovered', 'severe']


def read_categorical(group):
    if isinstance(group, h5py.Group) and 'categories' in group:
        cats  = group['categories'][:].astype(str)
        codes = group['codes'][:]
        return cats[codes]
    return group[:].astype(str)


def read_obs(h, col):
    if col in h['obs']:
        return read_categorical(h[f'obs/{col}'])
    return None


def load_filtered(path, obsm_key):
    with h5py.File(path, 'r') as h:
        X         = h[f'obsm/{obsm_key}'][:]
        disease   = read_obs(h, 'Disease')
        severity  = read_obs(h, 'covid_severity')
        sample_id = read_obs(h, 'Sample_ID')
        cell_id   = read_obs(h, 'Cell_ID')
    keep = np.isin(disease, KEEP_DISEASE) & np.isin(severity, KEEP_SEVERITY)
    return {
        'X':         X[keep],
        'disease':   disease[keep],
        'severity':  severity[keep],
        'sample_id': sample_id[keep],
        'cell_id':   (cell_id[keep] if cell_id is not None else np.arange(int(keep.sum())).astype(str)),
    }


def main():
    for fm, fname, key in EMBEDDINGS:
        path = os.path.join(ROOT, fname)
        out  = os.path.join(OUT_DIR, f'{fm}.csv')
        if os.path.exists(out):
            print(f'SKIP existing: {out}', flush=True)
            continue
        if not os.path.exists(path):
            print(f'SKIP missing input: {path}', flush=True)
            continue
        print(f'\n=== {fm} ===', flush=True)
        D = load_filtered(path, key)
        X = D['X']
        n, d = X.shape
        print(f'  {n} cells, {d} dims', flush=True)

        reducer = umap.UMAP(
            n_components=2,
            n_neighbors=30,
            min_dist=0.3,
            metric='cosine',
            random_state=42,
            verbose=True,
            n_jobs=-1,
        )
        emb2 = reducer.fit_transform(X)

        y = (D['severity'] == 'fatal').astype(int)
        df = pd.DataFrame({
            'cell_id':     D['cell_id'],
            'sample_id':   D['sample_id'],
            'disease':     D['disease'],
            'severity':    D['severity'],
            'y_mortality': y,
            'umap1':       emb2[:, 0],
            'umap2':       emb2[:, 1],
        })
        df.to_csv(out, index=False)
        print(f'  Wrote {out}', flush=True)

    print('\nALL DONE.', flush=True)


if __name__ == '__main__':
    main()
