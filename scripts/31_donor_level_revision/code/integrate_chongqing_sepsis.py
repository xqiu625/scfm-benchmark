#!/usr/bin/env python3
"""
Job S steps 2-3: integrate GSE279452 (Chongqing) adult-sepsis scRNA donors to
the option-C schema, then platelet-annotate with a cluster+marker procedure.

LIMITATION (flagged, not silently assumed): the original integration script
(scripts/data/integrate_platelet_data_harmony.R, referenced by
docs/cohort-sources.md) is not recoverable in this repo (same gap
Job F's Phase 1 inventory found) or on HPCC, so exact per-sample QC thresholds
and clustering parameters cannot be matched precisely. This uses standard,
documented scRNA-seq QC/clustering defaults instead (below) and reports every
threshold used so it can be audited or re-run with different values if the
original parameters ever surface.

Cohort: 91 donors (93 Adult-Sepsis donors in chongqing_sepsis_donor_mortality.csv,
minus 2 with missing mort28), per the 2026-09-03 coordinator decision to use
group=='Adult Sepsis' only (NOT the naive 108-donor cohort=adult filter, which
silently included 8 ICU-control + 7 healthy-control donors -- a real
cohort-contamination bug caught before this integration ran).

QC (per-sample, applied before concatenation):
  min_genes=200, min_cells=3 (gene-level, applied post-concat),
  pct_counts_mt < 20%

Platelet annotation (cluster + marker, not naive per-cell marker counting --
per the design doc's explicit instruction to match the existing integration's
approach): normalize_total -> log1p -> HVG(2000) -> scale -> PCA(50) ->
neighbors -> Leiden(resolution=1.0) -> score_genes on a 12-marker platelet
panel (PPBP, PF4, GP9, ITGA2B, TUBB1, ITGB3, GP1BA, GP6, PLEK, THBS1, SELP,
CLU) -> clusters with mean platelet score > 2 SD above the global mean score
are called platelet clusters.

Output:
  01_DATA/raw/chongqing_sepsis_platelet_raw_counts.h5ad  (platelet subset, raw counts, option-C obs schema)
  results_patient/chongqing_platelet_counts_per_donor.csv
"""
import gzip
import os
import shutil

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp
from anndata import AnnData

MANIFEST = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
            '31_donor_level_revision/results_patient/chongqing_manifest_sepsis93.csv')
EXTRACTED = '/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/chongqing_sepsis_GSE279452/extracted'
OUT_H5AD = ('/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/'
            'chongqing_sepsis_platelet_raw_counts.h5ad')
OUT_COUNTS = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
              '31_donor_level_revision/results_patient/chongqing_platelet_counts_per_donor.csv')

# real tar filename suffix differs from the mortality-table donor_id for these
# (GEO metadata/tar-filename mismatch caught 2026-09-03; see runbook Job S notes)
FILE_SUFFIX_OVERRIDE = {'Skin-S152': 'Skin-152'}

PLATELET_MARKERS = ['PPBP', 'PF4', 'GP9', 'ITGA2B', 'TUBB1', 'ITGB3',
                     'GP1BA', 'GP6', 'PLEK', 'THBS1', 'SELP', 'CLU']

MIN_GENES = 200
MIN_CELLS_GENE = 3
MAX_PCT_MT = 20.0


def load_one(donor_id, gsm):
    suffix = FILE_SUFFIX_OVERRIDE.get(donor_id, donor_id)
    prefix = os.path.join(EXTRACTED, f'{gsm}_{suffix}')
    mtx_path = f'{prefix}.matrix.mtx.gz'
    bc_path = f'{prefix}.barcodes.tsv.gz'
    feat_path = f'{prefix}.features.tsv.gz'
    if not (os.path.exists(mtx_path) and os.path.exists(bc_path) and os.path.exists(feat_path)):
        print(f'  MISSING files for {donor_id} ({gsm}) -- skipping', flush=True)
        return None

    with gzip.open(mtx_path, 'rb') as f:
        X = scipy.io.mmread(f).T.tocsr()  # genes x cells in file -> cells x genes
    with gzip.open(bc_path, 'rt') as f:
        barcodes = [line.strip() for line in f]
    with gzip.open(feat_path, 'rt') as f:
        features = [line.strip().split('\t') for line in f]
    gene_ids = [f[0] for f in features]  # Ensembl IDs -- unique by construction
    gene_names = [f[1] if len(f) > 1 else f[0] for f in features]

    a = AnnData(X=X)
    a.obs_names = [f'{donor_id}_{bc}' for bc in barcodes]
    a.var_names = gene_ids
    a.var['gene_symbol'] = gene_names
    return a


def main():
    manifest = pd.read_csv(MANIFEST)
    print(f'Integrating {len(manifest)} donors', flush=True)

    adatas = []
    for _, row in manifest.iterrows():
        donor_id, gsm, mort28 = row['donor_id'], row['gsm'], row['mort28']
        a = load_one(donor_id, gsm)
        if a is None:
            continue
        n_before = a.n_obs
        sc.pp.filter_cells(a, min_genes=MIN_GENES)
        mt_mask = a.var['gene_symbol'].str.startswith('MT-').to_numpy()
        if mt_mask.sum() > 0:
            pct_mt = np.asarray(a[:, mt_mask].X.sum(axis=1)).ravel() / (
                np.asarray(a.X.sum(axis=1)).ravel() + 1e-9) * 100
            a = a[pct_mt < MAX_PCT_MT].copy()
        a.obs['Sample_ID'] = donor_id
        a.obs['Disease'] = 'Sepsis'
        a.obs['Data_NO'] = 'chongqing2025'
        a.obs['covid_severity'] = 'fatal' if float(mort28) == 1.0 else 'severe'
        a.obs['mort14'] = row.get('mort14', np.nan)
        a.obs['mort28'] = mort28
        a.obs['mort90'] = row.get('mort90', np.nan)
        print(f'  {donor_id}: {n_before} -> {a.n_obs} cells post-QC '
              f'(mort28={mort28}, severity={a.obs["covid_severity"].iloc[0] if a.n_obs else "n/a"})',
              flush=True)
        adatas.append(a)

    print(f'\nConcatenating {len(adatas)} donors...', flush=True)
    full = sc.concat(adatas, join='outer', fill_value=0, merge='same')
    sc.pp.filter_genes(full, min_cells=MIN_CELLS_GENE)
    full.layers['counts'] = full.X.copy()
    print(f'  merged: {full.n_obs} cells x {full.n_vars} genes', flush=True)

    print('\nPlatelet annotation (cluster + marker)...', flush=True)
    work = full.copy()
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    sc.pp.highly_variable_genes(work, n_top_genes=2000)
    work_hvg = work[:, work.var['highly_variable']].copy()
    sc.pp.scale(work_hvg, max_value=10)
    sc.tl.pca(work_hvg, n_comps=50, svd_solver='arpack')
    sc.pp.neighbors(work_hvg, n_pcs=50)
    sc.tl.leiden(work_hvg, resolution=1.0)
    work.obs['leiden'] = work_hvg.obs['leiden']

    symbol_to_ensg = pd.Series(work.var_names.to_numpy(), index=work.var['gene_symbol']).to_dict()
    present_markers = [symbol_to_ensg[g] for g in PLATELET_MARKERS if g in symbol_to_ensg]
    present_marker_symbols = [g for g in PLATELET_MARKERS if g in symbol_to_ensg]
    print(f'  platelet markers present: {present_marker_symbols} '
          f'({len(present_markers)}/{len(PLATELET_MARKERS)})', flush=True)
    sc.tl.score_genes(work, present_markers, score_name='platelet_score')

    cluster_scores = work.obs.groupby('leiden')['platelet_score'].mean()
    global_mean, global_std = work.obs['platelet_score'].mean(), work.obs['platelet_score'].std()
    threshold = global_mean + 2 * global_std
    platelet_clusters = cluster_scores[cluster_scores > threshold].index.tolist()
    print(f'  cluster platelet scores:\n{cluster_scores.sort_values(ascending=False).to_string()}',
          flush=True)
    print(f'  threshold (global_mean + 2*SD): {threshold:.4f}', flush=True)
    print(f'  platelet clusters: {platelet_clusters}', flush=True)

    is_platelet = work.obs['leiden'].isin(platelet_clusters)
    full.obs['leiden'] = work.obs['leiden']
    full.obs['platelet_score'] = work.obs['platelet_score']
    full.obs['is_platelet'] = is_platelet.to_numpy()

    n_total = full.n_obs
    n_platelet = int(is_platelet.sum())
    print(f'\nTOTAL: {n_total} cells, PLATELET: {n_platelet} ({100*n_platelet/n_total:.1f}%)',
          flush=True)

    per_donor = full.obs.groupby('Sample_ID').agg(
        n_cells=('is_platelet', 'size'),
        n_platelet=('is_platelet', 'sum'),
    ).reset_index()
    per_donor['pct_platelet'] = 100 * per_donor['n_platelet'] / per_donor['n_cells']
    per_donor.to_csv(OUT_COUNTS, index=False)
    print(f'SAVED per-donor counts -> {OUT_COUNTS}', flush=True)

    platelet_subset = full[full.obs['is_platelet']].copy()
    platelet_subset.X = platelet_subset.layers['counts']
    del platelet_subset.layers['counts']
    platelet_subset.write(OUT_H5AD)
    print(f'SAVED platelet subset ({platelet_subset.n_obs} cells x '
          f'{platelet_subset.n_vars} genes) -> {OUT_H5AD}', flush=True)


if __name__ == '__main__':
    main()
