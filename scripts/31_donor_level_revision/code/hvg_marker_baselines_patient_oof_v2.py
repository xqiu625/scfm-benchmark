#!/usr/bin/env python3
"""
BIB revision memo item B3 (2026-09-13 memo, Tier 3): HVG-mean and
platelet-activation marker-score baselines -- "a one-number predictor that
competes with a 650M-parameter model is a strong, memorable result."

Same raw-count input, disease/severity filter, patient map, and CV harness
as pca_baseline_patient_oof_v2.py (this project's from-scratch baseline
template) -- only the per-fold feature construction changes:
  HVG-mean:    top-2000 highly-variable genes selected from TRAINING cells
               only (unsupervised, no label leakage), per-cell mean
               log1p-normalized expression across those genes = 1 feature.
  MarkerScore: fixed literature panel (PLATELET_ACTIVATION_GENES below, from
               this project's own scripts/24_option_c_canonical_input/code/
               score_activation_markers.py "core_activation" set: PF4, PPBP,
               SELP, ITGA2B, ITGB3, GP1BA, GP1BB, GP9) -- panel is fixed a
               priori, not data-derived, so no per-fold refitting needed for
               the gene LIST itself; per-cell mean log1p expression across
               the panel = 1 feature.
Both features are then fit with a single-feature logistic regression per
fold (same LR settings as every other baseline in this project), so the
comparison is fair: same CV protocol, same classifier family, only the
representation differs.

Output dir: 05_RESULTS/revision_donor_level_v2/ (same as every other
baseline, so donor_level_reanalysis.py --extra-methods can aggregate these
alongside STATE/UCE/Geneformer/scGPT/scPRINT-2/PCA-200 with zero code
changes).
"""
import argparse
import os

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.linear_model import LogisticRegression

from repeated_cv import DEFAULT_SEEDS, run_repeated_kfold

RAW = ('/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/'
       'human_platelet_covid_severity_raw_counts.h5ad')
PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2'

KEEP_DISEASE = ['COVID-19', 'Sepsis']
KEEP_SEVERITY = ['fatal', 'mild', 'moderate', 'recovered', 'severe']
N_HVGS = 2000

PLATELET_ACTIVATION_GENES = ['PF4', 'PPBP', 'SELP', 'ITGA2B', 'ITGB3', 'GP1BA', 'GP1BB', 'GP9']

SEVERITY_FIX = {
    'GSM4995436': 'moderate',
    'GSM4995442': 'severe',
    'GSM4995432': 'moderate',
}


def make_hvg_mean_fit_predict(var_names):
    def fit_predict(X_tr, y_tr, X_te):
        import anndata as ad
        a_tr = ad.AnnData(X=X_tr.copy())
        a_tr.var_names = var_names
        sc.pp.normalize_total(a_tr, target_sum=1e4)
        sc.pp.log1p(a_tr)
        sc.pp.highly_variable_genes(a_tr, n_top_genes=N_HVGS, flavor='seurat')
        hvg_mask = a_tr.var['highly_variable'].to_numpy()

        def hvg_score(Xraw):
            a = ad.AnnData(X=Xraw.copy())
            sc.pp.normalize_total(a, target_sum=1e4)
            sc.pp.log1p(a)
            Xn = a.X.toarray() if hasattr(a.X, 'toarray') else a.X
            return Xn[:, hvg_mask].mean(axis=1)

        feat_tr = hvg_score(X_tr).reshape(-1, 1)
        feat_te = hvg_score(X_te).reshape(-1, 1)
        clf = LogisticRegression(C=1.0, solver='lbfgs', max_iter=2000, class_weight='balanced')
        clf.fit(feat_tr, y_tr)
        p = clf.predict_proba(feat_te)[:, 1]
        return {'p_lr': p, 'p_rf': p}
    return fit_predict


def make_marker_score_fit_predict(var_names, marker_genes):
    marker_idx = [i for i, g in enumerate(var_names) if g in marker_genes]
    found = [var_names[i] for i in marker_idx]
    print(f'  Marker panel: {len(found)}/{len(marker_genes)} genes found in var_names: {found}', flush=True)
    assert len(marker_idx) >= 3, 'too few marker genes found -- check var_names gene symbol convention'

    def fit_predict(X_tr, y_tr, X_te):
        import anndata as ad

        def marker_score(Xraw):
            a = ad.AnnData(X=Xraw.copy())
            sc.pp.normalize_total(a, target_sum=1e4)
            sc.pp.log1p(a)
            Xn = a.X.toarray() if hasattr(a.X, 'toarray') else a.X
            return Xn[:, marker_idx].mean(axis=1)

        feat_tr = marker_score(X_tr).reshape(-1, 1)
        feat_te = marker_score(X_te).reshape(-1, 1)
        clf = LogisticRegression(C=1.0, solver='lbfgs', max_iter=2000, class_weight='balanced')
        clf.fit(feat_tr, y_tr)
        p = clf.predict_proba(feat_te)[:, 1]
        return {'p_lr': p, 'p_rf': p}
    return fit_predict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n-seeds', type=int, default=10)
    ap.add_argument('--which', choices=['hvg', 'marker', 'both'], default='both')
    args = ap.parse_args()
    seeds = DEFAULT_SEEDS[:args.n_seeds] if args.n_seeds > 1 else [42]

    os.makedirs(OUT_DIR, exist_ok=True)
    print('Loading raw counts...', flush=True)
    a = sc.read_h5ad(RAW)
    print(f'  {a.n_obs} cells x {a.n_vars} genes', flush=True)

    disease = a.obs['Disease'].astype(str).values
    severity = a.obs['covid_severity'].astype(str).values
    sample_id = a.obs['Sample_ID'].astype(str).values
    cell_id = a.obs_names.values.astype(str)
    for sid, new_sev in SEVERITY_FIX.items():
        severity[sample_id == sid] = new_sev

    keep = np.isin(disease, KEEP_DISEASE) & np.isin(severity, KEEP_SEVERITY)
    print(f'  Disease+severity filter: {keep.sum()} / {a.n_obs} cells', flush=True)

    pmap = pd.read_csv(PATIENT_MAP)[['sample_id', 'patient_id']]
    df = pd.DataFrame({'cell_id': cell_id[keep], 'sample_id': sample_id[keep], 'severity': severity[keep]})
    df['_row'] = np.arange(len(df))
    n_before = len(df)
    df = df.merge(pmap, on='sample_id', how='inner')
    print(f'  v2 map filter: {len(df)} / {n_before} cells', flush=True)

    X = a.X[keep][df['_row'].to_numpy()]
    var_names = a.var_names.astype(str).to_numpy()
    y = (df['severity'].to_numpy() == 'fatal').astype(int)
    groups = df['patient_id'].to_numpy()
    print(f'  {X.shape[0]} cells, {len(np.unique(groups))} patients, {int(y.sum())} fatal', flush=True)

    base = dict(cell_id=df['cell_id'].to_numpy(), sample_id=df['sample_id'].to_numpy(),
                patient_id=df['patient_id'].to_numpy(), severity=df['severity'].to_numpy(), y=y)

    jobs = []
    if args.which in ('hvg', 'both'):
        jobs.append(('HVG-mean', make_hvg_mean_fit_predict(var_names)))
    if args.which in ('marker', 'both'):
        jobs.append(('MarkerScore', make_marker_score_fit_predict(var_names, PLATELET_ACTIVATION_GENES)))

    for name, fit_predict_fn in jobs:
        print(f'\n=== {name} ===', flush=True)
        oof_by_seed = run_repeated_kfold(X, y, groups, fit_predict_fn=fit_predict_fn, seeds=seeds)
        for seed in seeds:
            out = pd.DataFrame({**base, **oof_by_seed[seed]})
            out_path = os.path.join(OUT_DIR, f'oof_patient_{name}_seed{seed}.csv')
            out.to_csv(out_path, index=False)
        print(f'SAVED {len(seeds)} seed files -> {OUT_DIR}/oof_patient_{name}_seed<42..{seeds[-1]}>.csv', flush=True)


if __name__ == '__main__':
    main()
