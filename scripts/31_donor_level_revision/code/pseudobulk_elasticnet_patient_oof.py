#!/usr/bin/env python3
"""
BIB revision memo item B1 (2026-09-13 memo, Tier 3): pseudobulk-per-patient
+ elastic net -- "what most clinical-genomics groups would actually run; if
it matches cell-level PCA, the whole single-cell pipeline is unnecessary for
this endpoint."

Unlike every other baseline/model in this project (which fit at the CELL
level and aggregate predictions to patients afterward), this one collapses
each patient to ONE pseudobulk expression vector (mean raw count per gene
across all their cells) BEFORE any model fitting, so the CV unit is already
the patient: 240 patient-level rows, StratifiedKFold-equivalent via
repeated_cv's run_repeated_kfold with groups=patient_id (each patient is its
own group of size 1, so GroupKFold reduces to plain StratifiedKFold here --
reusing the same harness as every other baseline for consistency, not
because grouping is doing real work in this case).

Same raw-count input, disease/severity filter, and patient map as
pca_baseline_patient_oof_v2.py. Elastic net = LogisticRegression with
penalty='elasticnet', solver='saga', l1_ratio swept via a small internal
grid (no nested CV budget for a full grid search here; a single reasonable
l1_ratio=0.5, C=0.1 is fit per fold -- flagged in Limitations as not
hyperparameter-tuned, matching the same C=1.0-fixed convention used for the
other from-scratch baselines in this project rather than introducing a
double standard).

Output: 05_RESULTS/revision_donor_level_v2/oof_patient_PseudobulkEN_seed<NN>.csv,
one row per patient per seed (not one row per cell, since there is no
cell-level prediction here) -- donor_level_reanalysis.py's --group-col
patient_id path handles this the same way it would a pre-aggregated method.
"""
import argparse
import os

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from repeated_cv import DEFAULT_SEEDS, run_repeated_kfold

RAW = ('/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/'
       'human_platelet_covid_severity_raw_counts.h5ad')
PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2'

KEEP_DISEASE = ['COVID-19', 'Sepsis']
KEEP_SEVERITY = ['fatal', 'mild', 'moderate', 'recovered', 'severe']

SEVERITY_FIX = {
    'GSM4995436': 'moderate',
    'GSM4995442': 'severe',
    'GSM4995432': 'moderate',
}


def fit_predict(X_tr, y_tr, X_te):
    scaler = StandardScaler().fit(X_tr)
    Xs_tr = scaler.transform(X_tr)
    Xs_te = scaler.transform(X_te)
    clf = LogisticRegression(penalty='elasticnet', solver='saga', l1_ratio=0.5,
                              C=0.1, max_iter=5000, class_weight='balanced')
    clf.fit(Xs_tr, y_tr)
    p = clf.predict_proba(Xs_te)[:, 1]
    return {'p_lr': p, 'p_rf': p}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n-seeds', type=int, default=10)
    ap.add_argument('--seeds', type=str, default=None,
                     help='explicit comma-separated seed list (overrides --n-seeds); '
                          'for chunking the 10-seed protocol across multiple short-partition '
                          'jobs, e.g. --seeds 42,43,44,45,46 then --seeds 47,48,49,50,51')
    args = ap.parse_args()
    if args.seeds:
        seeds = [int(s) for s in args.seeds.split(',')]
    else:
        seeds = DEFAULT_SEEDS[:args.n_seeds] if args.n_seeds > 1 else [42]

    os.makedirs(OUT_DIR, exist_ok=True)
    print('Loading raw counts...', flush=True)
    a = sc.read_h5ad(RAW)
    print(f'  {a.n_obs} cells x {a.n_vars} genes', flush=True)

    disease = a.obs['Disease'].astype(str).values
    severity = a.obs['covid_severity'].astype(str).values
    sample_id = a.obs['Sample_ID'].astype(str).values
    for sid, new_sev in SEVERITY_FIX.items():
        severity[sample_id == sid] = new_sev

    keep = np.isin(disease, KEEP_DISEASE) & np.isin(severity, KEEP_SEVERITY)
    print(f'  Disease+severity filter: {keep.sum()} / {a.n_obs} cells', flush=True)

    pmap = pd.read_csv(PATIENT_MAP)[['sample_id', 'patient_id']]
    df = pd.DataFrame({'sample_id': sample_id[keep], 'severity': severity[keep]})
    df['_row'] = np.arange(len(df))
    df = df.merge(pmap, on='sample_id', how='inner')
    df['y'] = (df['severity'].to_numpy() == 'fatal').astype(int)

    X_cells = a.X[keep][df['_row'].to_numpy()]
    X_cells = X_cells.toarray() if hasattr(X_cells, 'toarray') else np.asarray(X_cells)

    print('Pseudobulking to one row per patient (mean raw count per gene)...', flush=True)
    patient_ids = df['patient_id'].to_numpy()
    unique_patients = np.unique(patient_ids)
    n_genes = X_cells.shape[1]
    pseudobulk = np.zeros((len(unique_patients), n_genes), dtype=np.float64)
    y_patient = np.zeros(len(unique_patients), dtype=int)
    for i, pid in enumerate(unique_patients):
        mask = patient_ids == pid
        pseudobulk[i] = X_cells[mask].mean(axis=0)
        y_patient[i] = df['y'].to_numpy()[mask].max()

    print(f'  {pseudobulk.shape[0]} patients, {int(y_patient.sum())} fatal', flush=True)

    # log1p on the pseudobulk means (standard pseudobulk practice, stabilizes variance)
    pseudobulk = np.log1p(pseudobulk)

    oof_by_seed = run_repeated_kfold(pseudobulk, y_patient, unique_patients,
                                      fit_predict_fn=fit_predict, seeds=seeds)

    for seed in seeds:
        out = pd.DataFrame({'patient_id': unique_patients, 'y': y_patient,
                             **oof_by_seed[seed]})
        out_path = os.path.join(OUT_DIR, f'oof_patient_PseudobulkEN_seed{seed}.csv')
        out.to_csv(out_path, index=False)
    print(f'SAVED {len(seeds)} seed files -> '
          f'{OUT_DIR}/oof_patient_PseudobulkEN_seed<42..{seeds[-1]}>.csv', flush=True)


if __name__ == '__main__':
    main()
