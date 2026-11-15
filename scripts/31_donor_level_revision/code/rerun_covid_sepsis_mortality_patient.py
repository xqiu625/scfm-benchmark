#!/usr/bin/env python3
"""
Job H step 2: patient-level re-run of the COVID+Sepsis mortality CV.

Copy of scripts/27_covid_sepsis_only_rerun/code/rerun_covid_sepsis_mortality.py's
exact filters/classifiers/seed, restructured to (a) group StratifiedGroupKFold by
patient_id instead of Sample_ID (2026-09-02 multi-draw leakage fix; see
reports/manuscript/revision/patient_overlap_finding.md) and (b) emit per-cell
out-of-fold predictions (the oof_<FM>.csv schema donor_level_reanalysis.py and
null_model_cellcount.py consume: cell_id, sample_id, severity, y, p_lr, p_rf),
restricted to the 5 FMs used downstream (STATE, UCE, Geneformer, scPRINT-2, scGPT).

Filters/classifiers/seeds identical to the original script:
  Disease in {COVID-19, Sepsis}; covid_severity in
  {fatal, mild, moderate, recovered, severe}
  LR: C=1.0, lbfgs, class_weight=balanced, max_iter=1000
  RF: n_estimators=100, class_weight=balanced, n_jobs=-1, random_state=42
  StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42), groups=patient_id

Output: 05_RESULTS/revision_donor_level/oof_patient_<FM>.csv
"""
import h5py
import numpy as np
import pandas as pd
import os

from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedGroupKFold

ROOT = '/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c'
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level'
PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id.csv')

EMBEDDINGS = [
    ('STATE',      'state_covid_canonical.h5ad',      'X_state'),
    ('UCE',        'uce_covid_canonical.h5ad',        'X_uce'),
    ('Geneformer', 'geneformer_covid_canonical.h5ad', 'X_geneformer'),
    ('scPRINT-2',  'scprint2_covid_canonical.h5ad',   'X_scprint2'),
    ('scGPT',      'scgpt_covid_canonical.h5ad',      'X_scgpt'),
]

KEEP_DISEASE = ['COVID-19', 'Sepsis']
KEEP_SEVERITY = ['fatal', 'mild', 'moderate', 'recovered', 'severe']


def read_categorical(group):
    if isinstance(group, h5py.Group) and 'categories' in group:
        cats = group['categories'][:].astype(str)
        codes = group['codes'][:]
        return cats[codes]
    return group[:].astype(str)


def evaluate_fm(fm, fname, obsm_key, pmap):
    path = os.path.join(ROOT, fname)
    if not os.path.exists(path):
        print(f'SKIP: {path} (not found)', flush=True)
        return None

    print(f'\n=== {fm}: {fname} ===', flush=True)
    with h5py.File(path, 'r') as h:
        cell_id = h['obs/_index'][:].astype(str)
        obsm_keys = list(h.get('obsm', []))
        if obsm_key not in obsm_keys:
            print(f'  ERROR: obsm[{obsm_key}] missing -- skipping', flush=True)
            return None
        X = h[f'obsm/{obsm_key}'][:]
        disease = read_categorical(h['obs/Disease'])
        severity = read_categorical(h['obs/covid_severity'])
        sample_id = read_categorical(h['obs/Sample_ID'])

    keep = np.isin(disease, KEEP_DISEASE) & np.isin(severity, KEEP_SEVERITY)
    Xf = X[keep]
    cid_f = cell_id[keep]
    sev_f = severity[keep]
    sid_f = sample_id[keep]
    y = (sev_f == 'fatal').astype(int)

    df_map = pd.DataFrame({'sample_id': sid_f}).merge(pmap, on='sample_id', how='left')
    n_missing = int(df_map['patient_id'].isna().sum())
    assert n_missing == 0, f'{fm}: {n_missing} cells lack patient_id mapping'
    pid_f = df_map['patient_id'].to_numpy()

    n_donors = len(np.unique(sid_f))
    n_patients = len(np.unique(pid_f))
    n_fatal = int(y.sum())
    print(
        f'  filter: {Xf.shape[0]} cells | {n_donors} sample_ids | '
        f'{n_patients} patients | fatal_cells={n_fatal}',
        flush=True,
    )

    skf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)
    p_lr = np.full(Xf.shape[0], np.nan)
    p_rf = np.full(Xf.shape[0], np.nan)

    for fold, (tr, te) in enumerate(skf.split(Xf, y, groups=pid_f)):
        lr = LogisticRegression(C=1.0, solver='lbfgs', class_weight='balanced',
                                 max_iter=1000)
        lr.fit(Xf[tr], y[tr])
        p_lr[te] = lr.predict_proba(Xf[te])[:, 1]

        rf = RandomForestClassifier(n_estimators=100, class_weight='balanced',
                                     n_jobs=-1, random_state=42)
        rf.fit(Xf[tr], y[tr])
        p_rf[te] = rf.predict_proba(Xf[te])[:, 1]
        print(f'  fold {fold}: train={len(tr)} test={len(te)}', flush=True)

    assert not np.isnan(p_lr).any() and not np.isnan(p_rf).any(), \
        f'{fm}: some cells never held out (groups={n_patients} vs n_splits=5)'

    out = pd.DataFrame({
        'cell_id': cid_f, 'sample_id': sid_f, 'patient_id': pid_f,
        'severity': sev_f, 'y': y, 'p_lr': p_lr, 'p_rf': p_rf,
    })
    return out


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    pmap = pd.read_csv(PATIENT_MAP)[['sample_id', 'patient_id']]

    for fm, fname, obsm_key in EMBEDDINGS:
        out = evaluate_fm(fm, fname, obsm_key, pmap)
        if out is None:
            continue
        out_path = os.path.join(OUT_DIR, f'oof_patient_{fm}.csv')
        out.to_csv(out_path, index=False)
        print(f'  SAVED -> {out_path}', flush=True)


if __name__ == '__main__':
    main()
