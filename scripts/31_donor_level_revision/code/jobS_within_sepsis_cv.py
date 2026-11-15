#!/usr/bin/env python3
"""
Job S step 5a: within-sepsis mortality CV on the Chongqing cohort (91 adult
sepsis donors / 17 fatal, mort28 endpoint), patient-level, CV Protocol v2
(10 seeds x 5-fold StratifiedGroupKFold).

Chongqing has one sample per patient (Sample_ID IS the patient identifier;
no multi-draw structure like the platelet cohort's Stephenson/Ren/etc.), so
groups=Sample_ID directly -- no patient_id mapping file needed.

Output: 05_RESULTS/revision_donor_level_v2/oof_sepsis_<FM>_seed<NN>.csv
"""
import argparse
import os

import h5py
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression

from repeated_cv import DEFAULT_SEEDS, run_repeated_kfold

ROOT = '/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c_chongqing'
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2'

EMBEDDINGS = [
    ('STATE',      'state_chongqing_canonical.h5ad',      'X_state'),
    ('UCE',        'uce_chongqing_canonical.h5ad',        'X_uce'),
    ('Geneformer', 'geneformer_chongqing_canonical.h5ad', 'X_geneformer'),
    ('scPRINT-2',  'scprint2_chongqing_canonical.h5ad',   'X_scprint2'),
    ('scGPT',      'scgpt_chongqing_canonical.h5ad',      'X_scgpt'),
]


def read_categorical(group):
    if isinstance(group, h5py.Group) and 'categories' in group:
        cats = group['categories'][:].astype(str)
        codes = group['codes'][:]
        return cats[codes]
    return group[:].astype(str)


def load_fm(fm, fname, obsm_key):
    path = os.path.join(ROOT, fname)
    if not os.path.exists(path):
        print(f'SKIP: {path} (not found)', flush=True)
        return None
    print(f'\n=== {fm}: {fname} ===', flush=True)
    with h5py.File(path, 'r') as h:
        cell_id = h['obs/_index'][:].astype(str)
        if obsm_key not in h.get('obsm', []):
            print(f'  ERROR: obsm[{obsm_key}] missing -- skipping', flush=True)
            return None
        X = h[f'obsm/{obsm_key}'][:]
        severity = read_categorical(h['obs/covid_severity'])
        sample_id = read_categorical(h['obs/Sample_ID'])

    y = (severity == 'fatal').astype(int)
    n_patients = len(np.unique(sample_id))
    n_fatal = len(np.unique(sample_id[y == 1]))
    print(f'  {X.shape[0]} cells | {n_patients} patients | {n_fatal} fatal patients',
          flush=True)
    return dict(X=X, y=y, groups=sample_id, cell_id=cell_id, sample_id=sample_id,
                severity=severity)


def fit_predict(X_tr, y_tr, X_te):
    lr = LogisticRegression(C=1.0, solver='lbfgs', class_weight='balanced', max_iter=1000)
    lr.fit(X_tr, y_tr)
    rf = RandomForestClassifier(n_estimators=100, class_weight='balanced', n_jobs=-1,
                                 random_state=42)
    rf.fit(X_tr, y_tr)
    return {'p_lr': lr.predict_proba(X_te)[:, 1], 'p_rf': rf.predict_proba(X_te)[:, 1]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n-seeds', type=int, default=10)
    args = ap.parse_args()
    seeds = DEFAULT_SEEDS[:args.n_seeds] if args.n_seeds > 1 else [42]

    os.makedirs(OUT_DIR, exist_ok=True)
    for fm, fname, obsm_key in EMBEDDINGS:
        data = load_fm(fm, fname, obsm_key)
        if data is None:
            continue
        oof_by_seed = run_repeated_kfold(data['X'], data['y'], data['groups'],
                                          fit_predict_fn=fit_predict, seeds=seeds)
        base = dict(cell_id=data['cell_id'], sample_id=data['sample_id'],
                    patient_id=data['sample_id'], severity=data['severity'], y=data['y'])
        if args.n_seeds <= 1:
            out = pd.DataFrame({**base, **oof_by_seed[seeds[0]]})
            out.to_csv(os.path.join(OUT_DIR, f'oof_sepsis_{fm}.csv'), index=False)
        else:
            for seed in seeds:
                out = pd.DataFrame({**base, **oof_by_seed[seed]})
                out.to_csv(os.path.join(OUT_DIR, f'oof_sepsis_{fm}_seed{seed}.csv'),
                           index=False)
            print(f'  SAVED {len(seeds)} seed files -> oof_sepsis_{fm}_seed<NN>.csv', flush=True)


if __name__ == '__main__':
    main()
