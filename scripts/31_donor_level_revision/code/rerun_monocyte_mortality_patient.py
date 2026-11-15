#!/usr/bin/env python3
"""
Job F Phase 4: patient-level repeated-CV mortality evaluation on the monocyte
in-distribution canonical embeddings (STATE/UCE/scPRINT-2/scGPT/Geneformer).

Mirrors scripts/31_donor_level_revision/code/rerun_covid_sepsis_mortality_patient_v2.py
(Job H2 step 2's platelet-side script) exactly in classifier config and CV
protocol, so the two are directly comparable. Differences:
  - reads from 02_EMBEDDINGS/option_c_monocyte/ instead of option_c/
  - no Disease/severity filter or SEVERITY_FIX needed: build_monocyte_canonical.py
    already built exactly the target population (disease, y, patient_id already
    in obs) -- this script just reads obsm[<fm key>] + obs directly.
  - patient_id here is Job F's own (Ren/Qiu-sepsis: reused v1 donor_meta_with_
    patient_id.csv, validated unaffected by the v2 Stephenson/Combes/Wilk/Lee
    fixes; Stephenson: read directly off the source h5ad's own patient_id
    field, independently correct -- see Job F checkpoint-2 report). Not the
    platelet-side v2 mapping file (different cohort, different sample_ids).

CV Protocol v2: --n-seeds N (default 10) runs repeated StratifiedGroupKFold
across N seeds (42..42+N-1) via repeated_cv.py, one oof_monocyte_<FM>_seed<NN>.csv
per seed (or oof_monocyte_<FM>.csv when N=1, backward compatible).

Output dir: 05_RESULTS/revision_donor_level_v2/monocyte/ (new subdir; does not
touch the platelet-side files in the parent dir).
"""
import argparse
import os

import h5py
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression

from repeated_cv import DEFAULT_SEEDS, run_repeated_kfold

ROOT = '/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c_monocyte'
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2/monocyte'

EMBEDDINGS = [
    ('STATE',      'state_monocyte_canonical.h5ad',              'X_state'),
    ('UCE',        'uce_monocyte_canonical.h5ad',                'X_uce'),
    ('Geneformer', 'geneformer_monocyte_canonical.h5ad',         'X_geneformer'),
    ('scPRINT-2',  'scprint2_monocyte_canonical.h5ad_embedded.h5ad', 'scprint_emb'),
    ('scGPT',      'scgpt_monocyte_canonical.h5ad',              'X_scgpt'),
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
        obsm_keys = list(h.get('obsm', []))
        if obsm_key not in obsm_keys:
            print(f'  ERROR: obsm[{obsm_key}] missing (have {obsm_keys}) -- skipping', flush=True)
            return None
        X = h[f'obsm/{obsm_key}'][:]
        cell_id = h['obs/_index'][:].astype(str)
        sample_id = read_categorical(h['obs/Sample_ID'])
        patient_id = read_categorical(h['obs/patient_id'])
        y = h['obs/y'][:].astype(int)
        source_study = read_categorical(h['obs/source_study'])
        provenance = read_categorical(h['obs/monocyte_source_provenance'])

    n_donors = len(np.unique(sample_id))
    n_patients = len(np.unique(patient_id))
    n_fatal = int(y.sum())
    print(f'  {X.shape[0]} cells | {n_donors} sample_ids | {n_patients} patients | '
          f'fatal_cells={n_fatal}', flush=True)
    return dict(X=X, y=y, groups=patient_id, cell_id=cell_id, sample_id=sample_id,
                patient_id=patient_id, source_study=source_study, provenance=provenance)


def fit_predict(X_tr, y_tr, X_te):
    lr = LogisticRegression(C=1.0, solver='lbfgs', class_weight='balanced', max_iter=1000)
    lr.fit(X_tr, y_tr)
    rf = RandomForestClassifier(n_estimators=100, class_weight='balanced', n_jobs=-1,
                                 random_state=42)
    rf.fit(X_tr, y_tr)
    return {'p_lr': lr.predict_proba(X_te)[:, 1], 'p_rf': rf.predict_proba(X_te)[:, 1]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n-seeds', type=int, default=10,
                     help='CV Protocol v2: number of repeated 5-fold splits '
                          '(seeds 42..42+N-1). Default 10 per CV Protocol v2.')
    ap.add_argument('--fm', default=None, help='run only this FM (default: all)')
    args = ap.parse_args()
    seeds = DEFAULT_SEEDS[:args.n_seeds] if args.n_seeds > 1 else [42]

    os.makedirs(OUT_DIR, exist_ok=True)

    for fm, fname, obsm_key in EMBEDDINGS:
        if args.fm and fm != args.fm:
            continue
        data = load_fm(fm, fname, obsm_key)
        if data is None:
            continue

        oof_by_seed = run_repeated_kfold(data['X'], data['y'], data['groups'],
                                          fit_predict_fn=fit_predict, seeds=seeds)

        base = dict(cell_id=data['cell_id'], sample_id=data['sample_id'],
                    patient_id=data['patient_id'], y=data['y'],
                    source_study=data['source_study'], provenance=data['provenance'])

        if args.n_seeds <= 1:
            out = pd.DataFrame({**base, **oof_by_seed[seeds[0]]})
            out_path = os.path.join(OUT_DIR, f'oof_monocyte_{fm}.csv')
            out.to_csv(out_path, index=False)
            print(f'  SAVED -> {out_path}', flush=True)
        else:
            for seed in seeds:
                out = pd.DataFrame({**base, **oof_by_seed[seed]})
                out_path = os.path.join(OUT_DIR, f'oof_monocyte_{fm}_seed{seed}.csv')
                out.to_csv(out_path, index=False)
            print(f'  SAVED {len(seeds)} seed files -> '
                  f'{OUT_DIR}/oof_monocyte_{fm}_seed<42..{seeds[-1]}>.csv', flush=True)


if __name__ == '__main__':
    main()
