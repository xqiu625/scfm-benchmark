#!/usr/bin/env python3
"""
Job H2 step 2: patient-level re-run of the COVID+Sepsis mortality CV on the
v2 (pool-excluded) mapping.

Same as rerun_covid_sepsis_mortality_patient.py (Job H step 2) except:
  - patient_id comes from donor_meta_with_patient_id_v2.csv (Stephenson fix +
    Wilk C1 A/B merge + Lee nCoV pair merges + Combes pooled-GSM exclusion)
  - cells whose Sample_ID isn't in the v2 map (the 3 dropped Combes pooled
    GSMs: GSM4995425/26/30) are DROPPED from the cohort entirely, not just
    unmapped -- this is a data-validity exclusion, not a missing-mapping bug
  - output goes to 05_RESULTS/revision_donor_level_v2/ (new dir; v1's
    05_RESULTS/revision_donor_level/ is left untouched as historical record)

CV Protocol v2 (2026-09-04): --n-seeds N (default 1) runs repeated
StratifiedGroupKFold across N seeds (42, 43, ..., 42+N-1) via repeated_cv.py,
writing one oof_patient_<FM>_seed<NN>.csv per seed when N>1 (schema unchanged
per file), or the original single oof_patient_<FM>.csv when N=1 (backward
compatible with anything still reading the single-split file).

Filters/classifiers identical to the original stage-27 script:
  Disease in {COVID-19, Sepsis}; covid_severity in
  {fatal, mild, moderate, recovered, severe}
  LR: C=1.0, lbfgs, class_weight=balanced, max_iter=1000
  RF: n_estimators=100, class_weight=balanced, n_jobs=-1, random_state=42
  StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=<seed>), groups=patient_id
"""
import argparse
import os

import h5py
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression

from repeated_cv import DEFAULT_SEEDS, run_repeated_kfold

ROOT = '/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c'
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2'
PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')

EMBEDDINGS = [
    ('STATE',      'state_covid_canonical.h5ad',      'X_state'),
    ('UCE',        'uce_covid_canonical.h5ad',        'X_uce'),
    ('Geneformer', 'geneformer_covid_canonical.h5ad', 'X_geneformer'),
    ('scPRINT-2',  'scprint2_covid_canonical.h5ad',   'X_scprint2'),
    ('scGPT',      'scgpt_covid_canonical.h5ad',      'X_scgpt'),
]

KEEP_DISEASE = ['COVID-19', 'Sepsis']
KEEP_SEVERITY = ['fatal', 'mild', 'moderate', 'recovered', 'severe']

# Job H2 Step 0 (2026-09-03, Task-4 audit vs Combes 2021 Supp Table 1 p.3):
# three Combes covid_severity labels are wrong in the canonical h5ad obs.
# Applied at read-time here (not by mutating the shared h5ad files, which
# other concurrent jobs/agents may be reading) rather than in-place.
SEVERITY_FIX = {
    'GSM4995436': 'moderate',  # was fatal; paper: Discharged=Yes, survived
    'GSM4995442': 'severe',    # was fatal; paper: Discharged=Yes, survived
    'GSM4995432': 'moderate',  # was severe; paper: Mild/Moderate arm
}


def read_categorical(group):
    if isinstance(group, h5py.Group) and 'categories' in group:
        cats = group['categories'][:].astype(str)
        codes = group['codes'][:]
        return cats[codes]
    return group[:].astype(str)


def apply_severity_fix(sample_id, severity):
    severity = severity.copy()
    for sid, new_sev in SEVERITY_FIX.items():
        severity[sample_id == sid] = new_sev
    return severity


def load_fm(fm, fname, obsm_key, pmap):
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
    severity = apply_severity_fix(sample_id, severity)

    keep = np.isin(disease, KEEP_DISEASE) & np.isin(severity, KEEP_SEVERITY)
    Xf = X[keep]
    cid_f = cell_id[keep]
    sev_f = severity[keep]
    sid_f = sample_id[keep]

    df = pd.DataFrame({'cell_id': cid_f, 'sample_id': sid_f, 'severity': sev_f})
    df['_row'] = np.arange(len(df))
    n_before = len(df)
    df = df.merge(pmap, on='sample_id', how='inner')  # drops v2-excluded (Combes pooled) cells
    n_dropped = n_before - len(df)
    if n_dropped:
        print(f'  dropped {n_dropped} cells not in v2 map (Combes pooled GSMs)', flush=True)

    Xf = Xf[df['_row'].to_numpy()]
    cid_f = df['cell_id'].to_numpy()
    sid_f = df['sample_id'].to_numpy()
    sev_f = df['severity'].to_numpy()
    pid_f = df['patient_id'].to_numpy()
    y = (sev_f == 'fatal').astype(int)

    n_donors = len(np.unique(sid_f))
    n_patients = len(np.unique(pid_f))
    n_fatal = int(y.sum())
    print(
        f'  filter: {Xf.shape[0]} cells | {n_donors} sample_ids | '
        f'{n_patients} patients | fatal_cells={n_fatal}',
        flush=True,
    )
    return dict(X=Xf, y=y, groups=pid_f, cell_id=cid_f, sample_id=sid_f,
                patient_id=pid_f, severity=sev_f)


def fit_predict(X_tr, y_tr, X_te):
    lr = LogisticRegression(C=1.0, solver='lbfgs', class_weight='balanced', max_iter=1000)
    lr.fit(X_tr, y_tr)
    rf = RandomForestClassifier(n_estimators=100, class_weight='balanced', n_jobs=-1,
                                 random_state=42)
    rf.fit(X_tr, y_tr)
    return {'p_lr': lr.predict_proba(X_te)[:, 1], 'p_rf': rf.predict_proba(X_te)[:, 1]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n-seeds', type=int, default=1,
                     help='CV Protocol v2: number of repeated 5-fold splits '
                          '(seeds 42..42+N-1). Default 1 = single split, seed 42, '
                          'for backward compatibility.')
    args = ap.parse_args()
    seeds = DEFAULT_SEEDS[:args.n_seeds] if args.n_seeds > 1 else [42]

    os.makedirs(OUT_DIR, exist_ok=True)
    pmap = pd.read_csv(PATIENT_MAP)[['sample_id', 'patient_id']]

    for fm, fname, obsm_key in EMBEDDINGS:
        data = load_fm(fm, fname, obsm_key, pmap)
        if data is None:
            continue

        oof_by_seed = run_repeated_kfold(data['X'], data['y'], data['groups'],
                                          fit_predict_fn=fit_predict, seeds=seeds)

        base = dict(cell_id=data['cell_id'], sample_id=data['sample_id'],
                    patient_id=data['patient_id'], severity=data['severity'], y=data['y'])

        if args.n_seeds <= 1:
            out = pd.DataFrame({**base, **oof_by_seed[seeds[0]]})
            out_path = os.path.join(OUT_DIR, f'oof_patient_{fm}.csv')
            out.to_csv(out_path, index=False)
            print(f'  SAVED -> {out_path}', flush=True)
        else:
            for seed in seeds:
                out = pd.DataFrame({**base, **oof_by_seed[seed]})
                out_path = os.path.join(OUT_DIR, f'oof_patient_{fm}_seed{seed}.csv')
                out.to_csv(out_path, index=False)
            print(f'  SAVED {len(seeds)} seed files -> '
                  f'{OUT_DIR}/oof_patient_{fm}_seed<42..{seeds[-1]}>.csv', flush=True)


if __name__ == '__main__':
    main()
