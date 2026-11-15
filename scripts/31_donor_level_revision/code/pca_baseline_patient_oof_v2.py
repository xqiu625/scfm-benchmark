#!/usr/bin/env python3
"""
Job H2 step 3: PCA-200+LR baseline, patient-level, per-cell OOF (the
outstanding gap from the original Job H -- the existing
scripts/24_option_c_canonical_input/code/pca_baseline_mortality.py only
emits fold-level AUC, not per-cell OOF, so it can't drop into
donor_level_reanalysis.py's donor/patient aggregation).

Same raw-count input, PCA dim (200), and CV settings as that script, but:
  - explicit Disease in {COVID-19, Sepsis} filter added (the original script's
    valid mask is severity-only; without a Disease filter it would also pull
    in Flu/Lung_infect/Non_covid/SLE cells that carry a matching severity
    label -- not the same universe the FMs are evaluated on)
  - StratifiedGroupKFold groups = patient_id_v2 (not Sample_ID)
  - PCA/scaler fit on the train fold only, applied to test (no leakage)
  - emits per-cell OOF p_lr in the oof_<name>.csv schema
    (cell_id, sample_id, severity, y, p_lr, p_rf) so it drops straight into
    donor_level_reanalysis.py; p_rf is a duplicate of p_lr (this is an LR-only
    baseline -- no RF variant was requested or fit; flagged here so the
    downstream "PCA-200 RF" row in any table is understood as non-substantive)

CV Protocol v2 (2026-09-04): --n-seeds N (default 1) runs repeated
StratifiedGroupKFold across N seeds via repeated_cv.py, writing one
oof_PCA-200_seed<NN>.csv per seed when N>1, or oof_PCA-200.csv when N=1.

Output dir: 05_RESULTS/revision_donor_level_v2/
"""
import argparse
import os

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.decomposition import PCA
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
N_COMPONENTS = 200

# Job H2 Step 0 (2026-09-03, Task-4 audit vs Combes 2021 Supp Table 1 p.3):
# same fix as rerun_covid_sepsis_mortality_patient_v2.py -- applied at
# read-time, not by mutating the shared raw h5ad.
SEVERITY_FIX = {
    'GSM4995436': 'moderate',
    'GSM4995442': 'severe',
    'GSM4995432': 'moderate',
}


def fit_predict(X_tr, y_tr, X_te):
    scaler = StandardScaler(with_mean=False).fit(X_tr)
    Xs_tr = scaler.transform(X_tr)
    Xs_te = scaler.transform(X_te)
    pca = PCA(n_components=N_COMPONENTS, random_state=42)
    Z_tr = pca.fit_transform(Xs_tr.toarray() if hasattr(Xs_tr, 'toarray') else Xs_tr)
    Z_te = pca.transform(Xs_te.toarray() if hasattr(Xs_te, 'toarray') else Xs_te)
    clf = LogisticRegression(C=1.0, solver='lbfgs', max_iter=2000,
                              class_weight='balanced', n_jobs=4)
    clf.fit(Z_tr, y_tr)
    p = clf.predict_proba(Z_te)[:, 1]
    return {'p_lr': p, 'p_rf': p}  # LR-only baseline; p_rf duplicated for schema compat


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n-seeds', type=int, default=1,
                     help='CV Protocol v2: number of repeated 5-fold splits '
                          '(seeds 42..42+N-1). Default 1 = single split, seed 42.')
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
    df = pd.DataFrame({'cell_id': cell_id[keep], 'sample_id': sample_id[keep],
                        'severity': severity[keep]})
    df['_row'] = np.arange(len(df))
    n_before = len(df)
    df = df.merge(pmap, on='sample_id', how='inner')
    print(f'  v2 map filter: {len(df)} / {n_before} cells '
          f'(dropped {n_before - len(df)}, Combes pooled GSMs)', flush=True)

    X = a.X[keep][df['_row'].to_numpy()]
    y = (df['severity'].to_numpy() == 'fatal').astype(int)
    groups = df['patient_id'].to_numpy()
    n_patients = len(np.unique(groups))
    print(f'  {X.shape[0]} cells, {n_patients} patients, {int(y.sum())} fatal', flush=True)

    oof_by_seed = run_repeated_kfold(X, y, groups, fit_predict_fn=fit_predict, seeds=seeds)

    base = dict(cell_id=df['cell_id'].to_numpy(), sample_id=df['sample_id'].to_numpy(),
                severity=df['severity'].to_numpy(), y=y)

    if args.n_seeds <= 1:
        out = pd.DataFrame({**base, **oof_by_seed[seeds[0]]})
        out_path = os.path.join(OUT_DIR, 'oof_PCA-200.csv')
        out.to_csv(out_path, index=False)
        print(f'SAVED -> {out_path}', flush=True)
    else:
        for seed in seeds:
            out = pd.DataFrame({**base, **oof_by_seed[seed]})
            out_path = os.path.join(OUT_DIR, f'oof_PCA-200_seed{seed}.csv')
            out.to_csv(out_path, index=False)
        print(f'SAVED {len(seeds)} seed files -> '
              f'{OUT_DIR}/oof_PCA-200_seed<42..{seeds[-1]}>.csv', flush=True)


if __name__ == '__main__':
    main()
