#!/usr/bin/env python3
"""
Job S step 5: PCA-200+LR baseline on the Chongqing within-sepsis cohort,
patient-level, CV Protocol v2 (10 seeds x 5-fold StratifiedGroupKFold).
Same protocol as pca_baseline_patient_oof_v2.py (Job H2 step 3).

Output: 05_RESULTS/revision_donor_level_v2/oof_sepsis_PCA-200_seed<NN>.csv
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
       'chongqing_sepsis_platelet_raw_counts.h5ad')
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2'
N_COMPONENTS = 200


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
    return {'p_lr': p, 'p_rf': p}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n-seeds', type=int, default=10)
    args = ap.parse_args()
    seeds = DEFAULT_SEEDS[:args.n_seeds] if args.n_seeds > 1 else [42]

    os.makedirs(OUT_DIR, exist_ok=True)
    print('Loading raw counts...', flush=True)
    a = sc.read_h5ad(RAW)
    print(f'  {a.n_obs} cells x {a.n_vars} genes', flush=True)

    severity = a.obs['covid_severity'].astype(str).values
    sample_id = a.obs['Sample_ID'].astype(str).values
    cell_id = a.obs_names.values.astype(str)
    y = (severity == 'fatal').astype(int)
    n_patients = len(np.unique(sample_id))
    n_fatal = len(np.unique(sample_id[y == 1]))
    print(f'  {a.n_obs} cells, {n_patients} patients, {n_fatal} fatal', flush=True)

    oof_by_seed = run_repeated_kfold(a.X, y, sample_id, fit_predict_fn=fit_predict, seeds=seeds)
    base = dict(cell_id=cell_id, sample_id=sample_id, severity=severity, y=y)

    if args.n_seeds <= 1:
        out = pd.DataFrame({**base, **oof_by_seed[seeds[0]]})
        out.to_csv(os.path.join(OUT_DIR, 'oof_sepsis_PCA-200.csv'), index=False)
    else:
        for seed in seeds:
            out = pd.DataFrame({**base, **oof_by_seed[seed]})
            out.to_csv(os.path.join(OUT_DIR, f'oof_sepsis_PCA-200_seed{seed}.csv'), index=False)
        print(f'SAVED {len(seeds)} seed files -> oof_sepsis_PCA-200_seed<NN>.csv', flush=True)


if __name__ == '__main__':
    main()
