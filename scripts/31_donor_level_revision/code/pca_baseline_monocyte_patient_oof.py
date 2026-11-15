#!/usr/bin/env python3
"""
Job F Phase 4: PCA-200+LR baseline, patient-level, per-cell OOF, on the
monocyte canonical raw counts -- the in-distribution counterpart to
pca_baseline_patient_oof_v2.py (Job H2 step 3, platelet side). Same PCA dim,
classifier, and CV settings; no Disease/severity filter needed here since
build_monocyte_canonical.py already built exactly the target population.

CV Protocol v2: --n-seeds N (default 10) runs repeated StratifiedGroupKFold
across N seeds via repeated_cv.py.

Output dir: 05_RESULTS/revision_donor_level_v2/monocyte/
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

RAW = '/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c_monocyte/monocyte_canonical_raw.h5ad'
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2/monocyte'
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
    return {'p_lr': p, 'p_rf': p}  # LR-only baseline; p_rf duplicated for schema compat


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n-seeds', type=int, default=10)
    args = ap.parse_args()
    seeds = DEFAULT_SEEDS[:args.n_seeds] if args.n_seeds > 1 else [42]

    os.makedirs(OUT_DIR, exist_ok=True)
    print('Loading monocyte canonical raw counts...', flush=True)
    a = sc.read_h5ad(RAW)
    print(f'  {a.n_obs} cells x {a.n_vars} genes', flush=True)

    X = a.X
    y = a.obs['y'].to_numpy().astype(int)
    groups = a.obs['patient_id'].astype(str).to_numpy()
    cell_id = a.obs_names.values.astype(str)
    sample_id = a.obs['Sample_ID'].astype(str).to_numpy()
    source_study = a.obs['source_study'].astype(str).to_numpy()
    provenance = a.obs['monocyte_source_provenance'].astype(str).to_numpy()
    n_patients = len(np.unique(groups))
    print(f'  {X.shape[0]} cells, {n_patients} patients, {int(y.sum())} fatal', flush=True)

    oof_by_seed = run_repeated_kfold(X, y, groups, fit_predict_fn=fit_predict, seeds=seeds)

    base = dict(cell_id=cell_id, sample_id=sample_id, patient_id=groups, y=y,
                source_study=source_study, provenance=provenance)

    if args.n_seeds <= 1:
        out = pd.DataFrame({**base, **oof_by_seed[seeds[0]]})
        out_path = os.path.join(OUT_DIR, 'oof_monocyte_PCA-200.csv')
        out.to_csv(out_path, index=False)
        print(f'SAVED -> {out_path}', flush=True)
    else:
        for seed in seeds:
            out = pd.DataFrame({**base, **oof_by_seed[seed]})
            out_path = os.path.join(OUT_DIR, f'oof_monocyte_PCA-200_seed{seed}.csv')
            out.to_csv(out_path, index=False)
        print(f'SAVED {len(seeds)} seed files -> '
              f'{OUT_DIR}/oof_monocyte_PCA-200_seed<42..{seeds[-1]}>.csv', flush=True)


if __name__ == '__main__':
    main()
