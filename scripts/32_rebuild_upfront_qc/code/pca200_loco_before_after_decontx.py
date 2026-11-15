"""
Before-vs-after-decontamination robustness rerun, step 1 of 2 (step 2 is
Geneformer): PCA-200+LR baseline on the within-COVID validation set
(assemble_within_covid_validation.py - Liu+Ren+COMBAT, 2,394,993 cells,
283 donors, 24 events), run through the LOCO harness (repeated_cv.py's
run_loco/summarize_loco) on BOTH the raw and decontaminated matrices, so any
difference in headline AUC or per-cohort ranking is attributable to
decontamination alone -- same cells, same genes, same donors, same folds
(cohort-as-group is deterministic, no seed to vary between the two runs).

fit_predict follows scripts/31_donor_level_revision/code/pca_baseline_patient_oof_v2.py's
convention (StandardScaler(with_mean=False) + 200-component SVD + LogisticRegression,
fit on the train fold only), applied here to the LOCO protocol instead of
StratifiedGroupKFold. Uses TruncatedSVD, not that script's sklearn.decomposition.PCA
-- see the fit_predict docstring for why (PCA's required dense conversion
OOM'd a 400G job at this dataset's 2.4M-cell scale; TruncatedSVD on the
already-non-centered scaled data is the correctness-preserving fix, not an
approximation).

Usage:
    python pca200_loco_before_after_decontx.py --which raw
    python pca200_loco_before_after_decontx.py --which decontx
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.io import mmread
from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, '/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/code')
from repeated_cv import run_loco, summarize_loco

VALIDATION_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/within_covid_validation'
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/pca200_loco_validation'
N_COMPONENTS = 200


def fit_predict(X_tr, y_tr, X_te):
    """TruncatedSVD, not sklearn.decomposition.PCA -- this validation set is
    2.4M cells (LOCO training folds up to ~1.86M cells x 16,603 genes), and
    PCA requires a dense array (.toarray()), which OOM'd a 400G job here
    2026-09-16 (peak ~400G, right at the ceiling). TruncatedSVD operates
    directly on the sparse matrix, no densification, and is mathematically
    the right tool anyway: combined with StandardScaler(with_mean=False)
    (already non-centering, for the same sparse-compatibility reason),
    SVD-on-scaled-uncentered-data IS what PCA-without-centering computes --
    this is a correctness-preserving fix, not an approximation. The
    established single-cohort baseline (pca_baseline_patient_oof_v2.py) is
    untouched; this is specific to this validation script's larger scale."""
    scaler = StandardScaler(with_mean=False).fit(X_tr)
    Xs_tr = scaler.transform(X_tr)
    Xs_te = scaler.transform(X_te)
    svd = TruncatedSVD(n_components=N_COMPONENTS, random_state=42)
    Z_tr = svd.fit_transform(Xs_tr)
    Z_te = svd.transform(Xs_te)
    clf = LogisticRegression(C=1.0, solver='lbfgs', max_iter=2000,
                              class_weight='balanced', n_jobs=4)
    clf.fit(Z_tr, y_tr)
    p = clf.predict_proba(Z_te)[:, 1]
    return {'p_lr': p}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--which', required=True, choices=['raw', 'decontx'])
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)

    print(f'Loading {args.which} counts.mtx (genes x cells) + cells.csv...', flush=True)
    mat = mmread(f'{VALIDATION_DIR}/{args.which}/counts.mtx').tocsc()
    cells = pd.read_csv(f'{VALIDATION_DIR}/{args.which}/cells.csv')
    assert mat.shape[1] == len(cells)
    print(f'  {mat.shape[1]} cells x {mat.shape[0]} genes', flush=True)

    X = mat.T.tocsr()  # cells x genes, matching fit_predict's expected orientation
    y = cells['event'].to_numpy().astype(int)
    patient_id = cells['donor_id'].to_numpy()
    cohort = cells['cohort'].to_numpy()
    print(f'  {len(set(patient_id))} donors, {y.astype(bool).sum()} fatal cells-worth '
          f'(donor-level event count computed after aggregation)', flush=True)

    print('Running LOCO (PCA-200+LR, fit on train folds only)...', flush=True)
    oof, held_out = run_loco(X, y, cohort, fit_predict)
    print(f'  held out cohorts in order: {held_out}', flush=True)

    report = summarize_loco(oof, patient_id, y, cohort, score_key='p_lr')
    report.to_csv(f'{OUT_DIR}/loco_report_{args.which}.csv', index=False)
    print(f'\n=== {args.which} LOCO report ===')
    print(report.to_string(index=False))
    print(f'\nSAVED -> {OUT_DIR}/loco_report_{args.which}.csv')


if __name__ == '__main__':
    main()
