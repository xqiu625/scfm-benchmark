#!/usr/bin/env python3
"""
Round-6 item 1: cohort-only-predictor baseline (within-COVID, 240 patients/25 fatal).

Reviewer's point: the existing "null" (count+disease) has never included cohort
identity, and within COVID disease is constant so it is effectively a
cell-count-only null. LOBO's own fold rule (34 excluded patients all sit in
zero-death cohorts) implies a cohort-identity floor >=0.579 by a hand-computed
concordance argument. This script computes the REAL number: an out-of-fold
cohort-mortality-rate predictor under the identical repeated-CV protocol
(10 seeds x 5-fold patient-grouped StratifiedGroupKFold) used everywhere else
in this revision.

Predictor: for each fold, compute each training cohort's mean(y) (Laplace-
smoothed to avoid degenerate 0/1 rates for tiny cohorts: (sum_y+1)/(n+2));
score each test-fold patient by their own cohort's training-fold rate. This
is a legitimate out-of-fold estimate of "how predictive is cohort identity
alone", not a hand floor.

Sanity gate: n_patients=240, n_fatal=25.
Output: paired_inference_full/cohort_only_baseline.csv
"""
import os
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold

sys.path.insert(0, os.path.dirname(__file__))
from repeated_cv import DEFAULT_SEEDS, N_SPLITS

PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2/paired_inference_full'
os.makedirs(OUT_DIR, exist_ok=True)


def main():
    pm = pd.read_csv(PATIENT_MAP)
    covid = pm[pm['disease'] == 'COVID-19'].drop_duplicates('patient_id')
    n_patients, n_fatal = len(covid), int(covid['y'].sum())
    print(f'n_patients={n_patients}, n_fatal={n_fatal}')
    assert n_patients == 240, f'sanity gate failed: n_patients={n_patients} != 240'
    assert n_fatal == 25, f'sanity gate failed: n_fatal={n_fatal} != 25'

    patient_id = covid['patient_id'].to_numpy()
    y = covid['y'].to_numpy().astype(int)
    cohort = covid['Data_NO'].to_numpy()
    n_cohorts = len(np.unique(cohort))
    print(f'n_cohorts={n_cohorts}')

    aucs = []
    for seed in DEFAULT_SEEDS:
        skf = StratifiedGroupKFold(n_splits=N_SPLITS, shuffle=True, random_state=seed)
        oof = np.full(n_patients, np.nan)
        for tr_idx, te_idx in skf.split(np.zeros(n_patients), y, groups=patient_id):
            tr_cohort, tr_y = cohort[tr_idx], y[tr_idx]
            df_tr = pd.DataFrame({'cohort': tr_cohort, 'y': tr_y})
            rate = df_tr.groupby('cohort')['y'].agg(lambda s: (s.sum() + 1) / (len(s) + 2))
            overall_rate = (tr_y.sum() + 1) / (len(tr_y) + 2)
            te_cohort = cohort[te_idx]
            oof[te_idx] = [rate.get(c, overall_rate) for c in te_cohort]
        assert not np.isnan(oof).any()
        auc = roc_auc_score(y, oof)
        aucs.append(auc)
        print(f'seed={seed} cohort-only AUC={auc:.4f}')

    aucs = np.array(aucs)
    mean, sd = aucs.mean(), aucs.std(ddof=1)
    print(f'\nrepeat-mean={mean:.4f} sd={sd:.4f}')

    out = pd.DataFrame({
        'metric': ['cohort_only_baseline_auc'],
        'n_patients': [n_patients], 'n_fatal': [n_fatal], 'n_cohorts': [n_cohorts],
        'repeat_mean': [round(mean, 4)], 'repeat_sd': [round(sd, 4)],
        'n_seeds': [len(DEFAULT_SEEDS)],
        'method': ['out-of-fold per-cohort Laplace-smoothed mortality rate, '
                   '10-seed x 5-fold patient-grouped StratifiedGroupKFold'],
    })
    out.to_csv(os.path.join(OUT_DIR, 'cohort_only_baseline.csv'), index=False)
    per_seed = pd.DataFrame({'seed': DEFAULT_SEEDS, 'auc': aucs})
    per_seed.to_csv(os.path.join(OUT_DIR, 'cohort_only_baseline_per_seed.csv'), index=False)
    print(f'\nSAVED -> {OUT_DIR}/cohort_only_baseline.csv')


if __name__ == '__main__':
    main()
