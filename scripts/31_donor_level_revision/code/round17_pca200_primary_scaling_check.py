#!/usr/bin/env python3
"""
Round-17 item 1: does the pre-specified incremental-value primary test
(PCA-200 alone, nested LR, covariates-only vs covariates+PCA-200, C=0.3,
observed statistic 4.4772, permutation P<0.0001, round-10's corrected
construction) suffer from the same unstandardized-narrow-scale-feature
distortion round-16 found in the OOF harness?

This construction has NO cross-fold OOF pooling (round10_nested_lr_
permutation_null.py fits a single in-sample model on all 240 patients, both
for the observed statistic and for every permutation -- X never changes,
only y is permuted) -- so round-16's specific failure mode (fold-to-fold
intercept heterogeneity corrupting a POOLED probability ranking) has no
analogue here. But L2 at C=0.3 on unstandardized features could still
distort the FIT ITSELF even without pooling, and round-17 asks a related
question: does standardizing move the collapsed permutation null (round-10
found median 0.082 vs chi2(1)'s 0.455) substantially toward chi2(1),
confirming round-16's mechanism (L2-on-unstandardized-scale) as the same
phenomenon rather than "L2 in general"?

Runs both the ORIGINAL (unstandardized, matching round-10's exact
construction) and a STANDARDIZED variant (StandardScaler fit once on the
full 240-patient X -- X is fixed across permutations, only y is permuted,
so fitting the scaler once introduces no leakage), reports observed
statistic, empirical permutation null (2000 perms, cheaper than the full
10,000 since this is a diagnostic not a replacement primary test), and
compares each null's shape to chi2(1).
"""
import os
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from scipy import stats

OOF_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2'
QC_CSV = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
          '31_donor_level_revision/results_patient_v2/qc_baseline/patient_qc_distributions.csv')
PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')
SEEDS = list(range(42, 52))
NPERM = 2000
RNG_SEED = 20261017

qc = pd.read_csv(QC_CSV)
qc = qc[qc['disease'] == 'COVID-19'].copy()
qc['patient_id'] = qc['patient_id'].astype(str)
pm = pd.read_csv(PATIENT_MAP)
pm['patient_id'] = pm['patient_id'].astype(str)
cohort_map = pm.drop_duplicates('patient_id').set_index('patient_id')['Data_NO']
qc['cohort'] = qc['patient_id'].map(cohort_map)
qc = qc.dropna(subset=['cohort'])
qc = qc.sort_values('patient_id').reset_index(drop=True)
y_obs = qc['y'].to_numpy()
assert len(qc) == 240 and int(y_obs.sum()) == 25

qc['log_n_cells'] = np.log1p(qc['n_cells'])
qc['log_counts'] = np.log1p(qc['median_counts'])
qc['log_genes'] = np.log1p(qc['median_genes'])
X_qc_only = qc[['log_counts', 'log_genes', 'median_pct_mito', 'log_n_cells']].to_numpy()
cohort_dummies = pd.get_dummies(qc['cohort'], prefix='coh', drop_first=True).astype(float)
X_with_cohort = np.hstack([X_qc_only, cohort_dummies.to_numpy()])


def fit_ll(X, y, C=0.3):
    clf = LogisticRegression(C=C, solver='lbfgs', max_iter=2000)
    clf.fit(X, y)
    p = clf.predict_proba(X)[:, 1]
    eps = 1e-10
    return np.sum(y * np.log(p + eps) + (1 - y) * np.log(1 - p + eps))


def load_fm_ensemble_scores(fm, patient_ids):
    seed_scores = []
    for seed in SEEDS:
        df = pd.read_csv(os.path.join(OOF_DIR, f'oof_PCA-200_seed{seed}.csv'))
        df = df.merge(pm[['sample_id', 'patient_id']], on='sample_id', how='inner')
        df['patient_id'] = df['patient_id'].astype(str)
        df = df[df['patient_id'].isin(patient_ids)]
        seed_scores.append(df.groupby('patient_id')['p_lr'].mean())
    common = seed_scores[0].index
    for s in seed_scores[1:]:
        common = common.intersection(s.index)
    mat = np.column_stack([s.loc[common].to_numpy() for s in seed_scores])
    return pd.Series(mat.mean(axis=1), index=common)


patient_ids_set = set(qc['patient_id'])
fm_scores = load_fm_ensemble_scores('PCA-200', patient_ids_set)
common_idx = qc[qc['patient_id'].isin(fm_scores.index)].index.to_numpy()
fmvec = qc.loc[common_idx, 'patient_id'].map(fm_scores).to_numpy().reshape(-1, 1)
print(f'PCA-200 score scale: min={fmvec.min():.4f} max={fmvec.max():.4f} std={fmvec.std():.4f}')

Xb_full = X_with_cohort[common_idx]
yb_full = y_obs[common_idx]
X_aug_full = np.hstack([Xb_full, fmvec])

for label, scale in [('UNSTANDARDIZED (round-10 original)', False),
                      ('STANDARDIZED (StandardScaler, fit once on full X)', True)]:
    print(f'\n=== {label} ===')
    if scale:
        scaler_b = StandardScaler().fit(Xb_full)
        scaler_a = StandardScaler().fit(X_aug_full)
        Xb_use = scaler_b.transform(Xb_full)
        Xa_use = scaler_a.transform(X_aug_full)
    else:
        Xb_use = Xb_full
        Xa_use = X_aug_full

    ll_b_obs = fit_ll(Xb_use, yb_full, C=0.3)
    ll_a_obs = fit_ll(Xa_use, yb_full, C=0.3)
    lr_obs = 2 * (ll_a_obs - ll_b_obs)
    print(f'  Observed LR_stat = {lr_obs:.4f} (round-10 reported 4.4772 for unstandardized)')

    rng = np.random.default_rng(RNG_SEED)
    lr_perm = np.empty(NPERM)
    for b in range(NPERM):
        y_perm = rng.permutation(yb_full)
        if y_perm.sum() == 0 or y_perm.sum() == len(y_perm):
            lr_perm[b] = 0.0
            continue
        ll_b = fit_ll(Xb_use, y_perm, C=0.3)
        ll_a = fit_ll(Xa_use, y_perm, C=0.3)
        lr_perm[b] = 2 * (ll_a - ll_b)

    p_emp = (1 + np.sum(lr_perm >= lr_obs)) / (1 + NPERM)
    chi2_p = stats.chi2.sf(lr_obs, df=1)
    print(f'  Empirical permutation P ({NPERM} perms) = {p_emp:.5f}')
    print(f'  Null: median={np.median(lr_perm):.4f} (chi2(1) median=0.455), '
          f'95th={np.percentile(lr_perm,95):.4f} (chi2(1) 95th=3.841), '
          f'99th={np.percentile(lr_perm,99):.4f} (chi2(1) 99th=6.635)')
    print(f'  Asymptotic chi2(1) p for this observed stat = {chi2_p:.5f}')
