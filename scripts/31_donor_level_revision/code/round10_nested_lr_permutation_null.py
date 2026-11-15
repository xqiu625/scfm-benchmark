#!/usr/bin/env python3
"""
Round-10 item 1 fix: correctly-referenced nested-LR permutation test for the
QC increment-value family (6 representations: STATE, UCE, Geneformer,
scPRINT-2, scGPT, PCA-200).

Declared BEFORE running (per round-10's explicit requirement not to pick the
reference after seeing which one flips the conclusion): keep the exact same
as-fit model spec used throughout Sec 3a (covariates = QC 3 features +
log1p(n_cells) + 8 cohort dummies, L2 C=0.3, in-sample fit on the fixed,
already-cross-validated seed-ensemble OOF representation score) -- do NOT
drop the penalty. Replace only the reference distribution: instead of the
chi2(1) asymptotic reference (verified in round10_diagnose_null.py, on
STATE+PCA-200, to be a poor match for this penalized statistic -- the
unpenalized/C=1e6 null matches chi2(1) closely, the C=0.3 null does not),
use a proper label-permutation empirical null of the SAME penalized
statistic, with the representation score frozen (verified via object-id +
content-hash check in round10_diagnose_null.py to already be correctly
frozen -- this was never the round-9/round-8 bug's mechanism).

nperm=10,000 (matches the paper's other two permutation tests: primary
max-statistic and range-statistic, both rng_seed=20260908). One model per
array task (SLURM_ARRAY_TASK_ID 0-5) since a single model's 10k-perm loop
takes on the order of an hour serially on this hardware.
"""
import os
import sys
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

OOF_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2'
OUT_DIR = os.path.join(OOF_DIR, 'paired_inference_full')
QC_CSV = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
          '31_donor_level_revision/results_patient_v2/qc_baseline/patient_qc_distributions.csv')
PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')
SEEDS = list(range(42, 52))
NPERM = 10000
RNG_SEED = 20260908  # same convention as the primary/range-statistic tests
FMS_ALL = ['STATE', 'UCE', 'Geneformer', 'scPRINT-2', 'scGPT', 'PCA-200']

task_id = int(sys.argv[1]) if len(sys.argv) > 1 else int(os.environ.get('SLURM_ARRAY_TASK_ID', 0))
FM = FMS_ALL[task_id]
print(f'Running model: {FM} (task_id={task_id})')

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
        if fm == 'PCA-200':
            df = pd.read_csv(os.path.join(OOF_DIR, f'oof_PCA-200_seed{seed}.csv'))
            df = df.merge(pm[['sample_id', 'patient_id']], on='sample_id', how='inner')
        else:
            df = pd.read_csv(os.path.join(OOF_DIR, f'oof_patient_{fm}_seed{seed}.csv'))
        df['patient_id'] = df['patient_id'].astype(str)
        df = df[df['patient_id'].isin(patient_ids)]
        seed_scores.append(df.groupby('patient_id')['p_lr'].mean())
    common = seed_scores[0].index
    for s in seed_scores[1:]:
        common = common.intersection(s.index)
    mat = np.column_stack([s.loc[common].to_numpy() for s in seed_scores])
    return pd.Series(mat.mean(axis=1), index=common)


patient_ids_set = set(qc['patient_id'])
fm_scores = load_fm_ensemble_scores(FM, patient_ids_set)
common_idx = qc[qc['patient_id'].isin(fm_scores.index)].index.to_numpy()
fmvec = qc.loc[common_idx, 'patient_id'].map(fm_scores).to_numpy().reshape(-1, 1)  # FROZEN
print(f'{FM}: n_common={len(common_idx)}')

Xb_full = X_with_cohort[common_idx]
yb_full = y_obs[common_idx]
X_aug_full = np.hstack([Xb_full, fmvec])

ll_b_obs = fit_ll(Xb_full, yb_full, C=0.3)
ll_a_obs = fit_ll(X_aug_full, yb_full, C=0.3)
lr_obs = 2 * (ll_a_obs - ll_b_obs)
print(f'Observed LR_stat (C=0.3, in-sample, frozen representation) = {lr_obs:.4f}')

rng = np.random.default_rng(RNG_SEED + task_id)  # distinct stream per model, still reproducible
lr_perm = np.empty(NPERM)
for b in range(NPERM):
    y_perm = rng.permutation(yb_full)
    if y_perm.sum() == 0 or y_perm.sum() == len(y_perm):
        lr_perm[b] = 0.0
        continue
    ll_b = fit_ll(Xb_full, y_perm, C=0.3)
    ll_a = fit_ll(X_aug_full, y_perm, C=0.3)
    lr_perm[b] = 2 * (ll_a - ll_b)
    if (b + 1) % 1000 == 0:
        print(f'  {FM}: {b+1}/{NPERM} done', flush=True)

p_emp = (1 + np.sum(lr_perm >= lr_obs)) / (1 + NPERM)
print(f'\n{FM}: empirical permutation P = {p_emp:.5f} '
      f'(null median={np.median(lr_perm):.4f}, null 95th={np.percentile(lr_perm, 95):.4f})')

out = pd.DataFrame({
    'fm': [FM], 'lr_stat_obs': [lr_obs], 'nperm': [NPERM],
    'rng_seed': [RNG_SEED + task_id], 'empirical_P': [p_emp],
    'null_median': [np.median(lr_perm)], 'null_p95': [np.percentile(lr_perm, 95)],
    'null_max': [lr_perm.max()],
})
out_path = os.path.join(OUT_DIR, f'round10_nested_lr_permutation_{FM.replace("-", "")}.csv')
out.to_csv(out_path, index=False)
print(f'Saved {out_path}')
