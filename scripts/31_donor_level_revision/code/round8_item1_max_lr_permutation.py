#!/usr/bin/env python3
"""
Round-8 item 1: max-statistic permutation on the six nested-LR statistics
(QC+cells+cohort covariates-only vs covariates+representation), as the
secondary/exploratory "does ANY representation add value" readout, distinct
from the pre-specified PCA-200-only primary (nested LR p=0.0344, verified
against item3_nested_lr_test.csv, unchanged by this job).

Mirrors jobS2_round4_primary.py's max-statistic permutation logic, but the
per-model statistic here is the nested LR_stat (in-sample chi2, df=1) instead
of an out-of-fold AUC delta -- reuses round7_item3_item4.py's exact baseline
construction (X_with_cohort, fm ensemble scores) and fit_ll definition.

nperm=10,000, patient-level label permutation (240 patients, one row each,
no group structure to preserve -- already patient-aggregated).
"""
import os
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from scipy import stats

OOF_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2'
OUT_DIR = os.path.join(OOF_DIR, 'paired_inference_full')
QC_CSV = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
          '31_donor_level_revision/results_patient_v2/qc_baseline/patient_qc_distributions.csv')
PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')
SEEDS = list(range(42, 52))
FMS = ['STATE', 'UCE', 'Geneformer', 'scPRINT-2', 'scGPT']
NPERM = 10000
RNG_SEED = 20260908
os.makedirs(OUT_DIR, exist_ok=True)

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
assert len(qc) == 240 and int(y_obs.sum()) == 25, f'sanity gate failed: n={len(qc)}, fatal={int(y_obs.sum())}'
print(f'n_patients={len(qc)}, n_fatal={int(y_obs.sum())}')

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


# Pre-load each model's fixed representation vector, aligned to qc's patient order.
patient_ids_set = set(qc['patient_id'])
model_vecs = {}
model_idx = {}
for fm in FMS + ['PCA-200']:
    fm_scores = load_fm_ensemble_scores(fm, patient_ids_set)
    common_idx = qc[qc['patient_id'].isin(fm_scores.index)].index.to_numpy()
    fmvec = qc.loc[common_idx, 'patient_id'].map(fm_scores).to_numpy().reshape(-1, 1)
    model_vecs[fm] = fmvec
    model_idx[fm] = common_idx
    print(f'{fm}: n_common={len(common_idx)}')

# ---- Observed max LR_stat (cross-check against item3_nested_lr_test.csv) ----
obs_rows = []
for fm in FMS + ['PCA-200']:
    idx = model_idx[fm]
    Xb, yb = X_with_cohort[idx], y_obs[idx]
    X_aug = np.hstack([Xb, model_vecs[fm]])
    ll_b = fit_ll(Xb, yb)
    ll_a = fit_ll(X_aug, yb)
    lr_stat = 2 * (ll_a - ll_b)
    obs_rows.append({'fm': fm, 'lr_stat_obs': lr_stat})
obs_df = pd.DataFrame(obs_rows)
max_obs = obs_df['lr_stat_obs'].max()
print('\nObserved LR stats (cross-check vs item3_nested_lr_test.csv):')
print(obs_df.to_string(index=False))
print(f'max_obs = {max_obs:.4f} (argmax = {obs_df.loc[obs_df.lr_stat_obs.idxmax(), "fm"]})')

# ---- Permutation null on the max LR statistic ----
rng = np.random.default_rng(RNG_SEED)
max_perm = np.empty(NPERM)
n = len(y_obs)
for b in range(NPERM):
    y_perm = rng.permutation(y_obs)
    lr_stats_b = []
    for fm in FMS + ['PCA-200']:
        idx = model_idx[fm]
        Xb = X_with_cohort[idx]
        yb = y_perm[idx]
        if yb.sum() == 0 or yb.sum() == len(yb):
            lr_stats_b.append(0.0)
            continue
        X_aug = np.hstack([Xb, model_vecs[fm]])
        ll_b = fit_ll(Xb, yb)
        ll_a = fit_ll(X_aug, yb)
        lr_stats_b.append(2 * (ll_a - ll_b))
    max_perm[b] = max(lr_stats_b)
    if (b + 1) % 1000 == 0:
        print(f'  {b+1}/{NPERM} permutations done', flush=True)

p_val = (1 + np.sum(max_perm >= max_obs)) / (1 + NPERM)
print(f'\nmax-statistic permutation P (does ANY representation add value) = {p_val:.5f}')
print(f'nperm={NPERM}, rng_seed={RNG_SEED}')

out = pd.DataFrame({
    'metric': ['max_obs_lr_stat', 'argmax_model', 'nperm', 'rng_seed', 'permutation_P'],
    'value': [max_obs, obs_df.loc[obs_df.lr_stat_obs.idxmax(), 'fm'], NPERM, RNG_SEED, p_val],
})
out.to_csv(os.path.join(OUT_DIR, 'item1_max_lr_permutation.csv'), index=False)
obs_df.to_csv(os.path.join(OUT_DIR, 'item1_max_lr_permutation_obs.csv'), index=False)
print('Saved item1_max_lr_permutation.csv, item1_max_lr_permutation_obs.csv')
