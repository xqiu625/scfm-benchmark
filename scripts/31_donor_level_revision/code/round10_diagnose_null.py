#!/usr/bin/env python3
"""
Round-10 item 1 diagnosis: is the empirical permutation null of the
nested-LR statistic (in round8_item1_max_lr_permutation.py's construction)
a legitimate L2-shrunk null, or a collapsed artifact (same bug class as the
retired max-statistic P=0.0001)?

Reuses round8_item1_max_lr_permutation.py's exact data loading, X_with_cohort,
model_vecs (FIXED representation, computed once from true-label OOF ensemble),
and fit_ll (in-sample penalized LL, C=0.3) -- verbatim, no changes -- but:
  (a) saves the MARGINAL per-model null (not just the max) for STATE and PCA-200
  (b) also computes an UNPENALIZED (C=1e6, effectively no regularization) version
      of the same nested-LR statistic, to isolate whether the "collapse" is a
      regularization effect (expected, legitimate) vs a genuine implementation bug
  (c) sanity-checks that model_vecs[fm] truly never changes across permutations
      (prints array id + a hash before/after the loop)
"""
import os
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
NPERM = 2000
RNG_SEED = 20261010

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
model_vecs = {}
model_idx = {}
for fm in ['STATE', 'PCA-200']:
    fm_scores = load_fm_ensemble_scores(fm, patient_ids_set)
    common_idx = qc[qc['patient_id'].isin(fm_scores.index)].index.to_numpy()
    fmvec = qc.loc[common_idx, 'patient_id'].map(fm_scores).to_numpy().reshape(-1, 1)
    model_vecs[fm] = fmvec
    model_idx[fm] = common_idx
    print(f'{fm}: n_common={len(common_idx)}, representation vec id={id(fmvec)}, '
          f'hash={hash(fmvec.tobytes())}')

results = {fm: {'C0.3': [], 'C1e6': []} for fm in ['STATE', 'PCA-200']}
obs = {}
for fm in ['STATE', 'PCA-200']:
    idx = model_idx[fm]
    Xb, yb = X_with_cohort[idx], y_obs[idx]
    X_aug = np.hstack([Xb, model_vecs[fm]])
    for C, key in [(0.3, 'C0.3'), (1e6, 'C1e6')]:
        ll_b = fit_ll(Xb, yb, C=C)
        ll_a = fit_ll(X_aug, yb, C=C)
        obs[(fm, key)] = 2 * (ll_a - ll_b)

print('\nObserved LR stats (should match: STATE 5.7596 @ C=0.3, PCA-200 4.4772 @ C=0.3):')
for k, v in obs.items():
    print(f'  {k}: {v:.4f}')

rng = np.random.default_rng(RNG_SEED)
vec_id_before = {fm: id(model_vecs[fm]) for fm in ['STATE', 'PCA-200']}
vec_hash_before = {fm: hash(model_vecs[fm].tobytes()) for fm in ['STATE', 'PCA-200']}

for b in range(NPERM):
    y_perm = rng.permutation(y_obs)
    for fm in ['STATE', 'PCA-200']:
        idx = model_idx[fm]
        Xb = X_with_cohort[idx]
        yb = y_perm[idx]
        if yb.sum() == 0 or yb.sum() == len(yb):
            for key in ['C0.3', 'C1e6']:
                results[fm][key].append(0.0)
            continue
        X_aug = np.hstack([Xb, model_vecs[fm]])
        for C, key in [(0.3, 'C0.3'), (1e6, 'C1e6')]:
            ll_b = fit_ll(Xb, yb, C=C)
            ll_a = fit_ll(X_aug, yb, C=C)
            results[fm][key].append(2 * (ll_a - ll_b))
    if (b + 1) % 200 == 0:
        print(f'  {b+1}/{NPERM} done', flush=True)

print('\nSanity: representation vector unchanged across the loop?')
for fm in ['STATE', 'PCA-200']:
    same_id = id(model_vecs[fm]) == vec_id_before[fm]
    same_hash = hash(model_vecs[fm].tobytes()) == vec_hash_before[fm]
    print(f'  {fm}: same object id={same_id}, same content hash={same_hash}')

print('\nMarginal null summary:')
for fm in ['STATE', 'PCA-200']:
    for key in ['C0.3', 'C1e6']:
        arr = np.array(results[fm][key])
        p_emp = (1 + np.sum(arr >= obs[(fm, key)])) / (1 + NPERM)
        print(f'  {fm} [{key}]: obs={obs[(fm, key)]:.4f}  '
              f'null median={np.median(arr):.4f}  null 95th={np.percentile(arr, 95):.4f}  '
              f'null max={arr.max():.4f}  frac_zero={np.mean(arr == 0.0):.3f}  '
              f'empirical_p={p_emp:.4f}')

out_rows = []
for fm in ['STATE', 'PCA-200']:
    for key in ['C0.3', 'C1e6']:
        arr = np.array(results[fm][key])
        out_rows.append({
            'fm': fm, 'penalty': key, 'obs_lr_stat': obs[(fm, key)],
            'null_median': np.median(arr), 'null_p95': np.percentile(arr, 95),
            'null_max': arr.max(), 'frac_zero': np.mean(arr == 0.0),
            'empirical_p': (1 + np.sum(arr >= obs[(fm, key)])) / (1 + NPERM),
            'nperm': NPERM,
        })
pd.DataFrame(out_rows).to_csv(
    os.path.join(OUT_DIR, 'round10_marginal_null_diagnosis.csv'), index=False)
print('\nSaved round10_marginal_null_diagnosis.csv')
