#!/usr/bin/env python3
"""
Round-12 item 1: recompute the HEADLINE contrast (best FM vs PCA-200) under
the same penalized nested-LR / deviance statistic that makes the QC
incremental-value family (§3a) so sensitive, to check whether the negative
ΔAUC headline (permutation P=0.181) and this more powerful statistic
disagree on the identical 240-patient comparison.

Declared BEFORE running (per round-10/round-12's explicit requirement not to
choose machinery after seeing results): winner-corrected max-statistic
permutation, exactly mirroring jobS2_round4_primary.py's winner-correction
logic, but with the per-model statistic swapped from ΔAUC to the nested-LR
deviance statistic used throughout §3a:

  reduced model  = QC(3) + log1p(n_cells) + cohort one-hot(8) + PCA-200 score
                   (13 features, C=0.3 penalized, in-sample fit)
  full model     = reduced + one FM's frozen OOF ensemble score (14 features)
  per-FM stat    = 2*(ll_full - ll_reduced)
  winner stat    = max over the 5 FMs (STATE, UCE, Geneformer, scPRINT-2, scGPT)
  null           = patient-level label permutation, all representation scores
                    frozen (verified in round10_diagnose_null.py to already
                    be the correct construction), 10,000 permutations total,
                    split 4-way (2500/task) across a SLURM array for wall-clock.

ΔAUC (the existing headline test) remains the pre-specified primary
regardless of this result -- this is reported BESIDE it, not as a
replacement, matching how the two incremental-value readings (AUC-race vs
nested-LR) are already both reported in §3a.
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
FMS = ['STATE', 'UCE', 'Geneformer', 'scPRINT-2', 'scGPT']
RNG_SEED = 20260908
NPERM_TOTAL = 10000
N_TASKS = 4
NPERM_PER_TASK = NPERM_TOTAL // N_TASKS

task_id = int(sys.argv[1]) if len(sys.argv) > 1 else int(os.environ.get('SLURM_ARRAY_TASK_ID', 0))

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
for fm in ['PCA-200'] + FMS:
    fm_scores = load_fm_ensemble_scores(fm, patient_ids_set)
    common_idx = qc[qc['patient_id'].isin(fm_scores.index)].index.to_numpy()
    fmvec = qc.loc[common_idx, 'patient_id'].map(fm_scores).to_numpy().reshape(-1, 1)  # FROZEN
    model_vecs[fm] = fmvec
    model_idx[fm] = common_idx
    assert len(common_idx) == 240, f'{fm}: n_common={len(common_idx)} != 240'

# all five FMs + PCA-200 share the same 240-patient index (COVID-only, full overlap)
idx = model_idx['PCA-200']
X_reduced_true = np.hstack([X_with_cohort[idx], model_vecs['PCA-200']])  # 13 features

# ---- Observed winner statistic (task 0 only prints/saves this, to avoid duplication) ----
ll_reduced_obs = fit_ll(X_reduced_true, y_obs, C=0.3)
obs_stats = {}
for fm in FMS:
    X_full = np.hstack([X_reduced_true, model_vecs[fm]])
    ll_full_obs = fit_ll(X_full, y_obs, C=0.3)
    obs_stats[fm] = 2 * (ll_full_obs - ll_reduced_obs)
max_obs = max(obs_stats.values())
argmax_obs = max(obs_stats, key=obs_stats.get)
print(f'[task {task_id}] Observed per-FM deviance stats: {obs_stats}')
print(f'[task {task_id}] max_obs = {max_obs:.4f} (argmax={argmax_obs})')

# ---- Permutation null for this task's chunk ----
rng = np.random.default_rng(RNG_SEED + task_id)
max_perm = np.empty(NPERM_PER_TASK)
for b in range(NPERM_PER_TASK):
    y_perm = rng.permutation(y_obs)
    ll_reduced = fit_ll(X_reduced_true, y_perm, C=0.3)
    stats_b = []
    for fm in FMS:
        X_full = np.hstack([X_reduced_true, model_vecs[fm]])
        ll_full = fit_ll(X_full, y_perm, C=0.3)
        stats_b.append(2 * (ll_full - ll_reduced))
    max_perm[b] = max(stats_b)
    if (b + 1) % 250 == 0:
        print(f'[task {task_id}] {b+1}/{NPERM_PER_TASK} done', flush=True)

out = pd.DataFrame({'task_id': [task_id] * NPERM_PER_TASK, 'max_perm_stat': max_perm})
out_path = os.path.join(OUT_DIR, f'round12_headline_deviance_perm_task{task_id}.csv')
out.to_csv(out_path, index=False)

if task_id == 0:
    obs_df = pd.DataFrame({'fm': list(obs_stats.keys()), 'lr_stat_obs': list(obs_stats.values())})
    obs_df.to_csv(os.path.join(OUT_DIR, 'round12_headline_deviance_obs.csv'), index=False)
    print('Saved round12_headline_deviance_obs.csv')

print(f'[task {task_id}] Saved {out_path}')
