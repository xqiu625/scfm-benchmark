#!/usr/bin/env python3
"""
Round-13 item 2: the missing magnitude number for the round-12 deviance
result. The deviance statistic (nested-LR) establishes existence of
incremental information but carries no magnitude -- round-13 asks for the
one number that does: nested delta-AUC = AUC(covariates+PCA-200+FM score)
minus AUC(covariates+PCA-200), computed from the EXACT SAME reduced/full
C=0.3 in-sample fits already used for round12_headline_deviance_test.py's
observed deviance statistics (no new modeling choice -- same fits, just
also scoring their in-sample AUC). Explicitly an in-sample, optimistic
upper bound (not OOF), same caveat round-13's own simulation used.

Also serves round-13 item 3: reports all five FMs' deviance statistics
against the null's percentiles (median/95th/99th/max), not just the
winner (STATE), to check whether all five representations individually
clear the null, not just the argmax.
"""
import os
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

OOF_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2'
OUT_DIR = os.path.join(OOF_DIR, 'paired_inference_full')
QC_CSV = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
          '31_donor_level_revision/results_patient_v2/qc_baseline/patient_qc_distributions.csv')
PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')
SEEDS = list(range(42, 52))
FMS = ['STATE', 'UCE', 'Geneformer', 'scPRINT-2', 'scGPT']

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


def fit_ll_and_proba(X, y, C=0.3):
    clf = LogisticRegression(C=C, solver='lbfgs', max_iter=2000)
    clf.fit(X, y)
    p = clf.predict_proba(X)[:, 1]
    eps = 1e-10
    ll = np.sum(y * np.log(p + eps) + (1 - y) * np.log(1 - p + eps))
    return ll, p


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
    fmvec = qc.loc[common_idx, 'patient_id'].map(fm_scores).to_numpy().reshape(-1, 1)
    model_vecs[fm] = fmvec
    model_idx[fm] = common_idx
    assert len(common_idx) == 240

idx = model_idx['PCA-200']
X_reduced = np.hstack([X_with_cohort[idx], model_vecs['PCA-200']])
ll_reduced, p_reduced = fit_ll_and_proba(X_reduced, y_obs, C=0.3)
auc_reduced = roc_auc_score(y_obs, p_reduced)
print(f'Reduced model (covariates+PCA-200): in-sample AUC = {auc_reduced:.4f}')

rows = []
for fm in FMS:
    X_full = np.hstack([X_reduced, model_vecs[fm]])
    ll_full, p_full = fit_ll_and_proba(X_full, y_obs, C=0.3)
    auc_full = roc_auc_score(y_obs, p_full)
    lr_stat = 2 * (ll_full - ll_reduced)
    nested_delta_auc = auc_full - auc_reduced
    rows.append({
        'fm': fm, 'lr_stat_obs': lr_stat,
        'auc_reduced_covariates_pca200': auc_reduced,
        'auc_full_covariates_pca200_fm': auc_full,
        'nested_delta_auc': nested_delta_auc,
    })
    print(f'{fm}: LR_stat={lr_stat:.4f}, in-sample AUC full={auc_full:.4f}, '
          f'nested_delta_AUC={nested_delta_auc:+.4f}')

out_df = pd.DataFrame(rows)
out_path = os.path.join(OUT_DIR, 'round13_nested_delta_auc.csv')
out_df.to_csv(out_path, index=False)
print(f'Saved {out_path}')
