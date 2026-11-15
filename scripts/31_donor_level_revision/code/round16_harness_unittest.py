#!/usr/bin/env python3
"""
Round-16 item 1: unit-test the nested-comparison harness (used throughout
rounds 12-15) against round-16's specific diagnosis -- a single frozen
score refit through the 10-seed x 5-fold OOF protocol should reproduce its
own AUC to within ~0.01 (round-16's own simulation: median loss -0.0086,
0/200 replicates losing as much as the observed -0.068). If it instead
loses ~0.07, that confirms the harness itself is defective, most likely
from fitting LogisticRegression on unstandardized inputs under L2 (C=0.3)
-- a raw probability-scale feature (range ~0.05-0.4) needs a large
coefficient to matter, which is disproportionately punished by the penalty,
crushing the slope and letting fold-to-fold intercept noise dominate the
pooled ranking.

Test 1: reproduce the observed defect on the REAL PCA-200 score (this is
exactly what round12-15's harness did).
Test 2: same construction but with StandardScaler fit on each fold's TRAIN
split only (proper CV hygiene), to see if standardizing removes the loss.
Test 3: round-16's own simulation, reproduced locally, as an independent
check that the failure mode is real and not specific to this one score.
"""
import sys
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

sys.path.insert(0, '/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/code')
from repeated_cv import DEFAULT_SEEDS, run_repeated_kfold

OOF_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2'
QC_CSV = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
          '31_donor_level_revision/results_patient_v2/qc_baseline/patient_qc_distributions.csv')
PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')

qc = pd.read_csv(QC_CSV)
qc = qc[qc['disease'] == 'COVID-19'].copy()
qc['patient_id'] = qc['patient_id'].astype(str)
pm = pd.read_csv(PATIENT_MAP)
pm['patient_id'] = pm['patient_id'].astype(str)
qc = qc.sort_values('patient_id').reset_index(drop=True)
y_obs = qc['y'].to_numpy()
patient_id_arr = qc['patient_id'].to_numpy()


def load_fm_ensemble_scores(fm, patient_ids):
    seed_scores = []
    for seed in DEFAULT_SEEDS:
        df = pd.read_csv(f'{OOF_DIR}/oof_PCA-200_seed{seed}.csv')
        df = df.merge(pm[['sample_id', 'patient_id']], on='sample_id', how='inner')
        df['patient_id'] = df['patient_id'].astype(str)
        df = df[df['patient_id'].isin(patient_ids)]
        seed_scores.append(df.groupby('patient_id')['p_lr'].mean())
    common = seed_scores[0].index
    for s in seed_scores[1:]:
        common = common.intersection(s.index)
    mat = np.column_stack([s.loc[common].to_numpy() for s in seed_scores])
    return pd.Series(mat.mean(axis=1), index=common)


pca_scores = load_fm_ensemble_scores('PCA-200', set(qc['patient_id']))
pca_vec = qc['patient_id'].map(pca_scores).to_numpy().reshape(-1, 1)
raw_auc = roc_auc_score(y_obs, pca_vec.ravel())
print(f'Test 1 setup: raw PCA-200 seed-ensemble score AUC = {raw_auc:.4f} '
      f'(scale: min={pca_vec.min():.4f}, max={pca_vec.max():.4f}, std={pca_vec.std():.4f})\n')


def fit_predict_unscaled(X_tr, y_tr, X_te):
    clf = LogisticRegression(C=0.3, solver='lbfgs', max_iter=2000)
    clf.fit(X_tr, y_tr)
    return {'p_lr': clf.predict_proba(X_te)[:, 1]}


def fit_predict_scaled(X_tr, y_tr, X_te):
    scaler = StandardScaler().fit(X_tr)
    clf = LogisticRegression(C=0.3, solver='lbfgs', max_iter=2000)
    clf.fit(scaler.transform(X_tr), y_tr)
    return {'p_lr': clf.predict_proba(scaler.transform(X_te))[:, 1]}


print('=== Test 1: reproduce round12-15 harness on real PCA-200 score (unscaled, C=0.3) ===')
oof = run_repeated_kfold(pca_vec, y_obs, groups=patient_id_arr, fit_predict_fn=fit_predict_unscaled)
mat = np.column_stack([oof[s]['p_lr'] for s in sorted(oof.keys())])
oof_score = mat.mean(axis=1)
oof_auc = roc_auc_score(y_obs, oof_score)
print(f'  OOF AUC = {oof_auc:.4f}  (raw = {raw_auc:.4f}, loss = {raw_auc - oof_auc:+.4f})')
print(f'  This should reproduce the round-15 reported value of 0.6569 for reduced=PCA-200-alone.\n')

print('=== Test 2: same construction WITH StandardScaler fit per-fold on train only ===')
oof2 = run_repeated_kfold(pca_vec, y_obs, groups=patient_id_arr, fit_predict_fn=fit_predict_scaled)
mat2 = np.column_stack([oof2[s]['p_lr'] for s in sorted(oof2.keys())])
oof_score2 = mat2.mean(axis=1)
oof_auc2 = roc_auc_score(y_obs, oof_score2)
print(f'  OOF AUC = {oof_auc2:.4f}  (raw = {raw_auc:.4f}, loss = {raw_auc - oof_auc2:+.4f})')
print(f'  If standardizing fixes it, this loss should be << the Test 1 loss, close to ~0.01.\n')

print('=== Test 3: round-16-style simulation, reproduced independently ===')
rng = np.random.default_rng(42)
n_reps = 200
losses_unscaled = []
losses_scaled = []
for rep in range(n_reps):
    # Simulate a frozen score with true AUC ~0.725 on this exact cohort structure
    y_sim = y_obs.copy()
    rng.shuffle(y_sim)  # keep same class balance, break real association (we inject our own)
    # Build a score correlated with y_sim to hit ~0.725 AUC: latent + noise
    latent = y_sim * 1.1 + rng.normal(0, 1.0, len(y_sim))
    score_sim = 1 / (1 + np.exp(-latent))  # squashed to (0,1), probability-like scale
    achieved_auc = roc_auc_score(y_sim, score_sim)
    sv = score_sim.reshape(-1, 1)
    oof_u = run_repeated_kfold(sv, y_sim, groups=patient_id_arr, fit_predict_fn=fit_predict_unscaled,
                                seeds=[42])  # single seed per rep for speed, 200 reps
    p_u = oof_u[42]['p_lr']
    auc_oof_u = roc_auc_score(y_sim, p_u)
    losses_unscaled.append(achieved_auc - auc_oof_u)
losses_unscaled = np.array(losses_unscaled)
print(f'  Unscaled: median loss = {np.median(losses_unscaled):+.4f}, '
      f'5th-95th pct = [{np.percentile(losses_unscaled,5):+.4f}, {np.percentile(losses_unscaled,95):+.4f}], '
      f'max loss = {losses_unscaled.max():+.4f}')
n_extreme = int(np.sum(losses_unscaled >= (raw_auc - oof_auc)))
print(f'  Replicates losing as much as Test 1 ({raw_auc - oof_auc:.4f}): {n_extreme}/{n_reps}')
