#!/usr/bin/env python3
"""
Round-16 fix: round-15's harness (round15_delong_and_undamaged.py) has a
confirmed defect -- unit-tested in round16_harness_unittest.py -- that
inflates comparator AUC loss. A single frozen score (PCA-200's own OOF
ensemble score, raw AUC 0.7423) refit through the 10-seed x 5-fold
StratifiedGroupKFold protocol at C=0.3 with UNSTANDARDIZED features loses
0.0854 AUC (exactly reproducing round-15's reported reduced=PCA-200-alone
value of 0.6569); round-16's own simulation showed a healthy harness should
lose ~0.01. Diagnosis (round-16, confirmed): PCA-200's score sits on a
narrow probability scale (~0.04-0.76, std 0.146); L2 at C=0.3 applied to
UNSTANDARDIZED coefficients disproportionately shrinks a feature that needs
a large coefficient to matter, crushing the slope and letting fold-to-fold
intercept noise dominate the pooled OOF ranking.

Fix (round-16's suggested option 2, the standard remedy): StandardScaler
fit on each fold's TRAIN portion only, applied to train and test, before
the C=0.3 LogisticRegression fit -- proper CV hygiene, no leakage. Verified
in round16_harness_unittest.py: this drops the single-score loss from
0.0854 to 0.0182, close to the ~0.01 a healthy harness should show.

This script reruns the ENTIRE round-15 analysis (three reduced-comparator
configs x five FMs, exact DeLong) with the scaled harness, superseding
every AUC-level number round-15 reported. The round-15 DeLong significance
tests themselves are not necessarily wrong (round-16: a shared harness
distortion mostly cancels between reduced and full, making those tests
conservative rather than anti-conservative) but are rerun here anyway for
a clean, fully-corrected record rather than a partial patch.
"""
import os
import sys
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, '/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/code')
from repeated_cv import DEFAULT_SEEDS, run_repeated_kfold
from scipy.stats import norm, rankdata

OOF_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2'
OUT_DIR = os.path.join(OOF_DIR, 'paired_inference_full')
QC_CSV = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
          '31_donor_level_revision/results_patient_v2/qc_baseline/patient_qc_distributions.csv')
PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')
FMS = ['STATE', 'UCE', 'Geneformer', 'scPRINT-2', 'scGPT']
RNG_SEED = 20260915
N_BOOT = 5000

# ---------------------------------------------------------------------
# Exact DeLong -- verbatim copy of the function in
# scripts/24_option_c_canonical_input/code/primary_benchmark_option_c.py
# (this is the implementation already used elsewhere in the project for
# "exact DeLong"; reused here unmodified rather than reimplemented, to
# guarantee identical results to every other DeLong number in the paper --
# importing that module directly would pull in a scanpy dependency at
# module level, so the function body is copied instead).
# ---------------------------------------------------------------------

def delong_test(y_true, score_a, score_b):
    """DeLong test for two correlated ROC curves. Returns (auc_a, auc_b, z, p_two_sided)."""
    y = np.asarray(y_true).astype(int)
    pos = y == 1
    neg = y == 0
    n_pos = int(pos.sum())
    n_neg = int(neg.sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan"), float("nan"), float("nan"), float("nan"), float("nan")

    auc_a = roc_auc_score(y, score_a)
    auc_b = roc_auc_score(y, score_b)

    def fast_delong(scores):
        s_pos = scores[pos]
        s_neg = scores[neg]
        all_scores = np.concatenate([s_pos, s_neg])
        ranks = rankdata(all_scores, method="average")
        r_pos = ranks[:n_pos]
        r_neg = ranks[n_pos:]
        rank_pos_within = rankdata(s_pos, method="average")
        rank_neg_within = rankdata(s_neg, method="average")
        v10 = (r_pos - rank_pos_within) / n_neg
        v01 = 1 - (r_neg - rank_neg_within) / n_pos
        return v10, v01

    v10_a, v01_a = fast_delong(np.asarray(score_a))
    v10_b, v01_b = fast_delong(np.asarray(score_b))

    s10 = np.cov(np.vstack([v10_a, v10_b]), ddof=1) if n_pos > 1 else np.zeros((2, 2))
    s01 = np.cov(np.vstack([v01_a, v01_b]), ddof=1) if n_neg > 1 else np.zeros((2, 2))
    var = (s10[0, 0] + s10[1, 1] - 2 * s10[0, 1]) / n_pos + \
          (s01[0, 0] + s01[1, 1] - 2 * s01[0, 1]) / n_neg
    corr_num = s10[0, 1] / n_pos + s01[0, 1] / n_neg
    var_a = s10[0, 0] / n_pos + s01[0, 0] / n_neg
    var_b = s10[1, 1] / n_pos + s01[1, 1] / n_neg
    corr = corr_num / np.sqrt(var_a * var_b) if var_a > 0 and var_b > 0 else float('nan')
    if var <= 0:
        return auc_a, auc_b, float("nan"), float("nan"), corr
    z = (auc_a - auc_b) / np.sqrt(var)
    p = 2 * (1 - norm.cdf(abs(z)))
    return auc_a, auc_b, z, p, corr


def delong_paired_test(y, p_a, p_b):
    auc_a, auc_b, z, p, corr = delong_test(y, p_a, p_b)
    se = (auc_a - auc_b) / z if (z and not np.isnan(z) and z != 0) else float('nan')
    return {'auc_a': auc_a, 'auc_b': auc_b, 'delta': auc_a - auc_b,
            'se': abs(se) if not np.isnan(se) else float('nan'),
            'z': z, 'p': p, 'corr': corr}


def patient_bootstrap_se(y, p_a, p_b, n_boot=N_BOOT, seed=RNG_SEED):
    """Cross-check: patient-level paired bootstrap SE of AUC(p_a)-AUC(p_b)."""
    rng = np.random.default_rng(seed)
    y = np.asarray(y)
    p_a = np.asarray(p_a)
    p_b = np.asarray(p_b)
    n = len(y)
    deltas = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        yb = y[idx]
        if yb.sum() == 0 or yb.sum() == n:
            deltas[b] = np.nan
            continue
        deltas[b] = roc_auc_score(yb, p_a[idx]) - roc_auc_score(yb, p_b[idx])
    deltas = deltas[~np.isnan(deltas)]
    return float(np.std(deltas, ddof=1))


# ---------------------------------------------------------------------
# Sanity checks before trusting the implementation on real data
# ---------------------------------------------------------------------
def sanity_check():
    rng = np.random.default_rng(1)
    n = 240
    y = np.concatenate([np.ones(25), np.zeros(215)])
    p_a = y * 0.3 + rng.normal(0, 0.3, n)
    p_b = y * 0.15 + rng.normal(0, 0.3, n)
    auc_a_sklearn = roc_auc_score(y, p_a)
    auc_b_sklearn = roc_auc_score(y, p_b)
    res = delong_paired_test(y, p_a, p_b)
    assert abs(res['auc_a'] - auc_a_sklearn) < 1e-9, (res['auc_a'], auc_a_sklearn)
    assert abs(res['auc_b'] - auc_b_sklearn) < 1e-9, (res['auc_b'], auc_b_sklearn)
    boot_se = patient_bootstrap_se(y, p_a, p_b, n_boot=3000, seed=2)
    print(f'[sanity] AUC point estimates match sklearn exactly: '
          f'a={res["auc_a"]:.6f} vs {auc_a_sklearn:.6f}, '
          f'b={res["auc_b"]:.6f} vs {auc_b_sklearn:.6f}')
    print(f'[sanity] DeLong SE={res["se"]:.4f} vs bootstrap SE={boot_se:.4f} '
          f'(ratio={res["se"]/boot_se:.3f}, expect close to 1)')
    print(f'[sanity] DeLong p={res["p"]:.4f}, corr={res["corr"]:.3f}')
    assert 0.5 < res['se'] / boot_se < 2.0, 'DeLong SE and bootstrap SE disagree badly'
    print('[sanity] PASSED\n')


sanity_check()

# ---------------------------------------------------------------------
# Data loading (identical to round12/13/14)
# ---------------------------------------------------------------------
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
patient_id_arr = qc['patient_id'].to_numpy()
log_n_cells_vec = qc['log_n_cells'].to_numpy().reshape(-1, 1)


def load_fm_ensemble_scores(fm, patient_ids):
    seed_scores = []
    for seed in DEFAULT_SEEDS:
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
for fm in ['PCA-200'] + FMS:
    fm_scores = load_fm_ensemble_scores(fm, patient_ids_set)
    assert set(fm_scores.index) == patient_ids_set
    model_vecs[fm] = qc['patient_id'].map(fm_scores).to_numpy()
pca_vec = model_vecs['PCA-200'].reshape(-1, 1)


def make_fit_predict(X_full_feat):
    def fit_predict(X_tr, y_tr, X_te):
        scaler = StandardScaler().fit(X_tr)  # fit on TRAIN fold only, no leakage
        clf = LogisticRegression(C=0.3, solver='lbfgs', max_iter=2000)
        clf.fit(scaler.transform(X_tr), y_tr)
        return {'p_lr': clf.predict_proba(scaler.transform(X_te))[:, 1]}
    return fit_predict


def run_oof(X_feat):
    """Returns seed-ensemble-aggregated OOF prediction vector (one score per
    patient, mean over the 10 seeds' OOF predictions -- matches the paper's
    seed-ensemble OOF convention throughout)."""
    oof = run_repeated_kfold(X_feat, y_obs, groups=patient_id_arr,
                              fit_predict_fn=make_fit_predict(X_feat))
    mat = np.column_stack([oof[s]['p_lr'] for s in sorted(oof.keys())])
    return mat.mean(axis=1)


# --- Harness self-check (round-16 item 1's gating unit test) before trusting
# --- any downstream number: a single frozen score refit through this harness
# --- must reproduce its own AUC to within ~0.01-0.02, not lose 0.085 like the
# --- unscaled round-15 harness did.
raw_pca_auc = roc_auc_score(y_obs, pca_vec.ravel())
oof_pca_alone_check = run_oof(pca_vec)
oof_pca_alone_auc_check = roc_auc_score(y_obs, oof_pca_alone_check)
print(f'[harness self-check] raw PCA-200 AUC={raw_pca_auc:.4f}, scaled-harness OOF AUC='
      f'{oof_pca_alone_auc_check:.4f}, loss={raw_pca_auc - oof_pca_alone_auc_check:+.4f} '
      f'(round-15 unscaled harness lost +0.0854 on this exact test; round-16 unit test '
      f'showed scaling drops it to +0.0182)')
assert (raw_pca_auc - oof_pca_alone_auc_check) < 0.04, \
    'harness self-check FAILED: loss still too large even after scaling'
print('[harness self-check] PASSED\n')

REDUCED_CONFIGS = {
    'degraded_QC_cohort_PCA200': X_with_cohort,          # round-14's reduced (13 feat incl. 8 cohort dummies)
    'undamaged_PCA200_alone': np.zeros((240, 0)),         # placeholder, filled below
    'undamaged_PCA200_logNcells': log_n_cells_vec,        # PCA-200 + log1p(n_cells), 1 extra feat
}

results = []
for config_name, extra_covariates in REDUCED_CONFIGS.items():
    if config_name == 'degraded_QC_cohort_PCA200':
        X_reduced = np.hstack([X_with_cohort, pca_vec])
    elif config_name == 'undamaged_PCA200_alone':
        X_reduced = pca_vec.copy()
    elif config_name == 'undamaged_PCA200_logNcells':
        X_reduced = np.hstack([pca_vec, extra_covariates])
    print(f'\n=== reduced config: {config_name} (n_features={X_reduced.shape[1]}) ===')
    oof_reduced = run_oof(X_reduced)
    auc_reduced_sklearn = roc_auc_score(y_obs, oof_reduced)
    print(f'  reduced OOF AUC (seed-ensemble) = {auc_reduced_sklearn:.4f}')

    for fm in FMS:
        X_full = np.hstack([X_reduced, model_vecs[fm].reshape(-1, 1)])
        oof_full = run_oof(X_full)
        dl = delong_paired_test(y_obs, oof_full, oof_reduced)
        boot_se = patient_bootstrap_se(y_obs, oof_full, oof_reduced, n_boot=N_BOOT)
        print(f'  {fm}: full AUC={dl["auc_a"]:.4f}  delta={dl["delta"]:+.4f}  '
              f'DeLong SE={dl["se"]:.4f} (boot SE={boot_se:.4f})  '
              f'z={dl["z"]:.3f}  p={dl["p"]:.5f}  corr={dl["corr"]:.3f}')
        results.append({
            'reduced_config': config_name, 'fm': fm,
            'auc_reduced': dl['auc_b'], 'auc_full': dl['auc_a'],
            'delta_auc': dl['delta'], 'delong_se': dl['se'],
            'bootstrap_se': boot_se, 'z': dl['z'], 'delong_p': dl['p'],
            'pred_corr': dl['corr'],
        })

out_df = pd.DataFrame(results)
out_path = os.path.join(OUT_DIR, 'round15_delong_undamaged_comparator.csv')
out_df.to_csv(out_path, index=False)
print(f'\nSaved {out_path}')
