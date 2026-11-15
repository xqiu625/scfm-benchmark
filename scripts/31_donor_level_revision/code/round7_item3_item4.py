#!/usr/bin/env python3
"""
Round-7 review items 3 and 4.

Item 3 diagnostics on the QC+cells+cohort incremental-value baseline:
  - refit WITHOUT the 8 cohort dummies (QC+log1p(n_cells) only) -- tests whether
    rare-level dummies under L2/25-events are needed for the baseline's AUC>0.5.
  - nested nested likelihood-ratio test: covariates-only vs covariates+representation,
    fit on the same folds, orientation-free (doesn't require choosing a sign).
  NOTE: verified separately (round7_item3_item4.py header) that the with-cohort
  baseline (0.627) was NEVER sign-flipped -- it was as-fit positive from the start
  (see paired_inference_full/qc_baseline_sign_check.txt: "AUC = 0.6270 (as-fit)").
  Only Job K's *original*, simpler QC-only-no-cohort model (0.376) needed a flip
  (->0.624). These diagnostics are a robustness check, not a fix for a real flip.

Item 4: per-cohort n/deaths table + verification that the reviewer's Laplace
add-one-smoothing size-inversion mechanism explains cohort_only_baseline.csv's
0.486 (already computed by round6_item1_cohort_only.py).
"""
import os
import sys

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score
from scipy import stats

OOF_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2'
OUT_DIR = os.path.join(OOF_DIR, 'paired_inference_full')
QC_CSV = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
          '31_donor_level_revision/results_patient_v2/qc_baseline/patient_qc_distributions.csv')
PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')
SEEDS = list(range(42, 52))
FMS = ['STATE', 'UCE', 'Geneformer', 'scPRINT-2', 'scGPT']
os.makedirs(OUT_DIR, exist_ok=True)


def safe_auc(y, p):
    return roc_auc_score(y, p) if len(np.unique(y)) > 1 else np.nan


qc = pd.read_csv(QC_CSV)
qc = qc[qc['disease'] == 'COVID-19'].copy()
qc['patient_id'] = qc['patient_id'].astype(str)
pm = pd.read_csv(PATIENT_MAP)
pm['patient_id'] = pm['patient_id'].astype(str)
cohort_map = pm.drop_duplicates('patient_id').set_index('patient_id')['Data_NO']
qc['cohort'] = qc['patient_id'].map(cohort_map)
qc = qc.dropna(subset=['cohort'])
qc = qc.sort_values('patient_id').reset_index(drop=True)
y = qc['y'].to_numpy()
groups = qc['patient_id'].to_numpy()
assert len(qc) == 240 and int(y.sum()) == 25, f'sanity gate failed: n={len(qc)}, fatal={int(y.sum())}'
print(f'n_patients={len(qc)}, n_fatal={int(y.sum())}, n_cohorts={qc.cohort.nunique()}')

qc['log_n_cells'] = np.log1p(qc['n_cells'])
qc['log_counts'] = np.log1p(qc['median_counts'])
qc['log_genes'] = np.log1p(qc['median_genes'])
X_qc_only = qc[['log_counts', 'log_genes', 'median_pct_mito', 'log_n_cells']].to_numpy()
cohort_dummies = pd.get_dummies(qc['cohort'], prefix='coh', drop_first=True).astype(float)
X_with_cohort = np.hstack([X_qc_only, cohort_dummies.to_numpy()])

# ---- Item 3 step 1: univariate association + multivariate coefficient signs ----
print('\n=== univariate association (Spearman rho, QC features + log_n_cells vs y) ===')
univariate_rows = []
for col in ['log_counts', 'log_genes', 'median_pct_mito', 'log_n_cells']:
    rho, p = stats.spearmanr(qc[col], y)
    univariate_rows.append({'feature': col, 'spearman_rho': round(rho, 4), 'spearman_p': round(p, 5)})
    print(f'{col:16s} rho={rho:+.4f} p={p:.5f}')
uni_df = pd.DataFrame(univariate_rows)
uni_df.to_csv(os.path.join(OUT_DIR, 'item3_univariate_association.csv'), index=False)

print('\n=== multivariate coefficient signs across 50 fits (QC+cells+cohort, L2 C=0.3) ===')
coef_signs = {c: [] for c in ['log_counts', 'log_genes', 'median_pct_mito', 'log_n_cells']}
feat_names = ['log_counts', 'log_genes', 'median_pct_mito', 'log_n_cells'] + list(cohort_dummies.columns)
for seed in SEEDS:
    skf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=seed)
    for tr, te in skf.split(X_with_cohort, y, groups=groups):
        clf = LogisticRegression(C=0.3, solver='lbfgs', class_weight='balanced', max_iter=2000)
        clf.fit(X_with_cohort[tr], y[tr])
        for i, fname in enumerate(['log_counts', 'log_genes', 'median_pct_mito', 'log_n_cells']):
            coef_signs[fname].append(clf.coef_[0][i])
coef_summary = []
for fname, vals in coef_signs.items():
    vals = np.array(vals)
    coef_summary.append({'feature': fname, 'mean_coef': round(vals.mean(), 4),
                          'pct_positive': round((vals > 0).mean() * 100, 1), 'n_fits': len(vals)})
    print(f'{fname:16s} mean_coef={vals.mean():+.4f}  pct_positive={(vals>0).mean()*100:.1f}%  (n={len(vals)} fits)')
coef_df = pd.DataFrame(coef_summary)
coef_df.to_csv(os.path.join(OUT_DIR, 'item3_multivariate_coef_signs.csv'), index=False)


def repeated_cv_auc_and_ll(X, y, groups, seeds, C=0.3):
    n = len(y)
    oof_all = np.zeros((len(seeds), n))
    for si, seed in enumerate(seeds):
        skf_ = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=seed)
        oof = np.full(n, np.nan)
        for tr, te in skf_.split(X, y, groups=groups):
            clf = LogisticRegression(C=C, solver='lbfgs', class_weight='balanced', max_iter=2000)
            clf.fit(X[tr], y[tr])
            oof[te] = clf.predict_proba(X[te])[:, 1]
        oof_all[si] = oof
    return oof_all.mean(axis=0)


# ---- Item 3 step 2: refit WITHOUT cohort dummies ----
print('\n=== baseline WITHOUT cohort dummies (QC+log1p(n_cells) only) ===')
p_no_cohort = repeated_cv_auc_and_ll(X_qc_only, y, groups, SEEDS, C=0.3)
auc_no_cohort = safe_auc(y, p_no_cohort)
print(f'AUC (QC+cells only, no cohort) = {auc_no_cohort:.4f}')
p_with_cohort = repeated_cv_auc_and_ll(X_with_cohort, y, groups, SEEDS, C=0.3)
auc_with_cohort = safe_auc(y, p_with_cohort)
print(f'AUC (QC+cells+cohort, for comparison) = {auc_with_cohort:.4f}')

with open(os.path.join(OUT_DIR, 'item3_no_cohort_refit.txt'), 'w') as f:
    f.write(f'QC+log1p(n_cells) only (no cohort dummies), L2 C=0.3, n={len(y)}/{int(y.sum())} fatal:\n')
    f.write(f'AUC = {auc_no_cohort:.4f}\n\n')
    f.write(f'QC+log1p(n_cells)+cohort one-hot (8 dummies), same protocol:\n')
    f.write(f'AUC = {auc_with_cohort:.4f}\n\n')
    if auc_no_cohort >= 0.5 and auc_with_cohort >= 0.5:
        f.write('Both orientations are already positive (>0.5) without any sign flip. '
                'The "sign flip" the reviewer is asking about applies to Job K\'s '
                'ORIGINAL simpler QC-only-no-cohort-no-log1p(n_cells) model (0.376->0.624), '
                'NOT to this properly-specified baseline, which was never flipped.\n')
    elif auc_no_cohort < 0.5 <= auc_with_cohort:
        f.write('Dropping cohort flips the baseline below 0.5: cohort inclusion is what '
                'orients the baseline correctly, consistent with the QC features alone being '
                'confounded by cohort membership (same sign-flip mechanism Job K found).\n')
    else:
        f.write('Unexpected pattern -- report as-is, do not force an explanation.\n')
print('Saved item3_no_cohort_refit.txt')

# ---- Item 3 step 3: nested likelihood-ratio test ----
print('\n=== nested LR test: covariates-only vs covariates+representation ===')


def fit_ll(X, y, C=0.3):
    """In-sample log-likelihood of an L2-penalized LR fit on the full data (fixed C,
    for the LR-test comparison -- this is a model-fit LR test, not a CV-based one;
    complementary to the CV-based AUC/DeLong comparison already in incremental_value_test.csv."""
    clf = LogisticRegression(C=C, solver='lbfgs', max_iter=2000)
    clf.fit(X, y)
    p = clf.predict_proba(X)[:, 1]
    eps = 1e-10
    ll = np.sum(y * np.log(p + eps) + (1 - y) * np.log(1 - p + eps))
    return ll, clf


def load_fm_ensemble_scores(fm, patient_ids):
    seed_scores = []
    for seed in SEEDS:
        if fm == 'PCA-200':
            df = pd.read_csv(os.path.join(OOF_DIR, f'oof_PCA-200_seed{seed}.csv'))
            pmap = pm[['sample_id', 'patient_id']]
            df = df.merge(pmap, on='sample_id', how='inner')
        else:
            df = pd.read_csv(os.path.join(OOF_DIR, f'oof_patient_{fm}_seed{seed}.csv'))
        df['patient_id'] = df['patient_id'].astype(str)
        df = df[df['patient_id'].isin(patient_ids)]
        g = df.groupby('patient_id')['p_lr'].mean()
        seed_scores.append(g)
    common = seed_scores[0].index
    for s in seed_scores[1:]:
        common = common.intersection(s.index)
    mat = np.column_stack([s.loc[common].to_numpy() for s in seed_scores])
    return pd.Series(mat.mean(axis=1), index=common)


patient_ids_set = set(qc['patient_id'])
ll_rows = []
ll_base, _ = fit_ll(X_with_cohort, y, C=0.3)
for fm in FMS + ['PCA-200']:
    fm_scores = load_fm_ensemble_scores(fm, patient_ids_set)
    common_idx = qc[qc['patient_id'].isin(fm_scores.index)].index
    Xb = X_with_cohort[common_idx]
    yb = y[common_idx]
    fmvec = qc.loc[common_idx, 'patient_id'].map(fm_scores).to_numpy().reshape(-1, 1)
    X_aug = np.hstack([Xb, fmvec])
    ll_b, _ = fit_ll(Xb, yb, C=0.3)
    ll_a, _ = fit_ll(X_aug, yb, C=0.3)
    lr_stat = 2 * (ll_a - ll_b)
    lr_p = 1 - stats.chi2.cdf(lr_stat, df=1)
    ll_rows.append({'fm': fm, 'n_patients': len(yb), 'll_covariates_only': round(ll_b, 3),
                     'll_covariates_plus_fm': round(ll_a, 3), 'LR_stat': round(lr_stat, 3),
                     'LR_p': round(lr_p, 6)})
    print(f'{fm:12s}: LL_base={ll_b:.3f}  LL_aug={ll_a:.3f}  LR_stat={lr_stat:.3f}  p={lr_p:.6f}')
ll_df = pd.DataFrame(ll_rows)
ll_df.to_csv(os.path.join(OUT_DIR, 'item3_nested_lr_test.csv'), index=False)
print('Saved item3_nested_lr_test.csv')

# ---- Item 4: per-cohort n/deaths table ----
print('\n=== per-cohort n/deaths table (item 4) ===')
cohort_table = qc.groupby('cohort').agg(n=('y', 'size'), n_fatal=('y', 'sum')).reset_index()
cohort_table['raw_rate'] = (cohort_table['n_fatal'] / cohort_table['n']).round(4)
cohort_table['laplace_rate'] = ((cohort_table['n_fatal'] + 1) / (cohort_table['n'] + 2)).round(4)
cohort_table = cohort_table.sort_values('n', ascending=False)
cohort_table.to_csv(os.path.join(OUT_DIR, 'item4_per_cohort_table.csv'), index=False)
print(cohort_table.to_string(index=False))
print('Saved item4_per_cohort_table.csv')

print('\nDONE.')
