"""
Round-5 review item 1: incremental-value test.
Baseline: QC (median_counts/genes/pct_mito, log1p) + n_cells (log1p) + cohort (Data_NO one-hot),
L2-penalized LR, repeated 10-seed x 5-fold StratifiedGroupKFold (patient-grouped), within-COVID (240/25).
Augmented: baseline + each FM's own seed-ensemble OOF score (5 separate augmented models, one per FM,
plus one for PCA-200).
Also: QC baseline sign-oriented re-check, per-cohort QC-vs-mortality association, fold cohort balance.
"""
import numpy as np, pandas as pd, os, sys
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score
from scipy import stats

OOF_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2'
OUT_DIR = os.path.join(OOF_DIR, 'paired_inference_full')
QC_CSV = '/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/results_patient_v2/qc_baseline/patient_qc_distributions.csv'
PATIENT_MAP = '/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv'
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
print('cohort NA after map:', qc['cohort'].isna().sum())
qc = qc.dropna(subset=['cohort'])
print(f'n_patients={len(qc)}, n_fatal={int(qc.y.sum())}, n_cohorts={qc.cohort.nunique()}')

# ---- per-cohort QC-vs-mortality association (descriptive) ----
print('\n=== per-cohort QC medians by outcome (descriptive) ===')
rows = []
for coh, g in qc.groupby('cohort'):
    if g['y'].nunique() < 2 or len(g) < 4:
        rows.append({'cohort': coh, 'n': len(g), 'n_fatal': int(g.y.sum()),
                      'note': 'insufficient variation for within-cohort test'})
        continue
    fatal = g[g.y == 1]['median_counts']; surv = g[g.y == 0]['median_counts']
    try:
        stat, p = stats.mannwhitneyu(fatal, surv, alternative='two-sided')
    except Exception:
        p = np.nan
    rows.append({'cohort': coh, 'n': len(g), 'n_fatal': int(g.y.sum()),
                  'median_counts_fatal': fatal.median(), 'median_counts_surv': surv.median(),
                  'mwu_p': p})
cohort_qc_df = pd.DataFrame(rows)
cohort_qc_df.to_csv(os.path.join(OUT_DIR, 'qc_by_cohort_association.csv'), index=False)
print(cohort_qc_df.to_string(index=False))

# ---- fold cohort balance check (seed 42, illustrative) ----
qc_sorted = qc.sort_values('patient_id').reset_index(drop=True)
y = qc_sorted['y'].to_numpy()
groups = qc_sorted['patient_id'].to_numpy()
skf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)
fold_cohort_rows = []
for i, (tr, te) in enumerate(skf.split(qc_sorted, y, groups=groups)):
    te_cohorts = qc_sorted.iloc[te]['cohort'].value_counts().to_dict()
    fold_cohort_rows.append({'fold': i, 'n_test': len(te), 'n_fatal_test': int(y[te].sum()),
                              'cohorts_in_test': te_cohorts})
fold_df = pd.DataFrame(fold_cohort_rows)
fold_df.to_csv(os.path.join(OUT_DIR, 'fold_cohort_balance_seed42.csv'), index=False)
print('\n=== fold cohort balance (seed 42) ===')
print(fold_df.to_string(index=False))

# ---- build feature matrix: QC + n_cells (log1p) + cohort one-hot ----
qc_sorted['log_n_cells'] = np.log1p(qc_sorted['n_cells'])
qc_sorted['log_counts'] = np.log1p(qc_sorted['median_counts'])
qc_sorted['log_genes'] = np.log1p(qc_sorted['median_genes'])
cohort_dummies = pd.get_dummies(qc_sorted['cohort'], prefix='coh', drop_first=True).astype(float)
X_base = np.hstack([
    qc_sorted[['log_counts', 'log_genes', 'median_pct_mito', 'log_n_cells']].to_numpy(),
    cohort_dummies.to_numpy(),
])
print(f'\nBaseline feature matrix: {X_base.shape} (QC=3, n_cells=1, cohort_dummies={cohort_dummies.shape[1]})')

# ---- load each FM's seed-ensemble OOF score for the same patients (for augmented model) ----
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

patient_ids_set = set(qc_sorted['patient_id'])

def repeated_cv_auc(X, y, groups, seeds, C=1.0):
    """Returns per-patient ensemble OOF prob (averaged over seeds)."""
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

groups = qc_sorted['patient_id'].to_numpy()

# baseline model, regularized (C=0.3, since 13-14 covariates / 25 events is tight)
print('\n=== Fitting baseline (QC+cells+cohort), C=0.3 L2 ===')
p_base = repeated_cv_auc(X_base, y, groups, SEEDS, C=0.3)
auc_base = safe_auc(y, p_base)
print(f'baseline AUC (sign as-fit) = {auc_base:.4f}')
sign_note = 'as-fit' if auc_base >= 0.5 else 'BELOW 0.5, anti-predictive as originally reported'

# ---- item-1 augmented models: baseline + each FM score ----
results = []
for fm in FMS + ['PCA-200']:
    fm_scores = load_fm_ensemble_scores(fm, patient_ids_set)
    common_idx = qc_sorted[qc_sorted['patient_id'].isin(fm_scores.index)].index
    Xb = X_base[common_idx]
    yb = y[common_idx]
    gb = groups[common_idx]
    fmvec = qc_sorted.loc[common_idx, 'patient_id'].map(fm_scores).to_numpy().reshape(-1, 1)
    X_aug = np.hstack([Xb, fmvec])
    p_b = repeated_cv_auc(Xb, yb, gb, SEEDS, C=0.3)
    p_a = repeated_cv_auc(X_aug, yb, gb, SEEDS, C=0.3)
    auc_b = safe_auc(yb, p_b)
    auc_a = safe_auc(yb, p_a)
    # exact DeLong paired test between p_b and p_a on same patients
    pos = np.where(yb == 1)[0]; neg = np.where(yb == 0)[0]
    def placement(s):
        sx, sy = s[pos], s[neg]
        psi = (sx[:, None] > sy[None, :]).astype(float) + 0.5 * (sx[:, None] == sy[None, :])
        return psi.mean(axis=1), psi.mean(axis=0), psi.mean()
    v10_b, v01_b, auc_b2 = placement(p_b)
    v10_a, v01_a, auc_a2 = placement(p_a)
    S10 = np.cov(np.column_stack([v10_b, v10_a]), rowvar=False, ddof=1)
    S01 = np.cov(np.column_stack([v01_b, v01_a]), rowvar=False, ddof=1)
    Cov = S10/len(pos) + S01/len(neg)
    var_diff = Cov[0,0] + Cov[1,1] - 2*Cov[0,1]
    se = np.sqrt(max(var_diff, 0))
    delta = auc_a2 - auc_b2
    z = delta/se if se > 0 else np.nan
    from scipy.stats import norm
    p_delong = 2*(1-norm.cdf(abs(z))) if not np.isnan(z) else np.nan
    results.append({'fm': fm, 'n_patients': len(yb), 'n_fatal': int(yb.sum()),
                     'auc_baseline_QC_cohort': round(auc_b2,4), 'auc_augmented_with_FM': round(auc_a2,4),
                     'incremental_delta_auc': round(delta,4), 'delong_SE': round(se,5),
                     'delong_p': round(p_delong,5)})
    print(f'{fm:12s}: baseline={auc_b2:.4f}  +FM={auc_a2:.4f}  delta={delta:+.4f}  p={p_delong:.5f}')

res_df = pd.DataFrame(results)
res_df.to_csv(os.path.join(OUT_DIR, 'incremental_value_test.csv'), index=False)
print('\nSaved:', os.path.join(OUT_DIR, 'incremental_value_test.csv'))

# save baseline (sign-checked) summary
with open(os.path.join(OUT_DIR, 'qc_baseline_sign_check.txt'), 'w') as f:
    f.write(f'QC+cells+cohort baseline (L2 C=0.3, patient-level, within-COVID n={len(y)}/{int(y.sum())} fatal):\n')
    f.write(f'AUC = {auc_base:.4f} ({sign_note})\n')
    f.write(f'Note: this baseline (with cohort covariates) differs from Job K\'s original QC-only-no-cohort model (0.376/0.624 sign-flip).\n')
print('Saved sign-check note.')
