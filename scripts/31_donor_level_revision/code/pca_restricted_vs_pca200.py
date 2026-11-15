#!/usr/bin/env python3
"""
contamination_response_plan.md item 4, comparison step: exact DeLong +
patient-clustered bootstrap for PCA-200-restricted (non-platelet-lineage
markers removed before PCA) vs. the unrestricted PCA-200 baseline (0.7423,
Table 3), reusing baseline_ladder_delong.py's construction verbatim.

If PCA-200-restricted's AUC survives (indistinguishable from, or not
significantly below, unrestricted PCA-200), the headline transcriptomic
signal is not explained by cross-lineage contamination content.
"""
import os

import numpy as np
import pandas as pd
from scipy.stats import norm

OOF_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2'
PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')
OUT_DIR = ('/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/'
           'bib_revision_tier2/contamination_audit')
os.makedirs(OUT_DIR, exist_ok=True)

SEEDS = list(range(42, 52))
N_BOOT = 10000
RNG_SEED = 20260914
TABLE3_PCA200_ENSEMBLE_AUC = 0.7423  # sanity check target

m = pd.read_csv(PATIENT_MAP, dtype={'sample_id': str, 'patient_id': str})
covid_patients = set(m.loc[m['disease'] == 'COVID-19', 'patient_id'])
pmap = m[['sample_id', 'patient_id']]


def build_ensemble(file_prefix):
    seed_scores = []
    for seed in SEEDS:
        df = pd.read_csv(os.path.join(OOF_DIR, f'{file_prefix}_seed{seed}.csv'),
                          dtype={'sample_id': str})
        df = df.merge(pmap, on='sample_id', how='inner')
        df = df[df['patient_id'].isin(covid_patients)]
        g = df.groupby('patient_id').agg(y=('y', 'max'), score=('p_lr', 'mean'))
        seed_scores.append(g)
    common = seed_scores[0][[]]
    for s in seed_scores[1:]:
        common = common.join(s[[]], how='inner')
    y = seed_scores[0].loc[common.index, 'y']
    scores = np.mean([s.loc[common.index, 'score'].values for s in seed_scores], axis=0)
    return pd.DataFrame({'y': y.values, 'score': scores}, index=common.index)


ens_pca = build_ensemble('oof_PCA-200')
ens_restricted = build_ensemble('oof_PCA-200-restricted')

idx = ens_pca.index.intersection(ens_restricted.index)
ens_pca = ens_pca.loc[idx]
ens_restricted = ens_restricted.loc[idx]

y = ens_pca['y'].values
pos = np.where(y == 1)[0]
neg = np.where(y == 0)[0]
n_pos, n_neg = len(pos), len(neg)
print(f'n_patients={len(y)}, n_pos(fatal)={n_pos}, n_neg(survived)={n_neg}', flush=True)


def auc_of(scores):
    sx, sy = scores[pos], scores[neg]
    psi = (sx[:, None] > sy[None, :]).astype(float) + 0.5 * (sx[:, None] == sy[None, :])
    return psi, psi.mean(axis=1), psi.mean(axis=0), psi.mean()


psi_pca, V10_pca, V01_pca, auc_pca = auc_of(ens_pca['score'].values)
psi_res, V10_res, V01_res, auc_res = auc_of(ens_restricted['score'].values)

print(f'\nPCA-200 (unrestricted) ensemble-OOF AUC = {auc_pca:.4f} '
      f'(sanity check vs Table 3 0.7423: diff={abs(auc_pca - TABLE3_PCA200_ENSEMBLE_AUC):.4f})',
      flush=True)
print(f'PCA-200-restricted (non-platelet-lineage markers removed) ensemble-OOF AUC = {auc_res:.4f}',
      flush=True)
assert abs(auc_pca - TABLE3_PCA200_ENSEMBLE_AUC) < 0.01, "PCA-200 sanity check failed"

V10_mat = np.column_stack([V10_pca, V10_res])
V01_mat = np.column_stack([V01_pca, V01_res])
S10 = np.cov(V10_mat, rowvar=False, ddof=1)
S01 = np.cov(V01_mat, rowvar=False, ddof=1)
Cov = S10 / n_pos + S01 / n_neg

delta = auc_pca - auc_res
var_diff = Cov[0, 0] + Cov[1, 1] - 2 * Cov[0, 1]
se_exact = np.sqrt(max(var_diff, 0))
z = delta / se_exact if se_exact > 0 else np.nan
p = 2 * (1 - norm.cdf(abs(z))) if not np.isnan(z) else np.nan
rho_exact = Cov[0, 1] / np.sqrt(Cov[0, 0] * Cov[1, 1])

print(f'\nExact DeLong: PCA-200 - PCA-200-restricted = {delta:+.4f}, SE={se_exact:.5f}, '
      f'p={p:.5f}, rho(DeLong-implied)={rho_exact:.3f}', flush=True)

# Patient-clustered bootstrap for the restricted AUC, and paired bootstrap of the difference
rng = np.random.default_rng(RNG_SEED)
n = len(y)
boot_auc_res = []
boot_diff = []
s_pca = ens_pca['score'].values
s_res = ens_restricted['score'].values
for _ in range(N_BOOT):
    idx_boot = rng.integers(0, n, size=n)
    yb = y[idx_boot]
    if len(np.unique(yb)) < 2:
        continue
    pos_b = yb == 1
    sb_res = s_res[idx_boot]
    sb_pca = s_pca[idx_boot]
    auc_b_res = (sb_res[pos_b][:, None] > sb_res[~pos_b][None, :]).mean() + \
                0.5 * (sb_res[pos_b][:, None] == sb_res[~pos_b][None, :]).mean()
    auc_b_pca = (sb_pca[pos_b][:, None] > sb_pca[~pos_b][None, :]).mean() + \
                0.5 * (sb_pca[pos_b][:, None] == sb_pca[~pos_b][None, :]).mean()
    boot_auc_res.append(auc_b_res)
    boot_diff.append(auc_b_pca - auc_b_res)

boot_auc_res = np.array(boot_auc_res)
boot_diff = np.array(boot_diff)
ci_lo, ci_hi = np.percentile(boot_auc_res, [2.5, 97.5])
diff_lo, diff_hi = np.percentile(boot_diff, [2.5, 97.5])

print(f'PCA-200-restricted patient-clustered bootstrap 95% CI: [{ci_lo:.4f}, {ci_hi:.4f}] '
      f'({len(boot_auc_res)} resamples)', flush=True)
print(f'Paired bootstrap 95% CI of the difference (PCA-200 - restricted): [{diff_lo:.4f}, {diff_hi:.4f}]',
      flush=True)

out = pd.DataFrame([{
    'auc_pca200_unrestricted': auc_pca, 'auc_pca200_restricted': auc_res,
    'delta_auc': delta, 'delong_exact_SE': se_exact, 'delong_exact_p': p,
    'rho_exact_delong': rho_exact,
    'restricted_boot_ci_lo': ci_lo, 'restricted_boot_ci_hi': ci_hi,
    'paired_diff_boot_ci_lo': diff_lo, 'paired_diff_boot_ci_hi': diff_hi,
    'n_patients': len(y), 'n_pos': n_pos, 'n_neg': n_neg,
}])
out_path = os.path.join(OUT_DIR, 'pca200_restricted_vs_unrestricted.csv')
out.to_csv(out_path, index=False)
print(f'\nSAVED {out_path}', flush=True)
