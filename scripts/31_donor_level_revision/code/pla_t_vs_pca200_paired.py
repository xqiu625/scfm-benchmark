#!/usr/bin/env python3
"""
Proper paired comparison of the Platelet-to-T-cell ratio (AUC 0.7406, CI
[0.601, 0.860] on n=148/18-fatal) against the whole-transcriptome PCA-200
baseline (AUC 0.742 on the full n=240/25-fatal canonical cohort). These two
numbers are NOT directly comparable as reported -- different N, different
patient sets (Pla-T's 148 patients are a subset of PCA-200's 240, restricted
to Ren+Stephenson+Qiu-sepsis with a valid T-cell count). Comparing the two
point estimates directly ("0.740 is close to 0.742") would repeat the exact
overclaim pattern already caught twice in this thread (MarkerScore "above
chance", PseudobulkEN "equivalent") -- a valid comparison needs PCA-200's
score recomputed on the SAME 148 patients, then a genuinely paired exact
DeLong test + bootstrap CI of the difference.
"""
import os

import numpy as np
import pandas as pd
from scipy.stats import norm

OOF_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2'
PLA_T_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/bib_revision_tier2/pla_t_ratio'
PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/bib_revision_tier2/pla_t_ratio'

SEEDS = list(range(42, 52))
N_BOOT = 10000
RNG_SEED = 20260913

m = pd.read_csv(PATIENT_MAP, dtype={'sample_id': str, 'patient_id': str})
pmap = m[['sample_id', 'patient_id']]

# The Pla-T ratio's own patient set (148 patients) is the restriction -- load
# it first, then restrict PCA-200's ensemble to exactly these patients.
plat_seed_scores = []
for seed in SEEDS:
    df = pd.read_csv(os.path.join(PLA_T_DIR, f'oof_patient_PlaTRatio_seed{seed}.csv'),
                      dtype={'patient_id': str})
    df = df.set_index('patient_id')[['y', 'p_lr']].rename(columns={'p_lr': 'score'})
    plat_seed_scores.append(df)
plat_patients = plat_seed_scores[0].index
y_plat = plat_seed_scores[0]['y']
plat_ens = pd.DataFrame({'y': y_plat, 'score': np.mean([s['score'] for s in plat_seed_scores], axis=0)},
                         index=plat_patients)

pca_seed_scores = []
for seed in SEEDS:
    df = pd.read_csv(os.path.join(OOF_DIR, f'oof_PCA-200_seed{seed}.csv'), dtype={'sample_id': str})
    df = df.merge(pmap, on='sample_id', how='inner')
    df['patient_id'] = df['patient_id'].astype(str)
    df = df[df['patient_id'].isin(plat_patients)]
    g = df.groupby('patient_id').agg(y=('y', 'max'), score=('p_lr', 'mean'))
    pca_seed_scores.append(g)
common = pca_seed_scores[0].index
for s in pca_seed_scores[1:]:
    common = common.intersection(s.index)
print(f'Pla-T patients: {len(plat_patients)}; PCA-200-on-same-subset patients found: {len(common)}',
      flush=True)
missing = set(plat_patients) - set(common)
if missing:
    print(f'  {len(missing)} Pla-T patients have no PCA-200 score (not in canonical COVID cohort, '
          f'e.g. sepsis patients not evaluated by the within-COVID PCA-200 model): {sorted(missing)}',
          flush=True)

pca_ens_score = np.mean([s.loc[common, 'score'].to_numpy() for s in pca_seed_scores], axis=0)
pca_ens = pd.DataFrame({'y': pca_seed_scores[0].loc[common, 'y'], 'score': pca_ens_score}, index=common)

idx = plat_ens.index.intersection(pca_ens.index)
plat_ens, pca_ens = plat_ens.loc[idx], pca_ens.loc[idx]
y = plat_ens['y'].to_numpy()
s_plat = plat_ens['score'].to_numpy()
s_pca = pca_ens['score'].to_numpy()
n = len(y)
print(f'\nFinal paired set: n={n}, fatal={int(y.sum())}', flush=True)


def auc(yb, sb):
    pos = yb == 1
    return (sb[pos][:, None] > sb[~pos][None, :]).mean() + 0.5 * (sb[pos][:, None] == sb[~pos][None, :]).mean()


auc_plat, auc_pca = auc(y, s_plat), auc(y, s_pca)
print(f'On this SAME {n}-patient subset: Pla-T ratio AUC={auc_plat:.4f}, PCA-200 AUC={auc_pca:.4f}, '
      f'delta={auc_pca - auc_plat:+.4f}', flush=True)

# exact DeLong on this paired subset
pos, neg = np.where(y == 1)[0], np.where(y == 0)[0]
n_pos, n_neg = len(pos), len(neg)
V10, V01 = {}, {}
for name, s in [('PlaT', s_plat), ('PCA200', s_pca)]:
    sx, sy = s[pos], s[neg]
    psi = (sx[:, None] > sy[None, :]).astype(float) + 0.5 * (sx[:, None] == sy[None, :])
    V10[name] = psi.mean(axis=1)
    V01[name] = psi.mean(axis=0)
V10_mat = np.column_stack([V10['PlaT'], V10['PCA200']])
V01_mat = np.column_stack([V01['PlaT'], V01['PCA200']])
S10 = np.cov(V10_mat, rowvar=False, ddof=1)
S01 = np.cov(V01_mat, rowvar=False, ddof=1)
Cov = S10 / n_pos + S01 / n_neg
var_diff = Cov[0, 0] + Cov[1, 1] - 2 * Cov[0, 1]
se_exact = np.sqrt(max(var_diff, 0))
delta = auc_pca - auc_plat
z = delta / se_exact if se_exact > 0 else np.nan
p = 2 * (1 - norm.cdf(abs(z))) if not np.isnan(z) else np.nan
rho_exact = Cov[0, 1] / np.sqrt(Cov[0, 0] * Cov[1, 1])
print(f'Exact DeLong: delta={delta:+.4f}, SE={se_exact:.5f}, p={p:.5f}, rho={rho_exact:.3f}', flush=True)

# paired bootstrap of the difference
rng = np.random.default_rng(RNG_SEED)
diffs = []
for _ in range(N_BOOT):
    idx_b = rng.integers(0, n, size=n)
    yb = y[idx_b]
    if len(np.unique(yb)) < 2:
        continue
    diffs.append(auc(yb, s_pca[idx_b]) - auc(yb, s_plat[idx_b]))
diffs = np.array(diffs)
ci_lo, ci_hi = np.percentile(diffs, [2.5, 97.5])
print(f'Paired bootstrap 95% CI of (PCA-200 - Pla-T) difference, {len(diffs)} resamples: '
      f'[{ci_lo:+.4f}, {ci_hi:+.4f}]', flush=True)

pd.DataFrame({
    'n_patients': [n], 'n_fatal': [int(y.sum())],
    'auc_pla_t': [auc_plat], 'auc_pca200_same_subset': [auc_pca],
    'delta_pca_minus_plat': [delta], 'delong_SE': [se_exact], 'delong_p': [p], 'rho': [rho_exact],
    'boot_ci_lo': [ci_lo], 'boot_ci_hi': [ci_hi],
}).to_csv(os.path.join(OUT_DIR, 'pla_t_vs_pca200_paired.csv'), index=False)
print(f'\nSAVED {os.path.join(OUT_DIR, "pla_t_vs_pca200_paired.csv")}', flush=True)
