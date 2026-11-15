#!/usr/bin/env python3
"""
Follow-up to baseline_ladder_auc.py, prompted by an external review of that
result (BIB_revision_memo (1).md item B status box + baseline_ladder_results.tex,
2026-09-13): the review correctly flagged that reporting MarkerScore (AUC
0.585) as "above chance" was an overclaim without a proper interval, and
that the MORE VALUABLE, statistically supportable finding is that the
design resolves the baseline-vs-reduced-complexity-representation gap even
though it cannot resolve between-foundation-model gaps -- CONDITIONAL on the
actual score correlation (rho) between each reduced baseline and PCA-200,
which was not yet measured ("measure rho before leaning on that rung").

This script reuses round5_delong_exact.py's exact-DeLong construction
(fast DeLong via placement values, same method used for every other paired
contrast in this paper, e.g. Supplementary tab:se_check) verbatim, extended
to MarkerScore and HVG-mean vs PCA-200, computing:
  - the true patient-clustered bootstrap CI for each new baseline's AUC
    (10,000 resamples, matching this project's own convention e.g.
    Supplementary tab:winner_ensemble), not the Hanley-McNeil point-estimate
    approximation the memo used as a placeholder;
  - the exact DeLong paired SE/p for MarkerScore-vs-PCA-200 and
    HVG-mean-vs-PCA-200;
  - the EXACT rho implied by the DeLong covariance structure for each pair
    (rho = Cov[i,j] / sqrt(Cov[i,i]*Cov[j,j])), which is the precise version
    of the "measured rho" resolution_tool.py's resolvable_delta_auc() needs
    -- not a proxy from the raw score Pearson correlation, though that's
    also reported for cross-checking.
"""
import os

import numpy as np
import pandas as pd
from scipy.stats import norm

OOF_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2'
PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')
OUT_DIR = ('/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/'
           'bib_revision_tier2/clinical_utility')
os.makedirs(OUT_DIR, exist_ok=True)

SEEDS = list(range(42, 52))
N_BOOT = 10000
RNG_SEED = 20260913

m = pd.read_csv(PATIENT_MAP, dtype={'sample_id': str, 'patient_id': str})
covid_patients = set(m.loc[m['disease'] == 'COVID-19', 'patient_id'])
pmap = m[['sample_id', 'patient_id']]


def build_ensemble(name, file_prefix, needs_pmap_merge, already_patient_level=False):
    seed_scores = []
    for seed in SEEDS:
        df = pd.read_csv(os.path.join(OOF_DIR, f'{file_prefix}_seed{seed}.csv'),
                          dtype={'patient_id': str, 'sample_id': str})
        if needs_pmap_merge:
            df = df.merge(pmap, on='sample_id', how='inner')
        df = df[df['patient_id'].isin(covid_patients)]
        if already_patient_level:
            g = df.set_index('patient_id')[['y', 'p_lr']].rename(columns={'p_lr': 'score'})
        else:
            g = df.groupby('patient_id').agg(y=('y', 'max'), score=('p_lr', 'mean'))
        seed_scores.append(g)
    common = seed_scores[0][[]]
    for s in seed_scores[1:]:
        common = common.join(s[[]], how='inner')
    y = seed_scores[0].loc[common.index, 'y']
    scores = np.mean([s.loc[common.index, 'score'].values for s in seed_scores], axis=0)
    return pd.DataFrame({'y': y.values, 'score': scores}, index=common.index)


METHODS = {
    'PCA-200': ('oof_PCA-200', True, False),
    'MarkerScore': ('oof_patient_MarkerScore', False, False),
    'HVG-mean': ('oof_patient_HVG-mean', False, False),
    'PseudobulkEN': ('oof_patient_PseudobulkEN', False, True),
}

ensembled = {name: build_ensemble(name, prefix, needs_merge, patient_level)
             for name, (prefix, needs_merge, patient_level) in METHODS.items()}
idx = ensembled['PCA-200'].index
for name in METHODS:
    ensembled[name] = ensembled[name].loc[idx]

y = ensembled['PCA-200']['y'].values
pos = np.where(y == 1)[0]
neg = np.where(y == 0)[0]
n_pos, n_neg = len(pos), len(neg)
print(f'n_patients={len(y)}, n_pos(fatal)={n_pos}, n_neg(survived)={n_neg}', flush=True)

models = list(METHODS.keys())
V10, V01, AUC = {}, {}, {}
for name in models:
    s = ensembled[name]['score'].values
    sx, sy = s[pos], s[neg]
    psi = (sx[:, None] > sy[None, :]).astype(float) + 0.5 * (sx[:, None] == sy[None, :])
    V10[name] = psi.mean(axis=1)
    V01[name] = psi.mean(axis=0)
    AUC[name] = psi.mean()

V10_mat = np.column_stack([V10[fm] for fm in models])
V01_mat = np.column_stack([V01[fm] for fm in models])
S10 = np.cov(V10_mat, rowvar=False, ddof=1)
S01 = np.cov(V01_mat, rowvar=False, ddof=1)
Cov = S10 / n_pos + S01 / n_neg

print('\nAUCs:', {fm: round(AUC[fm], 4) for fm in models}, flush=True)

# ------------------------------------------------------------------
# Patient-clustered bootstrap CI for MarkerScore and HVG-mean (10,000
# resamples of patients with replacement, matching this project's own
# convention e.g. Supplementary tab:winner_ensemble).
# ------------------------------------------------------------------
rng = np.random.default_rng(RNG_SEED)
n = len(y)
boot_rows = []
for name in ['MarkerScore', 'HVG-mean', 'PseudobulkEN']:
    s = ensembled[name]['score'].values
    boot_aucs = []
    for _ in range(N_BOOT):
        idx_boot = rng.integers(0, n, size=n)
        yb, sb = y[idx_boot], s[idx_boot]
        if len(np.unique(yb)) < 2:
            continue
        pos_b = yb == 1
        auc_b = (sb[pos_b][:, None] > sb[~pos_b][None, :]).mean() + \
                0.5 * (sb[pos_b][:, None] == sb[~pos_b][None, :]).mean()
        boot_aucs.append(auc_b)
    boot_aucs = np.array(boot_aucs)
    ci_lo, ci_hi = np.percentile(boot_aucs, [2.5, 97.5])
    boot_rows.append({'method': name, 'auc': AUC[name], 'boot_ci_lo': ci_lo, 'boot_ci_hi': ci_hi,
                       'n_boot': len(boot_aucs)})
    print(f'{name}: AUC={AUC[name]:.4f}, patient-clustered bootstrap 95% CI '
          f'[{ci_lo:.4f}, {ci_hi:.4f}] ({len(boot_aucs)} resamples)', flush=True)

boot_df = pd.DataFrame(boot_rows)
boot_path = os.path.join(OUT_DIR, 'baseline_ladder_bootstrap_ci.csv')
boot_df.to_csv(boot_path, index=False)
print(f'SAVED {boot_path}', flush=True)

# ------------------------------------------------------------------
# Exact DeLong pairwise test + implied rho, vs PCA-200.
# ------------------------------------------------------------------
print('\nExact DeLong vs PCA-200, and the implied rho for each pair:', flush=True)
delong_rows = []
i_pca = models.index('PCA-200')
for name in ['MarkerScore', 'HVG-mean', 'PseudobulkEN']:
    j = models.index(name)
    var_diff = Cov[i_pca, i_pca] + Cov[j, j] - 2 * Cov[i_pca, j]
    se_exact = np.sqrt(max(var_diff, 0))
    delta = AUC['PCA-200'] - AUC[name]
    z = delta / se_exact if se_exact > 0 else np.nan
    p = 2 * (1 - norm.cdf(abs(z))) if not np.isnan(z) else np.nan
    rho_exact = Cov[i_pca, j] / np.sqrt(Cov[i_pca, i_pca] * Cov[j, j])
    rho_raw_score = np.corrcoef(ensembled['PCA-200']['score'], ensembled[name]['score'])[0, 1]
    delong_rows.append({'pair': f'PCA-200 - {name}', 'delta_auc': delta, 'delong_exact_SE': se_exact,
                         'delong_exact_p': p, 'rho_exact_delong': rho_exact,
                         'rho_raw_score_pearson': rho_raw_score})
    print(f'PCA-200 - {name}: delta={delta:+.4f}  DeLong_exact_SE={se_exact:.5f}  p={p:.5f}  '
          f'rho(DeLong-implied)={rho_exact:.3f}  rho(raw score Pearson)={rho_raw_score:.3f}', flush=True)

delong_df = pd.DataFrame(delong_rows)
delong_path = os.path.join(OUT_DIR, 'baseline_ladder_delong_vs_pca200.csv')
delong_df.to_csv(delong_path, index=False)
print(f'SAVED {delong_path}', flush=True)

# ------------------------------------------------------------------
# Resolution check using resolution_tool.py's exact formula, with the
# MEASURED rho (DeLong-implied) instead of an assumed value.
# ------------------------------------------------------------------
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from resolution_tool import resolvable_delta_auc

print('\nResolution check with the MEASURED rho (not assumed):', flush=True)
for row in delong_rows:
    name = row['pair'].split(' - ')[1]
    rho = row['rho_exact_delong']
    resolvable = resolvable_delta_auc(n_pos, n_neg, rho=max(rho, 0), auc=0.72)
    gap = abs(row['delta_auc'])
    verdict = 'RESOLVABLE' if gap > resolvable else 'NOT resolvable at 80% power'
    print(f'  {name}: observed gap={gap:.3f}, resolvable-at-measured-rho({rho:.3f})={resolvable:.3f} '
          f'-> {verdict}', flush=True)
