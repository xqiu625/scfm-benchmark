#!/usr/bin/env python3
"""
Follow-up to baseline_ladder_delong.py, prompted by a second external review
(baseline_ladder_results (2).tex, 2026-09-13): a non-significant exact DeLong
test (p=0.997 on a Delta=0.0002 gap) is NOT equivalence -- absence of
evidence of difference is not evidence of absence of difference. Reporting
PseudobulkEN and PCA-200 as "essentially a perfect match" from that p-value
alone was an overclaim, structurally the same mistake as the earlier
MarkerScore "above chance" overclaim this thread already corrected once.

Proper fix, per the review's own proposed method: report the PAIRED
bootstrap distribution of the difference (not two marginal CIs), and test
equivalence with TOST (two one-sided tests) at a margin that is PRE-SPECIFIED
before looking at the interval -- +/-0.05 patient-level AUC, justified as
"smaller than any difference that would change a deployment decision" (the
review's own suggested justification, used here since it's sound and I have
no better basis to pick a different number after already having seen the
point estimate).

Method: patient-clustered bootstrap (10,000 resamples, same patient index
draw applied to BOTH models' scores per resample so the difference is
genuinely paired), collect the per-resample AUC(PCA-200) - AUC(PseudobulkEN)
difference, take its 90% CI (5th/95th percentile -- the correct interval for
a TOST at alpha=0.05), and check whether that CI falls entirely within
[-0.05, +0.05]. Also reports the measured rho (already known from
baseline_ladder_delong.py: 0.531) for cross-reference against the review's
exact-threshold table.
"""
import os

import numpy as np
import pandas as pd

OOF_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2'
PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')
OUT_DIR = ('/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/'
           'bib_revision_tier2/clinical_utility')
os.makedirs(OUT_DIR, exist_ok=True)

SEEDS = list(range(42, 52))
N_BOOT = 10000
RNG_SEED = 20260913
MARGIN = 0.05  # pre-specified equivalence margin, patient-level AUC

m = pd.read_csv(PATIENT_MAP, dtype={'sample_id': str, 'patient_id': str})
covid_patients = set(m.loc[m['disease'] == 'COVID-19', 'patient_id'])
pmap = m[['sample_id', 'patient_id']]


def build_ensemble(file_prefix, needs_pmap_merge, already_patient_level):
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


pca = build_ensemble('oof_PCA-200', True, False)
pb = build_ensemble('oof_patient_PseudobulkEN', False, True)
idx = pca.index.intersection(pb.index)
pca, pb = pca.loc[idx], pb.loc[idx]
y = pca['y'].to_numpy()
s_pca = pca['score'].to_numpy()
s_pb = pb['score'].to_numpy()
n = len(y)
print(f'n_patients={n}, n_fatal={int(y.sum())}', flush=True)


def auc(yb, sb):
    pos_b = yb == 1
    return (sb[pos_b][:, None] > sb[~pos_b][None, :]).mean() + \
           0.5 * (sb[pos_b][:, None] == sb[~pos_b][None, :]).mean()


rng = np.random.default_rng(RNG_SEED)
point_auc_pca = auc(y, s_pca)
point_auc_pb = auc(y, s_pb)
point_delta = point_auc_pca - point_auc_pb
print(f'Point estimates: PCA-200={point_auc_pca:.4f}, PseudobulkEN={point_auc_pb:.4f}, '
      f'delta={point_delta:+.4f}', flush=True)

diffs = []
for _ in range(N_BOOT):
    idx_boot = rng.integers(0, n, size=n)
    yb = y[idx_boot]
    if len(np.unique(yb)) < 2:
        continue
    a1 = auc(yb, s_pca[idx_boot])
    a2 = auc(yb, s_pb[idx_boot])
    diffs.append(a1 - a2)
diffs = np.array(diffs)

ci90_lo, ci90_hi = np.percentile(diffs, [5, 95])
ci95_lo, ci95_hi = np.percentile(diffs, [2.5, 97.5])
print(f'\nPaired bootstrap difference (PCA-200 - PseudobulkEN), {len(diffs)} resamples:', flush=True)
print(f'  point delta = {point_delta:+.4f}', flush=True)
print(f'  90% CI (for TOST at alpha=0.05) = [{ci90_lo:+.4f}, {ci90_hi:+.4f}]', flush=True)
print(f'  95% CI (for reference) = [{ci95_lo:+.4f}, {ci95_hi:+.4f}]', flush=True)

tost_pass = (ci90_lo > -MARGIN) and (ci90_hi < MARGIN)
print(f'\nTOST at pre-specified margin +/-{MARGIN}: '
      f'{"PASS -- equivalence established" if tost_pass else "FAIL -- equivalence NOT established"}',
      flush=True)
if not tost_pass:
    # what margin WOULD the data support?
    supported_margin = max(abs(ci90_lo), abs(ci90_hi))
    print(f'  Margin the 90% CI actually supports: +/-{supported_margin:.4f}', flush=True)

pd.DataFrame({
    'point_delta': [point_delta],
    'boot_ci90_lo': [ci90_lo], 'boot_ci90_hi': [ci90_hi],
    'boot_ci95_lo': [ci95_lo], 'boot_ci95_hi': [ci95_hi],
    'prespecified_margin': [MARGIN], 'tost_pass': [tost_pass],
    'n_boot': [len(diffs)],
}).to_csv(os.path.join(OUT_DIR, 'pseudobulk_pca200_equivalence.csv'), index=False)
print(f'\nSAVED {os.path.join(OUT_DIR, "pseudobulk_pca200_equivalence.csv")}', flush=True)
