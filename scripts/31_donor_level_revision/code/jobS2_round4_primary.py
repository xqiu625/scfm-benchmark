#!/usr/bin/env python3
"""
Round-4 primary-test correction (HPCC_next_steps.md sec 12, mislabeled "Job S"
in the doc -- reusing an earlier job's name; referred to here as the round-4
primary-test correction to avoid confusion).

Fixes Job R's invalid max-statistic procedure (resample-then-max does not
construct a null distribution -- it re-estimates the winner's-curse-inflated
effect, not remove it; review simulation: Type I error ~0.36 under the true
null). Rebuilds the primary within-COVID FM-vs-PCA-200 contrast properly:

1. Permutation null (primary): shuffle mortality labels, nperm=10,000,
   P = (1+#{max_perm >= max_obs})/(1+nperm). Scoring-step null (frozen OOF
   scores, no retraining) -- weaker than a full retrain null, noted as such.
2. Recentered bootstrap (cross-check): B=10,000 patient resamples,
   max_j(delta_j^b - delta_j^obs), P = fraction(recentered max >= max_j delta_j^obs).
3. LOFO on the primary permutation test: 25 leave-one-fatal-out refits,
   nperm=10,000 each.
4. SE consistency check: 15 within-COVID pairs, bootstrap SE (recomputed,
   B=10,000, actual std of the delta distribution) vs DeLong-approximate
   paired SE (Hanley-McNeil per-model SE + empirical score correlation:
   SE_diff = sqrt(SE_A^2 + SE_B^2 - 2*r*SE_A*SE_B)). Flag bootstrap/DeLong < 0.8.

Reuses jobR_winner_test_ensemble.py's exact seed-ensemble construction
(same OOF_DIR/PATIENT_MAP/seeds/models).

Hard sanity gate: n_patients=240, n_fatal=25 on COVID-19-only cohort.

Outputs (OUT_DIR = .../paired_inference_full/):
  primary_permutation_null.csv
  primary_recentered_bootstrap.csv
  fragility_primary_lofo.csv
  se_check_delong.csv
"""
import os
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

OOF_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2'
OUT_DIR = os.path.join(OOF_DIR, 'paired_inference_full')
PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')
HANLEY_SE = os.path.join(OOF_DIR, 'paired_inference', 'hanley_mcneil_se.csv')
CORR_MAT = os.path.join(OOF_DIR, 'paired_inference_full', 'score_correlation_ensemble.csv')

SEEDS = list(range(42, 52))
FMS = ['STATE', 'UCE', 'Geneformer', 'scPRINT-2', 'scGPT']
NPERM = 10000
B = 10000
RNG_SEED = 20260908
N_PATIENTS_EXPECTED = 240
N_FATAL_EXPECTED = 25


def load_map():
    m = pd.read_csv(PATIENT_MAP)
    return set(m.loc[m['disease'] == 'COVID-19', 'patient_id'])


def load_seed_scores(name, seed, covid_patients, is_pca=False):
    if is_pca:
        path = os.path.join(OOF_DIR, f'oof_{name}_seed{seed}.csv')
        df = pd.read_csv(path)
        pmap = pd.read_csv(PATIENT_MAP)[['sample_id', 'patient_id']]
        df = df.merge(pmap, on='sample_id', how='inner')
    else:
        path = os.path.join(OOF_DIR, f'oof_patient_{name}_seed{seed}.csv')
        df = pd.read_csv(path)
    df = df[df['patient_id'].isin(covid_patients)]
    pat = df.groupby('patient_id').agg(y=('y', 'max'), score=('p_lr', 'mean'))
    return pat


def safe_auc(y, p):
    if len(np.unique(y)) < 2:
        return np.nan
    return roc_auc_score(y, p)


def build_ensemble(covid_patients):
    per_seed = {}
    for seed in SEEDS:
        per_seed[seed] = {}
        for fm in FMS:
            per_seed[seed][fm] = load_seed_scores(fm, seed, covid_patients)
        per_seed[seed]['PCA-200'] = load_seed_scores('PCA-200', seed, covid_patients, is_pca=True)
    models = FMS + ['PCA-200']
    ensembled = {}
    for m in models:
        allp = pd.concat([per_seed[s][m]['score'].rename(s) for s in SEEDS], axis=1)
        y_ = pd.concat([per_seed[s][m]['y'].rename(s) for s in SEEDS], axis=1).max(axis=1)
        ensembled[m] = pd.DataFrame({'y': y_, 'score': allp.mean(axis=1)})
    return ensembled, models


def max_delta_obs(y, pa_dict, pb):
    auc_pb = safe_auc(y, pb)
    return {fm: safe_auc(y, pa_dict[fm]) - auc_pb for fm in pa_dict}


def permutation_null(y, fm_scores, pca_score, rng, nperm=NPERM):
    """Shuffle y; recompute max_j(AUC_j - AUC_pca) each time."""
    n = len(y)
    obs_deltas = max_delta_obs(y, fm_scores, pca_score)
    obs_max = max(obs_deltas.values())
    max_perm = np.empty(nperm)
    for i in range(nperm):
        y_perm = rng.permutation(y)
        if len(np.unique(y_perm)) < 2:
            max_perm[i] = np.nan
            continue
        auc_pca = safe_auc(y_perm, pca_score)
        ds = [safe_auc(y_perm, fm_scores[fm]) - auc_pca for fm in fm_scores]
        max_perm[i] = max(ds)
    valid = ~np.isnan(max_perm)
    max_perm = max_perm[valid]
    p = (1 + np.sum(max_perm >= obs_max)) / (1 + len(max_perm))
    return obs_deltas, obs_max, p, max_perm


def recentered_bootstrap(y, fm_scores, pca_score, obs_deltas, obs_max, rng, b=B):
    """Recenter each resample's per-FM delta by the OBSERVED delta (not
    resample delta), then take max -- this is what makes it null-respecting,
    unlike Job R's flawed non-recentered version."""
    n = len(y)
    idx_all = np.arange(n)
    max_recentered = np.empty(b)
    for i in range(b):
        idx = rng.choice(idx_all, size=n, replace=True)
        yi = y[idx]
        if len(np.unique(yi)) < 2:
            max_recentered[i] = np.nan
            continue
        auc_pca_b = safe_auc(yi, pca_score[idx])
        ds = []
        for fm in fm_scores:
            delta_b = safe_auc(yi, fm_scores[fm][idx]) - auc_pca_b
            ds.append(delta_b - obs_deltas[fm])  # recenter
        max_recentered[i] = max(ds)
    valid = ~np.isnan(max_recentered)
    max_recentered = max_recentered[valid]
    p = float((max_recentered >= obs_max).mean())
    lo, hi = np.percentile(max_recentered, [2.5, 97.5])
    return p, lo, hi


def paired_bootstrap_sd(y, pa, pb, rng, b=B):
    """Same as jobL/jobR paired_bootstrap but also returns the SD of the
    delta distribution (not just percentile CI)."""
    n = len(y)
    idx_all = np.arange(n)
    deltas = np.empty(b)
    for i in range(b):
        idx = rng.choice(idx_all, size=n, replace=True)
        yi = y[idx]
        if len(np.unique(yi)) < 2:
            deltas[i] = np.nan
            continue
        deltas[i] = safe_auc(yi, pa[idx]) - safe_auc(yi, pb[idx])
    deltas = deltas[~np.isnan(deltas)]
    obs = safe_auc(y, pa) - safe_auc(y, pb)
    sd = float(np.std(deltas, ddof=1))
    lo, hi = np.percentile(deltas, [2.5, 97.5])
    p_two = min(2 * min((deltas <= 0).mean(), (deltas >= 0).mean()), 1.0)
    return obs, sd, lo, hi, p_two


def main():
    covid_patients = load_map()
    ensembled, models = build_ensemble(covid_patients)

    common_all = ensembled['STATE']
    for m in models[1:]:
        common_all = common_all.join(ensembled[m][[]], how='inner')
    n_patients = len(common_all)
    n_fatal = int(ensembled['STATE'].loc[common_all.index, 'y'].sum())
    print(f'Sanity gate: n_patients={n_patients}, n_fatal={n_fatal} '
          f'(expected {N_PATIENTS_EXPECTED}/{N_FATAL_EXPECTED})', flush=True)
    if n_patients != N_PATIENTS_EXPECTED or n_fatal != N_FATAL_EXPECTED:
        print('FATAL: sanity gate failed. Producing NO output files.', file=sys.stderr)
        sys.exit(1)

    os.makedirs(OUT_DIR, exist_ok=True)
    rng = np.random.default_rng(RNG_SEED)

    idx = common_all.index
    y = ensembled['STATE'].loc[idx, 'y'].to_numpy()
    scores = {m: ensembled[m].loc[idx, 'score'].to_numpy() for m in models}
    fm_scores = {fm: scores[fm] for fm in FMS}
    pca_score = scores['PCA-200']

    # ================= 1. permutation null (primary) =================
    print('\n=== (1) permutation null (primary) ===', flush=True)
    obs_deltas, obs_max, perm_p, max_perm = permutation_null(y, fm_scores, pca_score, rng, nperm=NPERM)
    print(f'observed per-model delta: {obs_deltas}')
    print(f'observed max delta = {obs_max:.4f}')
    print(f'permutation P = {perm_p:.5f} (nperm={NPERM})', flush=True)

    perm_rows = [{'row_type': 'per_model', 'model': fm,
                   'delta_auc_obs': round(obs_deltas[fm], 4)} for fm in FMS]
    perm_rows.append({'row_type': 'summary', 'model': 'max_over_5_FMs',
                       'delta_auc_obs': round(obs_max, 4),
                       'permutation_P': round(perm_p, 5), 'nperm': NPERM,
                       'rng_seed': RNG_SEED,
                       'note': 'scoring-step null: labels shuffled on frozen OOF scores, no retraining'})
    pd.DataFrame(perm_rows).to_csv(os.path.join(OUT_DIR, 'primary_permutation_null.csv'), index=False)

    # ================= 2. recentered bootstrap (cross-check) =================
    print('\n=== (2) recentered bootstrap (cross-check) ===', flush=True)
    recentered_p, rc_lo, rc_hi = recentered_bootstrap(y, fm_scores, pca_score, obs_deltas, obs_max, rng, b=B)
    print(f'recentered-bootstrap P = {recentered_p:.5f}, recentered-max 95% range [{rc_lo:.4f}, {rc_hi:.4f}]', flush=True)
    pd.DataFrame([{'obs_max_delta': round(obs_max, 4), 'recentered_bootstrap_P': round(recentered_p, 5),
                    'recentered_max_ci_lo': round(rc_lo, 4), 'recentered_max_ci_hi': round(rc_hi, 4),
                    'B': B, 'rng_seed': RNG_SEED}]).to_csv(
        os.path.join(OUT_DIR, 'primary_recentered_bootstrap.csv'), index=False)

    # ================= 3. LOFO on the primary permutation test =================
    print('\n=== (3) LOFO on primary permutation test (25 refits) ===', flush=True)
    fatal_patients = list(idx[y == 1])
    assert len(fatal_patients) == N_FATAL_EXPECTED
    lofo_rows = []
    p_list = []
    for drop_i, drop_patient in enumerate(fatal_patients):
        keep_mask = idx != drop_patient
        y_sub = y[keep_mask]
        fm_sub = {fm: fm_scores[fm][keep_mask] for fm in FMS}
        pca_sub = pca_score[keep_mask]
        _, obs_max_sub, p_sub, _ = permutation_null(y_sub, fm_sub, pca_sub, rng, nperm=NPERM)
        p_list.append(p_sub)
        lofo_rows.append({'row_type': 'refit', 'refit_idx': drop_i, 'dropped_patient': drop_patient,
                           'n_patients': int(keep_mask.sum()), 'obs_max_delta': round(obs_max_sub, 4),
                           'permutation_P': round(p_sub, 5)})
        print(f'  refit {drop_i}: drop={drop_patient}, obs_max={obs_max_sub:.4f}, P={p_sub:.5f}', flush=True)
    p_arr = np.array(p_list)
    lofo_rows.append({'row_type': 'summary', 'refit_idx': -1, 'dropped_patient': 'ALL_25',
                       'n_patients': N_PATIENTS_EXPECTED - 1,
                       'obs_max_delta': np.nan,
                       'permutation_P': f'range[{p_arr.min():.5f},{p_arr.max():.5f}]',
                       'n_p_gt_0.05': int((p_arr > 0.05).sum())})
    pd.DataFrame(lofo_rows).to_csv(os.path.join(OUT_DIR, 'fragility_primary_lofo.csv'), index=False)
    print(f'LOFO summary: P range [{p_arr.min():.5f}, {p_arr.max():.5f}], '
          f'{int((p_arr > 0.05).sum())}/25 refits with P>0.05', flush=True)

    # ================= 4. SE consistency check =================
    print('\n=== (4) SE consistency check (15 within-COVID pairs) ===', flush=True)
    hanley = pd.read_csv(HANLEY_SE).set_index('model')['hanley_mcneil_se'].to_dict()
    corr = pd.read_csv(CORR_MAT, index_col=0)
    all_models = FMS + ['PCA-200']
    se_rows = []
    for i in range(len(all_models)):
        for j in range(i + 1, len(all_models)):
            a, b_ = all_models[i], all_models[j]
            obs, boot_sd, lo, hi, p_two = paired_bootstrap_sd(y, scores[a], scores[b_], rng, b=B)
            se_a = hanley.get(a); se_b = hanley.get(b_)
            r = corr.loc[a, b_]
            if se_a is not None and se_b is not None:
                delong_se = float(np.sqrt(se_a**2 + se_b**2 - 2 * r * se_a * se_b))
            else:
                delong_se = np.nan
            ratio = boot_sd / delong_se if delong_se and not np.isnan(delong_se) and delong_se > 0 else np.nan
            se_rows.append({'pair': f'{a} - {b_}', 'delta_auc': round(obs, 4),
                             'bootstrap_SE': round(boot_sd, 5), 'delong_SE_approx': round(delong_se, 5)
                             if not np.isnan(delong_se) else np.nan,
                             'correlation': round(float(r), 4),
                             'bootstrap_over_delong_ratio': round(ratio, 3) if not np.isnan(ratio) else np.nan,
                             'flag_ratio_lt_0.8': bool(ratio < 0.8) if not np.isnan(ratio) else None})
            print(f'  {a} - {b_}: boot_SE={boot_sd:.5f}, DeLong_SE~={delong_se:.5f}, ratio={ratio:.3f}', flush=True)
    se_df = pd.DataFrame(se_rows)
    se_df.to_csv(os.path.join(OUT_DIR, 'se_check_delong.csv'), index=False)
    n_flagged = int(se_df['flag_ratio_lt_0.8'].sum())
    print(f'\n{n_flagged}/15 pairs flagged with bootstrap/DeLong SE ratio < 0.8', flush=True)

    print(f'\nSAVED -> {OUT_DIR}', flush=True)
    print(f'\n=== DECISION-RULE INPUT: permutation P = {perm_p:.5f} ({"< 0.05, headline stands" if perm_p < 0.05 else ">= 0.05, ESCALATE"}) ===')


if __name__ == '__main__':
    main()
