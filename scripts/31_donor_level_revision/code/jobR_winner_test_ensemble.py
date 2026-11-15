#!/usr/bin/env python3
"""
Job R -- seed-ensemble winner-corrected max-statistic test + fragility +
consensus predictor (external round-3 review decision item).

The primary within-COVID FM-vs-PCA-200 contrast must be a winner-corrected
max-statistic test (does ANY of the 5 FMs beat PCA-200, corrected for testing
5 models at once) rather than 5 separate per-FM tests. This also runs a
leave-one-fatal-out fragility check on the two headline pairs (Geneformer-PCA,
STATE-PCA) and a consensus-predictor (5-FM mean) sanity check.

Reuses jobL_paired_inference.py's exact seed-ensemble construction and
patient-level paired-bootstrap machinery (same OOF_DIR/PATIENT_MAP/seeds).

Hard sanity gate: n_patients must be 240 and n_fatal must be 25 on the
COVID-19-only filtered cohort, checked BEFORE any computation. If either
check fails, the script exits with an error and writes no output files.

Outputs (OUT_DIR = .../paired_inference_full/):
  winner_test_ensemble.csv    per-FM delta_auc vs PCA-200 (point estimate) +
                               one familywise-summary row (max-statistic P)
  fragility_ensemble.csv      25 leave-one-fatal-out refits x 2 pairs, plus a
                               summary row per pair (p range, n p>0.05)
  consensus_ensemble.csv      consensus (5-FM mean) AUC vs PCA-200
  score_correlation_ensemble.csv   6x6 Pearson corr (5 FM + PCA-200)
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

SEEDS = list(range(42, 52))
FMS = ['STATE', 'UCE', 'Geneformer', 'scPRINT-2', 'scGPT']
B = 10000
RNG_SEED = 20260908
N_PATIENTS_EXPECTED = 240
N_FATAL_EXPECTED = 25


def load_map():
    m = pd.read_csv(PATIENT_MAP)
    covid_patients = set(m.loc[m['disease'] == 'COVID-19', 'patient_id'])
    return covid_patients


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
        ens = pd.DataFrame({'y': y_, 'score': allp.mean(axis=1)})
        ensembled[m] = ens
    return ensembled, models


def paired_bootstrap(y, pa, pb, rng, b=B):
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
    lo, hi = np.percentile(deltas, [2.5, 97.5])
    p_two = 2 * min((deltas <= 0).mean(), (deltas >= 0).mean())
    p_two = min(p_two, 1.0)
    return obs, lo, hi, p_two


def max_statistic_familywise(y, fm_scores, pca_score, rng, b=B):
    """For each resample, compute max_fm(AUC_fm - AUC_pca); familywise P =
    fraction of resamples where that max is <= 0 (no FM beats PCA at all)."""
    n = len(y)
    idx_all = np.arange(n)
    max_deltas = np.empty(b)
    per_fm_deltas = {fm: np.empty(b) for fm in fm_scores}
    for i in range(b):
        idx = rng.choice(idx_all, size=n, replace=True)
        yi = y[idx]
        if len(np.unique(yi)) < 2:
            max_deltas[i] = np.nan
            for fm in fm_scores:
                per_fm_deltas[fm][i] = np.nan
            continue
        auc_pca = safe_auc(yi, pca_score[idx])
        ds = []
        for fm, sc in fm_scores.items():
            d = safe_auc(yi, sc[idx]) - auc_pca
            per_fm_deltas[fm][i] = d
            ds.append(d)
        max_deltas[i] = max(ds)
    valid = ~np.isnan(max_deltas)
    max_deltas = max_deltas[valid]
    familywise_p = float((max_deltas <= 0).mean())
    obs_max = max(safe_auc(y, fm_scores[fm]) - safe_auc(y, pca_score) for fm in fm_scores)
    lo, hi = np.percentile(max_deltas, [2.5, 97.5])
    return familywise_p, obs_max, lo, hi, per_fm_deltas


def main():
    covid_patients = load_map()
    ensembled, models = build_ensemble(covid_patients)

    # ---- hard sanity gate ----
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

    # align all models to the common patient set
    idx = common_all.index
    y = ensembled['STATE'].loc[idx, 'y'].to_numpy()
    scores = {m: ensembled[m].loc[idx, 'score'].to_numpy() for m in models}
    fm_scores = {fm: scores[fm] for fm in FMS}
    pca_score = scores['PCA-200']

    # ============================================================
    # 1. max-statistic familywise test
    # ============================================================
    familywise_p, obs_max, lo, hi, per_fm_deltas = max_statistic_familywise(
        y, fm_scores, pca_score, rng, b=B)
    print(f'\n=== (1) max-statistic familywise test ===')
    print(f'familywise P (no FM beats PCA) = {familywise_p:.5f}, obs max delta = '
          f'{obs_max:.4f} [{lo:.4f}, {hi:.4f}]', flush=True)

    wt_rows = []
    for fm in FMS:
        auc_fm = safe_auc(y, fm_scores[fm])
        auc_pca = safe_auc(y, pca_score)
        d = per_fm_deltas[fm]
        d = d[~np.isnan(d)]
        p_indiv = 2 * min((d <= 0).mean(), (d >= 0).mean())
        p_indiv = min(p_indiv, 1.0)
        dlo, dhi = np.percentile(d, [2.5, 97.5])
        wt_rows.append({
            'row_type': 'per_model', 'model': fm,
            'auc_model': round(auc_fm, 4), 'auc_pca': round(auc_pca, 4),
            'delta_auc': round(auc_fm - auc_pca, 4),
            'ci_lo': round(dlo, 4), 'ci_hi': round(dhi, 4),
            'p_indiv_uncorrected': round(p_indiv, 5),
            'B': B, 'rng_seed': RNG_SEED,
        })
    wt_rows.append({
        'row_type': 'familywise_summary', 'model': 'max_over_5_FMs',
        'auc_model': np.nan, 'auc_pca': round(safe_auc(y, pca_score), 4),
        'delta_auc': round(obs_max, 4),
        'ci_lo': round(lo, 4), 'ci_hi': round(hi, 4),
        'p_indiv_uncorrected': round(familywise_p, 5),
        'B': B, 'rng_seed': RNG_SEED,
    })
    wt_df = pd.DataFrame(wt_rows)
    wt_df.to_csv(os.path.join(OUT_DIR, 'winner_test_ensemble.csv'), index=False)
    print(wt_df.to_string(index=False))

    # ============================================================
    # 2. leave-one-fatal-out fragility (Geneformer-PCA, STATE-PCA)
    # ============================================================
    print(f'\n=== (2) leave-one-fatal-out fragility (25 refits) ===', flush=True)
    fatal_patients = list(idx[y == 1])
    assert len(fatal_patients) == N_FATAL_EXPECTED, \
        f'expected {N_FATAL_EXPECTED} fatal patients, got {len(fatal_patients)}'

    frag_pairs = [('Geneformer', 'PCA-200'), ('STATE', 'PCA-200')]
    frag_rows = []
    for a, b_ in frag_pairs:
        p_list = []
        for drop_i, drop_patient in enumerate(fatal_patients):
            keep_mask = idx != drop_patient
            y_sub = y[keep_mask]
            pa_sub = scores[a][keep_mask]
            pb_sub = scores[b_][keep_mask]
            obs, dlo, dhi, p = paired_bootstrap(y_sub, pa_sub, pb_sub, rng, b=B)
            p_list.append(p)
            frag_rows.append({
                'row_type': 'refit', 'pair': f'{a} - {b_}', 'refit_idx': drop_i,
                'dropped_patient': drop_patient, 'n_patients': int(keep_mask.sum()),
                'delta_auc': round(obs, 4), 'ci_lo': round(dlo, 4), 'ci_hi': round(dhi, 4),
                'p_boot': round(p, 5),
            })
        p_arr = np.array(p_list)
        frag_rows.append({
            'row_type': 'summary', 'pair': f'{a} - {b_}', 'refit_idx': -1,
            'dropped_patient': 'ALL_25', 'n_patients': N_PATIENTS_EXPECTED - 1,
            'delta_auc': np.nan, 'ci_lo': round(p_arr.min(), 5), 'ci_hi': round(p_arr.max(), 5),
            'p_boot': int((p_arr > 0.05).sum()),
        })
        print(f'{a} - {b_}: p range [{p_arr.min():.4f}, {p_arr.max():.4f}], '
              f'{int((p_arr > 0.05).sum())}/25 refits with p>0.05', flush=True)
    frag_df = pd.DataFrame(frag_rows)
    frag_df.to_csv(os.path.join(OUT_DIR, 'fragility_ensemble.csv'), index=False)

    # ============================================================
    # 3. consensus predictor + correlation matrix
    # ============================================================
    print(f'\n=== (3) consensus predictor + correlation matrix ===', flush=True)
    consensus_score = np.mean([scores[fm] for fm in FMS], axis=0)
    auc_consensus = safe_auc(y, consensus_score)
    auc_pca = safe_auc(y, pca_score)
    best_fm = max(FMS, key=lambda fm: safe_auc(y, scores[fm]))
    auc_best_fm = safe_auc(y, scores[best_fm])
    cons_df = pd.DataFrame([{
        'auc_consensus_5fm_mean': round(auc_consensus, 4),
        'auc_pca200': round(auc_pca, 4),
        'delta_consensus_minus_pca': round(auc_consensus - auc_pca, 4),
        'best_single_fm': best_fm,
        'auc_best_single_fm': round(auc_best_fm, 4),
        'delta_consensus_minus_best_fm': round(auc_consensus - auc_best_fm, 4),
        'n_patients': N_PATIENTS_EXPECTED, 'n_fatal': N_FATAL_EXPECTED,
    }])
    cons_df.to_csv(os.path.join(OUT_DIR, 'consensus_ensemble.csv'), index=False)
    print(cons_df.to_string(index=False))

    score_mat = pd.DataFrame({m: scores[m] for m in models}, index=idx)
    corr = score_mat.corr(method='pearson')
    corr.to_csv(os.path.join(OUT_DIR, 'score_correlation_ensemble.csv'))
    print('\nPearson correlation (6x6):')
    print(corr.round(3).to_string())

    print(f'\nSAVED -> {OUT_DIR}', flush=True)


if __name__ == '__main__':
    main()
