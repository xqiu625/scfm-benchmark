#!/usr/bin/env python3
"""
Regenerate and PERSIST the two raw 10,000-length permutation null vectors that
were computed in-memory by earlier rounds but never saved to disk:
  - the primary test's max-delta-vs-PCA-200 null (jobS2_round4_primary.py)
  - the range-statistic (leaders vs UCE) null (script no longer on disk --
    reconstructed here; P0_onepager.md states it uses "the same batch of
    label permutations" as the primary test, same rng_seed=20260908, which
    this script honors literally: ONE rng stream, ONE y_perm draw per
    iteration, both statistics computed from that same draw)

Deterministic given the same rng_seed and the same frozen OOF scores, so this
is a reproduction, not a new computation. Verified against both already-
reported summary numbers before trusting the saved arrays:
  primary:  observed max delta = 0.0644, permutation P = 0.181
  range:    observed range = 0.1302, permutation P = 0.0112
"""
import os

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

OOF_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2'
OUT_DIR = os.path.join(OOF_DIR, 'paired_inference_full')
PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')

SEEDS = list(range(42, 52))
FMS = ['STATE', 'UCE', 'Geneformer', 'scPRINT-2', 'scGPT']
NPERM = 10000
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


def main():
    covid_patients = load_map()
    ensembled, models = build_ensemble(covid_patients)
    common_all = ensembled['STATE']
    for m in models[1:]:
        common_all = common_all.join(ensembled[m][[]], how='inner')
    idx = common_all.index
    n_patients = len(idx)
    n_fatal = int(ensembled['STATE'].loc[idx, 'y'].sum())
    print(f'Sanity gate: n_patients={n_patients}, n_fatal={n_fatal} '
          f'(expected {N_PATIENTS_EXPECTED}/{N_FATAL_EXPECTED})')
    assert n_patients == N_PATIENTS_EXPECTED and n_fatal == N_FATAL_EXPECTED

    y = ensembled['STATE'].loc[idx, 'y'].to_numpy()
    scores = {m: ensembled[m].loc[idx, 'score'].to_numpy() for m in models}
    fm_scores = {fm: scores[fm] for fm in FMS}
    pca_score = scores['PCA-200']

    # observed statistics
    aucs_obs = {fm: safe_auc(y, fm_scores[fm]) for fm in FMS}
    auc_pca_obs = safe_auc(y, pca_score)
    obs_max_delta = max(aucs_obs[fm] - auc_pca_obs for fm in FMS)
    obs_range = max(aucs_obs.values()) - min(aucs_obs.values())
    print(f'Observed per-FM AUC: {aucs_obs}, PCA-200: {auc_pca_obs:.4f}')
    print(f'Observed max delta (primary) = {obs_max_delta:.4f} (expect 0.0644)')
    print(f'Observed range (leaders-vs-UCE) = {obs_range:.4f} (expect 0.1302)')

    rng = np.random.default_rng(RNG_SEED)
    max_perm = np.empty(NPERM)
    range_perm = np.empty(NPERM)
    for i in range(NPERM):
        y_perm = rng.permutation(y)
        if len(np.unique(y_perm)) < 2:
            max_perm[i] = np.nan
            range_perm[i] = np.nan
            continue
        auc_pca = safe_auc(y_perm, pca_score)
        aucs = {fm: safe_auc(y_perm, fm_scores[fm]) for fm in FMS}
        max_perm[i] = max(aucs[fm] - auc_pca for fm in FMS)
        range_perm[i] = max(aucs.values()) - min(aucs.values())
        if (i + 1) % 1000 == 0:
            print(f'  {i+1}/{NPERM} done', flush=True)

    valid = ~np.isnan(max_perm)
    max_perm_v = max_perm[valid]
    range_perm_v = range_perm[valid]
    p_primary = (1 + np.sum(max_perm_v >= obs_max_delta)) / (1 + len(max_perm_v))
    p_range = (1 + np.sum(range_perm_v >= obs_range)) / (1 + len(range_perm_v))
    print(f'\nRegenerated primary permutation P = {p_primary:.5f} (expect 0.181)')
    print(f'Regenerated range-statistic permutation P = {p_range:.5f} (expect 0.0112)')

    os.makedirs(OUT_DIR, exist_ok=True)
    pd.DataFrame({'perm_index': np.arange(NPERM), 'max_delta_vs_pca': max_perm,
                  'range_leaders_vs_uce': range_perm}).to_csv(
        os.path.join(OUT_DIR, 'primary_and_range_permutation_null_vectors.csv'), index=False)
    print(f'Saved {os.path.join(OUT_DIR, "primary_and_range_permutation_null_vectors.csv")}')


if __name__ == '__main__':
    main()
