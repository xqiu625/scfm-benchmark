#!/usr/bin/env python3
"""
Job Q -- paired bootstrap completion (external review, unblocks Figure 3/4/S2
significance layers). Same seed-ensemble patient-level paired bootstrap method
as Job L (B=10,000): average each patient's 10 per-seed scores, one bootstrap
on the ensembled scores. Runs 3 sub-tasks with all-pairs coverage (15 pairs for
6-model sets, 5 pairs for FM-vs-PCA-200-only monocyte).

Usage: python jobQ_paired_inference_full.py {within_covid|within_sepsis|monocyte}
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


def safe_auc(y, p):
    if len(np.unique(y)) < 2:
        return np.nan
    return roc_auc_score(y, p)


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


def load_within_covid():
    m = pd.read_csv(PATIENT_MAP)
    covid_patients = set(m.loc[m['disease'] == 'COVID-19', 'patient_id'])
    per_seed = {}
    for seed in SEEDS:
        per_seed[seed] = {}
        for fm in FMS:
            df = pd.read_csv(os.path.join(OOF_DIR, f'oof_patient_{fm}_seed{seed}.csv'))
            df = df[df['patient_id'].isin(covid_patients)]
            per_seed[seed][fm] = df.groupby('patient_id').agg(y=('y', 'max'), score=('p_lr', 'mean'))
        pca = pd.read_csv(os.path.join(OOF_DIR, f'oof_PCA-200_seed{seed}.csv'))
        pmap = m[['sample_id', 'patient_id']]
        pca = pca.merge(pmap, on='sample_id', how='inner')
        pca = pca[pca['patient_id'].isin(covid_patients)]
        per_seed[seed]['PCA-200'] = pca.groupby('patient_id').agg(y=('y', 'max'), score=('p_lr', 'mean'))
    return per_seed, FMS + ['PCA-200']


def load_within_sepsis():
    per_seed = {}
    for seed in SEEDS:
        per_seed[seed] = {}
        for fm in FMS + ['PCA-200']:
            df = pd.read_csv(os.path.join(OOF_DIR, f'oof_sepsis_{fm}_seed{seed}.csv'))
            if 'patient_id' not in df.columns:
                df['patient_id'] = df['sample_id']  # PCA-200 file lacks it; sample_id==patient_id here
            per_seed[seed][fm] = df.groupby('patient_id').agg(y=('y', 'max'), score=('p_lr', 'mean'))
    return per_seed, FMS + ['PCA-200']


def load_monocyte():
    # primary cut = COVID-only (Ren+Stephenson), matching Figure 4's primary framing
    per_seed = {}
    for seed in SEEDS:
        per_seed[seed] = {}
        for fm in FMS + ['PCA-200']:
            df = pd.read_csv(os.path.join(OOF_DIR, 'monocyte', f'oof_monocyte_{fm}_seed{seed}.csv'))
            df = df[df['source_study'] != 'Qiu_sepsis_GSE167363']
            per_seed[seed][fm] = df.groupby('patient_id').agg(y=('y', 'max'), score=('p_lr', 'mean'))
    return per_seed, FMS + ['PCA-200']


def main():
    mode = sys.argv[1]
    os.makedirs(OUT_DIR, exist_ok=True)
    rng = np.random.default_rng(RNG_SEED)

    if mode == 'within_covid':
        per_seed, models = load_within_covid()
        out_name = 'within_covid_15pairs.csv'
        pairs = [(models[i], models[j]) for i in range(len(models)) for j in range(i + 1, len(models))]
    elif mode == 'within_sepsis':
        per_seed, models = load_within_sepsis()
        out_name = 'within_sepsis_15pairs.csv'
        pairs = [(models[i], models[j]) for i in range(len(models)) for j in range(i + 1, len(models))]
    elif mode == 'monocyte':
        per_seed, models = load_monocyte()
        out_name = 'monocyte_vs_pca_5pairs.csv'
        pairs = [(fm, 'PCA-200') for fm in FMS]
    else:
        raise ValueError(mode)

    print(f'mode={mode}, models={models}, n_pairs={len(pairs)}', flush=True)
    for m in models:
        n_pat = len(per_seed[42][m])
        n_fat = int(per_seed[42][m]['y'].sum())
        print(f'  seed-42 view: {m} n_patients={n_pat} n_fatal={n_fat}', flush=True)

    # seed-ensemble: average each patient's 10 per-seed scores, one bootstrap
    ensembled = {}
    for m in models:
        allp = pd.concat([per_seed[s][m]['score'].rename(s) for s in SEEDS], axis=1)
        y_ = pd.concat([per_seed[s][m]['y'].rename(s) for s in SEEDS], axis=1).max(axis=1)
        ensembled[m] = pd.DataFrame({'y': y_, 'score': allp.mean(axis=1)})

    rows = []
    for a, b_ in pairs:
        common = ensembled[a].join(ensembled[b_], how='inner', lsuffix='_a', rsuffix='_b')
        y = common['y_a'].to_numpy()
        pa = common['score_a'].to_numpy()
        pb = common['score_b'].to_numpy()
        obs, lo, hi, p = paired_bootstrap(y, pa, pb, rng)
        rows.append({'pair': f'{a} - {b_}', 'n_patients': len(common),
                     'delta_auc': round(obs, 4), 'ci_lo': round(lo, 4), 'ci_hi': round(hi, 4),
                     'p_boot': round(p, 5)})
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT_DIR, out_name), index=False)
    print(f'\n=== {mode}: seed-ensemble paired bootstrap ({len(pairs)} pairs) ===')
    print(df.to_string(index=False))
    print(f'\nSAVED -> {os.path.join(OUT_DIR, out_name)}', flush=True)


if __name__ == '__main__':
    main()
