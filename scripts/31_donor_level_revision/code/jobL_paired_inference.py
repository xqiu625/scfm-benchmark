#!/usr/bin/env python3
"""
Job L -- definitive patient-level paired inference (external review A3).
Sets the actual headline wording for within-COVID FM-vs-PCA-200 and
STATE-vs-Geneformer comparisons, using all 10 seeds (not the seed-42 pilot).

For each of the 10 seeds' oof_patient_<FM>_seed<NN>.csv / oof_PCA-200_seed<NN>.csv
(05_RESULTS/revision_donor_level_v2/, within-COVID only): aggregate to
patient-level mean score, restrict to COVID-19-disease patients (via the v2
patient map), then:
  (a) per-seed paired patient bootstrap (B=10000) for FM-PCA200 and
      STATE-Geneformer
  (b) seed-ensemble version: average each patient's 10 per-seed scores first,
      then one paired bootstrap on the ensembled scores
  (c) Hanley-McNeil SE for each AUC
  (d) pairwise score correlation between models (Pearson + Spearman on
      per-patient mean scores, seed-42 view for a fixed reference point)
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

OOF_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2'
OUT_DIR = os.path.join(OOF_DIR, 'paired_inference')
PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')

SEEDS = list(range(42, 52))
FMS = ['STATE', 'UCE', 'Geneformer', 'scPRINT-2', 'scGPT']
B = 10000
RNG_SEED = 20260907


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


def hanley_mcneil_se(auc, n_pos, n_neg):
    q1 = auc / (2 - auc)
    q2 = 2 * auc ** 2 / (1 + auc)
    var = (auc * (1 - auc) + (n_pos - 1) * (q1 - auc ** 2) +
           (n_neg - 1) * (q2 - auc ** 2)) / (n_pos * n_neg)
    return float(np.sqrt(var))


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


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    covid_patients = load_map()
    rng = np.random.default_rng(RNG_SEED)

    # ---- load all seeds x all models (FMs + PCA-200) ----
    per_seed = {}  # seed -> model -> pd.DataFrame(patient_id-indexed: y, score)
    for seed in SEEDS:
        per_seed[seed] = {}
        for fm in FMS:
            per_seed[seed][fm] = load_seed_scores(fm, seed, covid_patients)
        per_seed[seed]['PCA-200'] = load_seed_scores('PCA-200', seed, covid_patients, is_pca=True)

    models = FMS + ['PCA-200']
    n_fatal = int(per_seed[42]['STATE']['y'].sum())
    n_surv = int((per_seed[42]['STATE']['y'] == 0).sum())
    print(f'n_fatal={n_fatal}, n_surv={n_surv} (seed-42 view; verify vs. spec assumption '
          f'of 25 fatal/215 survivors)', flush=True)

    # ---- (a) per-seed paired bootstrap for FM-PCA200 and STATE-Geneformer ----
    pairs = [(fm, 'PCA-200') for fm in FMS] + [('STATE', 'Geneformer')]
    rows = []
    for seed in SEEDS:
        for a, b_ in pairs:
            da = per_seed[seed][a]
            db = per_seed[seed][b_]
            common = da.join(db, how='inner', lsuffix='_a', rsuffix='_b')
            y = common['y_a'].to_numpy()
            pa = common['score_a'].to_numpy()
            pb = common['score_b'].to_numpy()
            obs, lo, hi, p = paired_bootstrap(y, pa, pb, rng)
            rows.append({'pair': f'{a} - {b_}', 'seed': seed, 'n_patients': len(common),
                         'delta_auc': round(obs, 4), 'ci_lo': round(lo, 4), 'ci_hi': round(hi, 4),
                         'p_boot': round(p, 5)})
    per_seed_df = pd.DataFrame(rows)
    per_seed_df.to_csv(os.path.join(OUT_DIR, 'per_seed_paired_bootstrap.csv'), index=False)
    print('\n=== (a) per-seed paired bootstrap, summary across 10 seeds ===')
    print(per_seed_df.groupby('pair')['delta_auc'].agg(['mean', 'std']).to_string())
    print('\nfraction of seeds with p<0.05:')
    print(per_seed_df.groupby('pair')['p_boot'].apply(lambda s: (s < 0.05).mean()).to_string())

    # ---- (b) seed-ensemble: average each patient's 10 per-seed scores, one bootstrap ----
    ens_rows = []
    ensembled = {}
    for m in models:
        allp = pd.concat([per_seed[s][m]['score'].rename(s) for s in SEEDS], axis=1)
        y_ = pd.concat([per_seed[s][m]['y'].rename(s) for s in SEEDS], axis=1).max(axis=1)
        ens = pd.DataFrame({'y': y_, 'score': allp.mean(axis=1)})
        ensembled[m] = ens
    for a, b_ in pairs:
        common = ensembled[a].join(ensembled[b_], how='inner', lsuffix='_a', rsuffix='_b')
        y = common['y_a'].to_numpy()
        pa = common['score_a'].to_numpy()
        pb = common['score_b'].to_numpy()
        obs, lo, hi, p = paired_bootstrap(y, pa, pb, rng)
        ens_rows.append({'pair': f'{a} - {b_}', 'n_patients': len(common),
                         'delta_auc': round(obs, 4), 'ci_lo': round(lo, 4), 'ci_hi': round(hi, 4),
                         'p_boot': round(p, 5)})
    ens_df = pd.DataFrame(ens_rows)
    ens_df.to_csv(os.path.join(OUT_DIR, 'seed_ensemble_paired_bootstrap.csv'), index=False)
    print('\n=== (b) seed-ensemble paired bootstrap (headline candidate) ===')
    print(ens_df.to_string(index=False))

    # ---- (c) Hanley-McNeil SE per model (seed-ensemble AUC) ----
    hm_rows = []
    for m in models:
        y = ensembled[m]['y'].to_numpy()
        p = ensembled[m]['score'].to_numpy()
        auc = safe_auc(y, p)
        n_pos = int(y.sum())
        n_neg = int((y == 0).sum())
        se = hanley_mcneil_se(auc, n_pos, n_neg)
        hm_rows.append({'model': m, 'auc_ensemble': round(auc, 4), 'n_fatal': n_pos,
                         'n_survived': n_neg, 'hanley_mcneil_se': round(se, 4)})
    hm_df = pd.DataFrame(hm_rows)
    hm_df.to_csv(os.path.join(OUT_DIR, 'hanley_mcneil_se.csv'), index=False)
    print('\n=== (c) Hanley-McNeil SE, seed-ensemble AUC ===')
    print(hm_df.to_string(index=False))

    # ---- (d) pairwise score correlation (seed-ensemble scores) ----
    score_mat = pd.concat([ensembled[m]['score'].rename(m) for m in models], axis=1).dropna()
    pearson = score_mat.corr(method='pearson')
    spearman = score_mat.corr(method='spearman')
    pearson.to_csv(os.path.join(OUT_DIR, 'score_correlation_pearson.csv'))
    spearman.to_csv(os.path.join(OUT_DIR, 'score_correlation_spearman.csv'))
    print('\n=== (d) pairwise score correlation (Pearson), seed-ensemble ===')
    print(pearson.round(3).to_string())

    print(f'\nSAVED -> {OUT_DIR}', flush=True)


if __name__ == '__main__':
    main()
