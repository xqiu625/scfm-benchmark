#!/usr/bin/env python3
"""
Ranking seed-stability (roadmap P2 / Devil's-Advocate alt-explanation #1):
"With per-fold SD 0.09-0.14 on 322 donors, a different CV seed than 42 would
plausibly reorder scPRINT-2/Geneformer/STATE entirely."

Re-runs the within-disease COVID+Sepsis mortality LR benchmark under N CV seeds
(loads each FM embedding once, loops seeds) and reports whether the per-FM ranking
is stable. Memory-safe (h5py-only, never loads X).

Output: 05_RESULTS/option_c_primary_benchmark/seed_stability.csv  (fm x seed x mean AUC)
        + printed ranking-stability summary (rank-1 frequency, Kendall tau).
"""
import h5py
import numpy as np
import pandas as pd
import os
from itertools import combinations
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score

ROOT = '/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c'
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/option_c_primary_benchmark'

EMBEDDINGS = [
    ('STATE',      'state_covid_canonical.h5ad',    'X_state'),
    ('scPRINT-2',  'scprint2_covid_canonical.h5ad', 'X_scprint2'),
    ('Geneformer', 'geneformer_covid_canonical.h5ad', 'X_geneformer'),
    ('scGPT',      'scgpt_covid_canonical.h5ad',    'X_scgpt'),
    ('UCE',        'uce_covid_canonical.h5ad',      'X_uce'),
]
KEEP_DISEASE = ['COVID-19', 'Sepsis']
KEEP_SEVERITY = ['fatal', 'mild', 'moderate', 'recovered', 'severe']
SEEDS = [42, 0, 1, 7, 13, 123, 2024, 31337]


def read_categorical(group):
    if isinstance(group, h5py.Group) and 'categories' in group:
        cats = group['categories'][:].astype(str)
        codes = group['codes'][:]
        return cats[codes]
    return group[:].astype(str)


def load_fm(fname, obsm_key):
    path = os.path.join(ROOT, fname)
    with h5py.File(path, 'r') as h:
        X = h[f'obsm/{obsm_key}'][:]
        disease = read_categorical(h['obs/Disease'])
        severity = read_categorical(h['obs/covid_severity'])
        sid = read_categorical(h['obs/Sample_ID'])
    keep = np.isin(disease, KEEP_DISEASE) & np.isin(severity, KEEP_SEVERITY)
    y = (severity[keep] == 'fatal').astype(int)
    return X[keep], y, sid[keep]


def mean_auc(X, y, sid, seed):
    skf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=seed)
    aucs = []
    for tr, te in skf.split(X, y, groups=sid):
        clf = LogisticRegression(C=1.0, solver='lbfgs',
                                 class_weight='balanced', max_iter=1000)
        clf.fit(X[tr], y[tr])
        p = clf.predict_proba(X[te])[:, 1]
        aucs.append(roc_auc_score(y[te], p))
    return float(np.mean(aucs))


def kendall_tau(r1, r2):
    keys = list(r1)
    n = len(keys)
    c = d = 0
    for i, j in combinations(range(n), 2):
        a = r1[keys[i]] - r1[keys[j]]
        b = r2[keys[i]] - r2[keys[j]]
        if a * b > 0:
            c += 1
        elif a * b < 0:
            d += 1
    return (c - d) / (c + d) if (c + d) else 1.0


def main():
    print(f'Seed-stability: {len(SEEDS)} seeds x {len(EMBEDDINGS)} FMs', flush=True)
    rows = []
    fm_data = {}
    for fm, fname, key in EMBEDDINGS:
        if not os.path.exists(os.path.join(ROOT, fname)):
            print(f'SKIP {fm}', flush=True)
            continue
        fm_data[fm] = load_fm(fname, key)
        print(f'loaded {fm}: {fm_data[fm][0].shape}', flush=True)

    per_seed_auc = {}  # seed -> {fm: auc}
    for seed in SEEDS:
        per_seed_auc[seed] = {}
        for fm in fm_data:
            X, y, sid = fm_data[fm]
            a = mean_auc(X, y, sid, seed)
            per_seed_auc[seed][fm] = a
            rows.append({'fm': fm, 'seed': seed, 'mean_auc': round(a, 4)})
        ranked = sorted(per_seed_auc[seed], key=per_seed_auc[seed].get, reverse=True)
        print(f'seed {seed}: ' + ' > '.join(f'{m}({per_seed_auc[seed][m]:.3f})' for m in ranked),
              flush=True)

    df = pd.DataFrame(rows)
    os.makedirs(OUT_DIR, exist_ok=True)
    df.to_csv(os.path.join(OUT_DIR, 'seed_stability.csv'), index=False)

    print('\n=== per-FM AUC across seeds ===', flush=True)
    summ = df.groupby('fm')['mean_auc'].agg(['mean', 'std', 'min', 'max'])
    print(summ.round(4).to_string(), flush=True)

    # rank-1 frequency
    print('\n=== rank-1 frequency ===', flush=True)
    rank1 = {}
    ranks = {}  # seed -> {fm: rank}
    for seed in SEEDS:
        ordered = sorted(per_seed_auc[seed], key=per_seed_auc[seed].get, reverse=True)
        rank1[ordered[0]] = rank1.get(ordered[0], 0) + 1
        ranks[seed] = {fm: i for i, fm in enumerate(ordered)}
    for fm, c in sorted(rank1.items(), key=lambda x: -x[1]):
        print(f'  {fm}: rank-1 in {c}/{len(SEEDS)} seeds', flush=True)

    # pairwise Kendall tau between seed rankings
    taus = [kendall_tau(ranks[s1], ranks[s2]) for s1, s2 in combinations(SEEDS, 2)]
    print(f'\n=== ranking stability: mean Kendall tau across seed-pairs = {np.mean(taus):.3f} '
          f'(1.0 = identical order) ===', flush=True)
    print(f'SAVED: {OUT_DIR}/seed_stability.csv', flush=True)


if __name__ == '__main__':
    main()
