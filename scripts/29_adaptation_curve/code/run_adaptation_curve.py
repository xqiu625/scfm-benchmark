#!/usr/bin/env python3
"""
Adaptation curve: how much sepsis training data is needed before a model is useful?

Decomposes the cross-disease question into three quantities per FM (+ PCA baseline):
  - Y-intercept = AUC at n_sepsis=0 (pure COVID->Sepsis zero-shot)
  - Early slope = (AUC at n=100) - (AUC at n=0), per added sepsis cell
  - Asymptote   = AUC at n=all (all available sepsis training cells)

Design:
  - Fixed patient-level sepsis test set (~30% of sepsis donors, stratified by fatal/survived)
  - Training pool = all COVID train cells + variable number of sepsis cells drawn from remaining sepsis donors
  - n_sepsis grid: 0, 10, 50, 100, 500, 1000, all
  - 10 bootstrap seeds for sepsis subsampling (test set fixed across seeds)

Outputs:
  adaptation_curve/<FM>.csv          per-seed per-n rows
  adaptation_curve/summary.csv       pooled, with mean and 95% interval

Memory-safe: h5py for FM embeddings; sparse-CSR load of raw X (from STATE h5ad) only for PCA baseline.
"""
import os
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.linear_model import LogisticRegression
from sklearn.decomposition import TruncatedSVD
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

ROOT    = '/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c'
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/option_c_primary_benchmark/covid_sepsis_only/adaptation_curve'
os.makedirs(OUT_DIR, exist_ok=True)

EMBEDDINGS = [
    ('STATE',            'state_covid_canonical.h5ad',            'X_state'),
    ('UCE',              'uce_covid_canonical.h5ad',              'X_uce'),
    ('TranscriptFormer', 'transcriptformer_covid_canonical.h5ad', 'X_transcriptformer'),
    ('scPRINT-2',        'scprint2_covid_canonical.h5ad',         'X_scprint2'),
    ('STACK',            'stack_covid_canonical.h5ad',            'X_stack'),
    ('Geneformer',       'geneformer_covid_canonical.h5ad',       'X_geneformer'),
]

KEEP_DISEASE  = ['COVID-19', 'Sepsis']
KEEP_SEVERITY = ['fatal', 'mild', 'moderate', 'recovered', 'severe']

N_SEPSIS_GRID = [0, 10, 50, 100, 500, 1000, -1]   # -1 = use all available
N_SEEDS       = 10
TEST_FRACTION = 1.0 / 3.0                          # ~33% of sepsis donors held out
RANDOM_STATE  = 42
PCA_N_COMP    = 200                                # PCA baseline dim


# ---------------- IO helpers (mirror rerun_covid_sepsis_all.py) -----------------
def read_categorical(group):
    if isinstance(group, h5py.Group) and 'categories' in group:
        cats  = group['categories'][:].astype(str)
        codes = group['codes'][:]
        return cats[codes]
    return group[:].astype(str)


def read_obs(h, col):
    if col in h['obs']:
        return read_categorical(h[f'obs/{col}'])
    return None


def load_embedding(path, obsm_key):
    with h5py.File(path, 'r') as h:
        X         = h[f'obsm/{obsm_key}'][:]
        disease   = read_obs(h, 'Disease')
        severity  = read_obs(h, 'covid_severity')
        sample_id = read_obs(h, 'Sample_ID')
    keep = np.isin(disease, KEEP_DISEASE) & np.isin(severity, KEEP_SEVERITY)
    return X[keep], disease[keep], severity[keep], sample_id[keep]


def load_sparse_X(path):
    with h5py.File(path, 'r') as h:
        if 'X' in h and isinstance(h['X'], h5py.Group):
            data    = h['X/data'][:]
            indices = h['X/indices'][:]
            indptr  = h['X/indptr'][:]
            shape   = tuple(h['X'].attrs.get('shape', [len(indptr) - 1, indices.max() + 1]))
            X = sparse.csr_matrix((data, indices, indptr), shape=shape)
        else:
            X = sparse.csr_matrix(h['X'][:])
        disease   = read_obs(h, 'Disease')
        severity  = read_obs(h, 'covid_severity')
        sample_id = read_obs(h, 'Sample_ID')
    keep = np.isin(disease, KEEP_DISEASE) & np.isin(severity, KEEP_SEVERITY)
    return X[keep], disease[keep], severity[keep], sample_id[keep]


# ---------------- core experiment -----------------
def make_sepsis_holdout(disease, severity, sample_id, seed=RANDOM_STATE):
    """Patient-level: hold out ~TEST_FRACTION of sepsis donors, stratified by fatal/survived."""
    is_sepsis = (disease == 'Sepsis')
    is_covid  = (disease == 'COVID-19')
    y = (severity == 'fatal').astype(int)

    sepsis_fatal_donors = np.unique(sample_id[is_sepsis & (y == 1)])
    sepsis_surv_donors  = np.unique(sample_id[is_sepsis & (y == 0)])

    rng = np.random.default_rng(seed)
    rng.shuffle(sepsis_fatal_donors)
    rng.shuffle(sepsis_surv_donors)
    n_test_fatal = max(1, int(round(len(sepsis_fatal_donors) * TEST_FRACTION)))
    n_test_surv  = max(1, int(round(len(sepsis_surv_donors)  * TEST_FRACTION)))
    test_donors = np.concatenate([sepsis_fatal_donors[:n_test_fatal],
                                  sepsis_surv_donors[:n_test_surv]])

    test_mask          = is_sepsis & np.isin(sample_id, test_donors)
    sepsis_train_mask  = is_sepsis & ~np.isin(sample_id, test_donors)
    covid_train_mask   = is_covid                                  # all COVID is train pool

    return test_mask, sepsis_train_mask, covid_train_mask, test_donors, y


def fit_eval(X_tr, y_tr, X_te, y_te, standardize=True):
    if standardize:
        sc = StandardScaler(with_mean=not sparse.issparse(X_tr))
        X_tr = sc.fit_transform(X_tr)
        X_te = sc.transform(X_te)
    clf = LogisticRegression(C=1.0, solver='lbfgs', class_weight='balanced',
                             max_iter=2000, n_jobs=1)
    clf.fit(X_tr, y_tr)
    p = clf.predict_proba(X_te)[:, 1]
    return roc_auc_score(y_te, p)


def adaptation_curve(name, X, y, test_mask, sepsis_train_mask, covid_train_mask,
                     n_sepsis_grid, n_seeds, transform=None):
    """transform: optional callable(X_train, X_test) -> (X_tr', X_te') applied per fit (e.g., PCA)."""
    sepsis_train_idx = np.where(sepsis_train_mask)[0]
    covid_train_idx  = np.where(covid_train_mask)[0]
    test_idx         = np.where(test_mask)[0]

    rows = []
    for n_req in n_sepsis_grid:
        actual = len(sepsis_train_idx) if n_req == -1 else min(n_req, len(sepsis_train_idx))
        for seed in range(n_seeds):
            rng = np.random.default_rng(seed)
            sampled = (rng.choice(sepsis_train_idx, size=actual, replace=False)
                       if actual > 0 else np.array([], dtype=int))
            train_idx = np.concatenate([covid_train_idx, sampled]).astype(int)

            X_tr = X[train_idx]
            y_tr = y[train_idx]
            X_te = X[test_idx]
            y_te = y[test_idx]

            if len(np.unique(y_tr)) < 2 or len(np.unique(y_te)) < 2:
                continue

            if transform is not None:
                X_tr, X_te = transform(X_tr, X_te, seed=seed)
                auc = fit_eval(X_tr, y_tr, X_te, y_te, standardize=False)
            else:
                auc = fit_eval(X_tr, y_tr, X_te, y_te, standardize=True)

            rows.append({
                'method':       name,
                'n_sepsis_req': n_req,
                'n_sepsis':     int(actual),
                'seed':         seed,
                'n_covid_train':  int(len(covid_train_idx)),
                'n_sepsis_test':  int(test_mask.sum()),
                'auc_sepsis':   float(auc),
            })
            print(f'  {name:18s} n={actual:>5} seed={seed} AUC={auc:.4f}', flush=True)

            # n=0 is identical across seeds; only run once
            if actual == 0:
                break
    return pd.DataFrame(rows)


def pca_transform_factory(n_comp):
    def transform(X_tr, X_te, seed=0):
        svd = TruncatedSVD(n_components=n_comp, random_state=seed)
        Z_tr = svd.fit_transform(X_tr)
        Z_te = svd.transform(X_te)
        return Z_tr, Z_te
    return transform


def summarize(df):
    g = df.groupby(['method', 'n_sepsis_req', 'n_sepsis'])['auc_sepsis']
    summary = g.agg(['mean', 'std', 'count']).reset_index()
    summary['ci_lo'] = summary['mean'] - 1.96 * summary['std'] / np.sqrt(summary['count'].clip(lower=1))
    summary['ci_hi'] = summary['mean'] + 1.96 * summary['std'] / np.sqrt(summary['count'].clip(lower=1))
    return summary


def main():
    # Establish the holdout against STATE's coordinate system (same filter, same donors across FMs).
    state_path = os.path.join(ROOT, 'state_covid_canonical.h5ad')
    X_state, disease, severity, sample_id = load_embedding(state_path, 'X_state')
    test_mask, sepsis_train_mask, covid_train_mask, test_donors, y = make_sepsis_holdout(
        disease, severity, sample_id, seed=RANDOM_STATE
    )
    print(f'Sepsis test donors:  {len(test_donors)}', flush=True)
    print(f'Sepsis test cells:   {int(test_mask.sum())}', flush=True)
    print(f'Sepsis train pool:   {int(sepsis_train_mask.sum())} cells '
          f'from {len(np.unique(sample_id[sepsis_train_mask]))} donors', flush=True)
    print(f'COVID train pool:    {int(covid_train_mask.sum())} cells', flush=True)
    print(f'Test class balance:  fatal={int((y[test_mask] == 1).sum())}, '
          f'surv={int((y[test_mask] == 0).sum())}', flush=True)

    all_dfs = []

    # ---- FM arms (use the same masks; reload per-FM embeddings; assume same row order via Cell_ID filter) ----
    for fm, fname, obsm_key in EMBEDDINGS:
        path = os.path.join(ROOT, fname)
        if not os.path.exists(path):
            print(f'SKIP missing input: {path}', flush=True)
            continue
        print(f'\n=== {fm} ===', flush=True)
        Xfm, dfm, sfm, idfm = load_embedding(path, obsm_key)
        # Each FM applies its own preprocessing and may drop a few cells; rebuild masks on its own ordering.
        tm, stm, ctm, _, yfm = make_sepsis_holdout(dfm, sfm, idfm, seed=RANDOM_STATE)
        df = adaptation_curve(fm, Xfm, yfm, tm, stm, ctm,
                              N_SEPSIS_GRID, N_SEEDS, transform=None)
        df.to_csv(os.path.join(OUT_DIR, f'{fm}.csv'), index=False)
        all_dfs.append(df)

    # ---- PCA-200 from-scratch baseline on raw counts ----
    print(f'\n=== PCA-{PCA_N_COMP} + LR (from-scratch baseline) ===', flush=True)
    Xraw, draw, sraw, idraw = load_sparse_X(state_path)
    tm, stm, ctm, _, yraw = make_sepsis_holdout(draw, sraw, idraw, seed=RANDOM_STATE)
    df_pca = adaptation_curve(f'PCA-{PCA_N_COMP}+LR', Xraw, yraw, tm, stm, ctm,
                              N_SEPSIS_GRID, N_SEEDS, transform=pca_transform_factory(PCA_N_COMP))
    df_pca.to_csv(os.path.join(OUT_DIR, f'PCA-{PCA_N_COMP}+LR.csv'), index=False)
    all_dfs.append(df_pca)

    # ---- pooled summary ----
    pooled  = pd.concat(all_dfs, ignore_index=True)
    pooled.to_csv(os.path.join(OUT_DIR, 'all_rows.csv'), index=False)
    summary = summarize(pooled)
    summary.to_csv(os.path.join(OUT_DIR, 'summary.csv'), index=False)
    print(f'\nWrote: all_rows.csv ({len(pooled)} rows), summary.csv ({len(summary)} rows)', flush=True)
    print('\nALL DONE.', flush=True)


if __name__ == '__main__':
    main()
