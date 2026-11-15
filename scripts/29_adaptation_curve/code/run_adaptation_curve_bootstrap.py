#!/usr/bin/env python3
"""Adaptation curve, bootstrap-corrected version.

Same training protocol as `run_adaptation_curve.py` (10 random training-
augmentation seeds at each n_sepsis), but adds a 200-iteration patient-
level bootstrap on the test predictions after each fit -- matching the
protocol used by `scripts/06_batch_generalization/code/cross_disease_fatal_vs_survived.py`.

Motivation (2026-05-22): the legacy 10-seed fixed-holdout protocol
reports a single point AUC per fit, with no uncertainty estimate. This
bootstrap-corrected version adds a patient-level bootstrap CI on the test
predictions after each fit, so both the point estimate and its spread are
available per (method, n_sepsis). NOTE: the manuscript adaptation-curve
table (tab:adaptation_curve, claim C005 in docs/claim-traceability.md)
uses the `point_auc_mean` column of this script's summary.csv, NOT
`boot_auc_mean`; the bootstrap columns provide the CIs.

Schema (per row in <FM>.csv):
  method, n_sepsis_req, n_sepsis, seed, n_covid_train, n_sepsis_test,
  auc_point, auc_boot_mean, auc_boot_ci_lo, auc_boot_ci_hi

Outputs:
  <OUT_DIR>/<FM>.csv         per-fit rows (10 training seeds * 7 n grid)
  <OUT_DIR>/all_rows.csv     pooled
  <OUT_DIR>/summary.csv      aggregate stats per (method, n)
"""
import argparse
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
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/option_c_primary_benchmark/covid_sepsis_only/adaptation_curve_bootstrap'
os.makedirs(OUT_DIR, exist_ok=True)

EMBEDDINGS = [
    ('STATE',            'state_covid_canonical.h5ad',            'X_state'),
    ('UCE',              'uce_covid_canonical.h5ad',              'X_uce'),
    ('TranscriptFormer', 'transcriptformer_covid_canonical.h5ad', 'X_transcriptformer'),
    ('scPRINT-2',        'scprint2_covid_canonical.h5ad',         'X_scprint2'),
    ('STACK',            'stack_covid_canonical.h5ad',            'X_stack'),
    ('Geneformer',       'geneformer_covid_canonical.h5ad',       'X_geneformer'),
    # scGPT Option C canonical embedding landed 2026-06-08, after the original
    # bootstrap run; added so `--models scGPT` (or a full rerun) fills its row.
    ('scGPT',            'scgpt_covid_canonical.h5ad',            'X_scgpt'),
]

KEEP_DISEASE  = ['COVID-19', 'Sepsis']
KEEP_SEVERITY = ['fatal', 'mild', 'moderate', 'recovered', 'severe']

N_SEPSIS_GRID = [0, 10, 50, 100, 500, 1000, -1]   # -1 = use all available
N_SEEDS       = 10
N_BOOTSTRAP   = 200                                # patient-level test bootstrap iterations
TEST_FRACTION = 1.0 / 3.0
RANDOM_STATE  = 42
PCA_N_COMP    = 200


# ---------------- IO helpers (verbatim from run_adaptation_curve.py) ----------------
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


# ---------------- holdout + bootstrap ----------------
def make_sepsis_holdout(disease, severity, sample_id, seed=RANDOM_STATE):
    """Patient-level: hold out ~TEST_FRACTION of sepsis donors, stratified by fatal/survived.

    Identical to the legacy protocol so curve shape is comparable across runs;
    bootstrap is layered on top of this fixed holdout.
    """
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
    covid_train_mask   = is_covid

    return test_mask, sepsis_train_mask, covid_train_mask, test_donors, y


def patient_bootstrap_auc(y_true, y_pred, sample_ids, n_iter=N_BOOTSTRAP, seed=0):
    """Patient-level bootstrap of test AUC.

    Resamples DONORS with replacement; for each draw, collects all cells
    from those donors and computes AUC over the resampled cohort.
    Returns (mean, 2.5%ile, 97.5%ile) across n_iter resamples.
    """
    rng = np.random.default_rng(seed)
    unique_donors = np.unique(sample_ids)

    # Pre-compute per-donor cell index lists for fast lookup
    donor_to_idx = {d: np.where(sample_ids == d)[0] for d in unique_donors}

    aucs = []
    for _ in range(n_iter):
        sampled_donors = rng.choice(unique_donors, size=len(unique_donors), replace=True)
        cell_idx = np.concatenate([donor_to_idx[d] for d in sampled_donors])
        y_t = y_true[cell_idx]
        y_p = y_pred[cell_idx]
        if len(np.unique(y_t)) < 2:
            continue
        aucs.append(roc_auc_score(y_t, y_p))

    if not aucs:
        return float('nan'), float('nan'), float('nan')
    aucs = np.asarray(aucs)
    return float(np.mean(aucs)), float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))


# ---------------- core experiment ----------------
def fit_predict(X_tr, y_tr, X_te, standardize=True):
    if standardize:
        sc = StandardScaler(with_mean=not sparse.issparse(X_tr))
        X_tr = sc.fit_transform(X_tr)
        X_te = sc.transform(X_te)
    clf = LogisticRegression(C=1.0, solver='lbfgs', class_weight='balanced',
                             max_iter=2000, n_jobs=1)
    clf.fit(X_tr, y_tr)
    return clf.predict_proba(X_te)[:, 1]


def adaptation_curve(name, X, y, test_mask, sepsis_train_mask, covid_train_mask,
                     sample_id, n_sepsis_grid, n_seeds, transform=None,
                     n_bootstrap=N_BOOTSTRAP):
    sepsis_train_idx = np.where(sepsis_train_mask)[0]
    covid_train_idx  = np.where(covid_train_mask)[0]
    test_idx         = np.where(test_mask)[0]
    te_sample_ids    = sample_id[test_idx]

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
                p = fit_predict(X_tr, y_tr, X_te, standardize=False)
            else:
                p = fit_predict(X_tr, y_tr, X_te, standardize=True)

            point_auc = float(roc_auc_score(y_te, p))
            boot_mean, ci_lo, ci_hi = patient_bootstrap_auc(
                y_te, p, te_sample_ids, n_iter=n_bootstrap, seed=seed
            )

            rows.append({
                'method':         name,
                'n_sepsis_req':   n_req,
                'n_sepsis':       int(actual),
                'seed':           seed,
                'n_covid_train':  int(len(covid_train_idx)),
                'n_sepsis_test':  int(test_mask.sum()),
                'auc_point':      point_auc,
                'auc_boot_mean':  boot_mean,
                'auc_boot_ci_lo': ci_lo,
                'auc_boot_ci_hi': ci_hi,
            })
            print(
                f'  {name:18s} n={actual:>5} seed={seed} '
                f'AUC_point={point_auc:.4f}  boot={boot_mean:.4f} [{ci_lo:.3f}, {ci_hi:.3f}]',
                flush=True,
            )

            if actual == 0:
                # n=0 training is deterministic; one fit per seed gives identical predictions.
                # Bootstrap below is what gives the CI; no need to repeat training seeds.
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
    g = df.groupby(['method', 'n_sepsis_req', 'n_sepsis'])
    summary = g.agg(
        point_auc_mean   =('auc_point',     'mean'),
        point_auc_std    =('auc_point',     'std'),
        boot_auc_mean    =('auc_boot_mean', 'mean'),
        boot_ci_lo_mean  =('auc_boot_ci_lo','mean'),
        boot_ci_hi_mean  =('auc_boot_ci_hi','mean'),
        n_training_seeds =('seed',          'count'),
    ).reset_index()
    return summary


def main():
    parser = argparse.ArgumentParser(description='Adaptation curve, bootstrap-corrected')
    parser.add_argument('--models', nargs='+', default=None,
                        help='Subset of methods to (re)compute, by display name '
                             '(e.g. --models scGPT). Default: all FMs + the '
                             'PCA-200+LR baseline (original behavior). When a '
                             'subset is given, the PCA baseline is skipped unless '
                             'PCA-200+LR is listed, and the new rows are merged '
                             'into the existing all_rows.csv / summary.csv instead '
                             'of overwriting them.')
    args = parser.parse_args()

    if args.models is None:
        selected = EMBEDDINGS
        run_pca_baseline = True
    else:
        selected = [e for e in EMBEDDINGS if e[0] in args.models]
        run_pca_baseline = f'PCA-{PCA_N_COMP}+LR' in args.models
        unknown = set(args.models) - {e[0] for e in EMBEDDINGS} - {f'PCA-{PCA_N_COMP}+LR'}
        if unknown:
            print(f'WARNING: unknown model names ignored: {sorted(unknown)}', flush=True)

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
    print(f'Bootstrap iters:     {N_BOOTSTRAP} (patient-level, with replacement)', flush=True)

    all_dfs = []

    for fm, fname, obsm_key in selected:
        path = os.path.join(ROOT, fname)
        if not os.path.exists(path):
            print(f'SKIP missing input: {path}', flush=True)
            continue
        print(f'\n=== {fm} ===', flush=True)
        Xfm, dfm, sfm, idfm = load_embedding(path, obsm_key)
        tm, stm, ctm, _, yfm = make_sepsis_holdout(dfm, sfm, idfm, seed=RANDOM_STATE)
        df = adaptation_curve(fm, Xfm, yfm, tm, stm, ctm, idfm,
                              N_SEPSIS_GRID, N_SEEDS, transform=None,
                              n_bootstrap=N_BOOTSTRAP)
        df.to_csv(os.path.join(OUT_DIR, f'{fm}.csv'), index=False)
        all_dfs.append(df)

    if run_pca_baseline:
        print(f'\n=== PCA-{PCA_N_COMP} + LR (from-scratch baseline) ===', flush=True)
        Xraw, draw, sraw, idraw = load_sparse_X(state_path)
        tm, stm, ctm, _, yraw = make_sepsis_holdout(draw, sraw, idraw, seed=RANDOM_STATE)
        df_pca = adaptation_curve(f'PCA-{PCA_N_COMP}+LR', Xraw, yraw, tm, stm, ctm, idraw,
                                  N_SEPSIS_GRID, N_SEEDS,
                                  transform=pca_transform_factory(PCA_N_COMP),
                                  n_bootstrap=N_BOOTSTRAP)
        df_pca.to_csv(os.path.join(OUT_DIR, f'PCA-{PCA_N_COMP}+LR.csv'), index=False)
        all_dfs.append(df_pca)

    pooled = pd.concat(all_dfs, ignore_index=True)

    # Subset runs merge into the existing pooled results (rows for methods
    # recomputed here are replaced; every other method is preserved) so a
    # single-model top-up does not clobber the other FMs' rows.
    all_rows_path = os.path.join(OUT_DIR, 'all_rows.csv')
    if args.models is not None and os.path.exists(all_rows_path):
        existing = pd.read_csv(all_rows_path)
        recomputed = set(pooled['method'].unique())
        kept = existing[~existing['method'].isin(recomputed)]
        pooled = pd.concat([kept, pooled], ignore_index=True)
        print(f'\nMerged into existing all_rows.csv: replaced methods {sorted(recomputed)}, '
              f'kept {kept["method"].nunique()} other methods', flush=True)
    elif args.models is not None:
        print('\nWARNING: no existing all_rows.csv found; pooled outputs will contain '
              'only the methods computed in this run.', flush=True)

    pooled.to_csv(all_rows_path, index=False)
    summary = summarize(pooled)
    summary.to_csv(os.path.join(OUT_DIR, 'summary.csv'), index=False)
    print(f'\nWrote: all_rows.csv ({len(pooled)} rows), summary.csv ({len(summary)} rows)', flush=True)
    print('\nALL DONE.', flush=True)


if __name__ == '__main__':
    main()
