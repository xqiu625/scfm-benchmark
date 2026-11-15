#!/usr/bin/env python3
"""Job H step 4: patient-level re-run of the adaptation curve (stage 29).

Copy of scripts/29_adaptation_curve/code/run_adaptation_curve_bootstrap.py
(the script backing the manuscript's tab:adaptation_curve, point_auc_mean
column) with the sepsis holdout split and the patient-bootstrap grouped by
patient_id instead of Sample_ID (2026-09-02 multi-draw leakage fix; see
reports/manuscript/revision/patient_overlap_finding.md). Protocol otherwise
identical: same grid, same 10 training seeds, same 200-iter bootstrap,
same classifiers/seeds.

Sepsis arm after the fix: 12 patients (2 fatal). The stratified 1/3 holdout
(max(1, round(n*1/3)) per class) yields test = 1 fatal + 3 survived patients,
train = 1 fatal + 7 survived patients -- >=1 fatal patient in both splits.

Output: 05_RESULTS/revision_donor_level/adaptation_curve_bootstrap_patient/
  <FM>.csv, all_rows.csv, summary.csv, split_documentation.txt
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

ROOT = '/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c'
OUT_DIR = ('/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/'
           'revision_donor_level/adaptation_curve_bootstrap_patient')
os.makedirs(OUT_DIR, exist_ok=True)
PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id.csv')

EMBEDDINGS = [
    ('STATE',            'state_covid_canonical.h5ad',            'X_state'),
    ('UCE',              'uce_covid_canonical.h5ad',              'X_uce'),
    ('TranscriptFormer', 'transcriptformer_covid_canonical.h5ad', 'X_transcriptformer'),
    ('scPRINT-2',        'scprint2_covid_canonical.h5ad',         'X_scprint2'),
    ('STACK',            'stack_covid_canonical.h5ad',            'X_stack'),
    ('Geneformer',       'geneformer_covid_canonical.h5ad',       'X_geneformer'),
    ('scGPT',            'scgpt_covid_canonical.h5ad',            'X_scgpt'),
]

KEEP_DISEASE = ['COVID-19', 'Sepsis']
KEEP_SEVERITY = ['fatal', 'mild', 'moderate', 'recovered', 'severe']

N_SEPSIS_GRID = [0, 10, 50, 100, 500, 1000, -1]
N_SEEDS = 10
N_BOOTSTRAP = 200
TEST_FRACTION = 1.0 / 3.0
RANDOM_STATE = 42
PCA_N_COMP = 200

PMAP = pd.read_csv(PATIENT_MAP).set_index('sample_id')['patient_id'].to_dict()


def to_patient_id(sample_id):
    pid = np.array([PMAP.get(s, None) for s in sample_id], dtype=object)
    n_missing = int(pd.isna(pid).sum())
    assert n_missing == 0, f'{n_missing} cells lack patient_id mapping'
    return pid.astype(str)


def read_categorical(group):
    if isinstance(group, h5py.Group) and 'categories' in group:
        cats = group['categories'][:].astype(str)
        codes = group['codes'][:]
        return cats[codes]
    return group[:].astype(str)


def read_obs(h, col):
    if col in h['obs']:
        return read_categorical(h[f'obs/{col}'])
    return None


def load_embedding(path, obsm_key):
    with h5py.File(path, 'r') as h:
        X = h[f'obsm/{obsm_key}'][:]
        disease = read_obs(h, 'Disease')
        severity = read_obs(h, 'covid_severity')
        sample_id = read_obs(h, 'Sample_ID')
    keep = np.isin(disease, KEEP_DISEASE) & np.isin(severity, KEEP_SEVERITY)
    sid = sample_id[keep]
    return X[keep], disease[keep], severity[keep], sid, to_patient_id(sid)


def load_sparse_X(path):
    with h5py.File(path, 'r') as h:
        if 'X' in h and isinstance(h['X'], h5py.Group):
            data = h['X/data'][:]
            indices = h['X/indices'][:]
            indptr = h['X/indptr'][:]
            shape = tuple(h['X'].attrs.get('shape', [len(indptr) - 1, indices.max() + 1]))
            X = sparse.csr_matrix((data, indices, indptr), shape=shape)
        else:
            X = sparse.csr_matrix(h['X'][:])
        disease = read_obs(h, 'Disease')
        severity = read_obs(h, 'covid_severity')
        sample_id = read_obs(h, 'Sample_ID')
    keep = np.isin(disease, KEEP_DISEASE) & np.isin(severity, KEEP_SEVERITY)
    sid = sample_id[keep]
    return X[keep], disease[keep], severity[keep], sid, to_patient_id(sid)


def make_sepsis_holdout(disease, severity, patient_id, seed=RANDOM_STATE):
    """Patient-level: hold out ~TEST_FRACTION of sepsis PATIENTS, stratified
    by fatal/survived (patient label = max(y) over the patient's draws)."""
    is_sepsis = (disease == 'Sepsis')
    is_covid = (disease == 'COVID-19')
    y = (severity == 'fatal').astype(int)

    df = pd.DataFrame({'patient_id': patient_id, 'is_sepsis': is_sepsis, 'y': y})
    pat_y = df[df['is_sepsis']].groupby('patient_id')['y'].max()
    sepsis_fatal_patients = pat_y[pat_y == 1].index.to_numpy()
    sepsis_surv_patients = pat_y[pat_y == 0].index.to_numpy()

    rng = np.random.default_rng(seed)
    rng.shuffle(sepsis_fatal_patients)
    rng.shuffle(sepsis_surv_patients)
    n_test_fatal = max(1, int(round(len(sepsis_fatal_patients) * TEST_FRACTION)))
    n_test_surv = max(1, int(round(len(sepsis_surv_patients) * TEST_FRACTION)))
    test_patients = np.concatenate([sepsis_fatal_patients[:n_test_fatal],
                                     sepsis_surv_patients[:n_test_surv]])

    test_mask = is_sepsis & np.isin(patient_id, test_patients)
    sepsis_train_mask = is_sepsis & ~np.isin(patient_id, test_patients)
    covid_train_mask = is_covid

    train_patients = np.setdiff1d(np.unique(patient_id[is_sepsis]), test_patients)
    return (test_mask, sepsis_train_mask, covid_train_mask, test_patients, y,
            {'sepsis_fatal_patients': sepsis_fatal_patients.tolist(),
             'sepsis_surv_patients': sepsis_surv_patients.tolist(),
             'test_patients': test_patients.tolist(),
             'train_patients': train_patients.tolist(),
             'n_test_fatal': n_test_fatal, 'n_test_surv': n_test_surv})


def patient_bootstrap_auc(y_true, y_pred, patient_ids, n_iter=N_BOOTSTRAP, seed=0):
    rng = np.random.default_rng(seed)
    unique_patients = np.unique(patient_ids)
    patient_to_idx = {p: np.where(patient_ids == p)[0] for p in unique_patients}
    aucs = []
    for _ in range(n_iter):
        sampled = rng.choice(unique_patients, size=len(unique_patients), replace=True)
        cell_idx = np.concatenate([patient_to_idx[p] for p in sampled])
        y_t = y_true[cell_idx]
        y_p = y_pred[cell_idx]
        if len(np.unique(y_t)) < 2:
            continue
        aucs.append(roc_auc_score(y_t, y_p))
    if not aucs:
        return float('nan'), float('nan'), float('nan')
    aucs = np.asarray(aucs)
    return float(np.mean(aucs)), float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))


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
                     patient_id, n_sepsis_grid, n_seeds, transform=None,
                     n_bootstrap=N_BOOTSTRAP):
    sepsis_train_idx = np.where(sepsis_train_mask)[0]
    covid_train_idx = np.where(covid_train_mask)[0]
    test_idx = np.where(test_mask)[0]
    te_patient_ids = patient_id[test_idx]

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
                y_te, p, te_patient_ids, n_iter=n_bootstrap, seed=seed
            )

            rows.append({
                'method': name, 'n_sepsis_req': n_req, 'n_sepsis': int(actual),
                'seed': seed, 'n_covid_train': int(len(covid_train_idx)),
                'n_sepsis_test': int(test_mask.sum()),
                'auc_point': point_auc, 'auc_boot_mean': boot_mean,
                'auc_boot_ci_lo': ci_lo, 'auc_boot_ci_hi': ci_hi,
            })
            print(f'  {name:18s} n={actual:>5} seed={seed} '
                  f'AUC_point={point_auc:.4f}  boot={boot_mean:.4f} [{ci_lo:.3f}, {ci_hi:.3f}]',
                  flush=True)
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
    g = df.groupby(['method', 'n_sepsis_req', 'n_sepsis'])
    return g.agg(
        point_auc_mean=('auc_point', 'mean'),
        point_auc_std=('auc_point', 'std'),
        boot_auc_mean=('auc_boot_mean', 'mean'),
        boot_ci_lo_mean=('auc_boot_ci_lo', 'mean'),
        boot_ci_hi_mean=('auc_boot_ci_hi', 'mean'),
        n_training_seeds=('seed', 'count'),
    ).reset_index()


def main():
    state_path = os.path.join(ROOT, 'state_covid_canonical.h5ad')
    X_state, disease, severity, sample_id, patient_id = load_embedding(state_path, 'X_state')
    test_mask, sepsis_train_mask, covid_train_mask, test_patients, y, split_doc = \
        make_sepsis_holdout(disease, severity, patient_id, seed=RANDOM_STATE)

    print(f'Sepsis test patients:  {len(test_patients)}', flush=True)
    print(f'Sepsis test cells:     {int(test_mask.sum())}', flush=True)
    print(f'Sepsis train pool:     {int(sepsis_train_mask.sum())} cells '
          f'from {len(np.unique(patient_id[sepsis_train_mask]))} patients', flush=True)
    print(f'COVID train pool:      {int(covid_train_mask.sum())} cells', flush=True)
    print(f'Split doc: {split_doc}', flush=True)
    with open(os.path.join(OUT_DIR, 'split_documentation.txt'), 'w') as f:
        f.write(str(split_doc) + '\n')

    all_dfs = []
    for fm, fname, obsm_key in EMBEDDINGS:
        path = os.path.join(ROOT, fname)
        if not os.path.exists(path):
            print(f'SKIP missing input: {path}', flush=True)
            continue
        print(f'\n=== {fm} ===', flush=True)
        Xfm, dfm, sfm, sidfm, pidfm = load_embedding(path, obsm_key)
        tm, stm, ctm, _, yfm, _ = make_sepsis_holdout(dfm, sfm, pidfm, seed=RANDOM_STATE)
        df = adaptation_curve(fm, Xfm, yfm, tm, stm, ctm, pidfm,
                              N_SEPSIS_GRID, N_SEEDS, transform=None,
                              n_bootstrap=N_BOOTSTRAP)
        df.to_csv(os.path.join(OUT_DIR, f'{fm}.csv'), index=False)
        all_dfs.append(df)

    print(f'\n=== PCA-{PCA_N_COMP} + LR (from-scratch baseline) ===', flush=True)
    Xraw, draw, sraw, sidraw, pidraw = load_sparse_X(state_path)
    tm, stm, ctm, _, yraw, _ = make_sepsis_holdout(draw, sraw, pidraw, seed=RANDOM_STATE)
    df_pca = adaptation_curve(f'PCA-{PCA_N_COMP}+LR', Xraw, yraw, tm, stm, ctm, pidraw,
                              N_SEPSIS_GRID, N_SEEDS,
                              transform=pca_transform_factory(PCA_N_COMP),
                              n_bootstrap=N_BOOTSTRAP)
    df_pca.to_csv(os.path.join(OUT_DIR, f'PCA-{PCA_N_COMP}+LR.csv'), index=False)
    all_dfs.append(df_pca)

    pooled = pd.concat(all_dfs, ignore_index=True)
    pooled.to_csv(os.path.join(OUT_DIR, 'all_rows.csv'), index=False)
    summary = summarize(pooled)
    summary.to_csv(os.path.join(OUT_DIR, 'summary.csv'), index=False)
    print(f'\nWrote: all_rows.csv ({len(pooled)} rows), summary.csv ({len(summary)} rows)', flush=True)
    print('\nALL DONE.', flush=True)


if __name__ == '__main__':
    main()
