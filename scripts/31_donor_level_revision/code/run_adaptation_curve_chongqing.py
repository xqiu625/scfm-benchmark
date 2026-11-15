#!/usr/bin/env python3
"""Job S step 5 (adaptation curve, Chongqing-expanded sepsis pool).

Copy of run_adaptation_curve_bootstrap_patient.py (Job H2, legacy 12-patient
sepsis arm) with two changes per the 2026-09-04 code review + coordinator
decision:
  (a) PATIENT_MAP switched from v1 to donor_meta_with_patient_id_v2.csv (the
      "adaptation curve unchanged" lock was about the 10-seed training-pool-
      subsampling DESIGN, not about freezing a since-corrected mapping --
      confirmed explicitly). Old v1-mapped adaptation curve numbers are
      SUPERSEDED by this run, not adjudicated.
  (b) the sepsis arm is expanded from 12 (Qiu+Jiang) to 103 patients (12 + 91
      Chongqing) by concatenating each FM's Chongqing embedding onto the
      platelet-canonical embedding before calling make_sepsis_holdout() --
      that function already operates generically on disease=='Sepsis', so no
      change to the holdout/subsampling logic itself is needed.

Chongqing cells: Sample_ID IS the patient identifier (no PMAP lookup needed,
unlike the platelet cohort's multi-draw patients). Only STATE/UCE/Geneformer/
scPRINT-2/scGPT have Chongqing embeddings -- TranscriptFormer/STACK do not,
so their sepsis pool stays at the legacy 12 (documented per-method, not
silently different without a note).

PCA-200+LR baseline: platelet raw counts use gene SYMBOLS as var_names,
Chongqing raw counts use Ensembl IDs (with a gene_symbol column) -- reindexed
to the intersection of symbols before concatenation.

Output: 05_RESULTS/revision_donor_level_v2/adaptation_curve_bootstrap_chongqing/
  <FM>.csv, all_rows.csv, summary.csv, split_documentation.txt
"""
import os

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

ROOT = '/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c'
CHONGQING_ROOT = '/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c_chongqing'
OUT_DIR = ('/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2/'
           'adaptation_curve_bootstrap_chongqing')
os.makedirs(OUT_DIR, exist_ok=True)
PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')
CHONGQING_RAW = ('/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/'
                 'chongqing_sepsis_platelet_raw_counts.h5ad')

EMBEDDINGS = [
    ('STATE',      'state_covid_canonical.h5ad',      'X_state',      'state_chongqing_canonical.h5ad', 'X_state'),
    ('UCE',        'uce_covid_canonical.h5ad',        'X_uce',        'uce_chongqing_canonical.h5ad', 'X_uce'),
    ('scPRINT-2',  'scprint2_covid_canonical.h5ad',   'X_scprint2',   'scprint2_chongqing_canonical.h5ad', 'X_scprint2'),
    ('Geneformer', 'geneformer_covid_canonical.h5ad', 'X_geneformer', 'geneformer_chongqing_canonical.h5ad', 'X_geneformer'),
    ('scGPT',      'scgpt_covid_canonical.h5ad',      'X_scgpt',      'scgpt_chongqing_canonical.h5ad', 'X_scgpt'),
]
# no Chongqing embedding for these -- sepsis pool stays at legacy 12
LEGACY_ONLY_EMBEDDINGS = [
    ('TranscriptFormer', 'transcriptformer_covid_canonical.h5ad', 'X_transcriptformer'),
    ('STACK',            'stack_covid_canonical.h5ad',            'X_stack'),
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
COMBES_DROP = {'GSM4995425', 'GSM4995426', 'GSM4995430'}  # pooled multi-patient libs
JIANG_DATA_NO = '302004data18'  # Job J (2026-09-07): no per-patient outcome record in
# the source paper (GSE151263) -- "survivor" labels were study-design inference, not
# verifiable. Excluded entirely from the mortality task (7 donors / 476 cells).


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


def load_platelet_embedding(path, obsm_key):
    with h5py.File(path, 'r') as h:
        X = h[f'obsm/{obsm_key}'][:]
        disease = read_obs(h, 'Disease')
        severity = read_obs(h, 'covid_severity')
        sample_id = read_obs(h, 'Sample_ID')
        data_no = read_obs(h, 'Data_NO')
    keep = (np.isin(disease, KEEP_DISEASE) & np.isin(severity, KEEP_SEVERITY) &
            ~np.isin(sample_id, list(COMBES_DROP)) & (data_no != JIANG_DATA_NO))
    sid = sample_id[keep]
    return X[keep], disease[keep], severity[keep], sid, to_patient_id(sid)


def load_chongqing_embedding(path, obsm_key):
    with h5py.File(path, 'r') as h:
        X = h[f'obsm/{obsm_key}'][:]
        severity = read_obs(h, 'covid_severity')
        sample_id = read_obs(h, 'Sample_ID')
    disease = np.full(len(severity), 'Sepsis', dtype=object)
    keep = np.isin(severity, KEEP_SEVERITY)
    sid = sample_id[keep]
    return X[keep], disease[keep], severity[keep], sid, sid.astype(str)  # patient_id == sample_id


def load_combined(fm, platelet_fname, platelet_key, chongqing_fname, chongqing_key):
    Xp, dp, sp_, sidp, pidp = load_platelet_embedding(os.path.join(ROOT, platelet_fname), platelet_key)
    cpath = os.path.join(CHONGQING_ROOT, chongqing_fname)
    if not os.path.exists(cpath):
        print(f'  [{fm}] no Chongqing embedding at {cpath} -- legacy 12-patient sepsis pool only',
              flush=True)
        return Xp, dp, sp_, sidp, pidp
    Xc, dc, sc_, sidc, pidc = load_chongqing_embedding(cpath, chongqing_key)
    assert Xp.shape[1] == Xc.shape[1], f'{fm}: dim mismatch platelet={Xp.shape[1]} chongqing={Xc.shape[1]}'
    X = np.concatenate([Xp, Xc], axis=0)
    disease = np.concatenate([dp, dc])
    severity = np.concatenate([sp_, sc_])
    sample_id = np.concatenate([sidp, sidc])
    patient_id = np.concatenate([pidp, pidc])
    print(f'  [{fm}] platelet {Xp.shape[0]} cells ({len(np.unique(pidp[dp=="Sepsis"]))} legacy sepsis '
          f'patients) + Chongqing {Xc.shape[0]} cells ({len(np.unique(pidc))} patients) = '
          f'{len(np.unique(patient_id[disease=="Sepsis"]))} total sepsis patients', flush=True)
    return X, disease, severity, sample_id, patient_id


def _first_occurrence_index(labels, common):
    """Map each label in `common` to the position of its first occurrence in
    `labels`. Used instead of pd.Series(...).reindex(common) because `labels`
    can carry duplicates (e.g. multiple Ensembl IDs sharing one gene symbol
    on the Chongqing side), which a plain reindex on a duplicate-indexed
    Series cannot handle (ValueError: cannot reindex on an axis with
    duplicate labels)."""
    seen = {}
    for i, lab in enumerate(labels):
        if lab not in seen:
            seen[lab] = i
    return np.array([seen[c] for c in common])


def load_pca_raw():
    """Platelet raw counts (gene SYMBOLS) + Chongqing raw counts (Ensembl IDs,
    gene_symbol column) -- reindexed to the intersection of symbols."""
    state_path = os.path.join(ROOT, 'state_covid_canonical.h5ad')
    with h5py.File(state_path, 'r') as h:
        if isinstance(h['X'], h5py.Group):
            data, indices, indptr = h['X/data'][:], h['X/indices'][:], h['X/indptr'][:]
            shape = tuple(h['X'].attrs.get('shape', [len(indptr) - 1, indices.max() + 1]))
            Xp_full = sparse.csr_matrix((data, indices, indptr), shape=shape)
        else:
            Xp_full = sparse.csr_matrix(h['X'][:])
        disease = read_obs(h, 'Disease')
        severity = read_obs(h, 'covid_severity')
        sample_id = read_obs(h, 'Sample_ID')
        data_no = read_obs(h, 'Data_NO')
        p_genes = h['var/_index'][:].astype(str)
    keep = (np.isin(disease, KEEP_DISEASE) & np.isin(severity, KEEP_SEVERITY) &
            ~np.isin(sample_id, list(COMBES_DROP)) & (data_no != JIANG_DATA_NO))
    Xp_full, disease, severity, sample_id = Xp_full[keep], disease[keep], severity[keep], sample_id[keep]
    patient_id_p = to_patient_id(sample_id)

    import scanpy as sc
    c = sc.read_h5ad(CHONGQING_RAW)
    c_symbols = c.var['gene_symbol'].to_numpy()
    c_severity = c.obs['covid_severity'].astype(str).to_numpy()
    c_sample_id = c.obs['Sample_ID'].astype(str).to_numpy()
    c_keep = np.isin(c_severity, KEEP_SEVERITY)
    Xc_full = c.X[c_keep]
    c_severity, c_sample_id = c_severity[c_keep], c_sample_id[c_keep]

    common = np.intersect1d(p_genes, c_symbols)
    print(f'  PCA gene alignment: platelet {len(p_genes)}, chongqing {len(c_symbols)}, '
          f'common {len(common)}', flush=True)
    p_idx = _first_occurrence_index(p_genes, common)
    c_idx = _first_occurrence_index(c_symbols, common)
    Xp_aligned = Xp_full[:, p_idx]
    Xc_aligned = Xc_full[:, c_idx]

    X = sparse.vstack([Xp_aligned, Xc_aligned]).tocsr()
    disease_all = np.concatenate([disease, np.full(len(c_severity), 'Sepsis', dtype=object)])
    severity_all = np.concatenate([severity, c_severity])
    sample_id_all = np.concatenate([sample_id, c_sample_id])
    patient_id_all = np.concatenate([patient_id_p, c_sample_id])  # Chongqing: sample_id == patient_id
    return X, disease_all, severity_all, sample_id_all, patient_id_all


def make_sepsis_holdout(disease, severity, patient_id, seed=RANDOM_STATE):
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
             'n_test_fatal': n_test_fatal, 'n_test_surv': n_test_surv,
             'n_total_sepsis_patients': int(is_sepsis.sum() and len(np.unique(patient_id[is_sepsis])))})


def patient_bootstrap_auc(y_true, y_pred, patient_ids, n_iter=N_BOOTSTRAP, seed=0):
    rng = np.random.default_rng(seed)
    unique_patients = np.unique(patient_ids)
    patient_to_idx = {p: np.where(patient_ids == p)[0] for p in unique_patients}
    aucs = []
    for _ in range(n_iter):
        sampled = rng.choice(unique_patients, size=len(unique_patients), replace=True)
        cell_idx = np.concatenate([patient_to_idx[p] for p in sampled])
        y_t, y_p = y_true[cell_idx], y_pred[cell_idx]
        if len(np.unique(y_t)) < 2:
            continue
        aucs.append(roc_auc_score(y_t, y_p))
    if not aucs:
        return float('nan'), float('nan'), float('nan')
    aucs = np.asarray(aucs)
    return float(np.mean(aucs)), float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))


def fit_predict(X_tr, y_tr, X_te, standardize=True):
    if standardize:
        scaler = StandardScaler(with_mean=not sparse.issparse(X_tr))
        X_tr = scaler.fit_transform(X_tr)
        X_te = scaler.transform(X_te)
    clf = LogisticRegression(C=1.0, solver='lbfgs', class_weight='balanced',
                              max_iter=2000, n_jobs=1)
    clf.fit(X_tr, y_tr)
    return clf.predict_proba(X_te)[:, 1]


def adaptation_curve(name, X, y, test_mask, sepsis_train_mask, covid_train_mask,
                      patient_id, n_sepsis_grid, n_seeds, transform=None,
                      n_bootstrap=N_BOOTSTRAP, test_set='expanded_96pt'):
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

            X_tr, y_tr = X[train_idx], y[train_idx]
            X_te, y_te = X[test_idx], y[test_idx]

            if len(np.unique(y_tr)) < 2 or len(np.unique(y_te)) < 2:
                continue

            if transform is not None:
                X_tr, X_te = transform(X_tr, X_te, seed=seed)
                p = fit_predict(X_tr, y_tr, X_te, standardize=False)
            else:
                p = fit_predict(X_tr, y_tr, X_te, standardize=True)

            point_auc = float(roc_auc_score(y_te, p))
            boot_mean, ci_lo, ci_hi = patient_bootstrap_auc(
                y_te, p, te_patient_ids, n_iter=n_bootstrap, seed=seed)

            rows.append({'method': name, 'test_set': test_set,
                         'n_sepsis_req': n_req, 'n_sepsis': int(actual),
                         'seed': seed, 'n_covid_train': int(len(covid_train_idx)),
                         'n_sepsis_test': int(test_mask.sum()),
                         'n_sepsis_test_patients': int(len(np.unique(te_patient_ids))),
                         'auc_point': point_auc, 'auc_boot_mean': boot_mean,
                         'auc_boot_ci_lo': ci_lo, 'auc_boot_ci_hi': ci_hi})
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
    g = df.groupby(['method', 'test_set', 'n_sepsis_req', 'n_sepsis'])
    return g.agg(
        point_auc_mean=('auc_point', 'mean'), point_auc_std=('auc_point', 'std'),
        boot_auc_mean=('auc_boot_mean', 'mean'), boot_ci_lo_mean=('auc_boot_ci_lo', 'mean'),
        boot_ci_hi_mean=('auc_boot_ci_hi', 'mean'), n_training_seeds=('seed', 'count'),
    ).reset_index()


def main():
    all_dfs = []
    split_doc_written = False

    for fm, platelet_fname, platelet_key, chongqing_fname, chongqing_key in EMBEDDINGS:
        ppath = os.path.join(ROOT, platelet_fname)
        if not os.path.exists(ppath):
            print(f'SKIP missing platelet input: {ppath}', flush=True)
            continue
        print(f'\n=== {fm} (Chongqing-expanded) ===', flush=True)
        X, disease, severity, sample_id, patient_id = load_combined(
            fm, platelet_fname, platelet_key, chongqing_fname, chongqing_key)
        tm, stm, ctm, test_patients, y, split_doc = make_sepsis_holdout(
            disease, severity, patient_id, seed=RANDOM_STATE)
        if not split_doc_written:
            print(f'Sepsis test patients:  {len(test_patients)}', flush=True)
            print(f'Sepsis test cells:     {int(tm.sum())}', flush=True)
            print(f'Sepsis train pool:     {int(stm.sum())} cells from '
                  f'{len(np.unique(patient_id[stm]))} patients', flush=True)
            with open(os.path.join(OUT_DIR, 'split_documentation.txt'), 'w') as f:
                f.write(str(split_doc) + '\n')
            split_doc_written = True
        df = adaptation_curve(fm, X, y, tm, stm, ctm, patient_id,
                               N_SEPSIS_GRID, N_SEEDS, transform=None, n_bootstrap=N_BOOTSTRAP)
        df.to_csv(os.path.join(OUT_DIR, f'{fm}.csv'), index=False)
        all_dfs.append(df)

    legacy_split_doc_written = False
    for fm, fname, obsm_key in LEGACY_ONLY_EMBEDDINGS:
        path = os.path.join(ROOT, fname)
        if not os.path.exists(path):
            print(f'SKIP missing input: {path}', flush=True)
            continue
        print(f'\n=== {fm} (LEGACY_ONLY test set -- 12-patient sepsis pool, no Chongqing '
              f'embedding -- NOT the same test set as the 6 expanded_96pt methods above) ===',
              flush=True)
        X, disease, severity, sample_id, patient_id = load_platelet_embedding(path, obsm_key)
        tm, stm, ctm, test_patients, y, split_doc = make_sepsis_holdout(
            disease, severity, patient_id, seed=RANDOM_STATE)
        if not legacy_split_doc_written:
            print(f'[legacy_5pt] Sepsis test patients:  {len(test_patients)} '
                  f'(fatal={split_doc["n_test_fatal"]}, survived={split_doc["n_test_surv"]})',
                  flush=True)
            print(f'[legacy_5pt] Sepsis test cells:     {int(tm.sum())}', flush=True)
            print(f'[legacy_5pt] Sepsis train pool:     {int(stm.sum())} cells from '
                  f'{len(np.unique(patient_id[stm]))} patients', flush=True)
            with open(os.path.join(OUT_DIR, 'split_documentation_legacy_only.txt'), 'w') as f:
                f.write(str(split_doc) + '\n')
            legacy_split_doc_written = True
        df = adaptation_curve(fm, X, y, tm, stm, ctm, patient_id,
                               N_SEPSIS_GRID, N_SEEDS, transform=None, n_bootstrap=N_BOOTSTRAP,
                               test_set='legacy_5pt')
        df.to_csv(os.path.join(OUT_DIR, f'{fm}.csv'), index=False)
        all_dfs.append(df)

    print(f'\n=== PCA-{PCA_N_COMP} + LR (from-scratch baseline, Chongqing-expanded) ===', flush=True)
    Xraw, draw, sraw, sidraw, pidraw = load_pca_raw()
    tm, stm, ctm, _, yraw, _ = make_sepsis_holdout(draw, sraw, pidraw, seed=RANDOM_STATE)
    df_pca = adaptation_curve(f'PCA-{PCA_N_COMP}+LR', Xraw, yraw, tm, stm, ctm, pidraw,
                               N_SEPSIS_GRID, N_SEEDS, transform=pca_transform_factory(PCA_N_COMP),
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
