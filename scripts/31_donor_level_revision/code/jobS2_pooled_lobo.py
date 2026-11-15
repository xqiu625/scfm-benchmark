#!/usr/bin/env python3
"""
Round-4 primary-test correction, item 5: concatenated-pooled LOBO.

Job N (jobN_lobo_v2.py) computed per-cohort held-out AUC but never persisted
the raw per-patient held-out predictions, so a pooled AUC can't be derived
post-hoc -- this reruns the identical LOBO loop (same filters/model/severity
fix) and additionally concatenates every computable cohort's held-out
(patient_id, y, score) into one set per FM, then computes ONE pooled AUC.

Output: paired_inference_full/lobo_pooled.csv
  (per-FM pooled AUC + n_patients_pooled + n_fatal_pooled, plus a
  per-cohort n/fatal breakdown re-saved alongside for reference)
"""
import os

import h5py
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

ROOT = '/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c'
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2/paired_inference_full'
PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')

EMBEDDINGS = [
    ('STATE',      'state_covid_canonical.h5ad',      'X_state'),
    ('UCE',        'uce_covid_canonical.h5ad',        'X_uce'),
    ('Geneformer', 'geneformer_covid_canonical.h5ad', 'X_geneformer'),
    ('scPRINT-2',  'scprint2_covid_canonical.h5ad',   'X_scprint2'),
    ('scGPT',      'scgpt_covid_canonical.h5ad',      'X_scgpt'),
]

KEEP_SEVERITY = ['fatal', 'mild', 'moderate', 'recovered', 'severe']
SEVERITY_FIX = {
    'GSM4995436': 'moderate', 'GSM4995442': 'severe', 'GSM4995432': 'moderate',
}


def read_categorical(group):
    if isinstance(group, h5py.Group) and 'categories' in group:
        cats = group['categories'][:].astype(str)
        codes = group['codes'][:]
        return cats[codes]
    return group[:].astype(str)


def safe_auc(y, p):
    if len(np.unique(y)) < 2:
        return np.nan
    return roc_auc_score(y, p)


def load_fm(fname, obsm_key, pmap):
    path = os.path.join(ROOT, fname)
    with h5py.File(path, 'r') as h:
        cell_id = h['obs/_index'][:].astype(str)
        X = h[f'obsm/{obsm_key}'][:]
        disease = read_categorical(h['obs/Disease'])
        severity = read_categorical(h['obs/covid_severity'])
        sample_id = read_categorical(h['obs/Sample_ID'])
    severity = severity.copy()
    for sid, new_sev in SEVERITY_FIX.items():
        severity[sample_id == sid] = new_sev

    keep = (disease == 'COVID-19') & np.isin(severity, KEEP_SEVERITY)
    df = pd.DataFrame({'cell_id': cell_id[keep], 'sample_id': sample_id[keep],
                        'severity': severity[keep]})
    df['_row'] = np.arange(keep.sum())
    Xk = X[keep]
    df = df.merge(pmap, on='sample_id', how='inner')
    Xf = Xk[df['_row'].to_numpy()]
    y = (df['severity'].to_numpy() == 'fatal').astype(int)
    return Xf, y, df['patient_id'].to_numpy(), df['Data_NO'].to_numpy()


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    pmap = pd.read_csv(PATIENT_MAP)
    pmap = pmap[pmap['disease'] == 'COVID-19'][['sample_id', 'patient_id', 'Data_NO']]

    pooled_rows = []
    cohort_rows = []
    for fm, fname, obsm_key in EMBEDDINGS:
        print(f'\n=== {fm} ===', flush=True)
        X, y, patient_id, data_no = load_fm(fname, obsm_key, pmap)
        cohort_df = pd.DataFrame({'patient_id': patient_id, 'data_no': data_no, 'y': y})
        cohort_summary = cohort_df.groupby(['data_no', 'patient_id'])['y'].max().reset_index()
        cohort_class_counts = cohort_summary.groupby('data_no')['y'].agg(['sum', 'count'])
        computable = cohort_class_counts[(cohort_class_counts['sum'] >= 1) &
                                          (cohort_class_counts['count'] - cohort_class_counts['sum'] >= 1)].index

        held_out_patients = []
        held_out_y = []
        held_out_score = []
        for cohort in computable:
            train_mask = data_no != cohort
            test_mask = data_no == cohort
            if train_mask.sum() == 0 or test_mask.sum() == 0:
                continue
            lr = LogisticRegression(C=1.0, solver='lbfgs', class_weight='balanced', max_iter=1000)
            lr.fit(X[train_mask], y[train_mask])
            p_test = lr.predict_proba(X[test_mask])[:, 1]
            test_df = pd.DataFrame({'patient_id': patient_id[test_mask], 'y': y[test_mask],
                                     'score': p_test})
            pat = test_df.groupby('patient_id').agg(y=('y', 'max'), score=('score', 'mean'))
            held_out_patients.extend(pat.index.tolist())
            held_out_y.extend(pat['y'].tolist())
            held_out_score.extend(pat['score'].tolist())
            n_pat = len(pat); n_fatal = int(pat['y'].sum())
            cohort_rows.append({'fm': fm, 'cohort_data_no': cohort,
                                 'n_patients_held_out': n_pat, 'n_fatal_held_out': n_fatal})
            print(f'    cohort {cohort}: n={n_pat} fatal={n_fatal} (added to pool)', flush=True)

        y_pooled = np.array(held_out_y)
        score_pooled = np.array(held_out_score)
        # a patient can appear in >1 held-out cohort only if data_no isn't
        # unique per patient -- verify uniqueness before pooling, since a
        # duplicate would corrupt an AUC computed over "patients"
        dup = pd.Series(held_out_patients).duplicated().sum()
        auc_pooled = safe_auc(y_pooled, score_pooled)
        pooled_rows.append({'fm': fm, 'n_patients_pooled': len(y_pooled),
                             'n_fatal_pooled': int(y_pooled.sum()),
                             'n_duplicate_patients_across_cohorts': int(dup),
                             'pooled_auc': round(auc_pooled, 4) if not np.isnan(auc_pooled) else np.nan,
                             'n_cohorts_pooled': len(computable)})
        print(f'  POOLED {fm}: n={len(y_pooled)} fatal={int(y_pooled.sum())} '
              f'dup={dup} AUC={auc_pooled:.4f}', flush=True)

    pooled_df = pd.DataFrame(pooled_rows)
    pooled_df.to_csv(os.path.join(OUT_DIR, 'lobo_pooled.csv'), index=False)
    print('\n=== pooled LOBO AUC per FM ===')
    print(pooled_df.to_string(index=False))

    cohort_df = pd.DataFrame(cohort_rows)
    cohort_df.to_csv(os.path.join(OUT_DIR, 'lobo_pooled_per_cohort_n.csv'), index=False)

    print(f'\nSAVED -> {OUT_DIR}', flush=True)


if __name__ == '__main__':
    main()
