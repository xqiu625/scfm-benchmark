#!/usr/bin/env python3
"""
Job N -- Leave-one-cohort-out (LOBO) within-COVID on v2 data (external
review C12, secondary). Only cohorts (Data_NO) with >=1 fatal AND >=1
survivor patient get a computable held-out AUC; other cohorts stay in the
training pool only (their cells are never scored).

Same filters/patient-map/severity-fix as rerun_covid_sepsis_mortality_patient_v2.py,
restricted to COVID-19 disease. LR only (class_weight=balanced), single
deterministic split per cohort (LOBO has no fold randomness to repeat).
"""
import os
import h5py
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

ROOT = '/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c'
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2/lobo_v2'
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

    all_rows = []
    for fm, fname, obsm_key in EMBEDDINGS:
        print(f'\n=== {fm} ===', flush=True)
        X, y, patient_id, data_no = load_fm(fname, obsm_key, pmap)
        cohort_df = pd.DataFrame({'patient_id': patient_id, 'data_no': data_no, 'y': y})
        cohort_summary = cohort_df.groupby(['data_no', 'patient_id'])['y'].max().reset_index()
        cohort_class_counts = cohort_summary.groupby('data_no')['y'].agg(['sum', 'count'])
        computable = cohort_class_counts[(cohort_class_counts['sum'] >= 1) &
                                          (cohort_class_counts['count'] - cohort_class_counts['sum'] >= 1)].index

        print(f'  cohorts total={cohort_class_counts.shape[0]}, computable (both classes)={len(computable)}',
              flush=True)

        for cohort in computable:
            train_mask = data_no != cohort
            test_mask = data_no == cohort
            if train_mask.sum() == 0 or test_mask.sum() == 0:
                continue
            lr = LogisticRegression(C=1.0, solver='lbfgs', class_weight='balanced', max_iter=1000)
            lr.fit(X[train_mask], y[train_mask])
            p_test = lr.predict_proba(X[test_mask])[:, 1]
            # aggregate to patient within held-out cohort
            test_df = pd.DataFrame({'patient_id': patient_id[test_mask], 'y': y[test_mask],
                                     'score': p_test})
            pat = test_df.groupby('patient_id').agg(y=('y', 'max'), score=('score', 'mean'))
            auc = safe_auc(pat['y'].to_numpy(), pat['score'].to_numpy())
            n_pat = len(pat)
            n_fatal = int(pat['y'].sum())
            all_rows.append({'fm': fm, 'cohort_data_no': cohort, 'held_out_auc': round(auc, 4)
                              if not np.isnan(auc) else np.nan,
                              'n_patients_held_out': n_pat, 'n_fatal_held_out': n_fatal})
            print(f'    LOBO {cohort}: n={n_pat} fatal={n_fatal} AUC={auc:.4f}', flush=True)

    lobo_df = pd.DataFrame(all_rows)
    lobo_df.to_csv(os.path.join(OUT_DIR, 'lobo_per_cohort_auc.csv'), index=False)

    print('\n=== Job N: LOBO summary, mean/sigma per FM across computable cohorts ===')
    summary = lobo_df.groupby('fm')['held_out_auc'].agg(['mean', 'std', 'count'])
    print(summary.to_string())
    summary.to_csv(os.path.join(OUT_DIR, 'lobo_fm_summary.csv'))

    print(f'\nSAVED -> {OUT_DIR}', flush=True)


if __name__ == '__main__':
    main()
