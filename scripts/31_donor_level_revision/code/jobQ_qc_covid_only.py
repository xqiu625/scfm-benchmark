#!/usr/bin/env python3
"""
Job Q (COVID-only variant of Job K) -- QC-metrics-only baseline restricted to
COVID-19 patients only (external review A1 cohort-scope follow-up).

Features per cell: log1p(nCount_RNA), log1p(nFeature_RNA), pct_counts_mt.
Same filters/patient-map/severity-fix/CV protocol as
rerun_covid_sepsis_mortality_patient_v2.py, but with these 3 QC features
in place of an FM embedding.

Decision rule (same as the null model): QC-only AUC >= 0.7 within-COVID ->
the main result would be a technical artifact -- flag immediately.
"""
import os
import h5py
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

from repeated_cv import (DEFAULT_SEEDS, run_repeated_kfold, repeated_patient_auc,
                          summarize_repeated)

RAW = '/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/human_platelet_covid_severity_raw_counts.h5ad'
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2/qc_baseline'
PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')

KEEP_DISEASE = ['COVID-19']  # Job Q: COVID-only variant, external review A1 cohort scope
KEEP_SEVERITY = ['fatal', 'mild', 'moderate', 'recovered', 'severe']

SEVERITY_FIX = {
    'GSM4995436': 'moderate',
    'GSM4995442': 'severe',
    'GSM4995432': 'moderate',
}


def read_categorical(group):
    if isinstance(group, h5py.Group) and 'categories' in group:
        cats = group['categories'][:].astype(str)
        codes = group['codes'][:]
        return cats[codes]
    return group[:].astype(str)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    pmap = pd.read_csv(PATIENT_MAP)[['sample_id', 'patient_id']]

    print('Loading raw counts h5ad...', flush=True)
    with h5py.File(RAW, 'r') as h:
        cell_id = h['obs/_index'][:].astype(str)
        disease = read_categorical(h['obs/Disease'])
        severity = read_categorical(h['obs/covid_severity'])
        sample_id = read_categorical(h['obs/Sample_ID'])
        n_counts = h['obs/nCount_RNA'][:].astype(float)
        n_genes = h['obs/nFeature_RNA'][:].astype(float)
        var_names = h['var/_index'][:].astype(str)
        mt_mask = np.array([g.startswith('MT-') for g in var_names])
        print(f'  {mt_mask.sum()} MT genes found', flush=True)
        X_all = h['X'][:]  # dense, (n_cells, n_genes)
        mt_counts = X_all[:, mt_mask].sum(axis=1).astype(float)
        total_counts_x = X_all.sum(axis=1).astype(float)
        del X_all

    pct_mito = np.divide(mt_counts, total_counts_x, out=np.zeros_like(mt_counts),
                          where=total_counts_x > 0) * 100.0

    severity = severity.copy()
    for sid, new_sev in SEVERITY_FIX.items():
        severity[sample_id == sid] = new_sev

    keep = np.isin(disease, KEEP_DISEASE) & np.isin(severity, KEEP_SEVERITY)
    print(f'Pre-filter: {len(cell_id)} cells; post-filter: {keep.sum()} cells', flush=True)

    df = pd.DataFrame({
        'cell_id': cell_id[keep], 'sample_id': sample_id[keep],
        'disease': disease[keep], 'severity': severity[keep],
        'n_counts': n_counts[keep], 'n_genes': n_genes[keep], 'pct_mito': pct_mito[keep],
    })
    n_before = len(df)
    df = df.merge(pmap, on='sample_id', how='inner')
    print(f'  dropped {n_before - len(df)} cells not in v2 map (Combes pooled GSMs)', flush=True)

    df['y'] = (df['severity'] == 'fatal').astype(int)
    df['log_counts'] = np.log1p(df['n_counts'])
    df['log_genes'] = np.log1p(df['n_genes'])

    n_patients = df['patient_id'].nunique()
    n_fatal_patients = df.groupby('patient_id')['y'].max().sum()
    print(f'Final: {len(df)} cells, {n_patients} patients, {int(n_fatal_patients)} fatal patients',
          flush=True)

    patient_qc = df.groupby('patient_id').agg(
        y=('y', 'max'), disease=('disease', 'first'),
        n_cells=('cell_id', 'count'),
        median_counts=('n_counts', 'median'), median_genes=('n_genes', 'median'),
        median_pct_mito=('pct_mito', 'median'),
    ).reset_index()
    patient_qc.to_csv(os.path.join(OUT_DIR, 'patient_qc_distributions_covid_only.csv'), index=False)
    print('\n=== per-patient QC, by disease x outcome (median of medians) ===')
    print(patient_qc.groupby(['disease', 'y'])[['median_counts', 'median_genes',
                                                 'median_pct_mito', 'n_cells']].median()
          .to_string())

    X = df[['log_counts', 'log_genes', 'pct_mito']].to_numpy()
    y = df['y'].to_numpy()
    groups = df['patient_id'].to_numpy()

    def fit_predict(X_tr, y_tr, X_te):
        lr = LogisticRegression(C=1.0, solver='lbfgs', class_weight='balanced', max_iter=1000)
        lr.fit(X_tr, y_tr)
        return {'p_lr': lr.predict_proba(X_te)[:, 1]}

    oof_by_seed = run_repeated_kfold(X, y, groups, fit_predict_fn=fit_predict,
                                      seeds=DEFAULT_SEEDS)
    aucs, seeds = repeated_patient_auc(oof_by_seed, groups, y, score_key='p_lr')
    headline = summarize_repeated(aucs)

    print('\n=== Job Q: QC-only baseline, COVID-only cohort, patient-level repeat-mean AUC ===')
    print(f'  mean={headline["mean"]:.4f}  sd={headline["sd"]:.4f}  n_seeds={headline["n"]}')
    print(f'  per-seed: {dict(zip(seeds, np.round(aucs, 4)))}')

    pd.DataFrame([{'analysis': 'qc_only_baseline', 'auc_repeat_mean': headline['mean'],
                   'auc_repeat_sd': headline['sd'], 'n_seeds': headline['n'],
                   'n_patients': n_patients, 'n_fatal_patients': int(n_fatal_patients)}]
                 ).to_csv(os.path.join(OUT_DIR, 'qc_baseline_auc_covid_only.csv'), index=False)
    pd.DataFrame({'seed': seeds, 'auc': aucs}).to_csv(
        os.path.join(OUT_DIR, 'qc_baseline_per_seed_covid_only.csv'), index=False)

    if headline['mean'] >= 0.7:
        print('\n*** DECISION RULE TRIGGERED: QC-only AUC >= 0.7 -- main result may be a '
              'technical artifact. STOP AND REPORT IMMEDIATELY. ***', flush=True)
    else:
        print(f'\nDecision rule NOT triggered ({headline["mean"]:.4f} < 0.7) -- '
              'supplementary null finding.', flush=True)

    print(f'\nSAVED -> {OUT_DIR}', flush=True)


if __name__ == '__main__':
    main()
