#!/usr/bin/env python3
"""
Re-run primary mortality benchmark on COVID + Sepsis subset only.

Filters every Option C FM embedding to Disease in {'COVID-19', 'Sepsis'}
(drops Flu, HC/healthy, Lung_infect, Non_covid, SLE -- 5,882 cells, 91 donors)
and re-runs LR + RF under 5-fold patient-level StratifiedGroupKFold on Sample_ID.

After filter: 42,095 cells / 322 donors / 8,850 fatal / 33,245 survived
(per-FM cell counts vary by ~13-868 cells depending on each FM's native preprocessing).

Memory-safe: h5py-only access (never loads X). Tested OK with 64G --mem.

Output: 05_RESULTS/option_c_primary_benchmark/covid_sepsis_only_mortality.csv
"""
import h5py
import numpy as np
import pandas as pd
import os

from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score, balanced_accuracy_score

ROOT = '/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c'
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/option_c_primary_benchmark'

EMBEDDINGS = [
    ('STATE',            'state_covid_canonical.h5ad',            'X_state'),
    ('UCE',              'uce_covid_canonical.h5ad',              'X_uce'),
    ('TranscriptFormer', 'transcriptformer_covid_canonical.h5ad', 'X_transcriptformer'),
    ('scPRINT-2',        'scprint2_covid_canonical.h5ad',         'X_scprint2'),
    ('STACK',            'stack_covid_canonical.h5ad',            'X_stack'),
    ('Geneformer',       'geneformer_covid_canonical.h5ad',       'X_geneformer'),
    ('scGPT',            'scgpt_covid_canonical.h5ad',             'X_scgpt'),
]

KEEP_DISEASE = ['COVID-19', 'Sepsis']
KEEP_SEVERITY = ['fatal', 'mild', 'moderate', 'recovered', 'severe']


def read_categorical(group):
    """Read an obs categorical column from an h5py group."""
    if isinstance(group, h5py.Group) and 'categories' in group:
        cats = group['categories'][:].astype(str)
        codes = group['codes'][:]
        return cats[codes]
    return group[:].astype(str)


def evaluate_fm(fm, fname, obsm_key):
    path = os.path.join(ROOT, fname)
    if not os.path.exists(path):
        print(f'SKIP: {path} (not found)', flush=True)
        return []

    print(f'\n=== {fm}: {fname} ===', flush=True)
    with h5py.File(path, 'r') as h:
        n_pre = h['obs/_index'].shape[0]
        if obsm_key not in list(h.get('obsm', [])):
            print(f'  ERROR: obsm[{obsm_key}] missing -- skipping', flush=True)
            return []
        X = h[f'obsm/{obsm_key}'][:]
        disease = read_categorical(h['obs/Disease'])
        severity = read_categorical(h['obs/covid_severity'])
        sample_id = read_categorical(h['obs/Sample_ID'])

    keep = np.isin(disease, KEEP_DISEASE) & np.isin(severity, KEEP_SEVERITY)
    Xf = X[keep]
    sev_f = severity[keep]
    sid_f = sample_id[keep]
    y = (sev_f == 'fatal').astype(int)

    n_donors = len(np.unique(sid_f))
    n_fatal = int(y.sum())
    n_surv = int((y == 0).sum())
    print(
        f'  filter: {n_pre} -> {Xf.shape[0]} cells | '
        f'{n_donors} donors | fatal={n_fatal}, survived={n_surv}',
        flush=True,
    )

    skf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)
    classifiers = [
        (
            'LR',
            lambda: LogisticRegression(
                C=1.0, solver='lbfgs', class_weight='balanced', max_iter=1000
            ),
        ),
        (
            'RF',
            lambda: RandomForestClassifier(
                n_estimators=100,
                class_weight='balanced',
                n_jobs=-1,
                random_state=42,
            ),
        ),
    ]

    rows = []
    for clf_name, factory in classifiers:
        fold_aucs, fold_baccs = [], []
        for fold, (tr, te) in enumerate(skf.split(Xf, y, groups=sid_f)):
            clf = factory()
            clf.fit(Xf[tr], y[tr])
            p = clf.predict_proba(Xf[te])[:, 1]
            yh = clf.predict(Xf[te])
            fold_aucs.append(roc_auc_score(y[te], p))
            fold_baccs.append(balanced_accuracy_score(y[te], yh))
        m_auc = float(np.mean(fold_aucs))
        s_auc = float(np.std(fold_aucs))
        m_bacc = float(np.mean(fold_baccs))
        print(
            f'  {clf_name}: AUC={m_auc:.3f} +/- {s_auc:.3f}   BalAcc={m_bacc:.3f}',
            flush=True,
        )
        rows.append({
            'fm': fm,
            'classifier': clf_name,
            'n_cells': int(Xf.shape[0]),
            'n_donors': n_donors,
            'n_fatal': n_fatal,
            'n_survived': n_surv,
            'auc_mean': m_auc,
            'auc_std': s_auc,
            'bacc_mean': m_bacc,
            'fold_aucs': ';'.join(f'{a:.3f}' for a in fold_aucs),
        })
    return rows


def main():
    print(f'COVID + Sepsis only mortality re-run', flush=True)
    print(f'  ROOT  = {ROOT}', flush=True)
    print(f'  Disease filter: {KEEP_DISEASE}', flush=True)
    print(f'  Severity filter: {KEEP_SEVERITY}', flush=True)

    all_rows = []
    for fm, fname, obsm_key in EMBEDDINGS:
        all_rows.extend(evaluate_fm(fm, fname, obsm_key))

    df = pd.DataFrame(all_rows)
    os.makedirs(OUT_DIR, exist_ok=True)
    out_path = os.path.join(OUT_DIR, 'covid_sepsis_only_mortality.csv')
    df.to_csv(out_path, index=False)
    print(f'\nSAVED: {out_path}', flush=True)

    print('\n=== Summary (LR only) ===', flush=True)
    lr_summary = df[df['classifier'] == 'LR'][[
        'fm', 'n_cells', 'n_donors', 'n_fatal', 'n_survived', 'auc_mean', 'auc_std'
    ]].sort_values('auc_mean', ascending=False)
    print(lr_summary.to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
