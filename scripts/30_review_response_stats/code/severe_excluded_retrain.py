#!/usr/bin/env python3
"""
Severe-excluded RETRAIN (roadmap F, definitive version).

Stage 30's severe sensitivity was eval-side (models trained on the full set, severe
negatives dropped only at scoring). This retrains from scratch with severe negatives
EXCLUDED from both train and test, applying the within-disease primary-task label rule
(Sec 2.5.1: negatives = mild/moderate/recovered; severe excluded) uniformly.

Confirms whether STATE's lead and the FM ranking are robust to the label definition.

Output: 05_RESULTS/option_c_primary_benchmark/severe_excluded_mortality.csv
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
    ('STATE',      'state_covid_canonical.h5ad',     'X_state'),
    ('UCE',        'uce_covid_canonical.h5ad',       'X_uce'),
    ('scPRINT-2',  'scprint2_covid_canonical.h5ad',  'X_scprint2'),
    ('Geneformer', 'geneformer_covid_canonical.h5ad', 'X_geneformer'),
    ('scGPT',      'scgpt_covid_canonical.h5ad',     'X_scgpt'),
]
KEEP_DISEASE = ['COVID-19', 'Sepsis']
# severe-EXCLUDED: negatives are mild/moderate/recovered only (primary-task rule)
KEEP_SEVERITY = ['fatal', 'mild', 'moderate', 'recovered']


def read_categorical(group):
    if isinstance(group, h5py.Group) and 'categories' in group:
        cats = group['categories'][:].astype(str)
        return cats[group['codes'][:]]
    return group[:].astype(str)


def evaluate(fm, fname, key):
    path = os.path.join(ROOT, fname)
    if not os.path.exists(path):
        print(f'SKIP {fm}', flush=True)
        return []
    with h5py.File(path, 'r') as h:
        X = h[f'obsm/{key}'][:]
        disease = read_categorical(h['obs/Disease'])
        severity = read_categorical(h['obs/covid_severity'])
        sid = read_categorical(h['obs/Sample_ID'])
    keep = np.isin(disease, KEEP_DISEASE) & np.isin(severity, KEEP_SEVERITY)
    Xf, sevf, sidf = X[keep], severity[keep], sid[keep]
    y = (sevf == 'fatal').astype(int)
    rows = []
    for clf_name, factory in [
        ('LR', lambda: LogisticRegression(C=1.0, solver='lbfgs',
                                          class_weight='balanced', max_iter=1000)),
        ('RF', lambda: RandomForestClassifier(n_estimators=100, class_weight='balanced',
                                              n_jobs=-1, random_state=42)),
    ]:
        skf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)
        aucs, baccs = [], []
        for tr, te in skf.split(Xf, y, groups=sidf):
            clf = factory(); clf.fit(Xf[tr], y[tr])
            p = clf.predict_proba(Xf[te])[:, 1]
            aucs.append(roc_auc_score(y[te], p))
            baccs.append(balanced_accuracy_score(y[te], clf.predict(Xf[te])))
        print(f'  {fm} {clf_name}: AUC={np.mean(aucs):.3f}+/-{np.std(aucs):.3f}', flush=True)
        rows.append({'fm': fm, 'classifier': clf_name, 'n_cells': int(Xf.shape[0]),
                     'n_donors': len(np.unique(sidf)), 'n_fatal': int(y.sum()),
                     'n_survived': int((y == 0).sum()),
                     'auc_mean': float(np.mean(aucs)), 'auc_std': float(np.std(aucs)),
                     'bacc_mean': float(np.mean(baccs)),
                     'fold_aucs': ';'.join(f'{a:.3f}' for a in aucs)})
    return rows


def main():
    print('Severe-EXCLUDED retrain (negatives = mild/moderate/recovered)', flush=True)
    allrows = []
    for fm, fn, k in EMBEDDINGS:
        allrows.extend(evaluate(fm, fn, k))
    df = pd.DataFrame(allrows)
    os.makedirs(OUT_DIR, exist_ok=True)
    out = os.path.join(OUT_DIR, 'severe_excluded_mortality.csv')
    df.to_csv(out, index=False)
    print(f'\nSAVED: {out}', flush=True)
    print(df[df.classifier == 'LR'][['fm', 'n_cells', 'n_fatal', 'n_survived',
                                     'auc_mean', 'auc_std']].sort_values(
        'auc_mean', ascending=False).to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
