#!/usr/bin/env python3
"""Job H step 3: patient-level re-run of cross-disease transfer (stage 06).

Reuses scripts/06_batch_generalization/code/cross_disease_fatal_vs_survived.py's
severe-excluded protocol (the one backing the manuscript's tab:crossdisease /
claim C002, Job 23761950: COVID-severe excluded from the survived class,
sepsis-survived = severe, 200-iter bootstrap, LR + RF), but with the bootstrap
resampling grouped by patient_id instead of Sample_ID (2026-09-02 multi-draw
leakage fix; see reports/manuscript/revision/patient_overlap_finding.md).

PROTOCOL NOTE: uses the canonical Option C embeddings (02_EMBEDDINGS/option_c/)
for all 5 primary FMs including scGPT, for consistency with the rest of Job H's
patient-level rerun. The manuscript's scGPT cross-disease row (forward 0.352)
came from a separate 2026-04-30 run on a legacy scGPT embedding path
(02_EMBEDDINGS/scgpt/covid_scgpt_clean.h5ad, not the canonical Option C file);
this rerun does NOT reproduce that legacy path and scGPT's number here may
therefore differ from 0.352 for reasons unrelated to the patient-id fix --
report both and flag the protocol difference, do not reconcile silently.

Models: STATE, UCE, scPRINT-2, Geneformer, scGPT (primary 5-FM set).
COVID survived tokens: mild,moderate,recovered (severe EXCLUDED).
Sepsis: fatal vs severe(=survived).
Directions: COVID->Sepsis and Sepsis->COVID.

Output: 05_RESULTS/revision_donor_level/cross_disease_patient_<timestamp>.csv
"""
import json
import logging
from datetime import datetime
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.preprocessing import StandardScaler

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
log = logging.getLogger('cross_disease_fatal_patient')

BASE = Path('/bigdata/godziklab/shared/Xinru/302006')
EMB = BASE / '02_EMBEDDINGS' / 'option_c'
RESULTS = BASE / '05_RESULTS' / 'revision_donor_level'
PATIENT_MAP = (BASE / '302006_updategit' / 'scripts' / '31_donor_level_revision' /
               'results' / 'donor_meta_with_patient_id.csv')

MODEL_PATHS = {
    'STATE': (EMB / 'state_covid_canonical.h5ad', 'X_state'),
    'UCE': (EMB / 'uce_covid_canonical.h5ad', 'X_uce'),
    'scPRINT-2': (EMB / 'scprint2_covid_canonical.h5ad', 'X_scprint2'),
    'Geneformer': (EMB / 'geneformer_covid_canonical.h5ad', 'X_geneformer'),
    'scGPT': (EMB / 'scgpt_covid_canonical.h5ad', 'X_scgpt'),
}

COVID_FATAL = {'fatal'}
COVID_SURVIVED = {'mild', 'moderate', 'recovered'}   # severe EXCLUDED
SEPSIS_FATAL = {'fatal'}
SEPSIS_SURVIVED = {'severe'}
N_BOOTSTRAP = 200
PMAP = pd.read_csv(PATIENT_MAP).set_index('sample_id')['patient_id'].to_dict()


def read_categorical(group):
    if isinstance(group, h5py.Group) and 'categories' in group:
        cats = group['categories'][:].astype(str)
        codes = group['codes'][:]
        return cats[codes]
    return group[:].astype(str)


def load(fm):
    path, obsm_key = MODEL_PATHS[fm]
    with h5py.File(path, 'r') as h:
        X = h[f'obsm/{obsm_key}'][:]
        disease = read_categorical(h['obs/Disease'])
        severity = read_categorical(h['obs/covid_severity'])
        sample_id = read_categorical(h['obs/Sample_ID'])
    pid = np.array([PMAP.get(s, None) for s in sample_id], dtype=object)
    return X, disease, severity, sample_id, pid


def build_labels(disease, severity):
    is_covid = disease == 'COVID-19'
    is_sepsis = disease == 'Sepsis'
    y = np.zeros(len(disease), dtype=np.int8)
    valid = np.zeros(len(disease), dtype=bool)
    y[is_covid & np.isin(severity, list(COVID_FATAL))] = 1
    valid[is_covid & np.isin(severity, list(COVID_FATAL | COVID_SURVIVED))] = True
    y[is_sepsis & np.isin(severity, list(SEPSIS_FATAL))] = 1
    valid[is_sepsis & np.isin(severity, list(SEPSIS_FATAL | SEPSIS_SURVIVED))] = True
    return y, valid, is_covid, is_sepsis


def patient_bootstrap_auc(y_true, y_prob, patient_ids, n_iter=N_BOOTSTRAP, seed=42):
    rng = np.random.default_rng(seed)
    unique_p = np.unique(patient_ids)
    p_to_idx = {p: np.where(patient_ids == p)[0] for p in unique_p}
    aucs = []
    for _ in range(n_iter):
        sampled = rng.choice(unique_p, size=len(unique_p), replace=True)
        idx = np.concatenate([p_to_idx[p] for p in sampled])
        if len(np.unique(y_true[idx])) < 2:
            continue
        try:
            aucs.append(roc_auc_score(y_true[idx], y_prob[idx]))
        except ValueError:
            continue
    if not aucs:
        return float('nan'), float('nan'), float('nan')
    arr = np.asarray(aucs)
    return float(arr.mean()), float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))


def transfer(X, y, valid, train_mask, test_mask, test_patient_ids):
    train_sel = train_mask & valid
    test_sel = test_mask & valid
    Xtr, ytr = X[train_sel], y[train_sel]
    Xte, yte = X[test_sel], y[test_sel]
    gte = test_patient_ids[test_sel]

    log.info('  Train n=%d (fatal=%d, survived=%d) | Test n=%d (fatal=%d, survived=%d, patients=%d)',
              len(ytr), int(ytr.sum()), int((1 - ytr).sum()),
              len(yte), int(yte.sum()), int((1 - yte).sum()), len(np.unique(gte)))

    if len(np.unique(ytr)) < 2 or len(np.unique(yte)) < 2:
        log.warning('  Skipping: insufficient class diversity')
        return {}

    scaler = StandardScaler()
    Xtr_s = scaler.fit_transform(Xtr)
    Xte_s = scaler.transform(Xte)

    out = {'n_train': int(len(ytr)), 'n_test': int(len(yte)),
           'fatal_train': int(ytr.sum()), 'fatal_test': int(yte.sum()),
           'n_test_patients': int(len(np.unique(gte))),
           'n_fatal_test_patients': int(pd.DataFrame({'p': gte, 'y': yte})
                                        .groupby('p')['y'].max().sum())}

    clf_lr = LogisticRegression(max_iter=2000, random_state=42, n_jobs=-1, class_weight='balanced')
    clf_lr.fit(Xtr_s, ytr)
    p_lr = clf_lr.predict_proba(Xte_s)[:, 1]
    auc_lr = float(roc_auc_score(yte, p_lr))
    bal_lr = float(balanced_accuracy_score(yte, clf_lr.predict(Xte_s)))
    mean_lr, lo_lr, hi_lr = patient_bootstrap_auc(yte, p_lr, gte)
    out['LogReg'] = {'auc': auc_lr, 'balanced_accuracy': bal_lr,
                     'boot_auc_mean': mean_lr, 'boot_ci_lo': lo_lr, 'boot_ci_hi': hi_lr}

    clf_rf = RandomForestClassifier(n_estimators=200, random_state=42, n_jobs=-1, class_weight='balanced')
    clf_rf.fit(Xtr_s, ytr)
    p_rf = clf_rf.predict_proba(Xte_s)[:, 1]
    auc_rf = float(roc_auc_score(yte, p_rf))
    bal_rf = float(balanced_accuracy_score(yte, clf_rf.predict(Xte_s)))
    mean_rf, lo_rf, hi_rf = patient_bootstrap_auc(yte, p_rf, gte)
    out['RandomForest'] = {'auc': auc_rf, 'balanced_accuracy': bal_rf,
                           'boot_auc_mean': mean_rf, 'boot_ci_lo': lo_rf, 'boot_ci_hi': hi_rf}
    return out


def main():
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    rows = []
    full = {}
    for fm in MODEL_PATHS:
        log.info('=' * 70)
        log.info('Model: %s', fm)
        X, disease, severity, sample_id, patient_id = load(fm)
        y, valid, cov_mask, sep_mask = build_labels(disease, severity)
        # patient_id is only needed for cells that are actually used (valid);
        # off-disease / excluded-severity cells (e.g. Flu, healthy, COVID-severe)
        # were never in the COVID+Sepsis-only patient map and are dropped by
        # `valid` before use, so only require mapping completeness there.
        n_missing = int((pd.isna(patient_id) & valid).sum())
        assert n_missing == 0, f'{fm}: {n_missing} valid cells lack patient_id mapping'
        patient_id = np.where(pd.isna(patient_id), 'NA', patient_id).astype(str)
        log.info('COVID valid: %d  Sepsis valid: %d',
                 int((cov_mask & valid).sum()), int((sep_mask & valid).sum()))

        full[fm] = {}
        for direction in ['COVID->Sepsis', 'Sepsis->COVID']:
            train, test = direction.split('->')
            tm = cov_mask if train == 'COVID' else sep_mask
            te = sep_mask if test == 'Sepsis' else cov_mask
            log.info('Direction: %s', direction)
            res = transfer(X, y, valid, tm, te, patient_id)
            full[fm][direction] = res
            for clf, m in res.items():
                if clf in ('n_train', 'n_test', 'fatal_train', 'fatal_test',
                           'n_test_patients', 'n_fatal_test_patients'):
                    continue
                rows.append({
                    'model': fm, 'direction': direction, 'classifier': clf,
                    'auc': m['auc'], 'balanced_accuracy': m['balanced_accuracy'],
                    'boot_auc_mean': m['boot_auc_mean'],
                    'boot_ci_lo': m['boot_ci_lo'], 'boot_ci_hi': m['boot_ci_hi'],
                    'n_train': res['n_train'], 'n_test': res['n_test'],
                    'fatal_train': res['fatal_train'], 'fatal_test': res['fatal_test'],
                    'n_test_patients': res['n_test_patients'],
                    'n_fatal_test_patients': res['n_fatal_test_patients'],
                })
                log.info('  %s: AUC=%.4f BootMean=%.4f [%.4f, %.4f] (n_test_patients=%d, fatal=%d)',
                         clf, m['auc'], m['boot_auc_mean'], m['boot_ci_lo'], m['boot_ci_hi'],
                         res['n_test_patients'], res['n_fatal_test_patients'])

    df = pd.DataFrame(rows)
    RESULTS.mkdir(parents=True, exist_ok=True)
    csv_path = RESULTS / f'cross_disease_patient_{timestamp}.csv'
    df.to_csv(csv_path, index=False)
    (RESULTS / f'cross_disease_patient_{timestamp}.json').write_text(json.dumps(full, indent=2))
    log.info('Saved: %s', csv_path)


if __name__ == '__main__':
    main()
