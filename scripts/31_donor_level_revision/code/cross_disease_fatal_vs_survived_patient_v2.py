#!/usr/bin/env python3
"""Job H step 3: patient-level re-run of cross-disease transfer (stage 06).

Reuses scripts/06_batch_generalization/code/cross_disease_fatal_vs_survived.py's
protocol backbone (200-iter bootstrap, LR + RF), but with:
  - bootstrap resampling grouped by patient_id (2026-09-02 multi-draw fix)
  - Fix-2 (code review 2026-09-04): PRIMARY AUC is patient-level -- each
    patient's cell predictions are averaged first, then AUC + the bootstrap
    are computed on that patient-level array. Cell-weighted AUC kept as a
    secondary column (`*_auc_cell_weighted`).
  - Fix-3 (code review 2026-09-04): severe=survived is now the PRIMARY COVID
    protocol (P0 decision); severe-excluded runs alongside as a sensitivity
    variant. Every stat is reported twice, column-prefixed `primary_` /
    `sens_severe_excluded_` -- no adjudication here, that's a manuscript call.

PROTOCOL NOTE: uses the canonical Option C embeddings (02_EMBEDDINGS/option_c/)
for all 5 primary FMs including scGPT, for consistency with the rest of Job H's
patient-level rerun. The manuscript's scGPT cross-disease row (forward 0.352)
came from a separate 2026-04-30 run on a legacy scGPT embedding path
(02_EMBEDDINGS/scgpt/covid_scgpt_clean.h5ad, not the canonical Option C file);
this rerun does NOT reproduce that legacy path and scGPT's number here may
therefore differ from 0.352 for reasons unrelated to the patient-id fix --
report both and flag the protocol difference, do not reconcile silently.

Models: STATE, UCE, scPRINT-2, Geneformer, scGPT (primary 5-FM set).
Sepsis: fatal vs severe(=survived).
Directions: COVID->Sepsis and Sepsis->COVID.

Output: 05_RESULTS/revision_donor_level_v2/cross_disease_patient_<timestamp>.csv
"""
import json
import logging
from datetime import datetime
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from transfer_utils import (COVID_FATAL, COVID_SURVIVED_PRIMARY, COVID_SURVIVED_SENS,
                             SEPSIS_FATAL, SEPSIS_SURVIVED, transfer_one)

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
log = logging.getLogger('cross_disease_fatal_patient')

BASE = Path('/bigdata/godziklab/shared/Xinru/302006')
EMB = BASE / '02_EMBEDDINGS' / 'option_c'
RESULTS = BASE / '05_RESULTS' / 'revision_donor_level_v2'
PATIENT_MAP = (BASE / '302006_updategit' / 'scripts' / '31_donor_level_revision' /
               'results' / 'donor_meta_with_patient_id_v2.csv')
COMBES_DROP = {'GSM4995425', 'GSM4995426', 'GSM4995430'}  # pooled multi-patient libs

MODEL_PATHS = {
    'STATE': (EMB / 'state_covid_canonical.h5ad', 'X_state'),
    'UCE': (EMB / 'uce_covid_canonical.h5ad', 'X_uce'),
    'scPRINT-2': (EMB / 'scprint2_covid_canonical.h5ad', 'X_scprint2'),
    'Geneformer': (EMB / 'geneformer_covid_canonical.h5ad', 'X_geneformer'),
    'scGPT': (EMB / 'scgpt_covid_canonical.h5ad', 'X_scgpt'),
}

PMAP = pd.read_csv(PATIENT_MAP).set_index('sample_id')['patient_id'].to_dict()
JIANG_DATA_NO = '302004data18'  # Job J (2026-09-07): GSE151263 has no per-patient
# outcome record in the source paper -- "survivor" labels were study-design
# inference, not verifiable. Excluded entirely from the mortality task
# (7 donors / 476 cells). Legacy sepsis arm: 12 -> 5 patients.

# Job H2 Step 0 (2026-09-03, Task-4 audit vs Combes 2021 Supp Table 1 p.3):
# GSM4995432's severe->moderate flip matters HERE specifically -- under the
# severe-excluded sensitivity protocol this changes which set the cell falls
# into (not just the fatal/y computation, which only cares about 'fatal').
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


def load(fm):
    path, obsm_key = MODEL_PATHS[fm]
    with h5py.File(path, 'r') as h:
        X = h[f'obsm/{obsm_key}'][:]
        disease = read_categorical(h['obs/Disease'])
        severity = read_categorical(h['obs/covid_severity'])
        sample_id = read_categorical(h['obs/Sample_ID'])
        data_no = read_categorical(h['obs/Data_NO'])
    for sid, new_sev in SEVERITY_FIX.items():
        severity[sample_id == sid] = new_sev
    pid = np.array([PMAP.get(s, None) for s in sample_id], dtype=object)
    return X, disease, severity, sample_id, pid, data_no


def build_labels(disease, severity, covid_survived_set):
    is_covid = disease == 'COVID-19'
    is_sepsis = disease == 'Sepsis'
    y = np.zeros(len(disease), dtype=np.int8)
    valid = np.zeros(len(disease), dtype=bool)
    y[is_covid & np.isin(severity, list(COVID_FATAL))] = 1
    valid[is_covid & np.isin(severity, list(COVID_FATAL | covid_survived_set))] = True
    y[is_sepsis & np.isin(severity, list(SEPSIS_FATAL))] = 1
    valid[is_sepsis & np.isin(severity, list(SEPSIS_FATAL | SEPSIS_SURVIVED))] = True
    return y, valid, is_covid, is_sepsis


def run_protocol(X, disease, severity, sample_id, patient_id, data_no, covid_survived_set, direction):
    """One protocol variant (primary or sensitivity), one direction. Returns
    the transfer_one() dict, or {} if insufficient class diversity."""
    y, valid, cov_mask, sep_mask = build_labels(disease, severity, covid_survived_set)
    valid = valid & ~np.isin(sample_id, list(COMBES_DROP)) & (data_no != JIANG_DATA_NO)
    train_mask = cov_mask if direction.startswith('COVID') else sep_mask
    test_mask = sep_mask if direction.endswith('Sepsis') else cov_mask
    train_sel = train_mask & valid
    test_sel = test_mask & valid
    Xtr, ytr = X[train_sel], y[train_sel]
    Xte, yte = X[test_sel], y[test_sel]
    gte = patient_id[test_sel]
    log.info('    train n=%d (fatal=%d) test n=%d (fatal=%d, patients=%d)',
              len(ytr), int(ytr.sum()), len(yte), int(yte.sum()), len(np.unique(gte)))
    return transfer_one(Xtr, ytr, Xte, yte, gte)


def main():
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    rows = []
    full = {}
    for fm in MODEL_PATHS:
        log.info('=' * 70)
        log.info('Model: %s', fm)
        X, disease, severity, sample_id, patient_id, data_no = load(fm)

        _, valid_primary, _, _ = build_labels(disease, severity, COVID_SURVIVED_PRIMARY)
        valid_primary = valid_primary & ~np.isin(sample_id, list(COMBES_DROP)) & (data_no != JIANG_DATA_NO)
        n_missing_primary = int((pd.isna(patient_id) & valid_primary).sum())
        assert n_missing_primary == 0, f'{fm}: {n_missing_primary} valid cells lack patient_id mapping'
        patient_id = np.where(pd.isna(patient_id), 'NA', patient_id).astype(str)

        full[fm] = {}
        for direction in ['COVID->Sepsis', 'Sepsis->COVID']:
            log.info('Direction: %s', direction)
            log.info('  primary (severe=survived):')
            res_primary = run_protocol(X, disease, severity, sample_id, patient_id, data_no,
                                        COVID_SURVIVED_PRIMARY, direction)
            log.info('  sensitivity (severe-excluded):')
            res_sens = run_protocol(X, disease, severity, sample_id, patient_id, data_no,
                                     COVID_SURVIVED_SENS, direction)
            full[fm][direction] = {'primary': res_primary, 'sens_severe_excluded': res_sens}

            for clf in ['LogReg', 'RandomForest']:
                if clf not in res_primary or clf not in res_sens:
                    continue
                mp, ms = res_primary[clf], res_sens[clf]
                row = {'model': fm, 'direction': direction, 'classifier': clf,
                       'primary_auc': mp['auc'],
                       'primary_boot_auc_mean': mp['boot_auc_mean'],
                       'primary_boot_ci_lo': mp['boot_ci_lo'],
                       'primary_boot_ci_hi': mp['boot_ci_hi'],
                       'primary_auc_cell_weighted': mp['auc_cell_weighted'],
                       'primary_balanced_accuracy': mp['balanced_accuracy'],
                       'primary_n_train': res_primary['n_train'],
                       'primary_n_test': res_primary['n_test'],
                       'primary_fatal_train': res_primary['fatal_train'],
                       'primary_fatal_test': res_primary['fatal_test'],
                       'primary_n_test_patients': res_primary['n_test_patients'],
                       'primary_n_fatal_test_patients': res_primary['n_fatal_test_patients'],
                       'sens_severe_excluded_auc': ms['auc'],
                       'sens_severe_excluded_boot_auc_mean': ms['boot_auc_mean'],
                       'sens_severe_excluded_boot_ci_lo': ms['boot_ci_lo'],
                       'sens_severe_excluded_boot_ci_hi': ms['boot_ci_hi'],
                       'sens_severe_excluded_auc_cell_weighted': ms['auc_cell_weighted'],
                       'sens_severe_excluded_balanced_accuracy': ms['balanced_accuracy'],
                       'sens_severe_excluded_n_train': res_sens['n_train'],
                       'sens_severe_excluded_n_test': res_sens['n_test'],
                       'sens_severe_excluded_fatal_train': res_sens['fatal_train'],
                       'sens_severe_excluded_fatal_test': res_sens['fatal_test'],
                       'sens_severe_excluded_n_test_patients': res_sens['n_test_patients'],
                       'sens_severe_excluded_n_fatal_test_patients': res_sens['n_fatal_test_patients']}
                rows.append(row)
                log.info('  %s: primary AUC=%.4f [%.4f,%.4f] | sens(severe-excl) AUC=%.4f [%.4f,%.4f]',
                         clf, mp['auc'], mp['boot_ci_lo'], mp['boot_ci_hi'],
                         ms['auc'], ms['boot_ci_lo'], ms['boot_ci_hi'])

    df = pd.DataFrame(rows)
    RESULTS.mkdir(parents=True, exist_ok=True)
    csv_path = RESULTS / f'cross_disease_patient_jobJ_{timestamp}.csv'
    df.to_csv(csv_path, index=False)
    (RESULTS / f'cross_disease_patient_jobJ_{timestamp}.json').write_text(json.dumps(full, indent=2))
    log.info('Saved: %s', csv_path)


if __name__ == '__main__':
    main()
