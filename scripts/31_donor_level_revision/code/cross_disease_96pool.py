#!/usr/bin/env python3
"""
Job J-supplement (HPCC_next_steps.md §11): reverse-transfer 96-patient pool
completion. Combines the legacy Qiu sepsis arm (5 patients, GSE167363, Jiang
excluded) with the Chongqing sepsis arm (91 patients, GSE279452) into one
96-patient sepsis-side pool for cross-disease transfer, matching the protocol
of jobS_cross_disease.py (primary/severe-excluded, LR+RF, 200-iter patient
bootstrap) and cross_disease_fatal_vs_survived_patient_v2.py (Jiang exclusion,
patient-level Fix-2/Fix-3 scoring via transfer_utils).

Output: 05_RESULTS/revision_donor_level_v2/cross_disease_patient_96pool_<ts>.csv/.json
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
log = logging.getLogger('cross_disease_96pool')

BASE = Path('/bigdata/godziklab/shared/Xinru/302006')
COVID_EMB = BASE / '02_EMBEDDINGS' / 'option_c'
CHONGQING_EMB = BASE / '02_EMBEDDINGS' / 'option_c_chongqing'
RESULTS = BASE / '05_RESULTS' / 'revision_donor_level_v2'
PATIENT_MAP = (BASE / '302006_updategit' / 'scripts' / '31_donor_level_revision' /
               'results' / 'donor_meta_with_patient_id_v2.csv')

COMBES_DROP = {'GSM4995425', 'GSM4995426', 'GSM4995430'}
JIANG_DATA_NO = '302004data18'
SEVERITY_FIX = {
    'GSM4995436': 'moderate',
    'GSM4995442': 'severe',
    'GSM4995432': 'moderate',
}

MODEL_PATHS = {
    'STATE': (COVID_EMB / 'state_covid_canonical.h5ad', 'X_state',
              CHONGQING_EMB / 'state_chongqing_canonical.h5ad', 'X_state'),
    'UCE': (COVID_EMB / 'uce_covid_canonical.h5ad', 'X_uce',
            CHONGQING_EMB / 'uce_chongqing_canonical.h5ad', 'X_uce'),
    'scPRINT-2': (COVID_EMB / 'scprint2_covid_canonical.h5ad', 'X_scprint2',
                  CHONGQING_EMB / 'scprint2_chongqing_canonical.h5ad', 'X_scprint2'),
    'Geneformer': (COVID_EMB / 'geneformer_covid_canonical.h5ad', 'X_geneformer',
                   CHONGQING_EMB / 'geneformer_chongqing_canonical.h5ad', 'X_geneformer'),
    'scGPT': (COVID_EMB / 'scgpt_covid_canonical.h5ad', 'X_scgpt',
              CHONGQING_EMB / 'scgpt_chongqing_canonical.h5ad', 'X_scgpt'),
}

PMAP = pd.read_csv(PATIENT_MAP).set_index('sample_id')['patient_id'].to_dict()


def read_categorical(group):
    if isinstance(group, h5py.Group) and 'categories' in group:
        cats = group['categories'][:].astype(str)
        codes = group['codes'][:]
        return cats[codes]
    return group[:].astype(str)


def load_covid_side(path, obsm_key, covid_survived_set):
    with h5py.File(path, 'r') as h:
        X = h[f'obsm/{obsm_key}'][:]
        disease = read_categorical(h['obs/Disease'])
        severity = read_categorical(h['obs/covid_severity'])
        sample_id = read_categorical(h['obs/Sample_ID'])
    for sid, new_sev in SEVERITY_FIX.items():
        severity[sample_id == sid] = new_sev
    is_covid = disease == 'COVID-19'
    y = np.zeros(len(disease), dtype=np.int8)
    valid = np.zeros(len(disease), dtype=bool)
    y[is_covid & np.isin(severity, list(COVID_FATAL))] = 1
    valid[is_covid & np.isin(severity, list(COVID_FATAL | covid_survived_set))] = True
    valid = valid & ~np.isin(sample_id, list(COMBES_DROP))
    pid = np.array([PMAP.get(s, None) for s in sample_id], dtype=object)
    n_missing = int(pd.isna(pid[valid]).sum())
    assert n_missing == 0, f'{n_missing} valid COVID cells lack patient_id mapping'
    pid = np.where(pd.isna(pid), 'NA', pid).astype(str)
    return X[valid], y[valid], pid[valid]


def load_legacy_qiu_sepsis(path, obsm_key):
    """Legacy sepsis side from the same canonical COVID h5ad: Disease=='Sepsis',
    Jiang (302004data18) excluded -> Qiu (GSE167363) only, 5 patients."""
    with h5py.File(path, 'r') as h:
        X = h[f'obsm/{obsm_key}'][:]
        disease = read_categorical(h['obs/Disease'])
        severity = read_categorical(h['obs/covid_severity'])
        sample_id = read_categorical(h['obs/Sample_ID'])
        data_no = read_categorical(h['obs/Data_NO'])
    is_sepsis = disease == 'Sepsis'
    valid = is_sepsis & np.isin(severity, list(SEPSIS_FATAL | SEPSIS_SURVIVED)) & (data_no != JIANG_DATA_NO)
    y = np.zeros(len(severity), dtype=np.int8)
    y[np.isin(severity, list(SEPSIS_FATAL))] = 1
    pid = np.array([PMAP.get(s, None) for s in sample_id], dtype=object)
    n_missing = int(pd.isna(pid[valid]).sum())
    assert n_missing == 0, f'{n_missing} valid legacy sepsis cells lack patient_id mapping'
    pid = np.where(pd.isna(pid), 'NA', pid).astype(str)
    return X[valid], y[valid], pid[valid]


def load_chongqing_sepsis(path, obsm_key):
    with h5py.File(path, 'r') as h:
        X = h[f'obsm/{obsm_key}'][:]
        severity = read_categorical(h['obs/covid_severity'])
        sample_id = read_categorical(h['obs/Sample_ID'])
    y = np.zeros(len(severity), dtype=np.int8)
    valid = np.isin(severity, list(SEPSIS_FATAL | SEPSIS_SURVIVED))
    y[np.isin(severity, list(SEPSIS_FATAL))] = 1
    return X[valid], y[valid], sample_id[valid]  # Chongqing: sample_id == patient_id


def main():
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    rows = []
    full = {}
    for fm, (covid_path, covid_key, chong_path, chong_key) in MODEL_PATHS.items():
        log.info('=' * 70)
        log.info('Model: %s', fm)
        if not chong_path.exists():
            log.warning('  %s missing -- skipping', chong_path)
            continue

        Xq, yq, gq = load_legacy_qiu_sepsis(covid_path, covid_key)
        Xc, yc_c, gc_c = load_chongqing_sepsis(chong_path, chong_key)
        Xs = np.concatenate([Xq, Xc], axis=0)
        ys = np.concatenate([yq, yc_c], axis=0)
        gs = np.concatenate([gq, gc_c], axis=0)
        n_pat_q = len(np.unique(gq))
        n_pat_c = len(np.unique(gc_c))
        n_pat_s = len(np.unique(gs))
        log.info('  Sepsis pool: Qiu %d cells/%d patients + Chongqing %d cells/%d patients '
                 '= %d cells/%d patients (expect 96)',
                 len(yq), n_pat_q, len(yc_c), n_pat_c, len(ys), n_pat_s)
        assert n_pat_s == n_pat_q + n_pat_c, 'patient id collision between Qiu and Chongqing pools'

        full[fm] = {}
        for direction in ['COVID->Sepsis', 'Sepsis->COVID']:
            log.info('Direction: %s', direction)
            proto_results = {}
            for proto_name, covid_set in [('primary', COVID_SURVIVED_PRIMARY),
                                           ('sens_severe_excluded', COVID_SURVIVED_SENS)]:
                Xcv, ycv, gcv = load_covid_side(covid_path, covid_key, covid_set)
                log.info('  [%s] COVID valid: %d cells, %d patients', proto_name,
                          len(ycv), len(np.unique(gcv)))
                if direction == 'COVID->Sepsis':
                    res = transfer_one(Xcv, ycv, Xs, ys, gs)
                else:
                    res = transfer_one(Xs, ys, Xcv, ycv, gcv)
                proto_results[proto_name] = res
            full[fm][direction] = proto_results

            res_primary = proto_results['primary']
            res_sens = proto_results['sens_severe_excluded']
            for clf in ['LogReg', 'RandomForest']:
                if clf not in res_primary or clf not in res_sens:
                    continue
                mp, ms = res_primary[clf], res_sens[clf]
                row = {'model': fm, 'direction': direction, 'classifier': clf,
                       'sepsis_pool': '96pool_qiu5_chongqing91',
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
    csv_path = RESULTS / f'cross_disease_patient_96pool_{timestamp}.csv'
    df.to_csv(csv_path, index=False)
    (RESULTS / f'cross_disease_patient_96pool_{timestamp}.json').write_text(json.dumps(full, indent=2))
    log.info('Saved: %s', csv_path)


if __name__ == '__main__':
    main()
