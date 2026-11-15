#!/usr/bin/env python3
"""
BIB revision memo item B (2026-09-13 memo, Tier 3): baseline-ladder headline
AUC table -- reconstructs the canonical within-COVID (240 patients, 25
fatal) repeat-mean and ensemble-OOF patient-level AUC for the new
from-scratch baselines (HVG-mean, MarkerScore, PseudobulkEN) using the exact
same construction as Table 3 (main.tex) and clinical_utility_and_consensus.py:
per seed, aggregate cells to patients (mean p_lr); ensemble-OOF = mean of
patient scores across the 10 seeds, scored once; repeat-mean = mean of the
10 per-seed patient-level AUCs. Restricted to the canonical COVID-only
cohort (the raw seed files pool COVID+sepsis, 252 patients -- see
clinical_utility_and_consensus.py's docstring for how this was caught).

PseudobulkEN's oof files are already one row per patient (no cell_id), a
different schema from the cell-level HVG-mean/MarkerScore/FM/PCA-200 files
-- handled separately below. Existing methods (the 5 FMs + PCA-200) are not
recomputed here; their Table 3 values are hard-coded for the ladder table
context (verified against clinical_utility_and_consensus.py's sanity check,
which passed to within 0.0005).

Gracefully skips any method whose seed files are not yet present (so this
can be re-run as new baselines finish landing) rather than failing outright.

Output: baseline_ladder_auc.csv (this dir, OUT_DIR).
"""
import glob
import os

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

IN_DIR = ('/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/'
          'revision_donor_level_v2')
PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')
PATIENT_LIST_CSV = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
                     '31_donor_level_revision/results_patient_v2/qc_baseline/'
                     'patient_qc_distributions.csv')
OUT_DIR = ('/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/'
           'bib_revision_tier2/clinical_utility')
os.makedirs(OUT_DIR, exist_ok=True)

SEEDS = list(range(42, 52))

_patient_list = pd.read_csv(PATIENT_LIST_CSV, dtype={'patient_id': str})
_patient_list = _patient_list[_patient_list['disease'] == 'COVID-19']
assert len(_patient_list) == 240 and int(_patient_list['y'].sum()) == 25
CANONICAL_COVID_PATIENTS = set(_patient_list['patient_id'])

_pmap = pd.read_csv(PATIENT_MAP, dtype={'sample_id': str, 'patient_id': str})[['sample_id', 'patient_id']]

# Existing methods, already verified against Table 3 (clinical_utility_and_consensus.py
# sanity check, max |diff| = 0.0005) -- included here for one complete ladder table.
EXISTING_TABLE3 = {
    'Geneformer': {'repeat_mean': 0.786, 'repeat_sd': 0.026, 'ensemble_oof': 0.807},
    'STATE': {'repeat_mean': 0.778, 'repeat_sd': 0.033, 'ensemble_oof': 0.795},
    'scPRINT-2': {'repeat_mean': 0.740, 'repeat_sd': 0.027, 'ensemble_oof': 0.756},
    'PCA-200': {'repeat_mean': 0.725, 'repeat_sd': 0.033, 'ensemble_oof': 0.742},
    'scGPT': {'repeat_mean': 0.720, 'repeat_sd': 0.047, 'ensemble_oof': 0.739},
    'UCE': {'repeat_mean': 0.661, 'repeat_sd': 0.046, 'ensemble_oof': 0.677},
}


def cell_level_scores(method_name, file_prefix):
    per_seed_patient_scores = []
    per_seed_auc = []
    y_ref, pid_ref = None, None
    for seed in SEEDS:
        path = os.path.join(IN_DIR, f'{file_prefix}_seed{seed}.csv')
        if not os.path.exists(path):
            return None
        df = pd.read_csv(path, dtype={'patient_id': str, 'sample_id': str})
        if 'patient_id' not in df.columns:
            df = df.merge(_pmap, on='sample_id', how='inner')
        df['patient_id'] = df['patient_id'].astype(str)
        df = df[df['patient_id'].isin(CANONICAL_COVID_PATIENTS)]
        pat = df.groupby('patient_id').agg(y=('y', 'max'), score=('p_lr', 'mean')).sort_index()
        if y_ref is None:
            y_ref, pid_ref = pat['y'].to_numpy(), pat.index.to_numpy()
        per_seed_patient_scores.append(pat['score'].to_numpy())
        per_seed_auc.append(roc_auc_score(pat['y'].to_numpy(), pat['score'].to_numpy()))
    ensemble_score = np.mean(per_seed_patient_scores, axis=0)
    return {
        'n_patients': len(pid_ref), 'n_fatal': int(y_ref.sum()),
        'repeat_mean': float(np.mean(per_seed_auc)), 'repeat_sd': float(np.std(per_seed_auc, ddof=1)),
        'ensemble_oof': float(roc_auc_score(y_ref, ensemble_score)),
    }


def patient_level_scores(method_name, file_prefix):
    """For PseudobulkEN-style files: already one row per patient, no cell aggregation."""
    per_seed_patient_scores = []
    per_seed_auc = []
    y_ref, pid_ref = None, None
    for seed in SEEDS:
        path = os.path.join(IN_DIR, f'{file_prefix}_seed{seed}.csv')
        if not os.path.exists(path):
            return None
        df = pd.read_csv(path, dtype={'patient_id': str})
        df = df[df['patient_id'].isin(CANONICAL_COVID_PATIENTS)].sort_values('patient_id')
        if y_ref is None:
            y_ref, pid_ref = df['y'].to_numpy(), df['patient_id'].to_numpy()
        per_seed_patient_scores.append(df['p_lr'].to_numpy())
        per_seed_auc.append(roc_auc_score(df['y'].to_numpy(), df['p_lr'].to_numpy()))
    ensemble_score = np.mean(per_seed_patient_scores, axis=0)
    return {
        'n_patients': len(pid_ref), 'n_fatal': int(y_ref.sum()),
        'repeat_mean': float(np.mean(per_seed_auc)), 'repeat_sd': float(np.std(per_seed_auc, ddof=1)),
        'ensemble_oof': float(roc_auc_score(y_ref, ensemble_score)),
    }


rows = []
for m, d in EXISTING_TABLE3.items():
    rows.append({'method': m, 'source': 'Table 3 (verified)', 'n_patients': 240, 'n_fatal': 25,
                 **d})

new_cell_level = {'HVG-mean': 'oof_patient_HVG-mean', 'MarkerScore': 'oof_patient_MarkerScore'}
for m, prefix in new_cell_level.items():
    result = cell_level_scores(m, prefix)
    if result is None:
        print(f'{m}: seed files not found yet, skipping', flush=True)
        continue
    rows.append({'method': m, 'source': 'this script (cell-level agg)', **result})
    print(f'{m}: repeat-mean={result["repeat_mean"]:.4f} +/- {result["repeat_sd"]:.4f}, '
          f'ensemble-OOF={result["ensemble_oof"]:.4f}', flush=True)

result = patient_level_scores('PseudobulkEN', 'oof_patient_PseudobulkEN')
if result is None:
    print('PseudobulkEN: seed files not found yet, skipping', flush=True)
else:
    rows.append({'method': 'PseudobulkEN', 'source': 'this script (patient-level)', **result})
    print(f'PseudobulkEN: repeat-mean={result["repeat_mean"]:.4f} +/- {result["repeat_sd"]:.4f}, '
          f'ensemble-OOF={result["ensemble_oof"]:.4f}', flush=True)

out_df = pd.DataFrame(rows).sort_values('ensemble_oof', ascending=False)
out_path = os.path.join(OUT_DIR, 'baseline_ladder_auc.csv')
out_df.to_csv(out_path, index=False)
print(f'\nSAVED {out_path}', flush=True)
print(out_df.to_string(index=False), flush=True)
