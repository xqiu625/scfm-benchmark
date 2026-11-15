#!/usr/bin/env python3
"""
BIB revision memo items C and F (2026-09-13 memo, Tier 2): clinical-utility
layer (PR-AUC / average precision, decision-curve net benefit) and the
consensus/stacked predictor -- both computed directly from the already-saved
canonical within-COVID out-of-fold patient predictions
(results_patient_v2/within_covid/donor_predictions_<Model>.csv), which
already carry the repeat-mean-aggregated 'mean' score per patient (the same
primary quantity reported in main.tex Table 3). No retraining, no new
embeddings -- this is a downstream re-analysis of existing saved scores.

CORRECTED 2026-09-13: the first version of this script used
results_patient_v2/within_covid/donor_predictions_<Model>.csv's 'mean'
column, which turns out to be the SEED-42-ONLY single-split OOF score (it
matches per_seed_auc_repeated.csv's seed=42 row exactly, e.g. STATE 0.7113),
NOT the 10-seed ensemble-OOF score the paper actually reports (STATE 0.795
ensemble-OOF / 0.778 repeat-mean, per Table 3) -- a ~0.04-0.12 AUC gap for
STATE/UCE/scGPT/PCA-200 gave this away immediately when cross-checked
against Table 3. Fixed by reconstructing the true ensemble-OOF score from
the raw per-seed, per-cell files (05_RESULTS/revision_donor_level_v2/
oof_patient_<Model>_seed<42..51>.csv): aggregate cells to patients (mean
p_lr) within each seed, then average that patient-level score across the
10 seeds -- exactly Methods' definition ("averages each patient's scores
across the ten repeats and scores them once").

Outputs (this dir, OUT_DIR):
  clinical_utility.csv   -- per-method: AUC, AP (average precision), and net
                            benefit at a grid of threshold probabilities,
                            compared against treat-all / treat-none.
  consensus_predictor.csv -- AUC/AP for two consensus scores (mean of the 5
                            FMs; mean of the 5 FMs + PCA-200) alongside every
                            individual method, so the table drops straight
                            into the paper's F item ("does averaging beat any
                            single model?").
  ensemble_auc_sanity_check.csv -- reconstructed ensemble-OOF AUC per method
                            vs. Table 3's reported value, to catch exactly
                            this kind of data-source mistake before trusting
                            the rest of the output.
"""
import os

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

IN_DIR = ('/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/'
          'revision_donor_level_v2')
PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')
PATIENT_LIST_CSV = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
                     '31_donor_level_revision/results_patient_v2/qc_baseline/'
                     'patient_qc_distributions.csv')
_patient_list = pd.read_csv(PATIENT_LIST_CSV)
_patient_list = _patient_list[_patient_list['disease'] == 'COVID-19']
assert len(_patient_list) == 240 and int(_patient_list['y'].sum()) == 25, \
    f'canonical cohort mismatch: n={len(_patient_list)}, fatal={int(_patient_list["y"].sum())}'
CANONICAL_COVID_PATIENTS = set(_patient_list['patient_id'].astype(str))
OUT_DIR = ('/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/'
           'bib_revision_tier2/clinical_utility')
os.makedirs(OUT_DIR, exist_ok=True)

METHODS = ['Geneformer', 'STATE', 'UCE', 'scGPT', 'scPRINT-2', 'PCA-200']
SEEDS = list(range(42, 52))
THRESHOLDS = np.round(np.arange(0.05, 0.51, 0.05), 2)  # 5%-50% risk thresholds

# Table 3 (main.tex) ensemble-OOF AUC, for the sanity check
TABLE3_ENSEMBLE_AUC = {'Geneformer': 0.807, 'STATE': 0.795, 'scPRINT-2': 0.756,
                        'PCA-200': 0.742, 'scGPT': 0.739, 'UCE': 0.677}


_pmap = pd.read_csv(PATIENT_MAP, dtype={'sample_id': str, 'patient_id': str})[['sample_id', 'patient_id']]


def load_ensemble_oof(method):
    """Reconstruct the ensemble-OOF patient-level score: per seed, aggregate
    cells to patients (mean p_lr); then average across the 10 seeds.
    PCA-200's oof files only carry sample_id (no patient_id column, unlike
    the 5 FM files) -- merge with the patient map for that one method."""
    fname_prefix = 'oof_PCA-200' if method == 'PCA-200' else f'oof_patient_{method}'
    per_seed_patient_scores = []
    y_ref, pid_ref = None, None
    for seed in SEEDS:
        path = os.path.join(IN_DIR, f'{fname_prefix}_seed{seed}.csv')
        df = pd.read_csv(path, dtype={'patient_id': str, 'sample_id': str})
        if 'patient_id' not in df.columns:
            df = df.merge(_pmap, on='sample_id', how='inner')
        df['patient_id'] = df['patient_id'].astype(str)
        df = df[df['patient_id'].isin(CANONICAL_COVID_PATIENTS)]
        pat = df.groupby('patient_id').agg(y=('y', 'max'), score=('p_lr', 'mean'))
        pat = pat.sort_index()
        if y_ref is None:
            y_ref, pid_ref = pat['y'].to_numpy(), pat.index.to_numpy()
        else:
            assert np.array_equal(pat['y'].to_numpy(), y_ref), f'{method} seed{seed}: label mismatch'
            assert np.array_equal(pat.index.to_numpy(), pid_ref), f'{method} seed{seed}: patient set mismatch'
        per_seed_patient_scores.append(pat['score'].to_numpy())
    ensemble_score = np.mean(per_seed_patient_scores, axis=0)
    return pid_ref, y_ref, ensemble_score


def net_benefit(y, p, pt):
    """Standard decision-curve net benefit at threshold probability pt."""
    n = len(y)
    pred_pos = p >= pt
    tp = int(np.sum(pred_pos & (y == 1)))
    fp = int(np.sum(pred_pos & (y == 0)))
    return (tp / n) - (fp / n) * (pt / (1 - pt))


def treat_all_net_benefit(y, pt):
    n = len(y)
    tp = int(np.sum(y == 1))
    fp = int(np.sum(y == 0))
    return (tp / n) - (fp / n) * (pt / (1 - pt))


scores = {}
labels = None
sanity_rows = []
for m in METHODS:
    pid, y, s = load_ensemble_oof(m)
    if labels is None:
        labels = y
        patient_ids = pid
    else:
        assert np.array_equal(y, labels), f'{m}: label mismatch vs first method'
        assert np.array_equal(pid, patient_ids), f'{m}: patient_id order mismatch'
    scores[m] = s
    recon_auc = roc_auc_score(y, s)
    sanity_rows.append({'method': m, 'reconstructed_ensemble_auc': recon_auc,
                         'table3_ensemble_auc': TABLE3_ENSEMBLE_AUC[m],
                         'diff': recon_auc - TABLE3_ENSEMBLE_AUC[m]})

print(f'{len(labels)} patients, {int(labels.sum())} fatal', flush=True)

sanity_df = pd.DataFrame(sanity_rows)
sanity_path = os.path.join(OUT_DIR, 'ensemble_auc_sanity_check.csv')
sanity_df.to_csv(sanity_path, index=False)
print(f'SAVED {sanity_path}', flush=True)
print(sanity_df.to_string(index=False), flush=True)
max_diff = sanity_df['diff'].abs().max()
if max_diff > 0.01:
    raise AssertionError(
        f'Reconstructed ensemble-OOF AUC diverges from Table 3 by up to {max_diff:.4f} -- '
        f'do not trust downstream clinical-utility/consensus numbers until this is resolved.')
print(f'Sanity check passed: max |diff| vs Table 3 = {max_diff:.4f}', flush=True)

# ---------------------------------------------------------------------
# C: clinical utility -- AUC, AP, and net benefit at a threshold grid
# ---------------------------------------------------------------------
rows = []
for m in METHODS:
    p = scores[m]
    auc = roc_auc_score(labels, p)
    ap = average_precision_score(labels, p)
    row = {'method': m, 'auc': auc, 'average_precision': ap}
    for pt in THRESHOLDS:
        row[f'net_benefit_pt{pt}'] = net_benefit(labels, p, pt)
    rows.append(row)

# reference strategies (same for every method, listed once each)
for pt in THRESHOLDS:
    pass  # computed below into a separate reference row block
ref_all = {'method': 'treat_all', 'auc': np.nan, 'average_precision': np.nan}
ref_none = {'method': 'treat_none', 'auc': np.nan, 'average_precision': np.nan}
for pt in THRESHOLDS:
    ref_all[f'net_benefit_pt{pt}'] = treat_all_net_benefit(labels, pt)
    ref_none[f'net_benefit_pt{pt}'] = 0.0
rows.append(ref_all)
rows.append(ref_none)

util_df = pd.DataFrame(rows)
util_path = os.path.join(OUT_DIR, 'clinical_utility.csv')
util_df.to_csv(util_path, index=False)
print(f'SAVED {util_path}', flush=True)
print(util_df.to_string(index=False), flush=True)

# ---------------------------------------------------------------------
# F: consensus / stacked predictor
# ---------------------------------------------------------------------
fm_names = ['Geneformer', 'STATE', 'UCE', 'scGPT', 'scPRINT-2']
consensus_fm = np.mean([scores[m] for m in fm_names], axis=0)
consensus_fm_pca = np.mean([scores[m] for m in fm_names + ['PCA-200']], axis=0)

consensus_rows = []
for m in METHODS:
    consensus_rows.append({'method': m, 'auc': roc_auc_score(labels, scores[m]),
                            'average_precision': average_precision_score(labels, scores[m])})
consensus_rows.append({'method': 'consensus_5FM_mean', 'auc': roc_auc_score(labels, consensus_fm),
                        'average_precision': average_precision_score(labels, consensus_fm)})
consensus_rows.append({'method': 'consensus_5FM_plus_PCA200_mean', 'auc': roc_auc_score(labels, consensus_fm_pca),
                        'average_precision': average_precision_score(labels, consensus_fm_pca)})

consensus_df = pd.DataFrame(consensus_rows).sort_values('auc', ascending=False)
consensus_path = os.path.join(OUT_DIR, 'consensus_predictor.csv')
consensus_df.to_csv(consensus_path, index=False)
print(f'SAVED {consensus_path}', flush=True)
print(consensus_df.to_string(index=False), flush=True)
