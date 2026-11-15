#!/usr/bin/env python3
"""
Follow-up to Job O: WHY does the exact patient-matched comparison flip the
official mono<platelet inversion for STATE/UCE/scGPT/PCA-200?

Three candidate explanations, tested one at a time by toggling each in turn
while holding the others fixed at the "official-comparison-like" setting:
  (1) single seed-42 split (noisy) vs. 10-seed ensemble average (matches
      how the official repeated-CV numbers are built)
  (2) all-disease patient pool (148, incl. 5 Qiu-sepsis) vs. COVID-only
      pool (143, Ren+Stephenson only, matching the official comparison's
      COVID-only framing)
  (3) implicitly tested by comparing against the ALREADY-COMPUTED official
      indistribution_vs_platelet.csv numbers, which use the full unmatched
      platelet n=240 and full unmatched monocyte n=148/143 -- if the
      ensemble+COVID-only matched version converges back to the official
      direction, the flip was (1)+(2) i.e. noise+composition, not patient
      selection per se. If it does NOT converge, something about the
      matched *subset itself* (which patients got excluded to build the
      match) is the real driver.
"""
import os
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')
MONO_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2/monocyte'
PLATELET_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2'
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2/monocyte_matched'

FMS = ['STATE', 'UCE', 'Geneformer', 'scPRINT-2', 'scGPT', 'PCA-200']
SEEDS = list(range(42, 52))
RESTRICT_DATA_NO = ['302004data15', '302004data17']  # Ren + Stephenson only (COVID-only)


def safe_auc(y, p):
    if len(np.unique(y)) < 2:
        return np.nan
    return roc_auc_score(y, p)


def ensemble_scores(fm, side, intersection, pmap=None):
    """side='mono' or 'plt'. Returns per-patient (y, score) averaged over 10 seeds."""
    frames = []
    for seed in SEEDS:
        if side == 'mono':
            path = os.path.join(MONO_DIR, f'oof_monocyte_{fm}_seed{seed}.csv')
        else:
            if fm == 'PCA-200':
                path = os.path.join(PLATELET_DIR, f'oof_PCA-200_seed{seed}.csv')
            else:
                path = os.path.join(PLATELET_DIR, f'oof_patient_{fm}_seed{seed}.csv')
        if not os.path.exists(path):
            return None
        df = pd.read_csv(path)
        if side == 'plt' and 'patient_id' not in df.columns:
            df = df.merge(pmap[['sample_id', 'patient_id']], on='sample_id', how='inner')
        df = df[df['patient_id'].isin(intersection)]
        pat = df.groupby('patient_id').agg(y=('y', 'max'), score=('p_lr', 'mean'))
        pat = pat.rename(columns={'score': f'score_{seed}'})
        frames.append(pat)
    if not frames:
        return None
    merged = pd.concat([f[[c for c in f.columns if c.startswith('score_')]] for f in frames], axis=1)
    y = frames[0]['y']
    ens_score = merged.mean(axis=1)
    return pd.DataFrame({'y': y, 'score': ens_score})


def main():
    pmap_full = pd.read_csv(PATIENT_MAP)
    mono_any = pd.read_csv(os.path.join(MONO_DIR, 'oof_monocyte_STATE_seed42.csv'))
    mono_patient_ids = set(mono_any['patient_id'].unique())
    platelet_all = set(pmap_full[pmap_full['Data_NO'].isin(
        ['302004data15', '302004data17', '302005data03'])]['patient_id'])
    intersection_all = platelet_all & mono_patient_ids

    platelet_covid_only = set(pmap_full[pmap_full['Data_NO'].isin(RESTRICT_DATA_NO)]['patient_id'])
    intersection_covidonly = platelet_covid_only & mono_patient_ids

    print(f'Intersection (all-disease, n=148 expected): {len(intersection_all)}')
    print(f'Intersection (COVID-only Ren+Stephenson): {len(intersection_covidonly)}')

    configs = [
        ('seed42_alldisease', 'seed42', intersection_all),
        ('ensemble_alldisease', 'ensemble', intersection_all),
        ('seed42_covidonly', 'seed42', intersection_covidonly),
        ('ensemble_covidonly', 'ensemble', intersection_covidonly),
    ]

    rows = []
    for cfg_name, mode, patients in configs:
        for fm in FMS:
            if mode == 'seed42':
                mono_path = os.path.join(MONO_DIR, f'oof_monocyte_{fm}_seed42.csv')
                plt_path = (os.path.join(PLATELET_DIR, 'oof_PCA-200_seed42.csv') if fm == 'PCA-200'
                            else os.path.join(PLATELET_DIR, f'oof_patient_{fm}_seed42.csv'))
                if not os.path.exists(mono_path) or not os.path.exists(plt_path):
                    continue
                mdf = pd.read_csv(mono_path)
                mdf = mdf[mdf['patient_id'].isin(patients)]
                mono_pat = mdf.groupby('patient_id').agg(y=('y', 'max'), score=('p_lr', 'mean'))
                pdf = pd.read_csv(plt_path)
                if 'patient_id' not in pdf.columns:
                    pdf = pdf.merge(pmap_full[['sample_id', 'patient_id']], on='sample_id', how='inner')
                pdf = pdf[pdf['patient_id'].isin(patients)]
                plt_pat = pdf.groupby('patient_id').agg(y=('y', 'max'), score=('p_lr', 'mean'))
            else:
                mono_pat = ensemble_scores(fm, 'mono', patients)
                plt_pat = ensemble_scores(fm, 'plt', patients, pmap=pmap_full)
                if mono_pat is None or plt_pat is None:
                    continue

            common = mono_pat.join(plt_pat, how='inner', lsuffix='_mono', rsuffix='_plt')
            if len(common) < 5 or len(np.unique(common['y_mono'])) < 2:
                continue
            auc_mono = safe_auc(common['y_mono'].to_numpy(), common['score_mono'].to_numpy())
            auc_plt = safe_auc(common['y_plt'].to_numpy(), common['score_plt'].to_numpy())
            rows.append({'config': cfg_name, 'fm': fm, 'n_matched': len(common),
                         'auc_monocyte': round(auc_mono, 4), 'auc_platelet': round(auc_plt, 4),
                         'delta_mono_minus_plt': round(auc_mono - auc_plt, 4)})

    result = pd.DataFrame(rows)
    result.to_csv(os.path.join(OUT_DIR, 'flip_investigation.csv'), index=False)
    print('\n=== Flip investigation: delta = monocyte AUC - platelet AUC (matched patients) ===')
    piv = result.pivot(index='fm', columns='config', values='delta_mono_minus_plt')
    piv = piv[['seed42_alldisease', 'ensemble_alldisease', 'seed42_covidonly', 'ensemble_covidonly']]
    print(piv.to_string())
    print('\n(positive = monocyte higher/flipped; negative = platelet higher/matches official direction)')
    print(f'\nSAVED -> {OUT_DIR}/flip_investigation.csv')


if __name__ == '__main__':
    main()
