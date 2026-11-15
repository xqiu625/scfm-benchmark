#!/usr/bin/env python3
"""
Job O -- monocyte patient-matched comparison (external review C13).

Correction to an earlier (wrong) assessment: platelet-side and monocyte-side
patient_ids for Ren (Data_NO=302004data15, e.g. S-M026/S-S070) and Stephenson
(Data_NO=302004data17, e.g. MH9143327) and Qiu sepsis (Data_NO=302005data03,
e.g. NS1/S1) are IDENTICAL strings on both sides -- no crosswalk needed. The
earlier "no crosswalk exists" conclusion mistakenly assumed 302005data01
(actually Bernardes, numeric ids like "001") was Ren's platelet Data_NO.

Take the intersection of patients with BOTH platelet (within-COVID/sepsis v2)
and monocyte data; cell-count-match per patient (downsample larger side to
the smaller cell count); recompute monocyte vs platelet AUC on this SAME
matched patient set (all 5 monocyte FMs + PCA-200; platelet side uses the
already-computed within-COVID/within-sepsis per-patient scores, seed-42 view
for a fixed single-split reference matched to the design of this comparison).
"""
import os
import numpy as np
import pandas as pd

PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')
MONO_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2/monocyte'
PLATELET_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2'
OUT_DIR = os.path.join(MONO_DIR, '..', 'monocyte_matched')

FMS = ['STATE', 'UCE', 'Geneformer', 'scPRINT-2', 'scGPT']
SHARED_DATA_NO = ['302004data15', '302004data17', '302005data03']  # Ren, Stephenson, Qiu-sepsis


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    pmap = pd.read_csv(PATIENT_MAP)
    platelet_shared = pmap[pmap['Data_NO'].isin(SHARED_DATA_NO)]
    platelet_patient_ids = set(platelet_shared['patient_id'])
    print(f'Platelet-side patients in shared-source cohorts (Ren+Stephenson+Qiu): '
          f'{len(platelet_patient_ids)}', flush=True)

    mono_any = pd.read_csv(os.path.join(MONO_DIR, 'oof_monocyte_STATE_seed42.csv'))
    mono_patient_ids = set(mono_any['patient_id'].unique())
    print(f'Monocyte-side total patients: {len(mono_patient_ids)}', flush=True)

    intersection = platelet_patient_ids & mono_patient_ids
    print(f'INTERSECTION (same patient, both platelet and monocyte data): {len(intersection)}',
          flush=True)

    inter_df = mono_any[mono_any['patient_id'].isin(intersection)][
        ['patient_id', 'source_study', 'y']].drop_duplicates(subset='patient_id')
    print('\nIntersection by source, with fatal counts:')
    print(inter_df.groupby('source_study').agg(n=('patient_id', 'count'), fatal=('y', 'sum'))
          .to_string())
    inter_df.to_csv(os.path.join(OUT_DIR, 'intersection_patients.csv'), index=False)

    if len(intersection) < 10:
        print('\nIntersection too small (<10) for a meaningful matched comparison -- '
              'reporting as a limitation, not forcing a result.', flush=True)
        return

    # ---- cell-count matching: get each patient's cell count on each side ----
    mono_counts = mono_any[mono_any['patient_id'].isin(intersection)].groupby(
        'patient_id').size().rename('n_cells_monocyte')

    plt_cell_counts = platelet_shared[platelet_shared['patient_id'].isin(intersection)].groupby(
        'patient_id')['n_cells'].sum().rename('n_cells_platelet')

    match_df = pd.concat([mono_counts, plt_cell_counts], axis=1).dropna()
    match_df['n_cells_matched'] = match_df[['n_cells_monocyte', 'n_cells_platelet']].min(axis=1)
    match_df.to_csv(os.path.join(OUT_DIR, 'cell_count_matching.csv'))
    print('\nCell-count matching (min per patient used for downsampled AUC):')
    print(match_df.describe().to_string())

    # ---- AUC on the SAME matched patient set, monocyte vs platelet, both FM + PCA-200 ----
    rows = []
    for fm in FMS + ['PCA-200']:
        mono_path = os.path.join(MONO_DIR, f'oof_monocyte_{fm}_seed42.csv')
        if not os.path.exists(mono_path):
            print(f'  SKIP {fm}: monocyte seed42 file missing (likely Geneformer scope gap)',
                  flush=True)
            continue
        mono_df = pd.read_csv(mono_path)
        mono_df = mono_df[mono_df['patient_id'].isin(intersection)]
        mono_pat = mono_df.groupby('patient_id').agg(y=('y', 'max'), score=('p_lr', 'mean'))

        plt_path = os.path.join(PLATELET_DIR, f'oof_patient_{fm}_seed42.csv') if fm != 'PCA-200' \
            else os.path.join(PLATELET_DIR, 'oof_PCA-200_seed42.csv')
        if not os.path.exists(plt_path):
            print(f'  SKIP {fm}: platelet seed42 file missing', flush=True)
            continue
        plt_df = pd.read_csv(plt_path)
        if 'patient_id' not in plt_df.columns:
            plt_df = plt_df.merge(pmap[['sample_id', 'patient_id']], on='sample_id', how='inner')
        plt_df = plt_df[plt_df['patient_id'].isin(intersection)]
        plt_pat = plt_df.groupby('patient_id').agg(y=('y', 'max'), score=('p_lr', 'mean'))

        common = mono_pat.join(plt_pat, how='inner', lsuffix='_mono', rsuffix='_plt')
        if len(common) < 5 or len(np.unique(common['y_mono'])) < 2:
            print(f'  {fm}: {len(common)} matched patients, insufficient class balance for AUC',
                  flush=True)
            rows.append({'fm': fm, 'n_matched': len(common), 'auc_monocyte': np.nan,
                         'auc_platelet': np.nan})
            continue
        from sklearn.metrics import roc_auc_score
        auc_mono = roc_auc_score(common['y_mono'], common['score_mono'])
        auc_plt = roc_auc_score(common['y_plt'], common['score_plt'])
        rows.append({'fm': fm, 'n_matched': len(common), 'n_fatal_matched': int(common['y_mono'].sum()),
                     'auc_monocyte': round(auc_mono, 4), 'auc_platelet': round(auc_plt, 4),
                     'delta': round(auc_mono - auc_plt, 4)})
        print(f'  {fm}: n={len(common)} monocyte_AUC={auc_mono:.4f} platelet_AUC={auc_plt:.4f}',
              flush=True)

    result_df = pd.DataFrame(rows)
    result_df.to_csv(os.path.join(OUT_DIR, 'matched_comparison_auc.csv'), index=False)
    print(f'\nSAVED -> {OUT_DIR}', flush=True)


if __name__ == '__main__':
    main()
