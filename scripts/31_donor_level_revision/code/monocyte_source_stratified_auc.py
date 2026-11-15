#!/usr/bin/env python3
"""
Job F follow-up: re-aggregate the existing monocyte repeated-CV OOF predictions
(no new CV -- same 10-seed StratifiedGroupKFold OOF used for the main
in-distribution table) under two additional cuts:
  1. COVID-only (Ren + Stephenson, dropping Qiu-sepsis) -- for direct disease-
     composition comparability with the platelet within-COVID (240 patients/
     25 fatal) comparison point.
  2. Per-source-study (Ren / Stephenson / Qiu-sepsis) -- tests whether the
     lowest-provenance arm (Qiu's marker-gene threshold call, not a curated
     atlas label like Ren/Stephenson) is dragging down the pooled aggregate.

Reads oof_monocyte_<FM>_seed<NN>.csv (source_study column already present from
build_monocyte_canonical.py / rerun_monocyte_mortality_patient.py) --
re-aggregation only, patient-level AUC computed within each cut's own subset
of patients each seed (StratifiedGroupKFold folds don't need to be redone;
each cell's OOF prediction, restricted to the cut's cells, aggregated to its
patient, then scored against that cut's patients only).
"""
import glob
import os

import numpy as np
import pandas as pd

from repeated_cv import DEFAULT_SEEDS, repeated_patient_auc, summarize_repeated

SRC_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2/monocyte'
OUT_DIR = SRC_DIR
METHODS = ['STATE', 'UCE', 'scPRINT-2', 'scGPT', 'PCA-200']
N_SEEDS = 10
SEEDS = DEFAULT_SEEDS[:N_SEEDS]

CUTS = {
    'ALL': None,
    'COVID_only_Ren_Stephenson': ['Ren_GSE158055', 'Stephenson_EMTAB10026'],
    'Ren_only': ['Ren_GSE158055'],
    'Stephenson_only': ['Stephenson_EMTAB10026'],
    'Qiu_sepsis_only': ['Qiu_sepsis_GSE167363'],
}


def load_seed_files(method):
    by_seed = {}
    for seed in SEEDS:
        path = os.path.join(SRC_DIR, f'oof_monocyte_{method}_seed{seed}.csv')
        by_seed[seed] = pd.read_csv(path)
    return by_seed


def main():
    rows = []
    for method in METHODS:
        by_seed = load_seed_files(method)
        for cut_name, studies in CUTS.items():
            oof_by_seed = {}
            pid_ref = None
            y_ref = None
            for seed in SEEDS:
                df = by_seed[seed]
                if studies is not None:
                    df = df[df['source_study'].isin(studies)]
                oof_by_seed[seed] = {'p_lr': df['p_lr'].to_numpy()}
                pid_this = df['patient_id'].to_numpy()
                y_this = df['y'].to_numpy()
                if pid_ref is None:
                    pid_ref, y_ref = pid_this, y_this
                elif len(pid_this) != len(pid_ref):
                    raise AssertionError(
                        f'{method}/{cut_name} seed{seed}: row count differs from seed{SEEDS[0]} '
                        f'({len(pid_this)} vs {len(pid_ref)}) -- OOF files should have identical '
                        f'cell sets across seeds (same FM embedding input each time)')
            aucs, seeds_out = repeated_patient_auc(oof_by_seed, pid_ref, y_ref, score_key='p_lr')
            summ = summarize_repeated(aucs)
            n_patients = len(np.unique(pid_ref))
            n_fatal = len(np.unique(pid_ref[y_ref == 1]))
            rows.append({
                'cut': cut_name, 'fm': method, 'n_patients': n_patients,
                'n_fatal_patients': n_fatal, 'n_seeds': summ['n'],
                'auc_repeat_mean': round(summ['mean'], 4) if not np.isnan(summ['mean']) else None,
                'auc_repeat_sd': round(summ['sd'], 4) if not np.isnan(summ['sd']) else None,
            })
            print(f"{method:10s} {cut_name:28s} n_pat={n_patients:4d} fatal={n_fatal:3d} "
                  f"AUC={summ['mean']:.4f}+/-{summ['sd']:.4f}" if not np.isnan(summ['mean'])
                  else f"{method:10s} {cut_name:28s} n_pat={n_patients:4d} fatal={n_fatal:3d} AUC=n/a (single class)",
                  flush=True)

    out = pd.DataFrame(rows)
    out_path = os.path.join(OUT_DIR, 'monocyte_source_stratified_auc.csv')
    out.to_csv(out_path, index=False)
    print(f'\nSAVED -> {out_path}', flush=True)


if __name__ == '__main__':
    main()
