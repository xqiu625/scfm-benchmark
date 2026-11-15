#!/usr/bin/env python3
"""Job S adaptation curve, PCA-200-only rerun + final compile.

run_adaptation_curve_chongqing.py (sbatch 28148282) FAILED at the PCA-200
step with `ValueError: cannot reindex on an axis with duplicate labels` in
load_pca_raw() -- the Chongqing raw-counts gene-symbol index has duplicate
labels (multiple Ensembl IDs sharing one symbol), which broke a plain
pd.Series.reindex(). Fixed in run_adaptation_curve_chongqing.py via
_first_occurrence_index() (dict-based first-occurrence lookup instead of
reindex). The 7 FM/method curves from that run already completed and are on
disk in OUT_DIR -- this script reruns ONLY PCA-200+LR with the fix, then
compiles all_rows.csv + summary.csv from all 8 methods' per-method CSVs
(reading from disk, not rerunning the 7 already-done methods).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_adaptation_curve_chongqing import (  # noqa: E402
    OUT_DIR, N_SEPSIS_GRID, N_SEEDS, N_BOOTSTRAP, PCA_N_COMP, RANDOM_STATE,
    load_pca_raw, make_sepsis_holdout, adaptation_curve, pca_transform_factory,
    summarize,
)
import pandas as pd  # noqa: E402

EXPECTED_METHODS = ['STATE', 'UCE', 'scPRINT-2', 'Geneformer', 'scGPT',
                    'STACK', 'TranscriptFormer', f'PCA-{PCA_N_COMP}+LR']


def main():
    pca_csv = os.path.join(OUT_DIR, f'PCA-{PCA_N_COMP}+LR.csv')
    if os.path.exists(pca_csv):
        print(f'PCA-{PCA_N_COMP}+LR.csv already exists, skipping rerun and '
              f'going straight to compile.', flush=True)
    else:
        print(f'=== PCA-{PCA_N_COMP} + LR (from-scratch baseline, Chongqing-expanded, '
              f'FIXED reindex) ===', flush=True)
        Xraw, draw, sraw, sidraw, pidraw = load_pca_raw()
        tm, stm, ctm, _, yraw, _ = make_sepsis_holdout(draw, sraw, pidraw, seed=RANDOM_STATE)
        df_pca = adaptation_curve(f'PCA-{PCA_N_COMP}+LR', Xraw, yraw, tm, stm, ctm, pidraw,
                                   N_SEPSIS_GRID, N_SEEDS, transform=pca_transform_factory(PCA_N_COMP),
                                   n_bootstrap=N_BOOTSTRAP)
        df_pca.to_csv(pca_csv, index=False)
        print(f'SAVED -> {pca_csv}', flush=True)

    # ---- compile all 8 methods from disk (7 already-done + fresh PCA-200) ----
    all_dfs = []
    missing = []
    for method in EXPECTED_METHODS:
        path = os.path.join(OUT_DIR, f'{method}.csv')
        if not os.path.exists(path):
            missing.append(method)
            continue
        df = pd.read_csv(path)
        n_nan = df[['auc_point', 'auc_boot_mean']].isna().sum().sum()
        if n_nan:
            print(f'  WARNING: {method}.csv has {n_nan} NaN values in auc columns', flush=True)
        all_dfs.append(df)
    if missing:
        print(f'MISSING method CSVs (not compiled): {missing}', flush=True)

    pooled = pd.concat(all_dfs, ignore_index=True)
    pooled.to_csv(os.path.join(OUT_DIR, 'all_rows.csv'), index=False)
    summary = summarize(pooled)
    summary.to_csv(os.path.join(OUT_DIR, 'summary.csv'), index=False)
    print(f'\nCompiled {len(all_dfs)}/{len(EXPECTED_METHODS)} methods '
          f'({len(pooled)} rows) -> all_rows.csv, summary.csv', flush=True)
    print('\n=== summary.csv ===')
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
