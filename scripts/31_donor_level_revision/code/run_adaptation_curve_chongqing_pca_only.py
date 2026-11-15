#!/usr/bin/env python3
"""Standalone rerun of just the PCA-200+LR piece of run_adaptation_curve_chongqing.py
-- the other 7 methods (5 combined-pool FMs + TranscriptFormer/STACK legacy-only)
already completed successfully in job 28148282; only PCA-200 failed (duplicate
gene symbols on reindex, since fixed in the shared module). Writes its own CSV;
merging into all_rows.csv/summary.csv is a separate cheap step."""
import os

import run_adaptation_curve_chongqing as m

print('Loading PCA raw data (combined platelet + Chongqing)...', flush=True)
Xraw, draw, sraw, sidraw, pidraw = m.load_pca_raw()
print(f'  loaded {Xraw.shape}', flush=True)

tm, stm, ctm, test_patients, yraw, split_doc = m.make_sepsis_holdout(
    draw, sraw, pidraw, seed=m.RANDOM_STATE)
print(f'  test patients: {len(test_patients)}, train pool cells: {int(stm.sum())}', flush=True)

df_pca = m.adaptation_curve(f'PCA-{m.PCA_N_COMP}+LR', Xraw, yraw, tm, stm, ctm, pidraw,
                             m.N_SEPSIS_GRID, m.N_SEEDS,
                             transform=m.pca_transform_factory(m.PCA_N_COMP),
                             n_bootstrap=m.N_BOOTSTRAP, test_set='expanded_103pt')
out_path = os.path.join(m.OUT_DIR, f'PCA-{m.PCA_N_COMP}+LR.csv')
df_pca.to_csv(out_path, index=False)
print(f'SAVED -> {out_path}', flush=True)
