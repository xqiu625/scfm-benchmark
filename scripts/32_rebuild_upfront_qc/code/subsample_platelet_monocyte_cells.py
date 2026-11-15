"""
Priority (b) step 1: build the ONE shared stratified subsample (cap 60
cells/donor, matching the Phase 3 pilot convention) of the final 8-cohort
Platelet + Monocyte_Myeloid dataset, used identically by all 5 FMs (Cindy's
2026-09-19 decision -- subsampled, not full ~1.2M-cell extraction). Doing
this once, shared, is deliberate: every FM must see the EXACT same cells for
the cross-model comparison to be fair, so this step -- not each FM's own
extraction script -- owns cell selection.

Excludes cells with undefined event (NaN, 53 donors/50,588 cells -- see
qc_caveats.md, caught 2026-09-19) before subsampling, since those donors
cannot contribute to the classification task at all regardless of which FM
embeds them.

Streams both raw and decontx counts.mtx (22G/50G) via read_mtx_columns
rather than mmread-then-subset, for the same reason as
assemble_platelet_monocyte_benchmark_input.py and
pca200_pooled_secondary_cv.py: files this size risk the memory/walltime
blowup Chongqing's DecontX write step hit.

Usage:
    python subsample_platelet_monocyte_cells.py --which raw|decontx
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.io import mmwrite

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from assemble_platelet_monocyte_benchmark_input import read_mtx_columns

DATA_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/platelet_monocyte_benchmark_input_v2'
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/fm_embeddings/subsampled_input'
CELLS_PER_DONOR = 60
RNG_SEED = 20260917


def stratified_subsample_idx(cells_df, cap, seed):
    """One row per (donor_id, harmonized_lineage) group, capped -- keeps the
    Platelet/Monocyte_Myeloid balance from drifting per-donor toward whichever
    lineage happens to be more numerous in that donor's raw data."""
    rng = np.random.default_rng(seed)
    keep = []
    for (_donor, _lineage), grp in cells_df.groupby(['donor_id', 'harmonized_lineage']):
        take = min(len(grp), cap)
        keep.extend(rng.choice(grp.index.to_numpy(), size=take, replace=False).tolist())
    return np.sort(np.array(keep))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--which', required=True, choices=['raw', 'decontx'])
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)

    print(f'Loading {args.which} cells.csv...', flush=True)
    cells = pd.read_csv(f'{DATA_DIR}/{args.which}/cells.csv', low_memory=False)
    n_total = len(cells)
    cells = cells[cells['event'].notna()].reset_index(drop=True)
    print(f'  {len(cells)} / {n_total} cells have a defined event', flush=True)

    sub_idx = stratified_subsample_idx(cells, CELLS_PER_DONOR, RNG_SEED)
    sub_cells = cells.iloc[sub_idx].reset_index(drop=True)
    print(f'  stratified subsample (cap={CELLS_PER_DONOR}/donor/lineage): {len(sub_cells)} cells, '
          f'{sub_cells.donor_id.nunique()} donors', flush=True)
    print(sub_cells.groupby(['cohort', 'harmonized_lineage'])['cell_id'].count(), flush=True)

    print(f'Streaming {args.which} counts.mtx for the {len(sub_idx)} subsampled cells...', flush=True)
    mat = read_mtx_columns(f'{DATA_DIR}/{args.which}/counts.mtx', sub_idx)  # genes x cells
    genes = pd.read_csv(f'{DATA_DIR}/genes.csv')['gene']

    mmwrite(f'{OUT_DIR}/{args.which}_counts.mtx', mat)
    genes.to_csv(f'{OUT_DIR}/genes.csv', index=False)
    sub_cells.to_csv(f'{OUT_DIR}/{args.which}_cells.csv', index=False)

    print(f'\nDone. {mat.shape[1]} cells x {mat.shape[0]} genes -> {OUT_DIR}/{args.which}_counts.mtx')


if __name__ == '__main__':
    main()
