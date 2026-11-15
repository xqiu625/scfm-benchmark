"""
Second half of the controlled comparison (see
ren_diagnostic_cellcount_vs_contam.py for the first, which ruled out small-N
instability). Combines 10 Ren channels spanning the full size range
(57-18,236 cells) into one small matrix, so run_decontx_per_channel.R's
ORIGINAL multi-batch decontX(batch=channel_id) call -- the exact approach
already validated on Liu -- can be run on them directly, for comparison
against the already-computed isolated-per-channel result for the same 10
channels.

Result: every one of the 10 channels' contamination estimates matched the
isolated-per-channel run to 6+ decimal places. Per-channel isolation (forced
by R's 32-bit sparse-matrix ceiling on Ren's full-cohort matrix, see
split_matrix_per_channel.py) is mathematically identical to celda's own
batch= splitting -- confirms Ren's elevated contamination (mean 52.7% vs
Liu's 29.4%) is genuine signal, not a methodological artifact. Full writeup
in project_302006_rebuild_architecture memory.

ren_comparison_channels.txt (the 10 channel IDs, one per line) is a small
sibling file in this directory.
"""
import pandas as pd
from scipy.io import mmread, mmwrite
from scipy import sparse
import os

channels = [l.strip() for l in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ren_comparison_channels.txt'))]
split_dir = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/ren_decontx_input/split_per_channel'
out_dir = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/ren_decontx_comparison'
os.makedirs(out_dir, exist_ok=True)

mats, cells_list = [], []
for ch in channels:
    m = mmread(os.path.join(split_dir, ch, 'counts.mtx')).tocsc()
    c = pd.read_csv(os.path.join(split_dir, ch, 'cells.csv'))
    mats.append(m)
    cells_list.append(c)
    print(ch, m.shape, len(c))

combined = sparse.hstack(mats, format='csc')
cells = pd.concat(cells_list, ignore_index=True)
assert combined.shape[1] == len(cells)
mmwrite(os.path.join(out_dir, 'counts.mtx'), combined)
cells.to_csv(os.path.join(out_dir, 'cells.csv'), index=False)
print(f'Combined: {combined.shape}, {len(cells)} cells, {cells.channel_id.nunique()} channels')
