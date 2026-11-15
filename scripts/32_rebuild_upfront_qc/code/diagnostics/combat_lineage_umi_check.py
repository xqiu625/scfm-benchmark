"""
Second COMBAT sanity check (2026-09-16), alongside combat_bcell_ambient_check.py:
tests whether B_cell's outlier contamination (mean 79.9%) is explained by B
cells having systematically lower native RNA content than other lineages (a
known mechanism - a roughly fixed per-droplet ambient contribution becomes a
larger fraction of a smaller denominator).

Result: refuted. B_cell median/mean raw UMI (3615/6143) is NOT lower than
other lineages - it's the second-highest mean of all 9 lineages, above
NK_cell, Monocyte_Myeloid, T_cell, DC, and Platelet. See
project_302006_rebuild_architecture memory and qc_caveats.md for the full
writeup and decision.
"""
import pandas as pd
import numpy as np
from scipy.io import mmread
import glob, os

split_dir = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/combat_decontx_input/split_per_channel'
test_channels = ['C4', 'K6', 'D4', 'F4', 'F3']

rows = []
for ch in test_channels:
    raw = mmread(f'{split_dir}/{ch}/counts.mtx').tocsc()
    cells = pd.read_csv(f'{split_dir}/{ch}/cells.csv')
    umi = np.asarray(raw.sum(axis=0)).flatten()
    cells['umi'] = umi
    cells['channel'] = ch
    rows.append(cells)

df = pd.concat(rows, ignore_index=True)
print('Median raw UMI per cell, by harmonized_lineage (pooled across the 5 test channels):')
print(df.groupby('harmonized_lineage')['umi'].agg(['median','mean','count']).sort_values('median'))
