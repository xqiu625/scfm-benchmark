"""
COMBAT gene-level sanity check (2026-09-16), requested before trusting COMBAT's
DecontX output: B_cell contamination is a striking outlier (mean 79.9% across
all 70 channels, uniform 52-89% range, ruled out as a per-channel-isolation or
small-N artifact by prior checks - see project_302006_rebuild_architecture
memory). Per-cell UMI count is also NOT lower for B_cell than other lineages
(combat_lineage_umi_check.py - median 3615/mean 6143, second-highest of all 9
lineages), ruling out a low-native-RNA-content denominator effect too.

This script checks the actual mechanism: for B-cell-labeled cells in 5 test
channels spanning the observed contamination range, which genes account for
the most REMOVED counts (raw - decontaminated)? If dominated by canonical
B-cell identity markers (CD79A/MS4A1/CD19/IGHM/IGHD/etc.), that would indicate
DecontX is stripping real B-cell signal, not ambient contamination. If instead
dominated by ubiquitous housekeeping/ribosomal genes shared across all cell
types, the "contamination" is not selectively destroying B-cell identity.

Result: top removed genes in all 5 channels are ACTB/B2M/EEF1A1/ribosomal
protein family (RPL*/RPS*)/FTL/FTH1/mitochondrial genes - universal
housekeeping genes, not B-cell-specific. Canonical B-cell markers account for
only 0.6-1.5% of removed counts. One partial exception: CD74 (MHC-II
invariant chain, real but not a core lineage-defining marker) shows
substantial removal, up to ~90% in some channels.

Per user decision 2026-09-16: this does not meet the "mainly removing
canonical B-cell genes" pause criterion, so proceed with the decontaminated
data. Record as an unresolved QC caveat rather than pursuing the mechanism
further - see qc_caveats.md.
"""
import pandas as pd
import numpy as np
from scipy.io import mmread

BCELL_MARKERS = {'CD79A','CD79B','MS4A1','CD19','IGHM','IGHD','CD74','TCL1A','BANK1',
                  'FCRL1','CD22','VPREB3','POU2AF1','BLK','FCER2','IGHG1','IGHG2','IGHG3',
                  'IGHG4','IGKC','IGLC1','IGLC2','IGLC3','IGLC7','JCHAIN'}
AMBIENT_CANDIDATES = {'HBB','HBA1','HBA2','PPBP','PF4','LYZ','S100A8','S100A9','S100A12',
                       'FCN1','VCAN','NKG7','GNLY','CD3D','CD3E'}

genes_csv = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/combat_decontx_input/genes.csv'
genes = pd.read_csv(genes_csv)['gene'].tolist()

split_dir = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/combat_decontx_input/split_per_channel'
out_dir = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/combat_decontx_output'

test_channels = ['C4', 'K6', 'D4', 'F4', 'F3']

for ch in test_channels:
    raw = mmread(f'{split_dir}/{ch}/counts.mtx').tocsc()
    cells = pd.read_csv(f'{split_dir}/{ch}/cells.csv')
    decont = mmread(f'{out_dir}/combat_{ch}_decontx_counts.mtx').tocsc()
    assert raw.shape == decont.shape
    assert len(cells) == raw.shape[1]

    b_idx = np.where(cells['harmonized_lineage'].to_numpy() == 'B_cell')[0]
    if len(b_idx) == 0:
        continue
    raw_b = raw[:, b_idx]
    decont_b = decont[:, b_idx]

    raw_sum = np.asarray(raw_b.sum(axis=1)).flatten()
    decont_sum = np.asarray(decont_b.sum(axis=1)).flatten()
    removed = raw_sum - decont_sum

    df = pd.DataFrame({'gene': genes, 'raw_total': raw_sum, 'decont_total': decont_sum, 'removed': removed})
    df = df.sort_values('removed', ascending=False)
    total_removed = removed.sum()
    total_raw = raw_sum.sum()

    print(f'=== {ch}: {len(b_idx)} B cells, total raw counts={total_raw:.0f}, '
          f'total removed={total_removed:.0f} ({100*total_removed/total_raw:.1f}%) ===')
    top30 = df.head(30).copy()
    top30['is_bcell_marker'] = top30['gene'].isin(BCELL_MARKERS)
    top30['is_ambient_candidate'] = top30['gene'].isin(AMBIENT_CANDIDATES)
    top30['pct_of_removed'] = 100 * top30['removed'] / total_removed
    print(top30[['gene','raw_total','decont_total','removed','pct_of_removed','is_bcell_marker','is_ambient_candidate']].to_string(index=False))

    n_bcell_marker_removed = df[df['gene'].isin(BCELL_MARKERS)]['removed'].sum()
    n_ambient_removed = df[df['gene'].isin(AMBIENT_CANDIDATES)]['removed'].sum()
    print(f'  Total removed from canonical B-cell markers: {n_bcell_marker_removed:.0f} '
          f'({100*n_bcell_marker_removed/total_removed:.1f}% of all removed counts)')
    print(f'  Total removed from non-B-cell ambient candidates: {n_ambient_removed:.0f} '
          f'({100*n_ambient_removed/total_removed:.1f}% of all removed counts)')
    print()
