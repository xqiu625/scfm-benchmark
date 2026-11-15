"""
Assemble the within-COVID validation set (Liu + Ren + COMBAT -- the 3
cohorts with DecontX done so far; Cillo not yet run) for the pre-specified
before-vs-after-decontamination robustness rerun (benchmark_design_decisions.md
Phases 3-4 validation order; user decision 2026-09-16 to prioritize this over
continuing Phase 1 to Cillo). Produces TWO matched combined matrices -- raw
and decontaminated -- over the identical cell/gene/donor set, so any
downstream difference between them is attributable to decontamination alone.

Labels: canonical_donor_labels.csv, filtered to disease=='COVID-19' and
cohort in {Liu 2021, Ren 2021, COMBAT 2022}. Genes: intersection of all three
cohorts' gene lists (16,603 genes, verified 2026-09-16), fixed sorted order,
identical for both output matrices.

Per-cohort channel/donor resolution (verified against each cohort's own
staged data, not assumed uniform):
  - Liu: single combined liu_decontx_input/{counts.mtx,cells.csv} +
    liu_decontx_input/liu_decontx_counts.mtx (same cell order as cells.csv,
    since run_decontx_per_channel.R sets colnames(sce) directly from
    cells$cell_id and decontXcounts(sce) preserves that order). donor_id is
    already a cells.csv column.
  - Ren: per-channel split (ren_decontx_input/split_per_channel/<channel>/
    {counts.mtx,cells.csv} for raw, ren_decontx_output/ren_<channel>_
    decontx_counts.mtx for decontaminated, same per-channel cell order).
    channel_id is the sampleID, NOT the donor -- donor comes from
    ren_gse158055_channel_inventory.csv's sampleID->donor column.
  - COMBAT: per-channel split (combat_decontx_input/split_per_channel/<channel>/
    {counts.mtx,cells.csv}, combat_decontx_output/combat_<channel>_
    decontx_counts.mtx). donor_id is already a cells.csv column (COMBAT_ID);
    donors can span multiple channels (COMBAT's pool structure), which this
    script handles by just processing every channel independently and letting
    donor_id do the grouping downstream, not by needing per-donor channel
    lookups the way Ren does.

Usage:
    python assemble_within_covid_validation.py --out-dir <output dir>

Output: <out_dir>/{raw,decontx}/{counts.mtx,cells.csv}, genes.csv (shared).
cells.csv columns: cell_id, cohort, channel_id, donor_id, event.
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread, mmwrite

RESULTS = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc'
LABELS_CSV = '/bigdata/godziklab/shared/Xinru/302006/302006_updategit/reports/manuscript/revision/literature_audit/canonical_donor_labels.csv'
REN_INVENTORY = '/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/32_rebuild_upfront_qc/channel_inventories/ren_gse158055_channel_inventory.csv'


def load_intersection_genes():
    liu = pd.read_csv(f'{RESULTS}/liu_decontx_input/genes.csv')['gene']
    ren = pd.read_csv(f'{RESULTS}/ren_decontx_input/genes.csv')['gene']
    combat = pd.read_csv(f'{RESULTS}/combat_decontx_input/genes.csv')['gene']
    inter = sorted(set(liu) & set(ren) & set(combat))
    return inter, {'liu': liu, 'ren': ren, 'combat': combat}


def subset_to_genes(mat, cohort_genes, target_genes):
    """mat: genes x cells sparse. Returns target_genes x cells, reordered/subset."""
    idx_map = {g: i for i, g in enumerate(cohort_genes)}
    row_idx = np.array([idx_map[g] for g in target_genes])
    return mat[row_idx, :]


def process_liu(target_genes, wanted_donors):
    """Liu's demuxlet-derived donor_id carries a timepoint/batch suffix for
    B2/B3-lane donors (e.g. HGR0000430_T1), which does NOT exact-match
    canonical_donor_labels.csv's bare donor_id (HGR0000430) -- verified
    directly 2026-09-16 (HGR0000430, one of the 4 confirmed-deceased donors,
    was silently dropped by an exact-match join). B1-lane donors have no
    suffix and match directly either way. Match on the base id (everything
    before the first underscore) instead of exact string equality."""
    genes = pd.read_csv(f'{RESULTS}/liu_decontx_input/genes.csv')['gene'].tolist()
    cells = pd.read_csv(f'{RESULTS}/liu_decontx_input/cells.csv')
    base_id = cells['donor_id'].str.split('_').str[0]
    keep = base_id.isin(wanted_donors)
    if keep.sum() == 0:
        return None
    raw = mmread(f'{RESULTS}/liu_decontx_input/counts.mtx').tocsc()
    decontx = mmread(f'{RESULTS}/liu_decontx_input/liu_decontx_counts.mtx').tocsc()
    idx = np.where(keep.to_numpy())[0]
    raw_sub = subset_to_genes(raw[:, idx], genes, target_genes)
    decontx_sub = subset_to_genes(decontx[:, idx], genes, target_genes)
    meta = cells.iloc[idx][['cell_id', 'channel_id']].copy()
    meta['donor_id'] = base_id.iloc[idx].to_numpy()  # normalized, matches canonical_donor_labels.csv
    meta['cohort'] = 'Liu 2021'
    print(f'  Liu: {len(idx)} cells from {meta.donor_id.nunique()} donors', flush=True)
    return raw_sub, decontx_sub, meta


def process_ren(target_genes, wanted_donors):
    genes = pd.read_csv(f'{RESULTS}/ren_decontx_input/genes.csv')['gene'].tolist()
    inv = pd.read_csv(REN_INVENTORY)
    sample_to_donor = dict(zip(inv['sampleID'], inv['donor']))
    wanted_samples = [s for s, d in sample_to_donor.items() if d in wanted_donors]

    raw_blocks, decontx_blocks, meta_rows = [], [], []
    for i, sample in enumerate(wanted_samples):
        ch_dir = f'{RESULTS}/ren_decontx_input/split_per_channel/{sample}'
        decontx_path = f'{RESULTS}/ren_decontx_output/ren_{sample}_decontx_counts.mtx'
        if not (os.path.exists(f'{ch_dir}/counts.mtx') and os.path.exists(decontx_path)):
            print(f'  Ren: WARNING missing files for sample {sample}, skipping', flush=True)
            continue
        cells = pd.read_csv(f'{ch_dir}/cells.csv')
        raw = mmread(f'{ch_dir}/counts.mtx').tocsc()
        decontx = mmread(decontx_path).tocsc()
        raw_blocks.append(subset_to_genes(raw, genes, target_genes))
        decontx_blocks.append(subset_to_genes(decontx, genes, target_genes))
        m = cells[['cell_id']].copy()
        m['channel_id'] = sample
        m['donor_id'] = sample_to_donor[sample]
        m['cohort'] = 'Ren 2021'
        meta_rows.append(m)
        if (i + 1) % 40 == 0:
            print(f'  Ren: {i+1}/{len(wanted_samples)} samples loaded', flush=True)
    if not raw_blocks:
        return None
    raw_all = sparse.hstack(raw_blocks, format='csc')
    decontx_all = sparse.hstack(decontx_blocks, format='csc')
    meta_all = pd.concat(meta_rows, ignore_index=True)
    print(f'  Ren: {raw_all.shape[1]} cells from {meta_all.donor_id.nunique()} donors, '
          f'{len(wanted_samples)} samples', flush=True)
    return raw_all, decontx_all, meta_all


def process_combat(target_genes, wanted_donors):
    genes = pd.read_csv(f'{RESULTS}/combat_decontx_input/genes.csv')['gene'].tolist()
    split_dir = f'{RESULTS}/combat_decontx_input/split_per_channel'
    channels = sorted(os.listdir(split_dir))

    raw_blocks, decontx_blocks, meta_rows = [], [], []
    for i, ch in enumerate(channels):
        ch_dir = f'{split_dir}/{ch}'
        if not os.path.isdir(ch_dir):
            continue
        cells = pd.read_csv(f'{ch_dir}/cells.csv')
        keep = cells['donor_id'].isin(wanted_donors)
        if keep.sum() == 0:
            continue
        decontx_path = f'{RESULTS}/combat_decontx_output/combat_{ch}_decontx_counts.mtx'
        if not os.path.exists(decontx_path):
            print(f'  COMBAT: WARNING missing decontx output for channel {ch}, skipping', flush=True)
            continue
        raw = mmread(f'{ch_dir}/counts.mtx').tocsc()
        decontx = mmread(decontx_path).tocsc()
        idx = np.where(keep.to_numpy())[0]
        raw_blocks.append(subset_to_genes(raw[:, idx], genes, target_genes))
        decontx_blocks.append(subset_to_genes(decontx[:, idx], genes, target_genes))
        m = cells.iloc[idx][['cell_id', 'channel_id', 'donor_id']].copy()
        m['cohort'] = 'COMBAT 2022'
        meta_rows.append(m)
        if (i + 1) % 20 == 0:
            print(f'  COMBAT: {i+1}/{len(channels)} channels scanned', flush=True)
    if not raw_blocks:
        return None
    raw_all = sparse.hstack(raw_blocks, format='csc')
    decontx_all = sparse.hstack(decontx_blocks, format='csc')
    meta_all = pd.concat(meta_rows, ignore_index=True)
    print(f'  COMBAT: {raw_all.shape[1]} cells from {meta_all.donor_id.nunique()} donors', flush=True)
    return raw_all, decontx_all, meta_all


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()
    os.makedirs(f'{args.out_dir}/raw', exist_ok=True)
    os.makedirs(f'{args.out_dir}/decontx', exist_ok=True)

    print('Loading labels and gene intersection...', flush=True)
    labels = pd.read_csv(LABELS_CSV)
    target_genes, _ = load_intersection_genes()
    print(f'  {len(target_genes)} shared genes', flush=True)

    raw_all_blocks, decontx_all_blocks, meta_all_blocks = [], [], []
    for cohort_name, proc_fn in [('Liu 2021', process_liu), ('Ren 2021', process_ren),
                                   ('COMBAT 2022', process_combat)]:
        sub = labels[(labels['cohort'] == cohort_name) & (labels['disease'] == 'COVID-19')]
        wanted_donors = set(sub['donor_id'])
        event_map = dict(zip(sub['donor_id'], sub['event']))
        print(f'\n{cohort_name}: {len(wanted_donors)} within-COVID donors in canonical_donor_labels.csv', flush=True)
        result = proc_fn(target_genes, wanted_donors)
        if result is None:
            print(f'  {cohort_name}: NO CELLS FOUND -- skipping entirely', flush=True)
            continue
        raw_sub, decontx_sub, meta = result
        meta['event'] = meta['donor_id'].map(event_map)
        raw_all_blocks.append(raw_sub)
        decontx_all_blocks.append(decontx_sub)
        meta_all_blocks.append(meta)

    print('\nCombining all cohorts...', flush=True)
    raw_combined = sparse.hstack(raw_all_blocks, format='csc')
    decontx_combined = sparse.hstack(decontx_all_blocks, format='csc')
    meta_combined = pd.concat(meta_all_blocks, ignore_index=True)
    assert raw_combined.shape == decontx_combined.shape
    assert raw_combined.shape[1] == len(meta_combined)

    pd.DataFrame({'gene': target_genes}).to_csv(f'{args.out_dir}/genes.csv', index=False)
    mmwrite(f'{args.out_dir}/raw/counts.mtx', raw_combined)
    mmwrite(f'{args.out_dir}/decontx/counts.mtx', decontx_combined)
    meta_combined.to_csv(f'{args.out_dir}/raw/cells.csv', index=False)
    meta_combined.to_csv(f'{args.out_dir}/decontx/cells.csv', index=False)

    print(f'\nDone. {raw_combined.shape[1]} cells x {raw_combined.shape[0]} genes.')
    print(f'Donors: {meta_combined.donor_id.nunique()}, '
          f'events: {meta_combined.groupby("donor_id")["event"].max().sum()}')
    print(meta_combined.groupby('cohort')['donor_id'].nunique())


if __name__ == '__main__':
    main()
