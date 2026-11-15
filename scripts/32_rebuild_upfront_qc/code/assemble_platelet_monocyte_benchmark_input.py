"""
Assemble the decontaminated Platelet + Monocyte_Myeloid cells across all
eight Phase-1-complete cohorts (Liu, Ren, COMBAT, Cillo, Qiu, Chongqing,
Stephenson, Wilk), richly labeled, ready to feed to FM embedding-extraction
scripts for the full benchmark. NOT restricted to within-COVID -- this
covers every donor with a canonical label in these cohorts (unlike
assemble_within_covid_validation.py, which is the disease-restricted
validation subset this script is modeled on). Schulte-Schrepping (EGA
controlled-access) and Bernardes (FastGenomics-gated) are not in Phase 1 yet
and so are not here either.

Qiu/Chongqing/Stephenson/Wilk all share Liu/Cillo's single-combined-matrix
build shape (counts.mtx + <tag>_decontx_counts.mtx + cells.csv + genes.csv
in one directory) -- process_simple_cohort() handles all four generically;
Liu keeps its own function only because of its donor_id base-id-stripping
special case (see process_liu's docstring below).

Cell-type filter: harmonized_lineage in {Platelet, Monocyte_Myeloid} (the
project's two thesis-relevant lineages -- Platelet is the endpoint itself,
Monocyte_Myeloid is the largest myeloid population and the lineage most
implicated in COVID/sepsis immunopathology).

Labels per cell: cohort, channel_id, donor_id, harmonized_lineage, disease,
event, event_definition (all from canonical_donor_labels.csv except
harmonized_lineage, which is cells.csv's own z-label), plus draw_time --
resolved per-cohort at DIFFERENT precision, verified directly against each
cohort's own metadata, not assumed uniform:
  - Ren: categorical sample_time (progression/convalescence/control) from
    ren_gse158055_channel_inventory.csv, per-channel (=per-sampleID) exact.
  - COMBAT: TimeSinceOnset, a per-cell numeric days-since-onset field already
    in the deposited h5ad's obs (joined back via cell_id, which IS the h5ad's
    own barcode index -- verified in build_combat_decontx_input.py). NaN for
    cells with no onset recorded (e.g. healthy controls).
  - Cillo: exact integer day, re-parsed from the h5ad's patient_id string
    (e.g. "COVID_11_day_5_PBMC" -> 5) via the same barcode+h5ad_suffix join
    build_cillo_decontx_input.py used -- NOT stored in Cillo's own cells.csv
    (which only kept the donor_id, stripped of its day token), so
    reconstructed here from the h5ad directly. NaN for the two ALIR-named
    patient_ids and Healthy donors, which carry no day token.
  - Liu: NOT reliably resolvable at cell level -- verified directly
    (T3_patient_sample_mapping/memo.md): GEO's own "timepoint" characteristic
    for Liu's pooled 10x lanes is literally "mixed", and demuxlet resolves
    donor GENOTYPE only, which is identical across a donor's repeat draws, so
    genotype-based demultiplexing cannot distinguish which of a donor's
    timepoints a given cell came from when more than one of that donor's
    draws is pooled in the same lane. A partial signal exists: demuxlet's
    donor_id carries a "_T<n>" suffix for some B2/B3-lane rows (the resolved
    donor happened to only have one draw in that lane) -- used where present,
    NaN elsewhere. draw_time_precision flags this explicitly so it is never
    silently treated as equally reliable as Ren/COMBAT/Cillo's.
  - Qiu: draw_time already computed at build time (exact T0/T6 sample tag,
    one GSM = one specific draw, no ambiguity) -- passed through as-is.
  - Chongqing: one GSM = one donor = one draw; there is no timepoint
    variation to record, so draw_time is NaN for every cell
    (draw_time_precision = "single_draw_no_variation", distinct from an
    unresolved case -- there is nothing to resolve).
  - Stephenson: draw_time already computed at build time (Days_from_onset,
    an exact per-cell numeric field in the deposited h5ad).
  - Wilk: draw_time already computed at build time (the "d<N>" day suffix on
    some donor_ids, e.g. "...d0"/"...d2"); NaN for the (majority) single,
    unspecified-day draws.

Usage:
    python assemble_platelet_monocyte_benchmark_input.py --out-dir <output dir>

Output: <out_dir>/{raw,decontx}/counts.mtx, <out_dir>/genes.csv (shared,
gene intersection across all 8 cohorts), <out_dir>/cells.csv (shared cell
metadata/labels, identical cell order for both raw and decontx matrices):
    cell_id, cohort, channel_id, donor_id, harmonized_lineage, disease,
    event, event_definition, draw_time, draw_time_precision
"""
import argparse
import os
import re

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread, mmwrite

RESULTS = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc'
UPDATEGIT = '/bigdata/godziklab/shared/Xinru/302006/302006_updategit'
LABELS_CSV = f'{UPDATEGIT}/reports/manuscript/revision/literature_audit/canonical_donor_labels.csv'
REN_INVENTORY = f'{UPDATEGIT}/scripts/32_rebuild_upfront_qc/channel_inventories/ren_gse158055_channel_inventory.csv'
CILLO_INVENTORY = f'{UPDATEGIT}/scripts/32_rebuild_upfront_qc/channel_inventories/cillo_gse180578_channel_inventory.csv'
COMBAT_H5AD = '/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/COMBAT/COMBAT-CITESeq-DATA.h5ad'
CILLO_H5AD = '/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/GSE180578_cillo/GSE180578_cillo_covid19_study_aggregrated_annotated_data.h5ad'
WANTED_LINEAGES = {'Platelet', 'Monocyte_Myeloid'}


COHORT_DIRS = {
    'Liu 2021': 'liu_decontx_input', 'Ren 2021': 'ren_decontx_input',
    'COMBAT 2022': 'combat_decontx_input', 'Cillo 2021': 'cillo_decontx_input',
    'Qiu 2021': 'qiu_decontx_input', 'Chongqing sepsis 2025': 'chongqing_decontx_input',
    'Stephenson 2021': 'stephenson_decontx_input', 'Wilk 2021 JEM': 'wilk_decontx_input',
}


def load_intersection_genes():
    gene_sets = [set(pd.read_csv(f'{RESULTS}/{d}/genes.csv')['gene']) for d in COHORT_DIRS.values()]
    inter = set.intersection(*gene_sets)
    return sorted(inter)


def subset_to_genes(mat, cohort_genes, target_genes):
    idx_map = {g: i for i, g in enumerate(cohort_genes)}
    row_idx = np.array([idx_map[g] for g in target_genes])
    return mat[row_idx, :]


def read_mtx_columns(path, wanted_col_idx, chunk_size=20_000_000):
    """Streams a MatrixMarket file in chunks, keeping only entries whose
    column falls in wanted_col_idx (0-indexed positions into the ORIGINAL
    matrix), remapped to their rank within wanted_col_idx (ascending) as the
    new column index -- never materializes the full matrix in memory.
    Critical at Chongqing/Stephenson's scale (20-44G decontx.mtx files):
    mmread-then-subset on files this size risks the same
    memory/walltime blowup that hit Chongqing's DecontX write step (see
    feedback_hpcc_checkpoint_short_partition memory) -- filtering down to
    the ~2-40k wanted (Platelet/Monocyte) columns while streaming keeps
    memory bounded by the OUTPUT size, not the input file size."""
    with open(path) as f:
        line = f.readline()
        assert line.startswith('%%MatrixMarket'), f'{path}: not a MatrixMarket file'
        n_header = 1
        line = f.readline()
        n_header += 1
        while line.startswith('%'):
            line = f.readline()
            n_header += 1
        n_row, n_col, _nnz = (int(x) for x in line.split())

    col_map = np.full(n_col, -1, dtype=np.int64)
    col_map[wanted_col_idx] = np.arange(len(wanted_col_idx))

    rows_out, cols_out, vals_out = [], [], []
    reader = pd.read_csv(path, sep=r'\s+', skiprows=n_header, header=None,
                          names=['i', 'j', 'v'], chunksize=chunk_size,
                          dtype={'i': np.int64, 'j': np.int64, 'v': np.float64})
    for chunk in reader:
        mapped = col_map[chunk['j'].to_numpy() - 1]
        mask = mapped >= 0
        if mask.any():
            rows_out.append(chunk['i'].to_numpy()[mask] - 1)
            cols_out.append(mapped[mask])
            vals_out.append(chunk['v'].to_numpy()[mask])

    if rows_out:
        rows = np.concatenate(rows_out)
        cols = np.concatenate(cols_out)
        vals = np.concatenate(vals_out)
    else:
        rows = np.array([], dtype=np.int64)
        cols = np.array([], dtype=np.int64)
        vals = np.array([], dtype=np.float64)
    return sparse.csc_matrix((vals, (rows, cols)), shape=(n_row, len(wanted_col_idx)))


def process_liu(target_genes, wanted_donors):
    genes = pd.read_csv(f'{RESULTS}/liu_decontx_input/genes.csv')['gene'].tolist()
    cells = pd.read_csv(f'{RESULTS}/liu_decontx_input/cells.csv')
    base_id = cells['donor_id'].str.split('_').str[0]
    lineage_ok = cells['harmonized_lineage'].isin(WANTED_LINEAGES)
    keep = lineage_ok & base_id.isin(wanted_donors)
    if keep.sum() == 0:
        return None
    raw = mmread(f'{RESULTS}/liu_decontx_input/counts.mtx').tocsc()
    decontx = mmread(f'{RESULTS}/liu_decontx_input/liu_decontx_counts.mtx').tocsc()
    idx = np.where(keep.to_numpy())[0]
    raw_sub = subset_to_genes(raw[:, idx], genes, target_genes)
    decontx_sub = subset_to_genes(decontx[:, idx], genes, target_genes)

    suffix = cells['donor_id'].str.extract(r'_T(\d+)$')[0]
    meta = cells.iloc[idx][['cell_id', 'channel_id', 'harmonized_lineage']].copy()
    meta['donor_id'] = base_id.iloc[idx].to_numpy()
    meta['cohort'] = 'Liu 2021'
    meta['draw_time'] = suffix.iloc[idx].apply(lambda t: f'T{t}' if pd.notna(t) else np.nan).to_numpy()
    meta['draw_time_precision'] = np.where(pd.notna(meta['draw_time']), 'partial_suffix', 'unresolved_mixed_pool')
    print(f'  Liu: {len(idx)} Platelet/Monocyte cells from {meta.donor_id.nunique()} donors '
          f'({meta.harmonized_lineage.value_counts().to_dict()})', flush=True)
    return raw_sub, decontx_sub, meta


def process_ren(target_genes, wanted_donors):
    genes = pd.read_csv(f'{RESULTS}/ren_decontx_input/genes.csv')['gene'].tolist()
    inv = pd.read_csv(REN_INVENTORY)
    sample_to_donor = dict(zip(inv['sampleID'], inv['donor']))
    sample_to_time = dict(zip(inv['sampleID'], inv['sample_time']))
    wanted_samples = [s for s, d in sample_to_donor.items() if d in wanted_donors]

    raw_blocks, decontx_blocks, meta_rows = [], [], []
    for i, sample in enumerate(wanted_samples):
        ch_dir = f'{RESULTS}/ren_decontx_input/split_per_channel/{sample}'
        decontx_path = f'{RESULTS}/ren_decontx_output/ren_{sample}_decontx_counts.mtx'
        if not (os.path.exists(f'{ch_dir}/counts.mtx') and os.path.exists(decontx_path)):
            continue
        cells = pd.read_csv(f'{ch_dir}/cells.csv')
        lineage_ok = cells['harmonized_lineage'].isin(WANTED_LINEAGES)
        if lineage_ok.sum() == 0:
            continue
        idx = np.where(lineage_ok.to_numpy())[0]
        raw = mmread(f'{ch_dir}/counts.mtx').tocsc()
        decontx = mmread(decontx_path).tocsc()
        raw_blocks.append(subset_to_genes(raw[:, idx], genes, target_genes))
        decontx_blocks.append(subset_to_genes(decontx[:, idx], genes, target_genes))
        m = cells.iloc[idx][['cell_id', 'harmonized_lineage']].copy()
        m['channel_id'] = sample
        m['donor_id'] = sample_to_donor[sample]
        m['cohort'] = 'Ren 2021'
        m['draw_time'] = sample_to_time.get(sample, np.nan)
        m['draw_time_precision'] = 'category'
        meta_rows.append(m)
        if (i + 1) % 40 == 0:
            print(f'  Ren: {i+1}/{len(wanted_samples)} samples scanned', flush=True)
    if not raw_blocks:
        return None
    raw_all = sparse.hstack(raw_blocks, format='csc')
    decontx_all = sparse.hstack(decontx_blocks, format='csc')
    meta_all = pd.concat(meta_rows, ignore_index=True)
    print(f'  Ren: {raw_all.shape[1]} Platelet/Monocyte cells from {meta_all.donor_id.nunique()} donors '
          f'({meta_all.harmonized_lineage.value_counts().to_dict()})', flush=True)
    return raw_all, decontx_all, meta_all


def load_combat_onset_map():
    f = h5py.File(COMBAT_H5AD, 'r')
    idx = np.array([x.decode() if isinstance(x, bytes) else x for x in f['obs']['_index'][:]])
    onset = f['obs']['TimeSinceOnset'][:]
    return dict(zip(idx, onset))


def process_combat(target_genes, wanted_donors, onset_map):
    genes = pd.read_csv(f'{RESULTS}/combat_decontx_input/genes.csv')['gene'].tolist()
    split_dir = f'{RESULTS}/combat_decontx_input/split_per_channel'
    channels = sorted(os.listdir(split_dir))

    raw_blocks, decontx_blocks, meta_rows = [], [], []
    for i, ch in enumerate(channels):
        ch_dir = f'{split_dir}/{ch}'
        if not os.path.isdir(ch_dir):
            continue
        cells = pd.read_csv(f'{ch_dir}/cells.csv')
        keep = cells['donor_id'].isin(wanted_donors) & cells['harmonized_lineage'].isin(WANTED_LINEAGES)
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
        m = cells.iloc[idx][['cell_id', 'channel_id', 'donor_id', 'harmonized_lineage']].copy()
        m['cohort'] = 'COMBAT 2022'
        m['draw_time'] = m['cell_id'].map(onset_map)
        m['draw_time_precision'] = 'exact_day'
        meta_rows.append(m)
        if (i + 1) % 20 == 0:
            print(f'  COMBAT: {i+1}/{len(channels)} channels scanned', flush=True)
    if not raw_blocks:
        return None
    raw_all = sparse.hstack(raw_blocks, format='csc')
    decontx_all = sparse.hstack(decontx_blocks, format='csc')
    meta_all = pd.concat(meta_rows, ignore_index=True)
    print(f'  COMBAT: {raw_all.shape[1]} Platelet/Monocyte cells from {meta_all.donor_id.nunique()} donors '
          f'({meta_all.harmonized_lineage.value_counts().to_dict()})', flush=True)
    return raw_all, decontx_all, meta_all


def load_cillo_day_map():
    """Re-derives per-cell draw day from the h5ad's patient_id (Cillo's own
    cells.csv only kept the day-stripped donor_id). Rebuilds the same
    barcode+h5ad_suffix join key build_cillo_decontx_input.py used, keyed by
    cell_id = f'{gsm}_{barcode}' (Cillo's cells.csv convention)."""
    f = h5py.File(CILLO_H5AD, 'r')
    idx = np.array([x.decode() if isinstance(x, bytes) else x for x in f['obs']['_index'][:]])
    pid = np.array([x.decode() if isinstance(x, bytes) else x for x in f['obs']['patient_id'][:]])
    h5ad_day = {}
    for k, p in zip(idx, pid):
        m = re.match(r'^.*_day_(\d+)_(PBMC|ETA)$', p)
        h5ad_day[k] = int(m.group(1)) if m else np.nan

    inv = pd.read_csv(CILLO_INVENTORY).set_index('gsm')
    cell_id_to_day = {}
    for gsm, suffix in inv['h5ad_suffix'].items():
        if pd.isna(suffix) or suffix == '':
            continue
        suffix = str(suffix)
        for k, day in h5ad_day.items():
            if k.endswith('_' + suffix):
                barcode = k.split('_')[0] + '-1'
                cell_id_to_day[f'{gsm}_{barcode}'] = day
    return cell_id_to_day


def process_simple_cohort(dir_name, decontx_tag, cohort_label, target_genes, wanted_donors,
                           draw_time_fn):
    """Shared loader for the single-combined-matrix cohorts (Liu/Cillo-shaped
    build): counts.mtx + <tag>_decontx_counts.mtx + cells.csv + genes.csv in
    one directory, donor_id already matching canonical_donor_labels.csv
    directly (no per-cohort id transform needed, unlike Liu).
    draw_time_fn(cells_df_subset) -> (draw_time array, draw_time_precision array).
    """
    genes = pd.read_csv(f'{RESULTS}/{dir_name}/genes.csv')['gene'].tolist()
    cells = pd.read_csv(f'{RESULTS}/{dir_name}/cells.csv', low_memory=False)
    keep = cells['donor_id'].isin(wanted_donors) & cells['harmonized_lineage'].isin(WANTED_LINEAGES)
    if keep.sum() == 0:
        return None
    idx = np.where(keep.to_numpy())[0]
    raw = read_mtx_columns(f'{RESULTS}/{dir_name}/counts.mtx', idx)
    decontx = read_mtx_columns(f'{RESULTS}/{dir_name}/{decontx_tag}_decontx_counts.mtx', idx)
    raw_sub = subset_to_genes(raw, genes, target_genes)
    decontx_sub = subset_to_genes(decontx, genes, target_genes)

    meta = cells.iloc[idx][['cell_id', 'channel_id', 'donor_id', 'harmonized_lineage']].copy()
    meta['cohort'] = cohort_label
    draw_time, draw_time_precision = draw_time_fn(cells.iloc[idx])
    meta['draw_time'] = draw_time
    meta['draw_time_precision'] = draw_time_precision
    print(f'  {cohort_label}: {len(idx)} Platelet/Monocyte cells from {meta.donor_id.nunique()} '
          f'donors ({meta.harmonized_lineage.value_counts().to_dict()})', flush=True)
    return raw_sub, decontx_sub, meta


def qiu_draw_time(cells_sub):
    # cells.csv already carries draw_time (T0/T6), exact per-GSM, no ambiguity.
    return cells_sub['draw_time'].to_numpy(), np.full(len(cells_sub), 'exact_draw_label')


def chongqing_draw_time(cells_sub):
    # One GSM = one donor = one draw; no timepoint variation exists to record.
    return np.full(len(cells_sub), np.nan), np.full(len(cells_sub), 'single_draw_no_variation')


def stephenson_draw_time(cells_sub):
    # cells.csv already carries draw_time (Days_from_onset, numeric), exact per-cell.
    dt = pd.to_numeric(cells_sub['draw_time'], errors='coerce').to_numpy()
    return dt, np.full(len(cells_sub), 'exact_day')


def wilk_draw_time(cells_sub):
    # cells.csv already carries draw_time (day suffix where present, else NaN).
    dt = pd.to_numeric(cells_sub['draw_time'], errors='coerce')
    precision = np.where(dt.notna(), 'exact_day', 'unresolved_no_day_token')
    return dt.to_numpy(), precision


def process_cillo(target_genes, wanted_donors, cillo_day_map):
    genes = pd.read_csv(f'{RESULTS}/cillo_decontx_input/genes.csv')['gene'].tolist()
    cells = pd.read_csv(f'{RESULTS}/cillo_decontx_input/cells.csv')
    keep = cells['donor_id'].isin(wanted_donors) & cells['harmonized_lineage'].isin(WANTED_LINEAGES)
    if keep.sum() == 0:
        return None
    raw = mmread(f'{RESULTS}/cillo_decontx_input/counts.mtx').tocsc()
    decontx = mmread(f'{RESULTS}/cillo_decontx_input/cillo_decontx_counts.mtx').tocsc()
    idx = np.where(keep.to_numpy())[0]
    raw_sub = subset_to_genes(raw[:, idx], genes, target_genes)
    decontx_sub = subset_to_genes(decontx[:, idx], genes, target_genes)

    meta = cells.iloc[idx][['cell_id', 'channel_id', 'donor_id', 'harmonized_lineage']].copy()
    meta['cohort'] = 'Cillo 2021'
    meta['draw_time'] = meta['cell_id'].map(cillo_day_map)
    meta['draw_time_precision'] = np.where(meta['draw_time'].notna(), 'exact_day', 'unresolved_no_day_token')
    print(f'  Cillo: {len(idx)} Platelet/Monocyte cells from {meta.donor_id.nunique()} donors '
          f'({meta.harmonized_lineage.value_counts().to_dict()})', flush=True)
    return raw_sub, decontx_sub, meta


COHORT_PROC_FNS = {
    'Liu 2021': lambda tg, wd, ctx: process_liu(tg, wd),
    'Ren 2021': lambda tg, wd, ctx: process_ren(tg, wd),
    'COMBAT 2022': lambda tg, wd, ctx: process_combat(tg, wd, ctx['onset_map']),
    'Cillo 2021': lambda tg, wd, ctx: process_cillo(tg, wd, ctx['cillo_day_map']),
    'Qiu 2021': lambda tg, wd, ctx: process_simple_cohort(
        'qiu_decontx_input', 'qiu', 'Qiu 2021', tg, wd, qiu_draw_time),
    'Chongqing sepsis 2025': lambda tg, wd, ctx: process_simple_cohort(
        'chongqing_decontx_input', 'chongqing', 'Chongqing sepsis 2025', tg, wd, chongqing_draw_time),
    'Stephenson 2021': lambda tg, wd, ctx: process_simple_cohort(
        'stephenson_decontx_input', 'stephenson', 'Stephenson 2021', tg, wd, stephenson_draw_time),
    'Wilk 2021 JEM': lambda tg, wd, ctx: process_simple_cohort(
        'wilk_decontx_input', 'wilk', 'Wilk 2021 JEM', tg, wd, wilk_draw_time),
}


def cohort_tag(cohort_name):
    return re.sub(r'[^A-Za-z0-9]+', '_', cohort_name).strip('_').lower()


def extract_one_cohort(cohort_name, target_genes, labels, ckpt_dir, ctx):
    """Checkpointed per-cohort extraction: writes raw.mtx/decontx.mtx/cells.csv
    to ckpt_dir/<tag>/ and skips entirely if already present. Each cohort is
    independent -- a timeout or failure on one (e.g. Chongqing/Stephenson's
    multi-GB matrices) costs only that cohort's redo, not the whole run, and
    each can be submitted as its own short/batch job rather than one
    monolithic multi-hour job (see feedback_hpcc_checkpoint_short_partition
    memory -- this is the direct fix for the lesson learned there)."""
    tag = cohort_tag(cohort_name)
    out_dir = f'{ckpt_dir}/{tag}'
    marker = f'{out_dir}/DONE'
    if os.path.exists(marker):
        print(f'{cohort_name}: checkpoint already exists, skipping', flush=True)
        return
    os.makedirs(out_dir, exist_ok=True)
    pd.DataFrame({'gene': target_genes}).to_csv(f'{out_dir}/genes.csv', index=False)

    sub = labels[labels['cohort'] == cohort_name]
    wanted_donors = set(sub['donor_id'])
    event_map = dict(zip(sub['donor_id'], sub['event']))
    edef_map = dict(zip(sub['donor_id'], sub['event_definition']))
    disease_map = dict(zip(sub['donor_id'], sub['disease']))
    print(f'\n{cohort_name}: {len(wanted_donors)} canonical donors', flush=True)

    result = COHORT_PROC_FNS[cohort_name](target_genes, wanted_donors, ctx)
    if result is None:
        print(f'  {cohort_name}: NO Platelet/Monocyte cells found', flush=True)
        mmwrite(f'{out_dir}/raw.mtx', sparse.csc_matrix((len(target_genes), 0)))
        mmwrite(f'{out_dir}/decontx.mtx', sparse.csc_matrix((len(target_genes), 0)))
        pd.DataFrame(columns=['cell_id', 'cohort', 'channel_id', 'donor_id', 'harmonized_lineage',
                               'disease', 'event', 'event_definition', 'draw_time',
                               'draw_time_precision']).to_csv(f'{out_dir}/cells.csv', index=False)
    else:
        raw_sub, decontx_sub, meta = result
        meta['event'] = meta['donor_id'].map(event_map)
        meta['event_definition'] = meta['donor_id'].map(edef_map)
        meta['disease'] = meta['donor_id'].map(disease_map)
        meta = meta[['cell_id', 'cohort', 'channel_id', 'donor_id', 'harmonized_lineage',
                     'disease', 'event', 'event_definition', 'draw_time', 'draw_time_precision']]
        # tmp names already end in ".mtx" -- scipy.io.mmwrite appends ".mtx"
        # to any target path that doesn't already have it, which would
        # otherwise silently write to a different filename than expected.
        tmp_raw, tmp_decontx, tmp_cells = (f'{out_dir}/_tmp_raw.mtx', f'{out_dir}/_tmp_decontx.mtx',
                                            f'{out_dir}/_tmp_cells.csv')
        mmwrite(tmp_raw, raw_sub)
        mmwrite(tmp_decontx, decontx_sub)
        meta.to_csv(tmp_cells, index=False)
        os.rename(tmp_raw, f'{out_dir}/raw.mtx')
        os.rename(tmp_decontx, f'{out_dir}/decontx.mtx')
        os.rename(tmp_cells, f'{out_dir}/cells.csv')
    with open(marker, 'w') as f:
        f.write('done\n')
    print(f'{cohort_name}: checkpoint written to {out_dir}', flush=True)


def combine_checkpoints(ckpt_dir, out_dir):
    os.makedirs(f'{out_dir}/raw', exist_ok=True)
    os.makedirs(f'{out_dir}/decontx', exist_ok=True)
    raw_blocks, decontx_blocks, meta_blocks, target_genes = [], [], [], None
    for cohort_name in COHORT_PROC_FNS:
        tag = cohort_tag(cohort_name)
        cdir = f'{ckpt_dir}/{tag}'
        if not os.path.exists(f'{cdir}/DONE'):
            raise SystemExit(f'{cohort_name}: checkpoint not done yet ({cdir}) -- run extraction first')
        genes = pd.read_csv(f'{cdir}/genes.csv')['gene'].tolist()
        if target_genes is None:
            target_genes = genes
        else:
            assert genes == target_genes, f'{cohort_name}: gene list mismatch across checkpoints'
        cells = pd.read_csv(f'{cdir}/cells.csv', low_memory=False)
        if len(cells) == 0:
            continue
        raw_blocks.append(mmread(f'{cdir}/raw.mtx').tocsc())
        decontx_blocks.append(mmread(f'{cdir}/decontx.mtx').tocsc())
        meta_blocks.append(cells)
        print(f'{cohort_name}: {len(cells)} cells loaded from checkpoint', flush=True)

    raw_combined = sparse.hstack(raw_blocks, format='csc')
    decontx_combined = sparse.hstack(decontx_blocks, format='csc')
    meta_combined = pd.concat(meta_blocks, ignore_index=True)
    assert raw_combined.shape == decontx_combined.shape
    assert raw_combined.shape[1] == len(meta_combined)

    pd.DataFrame({'gene': target_genes}).to_csv(f'{out_dir}/genes.csv', index=False)
    mmwrite(f'{out_dir}/raw/counts.mtx', raw_combined)
    mmwrite(f'{out_dir}/decontx/counts.mtx', decontx_combined)
    meta_combined.to_csv(f'{out_dir}/raw/cells.csv', index=False)
    meta_combined.to_csv(f'{out_dir}/decontx/cells.csv', index=False)
    return raw_combined, meta_combined


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ckpt-dir', required=True, help='per-cohort checkpoint directory')
    ap.add_argument('--out-dir', help='final combined output directory (only used with --combine)')
    ap.add_argument('--cohort', help='extract only this one cohort (checkpointed, skip if done)')
    ap.add_argument('--combine', action='store_true', help='combine all cohort checkpoints into --out-dir')
    args = ap.parse_args()
    os.makedirs(args.ckpt_dir, exist_ok=True)

    if args.combine:
        assert args.out_dir, '--out-dir is required with --combine'
        raw_combined, meta_combined = combine_checkpoints(args.ckpt_dir, args.out_dir)
        print(f'\nDone. {raw_combined.shape[1]} cells x {raw_combined.shape[0]} genes.')
        print(f'Donors: {meta_combined.donor_id.nunique()}, '
              f'events: {meta_combined.groupby("donor_id")["event"].max().sum()}')
        print('\nBy cohort x lineage:')
        print(meta_combined.groupby(['cohort', 'harmonized_lineage'])['cell_id'].count())
        print('\nBy cohort, donors and events:')
        print(meta_combined.groupby('cohort')['donor_id'].nunique())
        print(meta_combined.groupby('cohort').apply(
            lambda d: d.groupby('donor_id')['event'].max().sum(), include_groups=False))
        print('\ndraw_time_precision counts:')
        print(meta_combined.groupby(['cohort', 'draw_time_precision'])['cell_id'].count())
        return

    print('Loading labels and gene intersection...', flush=True)
    labels = pd.read_csv(LABELS_CSV)
    target_genes = load_intersection_genes()
    print(f'  {len(target_genes)} shared genes across all 8 cohorts', flush=True)

    ctx = {}
    if args.cohort in (None, 'COMBAT 2022'):
        print('Loading COMBAT TimeSinceOnset map...', flush=True)
        ctx['onset_map'] = load_combat_onset_map()
    if args.cohort in (None, 'Cillo 2021'):
        print('Loading Cillo per-cell draw-day map...', flush=True)
        ctx['cillo_day_map'] = load_cillo_day_map()

    cohorts = [args.cohort] if args.cohort else list(COHORT_PROC_FNS)
    for cohort_name in cohorts:
        assert cohort_name in COHORT_PROC_FNS, f'unknown cohort {cohort_name!r}'
        extract_one_cohort(cohort_name, target_genes, labels, args.ckpt_dir, ctx)


if __name__ == '__main__':
    main()
