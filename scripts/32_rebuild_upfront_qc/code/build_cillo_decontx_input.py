"""
Assemble Cillo/GSE180578's counts.mtx / genes.csv / cells.csv for
run_decontx_per_channel.R (the single-combined-matrix + batch= driver, same
one first validated on Liu -- Cillo's scale, ~150-200k cells x 33,539 genes,
is nowhere near R/Matrix's 32-bit ceiling that forced Ren's per-channel-split
workaround, so no need for that heavier path here).

Two data sources per GSM library, joined on cell barcode:
  1. RAW counts: GSM<n>_..._GEX.tar.gz -> <sample>/{barcodes,features,matrix}.tsv
     (uncompressed CellRanger 3.1.0 output, full 33,539-gene panel). This is
     what DecontX runs on -- NOT the deposited annotated h5ad's X, which is a
     downstream, HVG-subsetted (2,000 genes), non-raw object (verified
     directly 2026-09-17: only 2,000 var, no integer-count layer/raw slot
     usable at full gene resolution).
  2. DONOR LABEL + membership: GSE180578_cillo_covid19_study_aggregrated_annotated_data.h5ad's
     obs (patient_id, outcome) -- this is the authors' own QC'd/doublet-called
     cell set (104,962 cells total, fewer than the raw per-GSM filtered
     matrices sum to), so a raw barcode only gets a trustworthy donor call if
     it also appears here. cell barcode join key: h5ad obs._index is
     "<10x_barcode_without_-1_suffix>_<h5ad_suffix>" where h5ad_suffix is the
     per-GSM integer pair in cillo_gse180578_channel_inventory.csv's
     h5ad_suffix column (e.g. GSM5465125 -> suffix "1_1", confirmed directly:
     that GSM's first 5 raw barcodes' h5ad-style keys all matched real h5ad
     obs rows with patient_id starting "COVID_1_day_1_...").

Only cells with a h5ad match are kept (same discipline as Liu's demuxlet-SNG
filter -- an unlabeled cell cannot contribute to a downstream mortality
model regardless of what DecontX does with it, so there is no reason to carry
it through channel construction).

donor_id is derived from h5ad patient_id by stripping the "_day_<N>_<PBMC|ETA>"
suffix when present (e.g. "COVID_11_day_1_PBMC" -> "COVID_11"); the two ALIR-
named patient_ids ("ALIR_5097_10_PBMC_TEST", "ALIR_5310_10") have no such
suffix and are kept as-is -- both forms match canonical_donor_labels.csv's
Cillo donor_id values exactly, verified against a sample of both patterns.

channel_id = GSM accession (the physical GEM well) -- not the internal ALIR-
style tar subdirectory name, which is not a stable/meaningful id outside this
one file.

Only cillo_pbmc_manifest.csv's include_pbmc_only==True rows (42 of 48 GSMs)
are processed -- the 6 excluded rows mix PBMC with endotracheal-aspirate
material in the same GEM well (project_302006_rebuild_architecture: PBMC-only
by design decision, different platelet recovery than PBMC).

Usage:
    python build_cillo_decontx_input.py \
        --raw-dir <dir of GSM*_GEX.tar.gz + the annotated h5ad> \
        --h5ad <path to the annotated h5ad, decompressed> \
        --manifest <cillo_pbmc_manifest.csv> \
        --inventory <cillo_gse180578_channel_inventory.csv> \
        --work-dir <scratch dir for tar extraction> \
        --out-dir <output dir>

Output (matches run_decontx_per_channel.R's expected IO_DIR contract):
    <out_dir>/counts.mtx   -- genes (rows) x cells (cols), MatrixMarket
    <out_dir>/genes.csv    -- one column: gene
    <out_dir>/cells.csv    -- cell_id, channel_id, donor_id, harmonized_lineage
    <out_dir>/build_summary.csv -- per-GSM cell counts at each filter stage
"""
import argparse
import glob
import os
import re
import shutil
import tarfile

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy import sparse
from scipy.io import mmread, mmwrite

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from harmonized_lineage_labels import label_channel


def make_unique_names(names):
    counts = {}
    out = []
    for n in names:
        counts[n] = counts.get(n, 0) + 1
        out.append(n if counts[n] == 1 else f'{n}-{counts[n] - 1}')
    return out


def donor_id_from_patient_id(pid):
    m = re.match(r'^(.*)_day_\d+_(PBMC|ETA)$', pid)
    return m.group(1) if m else pid


def load_h5ad_donor_map(h5ad_path):
    f = h5py.File(h5ad_path, 'r')
    idx = np.array([x.decode() if isinstance(x, bytes) else x for x in f['obs']['_index'][:]])
    pid = np.array([x.decode() if isinstance(x, bytes) else x for x in f['obs']['patient_id'][:]])
    donor = np.array([donor_id_from_patient_id(p) for p in pid])
    return dict(zip(idx, donor))


def read_raw_10x(tar_path, work_dir):
    gsm = os.path.basename(tar_path).split('_')[0]
    extract_dir = os.path.join(work_dir, gsm)
    if os.path.exists(extract_dir):
        shutil.rmtree(extract_dir)
    os.makedirs(extract_dir, exist_ok=True)
    with tarfile.open(tar_path) as tf:
        tf.extractall(extract_dir)
    # exactly one sample subdirectory per tar
    subdirs = [d for d in glob.glob(os.path.join(extract_dir, '*')) if os.path.isdir(d)]
    assert len(subdirs) == 1, f'{tar_path}: expected 1 subdir, found {len(subdirs)}'
    sdir = subdirs[0]
    barcodes = pd.read_csv(os.path.join(sdir, 'barcodes.tsv'), header=None)[0].to_numpy()
    features = pd.read_csv(os.path.join(sdir, 'features.tsv'), sep='\t', header=None)
    ensembl = features[0].to_numpy()
    symbols = features[1].to_numpy()
    mat = mmread(os.path.join(sdir, 'matrix.mtx')).tocsc()  # genes x cells
    shutil.rmtree(extract_dir)
    return barcodes, ensembl, symbols, mat


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--raw-dir', required=True)
    ap.add_argument('--h5ad', required=True)
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--inventory', required=True)
    ap.add_argument('--work-dir', required=True)
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.work_dir, exist_ok=True)

    manifest = pd.read_csv(args.manifest)
    manifest = manifest[manifest['include_pbmc_only'] == True]
    inventory = pd.read_csv(args.inventory).set_index('gsm')

    print(f'Loading h5ad donor map from {args.h5ad}...', flush=True)
    donor_map = load_h5ad_donor_map(args.h5ad)
    print(f'  {len(donor_map)} h5ad cells with a donor label', flush=True)

    # Pass 1: per-channel extraction + h5ad-donor matching + z-labeling.
    # Gene panels are NOT assumed identical across GSMs -- verified false
    # 2026-09-17: GSM5465126's features.tsv has 33,538 rows vs GSM5465121-125's
    # 33,539 (minor 10x reference-version drift across processing batches).
    # Each channel keeps its OWN native gene axis for z-labeling (label_channel
    # only needs its marker panel to be present, not the full reference to
    # match byte-for-byte); alignment to one shared axis happens in pass 2,
    # keyed on Ensembl ID (stable across a 1-gene reference diff), not symbol.
    channel_blocks = []  # list of dict: gsm, ensembl, symbols, mat, donor_ids, barcodes, labels
    summary_rows = []

    for _, row in manifest.iterrows():
        gsm = row['gsm_id']
        tar_path = os.path.join(args.raw_dir, row['filename'])
        if not os.path.exists(tar_path):
            print(f'{gsm}: tar file not found at {tar_path}, skipping', flush=True)
            continue
        suffix = inventory.loc[gsm, 'h5ad_suffix']
        if pd.isna(suffix) or suffix == '':
            print(f'{gsm}: no h5ad_suffix in inventory (no h5ad cells) -- skipping', flush=True)
            continue

        barcodes, ensembl, symbols, mat = read_raw_10x(tar_path, args.work_dir)
        n_raw = len(barcodes)

        h5ad_keys = np.array([bc.split('-')[0] + '_' + str(suffix) for bc in barcodes])
        donor_ids = np.array([donor_map.get(k) for k in h5ad_keys])
        keep_mask = donor_ids != None  # noqa: E711
        n_matched = int(keep_mask.sum())

        if n_matched == 0:
            print(f'{gsm}: 0/{n_raw} barcodes matched the h5ad -- skipping', flush=True)
            summary_rows.append({'gsm': gsm, 'n_raw_barcodes': n_raw, 'n_h5ad_matched': 0})
            continue

        sub_mat = mat[:, keep_mask].tocsr()  # genes x cells, this channel's own gene axis
        sub_donor = donor_ids[keep_mask]
        sub_barcodes = barcodes[keep_mask]

        adata = AnnData(X=sub_mat.T.tocsr(), var=pd.DataFrame(index=symbols))  # cells x genes
        adata.var_names_make_unique()
        labels = label_channel(adata)

        channel_blocks.append({'gsm': gsm, 'ensembl': ensembl, 'symbols': symbols, 'mat': sub_mat,
                                'donor_ids': sub_donor, 'barcodes': sub_barcodes,
                                'labels': labels.to_numpy()})
        summary_rows.append({'gsm': gsm, 'n_raw_barcodes': n_raw, 'n_h5ad_matched': n_matched,
                              'n_genes_this_gsm': len(ensembl)})
        print(f'{gsm}: raw={n_raw}, h5ad-matched={n_matched}, genes={len(ensembl)}, '
              f'donors={len(set(sub_donor))}', flush=True)

    if not channel_blocks:
        raise SystemExit('no GSMs produced any cells -- nothing to write')

    # Pass 2: align every channel to the Ensembl-ID intersection across all
    # channels (order = first channel's order, filtered to the intersection),
    # then hstack.
    common_ensembl = set(channel_blocks[0]['ensembl'])
    for cb in channel_blocks[1:]:
        common_ensembl &= set(cb['ensembl'])
    ref_ensembl = [e for e in channel_blocks[0]['ensembl'] if e in common_ensembl]
    ref_symbols_raw = {e: s for e, s in zip(channel_blocks[0]['ensembl'], channel_blocks[0]['symbols'])}
    ref_genes = make_unique_names([ref_symbols_raw[e] for e in ref_ensembl])
    n_dropped = len(channel_blocks[0]['ensembl']) - len(ref_ensembl)
    print(f'\nCommon Ensembl-ID gene panel across {len(channel_blocks)} channels: '
          f'{len(ref_ensembl)} genes (dropped {n_dropped} from the first channel\'s '
          f'{len(channel_blocks[0]["ensembl"])}-gene panel not shared by every channel)', flush=True)

    X_blocks = []
    cells_rows = []
    for cb in channel_blocks:
        pos = {e: i for i, e in enumerate(cb['ensembl'])}
        row_idx = np.array([pos[e] for e in ref_ensembl])
        aligned = cb['mat'][row_idx, :]
        X_blocks.append(aligned)
        for bc, donor, lineage in zip(cb['barcodes'], cb['donor_ids'], cb['labels']):
            cells_rows.append({'cell_id': f"{cb['gsm']}_{bc}", 'channel_id': cb['gsm'],
                                'donor_id': donor, 'harmonized_lineage': lineage})

    combined = sparse.hstack(X_blocks, format='csr')
    cells_df = pd.DataFrame(cells_rows)
    assert combined.shape == (len(ref_genes), len(cells_df)), \
        f'shape mismatch: matrix {combined.shape} vs genes={len(ref_genes)} cells={len(cells_df)}'

    mmwrite(os.path.join(args.out_dir, 'counts.mtx'), combined)
    pd.DataFrame({'gene': ref_genes}).to_csv(os.path.join(args.out_dir, 'genes.csv'), index=False)
    cells_df.to_csv(os.path.join(args.out_dir, 'cells.csv'), index=False)
    pd.DataFrame(summary_rows).to_csv(os.path.join(args.out_dir, 'build_summary.csv'), index=False)

    print(f'\nDone. {combined.shape[1]} cells x {combined.shape[0]} genes across '
          f'{len(summary_rows)} GSMs, {cells_df.donor_id.nunique()} donors.')
    print('Harmonized lineage composition:')
    print(cells_df['harmonized_lineage'].value_counts())


if __name__ == '__main__':
    main()
