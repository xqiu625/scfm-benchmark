"""
Assemble Chongqing sepsis/GSE279452's counts.mtx / genes.csv / cells.csv for
run_decontx_per_channel.R. Structurally the simplest cohort so far: one GSM
= one library = one donor (no demultiplexing, no h5ad cross-referencing for
donor identity needed -- canonical_donor_labels.csv's donor_id IS the GSM
sample title suffix, e.g. "GSM8571091_Abd-S101.matrix.mtx.gz" -> donor_id
"Abd-S101", verified directly against canonical_donor_labels.csv's Chongqing
rows).

108 of 111 GSMs already extracted to raw gzipped CellRanger triples
(barcodes/features/matrix .tsv.gz, cellranger-6.1.1) at
01_DATA/raw/chongqing_sepsis_GSE279452/extracted/ -- 3 are missing from the
source archive (not an extraction bug, confirmed via extract_errors.log:
"Not found in archive"): GSM8571120/GSM8571121 (Abd-S8/Abd-S80, already
excluded per project_302006_rebuild_architecture -- ambiguous
metadata-to-matrix pairing per the CD8 companion preprint, both non-fatal)
and GSM8571198 (Skin-S152, a skin-tissue library, not in
canonical_donor_labels.csv anyway).

channel_id = GSM accession, donor_id = parsed sample title. Gene panels are
checked for identity across all GSMs (not assumed, per the Cillo precedent
of catching real per-GSM reference drift) via the same Ensembl-ID-
intersection two-pass approach as build_cillo_decontx_input.py.

Usage:
    python build_chongqing_decontx_input.py \
        --raw-dir <dir of GSM*.{barcodes,features,matrix}.tsv.gz> \
        --labels <canonical_donor_labels.csv> \
        --out-dir <output dir>

Output (matches run_decontx_per_channel.R's expected IO_DIR contract):
    <out_dir>/counts.mtx, genes.csv, cells.csv (cell_id, channel_id,
    donor_id, harmonized_lineage), build_summary.csv
"""
import argparse
import glob
import os
import re

import numpy as np
import pandas as pd
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--raw-dir', required=True)
    ap.add_argument('--labels', required=True)
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    labels = pd.read_csv(args.labels)
    wanted_donors = set(labels[labels['cohort'] == 'Chongqing sepsis 2025']['donor_id'])
    print(f'{len(wanted_donors)} canonical Chongqing donors', flush=True)

    matrix_files = sorted(glob.glob(os.path.join(args.raw_dir, '*.matrix.mtx.gz')))
    print(f'{len(matrix_files)} GSM matrix files found', flush=True)

    channel_blocks = []
    summary_rows = []
    for mtx_path in matrix_files:
        base = os.path.basename(mtx_path)[:-len('.matrix.mtx.gz')]
        m = re.match(r'^(GSM\d+)_(.+)$', base)
        assert m, f'unexpected filename: {mtx_path}'
        gsm, donor_id = m.group(1), m.group(2)

        if donor_id not in wanted_donors:
            summary_rows.append({'gsm': gsm, 'donor_id': donor_id, 'in_canonical_labels': False,
                                  'n_cells': 0})
            continue

        bc_path = os.path.join(args.raw_dir, f'{gsm}_{donor_id}.barcodes.tsv.gz')
        feat_path = os.path.join(args.raw_dir, f'{gsm}_{donor_id}.features.tsv.gz')
        barcodes = pd.read_csv(bc_path, header=None)[0].to_numpy()
        features = pd.read_csv(feat_path, sep='\t', header=None)
        ensembl = features[0].to_numpy()
        symbols = features[1].to_numpy()
        mat = mmread(mtx_path).tocsc()  # genes x cells
        n_cells = mat.shape[1]

        adata = AnnData(X=mat.T.tocsr(), var=pd.DataFrame(index=symbols))
        adata.var_names_make_unique()
        labels_z = label_channel(adata)

        channel_blocks.append({'gsm': gsm, 'donor_id': donor_id, 'ensembl': ensembl,
                                'symbols': symbols, 'mat': mat, 'barcodes': barcodes,
                                'labels': labels_z.to_numpy()})
        summary_rows.append({'gsm': gsm, 'donor_id': donor_id, 'in_canonical_labels': True,
                              'n_cells': n_cells, 'n_genes_this_gsm': len(ensembl)})
        if len(channel_blocks) % 20 == 0:
            print(f'  {len(channel_blocks)} channels loaded ({gsm} {donor_id}: {n_cells} cells)',
                  flush=True)

    if not channel_blocks:
        raise SystemExit('no GSMs matched canonical labels -- nothing to write')

    common_ensembl = set(channel_blocks[0]['ensembl'])
    for cb in channel_blocks[1:]:
        common_ensembl &= set(cb['ensembl'])
    ref_ensembl = [e for e in channel_blocks[0]['ensembl'] if e in common_ensembl]
    ref_symbols_raw = {e: s for e, s in zip(channel_blocks[0]['ensembl'], channel_blocks[0]['symbols'])}
    ref_genes = make_unique_names([ref_symbols_raw[e] for e in ref_ensembl])
    n_dropped = len(channel_blocks[0]['ensembl']) - len(ref_ensembl)
    print(f'\nCommon Ensembl-ID gene panel across {len(channel_blocks)} channels: '
          f'{len(ref_ensembl)} genes (dropped {n_dropped} not shared by every channel)', flush=True)

    X_blocks = []
    cells_rows = []
    for cb in channel_blocks:
        pos = {e: i for i, e in enumerate(cb['ensembl'])}
        row_idx = np.array([pos[e] for e in ref_ensembl])
        X_blocks.append(cb['mat'][row_idx, :])
        for bc, lineage in zip(cb['barcodes'], cb['labels']):
            cells_rows.append({'cell_id': f"{cb['gsm']}_{bc}", 'channel_id': cb['gsm'],
                                'donor_id': cb['donor_id'], 'harmonized_lineage': lineage})

    combined = sparse.hstack(X_blocks, format='csr')
    cells_df = pd.DataFrame(cells_rows)
    assert combined.shape == (len(ref_genes), len(cells_df))

    mmwrite(os.path.join(args.out_dir, 'counts.mtx'), combined)
    pd.DataFrame({'gene': ref_genes}).to_csv(os.path.join(args.out_dir, 'genes.csv'), index=False)
    cells_df.to_csv(os.path.join(args.out_dir, 'cells.csv'), index=False)
    pd.DataFrame(summary_rows).to_csv(os.path.join(args.out_dir, 'build_summary.csv'), index=False)

    print(f'\nDone. {combined.shape[1]} cells x {combined.shape[0]} genes across '
          f'{len(channel_blocks)} GSMs, {cells_df.donor_id.nunique()} donors '
          f'(of {len(wanted_donors)} canonical).')
    print('Harmonized lineage composition:')
    print(cells_df['harmonized_lineage'].value_counts())


if __name__ == '__main__':
    main()
