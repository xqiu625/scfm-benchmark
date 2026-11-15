"""
Assemble Qiu/GSE167363's counts.mtx / genes.csv / cells.csv for
run_decontx_per_channel.R. One GSM = one library = one (donor, timepoint)
draw -- no demultiplexing needed. 12 GSMs across 7 donors (2 healthy
controls, single draw each; 5 septic patients, T0 + T6 draws each -- "0 and
6 h post sepsis diagnosis", per benchmark_design_decisions.md).

donor_id resolution is NOT the tar filename's sample name (e.g. "P25_T0" --
an internal lab code that does not appear anywhere in
canonical_donor_labels.csv) -- it is GEO's own `donor`/`donor id` SAMPLE
CHARACTERISTIC field, parsed from the series SOFT metadata into
qiu_gse167363_channel_inventory.csv (e.g. GSM5102902 -> "Tissue donor 3";
GSM5511351 -> "donor 5" -- two different characteristic keys across the
cohort's two GEO submission batches, "donor:" vs "donor id:", verified
directly against the deposited family.soft.gz). draw_time (T0/T6) comes from
the same inventory, parsed from the `sample name` characteristic -- exact,
no pooling ambiguity (1 GSM = 1 specific draw).

Usage:
    python build_qiu_decontx_input.py \
        --raw-dir <dir of GSM*.{barcodes,features,matrix}.tsv.gz> \
        --inventory <qiu_gse167363_channel_inventory.csv> \
        --out-dir <output dir>

Output (matches run_decontx_per_channel.R's expected IO_DIR contract):
    <out_dir>/counts.mtx, genes.csv,
    cells.csv (cell_id, channel_id, donor_id, harmonized_lineage, draw_time),
    build_summary.csv
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
    ap.add_argument('--inventory', required=True)
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    inv = pd.read_csv(args.inventory).set_index('gsm')
    matrix_files = sorted(glob.glob(os.path.join(args.raw_dir, '*_matrix.mtx.gz')))
    print(f'{len(matrix_files)} GSM matrix files found', flush=True)

    channel_blocks = []
    summary_rows = []
    for mtx_path in matrix_files:
        base = os.path.basename(mtx_path)[:-len('_matrix.mtx.gz')]
        m = re.match(r'^(GSM\d+)_(.+)$', base)
        assert m, f'unexpected filename: {mtx_path}'
        gsm, sample_tag = m.group(1), m.group(2)
        if gsm not in inv.index:
            print(f'{gsm}: not in inventory, skipping', flush=True)
            continue
        donor_id = inv.loc[gsm, 'donor']
        draw_time = inv.loc[gsm, 'draw_time']

        bc_path = os.path.join(args.raw_dir, f'{gsm}_{sample_tag}_barcodes.tsv.gz')
        feat_path = os.path.join(args.raw_dir, f'{gsm}_{sample_tag}_features.tsv.gz')
        barcodes = pd.read_csv(bc_path, header=None)[0].to_numpy()
        features = pd.read_csv(feat_path, sep='\t', header=None)
        ensembl = features[0].to_numpy()
        symbols = features[1].to_numpy()
        mat = mmread(mtx_path).tocsc()  # genes x cells
        n_cells = mat.shape[1]

        adata = AnnData(X=mat.T.tocsr(), var=pd.DataFrame(index=symbols))
        adata.var_names_make_unique()
        labels_z = label_channel(adata)

        channel_blocks.append({'gsm': gsm, 'donor_id': donor_id, 'draw_time': draw_time,
                                'ensembl': ensembl, 'symbols': symbols, 'mat': mat,
                                'barcodes': barcodes, 'labels': labels_z.to_numpy()})
        summary_rows.append({'gsm': gsm, 'donor_id': donor_id, 'draw_time': draw_time,
                              'n_cells': n_cells, 'n_genes_this_gsm': len(ensembl)})
        print(f'  {gsm} ({donor_id}, {draw_time}): {n_cells} cells', flush=True)

    if not channel_blocks:
        raise SystemExit('no GSMs matched the inventory -- nothing to write')

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
                                'donor_id': cb['donor_id'], 'harmonized_lineage': lineage,
                                'draw_time': cb['draw_time']})

    combined = sparse.hstack(X_blocks, format='csr')
    cells_df = pd.DataFrame(cells_rows)
    assert combined.shape == (len(ref_genes), len(cells_df))

    mmwrite(os.path.join(args.out_dir, 'counts.mtx'), combined)
    pd.DataFrame({'gene': ref_genes}).to_csv(os.path.join(args.out_dir, 'genes.csv'), index=False)
    cells_df.to_csv(os.path.join(args.out_dir, 'cells.csv'), index=False)
    pd.DataFrame(summary_rows).to_csv(os.path.join(args.out_dir, 'build_summary.csv'), index=False)

    print(f'\nDone. {combined.shape[1]} cells x {combined.shape[0]} genes across '
          f'{len(channel_blocks)} GSMs, {cells_df.donor_id.nunique()} donors.')
    print('Harmonized lineage composition:')
    print(cells_df['harmonized_lineage'].value_counts())


if __name__ == '__main__':
    main()
