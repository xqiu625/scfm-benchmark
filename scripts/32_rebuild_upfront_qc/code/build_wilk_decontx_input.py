"""
Assemble Wilk/GSE174072's counts.mtx / genes.csv / cells.csv for
run_decontx_per_channel.R, from convert_wilk_rds_to_mtx.R's per-GSM
matrix.mtx/barcodes.tsv/features.tsv output (Seq-Well, EXON counts only,
gene SYMBOL only -- no Ensembl ID, unlike every other cohort, so gene-panel
alignment across channels falls back to symbol intersection directly).

Donor/sample resolution is NOT a clean 1-GSM-per-canonical-row mapping the
way Chongqing/Qiu are -- canonical_donor_labels.csv's `notes` column carries
the actual sample-tag list per donor (e.g. "scRNA samples: 568_ACK;568_PBMC"),
grouped by IDENTICAL notes text shared across one or more donor rows:
  - PBMC-only policy (project_302006_rebuild_architecture): where both an
    "_ACK" (whole-blood/ACK-lysed) and "_PBMC" tag exist for the same
    person, keep only "_PBMC"; where ONLY "_ACK" exists, that donor has NO
    usable sample and is excluded entirely (e.g. 55650-0132d0/d7, both
    non-fatal, share a single ACK-only tag -- both excluded, matches the
    file inventory: only one GSM was ever staged for that person).
  - Two-donor-rows-per-notes-string case (longitudinal draws, e.g.
    "28205-0564d0"/"28205-0564d2" both citing "564_A;564_B"): the notes
    text does not say which file is which day. INFERRED (not confirmed)
    by pairing sorted donor_id (d0 before d2) with sorted sample tag (564_A
    before 564_B) -- both rows are non-fatal, so this inference costs no
    event-relevant precision if wrong, flagged as inferred in draw_time.
  - Plain single-tag notes (e.g. "Imm_056") map directly, no ambiguity.

draw_time: the "d<N>" suffix on some donor_ids (e.g. "...d0", "...d7") when
present; NaN otherwise (most Wilk donors are single, unspecified-day draws).

Usage:
    python build_wilk_decontx_input.py \
        --converted-dir <dir of GSM*/matrix.mtx+barcodes.tsv+features.tsv> \
        --labels <canonical_donor_labels.csv> \
        --out-dir <output dir>
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


def parse_notes(notes):
    m = re.search(r'scRNA samples:\s*(.+)', notes)
    if not m:
        return []
    return [t.strip() for t in m.group(1).split(';') if t.strip()]


def build_donor_sample_map(labels):
    wilk = labels[labels['cohort'] == 'Wilk 2021 JEM']
    donor_to_sample = {}
    donor_to_day = {}
    for notes, grp in wilk.groupby('notes'):
        tags = parse_notes(notes)
        pbmc_tags = sorted(t for t in tags if t.endswith('_PBMC'))
        plain_tags = sorted(t for t in tags if not t.endswith('_PBMC') and not t.endswith('_ACK'))
        usable_tags = pbmc_tags + plain_tags
        rows = sorted(grp['donor_id'].tolist())
        if len(usable_tags) == 0:
            print(f'  EXCLUDED (no usable PBMC tag): {rows} (notes={notes!r})', flush=True)
            continue
        if len(usable_tags) != len(rows):
            print(f'  WARNING: {len(usable_tags)} usable tags vs {len(rows)} donor rows for '
                  f'notes={notes!r} -- skipping this group entirely', flush=True)
            continue
        for donor_id, tag in zip(rows, usable_tags):
            donor_to_sample[donor_id] = tag
            m = re.search(r'd(\d+)$', donor_id)
            donor_to_day[donor_id] = int(m.group(1)) if m else np.nan
    return donor_to_sample, donor_to_day


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--converted-dir', required=True)
    ap.add_argument('--labels', required=True)
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    labels = pd.read_csv(args.labels)
    print('Resolving donor -> sample-tag mapping from canonical_donor_labels.csv notes...', flush=True)
    donor_to_sample, donor_to_day = build_donor_sample_map(labels)
    print(f'{len(donor_to_sample)} donors with a usable PBMC sample tag', flush=True)

    sample_to_gsm_dir = {}
    for d in glob.glob(os.path.join(args.converted_dir, 'GSM*')):
        base = os.path.basename(d)
        m = re.match(r'^(GSM\d+)_(.+)$', base)
        if m:
            sample_to_gsm_dir[m.group(2)] = d

    channel_blocks = []
    summary_rows = []
    for donor_id, tag in donor_to_sample.items():
        gsm_dir = sample_to_gsm_dir.get(tag)
        if gsm_dir is None:
            print(f'{donor_id}: sample tag {tag!r} has no converted GSM directory, skipping', flush=True)
            continue
        mtx_path = os.path.join(gsm_dir, 'matrix.mtx')
        if not os.path.exists(mtx_path):
            print(f'{donor_id}: {mtx_path} missing (conversion may have failed), skipping', flush=True)
            continue
        barcodes = pd.read_csv(os.path.join(gsm_dir, 'barcodes.tsv'), header=None)[0].to_numpy()
        genes = pd.read_csv(os.path.join(gsm_dir, 'features.tsv'), header=None)[0].to_numpy()
        mat = mmread(mtx_path).tocsc()  # genes x cells
        n_cells = mat.shape[1]

        adata = AnnData(X=mat.T.tocsr(), var=pd.DataFrame(index=genes))
        adata.var_names_make_unique()
        labels_z = label_channel(adata)

        channel_id = os.path.basename(gsm_dir)
        channel_blocks.append({'channel_id': channel_id, 'donor_id': donor_id,
                                'draw_time': donor_to_day[donor_id], 'genes': genes,
                                'mat': mat, 'barcodes': barcodes, 'labels': labels_z.to_numpy()})
        summary_rows.append({'channel_id': channel_id, 'donor_id': donor_id, 'n_cells': n_cells,
                              'n_genes_this_gsm': len(genes)})
        print(f'  {channel_id} ({donor_id}): {n_cells} cells, {len(genes)} genes', flush=True)

    if not channel_blocks:
        raise SystemExit('no channels resolved -- nothing to write')

    common_genes = set(channel_blocks[0]['genes'])
    for cb in channel_blocks[1:]:
        common_genes &= set(cb['genes'])
    ref_genes_raw = [g for g in channel_blocks[0]['genes'] if g in common_genes]
    ref_genes = make_unique_names(ref_genes_raw)
    print(f'\nCommon symbol gene panel across {len(channel_blocks)} channels: {len(ref_genes)} genes',
          flush=True)

    X_blocks = []
    cells_rows = []
    for cb in channel_blocks:
        pos = {g: i for i, g in enumerate(cb['genes'])}
        row_idx = np.array([pos[g] for g in ref_genes_raw])
        X_blocks.append(cb['mat'][row_idx, :])
        for bc, lineage in zip(cb['barcodes'], cb['labels']):
            cells_rows.append({'cell_id': f"{cb['channel_id']}_{bc}", 'channel_id': cb['channel_id'],
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
          f'{len(channel_blocks)} channels, {cells_df.donor_id.nunique()} donors.')
    print('Harmonized lineage composition:')
    print(cells_df['harmonized_lineage'].value_counts())


if __name__ == '__main__':
    main()
