"""
Liu/GSE161918 demuxlet singlet filter -- Phase 1 preprocessing step, upstream of
run_decontx_per_channel.R (see project_302006_rebuild_architecture memory,
"Liu-specific correction"). demuxlet's .best files are NOT singlet-filtered:
each 30,000-barcode lane carries a real cross-donor doublet load (median 42.5%
DBL, 27-61% range across lanes, verified exact against
GSE161918_demuxletResults.tar.gz). Drop DBL/AMB on demuxlet's genotype
evidence here, before Scrublet ever sees the lane -- Scrublet only needs to
find same-donor doublets in what's left, not redo demuxlet's job on the full
~40% cross-donor rate from scratch.

Real .best column format (verified directly against
GSE161918_demuxletResults.tar.gz, B1_10xlane1_demuxed.best, 2026-09-15):
tab-separated, 20 columns, INT_ID/BARCODE/NUM.SNPS/NUM.READS/DROPLET.TYPE/
BEST.GUESS/.../SNG.BEST.GUESS/.../DBL.BEST.GUESS/.... The columns actually
used here:
  BARCODE          -- 10x cell barcode, e.g. "AAACCTGAGAAACCTA-1"
  DROPLET.TYPE     -- SNG / DBL / AMB
  SNG.BEST.GUESS   -- the singlet's donor ID, e.g. "HGR0000083" -- NOT
                      BEST.GUESS, which for a DBL row is a comma-separated
                      donor pair with a mixing fraction (e.g.
                      "HGR0000098,HGR0000083,0.50") and is meaningless for a
                      SNG row's donor identity in a form other scripts can use.

Usage:
    python filter_liu_demuxlet_singlets.py --best-dir <dir of *_demuxed.best files> \
        --out-dir <output dir>

Output, per lane <lane>: <out_dir>/<lane>_singlets.csv with columns
  barcode, cell_id, channel_id, donor_id
cell_id is channel-prefixed (f"{channel_id}:{barcode}") since raw 10x
barcodes repeat with the same -1 suffix across every lane and would collide
once lanes are combined into one matrix for run_decontx_per_channel.R.

Plus one combined summary CSV (<out_dir>/liu_demuxlet_filter_summary.csv)
across all lanes: lane, n_barcodes, n_sng, n_dbl, n_amb, sng_frac -- meant to
be diffed against
scripts/32_rebuild_upfront_qc/channel_inventories/liu_gse161918_channel_inventory.csv's
n_singlet/n_doublet/n_ambiguous columns as a correctness check (they should
match exactly; that inventory was itself derived from these same .best files).
"""
import argparse
import glob
import os

import pandas as pd


def parse_best_file(path):
    """Read one *_demuxed.best file, return the full DataFrame (all droplet
    types) with just the columns this pipeline needs, renamed to lowercase."""
    df = pd.read_csv(path, sep='\t', usecols=['BARCODE', 'DROPLET.TYPE', 'SNG.BEST.GUESS'])
    df = df.rename(columns={'BARCODE': 'barcode', 'DROPLET.TYPE': 'droplet_type',
                             'SNG.BEST.GUESS': 'donor_id'})
    return df


def lane_name_from_path(path):
    # "B1_10xlane1_demuxed.best" -> "B1_10xlane1"
    base = os.path.basename(path)
    suffix = '_demuxed.best'
    assert base.endswith(suffix), f'unexpected filename, expected *_demuxed.best: {base}'
    return base[:-len(suffix)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--best-dir', required=True, help='Directory of *_demuxed.best files (from GSE161918_demuxletResults.tar.gz)')
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    best_files = sorted(glob.glob(os.path.join(args.best_dir, '*_demuxed.best')))
    if not best_files:
        raise SystemExit(f'no *_demuxed.best files found under {args.best_dir}')

    summary_rows = []
    for path in best_files:
        lane = lane_name_from_path(path)
        df = parse_best_file(path)
        n_barcodes = len(df)
        counts = df['droplet_type'].value_counts()
        n_sng = int(counts.get('SNG', 0))
        n_dbl = int(counts.get('DBL', 0))
        n_amb = int(counts.get('AMB', 0))
        assert n_sng + n_dbl + n_amb == n_barcodes, \
            f'{lane}: droplet_type has unexpected values beyond SNG/DBL/AMB: {counts.to_dict()}'

        singlets = df[df['droplet_type'] == 'SNG'].copy()
        singlets['channel_id'] = lane
        singlets['cell_id'] = singlets['channel_id'] + ':' + singlets['barcode']
        singlets = singlets[['barcode', 'cell_id', 'channel_id', 'donor_id']]
        singlets.to_csv(os.path.join(args.out_dir, f'{lane}_singlets.csv'), index=False)

        summary_rows.append({'lane': lane, 'n_barcodes': n_barcodes, 'n_sng': n_sng,
                              'n_dbl': n_dbl, 'n_amb': n_amb,
                              'sng_frac': round(n_sng / n_barcodes, 4)})
        print(f'{lane}: {n_barcodes} barcodes -> {n_sng} SNG kept '
              f'({n_dbl} DBL + {n_amb} AMB dropped)', flush=True)

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(os.path.join(args.out_dir, 'liu_demuxlet_filter_summary.csv'), index=False)
    print(f'\n{len(best_files)} lanes processed. '
          f'Total SNG kept: {summary["n_sng"].sum()} / {summary["n_barcodes"].sum()} barcodes '
          f'({summary["n_sng"].sum() / summary["n_barcodes"].sum():.1%})')


if __name__ == '__main__':
    main()
