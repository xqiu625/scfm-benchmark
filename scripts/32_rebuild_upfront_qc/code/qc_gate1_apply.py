"""
GATE 1 application: compute the frozen basic-QC flags and APPLY them.

Rules frozen 2026-09-21 (qc_checklist.md §6c, scheme A) -- this script does not
choose thresholds, it executes the frozen ones and emits the same audit tables
the flags-only pass did, so every cohort has a traceable before/after record:
  min_genes >= 200, min_counts >= 500, pct_counts_mt < 20% (global, all lineages),
  MAD 5 unscaled per channel HIGH TAIL ONLY, min_cells per gene NOT applied,
  pct_counts_ribo / pct_counts_hb recorded but never filtered.

Handles both cohort layouts:
  --layout combined : one counts.mtx + cells.csv (Liu, Cillo, Qiu, Chongqing,
                      Stephenson, Wilk)
  --layout split    : split_per_channel/<channel>/{counts.mtx,cells.csv}
                      (Ren, COMBAT) -- MAD is naturally per channel here.

Output is the QC-retained FULL-LINEAGE cell set: Scrublet runs on this next
(§2 requires the channel's full population, not a pre-extracted lineage), and
DecontX is refit on it after that.

Usage:
    python qc_gate1_apply.py --cohort-dir <dir> --out-dir <dir> --layout combined|split
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread, mmwrite

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_compute_flags import (MIN_GENES, MIN_COUNTS, MT_MAX_PCT, N_MADS,
                               MAD_TAILS_THAT_DELETE, FLAG_COLS, mad_outlier_flags)


def gene_masks(genes):
    upper = np.asarray([g.upper() for g in genes], dtype=np.str_)
    mt = np.char.startswith(upper, 'MT-')
    ribo = np.char.startswith(upper, 'RPL') | np.char.startswith(upper, 'RPS')
    hb = np.array([g.startswith('HB') and not g.startswith('HBEGF') and not g.startswith('HBS1')
                    for g in upper])
    return mt, ribo, hb


def compute_flags(mat, cells, genes):
    """mat: genes x cells. Returns cells-with-flags dataframe. MAD is computed
    per channel_id within whatever cell set is passed in."""
    mt, ribo, hb = gene_masks(genes)
    total = np.asarray(mat.sum(axis=0)).ravel()
    n_genes = np.asarray((mat > 0).sum(axis=0)).ravel()
    denom = np.where(total > 0, total, np.nan)
    with np.errstate(divide='ignore', invalid='ignore'):
        pct_mt = 100.0 * np.asarray(mat[mt].sum(axis=0)).ravel() / denom
        pct_ribo = 100.0 * np.asarray(mat[ribo].sum(axis=0)).ravel() / denom
        pct_hb = 100.0 * np.asarray(mat[hb].sum(axis=0)).ravel() / denom

    qc = cells.copy()
    qc['total_counts'] = total
    qc['n_genes_by_counts'] = n_genes
    qc['pct_counts_mt'] = pct_mt
    qc['pct_counts_ribo'] = pct_ribo
    qc['pct_counts_hb'] = pct_hb
    qc['log1p_total_counts'] = np.log1p(total)
    qc['log1p_n_genes_by_counts'] = np.log1p(n_genes)
    qc['flag_min_genes'] = qc['n_genes_by_counts'] < MIN_GENES
    qc['flag_min_counts'] = qc['total_counts'] < MIN_COUNTS
    qc['flag_mito'] = qc['pct_counts_mt'] > MT_MAX_PCT
    for c in ['flag_mad_counts_low', 'flag_mad_counts_high',
              'flag_mad_genes_low', 'flag_mad_genes_high']:
        qc[c] = False
    qc['mad_status_counts'] = 'ok'
    qc['mad_status_genes'] = 'ok'
    mad_rows = []
    for ch, idx in qc.groupby('channel_id').groups.items():
        idx = np.asarray(idx)
        lo_c, hi_c, st_c, med_c, mad_c = mad_outlier_flags(
            qc.loc[idx, 'log1p_total_counts'].to_numpy(), N_MADS)
        lo_g, hi_g, st_g, med_g, mad_g = mad_outlier_flags(
            qc.loc[idx, 'log1p_n_genes_by_counts'].to_numpy(), N_MADS)
        qc.loc[idx, 'flag_mad_counts_low'] = lo_c
        qc.loc[idx, 'flag_mad_counts_high'] = hi_c
        qc.loc[idx, 'flag_mad_genes_low'] = lo_g
        qc.loc[idx, 'flag_mad_genes_high'] = hi_g
        qc.loc[idx, 'mad_status_counts'] = st_c
        qc.loc[idx, 'mad_status_genes'] = st_g
        mad_rows.append({'channel_id': ch, 'n_cells': len(idx),
                          'mad_status_counts': st_c, 'median_log1p_counts': med_c,
                          'mad_log1p_counts': mad_c,
                          'mad_status_genes': st_g, 'median_log1p_genes': med_g,
                          'mad_log1p_genes': mad_g})
    qc['n_flags'] = qc[FLAG_COLS].sum(axis=1)
    delete_cols = ['flag_min_genes', 'flag_min_counts', 'flag_mito'] + MAD_TAILS_THAT_DELETE
    qc['flag_would_delete'] = qc[delete_cols].any(axis=1)
    return qc, pd.DataFrame(mad_rows), delete_cols


def write_audit(qc, delete_cols, out_dir, tag=''):
    os.makedirs(out_dir, exist_ok=True)
    per_filter = pd.DataFrame({
        'filter': FLAG_COLS,
        'n_flagged': [int(qc[c].sum()) for c in FLAG_COLS],
        'pct_flagged': [round(100.0 * qc[c].mean(), 3) for c in FLAG_COLS]})
    per_filter.loc[len(per_filter)] = ['CUMULATIVE_would_delete',
                                        int(qc['flag_would_delete'].sum()),
                                        round(100.0 * qc['flag_would_delete'].mean(), 3)]
    per_filter.to_csv(f'{out_dir}/audit_per_filter{tag}.csv', index=False)
    by_lin = (qc.groupby('harmonized_lineage')
               .agg(n_in=('flag_would_delete', 'size'), n_removed=('flag_would_delete', 'sum')))
    by_lin['n_retained'] = by_lin.n_in - by_lin.n_removed
    by_lin['pct_removed'] = (100 * by_lin.n_removed / by_lin.n_in).round(2)
    by_lin.to_csv(f'{out_dir}/audit_by_lineage{tag}.csv')
    by_don = (qc.groupby(['donor_id', 'harmonized_lineage'])
               .agg(n_in=('flag_would_delete', 'size'), n_removed=('flag_would_delete', 'sum')))
    by_don['n_retained'] = by_don.n_in - by_don.n_removed
    by_don.reset_index().to_csv(f'{out_dir}/audit_by_donor_lineage{tag}.csv', index=False)
    return per_filter, by_lin


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cohort-dir', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--layout', required=True, choices=['combined', 'split'])
    ap.add_argument('--donor-map', default=None,
                     help='CSV with sampleID->donor, for cohorts whose per-channel cells.csv '
                          'carries no donor_id (Ren: channel_id is the sampleID, donor lives in '
                          'ren_gse158055_channel_inventory.csv). Donor-level auditing is required '
                          '(checklist 0.2), so this is joined rather than skipped.')
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    genes = pd.read_csv(f'{args.cohort_dir}/genes.csv')['gene'].astype(str).tolist()

    if args.layout == 'combined':
        if os.path.exists(f'{args.out_dir}/DONE'):
            print('already applied, skipping', flush=True)
            return
        print('Loading combined counts.mtx + cells.csv...', flush=True)
        mat = mmread(f'{args.cohort_dir}/counts.mtx').tocsc()
        cells = pd.read_csv(f'{args.cohort_dir}/cells.csv', low_memory=False)
        assert mat.shape[1] == len(cells)
        qc, mad_tbl, delete_cols = compute_flags(mat, cells, genes)
        per_filter, by_lin = write_audit(qc, delete_cols, args.out_dir)
        mad_tbl.to_csv(f'{args.out_dir}/audit_mad_per_channel.csv', index=False)
        keep = ~qc['flag_would_delete'].to_numpy()
        mmwrite(f'{args.out_dir}/counts.mtx', mat[:, keep])
        qc[keep].reset_index(drop=True).to_csv(f'{args.out_dir}/cells.csv', index=False)
        qc.to_csv(f'{args.out_dir}/qc_flags_per_cell.csv', index=False)
        pd.read_csv(f'{args.cohort_dir}/genes.csv').to_csv(f'{args.out_dir}/genes.csv', index=False)
        print(per_filter.to_string(index=False), flush=True)
        print(by_lin.to_string(), flush=True)
        print(f'\nRETAINED {int(keep.sum())} / {len(keep)} cells', flush=True)
    else:
        split_in = f'{args.cohort_dir}/split_per_channel'
        split_out = f'{args.out_dir}/split_per_channel'
        os.makedirs(split_out, exist_ok=True)
        channels = sorted(d for d in os.listdir(split_in) if os.path.isdir(f'{split_in}/{d}'))
        print(f'{len(channels)} channels', flush=True)
        all_qc = []
        for i, ch in enumerate(channels):
            co = f'{split_out}/{ch}'
            if os.path.exists(f'{co}/DONE'):
                all_qc.append(pd.read_csv(f'{co}/qc_flags_per_cell.csv', low_memory=False))
                continue
            os.makedirs(co, exist_ok=True)
            mat = mmread(f'{split_in}/{ch}/counts.mtx').tocsc()
            cells = pd.read_csv(f'{split_in}/{ch}/cells.csv', low_memory=False)
            qc, _mad, delete_cols = compute_flags(mat, cells, genes)
            keep = ~qc['flag_would_delete'].to_numpy()
            tmp = f'{co}/.counts.mtx'
            mmwrite(tmp, mat[:, keep])
            os.rename(tmp, f'{co}/counts.mtx')
            qc[keep].reset_index(drop=True).to_csv(f'{co}/cells.csv', index=False)
            qc.to_csv(f'{co}/qc_flags_per_cell.csv', index=False)
            with open(f'{co}/DONE', 'w') as f:
                f.write('done\n')
            all_qc.append(qc)
            if (i + 1) % 25 == 0 or i == len(channels) - 1:
                print(f'  {i+1}/{len(channels)} channels', flush=True)
        qc_all = pd.concat(all_qc, ignore_index=True)
        if 'donor_id' not in qc_all.columns:
            assert args.donor_map, ('per-channel cells.csv has no donor_id and no --donor-map '
                                     'given; donor-level auditing is required, not optional')
            inv = pd.read_csv(args.donor_map)
            m = dict(zip(inv['sampleID'].astype(str), inv['donor'].astype(str)))
            qc_all['donor_id'] = qc_all['channel_id'].astype(str).map(m)
            n_missing = int(qc_all['donor_id'].isna().sum())
            print(f'joined donor_id from {args.donor_map}: '
                  f'{qc_all.donor_id.nunique()} donors, {n_missing} cells unmapped', flush=True)
        delete_cols = ['flag_min_genes', 'flag_min_counts', 'flag_mito'] + MAD_TAILS_THAT_DELETE
        per_filter, by_lin = write_audit(qc_all, delete_cols, args.out_dir)
        pd.read_csv(f'{args.cohort_dir}/genes.csv').to_csv(f'{args.out_dir}/genes.csv', index=False)
        print(per_filter.to_string(index=False), flush=True)
        print(by_lin.to_string(), flush=True)
        print(f'\nRETAINED {int((~qc_all.flag_would_delete).sum())} / {len(qc_all)} cells', flush=True)

    with open(f'{args.out_dir}/DONE', 'w') as f:
        f.write('gate1 applied\n')


if __name__ == '__main__':
    main()
