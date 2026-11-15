"""
Flags-first QC gate (qc_checklist.md rev.3 §0). Computes EVERY per-cell QC flag
and emits the audit tables -- and deletes NOTHING. Deletion happens in a
separate step, only after the rules are reviewed and frozen.

Why flags-first: with several overlapping filters, per-filter flag counts and
cumulative removal counts diverge (a cell failing three filters is counted once
cumulatively but three times per-filter). Reporting only one of the two hides
which filter actually drives attrition. Both are emitted here, plus the overlap
structure.

Why donor-level auditing: the endpoint is PATIENT-level mortality, so what
matters is whether a filter removes patients, not only cells. A 10% overall
platelet loss concentrated in a few low-yield donors can eliminate those donors
entirely while the aggregate number looks harmless. Donors dropping to zero
retained platelets -- or below the prespecified per-donor minimum -- are listed
explicitly, not summarised away.

OUTCOME BLINDNESS: this script never reads event/mortality labels. Thresholds
are set blind to outcome and to predictive performance; the outcome-group
retention audit happens separately, AFTER the rules are frozen, as a check on
whether QC shifted patient composition -- never as a tuning signal.

Flags computed (none applied):
  flag_min_genes     n_genes_by_counts < MIN_GENES
  flag_min_counts    total_counts < MIN_COUNTS
  flag_mito          pct_counts_mt > MT_MAX_PCT
  flag_mad_counts    |log1p(total_counts) - median| > N_MADS * MAD      (per channel)
  flag_mad_genes     |log1p(n_genes_by_counts) - median| > N_MADS * MAD (per channel)
Recorded but never used as a filter (qc_checklist.md §1 #5/#6):
  pct_counts_ribo, pct_counts_hb

Usage:
    python qc_compute_flags.py --cohort-dir <dir with counts.mtx/cells.csv/genes.csv> \
        --out-dir <audit output dir> [--min-platelets-per-donor 20]
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import median_abs_deviation

# All PROVISIONAL candidates (qc_checklist.md rev.3 §6) -- frozen only after review.
MIN_GENES = 200
MIN_COUNTS = 500
MT_MAX_PCT = 20.0
N_MADS = 5

# Channels below this size do not get MAD flags at all: a median/MAD pair fitted
# on a handful of cells is not a meaningful distributional summary. Project rule,
# not a general boundary; such channels are recorded as mad_status='too_small'
# and their cells are NOT flagged (absence of evidence, not evidence of absence).
MIN_CELLS_FOR_MAD = 70

FLAG_COLS = ['flag_min_genes', 'flag_min_counts', 'flag_mito',
             'flag_mad_counts_low', 'flag_mad_counts_high',
             'flag_mad_genes_low', 'flag_mad_genes_high']

# Which MAD tails actually trigger REMOVAL once rules are frozen. Both tails are
# COMPUTED and reported; this states which ones are intended to delete. Low-tail
# count/gene outliers overlap heavily with the absolute floors (#2) and with
# genuinely small cells such as platelets, so the low tail is reported but NOT
# used for deletion by default -- the absolute floors handle that end explicitly
# and legibly. Revisit only with the audit in hand.
MAD_TAILS_THAT_DELETE = ['flag_mad_counts_high', 'flag_mad_genes_high']


def mad_outlier_flags(x, n_mads=N_MADS):
    """Two-sided MAD outlier flags on an already-log1p'd metric.

    Returns (low_flags, high_flags, status, med, mad).

    Explicit rules (qc_checklist.md rev.4 §2):
      * UNSCALED MAD -- scipy's `scale=1.0`, stated rather than left to a
        default, so the threshold means "5 raw MADs", not 5 normal-consistent
        sigma-equivalents (scale='normal' would divide by 0.67449 and silently
        make the same "5 MADs" ~1.48x tighter).
      * MAD == 0 -> NO cells flagged, status 'mad_zero'. This is the dangerous
        degenerate case: with MAD 0 the bounds collapse onto the median and
        every cell that differs at all would be flagged. A constant-ish metric
        means the statistic is uninformative here, not that the whole channel is
        an outlier.
      * NaN values -> never flagged (status records their count); a missing
        metric is not evidence of an outlier.
      * Tails kept separate so the audit can show which end drives attrition.
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    finite = np.isfinite(x)
    low = np.zeros(n, dtype=bool)
    high = np.zeros(n, dtype=bool)
    if finite.sum() < MIN_CELLS_FOR_MAD:
        return low, high, 'too_small', np.nan, np.nan
    med = float(np.median(x[finite]))
    mad = float(median_abs_deviation(x[finite], scale=1.0))
    if mad == 0:
        return low, high, 'mad_zero', med, mad
    low[finite] = x[finite] < med - n_mads * mad
    high[finite] = x[finite] > med + n_mads * mad
    return low, high, 'ok', med, mad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cohort-dir', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--min-platelets-per-donor', type=int, default=20,
                     help='prespecified per-donor platelet minimum used to define '
                          '"donor lost to QC" -- must be fixed BEFORE seeing the audit')
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print('Loading counts.mtx + cells.csv + genes.csv...', flush=True)
    mat = mmread(f'{args.cohort_dir}/counts.mtx').tocsc()  # genes x cells
    cells = pd.read_csv(f'{args.cohort_dir}/cells.csv', low_memory=False)
    genes = np.asarray(pd.read_csv(f'{args.cohort_dir}/genes.csv')['gene'].astype(str).tolist(),
                        dtype=np.str_)
    assert mat.shape[1] == len(cells) and mat.shape[0] == len(genes)
    print(f'  {mat.shape[1]} cells x {mat.shape[0]} genes, '
          f'{cells.channel_id.nunique()} channels, {cells.donor_id.nunique()} donors', flush=True)

    upper = np.char.upper(genes)  # genes is a true numpy str_ array, not object dtype
    mt_mask = np.char.startswith(upper, 'MT-')
    ribo_mask = np.char.startswith(upper, 'RPL') | np.char.startswith(upper, 'RPS')
    hb_mask = np.array([g.startswith('HB') and not g.startswith('HBEGF') and not g.startswith('HBS1')
                         for g in upper])
    print(f'  gene sets: {mt_mask.sum()} MT-, {ribo_mask.sum()} RPL/RPS, {hb_mask.sum()} HB*',
          flush=True)

    print('Computing per-cell QC metrics...', flush=True)
    total = np.asarray(mat.sum(axis=0)).ravel()
    n_genes = np.asarray((mat > 0).sum(axis=0)).ravel()
    with np.errstate(divide='ignore', invalid='ignore'):
        pct_mt = 100.0 * np.asarray(mat[mt_mask].sum(axis=0)).ravel() / np.where(total > 0, total, np.nan)
        pct_ribo = 100.0 * np.asarray(mat[ribo_mask].sum(axis=0)).ravel() / np.where(total > 0, total, np.nan)
        pct_hb = 100.0 * np.asarray(mat[hb_mask].sum(axis=0)).ravel() / np.where(total > 0, total, np.nan)

    qc = cells.copy()
    qc['total_counts'] = total
    qc['n_genes_by_counts'] = n_genes
    qc['pct_counts_mt'] = pct_mt
    qc['pct_counts_ribo'] = pct_ribo   # RECORD only, never a filter
    qc['pct_counts_hb'] = pct_hb       # RECORD only, never a filter
    qc['log1p_total_counts'] = np.log1p(total)
    qc['log1p_n_genes_by_counts'] = np.log1p(n_genes)

    qc['flag_min_genes'] = qc['n_genes_by_counts'] < MIN_GENES
    qc['flag_min_counts'] = qc['total_counts'] < MIN_COUNTS
    qc['flag_mito'] = qc['pct_counts_mt'] > MT_MAX_PCT

    # MAD is per channel by design -- it is a within-batch relative measure.
    for c in ['flag_mad_counts_low', 'flag_mad_counts_high',
              'flag_mad_genes_low', 'flag_mad_genes_high']:
        qc[c] = False
    qc['mad_status_counts'] = 'ok'
    qc['mad_status_genes'] = 'ok'
    mad_rows = []
    for ch, idx in qc.groupby('channel_id').groups.items():
        idx = np.asarray(idx)
        lo_c, hi_c, st_c, med_c, mad_c = mad_outlier_flags(
            qc.loc[idx, 'log1p_total_counts'].to_numpy())
        lo_g, hi_g, st_g, med_g, mad_g = mad_outlier_flags(
            qc.loc[idx, 'log1p_n_genes_by_counts'].to_numpy())
        qc.loc[idx, 'flag_mad_counts_low'] = lo_c
        qc.loc[idx, 'flag_mad_counts_high'] = hi_c
        qc.loc[idx, 'flag_mad_genes_low'] = lo_g
        qc.loc[idx, 'flag_mad_genes_high'] = hi_g
        qc.loc[idx, 'mad_status_counts'] = st_c
        qc.loc[idx, 'mad_status_genes'] = st_g
        mad_rows.append({'channel_id': ch, 'n_cells': len(idx),
                          'mad_status_counts': st_c, 'median_log1p_counts': med_c,
                          'mad_log1p_counts': mad_c, 'n_low_counts': int(lo_c.sum()),
                          'n_high_counts': int(hi_c.sum()),
                          'mad_status_genes': st_g, 'median_log1p_genes': med_g,
                          'mad_log1p_genes': mad_g, 'n_low_genes': int(lo_g.sum()),
                          'n_high_genes': int(hi_g.sum())})
    mad_per_channel = pd.DataFrame(mad_rows)

    qc['n_flags'] = qc[FLAG_COLS].sum(axis=1)
    qc['flag_any'] = qc['n_flags'] > 0
    # Deletion intent is explicit and separate from "was flagged": both MAD tails
    # are computed and audited, only the tails in MAD_TAILS_THAT_DELETE plus the
    # non-MAD filters are intended to remove cells once rules are frozen.
    delete_cols = ['flag_min_genes', 'flag_min_counts', 'flag_mito'] + MAD_TAILS_THAT_DELETE
    qc['flag_would_delete'] = qc[delete_cols].any(axis=1)

    # ---- audit tables (§0.1, §0.2) -------------------------------------------
    per_filter = pd.DataFrame({
        'filter': FLAG_COLS,
        'n_flagged': [int(qc[c].sum()) for c in FLAG_COLS],
        'pct_flagged': [round(100.0 * qc[c].mean(), 3) for c in FLAG_COLS],
    })
    per_filter.loc[len(per_filter)] = ['CUMULATIVE_any_flag', int(qc['flag_any'].sum()),
                                        round(100.0 * qc['flag_any'].mean(), 3)]
    per_filter.loc[len(per_filter)] = ['CUMULATIVE_would_delete',
                                        int(qc['flag_would_delete'].sum()),
                                        round(100.0 * qc['flag_would_delete'].mean(), 3)]
    print('\nPer-filter flags vs cumulative (overlap means these do not sum).', flush=True)
    print(f'Tails intended to DELETE: {delete_cols}', flush=True)
    print(per_filter.to_string(index=False), flush=True)

    overlap = qc['n_flags'].value_counts().sort_index().rename_axis('n_flags_failed').reset_index(
        name='n_cells')
    pair_rows = []
    for i, a in enumerate(FLAG_COLS):
        for b in FLAG_COLS[i + 1:]:
            pair_rows.append({'filter_a': a, 'filter_b': b,
                               'n_cells_failing_both': int((qc[a] & qc[b]).sum())})
    pairs = pd.DataFrame(pair_rows)

    by_lineage = (qc.groupby('harmonized_lineage')
                   .agg(n_in=('flag_would_delete', 'size'),
                        n_flagged=('flag_would_delete', 'sum'))
                   .assign(n_retained=lambda d: d.n_in - d.n_flagged,
                            pct_flagged=lambda d: (100 * d.n_flagged / d.n_in).round(2)))
    print('\nBy lineage (watch Platelet vs the rest):', flush=True)
    print(by_lineage.to_string(), flush=True)

    by_channel = (qc.groupby('channel_id')
                   .agg(n_in=('flag_would_delete', 'size'),
                        n_flagged=('flag_would_delete', 'sum'))
                   .assign(n_retained=lambda d: d.n_in - d.n_flagged,
                            pct_flagged=lambda d: (100 * d.n_flagged / d.n_in).round(2)))

    by_donor_lineage = (qc.groupby(['donor_id', 'harmonized_lineage'])
                         .agg(n_in=('flag_would_delete', 'size'),
                              n_flagged=('flag_would_delete', 'sum'))
                         .assign(n_retained=lambda d: d.n_in - d.n_flagged)
                         .reset_index())

    plt_ret = by_donor_lineage[by_donor_lineage['harmonized_lineage'] == 'Platelet'].copy()
    donors_all = set(qc['donor_id'].unique())
    donors_with_plt = set(plt_ret['donor_id'])
    zero_plt = sorted((donors_all - donors_with_plt) |
                       set(plt_ret.loc[plt_ret['n_retained'] == 0, 'donor_id']))
    below_min = sorted(plt_ret.loc[(plt_ret['n_retained'] > 0) &
                                    (plt_ret['n_retained'] < args.min_platelets_per_donor),
                                    'donor_id'])

    # CRITICAL distinction the raw counts above cannot make: a donor below the
    # minimum may have been below it BEFORE any filtering (a data-availability
    # limit of that cohort) or may have been pushed below it BY QC (attrition we
    # caused). Only the second is evidence against the thresholds. Cillo, for
    # instance, averages ~28 platelets/donor before any filter, so most of its
    # "below minimum" donors are pre-existing, not QC-induced -- reporting the
    # union of the two would overstate QC's damage substantially.
    m = args.min_platelets_per_donor
    plt_ret['was_below_before_qc'] = plt_ret['n_in'] < m
    plt_ret['is_below_after_qc'] = plt_ret['n_retained'] < m
    plt_ret['pushed_below_by_qc'] = (~plt_ret['was_below_before_qc']) & plt_ret['is_below_after_qc']
    plt_ret['zeroed_by_qc'] = (plt_ret['n_in'] > 0) & (plt_ret['n_retained'] == 0)
    pushed = sorted(plt_ret.loc[plt_ret['pushed_below_by_qc'], 'donor_id'])
    zeroed_by_qc = sorted(plt_ret.loc[plt_ret['zeroed_by_qc'], 'donor_id'])
    already_below = sorted(plt_ret.loc[plt_ret['was_below_before_qc'], 'donor_id'])
    print(f'\nDONOR-LEVEL IMPACT (endpoint is patient-level, so this is the number that matters):',
          flush=True)
    print(f'  ATTRIBUTABLE TO QC:', flush=True)
    print(f'    donors zeroed BY QC (had platelets, now none): {len(zeroed_by_qc)} {zeroed_by_qc[:20]}',
          flush=True)
    print(f'    donors pushed below min ({m}) BY QC (were at/above it before): '
          f'{len(pushed)} {pushed[:20]}', flush=True)
    print(f'  PRE-EXISTING (not caused by QC -- cohort data limit):', flush=True)
    print(f'    donors already below min before any filtering: {len(already_below)}', flush=True)
    print(f'  TOTALS after flags (union of both causes, do NOT read as QC damage):', flush=True)
    print(f'    zero retained platelets: {len(zero_plt)} | below min: {len(below_min)}', flush=True)

    qc.to_csv(f'{args.out_dir}/qc_flags_per_cell.csv', index=False)
    per_filter.to_csv(f'{args.out_dir}/audit_per_filter.csv', index=False)
    overlap.to_csv(f'{args.out_dir}/audit_flag_overlap.csv', index=False)
    pairs.to_csv(f'{args.out_dir}/audit_flag_pairs.csv', index=False)
    by_lineage.to_csv(f'{args.out_dir}/audit_by_lineage.csv')
    by_channel.to_csv(f'{args.out_dir}/audit_by_channel.csv')
    mad_per_channel.to_csv(f'{args.out_dir}/audit_mad_per_channel.csv', index=False)
    by_donor_lineage.to_csv(f'{args.out_dir}/audit_by_donor_lineage.csv', index=False)
    pd.DataFrame({'donor_id': zero_plt}).to_csv(
        f'{args.out_dir}/audit_donors_zero_platelets.csv', index=False)
    pd.DataFrame({'donor_id': below_min}).to_csv(
        f'{args.out_dir}/audit_donors_below_min_platelets.csv', index=False)
    plt_ret.to_csv(f'{args.out_dir}/audit_platelets_per_donor_cause.csv', index=False)
    pd.DataFrame([{'MIN_GENES': MIN_GENES, 'MIN_COUNTS': MIN_COUNTS, 'MT_MAX_PCT': MT_MAX_PCT,
                    'N_MADS': N_MADS,
                    'MIN_PLATELETS_PER_DONOR': args.min_platelets_per_donor}]).to_csv(
        f'{args.out_dir}/audit_thresholds_used.csv', index=False)

    print(f'\nFLAGS ONLY -- nothing deleted. Audit written to {args.out_dir}', flush=True)
    print('Next step is review and rule-freezing, not application.', flush=True)


if __name__ == '__main__':
    main()
