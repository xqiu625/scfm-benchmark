"""
Liu cell calling -- FALLBACK path (qc_checklist 6e).

Confirmed 2026-09-21: GSE161918 deposits 56 RNA + 56 CSP
*_filtered_feature_bc_matrix.h5 and NO raw/unfiltered droplet matrix (114
files, 0 raw matches). The lower ambient tail was already truncated upstream,
so EmptyDrops cannot be described here as standard raw-droplet cell calling.
Primary caller is therefore the per-lane barcode-rank curve; EmptyDrops, if run
later, is a sensitivity analysis only.

Explicitly forbidden (checklist): repackaging fixed 500-count / 200-gene floors
as cell calling.

Per lane: barcode-rank curve, knee and inflection, and a data-derived boundary
taken from the empirical valley between the two count components (the bimodal
split already visible in Q4: removed p90 = 530 vs retained p10 = 1,694).
Everything needed for the methods record is saved per lane.
"""
import numpy as np, pandas as pd, os
from scipy.io import mmread

R = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc'
OUT = f'{R}/liu_cellcall'
os.makedirs(OUT, exist_ok=True)


def knee_inflection(totals_sorted):
    """Knee / inflection on the log-log barcode-rank curve, as in the standard
    DropletUtils construction: knee = point of max curvature (steepest drop of
    the fitted signed curvature), inflection = point of steepest slope."""
    y = np.log10(totals_sorted + 1)
    x = np.log10(np.arange(1, len(y) + 1))
    dy = np.gradient(y, x)
    d2y = np.gradient(dy, x)
    inflection_i = int(np.argmin(dy))            # steepest descent
    knee_i = int(np.argmax(np.abs(d2y)))         # max curvature
    return knee_i, inflection_i


def valley_boundary(totals):
    """Data-derived boundary: the minimum-density point of the log10-count
    histogram BETWEEN the two modes -- i.e. the empirical valley, not a fixed
    threshold imported from elsewhere."""
    lv = np.log10(totals[totals > 0])
    hist, edges = np.histogram(lv, bins=120)
    centers = 0.5 * (edges[:-1] + edges[1:])
    peak1 = int(np.argmax(hist[:len(hist) // 2]))
    peak2 = len(hist) // 2 + int(np.argmax(hist[len(hist) // 2:]))
    if peak2 <= peak1:
        return np.nan, np.nan
    seg = hist[peak1:peak2]
    valley_i = peak1 + int(np.argmin(seg))
    return float(10 ** centers[valley_i]), float(centers[valley_i])


qc = pd.read_csv(f'{R}/gate1_qc/liu/qc_flags_per_cell.csv', low_memory=False)
rows = []
for lane, g in qc.groupby('channel_id'):
    tot = np.sort(g.total_counts.to_numpy())[::-1]
    knee_i, infl_i = knee_inflection(tot)
    v_count, v_log = valley_boundary(g.total_counts.to_numpy())
    rows.append({'lane': lane, 'n_barcodes': len(g),
                 'knee_rank': knee_i + 1, 'knee_counts': float(tot[knee_i]),
                 'inflection_rank': infl_i + 1, 'inflection_counts': float(tot[infl_i]),
                 'valley_counts': v_count,
                 'n_above_valley': int((g.total_counts > v_count).sum()) if np.isfinite(v_count) else np.nan,
                 'n_above_knee': int((g.total_counts > tot[knee_i]).sum()),
                 'n_above_inflection': int((g.total_counts > tot[infl_i]).sum())})
    np.save(f'{OUT}/rank_curve_{lane}.npy', tot)

t = pd.DataFrame(rows)
t['pct_above_valley'] = (100 * t.n_above_valley / t.n_barcodes).round(1)
t['pct_above_inflection'] = (100 * t.n_above_inflection / t.n_barcodes).round(1)
t.to_csv(f'{OUT}/liu_barcode_rank_per_lane.csv', index=False)
pd.set_option('display.width', 250)
print(t.head(12).to_string(index=False))
print('\n=== across 55 lanes ===')
for c in ['n_barcodes', 'knee_counts', 'inflection_counts', 'valley_counts',
          'pct_above_valley', 'pct_above_inflection']:
    print(f'  {c:22s} median={t[c].median():10.1f}  min={t[c].min():9.1f}  max={t[c].max():10.1f}')
print(f'\nSaved per-lane rank curves + table -> {OUT}')
