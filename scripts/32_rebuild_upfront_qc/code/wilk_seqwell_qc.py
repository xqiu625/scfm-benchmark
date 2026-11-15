"""Wilk Seq-Well-adapted QC, FROZEN rules (qc_checklist 6e):
 1. NO counts>=500 / genes>=200 -- kept as DIAGNOSTIC COLUMNS only.
 2. log1p counts & genes: LOWER-tail 5 MAD only (no two-sided deletion).
 3. thresholds estimated per Seq-Well x broad lineage, so naturally low-RNA
    platelets are not declared outliers by other cell types' distributions.
 4. within capture/sample where the stratum is large enough; small strata borrow
    the cohort-level same-lineage threshold -- a tiny sample never fits its own MAD.
 5. mito: NOT global 20%; lineage upper-tail 5 MAD, per-cell flags kept.
 6. donor-level platelet retention reported, never used to reverse-tune.
"""
import numpy as np, pandas as pd, os
from scipy.stats import median_abs_deviation

R='/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc'
OUT=f'{R}/wilk_seqwell_qc'; os.makedirs(OUT, exist_ok=True)
N_MADS=5; MIN_STRATUM=70

qc=pd.read_csv(f'{R}/gate1_qc/wilk/qc_flags_per_cell.csv', low_memory=False)
print(f'{len(qc)} cells, {qc.channel_id.nunique()} captures, {qc.harmonized_lineage.nunique()} lineages')

def lower_thr(v):
    v=v[np.isfinite(v)]
    if len(v)<MIN_STRATUM: return np.nan
    m=np.median(v); s=median_abs_deviation(v, scale=1.0)
    return np.nan if s==0 else m-N_MADS*s

def upper_thr(v):
    v=v[np.isfinite(v)]
    if len(v)<MIN_STRATUM: return np.nan
    m=np.median(v); s=median_abs_deviation(v, scale=1.0)
    return np.nan if s==0 else m+N_MADS*s

# cohort-level per-lineage fallbacks (rule 4)
fb={}
for lin,g in qc.groupby('harmonized_lineage'):
    fb[lin]={'counts':lower_thr(g.log1p_total_counts.to_numpy()),
             'genes':lower_thr(g.log1p_n_genes_by_counts.to_numpy()),
             'mito':upper_thr(g.pct_counts_mt.to_numpy())}

rows=[]; qc['thr_source']='capture_x_lineage'
for c in ['sw_thr_counts','sw_thr_genes','sw_thr_mito']: qc[c]=np.nan
for (ch,lin),g in qc.groupby(['channel_id','harmonized_lineage']):
    idx=g.index
    tc=lower_thr(g.log1p_total_counts.to_numpy())
    tg=lower_thr(g.log1p_n_genes_by_counts.to_numpy())
    tm=upper_thr(g.pct_counts_mt.to_numpy())
    src='capture_x_lineage'
    if not np.isfinite(tc) or not np.isfinite(tg) or not np.isfinite(tm):
        tc=fb[lin]['counts'] if not np.isfinite(tc) else tc
        tg=fb[lin]['genes'] if not np.isfinite(tg) else tg
        tm=fb[lin]['mito'] if not np.isfinite(tm) else tm
        src='borrowed_cohort_lineage'
    qc.loc[idx,['sw_thr_counts','sw_thr_genes','sw_thr_mito']]=[tc,tg,tm]
    qc.loc[idx,'thr_source']=src
    rows.append({'channel_id':ch,'lineage':lin,'n':len(g),'thr_source':src,
                 'thr_log1p_counts':tc,'thr_log1p_genes':tg,'thr_mito':tm})

qc['sw_flag_counts_low']=qc.log1p_total_counts < qc.sw_thr_counts
qc['sw_flag_genes_low'] =qc.log1p_n_genes_by_counts < qc.sw_thr_genes
qc['sw_flag_mito_high'] =qc.pct_counts_mt > qc.sw_thr_mito
qc['sw_would_delete']=qc[['sw_flag_counts_low','sw_flag_genes_low','sw_flag_mito_high']].any(axis=1)
# 10x floors retained as DIAGNOSTIC ONLY
qc['diag_10x_counts500']=qc.total_counts<500
qc['diag_10x_genes200']=qc.n_genes_by_counts<200

pd.DataFrame(rows).to_csv(f'{OUT}/wilk_thresholds_per_capture_lineage.csv', index=False)
qc.to_csv(f'{OUT}/wilk_cells_seqwell_flags.csv', index=False)

print('\n=== Seq-Well adaptive vs the rejected 10x floors ===')
print(f'  Seq-Well rules would delete : {int(qc.sw_would_delete.sum())} ({100*qc.sw_would_delete.mean():.1f}%)')
print(f'  old 10x floors flagged      : counts<500 {100*qc.diag_10x_counts500.mean():.1f}%  '
      f'genes<200 {100*qc.diag_10x_genes200.mean():.1f}%  (DIAGNOSTIC ONLY)')
print(f'  threshold source: ' + str(qc.thr_source.value_counts().to_dict()))
print('\n=== by lineage ===')
t=qc.groupby('harmonized_lineage').agg(n=('sw_would_delete','size'),
                                        removed=('sw_would_delete','sum'))
t['pct_removed']=(100*t.removed/t.n).round(2)
t['old_10x_pct']=qc.groupby('harmonized_lineage').apply(
    lambda d: round(100*(d.diag_10x_counts500|d.diag_10x_genes200).mean(),2), include_groups=False)
print(t.to_string())

plt=qc[qc.harmonized_lineage=='Platelet']
d=plt.groupby('donor_id').agg(before=('sw_would_delete','size'), removed=('sw_would_delete','sum'))
d['after']=d.before-d.removed
print(f'\n=== donor-level platelet retention (reported, never reverse-tuned) ===')
for m in [20,50]:
    qb=int((d.before>=m).sum()); qa=int((d.after>=m).sum())
    print(f'  >= {m} platelets: before QC {qb} donors -> after {qa} donors (lost {qb-qa})')
d.to_csv(f'{OUT}/wilk_platelet_per_donor.csv')
