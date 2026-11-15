"""Comparison table for the Gate-1 mito freeze decision:
  scheme A = current global 20% cap
  scheme B = per-lineage MAD on pct_counts_mt (high tail only, unscaled, 5 MADs)
Reports actual thresholds, platelet flag counts, the 3 affected donors, and
the effect on OTHER lineages (B changes their thresholds too)."""
import numpy as np, pandas as pd
from scipy.stats import median_abs_deviation

D='/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/qc_flags_trial/cillo'
qc=pd.read_csv(f'{D}/qc_flags_per_cell.csv', low_memory=False)
MIN_PLT, N_MADS = 20, 5

# scheme B thresholds: per-lineage, pooled within cohort.
# NOTE the tension this creates -- every other MAD rule in this pipeline is
# per-CHANNEL (a within-batch relative measure). Per lineage x per channel is
# not feasible for platelets here: 831 platelets / 35 channels ~= 24 per
# channel, far below the 70-cell floor the MAD spec requires.
rows=[]
for lin, g in qc.groupby('harmonized_lineage'):
    v=g.pct_counts_mt.dropna().to_numpy()
    med=np.median(v); mad=median_abs_deviation(v, scale=1.0)
    thr=med+N_MADS*mad if mad>0 else np.inf
    rows.append({'lineage':lin,'n':len(g),'median_mito':round(med,2),'MAD':round(mad,2),
                 'thr_B_per_lineage_MAD':round(thr,1),
                 'n_flag_A_global20':int((g.pct_counts_mt>20).sum()),
                 'n_flag_B':int((g.pct_counts_mt>thr).sum())})
t=pd.DataFrame(rows).set_index('lineage')
t['pct_A']=(100*t.n_flag_A_global20/t.n).round(2); t['pct_B']=(100*t.n_flag_B/t.n).round(2)
pd.set_option('display.width',200)
print('=== thresholds & mito-flag counts by lineage ===')
print(t[['n','median_mito','MAD','thr_B_per_lineage_MAD','n_flag_A_global20','pct_A','n_flag_B','pct_B']].to_string())

# donor-level, platelets, combining mito rule with the OTHER frozen filters
OTHER=['flag_min_genes','flag_min_counts','flag_mad_counts_high','flag_mad_genes_high']
plt=qc[qc.harmonized_lineage=='Platelet'].copy()
thrB=t.loc['Platelet','thr_B_per_lineage_MAD']
plt['del_A']=plt[OTHER].any(axis=1)|(plt.pct_counts_mt>20)
plt['del_B']=plt[OTHER].any(axis=1)|(plt.pct_counts_mt>thrB)
d=plt.groupby('donor_id').agg(before=('del_A','size'),delA=('del_A','sum'),delB=('del_B','sum'))
d['after_A']=d.before-d.delA; d['after_B']=d.before-d.delB
qb=d.before>=MIN_PLT
print(f'\n=== donor impact (platelet minimum = {MIN_PLT}, scheme B platelet threshold = {thrB}%) ===')
print(f'qualified before QC        : {int(qb.sum())}')
print(f'qualified after A (20%)    : {int((d.after_A>=MIN_PLT).sum())}')
print(f'qualified after B (per-lin): {int((d.after_B>=MIN_PLT).sum())}')
print(f'pushed below by A          : {int((qb&(d.after_A<MIN_PLT)).sum())}  {sorted(d[qb&(d.after_A<MIN_PLT)].index)}')
print(f'pushed below by B          : {int((qb&(d.after_B<MIN_PLT)).sum())}  {sorted(d[qb&(d.after_B<MIN_PLT)].index)}')
print('\nthe three donors in question:')
print(d.loc[[i for i in ['COVID_2','COVID_3','COVID_12'] if i in d.index],
            ['before','delA','after_A','delB','after_B']].to_string())
print(f'\ntotal platelets retained: A={int(len(plt)-plt.del_A.sum())}  B={int(len(plt)-plt.del_B.sum())}  (of {len(plt)})')
