import pandas as pd, numpy as np
D='/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/qc_flags_trial/cillo'
qc=pd.read_csv(f'{D}/qc_flags_per_cell.csv', low_memory=False)
plt=qc[qc.harmonized_lineage=='Platelet']; mono=qc[qc.harmonized_lineage=='Monocyte_Myeloid']
print('pct_counts_mt percentiles:')
qs=[50,75,90,95,97.5,99]
print(f'{"pct":>6} {"Platelet":>10} {"Monocyte":>10}')
for q in qs:
    print(f'{q:>6} {np.percentile(plt.pct_counts_mt.dropna(),q):>10.1f} {np.percentile(mono.pct_counts_mt.dropna(),q):>10.1f}')
print(f'\nplatelets above 20% mito: {(plt.pct_counts_mt>20).sum()} / {len(plt)} = {100*(plt.pct_counts_mt>20).mean():.1f}%')
print(f'monocytes above 20% mito: {(mono.pct_counts_mt>20).sum()} / {len(mono)} = {100*(mono.pct_counts_mt>20).mean():.1f}%')
# do the high-mito platelets look damaged (low complexity) or just platelet-like?
hi=plt[plt.pct_counts_mt>20]; lo=plt[plt.pct_counts_mt<=20]
print('\nhigh-mito(>20) vs low-mito platelets -- are the flagged ones actually degraded?')
for c in ['total_counts','n_genes_by_counts']:
    print(f'  {c:20s} high-mito median={hi[c].median():8.1f}   low-mito median={lo[c].median():8.1f}')
print(f'  also flagged by a count/gene filter: '
      f'{int((hi.flag_min_genes|hi.flag_min_counts|hi.flag_mad_counts_high|hi.flag_mad_genes_high).sum())} / {len(hi)}')
