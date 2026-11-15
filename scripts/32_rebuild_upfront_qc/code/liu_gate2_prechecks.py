"""Pre-Gate-2 checks for Liu (qc_checklist 6e):
 1. inspect the ONE lane whose inflection is itself an outlier
 2. donor-level target-lineage audit (outcome-blind)
 3. light Erythroid diagnostic -- documents why a PRE-CALL lineage label is not
    confirmed cell identity; never used to rescue barcodes or retune thresholds
"""
import numpy as np, pandas as pd, os
R='/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc'
OUT=f'{R}/liu_cellcall'
qc=pd.read_csv(f'{OUT}/liu_cells_with_callflag.csv', low_memory=False)
lane=pd.read_csv(f'{OUT}/liu_cellcall_per_lane.csv').set_index('channel_id')
rk=pd.read_csv(f'{OUT}/liu_barcode_rank_per_lane.csv').set_index('lane')
lane=lane.join(rk[['inflection_rank','knee_rank','n_barcodes']])

print('='*70); print('CHECK 1 -- lane(s) with an outlier inflection')
out=lane[lane.inflection_is_extreme]
print(f'flagged: {list(out.index)}')
infl=lane.inflection_counts
print(f'cohort inflection: median={infl.median():.0f} p5={np.percentile(infl,5):.0f} '
      f'p95={np.percentile(infl,95):.0f} min={infl.min():.0f} max={infl.max():.0f}')
for ln in out.index:
    r=lane.loc[ln]
    curve=np.load(f'{OUT}/rank_curve_{ln}.npy')
    ir=int(r.inflection_rank)
    print(f'\n  {ln}: inflection={r.inflection_counts:.0f} at rank {ir}/{len(curve)} '
          f'({100*ir/len(curve):.1f}% of barcodes); called={r.pct_called}%')
    print(f'    knee={r.knee_counts:.0f} valley={r.valley_counts:.0f}')
    # is there still a clear turn there? compare local slope around the point
    y=np.log10(curve+1); x=np.log10(np.arange(1,len(y)+1)); dy=np.gradient(y,x)
    w=max(50,len(curve)//100)
    lo=max(0,ir-w); hi=min(len(dy)-1,ir+w)
    print(f'    log-log slope: before={dy[lo]:.2f}  at={dy[ir-1]:.2f}  after={dy[hi]:.2f}'
          f'   (a clear turn shows a sharp local minimum)')
    print(f'    counts at rank 1/25%/50%/inflection/end: {curve[0]:.0f} / '
          f'{curve[len(curve)//4]:.0f} / {curve[len(curve)//2]:.0f} / {curve[ir-1]:.0f} / {curve[-1]:.0f}')
    batch=ln.split(chr(95))[0]
    peers=lane[[i.startswith(batch) for i in lane.index]]
    print(f'    same-batch peers ({batch}, n={len(peers)}): inflection median='
          f'{peers.inflection_counts.median():.0f}, called% median={peers.pct_called.median():.1f}')

print('\n'+'='*70); print('CHECK 2 -- donor-level target-lineage audit (outcome-blind)')
for lin in ['Platelet','Monocyte_Myeloid']:
    g=qc[qc.harmonized_lineage==lin]
    d=g.groupby('donor_id').agg(before=('called_cell','size'), called=('called_cell','sum'))
    d['pct']=(100*d.called/d.before).round(1)
    print(f'\n  {lin}: {len(d)} donors, called% median={d.pct.median():.1f} '
          f'min={d.pct.min():.1f} max={d.pct.max():.1f}')
    if lin=='Platelet':
        for m in [20,50]:
            qb=int((d.before>=m).sum()); qa=int((d.called>=m).sum())
            pushed=sorted(d[(d.before>=m)&(d.called<m)].index)
            print(f'    >={m} platelets: before {qb} donors -> after {qa} (lost {qb-qa}) {pushed[:10]}')
        print(f'    donors zeroed by calling: {int(((d.before>0)&(d.called==0)).sum())}')
        # concentration: is loss concentrated in few donors?
        loss=(d.before-d.called).sort_values(ascending=False)
        print(f'    loss concentration: top-5 donors hold '
              f'{100*loss.head(5).sum()/max(loss.sum(),1):.1f}% of all platelet loss')
    d.to_csv(f'{OUT}/liu_donor_{lin}_callcounts.csv')

print('\n'+'='*70); print('CHECK 3 -- Erythroid light diagnostic (documentation only)')
e=qc[qc.harmonized_lineage=='Erythroid']
for lbl,sub in [('called',e[e.called_cell]),('UNcalled',e[~e.called_cell])]:
    print(f'  Erythroid {lbl:8s} n={len(sub):6d}  counts med={sub.total_counts.median():7.0f}  '
          f'genes med={sub.n_genes_by_counts.median():6.0f}  mito med={sub.pct_counts_mt.median():5.1f}%')
u=qc[qc.harmonized_lineage=='Unclear']
for lbl,sub in [('called',u[u.called_cell]),('UNcalled',u[~u.called_cell])]:
    print(f'  Unclear   {lbl:8s} n={len(sub):6d}  counts med={sub.total_counts.median():7.0f}  '
          f'genes med={sub.n_genes_by_counts.median():6.0f}  mito med={sub.pct_counts_mt.median():5.1f}%')
print('\n  (pre-call lineage labels are assigned to CANDIDATE barcodes; ambient HBA/HBB is')
print('   abundant in blood, so empty droplets are readily labelled Erythroid, and background')
print('   barcodes without clear lineage signal fall into Unclear -- these two categories are')
print('   exactly where background barcodes concentrate, so a low called-rate is the expected')
print('   pattern, not evidence that real cells were deleted.)')
