"""Streaming min/max over ALL channels for Ren and COMBAT -- required before the
claim 'this cohort shows upstream truncation' can be written for the whole
cohort rather than a 40-channel sample (checklist 6d Q2 caveat)."""
import numpy as np, pandas as pd, os
R='/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/gate1_qc'
for c in ['ren','combat']:
    d=f'{R}/{c}/split_per_channel'
    chans=sorted(x for x in os.listdir(d) if os.path.isdir(f'{d}/{x}'))
    mn_c=mn_g=np.inf; mx_m=-np.inf; n=0; nch=0
    per=[]
    for ch in chans:
        p=f'{d}/{ch}/qc_flags_per_cell.csv'
        if not os.path.exists(p): continue
        q=pd.read_csv(p, low_memory=False,
                      usecols=['total_counts','n_genes_by_counts','pct_counts_mt'])
        if not len(q): continue
        a,b = q.total_counts.min(), q.n_genes_by_counts.min()
        m = q.pct_counts_mt.max()
        mn_c=min(mn_c,a); mn_g=min(mn_g,b); mx_m=max(mx_m,m); n+=len(q); nch+=1
        per.append({'channel':ch,'min_counts':a,'min_genes':b,'max_mito':m})
    t=pd.DataFrame(per)
    print(f'=== {c}: ALL {nch} channels, {n} cells ===')
    print(f'  cohort-wide min_counts={mn_c}  min_genes={mn_g}  max_mito={mx_m:.2f}')
    print(f'  per-channel max_mito: median={t.max_mito.median():.2f} '
          f'p95={np.percentile(t.max_mito,95):.2f} max={t.max_mito.max():.2f}')
    print(f'  channels whose max_mito exceeds 12%: {int((t.max_mito>12).sum())}/{nch}')
    print(f'  per-channel min_genes: median={t.min_genes.median():.0f} max={t.min_genes.max():.0f}')
    t.to_csv(f'{R}/{c}/audit_full_channel_truncation.csv', index=False)
