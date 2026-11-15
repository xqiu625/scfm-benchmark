"""Q1: are the QC metrics computed from RAW INTEGER counts (not normalized /
already-subset matrices)?  Q4: do Liu's removed cells look like low-quality or
empty barcodes, or like a normal population that simply sits lower?"""
import numpy as np, pandas as pd, os
from scipy.io import mmread
R='/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc'

print('=== Q1: are input matrices raw integer counts? ===')
for c,lay in [('qiu','c'),('cillo','c'),('liu','c'),('chongqing','c'),
              ('stephenson','c'),('wilk','c'),('combat','s'),('ren','s')]:
    try:
        if lay=='c':
            m=mmread(f'{R}/{c}_decontx_input/counts.mtx')
        else:
            d=f'{R}/{c}_decontx_input/split_per_channel'
            ch=sorted(os.listdir(d))[0]
            m=mmread(f'{d}/{ch}/counts.mtx')
        v=m.data[:200000]
        is_int=bool(np.all(np.equal(np.mod(v,1),0)))
        print(f'  {c:12s} dtype={m.dtype}  all_integer={is_int}  min={v.min():.3g} max={v.max():.3g} '
              f'shape={m.shape}')
    except Exception as e:
        print(f'  {c:12s} ERROR {e}')

print('\n=== Q4: Liu removed vs retained ===')
q=pd.read_csv(f'{R}/gate1_qc/liu/qc_flags_per_cell.csv', low_memory=False)
rm=q[q.flag_would_delete]; kp=q[~q.flag_would_delete]
for lbl,d in [('REMOVED',rm),('RETAINED',kp)]:
    print(f'  {lbl:9s} n={len(d):7d}  counts med={d.total_counts.median():7.0f} '
          f'p25={d.total_counts.quantile(.25):6.0f} p75={d.total_counts.quantile(.75):7.0f} | '
          f'genes med={d.n_genes_by_counts.median():5.0f} | mito med={d.pct_counts_mt.median():5.1f}')
# is the removed set a separable low tail, or overlapping?
print(f'  overlap check: removed p90 counts={rm.total_counts.quantile(.9):.0f}  '
      f'retained p10 counts={kp.total_counts.quantile(.10):.0f}')
print(f'  removed cells that are ONLY mito-flagged (otherwise fine): '
      f'{int((rm.flag_mito & ~rm.flag_min_counts & ~rm.flag_min_genes).sum())}')
print('\n  Liu barcode-layer accounting (candidate -> SNG -> QC-retained):')
print(f'    demuxlet-SNG barcodes entering Gate 1 : {len(q)}')
print(f'    QC-retained                            : {len(kp)} ({100*len(kp)/len(q):.1f}%)')

print('\n=== Q3 detail: Wilk whole-distribution vs tail ===')
w=pd.read_csv(f'{R}/gate1_qc/wilk/qc_flags_per_cell.csv', low_memory=False)
print(f'  Wilk counts percentiles: ' + ' '.join(
    f'p{p}={np.percentile(w.total_counts,p):.0f}' for p in [10,25,50,75,90,99]))
print(f'  Wilk genes  percentiles: ' + ' '.join(
    f'p{p}={np.percentile(w.n_genes_by_counts,p):.0f}' for p in [10,25,50,75,90,99]))
print(f'  fraction below the 10x floors: counts<500 {100*(w.total_counts<500).mean():.1f}%  '
      f'genes<200 {100*(w.n_genes_by_counts<200).mean():.1f}%')
