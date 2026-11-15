"""Do the cohorts ENTER this pipeline in the same QC state? A uniform absolute
threshold applied to non-uniformly-preprocessed inputs creates a systematic
cross-cohort difference in what counts as 'a cell' -- and cohort is the LOCO
split variable."""
import pandas as pd, numpy as np, os
R='/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/gate1_qc'
rows=[]
for c in ['qiu','cillo','liu','chongqing','stephenson','wilk','combat','ren']:
    f=f'{R}/{c}/qc_flags_per_cell.csv'
    if not os.path.exists(f):
        sub=[]
        d=f'{R}/{c}/split_per_channel'
        if os.path.isdir(d):
            for ch in sorted(os.listdir(d))[:40]:
                p=f'{d}/{ch}/qc_flags_per_cell.csv'
                if os.path.exists(p): sub.append(pd.read_csv(p, low_memory=False))
        if not sub: continue
        q=pd.concat(sub, ignore_index=True); note='(first 40 channels)'
    else:
        q=pd.read_csv(f, low_memory=False); note=''
    rows.append({'cohort':c+note,'n_cells':len(q),
                 'min_counts':int(q.total_counts.min()),
                 'p1_counts':int(np.percentile(q.total_counts,1)),
                 'median_counts':int(q.total_counts.median()),
                 'min_genes':int(q.n_genes_by_counts.min()),
                 'p1_genes':int(np.percentile(q.n_genes_by_counts,1)),
                 'median_genes':int(q.n_genes_by_counts.median()),
                 'max_mito':round(float(q.pct_counts_mt.max()),1),
                 'p99_mito':round(float(np.percentile(q.pct_counts_mt.dropna(),99)),1)})
t=pd.DataFrame(rows).set_index('cohort')
pd.set_option('display.width',250)
print(t.to_string())
