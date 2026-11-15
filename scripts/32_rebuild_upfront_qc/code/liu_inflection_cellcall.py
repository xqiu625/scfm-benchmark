"""Liu cell calling, FROZEN rule (qc_checklist 6e): called cell = total_counts >=
that LANE'S OWN inflection. Knee/valley are saved as diagnostics and never delete.
Not EmptyDrops -- GSE161918 deposits no raw-droplet matrix, so this is
barcode-rank-based cell calling with a lane-specific inflection threshold."""
import numpy as np, pandas as pd, os

R='/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc'
OUT=f'{R}/liu_cellcall'; os.makedirs(OUT, exist_ok=True)

qc=pd.read_csv(f'{R}/gate1_qc/liu/qc_flags_per_cell.csv', low_memory=False)
rank=pd.read_csv(f'{OUT}/liu_barcode_rank_per_lane.csv').set_index('lane')

# lane-level outlier flagging on the INFLECTION itself (not on knee/valley
# disagreement -- those never trigger a threshold change)
infl=rank['inflection_counts']
med=infl.median(); mad=np.median(np.abs(infl-med))
rank['inflection_is_extreme']= (np.abs(infl-med) > 5*mad) if mad>0 else False
rank['knee_vs_inflection_ratio']=(rank.knee_counts/rank.inflection_counts).round(2)
rank['valley_vs_inflection_ratio']=(rank.valley_counts/rank.inflection_counts).round(2)
rank['aux_detector_disagrees']=((rank.knee_vs_inflection_ratio<0.5)|
                                 (rank.knee_vs_inflection_ratio>2)|
                                 (rank.valley_vs_inflection_ratio<0.5)|
                                 (rank.valley_vs_inflection_ratio>2))

qc['lane_inflection']=qc.channel_id.map(rank['inflection_counts'])
qc['called_cell']=qc.total_counts >= qc.lane_inflection

# QC profile of barcodes within +-20% of each lane's threshold (methods record)
band=qc[(qc.total_counts >= 0.8*qc.lane_inflection) & (qc.total_counts <= 1.2*qc.lane_inflection)]
print(f'barcodes within +-20% of their lane threshold: {len(band)} '
      f'({100*len(band)/len(qc):.1f}%)  median genes={band.n_genes_by_counts.median():.0f} '
      f'median mito={band.pct_counts_mt.median():.1f}%')

per_lane=qc.groupby('channel_id').agg(n_candidate=('called_cell','size'),
                                       n_called=('called_cell','sum'))
per_lane['pct_called']=(100*per_lane.n_called/per_lane.n_candidate).round(1)
per_lane=per_lane.join(rank[['inflection_counts','knee_counts','valley_counts',
                              'knee_vs_inflection_ratio','valley_vs_inflection_ratio',
                              'inflection_is_extreme','aux_detector_disagrees']])
per_lane.to_csv(f'{OUT}/liu_cellcall_per_lane.csv')
qc.to_csv(f'{OUT}/liu_cells_with_callflag.csv', index=False)

print(f'\n=== Liu cell calling (lane-specific inflection) ===')
print(f'  candidate barcodes (demuxlet SNG): {len(qc)}')
print(f'  called cells                     : {int(qc.called_cell.sum())} '
      f'({100*qc.called_cell.mean():.1f}%)')
print(f'  lanes flagged: inflection extreme={int(per_lane.inflection_is_extreme.sum())}  '
      f'aux-detector disagreement (diagnostic only)={int(per_lane.aux_detector_disagrees.sum())}')
print(f'  pct_called per lane: median={per_lane.pct_called.median():.1f} '
      f'min={per_lane.pct_called.min():.1f} max={per_lane.pct_called.max():.1f}')
print('\n  by lineage (called vs not):')
print(qc.groupby('harmonized_lineage').called_cell.agg(
    n='size', n_called='sum').assign(pct=lambda d:(100*d.n_called/d.n).round(1)).to_string())
