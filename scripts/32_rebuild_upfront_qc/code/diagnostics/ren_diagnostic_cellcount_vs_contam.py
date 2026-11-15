"""
Diagnostic requested before trusting Ren's DecontX output (2026-09-16): Ren's
full-cohort mean contamination (52.7%) is much higher than Liu's (29.4%) and
structurally different (no lineage outlier, unlike Liu's Erythroid). First,
cheapest check - does contamination correlate with per-channel cell count,
consistent with small-channel EM instability (the per-channel isolation
workaround forced by R's 32-bit sparse-matrix ceiling)?

Result: no. Spearman rho=+0.223 (weak) with the WRONG SIGN for that
hypothesis - the smallest channels (<=100 cells) show the LOWEST
contamination, and the mean is flat from ~500 cells up to 18,000+ cells.
Ruled out small-N instability; see combine_comparison_channels.py +
run_decontx_per_channel.R for the follow-up controlled comparison that
confirms per-channel isolation gives identical results to celda's own
multi-batch call, and project_302006_rebuild_architecture memory for the
full writeup.
"""
import pandas as pd
from scipy import stats

df = pd.read_csv('/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/ren_decontx_output/ren_decontx_contamination_scores.csv')
per_ch = df.groupby('channel_id')['decontX_contamination'].agg(['mean','median']).reset_index()
counts = df.groupby('channel_id').size().rename('n_cells').reset_index()
per_ch = per_ch.merge(counts, on='channel_id')

rho, p = stats.spearmanr(per_ch['n_cells'], per_ch['mean'])
print(f'Spearman rho (n_cells vs mean contamination) = {rho:.4f}, p = {p:.3e}, n = {len(per_ch)}')

rho_log, p_log = stats.spearmanr(per_ch['n_cells'].apply(lambda x: x), per_ch['mean'])
print()
print('Cell-count bucket vs mean contamination:')
per_ch['bucket'] = pd.cut(per_ch['n_cells'], bins=[0,100,500,1000,3000,6000,10000,100000],
                            labels=['<=100','101-500','501-1000','1001-3000','3001-6000','6001-10000','>10000'])
print(per_ch.groupby('bucket', observed=True)['mean'].agg(['mean','median','count']))

print()
print('Channels with <=100 cells, sorted by contamination:')
print(per_ch[per_ch['n_cells']<=100].sort_values('mean', ascending=False)[['channel_id','n_cells','mean']].to_string(index=False))

per_ch.to_csv('/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/ren_decontx_output/ren_cellcount_vs_contamination.csv', index=False)
