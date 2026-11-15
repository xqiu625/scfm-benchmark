"""Targeted platelet review for the Gate-1 freeze decision (Cillo).
Reports BOTH denominators, the 3 QC-affected donors' before/after counts,
and platelet-specific flag attribution -- which filter actually did it."""
import pandas as pd
import numpy as np

D = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/qc_flags_trial/cillo'
qc = pd.read_csv(f'{D}/qc_flags_per_cell.csv', low_memory=False)
MIN = 20

plt = qc[qc.harmonized_lineage == 'Platelet']
per_donor = plt.groupby('donor_id').agg(
    n_before=('flag_would_delete', 'size'),
    n_deleted=('flag_would_delete', 'sum'))
per_donor['n_after'] = per_donor.n_before - per_donor.n_deleted
per_donor['qualified_before'] = per_donor.n_before >= MIN
per_donor['qualified_after'] = per_donor.n_after >= MIN

all_donors = qc.donor_id.nunique()
qual_before = int(per_donor.qualified_before.sum())
qual_after = int(per_donor.qualified_after.sum())
pushed = per_donor[per_donor.qualified_before & ~per_donor.qualified_after]

print('=== DENOMINATORS (both required) ===')
print(f'donors total (cohort)                 : {all_donors}')
print(f'donors with >=20 platelets BEFORE QC  : {qual_before}')
print(f'donors with >=20 platelets AFTER QC   : {qual_after}')
print(f'QC-newly-failing / ALL donors         : {len(pushed)}/{all_donors} = {100*len(pushed)/all_donors:.1f}%')
print(f'QC-newly-failing / PREVIOUSLY QUALIFYING: {len(pushed)}/{qual_before} = {100*len(pushed)/qual_before:.1f}%')

print('\n=== THE AFFECTED DONORS: before -> after ===')
print(pushed[['n_before', 'n_deleted', 'n_after']].to_string())

FLAGS = ['flag_min_genes', 'flag_min_counts', 'flag_mito',
         'flag_mad_counts_high', 'flag_mad_genes_high']
print('\n=== PLATELET-SPECIFIC FLAG ATTRIBUTION (which filter did it) ===')
print('--- cohort-wide platelets ---')
attr = pd.DataFrame({
    'n_platelets_flagged': [int(plt[f].sum()) for f in FLAGS],
    'pct_of_platelets': [round(100*plt[f].mean(), 2) for f in FLAGS]}, index=FLAGS)
mono = qc[qc.harmonized_lineage == 'Monocyte_Myeloid']
attr['pct_of_monocytes'] = [round(100*mono[f].mean(), 2) for f in FLAGS]
attr['platelet_vs_mono_ratio'] = (attr.pct_of_platelets /
                                   attr.pct_of_monocytes.replace(0, np.nan)).round(1)
print(attr.to_string())

print('\n--- the 3 affected donors only ---')
sub = plt[plt.donor_id.isin(pushed.index)]
print(pd.DataFrame({
    'n_flagged': [int(sub[f].sum()) for f in FLAGS],
    'pct': [round(100*sub[f].mean(), 2) for f in FLAGS]}, index=FLAGS).to_string())
print('\nper-donor x filter:')
print(sub.groupby('donor_id')[FLAGS].sum().to_string())

print('\n=== platelet QC metric distributions (context for the review) ===')
for col in ['total_counts', 'n_genes_by_counts', 'pct_counts_mt']:
    print(f'{col:22s} platelet median={plt[col].median():9.1f}  monocyte median={mono[col].median():9.1f}')
