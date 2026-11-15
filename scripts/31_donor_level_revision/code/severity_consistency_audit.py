#!/usr/bin/env python3
"""Severity/disease/mortality consistency audit across all patient_id draws
(corrected mapping, all 11 source cohorts), requested by coordinator 2026-09-03."""
import h5py
import numpy as np
import pandas as pd


def rc(g):
    if isinstance(g, h5py.Group) and 'categories' in g:
        return g['categories'][:].astype(str)[g['codes'][:]]
    return g[:].astype(str)


H5AD = '/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c/state_covid_canonical.h5ad'
MAP = '/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/results/donor_meta_with_patient_id.csv'
OUT = '/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/results_patient/severity_consistency_audit.csv'
OUT_SUMMARY = '/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/results_patient/severity_consistency_summary.txt'

KEEP_DISEASE = ['COVID-19', 'Sepsis']
KEEP_SEV = ['fatal', 'mild', 'moderate', 'recovered', 'severe']

with h5py.File(H5AD, 'r') as h:
    sid = rc(h['obs/Sample_ID'])
    dis = rc(h['obs/Disease'])
    sev = rc(h['obs/covid_severity'])
    dno = rc(h['obs/Data_NO'])

cell_df = pd.DataFrame({'sample_id': sid, 'Disease': dis, 'covid_severity': sev, 'Data_NO': dno})
keep = cell_df['Disease'].isin(KEEP_DISEASE) & cell_df['covid_severity'].isin(KEEP_SEV)
cell_df = cell_df[keep]

n_multi = cell_df.groupby('sample_id')['covid_severity'].nunique()
print('sample_ids with >1 severity value across own cells (should be 0):', int((n_multi > 1).sum()))
n_multi_d = cell_df.groupby('sample_id')['Disease'].nunique()
print('sample_ids with >1 Disease value across own cells (should be 0):', int((n_multi_d > 1).sum()))

sample_level = cell_df.drop_duplicates(subset='sample_id')
pmap = pd.read_csv(MAP)[['sample_id', 'patient_id', 'y', 'Data_NO']]
merged = sample_level.merge(pmap, on='sample_id', how='inner', suffixes=('', '_map'))
print('merged rows:', len(merged), 'sample_level rows:', len(sample_level), 'pmap rows:', len(pmap))
assert (merged['Data_NO'] == merged['Data_NO_map']).all(), 'Data_NO mismatch after merge'

n_consistent = 0
n_inconsistent = 0
inconsistent_rows = []

for pid, g in merged.groupby('patient_id'):
    n_draws = len(g)
    if n_draws == 1:
        n_consistent += 1
        continue
    sevs = sorted(g['covid_severity'].unique())
    diss = sorted(g['Disease'].unique())
    severity_consistent = len(sevs) == 1
    disease_consistent = len(diss) == 1
    if severity_consistent and disease_consistent:
        n_consistent += 1
    else:
        n_inconsistent += 1
        for _, row in g.sort_values('sample_id').iterrows():
            inconsistent_rows.append({
                'patient_id': pid, 'n_draws': n_draws, 'sample_id': row['sample_id'],
                'Data_NO': row['Data_NO'], 'Disease': row['Disease'],
                'covid_severity': row['covid_severity'], 'y': row['y'],
                'severity_consistent': severity_consistent,
                'disease_consistent': disease_consistent,
            })

n_total_patients = merged['patient_id'].nunique()
n_multi_draw = int((merged.groupby('patient_id').size() > 1).sum())

summary = []
summary.append(f'TOTAL patients: {n_total_patients}')
summary.append(f'patients with a single draw: {n_total_patients - n_multi_draw}')
summary.append(f'patients with >1 draw: {n_multi_draw}')
summary.append(f'  of those, severity+Disease consistent across draws: {n_multi_draw - n_inconsistent}')
summary.append(f'  of those, INCONSISTENT (severity and/or Disease differ across draws): {n_inconsistent}')
summary.append('')
summary.append('Note: mortality y (fatal/not) is derived FROM covid_severity (severity=="fatal" -> y=1),')
summary.append('so a severity-consistent patient is automatically y-consistent too; y is reported')
summary.append('per-draw below only as a cross-check, not as an independent consistency dimension.')

print('\n'.join(summary))
with open(OUT_SUMMARY, 'w') as f:
    f.write('\n'.join(summary) + '\n')

if inconsistent_rows:
    df_out = pd.DataFrame(inconsistent_rows)
    df_out.to_csv(OUT, index=False)
    print(f'\nSAVED inconsistent-patient detail -> {OUT}')
    print(df_out.to_string(index=False))
else:
    print('\nNo inconsistent patients found -- no detail file written.')
