#!/usr/bin/env python3
"""
Job H step 1: verify donor_meta_with_patient_id.csv against the canonical
STATE h5ad under the identical stage-27 filters (Disease in COVID-19/Sepsis,
covid_severity in fatal/mild/moderate/recovered/severe).

Checks: every post-filter Sample_ID appears in the mapping with exactly one
patient_id (no Sample_ID maps to >1 patient_id). Reports per-disease patient
and fatal-patient counts, using patient label = max(y) across a patient's draws.
"""
import h5py
import numpy as np
import pandas as pd
import sys

H5AD = '/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c/state_covid_canonical.h5ad'
MAP_CSV = '/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/results/donor_meta_with_patient_id.csv'

KEEP_DISEASE = ['COVID-19', 'Sepsis']
KEEP_SEVERITY = ['fatal', 'mild', 'moderate', 'recovered', 'severe']


def read_categorical(group):
    if isinstance(group, h5py.Group) and 'categories' in group:
        cats = group['categories'][:].astype(str)
        codes = group['codes'][:]
        return cats[codes]
    return group[:].astype(str)


def main():
    with h5py.File(H5AD, 'r') as h:
        disease = read_categorical(h['obs/Disease'])
        severity = read_categorical(h['obs/covid_severity'])
        sample_id = read_categorical(h['obs/Sample_ID'])

    keep = np.isin(disease, KEEP_DISEASE) & np.isin(severity, KEEP_SEVERITY)
    sid_f = sample_id[keep]
    dis_f = disease[keep]

    h5_sids = pd.DataFrame({'sample_id': sid_f, 'disease': dis_f}).drop_duplicates()
    print(f'post-filter unique Sample_IDs in h5ad: {h5_sids["sample_id"].nunique()}')

    m = pd.read_csv(MAP_CSV)
    print(f'mapping rows: {len(m)}  unique sample_id: {m["sample_id"].nunique()}  '
          f'unique patient_id: {m["patient_id"].nunique()}')

    # every post-filter sample_id must be in the mapping exactly once
    missing = set(h5_sids['sample_id']) - set(m['sample_id'])
    extra = set(m['sample_id']) - set(h5_sids['sample_id'])
    print(f'h5ad sample_ids missing from mapping: {len(missing)}'
          + (f'  e.g. {list(missing)[:10]}' if missing else ''))
    print(f'mapping sample_ids not in filtered h5ad: {len(extra)}'
          + (f'  e.g. {list(extra)[:10]}' if extra else ''))

    dup = m.groupby('sample_id')['patient_id'].nunique()
    bad = dup[dup > 1]
    print(f'sample_ids mapping to >1 patient_id: {len(bad)}'
          + (f'  e.g. {bad.index.tolist()[:10]}' if len(bad) else ''))

    # per-disease patient / fatal-patient counts (patient label = max(y))
    pat = m.groupby(['disease', 'patient_id'])['y'].max().reset_index()
    print('\nper-disease patient counts (patient label = max(y) over draws):')
    summ = pat.groupby('disease').agg(n_patients=('patient_id', 'nunique'),
                                       n_fatal=('y', 'sum'))
    print(summ.to_string())
    print(f'\nTOTAL patients: {pat["patient_id"].nunique()}  '
          f'TOTAL fatal patients: {int(pat["y"].max() if False else pat["y"].sum())}')

    ok = (len(missing) == 0) and (len(bad) == 0)
    print(f'\nVERIFICATION: {"PASS" if ok else "FAIL"}')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
