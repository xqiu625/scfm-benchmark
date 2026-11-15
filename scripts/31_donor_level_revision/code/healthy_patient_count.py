#!/usr/bin/env python3
"""
Job H step 6: apply the same patient_id collision rules (2026-09-02 finding,
reports/manuscript/revision/patient_overlap_finding.md) to the healthy-control
arm (covid_severity == 'healthy', excluded from the mortality task) and report
the healthy patient count.

Per-Data_NO regex rules (verbatim from the finding doc's evidence table):
  302005data03 (Qiu sepsis)         strip trailing _T\d+
  302004data16 (Liu)                strip trailing _T\d+
  302004data04 (Schulte-Schrepping) strip trailing _d\d+
  302004data17 (Stephenson)         CORRECTED 2026-09-03: BGCVnn is the sequencing
                                     POOL (E-MTAB-10026 SDRF Characteristics[sequence
                                     pool]), not the patient. Real patient = SDRF
                                     Characteristics[individual] = the CVxxxx suffix
                                     after the pool prefix (e.g. BGCV02_CV0052 ->
                                     CV0052). Non-BGCV Stephenson IDs (MH89191xx,
                                     newcastleNN) are already patient-level, unchanged.
  302005data01 (Bernardes)          strip trailing _T[A-Z]\d*
  302004data15 (Ren)                strip trailing -\d+ (conservative merge)
  all other Data_NO                 unchanged (Sample_ID == patient_id)
"""
import re
import h5py
import numpy as np
import pandas as pd

H5AD = '/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c/state_covid_canonical.h5ad'


def read_categorical(group):
    if isinstance(group, h5py.Group) and 'categories' in group:
        cats = group['categories'][:].astype(str)
        codes = group['codes'][:]
        return cats[codes]
    return group[:].astype(str)


def derive_patient_id(sample_id, data_no):
    if data_no in ('302005data03', '302004data16'):
        return re.sub(r'_T\d+$', '', sample_id)
    if data_no == '302004data04':
        return re.sub(r'_d\d+$', '', sample_id)
    if data_no == '302004data17':
        # BGCVnn_CVxxxx -> CVxxxx (real individual); non-BGCV IDs already patient-level
        if sample_id.startswith('BGCV') and '_' in sample_id:
            return sample_id.split('_', 1)[1]
        return sample_id
    if data_no == '302005data01':
        return re.sub(r'_T[A-Z]\d*$', '', sample_id)
    if data_no == '302004data15':
        return re.sub(r'-\d+$', '', sample_id)
    return sample_id


def main():
    with h5py.File(H5AD, 'r') as h:
        severity = read_categorical(h['obs/covid_severity'])
        sample_id = read_categorical(h['obs/Sample_ID'])
        disease = read_categorical(h['obs/Disease'])
        data_no = read_categorical(h['obs/Data_NO']) if 'Data_NO' in h['obs'] else None

    keep = severity == 'healthy'
    n_cells = int(keep.sum())
    sid = sample_id[keep]
    dno = data_no[keep] if data_no is not None else np.array(['unknown'] * n_cells)
    dis = disease[keep]

    df = pd.DataFrame({'sample_id': sid, 'Data_NO': dno, 'disease': dis})
    df['patient_id'] = [derive_patient_id(s, d) for s, d in zip(df['sample_id'], df['Data_NO'])]

    n_donors = df['sample_id'].nunique()
    n_patients = df['patient_id'].nunique()
    print(f'healthy arm: {n_cells} cells, {n_donors} Sample_IDs (donors) -> '
          f'{n_patients} patients')
    print('\nper-Data_NO breakdown:')
    g = df.groupby('Data_NO').agg(n_cells=('sample_id', 'size'),
                                   n_sample_ids=('sample_id', 'nunique'),
                                   n_patients=('patient_id', 'nunique'))
    print(g.to_string())

    out_path = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
                '31_donor_level_revision/results_patient/healthy_patient_mapping.csv')
    import os
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    df.drop_duplicates(subset='sample_id').to_csv(out_path, index=False)
    print(f'\nSAVED -> {out_path}')


if __name__ == '__main__':
    main()
