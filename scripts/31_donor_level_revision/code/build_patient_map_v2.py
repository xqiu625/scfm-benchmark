#!/usr/bin/env python3
"""
Job H2 step 1: rebuild donor_meta_with_patient_id_v2.csv on top of the v1
(Stephenson-corrected) mapping, applying three more corrections found by the
audit2 review (reports/manuscript/revision/literature_audit/T3_patient_sample_mapping/,
patient_overlap_finding.md Task-3 v2):

  Combes (302005data02): GSM4995425, GSM4995426, GSM4995430 are multi-patient
    pooled libraries (2,597 cells) -- DROP entirely, not merge.
  Wilk (302004data01): our Sample_IDs already carry the C-label directly
    ("C1 A", "C1 B", "C2", ...) -- the repeat pair is C1 A/C1 B, no external
    covid_5xx<->C-label mapping needed. Strip the trailing " A"/" B" token.
  Lee (302004data02): merge nCoV_3&4, nCoV_7&8, nCoV_9&10 (paper: C3/C6/C7,
    each sampled twice; pairing confirmed by matched age+sex, not by a
    patient field -- see literature_audit/T3_patient_sample_mapping/mapping_verdicts.csv).

Stephenson's v1 rule (BGCVnn pool prefix -> trailing CVxxxx = real patient) is
unchanged here; it was already verified against the live SDRF and carries
forward as-is.
"""
import re

import pandas as pd

V1 = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
      '31_donor_level_revision/results/donor_meta_with_patient_id.csv')
V2 = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
      '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')
DIFF = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
        '31_donor_level_revision/results_patient/patient_map_v1_to_v2_diff.csv')

COMBES_DROP = {'GSM4995425', 'GSM4995426', 'GSM4995430'}

# Job H2 Step 0 (2026-09-03, Task-4 audit vs Combes 2021 Supp Table 1 p.3):
# three Combes covid_severity labels are wrong -- GSM4995436/442 flip
# fatal->moderate/severe (y 1->0); GSM4995432 severe->moderate (y unchanged).
# Fix-1 (code review 2026-09-04): this was previously applied only at
# read-time in the CV/PCA/cross-disease scripts, never baked into the v2
# mapping CSV itself -- fixed here so anything reading the CSV directly
# (count reporting, null_model_cellcount.py) sees the corrected y.
SEVERITY_FIX_Y = {
    'GSM4995436': 0,  # was fatal (y=1); paper: Discharged=Yes, survived
    'GSM4995442': 0,  # was fatal (y=1); paper: Discharged=Yes, survived
    # GSM4995432: severe->moderate, y unchanged (already 0) -- no y edit needed
}


def main():
    m = pd.read_csv(V1)

    y_fix_rows = []
    for sid, new_y in SEVERITY_FIX_Y.items():
        mask = m['sample_id'] == sid
        old_y = m.loc[mask, 'y'].tolist()
        m.loc[mask, 'y'] = new_y
        y_fix_rows.append({'sample_id': sid, 'cohort': 'Combes',
                            'action': 'Task-4 label fix (covid_severity corrected)',
                            'v1_patient_id': sid, 'v2_patient_id': sid,
                            'n_cells': None, 'y_before': old_y, 'y_after': new_y})
    before_total = m['patient_id'].nunique()
    before_covid = m.loc[m['disease'] == 'COVID-19', 'patient_id'].nunique()
    before_sepsis = m.loc[m['disease'] == 'Sepsis', 'patient_id'].nunique()
    before_rows = len(m)

    diff_rows = list(y_fix_rows)

    # --- Combes: drop pooled GSMs entirely ---
    combes_mask = (m['Data_NO'] == '302005data02') & m['sample_id'].isin(COMBES_DROP)
    dropped = m[combes_mask]
    for _, r in dropped.iterrows():
        diff_rows.append({'sample_id': r['sample_id'], 'cohort': 'Combes',
                           'action': 'DROPPED (pooled library, multi-patient)',
                           'v1_patient_id': r['patient_id'], 'v2_patient_id': None,
                           'n_cells': r['n_cells']})
    m = m[~combes_mask].copy()

    # --- Wilk: "C1 A"/"C1 B" -> "C1" (repeat pair), others unchanged ---
    wilk_mask = m['Data_NO'] == '302004data01'
    def wilk_pid(sid):
        return re.sub(r' [AB]$', '', sid)
    for idx in m[wilk_mask].index:
        old = m.at[idx, 'patient_id']
        new = wilk_pid(m.at[idx, 'sample_id'])
        if new != old:
            diff_rows.append({'sample_id': m.at[idx, 'sample_id'], 'cohort': 'Wilk',
                               'action': 'MERGED (repeat draw, same patient)',
                               'v1_patient_id': old, 'v2_patient_id': new,
                               'n_cells': m.at[idx, 'n_cells']})
        m.at[idx, 'patient_id'] = new

    # --- Lee: merge nCoV_3&4, nCoV_7&8, nCoV_9&10 ---
    lee_mask = m['Data_NO'] == '302004data02'
    merge_pairs = {'nCoV_3': 'nCoV_3', 'nCoV_4': 'nCoV_3',
                   'nCoV_7': 'nCoV_7', 'nCoV_8': 'nCoV_7',
                   'nCoV_9': 'nCoV_9', 'nCoV_10': 'nCoV_9'}
    for idx in m[lee_mask].index:
        sid = m.at[idx, 'sample_id']
        if sid in merge_pairs:
            old = m.at[idx, 'patient_id']
            new = merge_pairs[sid]
            if new != old:
                diff_rows.append({'sample_id': sid, 'cohort': 'Lee',
                                   'action': 'MERGED (repeat draw, same patient, '
                                             'matched age+sex)',
                                   'v1_patient_id': old, 'v2_patient_id': new,
                                   'n_cells': m.at[idx, 'n_cells']})
            m.at[idx, 'patient_id'] = new

    m.to_csv(V2, index=False)

    after_total = m['patient_id'].nunique()
    after_covid = m.loc[m['disease'] == 'COVID-19', 'patient_id'].nunique()
    after_sepsis = m.loc[m['disease'] == 'Sepsis', 'patient_id'].nunique()
    after_rows = len(m)

    after_fatal_total = m.groupby('patient_id')['y'].max().sum()
    after_fatal_covid = (m[m['disease'] == 'COVID-19'].groupby('patient_id')['y']
                          .max().sum())
    after_fatal_sepsis = (m[m['disease'] == 'Sepsis'].groupby('patient_id')['y']
                           .max().sum())

    print(f'v1 -> v2 rows: {before_rows} -> {after_rows} '
          f'(dropped {before_rows - after_rows})')
    print(f'v1 -> v2 TOTAL patients: {before_total} -> {after_total}')
    print(f'v1 -> v2 COVID patients: {before_covid} -> {after_covid} '
          f'(fatal now {int(after_fatal_covid)})')
    print(f'v1 -> v2 Sepsis patients: {before_sepsis} -> {after_sepsis} '
          f'(fatal now {int(after_fatal_sepsis)})')
    print(f'v2 total fatal patients: {int(after_fatal_total)}')

    if diff_rows:
        pd.DataFrame(diff_rows).to_csv(DIFF, index=False)
        print(f'\nDIFF ({len(diff_rows)} rows) -> {DIFF}')
        print(pd.DataFrame(diff_rows).to_string(index=False))


if __name__ == '__main__':
    main()
