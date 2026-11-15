# Deliverable 4 / 4b — the 12 remaining fatal labels and the non-fatal severity grades

All four papers' own supplementary tables were retrieved and read. Three of the four came through Europe PMC's mirror at EBI (`/europepmc/webservices/rest/<PMCID>/supplementaryFiles`); Combes is marked not-open-access there, so its Supplementary Information PDF was fetched from the publisher after network access was approved. Nothing below is inferred from a severity score, and where a paper does not state an outcome the row says so.

Tables: [deliverable4_fatal_label_audit.csv]({{artifact:art_e8c76381-1ebc-41ed-b056-069324e1422f}}) (11 rows) and [deliverable4b_severity_audit.csv]({{artifact:art_b478a373-93e7-4ecf-8bba-d57ca5d4163b}}) (7 rows), with the full per-patient sources as [combes_comet_cohort.csv]({{artifact:art_6dc01aec-e2f5-4f31-a237-71a2113a7c1b}}), [liu_pbmc_subjects.csv]({{artifact:art_b9b03567-3ad1-4c3b-b77f-43ad5d9e5d7c}}), [schulte_cohort1_patients.csv]({{artifact:art_0a3dac0e-629c-4ca3-a0df-8177f16899b4}}), [bernardes_patients.csv]({{artifact:art_e4888be8-85f3-44d3-8a98-a70b8c6ee360}}) and [stephenson_sample_severity.csv]({{artifact:art_fa0a0667-e244-4df5-a800-2e764a98902a}}).

## Task 4 — fatal labels

**Combes 2021: both fatal labels are contradicted by the paper.** Supplementary Table 1 (Discovery Study Cohort) has no death or mortality column at all; its only outcome field is `Discharged`, and for both patients it reads `Yes`:

- COMET 1029 — `1029 POS Mild/Moderate 85 2 2 Viral No 0 Yes` (no ICU, 0 days mechanical ventilation, discharged)
- COMET 1047 — `1047 POS Severe 48.3 2 6 Viral + Bacterial Yes 28 Yes` (ICU, 28 days ventilated, discharged)

One caveat you should weigh before acting: the mapping from the GEO label `Patient N` to `COMET ID 100N` is **not stated** in either the paper or the GEO record. It rests on set equality — the 32 patient numbers in GSE163668 are exactly {COMET ID − 1000} for the 32 discovery-cohort patients in Supplementary Table 1, with no leftovers on either side. That is strong, but it is an inference; if it is wrong, these two rows are wrong with it.

**Liu 2021: all four fatal labels are confirmed.** Table S6B (PBMC Samples) lists 54 samples over 33 subjects — matching the 54 per-subject bulk records in GSE161918 — with an explicit `outcome` column. Exactly four rows read `Deceased`, one per subject:

| subject | sample_label | visit | severity | days from onset | WHO score |
|---|---|---|---|---|---|
| HGR0000071 | BS-EQ-0005 | T0 | Critical | 8 | 5 |
| HGR0000110 | BS-EQ-0013 | T0 | Critical | 3 | 5 |
| HGR0000141 | BS-EQ-0022 | T0 | Critical | 20 | 5 |
| HGR0000430 | BS-EQ-0044 | T1 | Critical | 20 | 5 |

All four are `Critical`, and the sheet has 25 Critical subjects — consistent with the paper's report of 4 deaths among its critical cases. Note HGR0000430's only PBMC row is T1; it has no T0 row in this table.

**Schulte-Schrepping 2020: confirmed, and the count works out exactly.** Table S1 (sheet `all donors`) carries an `Outcome` column. In the first block — scRNA-seq (10x), Cohort 1, the 18 patients / 27 samples that correspond to the EGA single-cell deposit — three patients read `deceased`: C19-CB-0013 (WHO 7, sampled days 8 and 15 after symptom onset), C19-CB-0016 (WHO 7, day 16) and C19-CB-0021 (WHO 6, day 18). That is 3 patients and **4 samples**, matching your count. This block's per-patient sampling-day lists also reproduce the EGA per-subject sample counts exactly, which retroactively settles the open question from Deliverable 3: the `_dN` suffix is the day after symptom onset.

**Bernardes 2020: confirmed, but it is 3 samples, not 2.** Table S1 (`full clinical table`) has a `Recovered` column; exactly two of the 15 patients read `no`, both with `WHO Scale (max) = 8` and `Max disease catergory = Critical (fatal)`:

- Patient 002 — visits TA and TA2, at days 19 and 21 since symptom onset
- Patient 003 — visit TA, at day 12 since symptom onset

So the fatal patient labels are right, but they cover three samples (`002_TA`, `002_TA2`, `003_TA`). If your benchmark carries two, one of patient 002's visits is missing.

## Task 4b — severity grades

**Combes 2021, the flagged mismatch: your `severe` for GSM4995432 disagrees with the paper.** That library is Patient 14 → COMET 1014, and Supplementary Table 1 reads `1014 POS Mild/Moderate 58.6 1 4 None Yes 0 Yes` — Mild/Moderate, ICU during the stay but zero days ventilated, discharged. GEO places the same library in the M/M arm, so both primary sources disagree with `severe`. Also note the paper's scale is binary (Mild/Moderate vs Severe): there is no separate `moderate` grade to map onto. For GSM4995436 (Patient 29 → COMET 1029) the M/M placement **is** consistent with the paper; it is only the fatal label that is not.

**Liu 2021.** Table S6B gives NIH severity per subject — 25 Critical, 5 Severe, 3 Moderate — and a `WHO_severity_score` that varies by visit (5 on 40 samples, 4 on 10, 3 on 2, and one each of 6 and 1). Severity is subject-level and constant across a subject's visits; the WHO score is the visit-level field.

**Schulte-Schrepping 2020.** In the cohort-1 scRNA-seq block, the `Condition` label and the WHO ordinal separate cleanly: every mild patient carries WHO 3 or 4, every severe patient WHO 5, 6 or 7, with no overlap — so the paper's own mild/severe split is exactly WHO ≤4 vs ≥5 here. For repeat-draw patients the WHO field is a per-visit list (C19-CB-0051 reads `4, 5, 5, 4, 4` across five draws), so a flat per-patient severity label discards information.

**Bernardes 2020 — a correction to the brief's premise.** Severity in Table S1 is **patient-level, not visit-level**: the columns are `WHO Scale (max)`, `Max disease catergory` (values such as `mild (recovering)/ asymptomatic`, `complicated (recovering)`, `complicated (hyperinflammatory)`, `Critical`, `Critical (fatal)`) and CRP / IL-6 / ferritin *ranges* across the stay. What is genuinely per-visit in that table is the sampling day — and it is anchored to symptom onset, not admission, which resolves the day-numbering conflict flagged in Deliverable 3:

![TA is recorded anywhere from day 4 to day 46 since symptom onset, so the visit codes' day numbers are nominal admission days]({{artifact:art_56542b98-c783-437b-b6b7-98ec4dbe071c}})

The figure legend's day 1 / 3 / 8 / 11 / 14 / 20 are nominal days after admission; Table S1's per-visit days since onset for the same code TA run from 4 to 46 (median 12). The two are not in conflict — they are different anchors, and only the Table S1 value is patient-specific.

**Stephenson 2021: not in the SDRF.** The SDRF's only clinical field is `Characteristics[disease]`; there is no WHO or severity column (full column list is in the audit table). The per-sample grades live in the published object's `obs`: `Status_on_day_collection` (nine values: Ward_noO2, Ward_O2, Ward_NIV, ITU_O2, ITU_NIV, ITU_intubated, Healthy, Non_covid, Staff screening — no LPS category; LPS belongs to the separate `Status` field, whose values are Covid 102, Healthy 24, LPS 12, Non_covid 5 samples), `Status_on_day_collection_summary` (Moderate 32, Mild 26, Healthy 24, Critical 17, Severe 15, Asymptomatic 12, LPS 12, Non_covid 5 samples) and `Worst_Clinical_Status`, which includes `Death` on 8 samples. The first is per draw and the last is per patient, so they differ by design for patients who deteriorated after their draw.

**Ren 2021: the scheme does not flatten cleanly.** GSE158055 gives three per-sample fields — `covid-19 severity` (mild/moderate, severe/critical, control), `sample time` (progression, convalescence, control) and `outcome` (discharged, deceased, control). Observed combinations: mild/moderate + convalescence + discharged 89; severe/critical + progression + discharged 68; severe/critical + convalescence + discharged 50; mild/moderate + progression + discharged 33; control 28; severe/critical + progression + deceased 15; severe/critical + convalescence + deceased 1. Mapping onto a flat mild / moderate / severe / recovered / fatal set therefore has three problems: mild and moderate cannot be separated (the paper's scale bins them together); `convalescence` is a disease phase, not an outcome, and one convalescent sample belongs to a patient who died, so mapping convalescence → `recovered` would mislabel it; and the deposit records 12 deceased patients across 16 samples, which is more than the 9 Ren fatal samples in the benchmark — worth checking which subset you kept.

## Open items

- The Combes `Patient N` ↔ `COMET 100N` correspondence is an inference from an exact set match, not a stated mapping. Only the authors (or a COMET data dictionary) can close it.
- Combes reports no deaths at all for the discovery cohort; two validation-cohort patients carry `Discharged = No`, but that is not a stated death and those patients are not in GSE163668.
- Bernardes' visit codes `TD` and `rec` remain undefined in the paper text; Table S1 does give their sampling days (patient 012 TD at day 50 since onset; patients 006 and 007 `rec` at days 44 and 45).
