# Deliverable 3 — patient ↔ sample mapping rules, verified against repository metadata

Scope: the nine cohorts named in the instruction file. For each one: the inferred rule as given, the verdict against authoritative repository metadata, what each draw means clinically, per-patient outcome where the deposit supports it, and the verbatim evidence. Every quote in `mapping_evidence.csv` was checked programmatically against the retrieved source text (15/15 exact after whitespace normalisation); nothing here is retyped from memory.

Tables: [mapping_verdicts.csv]({{artifact:art_d4b87b9b-6ab2-4332-8b69-69d919b8112c}}) (one row per cohort, 15 fields) and [patient_sample_map.csv]({{artifact:art_2ecd7f61-c55a-45f5-94f3-b4336aa3e0a1}}) (806 rows — every sample in the nine deposits with its patient key, draw label, draw time and outcome as recorded).

## Verdict summary

| Cohort | Accession | Verdict | Patients | Samples | Patients with >1 draw |
|---|---|---|---|---|---|
| Qiu 2021 (sepsis) | GSE167363 | CONFIRMED (structure); label strings corrected | 7 | 12 | 5 |
| Liu 2021 (COVID, N. Italy) | GSE161918 | CONFIRMED for bulk records; does **not** apply to the single-cell libraries | 33 | 54 | 17 |
| Schulte-Schrepping 2020 | EGAS00001004571 | CONFIRMED for grouping; day anchor inferred | 52 (+16 pooled entries) | 141 | 36 |
| Stephenson 2021 | E-MTAB-10026 | **CORRECTED — the rule is inverted** | 130 | 143 | 12 |
| Bernardes 2020 | GSE161777 | CONFIRMED as visit codes, **not** replicates | 14 | 101 | 11 |
| Combes 2021 | GSE163668 | **CORRECTED — 6 of 38 libraries are multi-patient pools** | 46 | 38 | 0 |
| Ren 2021 | GSE158055 | CONFIRMED for grouping; suffix semantics corrected | 196 | 284 | 31 |
| Wilk 2020 | GSE150728 | CONFIRMED for the repeat pair; label mapping and one sample unresolved | 12 | 13 | 1 |
| Lee 2020 | GSE149689 | SUPPORTED by metadata, not confirmable by label | 17 | 20 | 3 |

## The two rules that are wrong

**Stephenson 2021 — `BGCV02_CV0052` is pool_patient, not patient_sample.** The ArrayExpress sample table gives `Characteristics[individual]` = `CV0052` and `Characteristics[sequence pool]` = `BGCV02`; the published object's `obs` carries `patient_id` = `CV0052` and `sample_id` = `BGCV02_CV0052`. The first token is the 10x sequencing pool. There are 15 pools holding exactly 4 donors each, and 13 of the 15 contain more than one distinct `Outcome`. That is the whole explanation for the mixed labels you saw: pool `BGCV02` contains CV0052 (Death) together with CV0059, CV0068 and CV0902 (all Home).

![Pool BGCV02 and 12 other pools mix outcomes because each pool is four different patients]({{artifact:art_a2b48168-3a48-4b76-8de2-915bb385f1f6}})

Read the right way round, outcome is perfectly patient-consistent: across all 130 patients, zero have more than one distinct `Outcome` value (Home 108 samples, unknown 27, Death 8; 8 patients died — CV0052, CV0180, CV0198, MH8919329, MH9143320, MH9143325, MH9179821, newcastle004v2). Per-patient values are in [stephenson_patient_outcomes.csv]({{artifact:art_4f9e1311-2b5f-4362-b643-bbf63679ef4c}}), pool membership in [stephenson_pool_composition.csv]({{artifact:art_e16f4038-86ba-40cc-8f5d-b61e16332d81}}).

Two further cautions for this cohort. Repeat draws are marked by `Resample` (Initial/Resample) and `Collection_Day` (D0/D28): 12 patients have a D0 + D28 pair. The thirteenth multi-sample patient, CV0902, is a healthy donor whose two libraries sit in pools BGCV01 and BGCV02 on the same day — a duplicate library, not a second draw. And `Days_from_onset` is recorded once per patient at the initial draw; it carries the same value on the D28 row, so it is not the time of that draw.

**Combes 2021 — a bare GSM is not always one patient.** Six of the 38 libraries are pools: Patients 1+2+3, 5+6, 7+8+9, 10, 17+20+21, 50+51, each with `pooled = TRUE` and `covid_status`/`phenotype` set to `mixed (pooled sample)`. Donor identity inside them requires demultiplexing. The remaining 32 are single-subject (18 patients, 14 control donors) and no subject is sampled twice. One trap worth encoding in the loader: patient and donor numbers are separate namespaces, so `Patient 29` and `Donor 29` are different people.

## The rules that hold, with corrections

**Qiu 2021.** The two-draws-per-patient structure is right and the anchor is diagnosis, not admission: sequencing was performed on PBMCs from 5 sepsis patients "at 0 and 6 h post diagnosis", with three survivors discharged from ICU and two fatal courses. The label strings differ from the inferred ones — GEO titles are `NS ES_T0/T6`, `NS LS_T0/T6`, `S1 T0/T6`, `S2_T0/T6`, `S3_T0/T6`, `HC1`, `HC2`. There is no `NS1`/`NS2`, and the first two pairs separate the timepoint with a space rather than an underscore, which will break a strict `_T[06]$` parse.

**Liu 2021.** `HGR…_Tn` is correct — but only for the bulk RNA-seq records. Of 241 samples in GSE161918, 71 are per-subject bulk libraries titled `Batch N bulk RNAseq HGR0000071` with a `timepoint` characteristic of Timepoint 0–3 (54 COVID draws over 33 subjects: 16 with one draw, 14 with two, 2 with three, 1 with four), and 170 are 10x lanes titled `Batch N PBMC 10x Lane N …` whose `timepoint` value is literally `mixed`. The single-cell libraries are multi-subject pools carrying no HGR id, so if your pipeline reads the single-cell data, title-based patient assignment is not available at all. Clinically, T0 is "at or near the time of hospital admission" and the draws "covered a range of times post symptom onset" rather than fixed intervals.

**Schulte-Schrepping 2020.** Grouping by subject is confirmed, and by an exact-count argument rather than by reading the suffixes: the EGA study's `C19-CB-*` subject set contains exactly 27 samples across 18 subjects, matching the paper's cohort-1 statement of "27 samples from 18 COVID-19 patients … collected between day 3 and day 20 after symptom onset". EGA sample titles are null, so the `_dN` strings are not in the deposit; the day anchor is symptom onset by inference from the paper's reporting for both cohorts (cohort 2 "sampled between 2 and 25 days after symptom onset"), never stated as a naming convention. Sixteen additional entries have `subject_id` = `pooled_*` (e.g. `pooled_200416`) and cannot be attributed to a single subject.

**Bernardes 2020.** Visit codes, not replicates — the figure legend spells the scheme out: patient ID (001–014) with "day 1 (after admission TA), day 3 (TA2), day 8 (TB), day 11 (TB2), day 14 (TC), and day 20 (TE)". Two qualifications. The methods give a different scheme for the same visits — "day 0 (at admission), day 2, day 7, day 10, and day 13" — so treat the day numbers as nominal visit days with a ±1 day ambiguity, not measured intervals. And the whole-blood deposit GSE161777 (bulk RNA + BCR, with a `patient number` characteristic that confirms the grouping) contains two codes the paper never defines: `TD` (2 samples) and `rec` (patients 6 and 7), the latter presumably the discharge sample implied by "samples were taken at days 0, 2, 7, 10, 13 and/or at discharge". The single-cell data are on FastGenomics and were not retrieved, so per-cell labels were not checked against these records.

**Ren 2021.** Group by the `patient` characteristic (`P-M003`), which maps 1:1 onto the title base (`S-M003`) for all 284 samples — but the numeric suffix is a sample index, not a draw index. 41 patient-days carry more than one library with identical sampling day, sample type, platform and batch (P-M003 has `S-M003-2` and `S-M003-3` both at day 25 from fresh sputum), while other suffixes do mark distinct draws (`S-M003-1` at day 13). The authoritative draw key is patient + `sampling day (days after symptom onset)`: 165 patients have one draw, 26 have two, 5 have three. Outcome is per-sample in the record and unique within every one of the 196 patients (discharged 240, control 28, deceased 16). Sample types are heterogeneous — PBMC, sputum, sorted B cells — and should be filtered before any cross-sample comparison.

**Wilk 2020.** The repeat pair is real and clinically meaningful: "One patient (C1) was sampled twice: at nine days post-symptom onset while only requiring supplemental oxygen and at 11 days post-symptom onset following intubation" — so `covid_555_1` vs `_2` straddles intubation, and the cohort window is 2–16 days after onset. Two things do not reconcile. GEO holds 7 COVID libraries (`covid_555_1`, `covid_555_2`, `covid_556`–`559`, `covid_561`; note the gap at 560) = 6 patients, against the paper's "eight peripheral blood samples from seven hospitalized patients". And the paper labels patients C1–C7 while the deposit uses `covid_5xx`, with no mapping given in either — so per-patient clinical variables cannot be joined to the GSMs.

**Lee 2020.** The paper states that "PBMCs were obtained twice from three (the subject C3, C6, and C7) of the eight COVID-19 patients at different time points during hospitalization", and the repository is consistent with it in the only way available: exactly three adjacent duplicate age/sex pairs occur among the 11 nCoV samples — nCoV 3 & 4 (67 F), nCoV 7 & 8 (73 F), nCoV 9 & 10 (61 M). There is no patient field, so the pairing rests on that coincidence of demographics and record order, and which pair is C3, C6 or C7 cannot be determined. Note that `nCoV 9` is labelled severe and `nCoV 10` mild for what is almost certainly the same person: severity was assigned "according to the National Early Warning Score (NEWS; mild < 5, severe ≥ 5) evaluated on the day of whole blood sampling", so a within-patient label change is expected behaviour and is not evidence of mis-grouping.

## What remains open

- **Wilk 2020**: the eighth sample / seventh patient is absent from GSE150728, and the `covid_5xx` ↔ `C1–C7` mapping is unstated. Closing either needs the authors' clinical table keyed to GEO ids.
- **Lee 2020**: `nCoV n` ↔ `Cn` mapping unstated; no per-sample draw day.
- **Schulte-Schrepping 2020**: `_dN` suffix strings are not in EGA (titles null); the day anchor is inferred, and the per-sample day would need the paper's Table S1 keyed to EGA sample ids.
- **Liu 2021** and **Combes 2021**: donor identity for pooled single-cell libraries (168 lanes and 6 libraries respectively) requires demultiplexing metadata that is not in either deposit.
- **Bernardes 2020**: `TD` and `rec` are undefined in the paper; the ±1 day discrepancy between methods and figure legend is unresolved in the published text.
- **Combes 2021**: no draw timing in the repository record.
