# Literature Research Brief — Platelet scFM Benchmark Revision

You are doing literature research for a manuscript revision. The manuscript benchmarks
five single-cell foundation models (scFMs) on **human platelet** transcriptomes
(42,095 cells, 322 donors, COVID-19 + sepsis, in-hospital mortality task), claiming
platelets are "essentially absent" from scFM pretraining corpora (an out-of-distribution,
OOD, stress test). A senior reviewer pushed back on two things, which map to the two
tasks below.

## Hard rules

- Primary sources only: the papers themselves, GEO/ArrayExpress/EGA records, official
  model cards / GitHub / HuggingFace pages, CELLxGENE census. No secondary summaries.
- For every claim, give an **exact verbatim quote** + source URL/DOI/PMID.
- If a fact is not stated anywhere, write **"not reported"** — never infer or guess.
- If sources disagree, report both with quotes.
- Mark anything you could not verify as "UNVERIFIED" rather than smoothing it over.

---

## Task 1 — Quantify how OOD platelets really are for each scFM (revision item P1-E)

The reviewer says: "These models were not trained on platelets, not intended to use on
them" — but "never seen" and "0.3% of the corpus" are different claims. We need the
actual numbers.

For EACH of the five models below, find: (a) pretraining corpus name and total size
(cells); (b) whether **platelets** (anucleate, circulating) are in the corpus — count
and/or fraction, if reported; (c) whether **megakaryocytes** (the nucleated bone-marrow
progenitors of platelets) are in the corpus — count and/or fraction; (d) whether the
corpus is whole-atlas (many tissues) or restricted (e.g. perturbation screens, cell
lines); (e) exact quotes + links.

| Model | Checkpoint used in the benchmark | Paper |
|---|---|---|
| STATE | STATE (2,058-dim embeddings) | Adduri et al. 2025, bioRxiv doi:10.1101/2025.06.26.661135 (Arc Institute) |
| UCE | UCE (1,280-dim) | Rosen et al. 2026, Nature, doi:10.1038/s41586-026-10689-z |
| Geneformer | Geneformer-V2-316M (18 layers, 1,152-dim) | Theodoris et al. 2023, Nature 618:616–624 + any V2 update notes |
| scPRINT-2 | `small-v2.ckpt` (official HuggingFace scprint2 checkpoints) | Kalfon et al. 2025, Nat Commun 16:3607 |
| scGPT | whole-human pretrained model (512-dim, `best_model.pt`) | Cui et al. 2024, Nat Methods 21:1470–1480 |

Also answer, as a cross-check independent of any single paper:
1. In the current **CZ CELLxGENE census**, what fraction of cells are annotated
   "platelet" and what fraction "megakaryocyte" (cell_type ontology terms)? This bounds
   what any CELLxGENE-trained model could have seen. (cellxgene-census API or the
   CELLxGENE portal cell-type facet counts are both acceptable; cite what you used.)
2. For corpora that are NOT whole-atlas (e.g. if STATE is trained mainly on
   perturbation screens such as Tahoe-100M / Replogle): state explicitly what cell
   types/contexts dominate, with quotes.
3. Note explicitly the biological distinction: circulating platelets are anucleate
   cell fragments; megakaryocytes are nucleated progenitors. A corpus containing
   megakaryocytes is NOT the same as containing platelets — flag any source that
   conflates them.

**Deliverable 1:** a table (one row per model): corpus, size, platelet count/fraction,
megakaryocyte count/fraction, corpus type, evidence quotes + URLs. Plus a 3–5 sentence
verdict: for which models is "essentially absent" accurate, and for which must the claim
be softened (e.g. "absent except possibly rare megakaryocytes, ≤X%").

---

## Task 2 — Sampling windows and design of the 11 source cohorts (revision item P1-D)

The reviewer asks: "When were the samples taken during the disease? In sepsis this is
critical!" We are adding a cohort-characteristics table. For each study below, extract:

1. **Sampling window**: when blood was drawn relative to symptom onset / hospital
   admission / ICU admission / diagnosis — verbatim definition from the paper or GEO
   record.
2. **Design**: cross-sectional single draw vs longitudinal; if longitudinal, list the
   timepoints and which are in the public data.
3. **Platform**: scRNA-seq chemistry (e.g. 10x 3' v2/v3, 5').
4. **Severity/outcome definitions** used (e.g. WHO ordinal scale, qSOFA, ICU status).
5. Anything reported about **fatal cases specifically** (when they died relative to
   sampling), if stated.

| # | Study | Accession | Donors in benchmark | Notes |
|---|---|---|---|---|
| 1 | Qiu et al. 2021, J Leukoc Biol, PMID 34558746 | GSE167363 | 10 (sepsis, 4 fatal) | **ALREADY CONFIRMED — do not re-research:** ICU enrollment at sepsis recognition (qSOFA≥2), draws at T0 and T+6h. Use as the format example. |
| 2 | Jiang et al. 2020, JCI Insight, PMID 32554932 | GSE151263 | 7 (sepsis-arm, 0 fatal) | **PRIORITY** (the second sepsis-arm study). Framed as pneumonia/early-ARDS in a sepsis population — get its exact enrollment and draw timing. |
| 3 | Ren et al. 2021, Cell 184 | GSE158055 | 59 (COVID) | largest COVID atlas |
| 4 | Stephenson et al. 2021, Nat Med 27 | E-MTAB-10026 | 102 (COVID) | UK cohort, multi-omics |
| 5 | Liu et al. 2021, Cell 184(7):1836–1857 | (find GEO/SRA accession) | 52 (COVID) | "Time-resolved systems immunology reveals a late juncture linked to fatal COVID-19" — longitudinal; list ALL timepoints and which are public |
| 6 | Schulte-Schrepping et al. 2020, Cell 182 | EGAS00001004571 | 23 (COVID) | controlled access; use paper + ArrayExpress/EGA metadata |
| 7 | Bernardes et al. 2020, Immunity 53 | (FastGenomics; find accession) | 28 (COVID) | longitudinal megakaryocyte study |
| 8 | Combes et al. 2021, Nature 591 | GSE163668 | 16 (COVID) | |
| 9 | Wilk et al. 2020, Nat Med 26 | GSE150728 | 7 (COVID, 0 fatal) | |
| 10 | Lee et al. 2020, Sci Immunol 5 | GSE149689 | 11 (COVID, 0 fatal) | |
| 11 | Arunachalam et al. 2020, Science 369 | GSE155673 | 7 (COVID, 0 fatal) | |

**Deliverable 2:** a table with one row per study: accession, disease, donors,
sampling window (verbatim), design/timepoints, platform, severity definition,
evidence quotes + URLs. Studies that do not report timing get "not reported" in the
sampling-window column — that is an acceptable and expected outcome for several 2020
COVID studies; do not fill gaps with plausible guesses.

Context for prioritization: sepsis timing (studies 1–2) is the reviewer's critical
issue; the COVID studies (3–11) matter mainly to state honestly in Limitations that
sampling windows vary across cohorts and cannot be aligned across diseases.

---

## Task 3 — Verify patient↔sample mapping in the multi-draw cohorts (URGENT, blocks re-analysis)

We discovered that 6 cohorts encode multiple draws of one patient as separate sample
accessions, which our pipeline treated as independent donors. For each cohort below,
verify from the source metadata (GEO/ArrayExpress sample records, patient-supplement
tables) whether the grouping rule we inferred is CORRECT, and report the authoritative
patient↔sample mapping (patient identifier → list of sample accessions + draw time):

| Cohort | Our inferred rule (to verify) | Specific question |
|---|---|---|
| Qiu 2021 (GSE167363) | `NS1_T0`,`NS1_T6` = one patient at 0 h / 6 h | Confirm 5 patients × 2 draws; confirm NS1/NS2 are the fatal ones |
| Liu 2021 (Cell) | `HGR0000042_T0`,`_T1` = one patient, timepoints T0/T1/T2/T3 | Confirm HGR id = patient, suffix = timepoint; what do T0–T3 correspond to clinically (days since admission?) |
| Schulte-Schrepping 2020 (EGA) | `C19-CB-0001_d7`,`_d16` = one patient, day 7 / day 16 | Confirm `_dN` = day of sampling relative to what anchor |
| Stephenson 2021 (E-MTAB-10026) | `BGCV02_CV0052` = patient BGCV02, sample CV0052 | Confirm BGCVnn = patient; **and whether outcome (fatal/survived) is a patient-level attribute** — we see patients whose samples carry different labels (e.g. one draw "fatal", others not); report the authoritative per-patient outcome |
| Bernardes 2020 (FastGenomics) | `001_TA`,`001_TA2`,`001_TB`,`001_TC` = one patient, draws TA/TA2/TB/TC | Confirm; what do TA/TA2/TB/TC mean (visit codes? replicates?) |
| Combes 2021 (GSE163668) | unknown — 16 samples appear as bare GSM accessions | Map each GSM to a patient via GEO sample metadata/relations; are any GSMs draws of the same patient? |
| Ren 2021 (GSE158055) | `S-M026-1`,`S-M026-2` grouped as one patient | Are `-1/-2` technical replicates or distinct draws? Either way same patient? |
| Wilk 2020 (GSE150728) — NEW | 7 IDs (`C2`–`C7` + 2) treated as 7 distinct patients | The paper reports one donor sampled at 9 d AND 11 d post onset — do two of the 7 sample records map to one donor? |
| Lee 2020 (GSE149689) — NEW | 11 IDs (`nCoV_1`–`nCoV_11`) treated as 11 distinct patients | The paper reports 3 donors sampled twice — which samples are repeat draws, and are they among the 11 in the benchmark? |

**Deliverable 3:** per cohort: confirmed/corrected patient↔sample mapping rule,
the clinical meaning of each draw/timepoint, authoritative per-patient outcome where
available, quotes + URLs. If a cohort's metadata does not allow the mapping, say so
explicitly.

---

## Task 4 — Verify the 12 remaining fatal labels against paper supplements

Context: our benchmark's mortality label is `fatal` vs `survived` per patient.
A repository-metadata audit already verified 21 of our 33 fatal sample IDs
(Ren 9, Stephenson 8, Qiu 4 — 100% match). The remaining 12 cannot be checked
from repository records and need the papers' own supplementary tables. For each
case: does the paper/supplement record this patient as **died (in hospital)**?

| Cohort | Our fatal sample IDs | Where to look |
|---|---|---|
| Combes 2021, Nature 591 (GSE163668) | Patient 29 (GSM4995436), Patient 47 (GSM4995442) | Paper supplement patient table. NOTE: Patient 47 is marked "POSITIVE/SEVERE" in GEO but that is covid-status/severity, NOT outcome; also "Patient 29" and "Donor 29" are different namespaces — do not confuse |
| Liu 2021, Cell 184 (GSE161918) | 4 fatal samples (HGR-prefixed IDs) | Paper/supplement patient-level outcome table; the paper reports 4 deaths among 25 critical cases but does not time them to draws — we need the per-subject outcome, and ideally which subjects the single-cell libraries came from |
| Schulte-Schrepping 2020, Cell 182 (EGAS00001004571, controlled) | 4 fatal samples (C19-CB-xxxx_dN IDs) | Paper Table S1 (per-patient outcome); EGA itself is controlled-access |
| Bernardes 2020, Immunity 53 | 2 fatal samples (IDs like 00x_TA/TB) | Paper supplement per-patient outcomes; bulk companion dataset GSE161777 carries a patient-number characteristic |

**Deliverable 4:** one row per patient: our sample ID, paper's patient identifier,
outcome as stated in the paper (died / survived / not reported), exact quote +
page/table, URL. If a patient's outcome is not stated anywhere in the paper, write
"not reported" — that is an acceptable answer; do not infer from severity scores.
Same hard rules as before: primary sources only, verbatim quotes, no guessing.

### Task 4b — non-fatal severity grades (same rules: verbatim quotes or "not reported")

The repository cross-check already cleared Lee (11/11) and Combes (11/13; two
flagged below). For the remaining cohorts, verify our per-patient severity grade
(mild / moderate / severe / recovered) against the paper's own per-patient table
or supplement — using each paper's OWN severity scale (they differ: WHO ordinal,
NIH, NEWS, Berlin ARDS, ICU/ventilation, modified WHO+biomarkers):

| Cohort | What to verify |
|---|---|
| Combes 2021 (GSE163668) | **One mismatch to adjudicate:** our GSM4995432 is labeled "severe" but the GEO record groups it "POSITIVE/M/M" (mild/moderate). Also GSM4995436 (our "fatal") sits in the M/M arm — is that consistent with the paper's per-patient table? |
| Liu 2021 (GSE161918) | per-subject NIH grade (mild/moderate/severe/critical) for the subjects in our subset |
| Schulte-Schrepping 2020 (EGA) | per-patient WHO 2–4 vs 5–7 assignment (paper Table S1) |
| Bernardes 2020 (Immunity 53) | per-patient modified-WHO grade (with CRP/IL-6/ferritin) across visits — note severity is visit-level here |
| Stephenson 2021 (E-MTAB-10026) | per-sample WHO-style ordinal grade from the ArrayExpress SDRF (if present) |
| Ren 2021 (GSE158055) | the moderate/severe × progression/convalescence scheme — report how it maps onto a flat mild/moderate/severe/recovered/fatal label set |

**Deliverable 4b:** one row per cohort: our label, the paper's grade for the same
patient/visit, agree / disagree / not reported, exact quote + URL. Disagreements
get quoted, not explained away.

---

## Task 6 — Direction-framing evidence (positioning, not pipeline)

Context: our expanded sepsis arm (91 patients / 17 fatal, 28-day mortality,
draws at median 16 days post admission) shows a NULL within-sepsis result —
all five scFMs at 0.51–0.55 patient-level AUC, a PCA-200 baseline at 0.43 —
while within-COVID mortality is predictable (AUC ~0.78). Before we finalize the
paper's framing, answer these four questions from the literature (same hard
rules: primary sources, verbatim quotes, "not reported" where absent):

**6a. Is within-disease mortality/outcome prediction from platelet (or platelet-
containing PBMC) transcriptomics established for SEPSIS at all?** We know Qiu
2021 (fatal-sepsis platelet signatures at 0–6 h post diagnosis, tiny n) and the
COVID literature. What we need: any study claiming platelet transcriptome
predicts sepsis outcome/mortality, WITH the sampling window stated. The goal is
to learn whether outcome-predictive platelet signal in sepsis is documented
only in EARLY windows (hours–days) or also at >1–2 weeks post admission — i.e.,
whether our null at median-16-day draws is consistent with the field or an outlier.

**6b. Which clinical endpoints ARE established as predictable from platelet
transcriptomics (any disease)?** Build the list (COVID severity, sepsis mortality,
SLE, dengue, cancer therapy response, etc.) with cohort sizes and endpoints. We
need it to position which tasks are "signal-positive anchors" for a benchmarking
paper — and to check we are not claiming a first/negative result someone has
already published.

**6c. Does the Chongqing source paper itself (Nat Immunol 2025,
doi:10.1038/s41590-025-02345-x) report platelet or megakaryocyte findings in its
sepsis cohort?** If the original study found platelet-lineage signals, our null
is specifically about 28-day MORTALITY prediction, not about platelet biology in
this cohort generally — quote what they found (or state "no platelet analysis").

**6d. What do recent scFM evaluation/benchmark papers report for CLINICAL or
donor-level outcome prediction specifically?** (Not cell-type annotation or
batch integration.) We need the 3–6 closest evaluations: do foundation-model
embeddings beat simple baselines on patient-level clinical endpoints anywhere,
and under what evaluation protocols (patient-level splits or not)? This decides
whether our Discussion should claim "consistent with the emerging negative
literature" or something stronger/novel.

**Deliverable 6:** four short sections (one per question), each ending in a
one-sentence verdict we can cite in the Discussion, with verbatim quotes + URLs.
Explicitly flag any claim you could not verify.
