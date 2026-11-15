# Revision memo — scFM pretraining-corpus platelet content (P1-E) and source-cohort sampling windows

Scope: 5 single-cell foundation models, 11 source cohorts. Primary sources only (papers, GEO/ArrayExpress/EGA records, official model cards / GitHub / HuggingFace, CELLxGENE Census). Every quantitative claim below traces to a verbatim quote in the companion evidence tables; fields no source states are written **not reported**, and anything not closable from an accessible primary source is marked **UNVERIFIED**.

---

## 1. Task 1 — how out-of-distribution are platelets, really?

### 1.1 Verdict

The reviewer is right that "never trained on platelets" and "0.3 % of the corpus" are different claims — and the honest answer is that **no model reports either number.** For all five models, platelet and megakaryocyte content of the pretraining corpus is **not reported**: none of the five papers, and none of the official model or dataset cards, publishes per-cell-type counts or fractions for its pretraining corpus. Searching the full texts and cards for `platelet`, `thrombocyte`, `megakaryocyte` and `CL:0000233`/`CL:0000556` returns zero corpus-composition hits (the one `Megakaryocytes` occurrence in the scGPT text is a cell-type label in a downstream PBMC-10K annotation benchmark, not a corpus description).

So the claim **"essentially absent" is not supportable as written for any of the five models, and is not refutable either — it is unmeasured at the corpus level.** What *is* measurable is the upstream data source: four of the five models draw their corpus from CELLxGENE Census, and Census composition can be counted exactly. That converts the claim from an assertion into a bound.

**Recommended reframing:** replace "essentially absent" with a quantified, sourced statement about the underlying atlas — platelets are present but are a ~0.1 % minority of the Census releases these models were built from, and the models themselves never report their own cell-type composition. That concedes nothing the data does not support and is stronger against this reviewer than the original wording.

### 1.2 Per-model corpus facts

| Model | Checkpoint in benchmark | Corpus | Size (cells) | Corpus type | Platelets reported? | Megakaryocytes reported? |
|---|---|---|---|---|---|---|
| STATE | STATE (2,058-dim; SE component) | not reported (no corpus name in the bioRxiv record, GitHub README, or the SE-600M / SE-100M / Stack-Large cards) | SE: 167 M observational cells; ST: >100 M perturbed cells | MIXED / partly restricted — **UNVERIFIED** (SE observational; ST on perturbation screens — the repo names Tahoe-100M, i.e. 50 cancer cell lines × 1,100 small molecules, and Replogle) | not reported | not reported |
| UCE | UCE (1,280-dim) | Integrated Mega-scale Atlas (IMA) | >36 M across >300 datasets; 33.9 M / 285 datasets from CELLxGENE Census v2023-07-10; 2.3 M / 28 datasets across 8 species | whole-atlas, multi-species | not reported | not reported |
| Geneformer | Geneformer-V2-316M (18 layers, 1,152-dim) | Genecorpus-104M (V2); Genecorpus-30M (V1) | V2: 103,877,737 from 2,903 datasets. V1: 29,900,531 collected / 27,406,217 passing QC from 561 datasets | whole-atlas, non-cancer (high-mutational-burden cells excluded) | not reported | not reported |
| scPRINT-2 | `small-v2.ckpt` (official `jkobject/scPRINT`) | CELLxGENE-derived (no separate name); "a new database of 54,084,961 cells and 548 datasets" | >50 M; 54,084,961 / 548 datasets after preprocessing 1,200 CELLxGENE datasets | whole-atlas, human+mouse, primary cells only, no spatial | not reported | not reported |
| scGPT | whole-human pretrained (512-dim, `best_model.pt`) | whole-human (CELLxGENE-derived, no formal name) | 33 M normal (non-disease) human cells, 51 organs/tissues, 441 studies; Census release 15 May 2023 | whole-atlas, non-disease only | not reported | not reported |

Full field-by-field version with every verbatim quote and link: `fm_corpora.csv` and `fm_corpora_evidence.csv`.

### 1.3 The measurable bound — CELLxGENE Census composition

Counted directly from Census `summary_cell_counts` for `homo_sapiens`, across every LTS release from the era each model cites through the current release. Cell types by CL identifier: platelet `CL:0000233`; megakaryocyte `CL:0000556`; megakaryocyte progenitor `CL:0000553`; megakaryocyte-erythroid progenitor `CL:0000050`; lung megakaryocyte `CL:4033018`. "MK-lineage" = the four nucleated types summed (platelets excluded).

| Census release | Human cells | Platelets | Platelet % | MK lineage | MK-lineage % |
|---|---|---|---|---|---|
| 2023-05-15 | 53,794,728 | 67,363 | 0.125 % | 108,580 | 0.202 % |
| 2023-07-25 | 56,400,873 | 67,631 | 0.120 % | 109,322 | 0.194 % |
| 2023-12-15 | 62,998,417 | 67,657 | 0.107 % | 109,851 | 0.174 % |
| 2024-05-20 | 74,322,510 | 71,228 | 0.096 % | 109,865 | 0.148 % |
| 2025-01-30 | 109,085,698 | 71,996 | 0.066 % | 123,588 | 0.113 % |
| 2025-11-17 | 163,971,877 | 160,635 | 0.098 % | 128,286 | 0.078 % |

Mapping each model to its own stated release:

- **scGPT** cites the Census release of 15 May 2023 → platelets **0.125 %**, MK lineage 0.202 % of human cells in that release.
- **UCE** cites Census v2023-07-10; the nearest available release is 2023-07-25 → platelets **0.120 %**, MK lineage 0.194 %.
- **scPRINT** draws from the mid-2023 release (800+ datasets) → bracketed by the two rows above, **0.120–0.125 %**.
- **Geneformer** is *not* CELLxGENE-derived (Genecorpus is assembled from its own 2,903-dataset collection), so the Census numbers do not bound it. Its platelet content is **not reported** and cannot be bounded from any published source.
- **STATE** publishes no corpus name or source inventory, so **no bound is possible**; its perturbation component is explicitly cancer-cell-line data, in which platelets cannot occur.

Three caveats to state in the response, because a reviewer will find them otherwise. First, these are *Census* fractions, not corpus fractions: each model subsets Census (scGPT to non-disease; scPRINT to primary human/mouse non-spatial; UCE to 285 of 300+ datasets), and none reports the cell-type composition of its subset, so the Census figure is an upper bound on availability, not a measurement of what was trained on. Second, the fraction is not stable across releases — it drifts from 0.125 % to 0.066 % as the denominator grows, then jumps back to 0.098 % when the 2025-11 releases add ~89 k platelets — so any single number must be reported against a named release. Third, a Census "platelet" annotation is an annotation, not a validated anucleate platelet transcriptome; ambient-mRNA and platelet-monocyte-aggregate artifacts are not separable from these counts.

### 1.4 What this does to the OOD framing

The out-of-distribution argument survives, but its basis shifts from *absence* to *extreme minority plus annotation mismatch*: on the order of 1 cell in 1,000 in the source atlas, no reported representation in any model's own corpus description, and — for STATE — a corpus whose documented perturbation component consists of immortalized cell lines. That is defensible and quotable. "Essentially absent" is not.

---

## 2. Task 2 — source-cohort sampling windows

### 2.1 Verdict

Sampling is genuinely heterogeneous, and not only in width: the 11 cohorts do not share a **timing anchor**. Four different anchors are used (symptom onset, hospital admission/hospitalization, intubation, sepsis diagnosis), and among the cohorts that report a window, the widths span two orders of magnitude — from 6 hours (Qiu, sepsis) to 1–126 days after symptom onset (Ren, COVID). Two cohorts (Stephenson, Lee) report no draw timing at all, in the paper or in the repository record.

![Sampling windows across the 11 source cohorts]({{artifact:art_194b9ba4-b231-425b-9ca8-613ae40f7e21}})

### 2.2 Cohort table

Full 13-column table — accession, disease, donors in the benchmark, verbatim sampling window, design and timepoints with public availability, platform chemistry, severity/outcome definitions, fatal-case timing, per-row verification status, source links: `cohort_sampling_windows.csv`, with all 93 evidence rows in `cohort_sampling_evidence.csv`. Condensed:

| # | Study | Accession | Anchor | Sampling window | Fatal-case timing |
|---|---|---|---|---|---|
| 1 | Qiu 2021, J Leukoc Biol | GSE167363 | sepsis diagnosis | 0 and 6 h post diagnosis | fatality within 24 h (one non-survivor); the other within 30 days |
| 2 | Jiang 2020, JCI Insight | GSE151263 | intubation | within 24 h of intubation / initiation of mechanical ventilation | not reported |
| 3 | Ren 2021, Cell 184 | GSE158055 | symptom onset | paper: not reported as a window; **repository**: 1–126 d, median 23 (IQR 14–42), n = 255/284 | **repository**: 16 deceased-patient samples drawn on days 1, 5, 14, 14, 17, 18, 21, 22, 24, 28, 34, 34, 45, 60, 60, 64 |
| 4 | Stephenson 2021, Nat Med 27 | E-MTAB-10026 | not reported | not reported (only "sampling was performed before inclusion to any pharmacological interventional trials") | not reported |
| 5 | Liu 2021, Cell 184 | GSE161918 | symptom onset | T0–T3, T0 at/near admission, median 11 d post onset; full range not reported | deaths reported (4 of 25 critical) but not timed to draw |
| 6 | Schulte-Schrepping 2020, Cell 182 | EGAS00001004571 (controlled) | symptom onset | cohort 1 (10x): days 3–20; cohort 2 (BD Rhapsody): days 2–25; a figure legend instead gives days 4–13 for first draws | not reported |
| 7 | Bernardes 2020, Immunity 53 | scRNA-seq **not reported** ("to be confirmed"); bulk GSE161777 | hospital admission | days 0, 2, 7, 10, 13 and/or at discharge; sample-naming section instead gives days 1, 3, 8, 11, 14, 20 | second draw of non-survivors 4–35 d before death |
| 8 | Combes 2021, Nature 591 | GSE163668 | hospitalization | screened within 3 d of hospitalization; timing vs symptom onset not reported | not reported |
| 9 | Wilk 2020, Nat Med 26 | GSE150728 | symptom onset | 2–16 d post onset; one donor sampled at 9 d and 11 d | not reported (no fatal cases in subset) |
| 10 | Lee 2020, Sci Immunol 5 | GSE149689 | not reported | not reported in paper or GEO; 3 donors sampled twice, days unstated | not reported |
| 11 | Arunachalam 2020, Science 369 | GSE155673 | symptom onset | paper Methods inaccessible (supplementary-only); **repository**: 2, 8, 9, 9, 15, 15, 16 d (range 2–16, median 9) | not reported |

Row 1 follows the format you pre-confirmed with the authors and is the template for the rest.

### 2.3 Three findings worth putting in the response

**The priority row (Jiang 2020) is resolved.** PMC holds only metadata for this article — no `<body>`, open-access flags negative — so the enrollment and draw timing came from the publisher's free-access full text at insight.jci.org/articles/view/135678. It is a single-timepoint cross-sectional design, drawn within 24 h of intubation, with ARDS by the Berlin Definition and the sepsis-only arm serving as an at-risk comparator. Nothing about it is longitudinal, which matters if any benchmark statement implies otherwise.

**Repository metadata closes three gaps the papers leave open.** Ren's per-sample `sampling day (days after symptom onset)` characteristic in GSE158055 supplies both the cohort-wide window and — this is the reviewer's specific ask — the exact draw days for every sample from a patient recorded as deceased. Arunachalam's `days_since_symptom_onset` in GSE155673 supplies a window for a paper whose Methods are inaccessible. These are machine-read record values, not paper quotes, and are labeled as such in the evidence table; the distinction should be preserved in the response.

**Two cohorts contain internal timing disagreements, which we report rather than reconcile.** Schulte-Schrepping gives days 3–20 (cohort 1) and 2–25 (cohort 2) in the Results, but a figure legend gives days 4–13 as the first-draw window. Bernardes gives days 0, 2, 7, 10, 13 and/or discharge in the Methods — matched by the GSE161777 record — but the sample-naming section gives days 1, 3, 8, 11, 14, 20. Both versions are quoted in the evidence table.

### 2.4 Suggested sentence for the Methods/Limitations

> The 11 source cohorts were not sampled on a common schedule or against a common clinical anchor: draws are timed from sepsis diagnosis (0 and 6 h), from initiation of mechanical ventilation (within 24 h), from hospital admission or hospitalization (within 3 days; or days 0, 2, 7, 10 and 13), or from symptom onset (2–16, 2–25, or 1–126 days, median 23), and two cohorts report no draw timing at all. Because platelet transcriptomes are known to shift with time from illness onset, this heterogeneity is a source of between-cohort variance that our mortality task does not control for, and cross-cohort effect sizes should be read accordingly.

Adjust the emphasis to match how the manuscript already handles cohort as a covariate; every number in that sentence is sourced in the companion tables.

---

## 3. Gap audit — everything not closed, and why

**Not reported (the source genuinely does not state it):**

- Platelet and megakaryocyte counts/fractions in the pretraining corpus — all 5 models. No paper or card publishes per-cell-type composition. Closing this would require re-deriving composition from each corpus's own manifest: feasible for Geneformer (Genecorpus dataset card lists source datasets) and for the Census-derived models (re-run each paper's stated filters against the named Census release); impossible for STATE without a published source inventory.
- STATE corpus name and dataset inventory — absent from the bioRxiv record, GitHub README, and all three model cards.
- Stephenson 2021 — numeric COVID draw window. Days-from-onset appears only as a modeling covariate. Would need the E-MTAB-10026 sample table (SDRF), not retrieved.
- Lee 2020 — any draw timing. Absent from the full text including Methods, and GSE149689 sample characteristics carry only subject group, status, age, gender, cell type. Per-donor day-of-illness may exist in the supplementary tables. This is the expected gap you flagged in advance.
- Jiang 2020, Combes 2021, Schulte-Schrepping 2020, Wilk 2020, Lee 2020, Arunachalam 2020 — fatal-case timing relative to draw. Not stated. For Liu 2021, deaths are reported but not timed to the draw.
- Bernardes 2020 — scRNA-seq accession. The paper states the raw data is in EGA with the accession "to be confirmed"; processed data is on FastGenomics (`beta.fastgenomics.org/p/565003`) and the paired bulk data is GSE161777.

**UNVERIFIED (could not be verified from an accessible primary source):**

- STATE corpus type. The SE component is observational and the ST component is perturbational, but which data the 2,058-dim embedding checkpoint in the benchmark was trained on is not stated unambiguously across the paper and the three model cards. Blocking reason: incomplete model documentation.
- scPRINT `small-v2.ckpt` provenance. The file exists in the official repository listing, but neither the model card nor the Nature Communications paper documents a v2 pretraining run or corpus. The published corpus (54,084,961 cells / 548 datasets) is v1-era. Blocking reason: undocumented checkpoint. If the benchmark used `small-v2.ckpt`, the corpus attributed to it in the manuscript should carry this caveat.
- scGPT corpus wording. Verified against the bioRxiv preprint text; the published Nature Methods wording was not re-checked (paywall). Blocking reason: paywall. The 33 M / 51-organ / 441-study figures are from the preprint.
- Arunachalam 2020 — platform chemistry and severity definition. The PMC record of this Science paper contains no Materials and Methods section; Methods are in supplementary materials at `science.sciencemag.org/content/369/6508/1210/suppl/DC1`, which is not reachable from this environment. GSE155673 supplies sequencer (NovaSeq 6000) and per-sample `disease_severity` (Moderate n = 3, Severe n = 4) but not the 10x chemistry version or the clinical severity criteria. Blocking reason: supplementary-only Methods on a non-allowlisted domain.
- Schulte-Schrepping 2020 — per-sample timing. EGAS00001004571 is controlled-access; only the paper's aggregate windows are available. Blocking reason: controlled access.
- Stephenson 2021 — per-donor days from onset. Requires the ArrayExpress SDRF; the experiment-level record carries no timing. Blocking reason: not fetched (retrievable if wanted).

**Corrections applied during extraction, flagged so they are not reintroduced:** the timepoints "baseline / 90 min / 10 h" in Stephenson 2021 belong to a separate intravenous-LPS challenge arm in healthy volunteers, not to the COVID cohort; and the only fatal-case sentence in Lee 2020 is serological ("the last plasma sample in a lethal case") and gives no sampling-to-death interval, so that field is "not reported" rather than filled.

**Verification tally.** Deliverable 2: 93 evidence rows — 82 exact verbatim substring matches against the retrieved source text, 5 exact after whitespace normalization, 6 machine-read repository metadata (labeled separately, never presented as paper quotes), 0 failures. Census figures were computed from `summary_cell_counts` across 10 releases and are reproducible from `census_platelet_megakaryocyte_by_version.csv`. A file-level inventory of every source retrieved is in `fetch_log.csv` (68 files).
