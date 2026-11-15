# Benchmark design decisions (2026-09-15)

Goal that settles them: a BIB benchmark of single-cell foundation-model embeddings on uniformly
decontaminated cells, evaluated with no patient leakage. Two consequences run through everything
below. First, a benchmark compares models **on the same labels**, so label noise that is shared
across models shifts absolute AUC but largely preserves ranking — which makes an imperfect pooled
endpoint acceptable in a way it would not be for a clinical claim. Second, "no patient leakage"
in a pooled multiplexed dataset is a stronger requirement than patient-level splitting, because
patients who share a GEM well share ambient RNA.

## Decision 1 — endpoint: all-cause death during study follow-up

Seven of the ten labelled cohorts record only a binary outcome with no death date, so a 28-day
primary endpoint is not constructible without discarding them. Use all-cause death during study
follow-up as the primary label; carry each cohort's native window in `event_definition` and
report it as a supplementary table; run a 28-day sensitivity on the three cohorts that support it
(COMBAT `Death28`, GSE279452 `mort28`, GSE174072 `days_to_death`) and a 90-day one where
available (GSE279452, GSE180578).

Under the leave-one-cohort-out design below, the label is internally consistent within each
evaluation fold, which blunts the heterogeneity objection further: a reviewer's concern is that
window differences confound *comparisons between cohorts*, and LOCO never pools cohorts inside a
test set.

## Decision 2 — strata: two within-disease tasks, pooled as secondary

| task | donors | events |
|---|---|---|
| within-COVID (primary) | 497 | 53 |
| within-sepsis, incl. paediatric (primary) | 138 | 27 |
| pooled critical illness (secondary) | 721 | 86 |

Influenza (2 events) and non-COVID ARDS (1) cannot support their own AUC — keep them out of the
primary tasks and inside the pooled secondary. The 3 non-septic ICU-control deaths belong in the
pooled task only; folding them into the sepsis task would break the disease definition. Report
the pooled task as a generalization check, not as the headline, because a model can score on it
by learning disease rather than mortality.

## Decision 3 — platelet threshold: two tasks, both pre-registered

The scFM embeddings are computed per cell over all lineages, so the whole-PBMC benchmark can use
all 86 events. The platelet-restricted analysis — the project's biological thesis — should use
**≥20 platelet cells per donor: 75 of 86 events** (COVID 47, sepsis 21+2, ICU control 3,
influenza 1, ARDS 1), with sensitivity at 10 and 50. Fix both before looking at any model output;
the threshold is the real denominator of the platelet claim and is the obvious place for a
reviewer to suspect post-hoc selection.

## Leakage design — the part that needs building, not deciding

**Patients who share a GEM well share ambient RNA.** If donor A is in train and donor B is in
test and their cells came out of the same emulsion, B's cells carry A's transcripts. That is a
real information channel between folds, and it is exactly what upfront decontamination reduces —
so the decontamination step is not only QC here, it is a leakage-mitigation step, and the paper
should say so.

`canonical_donor_labels.csv` now carries two grouping columns computed from the channel
inventories: `well_block` (the pool/batch/library a donor was loaded in) and `cv_group`
(connected components of the "shares a well" graph — the minimum unit that can move between
folds without splitting a well).

- **546 of 721 donors are singleton `cv_group`s** — nothing else shares their well, so they are
  freely assignable. (A separate 536 donors carry `well_block == "singleton"`, meaning no well
  record was found for them at all; the two counts are not interchangeable.)
- **Liu** forms 3 blocks (14, 10, 9 donors; 4 events) — one per batch, since donors recur across
  that batch's lanes.
- **Cillo** forms 8 small blocks of 2–3 donors (8 events), from the hashtag-multiplexed
  libraries.
- **COMBAT collapses into a single 124-donor block carrying all 16 of its events**, because 13
  donors were loaded into more than one pool (up to 3), chaining all ten pools together. COMBAT
  therefore cannot be split internally under any well-respecting scheme.

That last point makes the choice for you:

**Primary evaluation: leave-one-cohort-out.** The group is the cohort, so every well stays intact
by construction, cohort-batch effects cannot be learned as a shortcut, and the claim being tested
— does an embedding transfer to a cohort the model has never seen — is the one a benchmark reader
cares about. Held-out event counts: GSE279452 22, COMBAT 16, Cillo 12, Ren 12, Stephenson 8,
Wilk 5, Liu 4, Schulte 3, Bernardes 2, Qiu 2.

**Secondary: grouped stratified 5-fold on the pooled set using `cv_group`**, which quantifies
how much of the LOCO performance drop is cohort shift rather than sample size. Note COMBAT's
block forces an imbalanced fold; report the fold composition rather than hiding it.

Three mechanical rules to enforce in code, since they are where benchmarks usually leak:

1. Every fitted transform — HVG selection, scaling, PCA-200, PseudobulkEN, any FM fine-tuning —
   is fit on training folds only. Zero-shot frozen embeddings are exempt; the baselines are not.
2. No integration (Harmony, scVI) fit across train and test together.
3. All samples and timepoints of a patient move as a unit, and so does the patient's whole
   `cv_group`.

## The leakage risk you cannot fix by splitting

Several of these cohorts are in the public corpora the scFMs were pretrained on — CELLxGENE
includes COMBAT, Ren, Stephenson and Teichlab among others — so test patients' cells may already
be in a model's pretraining data, and that advantage differs between models, which is precisely
what the benchmark measures. This is the standard objection to scFM benchmarks and a BIB reviewer
will raise it. It cannot be solved by fold design; it has to be documented per model × cohort,
and the project already has a corpus-verification plan artifact
(`plan_verify-scfm-pretraining-corpus-platelet_dd0417db.json`) to build on. Treat the overlap
table as a required figure, not an appendix note.

## Decision 4 — channel scope: 8 of 10 cohorts, Schulte-Schrepping and Bernardes excluded

Approved 2026-09-19 (Cindy). Phase 1 upfront decontamination covers 8 cohorts (Liu, Ren, COMBAT,
Cillo, Qiu, Chongqing, Stephenson, Wilk). Schulte-Schrepping (EGA controlled-access,
`EGAS00001004571`) and Bernardes (scRNA-seq accession unresolved in the source paper itself,
processed data gated behind a FastGenomics account) are **excluded from the benchmark**, not
deferred — the user will not file the EGA Data Access Request or register on FastGenomics.

Rationale: the canonical donor/event crosswalk (721 donors, 86 fatal events) was already fixed
before Phase 1 began, and the binding constraint on the benchmark's statistical power is the 86
events, not the channel count — Decision 3 above already notes a full 5-FM x 6-baseline x 3-task
DeLong matrix will not sustain many contrasts at this event count. Adding Schulte-Schrepping (3
events) and Bernardes (2 events) cannot expand the comparison budget enough to justify the access
friction. This must be stated explicitly in the methods as an exclusion with a reason, not left as
a silent gap.

## Decision 5 — embedding extraction: subsampled, not full-scale

Approved 2026-09-19 (Cindy). FM embedding extraction (Geneformer, scGPT, UCE, scPRINT-2, STATE)
uses the SAME stratified subsampling strategy validated in the Phase 3 pilot (cap 60 cells/donor,
per `geneformer_loco_before_after_decontx.py`'s convention) across the full 8-cohort Platelet +
Monocyte_Myeloid dataset — NOT full-scale extraction over all ~1.2M assembled cells.

Rationale: the Phase 3 before/after-decontamination pilot (PCA-200 + Geneformer, Liu+Ren+COMBAT
subset) already showed conclusions are stable under this subsampling — decontamination did not
flip the headline AUC for either baseline. Full-scale extraction across 5 models would burn
substantial GPU queue time without changing the direction of the answer this benchmark is built to
give. This decision should be revisited only if a specific downstream analysis (e.g. the
platelet-activation-axis variance-retention check, approved as a post-rerun mechanistic follow-up,
not a benchmark input) turns out to need full-scale embeddings.

## Milestone: mid-October 2026, September target cancelled

Approved 2026-09-19 (Cindy). The original September delivery target is cancelled. New milestone:
mid-October 2026, synchronized with the multiomics delivery. Proceed at normal pace — do not cut
corners (e.g. skipping subsampling validation, rushing statistical comparisons) to hit a September
date that no longer applies.

Priority order for remaining work: (a) wire `cv_group` into the actual CV callers (leftover TODO,
hours-scale); (b) 5-FM embedding extraction under subsampling; (c) full LOCO runs (5 FM x
baselines x 3 tasks: within-COVID / within-sepsis / pooled); (d) statistical comparisons (DeLong
etc., keeping the primary-comparison budget small per Decision 3); (e) results integration for the
manuscript (manuscript writing itself is the user's own work, not this session's).