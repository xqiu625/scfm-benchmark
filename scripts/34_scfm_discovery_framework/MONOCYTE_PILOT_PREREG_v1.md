# Monocyte severity pilot — preregistration draft v1

Per `review:20261014_204536` §5.A.3, revised per `review:20261014_221155` R2,
`review:20261014_231259` B, and now **`review:20261015_120031` and
`review:20261015_120317`** — Cindy's direct decisions in the automated review conversation,
recorded in `docs/decisions.md` (2026-10-15). Drafted from these decisions plus the real, dated
feasibility evidence in `results/monocyte_pilot_feasibility_v3.tsv` (supersedes `v2.tsv`'s
Stephenson row — the 90 raw `patient_id` strings were further found to OVERCOUNT true biological
donors; see §2). **This is a
draft for owner review, not an executed plan** — every field below is explicitly `confirmed` (the
owner's own recorded decision, or a fact independently verified from source this round),
`proposed` (this document's recommendation, not yet decided), or `pending_owner` (an open
question this document cannot resolve; §3 distinguishes `pending_owner` from `pending_evidence` —
a fact that is checkable and should be checked, not merely awaited). No real DE run, state
discovery, or pathway enrichment happens this round.

**The main comparison is now `confirmed` and LOCKED (2026-10-15, superseding the 2026-10-14
"severe vs. mild+moderate" / "critical pending" framing below wherever it conflicts):
comparison arm = severity in {Severe, Critical} AND confirmed ICU location at that SAME sample's
collection time; reference arm = severity in {Mild, Moderate}, regardless of location; direction
= comparison minus reference.** Original per-sample severity labels are always preserved
alongside the derived group. A sample whose severity is Severe/Critical but whose location is
non-ICU or unknown is EXCLUDED from both arms, never folded into the reference arm and never
guessed as ICU or non-ICU. See `results/stephenson_icu_severity_grouping.py` (the mapping, derived
only from the real `Status_on_day_collection` string, never from the severity grade itself) and
`results/stephenson_icu_severity_grouping_v1.tsv`/`_summary.json` (the per-sample/per-donor
evidence). Only `critical`'s merge-vs-separate handling and "ICU restricted or not" were ever
open — both are now closed; §2/§3/§7 below are updated accordingly. What remains open is **visit
selection** (§6) and the items in §2/§5/§8/§9 that this decision does not touch.

## 1. Question — `confirmed`

Within ONE COVID-19 cohort, do monocyte cell states/gene programs differ between **severe** and
**mild/moderate** disease. Patient/donor is the biological replicate. This selects the question;
it does not reinstate the former mortality benchmark, does not pool COVID-19 and Sepsis, and does
not imply any particular monocyte subtype scope by itself.

## 2. Cohort — `proposed`, pending the Stephenson verification below being accepted

Per `results/monocyte_pilot_feasibility_v2.tsv`, exactly two COVID-19 cohorts currently have a
curated, atlas-derived monocyte cell-type extraction in this pipeline's existing artifacts
(`option_c_monocyte/`): **Stephenson** and **Ren**. A third, Qiu-sepsis, is excluded by disease.
No other admitted COVID-19 cohort (Wilk, Lee, Arunachalam, Liu, Combes) has a monocyte-cell-type
extraction built yet — extracting one is real engineering work, out of this round's scope.

**Stephenson's numbers are now independently re-verified against the cohort's OWN raw object**
(`covid_portal_210320_with_raw.h5ad`, read-only, obs-only, `h5py`, 2026-10-14), per
`review:20261014_221155` R2 — superseding the earlier (incorrect) numbers. **A FURTHER
correction this round** (`review:20261014_231259` B): the 90 `patient_id` STRINGS found
in that extraction are not all genuinely distinct biological donors either.

- **CORRECTED AGAIN, true-donor level, n=86** (not 90): using Stephenson's OWN depositor-submitted
  SDRF (`E-MTAB-10026.sdrf.txt`, `Characteristics[individual]` — the SAME crosswalk resource
  `build_monocyte_canonical.py`'s own `load_stephenson_individual_map()` already reads) to group
  `patient_id` strings by their REAL biological individual, 3 confirmed hidden same-donor groups
  were found that no amount of `patient_id`-string matching could reveal: `{MH8919327, MH8919328,
  MH9143273}` (3 strings, 1 real donor, SDRF individual `C-8896`), `{MH9179821, newcastle004v2}`
  (2 strings, 1 real donor, `C-8821`), `{newcastle21, newcastle21v2}` (2 strings, 1 real donor,
  `C-8884`).
- **86 is the SDRF-CONFIRMED donor count — not a candidate pending further resolution.** Per
  `review:20261015_104226` §三, a prior handoff conflated three different quantities (85,
  87, 86) that measured different things; this round's rewrite of
  `stephenson_sdrf_identity_crosswalk.py` makes the roster and the summary JSON a single producer
  so they cannot diverge, and resolved the discrepancy with evidence rather than a wording choice.
  Two samples (`BGCV06_CV0201`, `BGCV13_CV0326`) have a Source Name string absent from the SDRF —
  the already-documented BGCV06/BGCV13 pool-membership conflict in
  `stephenson_identity_decision_sheet_v3.md`. Both resolve to a confirmed SDRF individual anyway,
  via the SAME patient_id's OTHER SDRF row(s) (`CV0201` via `BGCV08_CV0201`/`BGCV13_CV0201`;
  `CV0326` via `BGCV06_CV0326`), since every SDRF row sharing a given patient_id agrees on the
  individual regardless of pool label (checked across the FULL SDRF — 0 patient_id suffixes show
  internally conflicting individuals). All 102 cohort samples now resolve to a confirmed SDRF
  individual; there are zero remaining identity-unresolved samples. **The one thing that remains
  open is NOT the donor count** — it is WHICH POOL/batch `BGCV06_CV0201` and `BGCV13_CV0326`
  themselves belong to for dependency-grouping purposes (a question `stephenson_identity_decision_
  sheet_v3.md` already tracks and has not yet resolved). A narrow regression check
  (`results/check_stephenson_roster_consistency.py`) now fails if the roster and summary ever
  disagree about the identity-unit count again.
- **Corrected severity distribution, true-donor level, Initial-visit preferred, n=86**:
  Asymptomatic=9, Mild=23, Moderate=30, **Severe=13, Critical=11** (down from the patient_id-level
  Critical=15). The drop matters biologically, not just numerically: for 2 of the 3 collapsed
  groups, severity genuinely CHANGED between visits (e.g. one donor was `Severe` at their Initial
  draw and `Critical` only by a later resample) — the Resample-only "singleton" record that looked
  like an ordinary independent Critical patient was actually that SAME donor's LATER, deteriorated
  state, not their sampling-day severity. Using the Initial visit for every collapsed group (not
  the Resample-only record) is the policy applied here, consistent with §6's proposed MVP rule.
- **Candidate arithmetic, true-donor level** (per the review's own request, explicitly NOT a frozen
  final analysis N): severe alone = 13 vs. mild+moderate = 53; severe+critical = **24** vs.
  mild+moderate = 53 (this SUPERSEDES the earlier, patient_id-level "28 vs 53," which
  double/triple-counted the 3 collapsed groups as independent Critical patients).
- **Join uniqueness (cell-level), confirmed**: 128,503 of 128,503 existing monocyte-canonical
  Stephenson cells matched 1:1 to the fresh raw-file extraction via the EXACT crosswalk
  `build_monocyte_canonical.py` itself uses (`steph_<covid_index>`) — 0 unmatched, 0 duplicate
  crosswalk IDs on either side, 0 `sample_id`/`patient_id` conflicts.
- **Site is a real, confirmed 3-way confound**: Cambridge/Newcastle (`Ncl`)/Sanger — exact
  per-true-donor counts not recomputed this round; the `patient_id`-level counts (Cambridge=36,
  Ncl=43, Sanger=11) are a close approximation, not re-verified at the true-donor level.
- **Full evidence trail**: `results/stephenson_sample_level_real_fields.csv` (102 real sample-level
  rows), `results/stephenson_raw_metadata_extraction_summary.json`/`.py` (the raw-object
  extraction, hardened this round — see "Engineering fixes" note below), `results/
  stephenson_sdrf_identity_crosswalk_summary.json`/`.py` (the SDRF crosswalk), and `results/
  stephenson_donor_identity_roster_v1.tsv` — a per-sample DRAFT roster (original ID, confirmed/
  unresolved identity status, true donor group, visit, time, site, severity, inclusion rationale)
  for the owner's review, per the review's explicit request. **This roster does not itself
  authorize treating any row as an admitted independent-donor input, and does not change
  canonical identity labels or QC retention** — it is evidence for the owner's eventual decision,
  not a decision.
- **Engineering fix to the extraction itself, per `review:20261014_231259` A.2**: the
  metadata-extraction script's categorical decoder had a latent bug (the newer anndata-style
  group encoding did not handle a missing code (`-1`) the same safe way its legacy-encoding branch
  already did, silently wrapping a missing value to the LAST real category via numpy's negative
  indexing) — fixed, and confirmed this round that NONE of the actual columns this extraction
  reads ever hit a missing/out-of-range code, so this round's own real findings above are
  unaffected by the (now-fixed) latent bug. The script's duplicate/unmatched/identity-conflict
  checks and its internal-consistency check (now covering all 10 kept clinical fields, not 3) are
  now hard failures, not print-and-continue diagnostics.

**Ren's verification was ATTEMPTED this round and found genuinely blocked locally, not merely
skipped**: `monocyte_source/ren/`'s only files are the GEO single-cell matrix/barcode/annotation
files (`GSE158055_cell_annotation.csv.gz`'s own columns are exactly `cellName, sampleID, celltype,
majorType` — checked directly, 2026-10-14) — **no severity/clinical field exists anywhere in the
locally-present Ren raw data**. Ren's severity value must originate from the platelet-level
integration's own separate source (`scripts/data/integrate_platelet_data_harmony.R`, cited in
`docs/cohort-sources.md` but not itself present in this repository, and presumably from Ren et
al. 2021's own supplementary clinical table, not downloaded into `monocyte_source/ren/`). This is
reported as a genuine, checked blocker, not an unchecked assumption: Ren's severity numbers
in `monocyte_pilot_feasibility_v2.tsv` remain `UNVERIFIED`, and resolving this would require
either locating Ren et al.'s own supplementary metadata table (a new file, out of this round's
read-only-existing-files scope) or accepting the platelet-level join as the best available source
for Ren specifically (unlike Stephenson, where a better, native source was locally available and
used instead).

- **Proposed primary: Stephenson.** Now the better-evidenced candidate on two counts: the larger,
  three-tier-populated, independently-reverified severity distribution (true-donor level,
  Mild=23/Moderate=30/Severe=13/Critical=11), and a confirmed-clean cell-level identity join plus
  a genuine donor-identity crosswalk this round — not merely the larger raw patient count.
- **Proposed alternate: Ren**, with the explicit caveat that its severity numbers remain
  unverified, and — unlike Stephenson — cannot be verified from locally-present raw files at all.
- **`pending_owner`**: final cohort choice, informed by the asymmetry above.

## 3. Inclusion/exclusion and QC version — `pending_owner`

- Disease: COVID-19 only (no Sepsis, no influenza, no SLE/off-disease controls) — `confirmed` by
  the owner's question framing.
- Schulte-Schrepping and Bernardes remain **permanently excluded as data** (`docs/decisions.md`,
  2026-09-19) — `confirmed`, not reopened by this document.
- QC version: the broader platelet-level QC rebuild is unfinished per `docs/current-status.md`.
  This pilot's own input admissibility status is `provisional` for the SAME reason — `confirmed`
  as a constraint, with the consequence `pending_owner`: whether a provisional-QC pilot is
  acceptable for a feasibility/method pilot while the QC rebuild is pending, or whether the pilot
  must wait.
- **Severity-value handling, Stephenson, RESOLVED by Cindy's 2026-10-15 decision**:
  `Status_on_day_collection_summary` takes 5 real values for this cohort — `Asymptomatic`,
  `Mild`, `Moderate`, `Severe`, `Critical` (none literally named "fatal" or "recovered"; a
  SEPARATE collapsed `covid_severity` field, derived from the outcome `y`, is what previously
  supplied "fatal"/"survived" — never Stephenson's own per-sampling-day label, and still never
  used to design/select/tune this comparison). **`Critical` is now `confirmed` merged into the
  comparison arm together with `Severe`** (not a separate tier, not folded into "severe" as a
  renaming — both retain their original label alongside the derived group), **conditional on
  confirmed ICU location at sampling** (see the file header and §7). `Asymptomatic` patients
  remain `confirmed` EXCLUDED from both arms (never folded into "mild") — this was already this
  module's behavior before the ICU decision and is unchanged by it. Fatality/outcome is
  `confirmed` NOT used in this comparison.
- Ren's severity-value vocabulary was NOT independently checked this round (per §2, its raw files
  carry no severity field at all) — **`pending_owner`/`pending_evidence`** jointly: which
  vocabulary Ren even uses cannot be answered without finding its source metadata first.

## 4. Severity mapping and time window — `pending_owner`, now with real supporting fields identified

- **RESOLVED by Cindy's 2026-10-15 decision**: `Status_on_day_collection_summary` (sampling-day
  severity) together with `Status_on_day_collection` (the raw ICU/ward string used to derive
  confirmed ICU location at that same sampling day) is the `confirmed` field pair for this
  pilot's comparison — not `Worst_Clinical_Status` (peak-of-illness severity), which the raw
  object also carries but is explicitly NOT used to decide group membership. Later ICU
  admission, worst-ever clinical status, and ventilation status alone are each explicitly
  insufficient to prove ICU location at sampling; `care_setting_at_sampling` in
  `results/stephenson_icu_severity_grouping.py` is derived only from the literal
  `Status_on_day_collection` value recorded for that sample.
- **Confirmed this round**: `Days_from_onset` is a real, mostly-numeric field for this cohort (97
  of 102 samples numeric; the 5 non-numeric values are exactly `'Healthy'` (CV0231, both visits —
  consistent with being `Asymptomatic`), `'nan'` (CV0144, CV0198), and `'Not_known'` (newcastle49)
  — never blanket-cast to 0 or dropped silently). `Collection_Day` (`D0`/`D7`/`D9`/`D12`/`D13`/
  `D28`) is a SEPARATE sampling-SCHEDULE label, not days-from-onset — `D0` does NOT mean day 0 of
  illness. Whether a time-window restriction on `Days_from_onset` is needed for this pilot remains
  **`pending_owner`**, now with the real field identified and ready to filter on once decided.
- Whether `Status_on_day_collection_summary`/`Worst_Clinical_Status` reflects the SAME severity
  concept Stephenson et al. 2021's own Methods calls "severity" — not yet checked against the
  paper's text this round — remains **`pending_owner`** (a literature read, not a metadata pull).

## 5. Monocyte subtype scope — `proposed`, pending_owner confirmation

Both candidates have REAL subtype resolution available, but only in their PER-COHORT raw files
(`ren_monocyte_raw.h5ad`, `stephenson_monocyte_raw.h5ad`), not in the merged
`monocyte_canonical_raw.h5ad` currently in use (its `obs` does not carry `orig_celltype`):

- Stephenson `orig_celltype`: `CD14_mono` (58562 cells), `CD83_CD14_mono` (51735), `CD16_mono`
  (13619), `C1_CD16_mono` (3942), `Mono_prolif` (645).
- Ren `orig_celltype`: `Mono_c1-CD14-CCL3` (18605), `Mono_c2-CD14-HLA-DPB1` (28895),
  `Mono_c3-CD14-VCAN` (47599), `Mono_c4-CD14-CD16` (2936), `Mono_c5-CD16` (6043).

**Proposed, not decided**: a classical/CD14-dominant restriction (pooling the CD14-labeled
clusters in either cohort) is a real, available option, not a default — the owner must choose
between "all monocytes" and a classical/CD14 restriction; this document does not pick one, and
whichever is chosen must be re-joined from the per-cohort raw files (a small, real engineering
step, not yet done).

## 6. Repeated-visit handling and donor-as-biological-replicate — `confirmed` (method and incidence)

- **Scope clarification (2026-10-15, per `review:20261015_224111`): the 25/52/9
  independent-donor comparison view below (§7) is ONE SPECIFIC analysis view, never a global
  filter on this project's data.** `sdrf::C-8821`'s resolution into the comparison arm (using
  `newcastle004v2`, excluding `MH9179821` from THIS view) was decided under the explicit premise
  "if it can only enter one group" — it governs only the independent-donor cross-sectional
  design, and does not authorize deleting any other visit from the underlying `donor × visit ×
  monocyte subtype/state` index, nor does it forbid studying within-donor change over time. See
  `mono_LONGITUDINAL_FEASIBILITY_NOTE_v1.md` (new, this round) for what the actual repeated-visit
  structure supports: 86 true donors, 15 with >1 visit (14 pairs + 1 triple,
  `sdrf::C-8896`), of which only 2 (`sdrf::C-8821`, `sdrf::C-8896`) show ANY recorded CLINICAL
  severity-tier transition (no molecular/cell-state change has been measured for this question)
  and only 1 (`sdrf::C-8821`) crosses the comparison/reference arm boundary itself — far too few
  to support a population-level progression claim from this cohort alone, but sufficient for
  case-level description and for full-index scFM state/program discovery (which does not require
  picking one visit per donor at all). **Update 2026-10-16** (per `
  review:20261015_230801`, fixed in the inventory script and this note): the 12 same-string
  pairs' breakdown had undercounted by 2 (`Critical→Critical`); `sdrf::CV0231`'s identical
  `days_from_onset` value is the non-numeric string `"Healthy"`, not a resolved onset day; and all
  15 multi-visit donors (12 symptomatic + 3 Asymptomatic) now have a confirmed shared monocyte
  subtype across every visit (`results/stephenson_longitudinal_monocyte_coverage_v1.tsv`) — see
  `mono_LONGITUDINAL_FEASIBILITY_NOTE_v1.md` §1/§1a for the corrected numbers and evidence.
  **Second update, 2026-10-16** (per `review:20261016_102834`): "subtype present"
  is not "usable sample size" — e.g. `sdrf::C-8821`'s earlier visit has only 10 raw
  `CD83_CD14_mono` cells. Real per-subtype cell-count sensitivity and a donor×visit×subtype
  table now live in `results/stephenson_monocyte_pilot_design_readiness_table.py`'s outputs.
  See `mono_LONGITUDINAL_FEASIBILITY_NOTE_v1.md` §1b/§4(C)/§5a.
  **Third update, 2026-10-16** (per `review:20261016_120026`): the second update's
  "no QC-admitted input exists" was itself overstated (repository-only search, missed HPCC's
  actual QC directories). Corrected: Gate 1 QC (applied 2026-09-21) and DecontX (complete
  2026-09-18) both exist for Stephenson and join 100% by `cell_id` to all 128,503 monocyte-pilot
  cells — status is `qc_admission_not_verified_for_this_pilot`, not "absent." See
  `results/STEPHENSON_QC_SOURCE_CROSSCHECK_v1.md` for the full cross-check.
- **Confirmed method, carried from the existing framework interface
  (`scripts/34_scfm_discovery_framework/code/de_design_providers.py`)**: an independent-donor
  cross-sectional comparison (one sample per donor) uses the `independent_donors_unpaired`
  design (intercept + condition only) and must NOT add a per-donor fixed effect — a fixed effect
  per donor in a design with one sample per donor saturates the model (zero residual degrees of
  freedom), which `de_design_providers.py` already refuses explicitly (see its "full rank is
  never itself treated as adequate power" check). A donor fixed effect is only appropriate for a
  GENUINELY repeated/paired design (`paired_or_repeated_within_donor`), where it absorbs
  between-donor variance using within-donor pairs. Per `review:20261014_221155` R3's
  explicit caution: this pilot's actual target comparison is BETWEEN donors (severe vs.
  mild/moderate is a between-patient contrast, not a within-patient before/after), so even where
  repeated visits exist, a donor fixed effect is the WRONG tool here regardless — it would absorb
  the very between-donor severity effect this pilot wants to measure. Repeated visits are a
  SELECTION problem (pick one sample per donor), not a reason to switch design types.
- **CONFIRMED incidence this round, Stephenson, REVISED** (`results/
  stephenson_sample_level_real_fields.csv` + `results/stephenson_sdrf_identity_crosswalk_summary.json`):
  **"12 of 90 patients have a repeated visit" undercounted the true picture.** By `patient_id`
  STRING alone, 12 of 90 contribute exactly 2 samples under the SAME string (`Resample` design:
  `Collection_Day` `D0` "Initial" + a later `Resample`); severity was consistent across both
  visits for all 12 of these specific pairs. But per `review:20261014_231259` B, the
  SDRF crosswalk found 3 MORE real same-donor groups hidden behind DIFFERENT `patient_id` strings
  (see §2) — so the true number of real donors with more than one sample in this subset is **15**
  (12 same-string pairs + 3 cross-string groups), out of **86** true donors, not "12 of 90." Two
  of the 3 cross-string groups showed a REAL severity change between visits (unlike every one of
  the 12 same-string pairs, which were severity-consistent) — repeated-visit severity consistency
  is NOT a universal property of this cohort, only of the easily-matched same-string pairs. Ren's
  incidence was not checked this round (blocked by the same missing-severity-field issue in §2).
- **Proposed, result-independent MVP selection rule, now applied in the roster above**: always
  take the `Resample == 'Initial'` sample for any TRUE donor with more than one (using the SDRF
  identity crosswalk to find the true donor, not just matching `patient_id` strings) —
  deterministic, decided BEFORE looking at any severity/expression result, never chosen
  per-donor to favor an outcome. This is a proposal for the owner to confirm, not yet applied to
  any DE-ready artifact (only to the severity-distribution recount above).
- **Why this is not cosmetic, concrete evidence (2026-10-15, `results/
  stephenson_icu_severity_grouping_v1_summary.json`)**: applying the now-confirmed comparison
  (§7) to EVERY visit (not yet picking one per donor) found exactly one true donor,
  `sdrf::C-8821`, whose visits land in BOTH the comparison arm and the reference arm (one visit
  is Mild/Moderate, the other is Severe/Critical-and-ICU). Donor-level counts under "any visit
  counts" (before resolution): comparison=25, reference=53, excluded-other-severity=9, out of 86
  true donors.
- **`sdrf::C-8821`, RESOLVED (2026-10-15)**: per `docs/decisions.md`, Cindy directly resolved this
  ONE donor — "如果只能进入一组, 就进入目前分组少的那组" (if it can only enter one group, assign
  it to whichever group currently has fewer members) — into the comparison arm, via sample
  `newcastle004v2` (D12, Critical, ICU); its other sample, `MH9179821` (D0, Moderate, non-ICU), is
  excluded from this round's reference arm (raw data/identity/QC unaffected). This is a one-time
  decision under the current roster, NOT a general tie-break rule, NOT the authors' own analysis
  rule, and NOT derived from any expression result — see `results/
  stephenson_icu_severity_grouping_resolve_mixed_visits.py` and its `_resolved_v1_summary.json`
  for the exact hash-pinned before/after counts. **Resolved donor-level counts, each donor now in
  exactly one group: comparison=25, reference=52, excluded-other-severity=9** (86 total). This
  does NOT answer the separate, still-open visit-selection question below (which single sample
  represents a donor who has more than one qualifying visit WITHIN the same, already-resolved
  group) — `sdrf::C-8821` has only one ICU-qualifying sample, so is unaffected by that question
  either way.
- Cells and repeated samples/visits are never treated as independent donors, regardless of which
  design type is ultimately chosen — this is a hard constraint on the existing interface, not a
  pilot-specific choice.

## 6b. Site/center and known pool-identity confounds — `confirmed` (facts), `pending_owner` (handling)

- **Confirmed this round**: Stephenson's 90 monocyte-annotated patients span THREE real sites —
  Cambridge (36), Newcastle/`Ncl` (43), Sanger (11). Per `review:20261014_221155` R3's
  caution, "one admitted cohort" does not mean "no technical/clinical center confound" — Site is a
  real categorical variable available to check against (and, if needed, control for or stratify
  by), not an assumption either way about whether it actually confounds the severity contrast.
- **Confirmed this round**: donors `CV0201`/`CV0326` (the one open conflict in `scripts/
  33_batch_effect_audit/results/stephenson_identity_decision_sheet_v3.md` — which of two Cambridge
  sequencing pools, `BGCV06`/`BGCV13`, each is really drawn from) ARE among these 90 patients.
  Their severity/`patient_id` labels are unaffected by this conflict either way (0 donor-label
  mismatches found across every independently cross-checked source, per that decision sheet) — it
  only matters if/when a POOL-level batch or dependency grouping is built for this pilot, which
  has not been done this round.
- **`pending_owner`**: whether Site needs to be modeled, stratified on, or simply reported as a
  covariate availability note for this pilot's MVP scope; whether the `CV0201`/`CV0326` pool
  conflict needs resolving before any pool-level dependency grouping is built for this specific
  analysis (it does not block the severity/`patient_id`-only MVP design above).

## 7. Design matrix / contrast direction — `confirmed`, updated 2026-10-15

**(Severe + Critical, each confirmed at ICU location at sampling) minus (Mild + Moderate, pooled
as the single reference group) is this pilot's MAIN, confirmed comparison** — per Cindy's direct
decisions in `review:20261015_120031` (Severe+Critical pooling) and
`review:20261015_120317` (the further ICU-at-sampling restriction, her explicit choice
"B"), both recorded in `docs/decisions.md` 2026-10-15 and NOT reopened by this document:

```text
comparison_group = severity in {Severe, Critical} AND confirmed ICU location at sampling
reference_group  = severity in {Mild, Moderate}
direction         = comparison_group - reference_group
```

A Severe/Critical sample without confirmed ICU location at sampling is excluded from BOTH groups
(§3/§6) — it is never treated as non-ICU by default and never moved into the reference group.
This direction is recorded consistently in the design's `comparison_level`/`reference_level`/
`direction_definition` (per `gene_ranking_providers.ComparisonSpec`) and in any report — never
silently flipped; `results/stephenson_icu_severity_grouping.py` has a dedicated test
(`tests/test_stephenson_icu_severity_grouping.py`) asserting the comparison and reference labels
stay distinct and that ICU location is never re-derived from the severity label alone. A
SEPARATE, secondary/exploratory breakout of mild vs. moderate individually, or severe vs.
critical individually, MAY be reported alongside the main pooled contrast if useful, but never in
place of it, and never framed as reopening the grouping decision. The future sepsis-ICU vs.
COVID-19-ICU comparison (`docs/execution-plan.md` 2026-10-15) is a SEPARATE, later design stage —
it does not pool into this contrast and is not implied to be comparable by this decision alone.

## 8. Positive control, expression/model/PCA evidence, validation unit — see companion docs

- Positive control candidate and the literature conflict it surfaced: `POSITIVE_CONTROL_
  EVIDENCE_v1.md`. Direction not yet frozen — `pending_owner`.
- Model/checkpoint status: `MODEL_DECISION_v1.md`. Not yet confirmed — `pending_owner`.
- **Expression/PCA baseline vs. model-derived evidence must be kept SEPARATE and compared, never
  conflated** — per `review:20261014_204536` §2.1: a pseudobulk DE result computed from
  the SAME raw counts does not depend on which scFM is used, so it cannot by itself distinguish
  one model's contribution from another's; a model's contribution must show up in model-specific
  structure (state/embedding clusters, masking-based gene ranking), checked SEPARATELY against
  the expression-level evidence, not substituted for it.
- **Validation unit is the donor/patient, never the cell.** Any held-out split for validating a
  model-selection or masking-score decision must be at the donor level; cells from the same
  donor are never split across train/validation. If `ComparisonSpec`/the existing framework
  interfaces are used to build any such split, the split itself is a `pending_owner` decision,
  not implied by this document.

## 9. Pre-specified success / failure / indeterminate conditions — `proposed`, pending_owner

**Terminology kept identical to `POSITIVE_CONTROL_EVIDENCE_v1.md`'s own three-way outcome
language (per `review:20261014_231259` C), so the two documents cannot be read as
disagreeing about what counts as which:**

- **Success (feasibility, not biological claim)**: the pilot can build a real, full-rank,
  adequate-residual-df design through `de_design_providers.py` for the chosen cohort/subtype
  scope, and (once a real DE backend exists, per `DE_BACKEND_ENVIRONMENT_PLAN_v1.md`) the
  positive-control program's calibration-or-direction check (per `POSITIVE_CONTROL_EVIDENCE_v1.md`
  §4's fallback) is met.
- **Failure, distinguished explicitly from "no biological effect"**: an insufficient-residual-df
  design (per `de_design_providers.py`'s existing refusal), a positive control that fails even a
  calibration-only check (suggesting an input/pipeline problem, not necessarily absence of
  biology), or a donor-identity/QC problem surfaced during preparation.
- **Indeterminate** (same category `POSITIVE_CONTROL_EVIDENCE_v1.md` calls "a non-detected
  program" / "a null result"): a positive control that fails a DIRECTIONAL check but passes a
  calibration-only check, OR simply does not detect the program at all — given the real
  cross-cohort disagreement documented in `POSITIVE_CONTROL_EVIDENCE_v1.md` (never itself called a
  "direct contradiction," per that document's own R3 correction), this is a real, anticipated
  possible outcome, never automatically relabeled as pipeline Failure, and never upgraded to
  Success by relaxing the pre-specified direction after seeing the result.
- **Explicitly NOT a success/failure criterion**: sample size alone. Synthetic/small replicate
  counts used in this framework's own unit tests are never a stand-in for this pilot's real
  statistical power — per the review's own caution, "合成重复数从不是 Cindy 的真实样本量裁定."

## 10. What this document does NOT do

Does not run any real DE, state discovery, or pathway enrichment. Does not pick the cohort,
subtype scope, severity-value handling, or positive-control direction — each is marked
`pending_owner` above. Does not authorize a new SLURM job, install, or push.
