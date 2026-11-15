# Longitudinal feasibility note v1 — scope clarification + what the data actually supports

Per `review:20261015_224111`: the independent-donor DE view (25 comparison / 52
reference / 9 excluded, `sdrf::C-8821` resolved into the comparison arm) is **one specific
analysis view**, not a global filter on this project's data. Repeated visits are never deleted
from the underlying `donor × visit × monocyte subtype/state` index; other analyses — scFM-based
cell-state/program discovery, and a genuinely longitudinal within-donor comparison — can and
should use every qualifying visit. **This document does not decide a new primary/secondary
analysis** — it reports what the real repeated-visit structure actually is, so that choice can be
made on evidence rather than by continuing to ask the owner which visit to delete.

**2026-10-16 update, per `review:20261015_230801`:** that review found three real
defects in the first version of this note and its supporting inventory script (F1: the inventory
script and its own test deleted/overwrote the committed v1 result files in place — fixed, v2
output now lives alongside the untouched v1, see `results/stephenson_longitudinal_feasibility_
inventory.py`'s own module docstring; F2: a donor whose `days_from_onset` is the identical
non-numeric string `"Healthy"` was silently counted as a resolved "identical" onset day — fixed,
see §1 below; F4: only the first and last of a 3-visit donor's visits were ever compared — fixed,
the full trajectory is now recorded). It also found two real overstatements in this document
(§1's same-string-pair breakdown undercounted by 2, and §3 claimed a new design type was needed
where the existing one already suffices) and added real evidence this version now incorporates:
a per-sample monocyte-subtype cell-count check against the frozen per-cohort raw object
(`results/stephenson_sample_subtype_cell_counts_v1.tsv`/`_summary.json`) and an actionable
donor-level coverage table (`results/stephenson_longitudinal_monocyte_coverage_v1.tsv`/
`_summary.json`) — see §1a and §4. All corrections below are integrated in place; nothing from
the original version is removed, only corrected or extended.

## 1. What actually exists (per `results/stephenson_longitudinal_feasibility_inventory_v1.tsv`/
`_summary.json`, metadata-only — no expression matrix read)

- **86 true donors total; 71 single-visit; 15 multi-visit** (14 donors with 2 visits, 1 donor —
  `sdrf::C-8896` — with 3). This is the full repeated-visit INVENTORY, not an analysis-ready set:
  it has not been passed through this project's own QC/retention pipeline (still unfinished per
  `docs/current-status.md`) or admitted into any specific comparison.
- **Severity-tier transition and main-comparison-ARM transition are DIFFERENT questions and must
  not be conflated. "Transition" here means a change in the recorded CLINICAL severity label
  between visits — it is not a claim about any measured molecular/cell-state change; no
  expression data has been read for this question.**
  - 12 of 15 are the already-known same-`patient_id`-string Initial/Resample pairs — all
    severity-stable AND arm-stable (3 Asymptomatic→Asymptomatic, 3 Mild→Mild, 2 Moderate→Moderate,
    2 Severe→Severe, **2 Critical→Critical** — corrected 2026-10-16; the original version of this
    bullet omitted these 2 Critical→Critical pairs, undercounting 10 of the 12 — all 12 counted
    among the 13 "stable" pairs below).
  - Of the 3 SDRF cross-string groups: `sdrf::C-8884` (Critical→Critical) is stable in BOTH
    severity and care-setting/arm. `sdrf::C-8896` (Severe→Critical across 3 visits — the FULL
    ordered trajectory is D0 Severe → D7 Critical → D13 Critical) shows a REAL severity-TIER
    worsening but stays care-setting/arm-STABLE (ICU→ICU throughout — Severe and Critical are
    pooled into the same comparison arm, so this donor never crosses the comparison/reference
    boundary). Only `sdrf::C-8821` (Moderate/non-ICU → Critical/ICU) crosses the main-comparison
    ARM boundary itself.
  - **Severity-tier transitions: 13 stable, 1 `Moderate→Critical`, 1 `Severe→Critical`** (15
    total). **Care-setting/arm transitions: 8 `non-ICU→non-ICU`, 6 `ICU→ICU`, 1 `non-ICU→ICU`**
    (15 total) — only 1 donor crosses the arm boundary, confirming the independent-donor view's
    "exactly 1 cross-group donor" finding was already complete, not an undercount.
- **`Days_from_onset` is the IDENTICAL raw string across all visits for every one of the 15
  multi-visit donors — but "identical" and "numeric" are separate facts, corrected 2026-10-16.**
  14 of the 15 have an identical NUMERIC string (`identical_numeric`); the 15th, `sdrf::CV0231`,
  has the identical string `"Healthy"` on both visits (`identical_non_numeric`) — a real,
  non-numeric value, never a resolved onset day, and never to be folded into the 14 real cases.
  0 donors show differing values across visits. The original version of this bullet reported
  `0 unparseable`, which was wrong — CV0231 IS the unparseable case; it was merely counted under
  "identical" without being checked for numeric validity. For the 14 numeric-identical donors,
  the exact generation rule (baseline-only vs. re-measured per visit) remains **unconfirmed** —
  see `MONOCYTE_PILOT_PREREG_v1.md` §4 and the still-open author-code gap in `results/
  stephenson_author_code_investigation_v1.md` §2. A relative ORDERING from `Collection_Day`
  (`D0` < `D7`/`D9`/`D12`/`D13` < `D28`, WITHIN one donor only — never compared across donors)
  is available and was used to order each donor's own visits; absolute onset timing is not.
- **Same-monocyte-state/subtype availability across visits is NOW checked** (2026-10-16, see
  §1a) against the frozen per-cohort raw object's own `orig_celltype` annotation — the roster's
  `subtype_coverage_checked_this_round` column still correctly says "not re-checked this round"
  for the ROSTER table itself, since that check lives in a separate, newly-added artifact.

## 1a. Monocyte subtype coverage across repeated visits (NEW, 2026-10-16)

Per `results/stephenson_sample_subtype_cell_counts_v1.tsv`/`_summary.json` (a bounded,
`obs`-only read of `/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c_monocyte/
stephenson_monocyte_raw.h5ad`, sha256 `11eb7450ff5c2a87158add9b906521cac43db56b71561d75308788089f0f4f4b`
— only `obs['Sample_ID']`/`obs['orig_celltype']` read, never `X`/`layers`; 128,503 cells, 102
samples, 5 subtypes: `CD14_mono`, `CD83_CD14_mono`, `CD16_mono`, `C1_CD16_mono`, `Mono_prolif`)
joined against the roster in `results/stephenson_longitudinal_monocyte_coverage_v1.tsv`/
`_summary.json`:

- **All 102 roster `sample_id`s have an exact match in this object** (0 mapping gaps either
  direction) — confirmed by direct set comparison, not assumed.
- **All 12 symptomatic multi-visit donors (the 15 multi-visit donors minus the 3 Asymptomatic
  ones, listed separately below) share at least one monocyte subtype present (>0 cells) across
  EVERY one of their visits.** `CD83_CD14_mono` is the one subtype shared by all 12 without
  exception; several donors additionally share `CD14_mono` and/or `CD16_mono`; `sdrf::C-8821` and
  `sdrf::C-8896` (the two clinically-transitioning donors) share all 5 subtypes across every
  visit. 0 of the 12 have zero shared subtypes; 0 are `unknown`/incomplete.
- **The 3 Asymptomatic multi-visit donors** (`sdrf::CV0068`, `sdrf::CV0171`, `sdrf::CV0231`) are
  reported separately, never merged into the 12 — they fall outside the current Mild/Moderate vs.
  ICU-Severe/Critical comparison's severity scope. All 3 also have at least one shared subtype
  (`CD83_CD14_mono`, plus `CD14_mono` for 2 of the 3).
- **This is a presence check only** (count > 0 for this specific raw, not-yet-QC'd object) — it
  is not a statement about adequate cell COUNT for a given statistical test, and not a QC-admitted
  cell set. Several visits have low absolute counts (as low as single digits for the rarer
  subtypes); any actual comparison would need its own, separately-stated minimum-cell-count rule.

**2026-10-16 correction (per `review:20261016_102834`): "subtype coverage" and
"usable sample size" are different claims, and the earlier handoff's "subtype coverage is no
longer a limiting factor" language overstated the first as if it settled the second.** Presence
(>0 cells) is confirmed for all 15 multi-visit donors; whether any given visit has ENOUGH cells
of a subtype to support a comparison is a separate, still-open question — e.g. `sdrf::C-8821`'s
D0 visit has only **10** raw `CD83_CD14_mono` cells (`results/stephenson_monocyte_pilot_design_
readiness_donor_visit_subtype_v1.tsv`), next to its D12 visit's 1614.

## 1b. Real per-subtype cell-count sensitivity, NOT a new threshold (NEW, 2026-10-16)

Per `results/stephenson_monocyte_pilot_design_readiness_subtype_cutoff_sensitivity_v1.json`
(built from the SAME already-verified counts above — no new h5ad read): for each of the 5 real
subtypes, how many multi-visit donors have EVERY one of their visits at or above an illustrative
raw-count cutoff. **These are descriptive sensitivity numbers, not a new QC/inclusion rule, and no
single subtype is singled out as "the" subtype to use** — `CD83_CD14_mono` happens to be the only
subtype present in all 15 donors at cutoff 1, but that is a fact about this presence check, not a
decision to restrict the discovery or pilot question to it.

| Subtype | cutoff≥1 (symptomatic / asymptomatic) | cutoff≥20 | cutoff≥50 | cutoff≥100 |
|---|---|---|---|---|
| `CD83_CD14_mono` | 12 / 3 | 6 / 2 | 4 / 1 | 2 / 0 |
| `CD14_mono` | 6 / 2 | 3 / 0 | 3 / 0 | 3 / 0 |
| `CD16_mono` | 6 / 0 | 2 / 0 | 0 / 0 | 0 / 0 |
| `C1_CD16_mono` | 3 / 0 | 1 / 0 | 0 / 0 | 0 / 0 |
| `Mono_prolif` | 2 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |

(Symptomatic range here = the 12 symptomatic multi-visit donors — 11 with exactly 2 visits, plus
`sdrf::C-8896` with 3 — per §1a; asymptomatic = the 3 Asymptomatic multi-visit donors, reported
separately, never merged in.) At `CD83_CD14_mono`≥50, the 4 qualifying symptomatic donors are
`sdrf::C-8884`, `sdrf::C-8896`, `sdrf::CV0073`, `sdrf::CV0284` — `sdrf::C-8896` has 3 visits, so a
design restricted to "exactly 2 visits per donor" (§3, §4(C)) would still need a declared
2-of-3-visit selection rule or a genuinely new design type to include it; it cannot simply be
counted as a 4th 2-visit donor.

## 2. Which longitudinal questions have real support, right now

| Question | Supported? | Caveat |
|---|---|---|
| scFM cell-state/program discovery across ALL qualifying visits | **Yes, metadata + subtype-coverage supports it** | Donor-equal weighting is appropriate for a DESCRIPTIVE population summary (a 3-visit donor must not get 3x the weight of a 1-visit donor). It is NOT a universal fix for every statistical-dependency question — any INFERENTIAL step (a test, a p-value, a fitted trend) must instead use a paired/donor-correlated structure (e.g. the existing `paired_or_repeated_within_donor` design, or a mixed-effects model), never treat two visits of one donor as two independent samples. |
| Within-donor change for the 12 symptomatic + 3 Asymptomatic multi-visit donors, described individually | **Yes, and now with real subtype coverage (§1a)** | 13 of 15 show no CLINICAL severity-tier or arm change at all — only 2 donors (`C-8821`, `C-8896`) show ANY recorded clinical severity-tier transition, and only 1 (`C-8821`) crosses the main comparison/reference arm. All 15 have at least one shared monocyte subtype across every visit, so "which 2 donors are interesting" is a question about the CLINICAL label only, not about which donors have usable cells — a within-donor molecular comparison is not restricted to just these 2. A single donor's trajectory is a case description, not a population claim. |
| A population-level "does state/program X change as patients worsen" claim | **Not supported by this cohort alone** | Only 2 of 86 donors show a CLINICAL severity-tier transition (no molecular/cell-state change has been measured for this question); this is far too few paired observations to support a general population claim from Stephenson alone. Any such claim needs either a much larger within-cohort transition set (not available here) or a different data source. |
| Separating time-since-onset from severity in the current repeated-visit set | **Not supported** | `Days_from_onset` does not change across visits for any of the 15 donors in this dataset (identical raw string every time, including the 1 non-numeric case — see §1), so there is no actual within-donor time axis to separate from severity with the field as currently populated; `Collection_Day` gives only a within-donor schedule ordering, not a disease-time axis comparable across donors. |

## 3. Current code capability vs. this question (read directly, not assumed)

`scripts/34_scfm_discovery_framework/code/de_design_providers.py` lines 28–53 (verified by reading
the file, not from memory): the DE-design preparation module supports exactly two design types —
`independent_donors_unpaired` (one sample per donor) and `paired_or_repeated_within_donor`
(exactly one sample per donor in EACH of two compared groups, fully paired only). General
multi-visit donors (3 visits), partial pairing, and covariate values are explicitly REFUSED, not
silently pooled or dropped. **This is the current code's coverage boundary, not a scientific
reason to discard repeated visits** — the existing refusals stay exactly as they are (they
correctly prevent a repeated sample from being silently treated as an independent donor); per
`review:20261015_224111` item 3, this round does **not** weaken or remove any of those
guards to force repeated visits through the existing unpaired backend.

**Correction, 2026-10-16 (per `review:20261015_230801` F3.3): the original version
of this section incorrectly claimed a new THIRD design type was required even for the simplest
case. It is not.** `paired_or_repeated_within_donor` already builds exactly `~ donor + condition`
(donor fixed effects + one condition dummy) for "every donor contributes exactly ONE sample to
EACH of the two compared groups" — this is DESeq2's own standard paired-sample form. For the 14
donors with exactly 2 visits and no partial-overlap/covariate needs, `condition` can be defined as
`Initial` vs. `Follow-up` (it is not required to equal clinical severity) and the EXISTING design
type already covers this directly. No new capability is needed for that specific case.

A genuinely NEW (third) design type is needed only for what the existing two do not cover:
`sdrf::C-8896`'s 3 visits (an arbitrary, not-exactly-2, number of samples per donor), any donor
set with PARTIAL visit-count overlap (some donors with 2 visits, others with 3, mixed), and
covariate VALUES (still refused by both existing types, per `de_design_providers.py`'s own
documented scope). That third type would need to: build donor fixed effects from an arbitrary
number of visits per donor, handle partial overlap explicitly (refuse or require a declared
inclusion rule, never silently drop), and be estimability-checked the same way the existing two
types already are. This remains a feasibility statement, not an authorization to build it, and
not implemented this round.

## 4. Recommendation for Cindy (one concrete choice, not another "which visit to delete")

Given the above, the two near-term options that are actually well-supported by this data are:

- **(A) Proceed with scFM cell-state/program discovery across the FULL qualifying-visit index**
  (all 102 samples / 86 donor-weighted), reporting per-donor and per-visit structure, independent
  of the already-decided cross-sectional severe-vs-mild/moderate comparison. This directly serves
  the framework's stated primary goal (cell-state/program discovery) and needs no new DE-backend
  capability.
- **(B) A case-level longitudinal write-up of the 2 donors who show a real CLINICAL transition**
  (`sdrf::C-8821`, `sdrf::C-8896`), describing their own state/program change between visits —
  explicitly framed as case description, never a population claim, and never substituting for
  (A) or for the cross-sectional comparison's own result. Both donors share all 5 monocyte
  subtypes across every visit (§1a), so subtype coverage is not a limiting factor for this option.

**(C) — minimal EXECUTABLE within-donor comparison this data could actually support, named here
as a feasibility suggestion only (not fit, not run, not authorized by this note):** restricted to
the EXISTING `paired_or_repeated_within_donor` design (`~ donor + condition`, `condition` =
`Initial` vs. `Follow-up`, per §3's correction — no new backend capability needed). **Correction,
2026-10-16 (per `review:20261016_102834`): the prior version of this bullet said
"the 14 fully-paired 2-visit donors (or any subset)" and named `CD83_CD14_mono` as if it were
settled — both overstated what is actually in scope.** Of the 15 multi-visit donors with
2-visit-each structure (14, per §1) or exactly-2-visit structure within that set, **3 are
Asymptomatic and fall outside the current Mild/Moderate-vs-ICU-Severe/Critical comparison's own
severity scope** (§1a) — the population this design would actually need to draw from, if kept
consistent with that scope, is **11 symptomatic 2-visit donors**, with `sdrf::C-8896` (3 visits)
structurally excluded from a strict 2-visit paired design without a declared selection rule (see
§1b). No subtype is pre-selected here: §1b's real per-subtype sensitivity table shows which of
these 11 (plus C-8896, if a visit-selection rule is adopted) actually have adequate RAW cell
counts for which subtype at which illustrative cutoff — e.g. `sdrf::C-8821` has only 10 raw
`CD83_CD14_mono` cells at its earlier visit, so "shares CD83_CD14_mono across every visit" (true,
§1a) does not by itself mean that visit has enough cells for this donor to usefully enter a
`CD83_CD14_mono`-specific comparison. **This is named as a feasibility sketch of the smallest
concrete design the current evidence could support, not a recommendation to run it, not a subtype
choice, and not a final donor list — those remain open decisions.**

**Time-interval caveat, unresolved and must be stated alongside any such design**:
`Collection_Day` intervals are NOT comparable across donors (D28 is a schedule label, not a shared
calendar day — see `MONOCYTE_PILOT_PREREG_v1.md` §4), so a fitted "Initial vs. Follow-up" effect
would pool donors whose real elapsed time, treatment exposure, and clinical trajectory between
visits may differ substantially; `sdrf::C-8896`'s 3rd visit cannot enter this specific paired
design at all without either a new third design type or a declared, named 2-of-3 visit-selection
rule (never a silent drop). **QC-admission caveat, corrected 2026-10-16 (per `
review:20261016_120026`): the status is `qc_admission_not_verified_for_this_pilot`, NOT "no
QC input exists."** Real QC/decontamination artifacts for Stephenson DO exist — Gate 1 QC
(applied 2026-09-21) and DecontX (complete 2026-09-18, mean contamination 13.0%) both join 100%
by `cell_id` to every one of this object's 128,503 cells (0 mismatches; see
`results/STEPHENSON_QC_SOURCE_CROSSCHECK_v1.md` and `results/stephenson_monocyte_pilot_design_
readiness_table.py`'s `QC_SOURCE_EVIDENCE`, v2 outputs) — but no decision exists that either one
constitutes admission for THIS monocyte pilot specifically, and the decontaminated counts'
actual values were not read (only the per-cell contamination score). Any real run of this design
would first need an explicit owner decision on which of {raw, Gate-1-retained, DecontX-corrected}
counts to use for this pilot, informed by the real per-subtype flag-rate/contamination numbers
now in `STEPHENSON_QC_SOURCE_CROSSCHECK_v1.md` §3 — not an assumption from either object's name.

**Not recommended as a near-term choice**: a population-level "severity progression" DE/statistical
model from this specific repeated-visit set — the data does not currently support it (only 2
qualifying within-donor CLINICAL transitions in the whole cohort), and building new backend
capability for a 2-pair result is not a good use of the current implementation budget. This can be
revisited if a cohort with more repeated-visit transitions becomes available.

This note proposes (A)+(B) together as the recommended next increment, with (C) named as the
minimal executable longitudinal DE design the evidence supports if Cindy wants one; Cindy's
decision is which of (A), (B), (C), some combination, or neither to proceed with — not which
historical visit to discard, and not a request to wait further on time-distribution evidence that
has already been gathered above (§1, §1a).

## 5a. R1/R2 code fixes and the real QC-admission gap (2026-10-16, per `
review:20261016_102834`)

Two real code defects in the new scripts/inventory were found and fixed this round, neither of
which changed any already-committed real number (confirmed by re-running the mutation-killed
tests and re-hashing every v1/v2 output before/after):

- **R1**: `stephenson_sample_subtype_cell_counts.py` indexed categorical codes directly
  (`sample_cats[s_code]`); anndata's own on-disk format uses `-1` for "missing," and Python's
  negative indexing would have silently resolved that to the LAST category instead of raising.
  Reproduced on an isolated fixture; fixed with an explicit code-range validation that refuses
  (never silently remaps) any code outside `[0, len(categories))`. The real object's two `obs`
  fields have 0 missing/out-of-range codes, so this changes no committed count.
- **R2**: `stephenson_longitudinal_feasibility_inventory.py` correctly gated
  `severity_any_adjacent_transition` on `order_resolved`, but NOT the first-to-last fields or the
  summary's first-to-last transition-count dicts — an order-unresolved donor could still get a
  direction label computed from arbitrary input row order. Reproduced on an isolated fixture
  (an unparseable resample day): the unfixed code reported `severity_transition_first_to_last=
  "worsening"` and counted it under `{"Moderate->Critical": 1}` despite `order_resolved=False`.
  Fixed: an order-unresolved donor's first/last/transition fields now read `"unknown"`, and the
  donor is counted under a separate `"order_unresolved"` summary key. The real roster has 0
  order-unresolved multi-visit donors, so this changes no committed count either.
- **"No QC-admitted monocyte analysis input exists yet" (the claim just above, and the prior
  `qc_admitted_input_status` field) was itself overstated — corrected 2026-10-16, per
  `review:20261016_120026`.** That search only covered the repository, not HPCC's
  actual QC data directories. Real, bounded, read-only checks found: Gate 1 QC IS applied for
  Stephenson (2026-09-21) and DecontX IS complete for Stephenson (2026-09-18, mean contamination
  13.0%), and BOTH join 100% by `cell_id` to all 128,503 monocyte-pilot cells (0 mismatches) —
  see `results/STEPHENSON_QC_SOURCE_CROSSCHECK_v1.md` for the full command log and
  `results/stephenson_monocyte_pilot_design_readiness_table.py`'s `QC_SOURCE_EVIDENCE` (v2
  outputs; `qc_admitted_input_status` replaced by the dated, scoped `qc_source_evidence` record).
  **Corrected status: `qc_admission_not_verified_for_this_pilot`** — real QC/decontamination
  artifacts exist and are cell-identity-joinable, but no decision exists that either one
  constitutes admission for THIS monocyte pilot, and the decontaminated counts' actual values
  (not just the per-cell contamination score) were not read. This is a reported verification
  gap, not a claim that no such input could exist, and not a claim that it is already admitted.

## 5. What this note does NOT do

Does not implement any new DE-backend design type. Does not change `stephenson_icu_severity_
grouping_resolve_mixed_visits.py`'s existing, frozen-hash-gated C-8821 resolution (still the
correct view for the already-decided cross-sectional comparison). Does not authorize any new
SLURM job, install, model run, real DE, push, or statistics backend. Does not read the raw
expression matrix (§1a's subtype check reads only two `obs` columns of the per-cohort object, not
`X`). Does not assert the 12/15 donors' shared subtype is a QC-admitted or count-adequate cell
set — only that it is PRESENT (>0 cells) in this specific, not-yet-QC'd raw object; a usable
minimum-cell-count threshold for any real test remains undefined. Does not assert "same monocyte
SUBTYPE present" means "same scFM-discovered STATE" — `orig_celltype` is the authors' own coarse
published annotation, not a scFM output. Does not fit, run, or claim any within-donor molecular
comparison's result — option (C) in §4 is named, not executed.
