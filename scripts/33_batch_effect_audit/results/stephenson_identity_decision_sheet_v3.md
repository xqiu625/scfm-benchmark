# Stephenson (E-MTAB-10026) identity decision sheet — v3 (2026-10-11)

Supersedes `stephenson_identity_decision_sheet_v2.md` (2026-10-11, same
day), which a second review (`codex_output/20261011_124253.md`, R2/R5)
found overstated option 3's consequence and the TCR/multi-pool scope — see
the superseded-notice at the top of v2 for exactly what. Only §2's option-3
framing and two wording items in §1 are corrected here; everything else
carries over unchanged from v2. Produced by
`scripts/33_batch_effect_audit/code/audit_stephenson_source_identity.py`
and `scripts/33_batch_effect_audit/code/stephenson_disputed_donor_impact.py`
(commit `2ff4859` and this round's follow-up fixes). Read-only; no canonical
file modified.

## 1. What is confirmed (no owner decision needed)

- 13 of 15 Cambridge `BGCV*` pools have identical donor membership across
  SDRF, the local project table, the real deposited/processed GEX h5ad, BCR,
  and TCR (where TCR observes a pool at all).
- Sanger-site `S9–S16`/`S23–S24` multi-AP-patient sample membership is
  independently reproduced from the TCR file; these are **observed
  Sample-level groupings**, not confirmed physical captures (§4).
- **Corrected scope**: of the 12 Cambridge donors TCR itself observes in
  more than one pool, 11 are a real two-visit `Initial`(D0)/`Resample`(D28)
  pattern; 1 (`CV0902`) is the SAME day/flag in two different pools — a
  different phenomenon, not a resample. This count is scoped to what TCR
  happened to observe, not a complete Cambridge census — §2 below adds a
  13th real multi-visit donor (`CV0201`) that TCR's own coverage gap on
  `BGCV06` had hidden from this count; GEX's own fields see it directly.
- **Corrected scope**: zero BCR per-cell conflicts and zero TCR per-cell
  conflicts (same CellID, two different donors) were found anywhere either
  file observes a Cambridge cell — the check does not extend to either
  file's Sanger/Newcastle rows.
- `run_batch_diagnostics.py`/`step2_run_diagnostics.sh` use `cohort`/
  `patient_id`/`lineage`/`outcome` only — no `cv_group` reference anywhere
  in that call path (re-confirmed this round).
- Recomputing the proposed 21/22/4-donor Cambridge dependency components
  under both the SDRF-declared and the local/BCR/GEX-matching membership for
  all 15 pools changes **zero** donors' component.

## 2. The one genuine, unresolved source conflict

**`BGCV06`/`BGCV13`, donors `CV0201` and `CV0326`.** `BGCV08` is NOT itself
disputed — both SDRF and the local table agree `CV0201 ∈ BGCV08` — but it is
directly relevant evidence, so it is kept in view.

| Source | `BGCV06` | `BGCV08` | `BGCV13` |
|---|---|---|---|
| SDRF | `CV0037;CV0178;CV0234;CV0326` | `CV0073;CV0155;CV0201;CV0915` | `CV0073;CV0171;CV0201;CV0934` |
| Local table | `CV0037;CV0178;CV0201;CV0234` | `CV0073;CV0155;CV0201;CV0915` | `CV0073;CV0171;CV0326;CV0934` |
| Real GEX (`sample_id`, per-cell) | `CV0037;CV0178;CV0201;CV0234` | `CV0073;CV0155;CV0201;CV0915` | `CV0073;CV0171;CV0326;CV0934` |
| BCR (per-cell, observed) | includes `CV0201` (1,816 contig rows) | includes `CV0201` (645 contig rows) | includes `CV0326` (1,230 contig rows), no `CV0201` |
| TCR (per-cell, observed) | no `CV0201`/`CV0326` rows (coverage gap) | `CV0201`, `Initial`/`D0` (974 rows) | no `CV0201`/`CV0326` rows (coverage gap) |

(BCR/TCR row counts above are CONTIG/TSV ROW counts, not deduplicated cell
counts — `stephenson_disputed_donor_impact.py` v2 now reports both
separately; the distinct-cell counts are smaller. This table reports rows
because that is what each source file's granularity actually is.)

GEX's own `Resample`/`Collection_Day` fields show `CV0201` has a genuine,
internally coherent two-visit pattern: `Initial`/`D0` in `BGCV08` (824
cells), `Resample`/`D28` in `BGCV06` (933 cells) — the same D0/D28 resample
shape established for 11 other Cambridge donors (§1). SDRF's claim that
`CV0201`'s second pool is `BGCV13` has no corresponding Resample/Day
structure supporting it in any source checked. `CV0326` is single-pool
(`BGCV13`, `Initial`/`D0`) — not itself a multi-visit donor.

The real per-cell CellID join (GEX vs. BCR/TCR, Cambridge-wide) found a real
DISTINCT union of 94,457 shared cells (21,481 GEX∩BCR, 74,761 GEX∩TCR,
1,785 shared by both), 0 of which disagree with GEX across every source
that shares them.

This remains strong, convergent, but **not dispositive** evidence: GEX,
BCR, and TCR are all OBSERVED, processed, per-cell labels that may share
the same upstream genotype-demultiplexing lineage — their agreement is
real internal-consistency evidence, not proof of independent biological
ground truth. SDRF is the official depositor-submitted record, and the
possibility that `BGCV13_CV0201`'s SDRF registration reflects a real event
this project's copy of the processed object does not represent has not
been ruled out.

**Decision needed from the owner** — three options, with option 3's
consequence corrected this round:
1. Trust the real processed GEX/BCR/TCR convergence (and the new Resample/
   Day coherence) — `CV0201 ∈ {BGCV06, BGCV08}`, `CV0326 ∈ BGCV13`.
2. Trust the SDRF — `CV0201 ∈ {BGCV08, BGCV13}`, `CV0326 ∈ BGCV06`.
3. Mark both donors' `BGCV06`/`BGCV13` assignment `UNVERIFIED`.

**Consequence of each option — corrected this round (a 2026-10-11 review
found the prior wording overstated what was actually shown):**
- Options 1 vs. 2: the proposed 21/22/4-donor Cambridge components are
  identical either way — 0 donors change. This is a clean result: both full
  pool-membership rosters give the SAME partition.
- Option 3 is NOT a clean "no consequence" result, and this sheet does not
  claim it is. What was actually computed
  (`component_diff_excluding_disputed_donors.csv`): deleting `CV0201`/
  `CV0326` from the pool graph entirely leaves the 19 OTHER members of
  their former 21-donor component still connected to each other (component
  sizes become `[22, 19, 4]` — no further fragmentation among the
  retained donors). But this says nothing about what happens to `CV0201`/
  `CV0326` THEMSELVES — their own cells still exist (40 + 33 in Job 1's
  output) and deleting the edge information does not define how those
  cells should be treated in any future CV-dependency grouping: treating
  them as independent singletons, when they may in fact belong to one of
  these real multi-donor pools, would itself be a substantive (and not
  obviously safe) choice, not a neutral default.
- **Interim policy recommended, not a decision made for the owner**: keep
  BOTH source mappings and the conflict flag on record; do not
  auto-exclude `CV0201`/`CV0326` from any pool-dependent grouping; do not
  change `canonical_donor_labels.csv`. The clean options-1-vs-2 result
  above is sufficient to say the 21/22/4 proposal is robust to which of
  those two is eventually chosen; option 3 is not recommended as a
  default simply because it was computationally easy to test.
- These two donors' cells ARE present in Job 1's output (40 cells for
  `CV0201`, 33 for `CV0326`; Platelet + Monocyte_Myeloid) — not academic.
  A real full-CellID join between these exact 73 Job-1-selected cells and
  GEX's own cell_donor map (`job1_gex_cross_check`,
  `stephenson_disputed_donor_impact_v2_20261011_151807/`) found
  **0 unmatched, 0 donor-label mismatches** — every one of these cells'
  Job 1 donor label agrees with GEX's. Their outcome-eligibility (using
  `run_batch_diagnostics.py`'s own `bio_known_mask` definition: the `outcome`
  column IS loaded and inspected to compute this aggregate eligibility count
  — `column_loaded_and_inspected: true` in the real report — but no
  individual per-cell outcome value is printed or reported) is asymmetric: all
  40 of `CV0201`'s selected cells are outcome-`unknown` (0/40 known), while
  all 33 of `CV0326`'s are outcome-known (33/33). This means `CV0201`'s
  cells, as currently selected, would not actually participate in any
  outcome-stratified scoring in Job 2 regardless of this conflict's
  resolution — a real, checked fact, not an assumption either way.

## 3. Job 2 relevance — conditional, not unconditional

`run_batch_diagnostics.py`'s call path (`cohort`/`patient_id`/`lineage`/
`outcome`) never reads pool/`cv_group` identity, so a Job 2 run using
EXACTLY the current input cells, `patient_id`, `outcome`, `lineage`, and
parameters is unaffected by how this conflict is resolved. **This
invariance is conditional, stated precisely**: it holds only for changes
that are confined to the pool-only crosswalk (which pool a donor is
recorded as belonging to) while the input cells, their `patient_id`,
`outcome`, `lineage` values, and Job 2's parameters all remain exactly as
they are now. It does NOT mean the identity question is resolved, and
cross-file label agreement (§2) is not independent proof of identity — a
future donor-relabeling decision (if the conflict is ever resolved toward
option 2, SDRF) could still require different handling, just not of
anything Job 2's current call path reads. See
`scripts/33_batch_effect_audit/results/job2_preparation_packet_v2.md` for
the exact deployment identity, parameters, and acceptance criteria.

## 4. Continuing, non-blocking work

- Full genotype-level identity confirmation (demuxlet/souporcell-style) is
  the most direct further check, not the only one — the original sample
  sheets, a demultiplexing QC report, or a direct question to the
  depositors/ArrayExpress curators could also move this from "strong
  evidence" to "confirmed."
- `annotation_200112.csv` coverage: 128,083/130,644 (98.0%) of real
  Cambridge GEX cells have a matching CellID in the annotation file; 2,561
  do not. Not investigated further.

## 5. Sanger/Newcastle physical-unit caveat

`sanger_observed_membership.csv`'s `Sample` grouping is an **observed TCR
Sample-level co-occurrence**, not a confirmed physical 10x capture. Left
`UNKNOWN`, not asserted either way.

## 6. What this decision sheet is not

Not a recommendation on which source to trust generally, not a statistical
test of identity, not a request to contact the depositors, and not an
authorization to run Job 2, change `canonical_donor_labels.csv`, or promote
the proposed Cambridge grouping. Those remain separate owner decisions.
