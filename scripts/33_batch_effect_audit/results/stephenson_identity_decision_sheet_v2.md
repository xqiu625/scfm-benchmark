# Stephenson (E-MTAB-10026) identity decision sheet — v2 (2026-10-11)

**SUPERSEDED by `stephenson_identity_decision_sheet_v3.md` (2026-10-11,
same day).** A second review the same day
(`codex_output/20261011_124253.md`, R2/R5) found: §2's "option 3 has zero
consequence for everyone else" overstated what the component-size result
(`[22, 19, 4]`, `n_donors_with_changed_component: 21`) actually shows —
removing 2 members from a 21-member group does not fragment the remaining
19 further, but that is NOT the same claim as "every option is
consequence-free," since deleting the two donors leaves undefined how
their OWN cells should be treated in any future CV-dependency grouping.
§1 also said "zero TCR per-cell conflicts... anywhere either file observes
a cell," overstating scope — the check is Cambridge-only — and did not
flag that the "12 multi-pool/11 resample" count is itself scoped to what
TCR happened to observe, not a complete census (the new GEX evidence in
§2 found a 13th real multi-visit pattern, `CV0201`, that TCR's own
coverage gap had hidden from that count). v3 fixes both; §1's other items
and §2's core evidence table are NOT known to be wrong and carry over
with only wording/scope corrections, not new numbers.

Supersedes `stephenson_identity_decision_sheet_v1.md` (2026-10-10), which a
2026-10-11 review (`codex_output/20261011_113512.md`, F1) found overclaimed
in five places — see the superseded-notice at the top of v1 for exactly
what. This version corrects those five and adds real evidence v1 did not
have (GEX's own `Resample`/`Collection_Day` fields, §2). Produced by
`scripts/33_batch_effect_audit/code/audit_stephenson_source_identity.py`
(commit `fde3eea` and later) and the new
`scripts/33_batch_effect_audit/code/stephenson_disputed_donor_impact.py`.
Read-only; no canonical file modified. Full evidence in
`stephenson_identity_crosswalk_v2_20261011_123120/`,
`stephenson_disputed_donor_impact_20261011_123120/`, and
`job1_validation_addendum_v3_20261011_123120/`.

## 1. What is confirmed (no owner decision needed) — unchanged from v1

- 13 of 15 Cambridge `BGCV*` pools have identical donor membership across
  SDRF, the local project table, the real deposited/processed GEX h5ad, BCR,
  and TCR (where TCR observes a pool at all).
- Sanger-site `S9–S16`/`S23–S24` multi-AP-patient sample membership is
  independently reproduced from the TCR file; these are **observed
  Sample-level groupings**, not confirmed physical captures (§4).
- Of the 12 Cambridge donors seen in more than one pool, 11 are a real
  two-visit `Initial`(D0)/`Resample`(D28) pattern; 1 (`CV0902`) is the SAME
  day/flag in two different pools — a different phenomenon, not a resample.
- Zero BCR per-cell conflicts and zero TCR per-cell conflicts (same CellID,
  two different donors) were found anywhere either file observes a cell.
- `run_batch_diagnostics.py`/`step2_run_diagnostics.sh` use `cohort`/
  `patient_id`/`lineage`/`outcome` only — no `cv_group` reference anywhere
  in that call path (re-confirmed this round).
- Recomputing the proposed 21/22/4-donor Cambridge dependency components
  under both the SDRF-declared and the local/BCR/GEX-matching membership for
  all 15 pools changes **zero** donors' component (re-confirmed this round,
  `component_diff_sdrf_vs_local.csv`).

## 2. The one genuine, unresolved source conflict — corrected evidence

**`BGCV06`/`BGCV13`, donors `CV0201` and `CV0326`.** `BGCV08` is NOT itself
disputed — both SDRF and the local table agree `CV0201 ∈ BGCV08` — but it is
directly relevant evidence (below), so it is kept in view rather than
dropped.

| Source | `BGCV06` | `BGCV08` | `BGCV13` |
|---|---|---|---|
| SDRF | `CV0037;CV0178;CV0234;CV0326` | `CV0073;CV0155;CV0201;CV0915` | `CV0073;CV0171;CV0201;CV0934` |
| Local table | `CV0037;CV0178;CV0201;CV0234` | `CV0073;CV0155;CV0201;CV0915` | `CV0073;CV0171;CV0326;CV0934` |
| Real GEX (`sample_id`, per-cell) | `CV0037;CV0178;CV0201;CV0234` | `CV0073;CV0155;CV0201;CV0915` | `CV0073;CV0171;CV0326;CV0934` |
| BCR (per-cell, observed) | includes `CV0201` (1,816 cells) | includes `CV0201` (645 cells) | includes `CV0326` (1,230 cells), no `CV0201` |
| TCR (per-cell, observed) | no `CV0201`/`CV0326` rows (coverage gap) | `CV0201`, `Initial`/`D0` (974 rows) | no `CV0201`/`CV0326` rows (coverage gap) |

Three corrections to v1's framing of this evidence:

1. **GEX is not a fifth independent source.** It is one more OBSERVED,
   processed, per-cell label — internally consistent with BCR/TCR (which may
   share the same upstream genotype-demultiplexing lineage), not proof of
   independent biological ground truth. What IS now independently checked
   (this round): GEX's own `patient_id` obs field agrees with the
   `sample_id`-derived donor on **all 130,644** real Cambridge cells, 0
   mismatches (`stephenson_disputed_donor_impact_20261011_123120/
   disputed_donor_impact.json`) — a real internal-consistency check, not an
   assumption.
2. **TCR's silence on `BGCV06_CV0201`/`BGCV13_CV0326` is a coverage gap, not
   a vote for either side.** It does, however, positively OBSERVE `CV0201`
   in `BGCV08` (974 rows, `Initial`/`D0`) — consistent with every other
   source. BCR positively OBSERVES `CV0201` in BOTH `BGCV06` (1,816 cells)
   AND `BGCV08` (645 cells), and `CV0326` ONLY in `BGCV13` (1,230 cells,
   never `BGCV06`) — this is real, per-cell, observed evidence directly
   contradicting SDRF's specific claim for these two donors, not merely a
   subset/coverage-gap pattern.
3. **New evidence v1 did not have**: GEX's own `Resample`/`Collection_Day`
   fields (read directly this round, not assumed) show `CV0201` has a
   genuine, internally coherent two-visit pattern: `Initial`/`D0` in
   `BGCV08` (824 cells), `Resample`/`D28` in `BGCV06` (933 cells) — the SAME
   D0/D28 resample shape already established for 11 other Cambridge donors
   (§1). SDRF's claim that `CV0201`'s second pool is `BGCV13` (not `BGCV06`)
   has no corresponding Resample/Day structure supporting it in any source
   checked. `CV0326` shows only one pool (`BGCV13`, `Initial`/`D0`) in GEX —
   it is not itself a multi-visit donor.

The real per-cell CellID join (GEX vs. BCR/TCR, Cambridge-wide, not limited
to the disputed donors) found: 21,481 GEX∩BCR shared cells, 74,761
GEX∩TCR shared cells, with **1,785 cells shared by all three** — so the
real DISTINCT union is **94,457 cells**, not 96,242 (v1's unchecked sum of
the two intersection sizes). All 94,457 agree with every source that shares
them; 0 disagree (`stephenson_identity_crosswalk_v2_20261011_123120/
gex_cell_level_crosswalk.json`, `combined` block).

This remains strong, convergent, but **not dispositive** evidence: SDRF is
the official depositor-submitted record, and the possibility that
`BGCV13_CV0201`'s SDRF registration reflects a real, separate event this
project's copy of the processed object does not represent has not been
ruled out.

**Decision needed from the owner** — same three options as v1, restated:
1. Trust the real processed GEX/BCR/TCR convergence (and the new Resample/
   Day coherence) — `CV0201 ∈ {BGCV06, BGCV08}`, `CV0326 ∈ BGCV13`.
2. Trust the SDRF — `CV0201 ∈ {BGCV08, BGCV13}`, `CV0326 ∈ BGCV06`.
3. Mark both donors' `BGCV06`/`BGCV13` assignment `UNVERIFIED`, excluding
   them from any pool-dependent grouping decision.

**Consequence of each option, now actually computed (v1 asserted this
without testing option 3):**
- Options 1 vs. 2 (§1): the proposed 21/22/4-donor Cambridge components are
  identical either way — 0 donors change.
- Option 3, tested this round
  (`component_diff_excluding_disputed_donors.csv`/`_summary.json`): removing
  `CV0201`/`CV0326` entirely from the pool graph shrinks the 21-member
  component to 19 (component sizes `[22, 19, 4]`, still exactly 3
  components — no fragmentation). The remaining 19 donors stay together as
  ONE group; no donor outside `{CV0201, CV0326}` changes which other donors
  it is grouped with. So option 3 also has zero consequence for everyone
  else's grouping — the only real-world effect of any option is which
  donor→pool fact gets written for these two donors specifically.
- A fourth, newly concrete fact (not an "option," a measurement): these two
  donors' cells ARE present in Job 1's output (40 cells for `CV0201`, 33 for
  `CV0326`; Platelet + Monocyte_Myeloid;
  `stephenson_disputed_donor_impact_20261011_123120/
  disputed_donor_impact.json`). Job 2 (§3 below) does not read pool/`cv_group`
  identity, so this presence does not depend on resolving §2's conflict —
  but it means the conflict is not academic: real selected cells exist for
  both donors.

## 3. Job 2 relevance — see `job2_preparation_packet_v1.md`

That packet (new this round) states the exact deployment identity, checksums,
parameters, and the real fact that `run_batch_diagnostics.py`'s call path
never reads pool/`cv_group` identity, so Job 2's result (if run) does not
depend on how §2 is resolved. It also documents and fixes (in a new,
undeployed wrapper file) a `backed="r"` metadata-only-read defect inherited
from the original R1 finding.

## 4. Continuing, non-blocking work — unchanged from v1

- Full genotype-level identity confirmation (demuxlet/souporcell-style) was
  not attempted. **Correction to v1**: this is the MOST DIRECT further
  check, not the ONLY one — v1's original sample sheets, a demultiplexing
  QC report, or a direct question to the depositors/ArrayExpress curators
  could also move this from "strong evidence" to "confirmed," and are not
  ruled out as alternatives.
- `annotation_200112.csv` coverage: 128,083/130,644 (98.0%) of real
  Cambridge GEX cells have a matching CellID in the annotation file; 2,561
  do not. Not investigated further.

## 5. Sanger/Newcastle physical-unit caveat — unchanged from v1

`sanger_observed_membership.csv`'s `Sample` grouping is an **observed TCR
Sample-level co-occurrence**, not a confirmed physical 10x capture. Left
`UNKNOWN`, not asserted either way.

## 6. What this decision sheet is not

Not a recommendation on which source to trust generally, not a statistical
test of identity, not a request to contact the depositors, and not an
authorization to run Job 2, change `canonical_donor_labels.csv`, or promote
the proposed Cambridge grouping. Those remain separate owner decisions.
