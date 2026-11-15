# scripts/33_batch_effect_audit/results — manifest

v1 (2026-10-08), v2 (2026-10-09 morning), and v3 (2026-10-09) are all kept.
Each fixed real defects a review found in the previous version
(`codex_output/20261009_000948.md` for v2, `codex_output/20261009_101831.md`
for v3); see the module docstrings of `build_batch_metadata_inventory.py`,
`build_patient_batch_confounding_table.py`, and the new
`resolve_stephenson_grouping.py` for exactly what changed and why.
**Use v3.** v1/v2 are kept only so each fix is auditable against what it
replaced — this is an audit trail, not a record of correctness.

## Correction to this manifest's own v2 entry

The v2 entry below originally said the 75→58 row drop in Cillo's inventory
was "exactly the de-duplicated Cillo multi-draw rows plus the dropped
literal-nan rows." **That second clause was wrong**, caught by the
2026-10-09 review and confirmed directly (`awk -F',' '$1=="Cillo 2021" &&
$3==""' capture_batch_inventory_v2.csv | wc -l` → 4, same count as v1's
literal-`"nan"` rows). v2 did NOT drop the 4 rows with unknown donor
identity — it fixed their VALUE from the literal string `"nan"` to a
genuine missing value, and they are still present in the file, now as
empty/`NA`. The entire 75→58 reduction is from de-duplicating the Cillo
multi-draw rows alone.

Also corrected: the previous Claude Code handoff
(`claude-code_output/20261009_003147.md`) said "10 pairs resolved" for
Liu's alias normalization. The committed `donor_coverage_report_v2.csv`
correctly reports **19** pairs in its `n_matched_by_alias_normalization`
column; only the `note` column's *display* was truncated to the first 10
for readability. The handoff's prose, not the data, was wrong.

## v1 files (superseded, kept for audit trail only)

- `capture_batch_inventory_v1.csv` — has duplicate rows for any Cillo
  multi-draw library and 4 rows with the literal string `"nan"` as a donor
  ID. Do not use for counting.
- `donor_coverage_report_v1.csv` — reports COMBAT as "0 matched, 124
  only-in-canonical", which reads as a failed join; the real situation is
  that COMBAT's capture-level file has no donor identity column to join on.
- `cv_group_block_size_summary.csv`, `unsplittable_blocks.csv` — these only
  reported the 8 multi-donor `cv_group`s that mix fatal+survived donors, not
  the full set of 12 multi-donor components (4 homogeneous-outcome ones are
  missing), and used the terms "unsplittable block" / implied this
  demonstrates a batch-effect confound. A cv_group is a CV dependency
  component; mixing both outcomes shows outcome overlap within it, which is
  not the same claim as a demonstrated batch-biology association.

## v2 files (superseded by v3, kept for audit trail only)

- `capture_batch_inventory_v2.csv` — fixed the v1 Cillo duplicate-row and
  literal-`"nan"` defects (see correction above for what actually changed).
  **Still has the v2 defect v3 fixes**: Stephenson rows used a semicolon-
  joined, multi-pool `capture_id` with the PATIENT-level cell total written
  into `n_cells_capture` (confirmed: CV0068 shows 4,951 cells under a
  compound "BGCV02;BGCV12" capture_id, when the real per-sample counts are
  2,701 and 2,250 — `T4_outcome_severity_verification/stephenson_sample_severity.csv`).
- `donor_coverage_report_v2.csv` — adds `coverage_status`, exact vs.
  alias-normalized match counts. **Still has the v2 defect v3 fixes**: the
  underlying alias-matching dict comprehensions
  (`{normalizer(d): d for d in ...}`) would silently keep only one ID
  whenever two IDs shared a normalized key, with the survivor depending on
  dict/set iteration order — not observed to have actually dropped a match
  in this specific dataset, but not checked for either, and a real hazard
  in general.
- `cohort_x_event_table_v2.csv`, `cv_group_dependency_component_summary_v2.csv`,
  `all_multi_donor_dependency_components_v2.csv`,
  `mixed_known_outcome_dependency_components_v2.csv` — correct numbers
  (81 fatal / 546 survived / 61 NA across 688 donors; 12 multi-donor
  components, 8 mixed-known-outcome). **Still has the v2 defect v3 fixes**:
  `classify_event` used `np.isclose` (would accept 1.000001 as "fatal");
  donor-uniqueness was checked on bare `donor_id` globally rather than
  `(cohort, donor_id)`.
- `stephenson_pool_vs_cv_group_contradiction_v2.csv` — the finding itself
  (all 15 recorded pools contradict canonical_donor_labels.csv's singleton
  assignment) is confirmed and carries into v3's
  `_v3.csv` file unchanged. What v2 could NOT do — resolve whether a "pool"
  is a real physical capture — v3's new `resolve_stephenson_grouping.py`
  now answers for the Cambridge subset, with a primary-source citation.

## v3 files (current)

- `capture_batch_inventory_v3.csv` — Stephenson now contributes one row per
  (donor, sample_id), `capture_id = sample_id` (matching what
  `build_stephenson_decontx_input.py` itself calls `channel_id`), with the
  correct per-sample_id `n_cells_capture` from the T4 table, and a separate
  `pool_id` column. Duplicate `(cohort, capture_id, donor_id)` keys now make
  the script **raise**, not warn-and-publish.
- `alias_crosswalk_v3.csv` — **new.** Full, collision-aware crosswalk for
  every Liu and Chongqing donor ID (canonical AND capture-file side, not
  just the previously-unmatched subset), each row tagged `exact_match`,
  `normalized_match`, or `AMBIGUOUS_NOT_RESOLVED`. No collision was found in
  this run (0 ambiguous rows for Liu or Chongqing), but the check is now
  real rather than assumed.
- `donor_coverage_report_v3.csv` — now sourced from the crosswalk above;
  `not_assessable` rows (COMBAT) report `pd.NA`, not `0`, for match counts.
- `cohort_x_event_table_v3.csv`, `cv_group_dependency_component_summary_v3.csv`,
  `all_multi_donor_dependency_components_v3.csv`,
  `mixed_known_outcome_dependency_components_v3.csv` — same numbers as v2
  (81/546/61; 12 components, 8 mixed), now from exact-equality event
  classification and `(cohort, donor_id)`-scoped uniqueness, with fixed
  empty-output schemas.
- `stephenson_pool_vs_cv_group_contradiction_v3.csv` — same 15-pool finding
  as v2, re-emitted for the versioned set.
- `stephenson_cambridge_proposed_dependency_groups_v1.csv` — **new.** Built
  by `resolve_stephenson_grouping.py`. Primary-source-confirmed (Stephenson
  et al. 2021 *Nat Med*, Methods: Cambridge "pools of four samples were
  generated by combining 0.5 million live cells per individual" BEFORE
  antibody staining and 10x loading, with downstream genotype-based donor
  splitting) plus a local cross-check (every pooled donor's own `Site` field
  reads `"Cambridge"`): the 15 recorded pools ARE real physical multi-donor
  captures for these donors, and `canonical_donor_labels.csv`'s singleton
  `cv_group` for every one of them does not reflect that. This file proposes
  corrected connected-component groupings (union-find over the pool
  membership graph) for the Cambridge pool members **only** —
  `canonical_donor_labels.csv` itself is NOT modified; promoting this
  proposal is a separate decision. Non-Cambridge Stephenson donors (Sanger,
  Newcastle) have no pool evidence in local metadata at all and are left as
  unknown, not asserted independent.

  **Result:** the 47 Cambridge pool-member donors (currently 47 independent
  singletons in `canonical_donor_labels.csv`) collapse into only **3**
  connected components once cross-pool donor recurrence is followed
  transitively (sizes 21, 22, 4) — a Stephenson analog of COMBAT's single
  124-donor block, at smaller scale. Checked against `event`: the two large
  components **mix fatal and survived donors** (CV0025-component: 1 fatal /
  16 survived / 4 unknown; CV0050-component: 2 fatal / 13 survived / 7
  unknown). If this proposal is adopted, the project's mixed-known-outcome
  dependency-component count would rise from 8 to **10** (8 existing + the
  2 new large Cambridge components), and the multi-donor component count
  from 12 to **15** (12 existing + all 3 new Cambridge components, including
  the 4-donor one that stays single-outcome) — corrected 2026-10-09 from an
  arithmetic slip that said 14; a review caught it (12 + 3 = 15, not 14).
  Not applied here.

- **Addendum (2026-10-10, 2nd pass): the `BGCV06`/`BGCV13` conflict has now
  been recomputed both ways — the proposed components above are UNCHANGED
  either way.** `audit_stephenson_source_identity.py` cross-checked all 15
  Cambridge pools against the SDRF, the local table, the REAL deposited/
  processed GEX h5ad's own per-cell `sample_id`, BCR, and TCR. Exactly 2
  pools disagree: SDRF alone says `CV0201∈BGCV13`+`CV0326∈BGCV06`; the local
  table, BCR, AND the real GEX data all agree on the opposite. A full-CellID
  per-cell join between GEX and BCR/TCR found zero disagreements across
  96,242 shared Cambridge cells — not dispositive, but real evidence, not
  mere vote-counting. **`build_dependency_components`/`diff_component_mappings`
  recomputed the proposed grouping under BOTH the SDRF-declared and the
  local/BCR/GEX-matching membership for all 15 pools: 0 donors change
  component either way** (`BGCV06`, `BGCV08`, `BGCV13` are already bridged
  by other, unaffected donors — e.g. `CV0073`/`CV0171` sit in both `BGCV08`
  and `BGCV13` regardless of where `CV0201`/`CV0326` land). The 21/22/4
  proposed-component sizes above are therefore safe to read as-is regardless
  of how the `CV0201`/`CV0326` conflict is eventually resolved. **Job 1 did
  not use this grouping** and is unaffected either way. See
  `claude-code_output/` for exact result paths
  (`component_diff_sdrf_vs_local.csv`, `component_diff_summary.json`).

- **Addendum (2026-10-11, 3rd pass): a review
  (`codex_output/20261011_113512.md`) found the 2nd pass's "0/96,242
  disagreements" figure summed two overlapping cell-set intersections
  (21,481 GEX∩BCR + 74,761 GEX∩TCR) without checking their overlap — the
  real distinct union is 94,457 (1,785 cells shared by both), still 0
  disagreements.** It also found GEX had been treated as a 5th DECLARED
  source rather than one more OBSERVED, processed label (fixed:
  `gex_processed` moved to `OBSERVED_SOURCES`) and a legacy-encoding `-1`
  missing-category bug that could silently mislabel a cell as the LAST
  real category (fixed, with a regression test). "Option 3 (exclude the
  disputed donors) has zero consequence" was asserted in the 2nd pass
  without being tested — it is now actually tested
  (`component_diff_excluding_disputed_donors.csv`): excluding `CV0201`/
  `CV0326` shrinks the 21-component to 19 with NO fragmentation (sizes
  `[22, 19, 4]`), so the remaining 19 donors stay exactly one group either
  way. New evidence: GEX's own `Resample`/`Collection_Day` fields show
  `CV0201` has a real `Initial`/`BGCV08` -> `Resample`/`BGCV06` two-visit
  pattern, and BCR directly observes `CV0201` in `BGCV06` and `BGCV08`
  (never `BGCV13`) — see
  `scripts/33_batch_effect_audit/results/stephenson_identity_decision_sheet_v2.md`
  for the full corrected writeup. None of this changes the 21/22/4
  proposed-component sizes above.
