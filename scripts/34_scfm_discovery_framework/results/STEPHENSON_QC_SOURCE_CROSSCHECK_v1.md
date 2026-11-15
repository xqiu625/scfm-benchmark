# Stephenson monocyte pilot — QC source cross-check and minimal pilot readiness note

**Date:** 2026-10-16. **Per:** `codex_output/mono_20261016_120026.md`, items 2 and 4. **Scope:**
a bounded, read-only check of whether the current monocyte pilot's inputs can be joined by
cell/sample identity to the project's existing QC/decontamination work, and what that does and
does not settle. No large matrix was read; no model was fit; no job was submitted.

## 1. What was checked (commands, dates, scope)

All commands were run directly on the HPCC login node (`bluejay`), read-only, 2026-10-16:

| # | What | Size actually read | Result |
|---|---|---|---|
| 1 | `stephenson_monocyte_raw.h5ad` `obs/_index`, `obs/Sample_ID`, `obs/orig_celltype` | 3 small `obs` columns (~1.5 MB total), never `X` | 128,503 cells, index format `steph_<covid_index>` |
| 2 | `05_RESULTS/rebuild_upfront_qc/stephenson_decontx_input/cells.csv` | full file, 38 MB (not the 12 G `counts.mtx`) | 647,366 rows, `cell_id`/`channel_id`/`donor_id`/`harmonized_lineage`/`draw_time` |
| 3 | `.../stephenson_decontx_contamination_scores.csv` | full file, 44 MB (not the 25 G `stephenson_decontx_counts.mtx`) | per-cell `decontX_contamination` score |
| 4 | `05_RESULTS/rebuild_upfront_qc/gate1_qc/stephenson/qc_flags_per_cell.csv` | full file, 132 MB | per-cell `flag_would_delete`/`n_flags` and the individual QC flags it is built from |
| 5 | `gate2_scrublet/` (directory listing only) | `ls` | only `stephenson.INVALID_postqc_prior_20260921` exists — explicitly named INVALID |
| 6 | `scripts/31_donor_level_revision/code/build_monocyte_canonical.py`, `scripts/32_rebuild_upfront_qc/code/build_stephenson_decontx_input.py`, `qc_caveats.md`, `qc_checklist.md` | read, not executed | source/provenance of every file above |

**Corrected 2026-10-16 (per `codex_output/mono_20261016_130032.md` R3): the original version
of this line pointed to `claude_output/mono_20261016_*.md` as "the full command log" -- that
file is prose summarizing this check, not a reproducible implementation of it.** The actual
reproducible evidence is `results/stephenson_qc_source_crosscheck.py` (a committed, tested
script; `main()` re-derives every number in §2-3 below) and its real output,
`results/stephenson_qc_source_crosscheck_v1.json` (produced by running that script against
the real paths above, 2026-10-16) — re-running the script reproduces the same numbers, not
just restates them. `stephenson_monocyte_pilot_design_readiness_table.py`'s `QC_SOURCE_EVIDENCE`
constant is a dated summary of these same findings, not an independent check.

## 2. Cell-identity join — CONFIRMED, exact, 100%

Stripping the `steph_` prefix from `stephenson_monocyte_raw.h5ad`'s obs index gives the EXACT
`cell_id` string used in both `stephenson_decontx_input/cells.csv` and
`stephenson_decontx_contamination_scores.csv` (verified on the `AP1` channel directly, then for
all 128,503 cells): **0 of 128,503 monocyte-pilot cells failed to match**, and `channel_id`
agreed with the pilot's own `Sample_ID` for every matched cell (0 mismatches). The same join
against Gate 1's `qc_flags_per_cell.csv` also matched 100% (0 missing). This means every
monocyte-pilot cell already has: a Gate 1 per-cell QC flag record, and a DecontX per-cell
contamination score — **this is a real, checked fact, not an assumption.**

**Not read this round:** the decontaminated counts themselves (`stephenson_decontx_counts.mtx`,
25 G) and the pre-DecontX full counts (`counts.mtx`, 12 G). The row/column alignment between
`cells.csv` and those two matrices is documented directly in
`build_stephenson_decontx_input.py` (`cells_df`/`combined` built in the same per-channel loop,
written in the same order) — not independently re-verified by reading the matrices this round.

## 3. What this does and does not settle

**Settled, corrected 2026-10-16 (per `codex_output/mono_20261016_130032.md` R1 — the original
wording below overstated this):** HISTORICAL Gate 1 QC artifacts (dated 2026-09-21) and
HISTORICAL DecontX contamination scores (dated 2026-09-18, predating Gate 1) exist FOR
STEPHENSON and are cell-identity-joinable to every cell this pilot currently uses. This is NOT
the claim that a current, admitted, end-to-end QC pipeline is complete for this project's
present purposes — Gate 2 (Scrublet) has no valid Stephenson result (only a file explicitly
named `INVALID`), and whether any later DecontX refit or retention policy supersedes these two
historical stages is not addressed here. "Raw" in `stephenson_monocyte_raw.h5ad`'s own name
describes its PROVENANCE (`layers['raw']` from the authors' published atlas —
`build_monocyte_canonical.py:187`), not an absence of any QC; the authors' own atlas may carry
its own upstream cell-calling, which this round did not check.

**NOT settled (genuinely open):**
- Whether Gate 1's `flag_would_delete` and/or DecontX-corrected counts SHOULD be applied to this
  monocyte pilot specifically — no such decision exists anywhere in the repository. Gate 1's
  flag rate varies by monocyte subtype (real, this round's join): `CD83_CD14_mono` 0.8%,
  `CD14_mono` 0.6%, `CD16_mono` 0.8%, `C1_CD16_mono` 3.6%, `Mono_prolif` 1.1%.
- Whether raw or decontaminated counts are the right input for a DE test vs. an scFM
  preprocessing step — these have different input-requirement contracts (the DE backend's own
  interface contract, `de_design_providers.py`, should be consulted directly when this is
  decided; not assumed from either object's filename).
- DecontX contamination differs meaningfully by subtype (real, this round's join, mean %):
  `CD14_mono` 6.0%, `CD83_CD14_mono` 5.4%, `Mono_prolif` 9.4%, `CD16_mono` 29.7%,
  `C1_CD16_mono` 29.9% — this is descriptive evidence relevant to a future subtype choice, not a
  subtype choice itself.
- Gate 2 (Scrublet doublet detection) has no current, valid result for Stephenson — only a file
  explicitly marked `INVALID`. A mechanical doublet candidate is, per `CLAUDE.md`'s standing
  rule, never itself a deletion authorization even once a valid run exists.
- The broader platelet-level QC rebuild's overall completion status (`docs/current-status.md`)
  is a separate, larger question this note does not re-litigate or need to wait on — Gate 1 and
  DecontX are each independently complete FOR STEPHENSON specifically, which is what this pilot
  needs.

**Status, corrected from the prior round's overstated "no QC-admitted input exists":**
`qc_admission_not_verified_for_this_pilot`. Real QC/decontamination artifacts exist and are
cell-identity-joinable; no admission decision for THIS pilot has been made by anyone.

## 4. What can proceed now vs. what needs a decision first

**Can proceed now (engineering, no new scientific decision needed):** a bounded script that
joins Gate 1 flags and DecontX contamination scores onto the existing coverage/design-readiness
tables by the confirmed `cell_id` key, so any future pilot decision has real flag-rate/
contamination numbers per donor/visit/subtype in hand immediately, instead of needing a fresh
join each time. Not built this round (scope: this was a verification check, not a new table) —
named here as the next small, reversible, non-SLURM step if useful.

**Needs a decision first (scientific, Cindy's):**
1. Whether to accept the current provisional object (raw counts, no retention/decontamination
   applied) for a feasibility/method pilot, or require Gate 1 retention and/or DecontX correction
   to be applied first, specifically for monocyte work — independent of the broader QC rebuild's
   overall timeline.
2. Which of (A) full-index scFM discovery / (B) case-level description of the 2 transitioning
   donors / (C) a minimal paired DE design to pursue (per
   `mono_LONGITUDINAL_FEASIBILITY_NOTE_v1.md` §4) — unchanged by this round's QC findings.

**Not needed right now:** re-deciding repeated-visit selection for the independent-donor
cross-sectional view (Cindy has deferred it pending time evidence and the specific analysis
design — see `docs/execution-plan.md`'s 2026-10-16 third-follow-up entry), and no human web
lookup — this round's gaps were resolved from local/HPCC sources directly.
