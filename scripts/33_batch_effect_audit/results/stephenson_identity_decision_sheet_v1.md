# Stephenson (E-MTAB-10026) identity decision sheet — v1 (2026-10-10)

**SUPERSEDED by `stephenson_identity_decision_sheet_v2.md` (2026-10-11).**
Kept for the audit trail, not for use. A 2026-10-11 review
(`codex_output/20261011_113512.md`, F1) found this v1's §2 overclaimed in
five specific ways: "five independent sources"/"ground truth" language for
GEX (it is one more OBSERVED, processed label, not independent of BCR/TCR's
own demultiplexing lineage); it added two different cell-set intersection
sizes (21,481 + 74,761) and called the sum "96,242 distinct cells" without
checking their overlap; it implied TCR's silence on the disputed donors was
itself supporting evidence rather than a pure coverage gap; it asserted
"zero downstream consequence" for decision option 3 without actually testing
it; and it claimed only genotyping could resolve this further. v2 fixes all
five with real, computed evidence (including a new finding: GEX's own
`Resample`/`Collection_Day` fields show `CV0201` has a genuine D0/D28
two-visit pattern — `Initial`/`BGCV08` then `Resample`/`BGCV06` — which v1
did not have because it never read these GEX fields). The two real sections
below (§1 confirmed items, §4 Sanger caveat) are NOT known to be wrong and
are carried into v2 essentially unchanged; only §2's three "evidence"
claims, §3's `genotype-only` claim, and §5 are corrected.

Produced by `scripts/33_batch_effect_audit/code/audit_stephenson_source_identity.py`
(commit `0932cd2` and later). Read-only; no canonical file modified. Full evidence
in `claude-code_output/` and the timestamped `stephenson_identity_crosswalk_*` /
`job1_validation_addendum_v2_*` result directories this points to.

## 1. What is confirmed (no owner decision needed)

- 13 of 15 Cambridge `BGCV*` pools have identical donor membership across all five
  independent sources checked (SDRF, the local project table, the real deposited/
  processed GEX h5ad, BCR, TCR). No action needed on these.
- The Sanger-site `S9–S16`/`S23–S24` multi-AP-patient sample membership
  (`sanger_observed_membership.csv`) is independently reproduced from the TCR file.
  These are **observed Sample-level groupings**, not confirmed physical captures —
  see §4.
- Of the 12 Cambridge donors seen in more than one pool, 11 are a real two-visit
  `Initial`(D0)/`Resample`(D28) pattern; 1 (`CV0902`) is the SAME day/flag
  (`Initial`/D0) in two different pools — a different phenomenon, not a resample
  (`multi_pool_donor_classification.csv`).
- Zero BCR per-cell conflicts (one barcode, two donors) and zero TCR per-cell
  conflicts were found anywhere in either file.
- `run_batch_diagnostics.py`/`step2_run_diagnostics.sh` (Job 2) never reference
  `cv_group` — grepped directly. This identity question does not block Job 2 on
  that call path.
- Recomputing the existing proposed 21/22/4-donor Cambridge dependency components
  under BOTH possible resolutions of the one open conflict (§2) changes **zero**
  donors' component membership (`component_diff_sdrf_vs_local.csv`). The existing
  proposal (`stephenson_cambridge_proposed_dependency_groups_v1.csv`) does not need
  to wait on §2's resolution.

## 2. What is still a genuine, unresolved source conflict

**`BGCV06` / `BGCV13`, specifically donors `CV0201` and `CV0326`.**

| Source | `BGCV06` | `BGCV13` |
|---|---|---|
| SDRF (official ArrayExpress deposit) | `CV0037;CV0178;CV0234;CV0326` | `CV0073;CV0171;CV0201;CV0934` |
| Local project table (`stephenson_pool_composition.csv`) | `CV0037;CV0178;CV0201;CV0234` | `CV0073;CV0171;CV0326;CV0934` |
| Real GEX h5ad (`obs['sample_id']`, per-cell) | `CV0037;CV0178;CV0201;CV0234` | `CV0073;CV0171;CV0326;CV0934` |
| BCR (per-cell) | `CV0037;CV0178;CV0201;CV0234` | `CV0073;CV0171;CV0326;CV0934` |
| TCR (per-cell) | `CV0037;CV0178;CV0234` (no `CV0201` rows) | `CV0073;CV0171;CV0934` (no `CV0326` rows) |

Three sources derived independently of the SDRF text (the real processed GEX, BCR,
TCR) all agree with the local table and disagree with the SDRF specifically on
where `CV0201`/`CV0326` belong. A full-CellID join (never a bare barcode) between
the real GEX and BCR/TCR found **zero** per-cell disagreements across 96,242 shared
Cambridge cells (`gex_cell_level_crosswalk.json`) — every cell GEX calls `CV0201`
is also called `CV0201` everywhere else it is observed, and the same for `CV0326`.

This is strong evidence but **not dispositive**: the SDRF is the official
depositor-submitted record, and the possibility that the `BGCV13_CV0201` extract
named there was submitted/registered but genuinely excluded from, or never
represented a real cell in, the specific public processed object this project
built from has not been independently ruled out. TCR's silence on `BGCV06_CV0201`
and `BGCV13_CV0326` is a coverage gap (TCR libraries are not generated/retained for
every cell), not a vote for either side.

**What would resolve this with more certainty than is available now**: a response
from the depositors/original authors (Stephenson et al., EBI ArrayExpress curators)
confirming whether `BGCV13_CV0201`'s SDRF registration corresponds to a real,
separately-sequenced capture distinct from the processed object's own labeling, or
is itself a metadata entry error. This project does not contact external authors
without being asked to (per `AGENTS.md`).

**Decision needed from the owner**: pick one of —
1. **Trust the real processed GEX/BCR/TCR convergence** (treat `CV0201∈{BGCV06,
   BGCV08}`, `CV0326∈BGCV13`, matching the local table) — this is what every
   actual cell in this project's pipeline already reflects, since Job 1 and all
   prior Stage 32/33 work were built from this same GEX object.
2. **Trust the SDRF** (treat `CV0201∈{BGCV08,BGCV13}`, `CV0326∈BGCV06`) and accept
   that the real processed object may not fully match its own deposit's declared
   sample/pool assignment for these two donors specifically.
3. **Mark both donors' pool assignment `UNVERIFIED`** and exclude them from any
   pool-dependent grouping decision until external confirmation arrives, while
   still using their cells normally for anything that does not depend on pool
   identity (Job 1's per-donor sampling, for example, does not).

Consequence of each: since §1 already shows the proposed 21/22/4-donor components
are IDENTICAL either way, this decision currently has **no downstream effect** on
the Cambridge canonical-grouping proposal. Its only live consequence right now is
which donor->pool fact gets written into any future documentation or crosswalk
table that names `BGCV06`/`BGCV13`'s membership explicitly.

## 3. Continuing, non-blocking work

- Full genotype-level identity confirmation (demuxlet/souporcell-style) was not
  attempted and is out of this round's scope; it is the only evidence type that
  would move §2 from "strong convergent evidence" to "confirmed."
- `annotation_200112.csv` coverage: 128,083/130,644 (98.0%) of real Cambridge GEX
  cells have a matching CellID in the annotation file; 2,561 do not
  (`gex_cell_level_crosswalk.json`). Not investigated further this round — likely
  QC-filtering between the two releases, not itself evidence of an identity
  problem, but not confirmed as such either.

## 4. Sanger/Newcastle physical-unit caveat (unchanged from the prior round)

`sanger_observed_membership.csv`'s `Sample` grouping (e.g. `S9`/`S10` both showing
`AP3;AP8`) is an **observed TCR Sample-level co-occurrence**, not a confirmed
physical 10x capture. This project has not independently distinguished pool vs.
10x lane vs. sequencing library vs. analysis-sample for the Sanger/Newcastle sites
beyond what the IDF already states in prose. Left `UNKNOWN`, not asserted either
way.

## 5. What this decision sheet is not

Not a recommendation on which source to trust generally, not a statistical test of
identity, not a request to contact the depositors (that is a separate, explicit
owner decision if wanted), and not an authorization to run Job 2, change
`canonical_donor_labels.csv`, or promote the proposed Cambridge grouping. Those
remain separate owner decisions, now informed by the fact (§1) that this
conflict's resolution will not change the proposed component structure.
