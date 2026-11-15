# Capture provenance findings v1 — 2026-09-23

**Read-only.** Nothing was modified, recomputed, deleted or submitted. No builder, DecontX
output, Gate 2 use class or living document was touched. No parameter or rate was selected. No
outcome column was read.

Build: `extract_hpcc_membership.py` (raw HPCC → two extracts, `RAW_INPUT_MANIFEST.txt`) then
`build_capture_provenance.py` → `capture_provenance_long.csv` (**2,522 membership rows**, 8
cohorts), `cohort_summary.csv`, `donor_coverage_reconciliation.csv`,
`outcome_column_boundary.csv`, `INPUT_CHECKSUMS.txt`. All enumerated invariants held;
`test_capture_provenance.py` holds 15 negative tests.

**v2 corrections to this document and the table** (details in `AUDIT.md`): the 15 coverage-only
placeholder rows are gone from the membership table, Liu is `provisional` rather than `resolved`,
empty keys no longer count as a category, and every input is read with `usecols` so the outcome
boundary is real rather than asserted.

| cohort | donors | bio samples | physical captures | analysis channels | rows | capture status |
|---|---:|---:|---:|---:|---:|---|
| COMBAT 2022 | 124 | 124 | 70 | 70 | 979 | resolved 979 |
| Chongqing 2025 | 133 | 133 | 133 | 133 | 133 | resolved 133 |
| Cillo 2021 | 29 | 70 | 39 | 44 | 71 | resolved 66, unresolved 5 |
| Liu 2021 | 50 | 50 | 55 | 55 | 877 | **provisional 877** |
| Qiu 2021 | 7 | 12 | 12 | 12 | 12 | resolved 12 |
| Ren 2021 | 196 | 284 | 284 | 284 | 284 | resolved 284 |
| Stephenson 2021 | 130 | 143 | **0** | 143 | 143 | **unresolved 143** |
| Wilk 2021 JEM | 23 | 23 | 23 | 23 | 23 | provisional 23 |

15 canonical donors have no membership in any capture (Chongqing 1, Cillo 12, Wilk 2). They are
recorded in `donor_coverage_reconciliation.csv` **only** — they are coverage facts, not
memberships, and an invariant now refuses them in the long table.

---

## F1 — Stephenson: the analysis channel is a biological sample, not a physical capture

**Status: the physical capture is UNRESOLVED for all 143 rows.**

`build_stephenson_decontx_input.py:21` calls `obs/sample_id` "the GEM-well/library unit".
Measured against the deposit's own metadata, it is not.

The ArrayExpress SDRF resolves the 143 samples to **98 `sequence pool` values — 83 pools of one
sample and 15 pools of four** (`T7_doublet_prior_evidence/memo.md`, Stephenson section). This is
independently consistent with `stephenson_patient_outcomes.csv`: 47 patients / 60 samples carry a
`BGCV<nn>` pool prefix across 15 pools, and 83 patients / 83 samples carry none.

But a sequence pool is still **not** a GEM well. The same source states that the Cambridge pools
**each span two lanes, and which cells came from which lane is not recorded publicly**. So:

- the 60 pooled rows sit in 15 pools whose cells were spread over two sequencing lanes
  each, with the lane assignment unpublished. **This does not make them ~30 GEM wells.** A
  sequencing lane is a flow-cell partition; a GEM well is a droplet-generation channel. One
  library can be sequenced across several lanes, and the deposit says nothing about how many
  droplet-generation channels produced each pool. The pool-to-well relation is simply absent;
- the 83 unpooled rows have no pooling marker, and **absence of a marker is not evidence of a
  dedicated well** — the deposit never states one sample = one GEM well.

Every Stephenson row is therefore `unresolved`, not the `83 + 15 = 98` split an earlier note of
mine proposed, and not any well count derived from lane counts. That arithmetic happened to equal the SDRF's pool count, which made it look
confirmed; the SDRF is the evidence, and it stops one level above the GEM well.

**Corollary about circular evidence.** `gate2_source_crosswalk.csv` marks all 143 Stephenson rows
`matched`. That join's producer of record *is* the builder, so `matched` shows only that
`obs/sample_id` joins to the builder's own 143 channels. It cannot independently establish
physical-capture granularity, and it is not reused here as a capture status.

**Repository contradiction, left in place.** `docs/data-contracts.md:296` lists Stephenson's
capture unit as "GEM well / library = `obs/sample_id`", while `:359` and `:404` of the same file
record 143 captures as `unmapped_to_gem_well` with 98 sequence pools. Both sentences are in the
repo today. This table separates the two key levels rather than picking one; **the document is
not edited here.**

### Structural risk, registered only

Cross-donor doublets can only form within one physical capture, and ambient RNA is a property of
one physical capture. Stephenson's DecontX run batched on `channel_id = sample_id`. For the 60
pooled samples that is four batches per sequence pool, so whatever physical capture produced a
pool had its ambient profile estimated from a subset of itself four times over.

How large that effect is cannot be stated from here, because the pool-to-capture relation is
exactly what is missing. The risk is that the batching unit is known to be finer than the pool,
not that a specific number of wells was mis-split.

This is **registered, not acted on**. No existing result is declared invalid, nothing is
recomputed, and nothing is overwritten. Whether Stephenson needs rebuilding at a different
capture granularity is a decision for after this mapping evidence is itself audited.

**No impact on Gate 2 today:** Stephenson is `not_run` in the frozen design, so this neither
unblocks nor blocks the eleven-parameter review.

## F2 — "Wilk" names two different deposits

The benchmark cohort is **GSE174072** ("Wilk 2021 JEM", Seq-Well, donors `28205-0562`…).
The T3 mapping audit's Wilk rows are **GSE150728** ("Wilk 2020", 13 samples / 12 patients,
patients `HIP002`…). The donor id namespaces do not intersect.

Joining these two sources on the cohort label would silently mix studies. The table carries
`accession` on every row for this reason, and the two deposits are never joined.

Wilk's capture is a **Seq-Well array**, not a 10x GEM well, so the 10x multiplet tables are
inapplicable by construction and all 23 built rows are `provisional`: the deposit does not state
that one `channel_id` is exactly one array.

## F3 — Cillo: an empirically mapped Seurat suffix, and a source that does not cover the cohort

The h5ad join key is `<barcode without -1>_<h5ad_suffix>`, where the suffix is a nested Seurat
`merge()` double suffix (`"1_1"`, `"2_1"`). **It is not algorithmically derivable** — it is mapped
per GSM in `cillo_gse180578_channel_inventory.csv`, with `join_status`: 40 exact, 6
`PARTIAL — h5ad grouping differs from GEO library`, 2 `no h5ad cells`. Those 8 libraries stay
`unresolved`; their cells cannot be attributed to one GEM well through this key.

**40 `exact member-set match` libraries, 39 resolved captures.** Not a contradiction:
`GSM5465163` is an exact match of an *empty* member set (`n_donors_in_lib=0`, `members` blank), so
it has a suffix but contributes no membership row and therefore no capture.

Separately, GEO's `members` column names only 29 donors, while `canonical_donor_labels.csv` lists
40. The 12 not named are the 10 `Healthy_*` donors and the 2 `ALIR_*` samples; one donor
(`ARDS_5`) appears in GEO but not canonical. The 12 are kept as `unresolved` rows so coverage is
not overstated by dropping them.

## F4 — Liu: observed singlet recovery, not intended pool composition

877 membership rows over 55 lanes, from `DROPLET.TYPE == "SNG"` / `SNG.BEST.GUESS`. Every row is
`membership_kind=observed_singlet`, `evidence_grade=derived` and
**`physical_capture_status=provisional`**.

v1 marked these `resolved` while this same section admitted the lane-token reading is a
repository gloss — the status contradicted the prose. `SCHEMA.md` defines `resolved` as *the
deposit identifies the capture unit directly*, which a gloss is not.

This is **post-demultiplex recovery**. A donor loaded into a lane but yielding no singlet call
produces no row, so the absence of a row is not evidence of absence from the pool, and these
counts must not be reported as donor proportions or as pool composition.

Measured: 14–19 donors per lane, 7–24 lanes per donor, 50 distinct genotype ids. After stripping
the `_B<n>`/`_T<n>` suffix, the **33 HGR-form base ids match the 33 canonical Liu donors exactly**;
the remaining 17 are non-benchmark panel genotypes (`SHD1–8`, `AA*`, `CHI*`, `HDML`, `HDVO`, `IM*`).

**56 declared, 55 present.** `canonical_donor_labels.csv`'s `channel_unit` for Liu reads
"Batch x 10x Lane (56 wells; 14-19 donors/well)", but the deposit carries 55 lanes: B1 has 7
(lane 6 is absent), B2 and B3 have 24 each. The table holds the 55 that exist. The missing lane
is recorded, not reconciled -- nothing here establishes whether it was never deposited, never
run, or merged into another.

Liu's capture key is the lane token in the deposit's own GSM titles (`B1_10xlane1`), which is
also the unit demuxlet was run per. The gloss of that lane as a 10x chip well is this
repository's annotation, not a statement in the deposit.

Draw time stays `unresolved_mixed_pool`: GEO's own timepoint characteristic for these pooled
lanes is literally "mixed", and genotype demultiplexing cannot separate a donor's repeat draws.

## F4b — COMBAT: observed assignment, not intended loading

979 membership rows over 70 channels and 124 donors, from the deposited h5ad's own
`(obs/Channel_ID, obs/COMBAT_ID)` pairs across 836,148 cells. 13–14 donors per channel, 6–21
channels per donor. Independently equal to `sum(n_donors)=979` in `combat_channel_inventory.csv`,
which was produced separately.

`membership_kind='observed_assigned'`, **not** `intended`. These are donors observed among a
channel's **retained** cells after the depositors' QC. No pool-design artifact in this repository
lists the intended donors per channel — `combat_channel_inventory.csv` carries only an `n_donors`
count — so a donor absent from a channel here cannot be called unloaded. v2 marked these `intended`,
which claimed loading design from post-QC co-occurrence; the same limitation Liu's rows already
carried.

The capture itself stays `resolved`: `Channel_ID` is deposited metadata, and the T7 memo records
published per-channel loading ("COMBAT at 50,000 cells per channel", loading verdict VALIDATED), so
the channel is supported as the capture unit by more than a builder comment.

## F5 — Ren: 284 libraries, 196 donors

`sampleID` ↔ `gsm` is 1:1 across all 284 samples; `donor` (`P-*`) is 1:many over them. Confirmed
against `mapping_verdicts.csv` (patient key = the `patient` characteristic; title base `S-M003`
maps 1:1 onto `P-M003` for all 284).

`canonical_donor_labels.csv`'s `channel_unit` reads "GSM sample = one library" for Ren. That
describes the library side only; read alone it invites the conclusion that donor and capture are
1:1, which they are not.

## F6 — Donor coverage gaps, none dropped

| cohort | canonical | in table | normalisation applied | canonical without membership |
|---|---:|---:|---|---|
| Chongqing | 133 | 133 | remove spaces, upper-case (`pHC 01` → `PHC01`) | 1 (`Skin-S152`) |
| Wilk | 25 | 23 | none | 2 (`55650-0132d0`, `55650-0132d7`) |
| Cillo | 40 | 29 | none | 12 (10 `Healthy_*`, 2 `ALIR_*`) |
| Liu | 33 | 50 | strip `_B<n>`/`_T<n>` | 0 |
| Qiu / Ren / Stephenson / COMBAT | — | — | none | 0 |

Every uncovered canonical donor is recorded in `donor_coverage_reconciliation.csv`, and in that
table only. v1 also appended them to the membership table as `unresolved` rows, which inflated
Cillo to 83 rows / 82 samples and added a phantom analysis channel to three cohorts. A coverage
gap is not a membership.

Each normalisation rule is stated in the reconciliation table rather than applied silently.

---

## Not asserted

These are interpretations, deliberately outside the script's invariants:

- whether a `provisional` capture is in fact one physical capture;
- whether any existing DecontX or doublet result should be recomputed;
- whether Stephenson can be made runnable at 98, ~30, or any other granularity;
- whether the `docs/data-contracts.md` contradiction should be resolved toward either sentence.
