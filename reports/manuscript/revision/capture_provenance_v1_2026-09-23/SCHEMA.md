# Capture provenance long table — schema v2

**Read-only audit. Nothing here modifies a builder, a DecontX result, a Gate 2 use class, or any
existing artifact. No job was submitted. No parameter or rate was selected.**

One row = **one membership of a biological unit in a physical capture**. Not one donor, not one
channel. A donor loaded into two captures is two rows; four donors sharing one capture are four
rows.

**A donor with no membership in any capture is NOT a row here.** It belongs to
`donor_coverage_reconciliation.csv`. v1 appended 15 such placeholders (Chongqing 1, Cillo 12,
Wilk 2) to the long table, where they inflated `n_membership_rows`, `n_biological_samples` and the
capture-status distribution with rows recording no membership at all — Cillo read as 83 rows / 82
samples when 12 of those rows were coverage bookkeeping.

## Why the keys are separate columns

`biological_sample_key` and `physical_capture_key` coincide in some cohorts and differ in others.
Recording them in one field makes the coinciding cohorts look like a general rule, and it silently
converts a pooled capture into N independent captures.

They answer different questions. Cross-donor doublets can only form **within one physical
capture**; ambient RNA is a property **of one physical capture**. Outcome, disease and event
belong to the **biological** unit. A pipeline batched on the wrong one is not wrong by a constant
factor — it estimates a per-capture quantity from a fraction of the capture.

## Columns

| column | meaning |
|---|---|
| `cohort` | cohort label as used in this repository |
| `accession` | **the actual deposit.** Required: the cohort label is NOT a join key (see Wilk) |
| `donor_key` | the subject |
| `biological_sample_key` | the deposited sample/library record for that subject |
| `draw_key` | the timepoint, where the deposit resolves one |
| `physical_capture_key` | the GEM well / array / capture the cells were generated in |
| `current_analysis_channel_id` | what the builder currently uses as `channel_id` |
| `membership_kind` | `intended` / `observed_singlet` / `inferred` |
| `physical_capture_status` | `resolved` / `provisional` / `unresolved` |
| `analysis_matches_physical` | `true` / `false` / `unknown` |
| `evidence_source` | the artifact the row is asserted from |
| `evidence_locator` | file **plus the column/row rule**, never a filename alone |
| `evidence_grade` | `measured_metadata` / `published` / `derived` / `inferred` |
| `ambiguity_reason` | why the status is not `resolved`; blank only when it is |

`accession` is the one addition beyond the reviewed field list. It is load-bearing: "Wilk" names
two deposits in this repository (`GSE174072` in the benchmark, `GSE150728` in the T3 mapping
audit) with non-intersecting donor namespaces, so a future join on the cohort label alone would
silently mix two studies.

## Value definitions, and what each refuses to claim

### `membership_kind`

- `intended` — the deposit states which biological units went into this capture.
- `observed_singlet` — the unit was **recovered** as a singlet by genotype demultiplexing.
- `observed_assigned` — the unit appears in the deposit's own per-cell donor assignment among a
  capture's **retained** cells, after the depositors' QC. COMBAT is this: no pool-design artifact in
  this repository lists the intended donors per channel (`combat_channel_inventory.csv` carries only
  an `n_donors` count), so `intended` would have overclaimed.
- `inferred` — derived from structure rather than stated.

For **both** observed kinds, the absence of a row is **not** evidence of absence from the capture:
a unit loaded but not recovered, or filtered out by the depositors' QC, produces no row. Never
report either as donor proportions or pool composition. An invariant requires every
`observed_*` row to carry a stated limitation in `ambiguity_reason`.

`unresolved` is **not** a membership kind. A relation that cannot be determined is not a
membership, and an invariant refuses any such row in this table.

### `physical_capture_status`

- `resolved` — **the deposit identifies the capture unit directly.**
- `provisional` — a capture-level key exists, but its identity as one physical capture rests on
  something the deposit does not state: this repository's own annotation, or a key the builder
  derived.
- `unresolved` — no key identifies the physical capture.

Two things that do **not** earn `resolved`:

1. **A repository gloss.** Liu's capture key is the lane token in the deposit's GSM titles
   (`B1_10xlane1`). Reading that token as one 10x chip well is this repository's annotation, not a
   statement in the deposit, so all 877 Liu rows are `provisional`. v1 marked them `resolved` while
   the findings text simultaneously admitted the gloss — the status contradicted the prose.
2. **A builder-derived key.** Wilk's `channel_id` comes from the build summary rather than from the
   deposit, so it is `provisional`.

`resolved` is also never assigned because a sample *lacks* a pooling marker. Absence of a pooling
marker is not evidence of a dedicated capture.

### `analysis_matches_physical`

`false` means the builder's `channel_id` is a different granularity from the physical capture —
the finding, not an error in the row. `unknown` means the physical capture is not `resolved`, so
the comparison cannot be made.

## Outcome boundary — strict, at parse time

Four inputs carry outcome columns:

| input | outcome columns present |
|---|---|
| `stephenson_patient_outcomes.csv` | `Outcome`, `outcome_unique_within_patient` |
| `canonical_donor_labels.csv` | `event`, `event_definition`, `disease` |
| `ren_gse158055_channel_inventory.csv` | `outcome` |
| `qiu_gse167363_channel_inventory.csv` | `disease_state` |

v1 read every input with `csv.DictReader`, which materialises the whole row, and then claimed "no
outcome column was read". That was false: the values were parsed and held, merely never indexed.
Not referencing `r['Outcome']` is not the same as not reading it.

**The exact claim, and nothing stronger: no outcome column is parsed into a field, inspected,
propagated, exported or analysed.** Checksumming opens each input with `open(path, 'rb').read()`,
so the raw bytes of the whole file — an outcome field's bytes included — do pass through a hash.
Hashing bytes and reading a value are different acts, and conflating them would be the same kind of
overclaim as v1's.

Every input is now read through `read_columns()`, which uses `pandas.read_csv(usecols=...)` so the
other columns are **excluded at parse time**, refuses any request for a column matching the
outcome vocabulary, and records per input the columns read and the columns excluded. That record
is `outcome_column_boundary.csv` — the boundary is an artifact, not a sentence in a docstring.

The guard pattern is deliberately broader than the four known names, so a new input with a
differently spelled label column is refused rather than silently read.

The parser argument is verified by capturing the kwargs `read_columns()` actually passes, not by
inspecting the returned columns: the data read is followed by `df[list(usecols)]`, so dropping
`usecols=` still yields only the requested keys while pandas has parsed every column. A v2 test
checked the returned keys and was reported as killing that mutation — it did not, because that
mutation had also removed the projection, so the kill never isolated the parser argument.

Files **cited** as evidence but never parsed are checksummed too (the Stephenson rows rest on the
T7 memo's SDRF reading). An input whose checksum is absent can change without the bundle noticing.

## Invariants asserted by the build script

Only enumerable facts. Each is a count or key relation recomputable from the inputs:

1. every row has a non-empty `cohort`, `accession`, `donor_key`, `biological_sample_key` and
   `current_analysis_channel_id`;
2. no row carries `membership_kind='unresolved'` — coverage gaps are not memberships;
3. `physical_capture_status='resolved'` requires a non-empty `physical_capture_key`;
4. `ambiguity_reason` is non-empty whenever the status is not `resolved`;
5. `evidence_locator` names a column or row rule, not only a file;
6. `membership_kind='observed_singlet'` never carries `evidence_grade='measured_metadata'`;
7. an `unresolved` capture is never compared against the analysis channel;
8. no duplicated membership key within a cohort;
9. per-cohort donor / sample / row counts equal the counts in the cited source.

Summary counts exclude empty keys. v1 counted the empty `current_analysis_channel_id` of a
placeholder row as one channel, inflating Chongqing 133→134, Cillo 44→45 and Wilk 23→24 — and
because the summary was rebuilt with the same logic, "the summary is rebuildable" held while being
semantically wrong.

## What is NOT asserted

Interpretations, stated in `FINDINGS.md` and deliberately not encoded:

- whether a `provisional` capture is in fact one physical capture;
- whether any existing DecontX or doublet result should be recomputed;
- whether an unresolved cohort can be made runnable;
- how the `docs/data-contracts.md` contradiction should be resolved.
