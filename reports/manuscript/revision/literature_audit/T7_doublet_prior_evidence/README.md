# T7 — Doublet-prior evidence

**This is an evidence layer, not a configuration layer.** Nothing here selects a rate.
No Gate 2 runtime code may read any file in this directory. The boundary is enforced in
`docs/data-contracts.md`; only values that pass that contract may reach a prior configuration.

```
source audit      where each claim came from
capture evidence  mixed-level claims, each tagged with its own level
join key          identity mapping onto the 746 Gate 2 captures
verdicts          three orthogonal status fields per cohort
--------------- contract boundary; nothing above crosses it ---------------
prior config      only values that pass docs/data-contracts.md
sensitivity grid  rate grid, source, chemistry, model, limits
```

## Files

| File | Layer | Notes |
|---|---|---|
| `doublet_prior_source_audit.csv` | source audit | 23 sources, primary/secondary marked |
| `doublet_prior_capture_evidence.csv` | capture evidence | 605 rows, **mixed levels** — 603 `capture`, 2 `cohort-level` |
| `memo.md` | narrative | the survey this evidence was assembled from |
| `geo_capture_cell_counts.csv` | measured | repository-side barcode counts per deposit |
| `ren_gse158055_capture_qc.csv` | cohort evidence | Ren Table S1 per-sample QC ladder |
| `liu_demuxlet_per_lane.csv` | cohort evidence | **per-lane** SNG/DBL/AMB — see caveat below |
| `tenx_multiplet_rate_tables.csv` | vendor model | official 10x multiplet tables |
| `gate2_source_crosswalk.csv` | identity | 790 source artifacts → 746 captures |
| `capture_evidence_join_key.csv` | **derived** | built by `build_doublet_prior_evidence_layer.py` |
| `prior_evidence_verdicts.csv` | **derived** | built by the same script |

Regenerate the two derived tables with:

```bash
python3 scripts/32_rebuild_upfront_qc/code/build_doublet_prior_evidence_layer.py
```

## The join key is not optional

`capture_id` in the evidence table is **not a single key**. Verified per cohort:

| Cohorts | Resolves on | Coverage |
|---|---|---|
| Ren, COMBAT, Liu | `capture_id` | 409/409 |
| Qiu, Cillo, Chongqing | `repo_sample_id` | 171/171 |

Qiu's `capture_id` is the author's internal label (`HC1`), not the GEO accession
(`GSM5102900`). Joining on `capture_id` alone silently resolves 409 of 746 captures and
drops the remaining 337 without error.

## Coverage — 746 captures, three distinct kinds

| Kind | Captures | Meaning |
|---|---:|---|
| `mapped` | **580** | an evidence row is linked to this exact capture — **capture-specific evidence linkage** |
| `inherited` | **23** | Wilk: cohort protocol figure (15,000 cells/array, uniform) applied with `capture_specific=false` |
| `unmapped_to_gem_well` | **143** | Stephenson: deposit resolves to 98 sequence pools; GEM wells not identifiable |

**"Mapped" means linkage, not measurement.** A mapped evidence row may hold a directly measured
count, a published protocol value, a derived value, a proxy, or an upstream-processing record.
Identity mapping establishes *which* evidence belongs to a capture; it does not establish a
directly measured residual doublet prior for that capture, and no row in T7 supplies one.

The 23 and the 143 must additionally never be reported as capture-specific at all — a cohort
protocol statement broadcast across captures is still a cohort protocol statement.

## Two caveats that are easy to lose

**Ren.** The authors' `expected_doublet_rate = 0.08` was their Scrublet **input prior**, not a
measured removed fraction. Our Gate 2 input is demonstrably downstream of their filtering:
summing `Ncells_manual_removal` across the 284 samples gives 1,462,702, matching our
per-capture called-cell counts sample by sample. So 0.08 cannot be reused as our residual
prior, and the residual prevalence after their filtering is unknown.

**Liu.** The SNG/DBL/AMB figures are **per lane** and are kept that way in
`liu_demuxlet_per_lane.csv`. They are a genotype classification of a fixed 30,000-barcode
candidate whitelist that includes empty droplets — not a doublet rate among real cells. No
single cohort-wide number (such as the 36.3% from lane `B1_10xlane1`) may stand for Liu.

## Claude Science v3 — rate adjudication and scenario tables

Added 2026-09-21. These are **evidence and adjudication artifacts, not runtime configuration.**
They record which rate scenarios *could* be reasoned about and on what basis — not which rate
anything will use.

| File | Rows | Purpose |
|---|---:|---|
| `doublet-rate-adjudication-v3.md` | — | stage-alignment verdict, per-cohort reasoning, prospective manuscript template |
| `rate-sensitivity-scenarios-v3.csv` | 21 | candidate rate scenarios per stratum, each with its source stage, target stage and propagation assumption |
| `pooled-same-donor-pairing-scenarios-v3.csv` | 58 | modelled same-donor pairing probability Σpᵢ² per pooled capture |

### Runtime eligibility — nothing here is eligible

- **21 / 21** rows in `rate-sensitivity-scenarios-v3.csv` carry `runtime_eligible = false`.
- **58 / 58** rows in `pooled-same-donor-pairing-scenarios-v3.csv` carry `runtime_eligible = false`.
- **No Gate 2 code may read either table**, directly or through a wrapper. The same prohibition
  already covers every other file in T7.
- Turning any scenario into a runtime parameter requires a **separate, explicit, auditable
  decision** recorded in `docs/decisions.md`. It does not happen by a value being present in a
  table, by a column being readable, or by a scenario being labelled "middle".

### What the two CSVs are not

`Σpᵢ²` is a modelled probability that a randomly formed doublet joins two cells of the same
donor. It is **not** a measured residual doublet rate, **not** a demultiplexing efficiency, and
**not** evidence that some fraction of the original doublet burden survives. Its own
`propagation_assumption` column records that no relation has been established between it and the
doublet prevalence of the post-Gate-1 matrix.

The rate scenarios are likewise stage-mismatched by construction: vendor tables are indexed on
*recovered* cells while the available denominator is *called* cells, and each row states this
rather than papering over it.

### Manuscript wording

Section F of `doublet-rate-adjudication-v3.md` is a **prospective template**, marked as such in
the file. No grid has been run. Three passages were adjusted on intake so the text cannot be
mistaken for completed work:

1. the "Prospective manuscript template — do not use until the sensitivity analysis has been
   executed" warning was added before the block;
2. "Ambient and multiplet correction were applied to every capture using the same pipeline and
   the same parameter-selection rules" was replaced, because it overstated cross-cohort
   uniformity — cohort-specific upstream handling is real and must stay recorded;
3. "All results are reported across the full grid" and the clauses following it were put into
   the future tense.

A fourth, purely mechanical edit repointed the document's own companion-table references from
the received filenames to the names the files carry here. It also records that a third
companion named in the source, `v3_cohort_label_corrections.csv`, was not supplied.

### Source provenance

SHA-256 of the files as received:

| Received as | SHA-256 |
|---|---|
| `v3_doublet_rate_adjudication.md` | `a57f2d10681e1ab0775167801c5cf887a98055016f9e3791337df5142ce4a071` |
| `v3_rate_sensitivity_grid.csv` | `a70c6013f9481e65a61d5ab5ba4cc787d92d19df16fcb42679fe66c2c97d93c0` |
| `v3_pooled_same_donor_doublet_share.csv` | `1fa3b121c3e237a022aa64c387c376dcc11a630c2d2dfb4e75960bb9f1a4aff8` |

The two CSVs are committed byte-identical to the hashes above. The `.md` differs from its hash
by the three intake edits listed under "Manuscript wording"; nothing else in it was changed.

## Status fields

Three orthogonal fields, because one status cannot express a cohort whose loading is known
while its residual stays unknown. See `prior_evidence_verdicts.csv`.

`loading_evidence_status` · `upstream_doublet_handling` · `rate_model_status`

**No cohort is `DIRECT_RESIDUAL_ESTIMATE`.** The build script asserts this.
