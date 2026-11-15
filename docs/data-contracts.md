# Data Contracts

Explicit schemas for the major files this project depends on. Data contracts
prevent silent downstream bugs. Large binaries (`*.h5ad`, `*.pkl`, …) are
gitignored and live on HPCC under `project_root/data`; the contracts below
describe what those files must contain.

> Fill in the TODOs from the canonical HPCC files. Each contract should be
> validatable by `scripts/99_audit/code/check_data_contracts.py`.

---

## `data/processed/covid_sepsis_mortality.h5ad`  (primary mortality benchmark input)

### Required `obs` columns
| Column | Type | Description |
|---|---|---|
| donor_id | str | patient/donor identifier — CV grouping key (no donor split across folds) |
| cohort | str | source cohort (covid / sepsis) |
| disease | str | disease label |
| mortality | int | target: 1 = deceased, 0 = survived (severe **excluded** from y=0, see memory 302006 severe-excluded) |

### Required `var` columns
| Column | Type | Description |
|---|---|---|
| gene_symbol | str | gene name |

### Required layers / obsm
| Key | Description |
|---|---|
| `layers["counts"]` | raw integer counts |
| `obsm["X_<fm>"]` | per-FM embedding (geneformer, scgpt, uce, state, transcriptformer) |

### Invariants
- `.X` is log-normalized; `layers["counts"]` is raw integer counts.
- `obs_names` and `var_names` are unique.
- No donor appears in more than one CV fold.

### Provenance
- Original source: TODO (GEO accessions for COVID + sepsis cohorts)
- Download/access date: TODO
- Processing stage: `scripts/24_option_c_canonical_input/`
- Processing script: TODO `code/run_24_*.py`
- Manifest: `scripts/24_option_c_canonical_input/results/manifest.json` (to be generated)

### Integrity checks
- Expected cells: ~42,095 ; donors: ~322 (per `scripts/README.md` project state)
- Missing value policy: TODO
- Duplicates policy: unique obs_names enforced
- Known exclusions: severe cases excluded from y=0 (Job 23761950)

---

## `obsm["X_<fm>"]` embedding arrays

### Invariants
- One row per cell, aligned to `obs_names` order.
- No NaN/Inf (validated by `scripts/00_shared/code/embedding_sanity.py`).
- Dimensionality fixed per FM and recorded in the stage manifest.

---

## TODO — additional contracts to document
- [ ] Per-FM raw embedding files produced by `scripts/foundation_models/`
- [ ] Cross-disease train/test subsets used by `scripts/06_batch_generalization/`

---

## `gate2_input/<cohort>/prior_counts.csv`  (RETIRED — no longer produced or consumed)

> **Retired 2026-09-21.** `materialize_liu_wilk_for_gate2.py` no longer writes this file,
> and nothing reads it. It recorded cells *entering Gate 1* under the undefined name
> `n_pre_qc`, which was then fed to a loading-rate formula defined over droplets recovered
> from the instrument — two different levels of the count ladder. That ambiguity is what
> let the mismatch pass through six cohorts of invalid scoring. Existing copies under
> `gate2_input/{liu,wilk}/` are retained only as inputs to already-quarantined runs.
> The replacement is the runtime configuration below. The schema is kept as the record of
> what a valid prior would have had to contain.

**Status: no cohort currently yields a usable residual prior — but the reasons differ, and a
single verdict was too coarse.** Specified 2026-09-21 after an audit found the existing priors
count the wrong population. An earlier revision of this section said "UNSATISFIABLE — all eight
cohorts", which conflated *no instrument-level recovered count* with *no evidence at all*. Four
cohorts do have nominal loading evidence. Graded per cohort in
`reports/manuscript/revision/literature_audit/T7_doublet_prior_evidence/prior_evidence_verdicts.csv`.

**Contract boundary — enforced, not advisory.** The T7 evidence layer must never be read by
Gate 2 runtime code. `doublet_prior_capture_evidence.csv` mixes capture-level and cohort-level
claims; `capture_evidence_join_key.csv` is an identity mapping, not a prior configuration.
Only values that satisfy the columns and constraints below may reach a prior configuration.

The loading-rate formula `loading = min(0.008 × n_recovered / 1000, cap)` is defined over
**droplets recovered from the instrument**. Existing files were generated from
`flags.groupby('channel_id').size()` — cells *entering Gate 1*, i.e. an already-processed
candidate population — which is a different quantity. `n_pre_qc` was used as an undefined
catch-all name; that ambiguity is what let the mismatch pass unnoticed.

### Required columns

| Column | Type | Description |
|---|---|---|
| `channel_id` | str | capture unit, exactly as it appears in `cells.csv`; unique |
| `n_candidate_barcodes` | int | barcodes considered before cell calling |
| `n_called_pre_qc` | int | called cells before any QC |
| `n_post_demux_singlets` | int | after demultiplexing removed cross-donor doublets |
| `n_gate1_retained` | int | after frozen basic QC |
| `n_donors` | int | donors contributing to this capture |
| `donor_proportion_stage` | str | the stage at which donor shares were taken |
| `loading_prior_column` | str | **which column above feeds the loading rate** — required, no default |
| `loading_rate_source` | str | vendor figure or record the rate comes from |
| `is_proxy` | bool | true when the loading rate is borrowed from another platform |
| `source_artifact` | str | path the counts were derived from |
| `source_checksum` | str | checksum of that artifact |

### Constraints (all fail-loud)

- every `channel_id` in `cells.csv` present exactly once;
- `n_candidate_barcodes >= n_called_pre_qc >= n_post_demux_singlets >= n_gate1_retained`;
- the column named by `loading_prior_column` is a positive integer and **strictly greater** than
  the capture's scored cell count;
- `n_donors >= 1`; donor proportions sum to 1 within tolerance;
- `is_proxy` true requires a non-empty note naming the borrowed platform (Wilk / Seq-Well).

A missing file, a missing channel, a duplicate channel, or any violated constraint stops that
cohort. Falling back to a post-QC cell count is prohibited — it is the defect this contract
exists to prevent.

### Runtime configuration — the only route a rate reaches Scrublet

Implemented 2026-09-21 as `scripts/32_rebuild_upfront_qc/code/gate2_runtime_config.py`. This is a
**different artifact from `prior_counts.csv` above and from every T7 table**: the contract above
describes counts, T7 describes what could be reasoned about, this records what a human
explicitly authorized. `schema_version` = `gate2-runtime-config/1`.

**No authorized configuration exists.** `gate2_run.py` therefore exits 2 without one. That is the
intended state while rate selection is unresolved.

| Column | Must be |
|---|---|
| `schema_version` | exactly `gate2-runtime-config/1` |
| `cohort`, `channel_id` | match the run; every capture present exactly once, no extras |
| `platform` | `10x` or `seqwell`; unknown platforms raise rather than borrow a rate model |
| `gate2_input_stage` | exactly `post_gate1_scrublet_input` |
| `expected_doublet_rate` | finite, strictly inside `(0, 1)` |
| `rate_estimand` | must name `post_gate1_scrublet_input` — a rate estimated for one stage may not be applied to another |
| `rate_source`, `decision_reference` | non-empty; the decision reference points at `docs/decisions.md` |
| `runtime_eligible` | exactly `true`; blank, `1`, `yes` all raise |
| `value_basis` | `measured` / `modelled` / `extrapolated` / `arbitrary` — so a modelled number can never be read back as a measurement |
| `n_gate1_retained`, `n_donors` | positive integers |
| `n_candidate_barcodes`, `n_called_pre_qc`, `n_post_demux_singlets` | positive integers when known, else blank |
| `donor_proportion_stage`, `donor_proportions_checksum` | non-empty; the checksum binds the mixture into the checkpoint |
| `source_artifact`, `source_checksum`, `code_provenance` | non-empty |

Fail-closed, with no degraded mode: missing argument, missing file, missing channel, duplicate
channel, unknown platform, non-finite rate, rate outside `(0,1)`, stage/estimand mismatch, and
missing decision reference all raise `Gate2ConfigError` before any scoring starts.

**T7 paths are refused by name.** Any path containing `T7_doublet_prior_evidence`,
`literature_audit`, `rate-sensitivity-scenarios`, `pooled-same-donor-pairing`,
`doublet_prior_capture_evidence`, `prior_evidence_verdicts` or `capture_evidence_join_key` is
rejected as a runtime config, so the evidence layer cannot be wired in by accident or by a
wrapper.

Both historical fallbacks are removed, not repaired: the optional `--prior-counts` guarded by
`os.path.exists`, and `n_cells_recovered if ... is not None else n_cells`. The scripts that
contained them now exit with a refusal. Old `n_pre_qc` values are not reinterpreted.

### Gate 2 score-distribution diagnostics

Computed while both score arrays are still in memory — Scrublet's simulated scores live only
inside the call and cannot be recovered afterwards. Persisted in
`<cohort>_gate2_diagnostics.csv` and inside the checkpoint, carrying their own
`diagnostics_schema_version`.

| Field | Meaning |
|---|---|
| `n_observed` / `n_simulated` | sample size of each distribution |
| `observed_score_min/max`, `simulated_score_min/max` | observed range, **not clipped to [0,1]** |
| `observed_simulated_auc` | separability, **simulated is the POSITIVE class**, defined as `P(sim > obs) + 0.5·P(sim == obs)` |
| `simulated_recall_at_threshold` | `mean(sim >= threshold)` |
| `observed_tail_fraction_at_threshold` | `mean(obs >= threshold)` |
| `simulated_tail_fraction_at_threshold` | same as recall, stored separately so the two definitions cannot be confused |
| `*_percentiles_p0_to_p100` | 101 points, p0..p100 inclusive |
| `shared_histogram_bin_edges` + `observed/simulated_histogram_counts` | 50 bins over the **joint** finite range, so the two are directly comparable |
| `histogram_range_degenerate` | true when a distribution is constant and the range had to be widened |

**These are diagnostic evidence, not a reliability criterion.** Specifically:

- a high AUC does **not** authorize deletion, and nothing in the code turns it into a decision;
- `simulated_recall_at_threshold` is recall against **synthetic** doublets built by summing
  observed cells — it is **not** the true doublet recall, which is unmeasurable here;
- percentiles and histograms **summarize** the distributions, they do not reconstruct them; a
  claim needing the raw distribution needs the diagnostic-only dump;
- threshold-based fields are **null**, never `0`, when the threshold is not finite — `0` would
  read as a measured absence of tail mass rather than an unavailable measurement;
- distribution-shape fields are still computed for a `threshold_failed` capture, where they are
  the main evidence available;
- **no reliability threshold is set and `method_unreliable` is never assigned.**

Full simulated score vectors are **not** saved by default. `score_one_capture(...,
dump_simulated_scores_to=path)` is a diagnostic-only opt-in for validating these summaries
against the raw distribution in a future pilot; production runs must leave it off.

### Gate 2 call output — `<cohort>_gate2_calls.csv`

| Column | Meaning |
|---|---|
| `call_status` | **the authoritative state.** One of `predicted_doublet`, `scrublet_noncandidate`, `too_small`, `threshold_failed`, `method_unreliable`, `runtime_failed` |
| `candidate_doublet` | nullable integer: `1` doublet, `0` not-a-doublet, **blank = UNKNOWN** |
| `scrublet_score` | retained even when the capture failed, as evidence |
| `gate2_cell_key` | **DERIVED convenience column, not the identity.** Canonical JSON array `["cohort","channel_id","cell_id"]` — injective for any content including `|`, commas, quotes, whitespace and non-ASCII. An unescaped `a|b|c` join is *not* injective (`('A|B','c1')` and `('A','B|c1')` collide) |

**The authoritative identity is the three columns `(cohort, channel_id, cell_id)`.**
`cell_id` is unique only within a capture — barcodes legitimately repeat across captures — so
joining on `cell_id` alone can silently merge different cells. Join on the tuple.
`gate2_cell_key` is a derived convenience column; nothing, least of all deletion eligibility,
may depend on it.

**Two senses of "singlet", kept apart.** They are orthogonal and neither is an authorization:

| Term | Axis | Means |
|---|---|---|
| demuxlet / genotype **SNG** | donor assignment | the barcode carries ONE donor's genotype — says nothing about how many cells |
| `scrublet_noncandidate` | expression-based doublet detection | the capture scored successfully, the automatic threshold was valid, and this cell did not exceed it |

A Liu barcode can be a demuxlet singlet **and** a same-donor doublet simultaneously.
`scrublet_noncandidate` does **not** assert a true singlet, confers no retention authorization,
and does not exclude a homotypic doublet — Scrublet cannot detect those by construction. The
former name `auto_singlet` was withdrawn precisely because it collided with the demuxlet sense.

**Two layers, and only the first exists in this repository:**

```
mechanical candidate         call_status == "predicted_doublet"
final deletion authorization mechanical candidate
                             + cohort eligible for this scoring at all
                             + an approved runtime rate
                             + capture passed a reliability criterion
                             + a sensitivity analysis supports the call
```

`predicted_doublet` is **necessary but not sufficient** for removing a cell. **No
final-authorization artifact exists**, and none should be created — not even an all-false
per-cell column, which would later read as though a formal adjudication had already happened.
The accurate representation of the current state is that the artifact is absent.

A keep mask is not evidence of singlet status.

Four things this layer does NOT establish: an expected rate is not a deletion quota; top-`r`
selection is not an authorization rule; a successful automatic threshold does not mean the
capture is reliable; and passing a reliability criterion would still not delete anything on its
own. Outcome/mortality information plays no part in any QC rule or authorization decision.

**`candidate_doublet.fillna(False)` — and any equivalent — is prohibited.** Blank means the cell
was never successfully examined (`too_small`, `threshold_failed`, `runtime_failed`, or a future
`method_unreliable`). Coercing it to `False`/`0` asserts "examined and not a doublet", which is
the exact defect this column's tri-state exists to prevent. Read `call_status`; if you need a
boolean, derive it with `gate2_status.mechanical_candidate_mask` — which marks
mechanical candidates, not cells cleared for removal.

Order sensitivity: the layer guarantees *identity* (results join by `cell_id`, and a reordered
input refuses a stored checkpoint), **not numerical invariance**. Scrublet's neighbourhood
construction is position-dependent, and nothing here establishes that scores are invariant to
cell row order. Do not assume it.

### Provenance crosswalk (the authoritative record)

`05_RESULTS/rebuild_upfront_qc/gate2_prior_inventory/gate2_source_crosswalk.csv`, produced by
`scripts/32_rebuild_upfront_qc/code/gate2_source_crosswalk.py` (read-only). One row per source
artifact plus one per unmatched Gate 2 capture; 790 rows covering 746 Gate 2 captures.

Reconciliation — **every Gate 2 capture has a located source artifact** (`gate2_only` = 0):

| Cohort | matched | source_only | Gate 2 captures |
|---|---:|---:|---:|
| chongqing | 124 | 9 | 124 |
| cillo | 35 | 13 | 35 |
| combat | 70 | 0 | 70 |
| liu | 55 | 0 | 55 |
| qiu | 12 | 0 | 12 |
| ren | 284 | 0 | 284 |
| stephenson | 143 | 0 | 143 |
| wilk | 23 | 22 | 23 |
| **total** | **746** | **44** | **746** |

`source_only` rows are samples present in the deposit but not carried into the benchmark; they
are kept in the crosswalk rather than filtered out, so the difference is auditable.

### channel_id mapping — from each builder, which is the producer of record

| Cohort | `channel_id` is | obs field / source | Producer |
|---|---|---|---|
| qiu, chongqing | bare GSM accession | `barcodes.tsv.gz` filename | `build_{qiu,chongqing}_decontx_input.py` |
| cillo | GSM accession | GEX tar member | `build_cillo_decontx_input.py:220` |
| stephenson | GEM well / library | `obs/sample_id` | `build_stephenson_decontx_input.py:84` |
| combat | GEM well (not `Pool_ID`) | `obs/Channel_ID` | `build_combat_decontx_input.py:92` |
| wilk | converted dir basename | — | `convert_wilk_rds_to_mtx.R` → `build_wilk_decontx_input.py:131` |
| liu | lane id | h5 filename | `build_liu_decontx_input.py` |
| ren | sampleID | `cell_annotation.sampleID` | `build_ren_decontx_input.py` |

Three encoding traps that silently produced wrong or empty joins, fixed and recorded so they are
not reintroduced:

- **Stephenson** `obs/sample_id` is an `int16` dataset whose labels live in the legacy
  `obs/__categories/` group (AnnData encoding 0.1.0). Read as a modern codes/categories Group it
  returns raw integers — which yielded 130 fabricated "captures" instead of the real 143.
- **Qiu / Chongqing** deposit filenames carry a sample suffix (`GSM5102900_HC1`) while Gate 1
  keys on the bare accession (`GSM5102900`). Unstripped, every row failed to join and reported
  as zero retained.
- **Wilk** ships RDS, so a Python-only reader silently produced NA. The count is available from
  `converted/<channel>/barcodes.tsv`.

### What each cohort actually has

`n_input_barcodes` is **NA for all 746 captures**. No deposit contains an unfiltered droplet
matrix and none carries a Cell Ranger `metrics_summary`. That settles the *instrument-level*
question — but it does not settle the evidence question, because nominal loading is reported
for four cohorts and upstream doublet filtering has already happened in five.

### Evidence grading — three orthogonal fields

One status could not express a cohort whose loading is known while its residual stays unknown,
so the grading uses three independent axes. Full notes in `prior_evidence_verdicts.csv`.

| Cohort | Captures | `loading_evidence_status` | `upstream_doublet_handling` | `rate_model_status` |
|---|---:|---|---|---|
| qiu | 12 | LOADING_KNOWN | NONE_REPORTED | PLATFORM_MODEL_APPLICABLE_WITH_UNCERTAINTY |
| wilk | 23 | LOADING_KNOWN | AUTHOR_HEURISTIC_MULTIPLET_FILTER | PLATFORM_MODEL_APPLICABLE_WITH_UNCERTAINTY |
| combat | 70 | LOADING_KNOWN | GENOTYPE_DEMUX | PLATFORM_MODEL_EXTRAPOLATION_REQUIRED |
| ren | 284 | PARTIAL_OR_UNMAPPABLE_LOADING | AUTHOR_SCRUBLET_PLUS_MANUAL_FILTER | RESIDUAL_RATE_UNKNOWN |
| stephenson | 143 | PARTIAL_OR_UNMAPPABLE_LOADING | MIXED_OR_SITE_SPECIFIC | RESIDUAL_RATE_UNKNOWN |
| liu | 55 | LOADING_UNAVAILABLE | GENOTYPE_DEMUX | RESIDUAL_RATE_UNKNOWN |
| cillo | 35 | LOADING_UNAVAILABLE | HTO_FILTER | RESIDUAL_RATE_UNKNOWN |
| chongqing | 124 | LOADING_UNAVAILABLE | NONE_REPORTED | MODEL_UNAVAILABLE |

**No cohort is `DIRECT_RESIDUAL_ESTIMATE`.** `build_doublet_prior_evidence_layer.py` asserts it.

Three points that are easy to lose and change what may be claimed:

- **Ren** — the authors' `expected_doublet_rate = 0.08` was their Scrublet *input prior*, not a
  measured removed fraction. Our Gate 2 input is downstream of their filtering: summing
  `Ncells_manual_removal` over the 284 samples gives 1,462,702, matching our per-capture
  called-cell counts sample by sample. So 0.08 must not be reused as our residual prior, and
  the residual prevalence after their filtering is unknown. This is "upstream filtering already
  happened, residual unknown" — not "the prior has been spent".
- **COMBAT** — 50,000 cells loaded per channel is beyond the range of the official 10x tables.
  The vendor model may not be extended there without stating the extrapolation.
- **Wilk** — a Seq-Well loading-density model does exist (Hughes 2020), so this is not "no model
  available". The uncertainty is cross-protocol transfer plus an unknown residual after the
  authors' heuristic multiplet filter.

### Capture coverage — 746, in three distinct kinds

| Kind | Captures | Meaning |
|---|---:|---|
| `mapped` | **580** | capture-level evidence joined to this exact capture |
| `inherited` | **23** | Wilk cohort protocol (15,000 cells/array, uniform), `capture_specific=false` |
| `unmapped_to_gem_well` | **143** | Stephenson: 98 sequence pools, GEM wells not identifiable |

**580 Gate 2 captures have capture-specific evidence linkage.** That is linkage, not
measurement: a mapped row may hold a directly measured count, a published protocol value, a
derived value, a proxy, or an upstream-processing record. Identity mapping does not
establish a directly measured residual doublet prior, and no T7 row supplies one. The other
166 are not capture-specific at all.

The evidence table's `capture_id` is **not a single key**: Ren/COMBAT/Liu resolve on
`capture_id`, Qiu/Cillo/Chongqing only on `repo_sample_id`. Joining on `capture_id` alone
silently resolves 409 of 746 and drops 337 without error, which is why the resolved key is
materialized in `capture_evidence_join_key.csv` rather than re-derived per consumer.

### V3 scenario tables — documented assumptions, not a rate

`T7_doublet_prior_evidence/` additionally holds the v3 adjudication
(`doublet-rate-adjudication-v3.md`) and two scenario tables. They record which rate assumptions
*could* be reasoned about and where each one is stage-mismatched. They do not select a rate.

| Table | Rows | `runtime_eligible` |
|---|---:|---|
| `rate-sensitivity-scenarios-v3.csv` | 21 | **false on 21/21** |
| `pooled-same-donor-pairing-scenarios-v3.csv` | 58 | **false on 58/58** |

**No T7 evidence or scenario table is a runtime prior configuration.** No Gate 2 code may consume
them, directly or through a wrapper. The presence of a numeric value, a readable column, or a row
labelled `middle` does not authorize runtime use. Promoting any scenario to a runtime parameter
requires a separate, explicit, auditable decision recorded in `docs/decisions.md`.

What the scenarios are, stated at their own stage:

- **COMBAT** — 8 / 16 / 24% are `rho_capture_pre_filter` **extrapolation** scenarios; 50,000
  cells per channel lies outside the range of the official 10x tables. The 0.57 / 1.14 / 1.71%
  figures are hypothetical post-perfect-demultiplexing products under an equal-mixture
  assumption (Σp² = 1/14). None is an observed Gate 2 residual rate.
- **Wilk** — 6.8 / 8.5% are pre-filter Poisson occupancy scenarios; the table deliberately
  records **no central value**. Seq-Well is not governed by the 10x tables, and the authors'
  heuristic multiplet filter plus Gate 1 block propagation to a known residual rate.
- **Qiu / Cillo / Chongqing** — called-cell scenarios describe a **pre-Gate-1** stage. They are
  not post-Gate-1 runtime priors.
- **Liu** — Σp² is a modelled same-donor pairing probability, not a measured residual rate. The
  fixed 30,000-barcode whitelist is not a recovered-cell count, and SNG/DBL/AMB stay per lane.
- **Ren** — residual prevalence is **not identifiable** and must not be derived by subtracting
  the removed-cell fraction.
- **Stephenson** — loading cannot be mapped to GEM wells; 143 captures stay
  `unmapped_to_gem_well`.

Section F of the adjudication document is a **prospective manuscript template** and is labelled
as such in the file. No grid has been run; no scores or thresholds have been recomputed; no
results have been reported across any grid; correction has not been uniform across cohorts.

### Consequence for Gate 2

The 2026-09-24 adjudication is
`reports/manuscript/revision/methodology_review_2026-09-24/REVIEW.md`. The prior grid is
`0.005, 0.01, 0.02, 0.04, 0.08, 0.12, 0.16, 0.20`, labelled
`assumed_algorithmic_stress_setting`; it is not an interval of cohort prevalence.

- Scores at each grid point must use that prior, through verified recalibration from the fixed
  neighbour fractions or equivalent wrapper calls. Merely changing a threshold on scores
  calibrated with another prior is not the same experiment. Keep input, simulations, PCA and
  neighbours fixed within capture/seed.
- Report threshold/status/coverage and call-set sensitivity. Rank and observed-versus-simulated
  AUC stability are mathematical/implementation checks, not labelled accuracy or validation.
- Loading `r`, composition `s` and class-specific survival are explanatory uncertainties;
  observed post-demux shares and `e=1` do not generate a residual runtime prior.
- Cillo/Liu/COMBAT residual-input diagnostics must carry `residual_null_unvalidated`; their
  filtered-input simulations need not match actual remaining doublet classes.
- Never force top-`r` deletion. Unscored, failed or unmapped cells remain unknown/NA.
- Distinguish nominal protocol loading evidence (present for some cohorts) from actual
  capture-specific measured loading/recovery (not established in the supplied evidence).
  Wilk lacks an applicable 10x model; this is not proof that no Seq-Well model can exist.

Missing recovered counts mean "not available in the reviewed evidence", not "unobtainable in
principle". New primary evidence can change that conclusion; it cannot be silently substituted
into the approved stress experiment after its outputs are examined.
