# STALE — stage 31 contamination audit, invalid under old doublet semantics

**Registered 2026-09-22 (execution-plan P0b). These artifacts must not enter the manuscript
evidence chain, and must not be cited, until the code fix below has been re-run under a separate
authorization.**

Location on HPCC:
`/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/bib_revision_tier2/contamination_audit/`

| Artifact | Status |
|---|---|
| `doublet_scores_per_cell.csv` | **stale** — `predicted_doublet` column holds `False` for cells that were never examined |
| `doublet_rate_per_patient.csv` | **invalid** — patients with no examined cells appear as a measured 0% doublet rate |
| `contamination_audit_v2_doublet_excluded.csv` | **invalid** — computed on `a_clean`, in which unexamined cells were indistinguishable from examined non-candidates |
| doublet-rate × mortality Spearman (printed, and the rate column above) | **invalid** — computed over rates that include fabricated 0% values |

Nothing here has been overwritten, deleted or moved. The files remain exactly as produced, as the
record of what the old semantics generated.

## What was wrong

`doublet_detection_and_reaudit.py` wrote `predicted_doublet = False` for every cell in its
fallback branch, and left it `False` for samples with fewer than 30 cells and for samples whose
Scrublet call raised. Those cells were **never examined**, but `False` asserts *examined and not
a doublet*. The error propagated through four layers:

| Layer | Old behaviour |
|---|---|
| per-cell export | `False` exported for unexamined cells |
| per-patient aggregation | `groupby(...).mean()` counted them as 0, so an entirely unexamined patient scored a measured 0% |
| exclusion mask | `a[~obs["predicted_doublet"].fillna(False)]` — cells were kept, but `a_clean` then travelled downstream as a "doublet-cleaned" set with no way to tell unexamined from examined |
| statistics | Spearman over those rates, including the correlation against mortality |

The mortality correlation is the most consequential: unexamined samples entering as 0% is not
outcome-neutral.

## What was fixed (code only — nothing re-run)

The state is now a nullable boolean `doublet_call` (True / False / `<NA>`) plus an explicit
`examination_status` (`examined` / `too_small` / `threshold_failed` / `runtime_failed` /
`not_evaluated`), bound by an enforced invariant:

```
examined        <=>  doublet_call in {True, False}
unknown status  <=>  doublet_call is NA
```

A first pass closed only the explicit `fillna(False)` collapse and left a second path open: a
cell marked `examined` whose call is NA entered the rate denominator without contributing to
the numerator, diluting unknown into non-doublet just the same. Measured on that intermediate
code: one confirmed doublet among one examined cell plus three examined-but-NA cells reported
**0.250 instead of 1.000**. `validate_examination_state` now raises on both contradictions and
runs before export, aggregation and filtering.

- unexamined cells stay `<NA>`, never `False`;
- only `doublet_call == True` is removed, and retained unknown cells stay distinguishable in
  `a_clean.obs` — the result is explicitly **not** a doublet-cleaned set;
- per-patient rate uses **examined cells as the denominator** and is `NA`, not 0%, when a patient
  has none; `n_examined`, `n_unknown` and `examined_coverage` travel with it;
- the Spearman analyses exclude undefined rates explicitly, keep those patients visible in the
  exported table with `included_in_spearman=False`, and report the contributing patient count.

A sample is marked `examined` only after **three** things hold: a finite automatic threshold
read from the sample's own Scrublet metadata, strictly boolean (or 0/1 integer) calls of the
right length with no NA anywhere, and finite scores of the right length. Any failure leaves the
whole sample unknown.

Reading the threshold is the check the stage-32 failure actually needs. Scanpy *warns* rather
than raising when automatic thresholding fails, then returns finite scores and an all-False
prediction — output that passes every field-level test while resting on no threshold at all. An
all-False prediction is **not** itself evidence of failure (it can be legitimate), so the
threshold has to be read. On the batched path scanpy exposes per-batch metadata under
`uns['scrublet']['batches'][sample]`; if that is absent the run falls back to per-sample
execution rather than accepting samples on the strength of a pooled threshold.

`too_small` is decided **only** by the cell-count rule, never inferred from a NaN score, since a
NaN score in a large sample means the detector failed — a different fact.

The invariant requires the call column to be a **nullable boolean dtype**, not merely non-NA:
`"False"`, `0.0` and `2` are all non-NA and would otherwise pass while meaning nothing.

Covered by `tests/test_stage31_doublet_semantics.py` (188 tests), including that a genuinely
measured 0% and an undefined rate remain distinguishable, and that every contradictory
status/call combination is refused.

## The outcome boundary of `--integration-only`, stated exactly

The re-run can be validated technically without executing an outcome analysis. What that flag
does and does not guarantee, per input — the two are **not** the same guarantee, and claiming one
for both would be false:

| Input | Guarantee |
|---|---|
| `patient_qc_distributions.csv` | **strict** — `y` is excluded by `usecols` at parse time and is never constructed in memory. `usecols` is also fail-closed on a rename |
| `human_platelet_covid_severity_raw_counts.h5ad` | **incidental load, acknowledged** — `read_h5ad` cannot select obs columns, and this file's obs carries `Death` plus severity labels. They enter memory with the container. Only the columns the loader requires survive it; the rest are dropped before any use and recorded |

So the claim is not "outcome is never loaded". It is: excluded at parse time where the reader
allows it; where it does not, loaded incidentally and never inspected, propagated, exported or
analysed. The strict form for the h5ad would need a separate outcome-free technical input, which
does not exist and is not created by this run.

An earlier version claimed the strict guarantee for both while doing neither: it read the whole
CSV and then called `.drop(columns=["y"])`, and on the h5ad it dropped `"y"` — a column this file
does not have — while `Death` sat in obs untouched.

### Retention, not a deny-list

The first fix replaced that with a deny-list, `("y", "Death")`, and was wrong the same way: it let
`covid_severity`, `mild_vs_severe`, `severity_score`, `three_group` and `progression_group`
propagate into the AnnData, the slice and the temporary h5ad round trip, while the code and this
document said "outcome" throughout. A deny-list is only as good as the classification behind it,
and a clinical label added upstream would propagate by default.

So the direction is inverted. The loader needs exactly **one** column from the container —
`Sample_ID`, to map cells to patients; everything else it uses, it creates. `REQUIRED_OBS_COLUMNS`
is kept and all other columns are dropped, which makes the boundary independent of whether the
classification is complete: an unadjudicated clinical label cannot leak, because nothing but
`Sample_ID` survives. A missing required column is fatal.

`KNOWN_OUTCOME_OBS_COLUMNS` then serves only the *report*: the manifest lists
`known_outcome_dropped` separately from `unclassified_dropped`, so a column nobody has adjudicated
is visible as unadjudicated rather than quietly filed as technical. An omission there understates
the report; it cannot cause a leak.

### The manifest separates authorization, loading, inspection and execution

`outcome_analysis_requested` is known up front and is the only field that may be predicted.
Everything else starts at "has not happened" and is flipped by the code that does the thing, each
with an immediate atomic rewrite:

| Field | Under `--integration-only` |
|---|---|
| `known_outcome_values_loaded` / `_via` | **true** — the h5ad container carried `Death` and the severity labels into memory |
| `outcome_values_inspected` / `_via` | **false** — nothing read their values |
| `outcome_analyses_executed` | `[]` |
| `container_obs_columns` | what arrived, what was kept, what was dropped and under which classification |

The `known_` prefix is load-bearing: that flag speaks for `KNOWN_OUTCOME_OBS_COLUMNS` only. A
container carrying nothing but a new, unrecognised clinical label would leave it `false` while
`container_obs_columns.unclassified_dropped` names that label — which is why the scope limitation
is also written into the manifest record itself, under `classification_scope`. The anti-leak
guarantee never rests on that flag; retention does.

### The two records sit on opposite sides of the drop

Loading is recorded as soon as `read_h5ad` returns — the values are in memory at that moment, and
a job dying there did load them. The kept/dropped record is written only **after** the allowlist
has been applied and verified, because until then no column has actually been dropped. An earlier
version wrote `known_outcome_dropped` before applying the allowlist: a job dying in between left a
manifest asserting columns were gone while they were still on the object — the same
predicted-fact defect as `outcome_data_accessed`. `kept` is taken from `a.obs.columns` as observed,
not from the constant that was requested, so the record cannot be right while the object is wrong.

Loading and inspection are separate fields because one could not express the real state without
contradicting itself: the previous manifest recorded `outcome_data_accessed: false` beside
`outcome_incidentally_loaded: ["Death"]`, and "loaded but never looked at" is not "not accessed".
That field name is retired.

A job that completes the mortality correlation and then dies leaves `run_status: started` **and**
that analysis listed — both facts, neither inferred from the other.

### Verified by mutation, not by the suite passing

The tests run the real loader on fixture inputs and read the real manifest off disk. Eleven
mutations of the shipped code were each confirmed to fail the suite: parse-then-drop,
in-memory-only recording, predicted access, record-but-don't-drop, drop-`y`-ignore-`Death`,
deny-list instead of retention, incidental load not flipping the loaded flag, severity removed
from the classification, loading also claiming inspection, recording the drop before applying it,
and applying the drop without verifying it.

One mutation — taking `kept` from the constant rather than from the object — survives, and is
recorded here rather than papered over: on the success path the two agree, so no test can
distinguish them. What that mutation removes is the protection against a *failed* drop, and that
is covered instead by the explicit verification, whose own removal is caught.

Two earlier mutation attempts appeared to survive and did not: the string replacement missed
because a comment sat between the two statements it was joining, so the code was never mutated. A
mutation that "survives" has to be shown to have applied.

## What is still required before these artifacts can be replaced

1. A separate authorization to re-run the analysis.
2. Re-running writes to a **new versioned directory**; the files above are not overwritten.
3. A re-audit of the regenerated numbers before anything re-enters the evidence chain.

Re-running was deliberately **not** performed as part of the semantics fix: changing the rules
and the numbers in one step would leave no way to attribute a difference to either.


---

# v1 re-run (job 29010326) — technical boundary PASSED, dimension feasibility FAILED

`05_RESULTS/bib_revision_tier2/contamination_audit_v2_fixed_semantics_2026-09-22_v1/`
**must not enter the manuscript evidence chain.** It is retained as the record of the
boundary verification, not as a result. `--integration-only`; no outcome analysis ran.

## What the old semantics actually did, as a four-cell table

Read-only comparison of the two per-patient rate tables, 240 patients, all matched:

|  | new NA | new == 0 | new > 0 | old total |
|---|---|---|---|---|
| **old == 0** | **84** | **1** | **1** | 86 |
| **old > 0** | 1 | 0 | 153 | 154 |

- **84** — old reported a measured 0% for patients with no examined cell at all: fabricated zeros.
- **1** — the only genuinely measured 0% (60 examined cells, no doublet).
- **1** — CV0104: old reported 0%, actual rate **7.84%** (4 of 51 examined). A *wrong non-zero*,
  not a fabricated zero, and invisible to a comparison of marginal totals. The likely mechanism
  is the old code accepting an all-False prediction that rested on no threshold, but the v1-era
  run recorded no threshold or status diagnostics, so **this cannot be shown, only hypothesised**.
- **1** — old had a non-zero rate for a patient that has no validly examined cell now.

## Why v1 fails acceptance

`n_prin_comps` was **never passed**: both call sites took scanpy's implicit default, so the run
had a dimensionality nobody chose and nothing recorded. Four 30-cell samples then raised
`n_components=30 must be strictly less than min(n_samples, n_features)=30`, and the batched path
died with `min(n_samples, n_features)=17`. That message reports a *minimum* and does not say
which term was 17, or whether both were — an earlier reading of it as "seventeen usable
features" was an over-inference and is retracted. What it does establish is that feasibility is
not governed by the raw sample size or the 13,066 genes in the matrix, which is why both
`n_cells_effective` and `n_usable_features` are now measured per sample.

## What v2 must do differently

Fixed dimensionality, chosen on methodological grounds, **not** back-derived from the samples
that failed. Note the arithmetic, which an earlier version of this file got wrong: raising
*`n_prin_comps`* 30 → 31 fixes nothing (`31 < 30` is false as surely as `30 < 30`); raising
*`MIN_CELLS_PER_SAMPLE`* 30 → 31 only excludes those four samples one step earlier; and
*lowering* `n_prin_comps` is an algorithm configuration change, not a bug fix. See
`docs/scrublet_dimension_parameter.md`, which is open — and the `--dimension-audit` preflight
that must run before that decision, because requiring an approved value in order to measure the
dimensions the value is chosen from is circular. The preflight attempts a measurement on **every** sample, not
only those the current cell-count rule admits, and records failures as such, and reports both distributions separately; it
writes into a new isolated directory of its own (a manifest and `dimension_audit.csv`), leaving
the inputs and every existing artifact untouched. Its manifest states whether the result is even
usable: `run_status` is `completed` only if every sample was measured successfully,
`dimension_audit_sufficient_for_decision` is reported separately for the all-samples and
eligible-only distributions, and the process **exits 3** when neither is complete so the
scheduler cannot record a run that supports no decision as COMPLETED.

- **no default dimensionality**: `--n-prin-comps` is required and `scrublet_params()` fails
  closed, so the script cannot run under an unapproved value. The supplied value goes to the
  batched path and the fallback alike; no retry at a reduced dimension;
- `n_cells_effective` and `n_usable_features` measured by calling scanpy's own scrublet
  preprocessing in the same order — `filter_cells(min_genes=3)` changes the sample size and HVG
  selection sets the feature count, so neither is what the raw sample looks like;
- feasibility decided **before** Scrublet is called: `n_prin_comps < min(...)`, strict;
- infeasible samples marked `insufficient_dimension`, a distinct status, never folded into
  `runtime_failed`, calls and scores left `<NA>`;
- `per_sample_diagnostics.csv` per sample, reconciled against the cell-level counts by
  assertion, and distinguishing `scrublet_attempted` from a pre-call rejection:
  `n_prin_comps_used_for_attempt` is filled for every sample that entered the detector —
  including one later rejected on its threshold — and empty only for `too_small` and
  `insufficient_dimension`.

**A fixed dimensionality removes one source of incomparability. It does not make Scrublet rates
comparable across samples or biologically valid here** — sample sizes, score distributions and
automatic thresholds still differ, and 60 cells per sample makes any per-sample estimate
unstable. The v1 figure of 30.1% doublets among examined cells is **uncalibrated for this
material, subsampling and parameterisation**; it cannot be read as doublet prevalence and
authorises no deletion. Per-sample rates, thresholds, sizes and score distributions go to
methodological review; a generic 10x expectation does not adjudicate them.
