# Stage 31 `n_prin_comps` — deferred pending workflow redesign

**Current ruling, 2026-09-24:** the workflow redesign is conditionally approved in
`reports/manuscript/revision/methodology_review_2026-09-24/REVIEW.md`. No fixed dimension is
selected. First validate complete pre-subsample physical-capture RNA inputs and barcode mapping,
then audit dimensions on that new population. The old 301/161 table cannot select a parameter
for the new inputs. Coverage accounting must include every target unit; failure is unknown,
not zero, and no sample may receive an automatic lower dimension.

The follow-up Gate 2 eleven-parameter ruling fixes PC=30 **only for Gate 2's scoped scoring
pilot**. It does not close this Stage 31 decision. The owner's latest instruction delegates
authorization of a new dimension-only audit to the reviewer after raw-source provenance,
physical-capture identity and raw-barcode/canonical mapping have closed on the complete input.
That authorization has not been issued; it would permit measurements in a new output directory,
not doublet detection, outcome analysis or an advance choice of PCA dimension.

The material below is the historical technical rationale for the original Stage 31 audit.
It is not a pending instruction to choose a dimension on the old 60-cell analytic subset, nor
an adjudication of Gate 2's separate production parameter set. The seven-file review snapshot
remains unchanged.

**Status: OPEN. No value has methodological sign-off, and the script has no default: it
refuses to run without an explicit `--n-prin-comps`. This note exists so the value is chosen on
stated grounds rather than back-derived from whichever samples happened to fail — and so that
"not yet decided" is enforced by the code rather than asserted in a comment.**

Nothing in this note was computed by running the detector to see which value maximises
coverage. Doing that would tune a detector parameter against the data it is about to judge.

---

## Why this became a question

The v1 run (`contamination_audit_v2_fixed_semantics_2026-09-22_v1`, job 29010326) passed every
technical-boundary check and **failed dimension feasibility**. Two facts surfaced:

- `n_prin_comps` was **never passed**. Both call sites took scanpy's implicit default, so the
  run had a dimensionality nobody had chosen and nothing recorded.
- Four samples of exactly 30 cells raised
  `n_components=30 must be strictly less than min(n_samples, n_features)=30`, and the batched
  path died with `min(n_samples, n_features)=17`.

**What the second message does and does not say.** It reports a *minimum*. It does not identify
which term was 17 — `n_cells_effective`, `n_usable_features`, or both. An earlier version of this
note read it as "17 usable features"; that was an over-reading and is retracted. `filter_cells(min_genes=3)`
can reduce the cell count too, so either term could be the binding one. What the message does
establish is that **something other than the raw 13,066 genes and the raw sample size governs
feasibility**, which is why both quantities are now measured per sample instead of inferred.

## The constraint, stated correctly

scanpy 1.11.3 `pp.scrublet` preprocesses each batch as

```
filter_genes(min_cells=3) -> filter_cells(min_genes=3)
normalize_total -> log1p(temp layer) -> highly_variable_genes
subset to highly_variable -> ... -> PCA(n_components=n_prin_comps, svd_solver="arpack")
```

so the arpack requirement

```
n_prin_comps  <  min(n_cells_effective, n_usable_features)
```

involves two quantities that the raw sample size hides:

| Quantity | Why it is not what it looks like |
|---|---|
| `n_cells_effective` | after `filter_cells(min_genes=3)` — **not** the number of cells in the sample |
| `n_usable_features` | the **HVG count**, not the gene count of the matrix |

Both are now measured per sample by `scrublet_preprocessing_dimensions()`, which calls the same
scanpy functions in the same order rather than approximating them, and both appear in
`per_sample_diagnostics.csv`.

## What must NOT decide the value

- **Tuning either knob until the observed failures clear.** Two different moves get conflated
  here, and an earlier version of this note stated the arithmetic wrongly:

  | Move | What it actually does |
  |---|---|
  | `n_prin_comps` 30 → **31** | makes the constraint *harder*, not easier: `30 < 30` is false and `31 < 30` is false too. It fixes nothing. |
  | `MIN_CELLS_PER_SAMPLE` 30 → **31** | excludes the four 30-cell samples one step earlier. The failure disappears from the log; the samples are still not scored, and it mislabels a *dimensionality* constraint as a *cell-count* one. |
  | `n_prin_comps` **lowered** | can clear equality failures — and is a genuine algorithm configuration change, not a bug fix. It also does nothing for a batch whose minimum was 17 unless taken below that. |

  None of the three is a reason to pick a value.
- **Lowering it until coverage looks acceptable.** Coverage is an outcome of the choice, not an
  argument for it.
- **Per-sample adaptation.** Rejected: a sample-specific `n_prin_comps` gives each sample its
  own detector working point, and the per-patient rate aggregates across samples. The code has
  no fallback-to-smaller path, by design.

## What should decide it

1. **The published Scrublet default and its calibration basis.** Scrublet's own default is 30,
   established on typical nucleated-cell libraries with thousands of cells and hundreds to
   thousands of HVGs. Whether that calibration transfers to anucleate platelet fragments at
   ≤60 cells per sample is the actual question, and it is not answered by this repository.
2. **The variance actually available.** `min(n_cells_effective, n_usable_features)` fell to at
   least one value of 17 in v1, though which term that was is not recoverable from the message.
   If the binding term turns out to be the HVG count, a 30-dimensional embedding is asking for
   more components than the data has directions, and a defensible value may be well below 30.
   `per_sample_diagnostics.csv` from the next run reports both terms per sample, which is what
   this decision should rest on — not on the v1 error text.
3. **The consequence for eligibility, stated in advance.** Any value implies a minimum
   `min(n_cells_effective, n_usable_features)`. That implication should be written down and
   accepted *before* the run, so the resulting coverage is a prediction being tested rather than
   a result being rationalised.

## Breaking the circularity: the preflight

Requiring an approved dimensionality before the script will run, while saying the choice should
rest on measured per-sample dimensions, is circular — the distribution could only be obtained by
first picking a value arbitrarily, which is exactly what this note forbids.

`--dimension-audit` resolves it. It loads the cohort under the same outcome allowlist and the
same fixed subsample, runs scanpy's own scrublet preprocessing to **attempt** a measurement of `n_cells_effective`
and `n_usable_features` **for every sample** — successes and failures both recorded, since a
measurement can raise — writes `dimension_audit.csv`, prints percentiles,
and stops. It never calls Scrublet, produces no call, rate, contamination or detection-diagnostic
output, requires `--integration-only`, and records `scrublet_attempted: false` and
`scrublet_completed: false`.

**Measurement is not gated on eligibility.** An earlier version skipped the samples below
`MIN_CELLS_PER_SAMPLE`, which silently redefined "the per-sample dimension distribution" as "the
distribution over the samples the current cell-count rule already admits" — a detection-stage
decision leaking into the measurement meant to inform the next one. Every sample is *attempted*;
`meets_min_cells_rule`, `preflight_exclusion_reason` and `dimension_measurement_error` sit beside
the measurement, and the run
prints **two** distributions:

| Distribution | What it is |
|---|---|
| all measured samples | what the data can support, independent of the current rule |
| measured **and** meeting the current cell-count rule | what would be scored if `MIN_CELLS_PER_SAMPLE` stays as it is |

They are not interchangeable. **The parameter decision must state which one it used** — and the
preflight now states whether each is even usable.

**Acceptance.** Catching every per-sample measurement error and then writing `completed`
unconditionally meant a preflight in which all 301 samples failed would still exit 0 and look
usable, with the errors buried in a CSV column. The manifest therefore carries
`n_samples_total`, `n_samples_measured`, `n_measurement_errors`, `n_eligible_samples`,
`n_eligible_measured` and `dimension_audit_sufficient_for_decision`, the last split per
distribution:

| Decision uses | Requires |
|---|---|
| all-samples distribution | every sample measured **successfully** |
| eligible-only distribution | every sample meeting the cell-count rule measured **successfully** |

`run_status` is `completed` only when every sample was measured successfully; otherwise
`incomplete`.

**And the process exits non-zero (3) when neither distribution is complete.** Fixing the manifest
alone left two contradicting signals: SLURM recording COMPLETED while the manifest said the run
supported no decision. `write_manifest()` has already returned —
so `os.replace()` has succeeded — before the status is decided, and the exit that follows cannot
remove the replaced manifest. The claim is scoped to **this control flow**, not to signals or
crashes generally: an interruption before the replace leaves the previous final manifest in place
**if one exists** — during the first write, with the output guard requiring a new empty directory,
there may be no final manifest at all, and an uncatchable interruption can leave a `.inprogress`
file behind. And `write_manifest` does not fsync, so it provides **no guarantee** of survival across power
loss, node failure or storage failure — which is not the same as guaranteed loss: without fsync
the data may still survive, through OS writeback, the shared filesystem's own behaviour, or
simply when the failure lands. What is absent is the guarantee, not the possibility. The wording that said it was "fsynced" was carried over from `gate2_checkpoint.py`
without that module's fsync implementation. An incomplete all-samples distribution beside a
complete eligible-only one still exits 0 — it supports an explicitly narrowed decision, and the
methodological conclusion must then say it used the eligible-only distribution.

A sample that cannot be measured at all is a third fact, distinct from "measured and not enough"
(`insufficient_dimension`) and from "the detector broke" (`runtime_failed`): in detection mode it
is marked `dimension_measurement_failed` and reaches neither the feasibility verdict nor
Scrublet. An earlier version let such a sample fall through to
`dimension_feasibility(n_eff, ...)` with the names unbound, raising `UnboundLocalError` and
destroying the recorded reason along with it.

**Why the preflight refuses `--n-prin-comps`** — governance, not measurement. The value does not
enter the preprocessing, so supplying one would not change a single measured number. It would put
an unapproved value into the command line and the manifest of the run whose purpose is to inform
that choice, and so read as though it were already settled or applied. (An earlier version of this
note said the distribution would be measured "under" the value; that was wrong.)

**Order of operations:** preflight into its own **new, isolated output directory** — it writes a
manifest and `dimension_audit.csv`, so it is not a read-only run; what it leaves untouched is the
inputs and every existing artifact — under its own authorization → methodological decision on the
value, with the implied eligibility written down in advance → v2 in a new directory with that
value passed explicitly.

## What the code does in the meantime

- **no default at all.** `--n-prin-comps` is a required argument and `scrublet_params()` refuses
  a missing or invalid value, so the script cannot reach Scrublet under an unapproved parameter.
  A provisional constant would not have achieved this: recording it in the manifest still lets
  anyone produce results, and "provisional" in a comment binds nothing;
- the supplied value is passed to the batched path and the per-sample fallback alike, and
  recorded in `RUN_MANIFEST.json` under `params.N_PRIN_COMPS` and `params.scrublet_params`;
- feasibility decided **before** Scrublet is called, per sample;
- infeasible samples marked `insufficient_dimension` — a distinct status, never folded into
  `runtime_failed`, with `doublet_call` and `doublet_score` left `<NA>`;
- `scrublet_attempted` and `scrublet_any_call_completed` are separate manifest facts: attempt is
  recorded before the call, since a detector that raises at its entry point was still attempted,
  and completion only after a successful return. A single `scrublet_executed` written before the
  call claimed execution for a crash. The second is named "any call" deliberately — in the
  fallback the first sample that returns sets it while later samples may still fail, so it means
  "at least one call returned", **not** "the detection stage finished". Stage completion is
  `run_status`; per-sample outcomes are in `per_sample_diagnostics.csv`;
- no automatic retry at a reduced dimension.

## The limit of what fixing this buys

A fixed dimensionality removes **one** source of incomparability between samples. It does not
make Scrublet rates comparable across samples, and it does not make them biologically valid
here. Sample sizes, score distributions and automatic thresholds still differ per sample, and a
cap of 60 cells per sample makes any per-sample estimate unstable on its own terms.

**The re-run therefore remains a technical diagnostic. It is not evidence of doublet prevalence
and it does not authorise deleting any cell.**
