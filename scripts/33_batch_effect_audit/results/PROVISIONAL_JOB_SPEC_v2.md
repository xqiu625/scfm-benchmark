# Provisional batch-effect diagnostic job — specification for approval (v2)

**Status: specification only. No SLURM job has been submitted.** Per `CLAUDE.md`: "No SLURM
submission without explicit authorization, per job." v1 of this spec
(`PROVISIONAL_JOB_SPEC_v1.md`, kept for the audit trail, not deleted) had real defects a
2026-10-09 review found — a sampling bug that under-counted the Platelet lineage, an arithmetic
mismatch between its own table and its own acceptance text, and acceptance criteria that would
have passed a run with a silently-failed Harmony correction. All are fixed in the code this spec
now describes (commit `e93cc6a`); this document's own numbers are regenerated from that fixed
code, not copied forward from v1.

## What changed since v1, and why the numbers are different

The v1 estimate (665 donors → 13,232 cells) came from a sampler that capped cells **per donor
across both lineages combined** — a donor with cells in both Platelet and Monocyte_Myeloid could
lose most of the rarer lineage to the more abundant one in the random draw. Fixed: the cap now
applies **within each (cohort, donor, lineage) combination separately**. Re-running the identical
`--dry-run` command against the identical real input gives a materially different, and more
correct, estimate — nearly double the cell count, with the Platelet lineage especially recovering
cells the old sampler was losing (e.g. COMBAT's Platelet count: 43 → 1,775; Ren's: 481 → 3,444).

## Step 0 — a note on how this estimate was actually produced

The real `--dry-run` command below was run twice on the login node before succeeding: the first
two attempts were killed (`exit 137`, SIGKILL) while reading `decontx/cells.csv`, even though
`free -h` showed 441 GB free on that node at the time — consistent with this project's own
documented operational constraint that HPCC login-node process limits are not simply a function
of total node memory (`CLAUDE.md`: "Do not run the real 5 GB input on the login node... the
resulting exit code is a fact about the probe, not about the code"). The fix was a memory-pattern
change (read both files' `cell_id` column alone first, before loading the heavier 11-column
`raw/cells.csv`), not a logic change — the identity-check logic itself is unchanged and covered by
`tests/test_prepare_provisional_diagnostic_input.py::test_mismatched_raw_decontx_cell_order_is_rejected`.
After that fix, the dry-run succeeded on the login node in ~11 seconds. This is noted explicitly
because it is relevant to Step 1's resource request below: the real (non-dry-run) run also
performs this same two-file read before the much larger streamed `counts.mtx` reads, so the same
memory-pattern consideration applies, though a proper SLURM allocation (not a login-node shell)
removes the ambiguity entirely regardless of which exact mechanism caused the kill.

## Step 1 — input preparation (bounded, streaming; the non-streamed prefix verified live above)

**Command:**
```
python scripts/33_batch_effect_audit/code/prepare_provisional_diagnostic_input.py \
    --input-root /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/platelet_monocyte_benchmark_input_v2 \
    --out /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/batch_effect_audit/provisional_diagnostic_input_v2/input.h5ad \
    --n-per-donor 20 --seed 0
```

**Code identity:** commit `e93cc6a` (the commit immediately preceding this spec's own commit in
`git log`) — not a loose, unverified claim: `run_batch_diagnostics.py` (Step 2) refuses to run at
all unless `--code-revision` is supplied or `git rev-parse HEAD` resolves, and will fail the run
rather than publish `"UNKNOWN"` provenance.

**Real, verified input paths and sizes** (checked directly, unchanged from v1):
- `.../raw/cells.csv`, `.../decontx/cells.csv` — 208 MB each, 1,236,881 rows each, **confirmed
  identical `cell_id` ordering** by an actual run of the new preflight check (not assumed).
- `.../raw/counts.mtx` — **22.2 GB**. `.../decontx/counts.mtx` — **50.1 GB**.
- `.../genes.csv` — 11,709 genes.

**Real, verified estimated output size** — from an actual, successful (non-dry, metadata-only)
`--dry-run` execution against these exact files this session:

```
Patient-balanced cap (<=20 cells/donor/lineage, seed=0): 1,317 distinct (cohort,donor,lineage)
groups -> 24,841 cells total (11,843 Platelet + 12,998 Monocyte_Myeloid)
```

Per-cohort/lineage breakdown (the real numbers; sums independently verified: 11,843 + 12,998 =
24,841):

| Cohort | Platelet | Monocyte_Myeloid |
|---|---:|---:|
| COMBAT 2022 | 1,775 | 2,437 |
| Chongqing sepsis 2025 | 2,458 | 2,457 |
| Cillo 2021 | 385 | 580 |
| Liu 2021 | 660 | 660 |
| Qiu 2021 | 140 | 140 |
| Ren 2021 | 3,444 | 3,908 |
| Stephenson 2021 | 2,533 | 2,356 |
| Wilk 2021 JEM | 448 | 460 |

A `<out>.selection_manifest.csv` (cell_id, cohort, donor_id, lineage, orig_row_index for all
24,841 rows) is written automatically alongside the real (non-dry-run) output — the exact
selection is reviewable before or after the job runs, not just the aggregate counts above.

**Untested part, named exactly:** this dry-run verified everything through the patient-balanced
selection, including the raw/decontx cell-identity check, against the REAL files. It stops before
(does not perform) the two streamed `counts.mtx` reads. That exact streaming code path
(`read_mtx_columns`) is verified correct, including VALUE correctness (not just shape), against a
small real on-disk MatrixMarket fixture in
`tests/test_prepare_provisional_diagnostic_input.py::test_end_to_end_streaming_write_reload_on_real_small_mtx_fixture`
— and separately already used in production by `assemble_platelet_monocyte_benchmark_input.py`
against these same real 22 GB/50 GB files. It has not been re-run against these exact files at
this exact 24,841-cell selection in this session. That is what this job authorization would
exercise for real.

**Resources (Step 1):** `epyc` partition (default, 30-day cap — chosen over `short`'s 2h cap
because this is I/O-bound, sequentially scanning 72 GB total from shared storage once), 4 CPUs,
**32 GB memory** (bounded by the 24,841-column output plus the two full `cells.csv` loads, not by
the input file size — but given this session's login-node SIGKILL at a step using far less than
32 GB of genuinely necessary memory, this ceiling includes real margin, not a tight estimate),
**walltime ceiling 4:00:00**. Single job, no array, no checkpointing (one-shot, non-destructive —
a rerun is safe, nothing is deleted, and the new no-clobber check refuses to run into an existing
output path so a rerun cannot silently overwrite a prior attempt).

## Step 2 — diagnostic CLI on the prepared input

**Command:**
```
python scripts/33_batch_effect_audit/code/run_batch_diagnostics.py \
    --h5ad /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/batch_effect_audit/provisional_diagnostic_input_v2/input.h5ad \
    --out-dir /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/batch_effect_audit/provisional_diagnostic_input_v2/diagnostic_run \
    --population-stage "pre-Gate-2 decontx_input PROVISIONAL, not final (2026-10-09)" \
    --code-revision e93cc6a \
    --expression-source decontx_counts \
    --batch-col cohort --patient-col patient_id --lineage-col lineage --outcome-col outcome \
    --run-harmony --seed 0
```

(No `--embeddings`: this input has no FM embeddings — 5-FM full-scale extraction has not started.
This runs expression/PCA mode only, exercising the same code path the
`synthetic_fixture_v3` validation run exercised end to end this session, which completed with
`status: "complete"` — not `"complete_with_partial_results"` — for both lineages, both scopes.)

**Resources (Step 2):** same `epyc` partition, 4 CPUs, 16 GB memory, **walltime ceiling 1:00:00**
(the synthetic-fixture run, at ~1,500 cells/lineage after its own diagnostic cap, completed in
well under 10 minutes per lineage; this input's per-lineage counts above are of the same order of
magnitude before the diagnostic CLI's own `--diagnostic-n-per-patient` cap reduces them further).

**Acceptance criteria (strict — matches what the v3 synthetic run actually achieved, not an
aspiration):**
- Exit 0 AND `manifest.json`'s top-level `status` is exactly `"complete"` — **not**
  `"complete_with_partial_results"` and not `"failed"`. (v1 of this spec accepted
  `complete_with_partial_results`; a review correctly pointed out that is not an acceptable bar
  when Harmony was explicitly requested.)
- `manifest.json`'s `code_revision_verified` is `true`.
- `manifest.json`'s `input_sha256` matches a hash taken of Step 1's actual output file.
- Both lineages' `expression_pca` entries show `status: "complete"` for BOTH scopes (the Harmony
  `result` dict must contain real `batch_asw`/`neighborhood_mixing`/`bio_conservation` keys, not
  a `"blocked"`/`"failed"` placeholder).
- Per-lineage `n_diagnostic_cells` + `n_excluded_by_cap` + `n_excluded_outcome` reconciles against
  Step 1's own per-lineage counts above.
- `metrics_long.csv` has rows for both lineages × `expression_pca` × both scopes (4 combinations
  minimum), using the fixed `METRICS_COLS` schema even if a row is missing.
- Every saved `.npz` representation reloads with matching shape, cell-ID ordering, and values
  (enforced automatically by `save_representation()` — a failure here raises before the run can
  report success).

**Failure criteria (stop and report; do not retry automatically):**
- Non-zero exit, `manifest.json` absent, or `status` other than `"complete"`.
- Any `population_stage`/`expression_source` mismatch warning in the log (should not occur, since
  Step 1 declares both and Step 2 is given the matching values above — if one appears, something
  about the two steps' identity binding has drifted and must be resolved before interpreting
  results, not merely logged past).

## What this job does NOT authorize

- No cell deletion, no Gate 2 scoring, no change to any existing checkpoint or retained-cell set.
- No promotion of this run's numbers into the manuscript or `docs/current-status.md` as a final
  result — explicitly provisional, and the per-cohort imbalance above (Ren/Chongqing/Stephenson
  contribute far more cells than Qiu/Cillo) means any batch-effect reading from it must account
  for unequal cohort representation, not just interpret the aggregate number.
- No automatic retry if it fails; a failure gets reported, not resubmitted.
- No push and no change to `canonical_donor_labels.csv` or any other canonical artifact.

## One job or two?

Steps 1 and 2 are specified as **two separate, dependent SLURM allocations** (Step 2 cannot start
until Step 1's output file exists and is verified) rather than one combined job: they have
different resource profiles (Step 1 is I/O-heavy with modest CPU/memory needs; Step 2 is
compute-heavy on a small bounded input) and separating them means a Step 2 failure does not
require re-running Step 1's 72 GB scan. If authorized, **each needs its own per-job
authorization** per `CLAUDE.md`, requested separately (Step 2 only after Step 1's output and
acceptance checks are confirmed, not pre-authorized speculatively).
