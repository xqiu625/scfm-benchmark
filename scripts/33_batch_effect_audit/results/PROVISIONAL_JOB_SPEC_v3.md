# Provisional batch-effect diagnostic job — specification for approval (v3)

**Status: specification only. No SLURM job has been submitted.** Per `CLAUDE.md`: "No SLURM
submission without explicit authorization, per job." v2 of this spec
(`PROVISIONAL_JOB_SPEC_v2.md`, kept for the audit trail, not deleted) had a real defect a
2026-10-09 review found: its own Step 2 command's `--population-stage` value
(`"pre-Gate-2 decontx_input PROVISIONAL, not final (2026-10-09)"`) did not match the string
Step 1's own code actually writes to `uns["population_stage"]`
(`"pre-Gate-2 PROVISIONAL, not final (see docs/current-status.md)"`) — so following v2's
documented commands verbatim would itself have triggered the exact mismatch warning the spec's
own failure criteria said "should not occur". The review found four further execution-integrity
gaps (code-identity verification, concurrent-writer safety, selection/eligibility consistency,
corrected-UMAP status propagation); all five are fixed in the code this spec describes, each with
a passing regression test, not merely documented as a risk.

## What changed since v2, and why

1. **Stage identity is no longer two independently-authored free-text strings.** A new
   machine-checkable `--stage-id` (matched exactly against the input's own
   `uns["population_stage_id"]`, `scripts/00_shared/code/population_stage.py`) now HARD FAILS on
   mismatch, and Step 2's sbatch script (below) reads `--population-stage` directly from Step 1's
   actual output file at runtime rather than retyping it, so the two steps' free-text label
   cannot drift the way v2's did.
2. **`code_revision_verified` is a real check, not `code_revision is not None`.** v2 asserted
   `commit e93cc6a` without the code actually verifying that claim against anything; a review
   found the field would report "verified" for literally any string. Fixed
   (`_resolve_code_identity()`): verified now means a real git checkout, resolved `HEAD`, a clean
   working tree, and a match to the asserted value. `source_file_sha256` is always recorded as
   the honest substitute on a non-git scratch deployment.
3. **Atomic run ownership.** Both CLIs now claim an OS-atomic exclusive-create lock
   (`scripts/00_shared/code/run_ownership.py`) before writing, closing a real TOCTOU race the
   prior check-then-act guard had — proven with a genuine two-thread/two-process race test, not
   just reasoned about.
4. **Technical-integration eligibility is no longer coupled to outcome availability.** v2's
   `run_batch_diagnostics.py` excluded unknown-outcome cells from the ENTIRE diagnostic
   population (PCA, Harmony, neighbor-mixing), not just from the bio-labeled metrics that
   actually require a real outcome. Fixed: the full lineage-matched population (including
   unknown-outcome cells) is sampled and used for PCA/Harmony/neighbor-mixing;
   `bio_known_mask` restricts only `batch_asw` and `bio_conservation`. Both cell counts
   (`n_cells_technical`, `n_cells_bio_known`) are reported, not collapsed into one number.
5. **`patient_id` is cohort-qualified** (`f"{cohort}::{donor_id}"`), closing a real identity
   hazard: the prior bare `donor_id` could, in principle, merge two different patients from two
   different cohorts into one "patient" if their donor_id strings ever coincided.
6. **A failed corrected-UMAP no longer leaves `status` unchanged.** Status determination was
   extracted into a pure, directly-tested function (`determine_embedding_status()`) so the
   corrected-UMAP-failure case — asymmetric with the already-handled uncorrected case — is
   covered by a parametrized regression test, not left to a flaky real failure to exercise.
7. **The manifest is written as `status: "running"` before any fallible I/O** (code-identity
   resolution, input hashing), not constructed inline with the hash call outside the
   failure-handling `try`/`except` the way v2's was — closing a gap where a hashing failure left
   no persisted record of the run having started at all.
8. **A frozen, reviewable selection can now be bound to a real extraction.** `--dry-run` no
   longer deletes its own selection manifest (v2's did); a `.selection_manifest.status.json`
   sidecar distinguishes `"preflight_only"` from `"extracted"`, and a new `--selection-manifest
   <path>` lets a real run consume and re-validate a previously-reviewed selection instead of
   only recomputing deterministically from the same seed.

## Addendum (2026-10-10): both wrappers now bind to an isolated, verified deployment

A 2026-10-10 review found that both sbatch scripts below `cd`'d into
`/bigdata/godziklab/shared/Xinru/302006/302006_updategit` — the project's main working checkout
— and Step 2 resolved `--code-revision` from that same path. Checked directly: that checkout is
at commit `35131cf` (dozens of commits behind this spec's own `HEAD`) with five files showing
uncommitted local modifications (confirmed via `git -C
/bigdata/godziklab/shared/Xinru/302006/302006_updategit status --porcelain` on 2026-10-10). A job
submitted with the wrappers as they stood would have executed that old, dirty code and reported
its stale, dirty revision as the "resolved" `--code-revision` — silently failing this spec's own
acceptance criterion that the job run from a real, clean checkout (see below), without the wrapper
ever detecting it.

**Fixed.** A new isolated checkout was deployed specifically for this job, independent of the
main working copy (which is left exactly as found, uncommitted changes and all — not a
prerequisite for this job per "What this job does NOT authorize", below):

- Local `git clone --no-hardlinks` of this repo's own working tree at `HEAD` (`9ea4c6a`) into a
  fresh directory, confirmed clean (`git status --porcelain` empty — a plain `rsync` of the
  working tree would NOT have been clean, since four pre-existing untracked files live there).
- `rsync`'d to `/bigdata/godziklab/shared/Xinru/302006/302006_updategit_deploy_9ea4c6a` on HPCC.
- Verified independently ON HPCC (not just locally): `git rev-parse HEAD` ->
  `9ea4c6a43a4ec216c345c35a463943f7b550232d`; `git status --porcelain` -> empty.
- Both `step1_prepare_diagnostic_input.sh` and `step2_run_diagnostics.sh` now set
  `DEPLOY_ROOT=/bigdata/godziklab/shared/Xinru/302006/302006_updategit_deploy_9ea4c6a` and `cd` into
  it; both resolve a `CODE_REVISION` from `$DEPLOY_ROOT` and hard-fail (`exit 1`, before running
  any Python) if it does not resolve or does not equal `9ea4c6a43a4ec216c345c35a463943f7b550232d`
  — so a future edit to either the deployment or the wrapper that silently drifts them apart is
  caught at submission time, not discovered after the job ran.
- Checked (not assumed): `bash -n` on both scripts (syntax), the interpreter path
  (`/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean/bin/python`, exists and runs on HPCC,
  `Python 3.10.18`), Step 1's `--input-root` (exists), and that
  `05_RESULTS/batch_effect_audit/provisional_diagnostic_input_v3/` does not yet exist as a
  non-empty directory (first run, no collision with the no-clobber guard). The revision-guard
  logic itself was exercised against both the real new deployment (passes) and the real old
  checkout (correctly refuses) — not reasoned about in the abstract.
- The original `/bigdata/godziklab/shared/Xinru/302006/302006_updategit` checkout was not
  modified, reset, or cleaned — confirmed still at `35131cf` after the above.

Re-deploy to a new, differently-named `DEPLOY_ROOT` (never overwrite this one) if the code changes
again before this job is submitted.

## Addendum (2026-10-10, continued): Job 1 acceptance, preparer status timing, wrapper identity

Three further problems found and fixed, each confirmed against the real code/deployment first:

1. **Job 1's acceptance criteria (below, as originally written) checked `manifest.json` fields
   (`code_revision_verified`, `code_revision_verification_method`) that belong to Step 2's
   `run_batch_diagnostics.py` only.** Checked directly: `prepare_provisional_diagnostic_input.py`
   (Step 1) never writes a `manifest.json` at all — it writes the output `.h5ad`, a
   `.selection_manifest.csv`, a `.selection_manifest.status.json` sidecar, and a `.run_owner.json`.
   Step 1's acceptance criteria below are corrected to check these actual artifacts instead.
2. **The preparer published `status="extracted"` to that status sidecar BEFORE either
   `counts.mtx` file had been streamed or the h5ad written/reloaded.** Checked directly at
   `prepare_provisional_diagnostic_input.py:315-321` (pre-fix line numbers) — a process killed or
   raising during the streamed reads or the h5ad write left an `"extracted"` sidecar pointing at
   an h5ad that was never written, or only partially written. **Fixed**: `"selected"` is published
   first (atomically); `"extracted"` is published only after streaming, the h5ad write/reload
   check, and two reconciliation checks (cohort-qualified `patient_id`, per-lineage count
   agreement) all succeed — each RAISED on, not merely recorded, if violated. Any exception in
   between is caught, published as `"failed"` with the error recorded, and re-raised. Verified
   with a genuine mutation test, not just a new assertion: the new regression test
   (`test_extraction_failure_after_selection_does_not_report_extracted`, which corrupts
   `decontx/counts.mtx`'s MatrixMarket header for real, not a mock) was first run against the
   PRE-FIX code (via `git show`) and confirmed to FAIL exactly as the defect predicts
   (`status["status"] == "extracted"` after a real, injected failure); the fixed code was then
   confirmed to pass the same test, and the full file (9/9 tests) on HPCC.
3. **The wrapper's revision guard checked `HEAD` only, not working-tree cleanliness — and this
   round's own first deployment (`302006_updategit_deploy_9ea4c6a`) tripped exactly that gap.**
   That deployment was hand-patched in place (two wrapper files `rsync`'d in for a syntax check,
   without a matching commit) while preparing this spec, leaving it dirty despite still reporting
   `HEAD=9ea4c6a...`:
   ```
   $ ssh ... 'cd .../302006_updategit_deploy_9ea4c6a && git rev-parse HEAD && git status --porcelain'
   9ea4c6a43a4ec216c345c35a463943f7b550232d
    M scripts/33_batch_effect_audit/code/slurm/step1_prepare_diagnostic_input.sh
    M scripts/33_batch_effect_audit/code/slurm/step2_run_diagnostics.sh
   ```
   **Fixed**: both wrappers now also check `git status --porcelain` is empty, refusing (`exit 1`)
   if not, before any Python runs — verified directly against a deliberately-dirtied throwaway
   copy (revision matching, tree dirty) to confirm this specific check fires, not just the
   revision check. `302006_updategit_deploy_9ea4c6a` is **superseded** (left in place, untouched,
   as the record of this exact failure mode — not reused for anything). A **second-generation**
   deployment, `302006_updategit_deploy_aea9a61` (the preparer fix + wrapper cleanliness check,
   commit `aea9a6197539ddd5086e46ddd34561617ad250d7`), was deployed the same way as the first
   (`git clone --no-hardlinks` locally, confirmed clean, `rsync`'d, re-verified clean independently
   ON HPCC) and never hand-patched afterward.

   **Wrapper identity is now explicit and separate from application identity.** The wrapper
   `.sh` files committed alongside the filled-in `DEPLOY_ROOT`/`EXPECTED_REVISION` values
   (commit `a2b173859a0630933ba8720567e626207a75ec5e` — necessarily a commit AFTER `aea9a61`,
   since a commit cannot contain its own hash as a literal string) are staged at a path separate
   from `DEPLOY_ROOT`, so submitting them never dirties the application checkout the way gen-1
   was dirtied:
   ```
   /bigdata/godziklab/shared/Xinru/302006/stage33_wrappers/a2b1738/step1_prepare_diagnostic_input.sh
   /bigdata/godziklab/shared/Xinru/302006/stage33_wrappers/a2b1738/step2_run_diagnostics.sh
   ```
   Checksums (local == remote, confirmed on HPCC):
   ```
   91f5229d477b2cc29270a0d7a2ee04d71e1af956d2a5208b011f1702339ffdc4  step1_prepare_diagnostic_input.sh
   8a515c603dcc8761f03202cbf10ba5cd5ed24375c72a438b1decfbd01ed5d814  step2_run_diagnostics.sh
   ```
   `bash -n` syntax-checked on HPCC: both OK. The guard logic itself (not reasoned about
   abstractly) was exercised against all three real checkouts: `302006_updategit_deploy_aea9a61`
   (PROCEEDs), `302006_updategit_deploy_9ea4c6a` (REFUSEs — wrong revision), and the original
   `302006_updategit` (REFUSEs — wrong revision), plus a dedicated dirty-but-right-revision
   throwaway copy (REFUSEs — specifically the cleanliness check, deleted immediately after the
   check).

   **Application identity:** `DEPLOY_ROOT=/bigdata/godziklab/shared/Xinru/302006/302006_updategit_deploy_aea9a61`,
   `EXPECTED_REVISION=aea9a6197539ddd5086e46ddd34561617ad250d7`, independently confirmed clean on
   HPCC. **Wrapper identity:** the two staged files above, source revision
   `a2b173859a0630933ba8720567e626207a75ec5e`, checksums above. **Submission commands** (not yet
   run):
   ```
   sbatch /bigdata/godziklab/shared/Xinru/302006/stage33_wrappers/a2b1738/step1_prepare_diagnostic_input.sh
   sbatch /bigdata/godziklab/shared/Xinru/302006/stage33_wrappers/a2b1738/step2_run_diagnostics.sh   # separate authorization, after Step 1
   ```

## Step 0 — what this round's testing actually covers, named exactly

- `tests/test_execution_readiness.py` (new, 17 tests): the stage-id hard-fail, the
  code-identity verified/asserted distinction (including a real git-checkout positive case and a
  non-git-checkout negative case), the atomic ownership primitive under a genuine two-thread race
  AND a genuine two-process subprocess race against the real CLI, the
  `determine_embedding_status()` pure function across all 6 status-combination cases, and the
  unknown-outcome technical/bio-eligibility split against the real CLI. **15 passed, 2 skipped**
  (one skip is environment-only — this scratch deployment is not a git checkout, so the
  git-verified-identity positive test correctly does not apply here; the other is the two-process
  CLI race test skipping itself when the HPCC login node's own resource limits killed one
  competitor by signal before the race was fairly contested — CLAUDE.md's own documented stance
  that such an exit code is "a fact about the probe, not about the code". The same test **has**
  been observed to pass cleanly this round, confirming the mechanism works when both competitors
  get to run; the lighter two-thread race test against the ownership primitive itself, which does
  not depend on login-node resource availability, passed on every run).
- `tests/test_prepare_provisional_diagnostic_input.py` (8 tests, 3 new this round): **8 passed**,
  including the cohort-qualified-`patient_id` regression test, the dry-run-keeps-manifest test,
  and the selection-manifest-consume-and-detects-drift test (built against a real, small, on-disk
  MatrixMarket fixture, not a mock).
- `tests/test_preparer_to_diagnostic_integration.py` (new, 1 test): **passed** — chains the real
  preparer CLI's output into the real runner CLI using the exact argument contract this spec's
  Step 2 (below) uses, against a real small MatrixMarket fixture. This is the specific check a
  review asked be done before presenting this spec: "Run... a small preparer-to-diagnostic
  integration fixture with the exact argument contract intended for the real job."
- The other 6 pre-existing test files (`test_batch_asw`, `test_batch_metadata_inventory`,
  `test_batch_confounding`, `test_adapt_stephenson`, `test_evaluation_boundaries`,
  `test_stable_seed` — 52 tests total, confirmed passing in prior rounds and unmodified by this
  round's changes) were not re-run together as one combined `pytest tests/` invocation this round:
  two attempts were made, and both were killed by the HPCC login node under genuinely heavy
  concurrent load that evening (208 logged-in users, load average 6–7.5 at the time, confirmed via
  `uptime`) before completing. Each file **was** run individually and completely this round where
  touched (the two files above); the four files this round did not touch at all
  (`test_batch_asw`, `test_batch_confounding`, `test_adapt_stephenson`,
  `test_evaluation_boundaries`, `test_stable_seed`, `test_batch_metadata_inventory`) were not
  re-run this round and their passing status is carried forward from the prior round's
  confirmation, not re-measured here — named explicitly as `UNVERIFIED (this round)` rather than
  silently implied to be freshly confirmed.
- **The 24,841-cell estimate and its per-cohort breakdown (carried forward from v2, table below)
  were NOT regenerated this round.** This round's changes do not alter the sampling algorithm
  itself (only add identity/ownership/selection-freeze machinery around it); the number is
  reported as v2 measured it, not re-verified against the real 72 GB input this round.

## Step 1 — input preparation (bounded, streaming)

**Command** (staged at
`/bigdata/godziklab/shared/Xinru/302006/stage33_wrappers/a2b1738/step1_prepare_diagnostic_input.sh`
— see "Addendum (2026-10-10, continued)" above for its checksum and the application deployment it
runs from; submit with `sbatch /bigdata/godziklab/shared/Xinru/302006/stage33_wrappers/a2b1738/step1_prepare_diagnostic_input.sh`):
```
DEPLOY_ROOT=/bigdata/godziklab/shared/Xinru/302006/302006_updategit_deploy_aea9a61
cd "$DEPLOY_ROOT"
# (hard-fails first if this checkout's HEAD != aea9a6197539ddd5086e46ddd34561617ad250d7,
#  or if `git status --porcelain` is non-empty)

/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean/bin/python \
    scripts/33_batch_effect_audit/code/prepare_provisional_diagnostic_input.py \
    --input-root /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/platelet_monocyte_benchmark_input_v2 \
    --out /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/batch_effect_audit/provisional_diagnostic_input_v3/input.h5ad \
    --n-per-donor 20 --seed 0
```

**Resources:** `epyc` partition, 4 CPUs, 32 GB memory, walltime ceiling 4:00:00 (unchanged from
v2's reasoning — I/O-bound, scanning 72 GB total once; the ceiling includes real margin given
this session's own prior login-node SIGKILL at well under 32 GB of genuinely necessary memory).

**Estimated output size (carried forward from v2, NOT regenerated this round):**
```
Patient-balanced cap (<=20 cells/donor/lineage, seed=0): 1,317 distinct (cohort,donor,lineage)
groups -> 24,841 cells total (11,843 Platelet + 12,998 Monocyte_Myeloid)
```

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

A reviewer may optionally request a `--dry-run` first (same command plus `--dry-run`, safe on a
login node — metadata only, no counts.mtx streaming) to regenerate this table fresh and review the
exact selection before authorizing the real extraction; the resulting
`<out>.selection_manifest.csv` (status `"preflight_only"`) can then be passed to the real
extraction via `--selection-manifest <path>` to bind it to that specific reviewed list (v2 had no
such option — a dry-run's selection was deleted immediately after being printed).

**Job 1's own acceptance criteria (strict)** — corrected to check the artifacts this script
actually writes, not Step 2's `manifest.json` (see "Addendum (2026-10-10, continued)" above):
- Exit 0.
- `<out>.h5ad` exists and reloads (`anndata.read_h5ad`) with the expected shape.
- `<out>.h5ad.selection_manifest.status.json` has `"status": "extracted"` — NOT `"selected"` and
  NOT `"failed"` (those mean the run did not complete; see the addendum above for exactly what
  distinguishes them). Its `output_sha256` matches a hash taken of `<out>.h5ad` independently.
- The same status file's `patient_id_cohort_qualified` and `per_lineage_counts_reconciled` are
  both `true` (these are RAISED on, not merely recorded, if violated — see the addendum).
- The same status file's `stage_id` is `"pre_gate2_decontx_provisional"` and `expression_source`
  is `"decontx_counts"`.
- `<out>.h5ad.selection_manifest.csv` exists (the plan this run executed) and
  `<out>.h5ad.run_owner.json` exists (this run's atomic ownership claim).
- `adata.uns["population_stage_id"] == "pre_gate2_decontx_provisional"` on reload (the hard,
  machine-checkable stage identity — distinct from the free-text `population_stage` Step 2 reads
  dynamically from this file).

**Failure criteria (stop and report; do not retry automatically):** non-zero exit, the status
sidecar absent or showing `"failed"` (read its `error` field) or `"selected"` (extraction never
completed), or `<out>.h5ad` absent despite a `"extracted"` status (should not be possible given
the fix above, but would itself be a defect worth stopping on, not retrying past).

## Step 2 — diagnostic CLI on the prepared input

**Command** (staged at
`/bigdata/godziklab/shared/Xinru/302006/stage33_wrappers/a2b1738/step2_run_diagnostics.sh` — same
checksum/deployment notes as Step 1 above; submit with `sbatch /bigdata/godziklab/shared/Xinru/302006/stage33_wrappers/a2b1738/step2_run_diagnostics.sh`,
which reads `--population-stage` directly from Step 1's actual output file rather than retyping
it, and resolves `--code-revision` from the same isolated deployment Step 1 ran from, refusing to
guess or proceed if that fails or has drifted):
```
DEPLOY_ROOT=/bigdata/godziklab/shared/Xinru/302006/302006_updategit_deploy_aea9a61
cd "$DEPLOY_ROOT"
POP_STAGE=$(python -c "import anndata; print(anndata.read_h5ad('<Step 1 output>', backed='r').uns['population_stage'])")
CODE_REVISION=$(git rev-parse HEAD)
# (hard-fails here if CODE_REVISION != aea9a6197539ddd5086e46ddd34561617ad250d7,
#  or if `git status --porcelain` is non-empty)

/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean/bin/python \
    scripts/33_batch_effect_audit/code/run_batch_diagnostics.py \
    --h5ad /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/batch_effect_audit/provisional_diagnostic_input_v3/input.h5ad \
    --out-dir /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/batch_effect_audit/provisional_diagnostic_input_v3/diagnostic_run \
    --population-stage "$POP_STAGE" \
    --stage-id pre_gate2_decontx_provisional \
    --code-revision "$CODE_REVISION" \
    --expression-source decontx_counts \
    --batch-col cohort --patient-col patient_id --lineage-col lineage --outcome-col outcome \
    --run-harmony --seed 0
```

(No `--embeddings`: this input has no FM embeddings — 5-FM full-scale extraction has not started.
Expression/PCA mode only, exercising the same code path the small-fixture integration test just
confirmed end to end with this exact argument contract.)

**Resources:** `epyc` partition, 4 CPUs, 16 GB memory, walltime ceiling 1:00:00 (unchanged
reasoning from v2).

**Acceptance criteria (strict) — this is Step 2's `manifest.json`, distinct from Step 1's own
acceptance criteria above (which check the preparer's actual artifacts, not this file):**
- Exit 0 AND `manifest.json`'s top-level `status` is exactly `"complete"` — not
  `"complete_with_partial_results"` and not `"failed"`.
- `manifest.json`'s `stage_id` is `"pre_gate2_decontx_provisional"` and no stage-id mismatch error
  occurred (this is now a hard failure if it did, not something to merely notice in the log).
- `manifest.json`'s `code_revision_verified` is `true`, with `code_revision_verification_method`
  equal to `"git_head_clean_tree"` and `code_revision` equal to
  `aea9a6197539ddd5086e46ddd34561617ad250d7` (the isolated `302006_updategit_deploy_aea9a61`
  checkout this job submits from — see "Addendum (2026-10-10, continued)" above — IS a real, clean
  git checkout, independently confirmed on HPCC; the main `302006_updategit` working copy is NOT
  used by either wrapper and remains stale/dirty at `35131cf`. The scratch deployment this round's
  tests ran against has no `.git` at all and is correctly `false`/unverified by design — a
  different case from either checkout above).
- `manifest.json`'s `input_sha256` matches a hash taken of Step 1's actual output file.
- Both lineages' `expression_pca` entries show `status: "complete"` for BOTH scopes.
- Per-lineage `n_diagnostic_cells`, `n_diagnostic_cells_bio_known`, and
  `n_diagnostic_cells_bio_unknown` are all present and reconcile
  (`bio_known + bio_unknown == n_diagnostic_cells`).
- `metrics_long.csv` has rows for both lineages × `expression_pca` × both scopes.
- Every saved `.npz` representation reloads with matching shape, cell-ID ordering, and values
  (enforced automatically — a failure here raises before the run can report success).

**Failure criteria (stop and report; do not retry automatically):** non-zero exit, `manifest.json`
absent, `status` other than `"complete"`, or any `code_revision_verified: false` for this
particular submission (unlike the scratch-deployment tests, the real job runs from a proper git
checkout and should be verifiable).

## What this job does NOT authorize

Unchanged from v2: no cell deletion, no Gate 2 scoring, no change to any existing checkpoint or
retained-cell set; no promotion of this run's numbers into the manuscript or
`docs/current-status.md` as final; no automatic retry on failure; no push; no change to
`canonical_donor_labels.csv` or any other canonical artifact. Also not required as a
prerequisite for THIS provisional diagnostic job specifically: repairing the stale/dirty HPCC
`302006_updategit` checkout (an isolated, verified deployment is what matters here, not that
specific checkout), or resolving the Stephenson/Cambridge canonical-grouping promotion decision
(a separate decision, before any FINAL dependency-grouped evaluation, not this provisional one).

## One job or two?

Unchanged from v2: two separate, dependent SLURM allocations, each needing its own per-job
authorization, requested separately (Step 2 only after Step 1's output and acceptance checks are
confirmed, not pre-authorized speculatively).
