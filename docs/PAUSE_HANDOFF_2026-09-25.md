# 302006 two-week pause handoff

**Paused:** 2026-09-25  
**Suggested pickup date:** 2026-10-09  
**Repository:** `302006_updategit`  
**Authoritative commit at pause:** `07b5662` (`main == origin/main`)  
**Purpose:** allow the project to be resumed without reconstructing the September audit history.

## One-minute summary

The project remains in the upfront-QC rebuild. Gate 1 is complete. Gate 2 has not been run in
production and has deleted no cells. A scoring-only design, all eleven Scrublet parameters, one
complete Qiu capture (`GSM5102900_HC1`), and the permit-gated pilot implementation have passed the
scoped AI methodology/code review. The reviewed implementation was committed and pushed as
`07b5662`.

The immediate next action is **not** to run scoring directly. It is to create a tracked, scoped
pilot permit in a later commit, deploy that permit commit to a new clean HPCC worktree, and run the
wrapper in dry-run mode. Only after automated review checks that dry-run may the owner authorize the one
scoring pilot.

Stage 31 is a separate workstream and must not be used as a Gate 2 blocker or shortcut. Its old
30–60-cell/sample design is not approved for scientific use.

## State frozen at pause

| Item | State |
|---|---|
| Local and remote repository | `HEAD == origin/main == 07b5662` |
| Gate 2 methodology | conditionally approved for scoring-only sensitivity analysis |
| Gate 2 parameter ruling | 11/11 parameters approved |
| Audited capture | Qiu `GSM5102900_HC1`, complete 7,666-cell deposit |
| Score-free preflight | job `29073092`, completed; PC30 feasible |
| Pilot implementation | final audit passed; committed and pushed |
| Production permit | **does not exist** |
| Scoring pilot | **not submitted** |
| Production/cohort run | **not authorized** |
| Cell deletion | **not authorized; none performed** |
| Outcome analysis | outside this pilot and not authorized |
| Gate 2 outputs from earlier runs | quarantined; must not be reused |

Two unrelated/untracked paths remained at pause and were intentionally excluded from the Gate 2
commit:

- `reports/manuscript/revision/capture_provenance_v1_2026-09-23/`
- `reports/manuscript/revision/methodology_review_2026-09-24/parameter_ruling_tests.txt`

Do not sweep them into a future `git add -A` without a separate decision.

## Reviewed pilot identity

| Binding | Reviewed value |
|---|---|
| implementation commit | `07b5662` |
| capture | `GSM5102900_HC1` |
| complete scoring input | 7,666-cell GEO deposit matrix, not the 65-cell canonical subset |
| config SHA-256 | `0f02df7c7e2a81da1efa2f371d6890878ad8f747146986ae5bc5db7a94f1dd41` |
| implementation SHA-256 | `30831273d4d3030ef5f9b123ecefbb44af2ae37aa9b7a95aab18ffeb8de80d84` |
| expected output basename | `GSM5102900_HC1__cfg0f02df7c7e2a__impl30831273d4d3` |
| prior grid | `0.005, 0.01, 0.02, 0.04, 0.08, 0.12, 0.16, 0.20` |
| interpretation | algorithmic stress settings, not prevalence estimates or bounds |

The final review is:

`reports/manuscript/revision/gate2_capture_audit_2026-09-24/FINAL_REVIEW_2026-09-25.md`

The frozen evidence manifest is under:

`reports/manuscript/revision/gate2_capture_audit_2026-09-24/final_evidence_20260925T204838Z/`

The review passed focused tests, real Scanpy integration tests, the complete repository suite,
permit satisfiability checks, and the scoped mutation checks. Do not repeat the entire design
debate unless the reviewed bytes, dependencies, input hashes, scope, or scientific question
change.

## Exact next workflow after the pause

### 1. Re-establish facts without changing anything

Run read-only checks first:

```bash
git status --short
git rev-parse HEAD
git rev-parse origin/main
git log -3 --oneline --decorate
```

On HPCC, also check `squeue`/`sacct` before claiming that no job is active. Confirm that no new
Gate 2 output directory or permit appeared during the pause. Do not infer scheduler state from an
old report.

Expected repository identity is `07b5662`. If `main` has advanced, inspect the intervening diff
before using this handoff. Do not reset or clean a dirty checkout.

### 2. Create the scoped pilot permit

The permit must be a tracked repository artifact committed **after** `07b5662`. It must reference
`07b5662` as the reviewed implementation commit and exactly bind:

- one job;
- one capture: `GSM5102900_HC1`;
- all eight approved prior settings;
- scoring only;
- no outcome access;
- no cell deletion;
- no retry and no second job;
- the reviewed config, input and implementation hashes;
- the reviewed output root;
- resource ceiling: `short`, `02:00:00`, `48G`, 4 CPUs, maximum one job.

Do not fill or activate the global `PRODUCTION_PARAMETER_RULING`; the ruling remains scoped through
the validated permit.

### 3. Deploy and dry-run

Push the permit commit, create a **new clean detached HPCC worktree** at that commit, and use the
reviewed wrapper. The implementation files must remain blob-identical to `07b5662`; the permit
gate is designed to reject implementation drift, dirty/staged files, untracked permits, wrong
scope, wrong hashes, or a reused output directory.

Run only the wrapper's dry-run first. Preserve its repository commit, implementation/config/input
hashes, permit validation, resource ceiling and derived output directory. Dry-run must return
before loading the capture or scoring.

### 4. automated review dry-run review, then explicit owner authorization

automated review should review only the permit and dry-run evidence. If they match the frozen review, no new
architecture review is needed. The owner must then explicitly authorize the **single** scoring
pilot. A dry-run approval is not itself job authorization.

### 5. Post-run acceptance

After the one job completes, inspect before interpreting:

- SLURM `State`, `ExitCode`, `Elapsed`, `MaxRSS`;
- manifest is `complete`, not `incomplete` or `stopped`;
- all sidecars exist and their hashes, row counts and cell-order identity agree;
- all eight prior settings are present;
- graph construction count, `k_adj`, shape, non-empty checksums and cross-prior topology pass;
- threshold failure, if any, remains unknown/NA rather than False/0;
- identity mismatch publishes no score sidecar and cannot produce `complete`;
- no outcome fields were accessed or exported;
- no cells were deleted and no second job was submitted.

The pilot answers whether the scoring implementation behaves correctly on one complete capture.
It does **not** estimate biological doublet prevalence, validate the automatic threshold, authorize
deletion, or authorize the other cohorts.

## What remains after a successful pilot

1. Decide whether the single-capture pilot meets its prespecified technical acceptance criteria.
2. Obtain a separate decision before expanding to another capture or cohort.
3. For every additional conditional cohort, close its own C1–C3 evidence; do not inherit Qiu's
   acceptance.
4. Only after an approved Gate 2 retained set exists: refit DecontX on that retained set.
5. Then extract platelet/monocyte benchmark cells, run the five foundation-model benchmark,
   patient-grouped CV/LOCO analyses, statistics, figures and manuscript recompilation.

No pilot result alone authorizes cell removal. A separate reviewed deletion rule would be needed.

## Separate Stage 31 track

Stage 31's approved direction is full-capture scoring before mapping back to the biological subset.
Its raw source, physical capture and raw-barcode identity remain unresolved. Existing small-sample
dimension results do not choose a PCA dimension for the redesigned workflow. Do not run Stage 31
as part of resuming this Gate 2 pilot, and do not let it re-enter the Gate 2 permit conditions.

## Files to read first on pickup

Read these in order; archived memory-bank files are not current:

1. `docs/PAUSE_HANDOFF_2026-09-25.md` (this file)
2. `docs/current-status.md`
3. `docs/decisions.md`
4. `reports/manuscript/revision/gate2_capture_audit_2026-09-24/FINAL_REVIEW_2026-09-25.md`
5. `reports/manuscript/revision/gate2_capture_audit_2026-09-24/SCORING_CONFIG_AUDIT.md`
6. `reports/manuscript/revision/methodology_review_2026-09-24/REVIEW.md`
7. `reports/manuscript/revision/methodology_review_2026-09-24/GATE2_PARAMETER_RULING.md`

## Resume checklist

- [ ] Confirm date, repository and HPCC scheduler state.
- [ ] Confirm `07b5662` or audit every later commit.
- [ ] Preserve the two intentionally excluded untracked paths.
- [ ] Confirm no permit, scoring job or output directory already exists.
- [ ] Create one tracked scoped permit referencing `07b5662`.
- [ ] Push permit commit and deploy a new clean worktree.
- [ ] Run and archive dry-run only.
- [ ] Obtain automated review dry-run approval.
- [ ] Obtain explicit owner authorization for the unique scoring pilot.
- [ ] Submit exactly one job; do not retry automatically.
- [ ] Perform post-run acceptance before any scientific interpretation.

## Stop conditions

Stop and review rather than improvising if any of these occur:

- implementation/config/input hash differs from this handoff;
- dependency versions change;
- permit validation needs a bypass;
- output directory already exists or is non-empty;
- wrapper or implementation is dirty/staged/untracked;
- graph acceptance fails;
- manifest is not complete or sidecars disagree;
- outcome data appears in memory or output;
- a second job or retry appears necessary;
- anyone proposes converting the stress-grid settings into prevalence estimates or automatic
  deletion fractions.

The correct pause state is therefore: **reviewed implementation safely frozen and pushed; no
permit and no scoring job; resume at scoped permit creation, not at method redesign.**
