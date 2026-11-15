# Batch-effect audit — findings report v1

Covers stage 33 (`scripts/33_batch_effect_audit/`), commits `e9321a8`..`462e5be`
(2026-10-08 → 2026-10-09), seven commits across five review rounds with
`codex_output/*.md` → `claude-code_output/*.md` handoffs. Not pushed.

**Scope.** This audit asked three questions about the 302006 platelet/monocyte mortality
benchmark's 8-cohort pooled dataset, ahead of the paused Gate-2 rebuild: (1) is the declared
batch unit (`cv_group`) actually independent of outcome and of physical sample pooling, (2) does
the project's batch-mixing metric measure what it claims to, and (3) can a batch-correction
diagnostic actually be run end to end against this data's real shape, cheaply enough to iterate
on before committing to the full real run. No SLURM job was submitted at any point; no real
cohort data beyond metadata/membership was read in bulk; nothing in `canonical_donor_labels.csv`
or any other already-finalized data-selection artifact was modified.

## 1. The project's batch-mixing metric was not scIB-correct

`scripts/24_option_c_canonical_input/code/embedding_quality_option_c.py` computed batch
separation as one global, unstratified `silhouette_score(X, cohort)` — a single number that
conflates "cells of different cohorts are far apart" with "cells of different *cell types* happen
to also differ by cohort," which is exactly the confound a batch metric exists to rule out.
Replaced with `scib_batch_asw()`, now in `scripts/00_shared/code/batch_asw.py` so both this file
and the new diagnostic CLI share one implementation: per-cell silhouette stratified by a
biological label, per Luecken et al. 2022 (scIB). Hardened against three distinct degenerate
cases (fewer than 2 batches, one batch per cell, a batch below a minimum group size), each
reported as a named `UNDEFINED`/`POLICY` skip reason rather than a silently wrong number. A
second, older, non-scIB-correct caller in
`scripts/12_umap_visualization/code/selected_code/batch_effect_analysis.py` was retired to raise
`NotImplementedError` pointing at the shared module, rather than left callable with the old
behavior.

## 2. 12 multi-donor `cv_group` components exist; 8 of them mix known outcomes

**Corrected 2026-10-09** (a review caught this section's heading overclaiming "12 of its
multi-donor components" mix outcome — the real split is narrower and was already stated
correctly in `MANIFEST.md`, just not here). A patient-level inventory
(`build_patient_batch_confounding_table.py`) cross-tabulating `cv_group` against outcome found
**12 multi-donor components** across 3 cohorts (COMBAT 1, Cillo 8, Liu 3) — the first pass
over-claimed "8 unsplittable blocks" by only counting components that mix outcomes, silently
omitting 4 homogeneous-outcome multi-donor components that are no less a dependency structure,
just not an outcome-mixing one. Of the 12, **8 mix known (fatal/survived) outcomes**; the other
4 are outcome-homogeneous. A `cv_group` dependency component is not automatically a single
physical batch, and outcome mixing within one is not itself proof of confounding — it establishes
overlap, which is what makes it worth checking further (as §3 does for Stephenson).

## 3. The Stephenson cohort's declared batch unit does not match its own published pooling design

Fetched the primary source directly (Stephenson et al. 2021, *Nature Medicine* 27:904–916, via
PMC8121667) rather than inferring pooling structure from the data alone: Cambridge-site samples
were physically pooled — 0.5M live cells per individual, 4 donors per pool — **before** antibody
staining and 10x loading, then demultiplexed downstream by genotype, not by a TotalSeq-C
hashtag. Cross-checked against the local canonical labels: every pooled donor's `Site` field
reads "Cambridge", consistent with the paper.

`canonical_donor_labels.csv` currently treats each of these 47 Cambridge donors as an independent
singleton `cv_group`. A union-find over the actual pool-membership graph
(`resolve_stephenson_grouping.py`) shows they collapse into **3 real physical components** (sizes
21, 22, 4) — two of which mix fatal and survived donors. If adopted, this raises the
project-wide multi-donor-component count from 12 to **15** (12 existing + all 3 new Cambridge
components, not just the 2 mixed-outcome ones — an earlier draft of `MANIFEST.md` said "12→14,"
corrected after review), and the mixed-known-outcome component count (§2) from 8 to **10** (the 2
new large Cambridge components; the 4-donor one stays single-outcome). **This proposal is not
applied** — `canonical_donor_labels.csv` is
unmodified; promotion is an explicit owner decision (see §8).

A full readiness table now exists for the cohort's complete donor list, not just the 47 already
resolved: `stephenson_full_donor_readiness_v1.csv`, **130 rows** (every Stephenson donor, none
omitted — enforced by an assertion in the script, not just reported) — **47 `proposed_pooled` +
83 `unknown_no_pool_evidence`**. The 83 are not claimed to be unpooled; they are reported as
having no pool-membership evidence in the data consulted, which is a narrower and more honest
claim.

A related, separately-flagged finding: `build_stephenson_decontx_input.py` sets
`channel_id = sample_id`, which for these 47 Cambridge donors is a finer grain than the real
ambient-RNA-sharing physical pool (4 donors demultiplexed from one physical capture should share
one ambient-RNA channel, not four). This is reported as a methodological note about an
**already-completed** DecontX run — it was not rerun, and no DecontX output was touched.

A known, documented limitation, not silently resolved: the sample-to-pool pairing check used
here is membership-only and cannot detect a donor's pools being swapped in order when that donor
is a genuine member of two pools (e.g., CV0068 is in both BGCV02 and BGCV12). Characterized by
`test_swapped_pool_order_is_not_caught_known_limitation`, not claimed fixed.

## 4. The diagnostic CLI's own cell-sampling logic silently starved the rarer lineage by roughly half

`prepare_provisional_diagnostic_input.py` originally capped sampled cells **per donor across both
lineages combined**, grouping by bare `donor_id`. A donor contributing cells to both Platelet and
Monocyte_Myeloid could lose most of the rarer lineage's cells to the more abundant one in the
random draw. Fixed to cap within each `(cohort, donor_id, lineage)` combination independently.
**Measured impact**, from a real `--dry-run` against the actual (not synthetic) 302006 pre-Gate-2
input on HPCC: total selected cells roughly doubled, **13,232 → 24,841** (11,843 Platelet +
12,998 Monocyte_Myeloid), with Platelet recovering most dramatically per-cohort — COMBAT 43→1,775;
Ren 481→3,444; Cillo 14→385. A regression test
(`test_per_lineage_cap_does_not_starve_a_rarer_lineage`) constructs exactly this scenario and
checks the fix holds. These are the real numbers now carried in
`scripts/33_batch_effect_audit/results/PROVISIONAL_JOB_SPEC_v2.md`; no SLURM job built on them has
been submitted.

## 5. Both samplers' "reproducible" seeding was not actually reproducible

`run_batch_diagnostics.py` and `prepare_provisional_diagnostic_input.py` both derived their
per-patient sub-seed via Python's built-in `hash((seed, id))`. `hash()` on a string-containing
tuple is per-process-randomized (`PYTHONHASHSEED` varies by default across Python process
invocations), so a fixed `--seed` did not reproduce the same cell selection across two separate
runs, despite both tools' docstrings claiming determinism. Fixed with a new
`scripts/00_shared/code/stable_seed.py::derive_seed()` (SHA256-based, NUL-byte-separated parts).
Proven two ways: a positive test shows `derive_seed()` is stable across `PYTHONHASHSEED` values;
a negative-control test shows plain `hash()` really does vary in this environment, so the claim
that there was a real bug to fix is not just assumed.

The synthetic fixture generator had a related, independent instance of the same class of bug: its
outcome label was originally re-derived **per `(lineage, cohort, patient)`** with the lineage
folded into the seed, so the same `patient_id` could receive a different synthetic outcome in
different lineages. Fixed to assign outcome once per `(cohort, patient)` and reuse it across
lineages/draws, with a consistency check added at fixture-build time.

## 6. A run could report `status: "complete"` while required work had silently failed

Confirmed directly: if `--run-harmony` was requested but `harmonypy` failed (or was missing), or
an explicitly-requested `--embeddings` key did not exist in the input, the enclosing embedding's
status and the top-level `manifest.json` status could both still read `"complete"`, exit 0. Fixed:
required work failing now sets `status: "failed"` and a non-zero exit; optional, non-requested
misses report `"complete_with_partial_results"`, never silently `"complete"`. The same manifest
could separately carry `code_revision: "UNKNOWN (...)"` while still claiming `"complete"` — a
`--code-revision` flag now lets a caller supply it explicitly (required on the loose-file scratch
deployment used throughout this audit, where `git rev-parse HEAD` cannot resolve), and the run now
fails outright rather than publish unverified provenance. A failure raised before the final
manifest write previously left no persisted record at all; the manifest is now also written
(`status: "failed"`) from the exception handler that re-raises.

A harmonypy-specific bug surfaced while validating this: version 0.0.10's `sigma` parameter only
broadcasts a scalar to an array when `nclust > 1`, and its own default `nclust = round(N/30)`
returns 1 for any `N < 45`, crashing `init_cluster` on small inputs. Fixed by passing `nclust`
explicitly rather than relying on the library default.

## 7. Saved per-cell representations had no identity binding; the evaluation-boundary tests didn't prove what they claimed

UMAP coordinate CSVs carried `patient_id` (which repeats across cells within a patient) but no
`cell_id`; corrected/uncorrected PCA and FM arrays were saved as bare `.npy` with no identity at
all — nothing downstream could verify a given row of coordinates belonged to the cell it was
plotted or scored against. Fixed: every representation is now saved as `.npz` bundled with its
exact `cell_id` ordering (`save_representation()`), and immediately reloaded after writing to
verify shape, identity, and values survive the round trip, not just that the file exists. Neither
CLI had a no-clobber check before this round; both now refuse to run into an existing,
non-empty output path — this guard is what turned the login-node race described in §9 into a
clean, diagnosable `FileExistsError` instead of silent file corruption.

Separately, `tests/test_evaluation_boundaries.py` (train/test leakage and group-integrity checks
against the project's real `repeated_cv.py` harness) was found not to prove what it claimed: the
leakage test perturbed random rows then checked only an untouched column, which a genuinely leaky
implementation would also have passed; the group-integrity test checked a separately-constructed
splitter, not what the harness actually fed its callback. Both rewritten to instrument
`run_repeated_kfold`/`run_loco` directly via pandas `.iloc` row-index identity: group-integrity
now asserts disjointness from the real train/test row sets across all 15 real callback invocations
(3 seeds × 5 folds) for repeated k-fold and both for LOCO; the leakage test perturbs only one
real fold's actual held-out rows and shows a correct train-only `StandardScaler` fit is exactly
unchanged (`np.testing.assert_array_equal`) while a deliberately leaky one changes — both run
through the real harness, same seed, so the fold partition itself is shown not to depend on `X`'s
values.

A new `tests/test_prepare_provisional_diagnostic_input.py` builds a small, real, on-disk
MatrixMarket fixture (`scipy.io.mmwrite`, same genes×cells orientation as the real files) and runs
the actual preparer CLI as a subprocess end to end — streaming read via the real production
`read_mtx_columns()`, h5ad write, reload, and a value-correctness check (decontx built as exactly
0.7× raw, verified post-roundtrip). All 5 tests passed on first real run.

## 8. Synthetic validation: the diagnostic CLI does what it claims, end to end

A dedicated synthetic fixture (`make_synthetic_batch_fixture.py`) — 2 lineages, 4 cohorts, a known
injected batch shift, a known biological signal, repeated/seeded patient draws, a synthetic
FM-like embedding `X_synthFM` — was run through the full diagnostic CLI four times across this
audit (fixture v1→v3 as bugs above were found and fixed; a clean v4 diagnostic rerun as the final
step of this round). The latest clean run (`synthetic_fixture_v4/diagnostic_run`,
`--code-revision e93cc6a-plus-plots`) reports strict `status: "complete"` for both lineages × both
embeddings (not `"complete_with_partial_results"`), `code_revision_verified: true`.

Quantitatively: injected batch effect detected via `mean_neighbor_batch_mixing` (~0.08–0.20
pre-correction vs. ~0.75 expected under random mixing) and reduced by real Harmony correction
(~0.48–0.69) across all 8 lineage×embedding×scope combinations. Biological-conservation
(KMeans + ARI against the known synthetic signal) is often near zero pre-correction — the
synthetic batch effect dominates naive clustering — and rises sharply post-correction. This is
stated as a property of that specific metric on this specific synthetic design, not generalized
into "correction always preserves biology."

**Corrected 2026-10-09** (a review inspected all 4 rendered PNGs directly and found this
paragraph's blanket description too strong for 2 of the 4; the wording below reflects what is
actually visible in each, not a generalization from the cleanest case). Rendered before/after
UMAP plots (§9) show a clean textbook pattern for the two `X_synthFM` panels (Platelet and
Monocyte_Myeloid): 4 cohort-pure clusters collapse into thoroughly intermixed cohorts
post-Harmony, and the outcome signal — not visibly structured pre-correction, where cohort
dominates — becomes visually separable into the two remaining point-clouds post-correction. The
two `expression_pca` panels tell a messier, equally real story: cohort separation is already
imperfect pre-correction for 3 of the 4 synthetic cohorts in the Platelet panel (they partially
overlap even before any correction is applied), the outcome signal is visibly structured
pre-correction in both `expression_pca` panels (not absent, unlike the `X_synthFM` panels), and
post-correction BOTH `expression_pca` panels retain real residual cohort structure — Harmony
visibly reduces but does not eliminate it. These are measured differences between the 4
rendered panels, viewable at
`scripts/33_batch_effect_audit/results/synthetic_fixture_v4/plots/*.png`, not a single
"before vs. after" story generalized across all embeddings.

## 9. Before/after plots and the full Stephenson readiness table (this round)

Closed the two items explicitly left deferred at the end of the prior round:
`run_batch_diagnostics.py` previously saved UMAP coordinates only for the uncorrected scope, so no
before/after comparison was possible from its own output even though Harmony correction was
already being run. Fixed to save both scopes; a new `make_batch_effect_plots.py` renders the
4-panel before/after grids described above, validated against a clean rerun and visually confirmed
by direct inspection of the output PNGs (not committed — repo-wide `*.png` gitignore rule,
respected as existing policy; the images exist locally and on HPCC scratch, reproducible on
demand). The full 130-donor Stephenson readiness table from §3 was also completed and pulled back
this round. Full detail is in `claude-code_output/20261009_181353.md`, including a login-node
process race encountered while regenerating the diagnostic run — at the time attributed entirely
to an orphaned process (an operational/diagnosis issue, not a code defect). **Corrected by §10**:
the orphaned process explains why two writers existed that day, but the output-directory guard
they raced against was itself a real, fixable TOCTOU (check-then-act) defect, independent of how
the race started — fixed this round, not merely diagnosed.

## 10. Execution-readiness repair (2026-10-09/10, second review round)

A fresh review (`codex_output/20261009_230728.md`) inspected the provisional-job packet §9
produced and found it was not yet safe to execute as documented, independent of the scientific
audit above: following `PROVISIONAL_JOB_SPEC_v2.md`'s own two commands verbatim would have
triggered a stage-label mismatch its own acceptance criteria rejected; `code_revision_verified`
was `code_revision is not None`, so any hand-typed string — including v2's own
`e93cc6a` — reported as "verified" without being checked against anything; both CLIs' output
guards were check-then-act, not atomic, which is confirmed to be the real underlying cause of
this session's earlier login-node process race (§9 correction above); unknown-outcome cells were
withheld from the entire technical-integration population (PCA, Harmony, neighbor-mixing), not
just from the bio-labeled metrics that actually need a real outcome; and a failed corrected-UMAP
did not propagate to `status`, asymmetric with the already-handled uncorrected case.

All five were fixed, each with a real regression test, not merely documented:

1. **Stage identity**: new `scripts/00_shared/code/population_stage.py` gives both CLIs a small,
   fixed, machine-checkable `--stage-id`/`uns["population_stage_id"]` vocabulary, matched exactly
   and HARD-FAILING on mismatch — unlike the pre-existing free-text `--population-stage` check,
   kept unchanged as a descriptive-only warning.
2. **Code identity**: `_resolve_code_identity()` replaces the `is not None` check with a real one
   — verified now means a real git checkout, resolved `HEAD`, a clean working tree, and a match
   to the asserted value; `source_file_sha256` is always recorded as the honest substitute when
   none of that holds (e.g. on this session's loose-file HPCC scratch deployment, which is
   correctly *never* "verified" under the new check).
3. **Atomic run ownership**: new `scripts/00_shared/code/run_ownership.py::claim_ownership()`
   closes the TOCTOU race with an OS-atomic exclusive create, proven with a real two-thread race
   against the primitive itself and a real two-process race against the actual CLI (both tests
   pass; the CLI-level one is environment-sensitive to HPCC login-node resource contention and
   skips itself honestly rather than reporting a false result when that happens — see §10 testing
   note below).
4. **Eligibility decoupling**: `run_one_embedding()` now takes an explicit `bio_known_mask`;
   `neighborhood_batch_mixing` runs on every sampled cell regardless of outcome, while
   `batch_asw`/`bio_conservation` run only on the known-outcome subset — both cell counts
   (`n_cells_technical`, `n_cells_bio_known`) reported separately, not collapsed.
5. **Status propagation**: `determine_embedding_status()`, extracted as a pure function
   specifically so the corrected-UMAP-failure case could be regression-tested deterministically,
   now covers it symmetrically with the uncorrected case.

A sixth, lower-severity finding (no durable, reviewable selection artifact — a `--dry-run`
deleted its own selection manifest) was also fixed: the manifest is now kept with a
`preflight_only`/`extracted` status sidecar, and a new `--selection-manifest <path>` lets a real
extraction consume and re-validate a previously-reviewed selection instead of only recomputing
deterministically. A seventh, structural identity hazard (`patient_id = donor_id`, not
cohort-qualified, risking a cross-cohort collision) was fixed to `f"{cohort}::{donor_id}"`.

**Testing.** `tests/test_execution_readiness.py` (new, 17 tests): 15 passed, 2 skipped — one
inapplicable-by-design (this scratch deployment is not a git checkout), one an honest
environment-only skip when the HPCC login node killed a race competitor by signal before the race
was fairly contested (the same test has also been observed to pass cleanly this round; the
lighter ownership-primitive race test, independent of login-node load, passed on every run).
`tests/test_prepare_provisional_diagnostic_input.py` (8 tests, 3 new): all 8 passed.
`tests/test_preparer_to_diagnostic_integration.py` (new): a real end-to-end run chaining the
preparer's actual output into the runner using the exact argument contract
`PROVISIONAL_JOB_SPEC_v3.md` documents, against a real small MatrixMarket fixture — passed. The
six pre-existing, untouched-this-round test files (52 tests, confirmed passing in the prior
round) were not re-run together this round: two attempts to run the full `tests/` directory in
one combined invocation were killed by genuinely heavy login-node load that evening (208 users,
load average 6–7.5 via `uptime`) before completing, and are marked `UNVERIFIED (this round)`
rather than silently assumed still-passing.

Full detail, including every finding's exact code citation and the fix-by-fix reasoning, is in
`claude-code_output/20261010_004454.md`. Replaces `PROVISIONAL_JOB_SPEC_v2.md` with
`PROVISIONAL_JOB_SPEC_v3.md` (v1/v2 kept for the audit trail) and adds real, reviewable sbatch
scripts (`scripts/33_batch_effect_audit/code/slurm/step{1,2}_*.sh`) with exact interpreter paths
— neither submitted.

## What needs an owner decision next

1. **Promote or decline the Stephenson/Cambridge regrouping** (§3) — would change
   `canonical_donor_labels.csv`'s `cv_group` assignment for 47 donors and raise the multi-donor
   component count from 12 to 15 (mixed-known-outcome: 8 to 10). Not applied.
2. **Authorize or decline the two-step provisional diagnostic job** in
   `PROVISIONAL_JOB_SPEC_v3.md` (§10) — real ~24,841-cell input (carried forward from v2, not
   regenerated this round), Gate-2-pre data, two SLURM allocations with real sbatch scripts
   prepared, each needing its own per-job authorization. No SLURM job has been submitted at any
   point in this audit.
3. **The stale/dirty HPCC `302006_updategit` checkout** (commit `35131cf`, dozens behind local
   HEAD, uncommitted loose Gate-2 files at its root) — left untouched throughout; the Gate-2 pilot
   dry-run per `docs/PAUSE_HANDOFF_2026-09-25.md` has not been confirmed to have happened. Not a
   prerequisite for the provisional job above (an isolated, verified deployment is what matters
   there), but still open.

## Commits

`e9321a8` → `7285e56` → `e28710e` → `9c7767f` → `e93cc6a` → `02eb850` → `462e5be` → `62241ba` →
(this round's commit), all on `main`, none pushed. Round-by-round handoffs:
`claude-code_output/20261008_201608.md`, `20261009_003147.md`, `20261009_110215.md`,
`20261009_131641.md`, `20261009_143900.md`, `20261009_181353.md`, `20261010_004218.md`.
