# Gate 2 single-capture preflight — job package (NOT SUBMITTED)

**No SLURM submission is authorised. Nothing in this package has been run.**

## Exact command

```bash
cd /bigdata/godziklab/shared/Xinru/302006/worktrees/<commit>
mkdir -p logs        # REQUIRED: SLURM opens the relative --output path at submission
sbatch reports/manuscript/revision/gate2_capture_audit_2026-09-24/run_capture_preflight.sh
```

The wrapper resolves the repo from `SLURM_SUBMIT_DIR` and refuses to run unless **the script, the
wrapper itself, and every module it imports** (`gate2_parameter_ruling`, `gate2_stress_design`,
`gate2_stress_runner`) are tracked and unmodified. Validating only the script would let a changed
ruling reach the run unrecorded. It prints the commit and a SHA-256 per file, and writes to a
**new** versioned directory keyed by commit and timestamp, refusing if it exists.

## Output schema and manifest lifecycle

`CAPTURE_PREFLIGHT_MANIFEST.json` is the only output, written atomically at three points:

```
status=incomplete, stage=starting       before any heavy work
status=incomplete, stage=preprocessed   after filtering + HVG
status=complete                         only after every check passes
```

A kill during loading, HVG or PCA leaves an `incomplete` manifest naming the stage reached.

```
preflight, capture, cohort, channel_id, generated_utc, status, stage
produced_scores_calls_or_rates : false
outcome_columns_read           : []
production_parameter_ruling_written : false
provenance   { repo, commit, working_tree_clean,
               file_sha256 { preflight_script, wrapper, gate2_parameter_ruling,
                             gate2_stress_design, gate2_stress_runner } }
ruling       { version, reference, imported_not_copied, module }
inputs       { 3 deposit files + container: path, bytes, sha256, asserted }
dependency_identity, approved_scanpy_versions
seed         { value, derivation: "capture_seed(cohort, channel_id)" }
input_summary { shape_as_loaded, n_cells_effective, n_usable_features_hvg, n_simulated }
feasibility  { rule, n_prin_comps, limit, feasible }
parameters   { the eleven, n_neighbors resolved from N_eff }
k_adj        { value, rule }
prior_grid, prior_label, versions
manifold     { built_by, observed_shape, simulated_shape,
               observed_sha256, simulated_sha256, all_finite,
               rank_supports_n_prin_comps }
simulated_matrix_checksum
prior_independence { expected_doublet_rate_consumed_by, manifold_built_once,
                     evidence, graph_checksum: null, graph_checksum_reason }
```

## What it measures — and the limitation that must survive into the config audit

Measures: post-filter `N_eff`; `G_HVG`; the strict rule `n_prin_comps < min(N_eff, G_HVG)` for the
approved PC=30; whether the manifold is finite and has `n_prin_comps` columns; and SHA-256 of the
observed and simulated manifolds built by **scanpy's own** observed-fit pipeline.

**No neighbour graph is checked, and none is claimed.** The kNN graph is built inside
`calculate_doublet_scores` — it is what consumes `k_adj` — so checksumming it would cross the
score-free boundary. The manifest records `graph_checksum: null` with that reason and stores the
`k_adj` scanpy would use. **This limitation must be restated in the later config audit**: the graph
stays unverified until a scoring run is separately authorised.

**Prior independence is structural, not empirical.** `expected_doublet_rate` is consumed only by
`calculate_doublet_scores`, so every step up to the manifold is independent of it and the manifold
is built **once**. An earlier draft looped over eight priors with the loop variable unused, which
showed only that the script reruns deterministically; a test now asserts no such loop exists.

## Resources and the estimate behind them

| | |
|---|---|
| partition | `short` |
| walltime | `2:00:00` |
| memory | `48G` |
| cpus | 4 |

**48G is kept, and the earlier reduction to 24G was based on a false claim of mine.** I wrote that
scanpy's `zscore` branch "keeps these sparse". Measured from the pinned source
(`_scrublet/sparse_utils.py`):

```python
def sparse_zscore(E, *, gene_mean=None, gene_stdev=None):
    return sparse_multiply(np.asarray((E - gene_mean).T), 1 / gene_stdev).T
```

`np.asarray((E - gene_mean).T)` **densifies explicitly**, and `pipeline.mean_center` likewise
builds `sparse.csc_matrix(self._counts_obs_norm - gene_means)` from a dense intermediate. With the
approved `mean_center=True`, the path densifies — exactly as the ruling warns.

Bound: observed `7,666 × G_HVG` plus simulated `2 × 7,666 = 15,332 × G_HVG`, i.e. `22,998 × G_HVG`
float64. At `G_HVG ≈ 2,000` one copy is ~0.37 GB; at `G_HVG ≈ 5,000` it is ~0.92 GB, and the
subtract/transpose/multiply chain holds several at once. Add the `mmread` of a 9.9M-non-zero matrix
and the 33,538 × 7,666 CSR before subsetting. `G_HVG` is the quantity this job exists to measure,
so the bound cannot be tightened in advance.

`MaxRSS` from `sacct` should calibrate any later reduction. Until then the first run keeps headroom.

## Boundary

Preparing this package authorises nothing. Submission requires a separate explicit authorisation.
