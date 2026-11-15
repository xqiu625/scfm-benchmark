# Gate 2: adjudication packet for the eleven Scrublet runtime parameters

**For a methodological reviewer. This packet proposes no values and approves nothing.**

**Reviewer update, 2026-09-24:** all eleven choices are now adjudicated in
[`GATE2_PARAMETER_RULING.md`](../methodology_review_2026-09-24/GATE2_PARAMETER_RULING.md).
This packet remains the supporting request; its previous open-parameter status is historical.
The production slot remains inactive pending capture C1–C3, config audit and pilot authorization.

The 2026-09-24 ruling (`reports/manuscript/revision/methodology_review_2026-09-24/REVIEW.md`)
approved the eight-point prior grid without automatically approving the other eleven parameters.
The follow-up decision linked above now settles them. `PRODUCTION_PARAMETER_RULING` remains `{}`
as an inactive deployment slot; capture eligibility, config audit and pilot authorization remain.

**Filling a default here would be the defect the fail-closed design exists to prevent.** The
scanpy defaults are listed for orientation only: recording a default is provenance, not approval.
That distinction has already been violated once in this project — Stage 31 ran with an implicit
`n_prin_comps` because nothing forced the choice to be made.

## What a ruling must contain

For each parameter, per the ruling's own C4 ("parameters, software version, randomness, threshold
rule and effective-dimension condition fixed in advance"):

1. **value** — or the explicit option **"insufficient evidence, do not proceed"**;
2. **scope** — which cohorts / input stages it applies to, and whether one value serves all;
3. **basis** — what the value rests on (measurement here, published method, software default
   accepted deliberately, or judgement), and its evidence grade;
4. **invalidation** — what observation would retrospectively make the choice wrong.

A ruling without (4) is not usable: it cannot be falsified later.

## Environment this would apply to

`scanpy 1.11.3`, the only version in `APPROVED_SCANPY_VERSIONS`. Defaults measured from the
installed signature on 2026-09-24, not read from documentation.

## The eleven

Grouped by the interface each one reaches. The runner already separates these, because a
simulation parameter passed to the scoring call is silently ignored when `adata_sim` is supplied.

### Simulation (2) — `scrublet_simulate_doublets`

| parameter | scanpy default | what it controls | what makes it hard |
|---|---|---|---|
| `sim_doublet_ratio` | `2.0` | how many synthetic doublets per observed cell | Sets the size of the simulated set and thus `λ` in the score identity. Larger is not free: it changes the neighbour composition the score is read from. |
| `synthetic_doublet_umi_subsampling` | `1.0` | UMI subsampling when combining two cells | At 1.0 a synthetic doublet has the summed depth of its parents, which no real doublet need have. Whether to correct depth is a modelling choice, not a default. |

**Cross-cutting caveat the ruling already recorded:** for the upstream-filtered cohorts, pairing
surviving cells simulates cross-source combinations that upstream filtering already removed. No
value of these two parameters repairs that category mismatch.

### Preprocessing (1)

| parameter | scanpy default | what it controls | what makes it hard |
|---|---|---|---|
| `log_transform` | `False` | log1p before PCA | scanpy applies this **only** on the `adata_sim is None` branch. This path supplies `adata_sim`, so the runner applies it itself — the value must be ruled for the runner, and passing it to the scoring call would be silently ignored. |

### Scoring (8) — `sc.pp.scrublet`

| parameter | scanpy default | what it controls | what makes it hard |
|---|---|---|---|
| `n_prin_comps` | `30` | PCA dimensionality for the neighbour graph | Gate 2 requires its own deliberate choice and effective-dimension check. `docs/scrublet_dimension_parameter.md` concerns Stage 31's different population; its deferred decision does not automatically defer Gate 2. A fixed reference dimension for a scoring sensitivity experiment is not a biologically calibrated optimum. |
| `n_neighbors` | `None` → `round(0.5·√n_obs)` | kNN size for the score | Data-dependent by default, so it varies per capture unless fixed. Fixing it makes captures comparable in method but not in meaning; leaving it adaptive makes the neighbour graph differ across units. Either way the choice must be stated. |
| `stdev_doublet_rate` | `0.02` | prior-uncertainty contribution to score errors/Z | Tagged 1.11.3 core uses this in the error expression, not the point score or score-threshold call. It does not widen the score grid. The fixed-scenario interpretation and use of errors/Z need an explicit ruling. |
| `knn_dist_metric` | `'euclidean'` | neighbour distance metric | Changes the graph the score is computed on. Interacts with `normalize_variance` and `mean_center`. |
| `normalize_variance` | `True` | scale using observed gene variance before PCA | The Gate 2 input is a complete declared capture stage, not a platelet-only subset. Scaling must be shared between observed and simulated cells; zero/non-finite variance and downstream geometry need checking. |
| `mean_center` | `True` | centre before PCA | Centring a sparse count matrix densifies it and changes the leading components. |
| `use_approx_neighbors` | `None` | approximate vs exact kNN | An approximation changes results in a way that is not reproducible across library versions. Exact is defensible; the point is that it must be chosen, not inherited. |
| `get_doublet_neighbor_parents` | `False` | return parent identities of doublet neighbours | Diagnostic only. It does not change scores, but it does change what evidence is available afterward — worth ruling deliberately rather than leaving off. |

## What is already fixed and needs no ruling

| | |
|---|---|
| `expected_doublet_rate` | the approved eight-point grid, label `assumed_algorithmic_stress_setting` |
| `random_state` | `capture_seed(cohort, channel_id)` — fixed per capture, deliberately **not** a function of the rate, so the grid is a calibration sweep over common random numbers |
| `adata_sim` | one simulated set per capture, reused across the grid |
| `copy` | `True` |

## Interactions a reviewer should not have to rediscover

- `n_prin_comps`, `n_neighbors`, `knn_dist_metric`, `normalize_variance` and `mean_center` jointly
  determine one object: the neighbour graph. They cannot sensibly be ruled one at a time.
- With that graph fixed, moving the prior applies a monotone transform to the scores, so rank
  order and rank-based AUC are invariant across the grid. The runner treats that as an
  **implementation check only** — it is not evidence of robustness.
- The automatic histogram threshold is *not* invariant and may succeed at some grid points and
  fail at others. That is a result to report, not a failure to fix.

## Sequencing

Even a complete ruling does not make anything runnable. The per-cohort C1–C3 input preconditions
are unestablished for all six conditional cohorts, and `emit_runnable_config()` refuses on that
ground separately. Both must be satisfied.

## What this packet deliberately omits

No recommended values, no ranges, no "reasonable defaults", no ordering by likely acceptability.
Supplying any of those would convert an adjudication request into a rubber stamp.
