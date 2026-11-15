# Gate 2 scoring config audit — `GSM5102900_HC1` (AUDIT PASSED, NO PILOT PERMIT)

**No scoring job is authorised. The config is non-runnable by construction and the wrapper refuses
without a scoped pilot permit that does not exist.**

## What closed, and what did not

| condition | verdict | evidence |
|---|---|---|
| C1 un-normalised raw counts, source stated | **PASS** | MatrixMarket `integer` banner; 0 non-integer of 9,903,096; `feature_types = {Gene Expression}`; GEO `cellranger count` |
| C2 physical capture + raw-barcode namespace | **PASS** | 7,666 barcodes, all unique, all 16 nt, GEM suffix `-1` throughout |
| C3 complete declared-stage population | **PASS (deposit) / FAIL (container)** | deposit 7,666 cells; container keeps 65 = **0.85%**, join 65/65 clean |
| PC=30 strict feasibility | **PASS** | `N_eff=7,666`, `G_HVG=1,576`, 30 < 1,576 (job 29073092) |
| **kNN graph** | **UNVERIFIED** | built inside `calculate_doublet_scores`; unreachable from a score-free preflight |

`graph_preflight_verified = false` is carried in the config with its reason, and
`assert_permitted()` refuses if the reason is ever emptied — a bare `false` would read as an
oversight rather than a finding.

## What the config binds

| | |
|---|---|
| config version | `gate2-scoring-config-draft/GSM5102900_HC1/1` |
| config SHA-256 | `0f02df7c7e2a81da1efa2f371d6890878ad8f747146986ae5bc5db7a94f1dd41` |
| implementation digest | `30831273d4d3030ef5f9b123ecefbb44af2ae37aa9b7a95aab18ffeb8de80d84` |
| output directory | `GSM5102900_HC1__cfg0f02df7c7e2a__impl30831273d4d3` |
| capture | `GSM5102900_HC1` (Qiu), one complete capture |
| input | the three deposit files by full SHA-256; population is the **7,666-cell deposit matrix, not the 65-cell container subset** |
| score-free preflight commit | `9aab7befba0bd2a59ea1224aff95d261237611c6`; scoring implementation is bound separately by the final review and permit |
| preflight | job `29073092`, manifest SHA-256 `a78e5a01b4ff05f8…` |
| dependencies | scanpy 1.11.3, anndata 0.11.4, numpy 2.2.6, scipy 1.15.3, sklearn 1.7.1, python 3.10.18 |
| parameters | the eleven, **imported** from `gate2_parameter_ruling`, never copied |
| resolved | `n_neighbors = 44`, `k_adj = 132` — measured by the preflight on this capture |
| priors | the eight approved values, label `assumed_algorithmic_stress_setting` |
| seed | `10,137,777` from `capture_seed('qiu','GSM5102900')` |
| threshold | scanpy automatic; failure is unknown/NA, never 0/False/singlet |
| output | new directory required, reuse forbidden |

The config hash covers every field above, so **any** drift — a parameter, an input hash, the
capture — produces a different hash and invalidates a permit issued against the old one.

## The permit gate

`assert_permitted()` requires a JSON artifact at `$GATE2_PILOT_PERMIT` that is **inside the
deployed repository, tracked, and identical to HEAD**, carrying `permit_id`, `authorised_by`,
`authorised_utc`, `capture`, `config_sha256`, `input_sha256`, `implementation_sha256`,
`implementation_commit`, `output_root`, `resource_ceiling` and `scope`.

**The environment bypass was removed.** An earlier version honoured
`GATE2_PERMIT_TEST_BYPASS=1` in production code and called it test-only. It was not: anyone able
to set the job environment could switch the authenticity check off. Tests now inject a verifier
through `assert_permitted(..., repo_verifier=...)`, which the production call never passes, and
both the config and the wrapper **refuse outright if that variable is present at all** — presence,
not truth, so `=0` is refused too.

A docstring inside the permit check still read *"Tests set `GATE2_PERMIT_TEST_BYPASS` explicitly;
production cannot"* after the bypass had been deleted. Nothing executed it, but it documented a
mechanism that no longer exists, which is how a removed bypass gets reintroduced by someone
trusting the comment. It has been corrected, and the refusal is now covered by a test that sets
the variable alongside an **otherwise-valid** permit — a test using a broken permit would have
passed whether or not the branch honoured the flag.

`scope` and `resource_ceiling` are matched against exact values, not merely checked for presence:

```
scope            one job, one capture, eight approved priors, scoring only;
                 no outcome access, no cell deletion, no retry, no second job
resource_ceiling partition=short, time=2:00:00, mem=48G, cpus_per_task=4, max_jobs=1
```

`implementation_sha256` digests the **seven** files that are the implementation — entrypoint,
config, wrapper, `gate2_parameter_ruling`, `gate2_stress_design`, `gate2_stress_runner` and
**`gate2_diagnostics`**, which the runner imports and which an earlier list omitted, so a changed
diagnostic module could have altered output while the permit stayed valid. The permit additionally
names `implementation_commit`, which must be an ancestor of deployed `HEAD`, with all seven
implementation blobs unchanged and clean on disk. The permit is committed after the reviewed
implementation; see the satisfiability check below.

**The wrapper was guarding a smaller set than the permit digests.** `IMPLEMENTATION_FILES` carried
seven files; the wrapper's shell list carried six, omitting `gate2_diagnostics.py` — the same
transitive module the digest fix had just added. The permit gate would still have caught an edit
to it (the digest covers it), so this was not an open hole; what was missing is the wrapper's
tracked-and-clean check, which never looked at that file. The list is completed, and the wrapper
now **asserts its guarded set equals `IMPLEMENTATION_FILES`** and exits 4 on any difference,
rather than relying on two hand-maintained lists staying in step.

## Output destination

```
OUTPUT_ROOT/<capture>__cfg<config sha[:12]>__impl<implementation digest[:12]>
OUTPUT_ROOT = /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/gate2_capture_scoring
```

Any config or implementation change lands elsewhere, and a non-empty directory is refused.

**That derivation was false when first written, and is now true.** `config_sha256()` hashed the
full config *including* `output.directory`, while `output_directory()` had to hash the body alone
to avoid recursing through itself. The two produced different values: the permit bound
`8c7ac994…` while the directory displayed `fad6abb0264b`, a hash no published field equalled, so
an auditor holding a permit could not reproduce the output path. `config_sha256` now excludes
`output.directory` — the field derived *from* it — leaving one hash for both call sites.

The test that should have caught this asserted only that the substring `__cfg` appeared in the
path. It could not fail for this defect. It now asserts the segment **is** `config_sha256()[:12]`.

## Scoring entrypoint and its outputs

Calls `gate2_stress_runner.score_capture_across_grid` — the production path — rather than
reimplementing preprocessing. Loads only the three deposit files and refuses if `obs` carries
anything beyond cell identity.

**Per-cell scores are preserved.** The runner's return contract now carries a third sidecar of
per-prior observed doublet scores, retained **even when the threshold fails** — which is exactly
when the distribution is the remaining evidence. Scores are never inferred from calls; that is not
possible, and was not attempted.

Outputs, each with schema version, SHA-256, byte count and cell-order checksum:

```
SCORING_MANIFEST.json     gate2-scoring-manifest/1    rows, stability, graph, outputs
PER_PRIOR_CALLS.json      gate2-per-prior-calls/1     booleans in cell order
PER_PRIOR_SCORES.json     gate2-per-prior-scores/1    observed scores in cell order
CELL_ORDER.json           gate2-cell-order/1          the cell identifiers themselves
```

**Completion is recorded last.** Every sidecar is written to a `.inprogress` temp, validated for
length against the cell count and re-read byte-for-byte, and only then are all of them promoted;
the manifest is finalized `complete` after that. A previous version wrote `status=complete` first
and then a sidecar, so a failure on that second write left a manifest asserting a completion that
had not happened. On failure the manifest stays `incomplete`/`stopped` with the last stage reached.

**The atomicity claim has been narrowed.** `_promote_sidecars()` said "a failure part-way leaves no
half-published output set". That is false: promotion is a sequence of `os.replace` calls, not a
directory-level atomic operation, so a process kill between two replaces leaves a partial set. What
is true is narrower and is what the docstring now says — no sidecar appears under its final name
until *every* payload has been written, fsynced and re-read, so a **validation** failure publishes
nothing; and a partial set is always accompanied by a manifest that does not say `complete`. The
manifest is the commit marker; the sidecars are not self-certifying.

`json.dumps(..., default=str)` was also removed. Stringifying an unexpected object would have
serialized something other than what was computed, and passed.

## The permit/HEAD deadlock

Requiring `deployed HEAD == permit.implementation_commit` was a deadlock, not a safeguard. A
production permit must be tracked and identical to HEAD, so it can only arrive in a commit made
**after** the implementation commit — at which point HEAD is the permit commit, and a permit
cannot contain the hash of the commit that contains it. **No legitimate permit could have passed.**

The gate is now:

| requirement | why |
|---|---|
| `implementation_commit` is a commit in this repository | a permit naming nothing authorises nothing |
| it is an **ancestor** of deployed HEAD | the deployment must descend from what was reviewed |
| every `IMPLEMENTATION_FILES` path is **byte-identical** (git blob id) between that commit and HEAD | ancestry alone permits a later commit to have edited the implementation |
| those files are also clean against HEAD **on disk** | `implementation_sha256()` digests the working tree, so without this the digest could be satisfied by a local edit |
| the permit itself stays tracked and identical to HEAD | unchanged |

A two-commit fixture covers it behaviourally: commit A holds the implementation, commit B adds the
permit naming A, and B must pass — while an implementation file edited in B must fail.

## Graph evidence — the capture was attached to the wrong object

**The previous graph instrumentation captured nothing, and reported success.** It patched
`sklearn.neighbors.NearestNeighbors.kneighbors`. Measured on the pinned stack, that records
**zero** constructions, for two independent reasons:

- scanpy 1.11.3 does not reach `NearestNeighbors` here. `Scrublet._nearest_neighbor_classifier`
  goes through `scanpy.Neighbors.compute_neighbors`, whose backend is a `KNeighborsTransformer`.
- `NearestNeighbors.kneighbors` **is** the inherited `KNeighborsMixin.kneighbors`
  (`is` comparison: True), so assigning to the subclass shadows nothing the transformer resolves.

`pynndescent` is also installed in this environment, so the backend is not even guaranteed to be
sklearn. Any capture point below scanpy is contingent on a backend choice.

The capture point is now the value the doublet score **consumes** —
`_get_indices_distances_from_sparse_matrix(knn.distances, k_adj)` in
`scanpy/preprocessing/_scrublet/core.py:353` — which is backend-independent. Values pass through
untouched.

Measured on the pinned scanpy (300 observed / 600 simulated, `n_prin_comps=10`), and re-derived by
a `scanpy_integration` test so that a scanpy change breaks the test rather than silently changing
production:

| | measured |
|---|---|
| constructions per prior | **1** |
| `n_neighbors` argument | `k_adj` |
| returned index columns | **`k_adj − 1`** — *not* `k_adj`, despite the helper's docstring |
| first column is self | **False** — likewise contrary to the docstring |

Asserting self-inclusion or `k_adj` columns would have failed a **correct** run. For this capture
the expected shape is `(22,998, 131)` = `(7,666 + 15,332, 132 − 1)`.

### What is now enforced

Zero, missing, extra, wrong-`k`, wrong-shape, empty-checksum or topology-inconsistent captures all
write `status=stopped`, stage `graph_unacceptable`, and exit nonzero. The branch is on a **problem
list**, never on a tri-state field whose `None` means "nothing was captured" — that tri-state was
the defect: with zero captures it was `None` and with one capture it was `True`, so `is False`
never fired. `k_matches_expected` was computed and never read; the count, `k`, shape and checksum
emptiness are now all gates.

**Topology is enforced; distances are recorded only.** Distance invariance across priors is not
claimed, so it is not made a gate — `distances_enforced: false` states this in the record rather
than leaving a `*_verified` field that nothing branches on.

## Graph evidence — still not a preflight claim

`_GraphCapture` wraps the consumed-graph helper
`_get_indices_distances_from_sparse_matrix(knn.distances, k_adj)` in
`scanpy/preprocessing/_scrublet/core.py` for the duration of the run — **backend-independent**,
because it sits above whichever neighbour backend scanpy selects — recording the indices and
distances **on their way back to scanpy, passed through untouched**, so the graph recorded is the
graph the doublet score consumed. Per construction it records shape, the `n_neighbors` argument,
column count, whether the first column is self, and SHA-256 of both indices and distances; across
the eight priors it checks topology identity and compares the argument against `k_adj = 132`.

Zero, missing, extra, wrong-`k`, wrong-shape, empty-checksum or topology-inconsistent captures
**stop the run**. The record carries
`preflight_verified: false` with the reason — the score-free preflight could not construct this
graph, and graph verification remains a **post-run pilot acceptance item**, not something already
established.

## Scope

One capture. Nothing here extends to the other five conditional cohorts, whose C1–C3 remain
unestablished. The earlier 0-of-6 draft predates this capture's acceptance. Scoring the **deposit** rather than the canonical container is a condition of
the C3 pass, not a detail.

## Final-review identity correction (2026-09-25)

Identity-mismatched score vectors are no longer exported under the input cell order. The
entrypoint stops with `stage=drift_detected` before sidecars/complete when a runner row reports
`identity_verified=false`. Threshold failures with verified identity still retain scores and
unknown calls. Final byte-bound test evidence and decision are recorded in
`FINAL_REVIEW_2026-09-25.md`; no new production permit is created by this audit.
