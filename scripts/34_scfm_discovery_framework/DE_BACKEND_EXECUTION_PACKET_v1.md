# DE-backend execution packet v1 — consolidated for the next phase's review

Per `review:20261014_231259` C: consolidates the dependency-resolution evidence,
isolated-environment path, resource/execution commands, and synthetic-backend verification steps
already produced across this round and the prior one into a SINGLE document for reviewing
`ExpressionDERankingProvider`'s real-backend implementation phase. **This is a packaging step,
not new exploration** — no new dependency resolver run beyond what `DE_BACKEND_ENVIRONMENT_PLAN_v1.md`
already did; this document organizes that evidence for the next reviewer, and does not itself
authorize installing, training, fitting real data, or submitting SLURM.

**Revised per `review:20261015_104226` §四 item 5**: the prior version of this document
used elided (`.../`) paths and a version table that was not itself a locked, runnable command
sequence. This revision replaces every elided path with the actual absolute path, re-verifies the
three cited pip dry-run reports still exist with the stated hashes, adds an explicit install-failure
stop condition, a concrete new output directory for this phase, a draft SLURM resource proposal, and
restates environment-creation + native-3.11 verification + locked install + synthetic verification
as ONE scoped authorization request — while making explicit that writing this document is not that
authorization.

## 1. Target environment — absolute path, never the shared env

- **Proposed absolute path**: `/rhome/xinruqiu/bigdata/.conda/envs/302006_de_backend_v1`
- **Never** `302006_env_clean` — the shared environment this project's other stages depend on is
  not an install target for this backend. No command in this document touches it.
- **Primary Python/backend pin, fixed and unambiguous**: Python `3.11`, `pydeseq2==0.5.4`,
  `formulaic>=1.0.2`, `formulaic-contrasts>=0.2.0`, installed from the 40-package closure in §2
  below — not a looser "latest compatible" resolution at install time.
- **Fallback pin, same isolated env** (only if Python 3.11 provisioning is later found not to
  resolve cleanly — see §3's inconclusive conda check): Python `3.10`, `pydeseq2==0.5.2`, same
  `formulaic`/`formulaic-contrasts` floors.
- **New output directory for this phase's artifacts** (not reusing `claude_de_design_20261014_210331`,
  which belongs to the already-closed counts/decoder/roster round):
  `/bigdata/godziklab/shared/Xinru/302006/06_LOGS/gate2_capture_scoring/claude_de_backend_v1_<YYYYMMDD_HHMMSS>`,
  where `<YYYYMMDD_HHMMSS>` is the real wall-clock timestamp measured with `date +%Y%m%d_%H%M%S`
  immediately before the directory is created at actual execution time — never a value invented or
  estimated in this document, per `DEVELOPMENT.md` rule 4's timestamp discipline.

## 2. Dependency-resolution evidence — exact machine-readable reports, paths, and hashes

All three reports below are REAL, executed, read-only `pip install --dry-run` runs (never an
install) on `cluster.hpcc.ucr.edu`, from `/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean/bin/python
-m pip` (25.1.1). **Re-verified this round** (2026-10-15, via `sha256sum` over SSH against the live
files, not recalled from the prior write-up): all three files still exist at the stated paths and
all three hashes match exactly what the prior round recorded — no drift, no regeneration needed
this round.

| Report | HPCC path | sha256 (re-verified 2026-10-15) | What it proves |
|---|---|---|---|
| Full closure, `--ignore-installed` | `/bigdata/godziklab/shared/Xinru/302006/06_LOGS/gate2_capture_scoring/claude_de_design_20261014_210331/scratch_dir/pip_dry_run_report_py311_full.json` | `63206dc3a326b82739907da69890b6bea5e799fa9a960476df95a3cb417c086a` | REAL 40-package closure for `pydeseq2==0.5.4`/Python 3.11, not reusing `302006_env_clean`'s already-installed packages — the AUTHORITATIVE dependency-pin evidence |
| Incomplete, reuses installed packages (`py311`) | `/bigdata/godziklab/shared/Xinru/302006/06_LOGS/gate2_capture_scoring/claude_de_design_20261014_210331/scratch_dir/pip_dry_run_report_py311.json` | `163e99609e9154a6231fb55f6e79114f41d35351abe8fca5ef15e104257f8958` | Superseded by the full closure above; kept because it was the first pass that surfaced `narwhals`/`interface_meta`/`session_info`/`stdlib-list` |
| Incomplete, reuses installed packages (`py310`/`0.5.2` fallback) | `/bigdata/godziklab/shared/Xinru/302006/06_LOGS/gate2_capture_scoring/claude_de_design_20261014_210331/scratch_dir/pip_dry_run_report_py310_pydeseq052.json` | `d16c85e40e0852c5ac33b3bfb36852632865d08bb70baf68023fb73babdd2c77` | Same incompleteness caveat, for the fallback pin |

**Prior-round caveat, carried forward unchanged (per `review:20261015_104226` §四 item 6,
"protect existing evidence"):** these three files were previously accidentally deleted during an
over-broad cleanup command and regenerated from the same unmodified commands; the regeneration
produced an identical package set/size but this does NOT establish that the regenerated bytes are
byte-identical to whatever the very first run produced (pip's own JSON report formatting is not
guaranteed byte-stable across runs) — the hashes recorded above are of the CURRENTLY-LIVE,
regenerated files, re-verified again this round, not of some earlier, now-gone original.

**Residual, stated-not-resolved caveat** (per `DE_BACKEND_ENVIRONMENT_PLAN_v1.md`): pip's
cross-interpreter dry-run steers wheel-tag selection via `--python-version`/`--implementation`/
`--abi`/`--platform`, but the report's own `environment` block still reflects the RUNNING (3.10)
interpreter for PEP 508 marker evaluation, not the 3.11 target. The authoritative confirmation is
still a resolution run FROM an actual Python 3.11 interpreter once the isolated env exists (§3
below).

**Full 40-package closure** (from the `--ignore-installed` report's own `install` list):
`PyYAML-6.0.3 anndata-0.12.19 array-api-compat-1.15.0 cloudpickle-3.1.2 contourpy-1.3.2
cycler-0.12.1 donfig-0.8.1.post1 fonttools-4.67.0 formulaic-1.2.2 formulaic-contrasts-1.0.0
google-crc32c-1.9.0 h5py-3.14.0 interface_meta-2.0.1 joblib-1.6.0 kiwisolver-1.5.1
legacy-api-wrap-1.5 matplotlib-3.11.2 narwhals-2.27.1 natsort-8.4.0 numcodecs-0.16.5 numpy-2.2.6
packaging-26.3 pandas-2.3.2 pillow-12.2.0 pydeseq2-0.5.4 pyparsing-3.3.3
python-dateutil-2.9.0.post0 pytz-2026.5 scikit-learn-1.7.2 scipy-1.16.3 scverse-misc-0.0.3
session-info2-0.4.1 session_info-1.0.1 six-1.17.0 stdlib-list-0.12.0 threadpoolctl-3.7.0
typing_extensions-4.16.0 tzdata-2026.5 wrapt-2.5.1 zarr-3.1.6`.

**Install-failure stop condition (new this round):** if any package in this 40-package closure
fails to resolve or install into `302006_de_backend_v1` at the exact pinned version above, execution
stops immediately. The failing package name and the resolver's/pip's own error text are recorded
in the phase's output directory (§1) as `install_failure.log`; no unpinned or "nearest compatible"
substitute version is installed in its place, and no step downstream of the failure (native-3.11
dependency re-check, locked install into other packages, synthetic verification) proceeds until a
human has reviewed the failure and given a fresh, scoped decision — this is a genuine scope change
per §3's framing below, not a routine retry.

**Revised again per `review:20261015_142903` R4**: the 40-package closure above is now an
actual, committed, pip-constraints lock file — `results/de_backend_requirements.lock.txt` (40 lines,
one package per line, each pinned with `==`) — not only a prose version table. `pip install -c
results/de_backend_requirements.lock.txt pydeseq2==0.5.4 ...` fails loudly on any version this file
does not already pin, rather than silently resolving a "nearest compatible" substitute.

## 3. Environment-creation and dependency-check commands — one scoped authorization request

Per `review:20261015_104226` §四 item 5 and `review:20261015_142903` R4: environment
creation, native-Python-3.11 dependency verification, the locked install, and the synthetic
positive/null backend verification (§6) are bundled into ONE request a human can authorize in a
single decision, rather than four separate round-trips. **Writing this document is not that
authorization.** If anything about the scope changes during execution — a package outside the
40-package closure above, a different Python version, a different env path/prefix, a different
pinned `pydeseq2` version — execution stops and a fresh, scoped decision is requested; this
bundling does not pre-authorize scope changes.

**Rewritten per `review:20261015_150231` F1**, which found the prior version of
`results/de_backend_phase_wrapper_v1.sh` incomplete on 6 specific points: a lock-file path that
was never actually populated before use; an unbound deployment location ("or wherever this
phase's synced copy of the repo lives"); no plan for the test runner's (`pytest`'s) own
dependencies; treating `pip -c` as if it were a package whitelist (it only bounds the VERSION of
packages that get pulled in, never forbidding an unlisted extra, nor forcing every listed package
to be installed); an undetermined execution mode (an interactive `read -p` confirmation that a
non-interactive `sbatch` job has no way to answer); and treating a bare pytest exit code as proof
of a genuine pass when it cannot distinguish PASSED from SKIPPED. All six are now fixed in the
rewritten `results/de_backend_phase_wrapper_v1.sh` (`bash -n`-checked on HPCC; its own preflight
guards were run — safely, failing before `conda create` — against a wrong commit and a missing
SLURM allocation, both confirmed to stop the script at the right point for the right reason; NOT
otherwise run, NOT submitted to SLURM):

1. **Preflight, before anything is created or touched**: `REPO_CHECKOUT_DIR`/`EXPECTED_REPO_COMMIT`
   are REQUIRED inputs (no default) — the script runs `git -C "$REPO_CHECKOUT_DIR" rev-parse HEAD`
   and refuses outright on any mismatch, never "or wherever." `SLURM_JOB_ID`/`SLURM_CPUS_PER_TASK`
   are also REQUIRED — this decides the execution-mode question explicitly: the wrapper only runs
   inside a real SLURM allocation (`sbatch`/`salloc`), never interactively on the login node, which
   is also how it avoids retrying the 2026-10-14 login-node OOM without a separate manual
   reminder. The target `ENV_PREFIX` must not already exist. Only after all of these pass does the
   script create its output directory and copy both lock files plus the two new checker scripts
   (below) into it.
2. `conda create --prefix .../302006_de_backend_v1 -c conda-forge python=3.11 pip
   --repodata-fn current_repodata.json` — reduces the conda-forge metadata pulled into memory, per
   the already-identified mitigation; running inside an allocation (Step 1 above) is what actually
   keeps this off the login node.
3. A native dependency check from the NEW env's own `pip install --dry-run --report ... -c
   de_backend_requirements.lock.txt pydeseq2==0.5.4 ...` — a real resolution, not a cross-interpreter
   dry run, closing the residual caveat in §2.
4. **An AUTOMATED gate, replacing the prior `read -p` prompt**: `results/
   de_backend_verify_native_resolution.py` (new this round, functionally tested against the real
   dry-run reports already on record — both the matching-closure and a deliberately-mismatched
   case) parses the native report and the lock file and compares their package sets EXACTLY, in
   both directions (nothing extra, nothing missing, no version drift) — this is what actually
   enforces the closure, since `-c` alone does not. Any mismatch exits immediately with the exact
   differences; this is the genuine stop condition for a non-interactive job, not a human keypress.
5. TWO separate locked installs, from TWO separate lock files: `de_backend_requirements.lock.txt`
   (pydeseq2/formulaic — unchanged from §2) and the NEW `results/de_backend_test_tooling.lock.txt`
   (`pytest==9.1.1` and its own real transitive closure — `pluggy`, `iniconfig`, `exceptiongroup`,
   `Pygments`, `tomli`, `typing_extensions`, `packaging` — resolved via a REAL, executed
   `pip install --dry-run --report` run against `de_backend_requirements.lock.txt` as a constraint,
   2026-10-15, exit 0, no conflict). Keeping this separate makes clear pytest was never part of the
   original pydeseq2 dependency-pin evidence. On ANY failure in either install: `set -euo pipefail`
   stops immediately; stderr is captured to `install_stderr.log`; no unpinned substitute, no Step 6.
6. **The synthetic positive/null backend verification (§6)** — reuses the REAL, already-written test
   `tests/test_scfm_gene_ranking_providers.py::test_run_pydeseq2_from_prepared_design_known_effect_and_null_ranks_correctly`
   (currently `pytest.importorskip`-gated; this step is what makes it actually run instead of skip)
   — never a second, parallel synthetic-fixture script duplicating that test's logic. Output is
   REQUIRED to be structured (`--junitxml`), and `results/de_backend_verify_synthetic_test_passed.py`
   (new this round, functionally tested against both a real SKIPPED run and a real PASSED run of an
   unrelated test in this environment) parses that XML and refuses unless the named test genuinely
   PASSED — no `<skipped>`/`<failure>`/`<error>` child, and exactly one matching testcase. A bare
   pytest exit code or a log's last line is never treated as sufficient on its own.

A real, synced git checkout of this repository already exists on HPCC at
`/bigdata/godziklab/shared/Xinru/302006/302006_updategit` (confirmed via the preflight-guard test
above — it is currently at an OLDER commit than this round's work, so it would need re-syncing to
whatever commit is actually authorized before `REPO_CHECKOUT_DIR` points at it for a real run).

A `conda create --dry-run -c conda-forge python=3.11 pip` check (also read-only) was attempted on
the login node 2026-10-14 and **killed by OOM** before completing conda-forge metadata
resolution — reported as inconclusive, not negative, per `DEVELOPMENT.md`'s existing caution against
memory-heavy operations on a login node. **Neither this round nor Step 1 above retries that exact
solve on the login node.**

## 4. Resource/execution notes

- Interpreter once created: `/rhome/xinruqiu/bigdata/.conda/envs/302006_de_backend_v1/bin/python`.
- Standard project environment variables apply: `PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1
  OPENBLAS_NUM_THREADS=1`.
- All real fitting/testing happens on HPCC, never the login node for anything beyond trivial
  metadata reads (per this round's own OOM lesson, both for the conda solver above and for the
  Stephenson raw-metadata extraction in `results/stephenson_raw_metadata_extraction.py`).
- No GPU requirement identified for PyDESeq2 itself (a CPU-bound negative-binomial GLM fit); not
  independently re-verified against the package's own resource documentation this round.

**Draft SLURM resource proposal (new this round, NOT SUBMITTED, not self-authorizing):**

| Resource | Proposed value | Justification |
|---|---|---|
| Partition | `batch` (this project's existing default non-interactive partition; not re-verified against current cluster partition names this round — confirm before submitting) | Matches every other non-interactive job this project has run; no GPU/high-memory partition is implicated below |
| Nodes | 1 | PyDESeq2's GLM fit is single-process; this framework's design matrices are sample-count-scale (tens to low hundreds of donors), not a distributed workload |
| CPUs | 4 | Headroom for pandas/numpy's own internal BLAS threading on a sample x gene matrix of this project's scale (thousands of genes x tens of samples), while `OMP_NUM_THREADS=1`/`OPENBLAS_NUM_THREADS=1` are still set per this project's standing convention to keep the fit deterministic and single-threaded in practice |
| Memory | 8 GB | A dense counts/design DataFrame at this project's gene x sample scale is well under 1 GB; 8 GB leaves headroom for pandas/anndata overhead and the synthetic positive/null fixtures run in the same job, without approaching this project's documented OOM threshold on much larger inputs |
| Wall time | 30 minutes | A single negative-binomial GLM fit plus `lfc_shrink` over this project's scale is expected to run in minutes, not hours; 30 minutes leaves margin without parking a long reservation |

These are draft planning numbers only, carried into any future SLURM submission request, not
submitted by this document and not a substitute for the per-job authorization `DEVELOPMENT.md` already
requires. **Concrete `sbatch`/`salloc` directives matching the table above, and matching
`gene_ranking_providers._bounded_n_cpus()`'s own reading of `SLURM_CPUS_PER_TASK`** (per
`review:20261015_142903` R4, confirmed implemented): `--nodes=1 --ntasks=1
--cpus-per-task=4 --mem=8G --time=00:30:00 --partition=batch` (partition name not re-verified
against current cluster partition names this round — confirm before submitting).

## 5. Backend-call draft — verified against the ACTUAL recommended version's source

Per `DE_BACKEND_ENVIRONMENT_PLAN_v1.md`'s "Backend call draft" section (re-verified directly
against `pydeseq2` tag `v0.5.4`'s own source, confirmed identical to `v0.5.2`'s):

- `DeseqDataSet(counts=<pandas.DataFrame, samples x genes>, metadata=<pandas.DataFrame>,
  design=<pandas.DataFrame, samples x design columns>)` — accepts the design as a DataFrame
  directly, so `de_design_providers.PreparedDEDesign.design_matrix`/`.column_names` feed it
  without any re-derivation through a formula engine.
- `DeseqStats(dds, contrast=<numpy.ndarray>)` — a numeric contrast vector aligned to the design's
  own columns; `PreparedDEDesign.contrast_vector` is exactly this shape. **Per
  `review:20261015_142903` R2, now implemented**: the real code constructs this as
  `np.asarray(..., dtype=float)`, never a plain Python `list` — pydeseq2 v0.5.4's own constructor
  branches on this TYPE, and a list would be misread as a 3-item categorical contrast spec,
  raising `IndexError` for this project's 2-column design. Mutation-killed: reverting to a list is
  caught by `tests/test_scfm_gene_ranking_providers.py::
  test_run_pydeseq2_from_prepared_design_passes_an_ndarray_contrast_not_a_list`.
- `lfc_shrink(coeff=<verified column name>, adapt=True)` — `coeff` has NO usable default despite
  its own docstring text; the target coefficient column name must be read from
  `PreparedDEDesign.contrast_name` (already hash/recompute-bound), never left as `None`.
- **Per `review:20261015_142903` R4, now implemented**: both `DeseqDataSet(...)` and
  `DeseqStats(...)` are called with an explicit `n_cpus=` (`gene_ranking_providers._bounded_n_cpus()`),
  read from the job's own `SLURM_CPUS_PER_TASK`/`SLURM_CPUS_ON_NODE`, falling back to `1` (never
  `os.cpu_count()`) if neither is set. Verified directly against pydeseq2 v0.5.4's own source
  (`utils.get_num_processes`): `n_cpus=None` defaults to "the maximum number of available cpus" via
  `joblib.Parallel(n_jobs=...)` — a SEPARATE knob from this project's `OMP_NUM_THREADS`/
  `OPENBLAS_NUM_THREADS` convention, which bounds BLAS threads only, never joblib's own process pool.

## 6. Synthetic positive/null verification plan — now built, not yet run against real pydeseq2

Per `review:20261015_142903` R3: the fixture is now a 7-gene design (one large known
up-regulated effect, one large known down-regulated effect, one genuinely-null gene, and FOUR
stable background/normalization-anchor genes whose values vary per sample but whose group means
match within ~2%) — replacing a prior 3-gene fixture whose "null" gene was a hard-coded CONSTANT
across every sample, which a reviewer's own independent median-ratio diagnostic showed distorted
normalization (that null gene's own normalized log2FC came out ≈-2.68, not ≈0). This fixture's
own normalization was independently sanity-checked the same way (median-of-ratios against each
gene's geometric mean, by hand, not a pydeseq2 run): resulting size factors 0.98–1.02 across all 6
samples, normalized log2FC ≈+3.29/-3.28 for the two known effects, and |log2FC| < 0.02 for the
null and all four background genes. See `tests/test_scfm_gene_ranking_providers.py`'s
`_RB_X`/`_RB_GENE_NAMES` and the `test_run_pydeseq2_from_prepared_design_known_effect_and_null_ranks_correctly`
test's own docstring/assertions (including a raw-Wald-stat-sign-vs-shrunken-LFC-sign consistency
check, and a `contrast_name` consistency check). **Wording correction per
`review:20261015_150231` item 4**: the prior fixture's defect was not "a constant gene
can never be a valid normalization anchor" in general — with only 3 genes total and 2 of them
carrying a large injected shift, the single constant gene ended up as essentially the ONLY stable
reference the geometric-mean/median-of-ratios calculation had to work with, which skewed the
result; the 7-gene fixture's fix is having enough stable, varying background genes relative to
the genes that actually shift. **This 7-gene panel is a numeric SMOKE TEST only** — it verifies
the interface, direction, ranking, and write/read-back round-trip once pydeseq2 exists; it does
NOT establish genome-wide FDR control, statistical power, or dispersion-estimation calibration,
and must never be cited as if it does. Still gated on `pytest.importorskip("pydeseq2")` — not
executed against a real `pydeseq2` install this round; this is Step 4 in §3's single scoped
request, gated on Steps 1–3 succeeding first.

## 7. What this packet is NOT

Not an installation, not an environment creation, not a new dependency-resolver run beyond what
is cited above, not a decision about the real pilot's cohort/subtype/positive-control/checkpoint
(those remain `pending_owner`, tracked in `MONOCYTE_PILOT_PREREG_v1.md`/`MODEL_DECISION_v1.md`).
This packet exists so the NEXT phase (real `ExpressionDERankingProvider` backend implementation)
has one place to review the environment plan from, rather than re-deriving it. **Writing or revising
this document is not itself the "one scoped authorization" §3 describes** — that authorization is a
separate, explicit decision this document is organized to make easy to grant or withhold in one
pass, not a decision this document makes on its own.
