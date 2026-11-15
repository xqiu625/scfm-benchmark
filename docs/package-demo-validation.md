# Portable package validation — 2026-11-13

This is software validation on synthetic data, not acceptance of a biological
result, a research-cohort definition, or the old manuscript's performance claims.

## Executed boundary

Code was copied to a new HPCC staging tree and executed there without a SLURM
submission. The command below uses `STAGING` as a placeholder for that directory;
the original absolute path is retained in the local validation logs.

```bash
bash source/scripts/00_setup/verify_portable_package.sh \
  "$STAGING/run1"
```

The base interpreter was the approved `302006_env_clean` Python 3.10.18.
A new venv inherited that numerical stack; build tools were installed only in
the new venv. This proves the wheel works independently of the checkout on that
stack, not that every supported dependency combination has been tested.

Observed command outputs in `run1.log`:

- Ruff: `All checks passed!`
- `python -m build --no-isolation`: built both the sdist and a wheel from the sdist.
- `pip install --no-deps --ignore-installed ...whl`: installed `1.2.0.dev0`.
- Import from outside the checkout resolved to `run1/venv/lib/python3.10/site-packages/scfm_benchmark/`.
- `python -m pytest .../tests/test_portable_benchmark.py`: **11 passed in 19.10s**.
- Installed `scfm-bench demo`: completed with 320 synthetic cells and 40 donors.
- Installed `scfm-bench benchmark` using the saved CSV inputs: completed;
  `cmp` confirmed byte-identical metrics and donor predictions to the demo.

The first script revision ran Ruff/build/tests/demo; the subsequent formatting-only
test edit was checked separately with Black, installed with `>=24,<25`
(`1 file would be left unchanged`).
The checked-in verification script now includes that Black check for future runs.
No production Python source changed after the accepted build.

## Fault injection

Executed scratch harness `portable_benchmark_mutations.py` on HPCC. Each case
copied a hash-verified baseline into its own directory, asserted the targeted
replacement occurred exactly once, verified the imported module path pointed to
the mutant, and ran its targeted test in a separate process.

**11/11 mutations were killed** (pytest exit 1 with one failed test): cell-ID
alignment removed; fitting changed to all cells; missing-ID guard removed;
duplicate-ID guard removed; finite-feature guard removed; donor-label guard
removed; donor-group guard removed; label-feature guard removed; input-drift
guard removed; report escaping removed; synthetic data-kind label corrupted.

The baseline source hash remained unchanged. This proves those targeted tests
detect the injected faults; it is not a universal leakage or correctness proof.
Case logs and `mutations/mutations.json` remain in the staging directory.

## Delivered identities

| Artifact | SHA-256 |
|---|---|
| `src/scfm_benchmark/benchmark.py` | `0f42141c9cc6b81488405a9201a9cbc85efb803e7e96a60fd46d01c6233f2c23` |
| `dist/scfm_benchmark-1.2.0.dev0-py3-none-any.whl` | `3bd2e4be7e0d82a532d442686b3eadf5d875cb58f95af96fafd94cfc270e24d3` |
| `dist/scfm_benchmark-1.2.0.dev0.tar.gz` | `658bc3bebd00b39e07c6682c4cf915ca87899c00d42a2eb25147f3a7f37e364c` |

The example's `run.json` records NumPy 2.2.6, pandas 2.3.1 and scikit-learn 1.7.1,
the actual seeds/PCA setting, source/input hashes and generated-file hashes.
The three requested seeds produced identical metrics in this balanced synthetic
example. These are not three independent validation cohorts; across-repeat SD
must not be interpreted as a confidence interval.

## GitHub validation

The `Portable package demo` workflow passed for the initial demo commit
[`dcf9a40`](https://github.com/xqiu625/scfm-benchmark/actions/runs/43562145111)
and the repository-link update
[`8d79b47`](https://github.com/xqiu625/scfm-benchmark/actions/runs/43562963399).
These clean runners built the distribution, installed the wheel and its
dependencies, and ran the portable tests and synthetic demo outside the checkout.
The hashes above describe the original HPCC build, not the later CI artifacts.

## Not checked or performed

- No real-data reanalysis, GPU inference, clinical claims or atlas-scale profiling.
- The full research test suite and legacy HPCC dispatchers were not rerun.
- Browser visual inspection could not run because the computer-use browser host
  was unavailable. The generated HTML/SVG, links and data tables are available
  for review; no pixel-level browser QA is claimed.
- No Pages deployment, release tag or PyPI upload was performed.
