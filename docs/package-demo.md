# Portable demo and package boundary

The deliverable is an installable evaluation utility with a runnable synthetic
example. It salvages useful engineering from the research project without
promoting old rankings, unvalidated QC decisions or manuscript status.

## Public interfaces

- `scfm-bench demo --out NEW_DIRECTORY`
- `scfm-bench benchmark --metadata CSV --embedding NAME=CSV --out NEW_DIRECTORY`
- `scfm_benchmark.benchmark.evaluate_embeddings(...)` for in-memory DataFrames.
- `scfm_benchmark.benchmark.benchmark_from_csv(...)` for persisted results.

The new commands import code exclusively from `src/scfm_benchmark`. They require
no repository-relative scripts/configs. The existing research CLI remains for
compatibility and is explicitly not advertised as portable. The same Hatchling
build configuration is retained; no second package or framework is introduced.

`1.2.0.dev0` is an unreleased development version. Package dependencies retain
the prior CPU research stack to avoid breaking existing users; future dependency
trimming should first inventory actual callers. GPU/model-specific integrations
remain optional and outside this demo.

## Build and verify

Project maintainers follow the repository's HPCC-only execution rule. The commands
below also describe how an external user can build the package on their machine.
Do not install build tools into a shared research environment; use an isolated venv.

```bash
python -m pip install '.[package-test]'
python -m build
python -m pytest tests/test_portable_benchmark.py
```

Then install the built wheel in a fresh environment and run `scfm-bench demo`
**from outside the checkout**. This last check is essential: editable installs
can conceal dependencies on `scripts/`, `configs/`, or local fixture files.
The dedicated `package-demo.yml` workflow builds a wheel and sdist, tests the
installed wheel, runs the example in a temporary directory, and retains artifacts.
It does not publish to PyPI, GitHub Pages or a release.

## Evaluation design

This is a portable probe, not an adoption of a new real-data research design.
The split is `StratifiedGroupKFold` over unique donor rows, stratified by binary
label and grouped by the explicit `cv_group`. All representations reuse those
folds. Every donor must have exactly one label/group. Each representation is
aligned to metadata's cell IDs before fitting.

A `StandardScaler` → optional PCA → `LogisticRegression` pipeline is fitted on
training cells and predicts only held-out cells. The model is fitted anew for
each representation/repeat/fold. A `DummyClassifier(strategy='prior')` gives the
training-cell prevalence comparator. Training weights each cell equally, so
donors with more cells affect fitting more. Predictions are then averaged within
donor and metrics weight donors equally. This limitation is visible in the report.

Data-derived upstream preprocessing, pretrained-model training overlap,
donor/capture identity adjudication, cohort confounding, longitudinal label design,
LOCO validation, confidence intervals and clinical calibration claims are outside
this MVP. The user supplies preprocessed feature tables and frozen embeddings;
this package does not turn arbitrary raw counts into a justified PCA baseline.
CSV inputs are loaded as dense matrices in memory. Full-atlas memory use and
large-scale runtime have not been benchmarked; this release targets a small,
inspectable demo and moderate-sized precomputed representations.

## Releasing it

1. Review the changed files and the generated demo; resolve any real-data scope
   separately from software acceptance.
2. Push a scoped branch only after the owner's push decision. Existing research
   results/checkpoints must remain untouched.
3. Let the package workflow test the proposed revision; review its wheel/demo
   artifacts. Add a release tag only for a version actually reviewed.
4. PyPI is a separate optional release. Check ownership/name availability and use
   a scoped publishing setup; do not imply `pip install scfm-benchmark` retrieves
   this development version until a release actually exists.

The project can be described in a portfolio as an **installable donor-aware
embedding evaluation tool with ID alignment, grouped validation, reproducible
artifacts and an offline demonstration**. Avoid presenting it as a published
benchmark, a validated clinical model or evidence that one FM is superior.

References: [Python Packaging User Guide](https://packaging.python.org/en/latest/tutorials/packaging-projects/)
and [scikit-learn: avoiding data leakage](https://scikit-learn.org/stable/common_pitfalls.html).
