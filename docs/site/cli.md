# Portable CLI

Install the current checkout with `python -m pip install .` (Python 3.10+).

## Synthetic demonstration

```bash
scfm-bench demo --out demo-output
```

Produces `inputs/*.csv` and `results/{report.html,metrics.csv,donor_predictions.csv,folds.csv,run.json}`.
No external data or model downloads. Each output directory must be new.

## Evaluate frozen representations

```bash
scfm-bench benchmark --metadata metadata.csv \
  --embedding model_a=model_a.csv --embedding model_b=model_b.csv \
  --folds 5 --seeds 42 43 44 --out comparison
```

Metadata requires `cell_id,donor_id,cv_group,label`. IDs must be non-empty strings;
cell IDs must be globally unique. Each donor has one binary label (0/1) and one
CV dependency group. Several donors may share a group. Enough groups must contain
each label for the requested fold count; one-class train/test folds are refused.

Each representation is a CSV with `cell_id` and numeric features, containing
exactly the metadata's cells. File order may differ. IDs are aligned explicitly;
missing/extra/duplicated cells and NaN/Inf features are errors.

Optional `--pca NAME=DIMENSIONS` fits PCA on the training cells of that named
representation. To compare against expression features, supply a separate
preprocessed expression-feature CSV; raw-count normalization and feature selection
are outside this package. Repeated `--embedding` and `--pca` names must be unique.

All representations share donor/group folds. Scaling and probes fit training
cells only. Cell predictions are averaged within donor before metrics are
computed. Training weights cells equally; metrics weight donors equally.
Across-repeat SD is not a confidence interval. Frozen embeddings can still carry
upstream leakage or confounding; the tool does not validate their provenance.

## Historical commands

`embed`, `evaluate`, `diagnose` and the old `report` remain for compatibility.
`embed`/`evaluate`/`diagnose` depend on research scripts and HPCC paths outside the
wheel; `embed` may submit SLURM jobs. They are not part of the portable demo.
The old `report` is a file-index renderer; use the report automatically written
by `demo` or `benchmark` for the new evaluation outputs.
