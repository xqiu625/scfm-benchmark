"""Portable binary probes of frozen embeddings; no project paths or model downloads."""

from __future__ import annotations

import hashlib
import html
import json
import platform
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.dummy import DummyClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from scfm_benchmark import __version__

REQUIRED_COLUMNS = ["cell_id", "donor_id", "cv_group", "label"]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _validate_inputs(metadata, embeddings):
    if not isinstance(metadata, pd.DataFrame) or metadata.empty:
        raise ValueError("metadata must be a non-empty DataFrame")
    if not metadata.columns.is_unique or not set(REQUIRED_COLUMNS) <= set(metadata):
        raise ValueError(f"metadata needs unique columns including {REQUIRED_COLUMNS}")
    meta = metadata[REQUIRED_COLUMNS].copy()
    for column in REQUIRED_COLUMNS[:3]:
        if meta[column].isna().any() or not meta[column].map(
            lambda x: isinstance(x, str) and bool(x.strip())
        ).all():
            raise ValueError(f"{column} must contain non-empty string identifiers")
    if meta.cell_id.duplicated().any():
        raise ValueError("cell_id must be unique; include sample identity in barcodes")
    if meta.label.isna().any() or set(meta.label.unique()) != {0, 1}:
        raise ValueError("label must contain both 0 and 1, without missing values")
    meta["label"] = meta.label.astype(int)
    grouped = meta.groupby("donor_id", sort=True)
    if (grouped.label.nunique() != 1).any():
        raise ValueError("each donor must have exactly one binary label")
    if (grouped.cv_group.nunique() != 1).any():
        raise ValueError("a donor cannot belong to multiple cv_group values")
    if not embeddings or "prior" in embeddings:
        raise ValueError("provide at least one embedding; 'prior' is a reserved model name")
    aligned = {}
    for name, frame in embeddings.items():
        if not isinstance(name, str) or not name.strip():
            raise ValueError("model names must be non-empty strings")
        if not isinstance(frame, pd.DataFrame) or not frame.columns.is_unique:
            raise ValueError(f"{name}: expected a DataFrame with unique columns")
        if "cell_id" not in frame or frame.shape[1] < 2:
            raise ValueError(f"{name}: needs cell_id and at least one numeric feature")
        if frame.cell_id.isna().any() or frame.cell_id.duplicated().any():
            raise ValueError(f"{name}: cell_id must be non-missing and unique")
        if set(frame.cell_id) != set(meta.cell_id):
            raise ValueError(f"{name}: cell IDs must exactly match metadata; no silent intersection")
        features = frame.drop(columns="cell_id")
        if set(features) & set(REQUIRED_COLUMNS[1:]):
            raise ValueError(f"{name}: metadata/label columns cannot be features")
        if not all(pd.api.types.is_numeric_dtype(dtype) for dtype in features.dtypes):
            raise ValueError(f"{name}: all feature columns must be numeric")
        matrix = frame.set_index("cell_id").loc[meta.cell_id].to_numpy(dtype=float)
        if not np.isfinite(matrix).all():
            raise ValueError(f"{name}: features must be finite; imputation is not automatic")
        aligned[name] = matrix
    return meta.reset_index(drop=True), aligned


def evaluate_embeddings(metadata, embeddings, *, n_splits=5, seeds=(42, 43, 44), pca=None):
    """Return donor OOF predictions, per-repeat metrics and fold membership.

    Metadata requires cell_id, donor_id, cv_group, label (0/1). Each embedding
    DataFrame has cell_id plus numeric features. Rows are aligned by cell ID.
    ``pca`` optionally maps a representation name to a fold-fitted dimension.
    The split uses one row per donor; all donors sharing cv_group stay together.
    LR is fitted on cells, then held-out probabilities are averaged per donor.
    This lightweight probe is not the historical manuscript evaluation protocol.
    """
    meta, matrices = _validate_inputs(metadata, embeddings)
    seeds = list(seeds)
    if not seeds or any(type(s) is not int or not 0 <= s < 2**32 for s in seeds):
        raise ValueError("seeds must be non-empty integers in [0, 2**32)")
    if len(set(seeds)) != len(seeds):
        raise ValueError("seeds must be unique")
    if type(n_splits) is not int or n_splits < 2:
        raise ValueError("n_splits must be an integer >= 2")
    pca = dict(pca or {})
    if set(pca) - set(matrices) or any(type(n) is not int or n < 1 for n in pca.values()):
        raise ValueError("pca must map existing model names to positive integer dimensions")
    donors = meta.groupby("donor_id", sort=True).agg(
        label=("label", "first"), cv_group=("cv_group", "first"), n_cells=("cell_id", "size")
    )
    if donors.cv_group.nunique() < n_splits:
        raise ValueError("not enough independent cv_group values for n_splits")
    if donors.groupby("label").cv_group.nunique().min() < n_splits:
        raise ValueError("each label needs at least n_splits distinct cv_group values")

    plans = []
    fold_rows = []
    for seed in seeds:
        splitter = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        seen = np.zeros(len(meta), dtype=int)
        for fold, (train_d, test_d) in enumerate(
            splitter.split(donors, donors.label, donors.cv_group)
        ):
            train = np.flatnonzero(meta.donor_id.isin(donors.index[train_d]))
            test = np.flatnonzero(meta.donor_id.isin(donors.index[test_d]))
            if set(meta.cv_group.iloc[train]) & set(meta.cv_group.iloc[test]):
                raise ValueError("cv_group leakage across train and test")
            if any(meta.label.iloc[idx].nunique() != 2 for idx in (train, test)):
                raise ValueError(f"seed {seed}, fold {fold}: one-class fold; use fewer folds")
            for name, dimension in pca.items():
                if dimension > min(len(train), matrices[name].shape[1]):
                    raise ValueError(f"{name}: PCA dimension exceeds a training fold's size")
            seen[test] += 1
            plans.append((seed, fold, train, test))
            for donor in donors.index[test_d]:
                fold_rows.append({"seed": seed, "fold": fold, "donor_id": donor,
                                  "cv_group": donors.loc[donor, "cv_group"]})
        if not np.all(seen == 1):
            raise ValueError("each cell must be held out exactly once per repeat")

    prediction_parts = []
    for seed, fold, train, test in plans:
        y_train = meta.label.to_numpy()[train]
        for name in [*matrices, "prior"]:
            if name == "prior":
                x = np.zeros((len(meta), 1))
                model = DummyClassifier(strategy="prior")
            else:
                x = matrices[name]
                steps = [StandardScaler()]
                if name in pca:
                    steps.append(PCA(n_components=pca[name], svd_solver="full"))
                steps.append(LogisticRegression(max_iter=1000, random_state=seed))
                model = make_pipeline(*steps)
            model.fit(x[train], y_train)
            scores = model.predict_proba(x[test])[:, list(model.classes_).index(1)]
            if not np.isfinite(scores).all() or np.any((scores < 0) | (scores > 1)):
                raise ValueError("classifier produced invalid probabilities")
            held_out = meta.iloc[test].copy()
            held_out["probability"] = scores
            donor_predictions = held_out.groupby("donor_id", sort=True).agg(
                label=("label", "first"), probability=("probability", "mean"),
                n_cells=("cell_id", "size"), cv_group=("cv_group", "first")
            ).reset_index()
            donor_predictions["model"] = name
            donor_predictions["seed"] = seed
            donor_predictions["fold"] = fold
            prediction_parts.append(donor_predictions)
    predictions = pd.concat(prediction_parts, ignore_index=True)
    rows = []
    for (name, seed), group in predictions.groupby(["model", "seed"], sort=True):
        if len(group) != len(donors) or group.donor_id.duplicated().any():
            raise ValueError("expected one held-out prediction per donor, model and repeat")
        rows.append({"model": name, "seed": int(seed), "n_donors": len(group),
                     "roc_auc": roc_auc_score(group.label, group.probability),
                     "average_precision": average_precision_score(group.label, group.probability),
                     "brier_score": brier_score_loss(group.label, group.probability)})
    return {"predictions": predictions, "metrics": pd.DataFrame(rows),
            "folds": pd.DataFrame(fold_rows)}


def _write_report(metrics, out, data_kind):
    summary = metrics.groupby("model", sort=True).agg(
        auc_mean=("roc_auc", "mean"), auc_sd=("roc_auc", lambda s: s.std(ddof=1) if len(s) > 1 else 0.0),
        ap_mean=("average_precision", "mean"), brier_mean=("brier_score", "mean")
    )
    bars = []
    for i, (name, row) in enumerate(summary.iterrows()):
        y = 70 + i * 40
        bars.append(f'<text x="8" y="{y+16}"><title>{html.escape(name)}</title>{html.escape(name[:26])}</text>'
                    f'<rect x="240" y="{y}" width="{row.auc_mean*400:.1f}" height="22" fill="#296b87"/>'
                    f'<text x="655" y="{y+16}">{row.auc_mean:.3f}</text>')
    banner = ("SYNTHETIC DEMO — generated data, not biological results or model rankings."
              if data_kind == "synthetic-demo" else
              "USER-SUPPLIED DATA — exploratory probe; QC and upstream training are not audited.")
    chart_title = "Synthetic demo only — not foundation-model results" if data_kind == "synthetic-demo" else "Exploratory donor-level evaluation"
    chart = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 740 {len(summary)*40+100}" '
             'role="img" aria-label="Mean donor AUROC on a zero to one scale">'
             '<style>text{font:14px system-ui,sans-serif;fill:#18313b}</style>'
             f'<rect width="740" height="{len(summary)*40+100}" fill="white"/>'
             f'<text x="8" y="24">{chart_title}</text>'
             f'<line x1="440" x2="440" y1="55" y2="{len(summary)*40+75}" stroke="#7d8a90" stroke-dasharray="4"/>'
             '<text x="235" y="54">0</text><text x="430" y="54">0.5</text><text x="635" y="54">1</text>'
             + ''.join(bars) + '</svg>')
    out.with_name("preview.svg").write_text(chart, encoding="utf-8")
    page = f'''<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>scfm-benchmark | donor-level evaluation</title>
<style>body{{font:16px/1.6 system-ui,sans-serif;max-width:960px;margin:40px auto;padding:0 24px;color:#18313b}}
h1{{line-height:1.2}}.notice{{padding:16px;background:#fff0cb;border-left:4px solid #9b6400}}
table{{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}}th,td{{padding:10px;border-bottom:1px solid #ccd9de;text-align:right}}
th:first-child,td:first-child{{text-align:left}}svg{{width:100%;height:auto}}a{{color:#176480}}code{{background:#edf2f4;padding:2px 5px}}</style>
<h1>Frozen embeddings, auditable donor-level evaluation</h1><p class="notice">{banner}</p>
<p>{int(metrics.n_donors.iloc[0])} donors · {metrics.seed.nunique()} repeats · binary labels 0/1</p>
<h2>Mean pooled out-of-fold donor AUROC</h2>
{chart}
{summary.to_html(float_format=lambda x: f"{x:.4f}", border=0)}
<p>SD describes variation across the specified repeats; it is not a confidence interval.
AUROC/AP: higher is better. Brier score: lower is better. The prior baseline uses training-cell class prevalence.</p>
<h2>Evaluation contract</h2><ul>
<li>Cells align by cell ID. Donors have one label and one dependency group.</li>
<li>Splits use donor rows and keep entire dependency groups together; every representation uses identical folds.</li>
<li>Scaling, optional PCA and logistic regression fit training cells only. Training weights each cell equally.</li>
<li>Held-out cell probabilities are averaged within donor; metrics weight each donor equally.</li>
<li>This probes fixed input representations. It cannot rule out leakage during upstream pretraining, fine-tuning or feature selection.</li></ul>
<h2>Inspect the evidence</h2><p><a href="donor_predictions.csv">Donor predictions</a> ·
<a href="metrics.csv">Per-repeat metrics</a> · <a href="folds.csv">Donor folds</a> ·
<a href="run.json">Parameters, versions and hashes</a></p></html>'''
    out.write_text(page, encoding="utf-8")


def benchmark_from_csv(metadata_path, embedding_paths, output_dir, *, n_splits=5,
                       seeds=(42, 43, 44), pca=None, data_kind="user-supplied"):
    """Run the portable file interface into a new directory; write run.json last."""
    output = Path(output_dir)
    seeds = list(seeds)
    if output.exists():
        raise FileExistsError(f"output already exists; choose a new directory: {output}")
    metadata_path = Path(metadata_path)
    embedding_paths = {name: Path(path) for name, path in embedding_paths.items()}
    sources = {"metadata": metadata_path, **{f"embedding:{n}": p for n, p in embedding_paths.items()}}
    hashes = {name: _sha256(path) for name, path in sources.items()}
    metadata = pd.read_csv(metadata_path, dtype={c: str for c in REQUIRED_COLUMNS[:3]})
    embeddings = {name: pd.read_csv(path, dtype={"cell_id": str})
                  for name, path in embedding_paths.items()}
    result = evaluate_embeddings(metadata, embeddings, n_splits=n_splits, seeds=seeds, pca=pca)
    if hashes != {name: _sha256(path) for name, path in sources.items()}:
        raise ValueError("input changed during evaluation; refusing to publish results")
    output.mkdir(parents=True, exist_ok=False)
    for key, filename in [("predictions", "donor_predictions.csv"), ("metrics", "metrics.csv"),
                          ("folds", "folds.csv")]:
        result[key].to_csv(output / filename, index=False)
    _write_report(result["metrics"], output / "report.html", data_kind)
    manifest = {
        "status": "complete", "schema_version": 1, "data_kind": data_kind,
        "package_version": __version__, "python": platform.python_version(),
        "dependencies": {name: version(name) for name in ("numpy", "pandas", "scikit-learn")},
        "source_sha256": _sha256(Path(__file__)), "input_sha256": hashes,
        "parameters": {"n_splits": n_splits, "seeds": list(seeds), "pca": dict(pca or {}),
                       "probe": "StandardScaler + optional PCA + LogisticRegression(max_iter=1000)",
                       "training_unit": "cell", "metric_unit": "donor",
                       "donor_aggregation": "mean held-out probability", "split_group": "cv_group"},
        "n_cells": len(metadata), "n_donors": metadata.donor_id.nunique(),
        "artifacts": {path.name: _sha256(path) for path in sorted(output.iterdir())},
    }
    (output / "run.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return manifest
