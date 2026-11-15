"""Portable package boundaries, using synthetic cells and the real CSV/CLI path."""

import hashlib
import json
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
from sklearn.pipeline import Pipeline

import scfm_benchmark.benchmark as benchmark
from scfm_benchmark.demo import make_demo_data


def test_csv_cli_roundtrip_and_no_overwrite(tmp_path):
    out = tmp_path / "demo"
    command = [sys.executable, "-m", "scfm_benchmark.cli", "demo", "--out", str(out)]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    output = out / "results"
    manifest = json.loads((output / "run.json").read_text())
    assert manifest["status"] == "complete"
    assert manifest["data_kind"] == "synthetic-demo"
    assert manifest["n_cells"] == 320 and manifest["n_donors"] == 40
    assert manifest["parameters"]["pca"] == {"toy_expression_pca": 5}
    for name, digest in manifest["artifacts"].items():
        assert hashlib.sha256((output / name).read_bytes()).hexdigest() == digest
    predictions = pd.read_csv(output / "donor_predictions.csv")
    metrics = pd.read_csv(output / "metrics.csv")
    folds = pd.read_csv(output / "folds.csv")
    assert len(predictions) == 40 * 4 * 3
    assert len(metrics) == 4 * 3
    assert len(folds) == 40 * 3
    assert (folds.groupby(["seed", "cv_group"]).fold.nunique() == 1).all()
    assert (predictions.groupby(["seed", "donor_id"]).fold.nunique() == 1).all()
    assert not predictions.duplicated(["seed", "model", "donor_id"]).any()
    assert "SYNTHETIC DEMO" in (output / "report.html").read_text()
    before = (output / "run.json").read_bytes()
    repeated = subprocess.run(command, capture_output=True, text=True)
    assert repeated.returncode != 0
    assert (output / "run.json").read_bytes() == before


def test_alignment_uses_cell_identity_not_file_order():
    metadata, frames = make_demo_data()
    frames = {"toy_signal": frames["toy_signal"]}
    first = benchmark.evaluate_embeddings(metadata, frames, n_splits=3, seeds=[42])
    shuffled = {name: frame.sample(frac=1, random_state=4) for name, frame in frames.items()}
    second = benchmark.evaluate_embeddings(metadata, shuffled, n_splits=3, seeds=[42])
    pd.testing.assert_frame_equal(first["predictions"], second["predictions"])


def test_actual_training_state_excludes_held_out_cells(monkeypatch):
    metadata, frames = make_demo_data()
    frames = {"toy_signal": frames["toy_signal"]}
    states = []
    original_fit = Pipeline.fit

    def observed_fit(self, x, y, **kwargs):
        result = original_fit(self, x, y, **kwargs)
        states.append(
            (
                self.named_steps["standardscaler"].mean_.copy(),
                self.named_steps["pca"].components_.copy(),
                len(x),
            )
        )
        return result

    monkeypatch.setattr(Pipeline, "fit", observed_fit)
    first = benchmark.evaluate_embeddings(
        metadata, frames, n_splits=3, seeds=[42], pca={"toy_signal": 3}
    )
    original_state = states[0]
    test_donors = first["folds"].query("fold == 0").donor_id
    test_ids = set(metadata.loc[metadata.donor_id.isin(test_donors), "cell_id"])
    changed = frames["toy_signal"].copy()
    mask = changed.cell_id.isin(test_ids)
    changed.loc[mask, changed.columns != "cell_id"] += 1000
    assert original_state[2] == len(metadata) - len(test_ids)
    states.clear()
    second = benchmark.evaluate_embeddings(
        metadata, {"toy_signal": changed}, n_splits=3, seeds=[42], pca={"toy_signal": 3}
    )
    pd.testing.assert_frame_equal(first["folds"], second["folds"])
    np.testing.assert_array_equal(original_state[0], states[0][0])
    np.testing.assert_array_equal(original_state[1], states[0][1])


@pytest.mark.parametrize(
    "fault, message",
    [
        ("missing_id", "exactly match"),
        ("duplicate_id", "unique"),
        ("nonfinite", "finite"),
        ("donor_label", "one binary label"),
        ("donor_group", "multiple cv_group"),
        ("label_feature", "cannot be features"),
    ],
)
def test_invalid_inputs_fail_before_fitting(monkeypatch, fault, message):
    metadata, frames = make_demo_data()
    frames = {"toy_signal": frames["toy_signal"]}
    frame = frames["toy_signal"]
    if fault == "missing_id":
        frames["toy_signal"] = frame.iloc[:-1]
    elif fault == "duplicate_id":
        frame.loc[0, "cell_id"] = frame.loc[1, "cell_id"]
    elif fault == "nonfinite":
        frame.loc[0, "feature_000"] = np.inf
    elif fault == "donor_label":
        metadata.loc[0, "label"] = 1
    elif fault == "donor_group":
        metadata.loc[0, "cv_group"] = "different-group"
    else:
        frame["label"] = metadata.label

    def forbid_fit(*args, **kwargs):
        pytest.fail("invalid data reached fitting")

    monkeypatch.setattr(Pipeline, "fit", forbid_fit)
    with pytest.raises(ValueError, match=message):
        benchmark.evaluate_embeddings(metadata, frames)


def test_input_change_refuses_publication(tmp_path, monkeypatch):
    metadata, frames = make_demo_data()
    metadata_path = tmp_path / "metadata.csv"
    frame_path = tmp_path / "embedding.csv"
    metadata.to_csv(metadata_path, index=False)
    frames["toy_signal"].to_csv(frame_path, index=False)
    original = benchmark.evaluate_embeddings

    def changed_input(*args, **kwargs):
        result = original(*args, **kwargs)
        metadata_path.write_text(metadata_path.read_text() + "\n")
        return result

    monkeypatch.setattr(benchmark, "evaluate_embeddings", changed_input)
    with pytest.raises(ValueError, match="input changed"):
        benchmark.benchmark_from_csv(metadata_path, {"signal": frame_path}, tmp_path / "result")
    assert not (tmp_path / "result").exists()


def test_report_escapes_model_names(tmp_path):
    metrics = pd.DataFrame(
        [
            {
                "model": '<img src=x onerror="bad()">',
                "seed": 42,
                "n_donors": 20,
                "roc_auc": 0.5,
                "average_precision": 0.5,
                "brier_score": 0.25,
            }
        ]
    )
    benchmark._write_report(metrics, tmp_path / "report.html", "user-supplied")
    text = (tmp_path / "report.html").read_text()
    assert '<img src=x onerror="bad()">' not in text
    assert "&lt;img" in text
    assert "USER-SUPPLIED DATA" in text
