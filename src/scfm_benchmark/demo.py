"""Small, explicitly synthetic demo. No patient records, FM weights or network."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from scfm_benchmark.benchmark import benchmark_from_csv


def make_demo_data(seed=17):
    rng = np.random.default_rng(seed)
    n_donors, cells_per_donor = 40, 8
    donors = np.repeat(np.arange(n_donors), cells_per_donor)
    labels = (donors % 2).astype(int)
    cell_ids = [f"synthetic-cell-{i:04d}" for i in range(len(donors))]
    metadata = pd.DataFrame({
        "cell_id": cell_ids, "donor_id": [f"donor-{d:02d}" for d in donors],
        # Two opposite-label donors share a dependency group in this example.
        "cv_group": [f"group-{d//2:02d}" for d in donors], "label": labels,
    })
    donor_effect = rng.normal(size=(n_donors, 12))
    features = donor_effect[donors] + rng.normal(scale=0.8, size=(len(donors), 12))
    features[:, :3] += labels[:, None] * 0.8
    noise = rng.normal(size=(len(donors), 12))
    expression = features @ rng.normal(size=(12, 30)) + rng.normal(size=(len(donors), 30))
    frames = {}
    for name, matrix in [("toy_signal", features), ("toy_noise", noise), ("toy_expression_pca", expression)]:
        frame = pd.DataFrame(matrix, columns=[f"feature_{i:03d}" for i in range(matrix.shape[1])])
        frame.insert(0, "cell_id", cell_ids)
        # Exercise alignment: the noise file deliberately has another row order.
        frames[name] = frame.iloc[::-1].reset_index(drop=True) if name == "toy_noise" else frame
    return metadata, frames


def run_demo(output_dir):
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=False)
    inputs = output / "inputs"
    inputs.mkdir()
    metadata, embeddings = make_demo_data()
    metadata_path = inputs / "metadata.csv"
    metadata.to_csv(metadata_path, index=False)
    paths = {}
    for name, frame in embeddings.items():
        paths[name] = inputs / f"{name}.csv"
        frame.to_csv(paths[name], index=False)
    return benchmark_from_csv(metadata_path, paths, output / "results", n_splits=5,
                              seeds=(42, 43, 44), pca={"toy_expression_pca": 5},
                              data_kind="synthetic-demo")
