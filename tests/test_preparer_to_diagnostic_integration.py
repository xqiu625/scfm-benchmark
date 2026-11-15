"""End-to-end integration test chaining prepare_provisional_diagnostic_input.py's
real output into run_batch_diagnostics.py, using the EXACT argument contract
PROVISIONAL_JOB_SPEC_v3.md documents for the real job -- the piece a 2026-10-09
review asked be verified before presenting that spec for approval ("Run the
prescribed HPCC tests and a small preparer-to-diagnostic integration fixture with
the exact argument contract intended for the real job"). Not a mock of either
CLI: both run as real subprocesses against a small, real, on-disk MatrixMarket
fixture, exactly as each is tested in isolation elsewhere
(test_prepare_provisional_diagnostic_input.py, test_execution_readiness.py).
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

CODE_DIR = Path(__file__).parent.parent / "scripts/33_batch_effect_audit/code"


def _build_fixture(root: Path, n_cells=80, n_genes=15, seed=0):
    rng = np.random.RandomState(seed)
    (root / "raw").mkdir(parents=True)
    (root / "decontx").mkdir(parents=True)
    cohorts = ["CohortA", "CohortB", "CohortC"]
    lineages = ["Platelet", "Monocyte_Myeloid"]
    donors = [f"D{i:02d}" for i in range(10)]
    rows = []
    for i in range(n_cells):
        rows.append(
            dict(
                cell_id=f"cell_{i:04d}",
                cohort=cohorts[i % 3],
                donor_id=donors[i % len(donors)],
                harmonized_lineage=lineages[i % 2],
                event=(1.0 if i % 9 == 0 else (0.0 if i % 11 != 0 else np.nan)),
            )
        )
    cells = pd.DataFrame(rows)
    cells.to_csv(root / "raw" / "cells.csv", index=False)
    cells[["cell_id"]].to_csv(root / "decontx" / "cells.csv", index=False)
    genes = pd.DataFrame({"gene": [f"GENE{i:03d}" for i in range(n_genes)]})
    genes.to_csv(root / "genes.csv", index=False)
    raw_counts = sparse.random(
        n_genes,
        n_cells,
        density=0.4,
        random_state=seed,
        data_rvs=lambda n: rng.poisson(3, size=n).astype(np.float64) + 1,
    )
    decontx_counts = raw_counts.multiply(0.7)
    mmwrite(str(root / "raw" / "counts.mtx"), raw_counts.tocoo())
    mmwrite(str(root / "decontx" / "counts.mtx"), decontx_counts.tocoo())


def test_preparer_output_feeds_diagnostic_cli_with_documented_job_spec_arguments(tmp_path):
    input_root = tmp_path / "input"
    prepared_h5ad = tmp_path / "prepared" / "input.h5ad"
    diag_out = tmp_path / "diagnostic_run"
    _build_fixture(input_root)

    # Step 1, exactly as PROVISIONAL_JOB_SPEC_v3.md's Step 1 command (minus
    # cluster-specific paths).
    prep = subprocess.run(
        [
            sys.executable,
            str(CODE_DIR / "prepare_provisional_diagnostic_input.py"),
            "--input-root",
            str(input_root),
            "--out",
            str(prepared_h5ad),
            "--n-per-donor",
            "20",
            "--seed",
            "0",
        ],
        capture_output=True,
        text=True,
    )
    assert prep.returncode == 0, prep.stdout + prep.stderr
    assert prepared_h5ad.exists()

    import anndata

    adata = anndata.read_h5ad(prepared_h5ad)
    assert adata.uns["population_stage_id"] == "pre_gate2_decontx_provisional"
    population_stage_text = str(adata.uns["population_stage"])

    # Step 2, exactly as PROVISIONAL_JOB_SPEC_v3.md's Step 2 command -- the
    # --population-stage value here is taken FROM the preparer's own output
    # (not independently retyped), which is the actual fix for the stage-text
    # mismatch finding: the two steps' free-text label is bound to one source
    # instead of two independently-authored strings that can drift.
    diag = subprocess.run(
        [
            sys.executable,
            str(CODE_DIR / "run_batch_diagnostics.py"),
            "--h5ad",
            str(prepared_h5ad),
            "--out-dir",
            str(diag_out),
            "--population-stage",
            population_stage_text,
            "--stage-id",
            "pre_gate2_decontx_provisional",
            "--code-revision",
            "integration-test-fixture",
            "--expression-source",
            "decontx_counts",
            "--batch-col",
            "cohort",
            "--patient-col",
            "patient_id",
            "--lineage-col",
            "lineage",
            "--outcome-col",
            "outcome",
            "--diagnostic-n-per-patient",
            "20",
            "--min-group-size",
            "2",
        ],
        capture_output=True,
        text=True,
    )
    assert diag.returncode == 0, diag.stdout + diag.stderr

    manifest = json.loads((diag_out / "manifest.json").read_text())
    assert manifest["status"] in ("complete", "complete_with_partial_results"), manifest
    assert manifest["input_sha256"] is not None
    assert manifest["stage_id"] == "pre_gate2_decontx_provisional"
    assert len(manifest["per_lineage"]) >= 1
    # No population_stage mismatch warning should have fired, since Step 2 was
    # given exactly the text Step 1 wrote (not an independently retyped copy).
    assert "does not match the input file's own" not in diag.stdout + diag.stderr
