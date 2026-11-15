"""End-to-end test of prepare_provisional_diagnostic_input.py against a
small, real, on-disk MatrixMarket fixture -- exercising the actual streaming
reader (read_mtx_columns), h5ad write, and reload, exactly the path a
2026-10-09 review asked be verified before the real 72 GB scan is presented
for approval. Not a --dry-run check (which skips the streaming read
entirely) and not a mock of read_mtx_columns.
"""

import json
import subprocess
import sys
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.io import mmwrite

CODE_DIR = Path(__file__).parent.parent / "scripts/33_batch_effect_audit/code"
sys.path.insert(0, str(CODE_DIR))
from prepare_provisional_diagnostic_input import validate_event_column  # noqa: E402


def _build_fixture(root: Path, n_cells=40, n_genes=12, seed=0):
    """Genes x cells orientation (matches the real
    platelet_monocyte_benchmark_input_v2 files -- read_mtx_columns filters
    by COLUMN, i.e. cells)."""
    rng = np.random.RandomState(seed)
    (root / "raw").mkdir(parents=True)
    (root / "decontx").mkdir(parents=True)

    cohorts = ["CohortA", "CohortB"]
    lineages = ["Platelet", "Monocyte_Myeloid"]
    donors = [f"D{i:02d}" for i in range(8)]
    rows = []
    for i in range(n_cells):
        donor = donors[i % len(donors)]
        rows.append(
            dict(
                cell_id=f"cell_{i:04d}",
                cohort=cohorts[i % 2],
                donor_id=donor,
                harmonized_lineage=lineages[i % 2],
                event=(1.0 if i % 10 == 0 else (0.0 if i % 7 != 0 else np.nan)),
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
        density=0.3,
        random_state=seed,
        data_rvs=lambda n: rng.poisson(3, size=n).astype(np.float64) + 1,
    )
    decontx_counts = raw_counts.multiply(0.7)  # deterministic, distinguishable from raw
    mmwrite(str(root / "raw" / "counts.mtx"), raw_counts.tocoo())
    mmwrite(str(root / "decontx" / "counts.mtx"), decontx_counts.tocoo())
    return cells


def _run_cli(args_list):
    return subprocess.run(
        [sys.executable, str(CODE_DIR / "prepare_provisional_diagnostic_input.py")] + args_list,
        capture_output=True,
        text=True,
    )


def test_end_to_end_streaming_write_reload_on_real_small_mtx_fixture(tmp_path):
    input_root = tmp_path / "input"
    out_h5ad = tmp_path / "out" / "prepared.h5ad"
    _build_fixture(input_root)

    result = _run_cli(
        [
            "--input-root",
            str(input_root),
            "--out",
            str(out_h5ad),
            "--n-per-donor",
            "3",
            "--seed",
            "0",
        ]
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert out_h5ad.exists()

    manifest_path = out_h5ad.with_suffix(out_h5ad.suffix + ".selection_manifest.csv")
    assert manifest_path.exists()
    selection = pd.read_csv(manifest_path)
    assert set(selection.columns) >= {
        "cell_id",
        "cohort",
        "donor_id",
        "harmonized_lineage",
        "orig_row_index",
    }

    adata = anndata.read_h5ad(out_h5ad)
    assert adata.shape[0] == len(selection)
    assert set(adata.obs_names) == set(selection["cell_id"])
    assert adata.uns["population_stage"].startswith("pre-Gate-2 PROVISIONAL")
    assert adata.uns["expression_source"] == "decontx_counts"
    assert (
        "counts" in adata.layers
        and "counts_raw" in adata.layers
        and "counts_decontx" in adata.layers
    )

    # The decontx matrix was built as exactly 0.7x the raw matrix -- confirms
    # the streamed, column-filtered read preserved real VALUES, not just shape.
    raw_dense = np.asarray(adata.layers["counts_raw"].todense())
    decontx_dense = np.asarray(adata.layers["counts_decontx"].todense())
    np.testing.assert_allclose(decontx_dense, raw_dense * 0.7, rtol=1e-5)

    # Per-(cohort,donor,lineage) cap actually respected: no group exceeds 3.
    counts = selection.groupby(["cohort", "donor_id", "harmonized_lineage"]).size()
    assert counts.max() <= 3


def test_per_lineage_cap_does_not_starve_a_rarer_lineage(tmp_path):
    """Direct regression test for the review's finding: a donor with cells
    in both lineages must get up to n_per_donor of EACH lineage, not a
    combined cap split unpredictably between them."""
    input_root = tmp_path / "input"
    out_h5ad = tmp_path / "out" / "prepared.h5ad"
    # Build a fixture where one donor has lots of Monocyte_Myeloid cells and
    # very few Platelet cells.
    (input_root / "raw").mkdir(parents=True)
    (input_root / "decontx").mkdir(parents=True)
    rows = []
    for i in range(30):
        rows.append(
            dict(
                cell_id=f"cell_{i:04d}",
                cohort="CohortA",
                donor_id="D00",
                harmonized_lineage="Monocyte_Myeloid",
                event=0.0,
            )
        )
    for i in range(30, 35):
        rows.append(
            dict(
                cell_id=f"cell_{i:04d}",
                cohort="CohortA",
                donor_id="D00",
                harmonized_lineage="Platelet",
                event=0.0,
            )
        )
    cells = pd.DataFrame(rows)
    cells.to_csv(input_root / "raw" / "cells.csv", index=False)
    cells[["cell_id"]].to_csv(input_root / "decontx" / "cells.csv", index=False)
    pd.DataFrame({"gene": [f"GENE{i:03d}" for i in range(5)]}).to_csv(
        input_root / "genes.csv", index=False
    )
    rng = np.random.RandomState(0)
    raw_counts = sparse.random(
        5,
        35,
        density=0.5,
        random_state=0,
        data_rvs=lambda n: rng.poisson(3, size=n).astype(np.float64) + 1,
    )
    mmwrite(str(input_root / "raw" / "counts.mtx"), raw_counts.tocoo())
    mmwrite(str(input_root / "decontx" / "counts.mtx"), raw_counts.tocoo())

    result = _run_cli(
        [
            "--input-root",
            str(input_root),
            "--out",
            str(out_h5ad),
            "--n-per-donor",
            "20",
            "--seed",
            "0",
        ]
    )
    assert result.returncode == 0, result.stdout + result.stderr
    manifest_path = out_h5ad.with_suffix(out_h5ad.suffix + ".selection_manifest.csv")
    selection = pd.read_csv(manifest_path)
    n_platelet = (selection["harmonized_lineage"] == "Platelet").sum()
    n_monocyte = (selection["harmonized_lineage"] == "Monocyte_Myeloid").sum()
    assert (
        n_platelet == 5
    ), f"expected all 5 available Platelet cells kept (cap=20 > 5), got {n_platelet}"
    assert (
        n_monocyte == 20
    ), f"expected the cap (20) applied within Monocyte_Myeloid alone, got {n_monocyte}"


def test_mismatched_raw_decontx_cell_order_is_rejected(tmp_path):
    input_root = tmp_path / "input"
    out_h5ad = tmp_path / "out" / "prepared.h5ad"
    _build_fixture(input_root)
    # Corrupt decontx/cells.csv to have a different cell_id order than raw/cells.csv.
    decontx_cells = pd.read_csv(input_root / "decontx" / "cells.csv")
    decontx_cells = decontx_cells.iloc[::-1].reset_index(drop=True)
    decontx_cells.to_csv(input_root / "decontx" / "cells.csv", index=False)

    result = _run_cli(
        [
            "--input-root",
            str(input_root),
            "--out",
            str(out_h5ad),
            "--n-per-donor",
            "3",
            "--seed",
            "0",
        ]
    )
    assert result.returncode != 0
    assert "identical cell_id ordering" in result.stdout + result.stderr


def test_rejects_invalid_event_value():
    s = pd.Series([1.0, 0.0, np.nan, 2.0])
    with pytest.raises(ValueError, match="neither exactly 0"):
        validate_event_column(s)


def test_patient_id_is_cohort_qualified_not_bare_donor_id(tmp_path):
    """Direct regression test for the 2026-10-09 review's finding: the same
    donor_id string ("D00") reused across two DIFFERENT cohorts must not
    collide into one patient_id downstream."""
    input_root = tmp_path / "input"
    out_h5ad = tmp_path / "out" / "prepared.h5ad"
    (input_root / "raw").mkdir(parents=True)
    (input_root / "decontx").mkdir(parents=True)
    rows = []
    for i in range(10):
        rows.append(
            dict(
                cell_id=f"cellA_{i:03d}",
                cohort="CohortA",
                donor_id="D00",
                harmonized_lineage="Platelet",
                event=0.0,
            )
        )
    for i in range(10):
        rows.append(
            dict(
                cell_id=f"cellB_{i:03d}",
                cohort="CohortB",
                donor_id="D00",
                harmonized_lineage="Platelet",
                event=1.0,
            )
        )
    cells = pd.DataFrame(rows)
    cells.to_csv(input_root / "raw" / "cells.csv", index=False)
    cells[["cell_id"]].to_csv(input_root / "decontx" / "cells.csv", index=False)
    pd.DataFrame({"gene": [f"GENE{i:03d}" for i in range(5)]}).to_csv(
        input_root / "genes.csv", index=False
    )
    rng = np.random.RandomState(0)
    raw_counts = sparse.random(
        5,
        20,
        density=0.5,
        random_state=0,
        data_rvs=lambda n: rng.poisson(3, size=n).astype(np.float64) + 1,
    )
    mmwrite(str(input_root / "raw" / "counts.mtx"), raw_counts.tocoo())
    mmwrite(str(input_root / "decontx" / "counts.mtx"), raw_counts.tocoo())

    result = _run_cli(
        [
            "--input-root",
            str(input_root),
            "--out",
            str(out_h5ad),
            "--n-per-donor",
            "20",
            "--seed",
            "0",
        ]
    )
    assert result.returncode == 0, result.stdout + result.stderr
    adata = anndata.read_h5ad(out_h5ad)
    patient_ids = set(adata.obs["patient_id"])
    assert patient_ids == {"CohortA::D00", "CohortB::D00"}, (
        f"expected cohort-qualified patient_id values, got {patient_ids} -- a bare "
        f"'D00' would silently merge two different biological patients."
    )
    # The two "patients" must not share a single outcome -- if they had collided,
    # one of event=0.0/event=1.0 would have been lost/overwritten for 'D00'.
    assert set(adata.obs.loc[adata.obs["cohort"] == "CohortA", "outcome"]) == {"survived"}
    assert set(adata.obs.loc[adata.obs["cohort"] == "CohortB", "outcome"]) == {"fatal"}


def test_dry_run_keeps_selection_manifest_with_preflight_status(tmp_path):
    input_root = tmp_path / "input"
    out_h5ad = tmp_path / "out" / "prepared.h5ad"
    _build_fixture(input_root)

    result = _run_cli(
        [
            "--input-root",
            str(input_root),
            "--out",
            str(out_h5ad),
            "--n-per-donor",
            "3",
            "--seed",
            "0",
            "--dry-run",
        ]
    )
    assert result.returncode == 0, result.stdout + result.stderr
    manifest_path = out_h5ad.with_suffix(out_h5ad.suffix + ".selection_manifest.csv")
    status_path = out_h5ad.with_suffix(out_h5ad.suffix + ".selection_manifest.status.json")
    assert manifest_path.exists(), "dry-run must KEEP the selection manifest, not delete it"
    assert status_path.exists()
    status = json.loads(status_path.read_text())
    assert status["status"] == "preflight_only"
    assert not out_h5ad.exists(), "dry-run must not write the h5ad itself"


def test_selection_manifest_consume_matches_recompute_and_detects_drift(tmp_path):
    input_root = tmp_path / "input"
    dry_out = tmp_path / "dryrun" / "prepared.h5ad"
    real_out = tmp_path / "real" / "prepared.h5ad"
    _build_fixture(input_root)

    dry = _run_cli(
        [
            "--input-root",
            str(input_root),
            "--out",
            str(dry_out),
            "--n-per-donor",
            "3",
            "--seed",
            "0",
            "--dry-run",
        ]
    )
    assert dry.returncode == 0, dry.stdout + dry.stderr
    frozen_manifest = dry_out.with_suffix(dry_out.suffix + ".selection_manifest.csv")
    assert frozen_manifest.exists()
    frozen = pd.read_csv(frozen_manifest)

    consumed = _run_cli(
        [
            "--input-root",
            str(input_root),
            "--out",
            str(real_out),
            "--selection-manifest",
            str(frozen_manifest),
            "--seed",
            "0",
        ]
    )
    assert consumed.returncode == 0, consumed.stdout + consumed.stderr
    adata = anndata.read_h5ad(real_out)
    assert set(adata.obs["cell_id"]) == set(frozen["cell_id"])

    # Now corrupt the SOURCE cells.csv after the manifest was produced, and
    # confirm a consuming run detects the drift rather than extracting
    # against a stale selection.
    drifted = pd.read_csv(input_root / "raw" / "cells.csv")
    victim_row = frozen.iloc[0]["orig_row_index"]
    drifted.loc[victim_row, "donor_id"] = "SOMETHING_ELSE"
    drifted.to_csv(input_root / "raw" / "cells.csv", index=False)

    real_out2 = tmp_path / "real2" / "prepared.h5ad"
    drifted_run = _run_cli(
        [
            "--input-root",
            str(input_root),
            "--out",
            str(real_out2),
            "--selection-manifest",
            str(frozen_manifest),
            "--seed",
            "0",
        ]
    )
    assert drifted_run.returncode != 0
    assert "no longer match" in drifted_run.stdout + drifted_run.stderr


def test_extraction_failure_after_selection_does_not_report_extracted(tmp_path):
    """Direct regression test for the 2026-10-10 review's finding: the
    status sidecar must NOT say "extracted" if the process fails during the
    real streamed matrix read (after selection has already been published).
    Corrupts decontx/counts.mtx's MatrixMarket header for real -- not a
    mock of read_mtx_columns -- so the exact production code path that
    would run during a genuine failure (raw succeeds, decontx raises) is
    the one exercised here."""
    input_root = tmp_path / "input"
    out_h5ad = tmp_path / "out" / "prepared.h5ad"
    _build_fixture(input_root)

    decontx_mtx = input_root / "decontx" / "counts.mtx"
    original = decontx_mtx.read_text()
    assert original.startswith("%%MatrixMarket")
    decontx_mtx.write_text("NOT A MATRIXMARKET FILE\n" + original.split("\n", 1)[1])

    result = _run_cli(
        [
            "--input-root",
            str(input_root),
            "--out",
            str(out_h5ad),
            "--n-per-donor",
            "3",
            "--seed",
            "0",
        ]
    )
    assert result.returncode != 0, (
        "corrupting decontx/counts.mtx's header should make the real streaming "
        "read raise -- got a successful exit instead"
    )
    assert not out_h5ad.exists(), "a failed extraction must not leave a (partial) h5ad"

    status_path = out_h5ad.with_suffix(out_h5ad.suffix + ".selection_manifest.status.json")
    assert status_path.exists(), "the pre-extraction 'selected' status must have been published"
    status = json.loads(status_path.read_text())
    assert status["status"] == "failed", (
        f"expected status 'failed' after a real extraction failure, got {status['status']!r} "
        f"-- an interrupted/failed run must never look complete"
    )
    assert "error" in status and status["error"]
    manifest_path = out_h5ad.with_suffix(out_h5ad.suffix + ".selection_manifest.csv")
    assert manifest_path.exists(), "the selection manifest itself (the plan) is still durable"


def test_refuses_to_overwrite_existing_output(tmp_path):
    input_root = tmp_path / "input"
    out_h5ad = tmp_path / "out" / "prepared.h5ad"
    _build_fixture(input_root)
    out_h5ad.parent.mkdir(parents=True)
    out_h5ad.write_text("not a real h5ad, just occupying the path")

    result = _run_cli(
        [
            "--input-root",
            str(input_root),
            "--out",
            str(out_h5ad),
            "--n-per-donor",
            "3",
            "--seed",
            "0",
        ]
    )
    assert result.returncode != 0
    assert "already exists" in result.stdout + result.stderr
