"""Mutation-verified tests for `validate_job1_output.py` against a REAL
production h5ad (built by actually running
`prepare_provisional_diagnostic_input.py`'s CLI). v2 adds coverage for the
F4 defects `review:20261010_202340` found in v1: unenforced nnz
budgets, `--cells-per-group 0` producing a vacuous pass, tolerance-based
layer-copy checks, missing stage_id/expression_source/manifest-field
comparisons, and no CLI-entry-point test (only internal functions).
"""

import os
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

CODE_DIR = Path(__file__).parent.parent / "scripts/33_batch_effect_audit/code"
TESTS_DIR = Path(__file__).parent
sys.path.insert(0, str(CODE_DIR))
sys.path.insert(0, str(TESTS_DIR))

import validate_job1_output as vmod  # noqa: E402
from test_prepare_provisional_diagnostic_input import _build_fixture, _run_cli  # noqa: E402
from validate_job1_output import (  # noqa: E402
    ValidationError,
    _select_sample_rows,
    metadata_checks,
    sampled_numeric_checks,
)


def _produce_real_output(tmp_path, n_per_donor=3):
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
            str(n_per_donor),
            "--seed",
            "0",
        ]
    )
    assert result.returncode == 0, result.stdout + result.stderr
    manifest = out_h5ad.with_suffix(out_h5ad.suffix + ".selection_manifest.csv")
    genes_csv = input_root / "genes.csv"
    assert out_h5ad.exists() and manifest.exists() and genes_csv.exists()
    return out_h5ad, manifest, genes_csv


def _run_full_validation(h5ad_path, manifest_path, genes_csv, **kwargs):
    meta = metadata_checks(h5ad_path, manifest_path, genes_csv, **kwargs)
    obs_df = meta.pop("_obs_for_sampling")
    numeric = sampled_numeric_checks(h5ad_path, obs_df, cells_per_group=2, seed=0)
    return meta, numeric


def test_passes_on_real_unmutated_production_output(tmp_path):
    h5ad_path, manifest_path, genes_csv = _produce_real_output(tmp_path)
    meta, numeric = _run_full_validation(h5ad_path, manifest_path, genes_csv)

    assert meta["cell_id_unique"] is True
    assert meta["manifest_rows_match_n_obs"] is True
    assert meta["cell_id_order_matches_selection_manifest"] is True
    assert meta["manifest_cohort_donor_lineage_match_obs_independent"] is True
    assert meta["patient_id_cohort_qualified_independent"] is True
    assert meta["per_lineage_counts_reconciled_independent"] is True
    assert meta["var_order_matches_source_genes_csv"] is True
    assert meta["stage_id_matches_expected"] is True
    assert meta["expression_source_matches_expected"] is True

    assert numeric["sampled_indices_within_gene_bounds"] is True
    assert numeric["sampled_values_finite_and_nonnegative"] is True
    assert numeric["sampled_X_equals_counts_decontx_exact"] is True
    assert numeric["sampled_counts_layer_equals_counts_decontx_exact"] is True
    assert "NOT_VERIFIED" in numeric["source_mtx_full_numeric_equality"]


def test_detects_swapped_cell_order_against_manifest(tmp_path):
    h5ad_path, manifest_path, genes_csv = _produce_real_output(tmp_path)
    with h5py.File(h5ad_path, "r+") as f:
        ds = f["obs"]["cell_id"]
        before = ds[:2].copy()
        ds[0], ds[1] = ds[1], before[0]
        after = ds[:2].copy()
    assert not np.array_equal(before, after)

    meta, _numeric = _run_full_validation(h5ad_path, manifest_path, genes_csv)
    assert meta["cell_id_order_matches_selection_manifest"] is False


def test_detects_corrupted_sampled_layer_value(tmp_path):
    h5ad_path, manifest_path, genes_csv = _produce_real_output(tmp_path)
    meta = metadata_checks(h5ad_path, manifest_path, genes_csv)
    obs_df = meta.pop("_obs_for_sampling")
    sample = _select_sample_rows(obs_df, cells_per_group=2, seed=0)
    target_pos = next(iter(sample.values()))[0]

    with h5py.File(h5ad_path, "r+") as f:
        indptr = f["layers"]["counts_decontx"]["indptr"][:]
        lo, hi = int(indptr[target_pos]), int(indptr[target_pos + 1])
        assert hi > lo
        data_ds = f["layers"]["counts_decontx"]["data"]
        before = data_ds[lo]
        data_ds[lo] = before + 10_000.0
        after = data_ds[lo]
    assert after != before

    _meta, numeric = _run_full_validation(h5ad_path, manifest_path, genes_csv)
    assert (
        numeric["sampled_X_equals_counts_decontx_exact"] is False
        or numeric["sampled_counts_layer_equals_counts_decontx_exact"] is False
    )


def test_exact_equality_catches_a_tiny_difference_allclose_would_miss(tmp_path):
    """F4: v1 used np.allclose for the layer-copy checks -- a change too
    small for the default allclose tolerance would have passed. The fixed
    exact np.array_equal check must still catch it."""
    h5ad_path, manifest_path, genes_csv = _produce_real_output(tmp_path)
    meta = metadata_checks(h5ad_path, manifest_path, genes_csv)
    obs_df = meta.pop("_obs_for_sampling")
    sample = _select_sample_rows(obs_df, cells_per_group=2, seed=0)
    target_pos = next(iter(sample.values()))[0]

    with h5py.File(h5ad_path, "r+") as f:
        indptr = f["layers"]["counts_decontx"]["indptr"][:]
        lo, hi = int(indptr[target_pos]), int(indptr[target_pos + 1])
        assert hi > lo
        data_ds = f["layers"]["counts_decontx"]["data"]
        before = float(data_ds[lo])
        # Must be a handful of float32 ULPs (~1e-7 relative) so it actually
        # rounds to a DIFFERENT stored value, yet stay inside np.allclose's
        # default tolerance (atol=1e-8, rtol=1e-5) so allclose would miss it.
        tiny_delta = 5e-6
        data_ds[lo] = np.float32(before + tiny_delta)
        after = float(data_ds[lo])
    assert after != before
    assert np.allclose(before, after)  # confirms allclose WOULD have missed this

    _meta, numeric = _run_full_validation(h5ad_path, manifest_path, genes_csv)
    assert (
        numeric["sampled_X_equals_counts_decontx_exact"] is False
        or numeric["sampled_counts_layer_equals_counts_decontx_exact"] is False
    )


def test_detects_corrupted_categorical_decoding(tmp_path):
    h5ad_path, manifest_path, genes_csv = _produce_real_output(tmp_path)
    with h5py.File(h5ad_path, "r+") as f:
        cats_ds = f["obs"]["cohort"]["categories"]
        before = cats_ds[:].copy()
        assert len(before) >= 2
        swapped = before.copy()
        swapped[0], swapped[1] = before[1], before[0]
        cats_ds[:] = swapped
        after = cats_ds[:].copy()
    assert not np.array_equal(before, after)

    meta, _numeric = _run_full_validation(h5ad_path, manifest_path, genes_csv)
    assert meta["patient_id_cohort_qualified_independent"] is False


def test_refuses_above_max_obs_rows_ceiling(tmp_path, monkeypatch):
    h5ad_path, manifest_path, genes_csv = _produce_real_output(tmp_path)
    monkeypatch.setattr(vmod, "MAX_OBS_ROWS_FOR_METADATA_READ", 1)
    with pytest.raises(ValidationError):
        vmod.metadata_checks(h5ad_path, manifest_path, genes_csv)


def test_rejects_zero_cells_per_group_instead_of_vacuous_pass(tmp_path):
    """F4: v1 silently accepted --cells-per-group 0, producing an empty
    sample whose numeric checks were vacuously True."""
    h5ad_path, manifest_path, genes_csv = _produce_real_output(tmp_path)
    meta = metadata_checks(h5ad_path, manifest_path, genes_csv)
    obs_df = meta.pop("_obs_for_sampling")
    with pytest.raises(ValidationError):
        vmod.sampled_numeric_checks(h5ad_path, obs_df, cells_per_group=0, seed=0)


def test_enforces_max_nnz_per_row_budget_before_reading(tmp_path, monkeypatch):
    """F4: a per-row/cumulative nnz budget must be enforced BEFORE any
    data/indices slice, not just reported after the fact."""
    h5ad_path, manifest_path, genes_csv = _produce_real_output(tmp_path)
    meta = metadata_checks(h5ad_path, manifest_path, genes_csv)
    obs_df = meta.pop("_obs_for_sampling")
    monkeypatch.setattr(vmod, "MAX_NNZ_PER_SAMPLED_ROW", 0)
    with pytest.raises(ValidationError):
        vmod.sampled_numeric_checks(h5ad_path, obs_df, cells_per_group=2, seed=0)


def test_detects_stage_id_mismatch_against_expected_contract(tmp_path):
    h5ad_path, manifest_path, genes_csv = _produce_real_output(tmp_path)
    meta = metadata_checks(
        h5ad_path, manifest_path, genes_csv, expected_stage_id="some_other_stage_id_entirely"
    )
    assert meta["stage_id_matches_expected"] is False


def test_detects_manifest_field_mismatch_against_obs(tmp_path):
    """F4: v1 compared only cell_id order, never cohort/donor_id/
    harmonized_lineage row-by-row against obs."""
    h5ad_path, manifest_path, genes_csv = _produce_real_output(tmp_path)
    import pandas as pd

    manifest = pd.read_csv(manifest_path)
    before = manifest.loc[0, "cohort"]
    manifest.loc[0, "cohort"] = f"NOT_{before}"
    manifest.to_csv(manifest_path, index=False)
    after = pd.read_csv(manifest_path).loc[0, "cohort"]
    assert after != before

    meta = metadata_checks(h5ad_path, manifest_path, genes_csv)
    assert meta["manifest_cohort_donor_lineage_match_obs_independent"] is False


def test_cli_entry_point_exits_nonzero_on_a_real_failure(tmp_path):
    """F4: v1's tests covered only internal functions, never the CLI
    itself (sys.exit path, stderr reporting)."""
    h5ad_path, manifest_path, genes_csv = _produce_real_output(tmp_path)
    with h5py.File(h5ad_path, "r+") as f:
        ds = f["obs"]["cell_id"]
        ds[0], ds[1] = ds[1], ds[0]

    result = subprocess.run(
        [
            sys.executable,
            str(CODE_DIR / "validate_job1_output.py"),
            "--h5ad",
            str(h5ad_path),
            "--manifest",
            str(manifest_path),
            "--genes-csv",
            str(genes_csv),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "FAILED checks" in result.stderr


def test_cli_refuses_to_overwrite_existing_out_report(tmp_path):
    h5ad_path, manifest_path, genes_csv = _produce_real_output(tmp_path)
    out_path = tmp_path / "report.json"
    out_path.write_text("already here")
    result = subprocess.run(
        [
            sys.executable,
            str(CODE_DIR / "validate_job1_output.py"),
            "--h5ad",
            str(h5ad_path),
            "--manifest",
            str(manifest_path),
            "--genes-csv",
            str(genes_csv),
            "--out",
            str(out_path),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "already exists" in (result.stdout + result.stderr)


# ---- v3 (20261011_113512.md F4): deeper structural guarantees --------------


def test_detects_x_shape_attr_disagreeing_with_obs_index_length(tmp_path):
    """v2's `_n_obs_from_metadata_only` trusted `X.attrs['shape']` alone;
    a corrupted value that still looked plausible would pass silently."""
    h5ad_path, manifest_path, genes_csv = _produce_real_output(tmp_path)
    with h5py.File(h5ad_path, "r+") as f:
        shape = f["X"].attrs["shape"]
        f["X"].attrs["shape"] = (int(shape[0]) + 1, int(shape[1]))

    with pytest.raises(vmod.ValidationError):
        vmod.metadata_checks(h5ad_path, manifest_path, genes_csv)


def test_detects_corrupted_layer_shape_attribute(tmp_path):
    """A layer's own declared `shape` attribute (distinct from X's) must
    be checked against (n_obs, n_vars); v2 never read it at all."""
    h5ad_path, manifest_path, genes_csv = _produce_real_output(tmp_path)
    with h5py.File(h5ad_path, "r+") as f:
        grp = f["layers"]["counts_decontx"]
        shape = grp.attrs["shape"]
        grp.attrs["shape"] = (int(shape[0]), int(shape[1]) + 1)

    with pytest.raises(vmod.ValidationError):
        vmod.metadata_checks(h5ad_path, manifest_path, genes_csv)


def test_detects_non_integer_indptr_dtype(tmp_path):
    """v2 never checked indptr's dtype -- a float-typed indptr with the
    right VALUES would have passed every v2 structural check."""
    h5ad_path, manifest_path, genes_csv = _produce_real_output(tmp_path)
    with h5py.File(h5ad_path, "r+") as f:
        grp = f["layers"]["counts_decontx"]
        indptr = grp["indptr"][:]
        del grp["indptr"]
        grp.create_dataset("indptr", data=indptr.astype(np.float64))

    with pytest.raises(vmod.ValidationError):
        vmod.metadata_checks(h5ad_path, manifest_path, genes_csv)


def test_detects_wrong_rank_indices_array(tmp_path):
    """v2 never checked indices' rank. Reshaping to (-1, 1) deliberately
    keeps `.shape[0]` (and therefore the pre-existing indptr[-1]==
    len(indices) length check) UNCHANGED, so this isolates the new rank
    check specifically -- a mutation that only removed the rank check
    (confirmed: see review:20261011_113512 F4 mutation log)
    would make this test, and only this test among the three new
    structural tests, fail to raise."""
    h5ad_path, manifest_path, genes_csv = _produce_real_output(tmp_path)
    with h5py.File(h5ad_path, "r+") as f:
        grp = f["layers"]["counts_decontx"]
        indices = grp["indices"][:]
        del grp["indices"]
        grp.create_dataset("indices", data=indices.reshape(-1, 1))

    with pytest.raises(vmod.ValidationError):
        vmod.metadata_checks(h5ad_path, manifest_path, genes_csv)


def test_cli_out_creation_is_exclusive_not_check_then_write(tmp_path):
    """The fix replaces Path.exists()-then-write with os.O_CREAT|O_EXCL.
    This test exercises the real CLI path (not just the internal check)
    and confirms a pre-existing file is still refused."""
    h5ad_path, manifest_path, genes_csv = _produce_real_output(tmp_path)
    out_path = tmp_path / "exclusive_report.json"
    fd = os.open(str(out_path), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    os.close(fd)

    result = subprocess.run(
        [
            sys.executable,
            str(CODE_DIR / "validate_job1_output.py"),
            "--h5ad",
            str(h5ad_path),
            "--manifest",
            str(manifest_path),
            "--genes-csv",
            str(genes_csv),
            "--out",
            str(out_path),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "already exists" in (result.stdout + result.stderr)
