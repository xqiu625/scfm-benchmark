"""Tests for `pseudobulk_providers.PseudobulkAggregationProvider` -- the
first real increment toward `ExpressionDERankingProvider`'s intended
backend, per `review:20261013_222450` B: raw-count-sum
aggregation ONLY, no DE statistic, no p-value, no `GeneRanking`. Every
fixture here is a tiny, synthetic, local h5ad with HAND-COMPUTABLE sums --
no network, no real data.
"""

import os
import sys

import numpy as np
import pytest

CODE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "scripts",
    "34_scfm_discovery_framework",
    "code",
)
sys.path.insert(0, CODE)

import pseudobulk_providers as PB  # noqa: E402
from dataset_schema import NO_BATCH_SENTINEL  # noqa: E402


def _base_config(matrix_path, **overrides):
    config = {
        "dataset_id": "test_pseudobulk_pilot",
        "species": "human",
        "tissue_or_target_cell_type": "monocyte",
        "subject_id_column": "donor_id",
        "sample_id_column": "sample_id",
        "biological_replicate_unit": "donor_id",
        "technical_batch_column": NO_BATCH_SENTINEL,
        "expression_matrix_path": matrix_path,
        "counts_semantics": "raw_counts",
        "counts_slot": "X",
        "gene_id_namespace": "gene_symbol",
        "cell_id_column": "cell_barcode",
        "qc_status": "provisional",
        "qc_version_ref": "none (provisional)",
    }
    config.update(overrides)
    return config


def _fake_config_path(tmp_path):
    return str(tmp_path / "dataset_config.json")


#: Hand-computable PAIRED fixture: 2 donors, each contributing TWO samples
#: (one healthy/v1, one disease/v2) -- D0: S0+S1, D1: S2+S3.
_PAIRED_X = np.array(
    [
        [1, 2, 3, 4],  # CELL0 -> S0
        [5, 6, 7, 8],  # CELL1 -> S0      S0 sum = [6, 8, 10, 12], lib=36
        [1, 1, 1, 1],  # CELL2 -> S1
        [2, 2, 2, 2],  # CELL3 -> S1      S1 sum = [3, 3, 3, 3],   lib=12
        [0, 0, 1, 1],  # CELL4 -> S2
        [0, 0, 1, 1],  # CELL5 -> S2      S2 sum = [0, 0, 2, 2],   lib=4
        [10, 0, 0, 0],  # CELL6 -> S3
        [0, 10, 0, 0],  # CELL7 -> S3      S3 sum = [10, 10, 0, 0], lib=20
    ],
    dtype="float32",
)
_PAIRED_SAMPLE = ["S0", "S0", "S1", "S1", "S2", "S2", "S3", "S3"]
_PAIRED_DONOR = ["D0", "D0", "D0", "D0", "D1", "D1", "D1", "D1"]
_PAIRED_CONDITION = [
    "healthy",
    "healthy",
    "disease",
    "disease",
    "healthy",
    "healthy",
    "disease",
    "disease",
]
_PAIRED_VISIT = ["v1", "v1", "v2", "v2", "v1", "v1", "v2", "v2"]
_PAIRED_CELL_IDS = [f"CELL{i}" for i in range(8)]


def _make_h5ad(
    tmp_path,
    *,
    X,
    cell_ids,
    sample_ids,
    donor_ids,
    conditions=None,
    visits=None,
    name="pseudobulk.h5ad",
):
    import anndata as ad
    import pandas as pd

    obs = {"cell_barcode": cell_ids, "donor_id": donor_ids, "sample_id": sample_ids}
    if conditions is not None:
        obs["condition"] = conditions
    if visits is not None:
        obs["visit"] = visits
    obs_df = pd.DataFrame(obs, index=[f"idx{i}" for i in range(len(cell_ids))])
    var = pd.DataFrame(index=[f"GENE{i}" for i in range(X.shape[1])])
    a = ad.AnnData(X=np.asarray(X, dtype="float32"), obs=obs_df, var=var)
    path = tmp_path / name
    a.write_h5ad(path)
    return str(path)


def _paired_fixture(tmp_path, **config_overrides):
    h5ad = _make_h5ad(
        tmp_path,
        X=_PAIRED_X,
        cell_ids=_PAIRED_CELL_IDS,
        sample_ids=_PAIRED_SAMPLE,
        donor_ids=_PAIRED_DONOR,
        conditions=_PAIRED_CONDITION,
        visits=_PAIRED_VISIT,
    )
    config = _base_config(h5ad, **config_overrides)
    return h5ad, config


# ===========================================================================
# Core aggregation correctness -- hand-computable sums, conservation,
# sample-vs-donor count distinction, no fabricated zero rows.
# ===========================================================================


def test_pseudobulk_aggregates_paired_structure_with_hand_computable_sums(tmp_path):
    h5ad, config = _paired_fixture(tmp_path)
    provider = PB.PseudobulkAggregationProvider()
    result = provider.aggregate(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pb_out"),
        cell_ids_subset=tuple(_PAIRED_CELL_IDS),
        cell_type_or_state_label="monocyte",
        condition_column="condition",
        visit_or_state_column="visit",
    )

    assert result.status == "complete", result.reason
    assert result.n_pseudobulk_samples == 4
    assert result.n_distinct_donors == 2  # D0, D1 -- each contributes TWO samples

    rows = {r["sample_id"]: r for r in result.sample_rows}
    assert rows["S0"]["donor_or_animal_id"] == "D0"
    assert rows["S0"]["condition_value"] == "healthy"
    assert rows["S0"]["visit_or_state_value"] == "v1"
    assert rows["S0"]["n_cells"] == 2
    assert rows["S0"]["library_size"] == 36.0
    assert rows["S1"]["library_size"] == 12.0
    assert rows["S2"]["library_size"] == 4.0
    assert rows["S3"]["library_size"] == 20.0
    assert sorted(rows["S0"]["cell_ids"]) == ["CELL0", "CELL1"]

    assert result.total_counts_selected_cells == 72.0
    assert result.total_counts_aggregated == 72.0

    manifest = PB.read_pseudobulk_result(str(tmp_path / "pb_out" / "PSEUDOBULK_RESULT.json"))
    counts = PB.load_pseudobulk_counts(manifest)
    assert counts.shape == (4, 4)
    row_by_sid = {r["sample_id"]: i for i, r in enumerate(manifest["sample_rows"])}
    np.testing.assert_array_equal(counts[row_by_sid["S0"]], [6, 8, 10, 12])
    np.testing.assert_array_equal(counts[row_by_sid["S1"]], [3, 3, 3, 3])
    np.testing.assert_array_equal(counts[row_by_sid["S2"]], [0, 0, 2, 2])
    np.testing.assert_array_equal(counts[row_by_sid["S3"]], [10, 10, 0, 0])


def test_pseudobulk_unpaired_structure_reports_matching_sample_and_donor_counts(tmp_path):
    """Per B.4: the SAME interface handles an unpaired structure (one
    sample per donor) -- n_pseudobulk_samples == n_distinct_donors here,
    unlike the paired fixture above where they differ.
    """
    cell_ids = ["CELLA", "CELLB", "CELLC", "CELLD"]
    X = np.array([[1, 1], [1, 1], [2, 2], [2, 2]], dtype="float32")
    h5ad = _make_h5ad(
        tmp_path,
        X=X,
        cell_ids=cell_ids,
        sample_ids=["S4", "S4", "S5", "S5"],
        donor_ids=["D2", "D2", "D3", "D3"],
        conditions=["healthy", "healthy", "disease", "disease"],
    )
    config = _base_config(h5ad)
    provider = PB.PseudobulkAggregationProvider()
    result = provider.aggregate(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pb_out"),
        cell_ids_subset=tuple(cell_ids),
        cell_type_or_state_label="monocyte",
        condition_column="condition",
    )
    assert result.status == "complete", result.reason
    assert result.n_pseudobulk_samples == result.n_distinct_donors == 2


def test_pseudobulk_same_interface_supports_mouse_species(tmp_path):
    h5ad, config = _paired_fixture(tmp_path, species="mouse", dataset_id="test_mouse_pilot")
    provider = PB.PseudobulkAggregationProvider()
    result = provider.aggregate(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pb_out"),
        cell_ids_subset=tuple(_PAIRED_CELL_IDS),
        cell_type_or_state_label="monocyte",
        condition_column="condition",
    )
    assert result.status == "complete", result.reason
    assert result.dataset_identity["species"] == "mouse"


def test_pseudobulk_never_fabricates_a_zero_cell_sample_row(tmp_path):
    """Selecting cells for only S0/S1/S2 (excluding S3 entirely) must
    produce exactly 3 rows -- NOT a fabricated all-zero S3 row.
    """
    h5ad, config = _paired_fixture(tmp_path)
    provider = PB.PseudobulkAggregationProvider()
    subset = tuple(c for c, s in zip(_PAIRED_CELL_IDS, _PAIRED_SAMPLE) if s != "S3")
    result = provider.aggregate(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pb_out"),
        cell_ids_subset=subset,
        cell_type_or_state_label="monocyte",
        condition_column="condition",
    )
    assert result.status == "complete", result.reason
    sample_ids = {r["sample_id"] for r in result.sample_rows}
    assert sample_ids == {"S0", "S1", "S2"}


def test_pseudobulk_row_order_independent_of_cell_ids_subset_ordering(tmp_path):
    h5ad, config = _paired_fixture(tmp_path)
    provider = PB.PseudobulkAggregationProvider()
    shuffled = tuple(reversed(_PAIRED_CELL_IDS))
    r1 = provider.aggregate(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "run1"),
        cell_ids_subset=tuple(_PAIRED_CELL_IDS),
        cell_type_or_state_label="monocyte",
        condition_column="condition",
    )
    r2 = provider.aggregate(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "run2"),
        cell_ids_subset=shuffled,
        cell_type_or_state_label="monocyte",
        condition_column="condition",
    )
    assert r1.status == r2.status == "complete"
    c1 = PB.load_pseudobulk_counts(
        PB.read_pseudobulk_result(str(tmp_path / "run1" / "PSEUDOBULK_RESULT.json"))
    )
    c2 = PB.load_pseudobulk_counts(
        PB.read_pseudobulk_result(str(tmp_path / "run2" / "PSEUDOBULK_RESULT.json"))
    )
    np.testing.assert_array_equal(c1, c2)
    assert [r["sample_id"] for r in r1.sample_rows] == [r["sample_id"] for r in r2.sample_rows]


# ===========================================================================
# Structural/input refusals
# ===========================================================================


def test_pseudobulk_rejects_an_empty_cell_type_or_state_label(tmp_path):
    h5ad, config = _paired_fixture(tmp_path)
    provider = PB.PseudobulkAggregationProvider()
    result = provider.aggregate(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pb_out"),
        cell_ids_subset=tuple(_PAIRED_CELL_IDS),
        cell_type_or_state_label="",
    )
    assert result.status == "failed"
    assert "cell_type_or_state_label must be a non-empty string" in result.reason


def test_pseudobulk_rejects_an_empty_cell_ids_subset(tmp_path):
    h5ad, config = _paired_fixture(tmp_path)
    provider = PB.PseudobulkAggregationProvider()
    result = provider.aggregate(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pb_out"),
        cell_ids_subset=(),
        cell_type_or_state_label="monocyte",
    )
    assert result.status == "failed"
    assert "cell_ids_subset is empty" in result.reason


def test_pseudobulk_rejects_duplicate_cell_ids_in_subset(tmp_path):
    h5ad, config = _paired_fixture(tmp_path)
    provider = PB.PseudobulkAggregationProvider()
    result = provider.aggregate(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pb_out"),
        cell_ids_subset=("CELL0", "CELL0"),
        cell_type_or_state_label="monocyte",
    )
    assert result.status == "failed"
    assert "duplicate" in result.reason


def test_pseudobulk_rejects_non_raw_counts_semantics(tmp_path):
    h5ad, config = _paired_fixture(tmp_path, counts_semantics="normalized_counts")
    provider = PB.PseudobulkAggregationProvider()
    result = provider.aggregate(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pb_out"),
        cell_ids_subset=tuple(_PAIRED_CELL_IDS),
        cell_type_or_state_label="monocyte",
    )
    assert result.status == "failed"
    assert "requires 'raw_counts' explicitly" in result.reason


def test_pseudobulk_rejects_non_integer_counts(tmp_path):
    X = _PAIRED_X.copy()
    X[0, 0] = 1.5
    h5ad = _make_h5ad(
        tmp_path,
        X=X,
        cell_ids=_PAIRED_CELL_IDS,
        sample_ids=_PAIRED_SAMPLE,
        donor_ids=_PAIRED_DONOR,
        conditions=_PAIRED_CONDITION,
    )
    config = _base_config(h5ad)
    provider = PB.PseudobulkAggregationProvider()
    result = provider.aggregate(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pb_out"),
        cell_ids_subset=tuple(_PAIRED_CELL_IDS),
        cell_type_or_state_label="monocyte",
    )
    assert result.status == "failed"
    assert "non-integer value(s)" in result.reason


def test_pseudobulk_rejects_negative_counts(tmp_path):
    X = _PAIRED_X.copy()
    X[0, 0] = -1.0
    h5ad = _make_h5ad(
        tmp_path,
        X=X,
        cell_ids=_PAIRED_CELL_IDS,
        sample_ids=_PAIRED_SAMPLE,
        donor_ids=_PAIRED_DONOR,
        conditions=_PAIRED_CONDITION,
    )
    config = _base_config(h5ad)
    provider = PB.PseudobulkAggregationProvider()
    result = provider.aggregate(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pb_out"),
        cell_ids_subset=tuple(_PAIRED_CELL_IDS),
        cell_type_or_state_label="monocyte",
    )
    assert result.status == "failed"
    assert "negative value(s)" in result.reason


def test_pseudobulk_rejects_conflicting_donor_metadata_within_a_sample(tmp_path):
    donors = list(_PAIRED_DONOR)
    donors[1] = "D99"  # CELL1 shares sample_id S0 with CELL0, but a different donor
    h5ad = _make_h5ad(
        tmp_path,
        X=_PAIRED_X,
        cell_ids=_PAIRED_CELL_IDS,
        sample_ids=_PAIRED_SAMPLE,
        donor_ids=donors,
        conditions=_PAIRED_CONDITION,
    )
    config = _base_config(h5ad)
    provider = PB.PseudobulkAggregationProvider()
    result = provider.aggregate(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pb_out"),
        cell_ids_subset=tuple(_PAIRED_CELL_IDS),
        cell_type_or_state_label="monocyte",
        condition_column="condition",
    )
    assert result.status == "failed"
    assert "conflicting" in result.reason and "donor_id" in result.reason


def test_pseudobulk_rejects_conflicting_condition_metadata_within_a_sample(tmp_path):
    conditions = list(_PAIRED_CONDITION)
    conditions[1] = "disease"  # CELL1 shares sample_id S0 with CELL0, different condition
    h5ad = _make_h5ad(
        tmp_path,
        X=_PAIRED_X,
        cell_ids=_PAIRED_CELL_IDS,
        sample_ids=_PAIRED_SAMPLE,
        donor_ids=_PAIRED_DONOR,
        conditions=conditions,
    )
    config = _base_config(h5ad)
    provider = PB.PseudobulkAggregationProvider()
    result = provider.aggregate(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pb_out"),
        cell_ids_subset=tuple(_PAIRED_CELL_IDS),
        cell_type_or_state_label="monocyte",
        condition_column="condition",
    )
    assert result.status == "failed"
    assert "conflicting" in result.reason and "condition" in result.reason


# ===========================================================================
# F2 (review:20261014_100436) -- type/null validation must run
# BEFORE any computation/npz write, and a caller-provided-but-missing
# condition/visit value must be refused, never silently written as null or
# crash the conflict check's sorted() call.
# ===========================================================================


def test_pseudobulk_rejects_a_non_string_cell_type_or_state_label_without_orphan_npz(tmp_path):
    h5ad, config = _paired_fixture(tmp_path)
    provider = PB.PseudobulkAggregationProvider()
    out_dir = tmp_path / "pb_out"
    result = provider.aggregate(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(out_dir),
        cell_ids_subset=tuple(_PAIRED_CELL_IDS),
        cell_type_or_state_label=True,
        condition_column="condition",
    )
    assert result.status == "failed"
    assert "cell_type_or_state_label must be a non-empty string" in result.reason
    # A readable failed manifest must exist, and NO orphan npz.
    assert (out_dir / "PSEUDOBULK_RESULT.json").exists()
    assert not (out_dir / "pseudobulk_counts.npz").exists()
    manifest = PB.read_pseudobulk_result(str(out_dir / "PSEUDOBULK_RESULT.json"))
    assert manifest["status"] == "failed"


def test_pseudobulk_rejects_an_all_missing_condition_value_in_one_sample(tmp_path):
    """Per F2: `condition_column` explicitly provided, but every selected
    cell of one sample has a missing value -- must be refused outright, not
    silently written as `condition_value: null`.
    """
    conditions = list(_PAIRED_CONDITION)
    conditions[0] = None
    conditions[1] = None  # BOTH of S0's cells now missing their condition
    h5ad = _make_h5ad(
        tmp_path,
        X=_PAIRED_X,
        cell_ids=_PAIRED_CELL_IDS,
        sample_ids=_PAIRED_SAMPLE,
        donor_ids=_PAIRED_DONOR,
        conditions=conditions,
    )
    config = _base_config(h5ad)
    provider = PB.PseudobulkAggregationProvider()
    out_dir = tmp_path / "pb_out"
    result = provider.aggregate(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(out_dir),
        cell_ids_subset=tuple(_PAIRED_CELL_IDS),
        cell_type_or_state_label="monocyte",
        condition_column="condition",
    )
    assert result.status == "failed"
    assert "missing/non-string value" in result.reason
    assert not (out_dir / "pseudobulk_counts.npz").exists()


def test_pseudobulk_rejects_a_mixed_missing_and_present_condition_without_crashing(tmp_path):
    """Per F2: ONE cell missing, the other cell of the same sample present
    -- must be refused with a structured reason, never raise an uncaught
    `TypeError` from `sorted({None, 'healthy'})`.
    """
    conditions = list(_PAIRED_CONDITION)
    conditions[0] = None  # only CELL0; CELL1 (same sample S0) still 'healthy'
    h5ad = _make_h5ad(
        tmp_path,
        X=_PAIRED_X,
        cell_ids=_PAIRED_CELL_IDS,
        sample_ids=_PAIRED_SAMPLE,
        donor_ids=_PAIRED_DONOR,
        conditions=conditions,
    )
    config = _base_config(h5ad)
    provider = PB.PseudobulkAggregationProvider()
    out_dir = tmp_path / "pb_out"
    result = provider.aggregate(  # must not raise
        config,
        _fake_config_path(tmp_path),
        out_dir=str(out_dir),
        cell_ids_subset=tuple(_PAIRED_CELL_IDS),
        cell_type_or_state_label="monocyte",
        condition_column="condition",
    )
    assert result.status == "failed"
    assert "missing/non-string value" in result.reason
    manifest = PB.read_pseudobulk_result(str(out_dir / "PSEUDOBULK_RESULT.json"))
    assert manifest["status"] == "failed"


def test_pseudobulk_rejects_a_non_string_condition_column_argument(tmp_path):
    h5ad, config = _paired_fixture(tmp_path)
    provider = PB.PseudobulkAggregationProvider()
    result = provider.aggregate(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pb_out"),
        cell_ids_subset=tuple(_PAIRED_CELL_IDS),
        cell_type_or_state_label="monocyte",
        condition_column=True,
    )
    assert result.status == "failed"
    assert "condition_column must be None or a non-empty string" in result.reason


def test_mutation_kill_pseudobulk_write_pre_validation_catches_bad_label_before_npz(
    monkeypatch, tmp_path
):
    """Proves `_write_pseudobulk_result`'s pre-npz validation is itself a
    real guard: disabling `_validate_pseudobulk_result_core` entirely (as
    if the writer never checked anything before writing the npz) and
    forcing a structurally-invalid 'complete' result through `finish()`
    directly would leave an npz with no valid manifest -- the production
    entry point (`aggregate`) never reaches that state because of the
    entry-level type check proven above, so this isolates the writer's
    OWN pre-write guard as a second, independent line of defense.
    """
    import dataclasses

    h5ad, config = _paired_fixture(tmp_path)
    provider = PB.PseudobulkAggregationProvider()
    out_dir = tmp_path / "pb_out"
    result = provider.aggregate(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(out_dir),
        cell_ids_subset=tuple(_PAIRED_CELL_IDS),
        cell_type_or_state_label="monocyte",
        condition_column="condition",
    )
    assert result.status == "complete", result.reason
    counts = PB.load_pseudobulk_counts(
        PB.read_pseudobulk_result(str(out_dir / "PSEUDOBULK_RESULT.json"))
    )

    broken = dataclasses.replace(
        result, cell_type_or_state_label=True, aggregation_path="", aggregation_sha256=""
    )
    out_dir2 = tmp_path / "pb_out2"
    os.makedirs(out_dir2)
    with pytest.raises(ValueError, match="refusing to write an invalid"):
        PB._write_pseudobulk_result(str(out_dir2), broken, counts=counts)
    assert not (out_dir2 / "pseudobulk_counts.npz").exists()  # no orphan npz

    monkeypatch.setattr(PB, "_validate_pseudobulk_result_core", lambda payload: [])
    PB._write_pseudobulk_result(
        str(out_dir2), broken, counts=counts
    )  # with the mutant: an npz IS written
    assert (out_dir2 / "pseudobulk_counts.npz").exists()


# ===========================================================================
# Strict write/read contract -- tamper rejection
# ===========================================================================


def test_read_pseudobulk_result_rejects_a_tampered_conservation_claim(tmp_path):
    import json

    h5ad, config = _paired_fixture(tmp_path)
    provider = PB.PseudobulkAggregationProvider()
    result = provider.aggregate(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pb_out"),
        cell_ids_subset=tuple(_PAIRED_CELL_IDS),
        cell_type_or_state_label="monocyte",
        condition_column="condition",
    )
    assert result.status == "complete", result.reason
    path = tmp_path / "pb_out" / "PSEUDOBULK_RESULT.json"

    payload = json.loads(path.read_text())
    payload["total_counts_aggregated"] = payload["total_counts_aggregated"] + 1000.0
    path.write_text(json.dumps(payload))
    with pytest.raises(SystemExit, match="failed validation"):
        PB.read_pseudobulk_result(str(path))


def test_read_pseudobulk_result_rejects_a_tampered_sample_row_count(tmp_path):
    import json

    h5ad, config = _paired_fixture(tmp_path)
    provider = PB.PseudobulkAggregationProvider()
    result = provider.aggregate(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pb_out"),
        cell_ids_subset=tuple(_PAIRED_CELL_IDS),
        cell_type_or_state_label="monocyte",
        condition_column="condition",
    )
    assert result.status == "complete", result.reason
    path = tmp_path / "pb_out" / "PSEUDOBULK_RESULT.json"

    payload = json.loads(path.read_text())
    payload["n_pseudobulk_samples"] = 999
    path.write_text(json.dumps(payload))
    with pytest.raises(SystemExit, match="failed validation"):
        PB.read_pseudobulk_result(str(path))


def test_load_pseudobulk_counts_rejects_a_tampered_npz_file(tmp_path):
    h5ad, config = _paired_fixture(tmp_path)
    provider = PB.PseudobulkAggregationProvider()
    result = provider.aggregate(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pb_out"),
        cell_ids_subset=tuple(_PAIRED_CELL_IDS),
        cell_type_or_state_label="monocyte",
        condition_column="condition",
    )
    assert result.status == "complete", result.reason
    manifest = PB.read_pseudobulk_result(str(tmp_path / "pb_out" / "PSEUDOBULK_RESULT.json"))
    # overwrite the npz with different content after the manifest was written
    np.savez(manifest["aggregation_path"], counts=np.zeros((4, 4)))
    with pytest.raises(SystemExit, match="does not match"):
        PB.load_pseudobulk_counts(manifest)


# ===========================================================================
# F1 (review:20261014_100436) -- the manifest's LABELS (sample
# order, gene order, donor count, a row's library size, a condition value)
# must be bound to the actual npz content, not merely self-consistent JSON.
# Reproduces the review's own five independent probe cases.
# ===========================================================================


def _complete_manifest_dict(tmp_path):
    import json

    h5ad, config = _paired_fixture(tmp_path)
    provider = PB.PseudobulkAggregationProvider()
    out_dir = tmp_path / "pb_out"
    result = provider.aggregate(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(out_dir),
        cell_ids_subset=tuple(_PAIRED_CELL_IDS),
        cell_type_or_state_label="monocyte",
        condition_column="condition",
        visit_or_state_column="visit",
    )
    assert result.status == "complete", result.reason
    path = out_dir / "PSEUDOBULK_RESULT.json"
    return json.loads(path.read_text()), out_dir


def test_read_pseudobulk_result_rejects_reversed_sample_rows(tmp_path):
    import json

    manifest, out_dir = _complete_manifest_dict(tmp_path)
    manifest["sample_rows"] = list(reversed(manifest["sample_rows"]))
    tampered = out_dir / "tampered.json"
    tampered.write_text(json.dumps(manifest))
    with pytest.raises(SystemExit, match="failed validation"):
        PB.read_pseudobulk_result(str(tampered))


def test_read_pseudobulk_result_rejects_reversed_gene_ids(tmp_path):
    import json

    manifest, out_dir = _complete_manifest_dict(tmp_path)
    manifest["gene_ids"] = list(reversed(manifest["gene_ids"]))
    tampered = out_dir / "tampered.json"
    tampered.write_text(json.dumps(manifest))
    with pytest.raises(SystemExit, match="failed validation"):
        PB.read_pseudobulk_result(str(tampered))


def test_read_pseudobulk_result_rejects_a_wrong_distinct_donor_count(tmp_path):
    import json

    manifest, out_dir = _complete_manifest_dict(tmp_path)
    manifest["n_distinct_donors"] = 999
    tampered = out_dir / "tampered.json"
    tampered.write_text(json.dumps(manifest))
    with pytest.raises(SystemExit, match="failed validation"):
        PB.read_pseudobulk_result(str(tampered))


def test_read_pseudobulk_result_rejects_a_wrong_library_size(tmp_path):
    import json

    manifest, out_dir = _complete_manifest_dict(tmp_path)
    manifest["sample_rows"][0]["library_size"] = 999.0
    tampered = out_dir / "tampered.json"
    tampered.write_text(json.dumps(manifest))
    with pytest.raises(SystemExit, match="failed validation"):
        PB.read_pseudobulk_result(str(tampered))


def test_read_pseudobulk_result_rejects_a_flipped_condition_value(tmp_path):
    import json

    manifest, out_dir = _complete_manifest_dict(tmp_path)
    original = manifest["sample_rows"][0]["condition_value"]
    manifest["sample_rows"][0]["condition_value"] = (
        "disease" if original != "disease" else "healthy"
    )
    tampered = out_dir / "tampered.json"
    tampered.write_text(json.dumps(manifest))
    with pytest.raises(SystemExit, match="failed validation"):
        PB.read_pseudobulk_result(str(tampered))


def test_dataset_identity_carries_content_hashes_not_just_a_config_hash(tmp_path):
    manifest, _ = _complete_manifest_dict(tmp_path)
    ds = manifest["dataset_identity"]
    for key in (
        "selected_raw_counts_sha256",
        "selected_cell_metadata_sha256",
        "gene_ids_sha256",
        "dataset_config_sha256",
    ):
        assert isinstance(ds[key], str) and ds[key]
    # Two independently-named hash fields -- a config hash is NOT a stand-in
    # for a content hash of the actual expression matrix block that was read.
    assert ds["selected_raw_counts_sha256"] != ds["dataset_config_sha256"]


def test_mutation_kill_pseudobulk_identity_sha256_recompute(monkeypatch, tmp_path):
    """Proves `validate_pseudobulk_result`'s `identity_sha256` recompute is
    what actually catches a reordered `sample_rows` -- first shows the real
    reader rejects it, then disables ONLY the recompute (as if it always
    agreed with whatever hash is on disk) and confirms the SAME tampered
    file is then silently accepted.
    """
    import json

    manifest, out_dir = _complete_manifest_dict(tmp_path)
    manifest["sample_rows"] = list(reversed(manifest["sample_rows"]))
    tampered = out_dir / "tampered.json"
    tampered.write_text(json.dumps(manifest))

    with pytest.raises(SystemExit, match="failed validation"):
        PB.read_pseudobulk_result(str(tampered))  # real guard: rejects

    monkeypatch.setattr(
        PB, "_pseudobulk_identity_sha256", lambda *a, **k: manifest["identity_sha256"]
    )
    reread = PB.read_pseudobulk_result(str(tampered))  # with the mutant: silently accepted
    assert reread["sample_rows"] == manifest["sample_rows"]


def test_mutation_kill_load_pseudobulk_counts_library_size_cross_check(tmp_path):
    """Proves `load_pseudobulk_counts` itself (not just `validate_
    pseudobulk_result`) binds the manifest's row labels to the npz array's
    OWN stored identity: a manifest whose `sample_rows` order was swapped
    (and whose `identity_sha256`/`gene_total_counts` were ALSO recomputed to
    stay self-consistent with the new order -- a fully self-consistent JSON
    edit that `identity_sha256` alone cannot catch) is still rejected,
    because the relabeled rows no longer match the UNCHANGED npz's own
    `sample_ids` array (per `review:20261014_114357` R1 -- see
    `test_load_pseudobulk_counts_rejects_an_equal_marginals_row_swap` for
    the specific equal-library-size case R1 raised, where a library-size-
    only recompute cannot distinguish the swap at all).
    """
    manifest, out_dir = _complete_manifest_dict(tmp_path)
    swapped_rows = list(reversed(manifest["sample_rows"]))
    swapped_manifest = dict(manifest)
    swapped_manifest["sample_rows"] = swapped_rows
    swapped_manifest["identity_sha256"] = PB._pseudobulk_identity_sha256(
        manifest["gene_ids"],
        manifest["gene_total_counts"],
        swapped_rows,
        manifest["n_pseudobulk_samples"],
        manifest["n_distinct_donors"],
        manifest["total_counts_selected_cells"],
        manifest["total_counts_aggregated"],
    )
    # The npz file itself (and its hash) are UNCHANGED -- only the JSON
    # labels were reordered, and self-consistently rehashed.
    with pytest.raises(SystemExit, match="do not match the manifest's sample_rows order"):
        PB.load_pseudobulk_counts(swapped_manifest)


# ===========================================================================
# R1 (review:20261014_114357) -- a row/column-total recompute is a
# MANY-TO-ONE projection of content: an equal-marginals fixture ([[1, 9],
# [9, 1]], both row/column sums = 10) must still have `sample_rows`/
# `gene_ids` order and per-row metadata bound to the npz's OWN stored
# identity, not only to totals recomputed from `counts`.
# ===========================================================================


def _equal_marginals_fixture(tmp_path):
    """The review's own `[[1, 9], [9, 1]]` fixture: S0/S1 have IDENTICAL
    library sizes (10 each) and G0/G1 have IDENTICAL column totals (10
    each) -- a totals-only recompute cannot distinguish a row or column
    swap, or a condition relabel, from the real content at all.
    """
    import anndata as ad
    import pandas as pd

    X = np.array([[1, 9], [9, 1]], dtype="float32")
    obs = pd.DataFrame(
        dict(
            cell_barcode=["C0", "C1"],
            sample_id=["S0", "S1"],
            donor_id=["D0", "D1"],
            condition=["healthy", "disease"],
        ),
        index=["idx0", "idx1"],
    )
    var = pd.DataFrame(index=["G0", "G1"])
    path = tmp_path / "equal_marginals.h5ad"
    ad.AnnData(X=X, obs=obs, var=var).write_h5ad(path)
    config = _base_config(str(path))
    provider = PB.PseudobulkAggregationProvider()
    out_dir = tmp_path / "pb_out"
    result = provider.aggregate(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(out_dir),
        cell_ids_subset=("C0", "C1"),
        cell_type_or_state_label="synthetic",
        condition_column="condition",
    )
    assert result.status == "complete", result.reason
    import json

    manifest = json.loads((out_dir / "PSEUDOBULK_RESULT.json").read_text())
    return manifest, out_dir


def _rehash(manifest: dict) -> dict:
    """Self-consistently recomputes `identity_sha256` after a manifest edit
    -- the SAME thing a forger who read `_pseudobulk_identity_sha256`'s own
    definition could do; proves the fix does not merely rely on the forger
    not bothering to update the hash.
    """
    manifest = dict(manifest)
    manifest["identity_sha256"] = PB._pseudobulk_identity_sha256(
        manifest["gene_ids"],
        manifest["gene_total_counts"],
        manifest["sample_rows"],
        manifest["n_pseudobulk_samples"],
        manifest["n_distinct_donors"],
        manifest["total_counts_selected_cells"],
        manifest["total_counts_aggregated"],
    )
    return manifest


def test_load_pseudobulk_counts_rejects_an_equal_marginals_row_swap(tmp_path):
    manifest, out_dir = _equal_marginals_fixture(tmp_path)
    swapped = dict(manifest)
    swapped["sample_rows"] = list(reversed(manifest["sample_rows"]))
    swapped = _rehash(swapped)  # fully self-consistent: identity_sha256 updated too
    # The npz bytes are UNCHANGED -- row/column totals (10/10) agree either
    # way, so only the npz's own stored sample_ids can catch this.
    with pytest.raises(SystemExit, match="do not match the manifest's sample_rows order"):
        PB.load_pseudobulk_counts(swapped)


def test_load_pseudobulk_counts_rejects_an_equal_marginals_gene_swap(tmp_path):
    manifest, out_dir = _equal_marginals_fixture(tmp_path)
    swapped = dict(manifest)
    swapped["gene_ids"] = list(reversed(manifest["gene_ids"]))
    swapped["gene_total_counts"] = list(reversed(manifest["gene_total_counts"]))
    swapped = _rehash(swapped)
    with pytest.raises(SystemExit, match="gene axis is not the one actually bound"):
        PB.load_pseudobulk_counts(swapped)


def test_load_pseudobulk_counts_rejects_an_equal_marginals_condition_flip(tmp_path):
    manifest, out_dir = _equal_marginals_fixture(tmp_path)
    flipped = dict(manifest)
    flipped["sample_rows"] = [dict(r) for r in manifest["sample_rows"]]
    original = flipped["sample_rows"][0]["condition_value"]
    flipped["sample_rows"][0]["condition_value"] = "disease" if original != "disease" else "healthy"
    flipped = _rehash(flipped)
    # Row order/gene order/totals are all untouched -- only the condition
    # label changed, self-consistently rehashed.
    with pytest.raises(SystemExit, match="per-row metadata is not the metadata actually bound"):
        PB.load_pseudobulk_counts(flipped)


def test_mutation_kill_load_pseudobulk_counts_npz_identity_binding(monkeypatch, tmp_path):
    """Proves `_npz_identity_binding_problems` (the npz-embedded sample_ids/
    gene_ids/identity-snapshot comparison) is itself the load-bearing guard
    for the equal-marginals case: the real check rejects the row swap;
    disabling ONLY that comparison (as if `load_pseudobulk_counts` still
    recomputed row/column totals but never compared the npz's own stored
    identity arrays) lets the SAME swapped manifest load silently --
    because the swap's row/column totals (10/10 either way) agree with the
    totals-only checks that remain, reproducing R1's exact gap.
    """
    manifest, out_dir = _equal_marginals_fixture(tmp_path)
    swapped = dict(manifest)
    swapped["sample_rows"] = list(reversed(manifest["sample_rows"]))
    swapped = _rehash(swapped)

    with pytest.raises(SystemExit, match="do not match the manifest's sample_rows order"):
        PB.load_pseudobulk_counts(swapped)  # real guard: rejects

    monkeypatch.setattr(PB, "_npz_identity_binding_problems", lambda *a, **k: [])
    counts = PB.load_pseudobulk_counts(swapped)  # with the mutant: silently accepted
    # The first row is now labeled "S1" but is still the array's row 0
    # ([1, 9], the content that was actually D0/healthy's) -- the mutant
    # let a mislabeled row through.
    np.testing.assert_array_equal(counts[0], [1, 9])


def test_read_pseudobulk_result_rejects_a_non_int_n_pseudobulk_samples(tmp_path):
    """Per `review:20261014_114357` R2 (generalized to
    `n_pseudobulk_samples`, the same `isinstance(x, int) and ... != ...`
    shape): `read_pseudobulk_result` ALONE (no npz access) must reject a
    non-int `n_pseudobulk_samples`, not silently skip the consistency
    check the way the prior shape did.
    """
    manifest, out_dir = _complete_manifest_dict(tmp_path)
    tampered = dict(manifest)
    tampered["n_pseudobulk_samples"] = "two"
    import json

    path = out_dir / "tampered_n_samples.json"
    path.write_text(json.dumps(tampered))
    with pytest.raises(SystemExit, match="n_pseudobulk_samples must be an int"):
        PB.read_pseudobulk_result(str(path))


def test_read_pseudobulk_result_rejects_a_non_int_n_distinct_donors(tmp_path):
    """Same as above for `n_distinct_donors`, exercised through
    `read_pseudobulk_result` alone -- `load_pseudobulk_counts`'s own
    `isinstance` check already caught this (per the review's R2 note), but
    the READER entry point (no npz access) must not disagree with it.
    """
    manifest, out_dir = _complete_manifest_dict(tmp_path)
    tampered = dict(manifest)
    tampered["n_distinct_donors"] = None
    import json

    path = out_dir / "tampered_n_donors.json"
    path.write_text(json.dumps(tampered))
    with pytest.raises(SystemExit, match="n_distinct_donors must be an int"):
        PB.read_pseudobulk_result(str(path))


# ===========================================================================
# Mutation-kill: each factored guard function is proven to be load-bearing,
# not merely present and ignored (DEVELOPMENT.md rules 2 and 6).
# ===========================================================================


def test_mutation_kill_pseudobulk_raw_counts_validity_guard(monkeypatch):
    """Positive unit check plus a real reproduction: disabling the guard
    lets the non-integer-counts fixture (which must normally fail) report
    'complete' instead.
    """
    assert PB._raw_counts_validity_problem(np.array([[1.0, 2.0]])) is None
    assert "negative" in PB._raw_counts_validity_problem(np.array([[-1.0, 2.0]]))
    assert "non-integer" in PB._raw_counts_validity_problem(np.array([[1.5, 2.0]]))


def test_mutation_kill_pseudobulk_raw_counts_validity_guard_end_to_end(monkeypatch, tmp_path):
    monkeypatch.setattr(PB, "_raw_counts_validity_problem", lambda raw: None)
    X = _PAIRED_X.copy()
    X[0, 0] = 1.5
    h5ad = _make_h5ad(
        tmp_path,
        X=X,
        cell_ids=_PAIRED_CELL_IDS,
        sample_ids=_PAIRED_SAMPLE,
        donor_ids=_PAIRED_DONOR,
        conditions=_PAIRED_CONDITION,
    )
    config = _base_config(h5ad)
    provider = PB.PseudobulkAggregationProvider()
    result = provider.aggregate(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pb_out"),
        cell_ids_subset=tuple(_PAIRED_CELL_IDS),
        cell_type_or_state_label="monocyte",
    )
    assert result.status == "complete"  # with the mutant: non-integer counts no longer caught


def test_mutation_kill_pseudobulk_sample_metadata_conflict_guard_end_to_end(monkeypatch, tmp_path):
    monkeypatch.setattr(PB, "_sample_metadata_conflict_problem", lambda *a, **k: None)
    donors = list(_PAIRED_DONOR)
    donors[1] = "D99"
    h5ad = _make_h5ad(
        tmp_path,
        X=_PAIRED_X,
        cell_ids=_PAIRED_CELL_IDS,
        sample_ids=_PAIRED_SAMPLE,
        donor_ids=donors,
        conditions=_PAIRED_CONDITION,
    )
    config = _base_config(h5ad)
    provider = PB.PseudobulkAggregationProvider()
    result = provider.aggregate(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pb_out"),
        cell_ids_subset=tuple(_PAIRED_CELL_IDS),
        cell_type_or_state_label="monocyte",
        condition_column="condition",
    )
    assert result.status == "complete"  # with the mutant: conflicting donor IDs no longer caught


def test_mutation_kill_pseudobulk_counts_conservation_guard(monkeypatch, tmp_path):
    """The conservation mismatch cannot occur through real input (the
    grouping sums every selected row exactly once) -- proven instead by (a)
    a direct unit check of the pure function, and (b) proving its result
    actually controls `aggregate()`'s outcome: forcing it to always report
    a conflict makes an otherwise-valid aggregation fail, so the check is a
    real guard, not a computed-and-ignored field (DEVELOPMENT.md rule 6).
    """
    assert PB._counts_conservation_problem(10.0, 10.0) is None
    assert "conservation" in PB._counts_conservation_problem(10.0, 5.0)

    monkeypatch.setattr(
        PB, "_counts_conservation_problem", lambda sel, agg: "forced conservation failure"
    )
    h5ad, config = _paired_fixture(tmp_path)
    provider = PB.PseudobulkAggregationProvider()
    result = provider.aggregate(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pb_out"),
        cell_ids_subset=tuple(_PAIRED_CELL_IDS),
        cell_type_or_state_label="monocyte",
        condition_column="condition",
    )
    assert result.status == "failed"
    assert result.reason == "forced conservation failure"


def test_pseudobulk_result_to_dict_is_json_serializable():
    import json

    result = PB.PseudobulkResult(
        status="complete",
        cell_type_or_state_label="monocyte",
        gene_ids=("G1", "G2"),
        sample_rows=(
            {
                "sample_id": "S0",
                "donor_or_animal_id": "D0",
                "condition_value": "healthy",
                "visit_or_state_value": "",
                "n_cells": 2,
                "library_size": 4.0,
                "cell_ids": ["C0", "C1"],
            },
        ),
        n_pseudobulk_samples=1,
        n_distinct_donors=1,
        total_counts_selected_cells=4.0,
        total_counts_aggregated=4.0,
        dataset_identity={"dataset_id": "x"},
        backend="raw_count_sum_pseudobulk",
        backend_version="1",
        source_code_revision_or_hash="abc",
        dependency_versions={"python": "3.10"},
    )
    json.dumps(result.to_dict())  # must not raise
