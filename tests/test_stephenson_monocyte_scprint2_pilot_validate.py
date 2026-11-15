import os
import sys

import numpy as np
import pandas as pd
import pytest

RESULTS_DIR = os.path.join(
    os.path.dirname(__file__), "..", "scripts", "34_scfm_discovery_framework", "results"
)
sys.path.insert(0, RESULTS_DIR)

ad = pytest.importorskip("anndata")
import stephenson_monocyte_scprint2_pilot_validate as VAL  # noqa: E402


def _id_num(cell_id):
    """Extracts the trailing integer from a "cellN" id, falling back to 0 for ids
    outside that scheme (e.g. the "alien*" fixtures used for unrelated-id scenarios,
    where the derived identity/QC values are irrelevant to what those scenarios test).
    """
    try:
        return int(cell_id[len("cell") :])
    except ValueError:
        return 0


def _decontx_value_for_id(cell_id):
    """Derives a QC value from the cell id ITSELF (e.g. "cell2" -> 0.2), not from
    list position -- so a fixture with reordered obs_names still carries the
    correct-by-identity value for each cell, exactly like a real producer/consumer
    boundary must (never assume row order correspondence)."""
    return 0.1 * _id_num(cell_id)


def _patient_id_for_id(cell_id):
    return f"patient{_id_num(cell_id)}"


def _sample_id_for_id(cell_id):
    return f"sample{_id_num(cell_id)}"


def _make_expected_obs(n=3, required_columns=True):
    ids = [f"cell{i}" for i in range(n)]
    obs = pd.DataFrame(index=ids)
    if required_columns:
        obs["gate1_flag_would_delete"] = [False] * n
        obs["decontx_estimated_contamination"] = [_decontx_value_for_id(i) for i in ids]
        obs["patient_id"] = [_patient_id_for_id(i) for i in ids]
        obs["Sample_ID"] = [_sample_id_for_id(i) for i in ids]
    return obs


def _make_embedded_adata(obs_names, obsm=None, required_columns=True, column_values=None):
    """`column_values`, when given, overrides specific columns' VALUES (still one
    value per row of `obs_names`, in that order) while keeping all four required
    columns present -- used to simulate a present-but-misrouted annotation without
    having to also specify the other columns' correct-by-id values by hand."""
    n = len(obs_names)
    obs = pd.DataFrame(index=obs_names)
    if required_columns:
        overrides = column_values or {}
        obs["gate1_flag_would_delete"] = overrides.get("gate1_flag_would_delete", [False] * n)
        obs["decontx_estimated_contamination"] = overrides.get(
            "decontx_estimated_contamination", [_decontx_value_for_id(i) for i in obs_names]
        )
        obs["patient_id"] = overrides.get("patient_id", [_patient_id_for_id(i) for i in obs_names])
        obs["Sample_ID"] = overrides.get("Sample_ID", [_sample_id_for_id(i) for i in obs_names])
    var = pd.DataFrame(index=["g1", "g2"])
    adata = ad.AnnData(X=np.zeros((n, 2), dtype="float32"), obs=obs, var=var)
    if obsm:
        for k, v in obsm.items():
            adata.obsm[k] = v
    return adata


def test_find_embedding_keys_matches_known_prefixes():
    assert VAL.find_embedding_keys(["X_scprint_umap", "X_pca", "scprint_emb", "other"]) == [
        "X_scprint_umap",
        "scprint_emb",
    ]


def test_validate_embedding_output_passes_for_a_well_formed_result():
    ids = ["cell0", "cell1", "cell2"]
    adata = _make_embedded_adata(ids, obsm={"scprint_emb": np.ones((3, 4))})
    report, hard_failures = VAL.validate_embedding_output(adata, _make_expected_obs(3))
    assert hard_failures == []
    assert report["n_obs_difference"] == 0
    assert report["ids_match_as_sets"] is True
    assert report["order_matches_expected"] is True
    assert report["required_obs_columns"]["gate1_flag_would_delete"]["present_in_output"] is True
    assert report["required_obs_columns"]["patient_id"]["present_in_output"] is True
    assert report["required_obs_columns"]["Sample_ID"]["present_in_output"] is True


def test_validate_embedding_output_fails_hard_when_only_the_umap_key_is_present():
    """Per review:20261016_152556: the installed producer
    (scprint2/tasks/cell_emb.py:214,227-251) writes the UMAP to `X_scprint_umap` and
    the actual primary embedding to `scprint_emb` -- a file with only the former must
    not pass as if it had the real embedding.
    """
    ids = ["cell0", "cell1", "cell2"]
    adata = _make_embedded_adata(ids, obsm={"X_scprint_umap": np.ones((3, 2))})
    report, hard_failures = VAL.validate_embedding_output(adata, _make_expected_obs(3))
    assert any("scprint_emb" in f and "not found" in f for f in hard_failures)
    assert report["embedding_keys_found_all"] == ["X_scprint_umap"]


def test_validate_embedding_output_fails_hard_on_missing_cell_ids():
    """Fewer cells than expected: 1 of 3 (the review's 'one_of_three_cells' probe)."""
    adata = _make_embedded_adata(["cell0"], obsm={"scprint_emb": np.ones((1, 4))})
    report, hard_failures = VAL.validate_embedding_output(adata, _make_expected_obs(3))
    assert any("missing from the output" in f for f in hard_failures)
    assert report["n_missing_cell_ids"] == 2
    assert set(report["missing_cell_ids_sample"]) == {"cell1", "cell2"}


def test_validate_embedding_output_fails_hard_on_extra_cell_ids():
    """More cells than expected (the review's 'extra_cell' probe)."""
    ids = ["cell0", "cell1", "cell2", "cell3"]
    adata = _make_embedded_adata(ids, obsm={"scprint_emb": np.ones((4, 4))})
    report, hard_failures = VAL.validate_embedding_output(adata, _make_expected_obs(3))
    assert any("not in the expected" in f for f in hard_failures)
    assert report["n_extra_cell_ids"] == 1
    assert report["missing_cell_ids_sample"] == []


def test_validate_embedding_output_fails_hard_on_unrelated_ids_even_with_the_right_count():
    """The review's 'unrelated_ids' probe: same COUNT, completely different ids --
    v2's expected_n_obs-only contract passed this; the id-aware v3 contract must not.
    """
    adata = _make_embedded_adata(
        ["alien0", "alien1", "alien2"], obsm={"scprint_emb": np.ones((3, 4))}
    )
    report, hard_failures = VAL.validate_embedding_output(adata, _make_expected_obs(3))
    assert any("not in the expected" in f for f in hard_failures)
    assert any("missing from the output" in f for f in hard_failures)
    assert report["n_extra_cell_ids"] == 3
    assert report["n_missing_cell_ids"] == 3


def test_validate_embedding_output_fails_hard_on_duplicate_ids_even_with_the_right_count():
    """The review's 'duplicate_ids' probe: same COUNT, one id repeated 3x."""
    adata = _make_embedded_adata(["cell0", "cell0", "cell0"], obsm={"scprint_emb": np.ones((3, 4))})
    report, hard_failures = VAL.validate_embedding_output(adata, _make_expected_obs(3))
    assert any("duplicate cell id" in f for f in hard_failures)
    assert report["n_duplicate_cell_ids_in_output"] == 2


def test_validate_embedding_output_fails_hard_on_non_finite_embedding_values():
    emb = np.ones((3, 4))
    emb[0, 0] = float("nan")
    adata = _make_embedded_adata(["cell0", "cell1", "cell2"], obsm={"scprint_emb": emb})
    report, hard_failures = VAL.validate_embedding_output(adata, _make_expected_obs(3))
    assert any("non-finite" in f for f in hard_failures)
    assert report["non_finite_values_by_embedding_key"]["scprint_emb"] == 1


def test_validate_embedding_output_fails_hard_on_zero_columns():
    """The review's 'zero_dimensions' probe: an (N, 0) embedding must not pass the
    finiteness check by vacuously having nothing non-finite to find.
    """
    adata = _make_embedded_adata(["cell0", "cell1", "cell2"], obsm={"scprint_emb": np.ones((3, 0))})
    report, hard_failures = VAL.validate_embedding_output(adata, _make_expected_obs(3))
    assert any("invalid shape" in f for f in hard_failures)


def test_validate_embedding_output_fails_hard_on_a_string_embedding():
    """The review's 'string_embedding' probe."""
    emb = np.full((3, 4), "not_an_embedding")
    adata = _make_embedded_adata(["cell0", "cell1", "cell2"], obsm={"scprint_emb": emb})
    report, hard_failures = VAL.validate_embedding_output(adata, _make_expected_obs(3))
    assert any("non-numeric dtype" in f for f in hard_failures)


def test_validate_embedding_output_fails_hard_on_zero_cells():
    adata = _make_embedded_adata([], obsm={"scprint_emb": np.zeros((0, 4))})
    report, hard_failures = VAL.validate_embedding_output(adata, _make_expected_obs(3))
    assert any("0 cells" in f for f in hard_failures)


def test_validate_embedding_output_fails_hard_on_a_missing_required_column_even_when_ids_fully_match():
    """Per review:20261016_202903 'R2': the prior version reported a
    missing-but-"recoverable" column as a technical PASS without ever writing or
    persisting any recovered values anywhere -- the review's own `persisted_missing_qc`
    probe reopened the real output file and found the columns still absent. There is
    no recovery branch any more: missing is always a hard failure, period.
    """
    ids = ["cell0", "cell1", "cell2"]
    adata = _make_embedded_adata(ids, obsm={"scprint_emb": np.ones((3, 4))}, required_columns=False)
    report, hard_failures = VAL.validate_embedding_output(adata, _make_expected_obs(3))
    assert any(
        "required obs column 'gate1_flag_would_delete' is missing" in f for f in hard_failures
    )
    assert any("required obs column 'patient_id' is missing" in f for f in hard_failures)
    assert report["required_obs_columns"]["gate1_flag_would_delete"]["present_in_output"] is False
    assert (
        "recovered_from_expected_manifest"
        not in report["required_obs_columns"]["gate1_flag_would_delete"]
    )


def test_validate_embedding_output_fails_hard_on_a_missing_required_column_when_ids_do_not_fully_match():
    """Missing + id mismatch is still just "missing" -- the fail-closed rule does not
    depend on whether the ids otherwise match (see the no-recovery test above)."""
    adata = _make_embedded_adata(
        ["cell0"], obsm={"scprint_emb": np.ones((1, 4))}, required_columns=False
    )
    report, hard_failures = VAL.validate_embedding_output(adata, _make_expected_obs(3))
    assert any(
        "required obs column 'gate1_flag_would_delete' is missing" in f for f in hard_failures
    )


def test_validate_embedding_output_fails_hard_on_a_present_but_misrouted_qc_column():
    """Per review:20261016_152556: donor/sample/visit/QC correspondence
    must be checked BY ID, not position -- a present QC column whose values diverge
    from the frozen prepared-input manifest for a shared id is a misrouted annotation.
    """
    ids = ["cell0", "cell1", "cell2"]
    adata = _make_embedded_adata(
        ids,
        obsm={"scprint_emb": np.ones((3, 4))},
        column_values={"gate1_flag_would_delete": [True, True, True]},  # expected is all False
    )
    report, hard_failures = VAL.validate_embedding_output(adata, _make_expected_obs(3))
    assert any("diverges from the frozen prepared-input manifest" in f for f in hard_failures)
    assert (
        report["required_obs_columns"]["gate1_flag_would_delete"][
            "n_values_diverging_from_expected_manifest"
        ]
        == 3
    )


@pytest.mark.parametrize("col", ["patient_id", "Sample_ID"])
def test_validate_embedding_output_fails_hard_on_a_present_but_misrouted_identity_column(col):
    """Per review:20261016_202903 'R3': the old validator never checked
    donor/sample identity at all -- the review's independent `swapped_patient_id`/
    `swapped_Sample_ID` probes (same ids, same correct primary embedding, same correct
    QC values, only the identity column's VALUES swapped) both returned
    `hard_failures=[]`. `patient_id`/`Sample_ID` are real columns on
    `stephenson_monocyte_raw.h5ad` (confirmed via a bounded h5py read, 2026-10-16), so
    they are now checked by id exactly like the QC columns.
    """
    ids = ["cell0", "cell1", "cell2"]
    reversed_vals = (
        [_patient_id_for_id(i) for i in ids][::-1]
        if col == "patient_id"
        else [_sample_id_for_id(i) for i in ids][::-1]
    )
    adata = _make_embedded_adata(
        ids,
        obsm={"scprint_emb": np.ones((3, 4))},
        column_values={col: reversed_vals},
    )
    report, hard_failures = VAL.validate_embedding_output(adata, _make_expected_obs(3))
    assert any(
        f"obs[{col!r}] diverges from the frozen prepared-input manifest" in f for f in hard_failures
    )
    assert report["required_obs_columns"][col]["n_values_diverging_from_expected_manifest"] > 0


def test_validate_embedding_output_fails_hard_on_a_missing_identity_column():
    """A real run that somehow dropped `patient_id`/`Sample_ID` (e.g. a future prep-step
    regression) must fail closed exactly like a missing QC column -- these are not an
    optional, QC-only check."""
    ids = ["cell0", "cell1", "cell2"]
    adata = _make_embedded_adata(ids, obsm={"scprint_emb": np.ones((3, 4))})
    del adata.obs["patient_id"]
    report, hard_failures = VAL.validate_embedding_output(adata, _make_expected_obs(3))
    assert any("required obs column 'patient_id' is missing" in f for f in hard_failures)


def test_validate_embedding_output_records_order_change_without_failing_but_never_assumes_position():
    """A reordered-but-complete output must not fail on order alone -- but the QC
    cross-check above still has to catch a REAL misroute even when it would line up
    if checked positionally. This test only checks the order flag itself.
    """
    ids = ["cell2", "cell0", "cell1"]  # reordered, same set
    adata = _make_embedded_adata(ids, obsm={"scprint_emb": np.ones((3, 4))})
    report, hard_failures = VAL.validate_embedding_output(adata, _make_expected_obs(3))
    assert hard_failures == []
    assert report["ids_match_as_sets"] is True
    assert report["order_matches_expected"] is False


def test_main_writes_report_and_returns_nonzero_on_a_real_hard_failure(tmp_path):
    embedded_path = tmp_path / "embedded.h5ad"
    _make_embedded_adata(
        ["alien0", "alien1", "alien2"], obsm={"scprint_emb": np.ones((3, 4))}
    ).write_h5ad(embedded_path)
    expected_input_path = tmp_path / "expected_input.h5ad"
    ad.AnnData(
        X=np.zeros((3, 1), dtype="float32"),
        obs=_make_expected_obs(3),
        var=pd.DataFrame(index=["g1"]),
    ).write_h5ad(expected_input_path)
    out_json = tmp_path / "report.json"

    rc = VAL.main(str(embedded_path), str(expected_input_path), str(out_json))
    assert rc == 1
    assert out_json.exists()


def test_main_returns_zero_on_a_real_pass(tmp_path):
    ids = ["cell0", "cell1", "cell2"]
    embedded_path = tmp_path / "embedded.h5ad"
    _make_embedded_adata(ids, obsm={"scprint_emb": np.ones((3, 4))}).write_h5ad(embedded_path)
    expected_input_path = tmp_path / "expected_input.h5ad"
    ad.AnnData(
        X=np.zeros((3, 1), dtype="float32"),
        obs=_make_expected_obs(3),
        var=pd.DataFrame(index=["g1"]),
    ).write_h5ad(expected_input_path)
    out_json = tmp_path / "report.json"

    rc = VAL.main(str(embedded_path), str(expected_input_path), str(out_json))
    assert rc == 0


def test_main_fails_when_the_output_file_does_not_exist(tmp_path):
    expected_input_path = tmp_path / "expected_input.h5ad"
    ad.AnnData(
        X=np.zeros((3, 1), dtype="float32"),
        obs=_make_expected_obs(3),
        var=pd.DataFrame(index=["g1"]),
    ).write_h5ad(expected_input_path)
    rc = VAL.main(
        str(tmp_path / "missing.h5ad"), str(expected_input_path), str(tmp_path / "report.json")
    )
    assert rc == 1


def test_main_fails_when_the_expected_input_manifest_does_not_exist(tmp_path):
    embedded_path = tmp_path / "embedded.h5ad"
    _make_embedded_adata(["cell0"], obsm={"scprint_emb": np.ones((1, 4))}).write_h5ad(embedded_path)
    rc = VAL.main(
        str(embedded_path), str(tmp_path / "missing_expected.h5ad"), str(tmp_path / "report.json")
    )
    assert rc == 1
