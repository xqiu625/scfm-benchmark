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
import stephenson_monocyte_scprint2_pilot_prep as PREP  # noqa: E402


def _make_fixture_adata():
    obs = pd.DataFrame(index=["steph_AAA-S1", "steph_BBB-S1", "steph_CCC-S2"])
    var = pd.DataFrame(
        {"feature_types": ["Gene Expression", "Gene Expression", "Antibody Capture"]},
        index=["GENE1", "GENE2", "CD3_protein"],
    )
    X = np.array([[1, 2, 5], [3, 0, 7], [0, 4, 9]], dtype="float32")
    return ad.AnnData(X=X, obs=obs, var=var)


def _make_gate1_flags():
    return pd.DataFrame(
        {"flag_would_delete": [False, True, False], "n_flags": [0, 1, 0]},
        index=["AAA-S1", "BBB-S1", "CCC-S2"],
    )


def _make_decontx_scores():
    return pd.Series(
        [0.1, 0.2, 0.3], index=["AAA-S1", "BBB-S1", "CCC-S2"], name="decontX_contamination"
    )


def test_build_prepared_adata_drops_antibody_capture_and_keeps_all_cells():
    adata = _make_fixture_adata()
    prepared, report = PREP.build_prepared_adata(adata, _make_gate1_flags(), _make_decontx_scores())

    assert report["n_vars_in"] == 3
    assert report["n_vars_out"] == 2
    assert report["n_vars_dropped_non_gene_expression"] == 1
    assert list(prepared.var_names) == ["GENE1", "GENE2"]
    assert report["n_cells_in"] == report["n_cells_out"] == 3
    assert list(prepared.obs_names) == list(adata.obs_names)  # row order preserved


def test_build_prepared_adata_attaches_qc_annotations_without_filtering():
    adata = _make_fixture_adata()
    prepared, report = PREP.build_prepared_adata(adata, _make_gate1_flags(), _make_decontx_scores())

    assert list(prepared.obs["gate1_flag_would_delete"]) == [False, True, False]
    assert list(prepared.obs["gate1_n_flags"]) == [0, 1, 0]
    assert list(prepared.obs["decontx_estimated_contamination"]) == [0.1, 0.2, 0.3]
    # The flagged cell (BBB-S1) must still be PRESENT -- annotation only, never filtered.
    assert "steph_BBB-S1" in prepared.obs_names
    assert report["n_cells_missing_gate1_flag"] == 0
    assert report["n_cells_missing_decontx_score"] == 0


def test_build_prepared_adata_never_mutates_the_input_in_place():
    adata = _make_fixture_adata()
    n_vars_before = adata.n_vars
    PREP.build_prepared_adata(adata, _make_gate1_flags(), _make_decontx_scores())
    assert adata.n_vars == n_vars_before
    assert "gate1_flag_would_delete" not in adata.obs.columns


def test_build_prepared_adata_reports_missing_annotation_as_na_not_zero():
    adata = _make_fixture_adata()
    flags = _make_gate1_flags().drop(index="CCC-S2")  # CCC-S2 missing on purpose
    scores = _make_decontx_scores().drop(index="CCC-S2")

    prepared, report = PREP.build_prepared_adata(adata, flags, scores)
    assert report["n_cells_missing_gate1_flag"] == 1
    assert report["n_cells_missing_decontx_score"] == 1
    ccc_row = prepared.obs.loc["steph_CCC-S2"]
    assert pd.isna(ccc_row["gate1_flag_would_delete"])
    assert pd.isna(ccc_row["decontx_estimated_contamination"])


def test_build_prepared_adata_refuses_without_a_feature_types_column():
    obs = pd.DataFrame(index=["steph_AAA-S1"])
    var = pd.DataFrame(index=["GENE1"])  # no feature_types column
    adata = ad.AnnData(X=np.array([[1.0]], dtype="float32"), obs=obs, var=var)
    try:
        PREP.build_prepared_adata(adata, _make_gate1_flags(), _make_decontx_scores())
        assert False, "expected SystemExit without a feature_types column"
    except SystemExit:
        pass


def test_validate_prepared_adata_passes_for_a_well_formed_result():
    adata = _make_fixture_adata()
    prepared, report = PREP.build_prepared_adata(adata, _make_gate1_flags(), _make_decontx_scores())
    assert PREP.validate_prepared_adata(adata, prepared, report) is True


def test_mutation_kill_validate_catches_a_cell_count_change():
    """Per DEVELOPMENT.md Sec.2 and review:20261016_141242 R4: the previous
    version of this test wrapped BOTH the target call and its own `assert False` in
    the same `try`, with `except AssertionError: pass` catching either one -- so it
    passed whether or not the validator actually raised. `pytest.raises` wraps ONLY
    the target call, so a validator that stops raising makes this test fail, not pass.
    """
    import pytest

    adata = _make_fixture_adata()
    prepared, report = PREP.build_prepared_adata(adata, _make_gate1_flags(), _make_decontx_scores())

    mutated = prepared[:-1, :].copy()  # drop a cell, simulating an accidental cell-axis filter
    with pytest.raises(AssertionError):
        PREP.validate_prepared_adata(adata, mutated, report)


def test_mutation_kill_validate_catches_a_missing_qc_annotation():
    import pytest

    adata = _make_fixture_adata()
    flags = _make_gate1_flags().drop(index="CCC-S2")
    scores = _make_decontx_scores().drop(index="CCC-S2")
    prepared, report = PREP.build_prepared_adata(adata, flags, scores)

    with pytest.raises(AssertionError):
        PREP.validate_prepared_adata(adata, prepared, report)


def test_mutation_kill_a_no_op_validator_is_caught_by_these_tests(monkeypatch):
    """Per review:20261016_152556 ('两个小问题' #2): the prior version
    of this test only confirmed the no-op itself returns True without raising -- it
    never actually confirmed the OTHER mutation-kill tests in this file would FAIL
    under that mutation, so its name overclaimed what it proved. This version calls
    those tests directly, for real, with `validate_prepared_adata` monkeypatched to
    the no-op, and asserts each one actually fails (via `pytest.raises(AssertionError)`
    itself raising `Failed` when the expected exception never occurs) -- genuinely
    observing the target tests fail, not just inferring it from the no-op's behavior.
    """

    def noop(*a, **k):
        return True

    monkeypatch.setattr(PREP, "validate_prepared_adata", noop)

    with pytest.raises(BaseException):
        test_mutation_kill_validate_catches_a_cell_count_change()
    with pytest.raises(BaseException):
        test_mutation_kill_validate_catches_a_missing_qc_annotation()


def test_strip_prefix_refuses_an_unexpected_prefix():
    try:
        PREP._strip_prefix(["no_prefix_here"])
        assert False, "expected SystemExit on an unexpected prefix"
    except SystemExit:
        pass


def test_validate_rejects_reversed_cell_order_even_with_a_stale_report(tmp_path):
    """Reproduces the exact bypass an independent probe found in the prior round
    (review:20261016_141242 R4): reversing the prepared object's
    cell order while keeping the OLD (now-stale) report that still claims
    `obs_row_order_preserved=True`. The fixed validator must catch this by
    re-deriving order directly from the objects, never trusting that report field.
    """
    adata = _make_fixture_adata()
    prepared, report = PREP.build_prepared_adata(adata, _make_gate1_flags(), _make_decontx_scores())
    assert report["obs_row_order_preserved"] is True  # the stale, now-wrong claim

    reversed_prepared = prepared[::-1, :].copy()
    with pytest.raises(AssertionError):
        PREP.validate_prepared_adata(adata, reversed_prepared, report)


def test_validate_rejects_a_corrupted_count_value(tmp_path):
    """Reproduces the probe's 'change a count to 999 or 0.25' bypass."""
    adata = _make_fixture_adata()
    prepared, report = PREP.build_prepared_adata(adata, _make_gate1_flags(), _make_decontx_scores())

    corrupted = prepared.copy()
    corrupted.X[0, 0] = 999.0
    with pytest.raises(AssertionError):
        PREP.validate_prepared_adata(adata, corrupted, report)


def test_validate_rejects_an_empty_gex_object(tmp_path):
    """Reproduces the probe's 'zero GEX features' bypass -- an object with no
    Gene Expression features at all must never be accepted as prepared.
    """
    obs = pd.DataFrame(index=["steph_AAA-S1"])
    var_all_antibody = pd.DataFrame({"feature_types": ["Antibody Capture"]}, index=["CD3_protein"])
    original = ad.AnnData(X=np.array([[5.0]], dtype="float32"), obs=obs, var=var_all_antibody)

    empty_prepared = ad.AnnData(
        X=np.zeros((1, 0), dtype="float32"),
        obs=original.obs.copy(),
        var=pd.DataFrame(index=[]),
    )
    empty_prepared.obs["gate1_flag_would_delete"] = [False]
    empty_prepared.obs["gate1_n_flags"] = [0]
    empty_prepared.obs["decontx_estimated_contamination"] = [0.1]
    fake_report = {"obs_row_order_preserved": True}

    with pytest.raises(AssertionError):
        PREP.validate_prepared_adata(original, empty_prepared, fake_report)


def test_validate_rejects_a_non_finite_decontx_score():
    """Reproduces the probe's 'DecontX annotation contains Inf' bypass."""
    adata = _make_fixture_adata()
    prepared, report = PREP.build_prepared_adata(adata, _make_gate1_flags(), _make_decontx_scores())

    broken = prepared.copy()
    broken.obs["decontx_estimated_contamination"] = [0.1, float("inf"), 0.3]
    with pytest.raises(AssertionError):
        PREP.validate_prepared_adata(adata, broken, report)


def test_validate_rejects_a_corrupted_count_value_even_under_python_dash_O():
    """Per review:20261016_152556 ('两个小问题' #1): bare `assert`
    statements are stripped under `python -O`/`PYTHONOPTIMIZE` -- this runs the REAL
    validator in a separate `-O` interpreter (a process-level flag, not something a
    monkeypatch within this process can simulate) to confirm the `_require`-based
    rewrite still refuses a corrupted count under optimization, where the old bare
    `assert` would have silently passed.
    """
    import subprocess

    script = f"""
import sys
sys.path.insert(0, {RESULTS_DIR!r})
import numpy as np, pandas as pd
import anndata as ad
import stephenson_monocyte_scprint2_pilot_prep as PREP

obs = pd.DataFrame(index=["steph_AAA-S1", "steph_BBB-S1", "steph_CCC-S2"])
var = pd.DataFrame({{"feature_types": ["Gene Expression", "Gene Expression", "Antibody Capture"]}},
                    index=["GENE1", "GENE2", "CD3_protein"])
X = np.array([[1, 2, 5], [3, 0, 7], [0, 4, 9]], dtype="float32")
adata = ad.AnnData(X=X, obs=obs, var=var)
flags = pd.DataFrame({{"flag_would_delete": [False, True, False], "n_flags": [0, 1, 0]}},
                      index=["AAA-S1", "BBB-S1", "CCC-S2"])
scores = pd.Series([0.1, 0.2, 0.3], index=["AAA-S1", "BBB-S1", "CCC-S2"])

prepared, report = PREP.build_prepared_adata(adata, flags, scores)
prepared.X[0, 0] = 999.0
try:
    PREP.validate_prepared_adata(adata, prepared, report)
except AssertionError:
    sys.exit(42)
sys.exit(0)
"""
    result = subprocess.run([sys.executable, "-O", "-c", script], capture_output=True, text=True)
    assert result.returncode == 42, (
        "expected validate_prepared_adata to still raise AssertionError under -O; "
        f"got exit {result.returncode}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )


def _write_fixture_h5ad_file(path):
    _make_fixture_adata().write_h5ad(path)


def _write_gate1_csv(path):
    import csv

    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["cell_id", "channel_id", "n_flags", "flag_would_delete"])
        w.writeheader()
        w.writerows(
            [
                {
                    "cell_id": "AAA-S1",
                    "channel_id": "s",
                    "n_flags": 0,
                    "flag_would_delete": "False",
                },
                {"cell_id": "BBB-S1", "channel_id": "s", "n_flags": 1, "flag_would_delete": "True"},
                {
                    "cell_id": "CCC-S2",
                    "channel_id": "s",
                    "n_flags": 0,
                    "flag_would_delete": "False",
                },
            ]
        )


def _write_decontx_csv(path):
    import csv

    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["cell_id", "channel_id", "decontX_contamination"])
        w.writeheader()
        w.writerows(
            [
                {"cell_id": "AAA-S1", "channel_id": "s", "decontX_contamination": 0.1},
                {"cell_id": "BBB-S1", "channel_id": "s", "decontX_contamination": 0.2},
                {"cell_id": "CCC-S2", "channel_id": "s", "decontX_contamination": 0.3},
            ]
        )


def test_main_refuses_to_overwrite_an_existing_h5ad_output(tmp_path):
    h5ad_path = tmp_path / "raw.h5ad"
    _write_fixture_h5ad_file(h5ad_path)
    gate1_csv = tmp_path / "gate1.csv"
    _write_gate1_csv(gate1_csv)
    decontx_csv = tmp_path / "decontx.csv"
    _write_decontx_csv(decontx_csv)

    out_h5ad = tmp_path / "out.h5ad"
    out_h5ad.write_text("DO NOT TOUCH")

    try:
        PREP.main(str(h5ad_path), str(gate1_csv), str(decontx_csv), out_h5ad=str(out_h5ad))
        assert False, "expected SystemExit when out_h5ad already exists"
    except SystemExit:
        pass
    assert out_h5ad.read_text() == "DO NOT TOUCH"


def test_main_refuses_to_overwrite_an_existing_summary_even_if_h5ad_is_absent(tmp_path):
    """Reproduces the exact gap an independent fixture found (
    review:20261016_141242 R5): the prior guard only checked the h5ad path, so a
    pre-existing `_summary.json` with no h5ad yet was silently overwritten.
    """
    h5ad_path = tmp_path / "raw.h5ad"
    _write_fixture_h5ad_file(h5ad_path)
    gate1_csv = tmp_path / "gate1.csv"
    _write_gate1_csv(gate1_csv)
    decontx_csv = tmp_path / "decontx.csv"
    _write_decontx_csv(decontx_csv)

    out_h5ad = tmp_path / "out.h5ad"
    out_summary = tmp_path / "out_summary.json"
    out_summary.write_text("DO NOT TOUCH")
    assert not out_h5ad.exists()

    try:
        PREP.main(str(h5ad_path), str(gate1_csv), str(decontx_csv), out_h5ad=str(out_h5ad))
        assert False, "expected SystemExit when out_summary already exists, even without out_h5ad"
    except SystemExit:
        pass
    assert out_summary.read_text() == "DO NOT TOUCH"
    assert not out_h5ad.exists()


def test_mutation_kill_the_overwrite_guard_covers_both_paths(monkeypatch, tmp_path):
    h5ad_path = tmp_path / "raw.h5ad"
    _write_fixture_h5ad_file(h5ad_path)
    gate1_csv = tmp_path / "gate1.csv"
    _write_gate1_csv(gate1_csv)
    decontx_csv = tmp_path / "decontx.csv"
    _write_decontx_csv(decontx_csv)

    out_h5ad = tmp_path / "out.h5ad"
    out_summary = tmp_path / "out_summary.json"
    out_summary.write_text("DO NOT TOUCH")

    monkeypatch.setattr(PREP, "_require_output_does_not_already_exist", lambda *a, **k: None)
    PREP.main(str(h5ad_path), str(gate1_csv), str(decontx_csv), out_h5ad=str(out_h5ad))
    assert out_summary.read_text() != "DO NOT TOUCH"


def test_prep_counts_and_identity_survive_an_actual_h5ad_write_and_read_back(tmp_path):
    """Per review:20261016_141242 R4: exercise the REAL producer/
    consumer boundary (an actual h5ad write + `anndata.read_h5ad` round trip), not
    just in-memory objects -- a serializer bug (e.g. losing obs_names order, or
    silently casting a NaN annotation to 0) would not be caught by an in-memory-only
    test.
    """
    adata = _make_fixture_adata()
    prepared, report = PREP.build_prepared_adata(adata, _make_gate1_flags(), _make_decontx_scores())

    out_path = tmp_path / "prepared.h5ad"
    prepared.write_h5ad(out_path)
    reread = ad.read_h5ad(out_path)

    assert list(reread.obs_names) == list(prepared.obs_names) == list(adata.obs_names)
    assert list(reread.var_names) == list(prepared.var_names)
    assert np.array_equal(np.asarray(reread.X), np.asarray(prepared.X))
    assert list(reread.obs["gate1_flag_would_delete"]) == list(
        prepared.obs["gate1_flag_would_delete"]
    )
    assert PREP.validate_prepared_adata(adata, reread, report) is True
