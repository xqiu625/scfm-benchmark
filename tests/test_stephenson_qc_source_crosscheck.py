import os
import sys

import h5py
import numpy as np

RESULTS_DIR = os.path.join(
    os.path.dirname(__file__), "..", "scripts", "34_scfm_discovery_framework", "results"
)
sys.path.insert(0, RESULTS_DIR)

import stephenson_qc_source_crosscheck as XCHK  # noqa: E402


def _write_fixture_h5ad(path, obs_names, sample_codes, sample_cats, subtype_codes, subtype_cats):
    with h5py.File(path, "w") as f:
        obs = f.create_group("obs")
        obs.create_dataset(
            "_index", data=np.array(obs_names, dtype=object), dtype=h5py.string_dtype()
        )
        sg = obs.create_group("Sample_ID")
        sg.create_dataset(
            "categories", data=np.array(sample_cats, dtype=object), dtype=h5py.string_dtype()
        )
        sg.create_dataset("codes", data=np.array(sample_codes, dtype="int8"))
        tg = obs.create_group("orig_celltype")
        tg.create_dataset(
            "categories", data=np.array(subtype_cats, dtype=object), dtype=h5py.string_dtype()
        )
        tg.create_dataset("codes", data=np.array(subtype_codes, dtype="int8"))


def _write_gate1_csv(path, rows):
    import csv

    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["cell_id", "channel_id", "n_flags", "flag_would_delete"])
        w.writeheader()
        w.writerows(rows)


def _write_decontx_csv(path, rows):
    import csv

    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["cell_id", "channel_id", "decontX_contamination"])
        w.writeheader()
        w.writerows(rows)


def test_read_monocyte_obs_strips_prefix_and_reads_three_columns(tmp_path):
    h5ad_path = tmp_path / "fixture.h5ad"
    _write_fixture_h5ad(
        h5ad_path,
        obs_names=["steph_AAA-S1", "steph_BBB-S2"],
        sample_codes=[0, 1],
        sample_cats=["S1", "S2"],
        subtype_codes=[0, 1],
        subtype_cats=["CD14_mono", "CD16_mono"],
    )
    idx, stripped, sample_of_cell, subtype_of_cell = XCHK.read_monocyte_obs(str(h5ad_path))
    assert list(stripped) == ["AAA-S1", "BBB-S2"]
    assert list(sample_of_cell) == ["S1", "S2"]
    assert list(subtype_of_cell) == ["CD14_mono", "CD16_mono"]


def test_read_monocyte_obs_refuses_an_obs_name_without_the_expected_prefix(tmp_path):
    h5ad_path = tmp_path / "fixture.h5ad"
    _write_fixture_h5ad(
        h5ad_path,
        obs_names=["AAA-S1"],  # missing "steph_" prefix
        sample_codes=[0],
        sample_cats=["S1"],
        subtype_codes=[0],
        subtype_cats=["CD14_mono"],
    )
    try:
        XCHK.read_monocyte_obs(str(h5ad_path))
        assert False, "expected SystemExit on a missing prefix"
    except SystemExit:
        pass


def test_crosscheck_reports_zero_mismatches_for_a_fully_joinable_fixture(tmp_path):
    h5ad_path = tmp_path / "fixture.h5ad"
    _write_fixture_h5ad(
        h5ad_path,
        obs_names=["steph_AAA-S1", "steph_BBB-S1", "steph_CCC-S2"],
        sample_codes=[0, 0, 1],
        sample_cats=["S1", "S2"],
        subtype_codes=[0, 0, 1],
        subtype_cats=["CD14_mono", "CD16_mono"],
    )
    gate1_csv = tmp_path / "gate1.csv"
    _write_gate1_csv(
        gate1_csv,
        [
            {"cell_id": "AAA-S1", "channel_id": "S1", "n_flags": 0, "flag_would_delete": "False"},
            {"cell_id": "BBB-S1", "channel_id": "S1", "n_flags": 1, "flag_would_delete": "True"},
            {"cell_id": "CCC-S2", "channel_id": "S2", "n_flags": 0, "flag_would_delete": "False"},
        ],
    )
    decontx_csv = tmp_path / "decontx.csv"
    _write_decontx_csv(
        decontx_csv,
        [
            {"cell_id": "AAA-S1", "channel_id": "S1", "decontX_contamination": 0.1},
            {"cell_id": "BBB-S1", "channel_id": "S1", "decontX_contamination": 0.2},
            {"cell_id": "CCC-S2", "channel_id": "S2", "decontX_contamination": 0.3},
        ],
    )

    report = XCHK.crosscheck(str(h5ad_path), str(gate1_csv), str(decontx_csv))
    assert report["n_cells_total"] == 3
    assert report["n_gate1_flag_missing"] == 0
    assert report["n_gate1_channel_mismatch"] == 0
    assert report["n_decontx_score_missing"] == 0
    assert report["n_decontx_channel_mismatch"] == 0
    assert report["per_subtype"]["CD14_mono"]["n_cells"] == 2
    assert report["per_subtype"]["CD14_mono"]["n_gate1_flag_would_delete_true"] == 1
    assert (
        abs(report["per_subtype"]["CD14_mono"]["decontx_mean_estimated_contamination"] - 0.15)
        < 1e-9
    )
    assert report["per_subtype"]["CD16_mono"]["n_cells"] == 1


def test_crosscheck_reports_missing_and_channel_mismatch_not_silently(tmp_path):
    h5ad_path = tmp_path / "fixture.h5ad"
    _write_fixture_h5ad(
        h5ad_path,
        obs_names=["steph_AAA-S1", "steph_ZZZ-S1"],
        sample_codes=[0, 0],
        sample_cats=["S1"],
        subtype_codes=[0, 0],
        subtype_cats=["CD14_mono"],
    )
    gate1_csv = tmp_path / "gate1.csv"
    _write_gate1_csv(
        gate1_csv,
        [
            {
                "cell_id": "AAA-S1",
                "channel_id": "WRONG_CHANNEL",
                "n_flags": 0,
                "flag_would_delete": "False",
            },
            # ZZZ-S1 deliberately absent
        ],
    )
    decontx_csv = tmp_path / "decontx.csv"
    _write_decontx_csv(
        decontx_csv,
        [
            {"cell_id": "AAA-S1", "channel_id": "S1", "decontX_contamination": 0.1},
        ],
    )

    report = XCHK.crosscheck(str(h5ad_path), str(gate1_csv), str(decontx_csv))
    assert report["n_gate1_flag_missing"] == 1  # ZZZ-S1
    assert report["n_gate1_channel_mismatch"] == 1  # AAA-S1's channel disagrees
    assert report["n_decontx_score_missing"] == 1  # ZZZ-S1


def test_main_refuses_to_overwrite_an_existing_output_file(tmp_path):
    h5ad_path = tmp_path / "fixture.h5ad"
    _write_fixture_h5ad(
        h5ad_path,
        obs_names=["steph_AAA-S1"],
        sample_codes=[0],
        sample_cats=["S1"],
        subtype_codes=[0],
        subtype_cats=["CD14_mono"],
    )
    gate1_csv = tmp_path / "gate1.csv"
    _write_gate1_csv(
        gate1_csv,
        [{"cell_id": "AAA-S1", "channel_id": "S1", "n_flags": 0, "flag_would_delete": "False"}],
    )
    decontx_csv = tmp_path / "decontx.csv"
    _write_decontx_csv(
        decontx_csv, [{"cell_id": "AAA-S1", "channel_id": "S1", "decontX_contamination": 0.1}]
    )

    out_json = tmp_path / "out.json"
    out_json.write_text("DO NOT TOUCH")

    try:
        XCHK.main(str(h5ad_path), str(gate1_csv), str(decontx_csv), out_json=str(out_json))
        assert False, "expected SystemExit when out_json already exists"
    except SystemExit:
        pass
    assert out_json.read_text() == "DO NOT TOUCH"


def test_mutation_kill_the_overwrite_guard(tmp_path, monkeypatch):
    h5ad_path = tmp_path / "fixture.h5ad"
    _write_fixture_h5ad(
        h5ad_path,
        obs_names=["steph_AAA-S1"],
        sample_codes=[0],
        sample_cats=["S1"],
        subtype_codes=[0],
        subtype_cats=["CD14_mono"],
    )
    gate1_csv = tmp_path / "gate1.csv"
    _write_gate1_csv(
        gate1_csv,
        [{"cell_id": "AAA-S1", "channel_id": "S1", "n_flags": 0, "flag_would_delete": "False"}],
    )
    decontx_csv = tmp_path / "decontx.csv"
    _write_decontx_csv(
        decontx_csv, [{"cell_id": "AAA-S1", "channel_id": "S1", "decontX_contamination": 0.1}]
    )

    out_json = tmp_path / "out.json"
    out_json.write_text("DO NOT TOUCH")

    monkeypatch.setattr(XCHK, "_require_output_does_not_already_exist", lambda *a, **k: None)
    XCHK.main(str(h5ad_path), str(gate1_csv), str(decontx_csv), out_json=str(out_json))
    assert out_json.read_text() != "DO NOT TOUCH"


def test_crosscheck_against_real_committed_source_matches_known_numbers():
    """Calls crosscheck() directly against the REAL HPCC paths this module defaults to --
    only runs meaningfully on HPCC where those paths exist; skipped elsewhere.
    """
    if not os.path.exists(XCHK.RAW_H5AD):
        import pytest

        pytest.skip("real HPCC data paths not present in this environment")

    report = XCHK.crosscheck()
    assert report["n_cells_total"] == 128503
    assert report["n_gate1_flag_missing"] == 0
    assert report["n_gate1_channel_mismatch"] == 0
    assert report["n_decontx_score_missing"] == 0
    assert report["n_decontx_channel_mismatch"] == 0
    assert set(report["per_subtype"].keys()) == {
        "C1_CD16_mono",
        "CD14_mono",
        "CD16_mono",
        "CD83_CD14_mono",
        "Mono_prolif",
    }


def test_read_monocyte_obs_rejects_a_minus_one_sample_code(tmp_path):
    """R6 (review:20261016_141242): the same -1-wraps-to-last-category
    defect already fixed elsewhere, reintroduced in this NEW script."""
    h5ad_path = tmp_path / "fixture.h5ad"
    _write_fixture_h5ad(
        h5ad_path,
        obs_names=["steph_AAA-S1", "steph_BBB-S1"],
        sample_codes=[0, -1],
        sample_cats=["S1"],
        subtype_codes=[0, 0],
        subtype_cats=["CD14_mono"],
    )
    try:
        XCHK.read_monocyte_obs(str(h5ad_path))
        assert False, "expected SystemExit on a -1 (missing) Sample_ID code"
    except SystemExit:
        pass


def test_mutation_kill_code_validation_lets_minus_one_wrap_to_last_category(tmp_path):
    h5ad_path = tmp_path / "fixture.h5ad"
    _write_fixture_h5ad(
        h5ad_path,
        obs_names=["steph_AAA-S1", "steph_BBB-S1"],
        sample_codes=[0, -1],
        sample_cats=["S1"],
        subtype_codes=[0, 0],
        subtype_cats=["CD14_mono"],
    )
    with h5py.File(h5ad_path, "r") as f:
        sample_cats = [
            c.decode() if isinstance(c, bytes) else c
            for c in f["obs"]["Sample_ID"]["categories"][:]
        ]
        sample_codes = f["obs"]["Sample_ID"]["codes"][:]
    # Reproduce the unguarded indexing directly -- proves the defect is real.
    unguarded = np.array(sample_cats)[sample_codes]
    assert unguarded[1] == "S1"  # -1 silently wrapped to the (only, "last") category

    import pytest

    with pytest.raises(SystemExit):
        XCHK.read_monocyte_obs(str(h5ad_path))


def test_load_csv_rows_unique_by_id_refuses_duplicate_cell_ids(tmp_path):
    path = tmp_path / "dup.csv"
    _write_gate1_csv(
        path,
        [
            {"cell_id": "AAA-S1", "channel_id": "S1", "n_flags": 0, "flag_would_delete": "False"},
            {"cell_id": "AAA-S1", "channel_id": "S1", "n_flags": 1, "flag_would_delete": "True"},
        ],
    )
    try:
        XCHK._load_csv_rows_unique_by_id(
            str(path), ["cell_id", "channel_id", "n_flags", "flag_would_delete"]
        )
        assert False, "expected SystemExit on a duplicate cell_id"
    except SystemExit:
        pass


def test_parse_strict_bool_refuses_unrecognized_text():
    assert XCHK._parse_strict_bool("True", "x") is True
    assert XCHK._parse_strict_bool("false", "x") is False
    try:
        XCHK._parse_strict_bool("maybe", "x")
        assert False, "expected SystemExit on an unrecognized boolean string"
    except SystemExit:
        pass
