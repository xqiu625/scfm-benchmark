import json
import os
import sys

import h5py
import numpy as np

RESULTS_DIR = os.path.join(
    os.path.dirname(__file__), "..", "scripts", "34_scfm_discovery_framework", "results"
)
sys.path.insert(0, RESULTS_DIR)

import stephenson_sample_subtype_cell_counts as CNT  # noqa: E402


def _write_fixture_h5ad(path, sample_codes, sample_cats, subtype_codes, subtype_cats):
    """Hand-builds the exact minimal h5py structure the real anndata object
    uses for a categorical obs column (`categories` + `codes`), rather than
    depending on the anndata package being importable in the test env.
    """
    with h5py.File(path, "w") as f:
        obs = f.create_group("obs")
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


def test_read_sample_subtype_counts_tallies_a_small_fixture(tmp_path):
    h5ad_path = tmp_path / "fixture.h5ad"
    # 2 samples (s0, s1), 2 subtypes (A, B); 5 cells total.
    # s0: 2xA, 1xB ; s1: 1xA, 1xB
    _write_fixture_h5ad(
        h5ad_path,
        sample_codes=[0, 0, 0, 1, 1],
        sample_cats=["s0", "s1"],
        subtype_codes=[0, 0, 1, 0, 1],
        subtype_cats=["A", "B"],
    )

    counts, n_cells_total, sample_cats, subtype_cats = CNT.read_sample_subtype_counts(
        str(h5ad_path)
    )

    assert n_cells_total == 5
    assert counts["s0"] == {"A": 2, "B": 1}
    assert counts["s1"] == {"A": 1, "B": 1}


def test_read_sample_subtype_counts_rejects_missing_minus_one_sample_code(tmp_path):
    """R1 (review:20261016_102834): a `-1` ("missing", per
    anndata's own on-disk categorical format) sample code must never be
    silently resolved to the LAST category via Python negative indexing.
    """
    h5ad_path = tmp_path / "fixture.h5ad"
    _write_fixture_h5ad(
        h5ad_path,
        sample_codes=[0, -1],
        sample_cats=["s0", "s1"],
        subtype_codes=[-1, 0],
        subtype_cats=["A", "B"],
    )
    try:
        CNT.read_sample_subtype_counts(str(h5ad_path))
        assert False, "expected SystemExit on a -1 (missing) code"
    except SystemExit:
        pass


def test_mutation_kill_the_code_validation_guard_lets_minus_one_wrap_to_last_category(tmp_path):
    """Confirms the R1 defect actually reproduces when the guard is removed --
    the exact failure the regression check demonstrated: `sample_codes=[0, -1]`,
    `subtype_codes=[-1, 0]` silently counted as `s0/B=1, s1/A=1` with no error.
    """
    h5ad_path = tmp_path / "fixture.h5ad"
    _write_fixture_h5ad(
        h5ad_path,
        sample_codes=[0, -1],
        sample_cats=["s0", "s1"],
        subtype_codes=[-1, 0],
        subtype_cats=["A", "B"],
    )
    with h5py.File(h5ad_path, "r") as f:
        sample_group = f["obs"]["Sample_ID"]
        subtype_group = f["obs"]["orig_celltype"]
        sample_cats = [
            c.decode() if isinstance(c, bytes) else c for c in sample_group["categories"][:]
        ]
        subtype_cats = [
            c.decode() if isinstance(c, bytes) else c for c in subtype_group["categories"][:]
        ]
        sample_codes = sample_group["codes"][:]
        subtype_codes = subtype_group["codes"][:]

    # Reproduce the unfixed loop body directly (no call into the guarded
    # production function) to show the mutation this test targets is real.
    counts = {sid: {st: 0 for st in subtype_cats} for sid in sample_cats}
    for s_code, t_code in zip(sample_codes, subtype_codes):
        counts[sample_cats[s_code]][subtype_cats[t_code]] += 1
    assert counts["s0"]["B"] == 1  # the -1 subtype code silently wrapped to "B"
    assert counts["s1"]["A"] == 1  # the -1 sample code silently wrapped to "s1"

    # The guarded production path must refuse instead.
    try:
        CNT.read_sample_subtype_counts(str(h5ad_path))
        assert False, "expected SystemExit; the unguarded loop above proves the defect is real"
    except SystemExit:
        pass


def test_read_sample_subtype_counts_rejects_other_negative_and_out_of_range_codes(tmp_path):
    """R1: not just -1 -- any negative or any code >= len(categories) must be
    refused, never indexed (a more-negative code wraps further into the
    category list; an over-range code is an IndexError today but must never
    become a silent wraparound if the category list grows).
    """
    h5ad_path = tmp_path / "fixture.h5ad"
    _write_fixture_h5ad(
        h5ad_path,
        sample_codes=[0, -2],
        sample_cats=["s0", "s1"],
        subtype_codes=[0, 1],
        subtype_cats=["A", "B"],
    )
    try:
        CNT.read_sample_subtype_counts(str(h5ad_path))
        assert False, "expected SystemExit on a -2 (out-of-range negative) code"
    except SystemExit:
        pass

    h5ad_path2 = tmp_path / "fixture2.h5ad"
    _write_fixture_h5ad(
        h5ad_path2,
        sample_codes=[0, 1],
        sample_cats=["s0", "s1"],
        subtype_codes=[0, 5],  # only 2 subtypes exist (indices 0-1)
        subtype_cats=["A", "B"],
    )
    try:
        CNT.read_sample_subtype_counts(str(h5ad_path2))
        assert False, "expected SystemExit on an out-of-range positive code"
    except SystemExit:
        pass


def test_read_sample_subtype_counts_rejects_mismatched_code_array_lengths(tmp_path):
    h5ad_path = tmp_path / "fixture.h5ad"
    _write_fixture_h5ad(
        h5ad_path,
        sample_codes=[0, 0, 1],
        sample_cats=["s0", "s1"],
        subtype_codes=[0, 1],  # deliberately shorter
        subtype_cats=["A", "B"],
    )
    try:
        CNT.read_sample_subtype_counts(str(h5ad_path))
        assert False, "expected SystemExit on mismatched code-array lengths"
    except SystemExit:
        pass


def test_main_writes_expected_tsv_and_summary(tmp_path):
    h5ad_path = tmp_path / "fixture.h5ad"
    _write_fixture_h5ad(
        h5ad_path,
        sample_codes=[0, 0, 0, 1, 1],
        sample_cats=["s0", "s1"],
        subtype_codes=[0, 0, 1, 0, 1],
        subtype_cats=["A", "B"],
    )
    out_tsv = tmp_path / "out.tsv"
    out_summary = tmp_path / "out_summary.json"

    summary = CNT.main(h5ad_path=str(h5ad_path), out_tsv=str(out_tsv), out_summary=str(out_summary))

    assert summary["n_cells_total_this_object"] == 5
    assert summary["n_samples"] == 2
    assert sorted(summary["subtypes"]) == ["A", "B"]

    with open(out_tsv) as f:
        rows = list(f)
    assert rows[0].strip() == "sample_id\tn_cells_this_sample_total\tA\tB"
    body = {r.split("\t")[0]: r.strip().split("\t") for r in rows[1:]}
    assert body["s0"] == ["s0", "3", "2", "1"]
    assert body["s1"] == ["s1", "2", "1", "1"]

    with open(out_summary) as f:
        assert json.load(f)["n_samples"] == 2


def test_main_refuses_to_overwrite_an_existing_output_file(tmp_path):
    h5ad_path = tmp_path / "fixture.h5ad"
    _write_fixture_h5ad(
        h5ad_path,
        sample_codes=[0],
        sample_cats=["s0"],
        subtype_codes=[0],
        subtype_cats=["A"],
    )
    out_tsv = tmp_path / "out.tsv"
    out_summary = tmp_path / "out_summary.json"
    out_tsv.write_text("DO NOT TOUCH")

    try:
        CNT.main(h5ad_path=str(h5ad_path), out_tsv=str(out_tsv), out_summary=str(out_summary))
        assert False, "expected SystemExit when out_tsv already exists"
    except SystemExit:
        pass

    assert out_tsv.read_text() == "DO NOT TOUCH"
    assert not out_summary.exists()


def test_mutation_kill_the_overwrite_guard(tmp_path, monkeypatch):
    h5ad_path = tmp_path / "fixture.h5ad"
    _write_fixture_h5ad(
        h5ad_path,
        sample_codes=[0],
        sample_cats=["s0"],
        subtype_codes=[0],
        subtype_cats=["A"],
    )
    out_tsv = tmp_path / "out.tsv"
    out_summary = tmp_path / "out_summary.json"
    out_tsv.write_text("DO NOT TOUCH")

    monkeypatch.setattr(CNT, "_require_output_does_not_already_exist", lambda *a, **k: None)
    CNT.main(h5ad_path=str(h5ad_path), out_tsv=str(out_tsv), out_summary=str(out_summary))
    assert out_tsv.read_text() != "DO NOT TOUCH"
