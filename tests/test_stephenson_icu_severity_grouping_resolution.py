import json
import os
import sys

sys.path.insert(
    0,
    os.path.join(
        os.path.dirname(__file__), "..", "scripts", "34_scfm_discovery_framework", "results"
    ),
)

import stephenson_icu_severity_grouping_resolve_mixed_visits as RESOLVE  # noqa: E402


def test_accepts_the_frozen_approved_input_and_matches_known_counts(tmp_path):
    out_tsv = str(tmp_path / "resolved.tsv")
    out_summary = str(tmp_path / "resolved_summary.json")
    payload = RESOLVE.resolve(out_tsv=out_tsv, out_summary=out_summary)
    counts = payload["donor_level_group_counts_after_resolution_one_group_per_donor"]
    assert counts == {
        "comparison_severe_critical_icu_at_sampling": 25,
        "reference_mild_moderate": 52,
        "excluded_other_severity": 9,
    }
    assert os.path.exists(out_tsv)
    assert os.path.exists(out_summary)


def test_rejects_a_tsv_whose_hash_does_not_match_the_frozen_value(tmp_path, monkeypatch):
    """Reproduces review:20261015_150231 F2's probe: an UNRELATED
    donor's candidate_main_group changed, with the summary's own
    mixed-donor field left stale/unchanged. The frozen-hash check must
    refuse this outright, regardless of what the summary claims.
    """
    with open(RESOLVE.IN_TSV, newline="") as f:
        original_tsv_text = f.read()
    # Flip one character in an UNRELATED donor's row (not sdrf::C-8821's own
    # rows) -- e.g. the first row's candidate_main_group.
    tampered_tsv_text = original_tsv_text.replace(
        "reference_mild_moderate", "comparison_severe_critical_icu_at_sampling", 1
    )
    assert tampered_tsv_text != original_tsv_text
    tampered_tsv_path = tmp_path / "tampered.tsv"
    tampered_tsv_path.write_text(tampered_tsv_text)

    monkeypatch.setattr(RESOLVE, "IN_TSV", str(tampered_tsv_path))
    # The summary is left COMPLETELY UNCHANGED (still claims only
    # sdrf::C-8821 is mixed) -- exactly the probe's scenario.
    out_tsv = str(tmp_path / "resolved.tsv")
    out_summary = str(tmp_path / "resolved_summary.json")
    try:
        RESOLVE.resolve(out_tsv=out_tsv, out_summary=out_summary)
        assert False, "expected SystemExit for a tampered input TSV"
    except SystemExit as e:
        assert "does not match" in str(e) or "FROZEN" in str(e)
    assert not os.path.exists(out_tsv), "must not write any output on a rejected input"
    assert not os.path.exists(out_summary)


def test_rejects_a_summary_whose_hash_does_not_match_the_frozen_value(tmp_path, monkeypatch):
    with open(RESOLVE.IN_SUMMARY) as f:
        original_summary = json.load(f)
    tampered_summary = dict(original_summary)
    tampered_summary["donor_level_group_counts_any_visit"] = dict(
        tampered_summary["donor_level_group_counts_any_visit"]
    )
    tampered_summary["donor_level_group_counts_any_visit"]["reference_mild_moderate"] = 999
    tampered_summary_path = tmp_path / "tampered_summary.json"
    tampered_summary_path.write_text(json.dumps(tampered_summary))

    monkeypatch.setattr(RESOLVE, "IN_SUMMARY", str(tampered_summary_path))
    out_tsv = str(tmp_path / "resolved.tsv")
    out_summary = str(tmp_path / "resolved_summary.json")
    try:
        RESOLVE.resolve(out_tsv=out_tsv, out_summary=out_summary)
        assert False, "expected SystemExit for a tampered input summary"
    except SystemExit:
        pass
    assert not os.path.exists(out_tsv)
    assert not os.path.exists(out_summary)


def test_refuses_to_overwrite_an_existing_output_file(tmp_path):
    """Reproduces F2's sentinel probe: a pre-existing file at the output
    path must survive untouched, never silently overwritten.
    """
    out_tsv = tmp_path / "resolved.tsv"
    out_summary = str(tmp_path / "resolved_summary.json")
    sentinel = "SENTINEL -- must not be overwritten\n"
    out_tsv.write_text(sentinel)

    try:
        RESOLVE.resolve(out_tsv=str(out_tsv), out_summary=out_summary)
        assert False, "expected SystemExit when the output path already exists"
    except SystemExit as e:
        assert "overwrite" in str(e)
    assert out_tsv.read_text() == sentinel, "the existing output file was overwritten"
    assert not os.path.exists(out_summary)


def test_mutation_kill_the_overwrite_guard(tmp_path, monkeypatch):
    """Proves the overwrite guard is load-bearing: removing it lets the
    sentinel from the test above actually get clobbered, reproducing
    exactly the defect F2 found.
    """
    monkeypatch.setattr(RESOLVE, "_require_output_does_not_already_exist", lambda *a, **k: None)
    out_tsv = tmp_path / "resolved.tsv"
    out_summary = str(tmp_path / "resolved_summary.json")
    sentinel = "SENTINEL -- must not be overwritten\n"
    out_tsv.write_text(sentinel)

    RESOLVE.resolve(out_tsv=str(out_tsv), out_summary=out_summary)
    assert out_tsv.read_text() != sentinel, (
        "with the guard mutated away, the sentinel SHOULD be overwritten -- "
        "confirming the guard (not something else) is what protects it normally"
    )


def test_kept_and_excluded_samples_are_verified_against_their_real_rows(tmp_path, monkeypatch):
    """If the hardcoded KEPT/EXCLUDED sample ids somehow pointed at rows that
    are NOT actually in the expected groups, the script must refuse rather
    than silently resolve the wrong sample into the wrong arm. Simulated by
    swapping which constant is checked, without touching the frozen input.
    """
    monkeypatch.setattr(RESOLVE, "KEPT_SAMPLE_ID", RESOLVE.EXCLUDED_SAMPLE_ID)
    with open(RESOLVE.IN_TSV, newline="") as f:
        import csv

        rows = list(csv.DictReader(f, delimiter="\t"))
    try:
        RESOLVE._require_resolved_rows(rows)
        assert (
            False
        ), "expected SystemExit when KEPT_SAMPLE_ID is not actually in the comparison group"
    except SystemExit as e:
        assert "comparison" in str(e).lower()
