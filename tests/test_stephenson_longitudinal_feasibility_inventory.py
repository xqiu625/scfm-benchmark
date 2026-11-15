import os
import sys

RESULTS_DIR = os.path.join(
    os.path.dirname(__file__), "..", "scripts", "34_scfm_discovery_framework", "results"
)
sys.path.insert(0, RESULTS_DIR)

import stephenson_longitudinal_feasibility_inventory as INV  # noqa: E402


def test_transition_label_classifies_stable_worsening_improving_unknown():
    assert INV._transition_label("Mild", "Mild", INV.SEVERITY_ORDER) == "stable"
    assert INV._transition_label("Moderate", "Critical", INV.SEVERITY_ORDER) == "worsening"
    assert INV._transition_label("Critical", "Mild", INV.SEVERITY_ORDER) == "improving"
    assert INV._transition_label("Moderate", "Unknown_Value", INV.SEVERITY_ORDER) == "unknown"


def test_collection_day_sort_key_parses_d_prefixed_values():
    assert INV._collection_day_sort_key("D0") == 0
    assert INV._collection_day_sort_key("D28") == 28
    assert INV._collection_day_sort_key("not_a_day") is None


def test_classify_onset_separates_identical_non_numeric_from_identical_numeric():
    """Per F2 (review:20261015_230801): CV0231's two visits are
    both the literal string 'Healthy' -- identical, but not a day count. This
    must land in its own bucket, never merged with a real identical integer.
    """
    identical_numeric = INV._classify_onset(["9", "9"])
    assert identical_numeric[0] is True  # raw strings identical
    assert identical_numeric[1] is True  # all numeric
    assert identical_numeric[2] == "identical_numeric"

    identical_non_numeric = INV._classify_onset(["Healthy", "Healthy"])
    assert identical_non_numeric[0] is True  # raw strings identical
    assert identical_non_numeric[1] is False  # NOT numeric
    assert identical_non_numeric[2] == "identical_non_numeric"

    differs_numeric = INV._classify_onset(["9", "14"])
    assert differs_numeric[2] == "differs_numeric"

    differs_mixed = INV._classify_onset(["9", "Healthy"])
    assert differs_mixed[2] == "differs_mixed_or_non_numeric"


def test_any_adjacent_transition_catches_a_mid_trajectory_change_that_reverts():
    """A 3-visit donor can read 'stable' first-to-last while a middle visit
    genuinely changed (per F4) -- Mild -> Severe -> Mild. The any-adjacent
    check must catch what the first-to-last check structurally cannot.
    """
    sequence = ["Mild", "Severe", "Mild"]
    assert INV._transition_label(sequence[0], sequence[-1], INV.SEVERITY_ORDER) == "stable"
    assert INV._any_adjacent_transition(sequence, INV.SEVERITY_ORDER) is True

    stable_sequence = ["Mild", "Mild", "Mild"]
    assert INV._any_adjacent_transition(stable_sequence, INV.SEVERITY_ORDER) is False


def test_any_adjacent_transition_reports_unknown_not_false_for_an_unlabelable_step():
    """A sequence containing a value absent from the ordinal scale must report
    'unknown', never be silently read as 'no transition occurred'.
    """
    assert INV._any_adjacent_transition(["Mild", "Unmapped_Value"], INV.SEVERITY_ORDER) == "unknown"


def test_order_donor_visits_refuses_to_fabricate_an_order_from_an_unparseable_day():
    rows = [
        {"visit_resample_flag": "Initial", "collection_day": "D0", "sample_id": "s1"},
        {"visit_resample_flag": "Resample", "collection_day": "not_a_day", "sample_id": "s2"},
    ]
    ordered, order_resolved = INV._order_donor_visits(rows)
    assert order_resolved is False
    # Input row order preserved verbatim when order cannot be resolved -- not reordered.
    assert [r["sample_id"] for r in ordered] == ["s1", "s2"]


def test_order_donor_visits_refuses_to_fabricate_an_order_from_tied_resample_days():
    rows = [
        {"visit_resample_flag": "Initial", "collection_day": "D0", "sample_id": "s1"},
        {"visit_resample_flag": "Resample", "collection_day": "D7", "sample_id": "s2"},
        {"visit_resample_flag": "Resample", "collection_day": "D7", "sample_id": "s3"},
    ]
    _, order_resolved = INV._order_donor_visits(rows)
    assert order_resolved is False


def test_order_donor_visits_orders_three_real_visits_chronologically():
    """Reproduces sdrf::C-8896's own real row order (D13 listed before D7 in
    the roster) -- the sort must still land on D0, D7, D13.
    """
    rows = [
        {"visit_resample_flag": "Initial", "collection_day": "D0", "sample_id": "s1"},
        {"visit_resample_flag": "Resample", "collection_day": "D13", "sample_id": "s2"},
        {"visit_resample_flag": "Resample", "collection_day": "D7", "sample_id": "s3"},
    ]
    ordered, order_resolved = INV._order_donor_visits(rows)
    assert order_resolved is True
    assert [r["sample_id"] for r in ordered] == ["s1", "s3", "s2"]


def test_build_inventory_against_real_roster_matches_known_counts():
    """Calls the pure computation directly against the real, committed roster/
    grouping tables (never a synthetic stand-in, and never writing any file --
    `build_inventory` takes no output path at all).
    """
    donor_rows_out, summary = INV.build_inventory()

    assert summary["n_total_true_donors"] == 86
    assert summary["n_multi_visit_donors"] == 15
    assert summary["n_single_visit_donors"] == 71
    assert summary["visit_count_distribution_among_multi_visit_donors"] == {2: 14, 3: 1}
    assert summary["n_multi_visit_donors_with_order_unresolved"] == 0

    # F2: CV0231's Healthy/Healthy must be its own bucket, not folded into
    # the 14 real identical-integer donors.
    assert summary["days_from_onset_category_counts"] == {
        "identical_numeric": 14,
        "identical_non_numeric": 1,
    }

    worsening_pairs = {
        k: v
        for k, v in summary["severity_transition_counts_first_to_last_visit"].items()
        if "Moderate->Critical" in k
    }
    assert worsening_pairs == {"Moderate->Critical": 1}
    care_worsening = {
        k: v
        for k, v in summary["care_setting_transition_counts_first_to_last_visit"].items()
        if k == "non-ICU->ICU"
    }
    assert care_worsening == {"non-ICU->ICU": 1}

    # sdrf::C-8896: Severe -> Critical -> Critical. First-to-last is a real
    # worsening; the 3-visit trajectory itself must be present, not collapsed.
    c8896 = next(r for r in donor_rows_out if r["true_donor_group"] == "sdrf::C-8896")
    assert c8896["n_visits"] == 3
    assert c8896["order_resolved"] is True
    assert c8896["severity_sequence"] == "Severe,Critical,Critical"
    assert c8896["severity_transition_first_to_last"] == "worsening"
    assert c8896["severity_any_adjacent_transition"] is True

    cv0231 = next(r for r in donor_rows_out if r["true_donor_group"] == "sdrf::CV0231")
    assert cv0231["days_from_onset_category"] == "identical_non_numeric"


def test_build_inventory_does_not_report_a_direction_when_order_is_unresolved(tmp_path):
    """R2 (review:20261016_102834): the regression check reproduced this exact
    scenario (C-8821 with its resample date made unparseable) and got
    `order_resolved=False` together with `severity_transition_first_to_last=
    "worsening"` and the donor counted under `{"Moderate->Critical": 1}` in
    the summary -- a direction computed from arbitrary input row order. Fixed
    behaviour: the transition fields read "unknown" and the donor is counted
    under "order_unresolved", never under any directional label.
    """
    roster_path = tmp_path / "roster.tsv"
    grouping_path = tmp_path / "grouping.tsv"
    roster_path.write_text(
        "sample_id\ttrue_donor_group\tvisit_resample_flag\tcollection_day\tdays_from_onset\t"
        "own_severity_Status_on_day_collection_summary\tsubtype_coverage_checked_this_round\n"
        "s1\tdonor::A\tInitial\tD0\t9\tModerate\tnot re-checked this round\n"
        "s2\tdonor::A\tResample\tnot_a_day\t9\tCritical\tnot re-checked this round\n"
    )
    grouping_path.write_text("sample_id\tcare_setting_at_sampling\ns1\tnon-ICU\ns2\tICU\n")

    donor_rows_out, summary = INV.build_inventory(
        roster_path=str(roster_path), grouping_path=str(grouping_path)
    )

    row = donor_rows_out[0]
    assert row["order_resolved"] is False
    assert row["severity_first"] == "unknown"
    assert row["severity_last"] == "unknown"
    assert row["severity_transition_first_to_last"] == "unknown"
    assert row["severity_any_adjacent_transition"] == "unknown"
    assert row["care_setting_first"] == "unknown"
    assert row["care_setting_transition_first_to_last"] == "unknown"

    assert summary["severity_transition_counts_first_to_last_visit"] == {"order_unresolved": 1}
    assert summary["care_setting_transition_counts_first_to_last_visit"] == {"order_unresolved": 1}
    # The real defect the regression check demonstrated -- a directional label standing in
    # for the unresolved donor -- must not appear anywhere in the summary.
    assert "Moderate->Critical" not in summary["severity_transition_counts_first_to_last_visit"]
    assert "non-ICU->ICU" not in summary["care_setting_transition_counts_first_to_last_visit"]


def test_mutation_kill_the_order_resolved_gate_on_first_to_last_fields(tmp_path):
    """Proves the R2 fix is load-bearing: with the gate removed (mimicking
    the unfixed v2 code the review found), the same fixture as above
    DOES report a direction from arbitrary input-row order.
    """
    roster_path = tmp_path / "roster.tsv"
    grouping_path = tmp_path / "grouping.tsv"
    roster_path.write_text(
        "sample_id\ttrue_donor_group\tvisit_resample_flag\tcollection_day\tdays_from_onset\t"
        "own_severity_Status_on_day_collection_summary\tsubtype_coverage_checked_this_round\n"
        "s1\tdonor::A\tInitial\tD0\t9\tModerate\tnot re-checked this round\n"
        "s2\tdonor::A\tResample\tnot_a_day\t9\tCritical\tnot re-checked this round\n"
    )
    grouping_path.write_text("sample_id\tcare_setting_at_sampling\ns1\tnon-ICU\ns2\tICU\n")

    roster_rows = INV._load_tsv(str(roster_path))
    rows = roster_rows  # single donor, both rows
    ordered, order_resolved = INV._order_donor_visits(rows)
    assert order_resolved is False

    # Reproduce the pre-R2-fix code path directly: compute first-to-last
    # unconditionally, ignoring order_resolved.
    severity_sequence = [r["own_severity_Status_on_day_collection_summary"] for r in ordered]
    severity_first, severity_last = severity_sequence[0], severity_sequence[-1]
    unfixed_label = INV._transition_label(severity_first, severity_last, INV.SEVERITY_ORDER)
    assert unfixed_label == "worsening"  # the exact wrong answer the review reported

    # The guarded production path must not reproduce it.
    _, summary = INV.build_inventory(roster_path=str(roster_path), grouping_path=str(grouping_path))
    assert summary["severity_transition_counts_first_to_last_visit"] != {"Moderate->Critical": 1}


def test_build_inventory_with_no_multi_visit_donors_writes_header_only_tsv(tmp_path):
    """F1/F4: an input with zero repeated visits must not crash on
    `donor_rows_out[0]` and must still write a valid (empty) inventory.
    """
    roster_path = tmp_path / "roster.tsv"
    grouping_path = tmp_path / "grouping.tsv"
    roster_path.write_text(
        "sample_id\ttrue_donor_group\tvisit_resample_flag\tcollection_day\tdays_from_onset\t"
        "own_severity_Status_on_day_collection_summary\tsubtype_coverage_checked_this_round\n"
        "s1\tdonor::A\tInitial\tD0\t9\tMild\tnot re-checked this round\n"
    )
    grouping_path.write_text("sample_id\tcare_setting_at_sampling\ns1\tnon-ICU\n")

    out_tsv = tmp_path / "out.tsv"
    out_summary = tmp_path / "out_summary.json"
    summary = INV.main(
        roster_path=str(roster_path),
        grouping_path=str(grouping_path),
        out_tsv=str(out_tsv),
        out_summary=str(out_summary),
    )

    assert summary["n_multi_visit_donors"] == 0
    assert out_tsv.exists()
    with open(out_tsv) as f:
        assert f.readline().strip() == "\t".join(INV.FIELDNAMES)
        assert f.read() == ""


def test_main_refuses_to_overwrite_an_existing_output_file(tmp_path):
    roster_path = tmp_path / "roster.tsv"
    grouping_path = tmp_path / "grouping.tsv"
    roster_path.write_text(
        "sample_id\ttrue_donor_group\tvisit_resample_flag\tcollection_day\tdays_from_onset\t"
        "own_severity_Status_on_day_collection_summary\tsubtype_coverage_checked_this_round\n"
        "s1\tdonor::A\tInitial\tD0\t9\tMild\tnot re-checked this round\n"
    )
    grouping_path.write_text("sample_id\tcare_setting_at_sampling\ns1\tnon-ICU\n")

    out_tsv = tmp_path / "out.tsv"
    out_summary = tmp_path / "out_summary.json"
    sentinel = "DO NOT TOUCH"
    out_tsv.write_text(sentinel)

    try:
        INV.main(
            roster_path=str(roster_path),
            grouping_path=str(grouping_path),
            out_tsv=str(out_tsv),
            out_summary=str(out_summary),
        )
        assert False, "expected SystemExit when out_tsv already exists"
    except SystemExit:
        pass

    assert out_tsv.read_text() == sentinel
    assert not out_summary.exists()


def test_mutation_kill_the_overwrite_guard(tmp_path, monkeypatch):
    """Per DEVELOPMENT.md Sec.2/6: a guard the caller can switch off is not a
    guard. Confirms the sentinel from the previous test WOULD have been
    clobbered if `_require_output_does_not_already_exist` were a no-op.
    """
    roster_path = tmp_path / "roster.tsv"
    grouping_path = tmp_path / "grouping.tsv"
    roster_path.write_text(
        "sample_id\ttrue_donor_group\tvisit_resample_flag\tcollection_day\tdays_from_onset\t"
        "own_severity_Status_on_day_collection_summary\tsubtype_coverage_checked_this_round\n"
        "s1\tdonor::A\tInitial\tD0\t9\tMild\tnot re-checked this round\n"
    )
    grouping_path.write_text("sample_id\tcare_setting_at_sampling\ns1\tnon-ICU\n")

    out_tsv = tmp_path / "out.tsv"
    out_summary = tmp_path / "out_summary.json"
    out_tsv.write_text("DO NOT TOUCH")

    monkeypatch.setattr(INV, "_require_output_does_not_already_exist", lambda *a, **k: None)
    INV.main(
        roster_path=str(roster_path),
        grouping_path=str(grouping_path),
        out_tsv=str(out_tsv),
        out_summary=str(out_summary),
    )
    assert out_tsv.read_text() != "DO NOT TOUCH"


def test_main_never_touches_the_committed_v1_output_files(tmp_path):
    """F1: this module's default output path moved to v2; writing a fresh
    inventory into a tmp_path target must leave the committed v1 files'
    bytes exactly as they were.
    """
    before = {p: INV._sha256_file(p) for p in (INV.OUT_TSV_V1, INV.OUT_SUMMARY_V1)}

    roster_path = tmp_path / "roster.tsv"
    grouping_path = tmp_path / "grouping.tsv"
    roster_path.write_text(
        "sample_id\ttrue_donor_group\tvisit_resample_flag\tcollection_day\tdays_from_onset\t"
        "own_severity_Status_on_day_collection_summary\tsubtype_coverage_checked_this_round\n"
        "s1\tdonor::A\tInitial\tD0\t9\tMild\tnot re-checked this round\n"
    )
    grouping_path.write_text("sample_id\tcare_setting_at_sampling\ns1\tnon-ICU\n")
    INV.main(
        roster_path=str(roster_path),
        grouping_path=str(grouping_path),
        out_tsv=str(tmp_path / "out.tsv"),
        out_summary=str(tmp_path / "out_summary.json"),
    )

    after = {p: INV._sha256_file(p) for p in (INV.OUT_TSV_V1, INV.OUT_SUMMARY_V1)}
    assert before == after
