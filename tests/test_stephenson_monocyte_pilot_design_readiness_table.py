import json
import os
import sys

RESULTS_DIR = os.path.join(
    os.path.dirname(__file__), "..", "scripts", "34_scfm_discovery_framework", "results"
)
sys.path.insert(0, RESULTS_DIR)

import stephenson_monocyte_pilot_design_readiness_table as DRT  # noqa: E402


def _write_fixtures(tmp_path):
    roster_path = tmp_path / "roster.tsv"
    grouping_path = tmp_path / "grouping.tsv"
    counts_path = tmp_path / "counts.tsv"
    coverage_summary_path = tmp_path / "coverage_summary.json"

    roster_path.write_text(
        "sample_id\ttrue_donor_group\tvisit_resample_flag\tcollection_day\t"
        "own_severity_Status_on_day_collection_summary\n"
        "s1\tdonor::A\tInitial\tD0\tModerate\n"
        "s2\tdonor::A\tResample\tD28\tCritical\n"
        "s3\tdonor::B\tInitial\tD0\tAsymptomatic\n"
        "s4\tdonor::B\tResample\tD28\tAsymptomatic\n"
        "s5\tdonor::C\tInitial\tD0\tMild\n"
        "s6\tdonor::C\tResample\tD28\tMild\n"
    )
    grouping_path.write_text(
        "sample_id\tcare_setting_at_sampling\tstatus_on_day_collection_raw\n"
        "s1\tnon-ICU\tWard_noO2\n"
        "s2\tICU\tITU_NIV\n"
        "s3\tnon-ICU\tHome\n"
        "s4\tnon-ICU\tHome\n"
        "s5\tnon-ICU\tWard_noO2\n"
        "s6\tnon-ICU\tWard_noO2\n"
    )
    counts_path.write_text(
        "sample_id\tn_cells_this_sample_total\tA\tB\n"
        "s1\t30\t10\t20\n"
        "s2\t150\t60\t90\n"
        "s3\t5\t2\t3\n"
        "s4\t5\t2\t3\n"
        # s5/s6 (donor::C) deliberately missing from the counts table.
    )
    coverage_summary_path.write_text(
        json.dumps(
            {
                "symptomatic_donor_rows": [{"true_donor_group": "donor::A"}],
                "asymptomatic_donor_rows": [{"true_donor_group": "donor::B"}],
            }
        )
    )
    return str(roster_path), str(grouping_path), str(counts_path), str(coverage_summary_path)


def test_build_design_readiness_table_joins_counts_severity_and_icu_evidence(tmp_path):
    roster_path, grouping_path, counts_path, coverage_summary_path = _write_fixtures(tmp_path)

    table_rows, fieldnames, cutoff_sensitivity, summary = DRT.build_design_readiness_table(
        roster_path=roster_path,
        grouping_path=grouping_path,
        counts_path=counts_path,
        coverage_summary_path=coverage_summary_path,
    )

    assert summary["n_multi_visit_donors"] == 3
    assert summary["n_donor_visit_rows"] == 6

    s1 = next(r for r in table_rows if r["sample_id"] == "s1")
    assert s1["symptomatic_or_asymptomatic"] == "symptomatic"
    assert s1["own_severity"] == "Moderate"
    assert s1["care_setting_at_sampling"] == "non-ICU"
    assert s1["status_on_day_collection_raw_icu_evidence"] == "Ward_noO2"
    assert s1["n_A_cells_raw"] == "10"
    assert s1["n_B_cells_raw"] == "20"

    s3 = next(r for r in table_rows if r["sample_id"] == "s3")
    assert s3["symptomatic_or_asymptomatic"] == "asymptomatic"

    # donor::C's visits are missing from the counts table -- must report "unknown",
    # never fold into a 0 count.
    c_rows = [r for r in table_rows if r["true_donor_group"] == "donor::C"]
    assert all(r["n_monocyte_cells_this_sample_total"] == "unknown" for r in c_rows)
    assert all(r["n_A_cells_raw"] == "unknown" for r in c_rows)
    assert summary["n_sample_ids_missing_from_counts_table"] == 2
    assert sorted(summary["sample_ids_missing_from_counts_table"]) == ["s5", "s6"]


def test_cutoff_sensitivity_requires_every_visit_to_clear_the_threshold(tmp_path):
    roster_path, grouping_path, counts_path, coverage_summary_path = _write_fixtures(tmp_path)

    _, _, cutoff_sensitivity, _ = DRT.build_design_readiness_table(
        roster_path=roster_path,
        grouping_path=grouping_path,
        counts_path=counts_path,
        coverage_summary_path=coverage_summary_path,
    )

    # Subtype A: donor::A has [10, 60] -- clears cutoff 1, fails cutoff 20 (visit1=10).
    a_1 = cutoff_sensitivity["A"]["1"]
    a_1_donors = {r["true_donor_group"] for r in a_1["symptomatic"]}
    assert "donor::A" in a_1_donors

    a_20 = cutoff_sensitivity["A"]["20"]
    a_20_donors = {r["true_donor_group"] for r in a_20["symptomatic"]}
    assert "donor::A" not in a_20_donors  # visit1 (10) < 20

    # donor::B has [2, 2] for subtype A -- fails cutoff 1's "both visits >= 1"? No: 2>=1 passes.
    b_1 = cutoff_sensitivity["A"]["1"]
    assert "donor::B" in {r["true_donor_group"] for r in b_1["asymptomatic"]}

    # donor::C is missing from the counts table entirely -- must never appear
    # in any cutoff bucket (not silently treated as passing or failing).
    for cutoff_bucket in cutoff_sensitivity["A"].values():
        all_donors = (
            cutoff_bucket["symptomatic"] + cutoff_bucket["asymptomatic"] + cutoff_bucket["unknown"]
        )
        assert "donor::C" not in {r["true_donor_group"] for r in all_donors}


def test_main_refuses_to_overwrite_an_existing_output_file(tmp_path):
    roster_path, grouping_path, counts_path, coverage_summary_path = _write_fixtures(tmp_path)
    out_tsv = tmp_path / "out.tsv"
    out_json = tmp_path / "out.json"
    out_tsv.write_text("DO NOT TOUCH")

    try:
        DRT.main(
            roster_path=roster_path,
            grouping_path=grouping_path,
            counts_path=counts_path,
            coverage_summary_path=coverage_summary_path,
            out_tsv=str(out_tsv),
            out_json=str(out_json),
        )
        assert False, "expected SystemExit when out_tsv already exists"
    except SystemExit:
        pass

    assert out_tsv.read_text() == "DO NOT TOUCH"
    assert not out_json.exists()


def test_mutation_kill_the_overwrite_guard(tmp_path, monkeypatch):
    roster_path, grouping_path, counts_path, coverage_summary_path = _write_fixtures(tmp_path)
    out_tsv = tmp_path / "out.tsv"
    out_json = tmp_path / "out.json"
    out_tsv.write_text("DO NOT TOUCH")

    monkeypatch.setattr(DRT, "_require_output_does_not_already_exist", lambda *a, **k: None)
    DRT.main(
        roster_path=roster_path,
        grouping_path=grouping_path,
        counts_path=counts_path,
        coverage_summary_path=coverage_summary_path,
        out_tsv=str(out_tsv),
        out_json=str(out_json),
    )
    assert out_tsv.read_text() != "DO NOT TOUCH"


def test_build_design_readiness_table_against_real_committed_inputs():
    """Calls the pure computation directly against the real, already-committed
    roster/grouping/counts/coverage-summary tables -- never a synthetic
    stand-in, never writing any file.
    """
    table_rows, fieldnames, cutoff_sensitivity, summary = DRT.build_design_readiness_table()

    assert summary["n_multi_visit_donors"] == 15
    assert summary["n_donor_visit_rows"] == 31  # 14 x 2-visit + 1 x 3-visit
    assert summary["n_sample_ids_missing_from_counts_table"] == 0

    c8821_rows = [r for r in table_rows if r["true_donor_group"] == "sdrf::C-8821"]
    assert len(c8821_rows) == 2
    initial = next(r for r in c8821_rows if r["visit_resample_flag"] == "Initial")
    assert initial["n_CD83_CD14_mono_cells_raw"] == "10"
    assert initial["symptomatic_or_asymptomatic"] == "symptomatic"

    # CD83_CD14_mono at cutoff 1 must include all 12 symptomatic + 3 asymptomatic
    # multi-visit donors (per the already-reported universal coverage finding).
    cd83_1 = cutoff_sensitivity["CD83_CD14_mono"]["1"]
    assert len(cd83_1["symptomatic"]) == 12
    assert len(cd83_1["asymptomatic"]) == 3

    # At cutoff 50, C-8821 drops out (its D0 visit has only 10 CD83_CD14_mono cells).
    cd83_50 = cutoff_sensitivity["CD83_CD14_mono"]["50"]
    assert "sdrf::C-8821" not in {r["true_donor_group"] for r in cd83_50["symptomatic"]}


def test_qc_source_evidence_is_a_dated_scoped_record_not_a_blanket_claim():
    """Per review:20261016_120026 item 1: the status must be the
    accurate 'not verified for this pilot', never an implied global absence,
    and must carry a checked_date/checked_paths/findings -- not a bare string
    reprinted on every call as if freshly re-verified.
    """
    _, _, _, summary = DRT.build_design_readiness_table()
    evidence = summary["qc_source_evidence"]
    assert evidence["status"] == "qc_admission_not_verified_for_this_pilot"
    assert evidence["checked_date"] == "2026-10-16"
    assert evidence["checked_paths"]
    assert evidence["findings"]
    assert "gap" in evidence
    # Must not silently regress to the old overstated field name/claim.
    assert "qc_admitted_input_status" not in summary
