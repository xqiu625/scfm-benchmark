import os
import sys

RESULTS_DIR = os.path.join(
    os.path.dirname(__file__), "..", "scripts", "34_scfm_discovery_framework", "results"
)
sys.path.insert(0, RESULTS_DIR)

import stephenson_longitudinal_monocyte_coverage as COV  # noqa: E402

ROSTER_HEADER = (
    "sample_id\ttrue_donor_group\tvisit_resample_flag\tcollection_day\t"
    "own_severity_Status_on_day_collection_summary\n"
)
COUNTS_HEADER = "sample_id\tn_cells_this_sample_total\tA\tB\n"


def _write_roster(path, rows):
    with open(path, "w") as f:
        f.write(ROSTER_HEADER)
        for r in rows:
            f.write("\t".join(r) + "\n")


def _write_counts(path, rows):
    with open(path, "w") as f:
        f.write(COUNTS_HEADER)
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")


def test_build_coverage_separates_symptomatic_from_asymptomatic_and_finds_shared_subtype(tmp_path):
    roster_path = tmp_path / "roster.tsv"
    counts_path = tmp_path / "counts.tsv"

    _write_roster(
        roster_path,
        [
            ("sym1", "donor::sym", "Initial", "D0", "Moderate"),
            ("sym2", "donor::sym", "Resample", "D28", "Critical"),
            ("asym1", "donor::asym", "Initial", "D0", "Asymptomatic"),
            ("asym2", "donor::asym", "Resample", "D28", "Asymptomatic"),
            ("single1", "donor::single", "Initial", "D0", "Mild"),
        ],
    )
    _write_counts(
        counts_path,
        [
            # sym1: A and B present; sym2: A only -> shared = {A}
            ("sym1", 10, 6, 4),
            ("sym2", 5, 5, 0),
            ("asym1", 3, 1, 2),
            ("asym2", 3, 1, 2),
            ("single1", 7, 7, 0),
        ],
    )

    symptomatic, asymptomatic, summary = COV.build_coverage(str(roster_path), str(counts_path))

    assert len(symptomatic) == 1
    assert len(asymptomatic) == 1
    sym_row = symptomatic[0]
    assert sym_row["true_donor_group"] == "donor::sym"
    assert sym_row["coverage_status"] == "checked"
    assert sym_row["shared_subtypes_present_across_all_visits"] == ["A"]
    assert sym_row["any_visit_with_zero_monocyte_cells_in_this_object"] is False

    assert summary["n_symptomatic_multi_visit_donors"] == 1
    assert summary["n_asymptomatic_multi_visit_donors"] == 1
    assert summary["n_symptomatic_donors_with_at_least_one_shared_subtype"] == 1
    assert summary["n_symptomatic_donors_with_zero_shared_subtypes"] == 0


def test_build_coverage_reports_unknown_not_zero_for_a_visit_missing_from_counts_table(tmp_path):
    roster_path = tmp_path / "roster.tsv"
    counts_path = tmp_path / "counts.tsv"

    _write_roster(
        roster_path,
        [
            ("sym1", "donor::sym", "Initial", "D0", "Moderate"),
            ("sym2", "donor::sym", "Resample", "D28", "Critical"),
        ],
    )
    # sym2 deliberately absent from the counts table -- a mapping gap.
    _write_counts(
        counts_path,
        [
            ("sym1", 10, 6, 4),
        ],
    )

    symptomatic, _, summary = COV.build_coverage(str(roster_path), str(counts_path))

    row = symptomatic[0]
    assert row["coverage_status"] == "incomplete -- see per-visit rows"
    assert row["shared_subtypes_present_across_all_visits"] == (
        "unknown (at least one visit missing from the cell-count table)"
    )
    assert row["any_visit_with_zero_monocyte_cells_in_this_object"] == "unknown"

    sym2_visit = next(v for v in row["per_visit"] if v["sample_id"] == "sym2")
    assert sym2_visit["n_monocyte_cells_this_sample"] == "unknown"
    assert sym2_visit["subtypes_present_gt0"] == "unknown"

    assert summary["n_sample_visits_missing_from_cell_count_table"] == 1
    assert summary["n_symptomatic_donors_coverage_unknown"] == 1


def test_main_refuses_to_overwrite_an_existing_output_file(tmp_path):
    roster_path = tmp_path / "roster.tsv"
    counts_path = tmp_path / "counts.tsv"
    _write_roster(
        roster_path,
        [
            ("s1", "donor::a", "Initial", "D0", "Mild"),
            ("s2", "donor::a", "Resample", "D28", "Mild"),
        ],
    )
    _write_counts(counts_path, [("s1", 1, 1, 0), ("s2", 1, 1, 0)])

    out_tsv = tmp_path / "out.tsv"
    out_summary = tmp_path / "out_summary.json"
    out_tsv.write_text("DO NOT TOUCH")

    try:
        COV.main(
            roster_path=str(roster_path),
            counts_path=str(counts_path),
            out_tsv=str(out_tsv),
            out_summary=str(out_summary),
        )
        assert False, "expected SystemExit when out_tsv already exists"
    except SystemExit:
        pass
    assert out_tsv.read_text() == "DO NOT TOUCH"
    assert not out_summary.exists()


def test_mutation_kill_the_overwrite_guard(tmp_path, monkeypatch):
    roster_path = tmp_path / "roster.tsv"
    counts_path = tmp_path / "counts.tsv"
    _write_roster(
        roster_path,
        [
            ("s1", "donor::a", "Initial", "D0", "Mild"),
            ("s2", "donor::a", "Resample", "D28", "Mild"),
        ],
    )
    _write_counts(counts_path, [("s1", 1, 1, 0), ("s2", 1, 1, 0)])

    out_tsv = tmp_path / "out.tsv"
    out_summary = tmp_path / "out_summary.json"
    out_tsv.write_text("DO NOT TOUCH")

    monkeypatch.setattr(COV, "_require_output_does_not_already_exist", lambda *a, **k: None)
    COV.main(
        roster_path=str(roster_path),
        counts_path=str(counts_path),
        out_tsv=str(out_tsv),
        out_summary=str(out_summary),
    )
    assert out_tsv.read_text() != "DO NOT TOUCH"
