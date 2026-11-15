"""Narrow regression check (not a framework): fails if the committed
stephenson_donor_identity_roster_v1.tsv and stephenson_sdrf_identity_crosswalk_summary.json
ever disagree about the identity-unit count, or if the roster's row/column shape drifts.
Run this any time either file is hand-edited without regenerating the other via
stephenson_sdrf_identity_crosswalk.py.
"""
import csv
import json
import sys

SUMMARY = "stephenson_sdrf_identity_crosswalk_summary.json"
ROSTER = "stephenson_donor_identity_roster_v1.tsv"
EXPECTED_COLUMNS = 15
EXPECTED_ROWS = 102


def main() -> int:
    summary = json.load(open(SUMMARY))
    with open(ROSTER, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    problems = []
    if len(rows) != EXPECTED_ROWS:
        problems.append(f"roster has {len(rows)} data rows, expected {EXPECTED_ROWS}")
    if rows and len(rows[0]) != EXPECTED_COLUMNS:
        problems.append(f"roster has {len(rows[0])} columns, expected {EXPECTED_COLUMNS}")

    roster_groups = {r["true_donor_group"] for r in rows}
    expected = summary["sdrf_confirmed_identity_unit_count"]
    if len(roster_groups) != expected:
        problems.append(
            f"roster has {len(roster_groups)} distinct true_donor_group values, "
            f"summary.sdrf_confirmed_identity_unit_count={expected}"
        )

    n_roster_unresolved = sum(1 for r in rows if r["sdrf_individual_id"] == "n/a_unresolved")
    expected_unresolved = summary["n_samples_with_no_resolvable_sdrf_individual"]
    if n_roster_unresolved != expected_unresolved:
        problems.append(
            f"roster has {n_roster_unresolved} rows marked n/a_unresolved, "
            f"summary.n_samples_with_no_resolvable_sdrf_individual={expected_unresolved}"
        )

    if problems:
        print("FAIL:")
        for p in problems:
            print(f"  - {p}")
        return 1
    print(f"PASS: roster ({len(rows)} rows, {len(roster_groups)} donor groups) and summary "
         f"(sdrf_confirmed_identity_unit_count={expected}) agree.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
