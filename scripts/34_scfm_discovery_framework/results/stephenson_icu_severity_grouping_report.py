"""Builds the reviewable sample/donor list Cindy's ICU-at-sampling decision asked for
(review:20261015_120317 Sec: "给 implementation assistant 的直接接续指令", item 3).

Joins the real per-sample clinical fields (stephenson_sample_level_real_fields.csv,
`Status_on_day_collection` / `Status_on_day_collection_summary`) with the SDRF-confirmed
donor-identity roster (stephenson_donor_identity_roster_v1.tsv, `true_donor_group`,
`identity_status`, `visit_resample_flag`) via `stephenson_icu_severity_grouping.py`'s
pure mapping functions. Produces:

  - stephenson_icu_severity_grouping_v1.tsv: one row per SAMPLE, with raw severity,
    raw care-setting source value, derived care_setting_at_sampling, candidate group,
    and the reason -- never silently defaulting missing/conflicting evidence.
  - stephenson_icu_severity_grouping_v1_summary.json: SAMPLE-level counts per group, and
    DONOR-level counts (deduplicated via true_donor_group, using each donor's Initial-
    preferred sample where available) split into: eligible comparison, confirmed
    non-ICU-excluded, unknown-location-excluded, reference, excluded-other.

This is evidence preparation only -- it does NOT freeze repeated-visit handling (a donor
with both an eligible and an ineligible visit is reported, not resolved here; see the
"next pending" note in the summary) and does NOT itself admit any row as a final analysis
input.
"""
import csv
import json
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))
from stephenson_icu_severity_grouping import care_setting_at_sampling, candidate_main_group

HERE = os.path.dirname(__file__)
REAL_FIELDS_CSV = os.path.join(HERE, "stephenson_sample_level_real_fields.csv")
ROSTER_TSV = os.path.join(HERE, "stephenson_donor_identity_roster_v1.tsv")
OUT_TSV = os.path.join(HERE, "stephenson_icu_severity_grouping_v1.tsv")
OUT_SUMMARY = os.path.join(HERE, "stephenson_icu_severity_grouping_v1_summary.json")


def _load_real_fields():
    with open(REAL_FIELDS_CSV, newline="") as f:
        return {row["sample_id"]: row for row in csv.DictReader(f)}


def _load_roster():
    with open(ROSTER_TSV, newline="") as f:
        return {row["sample_id"]: row for row in csv.DictReader(f, delimiter="\t")}


def main():
    real_fields = _load_real_fields()
    roster = _load_roster()

    missing_in_roster = sorted(set(real_fields) - set(roster))
    missing_in_real_fields = sorted(set(roster) - set(real_fields))
    if missing_in_roster or missing_in_real_fields:
        raise SystemExit(
            f"refusing to produce a grouping report: sample_id sets disagree between "
            f"{REAL_FIELDS_CSV} and {ROSTER_TSV} -- "
            f"missing_in_roster={missing_in_roster} missing_in_real_fields={missing_in_real_fields}"
        )

    out_rows = []
    sample_group_counts = {}
    for sample_id in sorted(real_fields):
        rf = real_fields[sample_id]
        rr = roster[sample_id]
        raw_status = rf["Status_on_day_collection"]
        severity_summary = rf["Status_on_day_collection_summary"]
        setting = care_setting_at_sampling(raw_status)
        group, reason = candidate_main_group(severity_summary, setting)
        sample_group_counts[group] = sample_group_counts.get(group, 0) + 1
        out_rows.append(
            {
                "sample_id": sample_id,
                "true_donor_group": rr["true_donor_group"],
                "identity_status": rr["identity_status"],
                "visit_resample_flag": rr["visit_resample_flag"],
                "collection_day": rr["collection_day"],
                "severity_at_sampling_original": severity_summary,
                "status_on_day_collection_raw": raw_status,
                "care_setting_at_sampling": setting,
                "candidate_main_group": group,
                "inclusion_exclusion_reason": reason,
            }
        )

    fieldnames = list(out_rows[0].keys())
    with open(OUT_TSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    # Donor-level dedup: a donor is eligible for the comparison arm only if AT LEAST ONE
    # of their samples lands in the comparison group; flag donors whose samples disagree
    # across visits (e.g. Initial Moderate -> later ICU Critical) rather than resolving
    # the visit-selection question here.
    donor_groups = {}
    for row in out_rows:
        donor_groups.setdefault(row["true_donor_group"], []).append(row["candidate_main_group"])

    donor_level_counts = {}
    donors_with_mixed_group_across_visits = []
    for donor, groups in donor_groups.items():
        distinct_groups = sorted(set(groups))
        if len(distinct_groups) > 1:
            donors_with_mixed_group_across_visits.append(
                {"true_donor_group": donor, "groups_across_visits": distinct_groups}
            )
        # count each donor once per distinct group they ever appear under this round,
        # explicitly NOT collapsing to one row -- visit selection is a separate open decision
        for g in distinct_groups:
            donor_level_counts[g] = donor_level_counts.get(g, 0) + 1

    summary = {
        "source_files": {"real_fields_csv": REAL_FIELDS_CSV, "roster_tsv": ROSTER_TSV},
        "n_samples": len(out_rows),
        "n_distinct_true_donor_groups": len(donor_groups),
        "sample_level_group_counts": sample_group_counts,
        "donor_level_group_counts_any_visit": donor_level_counts,
        "n_donors_with_mixed_group_across_visits": len(donors_with_mixed_group_across_visits),
        "donors_with_mixed_group_across_visits": donors_with_mixed_group_across_visits,
        "note": (
            "donor_level_group_counts_any_visit counts a donor once per distinct candidate "
            "group they have ANY sample in -- it is NOT a final per-donor analysis N, because "
            "repeated-visit handling (which single visit represents a donor) is still "
            "pending_owner. See docs/decisions.md 2026-10-15 entries and "
            "MONOCYTE_PILOT_PREREG_v1.md Sec.6 for that open question."
        ),
    }
    with open(OUT_SUMMARY, "w") as f:
        json.dump(summary, f, indent=2, sort_keys=True)

    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
