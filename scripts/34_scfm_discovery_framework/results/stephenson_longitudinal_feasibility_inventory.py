"""Per review:20261015_224111 item 2, hardened per
review:20261015_230801 (F1/F2/F4): a SMALL, metadata-only longitudinal
feasibility inventory -- never a read of the full expression matrix, never a
new statistical analysis. This answers "what repeated-visit structure
actually exists in this cohort" so the independent-donor DE view's
one-visit-per-donor simplification (already decided, already implemented in
`stephenson_icu_severity_grouping_resolve_mixed_visits.py`) is never mistaken
for "every other visit is deleted from this project." It does NOT decide
which longitudinal analysis to run -- see the companion feasibility note this
script's own summary feeds.

Joins the already-produced, already-reviewed `stephenson_donor_identity_
roster_v1.tsv` (true_donor_group, visit_resample_flag, collection_day,
days_from_onset, own severity) with `stephenson_icu_severity_grouping_v1.tsv`
(care_setting_at_sampling, candidate_main_group) by sample_id. Only reads
these two small, already-committed tables -- no raw h5ad, no expression
values.

v2 fixes three real defects the review found in v1 (`_v1.tsv`/`_v1_summary.
json`, left on disk UNCHANGED as the historical record of what that code
produced -- never deleted, never silently regenerated in place):

F1 -- v1's own test deleted and rewrote the committed v1 result files in
place, and `main()` unconditionally `open(..., "w")`-overwrote whatever was
at its hardcoded output path. This module now takes explicit input/output
paths (`build_inventory(roster_path=..., grouping_path=..., out_tsv=...,
out_summary=...)`), refuses outright (`SystemExit`) if either output path
already exists, and the CLI defaults to a NEW path
(`..._v2.tsv`/`..._v2_summary.json`), never the v1 path.

F2 -- v1 only checked whether the raw `days_from_onset` strings were
IDENTICAL across a donor's visits, then called that "identical" without
checking they were numeric. `sdrf::CV0231`'s two visits are both the literal
string `"Healthy"` -- identical, but not a day count -- and v1's logic
silently filed it under "identical" together with the 14 donors that have a
real identical integer, reporting `days_from_onset_unparseable=0` (false; it
is 1). This version keeps three independent booleans per donor (raw-string
identity, numeric validity, and the still-open question of whether the field
is re-measured per visit or copied from a single baseline) and classifies
into four disjoint buckets instead of three.

F4 -- v1 only ever compared the first and last visit, discarding whatever
happened at `sdrf::C-8896`'s middle visit, and would silently claim a
day-based order for any future donor whose `collection_day` cannot be
parsed (by sorting it to the end rather than refusing to order it). This
version records the full ordered trajectory (when order can actually be
resolved) and reports BOTH "first-to-last" and "any adjacent-step"
transition labels, and an explicit `order_resolved` flag per donor rather
than inferring a chronological position from an unparseable value or from
input row order. An empty multi-visit set now writes an empty TSV instead of
crashing on `donor_rows_out[0]`.

R2 (review:20261016_102834) -- the F4 fix above correctly gated
`severity_any_adjacent_transition`/`care_setting_any_adjacent_transition` on
`order_resolved`, but NOT the first-to-last fields or the summary's
first-to-last transition-count dicts: those still computed a direction from
whatever order the input rows happened to be in when `order_resolved` is
False. Reproduced on an isolated fixture (C-8821 with an unparseable
resample day): `order_resolved=False` yet
`severity_transition_first_to_last="worsening"` and the summary counted it
under `"Moderate->Critical": 1` as if the order were known. Fixed: when
`order_resolved` is False, `severity_first`/`severity_last`/
`severity_transition_first_to_last` (and the care-setting equivalents) are
now `"unknown"`, and the donor is counted under a separate
`"order_unresolved"` key in the summary's transition-count dicts rather than
under any directional label. This does not change any already-resolved real
donor's output -- the real roster has 0 order-unresolved multi-visit donors
(`n_multi_visit_donors_with_order_unresolved=0`), confirmed again after this
fix.
"""
import csv
import hashlib
import json
import os
import subprocess

HERE = os.path.dirname(__file__)
ROSTER_TSV = os.path.join(HERE, "stephenson_donor_identity_roster_v1.tsv")
GROUPING_TSV = os.path.join(HERE, "stephenson_icu_severity_grouping_v1.tsv")

#: v1 outputs -- historical record of the earlier (F1/F2/F4-affected) code's
#: result. Left on disk, never written to by this module again.
OUT_TSV_V1 = os.path.join(HERE, "stephenson_longitudinal_feasibility_inventory_v1.tsv")
OUT_SUMMARY_V1 = os.path.join(HERE, "stephenson_longitudinal_feasibility_inventory_v1_summary.json")

#: v2 outputs -- this module's own default target. Never overwrites v1, and
#: refuses to overwrite itself if already present (see `_require_output_does_
#: not_already_exist`).
OUT_TSV_V2 = os.path.join(HERE, "stephenson_longitudinal_feasibility_inventory_v2.tsv")
OUT_SUMMARY_V2 = os.path.join(HERE, "stephenson_longitudinal_feasibility_inventory_v2_summary.json")

#: Explicit ordinal scale, per Stephenson et al. 2021 Methods (Clinical status
#: assignment, PDF p.14) -- used ONLY to classify a transition's direction
#: (stable/worsening/improving), never to recompute group membership (that remains
#: `stephenson_icu_severity_grouping.py`'s job).
SEVERITY_ORDER = {"Asymptomatic": 0, "Mild": 1, "Moderate": 2, "Severe": 3, "Critical": 4}
CARE_SETTING_ORDER = {"non-ICU": 0, "ICU": 1}

FIELDNAMES = [
    "true_donor_group",
    "n_visits",
    "order_resolved",
    "sample_ids_in_order",
    "collection_days_in_order",
    "severity_sequence",
    "care_setting_sequence",
    "severity_first",
    "severity_last",
    "severity_transition_first_to_last",
    "severity_any_adjacent_transition",
    "care_setting_first",
    "care_setting_last",
    "care_setting_transition_first_to_last",
    "care_setting_any_adjacent_transition",
    "days_from_onset_raw_values",
    "days_from_onset_raw_strings_identical",
    "days_from_onset_all_numeric",
    "days_from_onset_category",
    "days_from_onset_note",
    "same_monocyte_state_availability_checked_this_round",
]


def _collection_day_sort_key(collection_day):
    """'D0'/'D7'/'D9'/'D12'/'D13'/'D28' -> an int for chronological sort within one
    donor. Never compared ACROSS donors (D0 is a schedule label, not a shared
    calendar day) -- see MONOCYTE_PILOT_PREREG_v1.md Sec.4.
    """
    if collection_day.startswith("D") and collection_day[1:].isdigit():
        return int(collection_day[1:])
    return None


def _load_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _git_head_or_unknown():
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=HERE, capture_output=True, text=True, timeout=10,
        )
        if out.returncode == 0:
            return out.stdout.strip()
    except OSError:
        pass
    return "unknown"


def _transition_label(first_value, last_value, order):
    if first_value not in order or last_value not in order:
        return "unknown"
    if order[last_value] > order[first_value]:
        return "worsening"
    if order[last_value] < order[first_value]:
        return "improving"
    return "stable"


def _any_adjacent_transition(sequence, order):
    """True if any CONSECUTIVE pair in an ORDER-RESOLVED sequence differs,
    even if the first and last entries happen to match (a donor can go
    Mild -> Severe -> Mild and still read 'stable' first-to-last). Returns
    'unknown' (not True/False) when the order itself is not resolved, so an
    unresolved sequence is never silently read as a stable one.
    """
    labels = [_transition_label(a, b, order) for a, b in zip(sequence[:-1], sequence[1:])]
    if any(label == "unknown" for label in labels):
        return "unknown"
    return any(label != "stable" for label in labels)


def _order_donor_visits(rows):
    """Order one donor's visit rows chronologically when that is actually
    resolvable, and say explicitly when it is not -- never silently fall
    back to input row order and call it a timeline.

    Resolved iff: exactly one row is flagged 'Initial', and every other
    ('Resample') row has a `collection_day` this module can parse into a
    distinct day number (ties or unparseable days make relative order among
    the Resamples unknown).
    """
    initials = [r for r in rows if r["visit_resample_flag"] == "Initial"]
    resamples = [r for r in rows if r["visit_resample_flag"] != "Initial"]
    if len(initials) != 1:
        return list(rows), False

    resample_keys = [_collection_day_sort_key(r["collection_day"]) for r in resamples]
    if any(k is None for k in resample_keys) or len(set(resample_keys)) != len(resample_keys):
        return list(rows), False

    ordered_resamples = [r for _, r in sorted(zip(resample_keys, resamples), key=lambda t: t[0])]
    return initials + ordered_resamples, True


def _classify_onset(raw_values):
    """Three independent facts, never collapsed into one: are the raw
    strings identical; are they all numeric; and (never answered here --
    named as an open gap) whether the field is re-measured per visit.
    """
    raw_identical = len(set(raw_values)) == 1
    all_numeric = all(v.lstrip("-").isdigit() for v in raw_values)

    if raw_identical and all_numeric:
        category = "identical_numeric"
        note = (
            "identical numeric string across all this donor's visits -- consistent with (not "
            "proof of) a baseline-only field; per-visit re-measurement semantics remain "
            "unconfirmed (see MONOCYTE_PILOT_PREREG_v1.md Sec.4)"
        )
    elif raw_identical and not all_numeric:
        category = "identical_non_numeric"
        note = (
            f"identical NON-NUMERIC string across all this donor's visits ({sorted(raw_values)[0]!r}) "
            "-- not a day count; must not be counted as a resolved onset day, and this donor's "
            "'identical' status says nothing about disease-onset timing"
        )
    elif not raw_identical and all_numeric:
        category = "differs_numeric"
        note = f"numeric strings differ across visits: {sorted(raw_values)}"
    else:
        category = "differs_mixed_or_non_numeric"
        note = f"differing, at least one non-numeric value present: {sorted(raw_values)}"

    return raw_identical, all_numeric, category, note


def build_inventory(roster_path=ROSTER_TSV, grouping_path=GROUPING_TSV):
    """Pure computation (no file writes): returns (donor_rows_out, summary_without_
    provenance). Separated from `main()` so tests can exercise it against
    arbitrary fixture paths without touching any committed result file.
    """
    roster_rows = _load_tsv(roster_path)
    grouping_by_sample = {r["sample_id"]: r for r in _load_tsv(grouping_path)}

    missing_in_grouping = sorted(
        r["sample_id"] for r in roster_rows if r["sample_id"] not in grouping_by_sample)
    if missing_in_grouping:
        raise SystemExit(
            f"refusing to build the inventory: {len(missing_in_grouping)} roster sample_id(s) "
            f"have no row in {grouping_path}: {missing_in_grouping[:5]}"
        )

    by_donor = {}
    for r in roster_rows:
        by_donor.setdefault(r["true_donor_group"], []).append(r)

    n_total_donors = len(by_donor)
    multi_visit = {d: rs for d, rs in by_donor.items() if len(rs) > 1}
    n_single_visit = n_total_donors - len(multi_visit)

    donor_rows_out = []
    severity_transition_counts = {}
    care_setting_transition_counts = {}
    onset_category_counts = {}
    visit_count_distribution = {}
    n_order_unresolved = 0

    for donor, rows in sorted(multi_visit.items()):
        visit_count_distribution[len(rows)] = visit_count_distribution.get(len(rows), 0) + 1

        ordered, order_resolved = _order_donor_visits(rows)
        if not order_resolved:
            n_order_unresolved += 1

        # `ordered` is the donor's rows in input-row order when `order_resolved`
        # is False (see `_order_donor_visits`) -- these sequences are raw,
        # input-order values in that case, NOT a chronology, and must never be
        # used to compute a direction (R2, review:20261016_102834).
        severity_sequence = [r["own_severity_Status_on_day_collection_summary"] for r in ordered]
        care_sequence = [grouping_by_sample[r["sample_id"]]["care_setting_at_sampling"] for r in ordered]

        if order_resolved:
            severity_first, severity_last = severity_sequence[0], severity_sequence[-1]
            severity_label = _transition_label(severity_first, severity_last, SEVERITY_ORDER)
            severity_transition_counts[f"{severity_first}->{severity_last}"] = (
                severity_transition_counts.get(f"{severity_first}->{severity_last}", 0) + 1)

            care_first, care_last = care_sequence[0], care_sequence[-1]
            care_label = _transition_label(care_first, care_last, CARE_SETTING_ORDER)
            care_setting_transition_counts[f"{care_first}->{care_last}"] = (
                care_setting_transition_counts.get(f"{care_first}->{care_last}", 0) + 1)

            severity_any = _any_adjacent_transition(severity_sequence, SEVERITY_ORDER)
            care_any = _any_adjacent_transition(care_sequence, CARE_SETTING_ORDER)
        else:
            # Order unresolved: a direction (worsening/improving/stable) requires
            # knowing which visit came first, which this donor's data does not
            # resolve. Report "unknown" rather than a direction computed from
            # arbitrary input row order, and count this donor separately from
            # any real, resolved-order transition (R2 fix).
            severity_first = severity_last = "unknown"
            severity_label = "unknown"
            care_first = care_last = "unknown"
            care_label = "unknown"
            severity_any = "unknown"
            care_any = "unknown"
            severity_transition_counts["order_unresolved"] = (
                severity_transition_counts.get("order_unresolved", 0) + 1)
            care_setting_transition_counts["order_unresolved"] = (
                care_setting_transition_counts.get("order_unresolved", 0) + 1)

        raw_onset_values = [r["days_from_onset"] for r in rows]
        onset_identical, onset_numeric, onset_category, onset_note = _classify_onset(raw_onset_values)
        onset_category_counts[onset_category] = onset_category_counts.get(onset_category, 0) + 1

        subtype_checked = not all(
            r["subtype_coverage_checked_this_round"].strip().lower().startswith("not re-checked")
            for r in rows)

        donor_rows_out.append({
            "true_donor_group": donor,
            "n_visits": len(rows),
            "order_resolved": order_resolved,
            "sample_ids_in_order": ",".join(r["sample_id"] for r in ordered),
            "collection_days_in_order": ",".join(r["collection_day"] for r in ordered),
            "severity_sequence": ",".join(severity_sequence),
            "care_setting_sequence": ",".join(care_sequence),
            "severity_first": severity_first,
            "severity_last": severity_last,
            "severity_transition_first_to_last": severity_label,
            "severity_any_adjacent_transition": severity_any,
            "care_setting_first": care_first,
            "care_setting_last": care_last,
            "care_setting_transition_first_to_last": care_label,
            "care_setting_any_adjacent_transition": care_any,
            "days_from_onset_raw_values": ",".join(raw_onset_values),
            "days_from_onset_raw_strings_identical": onset_identical,
            "days_from_onset_all_numeric": onset_numeric,
            "days_from_onset_category": onset_category,
            "days_from_onset_note": onset_note,
            "same_monocyte_state_availability_checked_this_round": subtype_checked,
        })

    summary = {
        "n_total_true_donors": n_total_donors,
        "n_single_visit_donors": n_single_visit,
        "n_multi_visit_donors": len(multi_visit),
        "n_multi_visit_donors_with_order_unresolved": n_order_unresolved,
        "visit_count_distribution_among_multi_visit_donors": visit_count_distribution,
        "severity_transition_counts_first_to_last_visit": severity_transition_counts,
        "care_setting_transition_counts_first_to_last_visit": care_setting_transition_counts,
        "days_from_onset_category_counts": onset_category_counts,
        "same_monocyte_state_subtype_availability_across_visits": (
            "NOT independently checked this round for any multi-visit donor -- "
            "the roster's own subtype_coverage_checked_this_round column already "
            "says 'not re-checked this round' for every row; this is a real, named "
            "gap for the longitudinal feasibility question, not an assumed yes/no"
        ),
        "note": (
            "Severity-tier transition and main-comparison-ARM transition are DIFFERENT questions "
            "and must not be conflated -- see per-donor rows for the full trajectory, not just "
            "first-to-last. 'days_from_onset_category_counts' separates raw-string identity from "
            "numeric validity: 'identical_non_numeric' (e.g. CV0231's Healthy/Healthy) must never "
            "be read as a resolved onset day, and an identical NUMERIC string is still only "
            "consistent with, not proof of, a baseline-only (non-re-measured) field. "
            "'order_resolved'=False means this donor's visit order could not be determined from "
            "visit_resample_flag + collection_day alone; its 'sequence' fields reflect raw "
            "input row order ONLY and must not be read as chronological, and its "
            "'severity_first'/'severity_last'/'severity_transition_first_to_last' (and the "
            "care_setting equivalents) are reported as 'unknown' rather than computed from that "
            "input order -- such a donor is counted under the 'order_unresolved' key in the "
            "transition-count dicts below, never under a directional label (R2 fix, "
            "codex_output/mono_20261016_102834.md)."
        ),
    }
    return donor_rows_out, summary


def _require_output_does_not_already_exist(out_tsv, out_summary):
    existing = [p for p in (out_tsv, out_summary) if os.path.exists(p)]
    if existing:
        raise SystemExit(
            "refusing to overwrite an existing inventory output -- rerun against a new "
            f"out_tsv/out_summary path (existing: {existing})"
        )


def main(roster_path=ROSTER_TSV, grouping_path=GROUPING_TSV, out_tsv=OUT_TSV_V2, out_summary=OUT_SUMMARY_V2):
    _require_output_does_not_already_exist(out_tsv, out_summary)

    donor_rows_out, summary = build_inventory(roster_path=roster_path, grouping_path=grouping_path)

    with open(out_tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(donor_rows_out)

    summary["source_files"] = {
        "roster_tsv": roster_path,
        "roster_tsv_sha256": _sha256_file(roster_path),
        "grouping_tsv": grouping_path,
        "grouping_tsv_sha256": _sha256_file(grouping_path),
    }
    summary["code_revision"] = {
        "git_head": _git_head_or_unknown(),
        "script_sha256": _sha256_file(os.path.abspath(__file__)),
    }
    summary["output_files"] = {"out_tsv": out_tsv, "out_summary": out_summary}

    with open(out_summary, "w") as f:
        json.dump(summary, f, indent=2, sort_keys=True)
        f.write("\n")

    print(json.dumps(summary, indent=2, sort_keys=True))
    return summary


if __name__ == "__main__":
    main()
