"""Per review:20261015_230801 increment 2: turn the metadata-only
repeated-visit inventory (`stephenson_longitudinal_feasibility_inventory.py`)
into an ACTIONABLE coverage table -- does each multi-visit donor's repeated
visit actually share a usable monocyte subtype, using the real per-sample
subtype cell counts (`stephenson_sample_subtype_cell_counts_v1.tsv`, itself a
bounded `obs`-only read of the frozen per-cohort raw object, never `X`).

Reports the 12 SYMPTOMATIC multi-visit donors (excludes the 3 Asymptomatic
multi-visit donors, which are out of scope for the current Mild/Moderate vs.
ICU-Severe/Critical comparison and are listed separately, never silently
merged in or silently dropped) plus, separately, the 3 Asymptomatic ones.

A sample absent from the cell-count table is `unknown` -- a mapping gap, never
read as 0 cells. A sample present with a genuine 0 count for a subtype is a
real, confirmed zero FOR THIS SPECIFIC FILE/VERSION (qc_status: provisional,
same caveat the cell-count script and `monocyte_pilot_feasibility_v3.tsv`'s
Stephenson row already carry) -- not a final QC-admitted count, and not proof
the subtype is biologically absent from the patient.

This script does NOT decide which longitudinal comparison to run, does NOT
fit anything, and does NOT touch the raw expression matrix.
"""
import csv
import json
import os

HERE = os.path.dirname(__file__)
ROSTER_TSV = os.path.join(HERE, "stephenson_donor_identity_roster_v1.tsv")
COUNTS_TSV = os.path.join(HERE, "stephenson_sample_subtype_cell_counts_v1.tsv")
OUT_TSV = os.path.join(HERE, "stephenson_longitudinal_monocyte_coverage_v1.tsv")
OUT_SUMMARY = os.path.join(HERE, "stephenson_longitudinal_monocyte_coverage_v1_summary.json")

#: Same scale as stephenson_longitudinal_feasibility_inventory.py -- reused by
#: value, not reimported, so this module has no runtime dependency on that
#: module's own output format.
ASYMPTOMATIC_LABEL = "Asymptomatic"


def _load_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _load_counts(counts_path):
    rows = _load_tsv(counts_path)
    subtype_fields = [k for k in rows[0].keys() if k not in ("sample_id", "n_cells_this_sample_total")]
    by_sample = {}
    for r in rows:
        by_sample[r["sample_id"]] = {
            "n_cells_total": int(r["n_cells_this_sample_total"]),
            "subtype_counts": {st: int(r[st]) for st in subtype_fields},
        }
    return by_sample, subtype_fields


def _present_subtypes(sample_counts, subtype_fields):
    """Subtypes with count > 0 for this one sample. Returns None (never an
    empty-set substitute) if the sample has no row in the count table at
    all -- a mapping gap must never be read as 'zero of everything'."""
    if sample_counts is None:
        return None
    return {st for st in subtype_fields if sample_counts["subtype_counts"][st] > 0}


def build_coverage(roster_path=ROSTER_TSV, counts_path=COUNTS_TSV):
    roster_rows = _load_tsv(roster_path)
    counts_by_sample, subtype_fields = _load_counts(counts_path)

    by_donor = {}
    for r in roster_rows:
        by_donor.setdefault(r["true_donor_group"], []).append(r)
    multi_visit = {d: rs for d, rs in by_donor.items() if len(rs) > 1}

    symptomatic_rows = []
    asymptomatic_rows = []
    n_samples_missing_from_counts_table = 0

    for donor, rows in sorted(multi_visit.items()):
        is_asymptomatic = all(
            r["own_severity_Status_on_day_collection_summary"] == ASYMPTOMATIC_LABEL for r in rows)

        per_visit = []
        present_sets = []
        for r in rows:
            sc = counts_by_sample.get(r["sample_id"])
            if sc is None:
                n_samples_missing_from_counts_table += 1
            present = _present_subtypes(sc, subtype_fields)
            per_visit.append({
                "sample_id": r["sample_id"],
                "visit_resample_flag": r["visit_resample_flag"],
                "collection_day": r["collection_day"],
                "n_monocyte_cells_this_sample": "unknown" if sc is None else sc["n_cells_total"],
                "subtypes_present_gt0": "unknown" if present is None else sorted(present),
            })
            if present is not None:
                present_sets.append(present)

        if len(present_sets) == len(rows):
            shared_across_all_visits = sorted(set.intersection(*present_sets)) if present_sets else []
            any_visit_zero_monocyte_cells = any(v["n_monocyte_cells_this_sample"] == 0 for v in per_visit)
            coverage_status = "checked"
        else:
            shared_across_all_visits = "unknown (at least one visit missing from the cell-count table)"
            any_visit_zero_monocyte_cells = "unknown"
            coverage_status = "incomplete -- see per-visit rows"

        out_row = {
            "true_donor_group": donor,
            "n_visits": len(rows),
            "coverage_status": coverage_status,
            "shared_subtypes_present_across_all_visits": shared_across_all_visits,
            "any_visit_with_zero_monocyte_cells_in_this_object": any_visit_zero_monocyte_cells,
            "per_visit": per_visit,
        }
        (asymptomatic_rows if is_asymptomatic else symptomatic_rows).append(out_row)

    summary = {
        "n_symptomatic_multi_visit_donors": len(symptomatic_rows),
        "n_asymptomatic_multi_visit_donors": len(asymptomatic_rows),
        "n_sample_visits_missing_from_cell_count_table": n_samples_missing_from_counts_table,
        "subtypes_in_source_object": subtype_fields,
        "n_symptomatic_donors_with_at_least_one_shared_subtype": sum(
            1 for r in symptomatic_rows
            if isinstance(r["shared_subtypes_present_across_all_visits"], list)
            and len(r["shared_subtypes_present_across_all_visits"]) > 0
        ),
        "n_symptomatic_donors_with_zero_shared_subtypes": sum(
            1 for r in symptomatic_rows
            if isinstance(r["shared_subtypes_present_across_all_visits"], list)
            and len(r["shared_subtypes_present_across_all_visits"]) == 0
        ),
        "n_symptomatic_donors_coverage_unknown": sum(
            1 for r in symptomatic_rows if r["coverage_status"] != "checked"
        ),
        "provenance_caveat": (
            "subtype presence (>0 cells) is read from the per-cohort RAW object's own "
            "obs['orig_celltype'] (curated_atlas_label), qc_status: provisional -- this is a "
            "feasibility coverage check, not a QC-admitted analysis input. A sample absent from "
            "the cell-count table is reported 'unknown', never folded into a 0 count."
        ),
    }
    return symptomatic_rows, asymptomatic_rows, summary


def _require_output_does_not_already_exist(out_tsv, out_summary):
    existing = [p for p in (out_tsv, out_summary) if os.path.exists(p)]
    if existing:
        raise SystemExit(f"refusing to overwrite an existing coverage output: {existing}")


def main(roster_path=ROSTER_TSV, counts_path=COUNTS_TSV, out_tsv=OUT_TSV, out_summary=OUT_SUMMARY):
    _require_output_does_not_already_exist(out_tsv, out_summary)

    symptomatic_rows, asymptomatic_rows, summary = build_coverage(roster_path, counts_path)

    fieldnames = [
        "donor_category", "true_donor_group", "n_visits", "coverage_status",
        "shared_subtypes_present_across_all_visits",
        "any_visit_with_zero_monocyte_cells_in_this_object", "per_visit_json",
    ]
    with open(out_tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for category, rows in (("symptomatic", symptomatic_rows), ("asymptomatic", asymptomatic_rows)):
            for r in rows:
                w.writerow({
                    "donor_category": category,
                    "true_donor_group": r["true_donor_group"],
                    "n_visits": r["n_visits"],
                    "coverage_status": r["coverage_status"],
                    "shared_subtypes_present_across_all_visits": json.dumps(r["shared_subtypes_present_across_all_visits"]),
                    "any_visit_with_zero_monocyte_cells_in_this_object": r["any_visit_with_zero_monocyte_cells_in_this_object"],
                    "per_visit_json": json.dumps(r["per_visit"]),
                })

    summary["symptomatic_donor_rows"] = symptomatic_rows
    summary["asymptomatic_donor_rows"] = asymptomatic_rows
    with open(out_summary, "w") as f:
        json.dump(summary, f, indent=2, sort_keys=True)
        f.write("\n")

    print(json.dumps(summary, indent=2, sort_keys=True))
    return summary


if __name__ == "__main__":
    main()
