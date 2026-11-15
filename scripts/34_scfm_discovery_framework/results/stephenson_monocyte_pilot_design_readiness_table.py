"""Per review:20261016_102834 item 6 (implementation assistant instruction,
items 1-3): a minimal "monocyte real-pilot input and design-readiness table" --
NOT a new A/B/C menu, NOT a model run, NOT a new h5ad read.

Reuses only the three already-committed, already-verified small tables:
  - stephenson_donor_identity_roster_v1.tsv (true_donor_group, visit labels,
    severity, collection_day)
  - stephenson_icu_severity_grouping_v1.tsv (care_setting_at_sampling, and the
    raw ICU/ward string `status_on_day_collection_raw` that evidence is drawn
    from -- not a new ICU judgment)
  - stephenson_sample_subtype_cell_counts_v1.tsv (the already-produced,
    already-hash-verified per-sample x per-subtype RAW cell counts from
    `stephenson_monocyte_raw.h5ad`)

Produces two things (v1 -- unchanged on disk, superseded as the default target
only -- was the first version; v2 is this module's current default):

1. A donor x visit x subtype table (`..._donor_visit_subtype_v2.tsv`): one row
   per (donor, visit), with every subtype's raw count as its own column, plus
   severity, care-setting/ICU evidence, visit label, true_donor_group, and
   `symptomatic_or_asymptomatic` (from the already-produced coverage summary's
   own classification -- not recomputed here).

2. A per-subtype sensitivity cutoff listing (`..._subtype_cutoff_sensitivity_
   v2.json`): for each of the 5 real subtypes and each of the cutoffs
   {1, 20, 50, 100}, which donors have EVERY one of their visits at or above
   that raw count for that subtype -- split into symptomatic/asymptomatic.
   This is an explicitly-labelled SENSITIVITY DESCRIPTION, never a new QC/
   inclusion threshold (see `provenance_caveat` in the output and
   mono_LONGITUDINAL_FEASIBILITY_NOTE_v1.md Sec.1).

Scope boundary (explicit, per the review): this script does NOT touch the raw
h5ad's X, does NOT decide a QC-admitted monocyte input exists, and does NOT
pick a subtype or donor set for any DE run.

v2 (per review:20261016_120026 item 1): v1's `QC_ADMITTED_
INPUT_STATUS` was a fixed string, written into every output as if `build_
design_readiness_table()` had freshly re-verified it -- it had not; the
function never touches any QC artifact. It also understated the real
evidence: a repository-only search misses the HPCC data directories where
the actual QC/decontamination work lives. `QC_SOURCE_EVIDENCE` below
replaces it with a dated, scoped record of the ACTUAL bounded checks run on
2026-10-16 (commands in the module docstring of the companion review, see
`docs/current-status.md`'s entry for this round) -- not a live scan this
module repeats on every call. Status renamed from an implied "confirmed
absent" to the accurate `qc_admission_not_verified_for_this_pilot`: real
Gate 1 QC flags and DecontX contamination scores DO exist and DO join
100% by cell_id to every one of this pilot's 128,503 monocyte cells, but no
decision exists that any of this constitutes admission for the monocyte
pilot specifically, and the decontaminated COUNTS matrix's actual values
(as opposed to the per-cell contamination score) were not read this round.
"""
import csv
import hashlib
import json
import os

HERE = os.path.dirname(__file__)
ROSTER_TSV = os.path.join(HERE, "stephenson_donor_identity_roster_v1.tsv")
GROUPING_TSV = os.path.join(HERE, "stephenson_icu_severity_grouping_v1.tsv")
COUNTS_TSV = os.path.join(HERE, "stephenson_sample_subtype_cell_counts_v1.tsv")
COVERAGE_SUMMARY_JSON = os.path.join(HERE, "stephenson_longitudinal_monocyte_coverage_v1_summary.json")

OUT_DONOR_VISIT_SUBTYPE_TSV = os.path.join(
    HERE, "stephenson_monocyte_pilot_design_readiness_donor_visit_subtype_v2.tsv")
OUT_CUTOFF_SENSITIVITY_JSON = os.path.join(
    HERE, "stephenson_monocyte_pilot_design_readiness_subtype_cutoff_sensitivity_v2.json")

CUTOFFS = (1, 20, 50, 100)

#: A DATED, SCOPED EVIDENCE RECORD -- not a live check this module re-runs.
#: Every number here was produced by `stephenson_qc_source_crosscheck.py`'s `main()`
#: (a committed, tested, re-runnable script -- NOT prose) run directly on HPCC on
#: 2026-10-16; see its own output `stephenson_qc_source_crosscheck_v1.json` and
#: `STEPHENSON_QC_SOURCE_CROSSCHECK_v1.md` for the narrative read of the same numbers.
#: Re-running `main()` HERE does NOT re-verify any of this -- it reports this
#: recorded finding, exactly as `source_files` below reports hashes computed
#: once, not on every read.
QC_SOURCE_EVIDENCE = {
    "status": "qc_admission_not_verified_for_this_pilot",
    "checked_date": "2026-10-16",
    "checked_paths": [
        "/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c_monocyte/stephenson_monocyte_raw.h5ad",
        "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/stephenson_decontx_input/"
        "{cells.csv,genes.csv,stephenson_decontx_contamination_scores.csv} (headers/full small CSVs; "
        "counts.mtx and stephenson_decontx_counts.mtx NOT read -- 12G/25G, bounded-read only)",
        "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/gate1_qc/stephenson/"
        "qc_flags_per_cell.csv (full read, 132MB, cell-level QC flags)",
        "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/gate2_scrublet/ (ls only)",
        "scripts/31_donor_level_revision/code/build_monocyte_canonical.py (read, not executed)",
        "scripts/32_rebuild_upfront_qc/code/build_stephenson_decontx_input.py (read, not executed)",
        "scripts/32_rebuild_upfront_qc/qc_caveats.md, qc_checklist.md (read)",
    ],
    "findings": [
        "Gate 1 QC is APPLIED for Stephenson (gate1_qc/stephenson/DONE = 'gate1 applied', dated "
        "2026-09-21) with a per-cell qc_flags_per_cell.csv carrying flag_would_delete/n_flags.",
        "DecontX is COMPLETE for Stephenson (qc_caveats.md: 'DecontX complete: mean contamination "
        "13.0%, median 2.3%', dated 2026-09-18) with a per-cell contamination score file.",
        "cell_id join CONFIRMED 100%: stripping the 'steph_' prefix from stephenson_monocyte_raw."
        "h5ad's obs index gives an EXACT cell_id match, with 0 mismatches, for ALL 128,503 "
        "monocyte-pilot cells against BOTH gate1_qc/stephenson/qc_flags_per_cell.csv and "
        "stephenson_decontx_input/stephenson_decontx_contamination_scores.csv (channel_id/sample_id "
        "also verified consistent for every matched cell).",
        "Gate 1 flag_would_delete rate by monocyte subtype (real, this join): C1_CD16_mono 140/3942 "
        "(3.6%), CD14_mono 349/58562 (0.6%), CD16_mono 114/13619 (0.8%), CD83_CD14_mono 425/51735 "
        "(0.8%), Mono_prolif 7/645 (1.1%).",
        "DecontX contamination by monocyte subtype (real, this join, MEAN ESTIMATED score under "
        "the historical DecontX fit -- not an independently confirmed contamination level): "
        "C1_CD16_mono 29.9%, CD14_mono 6.0%, CD16_mono 29.7%, CD83_CD14_mono 5.4%, Mono_prolif "
        "9.4%. Corrected 2026-10-16 (per codex_output/mono_20261016_130032.md R2 -- the prior "
        "wording called C1_CD16_mono/CD16_mono 'coarser/less-specific' and stated they 'carry' "
        "higher contamination, neither of which this check supports): this is a difference in "
        "the historical model's ESTIMATED score only, not a validated biological or technical "
        "property of these subtypes, and is not grounds to exclude CD16_mono/C1_CD16_mono, "
        "declare their signal unreliable, or prefer CD83_CD14_mono/CD14_mono for that reason.",
        "Gate 2 Scrublet for Stephenson exists only as 'gate2_scrublet/stephenson."
        "INVALID_postqc_prior_20260921' -- explicitly marked INVALID by its own directory name, "
        "not a current, usable doublet result.",
        "'raw' in stephenson_monocyte_raw.h5ad names its PROVENANCE (layers['raw'] from the "
        "authors' own published atlas, per build_monocyte_canonical.py:187), not an absence of "
        "any QC -- the authors' own atlas may already reflect their own cell-calling/filtering; "
        "this was not independently checked.",
    ],
    "gap": (
        "No decision or artifact exists that says Gate 1 flags and/or DecontX-corrected counts "
        "ARE (or are not) the admitted input for THIS monocyte pilot -- that decision has not been "
        "made, by anyone, for this specific use. The decontaminated counts matrix's actual VALUES "
        "were not read this round (only the per-cell contamination SCORE, which is a small separate "
        "file) -- whether decontaminated counts should replace raw counts for a monocyte DE/scFM "
        "pilot, and whether Gate 1's flag_would_delete cells should be excluded, are open scientific "
        "choices, not resolved by this join existing. The broader platelet-level QC rebuild overall "
        "completion status (docs/current-status.md) is a separate, larger question this round does "
        "not re-litigate."
    ),
}


def _sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _load_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _load_json(path):
    with open(path) as f:
        return json.load(f)


def _donor_category_by_donor(coverage_summary_path=COVERAGE_SUMMARY_JSON):
    """Re-reads the ALREADY-COMPUTED symptomatic/asymptomatic split from the
    coverage summary's own per-donor rows -- does not recompute the
    Asymptomatic-vs-symptomatic classification here."""
    summary = _load_json(coverage_summary_path)
    by_donor = {}
    for row in summary.get("symptomatic_donor_rows", []):
        by_donor[row["true_donor_group"]] = "symptomatic"
    for row in summary.get("asymptomatic_donor_rows", []):
        by_donor[row["true_donor_group"]] = "asymptomatic"
    return by_donor, summary


def build_design_readiness_table(roster_path=ROSTER_TSV, grouping_path=GROUPING_TSV,
                                  counts_path=COUNTS_TSV, coverage_summary_path=COVERAGE_SUMMARY_JSON):
    roster_rows = _load_tsv(roster_path)
    grouping_by_sample = {r["sample_id"]: r for r in _load_tsv(grouping_path)}
    counts_rows = _load_tsv(counts_path)
    counts_by_sample = {r["sample_id"]: r for r in counts_rows}
    subtype_fields = [c for c in counts_rows[0].keys() if c not in ("sample_id", "n_cells_this_sample_total")]

    donor_category, coverage_summary_raw = _donor_category_by_donor(coverage_summary_path)

    by_donor = {}
    for r in roster_rows:
        by_donor.setdefault(r["true_donor_group"], []).append(r)
    multi_visit_donors = {d: rs for d, rs in by_donor.items() if len(rs) > 1}

    table_rows = []
    missing_from_counts = []
    for donor, rows in sorted(multi_visit_donors.items()):
        category = donor_category.get(donor, "unknown (not classified by coverage summary)")
        for r in rows:
            sid = r["sample_id"]
            grp = grouping_by_sample.get(sid, {})
            cnt = counts_by_sample.get(sid)
            row_out = {
                "true_donor_group": donor,
                "symptomatic_or_asymptomatic": category,
                "sample_id": sid,
                "visit_resample_flag": r["visit_resample_flag"],
                "collection_day": r["collection_day"],
                "own_severity": r["own_severity_Status_on_day_collection_summary"],
                "care_setting_at_sampling": grp.get("care_setting_at_sampling", "unknown (sample missing from grouping table)"),
                "status_on_day_collection_raw_icu_evidence": grp.get("status_on_day_collection_raw", "unknown (sample missing from grouping table)"),
                "n_monocyte_cells_this_sample_total": cnt["n_cells_this_sample_total"] if cnt else "unknown",
            }
            for st in subtype_fields:
                row_out[f"n_{st}_cells_raw"] = cnt[st] if cnt else "unknown"
            if cnt is None:
                missing_from_counts.append(sid)
            table_rows.append(row_out)

    fieldnames = list(table_rows[0].keys()) if table_rows else []

    cutoff_sensitivity = {}
    for st in subtype_fields:
        cutoff_sensitivity[st] = {}
        for cutoff in CUTOFFS:
            qualifying = {"symptomatic": [], "asymptomatic": [], "unknown": []}
            for donor, rows in sorted(multi_visit_donors.items()):
                cat = donor_category.get(donor, "unknown")
                per_visit_counts = []
                any_missing = False
                for r in rows:
                    cnt = counts_by_sample.get(r["sample_id"])
                    if cnt is None:
                        any_missing = True
                        break
                    per_visit_counts.append(int(cnt[st]))
                if any_missing:
                    continue
                if all(c >= cutoff for c in per_visit_counts):
                    bucket = qualifying.get(cat, qualifying["unknown"])
                    bucket.append({"true_donor_group": donor, "n_visits": len(rows), "per_visit_counts": per_visit_counts})
            cutoff_sensitivity[st][str(cutoff)] = qualifying

    summary = {
        "source_files": {
            "roster_tsv": roster_path,
            "roster_tsv_sha256": _sha256_file(roster_path),
            "grouping_tsv": grouping_path,
            "grouping_tsv_sha256": _sha256_file(grouping_path),
            "counts_tsv": counts_path,
            "counts_tsv_sha256": _sha256_file(counts_path),
            "coverage_summary_json": coverage_summary_path,
            "coverage_summary_json_sha256": _sha256_file(coverage_summary_path),
        },
        "n_multi_visit_donors": len(multi_visit_donors),
        "n_donor_visit_rows": len(table_rows),
        "n_sample_ids_missing_from_counts_table": len(missing_from_counts),
        "sample_ids_missing_from_counts_table": missing_from_counts,
        "subtype_fields": subtype_fields,
        "cutoffs_checked": list(CUTOFFS),
        "qc_source_evidence": QC_SOURCE_EVIDENCE,
        "provenance_caveat": (
            "All counts here are the SAME already-committed, already-verified raw counts from "
            "stephenson_sample_subtype_cell_counts_v1.tsv (source: the per-cohort RAW, provisional "
            "stephenson_monocyte_raw.h5ad) -- no new h5ad read was performed to build this table. "
            "The per-subtype cutoff sensitivity numbers are a DESCRIPTIVE exploration of how many "
            "donors would qualify at illustrative thresholds, NOT a new inclusion/QC rule, NOT a "
            "subtype choice, and NOT a claim that the surviving donors are ready for a real DE run."
        ),
    }
    return table_rows, fieldnames, cutoff_sensitivity, summary


def _require_output_does_not_already_exist(*paths):
    existing = [p for p in paths if os.path.exists(p)]
    if existing:
        raise SystemExit(f"refusing to overwrite an existing design-readiness output: {existing}")


def main(roster_path=ROSTER_TSV, grouping_path=GROUPING_TSV, counts_path=COUNTS_TSV,
         coverage_summary_path=COVERAGE_SUMMARY_JSON,
         out_tsv=OUT_DONOR_VISIT_SUBTYPE_TSV, out_json=OUT_CUTOFF_SENSITIVITY_JSON):
    _require_output_does_not_already_exist(out_tsv, out_json)

    table_rows, fieldnames, cutoff_sensitivity, summary = build_design_readiness_table(
        roster_path=roster_path, grouping_path=grouping_path, counts_path=counts_path,
        coverage_summary_path=coverage_summary_path,
    )

    with open(out_tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(table_rows)

    out = dict(summary)
    out["per_subtype_cutoff_sensitivity"] = cutoff_sensitivity
    out["output_files"] = {"donor_visit_subtype_tsv": out_tsv, "cutoff_sensitivity_json": out_json}
    with open(out_json, "w") as f:
        json.dump(out, f, indent=2, sort_keys=True)
        f.write("\n")

    print(json.dumps({k: v for k, v in out.items() if k != "per_subtype_cutoff_sensitivity"}, indent=2, sort_keys=True))
    return out


if __name__ == "__main__":
    main()
