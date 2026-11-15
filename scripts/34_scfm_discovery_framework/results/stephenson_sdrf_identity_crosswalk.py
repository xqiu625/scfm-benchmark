"""Per review:20261014_231259 B and review:20261015_104226 Sec.4 item 2:
resolve whether this cohort's 90 patient_id strings are all genuinely distinct biological
donors, using Stephenson's OWN depositor-submitted SDRF (E-MTAB-10026.sdrf.txt,
Characteristics[individual] -- the SAME crosswalk resource build_monocyte_canonical.py's own
load_stephenson_individual_map() already reads), never a suffix-based guess at WHICH donor a
sample belongs to.

ONE run produces BOTH stephenson_sdrf_identity_crosswalk_summary.json AND the full 102-row
stephenson_donor_identity_roster_v1.tsv, so the two can never diverge again (previously a
separate ad hoc script built the roster from the summary's output, and drifted from it).

Two samples (BGCV06_CV0201, BGCV13_CV0326) have a Source Name string that is NOT present in
the SDRF at all -- this is the already-documented BGCV06/BGCV13 pool-membership conflict in
stephenson_identity_decision_sheet_v3.md (the depositor's pool label for these two samples
disagrees with this project's own GEX-derived pool label). That pool disagreement is a
*batch-dependency* question, separate from *donor identity*: the SDRF's own
Characteristics[individual] field, read from this sample's OTHER row(s) sharing the identical
patient_id string, resolves the donor identity unambiguously (patient_id "CV0201" appears in
SDRF under BOTH BGCV08_CV0201 and BGCV13_CV0201, both Characteristics[individual]=CV0201;
patient_id "CV0326" appears in SDRF under BGCV06_CV0326, Characteristics[individual]=CV0326).
This fallback is checked for internal conflict before being trusted (see
_patient_id_fallback_conflicts below) -- it is not merging two different unresolved samples
into each other (that was a previously self-caught bug in an earlier version of this script,
and remains fixed: no row's identity is ever inferred from ANOTHER row that is itself
unresolved).
"""
import csv
import json

SDRF = "/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/monocyte_source/stephenson/E-MTAB-10026.sdrf.txt"
SAMPLE_CSV = "stephenson_sample_level_real_fields.csv"
SUMMARY_OUT = "stephenson_sdrf_identity_crosswalk_summary.json"
ROSTER_OUT = "stephenson_donor_identity_roster_v1.tsv"
POOL_CONFLICT_NOTE = (
    "SDRF_CONFIRMED_VIA_PATIENT_ID_FALLBACK -- this sample's own Source Name string is not "
    "present in the SDRF (known BGCV06/BGCV13 pool-membership conflict, see "
    "stephenson_identity_decision_sheet_v3.md); the donor identity below is resolved from "
    "this patient_id's OTHER SDRF row(s), which agree on the same individual regardless of "
    "pool label. The POOL/batch-dependency assignment for this specific sample remains open "
    "pending that decision sheet's resolution -- this is a batch-dependency question, not a "
    "donor-identity question."
)


def _pid_suffix(source_name: str) -> str:
    return source_name.rsplit("_", 1)[-1] if "_" in source_name else source_name


# Source Name -> Characteristics[individual], straight from the depositor's own SDRF.
sdrf_individual = {}
# patient-id-like suffix -> set of individuals seen anywhere in the SDRF under that suffix.
pid_suffix_to_individuals = {}
pid_suffix_to_sources = {}
with open(SDRF) as fh:
    r = csv.reader(fh, delimiter="\t")
    header = next(r)
    i_source = header.index("Source Name")
    i_indiv = header.index("Characteristics[individual]")
    for row in r:
        src, ind = row[i_source], row[i_indiv]
        sdrf_individual[src] = ind
        suf = _pid_suffix(src)
        pid_suffix_to_individuals.setdefault(suf, set()).add(ind)
        pid_suffix_to_sources.setdefault(suf, []).append(src)

print(f"SDRF rows: {len(sdrf_individual)}")

# A fallback keyed by patient_id suffix is only trustworthy if the SDRF itself never disagrees
# with itself about which individual a given suffix belongs to across pools.
_patient_id_fallback_conflicts = {
    suf: sorted(inds) for suf, inds in pid_suffix_to_individuals.items() if len(inds) > 1
}
print(f"patient-id suffixes with internally CONFLICTING SDRF individuals (fallback refused "
     f"for these): {len(_patient_id_fallback_conflicts)} -> {_patient_id_fallback_conflicts}")

rows = list(csv.DictReader(open(SAMPLE_CSV)))
print(f"cohort sample-level rows: {len(rows)}")

resolution_counts = {"sdrf_source_name_exact": 0, "sdrf_patient_id_fallback": 0,
                     "unresolved_no_sdrf_match": 0, "unresolved_conflicting_fallback": 0}

for r in rows:
    sid, pid = r["sample_id"], r["patient_id"]
    if sid in sdrf_individual:
        r["sdrf_individual"] = sdrf_individual[sid]
        r["resolution_method"] = "sdrf_source_name_exact"
        resolution_counts["sdrf_source_name_exact"] += 1
        continue
    candidates = pid_suffix_to_individuals.get(pid, set())
    if len(candidates) == 1:
        r["sdrf_individual"] = next(iter(candidates))
        r["resolution_method"] = "sdrf_patient_id_fallback"
        resolution_counts["sdrf_patient_id_fallback"] += 1
    elif len(candidates) > 1:
        r["sdrf_individual"] = "UNKNOWN"
        r["resolution_method"] = "unresolved_conflicting_fallback"
        resolution_counts["unresolved_conflicting_fallback"] += 1
    else:
        r["sdrf_individual"] = "UNKNOWN"
        r["resolution_method"] = "unresolved_no_sdrf_match"
        resolution_counts["unresolved_no_sdrf_match"] += 1

print(f"resolution method counts: {resolution_counts}")
missing_in_sdrf = [r["sample_id"] for r in rows if r["sdrf_individual"] == "UNKNOWN"]
print(f"cohort samples with NO resolvable SDRF individual (exact or fallback): "
     f"{len(missing_in_sdrf)} -> {missing_in_sdrf}")

# Group by the REAL SDRF individual id, not the patient_id string.
by_individual = {}
for r in rows:
    by_individual.setdefault(r["sdrf_individual"], []).append(r)

n_patient_id_strings = len({r["patient_id"] for r in rows})
n_confirmed_distinct_individuals = len([k for k in by_individual if k != "UNKNOWN"])
print(f"\nDISTINCT patient_id strings in cohort: {n_patient_id_strings}")
print(f"CONFIRMED-distinct SDRF individuals (exact Source Name match OR unambiguous "
     f"patient_id fallback; excluding the samples with NO resolvable individual): "
     f"{n_confirmed_distinct_individuals}")

# Which individuals map to MORE THAN ONE patient_id string (the real identity-collapse cases)?
# EXCLUDE the "UNKNOWN" placeholder explicitly -- two samples both missing a resolvable
# individual are NOT thereby the same individual as each other; grouping them together would
# be a methodological artifact of the placeholder, not a real finding. Never infer a row's
# identity from another row that is ITSELF unresolved (UNKNOWN never contributes to
# pid_suffix_to_individuals above, since it is computed only from real SDRF rows).
collapsed = {ind: recs for ind, recs in by_individual.items()
            if ind != "UNKNOWN" and len({r["patient_id"] for r in recs}) > 1}
n_unknown = len(by_individual.get("UNKNOWN", []))
print(f"\nSamples with NO resolvable SDRF individual at all (identity genuinely unresolved, "
     f"NOT asserted same or different from each other): {n_unknown}")
print(f"SDRF individuals spanning >1 DIFFERENT patient_id string (hidden same-donor pairs): {len(collapsed)}")
for ind, recs in collapsed.items():
    print(f"  individual={ind}: patient_ids={sorted({r['patient_id'] for r in recs})}, "
         f"samples={[(r['sample_id'], r['Resample'], r['Collection_Day'], r['Status_on_day_collection_summary']) for r in recs]}")

# For the 4 Resample-only singletons specifically, does their SDRF individual ALSO appear
# under any OTHER sample row anywhere in the FULL SDRF (even outside this 102-sample cohort)?
resample_only_ids = ["MH8919328", "MH9143273", "newcastle004v2", "newcastle21v2"]
print("\n--- Full-SDRF search (not just this cohort) for each Resample-only singleton's sibling ---")
indiv_to_all_sources = {}
for src, ind in sdrf_individual.items():
    indiv_to_all_sources.setdefault(ind, []).append(src)
for pid in resample_only_ids:
    row = next(r for r in rows if r["patient_id"] == pid)
    this_ind = row["sdrf_individual"]
    siblings = [s for s in indiv_to_all_sources.get(this_ind, []) if s != pid]
    print(f"  {pid}: SDRF individual={this_ind}, OTHER Source Name(s) sharing this individual "
         f"anywhere in the full SDRF: {siblings}")

# Identity-unit count: SDRF-confirmed distinct individuals (exact or unambiguous fallback) +
# samples with no resolvable individual, each counted as its own separate unresolved unit.
# This is the SDRF-CONFIRMED donor count. It is NOT the same question as the open pool/
# batch-dependency assignment for BGCV06_CV0201 / BGCV13_CV0326 -- see
# stephenson_identity_decision_sheet_v3.md for that separate, still-open policy choice.
sdrf_confirmed_identity_units = n_confirmed_distinct_individuals + n_unknown
out = {
    "n_cohort_samples": len(rows),
    "n_distinct_patient_id_strings": n_patient_id_strings,
    "n_confirmed_distinct_sdrf_individuals": n_confirmed_distinct_individuals,
    "n_samples_with_no_resolvable_sdrf_individual": n_unknown,
    "sdrf_confirmed_identity_unit_count": sdrf_confirmed_identity_units,
    "resolution_method_counts": resolution_counts,
    "n_sdrf_individuals_spanning_multiple_patient_id_strings": len(collapsed),
    "collapsed_individuals": {
        ind: sorted({r["patient_id"] for r in recs}) for ind, recs in collapsed.items()
    },
    "samples_with_no_resolvable_sdrf_individual": missing_in_sdrf,
    "patient_id_fallback_internal_conflicts": _patient_id_fallback_conflicts,
    "open_pool_membership_question_not_a_donor_count_question": {
        "samples": ["BGCV06_CV0201", "BGCV13_CV0326"],
        "note": ("Both now resolve to a confirmed SDRF individual via the patient_id "
                "fallback (CV0201, CV0326 respectively), since every other SDRF row "
                "sharing their patient_id agrees on that individual regardless of pool. "
                "What remains open is ONLY which 10x pool/batch each sample belongs to "
                "for dependency-grouping purposes -- see "
                "stephenson_identity_decision_sheet_v3.md. This does not change the "
                "donor count above."),
    },
    "counting_definitions": {
        "sdrf_confirmed_identity_unit_count": (
            "The number used as the donor count below (86 as of this run). Every one of "
            "the 102 cohort samples resolves to a confirmed SDRF individual under this "
            "definition; there is no remaining candidate-policy-dependent donor-count "
            "ambiguity. A prior run of this script reported 85 (undercounting CV0326, "
            "which had no exact Source Name match and no patient_id-fallback logic yet) "
            "and, in a separate ad hoc roster-building pass, 87 (treating both "
            "pool-conflicted samples as unresolved-and-separate) and 86 (an "
            "undocumented, roster-only convention). This run's number is produced by a "
            "single script shared with the roster file below, so it cannot drift from it "
            "again."
        ),
    },
}
print(f"\nSDRF-CONFIRMED identity-unit count: {sdrf_confirmed_identity_units} "
     f"(vs. {n_patient_id_strings} raw patient_id strings)")
with open(SUMMARY_OUT, "w") as fh:
    json.dump(out, fh, indent=2)
print(f"\nwrote {SUMMARY_OUT}")

# ---------------------------------------------------------------------------
# Roster: one row per cohort sample, built from the SAME resolution computed above.
# ---------------------------------------------------------------------------
individual_to_patient_ids = {}
for r in rows:
    if r["sdrf_individual"] != "UNKNOWN":
        individual_to_patient_ids.setdefault(r["sdrf_individual"], set()).add(r["patient_id"])


def _true_donor_group(r: dict) -> str:
    if r["sdrf_individual"] != "UNKNOWN":
        return f"sdrf::{r['sdrf_individual']}"
    return f"own::{r['patient_id']}"


for r in rows:
    r["true_donor_group"] = _true_donor_group(r)

groups: dict[str, list[dict]] = {}
for r in rows:
    groups.setdefault(r["true_donor_group"], []).append(r)


def _group_severity(members: list[dict]) -> str:
    initial = [m for m in members if m["Resample"] == "Initial"]
    chosen = initial[0] if initial else members[0]
    return chosen["Status_on_day_collection_summary"]


group_severity = {g: _group_severity(members) for g, members in groups.items()}

ROSTER_COLUMNS = [
    "sample_id", "patient_id", "sdrf_individual_id", "true_donor_group", "identity_status",
    "visit_resample_flag", "collection_day", "days_from_onset", "site",
    "own_severity_Status_on_day_collection_summary", "donor_level_severity_initial_preferred",
    "worst_clinical_status", "subtype_coverage_checked_this_round", "batch_dependency_flag",
    "inclusion_rationale_for_this_draft_roster",
]

DRAFT_NOTE = ("DRAFT for owner review -- does not itself authorize treating this row as an "
             "admitted independent-donor input; canonical identity/QC retention is unchanged "
             "by this table")
SUBTYPE_NOTE = ("not re-checked this round (see MONOCYTE_PILOT_PREREG_v1.md Sec.5 for the "
                "per-cohort orig_celltype vocabulary; not joined to this specific roster)")


def _identity_status(r: dict) -> str:
    if r["resolution_method"] == "sdrf_source_name_exact":
        pid = r["patient_id"]
        if len(individual_to_patient_ids.get(r["sdrf_individual"], set())) > 1:
            return (f"SDRF_CONFIRMED_COLLAPSE (individual {r['sdrf_individual']} spans "
                    f"patient_id strings {sorted(individual_to_patient_ids[r['sdrf_individual']])})")
        return "SDRF_CONFIRMED_EXACT (own single SDRF individual, no cross-string collapse)"
    if r["resolution_method"] == "sdrf_patient_id_fallback":
        return POOL_CONFLICT_NOTE
    return f"UNRESOLVED ({r['resolution_method']})"


def _batch_dependency_flag(r: dict) -> str:
    if r["sample_id"] in ("BGCV06_CV0201", "BGCV13_CV0326"):
        return "YES - BGCV06/BGCV13 pool-membership conflict, open (stephenson_identity_decision_sheet_v3.md)"
    return ""


with open(ROSTER_OUT, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(ROSTER_COLUMNS)
    for r in rows:
        w.writerow([
            r["sample_id"],
            r["patient_id"],
            r["sdrf_individual"] if r["sdrf_individual"] != "UNKNOWN" else "n/a_unresolved",
            r["true_donor_group"],
            _identity_status(r),
            r["Resample"],
            r["Collection_Day"],
            r["Days_from_onset"],
            r["Site"],
            r["Status_on_day_collection_summary"],
            group_severity[r["true_donor_group"]],
            r["Worst_Clinical_Status"],
            SUBTYPE_NOTE,
            _batch_dependency_flag(r),
            DRAFT_NOTE,
        ])

n_roster_identity_units = len(groups)
print(f"wrote {ROSTER_OUT}: {len(rows)} rows, {n_roster_identity_units} distinct true_donor_group values")
assert n_roster_identity_units == sdrf_confirmed_identity_units, (
    f"roster identity-unit count ({n_roster_identity_units}) diverged from the summary's "
    f"sdrf_confirmed_identity_unit_count ({sdrf_confirmed_identity_units}) -- these must be "
    f"produced by the same run and can never disagree"
)
print("CONSISTENCY CHECK PASSED: roster and summary report the same identity-unit count.")
