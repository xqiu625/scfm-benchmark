"""Applies Cindy's one-off cross-group resolution for the single mixed-visit
donor `stephenson_icu_severity_grouping_v1_summary.json` flagged, per her
direct reply in the automated review conversation (`review:20261015_142903`
Sec.1, quoting: "如果只能进入一组, 就进入目前分组少的那组" -- "if it can only
enter one group, put it in whichever group currently has fewer members").

This is explicitly a ONE-TIME, documented, hash-pinned decision about ONE
named donor (`sdrf::C-8821`) under the EXACT, frozen candidate roster this
decision was made against -- NOT a general tie-break rule, and NOT a script
that re-derives "the right answer" from whatever input it is handed.

Per `review:20261015_150231` F2: an earlier version of this script
only checked that the INPUT's own self-reported `donors_with_mixed_group_
across_visits` field equalled `[RESOLVED_DONOR]` -- it never verified the
input was the SPECIFIC, approved roster this decision was actually made
against. A probe showed that changing an UNRELATED donor's
`candidate_main_group` (changing the overall counts, e.g. reference 53->51)
while leaving the stale `donors_with_mixed_group_across_visits` field
untouched still passed every check and silently overwrote the output,
including overwriting a pre-existing output file (a sentinel placed at the
output path was clobbered without complaint).

Fixed here: (1) the input TSV/summary are checked against FROZEN, hardcoded
sha256 hashes of the specific, approved roster this decision was made
against -- never a hash recomputed from "whatever was just read" and
declared to be the expected value; a single changed byte anywhere in either
input file, related or not, refuses the whole run. (2) The two named sample
rows are independently verified (exactly one row each, correct donor, and
the KEPT sample is actually in the comparison group / the EXCLUDED sample is
actually in the reference group) before anything is computed or written.
(3) ALL conflict/count computation happens before any file is opened for
writing. (4) The script refuses outright if its own output path already
exists -- it never overwrites a prior run's result; a caller who wants to
regenerate must pass a new, not-yet-existing output path explicitly.
"""
import csv
import hashlib
import json
import os

HERE = os.path.dirname(__file__)
IN_TSV = os.path.join(HERE, "stephenson_icu_severity_grouping_v1.tsv")
IN_SUMMARY = os.path.join(HERE, "stephenson_icu_severity_grouping_v1_summary.json")
DEFAULT_OUT_TSV = os.path.join(HERE, "stephenson_icu_severity_grouping_resolved_v1.tsv")
DEFAULT_OUT_SUMMARY = os.path.join(HERE, "stephenson_icu_severity_grouping_resolved_v1_summary.json")

#: The EXACT, approved roster this decision was made against. A caller
#: supplying a different `stephenson_icu_severity_grouping_v1.tsv`/
#: `_v1_summary.json` -- for ANY reason, including an edit unrelated to the
#: one resolved donor -- is refused outright. This is deliberately NOT
#: regenerated from whatever is on disk; it is the one, frozen, reviewed
#: value, re-verified by `review:20261015_150231`'s own
#: independent re-hash before being hardcoded here.
FROZEN_IN_TSV_SHA256 = "1bd3e72e503d50e24205488f3a0d9eadc80bf797908bc8fee7ac82c813b85090"
FROZEN_IN_SUMMARY_SHA256 = "e758edf6a6edcc7d3bd6ea37f49fc5bc540c67cf7755422bc10495fd49e16ee2"

#: The ONE donor this decision resolves, and nothing else.
RESOLVED_DONOR = "sdrf::C-8821"
GROUP_COMPARISON = "comparison_severe_critical_icu_at_sampling"
GROUP_REFERENCE = "reference_mild_moderate"
#: Per Cindy's own citation of the specific candidate sample
#: (review:20261015_142903 Sec.1): D12, Critical, ITU_intubated --
#: the sample kept in the comparison (ICU severe/critical) arm.
KEPT_SAMPLE_ID = "newcastle004v2"
#: D0, Moderate, Ward_O2 -- removed from the reference (mild/moderate) arm
#: by this decision; the underlying sample/row is NOT deleted from the
#: roster, only excluded from this round's main-comparison candidate set.
EXCLUDED_SAMPLE_ID = "MH9179821"
DECISION_CITATION = (
    "docs/decisions.md 2026-10-15 (via codex_output/20261015_142903.md Sec.1): "
    "Cindy's direct reply 'ru guo zhi neng jin ru yi zu, jiu jin ru mu qian "
    "fen zu shao de na zu' (if it can only enter one group, assign it to "
    "whichever group currently has FEWER members) -- a one-time decision "
    "under the CURRENT candidate roster, not the authors' original analysis "
    "rule, and not derived from any expression/pathway result.")


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def _require_frozen_input():
    actual_tsv = _sha256(IN_TSV)
    actual_summary = _sha256(IN_SUMMARY)
    problems = []
    if actual_tsv != FROZEN_IN_TSV_SHA256:
        problems.append(
            f"{IN_TSV} has sha256 {actual_tsv!r}, expected the FROZEN "
            f"{FROZEN_IN_TSV_SHA256!r} -- this script resolves ONE specific, "
            "approved roster, never whatever the file currently on disk happens "
            "to contain")
    if actual_summary != FROZEN_IN_SUMMARY_SHA256:
        problems.append(
            f"{IN_SUMMARY} has sha256 {actual_summary!r}, expected the FROZEN "
            f"{FROZEN_IN_SUMMARY_SHA256!r}")
    if problems:
        raise SystemExit(
            "refusing to apply this one-off resolution -- the input does not match the "
            f"exact, approved roster this decision was made against: {problems}"
        )


def _require_resolved_rows(rows_by_sample):
    """Verifies the two named samples independently -- never trusts the
    hardcoded sample_id constants alone without checking what they actually
    are in THIS input.
    """
    for sample_id in (KEPT_SAMPLE_ID, EXCLUDED_SAMPLE_ID):
        matches = [r for r in rows_by_sample if r["sample_id"] == sample_id]
        if len(matches) != 1:
            raise SystemExit(
                f"expected exactly 1 row for sample_id={sample_id!r}, found {len(matches)}")
        (row,) = matches
        if row["true_donor_group"] != RESOLVED_DONOR:
            raise SystemExit(
                f"sample_id={sample_id!r} belongs to true_donor_group="
                f"{row['true_donor_group']!r}, expected {RESOLVED_DONOR!r}")

    kept = next(r for r in rows_by_sample if r["sample_id"] == KEPT_SAMPLE_ID)
    excluded = next(r for r in rows_by_sample if r["sample_id"] == EXCLUDED_SAMPLE_ID)
    if kept["candidate_main_group"] != GROUP_COMPARISON:
        raise SystemExit(
            f"{KEPT_SAMPLE_ID!r} (the sample this decision KEEPS) has "
            f"candidate_main_group={kept['candidate_main_group']!r}, expected "
            f"{GROUP_COMPARISON!r} -- refusing to keep a sample that is not actually "
            "in the comparison arm")
    if excluded["candidate_main_group"] != GROUP_REFERENCE:
        raise SystemExit(
            f"{EXCLUDED_SAMPLE_ID!r} (the sample this decision EXCLUDES) has "
            f"candidate_main_group={excluded['candidate_main_group']!r}, expected "
            f"{GROUP_REFERENCE!r} -- refusing to exclude a sample that is not actually "
            "in the reference arm")


def _require_output_does_not_already_exist(out_tsv, out_summary):
    existing = [p for p in (out_tsv, out_summary) if os.path.exists(p)]
    if existing:
        raise SystemExit(
            f"refusing to overwrite existing output file(s): {existing} -- this script "
            "never overwrites a prior run's result; pass a new, not-yet-existing output "
            "path to regenerate"
        )


def resolve(out_tsv=DEFAULT_OUT_TSV, out_summary=DEFAULT_OUT_SUMMARY):
    _require_frozen_input()
    _require_output_does_not_already_exist(out_tsv, out_summary)

    with open(IN_SUMMARY) as f:
        summary_in = json.load(f)
    # Defense in depth: also recompute the mixed-donor list directly from the
    # TSV (never only trust the summary file's own self-reported field),
    # though the frozen-hash check above already pins both inputs exactly.
    mixed_donor_ids_claimed = sorted(
        d["true_donor_group"] for d in summary_in["donors_with_mixed_group_across_visits"])
    if mixed_donor_ids_claimed != [RESOLVED_DONOR]:
        raise SystemExit(
            f"the frozen summary's own mixed-donor list is {mixed_donor_ids_claimed!r}, "
            f"not exactly [{RESOLVED_DONOR!r}] -- refusing (this should be unreachable "
            "once the frozen-hash check above has passed; kept as a second, independent "
            "guard)"
        )

    with open(IN_TSV, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    _require_resolved_rows(rows)

    # ALL conflict/count computation happens here, BEFORE any file is opened
    # for writing -- a refusal from this point on never leaves a partially
    # written output.
    resolved_rows = []
    for row in rows:
        r = dict(row)
        if r["true_donor_group"] == RESOLVED_DONOR:
            if r["sample_id"] == KEPT_SAMPLE_ID:
                r["resolved_main_group"] = r["candidate_main_group"]
                r["resolution_note"] = (
                    f"kept -- fewer-populated-group decision, see {DECISION_CITATION}")
            else:
                r["resolved_main_group"] = "excluded_by_cross_group_resolution"
                r["resolution_note"] = (
                    "excluded from this round's main-comparison candidate set by the "
                    f"fewer-populated-group decision (raw data/QC/identity unaffected): "
                    f"{DECISION_CITATION}")
        else:
            r["resolved_main_group"] = r["candidate_main_group"]
            r["resolution_note"] = ""
        resolved_rows.append(r)

    donor_resolved_group = {}
    for r in resolved_rows:
        g = r["resolved_main_group"]
        if g == "excluded_by_cross_group_resolution":
            continue
        donor_resolved_group.setdefault(r["true_donor_group"], set()).add(g)

    conflicts = {d: sorted(gs) for d, gs in donor_resolved_group.items() if len(gs) > 1}
    if conflicts:
        raise SystemExit(f"resolution still leaves donor(s) in multiple groups: {conflicts}")

    resolved_counts = {}
    for gs in donor_resolved_group.values():
        (g,) = gs
        resolved_counts[g] = resolved_counts.get(g, 0) + 1

    out_summary_payload = {
        "decision_citation": DECISION_CITATION,
        "resolved_donor": RESOLVED_DONOR,
        "kept_sample_id": KEPT_SAMPLE_ID,
        "excluded_sample_id": EXCLUDED_SAMPLE_ID,
        "frozen_input_hashes_this_resolution_was_applied_against": {
            os.path.basename(IN_TSV): FROZEN_IN_TSV_SHA256,
            os.path.basename(IN_SUMMARY): FROZEN_IN_SUMMARY_SHA256,
        },
        "donor_level_group_counts_before_resolution_any_visit": summary_in[
            "donor_level_group_counts_any_visit"],
        "donor_level_group_counts_after_resolution_one_group_per_donor": resolved_counts,
        "note": (
            "This resolves ONLY the one named cross-group donor above, as a one-time, "
            "hash-pinned decision against a FROZEN, exact input -- it is NOT a general "
            "tie-break rule, and this script refuses to run against any input that is "
            "not byte-identical to the approved roster. A future round with additional "
            "mixed-group donors needs a fresh, explicit counting-basis decision. "
            "Repeated-visit selection WITHIN an already-resolved group (e.g. a donor "
            "with two Mild visits) remains a separate, still-open question -- see "
            "MONOCYTE_PILOT_PREREG_v1.md Sec.6."
        ),
    }

    # Only now, after every check above has passed, is anything written.
    fieldnames = list(resolved_rows[0].keys())
    with open(out_tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(resolved_rows)
    with open(out_summary, "w") as f:
        json.dump(out_summary_payload, f, indent=2, sort_keys=True)

    return out_summary_payload


def main():
    payload = resolve()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
