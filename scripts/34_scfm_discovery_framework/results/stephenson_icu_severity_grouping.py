"""Pure mapping functions for the owner-approved monocyte-pilot comparison.

Per Cindy's decisions recorded in `docs/decisions.md` (2026-10-15, via
`review:20261015_120031` and `review:20261015_120317`):

    comparison_group = severity in {Severe, Critical} AND confirmed ICU location
                        at the SAME sample's collection time
    reference_group  = severity in {Mild, Moderate}
    direction         = comparison_group - reference_group

Two rules this module exists to enforce, because both were explicitly
prohibited by the owner's decision text:

1. ICU location must come from the real, source-recorded clinical-status
   string (Stephenson's own `Status_on_day_collection` vocabulary), never
   inferred from the severity grade itself. A `Severe`/`Critical` label is
   not treated as proof of ICU location.
2. A sample whose severity is Severe/Critical but whose location is
   non-ICU or unknown is EXCLUDED from both groups -- it must never be
   silently folded into the Mild/Moderate reference group, and never
   silently relabeled "non-ICU" when the source value is actually missing
   or in conflict.
"""
from __future__ import annotations

ICU_STATUS_VALUES = frozenset({"ITU_NIV", "ITU_O2", "ITU_intubated"})
NON_ICU_STATUS_VALUES = frozenset({"Ward_NIV", "Ward_O2", "Ward_noO2", "Staff screening"})

SEVERE_CRITICAL_SUMMARY_VALUES = frozenset({"Severe", "Critical"})
MILD_MODERATE_SUMMARY_VALUES = frozenset({"Mild", "Moderate"})

CARE_SETTING_ICU = "ICU"
CARE_SETTING_NON_ICU = "non-ICU"
CARE_SETTING_UNKNOWN = "unknown"

GROUP_COMPARISON = "comparison_severe_critical_icu_at_sampling"
GROUP_REFERENCE = "reference_mild_moderate"
GROUP_EXCLUDED_PENDING_ICU_EVIDENCE = "excluded_severe_critical_without_confirmed_icu"
GROUP_EXCLUDED_OTHER = "excluded_other_severity"


def care_setting_at_sampling(raw_status_on_day_collection: str | None) -> str:
    """Map Stephenson's own `Status_on_day_collection` string to ICU/non-ICU/unknown.

    This is a literal, value-by-value lookup against the source vocabulary --
    it never looks at the severity summary field. An unrecognized or missing
    value returns "unknown"; it is never guessed as non-ICU.
    """
    if raw_status_on_day_collection is None:
        return CARE_SETTING_UNKNOWN
    value = raw_status_on_day_collection.strip()
    if not value or value.lower() == "nan":
        return CARE_SETTING_UNKNOWN
    if value in ICU_STATUS_VALUES:
        return CARE_SETTING_ICU
    if value in NON_ICU_STATUS_VALUES:
        return CARE_SETTING_NON_ICU
    return CARE_SETTING_UNKNOWN


def candidate_main_group(severity_summary: str | None, care_setting: str) -> tuple[str, str]:
    """Return (group_label, reason) for the owner-approved comparison.

    `severity_summary` is the sampling-day severity category
    (`Status_on_day_collection_summary`: Asymptomatic/Mild/Moderate/Severe/Critical).
    `care_setting` must already be the output of `care_setting_at_sampling`
    (ICU/non-ICU/unknown) -- this function does not re-derive it from severity.
    """
    summary = (severity_summary or "").strip()
    if summary in SEVERE_CRITICAL_SUMMARY_VALUES:
        if care_setting == CARE_SETTING_ICU:
            return GROUP_COMPARISON, "severity Severe/Critical, confirmed ICU at sampling"
        if care_setting == CARE_SETTING_NON_ICU:
            return (
                GROUP_EXCLUDED_PENDING_ICU_EVIDENCE,
                "severity Severe/Critical but confirmed NON-ICU at sampling -- "
                "excluded from the comparison arm; not moved to the reference arm",
            )
        return (
            GROUP_EXCLUDED_PENDING_ICU_EVIDENCE,
            "severity Severe/Critical but ICU location at sampling is unknown/unrecorded -- "
            "excluded pending evidence; not assumed ICU, not assumed non-ICU, "
            "not moved to the reference arm",
        )
    if summary in MILD_MODERATE_SUMMARY_VALUES:
        return GROUP_REFERENCE, "severity Mild/Moderate"
    return (
        GROUP_EXCLUDED_OTHER,
        f"severity_summary={summary!r} is neither Severe/Critical nor Mild/Moderate "
        "(e.g. Asymptomatic, missing) -- never defaulted into the reference group",
    )
