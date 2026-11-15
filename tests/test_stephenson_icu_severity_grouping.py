import os
import sys

sys.path.insert(
    0,
    os.path.join(
        os.path.dirname(__file__), "..", "scripts", "34_scfm_discovery_framework", "results"
    ),
)

from stephenson_icu_severity_grouping import (
    CARE_SETTING_ICU,
    CARE_SETTING_NON_ICU,
    CARE_SETTING_UNKNOWN,
    GROUP_COMPARISON,
    GROUP_EXCLUDED_OTHER,
    GROUP_EXCLUDED_PENDING_ICU_EVIDENCE,
    GROUP_REFERENCE,
    candidate_main_group,
    care_setting_at_sampling,
)


def test_all_real_icu_status_values_map_to_icu():
    for raw in ("ITU_NIV", "ITU_O2", "ITU_intubated"):
        assert care_setting_at_sampling(raw) == CARE_SETTING_ICU


def test_all_real_non_icu_status_values_map_to_non_icu():
    for raw in ("Ward_NIV", "Ward_O2", "Ward_noO2", "Staff screening"):
        assert care_setting_at_sampling(raw) == CARE_SETTING_NON_ICU


def test_missing_or_unrecognized_status_maps_to_unknown_not_non_icu():
    for raw in (None, "", "nan", "NaN", "some_future_status_code"):
        assert care_setting_at_sampling(raw) == CARE_SETTING_UNKNOWN


def test_severe_critical_with_icu_is_the_comparison_group():
    for severity in ("Severe", "Critical"):
        group, reason = candidate_main_group(severity, CARE_SETTING_ICU)
        assert group == GROUP_COMPARISON
        assert "ICU" in reason


def test_severe_critical_without_confirmed_icu_is_excluded_not_reference():
    for severity in ("Severe", "Critical"):
        for setting in (CARE_SETTING_NON_ICU, CARE_SETTING_UNKNOWN):
            group, reason = candidate_main_group(severity, setting)
            assert group == GROUP_EXCLUDED_PENDING_ICU_EVIDENCE
            assert group != GROUP_REFERENCE


def test_mild_moderate_is_the_reference_group_regardless_of_location():
    for severity in ("Mild", "Moderate"):
        for setting in (CARE_SETTING_ICU, CARE_SETTING_NON_ICU, CARE_SETTING_UNKNOWN):
            group, _ = candidate_main_group(severity, setting)
            assert group == GROUP_REFERENCE


def test_asymptomatic_and_missing_severity_are_excluded_never_defaulted_into_reference():
    for severity in ("Asymptomatic", None, "", "Death"):
        group, reason = candidate_main_group(severity, CARE_SETTING_ICU)
        assert group == GROUP_EXCLUDED_OTHER
        assert group != GROUP_REFERENCE
        assert group != GROUP_COMPARISON


def test_location_is_never_rederived_from_severity_label_alone():
    """Mutation-kill target: a defect that inferred ICU from the severity grade
    itself (ignoring the real Status_on_day_collection string) would make a
    Severe/non-ICU sample wrongly land in the comparison group. This asserts
    the real, source-grounded behaviour instead."""
    group, reason = candidate_main_group("Severe", CARE_SETTING_NON_ICU)
    assert group != GROUP_COMPARISON
    group, reason = candidate_main_group("Critical", CARE_SETTING_UNKNOWN)
    assert group != GROUP_COMPARISON


def test_direction_is_comparison_minus_reference_not_reversed():
    """Guards against a flipped-sign defect: comparison and reference must
    remain distinct labels that a downstream contrast-direction check can
    tell apart, never collapsing to the same group."""
    comparison_group, _ = candidate_main_group("Severe", CARE_SETTING_ICU)
    reference_group, _ = candidate_main_group("Mild", CARE_SETTING_NON_ICU)
    assert comparison_group != reference_group
    assert comparison_group == GROUP_COMPARISON
    assert reference_group == GROUP_REFERENCE
