"""The Gate 2 frozen stress design: the grid, the seed, and who may be scored.

These tests exist to make the ruling mechanically enforceable rather than
documented. The recurring failure in this repository is a constraint that lives
only in prose, so each test below pins a specific way the design could be
loosened without anyone noticing.

Nothing here runs Scrublet, reads cohort data, or touches a result file. No
outcome value is used.
"""

import os
import sys

import pytest

CODE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "scripts",
    "32_rebuild_upfront_qc",
    "code",
)
sys.path.insert(0, CODE)

import gate2_stress_design as D  # noqa: E402

# ---- the grid ------------------------------------------------------------


def test_grid_is_exactly_the_ruled_values():
    assert D.RESIDUAL_RATE_GRID == (0.005, 0.01, 0.02, 0.04, 0.08, 0.12, 0.16, 0.20)


def test_grid_is_an_immutable_tuple():
    """A list could be appended to at import time by anything that imports it."""
    assert isinstance(D.RESIDUAL_RATE_GRID, tuple)
    with pytest.raises((TypeError, AttributeError)):
        D.RESIDUAL_RATE_GRID.append(0.5)


@pytest.mark.parametrize(
    "bad",
    [
        (0.005, 0.01, 0.02, 0.04, 0.08, 0.12, 0.16),  # dropped a point
        (0.005, 0.01, 0.02, 0.04, 0.08, 0.12, 0.16, 0.20, 0.3),  # added one
        (0.01, 0.005, 0.02, 0.04, 0.08, 0.12, 0.16, 0.20),  # reordered
        (0.005, 0.01, 0.02, 0.05, 0.08, 0.12, 0.16, 0.20),  # substituted
    ],
)
def test_a_modified_grid_is_refused(bad):
    """Elementwise, so reordering and substitution are caught as well as
    length changes -- a set or length comparison would pass three of these."""
    with pytest.raises(D.Gate2DesignError, match="frozen grid"):
        D.assert_grid_is_frozen(bad)


def test_the_frozen_grid_itself_passes():
    assert D.assert_grid_is_frozen(D.RESIDUAL_RATE_GRID) == D.RESIDUAL_RATE_GRID
    assert D.assert_grid_is_frozen(list(D.RESIDUAL_RATE_GRID)) == D.RESIDUAL_RATE_GRID


@pytest.mark.parametrize("rate", list(D.RESIDUAL_RATE_GRID))
def test_every_grid_value_is_accepted_as_a_runtime_rate(rate):
    assert D.assert_rate_on_grid(rate) == float(rate)


@pytest.mark.parametrize("bad", [None, "0.02", [0.02], float("nan"), float("inf")])
def test_a_non_numeric_or_non_finite_rate_raises_the_design_error(bad):
    """All design violations produce Gate2DesignError, so a runner can fail
    closed on one exception type rather than on a mix of TypeError/ValueError."""
    with pytest.raises(D.Gate2DesignError):
        D.assert_rate_on_grid(bad)


@pytest.mark.parametrize("flag", [True, False])
def test_a_bool_rate_is_rejected_AS_A_BOOL(flag):
    """bool is a subclass of int, so True reaches the grid check and is refused
    for being off-grid -- which means a test that only asserts the exception
    TYPE cannot tell whether the bool guard exists. Asserting the message makes
    the guard load-bearing: if a grid value were ever 1.0, True would otherwise
    arrive as a valid rate."""
    with pytest.raises(D.Gate2DesignError, match="bool"):
        D.assert_rate_on_grid(flag)


@pytest.mark.parametrize("off", [0.0, 0.05, 0.1, 0.006, 0.19999, 0.25, 1.0])
def test_a_rate_off_the_grid_is_refused(off):
    """Exact membership, not a tolerance: 0.19999 came from somewhere other
    than the grid, and the point of the design is that nothing else supplies a
    rate."""
    with pytest.raises(D.Gate2DesignError, match="not on the frozen grid"):
        D.assert_rate_on_grid(off)


# ---- per-cohort permitted use -------------------------------------------


@pytest.mark.parametrize(
    "cohort,expected",
    [
        ("qiu", D.USE_STANDARD_SCORING),
        ("chongqing", D.USE_STANDARD_SCORING),
        ("cillo", D.USE_RESIDUAL_DIAGNOSTIC),
        ("liu", D.USE_RESIDUAL_DIAGNOSTIC),
        ("combat", D.USE_RESIDUAL_DIAGNOSTIC),
        ("ren", D.USE_RESIDUAL_SCORING_ONLY),
        ("stephenson", D.USE_NOT_RUN),
        ("wilk", D.USE_NOT_RUN),
    ],
)
def test_each_ruled_cohort_has_its_ruled_use(cohort, expected):
    assert D.use_class_for(cohort) == expected


def test_all_eight_cohorts_are_ruled_and_no_more():
    assert set(D.COHORT_USE_CLASS) == {
        "qiu",
        "chongqing",
        "cillo",
        "liu",
        "combat",
        "ren",
        "stephenson",
        "wilk",
    }


@pytest.mark.parametrize("cohort", ["stephenson", "wilk"])
def test_a_not_run_cohort_cannot_be_scored(cohort):
    with pytest.raises(D.Gate2DesignError, match="not_run"):
        D.assert_may_score(cohort)


@pytest.mark.parametrize("cohort", ["qiu", "chongqing", "cillo", "liu", "combat", "ren"])
def test_every_other_cohort_may_be_scored(cohort):
    assert D.assert_may_score(cohort) in D.PERMITTED_USES


def test_an_unruled_cohort_raises_rather_than_defaulting():
    """Absence is not permission. A new cohort arriving without a ruling is the
    case that must stop the run."""
    with pytest.raises(D.Gate2DesignError, match="no recorded use class"):
        D.use_class_for("schulte_schrepping")


def test_cohort_lookup_is_case_and_whitespace_insensitive():
    assert D.use_class_for("  QIU ") == D.USE_STANDARD_SCORING


def test_every_use_class_has_a_stated_meaning():
    """A use class with no meaning is a label, and a label does not travel with
    the number into a figure caption."""
    assert set(D.USE_CLASS_MEANING) == D.PERMITTED_USES
    for use, text in D.USE_CLASS_MEANING.items():
        assert len(text) > 40, use


def test_ren_is_never_pooled():
    assert "ren" in D.NEVER_POOLED
    assert D.use_class_for("ren") == D.USE_RESIDUAL_SCORING_ONLY


# ---- seeds ---------------------------------------------------------------


def test_capture_seed_is_deterministic():
    assert D.capture_seed("qiu", "ch1") == D.capture_seed("qiu", "ch1")


def test_one_capture_uses_ONE_seed_across_the_whole_grid():
    """Common random numbers. The first version of this suite asserted the
    opposite -- eight distinct seeds for one capture's eight grid points -- which
    would have made every call-set difference across the grid carry both the
    calibration effect and Monte Carlo variation, with no way to separate them.

    The rate is not an input to the seed, so this is now structural: there is no
    parameter through which it could vary."""
    import inspect

    params = list(inspect.signature(D.capture_seed).parameters)
    assert params == ["cohort", "channel_id"], f"rate must not be a seed input: {params}"

    seed = D.capture_seed("qiu", "ch1")
    for r in D.RESIDUAL_RATE_GRID:
        D.assert_rate_on_grid(r)  # the rate is validated elsewhere
        assert D.capture_seed("qiu", "ch1") == seed


def test_different_captures_get_different_seeds():
    seeds = {
        D.capture_seed(c, ch)
        for c in ("qiu", "chongqing", "combat")
        for ch in ("ch1", "ch2", "pool7")
    }
    assert len(seeds) == 9, "two captures shared a stream"


def test_capture_seed_is_stable_across_reruns_and_case():
    a = D.capture_seed("qiu", "ch1")
    assert D.capture_seed("  QIU  ", "ch1") == a


@pytest.mark.parametrize(
    "cohort,channel",
    [
        ("", "ch1"),
        ("qiu", ""),
        ("  ", "ch1"),
        ("qiu", "   "),
        (None, "ch1"),
        ("qiu", None),
    ],
)
def test_capture_seed_requires_both_identifiers(cohort, channel):
    """None used to slip through: str(None) == 'None' is truthy, and the first
    version only tested empty strings."""
    with pytest.raises(D.Gate2DesignError, match="required"):
        D.capture_seed(cohort, channel)


def test_capture_seed_is_in_numpy_range():
    for ch in ("ch1", "pool7", "B2_10xlane01"):
        s = D.capture_seed("combat", ch)
        assert 0 <= s < 2**32 - 1


# ---- the rulings cannot be edited at runtime -----------------------------


@pytest.mark.parametrize("mapping_name", ["COHORT_USE_CLASS", "USE_CLASS_MEANING"])
def test_the_rulings_are_read_only_at_runtime(mapping_name):
    """A plain dict could be reassigned after import, and assert_may_score()
    would then let a not_run cohort through -- invisible to any source-level
    mutation test."""
    m = getattr(D, mapping_name)
    with pytest.raises(TypeError):
        m["stephenson"] = D.USE_STANDARD_SCORING
    with pytest.raises(TypeError):
        del m["qiu"]


def test_a_runtime_edit_attempt_leaves_the_ruling_intact():
    try:
        D.COHORT_USE_CLASS["stephenson"] = D.USE_STANDARD_SCORING
    except TypeError:
        pass
    assert D.use_class_for("stephenson") == D.USE_NOT_RUN
    with pytest.raises(D.Gate2DesignError):
        D.assert_may_score("stephenson")


def test_the_manifest_is_a_plain_serialisable_dict():
    import json as _json

    m = D.design_manifest()
    assert isinstance(m["cohort_use_class"], dict)
    _json.dumps(m)


# ---- fail-closed rules ---------------------------------------------------


def test_an_all_false_prediction_is_not_itself_a_failure_condition():
    """It can be legitimate. The threshold is what must be read -- scanpy warns
    rather than raising when auto-thresholding fails, then returns finite scores
    and an all-False prediction, which passes every field-level check."""
    joined = " ".join(D.UNKNOWN_CONDITIONS)
    assert "all_false" not in joined
    assert "threshold_absent" in D.UNKNOWN_CONDITIONS
    assert "threshold_non_finite" in D.UNKNOWN_CONDITIONS


def test_unknown_is_never_filled_with_false():
    assert D.UNKNOWN_FILL_VALUE is None
    assert D.UNKNOWN_FILL_VALUE is not False


def test_modality_is_a_reliability_diagnostic_not_a_mechanical_failure():
    """A first version listed `simulation_unimodal` among the mechanical UNKNOWN
    conditions while defining no detector, statistic, threshold or
    implementation -- which would have let a runner decide case by case when to
    reject a result, bypassing the method_unreliable layer whose cutoff has not
    been adjudicated."""
    assert "simulation_unimodal" not in D.UNKNOWN_CONDITIONS
    joined = " ".join(D.RELIABILITY_DIAGNOSTICS_REQUIRED)
    assert "histogram" in joined and "seed_stability" in joined
    assert "none adjudicated" in D.design_manifest()["reliability_cutoffs"]


def test_reliability_diagnostics_cannot_fail_a_capture():
    overlap = set(D.RELIABILITY_DIAGNOSTICS_REQUIRED) & set(D.UNKNOWN_CONDITIONS)
    assert not overlap, f"a reliability diagnostic is also a failure condition: {overlap}"


def test_the_failure_conditions_cover_the_known_silent_modes():
    for cond in (
        "threshold_absent",
        "threshold_non_finite",
        "scrublet_raised",
        "calls_not_strict_boolean",
        "calls_contain_na",
        "scores_not_finite",
        "scores_wrong_shape",
    ):
        assert cond in D.UNKNOWN_CONDITIONS


# ---- what the design claims ---------------------------------------------


def test_the_design_authorises_scoring_only():
    m = D.design_manifest()
    assert "scoring only" in m["authorises"]
    assert "separate ruling" in m["authorises"]


def test_the_grid_is_labelled_as_a_stress_axis_not_a_prevalence_estimate():
    m = D.design_manifest()
    assert "NOT a prevalence estimate" in m["grid_semantics"]
    assert "stress" in m["grid_semantics"]


def test_the_manifest_records_everything_a_reader_would_have_to_trust():
    m = D.design_manifest()
    for k in (
        "design_version",
        "decision_reference",
        "residual_rate_grid",
        "grid_semantics",
        "design_seed",
        "cohort_use_class",
        "never_pooled",
        "unknown_conditions",
        "reliability_diagnostics_required",
        "reliability_cutoffs",
        "seeding",
        "authorises",
        "frozen_before_results",
    ):
        assert k in m, k
    assert m["residual_rate_grid"] == list(D.RESIDUAL_RATE_GRID)
    assert m["cohort_use_class"]["stephenson"] == D.USE_NOT_RUN


def test_the_pre_result_freeze_list_names_the_no_selection_rule():
    joined = " ".join(D.freeze_before_results()).lower()
    assert "before any result is inspected" in joined
    assert "no grid point is chosen because" in joined
    assert "method is unstable" in joined


def test_the_runtime_axis_is_one_dimensional():
    """The `base rate x removal efficiency` design was rejected as
    unidentifiable. Asserted structurally rather than by scanning the docstring:
    a first version of this test searched for a phrase and failed only because
    the sentence wrapped across two lines, which is what a prose scan is worth."""
    for r in D.RESIDUAL_RATE_GRID:
        assert isinstance(r, float), "a grid point is a scalar rate, not a (rate, e) pair"
    assert isinstance(D.design_manifest()["residual_rate_grid"][0], float)


def test_no_removal_efficiency_parameter_is_exposed_at_runtime():
    """`e` may appear in the reporting layer; it must not be a runtime axis, or
    the researcher degrees of freedom this design removes come straight back."""
    # Token-level, not substring: a first version banned 'e_grid', which matches
    # the legitimate RESIDUAL_RAT-E_GRID and failed on clean code -- while also
    # appearing in every mutation result and inflating the kill counts.
    import re

    public = [n for n in dir(D) if not n.startswith("_")]
    banned = [n for n in public if re.search(r"(^|_)(removal|efficiency|eff)(_|$)", n.lower())]
    assert not banned, f"runtime removal-efficiency axis exposed: {banned}"
    # and the manifest declares exactly one runtime axis
    axes = [k for k in D.design_manifest() if k.endswith("_grid")]
    assert axes == ["residual_rate_grid"], f"more than one runtime axis: {axes}"


def test_the_grid_is_the_only_source_of_a_runtime_rate():
    """Nothing converts a vendor recovered-stage rate into a runtime rate. If
    such a helper existed, the unidentifiable mapping would be one call away."""
    public = [n for n in dir(D) if not n.startswith("_") and callable(getattr(D, n))]
    converters = [
        n
        for n in public
        if any(
            k in n.lower()
            for k in ("vendor", "convert", "from_recovered", "to_residual", "estimate_rate")
        )
    ]
    assert not converters, f"a vendor->runtime rate conversion is exposed: {converters}"
    with pytest.raises(D.Gate2DesignError):
        D.assert_rate_on_grid(0.076)  # a plausible vendor-derived number
