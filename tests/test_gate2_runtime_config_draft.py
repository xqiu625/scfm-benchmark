"""Parameter approval must not bypass capture/config/pilot execution gates."""

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

import gate2_runtime_config_draft as C  # noqa: E402
import gate2_stress_design as D  # noqa: E402
import gate2_stress_runner as R  # noqa: E402

# ---- the production runner still refuses ---------------------------------


def test_the_production_runner_still_refuses_to_run():
    """The whole point of the draft: it does not unlock the runner."""
    assert R.PRODUCTION_PARAMETER_RULING == {}
    with pytest.raises(R.Gate2RunnerError, match="production configuration is not activated"):
        R.validate_parameter_ruling(None)


def test_emitting_a_runnable_config_is_refused_and_names_the_gap():
    with pytest.raises(R.Gate2RunnerError, match="per-cohort eligibility artifact"):
        C.emit_runnable_config()


def test_no_unit_is_runtime_eligible_in_the_draft():
    cfg = C.config_draft()
    assert cfg["is_runnable_config"] is False
    assert cfg["any_unit_runtime_eligible"] is False
    for name, c in cfg["cohorts"].items():
        assert c["runtime_eligible"] is False, name
        assert c["runtime_eligible_blockers"], name


def test_all_parameters_are_approved_but_production_is_not_activated():
    cfg = C.config_draft()
    assert cfg["unresolved_parameters"] == []
    assert cfg["n_unresolved_parameters"] == 0
    assert set(cfg["ruled_parameters"]) == set(R.REQUIRED_PARAMETER_NAMES)
    assert cfg["production_configuration_activated"] is False
    for cohort in cfg["cohorts"].values():
        assert not any("parameter(s) unruled" in b for b in cohort["runtime_eligible_blockers"])


# ---- the ruling's own constraints ----------------------------------------


def test_the_grid_is_the_full_eight_points_and_is_frozen():
    cfg = C.config_draft()
    assert cfg["grid"] == [0.005, 0.01, 0.02, 0.04, 0.08, 0.12, 0.16, 0.20]
    assert cfg["grid"] == list(D.RESIDUAL_RATE_GRID)


def test_the_only_permitted_label_is_the_stress_setting():
    assert C.STRESS_SETTING_LABEL == "assumed_algorithmic_stress_setting"
    assert C.assert_label_permitted(C.STRESS_SETTING_LABEL)
    for bad in C.FORBIDDEN_LABELS:
        with pytest.raises(D.Gate2DesignError, match="refused"):
            C.assert_label_permitted(bad)


def test_the_removal_efficiency_axis_cannot_become_a_runtime_axis():
    """`e`, `r` and `s` stay in the interpretation layer: several
    unidentifiable combinations give the same residual."""
    assert C.assert_no_runtime_axis_beyond_the_prior(["expected_doublet_rate"])
    for axis in C.INTERPRETATION_ONLY_AXES:
        with pytest.raises(D.Gate2DesignError, match="interpretation-layer"):
            C.assert_no_runtime_axis_beyond_the_prior(["expected_doublet_rate", axis])


def test_only_the_prior_is_a_runtime_axis_in_the_draft():
    cfg = C.config_draft()
    assert cfg["runtime_axes"] == ["expected_doublet_rate"]
    for axis in C.INTERPRETATION_ONLY_AXES:
        assert axis not in cfg["runtime_axes"]


def test_rank_and_auc_invariance_are_labelled_implementation_checks():
    inv = C.config_draft()["invariance_checks"]
    assert inv["rank_invariance_across_grid"] == "implementation_check_only"
    assert inv["observed_vs_simulated_auc_invariance"] == "implementation_check_only"
    assert "NOT diagnostic accuracy" in inv["auc_interpretation"]


def test_the_per_capture_fixed_set_covers_the_neighbour_graph():
    fixed = C.config_draft()["fixed_per_capture_seed"]
    for item in (
        "simulated_parent_pairs",
        "pca_representation",
        "neighbour_graph",
        "hvg_selection",
        "input_cell_order",
    ):
        assert item in fixed


# ---- cohort classes match the frozen design and the ruling ---------------


def test_cohort_classes_agree_with_the_frozen_design():
    for name, r in C.COHORT_RULING.items():
        assert r["use"] == D.COHORT_USE_CLASS[name], name


def test_not_run_cohorts_get_no_grid_and_no_label():
    cfg = C.config_draft()
    for name in ("stephenson", "wilk"):
        assert cfg["cohorts"][name]["can_design_scoring_experiment"] == "no"
        assert cfg["cohorts"][name]["grid"] == []
        assert cfg["cohorts"][name]["grid_label"] is None


def test_the_three_upstream_filtered_cohorts_carry_the_residual_null_label():
    cfg = C.config_draft()
    for name in ("cillo", "liu", "combat"):
        assert "residual_null_unvalidated" in cfg["cohorts"][name]["evidence_labels"], name
        assert cfg["cohorts"][name]["use_class"] == D.USE_RESIDUAL_DIAGNOSTIC


def test_ren_is_residual_scoring_only_and_never_pooled():
    cfg = C.config_draft()["cohorts"]["ren"]
    assert cfg["use_class"] == D.USE_RESIDUAL_SCORING_ONLY
    assert cfg["never_pooled"] is True


def test_every_conditional_cohort_lists_stop_conditions():
    cfg = C.config_draft()
    for name, c in cfg["cohorts"].items():
        assert c["stop_conditions"], name
        assert c["evidence_labels"], name


def test_the_common_conditions_are_all_present():
    cond = C.config_draft()["common_conditions"]
    assert sorted(cond) == [
        "C1_input",
        "C2_unit",
        "C3_support",
        "C4_fixed_implementation",
        "C5_interpretation",
    ]
    assert "never False, 0 or singlet" in cond["C4_fixed_implementation"]


def test_the_draft_is_json_serialisable():
    import json

    json.dumps(C.config_draft(), allow_nan=False)


# ---- the two errors the 2026-09-24 scoped audit caught --------------------


def test_the_draft_never_cites_the_stage_31_preflight_as_a_gate_2_blocker():
    """Stage 31's 302 Sample_IDs are a donor-level platelet analysis, not the
    input evidence for these eight cohorts. An earlier draft listed
    "Stage 31 preflight, 302/302 unknown" as a blocker on all of them, which let
    one work line's input audit silently gate seven unrelated cohorts."""
    cfg = C.config_draft()
    for name, c in cfg["cohorts"].items():
        for b in c["runtime_eligible_blockers"]:
            assert "302/302" not in b, f"{name}: {b}"
            assert "Stage 31" not in b, f"{name}: {b}"
    import json

    assert "302/302" not in json.dumps(cfg)


def test_each_conditional_cohort_states_its_own_preconditions():
    cfg = C.config_draft()
    for name, c in cfg["cohorts"].items():
        if c["can_design_scoring_experiment"] != "conditional":
            continue
        conds = [b.split(":")[0] for b in c["runtime_eligible_blockers"]]
        assert "global" in conds, name
        for needed in ("C1", "C2", "C3"):
            assert needed in conds, f"{name} does not state {needed}"


def test_a_complete_ruling_is_refused_for_missing_per_cohort_eligibility():
    """The refusal must name the per-cohort input preconditions, not Stage 31."""
    complete = {name: "test-only" for name in R.REQUIRED_PARAMETER_NAMES}
    with pytest.raises(R.Gate2RunnerError) as exc:
        C.emit_runnable_config(complete)
    msg = str(exc.value)
    assert "per-cohort eligibility artifact" in msg
    assert "Stage 31" not in msg and "302" not in msg
