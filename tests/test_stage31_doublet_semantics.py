"""Stage 31 doublet-call semantics: unknown must never collapse to False/0%.

Covers the entry point and three of the four layers the P0b defect ran through
-- filtering, aggregation and statistics -- by testing the pure functions the
script now uses, plus a round-trip check of the exported column semantics.

The EXPORT layer is covered only at the column level: the test below writes and
re-reads the two columns and checks that a blank does not come back as False.
It does NOT exercise the script's own `to_csv` call, which needs a real AnnData
and is therefore left to the re-run audit.

Nothing here runs Scrublet, reads cohort data, or touches any result file. No
outcome/mortality value is used to design a rule or a threshold: `y` appears
only as an opaque column that must survive the row filtering.
"""

import os
import sys

import numpy as np
import pandas as pd
import pytest

CODE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "scripts",
    "31_donor_level_revision",
    "code",
)
sys.path.insert(0, CODE)


def _load_module():
    """Import the semantics without executing the script body.

    The script reads HPCC paths, parses argv and runs the analysis at import
    time, so it cannot simply be imported. Static classification of which
    top-level statements are "configuration" kept getting it wrong -- `X =
    a_clean.X` is uppercase, `_MANIFEST = {...}` looks like a literal but
    references runtime names, and `EXAM_UNKNOWN_STATUSES` is a tuple of Name
    references that literal_eval cannot evaluate. So: take imports, classes and
    functions, then try each top-level assignment in order and keep the ones
    that actually evaluate against what has been defined so far.
    """
    import ast
    import types

    src = open(os.path.join(CODE, "doublet_detection_and_reaudit.py")).read()
    tree = ast.parse(src)
    mod = types.ModuleType("stage31_semantics")
    mod.__dict__["__file__"] = os.path.join(CODE, "doublet_detection_and_reaudit.py")

    for node in tree.body:
        if isinstance(node, (ast.Import, ast.ImportFrom, ast.FunctionDef, ast.ClassDef)):
            block = ast.Module(body=[node], type_ignores=[])
        elif isinstance(node, ast.Assign):
            block = ast.Module(body=[node], type_ignores=[])
        else:
            continue  # if/try/expression statements are script body
        try:
            exec(compile(block, "<stage31>", "exec"), mod.__dict__)
        except BaseException:
            # BaseException, not Exception: `_args = _ap.parse_args()` is an
            # Assign whose evaluation calls sys.exit when argv lacks --out-dir,
            # and SystemExit does not inherit from Exception.
            pass
    return mod


s31 = _load_module()


def make_obs(spec):
    """spec: list of (patient_id, status, call) with call in {True, False, None}."""
    rows = []
    for i, (pid, status, call) in enumerate(spec):
        rows.append(
            {
                "patient_id": pid,
                "Sample_ID": f"{pid}_s",
                "y": 0,
                "doublet_score": np.nan if call is None else 0.5,
                "doublet_call": pd.NA if call is None else call,
                "examination_status": status,
            }
        )
    df = pd.DataFrame(rows)
    df["doublet_call"] = df["doublet_call"].astype("boolean")
    return df


# --------------------------------------------------------------------------
# aggregation layer
# --------------------------------------------------------------------------


def test_rate_denominator_is_examined_cells_only():
    obs = make_obs(
        [
            ("p1", s31.EXAM_EXAMINED, True),
            ("p1", s31.EXAM_EXAMINED, False),
            ("p1", s31.EXAM_TOO_SMALL, None),
            ("p1", s31.EXAM_TOO_SMALL, None),
        ]
    )
    r = s31.doublet_rate_per_patient(obs)
    # 1 doublet of 2 EXAMINED, not 1 of 4
    assert r.loc["p1", "doublet_rate"] == 0.5
    assert r.loc["p1", "n_examined"] == 2
    assert r.loc["p1", "n_unknown"] == 2
    assert r.loc["p1", "examined_coverage"] == 0.5


def test_patient_with_zero_examined_cells_is_na_not_zero():
    """The defect in one line: every cell unexamined used to give a measured 0%."""
    obs = make_obs([("p1", s31.EXAM_TOO_SMALL, None), ("p1", s31.EXAM_RUNTIME_FAILED, None)])
    r = s31.doublet_rate_per_patient(obs)
    assert pd.isna(r.loc["p1", "doublet_rate"]), "must be NA, never 0.0"
    assert r.loc["p1", "doublet_rate"] != 0.0
    assert r.loc["p1", "n_examined"] == 0


@pytest.mark.parametrize(
    "status", [s31.EXAM_TOO_SMALL, s31.EXAM_RUNTIME_FAILED, s31.EXAM_NOT_EVALUATED]
)
def test_every_unknown_status_is_excluded_from_the_denominator(status):
    obs = make_obs([("p1", s31.EXAM_EXAMINED, True)] + [("p1", status, None)] * 9)
    r = s31.doublet_rate_per_patient(obs)
    assert r.loc["p1", "doublet_rate"] == 1.0, f"{status} must not dilute the rate"
    assert r.loc["p1", "n_examined"] == 1


def test_mixed_patients_each_get_their_own_denominator():
    obs = make_obs(
        [
            ("fully_examined", s31.EXAM_EXAMINED, True),
            ("fully_examined", s31.EXAM_EXAMINED, False),
            ("fully_examined", s31.EXAM_EXAMINED, False),
            ("partly", s31.EXAM_EXAMINED, False),
            ("partly", s31.EXAM_TOO_SMALL, None),
            ("none_examined", s31.EXAM_RUNTIME_FAILED, None),
        ]
    )
    r = s31.doublet_rate_per_patient(obs)
    assert r.loc["fully_examined", "doublet_rate"] == pytest.approx(1 / 3)
    assert r.loc["partly", "doublet_rate"] == 0.0  # genuinely measured 0 of 1
    assert pd.isna(r.loc["none_examined", "doublet_rate"])
    assert r.loc["partly", "examined_coverage"] == 0.5


def test_a_real_zero_and_an_undefined_rate_are_distinguishable():
    """0.0 must mean 'examined, none were doublets'. NaN must mean 'not measured'."""
    obs = make_obs(
        [
            ("measured_zero", s31.EXAM_EXAMINED, False),
            ("measured_zero", s31.EXAM_EXAMINED, False),
            ("undefined", s31.EXAM_TOO_SMALL, None),
            ("undefined", s31.EXAM_TOO_SMALL, None),
        ]
    )
    r = s31.doublet_rate_per_patient(obs)
    assert r.loc["measured_zero", "doublet_rate"] == 0.0
    assert r.loc["measured_zero", "n_examined"] == 2
    assert pd.isna(r.loc["undefined", "doublet_rate"])
    assert r.loc["undefined", "n_examined"] == 0


# --------------------------------------------------------------------------
# summary layer
# --------------------------------------------------------------------------


def test_summary_never_counts_unknown_as_non_doublet():
    obs = make_obs(
        [
            ("p1", s31.EXAM_EXAMINED, True),
            ("p1", s31.EXAM_EXAMINED, False),
            ("p1", s31.EXAM_TOO_SMALL, None),
            ("p2", s31.EXAM_RUNTIME_FAILED, None),
        ]
    )
    c = s31.summarize_examination(obs)
    assert c["n_cells"] == 4
    assert c["n_examined"] == 2
    assert c["n_unknown"] == 2
    assert c["n_doublet"] == 1
    assert c["n_non_doublet"] == 1, "the 2 unknown cells must NOT be counted here"
    assert c["n_doublet"] + c["n_non_doublet"] == c["n_examined"]
    assert c["coverage"] == 0.5


def test_summary_on_an_all_unknown_cohort():
    obs = make_obs([("p1", s31.EXAM_TOO_SMALL, None)] * 5)
    c = s31.summarize_examination(obs)
    assert c["n_examined"] == 0 and c["n_unknown"] == 5
    assert c["n_doublet"] == 0 and c["n_non_doublet"] == 0
    assert c["coverage"] == 0.0


# --------------------------------------------------------------------------
# filtering layer
# --------------------------------------------------------------------------


def test_exclusion_removes_only_confirmed_doublets_and_keeps_unknown_visible():
    obs = make_obs(
        [
            ("p1", s31.EXAM_EXAMINED, True),
            ("p1", s31.EXAM_EXAMINED, False),
            ("p1", s31.EXAM_TOO_SMALL, None),
        ]
    )
    drop = (obs["doublet_call"] == True).to_numpy(dtype=bool, na_value=False)  # noqa: E712
    kept = obs[~drop]
    assert drop.sum() == 1
    assert len(kept) == 2
    # the unknown cell survives AND is still identifiable as unknown
    unknown = kept[kept["examination_status"] != s31.EXAM_EXAMINED]
    assert len(unknown) == 1
    assert unknown["doublet_call"].isna().all(), "must not have become False"


def test_na_is_not_dropped_by_the_equality_mask():
    """`obs['doublet_call'] == True` on a nullable column must treat NA as
    not-True without raising and without silently dropping the row."""
    obs = make_obs([("p1", s31.EXAM_TOO_SMALL, None)] * 3)
    drop = (obs["doublet_call"] == True).to_numpy(dtype=bool, na_value=False)  # noqa: E712
    assert drop.sum() == 0
    assert len(obs[~drop]) == 3


# --------------------------------------------------------------------------
# statistics layer
# --------------------------------------------------------------------------


def test_spearman_inclusion_is_explicit_and_counts_are_reported():
    obs = make_obs(
        [
            ("p1", s31.EXAM_EXAMINED, True),
            ("p2", s31.EXAM_EXAMINED, False),
            ("p3", s31.EXAM_TOO_SMALL, None),
        ]
    )
    r = s31.doublet_rate_per_patient(obs)
    valid = r["doublet_rate"].notna()
    assert int(valid.sum()) == 2, "p3 has no examined cells and must be excluded"
    assert not valid.loc["p3"]
    # and the excluded patient is still present in the exported table
    assert "p3" in r.index


def test_undefined_rates_are_neither_dropped_silently_nor_filled():
    obs = make_obs([("p1", s31.EXAM_EXAMINED, True), ("p2", s31.EXAM_TOO_SMALL, None)])
    r = s31.doublet_rate_per_patient(obs)
    assert len(r) == 2, "the undefined patient stays in the table for transparency"
    assert pd.isna(r.loc["p2", "doublet_rate"])
    assert not (
        r["doublet_rate"].fillna(0) == r["doublet_rate"]
    ).all(), "filling with 0 must change the values -- i.e. nothing pre-filled them"


# --------------------------------------------------------------------------
# the defect must not reappear
# --------------------------------------------------------------------------


def test_source_has_no_live_false_collapse():
    """Scan CODE only. Docstrings and comments quote the old expressions to
    explain why they were removed; matching those would forbid documenting the
    defect."""
    import ast

    src = open(os.path.join(CODE, "doublet_detection_and_reaudit.py")).read()
    tree = ast.parse(src)
    # drop every docstring, then unparse what remains
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.ClassDef)) and node.body:
            first = node.body[0]
            if (
                isinstance(first, ast.Expr)
                and isinstance(first.value, ast.Constant)
                and isinstance(first.value.value, str)
            ):
                node.body = node.body[1:] or [ast.Pass()]
    code = ast.unparse(tree)  # comments are not preserved by unparse
    for forbidden in (
        "fillna(False)",
        "predicted_doublet'] = False",
        'predicted_doublet"] = False',
    ):
        assert forbidden not in code, f"collapse reintroduced in code: {forbidden}"
    # and the nullable dtype must still be in use
    assert "dtype='boolean'" in code or 'dtype="boolean"' in code


def test_status_constants_are_distinct_and_exhaustive():
    all_statuses = {s31.EXAM_EXAMINED, *s31.EXAM_UNKNOWN_STATUSES}
    assert len(all_statuses) == len(s31.EXAM_UNKNOWN_STATUSES) + 1, "statuses must be distinct"
    assert s31.EXAM_EXAMINED not in s31.EXAM_UNKNOWN_STATUSES
    # threshold_failed was added after the audit: scrublet returning without
    # raising is not proof of a usable result
    assert s31.EXAM_THRESHOLD_FAILED in s31.EXAM_UNKNOWN_STATUSES


# ==========================================================================
# Round 2 — the second collapse path: `examined` + NA call
# ==========================================================================
#
# Closing `fillna(False)` removed the explicit collapse but left an implicit
# one. A cell marked `examined` whose call is NA enters the rate denominator
# while contributing nothing to the numerator -- diluting unknown into
# non-doublet exactly as before, just through a different door. Measured on the
# pre-fix code: one confirmed doublet among one examined cell plus three
# examined-but-NA cells reported 0.250 instead of 1.000.


def contradictory_obs(status, call, n=4):
    return make_obs([("p1", status, call)] * n)


def test_r2_examined_with_na_call_is_refused():
    obs = contradictory_obs(s31.EXAM_EXAMINED, None)
    with pytest.raises(s31.ExaminationStateError, match="NA doublet_call"):
        s31.validate_examination_state(obs)


def test_r2_the_dilution_path_is_closed_in_the_rate():
    """The exact shape that reported 0.250 instead of 1.000."""
    obs = make_obs(
        [
            ("p1", s31.EXAM_EXAMINED, True),
            ("p1", s31.EXAM_EXAMINED, None),
            ("p1", s31.EXAM_EXAMINED, None),
            ("p1", s31.EXAM_EXAMINED, None),
        ]
    )
    with pytest.raises(s31.ExaminationStateError):
        s31.doublet_rate_per_patient(obs)


def test_r2_the_dilution_path_is_closed_in_the_summary():
    obs = make_obs([("p1", s31.EXAM_EXAMINED, True), ("p1", s31.EXAM_EXAMINED, None)])
    with pytest.raises(s31.ExaminationStateError):
        s31.summarize_examination(obs)


@pytest.mark.parametrize("status", list(s31.EXAM_UNKNOWN_STATUSES))
@pytest.mark.parametrize("call", [True, False])
def test_r2_unknown_status_carrying_a_call_is_refused(status, call):
    obs = contradictory_obs(status, call)
    with pytest.raises(s31.ExaminationStateError, match="carry a doublet_call"):
        s31.validate_examination_state(obs)


def test_r2_unrecognized_status_is_refused():
    obs = contradictory_obs("scored_probably", None)
    with pytest.raises(s31.ExaminationStateError, match="unknown examination_status"):
        s31.validate_examination_state(obs)


def test_r2_legal_states_pass_validation():
    obs = make_obs(
        [
            ("p1", s31.EXAM_EXAMINED, True),
            ("p1", s31.EXAM_EXAMINED, False),
            ("p1", s31.EXAM_TOO_SMALL, None),
            ("p1", s31.EXAM_THRESHOLD_FAILED, None),
            ("p1", s31.EXAM_RUNTIME_FAILED, None),
            ("p1", s31.EXAM_NOT_EVALUATED, None),
        ]
    )
    assert s31.validate_examination_state(obs) is obs
    c = s31.summarize_examination(obs)
    assert c["n_doublet"] + c["n_non_doublet"] == c["n_examined"] == 2


def test_r2_summary_arithmetic_identity_is_enforced():
    obs = make_obs(
        [
            ("p1", s31.EXAM_EXAMINED, True),
            ("p1", s31.EXAM_EXAMINED, False),
            ("p1", s31.EXAM_TOO_SMALL, None),
        ]
    )
    c = s31.summarize_examination(obs)
    assert c["n_doublet"] + c["n_non_doublet"] == c["n_examined"]


# ---- strict boolean extraction at the entry point -------------------------


@pytest.mark.parametrize(
    "values,why",
    [
        (np.array([0.5, 0.5, 0.5]), "strict boolean"),
        (np.array([0, 1, 2]), "strict boolean"),
        (np.array(["True", "False", "True"]), "strict boolean"),
        (np.array([[True, False], [True, False], [True, False]]), "not 1-D"),
        (np.array([True, False]), "length"),
        (None, "not 1-D"),  # np.asarray(None) is 0-d, not unreadable
    ],
)
def test_r2_non_boolean_scrublet_output_is_refused(values, why):
    calls, reason = s31.strict_boolean_calls(values, 3)
    assert calls is None, f"{values!r} must not be accepted as calls"
    assert why in reason, reason


@pytest.mark.parametrize(
    "values",
    [
        np.array([True, False, True]),
        np.array([1, 0, 1], dtype=np.int8),
        np.array([1, 0, 1], dtype=np.int64),
    ],
)
def test_r2_real_boolean_output_is_accepted(values):
    calls, reason = s31.strict_boolean_calls(values, 3)
    assert reason is None and calls is not None
    assert list(calls) == [True, False, True]


def test_r2_an_na_anywhere_disqualifies_the_whole_sample():
    """A per-cell NA means the detector produced no usable result for that
    sample; keeping the rest would leave `examined` cells with no call."""
    calls, reason = s31.strict_boolean_calls(np.array([True, None, False], dtype=object), 3)
    assert calls is None and "NA" in reason


def test_r2_numpy_astype_bool_would_have_accepted_what_this_refuses():
    """Records which conversion is the permissive one.

    pandas' own `pd.array(..., dtype="boolean")` actually REJECTS 0.5 and 2.0 --
    an earlier version of this comment claimed otherwise and was wrong. The
    permissive path is numpy's `.astype(bool)`, which maps every non-zero value
    to True and would invent calls from whatever the column held.
    """
    with pytest.raises(TypeError):
        pd.array(np.array([0.5, 2.0, 0.0]), dtype="boolean")
    assert list(np.array([0.5, 2.0, 0.0]).astype(bool)) == [
        True,
        True,
        False,
    ], "numpy is the permissive one"
    strict, reason = s31.strict_boolean_calls(np.array([0.5, 2.0, 0.0]), 3)
    assert strict is None and "strict boolean" in reason


def test_r2_too_small_is_never_inferred_from_a_nan_score():
    """too_small must come from the cell-count rule only. A NaN score in a large
    sample means the detector failed, which is a different fact."""
    import ast

    src = open(os.path.join(CODE, "doublet_detection_and_reaudit.py")).read()
    tree = ast.parse(src)
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and "EXAM_TOO_SMALL" in ast.unparse(node.value):
            ctx = ast.unparse(node)
            assert (
                "isna" not in ctx and "doublet_score" not in ctx
            ), f"too_small inferred from a score: {ctx[:90]}"
    code = ast.unparse(tree)
    assert "too_small_set" in code, "the cell-count rule must drive too_small"


# ==========================================================================
# Round 3 — threshold metadata, fail-closed score extraction, strict dtype
# ==========================================================================

# ---- the stage-32 failure mode: finite scores, all-False calls, NaN threshold


def test_r3_nan_threshold_is_refused_even_when_every_field_looks_valid():
    """The exact shape that passed the previous version: finite scores, a
    strictly boolean all-False prediction, and a threshold that is NaN. Every
    field-level check succeeds; only reading the threshold catches it."""
    thr, why = s31._finite_threshold({"threshold": float("nan")})
    assert thr is None and "not finite" in why


@pytest.mark.parametrize(
    "uns,why",
    [
        (None, "no scrublet metadata"),
        ({}, "carries no threshold"),
        ({"threshold": None}, "not a number"),
        ({"threshold": "auto"}, "not a number"),
        ({"threshold": float("nan")}, "not finite"),
        ({"threshold": float("inf")}, "not finite"),
    ],
)
def test_r3_unusable_threshold_metadata_is_refused(uns, why):
    thr, reason = s31._finite_threshold(uns)
    assert thr is None, f"{uns!r} must not yield a usable threshold"
    assert why in reason, reason


def test_r3_a_finite_threshold_is_accepted():
    thr, reason = s31._finite_threshold({"threshold": 0.32})
    assert reason is None and thr == 0.32


def test_r3_all_false_predictions_are_not_themselves_a_failure():
    """An all-False prediction can be legitimate. The threshold is what decides,
    so a valid threshold with no doublets must still be accepted."""
    thr, reason = s31._finite_threshold({"threshold": 0.9})
    assert reason is None
    calls, why = s31.strict_boolean_calls(np.zeros(5, dtype=bool), 5)
    assert why is None and not any(calls)


# ---- score extraction must fail closed, not raise


def test_r3_missing_score_column_does_not_raise_from_numpy():
    """np.asarray(None, dtype=float) is 0-d, so `.shape[0]` raised IndexError --
    from outside the per-sample try/except, killing the whole run instead of
    marking one sample unknown."""
    with pytest.raises(IndexError):
        np.asarray(None, dtype=float).shape[0]
    # the fixed code never reaches that expression: it checks for None first
    import ast

    src = open(os.path.join(CODE, "doublet_detection_and_reaudit.py")).read()
    fn = next(
        n
        for n in ast.walk(ast.parse(src))
        if isinstance(n, ast.FunctionDef) and n.name == "evaluate_sample"
    )
    body = ast.unparse(fn)
    assert "raw_scores is None" in body, "a missing column must be handled before conversion"
    # `scores.shape[0]` is present and SAFE: the ndim check returns before it,
    # so a 0-d array can never reach the subscript. Banning the string was the
    # wrong property -- ORDER is what makes it safe, so assert the order.
    ndim_at = body.index("scores.ndim != 1")
    assert (
        "scores.shape[0]" not in body[:ndim_at]
    ), "the length check must come AFTER the dimensionality check"


def test_r3_zero_dimensional_score_is_refused_before_any_indexing():
    # Behavioural counterpart: a 0-d score returns a reason, never IndexError.
    status, calls, scores, why = s31.evaluate_sample(
        {"predicted_doublet": np.zeros(3, dtype=bool), "doublet_score": np.float64(0.5)},
        "s1",
        3,
        {"threshold": 0.3},
    )
    assert status == s31.EXAM_RUNTIME_FAILED
    assert scores is None
    assert "not 1-D" in why


# ---- the invariant must require a strict boolean dtype, not merely non-NA


@pytest.mark.parametrize("bad", ["False", "True", 0.0, 1.0, 2, 0.5])
def test_r3_non_boolean_call_column_is_refused(bad):
    """'False', 0.0 and 2 are all non-NA, so an NA-only check let them through."""
    obs = pd.DataFrame(
        {"patient_id": ["p1"], "doublet_call": [bad], "examination_status": [s31.EXAM_EXAMINED]}
    )
    with pytest.raises(s31.ExaminationStateError, match="nullable boolean"):
        s31.validate_examination_state(obs)


@pytest.mark.parametrize("bad", ["False", 0.0, 2])
def test_r3_non_boolean_call_on_an_unknown_status_is_also_refused(bad):
    obs = pd.DataFrame(
        {"patient_id": ["p1"], "doublet_call": [bad], "examination_status": [s31.EXAM_TOO_SMALL]}
    )
    with pytest.raises(s31.ExaminationStateError):
        s31.validate_examination_state(obs)


@pytest.mark.parametrize("fn_name", ["doublet_rate_per_patient", "summarize_examination"])
def test_r3_non_boolean_column_fails_closed_before_rate_or_summary(fn_name):
    obs = pd.DataFrame(
        {
            "patient_id": ["p1", "p1"],
            "doublet_call": ["True", "False"],
            "examination_status": [s31.EXAM_EXAMINED] * 2,
        }
    )
    with pytest.raises(s31.ExaminationStateError, match="nullable boolean"):
        getattr(s31, fn_name)(obs)


def test_r3_a_proper_nullable_boolean_column_still_passes():
    obs = make_obs([("p1", s31.EXAM_EXAMINED, True), ("p1", s31.EXAM_TOO_SMALL, None)])
    assert s31.validate_examination_state(obs) is obs


# ==========================================================================
# Round 4 — no partial writes, strict score dimensionality, export round-trip
# ==========================================================================


def _uns(thr=0.3):
    return {"threshold": thr}


def _obs_for(n, calls=None, scores=None):
    d = {}
    if calls is not None:
        d["predicted_doublet"] = calls
    if scores is not None:
        d["doublet_score"] = scores
    return d


# ---- evaluate_sample is pure: a failure carries no values -----------------


def test_r4_evaluate_sample_is_pure_and_returns_nothing_on_failure():
    """Every failure returns (status, None, None, why). A caller cannot keep a
    stale call beside a failure, because there is nothing to keep."""
    status, calls, scores, why = s31.evaluate_sample(
        _obs_for(3, calls=np.array([True, False, True]), scores=np.array([0.1, 0.2, 0.3])),
        "s1",
        3,
        _uns(float("nan")),
    )
    assert status == s31.EXAM_THRESHOLD_FAILED
    assert calls is None and scores is None
    assert "not finite" in why


def test_r4_evaluate_sample_success_returns_both_values():
    status, calls, scores, why = s31.evaluate_sample(
        _obs_for(3, calls=np.array([True, False, False]), scores=np.array([0.9, 0.1, 0.2])),
        "s1",
        3,
        _uns(),
    )
    assert status == s31.EXAM_EXAMINED and why is None
    assert list(calls) == [True, False, False]
    assert list(scores) == [0.9, 0.1, 0.2]


def test_r4_evaluate_sample_does_not_mutate_its_input():
    obs = _obs_for(3, calls=np.array([True, False, True]), scores=np.array([0.1, 0.2, 0.3]))
    before = {k: np.array(v, copy=True) for k, v in obs.items()}
    s31.evaluate_sample(obs, "s1", 3, _uns())
    for k, v in before.items():
        assert np.array_equal(np.asarray(obs[k]), v), f"{k} was mutated"


# ---- 2-D scores must be refused, not flattened ---------------------------


@pytest.mark.parametrize("shape", [(3, 1), (1, 3), (2, 2)])
def test_r4_two_dimensional_scores_are_refused_not_flattened(shape):
    """`.ravel()` accepted (n,1), (1,n) and any shape whose size matched. A 2-D
    return means the schema is wrong; flattening would guess an element order."""
    n = int(np.prod(shape))
    status, calls, scores, why = s31.evaluate_sample(
        _obs_for(n, calls=np.zeros(n, dtype=bool), scores=np.zeros(shape)), "s1", n, _uns()
    )
    assert status == s31.EXAM_RUNTIME_FAILED
    assert scores is None
    assert "not 1-D" in why


def test_r4_one_dimensional_scores_of_the_right_length_still_pass():
    status, _, scores, why = s31.evaluate_sample(
        _obs_for(3, calls=np.zeros(3, dtype=bool), scores=np.array([0.1, 0.2, 0.3])),
        "s1",
        3,
        _uns(),
    )
    assert status == s31.EXAM_EXAMINED and why is None and scores.ndim == 1


# ---- reset before fallback ------------------------------------------------


def test_r4_reset_clears_every_detection_field():
    obs = make_obs([("p1", s31.EXAM_EXAMINED, True), ("p1", s31.EXAM_EXAMINED, False)])
    obs["doublet_score"] = [0.9, 0.1]
    s31.reset_detection_state(obs)
    assert obs["doublet_call"].isna().all()
    assert obs["doublet_score"].isna().all()
    assert (obs["examination_status"] == s31.EXAM_NOT_EVALUATED).all()
    assert s31.validate_examination_state(obs) is obs


def test_r4_partial_batched_write_then_fallback_leaves_no_contradiction():
    """The scenario the audit named: sample A accepted by the batched loop, a
    later sample raises, the fallback runs, and A then fails there. Without the
    reset, A would keep its call beside runtime_failed."""
    obs = make_obs(
        [
            ("pA", s31.EXAM_EXAMINED, True),  # accepted by batched
            ("pB", s31.EXAM_NOT_EVALUATED, None),
        ]
    )
    obs["doublet_score"] = [0.9, np.nan]
    # fallback begins
    s31.reset_detection_state(obs)
    # and now A fails in the fallback -- status only, no values to conflict with
    obs.loc[obs["patient_id"] == "pA", "examination_status"] = s31.EXAM_RUNTIME_FAILED
    s31.validate_examination_state(obs)  # must not raise
    assert obs["doublet_call"].isna().all()


def test_r4_without_the_reset_the_contradiction_is_detectable():
    """Confirms the reset is load-bearing rather than decorative."""
    obs = make_obs([("pA", s31.EXAM_EXAMINED, True)])
    obs.loc[0, "examination_status"] = s31.EXAM_RUNTIME_FAILED  # no reset
    with pytest.raises(s31.ExaminationStateError, match="carry a doublet_call"):
        s31.validate_examination_state(obs)


def test_r4_commit_happens_after_the_whole_loop():
    """a.obs must not be touched until every sample has been evaluated."""
    import ast

    src = open(os.path.join(CODE, "doublet_detection_and_reaudit.py")).read()
    code = ast.unparse(ast.parse(src))
    assert "staged.append" in code, "results must be staged, not written per sample"
    assert "_commit(staged)" in code, "a single commit after the loop"
    assert "reset_detection_state(a.obs)" in code, "fallback must reset first"


# ---- export column semantics (column level only) --------------------------


def test_r4_exported_blank_does_not_read_back_as_false(tmp_path):
    """A nullable boolean writes <NA> as an empty field; re-reading must not
    turn it into False. Column level only -- the script's own to_csv needs a
    real AnnData and is left to the re-run audit."""
    obs = make_obs(
        [
            ("p1", s31.EXAM_EXAMINED, True),
            ("p1", s31.EXAM_EXAMINED, False),
            ("p1", s31.EXAM_TOO_SMALL, None),
        ]
    )
    path = tmp_path / "calls.csv"
    obs[["patient_id", "doublet_call", "examination_status"]].to_csv(path, index=False)
    back = pd.read_csv(path)
    assert back["doublet_call"].isna().sum() == 1, "the unknown cell must come back blank"
    assert list(back["doublet_call"].dropna().astype(str)) == ["True", "False"]
    # and the status column still explains why
    assert back.loc[back["doublet_call"].isna(), "examination_status"].iloc[0] == s31.EXAM_TOO_SMALL


# ==========================================================================
# Round 5 — rerun prerequisites: outcome access, manifest lifecycle, round trip
# ==========================================================================

import ast as _ast

SCRIPT = os.path.join(CODE, "doublet_detection_and_reaudit.py")


def _src():
    return open(SCRIPT).read()


def _fn(name):
    return next(
        n
        for n in _ast.walk(_ast.parse(_src()))
        if isinstance(n, _ast.FunctionDef) and n.name == name
    )


# ---- integration-only must not ACCESS outcome, not merely skip computing it
#
# These are SOURCE-level and weak by construction. An earlier version of this
# block asserted `drop(columns=["y"])` was present, which made "parse the whole
# file, then delete the column" the passing condition -- exactly the defect.
# The real proof of the boundary is the behavioural block (round 6) below;
# these remain only to pin the shape of the branch.


def test_r5_loader_excludes_y_at_parse_time_not_by_dropping():
    body = _ast.unparse(_fn("load_canonical_cohort"))
    assert (
        "usecols" in body
    ), "the patient list must be parsed with usecols; dropping y afterwards is not exclusion"
    assert (
        "drop(columns=['y'])" not in body and 'drop(columns=["y"])' not in body
    ), "a post-parse drop of y is the defect this replaced"
    i_guard = body.index("INTEGRATION_ONLY")
    i_ysum = body.index("patient_list['y'].sum()")
    assert i_guard < i_ysum, "the event-count assertion must sit inside the full-mode branch"


def test_r5_event_count_assertion_is_not_unconditional():
    """`assert ... y.sum() == 25` reads outcome. It belongs in full mode only."""
    body = _ast.unparse(_fn("load_canonical_cohort"))
    for line in body.splitlines():
        if "sum() == 25" in line:
            assert line.startswith(" "), "the event-count assert must be indented under a branch"


def test_r5_obs_never_gains_y_under_integration_only():
    body = _ast.unparse(_fn("load_canonical_cohort"))
    assert (
        "if y_by_patient is not None" in body
    ), "writing y into obs must be conditional on having read it"


def test_r5_downstream_y_is_built_inside_the_full_mode_branch():
    src = _src()
    # the per-patient y Series must not be constructed at module level
    for line in src.splitlines():
        if "y_by_patient[p] for p in per_patient.index" in line:
            assert line.startswith("    "), "y must be built inside the full-mode branch"


# ---- manifest lifecycle


def test_r5_manifest_starts_as_started_and_claims_nothing_executed():
    src = _src()
    assert '"run_status": "started"' in src
    assert (
        '"outcome_analyses_executed": []' in src
    ), "the initial manifest must claim no analysis has run yet"
    assert '"outcome_analysis_requested"' in src, "authorization is separate from what happened"
    # and the facts must start at "has not happened", never be predicted from
    # the mode -- `"outcome_data_accessed": not INTEGRATION_ONLY` claimed an
    # access for any full-mode job that died before the read.
    for field in ('"known_outcome_values_loaded": False', '"outcome_values_inspected": False'):
        assert field in src, 'each fact starts at "has not happened", never predicted from the flag'
    assert (
        '"outcome_data_accessed"' not in src
    ), "the ambiguous field name is retired: loaded and inspected are separate facts"
    assert 'not INTEGRATION_ONLY,\n        "outcome_values' not in src


def test_r5_manifest_marks_completed_only_at_the_end():
    src = _src()
    i_started = src.index('"run_status": "started"')
    # There are TWO completion sites -- the preflight's and the detection
    # path's -- and each must be stamped after ITS OWN outputs exist. Checking
    # only the first against the last CSV compared a preflight completion with
    # a detection output and mixed the two modes up.
    stamps = [
        i for i in range(len(src)) if src.startswith('_MANIFEST["run_status"] = "completed"', i)
    ]
    assert len(stamps) == 2, f"expected a preflight and a detection completion, got {len(stamps)}"
    for i_completed in stamps:
        assert i_started < i_completed
        csvs_before = [i for i in range(i_completed) if src.startswith("write_csv_atomic(", i)]
        assert csvs_before, "completion must be stamped after the outputs exist"

    # the preflight completes after its audit table and before the detector
    i_preflight = stamps[0]
    assert src.index("dimension_audit.csv") < i_preflight
    assert i_preflight < src.index("sc.pp.scrublet(scored")
    # the detection path completes after the last CSV of all
    assert src.rindex("write_csv_atomic(") < stamps[1]


def test_r5_executed_list_is_populated_only_where_the_analysis_runs():
    src = _src()
    for marker in ("doublet_rate_vs_mortality", "flagged_gene_vs_mortality"):
        idx = src.index(f'record_outcome_analysis("{marker}")')
        assert "INTEGRATION_ONLY" in src[:idx]
    # recording must go through the flushing helper: a bare append keeps the
    # fact in memory until the end, so a crash in between loses it
    assert (
        'OUTCOME_ANALYSIS_RAN.append("' not in src
    ), "appending without flushing loses the fact if the job dies later"


def test_r5_manifest_is_written_atomically():
    body = _ast.unparse(_fn("write_manifest"))
    assert (
        "os.replace" in body and ".inprogress" in body
    ), "the manifest needs the same durability as the CSVs"


# ---- AnnData nullable round trip


def test_r5_a_clean_is_revalidated_after_slicing():
    src = _src()
    assert (
        "validate_examination_state(a_clean.obs" in src
    ), "the nullable state must be re-checked after the AnnData slice"
    assert (
        "is_bool_dtype(a_clean.obs['doublet_call'].dtype)" in src
        or 'is_bool_dtype(a_clean.obs["doublet_call"].dtype)' in src
    )


def test_r5_h5ad_round_trip_is_checked_and_the_file_is_temporary():
    src = _src()
    assert (
        "write_h5ad(_rt)" in src and "read_h5ad(_rt)" in src
    ), "a disk round trip must actually happen"
    assert "_nullable_roundtrip_check.h5ad" in src
    i_write = src.index("write_h5ad(_rt)")
    assert "os.remove(_rt)" in src[i_write:], "the check file must be deleted, not kept"
    assert (
        "finally:" in src[i_write - 400 : i_write + 900]
    ), "deletion must happen even if the assertions fail"


def test_r5_round_trip_checks_dtype_and_na_count():
    src = _src()
    seg = src[src.index("write_h5ad(_rt)") : src.index("n_unknown_kept")]
    assert "is_bool_dtype" in seg, "dtype must survive the round trip"
    assert "isna().sum()" in seg, "the NA count must survive the round trip"
    assert "validate_examination_state(_back.obs" in seg


# ---- the output guard, as behaviour rather than prose


def test_r5_guard_refuses_legacy_and_nonempty(tmp_path):
    src = _src()
    assert "LEGACY_OUT_DIR" in src
    assert (
        "os.path.realpath(OUT_DIR) == os.path.realpath(LEGACY_OUT_DIR)" in src
    ), "both sides must be resolved, or a symlink slips through"
    assert "is not empty" in src, "a non-empty directory must be refused"


# ==========================================================================
# Round 6 — the outcome boundary as BEHAVIOUR, not as source strings
#
# Round 5's checks were source-level, and one of them made the defect the
# passing condition: it asserted `drop(columns=["y"])` was present, which is
# satisfied by parsing the whole file and then deleting the column. Nothing in
# it could distinguish "y was never read" from "y was read and discarded".
#
# These tests run the real loader against real fixture files and inspect the
# real manifest on disk. They also encode the POLICY for an input container
# that carries outcome: read_h5ad cannot select obs columns, so the guarantee
# for that input is incidental-load-then-drop, recorded in the manifest -- not
# "never loaded". A test that asserted the latter would be asserting something
# the reader cannot deliver.
# ==========================================================================

import json as _json

import anndata as _ad

#: Stand-in dimensionality for the fixtures and the feasibility tests. There is
#: no default in the script, so the tests must supply one exactly as a real run
#: does -- and must not depend on its particular value.
NPC = 7


def _fixture_inputs(
    tmp_path, *, obs_outcome_cols=("Death",), n_patients=240, n_events=25, cells_per_patient=4
):
    """Real CSVs and a real h5ad, small but the right shape.

    `obs_outcome_cols` is what the container carries. Their values are
    deliberately WRONG (all 9) so that any leak into the analysis frame is
    visible rather than coincidentally matching the patient list.
    """
    pids = [f"P{i:03d}" for i in range(n_patients)]
    y = [1] * n_events + [0] * (n_patients - n_events)
    pl = pd.DataFrame(
        {
            "patient_id": pids,
            "y": y,
            "disease": "COVID-19",
            "n_cells": 100,
            "median_counts": 1000.0,
            "median_genes": 500.0,
            "median_pct_mito": 5.0,
        }
    )
    pl_path = tmp_path / "patient_list.csv"
    pl.to_csv(pl_path, index=False)

    sm = pd.DataFrame({"sample_id": [f"{p}_s" for p in pids], "patient_id": pids})
    sm_path = tmp_path / "sample_map.csv"
    sm.to_csv(sm_path, index=False)

    obs = pd.DataFrame({"Sample_ID": np.repeat([f"{p}_s" for p in pids], cells_per_patient)})
    for col in obs_outcome_cols:
        obs[col] = 9
    obs.index = [f"c{i}" for i in range(len(obs))]
    h5_path = tmp_path / "raw.h5ad"
    _ad.AnnData(X=np.ones((len(obs), 3), dtype=np.int64), obs=obs).write_h5ad(h5_path)
    return str(pl_path), str(sm_path), str(h5_path)


def _wire(monkeypatch, tmp_path, *, integration_only, obs_outcome_cols=("Death",)):
    """Point the module at the fixtures and give it a fresh, real manifest.

    The manifest helpers are NOT stubbed: they run for real and write real
    JSON, so these tests exercise the atomic write too.
    """
    pl, sm, h5 = _fixture_inputs(tmp_path, obs_outcome_cols=obs_outcome_cols)
    monkeypatch.setattr(s31, "PATIENT_LIST_CSV", pl, raising=False)
    monkeypatch.setattr(s31, "SAMPLE_MAP_CSV", sm, raising=False)
    monkeypatch.setattr(s31, "RAW_H5AD", h5, raising=False)
    monkeypatch.setattr(s31, "INTEGRATION_ONLY", integration_only, raising=False)
    monkeypatch.setattr(s31, "CELLS_PER_PATIENT", 2, raising=False)
    monkeypatch.setattr(s31, "OUT_DIR", str(tmp_path), raising=False)
    monkeypatch.setattr(s31, "_commit", "test", raising=False)
    # Everything derived from argv, set together. The test loader cannot
    # evaluate `X = _args.…` (parse_args exits without a command line), so any
    # such name is simply absent from the module and the first function to
    # reference it raises NameError. Twice now that has surfaced as a wave of
    # unrelated-looking failures -- once for N_PRIN_COMPS, once for
    # DIMENSION_AUDIT -- so they are listed in one place.
    monkeypatch.setattr(s31, "INTEGRATION_ONLY", integration_only, raising=False)
    monkeypatch.setattr(s31, "DIMENSION_AUDIT", False, raising=False)
    monkeypatch.setattr(s31, "N_PRIN_COMPS", NPC, raising=False)
    # through scrublet_params, not as a literal dict, so an invalid value would
    # fail here exactly as it would in a real run
    monkeypatch.setattr(s31, "SCRUBLET_PARAMS", s31.scrublet_params(NPC), raising=False)
    # the SCRIPT's own starting values, not a dict written here: a fixture that
    # supplies its own manifest cannot detect a wrong starting value, which is
    # exactly how a predicted `outcome_data_accessed` would survive this suite
    manifest = s31.initial_manifest()
    monkeypatch.setattr(s31, "_MANIFEST", manifest, raising=False)
    monkeypatch.setattr(s31, "MANIFEST_PATH", str(tmp_path / "RUN_MANIFEST.json"), raising=False)
    monkeypatch.setattr(s31, "OUTCOME_ANALYSIS_RAN", [], raising=False)
    return manifest


def _on_disk(tmp_path):
    with open(tmp_path / "RUN_MANIFEST.json") as fh:
        return _json.load(fh)


# ---- the patient list: exclusion happens at parse time --------------------


def test_r6_integration_only_passes_usecols_without_y(monkeypatch, tmp_path):
    """The discriminating test: the previous 'read all, then drop' version
    passed NO usecols at all and would fail here."""
    _wire(monkeypatch, tmp_path, integration_only=True)
    seen = []
    real = pd.read_csv

    def spy(path, **kw):
        seen.append((os.path.basename(str(path)), kw.get("usecols")))
        return real(path, **kw)

    monkeypatch.setattr(s31.pd, "read_csv", spy)
    a, y_by_patient, _ = s31.load_canonical_cohort()

    usecols = dict(seen)["patient_list.csv"]
    assert usecols is not None, "the patient list must be parsed with usecols"
    assert "y" not in usecols, "y must be excluded at parse time"
    assert y_by_patient is None
    assert "y" not in a.obs.columns


def test_r6_full_mode_passes_usecols_including_y(monkeypatch, tmp_path):
    _wire(monkeypatch, tmp_path, integration_only=False)
    seen = []
    real = pd.read_csv

    def spy(path, **kw):
        seen.append((os.path.basename(str(path)), kw.get("usecols")))
        return real(path, **kw)

    monkeypatch.setattr(s31.pd, "read_csv", spy)
    a, y_by_patient, _ = s31.load_canonical_cohort()

    assert "y" in dict(seen)["patient_list.csv"]
    assert y_by_patient is not None and sum(y_by_patient.values()) == 25
    assert "y" in a.obs.columns


def test_r6_a_renamed_outcome_column_fails_the_run(monkeypatch, tmp_path):
    """usecols is fail-closed: a missing name raises rather than yielding a
    frame silently short of a column."""
    _wire(monkeypatch, tmp_path, integration_only=False)
    pl = pd.read_csv(s31.PATIENT_LIST_CSV)
    pl.rename(columns={"y": "outcome"}).to_csv(s31.PATIENT_LIST_CSV, index=False)
    with pytest.raises(ValueError):
        s31.load_canonical_cohort()


# ---- the h5ad container: retention, not a deny-list ----------------------


def test_r6_only_required_obs_columns_survive_the_loader(monkeypatch, tmp_path):
    """The boundary is enforced by RETENTION. A deny-list is only as good as
    the classification behind it -- the first version listed ("y", "Death") and
    let covid_severity through -- so the loader keeps what it needs and drops
    everything else, classified or not."""
    _wire(
        monkeypatch,
        tmp_path,
        integration_only=True,
        obs_outcome_cols=(
            "Death",
            "covid_severity",
            "mild_vs_severe",
            "severity_score",
            "three_group",
            "progression_group",
        ),
    )
    a, _, _ = s31.load_canonical_cohort()

    # patient_id is created by the loader; nothing else from the container remains
    assert set(a.obs.columns) == {"Sample_ID", "patient_id"}


def test_r6_an_unclassified_clinical_label_cannot_leak(monkeypatch, tmp_path):
    """The point of retention: a label nobody has adjudicated is dropped
    because it is not required, not because it was recognised as clinical."""
    _wire(
        monkeypatch,
        tmp_path,
        integration_only=True,
        obs_outcome_cols=("Death", "who_ordinal_scale_day28"),
    )
    a, _, _ = s31.load_canonical_cohort()

    assert "who_ordinal_scale_day28" not in a.obs.columns
    rec = _on_disk(tmp_path)["container_obs_columns"]
    assert (
        "who_ordinal_scale_day28" in rec["unclassified_dropped"]
    ), "an unadjudicated column must be reported as unclassified, not filed as technical"
    assert rec["known_outcome_dropped"] == ["Death"]


def test_r6_severity_columns_are_classified_as_outcome_bearing(monkeypatch, tmp_path):
    """covid_severity and its relatives are clinical outcome. The code and the
    docs say 'outcome' throughout, so the report must match that scope."""
    sev = ("covid_severity", "mild_vs_severe", "severity_score", "three_group", "progression_group")
    for col in sev:
        assert col in s31.KNOWN_OUTCOME_OBS_COLUMNS, f"{col} is a severity label"
    _wire(monkeypatch, tmp_path, integration_only=True, obs_outcome_cols=("Death",) + sev)
    s31.load_canonical_cohort()
    rec = _on_disk(tmp_path)["container_obs_columns"]
    assert set(rec["known_outcome_dropped"]) == {"Death", *sev}
    assert rec["unclassified_dropped"] == []


def test_r6_the_container_load_flips_loaded_but_never_inspected(monkeypatch, tmp_path):
    """The contradiction this split exists to remove: an earlier manifest said
    outcome_data_accessed=false beside outcome_incidentally_loaded=["Death"].
    Loaded and inspected are now separate facts and both are stated."""
    _wire(
        monkeypatch, tmp_path, integration_only=True, obs_outcome_cols=("Death", "covid_severity")
    )
    s31.load_canonical_cohort()
    disk = _on_disk(tmp_path)

    assert (
        disk["known_outcome_values_loaded"] is True
    ), "read_h5ad brought Death into memory; saying otherwise is false"
    assert set(disk["known_outcome_values_loaded_via"]) == {
        "raw_h5ad_obs:Death",
        "raw_h5ad_obs:covid_severity",
    }
    assert disk["outcome_values_inspected"] is False
    assert disk["outcome_values_inspected_via"] == []
    assert disk["outcome_analyses_executed"] == []


def test_r6_a_container_y_is_dropped_too_and_never_reaches_the_analysis(monkeypatch, tmp_path):
    """The original code dropped "y" -- a column this h5ad does not have --
    while `Death` sat in obs untouched."""
    _wire(monkeypatch, tmp_path, integration_only=True, obs_outcome_cols=("Death", "y"))
    a, _, _ = s31.load_canonical_cohort()

    assert "y" not in a.obs.columns and "Death" not in a.obs.columns
    assert set(_on_disk(tmp_path)["container_obs_columns"]["known_outcome_dropped"]) == {
        "y",
        "Death",
    }


def test_r6_full_mode_y_comes_from_the_patient_list_not_the_container(monkeypatch, tmp_path):
    """The container's y is 9 everywhere. If any of it survived, obs would
    carry a 9."""
    _wire(monkeypatch, tmp_path, integration_only=False, obs_outcome_cols=("Death", "y"))
    a, y_by_patient, _ = s31.load_canonical_cohort()

    assert set(a.obs["y"].unique()) <= {0, 1}, "the container value leaked into obs"
    for pid, grp in a.obs.groupby("patient_id", observed=True):
        assert set(grp["y"]) == {y_by_patient[pid]}


def test_r6_recording_the_empty_case_distinguishes_it_from_not_reached(monkeypatch, tmp_path):
    """`null` must mean the load was never reached; empty lists must mean it
    was reached and the container carried nothing clinical."""
    manifest = _wire(monkeypatch, tmp_path, integration_only=True, obs_outcome_cols=())
    assert manifest["container_obs_columns"] is None
    s31.load_canonical_cohort()
    disk = _on_disk(tmp_path)
    assert disk["container_obs_columns"]["known_outcome_dropped"] == []
    assert (
        disk["known_outcome_values_loaded"] is False
    ), "nothing clinical arrived, so nothing was loaded"


def test_r6_a_missing_required_obs_column_is_fatal(monkeypatch, tmp_path):
    _wire(monkeypatch, tmp_path, integration_only=True)
    import anndata as _a

    obs = pd.DataFrame({"Death": [0, 1]}, index=["c0", "c1"])
    _a.AnnData(X=np.ones((2, 3), dtype=np.int64), obs=obs).write_h5ad(s31.RAW_H5AD)
    with pytest.raises(SystemExit, match="missing required column"):
        s31.load_canonical_cohort()


# ---- inspection is a fact, recorded when it happens ----------------------


def test_r6_accessed_stays_false_through_an_integration_only_load(monkeypatch, tmp_path):
    _wire(monkeypatch, tmp_path, integration_only=True)
    s31.load_canonical_cohort()
    disk = _on_disk(tmp_path)
    assert disk["outcome_values_inspected"] is False
    assert disk["outcome_values_inspected_via"] == []
    assert disk["outcome_analysis_requested"] is False


def test_r6_accessed_flips_only_when_the_column_is_really_read(monkeypatch, tmp_path):
    manifest = _wire(monkeypatch, tmp_path, integration_only=False)
    assert manifest["outcome_values_inspected"] is False, 'must start at "has not happened"'
    s31.load_canonical_cohort()
    disk = _on_disk(tmp_path)
    assert disk["outcome_values_inspected"] is True
    assert disk["outcome_values_inspected_via"] == ["patient_list_csv:y"]


def test_r6_a_full_mode_job_that_dies_before_the_read_claims_no_access(monkeypatch, tmp_path):
    """The old `"outcome_data_accessed": not INTEGRATION_ONLY` recorded an
    access for a job that never got there."""
    _wire(monkeypatch, tmp_path, integration_only=False)
    s31.write_manifest()  # the up-front write, then the job dies
    disk = _on_disk(tmp_path)
    assert disk["outcome_analysis_requested"] is True, "authorization is known up front"
    assert disk["outcome_values_inspected"] is False, "but the access has not happened"


def test_r6_access_is_recorded_once_not_appended_per_call(monkeypatch, tmp_path):
    _wire(monkeypatch, tmp_path, integration_only=False)
    s31.mark_outcome_values_inspected("patient_list_csv:y")
    s31.mark_outcome_values_inspected("patient_list_csv:y")
    assert _on_disk(tmp_path)["outcome_values_inspected_via"] == ["patient_list_csv:y"]


# ---- partial execution must reach disk immediately -----------------------


def test_r6_an_analysis_that_ran_survives_a_later_crash(monkeypatch, tmp_path):
    """Run one outcome analysis, then let the job die. The manifest must say
    BOTH that the run did not finish AND that this analysis really ran. The
    in-memory-only version left `[]` here."""
    _wire(monkeypatch, tmp_path, integration_only=False)
    s31.write_manifest()
    s31.record_outcome_analysis("doublet_rate_vs_mortality")
    # ... the gene re-audit now raises and nothing further touches the manifest

    disk = _on_disk(tmp_path)
    assert disk["run_status"] == "started", "an unfinished run must not say completed"
    assert disk["outcome_analyses_executed"] == ["doublet_rate_vs_mortality"]
    assert disk["finished"] is None


def test_r6_each_analysis_reaches_disk_as_it_runs(monkeypatch, tmp_path):
    _wire(monkeypatch, tmp_path, integration_only=False)
    s31.record_outcome_analysis("doublet_rate_vs_mortality")
    assert _on_disk(tmp_path)["outcome_analyses_executed"] == ["doublet_rate_vs_mortality"]
    s31.record_outcome_analysis("flagged_gene_vs_mortality")
    assert _on_disk(tmp_path)["outcome_analyses_executed"] == [
        "doublet_rate_vs_mortality",
        "flagged_gene_vs_mortality",
    ]


def test_r6_repeated_recording_does_not_duplicate(monkeypatch, tmp_path):
    """The gene loop calls this once per gene."""
    _wire(monkeypatch, tmp_path, integration_only=False)
    for _ in range(8):
        s31.record_outcome_analysis("flagged_gene_vs_mortality")
    assert _on_disk(tmp_path)["outcome_analyses_executed"] == ["flagged_gene_vs_mortality"]


def test_r6_manifest_write_leaves_no_temp_behind(monkeypatch, tmp_path):
    _wire(monkeypatch, tmp_path, integration_only=True)
    s31.write_manifest()
    assert not [p for p in os.listdir(tmp_path) if p.endswith(".inprogress")]


# ---- the drop record is an accomplished fact, not an intention -------------


def test_r6_container_record_reflects_the_object_after_the_drop(monkeypatch, tmp_path):
    """`kept` must be what is observed on the object, not the constant that was
    requested. A version that recorded before applying the allowlist would put
    every arrived column in `kept` and fail here."""
    _wire(
        monkeypatch, tmp_path, integration_only=True, obs_outcome_cols=("Death", "covid_severity")
    )
    s31.load_canonical_cohort()
    rec = _on_disk(tmp_path)["container_obs_columns"]

    assert rec["kept"] == [
        "Sample_ID"
    ], "kept must be the post-drop columns, so the record cannot be right while the object is wrong"
    assert set(rec["arrived"]) == {"Sample_ID", "Death", "covid_severity"}
    assert "Death" not in rec["kept"] and "covid_severity" not in rec["kept"]


def test_r6_loading_is_recorded_before_the_drop_record(monkeypatch, tmp_path):
    """Two records on opposite sides of the drop. If the process dies between
    them, the manifest must say the values were loaded and must NOT yet claim
    any column was dropped."""
    _wire(monkeypatch, tmp_path, integration_only=True, obs_outcome_cols=("Death",))

    def boom(**kwargs):
        raise RuntimeError("died after the drop, before the record")

    monkeypatch.setattr(s31, "record_container_obs_handling", boom)
    with pytest.raises(RuntimeError):
        s31.load_canonical_cohort()

    disk = _on_disk(tmp_path)
    assert (
        disk["known_outcome_values_loaded"] is True
    ), "read_h5ad had returned, so the values were in memory"
    assert (
        disk["container_obs_columns"] is None
    ), "nothing may claim a drop that has not been recorded as done"


def test_r6_verification_sits_between_the_drop_and_the_record():
    """Order, asserted on the source because selection and verification use the
    same constant and so cannot be made to disagree from a test. It guards a
    future refactor in which they diverge: drop, then verify, then record."""
    # _ast.unparse, not inspect.getsource: the module is exec'd from compiled
    # AST blocks, so it has no source file to read back.
    src = _ast.unparse(_fn("load_canonical_cohort"))
    assert "allowlist did not apply" in src, "the drop must be verified, not assumed"
    i_drop = src.index("a.obs = a.obs[list(REQUIRED_OBS_COLUMNS)].copy()")
    i_verify = src.index("allowlist did not apply")
    i_record = src.index("record_container_obs_handling(")
    assert i_drop < i_verify < i_record, "verify after dropping, record after verifying"


def test_r6_known_flag_speaks_only_for_the_classified_list(monkeypatch, tmp_path):
    """A container carrying ONLY an unrecognised clinical label leaves
    known_outcome_values_loaded false. That is why "known_" is in the name --
    the field must not be read as "no outcome value was ever in memory"."""
    _wire(
        monkeypatch, tmp_path, integration_only=True, obs_outcome_cols=("who_ordinal_scale_day28",)
    )
    a, _, _ = s31.load_canonical_cohort()
    disk = _on_disk(tmp_path)

    assert disk["known_outcome_values_loaded"] is False
    assert disk["container_obs_columns"]["unclassified_dropped"] == ["who_ordinal_scale_day28"]
    assert (
        "who_ordinal_scale_day28" not in a.obs.columns
    ), "the anti-leak guarantee comes from retention, not from the flag"
    assert (
        "classification_scope" in disk["container_obs_columns"]
    ), "the scope limitation must be stated in the manifest itself"


# ==========================================================================
# Round 7 — fixed dimensionality: feasibility, not a caught crash
#
# The v1 run failed on n_components=30 vs min(n_samples, n_features). The
# constraint is a DIMENSION fact knowable before Scrublet is called, and the
# binding term is not the raw cell count -- one v1 batch had
# min(n_cells_effective, n_usable_features)=17; which term that was is not
# recorded, and reading it as "17 usable features" was an over-inference.
# ==========================================================================


@pytest.mark.parametrize(
    "n_cells,n_feat,ok",
    [
        (NPC, NPC + 5, False),  # n_cells == n_prin_comps
        (NPC + 1, NPC + 5, True),  # n_cells == n_prin_comps + 1
        (NPC + 5, NPC, False),  # n_usable_features == n_prin_comps
        (NPC + 5, NPC + 1, True),  # n_usable_features == n_prin_comps + 1
    ],
)
def test_r7_feasibility_is_strict_on_both_terms(n_cells, n_feat, ok):
    """arpack needs n_components STRICTLY less than min(n_samples, n_features);
    equality raises, so equality must be refused."""
    got, why = s31.dimension_feasibility(n_cells, n_feat, NPC)
    assert got is ok, why
    if not ok:
        assert str(NPC) in why and "strictly greater" in why


def test_r7_the_binding_term_can_be_features_not_cells():
    """A sample with plenty of cells and almost no HVGs is infeasible. Judging
    by cell count alone -- which is what raising MIN_CELLS_PER_SAMPLE does --
    would have passed this."""
    ok, why = s31.dimension_feasibility(600, NPC - 1, NPC)
    assert ok is False
    assert f"usable_features={NPC - 1}" in why


def test_r7_insufficient_dimension_is_its_own_unknown_status():
    assert s31.EXAM_INSUFFICIENT_DIMENSION in s31.EXAM_UNKNOWN_STATUSES
    assert s31.EXAM_INSUFFICIENT_DIMENSION != s31.EXAM_RUNTIME_FAILED
    assert s31.EXAM_INSUFFICIENT_DIMENSION not in (
        s31.EXAM_EXAMINED,
        s31.EXAM_TOO_SMALL,
        s31.EXAM_THRESHOLD_FAILED,
        s31.EXAM_NOT_EVALUATED,
    )


def test_r7_insufficient_dimension_cells_stay_na():
    """It is an unknown status, so the invariant must force NA -- a dimension
    failure may not leave a call behind."""
    obs = make_obs([("p1", s31.EXAM_INSUFFICIENT_DIMENSION, None), ("p1", s31.EXAM_EXAMINED, True)])
    s31.validate_examination_state(obs)
    r = s31.doublet_rate_per_patient(obs)
    assert r.loc["p1", "doublet_rate"] == 1.0, "it must not dilute the denominator"
    assert r.loc["p1", "n_unknown"] == 1


@pytest.mark.parametrize("call", [True, False])
def test_r7_a_dimension_failure_carrying_a_call_is_refused(call):
    obs = make_obs([("p1", s31.EXAM_INSUFFICIENT_DIMENSION, call)])
    with pytest.raises(s31.ExaminationStateError):
        s31.validate_examination_state(obs)


def test_r7_dimension_failure_is_not_a_runtime_failure():
    """The two must stay distinguishable: 'the detector broke' and 'we never
    had the dimensions to ask' are different facts about a sample."""
    src = _src()
    i = src.index("if sid in insufficient:")
    seg = src[i : i + 400]
    assert "EXAM_INSUFFICIENT_DIMENSION" in seg
    assert "EXAM_RUNTIME_FAILED" not in seg
    # and the feasibility verdict is reached before scrublet is called
    assert src.index("dimension_feasibility(") < src.index("sc.pp.scrublet(")


def test_r7_one_fixed_dimension_reaches_every_call_site():
    """No dynamic reduction, and no call site silently taking scanpy's default
    -- which is what v1 did at both of them."""
    src = _src()
    calls = [line for line in src.splitlines() if "sc.pp.scrublet(" in line]
    assert len(calls) >= 2, "expected the batched path and the fallback"
    for line in calls:
        assert "**SCRUBLET_PARAMS" in line, f"unbound dimensionality: {line.strip()}"
    assert s31.scrublet_params(NPC) == {"n_prin_comps": NPC}
    # A retry at a smaller dimension would defeat cross-sample comparability, so
    # no literal dimensionality may appear anywhere -- the constant is the only
    # source. (The f-string `n_prin_comps={...}` is not a literal.)
    import re

    assert not re.search(
        r"n_prin_comps\s*=\s*\d", src
    ), "a hard-coded dimensionality bypasses the required argument"
    assert (
        src.count("N_PRIN_COMPS = _args.n_prin_comps") == 1
    ), "the value is bound once, from the command that was run"


def test_r7_feasibility_is_computed_once_and_shared_by_both_paths():
    src = _src()
    # one definition, one call site -- two evaluations could disagree about
    # which samples were askable
    assert src.count("def dimension_feasibility(") == 1
    assert (
        src.count("= dimension_feasibility(") == 1
    ), "two evaluations could disagree about which samples were askable"
    i_feas = src.index("= dimension_feasibility(")
    assert src.index("if sid in insufficient:") > i_feas
    assert (
        src.count("if sid in insufficient:") == 2
    ), "both the batched loop and the fallback must honour the same verdict"


def test_r7_usable_features_uses_scanpys_own_preprocessing():
    """Not an approximation of the pipeline -- the same calls in the same
    order, because filter_cells changes n_cells and HVG selection sets
    n_features."""
    # AST call nodes, NOT a substring of the unparsed source: this function's
    # docstring spells out the pipeline, so a substring check passes even when
    # every call has been deleted. (Found by mutation -- replacing the HVG
    # count with `ad.n_vars` survived the substring version of this test.)
    fn = _fn("scrublet_preprocessing_dimensions")
    called = {
        n.func.attr
        for n in _ast.walk(fn)
        if isinstance(n, _ast.Call) and isinstance(n.func, _ast.Attribute)
    }
    for call in (
        "filter_genes",
        "filter_cells",
        "normalize_total",
        "log1p",
        "highly_variable_genes",
    ):
        assert call in called, f"{call} is part of scanpy scrublet preprocessing"

    kwargs = {
        (k.arg, getattr(k.value, "value", None))
        for n in _ast.walk(fn)
        if isinstance(n, _ast.Call)
        for k in n.keywords
    }
    assert ("min_cells", 3) in kwargs and (
        "min_genes",
        3,
    ) in kwargs, "the gene/cell filters must match scanpy exactly"

    # and the returned feature count must come from the HVG mask, not n_vars
    returns = [_ast.unparse(n) for n in _ast.walk(fn) if isinstance(n, _ast.Return)]
    assert any(
        "highly_variable" in r for r in returns
    ), "usable_features is the HVG count; ad.n_vars is the pre-selection count"
    assert not any("n_vars" in r for r in returns)


def test_r7_diagnostics_table_has_the_required_columns():
    src = _src()
    for col in (
        "Sample_ID",
        "n_cells",
        "n_input_features",
        "n_usable_features",
        "n_prin_comps_requested",
        "scrublet_attempted",
        "n_prin_comps_used_for_attempt",
        "threshold",
        "examination_status",
        "failure_reason",
        "n_examined",
        "n_doublet",
    ):
        assert f'"{col}"' in src, f"diagnostics must report {col}"
    assert "per_sample_diagnostics.csv" in src


def test_r7_diagnostics_are_reconciled_against_the_cell_level_counts():
    """A diagnostics table that disagrees with the cells is describing a run
    that did not happen, so the disagreement must be fatal."""
    src = _src()
    # anchored on the statements, not on a filename that also appears in a
    # comment 500 lines earlier -- index() found that one and sliced backwards
    seg = src[src.index("_by_status = ") : src.index("write_csv_atomic(_diag")]
    assert "assert" in seg
    assert "counts['by_status']" in seg or 'counts["by_status"]' in seg
    assert "counts['n_examined']" in seg or 'counts["n_examined"]' in seg
    assert "counts['n_doublet']" in seg or 'counts["n_doublet"]' in seg


def test_r7_the_minimum_cell_rule_has_one_definition():
    """v1 had MIN_CELLS_PER_SAMPLE = 30 and a separate literal `< 30`."""
    src = _src()
    assert "sample_sizes < MIN_CELLS_PER_SAMPLE" in src
    assert "sample_sizes < 30" not in src


def test_r7_the_dimension_value_is_not_justified_in_code():
    """The value is a methodological choice; the code must point at the note
    rather than assert a rationale of its own."""
    src = _src()
    assert "docs/scrublet_dimension_parameter.md" in src
    assert "N_PRIN_COMPS = 30" not in src, "a default value makes an unapproved parameter runnable"


# ---- the dimensionality must not be runnable before it is approved ---------


def test_r8_there_is_no_default_dimensionality():
    """docs/scrublet_dimension_parameter.md is OPEN. A provisional constant
    recorded in the manifest still lets anyone produce results under an
    unapproved parameter -- "provisional" in a comment binds nothing."""
    assert s31.N_PRIN_COMPS_UNSET is None
    with pytest.raises(SystemExit, match="required"):
        s31.scrublet_params(None)


@pytest.mark.parametrize("bad", [0, -1, 2.5, "30", True, None])
def test_r8_an_invalid_dimensionality_is_refused(bad):
    with pytest.raises(SystemExit):
        s31.scrublet_params(bad)


def test_r8_a_valid_dimensionality_binds_exactly_one_parameter():
    assert s31.scrublet_params(12) == {"n_prin_comps": 12}


def test_r8_the_cli_refuses_to_detect_without_a_dimensionality():
    """Not argparse `required=True` any more -- the preflight legitimately runs
    without one -- so the refusal must be explicit and mode-aware."""
    src = _src()
    i = src.index('"--n-prin-comps"')
    assert "default=None" in src[i : i + 200], "no default value may be supplied"
    # detection mode calls the fail-closed validator; the preflight does not
    assert "None if DIMENSION_AUDIT else scrublet_params(N_PRIN_COMPS)" in src


def test_r8_the_wrapper_requires_it_too():
    """A wrapper with a baked-in value would reintroduce the default the script
    just removed."""
    sh = open(os.path.join(CODE, "run_doublet_detection.sh")).read()
    assert "--n-prin-comps" in sh
    assert 'N_PRIN_COMPS="$1"' in sh, "it comes from the caller, not from the wrapper"
    assert "[ $# -lt 2 ]" in sh
    import re as _re

    assert not _re.search(r"--n-prin-comps\s+\d", sh), "no hard-coded value in the wrapper"


# ---- attempted vs rejected before the call --------------------------------


def test_r8_a_threshold_failure_still_records_the_dimension_it_ran_at():
    """Requested by audit: a sample rejected on its threshold DID enter the
    detector at the fixed dimension. Reporting an empty dimension there
    described a run that did not happen."""
    src = _src()
    seg = src[src.index('row["scrublet_attempted"]') : src.index("_diag_cols")]
    assert "if sid in SAMPLE_ATTEMPTED:" in seg
    assert 'n_prin_comps_used_for_attempt"] = N_PRIN_COMPS' in seg
    # the dimension is filled by ATTEMPT, not by success
    i_attempt = seg.index("if sid in SAMPLE_ATTEMPTED:")
    i_examined = seg.index("if st == EXAM_EXAMINED:")
    assert i_attempt < i_examined
    assert "EXAM_EXAMINED" not in seg[:i_attempt]


def test_r8_only_pre_call_rejections_may_be_unattempted():
    """too_small and insufficient_dimension are decided before the call and
    have no dimension; any other unexamined status must have been attempted,
    and the script asserts it rather than trusting it."""
    src = _src()
    seg = src[src.index('row["scrublet_attempted"]') : src.index("_diag_cols")]
    assert "assert st in (EXAM_TOO_SMALL, EXAM_INSUFFICIENT_DIMENSION," in seg
    assert "EXAM_DIMENSION_MEASUREMENT_FAILED), (" in seg


def test_r8_attempt_is_recorded_before_the_call_not_after():
    """A call that raises still attempted the sample. Recording afterwards
    would report a crashed sample as never attempted."""
    src = _src()
    for call in ("sc.pp.scrublet(scored", "sc.pp.scrublet(sub"):
        i = src.index(call)
        before = src[max(0, i - 600) : i]
        assert "SAMPLE_ATTEMPTED" in before, f"attempt not recorded before {call}"


def test_r8_falling_back_clears_the_batched_attempt_records():
    body = _src()
    seg = body[body.index("reset_detection_state(a.obs)") :]
    for name in ("SAMPLE_THRESHOLD.clear()", "SAMPLE_FAILURE.clear()", "SAMPLE_ATTEMPTED.clear()"):
        assert name in seg[:600], f"{name} missing: stale batched state would be reported"


# ==========================================================================
# Round 9 — the preflight, because the parameter rule was circular
#
# "Refuse to run without an approved dimensionality" and "choose it from this
# run's measured per-sample dimensions" cannot both hold: the distribution was
# only obtainable by first picking a value arbitrarily. --dimension-audit
# measures without detecting.
# ==========================================================================


def test_r9_the_preflight_needs_no_dimensionality():
    src = _src()
    i = src.index('"--dimension-audit"')
    assert 'action="store_true"' in src[i : i + 200]
    # and it does not construct scrublet params at all
    assert "SCRUBLET_PARAMS = None if DIMENSION_AUDIT else" in src


def test_r9_the_preflight_refuses_a_dimensionality():
    """Accepting one would measure the distribution under the very value it is
    supposed to inform."""
    src = _src()
    assert "if DIMENSION_AUDIT and _args.n_prin_comps is not None:" in src
    i = src.index("if DIMENSION_AUDIT and _args.n_prin_comps is not None:")
    assert "raise SystemExit" in src[i : i + 400]


def test_r9_the_preflight_requires_integration_only():
    src = _src()
    assert "if DIMENSION_AUDIT and not _args.integration_only:" in src
    i = src.index("if DIMENSION_AUDIT and not _args.integration_only:")
    assert "raise SystemExit" in src[i : i + 400]


def test_r9_the_preflight_stops_before_any_detector_call():
    src = _src()
    i_exit = src.index("if DIMENSION_AUDIT:\n    # ----")
    i_first_scrublet = src.index("sc.pp.scrublet(scored")
    assert i_exit < i_first_scrublet, "the preflight must exit before the detector"
    seg = src[i_exit:i_first_scrublet]
    assert "raise SystemExit(_exit)" in seg
    assert "sc.pp.scrublet" not in src[i_exit : src.index("raise SystemExit(_exit)")]


def test_r9_the_preflight_produces_no_call_rate_or_contamination_output():
    src = _src()
    seg = src[src.index("if DIMENSION_AUDIT:\n    # ----") : src.index("raise SystemExit(_exit)")]
    assert "dimension_audit.csv" in seg
    for forbidden in (
        "doublet_scores_per_cell.csv",
        "doublet_rate_per_patient.csv",
        "contamination_audit_v2_doublet_excluded.csv",
        "per_sample_diagnostics.csv",
    ):
        assert forbidden not in seg, f"the preflight must not write {forbidden}"


def test_r9_the_preflight_reports_only_measured_dimensions():
    src = _src()
    seg = src[src.index("_audit_cols = [") : src.index("raise SystemExit(_exit)")]
    for col in (
        "Sample_ID",
        "n_cells",
        "n_input_features",
        "n_cells_effective",
        "n_usable_features",
    ):
        assert f'"{col}"' in seg
    # no verdict: there is no approved value to judge against
    assert "insufficient" not in seg and "examination_status" not in seg


def test_r9_no_feasibility_verdict_is_reached_in_the_preflight():
    src = _src()
    i = src.index("= dimension_feasibility(")
    assert (
        "if not DIMENSION_AUDIT:" in src[max(0, i - 400) : i]
    ), "judging feasibility needs a value, and the preflight has none"


def test_r9_attempted_and_completed_are_separate_manifest_facts():
    """A single `scrublet_executed` set before the call claimed execution for a
    detector that raised at its entry point. Attempt precedes the call because
    a crash did attempt it; completion follows a successful return."""
    src = _src()
    # the field is retired as a FIELD; the comment explaining why it was
    # retired legitimately names it
    assert '"scrublet_executed"' not in src, "the conflated field is retired"
    assert '_MANIFEST["scrublet_executed"]' not in src
    assert '"scrublet_attempted": False,' in src and '"scrublet_any_call_completed": False,' in src

    for call in ("sc.pp.scrublet(scored", "sc.pp.scrublet(sub"):
        i = src.index(call)
        assert (
            '_MANIFEST["scrublet_attempted"] = True' in src[max(0, i - 700) : i]
        ), f"attempt must be recorded BEFORE {call}"
        assert (
            '_MANIFEST["scrublet_any_call_completed"] = True' not in src[max(0, i - 700) : i]
        ), f"completion must not be claimed before {call}"
        assert (
            '_MANIFEST["scrublet_any_call_completed"] = True' in src[i : i + 700]
        ), f"completion must follow a successful {call}"

    seg = src[src.index("if DIMENSION_AUDIT:\n    # ----") : src.index("raise SystemExit(_exit)")]
    assert '_MANIFEST["scrublet_attempted"] = False' in seg
    assert '_MANIFEST["scrublet_any_call_completed"] = False' in seg


def test_r9_every_sample_is_measured_not_only_the_eligible_ones():
    """Gating measurement on the current cell-count rule would redefine "the
    per-sample dimension distribution" as "the distribution over the samples
    that rule already admits" -- the decision leaking into its own evidence."""
    src = _src()
    i_meas = src.index("scrublet_preprocessing_dimensions(a[idx])")
    before = src[src.index("for sid, idx in sample_groups:") : i_meas]
    assert "continue" not in before, "no sample may be skipped before measurement"
    # eligibility is recorded beside the measurement, not used to skip it
    assert 'row["meets_min_cells_rule"] = sid not in too_small_set' in src
    assert '"preflight_exclusion_reason"' in src


def test_r9_a_measurement_failure_is_recorded_not_fatal():
    src = _src()
    i = src.index("scrublet_preprocessing_dimensions(a[idx])")
    seg = src[i : i + 600]
    assert "except Exception" in seg
    assert "dimension_measurement_error" in seg


def test_r9_both_distributions_are_reported_and_labelled():
    """They are not interchangeable, so neither may be printed as "the"
    distribution."""
    src = _src()
    seg = src[src.index("_audit_cols = [") : src.index("raise SystemExit(_exit)")]
    assert "ALL measured samples" in seg
    assert "meeting the current cell-count rule" in seg
    for col in (
        "meets_min_cells_rule",
        "preflight_exclusion_reason",
        "dimension_measurement_error",
    ):
        assert f'"{col}"' in seg, f"{col} must travel with the measurements"


def test_r9_the_parameter_is_refused_for_governance_not_measurement():
    """n_prin_comps does not enter the preprocessing, so it cannot change a
    measured number. Saying it would was wrong, and a wrong reason invites
    someone to "fix" the guard later."""
    src = _src()
    i = src.index("if DIMENSION_AUDIT and _args.n_prin_comps is not None:")
    seg = src[max(0, i - 900) : i + 500]
    assert "does not affect what is measured" in seg or "does not enter the preprocessing" in seg
    assert "measures the dimensions available, which is what the value is chosen from" not in src


def test_r9_the_wrapper_supports_the_preflight_without_a_dimensionality():
    sh = open(os.path.join(CODE, "run_doublet_detection.sh")).read()
    assert "--dimension-audit" in sh
    assert 'N_PRIN_COMPS=""' in sh, "the preflight form takes no dimensionality"
    assert (
        "${N_PRIN_COMPS:+--n-prin-comps" in sh
    ), "the flag must be omitted entirely, not passed empty"


# ==========================================================================
# Round 10 — the measurement failure path, and what makes a preflight usable
# ==========================================================================


def test_r10_a_failed_measurement_cannot_reach_feasibility_unbound():
    """The names must be bound before the try. Without that, an ELIGIBLE sample
    whose measurement raised fell through to dimension_feasibility(n_eff, ...)
    with n_eff unbound -- an UnboundLocalError that also destroyed the recorded
    reason."""
    src = _src()
    i_try = src.index("n_eff, n_usable = scrublet_preprocessing_dimensions(a[idx])")
    assert (
        "n_eff = n_usable = None" in src[max(0, i_try - 500) : i_try]
    ), "bind before the try, or the except path leaves them unbound"
    # and an unmeasured eligible sample must leave the loop before the verdict
    i_guard = src.index("if n_eff is None:")
    i_verdict = src.index("= dimension_feasibility(")
    assert i_guard < i_verdict
    assert "continue" in src[i_guard:i_verdict]


def test_r10_measurement_failure_is_its_own_status():
    """Three different facts: we could not measure it, we measured and it is
    not enough, and the detector broke."""
    assert s31.EXAM_DIMENSION_MEASUREMENT_FAILED in s31.EXAM_UNKNOWN_STATUSES
    for other in (
        s31.EXAM_INSUFFICIENT_DIMENSION,
        s31.EXAM_RUNTIME_FAILED,
        s31.EXAM_TOO_SMALL,
        s31.EXAM_EXAMINED,
    ):
        assert s31.EXAM_DIMENSION_MEASUREMENT_FAILED != other


def test_r10_unmeasurable_samples_stay_na_and_do_not_dilute():
    obs = make_obs(
        [("p1", s31.EXAM_DIMENSION_MEASUREMENT_FAILED, None), ("p1", s31.EXAM_EXAMINED, True)]
    )
    s31.validate_examination_state(obs)
    r = s31.doublet_rate_per_patient(obs)
    assert r.loc["p1", "doublet_rate"] == 1.0
    assert r.loc["p1", "n_unknown"] == 1


@pytest.mark.parametrize("call", [True, False])
def test_r10_an_unmeasured_sample_may_not_carry_a_call(call):
    obs = make_obs([("p1", s31.EXAM_DIMENSION_MEASUREMENT_FAILED, call)])
    with pytest.raises(s31.ExaminationStateError):
        s31.validate_examination_state(obs)


def test_r10_unmeasured_samples_are_withheld_from_both_detection_paths():
    src = _src()
    assert "and sid not in measurement_failed" in src, "excluded from askable_idx"
    assert (
        src.count("if sid in measurement_failed:") == 2
    ), "the batched loop and the fallback must both honour it"
    # and they are never attempted, so the diagnostics assertion must allow it
    assert "EXAM_DIMENSION_MEASUREMENT_FAILED), (" in src


def test_r10_preflight_completion_requires_every_sample_measured():
    """Catching every per-sample error and then writing `completed` meant a run
    in which all 301 samples failed still exited 0 looking usable."""
    src = _src()
    assert '_MANIFEST["run_status"] = "completed" if _suff_all else "incomplete"' in src


def test_r10_sufficiency_is_reported_per_distribution():
    """The two distributions have different requirements, so one boolean cannot
    speak for both."""
    src = _src()
    seg = src[src.index('_MANIFEST["dimension_audit"] = {') : src.index("raise SystemExit(_exit)")]
    for field in (
        "n_samples_total",
        "n_samples_measured",
        "n_measurement_errors",
        "n_eligible_samples",
        "n_eligible_measured",
        "dimension_audit_sufficient_for_decision",
    ):
        assert f'"{field}"' in seg, f"the acceptance record must report {field}"
    assert '"all_samples_distribution"' in seg and '"eligible_only_distribution"' in seg


def test_r10_a_preflight_with_no_usable_distribution_says_so():
    src = _src()
    assert "NEITHER distribution is complete" in src


def test_r10_any_call_completed_is_named_for_what_it_means():
    """In the fallback the first sample that returns sets it while later
    samples may still fail, so it is not stage completion."""
    src = _src()
    assert '"scrublet_completed"' not in src, "the over-broad name is retired"
    assert '"scrublet_any_call_completed": False,' in src
    i = src.index('"scrublet_any_call_completed": False,')
    assert (
        "at least one" in src[max(0, i - 600) : i]
    ), "the manifest must say what the field actually covers"


# ==========================================================================
# Round 11 — the exit status, because SLURM reads nothing else
#
# Fixing run_status alone left two contradicting signals: the scheduler
# recording COMPLETED while the manifest said the run supported no decision.
# ==========================================================================


@pytest.mark.parametrize(
    "suff_all,suff_elig,expected",
    [
        (True, True, 0),  # both complete
        (True, False, 0),  # cannot actually occur, but must not be a failure
        (False, True, 0),  # supports an explicitly narrowed decision
        (False, False, 3),  # supports nothing
    ],
)
def test_r11_exit_code_is_zero_iff_some_distribution_is_usable(suff_all, suff_elig, expected):
    """Behavioural, over every combination -- not a search for a message."""
    assert s31.preflight_exit_code(suff_all, suff_elig) == expected


def test_r11_the_failure_code_is_distinct_from_the_guard_codes():
    """argparse uses 2 and the directory guards raise their own; 3 has to be
    readable as 'the audit ran and produced nothing usable'."""
    assert s31.PREFLIGHT_EXIT_OK == 0
    assert s31.PREFLIGHT_EXIT_INSUFFICIENT not in (0, 1, 2)


def test_r11_the_manifest_is_written_before_the_status_is_decided():
    """A non-zero exit must never cost the record of why."""
    src = _src()
    i_write = src.index('_MANIFEST["scrublet_any_call_completed"] = False')
    i_decide = src.index("_exit = preflight_exit_code(")
    assert i_write < i_decide
    assert "write_manifest()" in src[i_write:i_decide]


def test_r11_the_preflight_exits_with_the_computed_code():
    src = _src()
    assert "raise SystemExit(_exit)" in src
    assert "raise SystemExit(0)" not in src, "a hard-coded 0 defeats the whole point"


def test_r11_the_failure_message_names_where_to_look():
    src = _src()
    i = src.index("if _exit != PREFLIGHT_EXIT_OK:")
    seg = src[i : i + 1200]
    assert "file=sys.stderr" in seg, "a job failure belongs on stderr"
    assert "dimension_audit.csv" in seg and "MANIFEST_PATH" in seg
    assert "dimension_measurement_error" in seg


def test_r11_a_narrowed_decision_is_announced():
    """Exit 0 on eligible-only must not read as 'everything was fine'."""
    src = _src()
    assert "Only the ELIGIBLE-ONLY distribution is complete" in src
    assert "must state that it used that distribution" in src


def test_r11_attempted_is_not_reported_as_measured():
    """A measurement can raise, so 'every sample is measured' overstates it."""
    src = _src()
    # the counters keep the two apart
    assert '"n_samples_measured"' in src and '"n_measurement_errors"' in src
    i = src.index("_n_measured = len(_measured)")
    assert (
        '_measured = _audit[_audit["n_cells_effective"].notna()]' in src[:i]
    ), "measured means a value came back, not that it was tried"


# ---- durability claims must match the implementation ----------------------


def _calls_in(fn):
    """Fully-qualified call targets inside one function, from the AST."""
    return {_ast.unparse(n.func) for n in _ast.walk(fn) if isinstance(n, _ast.Call)}


def test_r11_write_manifest_replaces_and_does_not_fsync():
    """Bound to the FUNCTION, not to the file. The first version of this test
    used `'os.fsync(' in src`, so one unrelated fsync anywhere in the script
    would have licensed every durability claim in it."""
    fn = _fn("write_manifest")
    calls = _calls_in(fn)
    assert "os.replace" in calls, "the atomic rename is the whole guarantee"
    assert "os.fsync" not in calls, "this design deliberately does not fsync"


def test_r11_the_module_contains_no_fsync_call_at_all():
    """The chosen option was "do not implement fsync", so there is nothing in
    this module that could back a positive claim."""
    tree = _ast.parse(_src())
    calls = {_ast.unparse(n.func) for n in _ast.walk(tree) if isinstance(n, _ast.Call)}
    assert "os.fsync" not in calls


def test_r11_every_fsync_mention_denies_it_for_this_module():
    """A LINT-LEVEL TRIPWIRE, not a semantic proof.

    It is an English keyword scan; "without" and friends are broad and cannot
    establish that a sentence is really a denial. The tests that actually
    establish the facts are the two AST ones above -- write_manifest calls
    os.replace and not os.fsync, and the module contains no os.fsync call. This
    one exists only to catch a new positive-sounding sentence being added, and
    should not be described as more than that.

    It still beats its predecessor, which accepted any line containing
    "gate2_checkpoint" as a denial, so "the manifest is fsynced exactly like
    gate2_checkpoint" walked straight through the guard built to catch it.
    """
    DENIALS = (
        "not fsync",
        "does not fsync",
        "no fsync",
        "without",
        "deliberately does not",
        "claimed it did",
        "said it did",
    )
    offenders = []
    for n, line in enumerate(_src().splitlines(), 1):
        low = line.lower()
        if "fsync" not in low:
            continue
        if not any(d in low for d in DENIALS):
            offenders.append(f"{n}: {line.strip()}")
    assert not offenders, "positive fsync claim with nothing to back it:\n" + "\n".join(offenders)


def test_r11_the_docstring_states_the_non_guarantees_explicitly():
    doc = _ast.get_docstring(_fn("write_manifest"))
    assert "os.replace()" in doc, "the guarantee must be pinned to the replace returning"
    for absent in ("does not fsync", "power loss", "BEFORE the replace"):
        assert absent in doc, f"the docstring must say it does not cover: {absent}"


def test_r11_the_guarantee_is_scoped_to_the_replace_returning():
    """ "already-written" was not pinned to "os.replace has returned", which let
    it be read as covering an interruption before the replace."""
    doc = _ast.get_docstring(_fn("write_manifest"))
    assert "After `os.replace()` returns" in doc
    assert "signal, non-zero status or crash -- can lose an already-written one" not in doc


# ---- what a pre-replace failure actually leaves behind ---------------------
#
# SCOPE, stated so it is not read as more: these force os.replace to raise a
# CATCHABLE OSError. They do not simulate SIGKILL or a node failure, and cannot
# -- the finally block runs here and would not run there. That an uncatchable
# interruption may orphan a .inprogress file is an inference from the code, not
# something these tests observe.


def _manifest_at(monkeypatch, tmp_path, payload):
    monkeypatch.setattr(s31, "MANIFEST_PATH", str(tmp_path / "RUN_MANIFEST.json"), raising=False)
    monkeypatch.setattr(s31, "_MANIFEST", payload, raising=False)
    return tmp_path / "RUN_MANIFEST.json"


def test_r11_a_failure_before_replace_leaves_an_existing_manifest_untouched(monkeypatch, tmp_path):
    path = _manifest_at(monkeypatch, tmp_path, {"run_status": "started"})
    s31.write_manifest()
    before = path.read_text()

    monkeypatch.setattr(s31, "_MANIFEST", {"run_status": "completed"}, raising=False)
    monkeypatch.setattr(
        s31.os, "replace", lambda *a, **k: (_ for _ in ()).throw(OSError("disk full"))
    )
    with pytest.raises(OSError):
        s31.write_manifest()

    assert path.read_text() == before, "the previous final manifest must survive"
    assert not list(tmp_path.glob("*.inprogress")), "the catchable path cleans up"


def test_r11_a_failure_on_the_FIRST_write_leaves_no_final_manifest(monkeypatch, tmp_path):
    """The output guard requires a NEW EMPTY directory, so on the first write
    there is no previous manifest to fall back to. "leaves the previous
    manifest" is false here -- the real outcome is no final manifest at all.

    Catchable failure only; see the scope note above the helper.
    """
    path = _manifest_at(monkeypatch, tmp_path, {"run_status": "started"})
    assert not path.exists()

    monkeypatch.setattr(
        s31.os, "replace", lambda *a, **k: (_ for _ in ()).throw(OSError("node died"))
    )
    with pytest.raises(OSError):
        s31.write_manifest()

    assert not path.exists(), "no final manifest may be claimed to exist"
    assert not list(tmp_path.glob("*.inprogress"))


def test_r11_absence_of_a_guarantee_is_not_a_guarantee_of_loss():
    """The mirror image of the fsync overclaim: "nothing survives power loss"
    asserts certain loss, which is as unfounded as asserting durability.
    Without fsync the data may still survive -- OS writeback, the filesystem's
    behaviour, when the failure lands. What is absent is the guarantee."""
    doc = _ast.get_docstring(_fn("write_manifest"))
    assert "GUARANTEE" in doc or "guarantee" in doc
    for forbidden in ("nothing survives", "cannot survive", "will be lost"):
        assert forbidden not in doc.lower().replace("guarantee of loss", "")
    doc_md = open(
        os.path.join(os.path.dirname(CODE), "..", "..", "docs", "scrublet_dimension_parameter.md")
    ).read()
    assert "nothing here survives" not in doc_md
    assert "no guarantee" in doc_md


def test_r11_the_docstring_covers_the_first_write_case():
    doc = _ast.get_docstring(_fn("write_manifest"))
    # the claim must be conditional, not absolute
    assert "IF ONE" in doc and "EXISTS" in doc, "the previous-manifest claim must be conditional"
    assert "first write" in doc
    assert "inprogress" in doc, "an uncatchable interruption can leave the temp behind"
