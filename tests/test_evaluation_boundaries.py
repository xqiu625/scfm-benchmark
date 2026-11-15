"""Evaluation-boundary tests requested by the 2026-10-09 review: primary
LOCO, secondary validated-dependency-component CV, train-only fitted
preprocessing, and the "no out-of-sample transform" contract for jointly
fit Harmony correction.

`run_repeated_kfold`/`run_loco` (scripts/31_donor_level_revision/code/
repeated_cv.py) already existed and are already the project's real
evaluation harness -- this file is the first to test them directly, rather
than writing a parallel implementation. A caller's `fit_predict_fn(X_train,
y_train, X_test)` controls what is used for fitting; these tests prove the
harness actually enforces group integrity and that held-out rows cannot
leak into a fitted preprocessing step when the caller's own fit_predict_fn
is written correctly (and that a BADLY written one, fit on all rows, is
what a leakage check must be able to catch -- tested below too).
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

sys.path.insert(
    0,
    str(Path(__file__).parent.parent / "scripts/31_donor_level_revision/code"),
)
from repeated_cv import run_loco, run_repeated_kfold


def _synthetic_grouped_data(seed=0, n_groups=30, cells_per_group=6, n_features=10):
    rng = np.random.RandomState(seed)
    groups = np.repeat(np.arange(n_groups), cells_per_group)
    y = (groups % 2 == 0).astype(int)  # group-level label (like patient-level outcome)
    y = np.repeat((rng.rand(n_groups) < 0.4).astype(int), cells_per_group)
    X = rng.normal(size=(len(groups), n_features))
    X[:, 0] += y * 2.0  # a real, learnable signal
    return X, y, groups


def test_run_repeated_kfold_never_splits_a_group_across_train_and_test():
    """The exact property 'validated dependency components stay together'
    requires: StratifiedGroupKFold must place a whole group id in train XOR
    test, never both, within a single fold -- asserted here from what the
    REAL harness actually fed the production callback (via row identity
    carried through pandas .iloc), not from a separately re-constructed
    splitter instance that merely uses the same parameters (a 2026-10-09
    review's point: a harness bug that split groups could still pass a
    test that only checks an independently-built splitter)."""
    X, y, groups = _synthetic_grouped_data()
    X_df = pd.DataFrame(X)  # default RangeIndex -- .iloc[tr]/.iloc[te] carry it through
    violations = []
    n_calls = [0]

    def fit_predict_fn(X_train, y_train, X_test):
        n_calls[0] += 1
        train_groups = set(groups[X_train.index.to_numpy()])
        test_groups = set(groups[X_test.index.to_numpy()])
        if not train_groups.isdisjoint(test_groups):
            violations.append(train_groups & test_groups)
        clf = LogisticRegression().fit(X_train.to_numpy(), y_train)
        return {"p": clf.predict_proba(X_test.to_numpy())[:, 1]}

    result = run_repeated_kfold(X_df, y, groups, fit_predict_fn, seeds=[0, 1, 2])
    assert (
        n_calls[0] == 3 * 5
    ), f"expected 3 seeds x 5 folds = 15 real callback invocations, got {n_calls[0]}"
    assert (
        violations == []
    ), f"group(s) split across train/test in the ACTUAL harness calls: {violations}"
    assert set(result.keys()) == {0, 1, 2}
    for seed, preds in result.items():
        assert not np.isnan(
            preds["p"]
        ).any(), f"seed {seed}: every row must be held out exactly once"


def test_run_loco_holds_out_one_whole_cohort_per_fold():
    """Same instrumentation approach as above, applied to run_loco: asserts
    the ACTUAL train/test row identity fed to the callback leaves each
    cohort's rows entirely on one side."""
    X, y, groups = _synthetic_grouped_data(n_groups=20, cells_per_group=5)
    cohort = np.where(groups < 10, "CohortA", "CohortB")
    X_df = pd.DataFrame(X)
    violations = []
    n_calls = [0]

    def fit_predict_fn(X_train, y_train, X_test):
        n_calls[0] += 1
        train_cohorts = set(cohort[X_train.index.to_numpy()])
        test_cohorts = set(cohort[X_test.index.to_numpy()])
        if not train_cohorts.isdisjoint(test_cohorts):
            violations.append(train_cohorts & test_cohorts)
        clf = LogisticRegression().fit(X_train.to_numpy(), y_train)
        return {"p": clf.predict_proba(X_test.to_numpy())[:, 1]}

    oof, held_out = run_loco(X_df, y, cohort, fit_predict_fn)
    assert n_calls[0] == 2, f"expected exactly 2 LOCO folds (CohortA, CohortB), got {n_calls[0]}"
    assert (
        violations == []
    ), f"cohort split across train/test in the ACTUAL harness calls: {violations}"
    assert set(held_out) == {"CohortA", "CohortB"}
    assert not np.isnan(oof["p"]).any()


def _fixed_fold_perturbation_check(fit_fn_factory, seed=1):
    """Shared machinery for the two tests below. Runs `run_repeated_kfold`
    once on the original data, capturing (via pandas .iloc index identity)
    the FIRST fold's actual train/test row positions as the real harness
    produced them. Then builds X_mutated, which is BYTE-IDENTICAL to X on
    every row EXCEPT that fold's held-out (test) rows, and reruns with the
    same seed -- which reproduces the identical fold partition (the split
    depends only on y/groups/random_state, never on X's values, confirmed
    by the group-integrity test above). Returns (state_run1, state_run2)
    captured by `fit_fn_factory()` on the first call only, so a correct
    (train-only) implementation must return IDENTICAL state and a leaky one
    must not -- using the REAL harness and REAL row identity, not a
    reimplemented split (the 2026-10-09 review's objection to the prior
    version of this test, which only compared an untouched column)."""
    X, y, groups = _synthetic_grouped_data(seed=seed)
    X_df = pd.DataFrame(X)
    captured = {}
    held_out_idx = {}

    def make_capturing_fn(real_fit_fn):
        def fn(X_train, y_train, X_test):
            if "state" not in captured:
                captured["state"] = real_fit_fn(X_train, y_train, X_test)
                held_out_idx["test"] = X_test.index.to_numpy().copy()
                held_out_idx["train"] = X_train.index.to_numpy().copy()
            return {"p": np.zeros(len(X_test))}

        return fn

    run_repeated_kfold(X_df, y, groups, make_capturing_fn(fit_fn_factory), seeds=[seed])
    state_run1 = captured["state"]
    test_idx, train_idx = held_out_idx["test"], held_out_idx["train"]
    assert len(test_idx) > 0 and len(train_idx) > 0

    X_mutated = X.copy()
    X_mutated[test_idx, :] += 1000.0  # perturb ONLY this fold's real held-out rows
    assert np.array_equal(
        X_mutated[train_idx], X[train_idx]
    ), "test setup error: train rows must be untouched by this perturbation"
    X_mutated_df = pd.DataFrame(X_mutated)

    captured.clear()
    held_out_idx.clear()
    run_repeated_kfold(X_mutated_df, y, groups, make_capturing_fn(fit_fn_factory), seeds=[seed])
    state_run2 = captured["state"]
    assert np.array_equal(held_out_idx["train"], train_idx), (
        "test setup error: the fold partition changed between runs -- the premise "
        "that X's values don't influence the split no longer holds"
    )
    return state_run1, state_run2


def test_correctly_written_fit_predict_fn_is_unaffected_by_held_out_row_values():
    """Train-only fitting, proven against the REAL fold the harness drew and
    REAL held-out-only row perturbation (see `_fixed_fold_perturbation_check`
    docstring) -- a correct fit_predict_fn's captured fitted state must be
    EXACTLY unchanged when only held-out rows are perturbed."""

    def correct_fn(X_train, y_train, X_test):
        scaler = StandardScaler().fit(X_train.to_numpy())
        return scaler.mean_.copy()

    mean1, mean2 = _fixed_fold_perturbation_check(correct_fn)
    np.testing.assert_array_equal(
        mean1,
        mean2,
        err_msg=(
            "a correct train-only StandardScaler's fitted mean changed when ONLY "
            "this fold's held-out rows were perturbed -- held-out data influenced "
            "a quantity it must not have."
        ),
    )


def test_leakage_is_actually_detectable_if_fit_predict_fn_is_written_wrong():
    """Negative control using the SAME real-fold-perturbation machinery: a
    deliberately leaky fit_predict_fn (fits on train+test together) MUST
    show a different fitted mean under the identical held-out-only
    perturbation -- proving the invariance check above would actually catch
    this exact leakage bug through the real harness, not pass vacuously."""

    def leaky_fn(X_train, y_train, X_test):
        X_all = np.vstack([X_train.to_numpy(), X_test.to_numpy()])  # LEAK
        return StandardScaler().fit(X_all).mean_.copy()

    mean1, mean2 = _fixed_fold_perturbation_check(leaky_fn, seed=2)
    assert not np.allclose(mean1, mean2), (
        "expected the leaky (train+test) fitted mean to change when held-out "
        "rows were perturbed -- if it didn't, this perturbation check's "
        "ability to catch real leakage is itself broken."
    )


def test_harmony_correction_result_depends_only_on_its_own_call_not_prior_calls():
    """Behavioral enforcement of the 'no out-of-sample transform' contract:
    run_harmony_correction() must jointly fit and return corrected
    coordinates for exactly the cells given each call, with no state shared
    across calls. Proven here by calling it on dataset A, then on a
    DIFFERENT dataset B, then on dataset A again -- the two A results must
    be identical (same seed), which would NOT be guaranteed if any fitted
    object persisted across calls and silently influenced a later one on
    different data. This is the structural reason the function has no
    .fit()/.transform() split and must be labeled transductive/exploratory,
    never a reusable correction model applicable to an unseen cohort."""
    sys.path.insert(0, str(Path(__file__).parent.parent / "scripts/33_batch_effect_audit/code"))
    from run_batch_diagnostics import run_harmony_correction

    rng = np.random.RandomState(0)
    X_a = rng.normal(size=(60, 5)).astype(np.float64)
    batch_a = np.array(["A1", "A2"] * 30)
    X_b = rng.normal(size=(40, 5)).astype(np.float64) + 50.0  # very different data
    batch_b = np.array(["B1", "B2"] * 20)

    result_a1 = run_harmony_correction(X_a, batch_a, seed=0)
    _ = run_harmony_correction(X_b, batch_b, seed=0)  # different call, in between
    result_a2 = run_harmony_correction(X_a, batch_a, seed=0)

    np.testing.assert_allclose(
        result_a1,
        result_a2,
        err_msg=(
            "calling run_harmony_correction on unrelated data B between two "
            "identical calls on A changed A's result -- this would mean fitted "
            "state leaks across calls, contradicting the transductive-only, "
            "no-persisted-model contract this function is documented to have."
        ),
    )
    assert result_a1.shape == (60, X_a.shape[1]) or result_a1.shape[0] == 60
