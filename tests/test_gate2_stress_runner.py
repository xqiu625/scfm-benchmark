"""Runner draft: the design must reach the run, and failures must stay unknown.

SCOPE OF THESE TESTS -- read before trusting one.
-------------------------------------------------
Everything below the integration section uses a DOUBLE for scanpy. A double can
establish this runner's CONTROL FLOW: what it calls, in what order, with which
arguments, and what it records. It cannot establish scanpy's field names, dtypes
or numerics -- and in this repository it actively misled: an earlier double put
`_calls`/`_scores` into `uns`, keys scanpy never emits, and the runner was
written to match the double. Every test passed while a real successful run would
have been classified `calls_not_strict_boolean`.

The real schema is therefore pinned by `test_integration_*`, which runs the
installed scanpy on a small synthetic AnnData. Those are the only tests here
that prove anything about the producer-consumer interface.

Synthetic fixtures only. No cohort data, no HPCC path, no outcome column, no
Scrublet on real counts. scanpy is injected as a fake wherever the real one
would read a matrix, so these tests check the runner's control flow rather than
scanpy's numerics.
"""

import os
import sys

import numpy as np
import pytest

CODE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "scripts",
    "32_rebuild_upfront_qc",
    "code",
)
sys.path.insert(0, CODE)

import gate2_stress_design as D  # noqa: E402
import gate2_stress_runner as R  # noqa: E402

# ---- fakes ---------------------------------------------------------------


class _FakeObs(dict):
    """Just enough of a DataFrame: `columns` and item access."""

    @property
    def columns(self):
        return list(self.keys())


class FakeAD:
    """Minimal AnnData stand-in: enough shape and layers for the control flow."""

    def __init__(self, n_obs=40, n_vars=30, tag="obs"):
        rng = np.random.default_rng(0)
        self.X = rng.poisson(3, size=(n_obs, n_vars)).astype(np.int32)
        self.layers = {}
        self.var = {"highly_variable": np.ones(n_vars, dtype=bool)}
        self.uns = {}
        self.obs = _FakeObs()
        self.obs_names = [f"c{i}" for i in range(n_obs)]
        self.var_names = [f"g{i}" for i in range(n_vars)]
        self.tag = tag

    @property
    def n_obs(self):
        return self.X.shape[0]

    @property
    def n_vars(self):
        return self.X.shape[1]

    def copy(self):
        new = FakeAD(self.n_obs, self.n_vars, self.tag)
        new.X = self.X.copy()
        new.layers = dict(self.layers)
        new.var = dict(self.var)
        return new

    def __getitem__(self, key):
        return self.copy()


class FakeSC:
    """Records the seed the simulator was called with, and makes the simulated
    matrix depend on it -- so a seed that varied per rate would show up.

    It declares `gate2_test_identity` and deliberately has NO `__version__`, so
    every row produced under it records `scanpy_version=None` and
    `scanpy_identity_source='test_double'`. A fake cannot leave behind a row that
    claims a production library version.
    """

    gate2_test_identity = "tests.FakeSC"

    def __init__(self):
        self.sim_seeds = []
        self.sim_kwargs = []
        self.log1p_calls = []
        self.pp = self

    def filter_genes(self, ad, min_cells=None):
        pass

    def filter_cells(self, ad, min_genes=None):
        pass

    def normalize_total(self, ad, target_sum=None):
        pass

    def log1p(self, ad, layer=None):
        self.log1p_calls.append(layer)

    def highly_variable_genes(self, ad, layer=None):
        ad.var["highly_variable"] = np.ones(ad.n_vars, dtype=bool)

    def scrublet_simulate_doublets(self, ad, layer=None, random_seed=None, **kwargs):
        self.sim_seeds.append(random_seed)
        # Recorded as **kwargs, not as named defaults: a named default would
        # silently supply a value the ruling failed to pass.
        self.sim_kwargs.append(dict(kwargs))
        sim = FakeAD(ad.n_obs * 2, ad.n_vars, tag="sim")
        sim.X = (np.ones_like(sim.X) * (int(random_seed) % 97 + 1)).astype(np.int32)
        return sim


def make_scrublet(behaviour="ok", log=None):
    """A scrublet stand-in with selectable failure modes.

    It receives a REAL AnnData, because minimise_input builds one at the runner
    boundary, and writes to the REAL scanpy locations. It can still only prove
    control flow -- see the module docstring.
    """

    def fn(
        ad_obs, adata_sim=None, expected_doublet_rate=None, random_state=None, copy=True, **kwargs
    ):
        if log is not None:
            log.append(
                {
                    "random_state": random_state,
                    "sim_id": id(adata_sim),
                    "rate": expected_doublet_rate,
                    "kwargs": dict(kwargs),
                }
            )
        if behaviour == "raise":
            raise ValueError("boom")  # a DATA failure scanpy really raises
        if behaviour == "raise_bug":
            raise TypeError("unexpected keyword")  # a CODE defect: must propagate
        n = ad_obs.n_obs
        out = ad_obs.copy()
        rng = np.random.default_rng(int(random_state))
        scores = rng.random(n)
        calls = scores > (1.0 - float(expected_doublet_rate))
        if behaviour == "all_false":
            calls = np.zeros(n, dtype=bool)
        if behaviour == "bad_calls":
            calls = np.array(["x"] * n, dtype=object)
        if behaviour == "nan_scores":
            scores = np.full(n, np.nan)
        out.obs["predicted_doublet"] = calls
        out.obs["doublet_score"] = scores
        uns = {
            "threshold": 0.5,
            "doublet_scores_sim": rng.random(n * 2),
            "parameters": {
                "expected_doublet_rate": float(expected_doublet_rate),
                "random_state": int(random_state),
            },
        }
        if behaviour == "wrong_echo":
            uns["parameters"]["expected_doublet_rate"] = 0.999
        if behaviour in ("reordered_identity", "no_threshold_reordered"):
            out.obs_names = list(reversed(list(out.obs_names)))
        if behaviour == "no_threshold_reordered":
            uns.pop("threshold", None)  # `uns` is assigned to out AFTER this
        if behaviour == "no_threshold":
            uns.pop("threshold")
        if behaviour == "nan_threshold":
            uns["threshold"] = float("nan")
        out.uns["scrublet"] = uns
        return out

    return fn


#: A COMPLETE, EXPLICIT, TEST-ONLY ruling. Every value differs from scanpy's
#: default, so a test cannot pass by the ruling and the default coinciding.
#:
#: A previous draft passed `unadjudicated={}` instead -- the caller declaring
#: that nothing was open. That is the defect this structure exists to prevent,
#: so there is no longer any argument that can say it. The only way to reach the
#: scoring path is to supply a ruling covering EVERY required parameter.
TEST_RULING = {
    "sim_doublet_ratio": 3.0,  # scanpy default 2.0
    "synthetic_doublet_umi_subsampling": 0.75,  # scanpy default 1.0
    "log_transform": True,  # scanpy default False
    "n_prin_comps": 11,  # scanpy default 30
    "n_neighbors": 7,  # scanpy default None
    # A NumPy scalar on purpose: it lands in `scrublet_parameters` on EVERY
    # row including the failure paths, so a row that skipped the finalizer is
    # not natively serialisable and the test can see it.
    "stdev_doublet_rate": np.float64(0.03),  # scanpy default 0.02
    "knn_dist_metric": "manhattan",  # scanpy default 'euclidean'
    "normalize_variance": False,  # scanpy default True
    "mean_center": False,  # scanpy default True
    "use_approx_neighbors": False,  # scanpy default None
    "get_doublet_neighbor_parents": True,  # scanpy default False
}
assert set(TEST_RULING) == set(R.REQUIRED_PARAMETER_NAMES)


def run_one(behaviour="ok", cohort="qiu", channel="ch1", log=None, adata=None):
    sc = FakeSC()
    rows, side, scores = R.score_capture_across_grid(
        cohort,
        channel,
        adata if adata is not None else FakeAD(),
        scrublet_fn=make_scrublet(behaviour, log),
        sc=sc,
        ruling=TEST_RULING,
    )
    return rows, side, sc


# ---- the design must reach the run --------------------------------------


def test_no_constant_is_copied_from_the_design_module():
    """A copied grid can drift from the ruling, and a drifted copy looks like a
    deliberate change."""
    src = open(os.path.join(CODE, "gate2_stress_runner.py")).read()
    assert "0.005" not in src, "the grid is imported, never restated"
    assert "from gate2_stress_design import" in src


def test_the_grid_is_validated_against_the_design_before_use():
    sc = FakeSC()
    with pytest.raises(D.Gate2DesignError, match="frozen grid"):
        R.score_capture_across_grid(
            "qiu",
            "ch1",
            FakeAD(),
            scrublet_fn=make_scrublet(),
            sc=sc,
            grid=(0.01, 0.02),
            ruling=TEST_RULING,
        )


def test_every_row_carries_the_design_provenance():
    rows, _side, _sc = run_one()
    assert len(rows) == len(D.RESIDUAL_RATE_GRID)
    for r in rows:
        for k in (
            "design_version",
            "decision_reference",
            "grid_semantics",
            "use_class",
            "capture_seed",
            "expected_doublet_rate",
            "mechanical_status",
            "simulated_checksum",
            "never_pooled",
            "seed_stability",
        ):
            assert k in r, k
        assert r["design_version"] == D.DESIGN_VERSION
        assert "NOT a prevalence estimate" in r["grid_semantics"]
    assert [r["expected_doublet_rate"] for r in rows] == list(D.RESIDUAL_RATE_GRID)


def test_every_reliability_diagnostic_has_a_slot_on_every_row():
    rows, _side, _sc = run_one()
    for r in rows:
        for diag in D.RELIABILITY_DIAGNOSTICS_REQUIRED:
            assert diag in r, diag


# ---- not_run cohorts stop before any I/O --------------------------------


@pytest.mark.parametrize("cohort", ["stephenson", "wilk"])
def test_a_not_run_cohort_is_refused_before_the_matrix_is_touched(cohort):
    sc = FakeSC()
    with pytest.raises(D.Gate2DesignError, match="not_run"):
        R.score_capture_across_grid(
            cohort, "ch1", FakeAD(), scrublet_fn=make_scrublet(), sc=sc, ruling=TEST_RULING
        )
    assert sc.sim_seeds == [], "the capture was simulated before the refusal"


def test_an_unruled_cohort_is_refused_too():
    with pytest.raises(D.Gate2DesignError, match="no recorded use class"):
        R.score_capture_across_grid(
            "schulte_schrepping",
            "ch1",
            FakeAD(),
            scrublet_fn=make_scrublet(),
            sc=FakeSC(),
            ruling=TEST_RULING,
        )


# ---- common random numbers ----------------------------------------------


def test_the_capture_is_simulated_exactly_once_for_the_whole_grid():
    rows, _side, sc = run_one()
    assert len(sc.sim_seeds) == 1, f"simulated {len(sc.sim_seeds)} times, not once"
    assert sc.sim_seeds[0] == D.capture_seed("qiu", "ch1")


def test_the_simulated_checksum_is_identical_across_the_grid():
    """The empirical check, not just the interface: the simulated dataset is
    hashed and the hash must be the same at every grid point."""
    rows, _side, _sc = run_one()
    assert len({r["simulated_checksum"] for r in rows}) == 1
    assert R.grid_stability(rows, _side)["common_random_numbers_held"] is True


def test_the_same_random_state_and_the_same_sim_object_reach_every_grid_point():
    """Checksum equality is necessary but not sufficient: the simulated dataset
    can be built once while a per-rate random_state still perturbs the KNN and
    PCA at call time. A mutation that passed `seed + int(rate*1000)` survived
    the checksum test, so this inspects what the runner actually passed."""
    log = []
    rows, _side, _sc = run_one(log=log)
    assert len(log) == len(D.RESIDUAL_RATE_GRID)
    seed = D.capture_seed("qiu", "ch1")
    assert {e["random_state"] for e in log} == {seed}, "random_state varied across the grid"
    assert (
        len({e["sim_id"] for e in log}) == 1
    ), "a different simulated object reached some grid point"
    assert [e["rate"] for e in log] == list(D.RESIDUAL_RATE_GRID)


def test_two_captures_get_different_simulated_datasets():
    a, _sa, _ = run_one(channel="ch1")
    b, _sb, _ = run_one(channel="ch2")
    assert a[0]["simulated_checksum"] != b[0]["simulated_checksum"]


def test_grid_stability_reports_when_common_random_numbers_did_not_hold():
    rows, _side, _sc = run_one()
    rows[3]["simulated_checksum"] = "different"
    st = R.grid_stability(rows, _side)
    assert st["common_random_numbers_held"] is False
    assert "cannot be separated" in st["common_random_numbers_detail"]


# ---- mechanical failures stay unknown, per grid point -------------------


@pytest.mark.parametrize(
    "behaviour,status",
    [
        ("raise", "scrublet_raised"),
        ("no_threshold", "threshold_absent"),
        ("nan_threshold", "threshold_non_finite"),
        ("bad_calls", "calls_not_strict_boolean"),
        ("nan_scores", "scores_not_finite"),
    ],
)
def test_each_failure_mode_is_recorded_not_raised(behaviour, status):
    rows, _side, _sc = run_one(behaviour)
    assert {r["mechanical_status"] for r in rows} == {status}
    for r in rows:
        assert r["n_called"] is D.UNKNOWN_FILL_VALUE
        assert r["n_called"] is not False, "unknown was filled with False"


def test_an_all_false_prediction_is_examined_not_a_failure():
    """It can be legitimate; the threshold is what decides."""
    rows, _side, _sc = run_one("all_false")
    assert {r["mechanical_status"] for r in rows} == {"examined"}
    assert {r["n_called"] for r in rows} == {0}


def test_a_failure_at_one_grid_point_does_not_substitute_another_rate():
    rows, _side, _sc = run_one("no_threshold")
    assert [r["expected_doublet_rate"] for r in rows] == list(D.RESIDUAL_RATE_GRID)
    assert len(rows) == len(D.RESIDUAL_RATE_GRID)


def test_every_unknown_status_is_a_declared_mechanical_condition():
    for behaviour in ("raise", "no_threshold", "nan_threshold", "bad_calls", "nan_scores"):
        rows, _side, _sc = run_one(behaviour)
        for r in rows:
            assert r["mechanical_status"] in D.UNKNOWN_CONDITIONS


# ---- no best rate, no deletion ------------------------------------------


def test_grid_stability_offers_no_recommendation():
    rows, _side, _sc = run_one()
    st = R.grid_stability(rows, _side)
    joined = " ".join(k.lower() for k in st)
    for banned in ("best", "optimal", "recommended", "selected", "chosen", "argmax"):
        assert banned not in joined, f"{banned} appears as a field name"
    assert "does NOT indicate which rate" in st["interpretation_note"]


def test_the_runner_produces_no_deletion_mask_or_retained_matrix():
    src = open(os.path.join(CODE, "gate2_stress_runner.py")).read()
    for banned in (
        "deletion_mask",
        "retained_matrix",
        "write_h5ad",
        "to_csv",
        "remove_cells",
        "drop_cells",
    ):
        assert banned not in src, f"{banned} appears in a scoring-only runner"


def _runner_ast():
    import ast

    return ast.parse(open(os.path.join(CODE, "gate2_stress_runner.py")).read())


def _code_strings_and_names(tree):
    """String literals and identifiers from the CODE, with docstrings excluded.

    Three source-scanning tests in this repository have now failed or passed
    wrongly because they matched the module's own prose -- a docstring that says
    "no __main__" contains "__main__", and a sentence promising that no outcome
    column is loaded contains "outcome". The AST is what the runner does; the
    docstring is what it says about itself.
    """
    import ast

    docstrings = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            d = ast.get_docstring(node, clean=False)
            if d:
                docstrings.add(d)
    out = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            if node.value not in docstrings:
                out.append(node.value)
        elif isinstance(node, ast.Name):
            out.append(node.id)
        elif isinstance(node, ast.Attribute):
            out.append(node.attr)
    return [x.lower() for x in out]


def test_the_runner_has_no_entry_point_and_no_paths():
    """It cannot be pointed at cohort data as it stands."""
    import ast

    tree = _runner_ast()
    imports = {
        a.name.split(".")[0] for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names
    }
    imports |= {
        n.module.split(".")[0] for n in ast.walk(tree) if isinstance(n, ast.ImportFrom) and n.module
    }
    assert "argparse" not in imports, "a CLI would make this runnable"

    guards = [
        n
        for n in ast.walk(tree)
        if isinstance(n, ast.If)
        and any(isinstance(c, ast.Constant) and c.value == "__main__" for c in ast.walk(n.test))
    ]
    assert not guards, "an entry-point guard would make this runnable"

    # An absolute PATH, not merely a string starting with "/": the separator in
    # an f-string like f"{cohort}/{channel_id}: ..." is itself the literal "/".
    import re as _re

    for lit in _code_strings_and_names(tree):
        assert not _re.match(
            r"^/[A-Za-z0-9_.]", lit
        ), f"absolute path literal in the runner: {lit!r}"


def test_no_outcome_column_is_referenced_anywhere():
    """AST-level: docstrings may say the runner never touches outcome; the code
    must not name it."""
    for lit in _code_strings_and_names(_runner_ast()):
        for banned in ("mortality", "death", "outcome", "patient_list", "covid_severity"):
            assert banned not in lit, f"{banned!r} appears in runner code as {lit!r}"


def test_no_vendor_or_removal_efficiency_conversion_is_called():
    for lit in _code_strings_and_names(_runner_ast()):
        for banned in ("vendor", "removal_efficiency", "base_rate"):
            assert banned not in lit, f"{banned!r} appears in runner code as {lit!r}"


# ---- Ren is never pooled -------------------------------------------------


def test_ren_rows_are_flagged_never_pooled():
    rows, _side, _sc = run_one(cohort="ren")
    assert all(r["never_pooled"] for r in rows)
    assert all(r["use_class"] == D.USE_RESIDUAL_SCORING_ONLY for r in rows)


def test_aggregating_ren_with_another_cohort_raises():
    ren, _sr, _ = run_one(cohort="ren")
    qiu, _sq, _ = run_one(cohort="qiu")
    with pytest.raises(R.Gate2RunnerError, match="never pooled"):
        R.assert_not_pooled_with_others(ren + qiu)


def test_ren_alone_aggregates_fine():
    ren, _sr, _ = run_one(cohort="ren")
    assert R.assert_not_pooled_with_others(ren) == ["ren"]


def test_other_cohorts_aggregate_together():
    a, _sa, _ = run_one(cohort="qiu")
    b, _sb, _ = run_one(cohort="chongqing")
    assert R.assert_not_pooled_with_others(a + b) == ["chongqing", "qiu"]


# ---- seed stability is not faked ----------------------------------------


def test_seed_stability_is_reported_as_not_computed():
    """One capture seed cannot support a multi-seed stability estimate, and
    re-running the same seed would measure determinism instead."""
    rows, _side, _sc = run_one()
    assert all(r["seed_stability"] == R.SEED_STABILITY_STATUS for r in rows)
    assert "pending_secondary_seed_design" in R.SEED_STABILITY_STATUS
    assert R.grid_stability(rows, _side)["seed_stability"] == R.SEED_STABILITY_STATUS


def test_the_preprocessing_order_is_recorded_for_comparison():
    """A scanpy upgrade that reorders its own preprocessing should be a visible
    mismatch, not a silent divergence from what pp.scrublet does internally."""
    assert R.SCANPY_PREPROCESSING_ORDER[0].startswith("filter_genes")
    assert "highly_variable_genes" in " ".join(R.SCANPY_PREPROCESSING_ORDER)
    assert "scrublet_simulate_doublets" in " ".join(R.SCANPY_PREPROCESSING_ORDER)


# ==========================================================================
# Integration: the REAL scanpy, on a small synthetic AnnData.
#
# These are the only tests here that can establish the schema. Everything above
# uses a double and proves control flow only.
# ==========================================================================


# Module-level importorskip would skip THIS WHOLE FILE when scanpy is absent,
# taking every control-flow unit test with it. Scoped to the fixture instead.
@pytest.fixture(scope="module")
def sc():
    return pytest.importorskip("scanpy")


@pytest.fixture(scope="module")
def anndata_mod():
    return pytest.importorskip("anndata")


def _real_capture(anndata, n_obs=400, n_vars=300, seed=0):
    rng = np.random.default_rng(seed)
    X = rng.poisson(2, size=(n_obs, n_vars)).astype(np.float32)
    X[: n_obs // 10] += rng.poisson(6, size=(n_obs // 10, n_vars))
    a = anndata.AnnData(X=X)
    a.obs_names = [f"c{i}" for i in range(n_obs)]
    a.var_names = [f"g{i}" for i in range(n_vars)]
    return a


@pytest.mark.scanpy_integration
def test_integration_real_scanpy_emits_the_fields_the_runner_reads(sc, anndata_mod):
    """Pins the producer-consumer interface. If scanpy moves a field, this
    fails here rather than silently turning every success into an UNKNOWN."""
    out = sc.pp.scrublet(
        _real_capture(anndata_mod),
        expected_doublet_rate=0.08,
        random_state=7,
        copy=True,
        verbose=False,
    )
    assert R.SCANPY_OBS_CALLS in out.obs.columns
    assert R.SCANPY_OBS_SCORES in out.obs.columns
    assert out.obs[R.SCANPY_OBS_CALLS].dtype == bool
    uns = out.uns[R.SCANPY_UNS_KEY]
    assert R.SCANPY_UNS_THRESHOLD in uns
    assert R.SCANPY_UNS_SIM_SCORES in uns
    assert R.SCANPY_UNS_PARAMETERS in uns
    # and the keys the earlier double invented are absent
    assert "_calls" not in uns and "_scores" not in uns


@pytest.mark.scanpy_integration
def test_integration_mechanical_status_reads_a_real_success(sc, anndata_mod):
    out = sc.pp.scrublet(
        _real_capture(anndata_mod),
        expected_doublet_rate=0.08,
        random_state=7,
        copy=True,
        verbose=False,
    )
    status, calls, obs_scores, sim_scores, thr = R._mechanical_status(out, out.n_obs)
    assert status == "examined", status
    assert calls.dtype == bool and calls.shape == (out.n_obs,)
    assert obs_scores.shape == (out.n_obs,) and np.all(np.isfinite(obs_scores))
    assert sim_scores is not None and sim_scores.size > out.n_obs
    assert np.isfinite(thr)


@pytest.mark.scanpy_integration
def test_integration_observed_and_simulated_scores_are_not_the_same_vector(sc, anndata_mod):
    """Swapping them would still produce plausible diagnostics, so their
    distinctness is asserted rather than assumed."""
    out = sc.pp.scrublet(
        _real_capture(anndata_mod),
        expected_doublet_rate=0.08,
        random_state=7,
        copy=True,
        verbose=False,
    )
    _, _, obs_scores, sim_scores, _ = R._mechanical_status(out, out.n_obs)
    assert obs_scores.size != sim_scores.size
    assert obs_scores.size == out.n_obs
    assert sim_scores.size == out.n_obs * 2  # sim_doublet_ratio default 2.0


@pytest.mark.scanpy_integration
def test_integration_a_missing_threshold_is_detected_on_a_real_object(sc, anndata_mod):
    out = sc.pp.scrublet(
        _real_capture(anndata_mod),
        expected_doublet_rate=0.08,
        random_state=7,
        copy=True,
        verbose=False,
    )
    del out.uns[R.SCANPY_UNS_KEY][R.SCANPY_UNS_THRESHOLD]
    status, calls, _, sim, thr = R._mechanical_status(out, out.n_obs)
    assert status == "threshold_absent"
    assert calls is None and thr is None
    assert sim is not None, "the simulated scores survive a threshold failure"


@pytest.mark.scanpy_integration
def test_integration_reliability_block_is_actually_computed(sc, anndata_mod):
    """The five slots must hold values, not permanent None."""
    out = sc.pp.scrublet(
        _real_capture(anndata_mod),
        expected_doublet_rate=0.08,
        random_state=7,
        copy=True,
        verbose=False,
    )
    _, _, obs_scores, sim_scores, thr = R._mechanical_status(out, out.n_obs)
    block = R._reliability_block(obs_scores, sim_scores, thr)
    for slot in D.RELIABILITY_DIAGNOSTICS_REQUIRED:
        assert slot in block, slot
    assert block["observed_score_histogram"]["counts"] is not None
    assert block["simulated_score_histogram"]["counts"] is not None
    assert len(block["observed_score_histogram"]["percentiles_p0_to_p100"]) == 101
    assert block["tail_enrichment"]["observed_simulated_auc"] is not None
    assert block["seed_stability"] == R.SEED_STABILITY_STATUS


@pytest.mark.scanpy_integration
def test_integration_threshold_failure_keeps_BOTH_distributions(sc, anndata_mod):
    """Through the FULL path, not by feeding threshold=None to the helper.

    The earlier version of this test called _mechanical_status on a SUCCESSFUL
    object to obtain observed scores, then bypassed it and handed threshold=None
    straight to _reliability_block. That proved the helper tolerates a missing
    threshold; it could not see that _mechanical_status was returning
    observed=None whenever the threshold failed, which is what it did.
    """
    out = sc.pp.scrublet(
        _real_capture(anndata_mod),
        expected_doublet_rate=0.08,
        random_state=7,
        copy=True,
        verbose=False,
    )
    del out.uns[R.SCANPY_UNS_KEY][R.SCANPY_UNS_THRESHOLD]

    status, calls, obs_scores, sim_scores, thr = R._mechanical_status(out, out.n_obs)
    assert status == "threshold_absent"
    assert calls is None and thr is None
    assert obs_scores is not None, "observed scores were discarded by the threshold failure"
    assert sim_scores is not None

    block = R._reliability_block(obs_scores, sim_scores, thr)
    assert block["observed_score_histogram"]["counts"] is not None
    assert block["simulated_score_histogram"]["counts"] is not None
    assert len(block["observed_score_histogram"]["percentiles_p0_to_p100"]) == 101
    assert len(block["simulated_score_histogram"]["percentiles_p0_to_p100"]) == 101
    assert block["tail_enrichment"]["observed_simulated_auc"] is not None
    # threshold-derived fields are null, not zero
    summary = block["distribution_summary"]
    thr_fields = [k for k in summary if "threshold" in k and k != "diagnostics_schema_version"]
    assert thr_fields, "expected threshold-derived fields to exist as keys"
    assert all(
        summary[k] is None for k in thr_fields
    ), "a threshold-derived field was filled despite no threshold"


# ---- parameters are adjudicated, or nothing runs ------------------------


def test_production_deployment_hold_stops_before_any_matrix_work():
    """Methodological approval does not activate a production configuration."""
    sc = FakeSC()
    calls = []
    with pytest.raises(R.Gate2RunnerError, match="production configuration is not activated"):
        R.score_capture_across_grid(
            "qiu", "ch1", FakeAD(), scrublet_fn=lambda *a, **k: calls.append(1), sc=sc
        )  # production open set
    assert sc.sim_seeds == [], "the capture was prepared despite open parameters"
    assert calls == [], "scrublet was called despite open parameters"


def test_methodological_ruling_is_complete_without_activating_production():
    from gate2_parameter_ruling import METHODOLOGY_PARAMETER_RULING

    assert set(METHODOLOGY_PARAMETER_RULING) == R.REQUIRED_PARAMETER_NAMES
    assert R.PRODUCTION_PARAMETER_RULING == {}


def test_default_neighbor_rule_is_resolved_before_library_echo_comparison():
    from gate2_parameter_ruling import METHODOLOGY_PARAMETER_RULING

    ruling = dict(METHODOLOGY_PARAMETER_RULING)
    for n_obs, expected in [(36, 3), (100, 5), (400, 10)]:
        for rate in D.RESIDUAL_RATE_GRID:
            params = R.scrublet_call_parameters(rate, 123, ruling, n_obs=n_obs)
            assert params["n_neighbors"] == expected
            assert R._echo_agrees(params, {"n_neighbors": expected})["agree"]
            assert not R._echo_agrees(params, {"n_neighbors": expected + 1})["agree"]
    assert ruling["n_neighbors"] is None
    with pytest.raises(R.Gate2RunnerError, match="post-filter"):
        R.scrublet_call_parameters(0.04, 123, ruling)


def test_the_groups_partition_the_required_parameters():
    """Not a restatement of the module assertion: this fails if a later edit
    adds a required name without assigning it to an interface."""
    grouped = [k for names in R.PARAMETER_GROUPS.values() for k in names]
    assert sorted(grouped) == sorted(R.REQUIRED_PARAMETER_NAMES)
    assert len(grouped) == len(set(grouped)), "a parameter is in two groups"


def test_an_incomplete_ruling_is_refused_and_names_what_is_missing():
    partial = dict(TEST_RULING)
    del partial["knn_dist_metric"]
    with pytest.raises(R.Gate2RunnerError, match="knn_dist_metric"):
        R.validate_parameter_ruling(partial)


def test_a_ruling_with_an_unknown_parameter_is_refused():
    with pytest.raises(R.Gate2RunnerError, match="does not use"):
        R.validate_parameter_ruling(dict(TEST_RULING, invented_parameter=1))


def test_an_empty_ruling_cannot_declare_the_problem_away():
    """The exact shape of the previous defect: a caller asserting that nothing
    is open. There is now no argument that says that -- an empty ruling is an
    incomplete ruling."""
    sc = FakeSC()
    with pytest.raises(R.Gate2RunnerError, match="no methodological ruling"):
        R.score_capture_across_grid(
            "qiu", "ch1", FakeAD(), scrublet_fn=make_scrublet(), sc=sc, ruling={}
        )
    assert sc.sim_seeds == []


@pytest.mark.parametrize("dropped", sorted(TEST_RULING))
def test_dropping_any_single_parameter_stops_the_run_before_matrix_work(dropped):
    """The mutation the audit asked for, as a test: remove one key from the
    complete test ruling and the entry point must refuse before preparing."""
    sc = FakeSC()
    partial = {k: v for k, v in TEST_RULING.items() if k != dropped}
    with pytest.raises(R.Gate2RunnerError, match=dropped):
        R.score_capture_across_grid(
            "qiu", "ch1", FakeAD(), scrublet_fn=make_scrublet(), sc=sc, ruling=partial
        )
    assert sc.sim_seeds == [], "the capture was prepared with an open parameter"


def test_every_ruled_value_reaches_its_own_interface():
    """Approval that does not reach the call is not approval. The simulation
    values are compared against what the SIMULATOR received, not against the
    row -- the previous draft's row and simulator were free to disagree."""
    log = []
    sc = FakeSC()
    rows, _side, _sc = R.score_capture_across_grid(
        "qiu", "ch1", FakeAD(), scrublet_fn=make_scrublet(log=log), sc=sc, ruling=TEST_RULING
    )

    # simulation: what the simulator actually got
    assert sc.sim_kwargs, "the simulator was never called"
    for got in sc.sim_kwargs:
        assert got == {"sim_doublet_ratio": 3.0, "synthetic_doublet_umi_subsampling": 0.75}
    # preprocessing: log_transform=True means the observed and simulated objects
    # are both log1p'd with no layer, on top of the layer call
    assert sc.log1p_calls.count(None) == 2, sc.log1p_calls
    # scoring: what scrublet actually got
    for entry in log:
        for k in R.PARAMETER_GROUPS["scoring"]:
            assert entry["kwargs"][k] == TEST_RULING[k], k
    # and the row says the same thing as the call
    for row in rows:
        assert row["simulation_parameters"] == sc.sim_kwargs[0]
        assert row["preprocessing_parameters"] == {"log_transform": True}
        for k in R.PARAMETER_GROUPS["scoring"]:
            assert row["scrublet_parameters"][k] == TEST_RULING[k], k


# ---- a recorded check that changes nothing is not a check ----------------


def test_an_identity_mismatch_cannot_be_examined_or_enter_the_jaccard():
    """The Jaccard is positional. A previous draft recorded
    identity_verified=False and carried on to `examined`."""
    rows, side, _sc = run_one("reordered_identity")
    assert all(r["identity_verified"] is False for r in rows)
    assert {r["mechanical_status"] for r in rows} == {"identity_mismatch"}
    assert side == {}, "a misaligned grid point reached the Jaccard"
    assert all(r["n_called"] is D.UNKNOWN_FILL_VALUE for r in rows)
    D.assert_status_is_frozen("identity_mismatch")
    st = R.grid_stability(rows, side)
    assert st["n_examined"] == 0


def test_a_scanpy_level_failure_is_not_masked_by_an_identity_mismatch():
    """Both exclude the point; the more specific diagnosis is the one worth
    keeping, and identity_verified stays on the row either way."""
    rows, side, _sc = run_one("no_threshold_reordered")
    assert {r["mechanical_status"] for r in rows} == {"threshold_absent"}
    assert all(r["identity_verified"] is False for r in rows)
    assert side == {}


def test_a_parameter_echo_mismatch_cannot_be_examined():
    """A field that exists to verify execution has to be able to fail."""
    rows, side, _sc = run_one("wrong_echo")
    assert {r["mechanical_status"] for r in rows} == {"parameter_echo_mismatch"}
    assert all(r["parameters_echo_agree"]["agree"] is False for r in rows)
    assert all("expected_doublet_rate" in r["failure_detail"] for r in rows)
    assert side == {}
    D.assert_status_is_frozen("parameter_echo_mismatch")


def test_an_agreeing_echo_still_names_what_scanpy_does_not_report():
    rows, _side, _sc = run_one()
    echo = rows[0]["parameters_echo_agree"]
    assert echo["agree"] is True
    assert set(echo["compared"]) == {"expected_doublet_rate", "random_state"}
    # the eight scoring parameters are requested but not echoed back
    assert set(echo["not_echoed_by_scanpy"]) == set(R.PARAMETER_GROUPS["scoring"])


def test_the_recorded_parameters_are_the_ones_the_call_received():
    """Same object, not a parallel record. The first draft logged eight
    parameters while the call passed four unrelated ones."""
    log = []
    rows, _side, _sc = run_one(log=log)
    assert len(log) == len(D.RESIDUAL_RATE_GRID)
    for r, entry in zip(rows, log):
        p = r["scrublet_parameters"]
        assert p["expected_doublet_rate"] == entry["rate"]
        assert p["random_state"] == entry["random_state"]
        assert p["random_state"] == r["capture_seed"]
        for k, v in p.items():
            if k in ("expected_doublet_rate", "random_state"):
                continue
            assert entry["kwargs"].get(k) == v, f"{k} recorded but not passed"


# ---- call-set stability, not count stability ----------------------------


def test_stability_uses_call_sets_not_counts():
    """Two rates can call the same NUMBER of cells and share none of them."""
    rows, _side, _sc = run_one()
    st = R.grid_stability(rows, _side)
    assert "call_set_jaccard_pairwise" in st
    assert st["call_set_jaccard_min"] is not None
    n = st["n_examined"]
    assert len(st["call_set_jaccard_pairwise"]) == n * (n - 1) // 2


def test_equal_counts_with_disjoint_sets_are_reported_as_unstable():
    rows, side, _sc = run_one()
    n = rows[0]["n_cells"]
    a = np.zeros(n, dtype=bool)
    a[:5] = True
    b = np.zeros(n, dtype=bool)
    b[5:10] = True
    side = {rows[0]["expected_doublet_rate"]: a, rows[1]["expected_doublet_rate"]: b}
    rows[0]["n_called"], rows[1]["n_called"] = 5, 5
    st = R.grid_stability(rows[:2], side)
    assert st["n_called_span"] == 0, "counts agree"
    assert st["call_set_jaccard_min"] == 0.0, "sets are disjoint"


# ---- capture-preparation failure is not a per-grid failure --------------


def test_a_preparation_failure_yields_one_unknown_row_per_grid_point():
    class EmptySC(FakeSC):
        def filter_genes(self, ad, min_cells=None):
            ad.X = ad.X[:, :0]

    rows, _side, _sc = R.score_capture_across_grid(
        "qiu", "ch1", FakeAD(), scrublet_fn=make_scrublet(), sc=EmptySC(), ruling=TEST_RULING
    )
    assert len(rows) == len(D.RESIDUAL_RATE_GRID)
    assert {r["mechanical_status"] for r in rows} == {D.CAPTURE_PREPARATION_FAILED}
    assert all(r["n_called"] is D.UNKNOWN_FILL_VALUE for r in rows)
    assert all(r["failure_detail"] for r in rows), "the reason must be recorded"
    # The status is frozen in the design module, not invented here
    D.assert_status_is_frozen(rows[0]["mechanical_status"])


@pytest.mark.parametrize("exc_type", [ValueError, IndexError, ZeroDivisionError])
def test_a_library_data_exception_becomes_a_capture_level_unknown(exc_type):
    """Scanpy, numpy, anndata and scipy raise these for degenerate captures.
    Caught at the minimal call boundary and chained, with the STAGE recorded."""

    class RaisingSC(FakeSC):
        def filter_genes(self, ad, min_cells=None):
            raise exc_type("degenerate capture")

    rows, side, _sc = R.score_capture_across_grid(
        "qiu", "ch1", FakeAD(), scrublet_fn=make_scrublet(), sc=RaisingSC(), ruling=TEST_RULING
    )
    assert len(rows) == len(D.RESIDUAL_RATE_GRID)
    assert {r["mechanical_status"] for r in rows} == {D.CAPTURE_PREPARATION_FAILED}
    assert side == {}
    assert all(exc_type.__name__ in r["failure_detail"] for r in rows)
    assert all(
        "filter_genes" in r["failure_detail"] for r in rows
    ), "the failing stage must be recorded"


@pytest.mark.parametrize(
    "exc_type", [TypeError, AttributeError, NameError, KeyError, AssertionError]
)
def test_a_programming_error_propagates_instead_of_posing_as_data_quality(exc_type):
    """An unexpected keyword or a typo is a defect in the runner. A previous
    draft caught every Exception bar four, and a previous TEST asserted that an
    arbitrary TypeError was data degeneracy -- which is how a bug becomes a
    permanent finding about a cohort."""

    class BuggySC(FakeSC):
        def filter_genes(self, ad, min_cells=None):
            raise exc_type("this is a bug in the runner")

    with pytest.raises(exc_type):
        R.score_capture_across_grid(
            "qiu", "ch1", FakeAD(), scrublet_fn=make_scrublet(), sc=BuggySC(), ruling=TEST_RULING
        )


def test_a_programming_error_in_the_scoring_call_also_propagates():
    with pytest.raises(TypeError):
        R.score_capture_across_grid(
            "qiu",
            "ch1",
            FakeAD(),
            scrublet_fn=make_scrublet("raise_bug"),
            sc=FakeSC(),
            ruling=TEST_RULING,
        )


@pytest.mark.parametrize("exc_type", [KeyboardInterrupt, SystemExit, MemoryError])
def test_fatal_exceptions_are_never_swallowed(exc_type):
    """A design violation or an interrupt must not become a tidy
    'this capture could not be prepared'."""

    class FatalSC(FakeSC):
        def filter_genes(self, ad, min_cells=None):
            raise exc_type("stop")

    with pytest.raises(exc_type):
        R.score_capture_across_grid(
            "qiu", "ch1", FakeAD(), scrublet_fn=make_scrublet(), sc=FatalSC(), ruling=TEST_RULING
        )


def test_a_design_violation_during_preparation_propagates():
    class BadCohortSC(FakeSC):
        def filter_genes(self, ad, min_cells=None):
            raise D.Gate2DesignError("ruling violated")

    with pytest.raises(D.Gate2DesignError):
        R.score_capture_across_grid(
            "qiu",
            "ch1",
            FakeAD(),
            scrublet_fn=make_scrublet(),
            sc=BadCohortSC(),
            ruling=TEST_RULING,
        )


# ---- metadata minimisation ----------------------------------------------


def test_unapproved_obs_columns_do_not_reach_scanpy():
    """The runner cannot control what the CALLER loaded -- it is handed a
    constructed AnnData. What it controls is what it passes on. Fictional column
    names, and no value is read."""
    import anndata as _ad
    import pandas as _pd

    rng = np.random.default_rng(0)
    X = rng.poisson(3, size=(40, 30)).astype(np.float32)
    a = _ad.AnnData(X=X)
    a.obs_names = [f"c{i}" for i in range(40)]
    a.var_names = [f"g{i}" for i in range(30)]
    a.obs = _pd.DataFrame(
        {"fictional_clinical_label": ["x"] * 40, "fictional_followup_flag": [1] * 40},
        index=a.obs_names,
    )

    seen = {}

    def spy(ad_obs, **kwargs):
        seen["cols"] = list(ad_obs.obs.columns)
        return make_scrublet()(ad_obs, **kwargs)

    rows, _, _sc = R.score_capture_across_grid(
        "qiu", "ch1", a, scrublet_fn=spy, sc=FakeSC(), ruling=TEST_RULING
    )
    assert "fictional_clinical_label" not in seen["cols"]
    assert "fictional_followup_flag" not in seen["cols"]
    assert set(rows[0]["obs_columns_dropped"]) == {
        "fictional_clinical_label",
        "fictional_followup_flag",
    }
    assert rows[0]["obs_columns_retained"] == []


def _assert_native(obj, path="row"):
    """Every leaf is a built-in type -- no numpy scalars, no arrays."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            assert isinstance(k, str), f"{path}: non-string key {k!r}"
            _assert_native(v, f"{path}.{k}")
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            _assert_native(v, f"{path}[{i}]")
    else:
        # `type(obj) in (...)`, never isinstance: np.float64 subclasses float and
        # np.bool_ does not subclass bool, so isinstance would wave a NumPy scalar
        # straight through. That is the bug this check found in _to_native.
        assert obj is None or type(obj) in (
            str,
            bool,
            int,
            float,
        ), f"{path}: {type(obj).__name__} is not a native JSON type"


def test_the_row_round_trips_through_json_without_a_default_hook():
    """`json.dumps(rows, default=str)` -- what a previous test used -- stringifies
    anything it cannot serialise, so it would have passed on a row full of
    ndarrays and np.float64 and proved nothing about the checkpoint write."""
    import json as _json

    rows, side, _ = run_one()
    text = _json.dumps(rows)  # NO default= hook
    restored = _json.loads(text)
    assert restored == rows, "the row did not survive a JSON round trip unchanged"
    for r in rows:
        _assert_native(r)
        assert "_calls_vector" not in r
    assert all(isinstance(v, np.ndarray) for v in side.values())


def test_an_unconvertible_field_is_refused_rather_than_stringified():
    with pytest.raises(R.Gate2RunnerError, match="not JSON-serialisable"):
        R._to_native({"x": object()})


@pytest.mark.scanpy_integration
def test_integration_every_required_parameter_is_a_real_scanpy_parameter(sc):
    """A ruling over invented names would be a ruling over nothing. Measured
    against the installed signature, not read from documentation."""
    import inspect

    score_params = set(inspect.signature(sc.pp.scrublet).parameters)
    sim_params = set(inspect.signature(sc.pp.scrublet_simulate_doublets).parameters)

    for name in R.PARAMETER_GROUPS["scoring"]:
        assert name in score_params, f"{name} is not a sc.pp.scrublet parameter"
    for name in R.PARAMETER_GROUPS["simulation"]:
        assert name in sim_params, f"{name} is not a simulate_doublets parameter"

    # log_transform is grouped as PREPROCESSING, not scoring, even though
    # sc.pp.scrublet accepts it: scanpy applies it only on the `adata_sim is
    # None` branch, so passing it alongside adata_sim would be silently ignored
    # while the row claimed it had been applied. The runner applies it itself.
    assert "log_transform" in score_params, "scanpy does accept the kwarg"
    assert "log_transform" not in R.PARAMETER_GROUPS["scoring"]
    assert "log_transform" in R.PARAMETER_GROUPS["preprocessing"]


@pytest.mark.scanpy_integration
def test_integration_a_complete_ruling_drives_a_real_end_to_end_run(sc, anndata_mod):
    """The whole path against real scanpy: every ruled value is accepted by the
    real interfaces, the grid runs, and the rows are natively serialisable.

    Real scanpy raises TypeError on an unknown keyword, so this also fails if a
    ruling name drifts from the installed signature."""
    rng = np.random.default_rng(7)
    X = rng.poisson(2.0, size=(300, 200)).astype(np.float32)
    X[:, :20] += rng.poisson(8.0, size=(300, 20))  # some variable genes
    a = anndata_mod.AnnData(X=X)
    a.obs_names = [f"cell{i}" for i in range(300)]
    a.var_names = [f"gene{i}" for i in range(200)]
    a.obs["fictional_clinical_label"] = ["x"] * 300

    rows, side, _sc = R.score_capture_across_grid(
        "qiu", "ch1", a, scrublet_fn=sc.pp.scrublet, sc=sc, ruling=TEST_RULING
    )

    assert len(rows) == len(D.RESIDUAL_RATE_GRID)
    # "every status is frozen" would also pass if every point had FAILED, which
    # would prove nothing about the approved path. Real scanpy must actually
    # have accepted the ruling and run the grid to completion.
    statuses = {r["mechanical_status"] for r in rows}
    assert D.CAPTURE_PREPARATION_FAILED not in statuses, statuses
    assert "scrublet_raised" not in statuses, statuses
    assert "examined" in statuses, statuses
    for r in rows:
        D.assert_status_is_frozen(r["mechanical_status"])
        assert r["obs_columns_dropped"] == ["fictional_clinical_label"]
        assert r["obs_columns_retained"] == []
        assert r["scrublet_parameters"]["n_prin_comps"] == TEST_RULING["n_prin_comps"]
        assert r["simulation_parameters"]["sim_doublet_ratio"] == 3.0
        _assert_native(r)
    # one simulated dataset for the whole grid -- the common-random-numbers claim
    assert len({r["simulated_checksum"] for r in rows}) == 1
    import json as _json

    assert _json.loads(_json.dumps(rows)) == rows  # no default= hook
    assert set(side) <= {float(x) for x in D.RESIDUAL_RATE_GRID}


@pytest.mark.scanpy_integration
def test_integration_a_reordered_real_result_cannot_be_examined(sc, anndata_mod):
    """Identity is verified against the real object scanpy returns, not a fake."""
    rng = np.random.default_rng(11)
    X = rng.poisson(2.0, size=(300, 200)).astype(np.float32)
    X[:, :20] += rng.poisson(8.0, size=(300, 20))
    a = anndata_mod.AnnData(X=X)
    a.obs_names = [f"cell{i}" for i in range(300)]
    a.var_names = [f"gene{i}" for i in range(200)]

    def reordering(ad_obs, **kwargs):
        out = sc.pp.scrublet(ad_obs, **kwargs)
        out.obs = out.obs.iloc[::-1]  # same cells, different order
        return out

    rows, side, _sc = R.score_capture_across_grid(
        "qiu", "ch1", a, scrublet_fn=reordering, sc=sc, ruling=TEST_RULING
    )
    assert {r["mechanical_status"] for r in rows} == {"identity_mismatch"}
    assert _sc == {}, "misaligned scores must not be exported under input cell identities"
    assert side == {}, "a misaligned real result reached the Jaccard"


# ---- every exit satisfies the same contract -----------------------------


def test_every_exit_path_produces_the_same_schema():
    """Success, threshold failure, mismatch, scrublet exception and capture
    preparation failure are one schema, not two."""

    class EmptySC(FakeSC):
        def filter_genes(self, ad, min_cells=None):
            ad.X = ad.X[:, :0]

    batches = [
        run_one("ok")[0],
        run_one("no_threshold")[0],
        run_one("raise")[0],
        run_one("reordered_identity")[0],
        run_one("wrong_echo")[0],
        R.score_capture_across_grid(
            "qiu", "ch1", FakeAD(), scrublet_fn=make_scrublet(), sc=EmptySC(), ruling=TEST_RULING
        )[0],
    ]
    for rows in batches:
        for r in rows:
            assert set(r) == set(R.ROW_SCHEMA), set(r) ^ set(R.ROW_SCHEMA)


def test_the_scrublet_failure_path_is_not_exempt_from_the_json_contract():
    """A previous draft ran `_to_native` on the success path and appended the raw
    dict when scrublet raised -- the rows most likely to carry an odd value were
    exactly the ones skipping the conversion."""
    import json as _json

    for behaviour in ("ok", "raise", "no_threshold", "nan_scores", "nan_threshold"):
        rows, _side, _sc = run_one(behaviour)
        text = _json.dumps(rows, allow_nan=False)  # no default=, no NaN
        assert _json.loads(text) == rows, behaviour
        for r in rows:
            _assert_native(r)


def test_the_strict_dump_is_an_independent_backstop_to_the_converter():
    """`_to_native` maps NaN to None, so in normal operation the `allow_nan=False`
    dump can never fire. Tested here with the converter stubbed out, so the
    backstop is exercised on its own rather than assumed from the pair."""
    rows, _side, _sc = run_one()
    original = R._to_native
    R._to_native = lambda obj: obj  # the converter stops converting
    try:
        with pytest.raises(R.Gate2RunnerError, match="JSON contract"):
            R._finalize_row(dict(rows[0], threshold=float("nan")))
    finally:
        R._to_native = original


def test_a_non_finite_float_never_reaches_the_json_contract():
    """`allow_nan=False` is what a checkpoint writer would use; a NaN surfacing
    there would be found long after the run that produced it."""
    assert R._to_native(float("nan")) is None
    assert R._to_native(float("inf")) is None
    assert R._to_native(np.float64("nan")) is None


def test_a_row_missing_a_schema_field_is_refused():
    rows, _side, _sc = run_one()
    broken = dict(rows[0])
    del broken["scanpy_version"]
    with pytest.raises(R.Gate2RunnerError, match="missing required fields"):
        R._finalize_row(broken)


def test_a_row_with_an_undeclared_field_is_refused():
    """An added field that nothing declared is a schema drift, not a bonus."""
    rows, _side, _sc = run_one()
    with pytest.raises(R.Gate2RunnerError, match="undeclared fields"):
        R._finalize_row(dict(rows[0], invented_field=1))


# ---- the dependency that actually ran is recorded -----------------------


def test_a_module_that_declares_no_identity_is_refused_before_matrix_work():
    class Anonymous(FakeSC):
        gate2_test_identity = None

    sc = Anonymous()
    with pytest.raises(R.Gate2RunnerError, match="gate2_test_identity"):
        R.score_capture_across_grid(
            "qiu", "ch1", FakeAD(), scrublet_fn=make_scrublet(), sc=sc, ruling=TEST_RULING
        )
    assert sc.sim_seeds == []


def test_an_unapproved_scanpy_version_fails_closed():
    class WrongVersion(FakeSC):
        __version__ = "1.9.0"

    sc = WrongVersion()
    with pytest.raises(R.Gate2RunnerError, match="not in the approved set"):
        R.score_capture_across_grid(
            "qiu", "ch1", FakeAD(), scrublet_fn=make_scrublet(), sc=sc, ruling=TEST_RULING
        )
    assert sc.sim_seeds == [], "the capture was prepared on an unreviewed version"


def test_a_fake_cannot_leave_a_row_claiming_a_production_version():
    rows, _side, _sc = run_one()
    for r in rows:
        assert r["scanpy_identity_source"] == "test_double"
        assert r["scanpy_version"] is None
        assert r["preprocessing_contract_version"] == R.PREPROCESSING_CONTRACT_VERSION


def test_a_preparation_failure_row_still_carries_the_dependency():
    """The dependency is resolved before the matrix work, so a failure row knows
    which module failed."""

    class EmptySC(FakeSC):
        def filter_genes(self, ad, min_cells=None):
            ad.X = ad.X[:, :0]

    rows, _side, _sc = R.score_capture_across_grid(
        "qiu", "ch1", FakeAD(), scrublet_fn=make_scrublet(), sc=EmptySC(), ruling=TEST_RULING
    )
    for r in rows:
        assert r["scanpy_identity_source"] == "test_double"
        assert r["preprocessing_contract_version"] == R.PREPROCESSING_CONTRACT_VERSION


def test_the_preprocessing_order_tuple_is_documentation_not_detection():
    """It records THIS runner's expected order. Nothing reads scanpy's source, so
    it cannot notice an upstream change -- the version guard is what does.

    Checked over the AST, not the text: a text scan matches this module's own
    docstring, which cites `inspect.getsource` as the measurement behind a
    comment. That is the failure mode DEVELOPMENT.md rule 2 names, and it was hit
    again writing this very test."""
    import ast as _ast

    tree = _ast.parse(open(os.path.join(CODE, "gate2_stress_runner.py")).read())
    called = set()
    for node in _ast.walk(tree):
        if isinstance(node, _ast.Call):
            fn = node.func
            called.add(fn.attr if isinstance(fn, _ast.Attribute) else getattr(fn, "id", None))
    assert "getsource" not in called, "nothing here inspects scanpy source"
    imported = {
        a.name.split(".")[0]
        for node in _ast.walk(tree)
        if isinstance(node, _ast.Import)
        for a in node.names
    }
    imported |= {
        node.module.split(".")[0]
        for node in _ast.walk(tree)
        if isinstance(node, _ast.ImportFrom) and node.module
    }
    assert "inspect" not in imported
    assert R.APPROVED_SCANPY_VERSIONS, "the version guard is the actual detector"


@pytest.mark.scanpy_integration
def test_integration_the_recorded_version_is_the_running_scanpy(sc, anndata_mod):
    """Provenance resolved from the module that did the work."""
    rng = np.random.default_rng(3)
    X = rng.poisson(2.0, size=(300, 200)).astype(np.float32)
    X[:, :20] += rng.poisson(8.0, size=(300, 20))
    a = anndata_mod.AnnData(X=X)
    a.obs_names = [f"cell{i}" for i in range(300)]
    a.var_names = [f"gene{i}" for i in range(200)]

    rows, _side, _sc = R.score_capture_across_grid(
        "qiu", "ch1", a, scrublet_fn=sc.pp.scrublet, sc=sc, ruling=TEST_RULING
    )
    for r in rows:
        assert r["scanpy_identity_source"] == "module.__version__"
        assert r["scanpy_version"] == sc.__version__
        assert r["scanpy_version"] in R.APPROVED_SCANPY_VERSIONS
        assert r["preprocessing_contract_version"] == R.PREPROCESSING_CONTRACT_VERSION
    import json as _json

    assert _json.loads(_json.dumps(rows, allow_nan=False)) == rows
