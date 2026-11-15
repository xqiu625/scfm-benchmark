"""Mutation-killing tests for the single-capture preflight.

Each test exists because the 2026-09-24 preflight audit named that exact defect.
None of them runs the preflight against real data; they check the contract the
script must satisfy, at the level a mutation would break.
"""

import ast
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE = os.path.join(ROOT, "scripts", "32_rebuild_upfront_qc", "code")
AUDIT = os.path.join(ROOT, "reports", "manuscript", "revision", "gate2_capture_audit_2026-09-24")
PREFLIGHT = os.path.join(AUDIT, "capture_preflight_GSM5102900_HC1.py")
sys.path.insert(0, CODE)
sys.path.insert(0, AUDIT)

import gate2_parameter_ruling as RULING  # noqa: E402
import gate2_stress_runner as R  # noqa: E402

SRC = open(PREFLIGHT).read()
TREE = ast.parse(SRC)


def _calls(tree):
    out = set()
    for n in ast.walk(tree):
        if isinstance(n, ast.Call):
            f = n.func
            out.add(f.attr if isinstance(f, ast.Attribute) else getattr(f, "id", None))
    return out


CALLS = _calls(TREE)
IMPORTED = {a.name for n in ast.walk(TREE) if isinstance(n, ast.ImportFrom) for a in n.names}


# ---- 1. input identity is fail-closed ------------------------------------


def test_every_input_has_a_full_64_character_expected_sha():
    import capture_preflight_GSM5102900_HC1 as P  # noqa: E402

    assert len(P.EXPECTED_SHA256) == 3
    for name, digest in P.EXPECTED_SHA256.items():
        assert len(digest) == 64, f"{name}: {len(digest)}-char digest is not a full sha256"
        assert set(digest) <= set("0123456789abcdef"), name


def _guarded_comparisons(fn):
    """Every `if <a> != <b>: raise` in a function, as (left, right) source pairs.

    AST, not substrings: replacing `if got != EXPECTED_SHA256[name]:` with
    `if False:` leaves the word EXPECTED_SHA256 and a literal `!=` inside the
    f-string of the raise, so a text scan still passes while nothing is compared.
    Both mutations survived exactly that way."""
    out = []
    for node in ast.walk(fn):
        if not isinstance(node, ast.If):
            continue
        raises = any(isinstance(b, ast.Raise) for b in node.body)
        if not raises or not isinstance(node.test, ast.Compare):
            continue
        if not any(isinstance(o, ast.NotEq) for o in node.test.ops):
            continue
        out.append((ast.unparse(node.test.left), ast.unparse(node.test.comparators[0])))
    return out


def test_the_script_asserts_hashes_rather_than_only_recording_them():
    """A recomputed hash written into a manifest is an observation, not a check."""
    fn = next(
        n for n in ast.walk(TREE) if isinstance(n, ast.FunctionDef) and n.name == "assert_inputs"
    )
    cmps = _guarded_comparisons(fn)
    rights = [r for _, r in cmps]
    assert any(
        "EXPECTED_SHA256" in r for r in rights
    ), f"no guarded comparison against EXPECTED_SHA256; found {cmps}"
    assert any(
        "EXPECTED_BYTES" in r for r in rights
    ), f"no guarded comparison against EXPECTED_BYTES; found {cmps}"
    all_rights = [r for _, r in _guarded_comparisons(TREE)]
    assert any(
        "EXPECTED_CONTAINER_SHA256" in r for r in all_rights
    ), f"the container sha256 is not compared; found {all_rights}"
    assert any(
        "EXPECTED_CONTAINER_BYTES" in r for r in all_rights
    ), f"the container byte size is not compared; found {all_rights}"


def test_inputs_are_asserted_before_the_matrix_is_read():
    main = next(n for n in ast.walk(TREE) if isinstance(n, ast.FunctionDef) and n.name == "main")
    order = [n for n in ast.walk(main) if isinstance(n, ast.Call)]
    names = [getattr(n.func, "id", getattr(n.func, "attr", "")) for n in order]
    assert "assert_inputs" in names
    assert names.index("assert_inputs") < names.index(
        "mmread"
    ), "the matrix is read before input identity is asserted"


def test_the_container_identity_is_a_full_sha_not_path_plus_size():
    """The defect: `EXPECTED_CONTAINER_IDENTITY = '5027039977:'` with a comment
    promising the digest later, so only the size was ever compared."""
    import capture_preflight_GSM5102900_HC1 as P

    assert len(P.EXPECTED_CONTAINER_SHA256) == 64
    assert set(P.EXPECTED_CONTAINER_SHA256) <= set("0123456789abcdef")
    assert isinstance(P.EXPECTED_CONTAINER_BYTES, int)
    assert "EXPECTED_CONTAINER_IDENTITY" not in SRC, "the path+size stand-in is back"


# ---- 2. simulated matrix is normalised too -------------------------------


def test_both_observed_and_simulated_are_normalised():
    """The defect: normalising only the observed matrix, leaving the simulated one
    on a different scale before stacking."""
    norm_calls = [
        ast.unparse(n)
        for n in ast.walk(TREE)
        if isinstance(n, ast.Call) and getattr(n.func, "attr", "") == "normalize_total"
    ]
    targeted = [c for c in norm_calls if "target_sum" in c]
    assert len(targeted) == 2, f"expected 2 target_sum normalisations, got {targeted}"
    assert any("(a," in c or "(a)" in c for c in targeted), "observed not normalised"
    assert any("sim" in c for c in targeted), "simulated not normalised"


# ---- 3. observed-fit PCA, not a joint fit --------------------------------


def test_no_joint_svd_is_used():
    """scanpy fits TruncatedSVD on the OBSERVED matrix and transforms both. A joint
    SVD over stacked observed+simulated is a different pipeline.

    Checked over the AST, not the text: this module's own docstring explains the
    rejected approach and names `np.linalg.svd`, so a substring scan matches the
    explanation rather than the code. That is the failure mode DEVELOPMENT.md rule 2
    names, and it caught this test on the first run."""
    forbidden = {"svd", "vstack", "hstack"}
    for n in ast.walk(TREE):
        if isinstance(n, ast.Call):
            name = getattr(n.func, "attr", getattr(n.func, "id", ""))
            if name in forbidden:
                raise AssertionError(f"hand-rolled {name}() call is back: {ast.unparse(n)}")


def test_the_manifold_comes_from_scanpys_own_pipeline():
    assert "_scrublet" in SRC and "pipeline" in SRC
    assert any(
        c in CALLS for c in ("pca", "truncated_svd")
    ), "scanpy's own projection is not called"
    assert "observed-fit" in SRC


# ---- 4. the graph boundary is stated, not faked ---------------------------


def test_no_sklearn_nearest_neighbors_graph_is_built():
    """The previous draft built a graph with base k via sklearn and checksummed it as
    if it were the production graph."""
    assert "NearestNeighbors" not in SRC
    assert "kneighbors" not in SRC


def test_k_adj_is_recorded_with_the_scanpy_rule():
    assert "k_adj" in SRC
    assert "1 + n_sim / n_eff" in SRC or "1 + n_sim/n_eff" in SRC


def test_the_graph_checksum_is_null_with_a_stated_reason():
    assert "'graph_checksum': None" in SRC
    assert "graph_checksum_reason" in SRC


# ---- 5. seed comes from the repository -----------------------------------


def test_the_seed_is_derived_not_invented():
    assert "capture_seed" in IMPORTED, "capture_seed is not imported"
    assert "capture_seed(COHORT, CHANNEL)" in SRC
    assert "SEED = 2026" not in SRC, "a hand-written seed constant is back"


# ---- 6. the ruling is imported, not copied -------------------------------


def test_the_eleven_values_are_not_duplicated_in_the_preflight():
    """The defect: copying all eleven VALUES into the script, so two lists can drift.

    A parameter NAME may appear -- the manifest records what was used, and reading a
    value out of the imported ruling requires naming it. What must not appear is a
    dict literal that binds the names to values, i.e. a second copy of the ruling."""
    assert "approved_for" in IMPORTED
    names = set(RULING.METHODOLOGY_PARAMETER_RULING)
    for n in ast.walk(TREE):
        if isinstance(n, ast.Dict):
            keys = {
                k.value for k in n.keys if isinstance(k, ast.Constant) and isinstance(k.value, str)
            }
            overlap = keys & names
            if len(overlap) >= 3:
                raise AssertionError(
                    f"a dict literal re-binds {len(overlap)} ruling parameters "
                    f"({sorted(overlap)[:4]}...): the ruling has been copied"
                )


def test_the_ruling_version_is_checked():
    """Importing the checker is not calling it. Replacing the call with `pass` left
    the import in place and this test passed -- so the call itself is asserted."""
    assert "assert_ruling_version" in IMPORTED
    called = [
        ast.unparse(n)
        for n in ast.walk(TREE)
        if isinstance(n, ast.Call) and getattr(n.func, "id", "") == "assert_ruling_version"
    ]
    assert called, "assert_ruling_version is imported but never called"
    assert any("EXPECT_RULING_VERSION" in c for c in called), called
    with pytest.raises(SystemExit, match="ruling version mismatch"):
        RULING.assert_ruling_version("gate2-parameter-ruling/not-this-one")


def test_the_ruling_covers_exactly_the_required_parameters():
    got = RULING.approved_for(R.REQUIRED_PARAMETER_NAMES)
    assert set(got) == set(R.REQUIRED_PARAMETER_NAMES)
    with pytest.raises(SystemExit, match="does not cover"):
        RULING.approved_for(set(R.REQUIRED_PARAMETER_NAMES) | {"invented_parameter"})


def test_importing_the_ruling_does_not_unlock_the_runner():
    """The runner's refusal message is read from the runner, not hard-coded here:
    a previous version asserted an obsolete phrase and passed only because it ran
    against a stale worktree."""
    assert R.PRODUCTION_PARAMETER_RULING == {}
    with pytest.raises(R.Gate2RunnerError) as exc:
        R.validate_parameter_ruling(None)
    msg = str(exc.value)
    assert "not activated" in msg or "no methodological ruling" in msg, msg


# ---- 7. no artificial prior loop -----------------------------------------


def test_there_is_no_loop_over_the_prior_grid():
    """An unused `rho` loop shows only that the script reruns deterministically."""
    for n in ast.walk(TREE):
        if isinstance(n, ast.For):
            it = ast.unparse(n.iter)
            assert (
                "RESIDUAL_RATE_GRID" not in it and "PRIOR_GRID" not in it
            ), "the preflight loops over the prior grid again"


def test_prior_independence_is_claimed_structurally():
    assert "'manifold_built_once': True" in SRC
    assert "structural" in SRC
    assert "calculate_doublet_scores only" in SRC


# ---- 8. provenance and write safety --------------------------------------


def test_the_manifest_is_written_atomically():
    assert "def write_atomic" in SRC
    assert "os.replace" in SRC and "fsync" in SRC
    assert "json.dump(payload" in SRC


def test_the_manifest_records_provenance_and_a_lifecycle():
    for field in ("'status'", "dependency_identity", "'ruling'", "'seed'", "'inputs'", "versions"):
        assert field in SRC, f"manifest lacks {field}"
    assert "'status': 'incomplete'" in SRC and "'status'] = 'complete'" in SRC


def test_no_score_call_or_rate_is_produced():
    assert "calculate_doublet_scores" not in CALLS
    assert "call_doublets" not in CALLS
    assert "sc.pp.scrublet(" not in SRC
    assert "'produced_scores_calls_or_rates': False" in SRC


def test_no_outcome_column_is_read():
    """`covid_severity` occurs inside the container FILENAME
    (human_platelet_covid_severity_raw_counts.h5ad), so a substring scan flags a path
    that is never even opened here. What matters is whether an outcome column is
    used as a key, so string CONSTANTS are checked, with the known path excluded."""
    outcome = {
        "Death",
        "covid_severity",
        "mild_vs_severe",
        "severity_score",
        "progression_group",
        "three_group",
        "y",
    }
    for n in ast.walk(TREE):
        if isinstance(n, ast.Constant) and isinstance(n.value, str):
            v = n.value
            if v.endswith(".h5ad") or "/" in v:
                continue
            assert v not in outcome, f"outcome column named as a constant: {v!r}"
    assert "'outcome_columns_read': []" in SRC


def test_exactly_one_pre_simulation_raw_snapshot():
    """`layers['raw']` is what `scrublet_simulate_doublets` reads. Assigning it
    twice would silently change which matrix is simulated from -- the second
    snapshot would be taken after normalisation. Counted over the AST, so a
    duplicate cannot hide behind different whitespace."""
    assigns = []
    for n in ast.walk(TREE):
        if isinstance(n, ast.Assign):
            for t in n.targets:
                if (
                    isinstance(t, ast.Subscript)
                    and isinstance(t.value, ast.Attribute)
                    and t.value.attr == "layers"
                    and isinstance(t.slice, ast.Constant)
                    and t.slice.value == "raw"
                ):
                    assigns.append(ast.unparse(n))
    assert len(assigns) == 1, f"expected exactly one raw snapshot, found {assigns}"


def test_the_raw_snapshot_precedes_normalisation():
    """It must be taken before normalize_total, or the simulation reads
    normalised counts."""
    lines = SRC.splitlines()
    raw_at = next(i for i, line in enumerate(lines) if "layers['raw'] = " in line)
    norm_at = next(i for i, line in enumerate(lines) if "sc.pp.normalize_total(a)" in line)
    assert raw_at < norm_at, "the raw snapshot is taken after normalisation"
