"""Integration test for the scoring wrapper's dry-run control flow.

The wrapper runs under `set -e`. When the embedded validation exits 3 for an
absent permit, errexit would kill the shell before `rc=$?` could be read, making
the documented dry-run success path unreachable. This test runs the real wrapper
from a clean tracked fixture and proves the two exit codes differ.
"""

import hashlib
import os
import re
import shutil
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WRAPPER_REL = "reports/manuscript/revision/gate2_capture_audit_2026-09-24/" "run_capture_scoring.sh"
WRAPPER = os.path.join(ROOT, WRAPPER_REL)


def _fixture(tmp_path):
    """A clean tracked git repo holding the wrapper and the modules it validates."""
    if shutil.which("git") is None:
        pytest.skip("git unavailable")
    repo = tmp_path / "repo"
    files = [
        WRAPPER_REL,
        "scripts/32_rebuild_upfront_qc/code/gate2_scoring_config_GSM5102900_HC1.py",
        "scripts/32_rebuild_upfront_qc/code/gate2_score_capture_GSM5102900_HC1.py",
        "scripts/32_rebuild_upfront_qc/code/gate2_parameter_ruling.py",
        "scripts/32_rebuild_upfront_qc/code/gate2_stress_design.py",
        "scripts/32_rebuild_upfront_qc/code/gate2_stress_runner.py",
        "scripts/32_rebuild_upfront_qc/code/gate2_diagnostics.py",
    ]
    for rel in files:
        src = os.path.join(ROOT, rel)
        if not os.path.exists(src):
            pytest.skip(f"missing {rel}")
        dst = repo / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dst)
    # The wrapper validates bytes and hashes; CI cannot read the HPCC deposit.
    deposit = tmp_path / "deposit"
    deposit.mkdir()
    hashes = {}
    for suffix in ("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"):
        source = deposit / f"GSM5102900_HC1_{suffix}"
        source.write_bytes(f"synthetic {suffix}\n".encode())
        hashes[source.name] = hashlib.sha256(source.read_bytes()).hexdigest()
    config = repo / "scripts/32_rebuild_upfront_qc/code/gate2_scoring_config_GSM5102900_HC1.py"
    config.write_text(
        config.read_text() + f"\nDEPOSIT_DIR = {str(deposit)!r}\nINPUT_SHA256 = {hashes!r}\n"
    )
    subprocess.run(["git", "init", "-q"], cwd=repo, check=True)
    subprocess.run(["git", "add", "-A"], cwd=repo, check=True)
    subprocess.run(
        ["git", "-c", "user.email=t@t", "-c", "user.name=t", "commit", "-qm", "fixture"],
        cwd=repo,
        check=True,
    )
    return repo


def _run(repo, *args):
    env = dict(os.environ)
    env.pop("GATE2_PILOT_PERMIT", None)
    env.pop("GATE2_PERMIT_TEST_BYPASS", None)
    env["PATH"] = os.path.dirname(sys.executable) + os.pathsep + env["PATH"]
    return subprocess.run(
        ["bash", WRAPPER_REL, *args], cwd=repo, capture_output=True, text=True, env=env
    )


def test_dry_run_exits_zero_when_the_permit_is_absent(tmp_path):
    """The whole point of the fix: the validation's exit 3 must be caught, not
    kill the shell."""
    repo = _fixture(tmp_path)
    p = _run(repo, "--dry-run")
    assert p.returncode == 0, f"stdout:\n{p.stdout}\nstderr:\n{p.stderr}"
    assert "DRY RUN COMPLETE" in p.stdout
    assert "permit  : ABSENT" in p.stdout


def test_normal_mode_exits_nonzero_when_the_permit_is_absent(tmp_path):
    repo = _fixture(tmp_path)
    p = _run(repo)
    assert p.returncode != 0, f"stdout:\n{p.stdout}"
    assert p.returncode == 3, p.returncode
    assert "REFUSING to score without a valid scoped pilot permit" in p.stderr


@pytest.mark.parametrize("missing", [False, True])
def test_dry_run_refuses_changed_or_missing_input(tmp_path, missing):
    repo = _fixture(tmp_path)
    source = tmp_path / "deposit" / "GSM5102900_HC1_barcodes.tsv.gz"
    if missing:
        source.unlink()
    else:
        source.write_bytes(b"changed input")
    p = _run(repo, "--dry-run")
    assert p.returncode == 1, f"stdout:\n{p.stdout}\nstderr:\n{p.stderr}"
    expected = "input missing" if missing else "does not describe this file"
    assert expected in p.stderr
    assert "DRY RUN COMPLETE" not in p.stdout


def test_the_wrapper_refuses_modified_files(tmp_path):
    repo = _fixture(tmp_path)
    target = repo / "scripts/32_rebuild_upfront_qc/code/gate2_parameter_ruling.py"
    target.write_text(target.read_text() + "\n# local edit\n")
    p = _run(repo, "--dry-run")
    assert p.returncode == 2
    assert "unstaged changes" in p.stderr


def test_the_wrapper_guards_exactly_what_the_permit_digests(tmp_path):
    """The shell's guard list and IMPLEMENTATION_FILES must be the same set.

    They were not: the wrapper guarded six files while the config digested seven,
    so `gate2_diagnostics.py` could be edited locally without the wrapper's
    tracked/clean check ever looking at it.
    """
    sys.path.insert(0, os.path.join(ROOT, "scripts/32_rebuild_upfront_qc/code"))
    import gate2_scoring_config_GSM5102900_HC1 as C

    src = open(WRAPPER).read()
    guarded = set(re.findall(r'"(scripts/[^"]+\.py|reports/[^"]+\.sh)"', src))
    assert guarded == set(C.IMPLEMENTATION_FILES), (
        f"guarded not digested: {sorted(guarded - set(C.IMPLEMENTATION_FILES))}; "
        f"digested not guarded: {sorted(set(C.IMPLEMENTATION_FILES) - guarded)}"
    )


def test_a_divergent_guard_list_is_refused_by_the_running_wrapper(tmp_path):
    """Behavioural, not textual: drop a file from the shell list, commit it so the
    tracked/clean check passes, and the wrapper must refuse at the validation."""
    repo = _fixture(tmp_path)
    w = repo / WRAPPER_REL
    src = w.read_text()
    drop = '  "scripts/32_rebuild_upfront_qc/code/gate2_diagnostics.py"\n'
    assert drop in src, "anchor missing; the mutation would be a no-op"
    w.write_text(src.replace(drop, "", 1))
    subprocess.run(["git", "add", "-A"], cwd=repo, check=True)
    subprocess.run(
        [
            "git",
            "-c",
            "user.email=t@t",
            "-c",
            "user.name=t",
            "commit",
            "-qm",
            "drop one guarded file",
        ],
        cwd=repo,
        check=True,
    )
    p = _run(repo, "--dry-run")
    assert p.returncode == 4, f"rc={p.returncode}\nstdout:\n{p.stdout}\nstderr:\n{p.stderr}"
    assert "guards a different file set than the permit digests" in p.stderr
    assert "gate2_diagnostics.py" in p.stderr


def test_the_wrapper_refuses_staged_but_uncommitted_files(tmp_path):
    """`git diff` alone reports clean when a change is staged."""
    repo = _fixture(tmp_path)
    target = repo / "scripts/32_rebuild_upfront_qc/code/gate2_stress_design.py"
    target.write_text(target.read_text() + "\n# staged edit\n")
    subprocess.run(["git", "add", str(target)], cwd=repo, check=True)
    p = _run(repo, "--dry-run")
    assert p.returncode == 2
    assert "staged changes" in p.stderr


def test_errexit_is_suspended_only_around_the_validation():
    """Ordering checked over CODE lines only: the comment above the fix quotes
    `rc=$?` while explaining it, so scanning the raw text finds the explanation
    before the code and reports a false inversion. That is the same substring
    trap DEVELOPMENT.md rule 2 names, and it caught this test on first run."""
    code = [
        ln.strip()
        for ln in open(WRAPPER).read().splitlines()
        if ln.strip() and not ln.strip().startswith("#")
    ]
    assert "set +e" in code, "errexit is never suspended"
    off = code.index("set +e")
    rc = next(i for i, ln in enumerate(code) if ln.startswith("rc=$?"))
    on = next(i for i, ln in enumerate(code) if ln == "set -e" and i > off)
    assert off < rc < on, f"suspend/read/restore out of order: {off}, {rc}, {on}"


# ---- entrypoint output contract ------------------------------------------

import ast as _ast  # noqa: E402

ENTRY = os.path.join(
    ROOT, "scripts", "32_rebuild_upfront_qc", "code", "gate2_score_capture_GSM5102900_HC1.py"
)
ENTRY_SRC = open(ENTRY).read()
ENTRY_TREE = _ast.parse(ENTRY_SRC)


def _stmt_order(name_a, name_b):
    """Line numbers of the first call to each name, over the AST."""
    pos = {}
    for n in _ast.walk(ENTRY_TREE):
        if isinstance(n, _ast.Call):
            fn = getattr(n.func, "id", getattr(n.func, "attr", ""))
            if fn in (name_a, name_b) and fn not in pos:
                pos[fn] = n.lineno
    return pos


def test_sidecars_are_promoted_before_the_manifest_is_finalized():
    """A previous version wrote status=complete and THEN wrote a sidecar, so a
    failure on that write left a manifest asserting a completion that had not
    happened."""
    assert "_promote_sidecars" in ENTRY_SRC
    src = ENTRY_SRC
    promote = src.index("outputs = _promote_sidecars(")
    complete = src.index("status='complete'")
    assert promote < complete, "the manifest is finalized before sidecars are promoted"


def test_every_sidecar_is_written_to_a_temp_and_replaced():
    """Ordering compared over AST nodes, not over the unparsed text.

    The text form failed the moment the docstring explained the limitation: the
    docstring says "a sequence of `os.replace` calls", so `body.index("os.replace")`
    found the PROSE at offset 398 and reported the promotion as happening before
    the validation at 1698. Same substring trap as DEVELOPMENT.md rule 2.
    """
    fn = next(
        n
        for n in _ast.walk(ENTRY_TREE)
        if isinstance(n, _ast.FunctionDef) and n.name == "_promote_sidecars"
    )
    stmts = [n for n in _ast.walk(fn) if not isinstance(n, _ast.Constant)]

    replaces = [
        n.lineno
        for n in stmts
        if isinstance(n, _ast.Call)
        and getattr(n.func, "attr", None) == "replace"
        and getattr(getattr(n.func, "value", None), "id", None) == "os"
    ]
    fsyncs = [
        n.lineno
        for n in stmts
        if isinstance(n, _ast.Call) and getattr(n.func, "attr", None) == "fsync"
    ]
    validations = [
        n.lineno
        for n in stmts
        if isinstance(n, _ast.Raise) and "does not match its payload" in _ast.unparse(n)
    ]

    assert replaces, "no os.replace call"
    assert fsyncs, "no fsync call"
    assert validations, "no byte-for-byte validation"
    assert max(validations) < min(replaces), "temps must be validated before any is promoted"
    assert min(fsyncs) < min(replaces), "temps must be fsynced before promotion"


def test_per_cell_scores_are_written_not_inferred():
    """Aggregate summaries cannot support cell-level audit or threshold
    reconstruction, and scores cannot be recovered from calls."""
    assert "PER_PRIOR_SCORES.json" in ENTRY_SRC
    assert "score_sidecars" in ENTRY_SRC
    assert "never inferred from calls" in ENTRY_SRC


def test_the_manifest_lists_every_output_with_hashes_and_schema():
    assert "'outputs': " in ENTRY_SRC or "outputs=outputs" in ENTRY_SRC
    fn = next(
        n
        for n in _ast.walk(ENTRY_TREE)
        if isinstance(n, _ast.FunctionDef) and n.name == "_promote_sidecars"
    )
    body = _ast.unparse(fn)
    for field in ("'sha256'", "'bytes'", "'schema_version'"):
        assert field in body, field


# ---- graph acceptance: behavioural, not source-scanning -------------------
#
# The three tests that stood here scanned ENTRY_SRC for field names. All three
# passed while the capture recorded ZERO constructions, because a substring
# cannot tell you whether the instrumentation is attached to the right object.

sys.path.insert(0, os.path.join(ROOT, "scripts/32_rebuild_upfront_qc/code"))
import gate2_score_capture_GSM5102900_HC1 as E  # noqa: E402
import gate2_scoring_config_GSM5102900_HC1 as C2  # noqa: E402
import numpy as _np  # noqa: E402

K_ADJ = C2.RESOLVED["k_adj"]
SHAPE = C2.expected_graph_shape()
N_PRIORS = len(C2.RESIDUAL_RATE_GRID)


class _FakeCapture:
    """Stands in for _GraphCapture's digest() only -- see the scanpy_integration
    test for the real attachment, which is what a double cannot establish."""

    def __init__(self, records):
        self._records = records

    def digest(self):
        return self._records


def _rec(n=None, shape=None, idx_sha="a" * 64, dist_sha="b" * 64):
    return {
        "shape": list(shape or SHAPE),
        "n_neighbors_arg": K_ADJ if n is None else n,
        "n_queries": (shape or SHAPE)[0],
        "columns": (shape or SHAPE)[1],
        "indices_sha256": idx_sha,
        "distances_sha256": dist_sha,
        "self_first_column": False,
    }


def _evidence(records):
    return E._graph_evidence(_FakeCapture(records), N_PRIORS, K_ADJ, SHAPE)


def test_a_correct_capture_set_is_accepted():
    """Or every other test here proves only that the gate is stuck shut."""
    ev, problems = _evidence([_rec() for _ in range(N_PRIORS)])
    assert problems == [], problems
    assert ev["topology_identical_across_priors"] is True
    assert ev["n_graph_constructions"] == N_PRIORS


def test_zero_captures_is_a_failure_not_a_pass():
    """THE defect: with no captures the old field was None, and `is False` was
    never true, so a run that built no graph at all finalized complete."""
    ev, problems = _evidence([])
    assert problems, "zero graph constructions must not be acceptable"
    assert any("graph constructions" in p for p in problems)
    assert ev["topology_identical_across_priors"] is False


def test_a_single_capture_is_a_failure_not_a_pass():
    """With one capture the old identity check was trivially True."""
    _, problems = _evidence([_rec()])
    assert any("graph constructions" in p for p in problems)


def test_extra_constructions_are_a_failure():
    _, problems = _evidence([_rec() for _ in range(N_PRIORS + 1)])
    assert any("graph constructions" in p for p in problems)


def test_wrong_k_is_enforced_not_merely_recorded():
    """`k_matches_expected` was computed and never read."""
    recs = [_rec() for _ in range(N_PRIORS)]
    recs[3] = _rec(n=K_ADJ - 1)
    _, problems = _evidence(recs)
    assert any("n_neighbors argument" in p for p in problems), problems


def test_wrong_shape_is_enforced():
    recs = [_rec() for _ in range(N_PRIORS)]
    recs[5] = _rec(shape=(SHAPE[0] - 1, SHAPE[1]))
    _, problems = _evidence(recs)
    assert any("shape" in p for p in problems), problems


def test_an_empty_checksum_is_enforced():
    recs = [_rec() for _ in range(N_PRIORS)]
    recs[0] = _rec(idx_sha="")
    _, problems = _evidence(recs)
    assert any("empty indices checksum" in p for p in problems), problems


def test_topology_differing_across_priors_is_enforced():
    recs = [_rec() for _ in range(N_PRIORS)]
    recs[7] = _rec(idx_sha="c" * 64)
    _, problems = _evidence(recs)
    assert any("topology differed" in p for p in problems), problems


def test_distances_are_recorded_but_not_claimed_invariant():
    """Differing distances must NOT stop the run, because invariance is not
    claimed. The observation is still recorded, and the field says it is not
    enforced -- rather than a *_verified field that nothing branches on."""
    recs = [_rec() for _ in range(N_PRIORS)]
    recs[2] = _rec(dist_sha="d" * 64)
    ev, problems = _evidence(recs)
    assert problems == [], problems
    assert ev["distances_identical_across_priors_observed"] is False
    assert ev["distances_enforced"] is False


def test_the_run_stops_on_any_acceptance_problem():
    """The branch must be on the problem list, never on a tri-state field."""
    main = next(
        n for n in _ast.walk(ENTRY_TREE) if isinstance(n, _ast.FunctionDef) and n.name == "main"
    )
    src = _ast.unparse(main)
    assert "graph_problems" in src
    assert "graph_unacceptable" in src
    assert "identical_across_priors'] is False" not in src


def test_the_graph_capture_passes_values_through_untouched():
    """Instrumentation must observe, not alter, the graph scanpy uses."""
    fn = next(
        n for n in _ast.walk(ENTRY_TREE) if isinstance(n, _ast.FunctionDef) and n.name == "patched"
    )
    body = _ast.unparse(fn)
    assert body.rstrip().endswith("return (idx, dist)"), "the captured value is not returned as-is"


@pytest.mark.scanpy_integration
def test_the_capture_actually_attaches_to_the_pinned_scanpy_path():
    """The one test a double cannot replace.

    A prior version patched `sklearn.neighbors.NearestNeighbors.kneighbors` and
    captured nothing: scanpy reaches a `KNeighborsTransformer`, and
    `NearestNeighbors.kneighbors` IS the inherited mixin method, so assigning to
    the subclass shadows nothing the transformer resolves.

    This also DERIVES the reviewed expectations rather than restating them, so a
    scanpy change breaks this test instead of silently changing production.
    """
    import anndata as ad
    import scanpy as sc
    import scipy.sparse as sp

    rng = _np.random.default_rng(0)
    X = sp.csr_matrix(rng.poisson(0.4, size=(300, 200)).astype(_np.float32))
    a = ad.AnnData(X=X)

    per_prior = []
    for rate in (0.05, 0.10):
        b = a.copy()
        with E._GraphCapture() as gcap:
            sc.pp.scrublet(
                b, expected_doublet_rate=rate, n_prin_comps=10, random_state=0, verbose=False
            )
        per_prior.append(gcap.digest())

    for d in per_prior:
        assert len(d) == C2.GRAPH_CONSTRUCTIONS_PER_PRIOR, (
            f"captured {len(d)} constructions, expected "
            f"{C2.GRAPH_CONSTRUCTIONS_PER_PRIOR} per prior"
        )
        g = d[0]
        assert g["indices_sha256"] and g["distances_sha256"]
        # the columns rule production relies on, re-derived here
        assert g["columns"] == g["n_neighbors_arg"] - 1, C2.GRAPH_COLUMNS_RULE
        assert g["shape"][0] == 900, "300 observed + 600 simulated"

    # the manifold is built once per prior from the same seed, so topology holds
    assert per_prior[0][0]["indices_sha256"] == per_prior[1][0]["indices_sha256"]


@pytest.mark.scanpy_integration
def test_patching_nearest_neighbors_would_have_captured_nothing():
    """Pins the reason the old capture point was wrong, so it is not reinstated."""
    import anndata as ad
    import scanpy as sc
    import scipy.sparse as sp
    from sklearn.neighbors import NearestNeighbors

    seen = []
    orig = NearestNeighbors.kneighbors

    def patched(self, X=None, n_neighbors=None, return_distance=True):
        out = orig(self, X, n_neighbors, return_distance)
        seen.append(1)
        return out

    NearestNeighbors.kneighbors = patched
    try:
        rng = _np.random.default_rng(0)
        X = sp.csr_matrix(rng.poisson(0.4, size=(300, 200)).astype(_np.float32))
        b = ad.AnnData(X=X)
        sc.pp.scrublet(
            b, expected_doublet_rate=0.05, n_prin_comps=10, random_state=0, verbose=False
        )
    finally:
        NearestNeighbors.kneighbors = orig
    assert seen == [], (
        "NearestNeighbors.kneighbors now fires; the recorded reason for moving "
        "the capture point no longer holds and must be re-reviewed"
    )


def test_no_scoring_happens_without_the_permit_gate():
    main = next(
        n for n in _ast.walk(ENTRY_TREE) if isinstance(n, _ast.FunctionDef) and n.name == "main"
    )
    order = [
        n.lineno
        for n in _ast.walk(main)
        if isinstance(n, _ast.Call)
        and getattr(n.func, "attr", getattr(n.func, "id", "")) == "assert_permitted"
    ]
    load = [
        n.lineno
        for n in _ast.walk(main)
        if isinstance(n, _ast.Call) and getattr(n.func, "id", "") == "_load_deposit"
    ]
    assert order and load and min(order) < min(load), "the deposit is loaded before the permit gate"


def test_identity_mismatch_cannot_publish_a_complete_manifest(tmp_path, monkeypatch):
    """Exercise the entrypoint's exit and on-disk manifest, with scoring stubbed."""
    import contextlib
    import json
    import types

    out = tmp_path / "output"
    monkeypatch.setattr(sys, "argv", ["pilot", "--out-dir", str(out)])
    monkeypatch.setattr(C2, "assert_permitted", lambda *a: {"permit_id": "TEST-ONLY"})
    monkeypatch.setattr(C2, "scoring_config", lambda: {"output": {"directory": str(out)}})
    monkeypatch.setattr(C2, "config_sha256", lambda: "test-config")
    versions = dict(C2.DEPENDENCY_VERSIONS, python=sys.version.split()[0])
    monkeypatch.setattr(C2, "DEPENDENCY_VERSIONS", versions)
    for name in ("scanpy", "anndata", "numpy", "scipy", "sklearn"):
        monkeypatch.setitem(sys.modules, name, types.SimpleNamespace(__version__=versions[name]))
    monkeypatch.setattr(E, "_load_deposit", lambda: types.SimpleNamespace(shape=(7666, 1576)))
    rows = [
        {
            "expected_doublet_rate": rate,
            "n_cells": 7666,
            "identity_verified": False,
            "parameters_echo_agree": {"agree": True},
            "simulated_checksum": "same",
        }
        for rate in C2.RESIDUAL_RATE_GRID
    ]
    monkeypatch.setattr(E, "score_capture", lambda *a: (rows, {}, {}))

    @contextlib.contextmanager
    def graph():
        yield _FakeCapture([_rec() for _ in range(N_PRIORS)])

    monkeypatch.setattr(E, "_GraphCapture", graph)
    with pytest.raises(SystemExit, match="identity mismatch"):
        E.main()
    manifest = json.loads((out / "SCORING_MANIFEST.json").read_text())
    assert manifest["status"] == "stopped"
    assert manifest["stage"] == "drift_detected"
    assert not (out / "PER_PRIOR_SCORES.json").exists()
