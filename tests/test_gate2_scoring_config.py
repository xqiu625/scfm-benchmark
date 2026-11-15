"""The scoring config must fail closed on every substitution the audit named."""

import hashlib
import json
import os
import sys
import tempfile
import textwrap

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE = os.path.join(ROOT, "scripts", "32_rebuild_upfront_qc", "code")
sys.path.insert(0, CODE)

import gate2_scoring_config_GSM5102900_HC1 as C  # noqa: E402
import gate2_stress_runner as R  # noqa: E402


def _permit(**over):
    p = {
        "permit_id": "TEST-ONLY",
        "authorised_by": "test",
        "authorised_utc": "2026-09-24T00:00:00Z",
        "capture": C.CAPTURE,
        "config_sha256": C.config_sha256(),
        "input_sha256": dict(C.INPUT_SHA256),
        "implementation_sha256": C.implementation_sha256(),
        "implementation_commit": C._head_commit(),
        "output_root": C.OUTPUT_ROOT,
        "resource_ceiling": dict(C.REQUIRED_RESOURCE_CEILING),
        "scope": C.REQUIRED_SCOPE,
    }
    p.update(over)
    return p


def _ok(path):
    """Stub repository verifier. Injected through the signature, which production
    never does -- there is no environment flag that weakens the real check."""
    return "test-injected"


def _two_commit_repo(tmp_path, edit_implementation_in_b=False):
    """Commit A = the implementation. Commit B = the permit naming A.

    This is the real shape a production permit has to take, and the shape the old
    HEAD-equality check made impossible.
    """
    import shutil
    import subprocess

    if shutil.which("git") is None:
        pytest.skip("git unavailable")
    repo = tmp_path / "permitrepo"
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    for rel in C.IMPLEMENTATION_FILES:
        src = os.path.join(root, rel)
        if not os.path.exists(src):
            pytest.skip(f"missing {rel}")
        dst = repo / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dst)

    def git(*a, **kw):
        return subprocess.run(
            ["git", *a], cwd=repo, check=kw.pop("check", True), capture_output=True, text=True
        )

    git("init", "-q")
    git("add", "-A")
    git("-c", "user.email=t@t", "-c", "user.name=t", "commit", "-qm", "A: implementation")
    commit_a = git("rev-parse", "HEAD").stdout.strip()

    if edit_implementation_in_b:
        target = repo / C.IMPLEMENTATION_FILES[2]
        target.write_text(target.read_text() + "\n# edited after review\n")
    (repo / "permit.json").write_text(json.dumps({"implementation_commit": commit_a}))
    git("add", "-A")
    git("-c", "user.email=t@t", "-c", "user.name=t", "commit", "-qm", "B: permit")
    return repo, commit_a


def _write(permit):
    fh = tempfile.NamedTemporaryFile("w", suffix=".json", delete=False)
    json.dump(permit, fh)
    fh.close()
    return fh.name


# ---- the config is not runnable -----------------------------------------


def test_the_config_is_not_runnable():
    cfg = C.scoring_config()
    assert cfg["runnable"] is False
    assert cfg["permit"]["required"] is True and cfg["permit"]["present"] is False


def test_the_production_ruling_is_still_empty():
    assert R.PRODUCTION_PARAMETER_RULING == {}


# ---- missing permit ------------------------------------------------------


def test_scoring_is_refused_without_a_permit(monkeypatch):
    monkeypatch.delenv(C.PERMIT_PATH_ENV, raising=False)
    with pytest.raises(R.Gate2RunnerError, match="no scoped pilot permit"):
        C.assert_permitted()


def test_a_missing_permit_file_is_refused():
    with pytest.raises(R.Gate2RunnerError, match="permit not found"):
        C.assert_permitted("/nonexistent/permit.json")


def test_a_permit_missing_required_fields_is_refused():
    p = _permit()
    del p["authorised_by"]
    with pytest.raises(R.Gate2RunnerError, match="missing required field"):
        C.assert_permitted(_write(p), repo_verifier=_ok)


def test_a_well_formed_matching_permit_is_accepted(tmp_path, monkeypatch):
    """The gate must be openable by a correct permit, or the other tests prove
    only that it is stuck shut. The repository-artifact rule is bypassed with the
    explicit keyword, which production never passes.

    This runs against a two-commit fixture, not the live repo: the scoring
    entrypoint is not committed here, so `HEAD` is NOT a commit at which the
    implementation exists. Naming it would be the very thing the gate refuses --
    and it did refuse, which is how this test found its own premise was wrong.
    """
    repo, commit_a = _two_commit_repo(tmp_path)
    monkeypatch.setattr(C, "_repo_root", lambda: str(repo))
    p = _permit(implementation_commit=commit_a)
    assert C.assert_permitted(_write(p), repo_verifier=_ok)["permit_id"] == "TEST-ONLY"


def test_the_permit_commit_gate_is_satisfiable_at_all(tmp_path, monkeypatch):
    """Two commits: A holds the implementation, B adds the permit naming A.

    Requiring `deployed HEAD == implementation_commit` was a DEADLOCK. A permit
    must be tracked and identical to HEAD, so it can only arrive in a commit after
    the implementation; HEAD is then the permit commit, and a permit cannot carry
    the hash of the commit that contains it. No legitimate permit could pass.
    """
    repo, commit_a = _two_commit_repo(tmp_path)
    monkeypatch.setattr(C, "_repo_root", lambda: str(repo))
    assert C._assert_implementation_commit(commit_a)


def test_an_implementation_edited_after_the_reviewed_commit_is_refused(tmp_path, monkeypatch):
    """Ancestry alone is too weak: B could have edited the implementation."""
    repo, commit_a = _two_commit_repo(tmp_path, edit_implementation_in_b=True)
    monkeypatch.setattr(C, "_repo_root", lambda: str(repo))
    with pytest.raises(R.Gate2RunnerError, match="implementation changed between"):
        C._assert_implementation_commit(commit_a)


def test_a_commit_that_is_not_an_ancestor_is_refused(tmp_path, monkeypatch):
    repo, _ = _two_commit_repo(tmp_path)
    monkeypatch.setattr(C, "_repo_root", lambda: str(repo))
    import subprocess

    subprocess.run(["git", "checkout", "-q", "-b", "side", "HEAD~1"], cwd=repo, check=True)
    (repo / "unrelated.txt").write_text("x")
    subprocess.run(["git", "add", "-A"], cwd=repo, check=True)
    subprocess.run(
        ["git", "-c", "user.email=t@t", "-c", "user.name=t", "commit", "-qm", "side"],
        cwd=repo,
        check=True,
    )
    side = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=repo, capture_output=True, text=True
    ).stdout.strip()
    subprocess.run(["git", "checkout", "-q", "master"], cwd=repo, capture_output=True)
    subprocess.run(["git", "checkout", "-q", "main"], cwd=repo, capture_output=True)
    with pytest.raises(R.Gate2RunnerError, match="is not an ancestor"):
        C._assert_implementation_commit(side)


def test_an_implementation_dirty_on_disk_is_refused(tmp_path, monkeypatch):
    """Both commits can agree while the working tree does not.

    `implementation_sha256()` digests the files ON DISK, so without this check the
    digest could be satisfied by a local edit that neither commit contains. Found
    by mutation M11, which survived because nothing exercised this branch.
    """
    import subprocess

    repo, commit_a = _two_commit_repo(tmp_path)
    monkeypatch.setattr(C, "_repo_root", lambda: str(repo))
    assert C._assert_implementation_commit(commit_a)  # clean first
    target = repo / C.IMPLEMENTATION_FILES[4]
    target.write_text(target.read_text() + "\n# uncommitted local edit\n")
    with pytest.raises(R.Gate2RunnerError, match="differ from deployed HEAD on disk"):
        C._assert_implementation_commit(commit_a)
    # and staged-but-uncommitted counts too: `git diff HEAD` alone would miss it
    subprocess.run(["git", "add", str(target)], cwd=repo, check=True)
    with pytest.raises(R.Gate2RunnerError, match="differ from deployed HEAD on disk"):
        C._assert_implementation_commit(commit_a)


def test_a_commit_predating_the_implementation_is_refused(tmp_path, monkeypatch):
    """Naming an ancestor at which the implementation does not yet exist.

    This is not hypothetical: the live repo's HEAD (the preflight commit) is such
    a commit, and the first run of the new gate refused it -- which is how the
    positive tests were found to be asserting an impossible premise.
    """
    import shutil
    import subprocess

    if shutil.which("git") is None:
        pytest.skip("git unavailable")
    repo = tmp_path / "early"
    repo.mkdir()
    (repo / "unrelated.txt").write_text("before the implementation existed\n")

    def git(*a):
        return subprocess.run(["git", *a], cwd=repo, check=True, capture_output=True, text=True)

    git("init", "-q")
    git("add", "-A")
    git("-c", "user.email=t@t", "-c", "user.name=t", "commit", "-qm", "early")
    early = git("rev-parse", "HEAD").stdout.strip()
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    for rel in C.IMPLEMENTATION_FILES:
        dst = repo / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(os.path.join(root, rel), dst)
    git("add", "-A")
    git("-c", "user.email=t@t", "-c", "user.name=t", "commit", "-qm", "implementation")

    monkeypatch.setattr(C, "_repo_root", lambda: str(repo))
    with pytest.raises(R.Gate2RunnerError, match="predates the implementation"):
        C._assert_implementation_commit(early)


def test_an_unknown_commit_is_refused(tmp_path, monkeypatch):
    repo, _ = _two_commit_repo(tmp_path)
    monkeypatch.setattr(C, "_repo_root", lambda: str(repo))
    with pytest.raises(R.Gate2RunnerError, match="not a commit in this repository"):
        C._assert_implementation_commit("0" * 40)


def test_the_removed_bypass_variable_cannot_open_the_gate(monkeypatch):
    """The strongest form: an OTHERWISE-VALID permit must still be refused.

    Asserting that the variable is refused when everything else is broken would
    pass even if the branch honoured it, because the run would fail anyway.
    """
    monkeypatch.setenv(C.REMOVED_BYPASS_ENV, "1")
    with pytest.raises(R.Gate2RunnerError, match="bypass was removed"):
        C.assert_permitted(_write(_permit()), repo_verifier=_ok)


def test_the_removed_bypass_is_refused_even_when_set_to_a_falsy_string(monkeypatch):
    """Presence, not truth: `=0` is still an attempt to address the removed flag."""
    monkeypatch.setenv(C.REMOVED_BYPASS_ENV, "0")
    with pytest.raises(R.Gate2RunnerError, match="bypass was removed"):
        C.assert_permitted(_write(_permit()), repo_verifier=_ok)


def test_no_source_branch_skips_verification_on_the_bypass_variable():
    """Lint-level tripwire, not a semantic proof: the only mentions of the removed
    name must be the constant and the refusal message -- never a call to
    os.environ.get(...) whose truth causes an early `return`."""
    import ast
    import inspect

    src = inspect.getsource(C.assert_permitted)
    tree = ast.parse(textwrap.dedent(src))
    for node in ast.walk(tree):
        if isinstance(node, ast.If):
            body = getattr(node, "body", [])
            mentions = C.REMOVED_BYPASS_ENV in ast.dump(node.test) or any(
                isinstance(n, ast.Name) and n.id == "REMOVED_BYPASS_ENV"
                for n in ast.walk(node.test)
            )
            if mentions:
                assert any(
                    isinstance(s, ast.Raise) for s in body
                ), "a branch on the removed bypass must raise, never return early"


# ---- capture / input substitution ---------------------------------------


def test_a_permit_for_another_capture_is_refused():
    with pytest.raises(R.Gate2RunnerError, match="permit is for capture"):
        C.assert_permitted(_write(_permit(capture="GSM5102901_HC2")), repo_verifier=_ok)


def test_a_permit_with_substituted_inputs_is_refused():
    bad = dict(C.INPUT_SHA256)
    k = sorted(bad)[0]
    bad[k] = "0" * 64
    with pytest.raises(R.Gate2RunnerError, match="input_sha256 does not match"):
        C.assert_permitted(_write(_permit(input_sha256=bad)), repo_verifier=_ok)


# ---- config / parameter drift -------------------------------------------


def test_a_permit_for_a_different_config_is_refused():
    with pytest.raises(R.Gate2RunnerError, match="config_sha256"):
        C.assert_permitted(_write(_permit(config_sha256="0" * 64)), repo_verifier=_ok)


def test_parameter_drift_changes_the_config_hash():
    """A permit must not survive a parameter change: the hash is over the config,
    so altering a parameter invalidates it."""
    cfg = C.scoring_config()
    before = C.config_sha256(cfg)
    drifted = json.loads(json.dumps(cfg))
    drifted["parameters"]["n_prin_comps"] = 29
    assert C.config_sha256(drifted) != before


def test_resolved_values_match_the_preflight_measurements():
    assert C.RESOLVED["n_cells_effective"] == 7666
    assert C.RESOLVED["n_usable_features_hvg"] == 1576
    assert C.RESOLVED["n_neighbors"] == 44
    assert C.RESOLVED["k_adj"] == 132
    cfg = C.scoring_config()
    assert cfg["parameters"]["n_neighbors"] == 44


def test_the_eleven_parameters_come_from_the_ruling():
    cfg = C.scoring_config()
    import gate2_parameter_ruling as RULING

    for k, v in RULING.METHODOLOGY_PARAMETER_RULING.items():
        if k == "n_neighbors":
            continue
        assert cfg["parameters"][k] == v, k
    assert cfg["ruling"]["imported_not_copied"] is True


# ---- input hash drift on disk -------------------------------------------


def test_input_hash_drift_on_disk_fails_closed(tmp_path, monkeypatch):
    """If a deposit file changes, the C1-C3 evidence no longer describes it."""
    source = tmp_path / sorted(C.INPUT_SHA256)[0]
    source.write_bytes(b"original input")
    monkeypatch.setattr(C, "DEPOSIT_DIR", str(tmp_path))
    monkeypatch.setattr(
        C, "INPUT_SHA256", {source.name: hashlib.sha256(source.read_bytes()).hexdigest()}
    )
    assert C.assert_inputs_unchanged()
    source.write_bytes(b"changed input")
    with pytest.raises(R.Gate2RunnerError, match="does not describe this file"):
        C.assert_inputs_unchanged()


def test_a_missing_input_fails_closed(monkeypatch):
    monkeypatch.setattr(C, "DEPOSIT_DIR", "/nonexistent/dir")
    with pytest.raises(R.Gate2RunnerError, match="input missing"):
        C.assert_inputs_unchanged()


# ---- the graph limitation must stay explicit -----------------------------


def test_the_graph_limitation_is_present_and_false():
    cfg = C.scoring_config()
    assert cfg["graph_preflight_verified"] is False
    assert "calculate_doublet_scores" in cfg["graph_preflight_reason"]


def test_dropping_the_graph_reason_fails_closed(monkeypatch):
    """`false` without a stated reason reads as an oversight rather than a finding."""

    def no_reason():
        cfg = dict(original())
        cfg["graph_preflight_reason"] = ""
        return cfg

    original = C.scoring_config
    monkeypatch.setattr(C, "scoring_config", no_reason)
    with pytest.raises(R.Gate2RunnerError, match="graph limitation must be stated"):
        C.assert_permitted("/whatever")


# ---- outcome fields and output reuse ------------------------------------


def test_no_outcome_field_appears_in_the_config():
    blob = json.dumps(C.scoring_config())
    for bad in C.FORBIDDEN_OUTPUT_FIELDS:
        assert f'"{bad}"' not in blob, bad


def test_output_reuse_is_forbidden():
    cfg = C.scoring_config()
    assert cfg["output"]["new_directory_required"] is True
    assert cfg["output"]["reuse_forbidden"] is True


def test_the_input_population_is_the_deposit_not_the_container_subset():
    cfg = C.scoring_config()
    assert "7,666" in cfg["input"]["population"]
    assert "65-cell" in cfg["input"]["population"]


# ---- permit authenticity: a tmp file is not a permit ---------------------


def test_a_permit_outside_the_repository_is_refused_without_the_bypass(monkeypatch):
    """PRODUCTION calls assert_permitted() with no verifier argument, so it always
    gets the real check. $GATE2_PILOT_PERMIT pointing at any writable path would
    otherwise be a permit factory."""
    with pytest.raises(R.Gate2RunnerError, match="outside the deployed repository"):
        C.assert_permitted(_write(_permit()))


def test_an_untracked_permit_inside_the_repository_is_refused(tmp_path):
    root = C._repo_root()
    p = os.path.join(root, "_test_untracked_permit.json")
    try:
        with open(p, "w") as fh:
            json.dump(_permit(), fh)
        with pytest.raises(R.Gate2RunnerError, match="not tracked"):
            C.assert_permitted(p)
    finally:
        os.remove(p)


# ---- implementation, output root, resource ceiling -----------------------


def test_a_permit_for_a_different_implementation_is_refused(monkeypatch):
    with pytest.raises(R.Gate2RunnerError, match="implementation_sha256"):
        C.assert_permitted(_write(_permit(implementation_sha256="0" * 64)), repo_verifier=_ok)


def test_a_permit_with_a_different_output_root_is_refused(monkeypatch):
    with pytest.raises(R.Gate2RunnerError, match="output_root"):
        C.assert_permitted(_write(_permit(output_root="/tmp/elsewhere")), repo_verifier=_ok)


def test_a_permit_without_a_resource_ceiling_is_refused(monkeypatch):
    with pytest.raises(R.Gate2RunnerError, match="resource_ceiling"):
        C.assert_permitted(_write(_permit(resource_ceiling=None)), repo_verifier=_ok)


def test_editing_any_implementation_file_changes_the_digest(tmp_path):
    """The digest is over file contents, so a permit cannot survive a code edit."""
    before = C.implementation_sha256()
    root = C._repo_root()
    target = os.path.join(root, C.IMPLEMENTATION_FILES[0])
    original = open(target, "rb").read()
    try:
        with open(target, "ab") as fh:
            fh.write(b"\n# mutated\n")
        assert C.implementation_sha256() != before
    finally:
        with open(target, "wb") as fh:
            fh.write(original)
    assert C.implementation_sha256() == before


# ---- output directory is derived, not merely "required" ------------------


def test_the_output_directory_is_derived_and_under_the_bound_root():
    cfg = C.scoring_config()
    d = cfg["output"]["directory"]
    assert d.startswith(C.OUTPUT_ROOT + os.sep)
    assert C.CAPTURE in d and "__cfg" in d and "__impl" in d


def test_the_graph_acceptance_criteria_are_bound_by_the_config_hash():
    """Relaxing an acceptance criterion must invalidate a permit.

    The criteria are what makes the graph record evidence rather than decoration;
    if they sat outside the hash, a permit would survive someone lowering them.
    """
    cfg = C.scoring_config()
    ga = cfg["graph_acceptance"]
    assert ga["constructions_per_prior"] == C.GRAPH_CONSTRUCTIONS_PER_PRIOR
    assert ga["expected_shape"] == list(C.expected_graph_shape())
    before = C.config_sha256()
    body = C._config_body()
    body["graph_acceptance"] = dict(body["graph_acceptance"], constructions_per_prior=0)
    assert C.config_sha256(body) != before


def test_there_is_exactly_one_config_hash():
    """The body and the full config must hash alike.

    `output_directory()` has to hash the body (hashing the full config would
    recurse through it), while the permit check hashes the full config. If those
    two differ, the permit binds one value and the directory displays another.
    """
    assert C.config_sha256() == C.config_sha256(C._config_body())


def test_the_directory_cfg_segment_IS_the_published_config_sha256():
    """Not merely that a `__cfg` segment exists -- that it is THE bound hash.

    An auditor holding a permit must be able to reproduce the output path from
    its `config_sha256`. Asserting only that the substring `__cfg` appears is a
    check whose result cannot fail for the defect it is supposed to catch.
    """
    d = os.path.basename(C.output_directory())
    cfg_seg = d.split("__cfg")[1].split("__impl")[0]
    assert cfg_seg == C.config_sha256()[:12]
    impl_seg = d.split("__impl")[1]
    assert impl_seg == C.implementation_sha256()[:12]


def test_excluding_the_directory_from_the_hash_does_not_excuse_the_rest_of_output():
    """Only `directory` is excluded; the rest of the output block is still bound."""
    before = C.config_sha256()
    body = C._config_body()
    body["output"] = dict(body["output"], reuse_forbidden=False)
    assert C.config_sha256(body) != before


def test_a_config_change_changes_the_output_directory(monkeypatch):
    before = C.output_directory()
    before_hash = C.config_sha256()
    monkeypatch.setattr(C, "PRIOR_LABEL", "something_else")
    assert C.output_directory() != before
    # the two must move together, or the directory stops identifying the config
    assert C.config_sha256() != before_hash
    assert C.config_sha256()[:12] in os.path.basename(C.output_directory())


# ---- the scoped ruling is never the global one ---------------------------


def test_the_scoped_ruling_requires_a_permit():
    with pytest.raises(R.Gate2RunnerError, match="without a validated permit"):
        C.scoped_ruling(None)


def test_the_scoped_ruling_does_not_activate_the_global_one(tmp_path, monkeypatch):
    repo, commit_a = _two_commit_repo(tmp_path)
    monkeypatch.setattr(C, "_repo_root", lambda: str(repo))
    permit = C.assert_permitted(_write(_permit(implementation_commit=commit_a)), repo_verifier=_ok)
    scoped = C.scoped_ruling(permit)
    import gate2_parameter_ruling as RULING

    assert scoped == dict(RULING.METHODOLOGY_PARAMETER_RULING)
    assert R.PRODUCTION_PARAMETER_RULING == {}


def test_the_preflight_commit_is_not_called_a_deploy_commit():
    """9aab7be predates the scoring entrypoint; calling it the deployment commit
    would claim the scoring code was deployed then."""
    assert hasattr(C, "PREFLIGHT_COMMIT")
    assert not hasattr(C, "DEPLOY_COMMIT")
    assert "preflight_commit" in C.scoring_config()["provenance"]


# ---- scope and ceiling are exact values, not merely present --------------


def test_a_permit_with_a_vague_scope_is_refused():
    with pytest.raises(R.Gate2RunnerError, match="the only authorised scope"):
        C.assert_permitted(_write(_permit(scope="scoring")), repo_verifier=_ok)


def test_a_permit_with_a_larger_resource_ceiling_is_refused():
    bigger = dict(C.REQUIRED_RESOURCE_CEILING, mem="96G")
    with pytest.raises(R.Gate2RunnerError, match="resource_ceiling"):
        C.assert_permitted(_write(_permit(resource_ceiling=bigger)), repo_verifier=_ok)


def test_a_permit_for_another_commit_is_refused():
    with pytest.raises(R.Gate2RunnerError, match="implementation_commit"):
        C.assert_permitted(_write(_permit(implementation_commit="0" * 40)), repo_verifier=_ok)


def test_the_implementation_digest_covers_the_transitive_diagnostics_module():
    """gate2_stress_runner imports gate2_diagnostics; omitting it let a changed
    diagnostic module alter output while the permit stayed valid."""
    assert any("gate2_diagnostics.py" in f for f in C.IMPLEMENTATION_FILES)
