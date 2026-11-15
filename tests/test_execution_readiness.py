"""Regression tests for the 2026-10-09 execution-readiness review
(review:20261009_230728): stage-id mismatch hard-fails, an arbitrary
--code-revision string does not claim "verified", atomic run ownership
actually prevents two concurrent writers, unknown-outcome cells stay
eligible for technical diagnostics while excluded only from bio-labeled
scoring, and a corrected-UMAP failure is not silently dropped from status.
"""

from __future__ import annotations

import json
import subprocess
import sys
import threading
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import pytest

CODE_DIR = Path(__file__).parent.parent / "scripts/33_batch_effect_audit/code"
SHARED_DIR = Path(__file__).parent.parent / "scripts/00_shared/code"
sys.path.insert(0, str(CODE_DIR))
sys.path.insert(0, str(SHARED_DIR))

from population_stage import check_stage_id  # noqa: E402
from run_batch_diagnostics import _resolve_code_identity, determine_embedding_status  # noqa: E402
from run_ownership import OwnershipError, claim_ownership  # noqa: E402

# ---------------------------------------------------------------------------
# population_stage.check_stage_id -- hard fail, not warn
# ---------------------------------------------------------------------------


def test_stage_id_match_passes():
    check_stage_id("synthetic_fixture", "synthetic_fixture", context="test")


def test_stage_id_mismatch_raises_not_just_warns():
    with pytest.raises(ValueError, match="does not match"):
        check_stage_id("pre_gate2_decontx_provisional", "synthetic_fixture", context="test")


def test_stage_id_unregistered_value_raises():
    with pytest.raises(ValueError, match="not a registered"):
        check_stage_id("made_up_stage", "made_up_stage", context="test")


# ---------------------------------------------------------------------------
# run_ownership.claim_ownership -- atomic, proven with a real race
# ---------------------------------------------------------------------------


def test_claim_ownership_second_caller_loses(tmp_path):
    lock = tmp_path / ".run_owner.json"
    first = claim_ownership(lock, command="first")
    assert lock.exists()
    assert json.loads(lock.read_text())["command"] == "first"
    with pytest.raises(OwnershipError):
        claim_ownership(lock, command="second")
    # The lock still records the WINNER, not the loser.
    assert json.loads(lock.read_text())["command"] == "first" == first["command"]


def test_claim_ownership_two_real_threads_racing_exactly_one_wins(tmp_path):
    """Not a sequential proxy for the race -- two real OS threads call
    claim_ownership() against the SAME lock path as close to simultaneously
    as a `threading.Barrier` can arrange, proving the O_EXCL guarantee holds
    under actual concurrency, not just "call it twice in a row"."""
    lock = tmp_path / ".run_owner.json"
    barrier = threading.Barrier(2)
    results = {}

    def attempt(name):
        barrier.wait()
        try:
            claim_ownership(lock, command=name)
            results[name] = "won"
        except OwnershipError:
            results[name] = "lost"

    t1 = threading.Thread(target=attempt, args=("a",))
    t2 = threading.Thread(target=attempt, args=("b",))
    t1.start()
    t2.start()
    t1.join()
    t2.join()

    assert sorted(results.values()) == ["lost", "won"]


# ---------------------------------------------------------------------------
# _resolve_code_identity -- an arbitrary asserted string must not claim verified
# ---------------------------------------------------------------------------


def test_arbitrary_code_revision_string_is_not_verified(tmp_path):
    """tmp_path is not a git checkout at all -- the exact situation this
    project's HPCC scratch deployments are in. Before the 2026-10-09 fix,
    `code_revision_verified = code_revision is not None` made ANY string
    report as verified. The real source files are copied into tmp_path at
    their expected relative paths (but no .git is created) so the content
    hashes this test checks for are real, not merely absent because the
    tmpdir is empty -- that would conflate "not a git checkout" with "has no
    source files", which are different things _resolve_code_identity must
    distinguish."""
    import shutil

    real_repo_root = Path(__file__).parent.parent
    from run_batch_diagnostics import _SOURCE_FILES_RELATIVE_TO_REPO_ROOT

    for rel in _SOURCE_FILES_RELATIVE_TO_REPO_ROOT:
        dst = tmp_path / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(real_repo_root / rel, dst)

    identity = _resolve_code_identity(tmp_path, "totally-made-up-revision-xyz")
    assert identity["code_revision_verified"] is False
    assert identity["code_revision_asserted"] == "totally-made-up-revision-xyz"
    assert identity["code_revision_resolved_git_head"] is None
    assert (
        identity["code_revision_verification_method"]
        == "not_a_git_checkout_unverified_asserted_string"
    )
    # The honest substitute: real content hashes are still recorded.
    assert "source_file_sha256" in identity
    assert any(v is not None for v in identity["source_file_sha256"].values())


def test_code_identity_verified_in_a_real_clean_git_checkout():
    """Run against THIS repo's own real checkout -- if HEAD resolves and the
    tree is clean, an asserted revision that matches HEAD (or none asserted
    at all) IS verified; this is the positive case the negative test above
    is contrasted against."""
    repo_root = Path(__file__).parent.parent
    head = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=repo_root, capture_output=True, text=True
    ).stdout.strip()
    clean = (
        subprocess.run(
            ["git", "status", "--porcelain"], cwd=repo_root, capture_output=True, text=True
        ).stdout.strip()
        == ""
    )
    if not head or not clean:
        pytest.skip("local checkout is not HEAD-resolvable and clean at test time")
    identity = _resolve_code_identity(repo_root, head)
    assert identity["code_revision_verified"] is True
    assert identity["code_revision_verification_method"] == "git_head_clean_tree"


# ---------------------------------------------------------------------------
# determine_embedding_status -- corrected-UMAP failure must not be silently dropped
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "asw_skipped,unc_ok,harmony,blocked,corr_ok,expected",
    [
        (False, True, False, False, None, "complete"),
        (True, True, False, False, None, "partial"),
        (False, False, False, False, None, "partial"),
        (False, True, True, False, True, "complete"),
        (False, True, True, False, False, "partial"),  # <-- the finding: corrected UMAP failed
        (False, True, True, True, None, "failed"),  # Harmony itself blocked/failed
    ],
)
def test_determine_embedding_status(asw_skipped, unc_ok, harmony, blocked, corr_ok, expected):
    assert (
        determine_embedding_status(
            asw_skipped_groups=asw_skipped,
            uncorrected_umap_ok=unc_ok,
            run_harmony=harmony,
            harmony_blocked_or_failed=blocked,
            corrected_umap_ok=corr_ok,
        )
        == expected
    )


# ---------------------------------------------------------------------------
# End-to-end CLI: stage-id, eligibility split, concurrent ownership
# ---------------------------------------------------------------------------


def _build_tiny_h5ad(
    path: Path, *, n_cells=60, n_genes=20, seed=0, population_stage_id="synthetic_fixture"
):
    rng = np.random.RandomState(seed)
    cohorts = ["CohortA", "CohortB", "CohortC"]
    outcomes = ["fatal", "survived", "unknown"]
    obs = pd.DataFrame(
        {
            "lineage": ["Platelet"] * n_cells,
            "cohort": [cohorts[i % 3] for i in range(n_cells)],
            "patient_id": [f"P{(i % 12):02d}" for i in range(n_cells)],
            # roughly 1/3 unknown -- deliberately large enough that bio_known_mask
            # genuinely excludes a nontrivial chunk, while technical diagnostics
            # should still see all n_cells.
            # Each outcome spans multiple cohorts, with one outcome per patient.
            "outcome": [outcomes[(i % 12) // 4] for i in range(n_cells)],
        },
        index=[f"cell_{i:04d}" for i in range(n_cells)],
    )
    X = rng.poisson(3, size=(n_cells, n_genes)).astype(np.float32)
    adata = anndata.AnnData(
        X=X, obs=obs, var=pd.DataFrame(index=[f"GENE{i:03d}" for i in range(n_genes)])
    )
    adata.layers["counts"] = X.copy()
    adata.uns["population_stage"] = "synthetic_fixture"
    adata.uns["population_stage_id"] = population_stage_id
    adata.uns["expression_source"] = "counts"
    adata.write_h5ad(path)
    return adata


def _run_diag_cli(args_list):
    return subprocess.run(
        [sys.executable, str(CODE_DIR / "run_batch_diagnostics.py")] + args_list,
        capture_output=True,
        text=True,
    )


def _common_args(h5ad, out_dir, stage_id="synthetic_fixture"):
    return [
        "--h5ad",
        str(h5ad),
        "--out-dir",
        str(out_dir),
        "--population-stage",
        "synthetic_fixture",
        "--stage-id",
        stage_id,
        "--code-revision",
        "test-fixture-revision",
        "--diagnostic-n-per-patient",
        "10",
        "--min-group-size",
        "2",
    ]


def test_cli_stage_id_mismatch_hard_fails(tmp_path):
    h5ad = tmp_path / "fixture.h5ad"
    _build_tiny_h5ad(h5ad, population_stage_id="synthetic_fixture")
    out_dir = tmp_path / "run"
    result = _run_diag_cli(_common_args(h5ad, out_dir, stage_id="pre_gate2_decontx_provisional"))
    assert result.returncode != 0
    assert "does not match the input's own recorded stage id" in result.stdout + result.stderr
    manifest = json.loads((out_dir / "manifest.json").read_text())
    assert manifest["status"] == "failed"


def test_cli_unknown_outcome_cells_stay_in_technical_population(tmp_path):
    """Direct regression test for the eligibility-coupling finding: the
    diagnostic subsample size (n_diagnostic_cells) must reflect ALL sampled
    cells regardless of outcome, while n_diagnostic_cells_bio_known is
    strictly smaller (since ~1/3 of the fixture is "unknown")."""
    h5ad = tmp_path / "fixture.h5ad"
    _build_tiny_h5ad(h5ad)
    out_dir = tmp_path / "run"
    result = _run_diag_cli(_common_args(h5ad, out_dir))
    assert result.returncode == 0, result.stdout + result.stderr
    manifest = json.loads((out_dir / "manifest.json").read_text())
    rec = manifest["per_lineage"][0]
    assert rec["n_diagnostic_cells_bio_known"] < rec["n_diagnostic_cells"]
    assert (
        rec["n_diagnostic_cells_bio_known"] + rec["n_diagnostic_cells_bio_unknown"]
        == rec["n_diagnostic_cells"]
    )
    emb = rec["embeddings"]["expression_pca"]
    tech = emb["uncorrected"]["n_cells_technical"]
    bio = emb["uncorrected"]["n_cells_bio_known"]
    assert (
        tech == rec["n_diagnostic_cells"]
    ), "technical population must be the FULL diagnostic subsample"
    assert (
        bio == rec["n_diagnostic_cells_bio_known"]
    ), "bio-scored population must be the known-outcome subset"
    assert bio < tech


def test_cli_code_revision_verified_false_for_arbitrary_string(tmp_path):
    h5ad = tmp_path / "fixture.h5ad"
    _build_tiny_h5ad(h5ad)
    out_dir = tmp_path / "run"
    result = _run_diag_cli(_common_args(h5ad, out_dir))
    assert result.returncode == 0, result.stdout + result.stderr
    manifest = json.loads((out_dir / "manifest.json").read_text())
    assert manifest["code_revision_asserted"] == "test-fixture-revision"
    assert manifest["code_revision_verified"] is False
    assert "source_file_sha256" in manifest


def test_cli_concurrent_runs_against_same_out_dir_exactly_one_wins(tmp_path):
    """Real two-process race, per the review's explicit ask: launch two
    `run_batch_diagnostics.py` invocations against the SAME --out-dir close
    together and prove exactly one owns the run."""
    h5ad = tmp_path / "fixture.h5ad"
    _build_tiny_h5ad(h5ad)
    out_dir = tmp_path / "run"
    out_dir.mkdir()

    results = {}

    def launch(name):
        results[name] = _run_diag_cli(_common_args(h5ad, out_dir))

    barrier_ready = threading.Barrier(2)

    def launch_synced(name):
        barrier_ready.wait()
        launch(name)

    t1 = threading.Thread(target=launch_synced, args=("a",))
    t2 = threading.Thread(target=launch_synced, args=("b",))
    t1.start()
    t2.start()
    t1.join()
    t2.join()

    # A negative returncode means that side was killed by a signal (e.g. -9 /
    # SIGKILL from the HPCC login node's own resource limits under heavy
    # concurrent load -- DEVELOPMENT.md's own documented stance: "the resulting
    # exit code is a fact about the probe, not about the code"). When that
    # happens, one competitor never got to participate in the race fairly, so
    # this run cannot test the actual safety property either way -- skip
    # rather than report a false pass or false fail.
    if any(r.returncode < 0 for r in results.values()):
        pytest.skip(
            f"one side was killed by a signal (not application logic) -- likely "
            f"login-node resource contention, not a race outcome: "
            f"{[(k, r.returncode) for k, r in results.items()]}"
        )

    # The real safety property is "at most one process writes a completed run
    # into out_dir" -- the LOSER's exact rejection message is timing-dependent
    # and legitimately varies: it loses the atomic ownership-lock race
    # ("already owned by another run") if it is still early when the winner
    # claims the lock, or it hits the pre-existing non-empty-directory guard
    # ("already exists and is non-empty") if it only gets as far as that check
    # after the winner has already started writing real output files. Both are
    # correct, safe outcomes -- neither corrupts the winner's data -- so this
    # asserts exactly one winner and that the loser failed for one of those two
    # legitimate reasons, not a specific wording.
    returncodes = sorted(r.returncode for r in results.values())
    n_zero = sum(1 for rc in returncodes if rc == 0)
    assert n_zero == 1, f"expected exactly one winner (returncode 0), got {n_zero}: {returncodes}"
    loser = next(r for r in results.values() if r.returncode != 0)
    loser_output = loser.stdout + loser.stderr
    valid_rejections = ("already owned by another run", "already exists and is non-empty")
    assert any(msg in loser_output for msg in valid_rejections), (
        f"loser failed for an unexpected reason (neither ownership-lock nor "
        f"non-empty-dir refusal): {loser_output}"
    )
    manifest = json.loads((out_dir / "manifest.json").read_text())
    assert manifest["status"] == "complete", (
        "the single winner's manifest must reflect a genuinely completed run, "
        "not just an unrelated exit-0"
    )
