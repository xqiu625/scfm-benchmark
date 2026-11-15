"""Real-fixture tests for `stephenson_monocyte_scprint2_pilot_row_order_verify.py`.

These tests prove the VERIFIER's own sensitivity using real small h5ad files
written through the real anndata writer (so the h5py-reading code in the
verifier is exercised against the real on-disk format) with SYNTHETIC
embedding arrays standing in for two real `scprint2 embed` runs. This is a
stub at the level DEVELOPMENT.md calls out explicitly: it proves the CHECKER
correctly distinguishes aligned/misaligned embeddings, not that the real
model's forward pass is permutation-invariant under `how="random expr"` --
that remains open pending the real GPU/CPU smoke test
(`run_scprint2_row_order_smoke_test.sh`, Packet B, not submitted this round).

`test_cross_run_consistency_alone_misses_a_shared_systematic_bug` and
`test_duplicate_row_self_consistency_catches_the_same_bug` together
reproduce review:20261018_000122's independent probe: two runs
sharing the exact same systematic reversal pass the cross-run check with
`max_abs_diff=0.0` even though every non-center row is wrong, and the
duplicate-row self-consistency check catches that same bug within a single
run. Documented as a real, asserted scope limit -- not silently left open.
"""
import importlib.util
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

THIS_DIR = Path(__file__).resolve().parent


def _load_module(name, filename):
    spec = importlib.util.spec_from_file_location(name, THIS_DIR / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


verify = _load_module("row_order_verify_mod", "stephenson_monocyte_scprint2_pilot_row_order_verify.py")

N = 9  # deliberately not a power of 2 / not evenly divisible by a batch size -- a non-full tail batch
IDS = [f"steph_BARCODE{i}-DONOR1" for i in range(N)]


def _write_recovered_run(path, ids, emb, obs_extra=None):
    obs = pd.DataFrame({"placeholder": [0] * len(ids)}, index=ids)
    if obs_extra:
        for k, v in obs_extra.items():
            obs[k] = v
    adata = ad.AnnData(X=np.zeros((len(ids), 2), dtype=np.float32), obs=obs)
    adata.obsm["scprint_emb"] = emb
    adata.write_h5ad(path)
    return str(path)


def _make_base_embeddings(seed):
    rng = np.random.default_rng(seed)
    return rng.random((N, 6)).astype(np.float32)


def test_identical_runs_confirm_row_order(tmp_path):
    emb = _make_base_embeddings(0)
    run_a = _write_recovered_run(tmp_path / "run_a.h5ad", IDS, emb)
    run_b = _write_recovered_run(tmp_path / "run_b.h5ad", IDS, emb.copy())

    result = verify.compare_two_runs_by_recovered_id(run_a, run_b, atol=1e-4, rtol=1e-3)
    assert result["status"] == "cross_run_consistency_confirmed_for_tested_cells"
    assert result["n_common_ids"] == N
    assert result["max_abs_diff_among_common_ids"] == 0.0


def test_permuted_but_correctly_matched_runs_still_confirm(tmp_path):
    """Mimics the real smoke test's own structure: run B processed cells in a
    DIFFERENT order, but each cell's OWN embedding (by real-world content) is
    unchanged -- only the row POSITION differs, not which value belongs to
    which id. The verifier must match by id, not position, and still pass."""
    emb_a = _make_base_embeddings(1)
    perm = list(reversed(range(N)))
    ids_b = [IDS[i] for i in perm]
    emb_b = emb_a[perm]

    run_a = _write_recovered_run(tmp_path / "run_a.h5ad", IDS, emb_a)
    run_b = _write_recovered_run(tmp_path / "run_b.h5ad", ids_b, emb_b)

    result = verify.compare_two_runs_by_recovered_id(run_a, run_b, atol=1e-4, rtol=1e-3)
    assert result["status"] == "cross_run_consistency_confirmed_for_tested_cells"
    assert result["n_common_ids"] == N


def test_fault_injection_swap_is_caught(tmp_path):
    """The exact scenario review:20261017_223541 requires be
    covered: two same-sample rows' embedding-only values swapped between the
    two runs, with ids/obs otherwise intact. Mutation-kill target: this test
    must FAIL if compare_two_runs_by_recovered_id's value-mismatch raise is
    removed or its tolerance check is broken."""
    emb_a = _make_base_embeddings(2)
    emb_b = emb_a.copy()
    emb_b[[0, 1]] = emb_b[[1, 0]]  # swap rows 0 and 1's embedding values only

    run_a = _write_recovered_run(tmp_path / "run_a.h5ad", IDS, emb_a)
    run_b = _write_recovered_run(tmp_path / "run_b.h5ad", IDS, emb_b)

    with pytest.raises(verify.RowOrderMismatch, match="value_mismatch|value mismatch"):
        verify.compare_two_runs_by_recovered_id(run_a, run_b, atol=1e-4, rtol=1e-3)


def test_refuses_on_id_set_mismatch(tmp_path):
    emb_a = _make_base_embeddings(3)
    emb_b = _make_base_embeddings(3)
    ids_b = IDS[:-1] + ["steph_NOT_IN_A-DONOR9"]  # drop last id, add an unrelated one

    run_a = _write_recovered_run(tmp_path / "run_a.h5ad", IDS, emb_a)
    run_b = _write_recovered_run(tmp_path / "run_b.h5ad", ids_b, emb_b)

    with pytest.raises(verify.RowOrderMismatch, match="id sets differ"):
        verify.compare_two_runs_by_recovered_id(run_a, run_b, atol=1e-4, rtol=1e-3)


def test_refuses_on_duplicate_ids_within_a_run(tmp_path):
    emb_a = _make_base_embeddings(4)
    dup_ids = [IDS[0]] + IDS[1:]  # make id 0 unique, then duplicate it onto index 1 too
    dup_ids[1] = IDS[0]

    run_a_path = tmp_path / "run_a.h5ad"
    _write_recovered_run(run_a_path, dup_ids, emb_a)
    run_b = _write_recovered_run(tmp_path / "run_b.h5ad", IDS, _make_base_embeddings(4))

    with pytest.raises(verify.RowOrderMismatch, match="not unique"):
        verify.compare_two_runs_by_recovered_id(str(run_a_path), run_b, atol=1e-4, rtol=1e-3)


def test_tolerance_is_respected_both_directions(tmp_path):
    """A small, pre-declared-tolerance-sized perturbation must PASS; a
    perturbation just beyond it must FAIL -- proves the declared (atol, rtol)
    actually governs the outcome rather than some hidden default."""
    emb_a = _make_base_embeddings(5)

    within_tol = emb_a.copy()
    within_tol[0, 0] += 5e-5  # smaller than atol=1e-4
    run_a = _write_recovered_run(tmp_path / "run_a.h5ad", IDS, emb_a)
    run_b_ok = _write_recovered_run(tmp_path / "run_b_ok.h5ad", IDS, within_tol)
    result = verify.compare_two_runs_by_recovered_id(run_a, run_b_ok, atol=1e-4, rtol=1e-3)
    assert result["status"] == "cross_run_consistency_confirmed_for_tested_cells"

    beyond_tol = emb_a.copy()
    beyond_tol[0, 0] += 1.0  # far beyond atol/rtol
    run_b_bad = _write_recovered_run(tmp_path / "run_b_bad.h5ad", IDS, beyond_tol)
    with pytest.raises(verify.RowOrderMismatch, match="value_mismatch|value mismatch"):
        verify.compare_two_runs_by_recovered_id(run_a, run_b_bad, atol=1e-4, rtol=1e-3)


def test_refuses_on_matching_infinity(tmp_path):
    """Mutation-kill target: review:20261018_000122's independent
    probe found two runs that were BOTH Inf at the same position returned a
    clean PASS (np.allclose(inf, inf, equal_nan=False) is True). Must now
    refuse outright rather than treat Inf==Inf as agreement."""
    emb_a = _make_base_embeddings(6)
    emb_a[0, 0] = float("inf")
    emb_b = emb_a.copy()  # same Inf, same position -- the exact probe shape

    run_a = _write_recovered_run(tmp_path / "run_a.h5ad", IDS, emb_a)
    run_b = _write_recovered_run(tmp_path / "run_b.h5ad", IDS, emb_b)

    with pytest.raises(verify.RowOrderMismatch, match="non-finite"):
        verify.compare_two_runs_by_recovered_id(run_a, run_b, atol=1e-4, rtol=1e-3)


def test_refuses_on_nan(tmp_path):
    emb_a = _make_base_embeddings(7)
    emb_a[2, 3] = float("nan")
    run_a = _write_recovered_run(tmp_path / "run_a.h5ad", IDS, emb_a)
    run_b = _write_recovered_run(tmp_path / "run_b.h5ad", IDS, _make_base_embeddings(7))

    with pytest.raises(verify.RowOrderMismatch, match="non-finite"):
        verify.compare_two_runs_by_recovered_id(run_a, run_b, atol=1e-4, rtol=1e-3)


def test_refuses_on_embedding_row_count_mismatch(tmp_path):
    """anndata itself refuses to let obsm/X row counts disagree at write time,
    so this simulates a corrupted/truncated embedding array the way it could
    actually occur on disk: writing a normal, consistent h5ad first, then
    truncating obsm/scprint_emb in place via h5py (bypassing anndata's own
    write-time validation, the way file corruption or a partial write would)."""
    import h5py as _h5py

    run_a_path = tmp_path / "run_a.h5ad"
    _write_recovered_run(run_a_path, IDS, _make_base_embeddings(8))
    with _h5py.File(run_a_path, "r+") as f:
        emb_ds = f["obsm"]["scprint_emb"]
        truncated = emb_ds[:-1]
        del f["obsm"]["scprint_emb"]
        f["obsm"].create_dataset("scprint_emb", data=truncated)

    run_b = _write_recovered_run(tmp_path / "run_b.h5ad", IDS, _make_base_embeddings(8))

    with pytest.raises(verify.RowOrderMismatch, match="does not match"):
        verify.compare_two_runs_by_recovered_id(str(run_a_path), run_b, atol=1e-4, rtol=1e-3)


@pytest.mark.parametrize("bad_atol,bad_rtol", [(-1e-4, 1e-3), (1e-4, -1e-3), (float("nan"), 1e-3), (1e-4, float("inf"))])
def test_refuses_on_invalid_tolerance(tmp_path, bad_atol, bad_rtol):
    emb = _make_base_embeddings(9)
    run_a = _write_recovered_run(tmp_path / "run_a.h5ad", IDS, emb)
    run_b = _write_recovered_run(tmp_path / "run_b.h5ad", IDS, emb.copy())
    with pytest.raises(verify.RowOrderMismatch, match="tolerance"):
        verify.compare_two_runs_by_recovered_id(run_a, run_b, atol=bad_atol, rtol=bad_rtol)


def test_duplicate_row_self_consistency_happy_path(tmp_path):
    """A real cell's content duplicated at two different positions within ONE
    run must get the same embedding -- the independent, position-independent
    oracle this module adds alongside cross-run comparison."""
    emb = _make_base_embeddings(10)
    ids = list(IDS)
    dup_source_pos, dup_target_pos = 0, N - 1
    dup_id = "steph_DUPLICATE_OF_0-DONOR1"
    ids_with_dup = ids + [dup_id]
    emb_with_dup = np.vstack([emb, emb[dup_source_pos : dup_source_pos + 1]])  # identical content

    run_path = _write_recovered_run(tmp_path / "run.h5ad", ids_with_dup, emb_with_dup)
    result = verify.check_duplicate_row_self_consistency(
        run_path, [(ids[dup_source_pos], dup_id)], atol=1e-4, rtol=1e-3
    )
    assert result["status"] == "content_identity_self_consistency_confirmed"
    assert result["max_abs_diff_among_pairs"] == 0.0


def test_cross_run_consistency_alone_misses_a_shared_systematic_bug(tmp_path):
    """Reproduces review:20261018_000122's independent probe
    exactly: both runs apply the SAME systematic reversal to their embedding
    rows (independent of input order), so matching by recovered id lines the
    wrong values back up with each other. This test asserts the known
    scope limit rather than silently hiding it -- see
    test_duplicate_row_self_consistency_catches_the_same_bug for the check
    that DOES catch this."""
    true_emb = _make_base_embeddings(11)
    # Run A: labeled in natural order, but content is the reverse of that order
    # (the "producer always reverses regardless of input order" bug).
    emb_a = true_emb[::-1].copy()
    # Run B: labeled in reversed order, content is the reverse of THAT order,
    # i.e. back to natural order -- the exact pairing that cancels out.
    ids_b = IDS[::-1]
    emb_b = true_emb.copy()

    run_a = _write_recovered_run(tmp_path / "run_a.h5ad", IDS, emb_a)
    run_b = _write_recovered_run(tmp_path / "run_b.h5ad", ids_b, emb_b)

    result = verify.compare_two_runs_by_recovered_id(run_a, run_b, atol=1e-4, rtol=1e-3)
    assert result["status"] == "cross_run_consistency_confirmed_for_tested_cells"
    assert result["max_abs_diff_among_common_ids"] == 0.0
    # every row except the exact center is actually wrong relative to true_emb,
    # despite the clean PASS above -- this is the documented scope limit.
    n_actually_correct = sum(
        np.allclose(emb_a[i], true_emb[i]) for i in range(N)
    )
    assert n_actually_correct < N


def test_duplicate_row_self_consistency_catches_the_same_bug(tmp_path):
    """The same systematic-reversal bug as the test above, but checked via
    the independent, within-run oracle: a duplicated cell's two occurrences
    (inserted at two different positions, both subject to the same fixed
    reversal) get DIFFERENT (wrong) embeddings, which this function must
    catch even though the cross-run check above could not."""
    true_emb = _make_base_embeddings(12)
    dup_id = "steph_DUPLICATE_OF_2-DONOR1"
    # Duplicate of cell 2 (NOT 0): with 0 the appended duplicate's post-reversal
    # position coincidentally mirrors cell 0's own position, making the two
    # (both-wrong) values accidentally equal -- the same class of coincidence
    # as the cross-run test above. Position 2 avoids that coincidence.
    ids_with_dup = list(IDS) + [dup_id]
    content_with_dup = np.vstack([true_emb, true_emb[2:3]])
    # Apply the SAME fixed, position-dependent bug (reverse the whole array)
    # that a buggy producer would apply regardless of input content/order.
    emb_with_dup_buggy = content_with_dup[::-1].copy()

    run_path = _write_recovered_run(tmp_path / "run.h5ad", ids_with_dup, emb_with_dup_buggy)
    with pytest.raises(verify.RowOrderMismatch, match="position-dependent"):
        verify.check_duplicate_row_self_consistency(run_path, [(IDS[2], dup_id)], atol=1e-4, rtol=1e-3)
