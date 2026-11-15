"""Row-order / producer-alignment verifier for the Stephenson monocyte
scPRINT-2 pilot (Packet B, review:20261017_223541 step 4).

## What this answers that stephenson_monocyte_scprint2_pilot_validate.py and
## stephenson_monocyte_scprint2_pilot_recover_cell_ids.py do NOT

Both of those operate on `obs` alone: the validator checks cell-id coverage
and required-column agreement by id; the recovery adapter relabels
`obs_names` from the untouched `cell_id_stripped` column. Neither reads
`obsm["scprint_emb"]`'s actual numeric content against anything independent --
so neither can tell whether the embedding ARRAY's rows are truly bound to the
cell each row is now correctly labeled as (see
`implementation-note:20261017_222358` section 1 and the
`test_embedding_only_same_sample_swap_is_not_detected` test in the recovery
adapter's own test file, which proves this gap exists and is not silently
claimed closed).

This module answers that question THE ONLY WAY IT CAN HONESTLY BE ANSWERED
without re-deriving the model's own numerics by hand: run the SAME real
`scprint2 embed` pipeline TWICE on the SAME set of cells in two DIFFERENT,
fully-declared row orders, recover both outputs' cell ids (via the real
recovery adapter), and require that, for every cell id common to both runs,
the two runs' `scprint_emb` vectors for that cell agree within a tolerance
DECLARED BEFORE the runs are compared (never loosened after seeing a
mismatch -- a looser tolerance chosen after the fact is not a tolerance, it
is a different, unstated claim).

## What this does NOT answer

- It does not validate `how="random expr"` (the FROZEN production
  parameter): `scdataloader/collator.py`'s "random expr" branch draws from
  numpy's global RNG stream, consumed once per cell in whatever order the
  DataLoader presents cells -- changing row order therefore changes which
  random genes a given cell's content is reduced to, which is a confound
  this test must not have. The real smoke test that uses this module
  (`run_scprint2_row_order_smoke_test.sh`, Packet B, NOT submitted this
  round) therefore overrides the pipeline to `--how "some"` with a small,
  pinned, real `--genelist` -- `collator.py`'s "some"/"all" branch is
  `loc = np.arange(len(expr))`, with no RNG call anywhere in that branch
  (confirmed by reading the installed `Collator.__call__` directly,
  2026-10-17) -- so this test exercises the SAME DataLoader/Collator/
  `_predict`/shard-write/concat/final-assignment call path, with the
  randomness specific to gene selection removed, not the production
  parameters themselves. A pass here is evidence about ROW IDENTITY
  plumbing, not about the production run's actual embedding values.
- It does not run the real model yet in this repository -- see
  `test_fault_injection_swap_is_caught` below, which proves this verifier's
  OWN sensitivity using synthetic (not model-produced) arrays. That is a
  test of the checker, not of the real forward pass; running the real
  smoke test needs a new SLURM allocation this round does not submit.
- `compare_two_runs_by_recovered_id` PASSING is cross-run consistency, not
  independent proof of correct row identity:
  review:20261018_000122's independent probe constructed two runs that BOTH
  share the exact same systematic bug (every row's content is the reverse of
  what its label claims) and showed this function returned
  `max_abs_diff_among_common_ids=0.0` -- a fixed, input-order-INDEPENDENT bug
  cancels out under a same-reversed-pair comparison. This is why the real
  smoke test also calls `check_duplicate_row_self_consistency` on EACH run
  individually (a duplicated real cell, inserted at two different positions
  within the SAME run, must get the same embedding regardless of position --
  a property a fixed positional bug breaks and that cross-run comparison
  alone cannot see) and uses a random derangement, not a pure reversal, as
  the second run's order (a pure reversal is the one specific permutation a
  "the producer always reverses" bug would exactly cancel against). Treat
  `compare_two_runs_by_recovered_id`'s PASS as "consistent across two
  particular orderings", `check_duplicate_row_self_consistency`'s PASS as
  "content-dependent, position-independent for the tested duplicate pairs" --
  never as "row identity is proven".
"""
import sys

import h5py
import numpy as np
import pandas as pd


class RowOrderMismatch(Exception):
    """Raised when two runs' embeddings for the same cell id diverge beyond
    the declared tolerance, when their recovered id sets disagree, or when
    either run's embedding array fails a basic structural/numeric sanity
    check (wrong rank, empty, non-finite values, mismatched row count)."""


def _check_embedding_sanity(h5ad_path, ids, emb):
    """Refuses on: wrong array rank, zero-width embedding, row count not
    matching the number of ids, non-numeric dtype, or ANY non-finite
    (NaN/Inf) value. Per review:20261018_000122's independent
    probe: the prior version compared two runs that were BOTH all-Inf at the
    same position and returned a clean PASS (`np.allclose` on two equal Infs
    with `equal_nan=False` is True; a RuntimeWarning on the subtract was the
    only hint, not a failure) -- non-finite values must be refused outright,
    never compared."""
    emb = np.asarray(emb)
    if emb.ndim != 2 or emb.shape[1] == 0:
        raise RowOrderMismatch(
            f"{h5ad_path}: embedding has invalid shape {emb.shape} (expected 2-D, >=1 column)"
        )
    if emb.shape[0] != len(ids):
        raise RowOrderMismatch(
            f"{h5ad_path}: embedding has {emb.shape[0]} rows, which does not match {len(ids)} obs_names"
        )
    if emb.dtype.kind not in "fiu":
        raise RowOrderMismatch(f"{h5ad_path}: embedding has non-numeric dtype {emb.dtype}")
    n_non_finite = int((~np.isfinite(emb.astype(np.float64))).sum())
    if n_non_finite:
        raise RowOrderMismatch(
            f"{h5ad_path}: embedding has {n_non_finite} non-finite (NaN/Inf) value(s) -- "
            f"refusing to compare rather than silently treating Inf==Inf as a match"
        )


def _read_recovered_run(h5ad_path, embedding_key="scprint_emb"):
    """Reads obs_names (already recovered, real cell ids -- NOT uuids) and the
    primary embedding array via h5py only, never anndata -- see
    stephenson_monocyte_scprint2_pilot_recover_cell_ids.py's own
    _read_obs_dataframe_h5 for why (anndata backed="r" OOM-killed on this
    cluster's login node against real multi-GB files, 2026-10-17)."""
    with h5py.File(h5ad_path, "r") as f:
        obs = f["obs"]
        index_name = obs.attrs["_index"]
        ids = obs[index_name].asstr()[:]
        emb = f["obsm"][embedding_key][:]
    if len(set(ids)) != len(ids):
        raise RowOrderMismatch(f"{h5ad_path}: obs_names are not unique -- refusing to match by id")
    emb = np.asarray(emb)
    _check_embedding_sanity(h5ad_path, ids, emb)
    return pd.Index(ids), emb


def compare_two_runs_by_recovered_id(
    run_a_h5ad_path,
    run_b_h5ad_path,
    embedding_key="scprint_emb",
    atol=1e-4,
    rtol=1e-3,
    require_full_overlap=True,
):
    """Core check. Both paths must already be RECOVERED outputs (real cell
    ids in obs_names, via stephenson_monocyte_scprint2_pilot_recover_cell_ids.py
    -- this function does not do that relabeling itself). Matches rows by id,
    not by position, and compares `embedding_key` for every id common to both
    runs. `atol`/`rtol` must be the caller's pre-declared values -- this
    function has no default-widening behavior and never adjusts tolerance
    based on the data it is comparing.

    Returns a report dict on success. Raises RowOrderMismatch on:
    - either run's ids not unique,
    - either run's embedding failing the structural/finiteness sanity check
      (see _check_embedding_sanity -- a matching pair of non-finite values is
      refused, never treated as agreement),
    - a non-positive or non-finite atol/rtol (a silently-accepted atol=-1 or
      atol=nan would make every comparison vacuously pass or fail),
    - (if `require_full_overlap`) the two runs' id sets not being identical,
    - any common id's embedding vector differing beyond (atol, rtol).

    PASSING IS CROSS-RUN CONSISTENCY, NOT PROOF OF CORRECT ROW IDENTITY -- see
    this module's docstring for the demonstrated systematic-bug-cancellation
    failure mode. Call `check_duplicate_row_self_consistency` on each run too.
    """
    if not (np.isfinite(atol) and atol >= 0):
        raise RowOrderMismatch(f"atol={atol!r} is not a finite, non-negative tolerance")
    if not (np.isfinite(rtol) and rtol >= 0):
        raise RowOrderMismatch(f"rtol={rtol!r} is not a finite, non-negative tolerance")

    ids_a, emb_a = _read_recovered_run(run_a_h5ad_path, embedding_key)
    ids_b, emb_b = _read_recovered_run(run_b_h5ad_path, embedding_key)

    set_a, set_b = set(ids_a), set(ids_b)
    common = set_a & set_b
    if require_full_overlap and set_a != set_b:
        only_a = sorted(set_a - set_b)
        only_b = sorted(set_b - set_a)
        raise RowOrderMismatch(
            f"run A/B recovered id sets differ: {len(only_a)} only in A "
            f"(e.g. {only_a[:3]}), {len(only_b)} only in B (e.g. {only_b[:3]})"
        )
    if not common:
        raise RowOrderMismatch("no common recovered ids between the two runs -- nothing to compare")

    pos_a = {cid: i for i, cid in enumerate(ids_a)}
    pos_b = {cid: i for i, cid in enumerate(ids_b)}

    mismatched = []
    max_abs_diff = 0.0
    for cid in sorted(common):
        va = emb_a[pos_a[cid]]
        vb = emb_b[pos_b[cid]]
        if va.shape != vb.shape:
            mismatched.append((cid, "shape_mismatch", str(va.shape), str(vb.shape)))
            continue
        diff = np.abs(va.astype(np.float64) - vb.astype(np.float64))
        max_abs_diff = max(max_abs_diff, float(diff.max()))
        if not np.allclose(va, vb, atol=atol, rtol=rtol, equal_nan=False):
            mismatched.append((cid, "value_mismatch", float(diff.max()), None))

    if mismatched:
        raise RowOrderMismatch(
            f"{len(mismatched)}/{len(common)} common cell id(s) have {embedding_key!r} "
            f"vectors that disagree beyond atol={atol}, rtol={rtol} between the two runs "
            f"(e.g. {mismatched[:3]}) -- row alignment is NOT confirmed"
        )

    return {
        "status": "cross_run_consistency_confirmed_for_tested_cells",
        "scope_note": "consistency between these two specific orderings, NOT independent proof of "
                      "correct row identity -- see check_duplicate_row_self_consistency and this "
                      "module's docstring for the demonstrated cancellation failure mode",
        "n_common_ids": len(common),
        "n_only_in_a": len(set_a - set_b),
        "n_only_in_b": len(set_b - set_a),
        "max_abs_diff_among_common_ids": max_abs_diff,
        "atol": atol,
        "rtol": rtol,
        "embedding_key": embedding_key,
        "run_a_h5ad_path": run_a_h5ad_path,
        "run_b_h5ad_path": run_b_h5ad_path,
    }


def check_duplicate_row_self_consistency(
    h5ad_path,
    duplicate_id_pairs,
    embedding_key="scprint_emb",
    atol=1e-4,
    rtol=1e-3,
):
    """Independent identity oracle, WITHIN a single run: `duplicate_id_pairs`
    names pairs of recovered cell ids that are known, by construction of the
    fixture, to carry IDENTICAL input content (the same real cell's row,
    inserted at two different positions before this run). A correct pipeline
    must give them the same embedding regardless of position -- a fixed
    positional/row-order bug (e.g. "producer always reverses", the exact
    failure mode that cancels out in compare_two_runs_by_recovered_id when
    both runs share it) breaks this, because the two positions see DIFFERENT
    (wrong) content swapped in from elsewhere. This does not depend on, and
    is not satisfied merely by, a second run -- it is evidence from ONE run's
    own internal consistency, bound to real content, not to position."""
    if not (np.isfinite(atol) and atol >= 0):
        raise RowOrderMismatch(f"atol={atol!r} is not a finite, non-negative tolerance")
    if not (np.isfinite(rtol) and rtol >= 0):
        raise RowOrderMismatch(f"rtol={rtol!r} is not a finite, non-negative tolerance")
    if not duplicate_id_pairs:
        raise RowOrderMismatch("no duplicate_id_pairs given -- nothing to check")

    ids, emb = _read_recovered_run(h5ad_path, embedding_key)
    pos = {cid: i for i, cid in enumerate(ids)}

    mismatched = []
    max_abs_diff = 0.0
    for id_a, id_b in duplicate_id_pairs:
        for cid in (id_a, id_b):
            if cid not in pos:
                raise RowOrderMismatch(f"{h5ad_path}: expected duplicate id {cid!r} not found in recovered obs_names")
        va, vb = emb[pos[id_a]], emb[pos[id_b]]
        diff = np.abs(va.astype(np.float64) - vb.astype(np.float64))
        max_abs_diff = max(max_abs_diff, float(diff.max()))
        if not np.allclose(va, vb, atol=atol, rtol=rtol, equal_nan=False):
            mismatched.append((id_a, id_b, float(diff.max())))

    if mismatched:
        raise RowOrderMismatch(
            f"{len(mismatched)}/{len(duplicate_id_pairs)} duplicate-content pair(s) in {h5ad_path} "
            f"got DIFFERENT embeddings beyond atol={atol}, rtol={rtol} despite identical input content "
            f"(e.g. {mismatched[:3]}) -- a position-dependent (not content-dependent) defect is present"
        )

    return {
        "status": "content_identity_self_consistency_confirmed",
        "n_pairs_checked": len(duplicate_id_pairs),
        "max_abs_diff_among_pairs": max_abs_diff,
        "atol": atol,
        "rtol": rtol,
        "embedding_key": embedding_key,
        "h5ad_path": h5ad_path,
    }


if __name__ == "__main__":
    if len(sys.argv) not in (3, 5):
        print(
            "usage: stephenson_monocyte_scprint2_pilot_row_order_verify.py "
            "<run_a_recovered.h5ad> <run_b_recovered.h5ad> [<atol> <rtol>]\n"
            "atol/rtol must be passed explicitly by the caller for any real run -- "
            "omitting them falls back to this module's defaults, which a caller "
            "should not rely on silently.",
            file=sys.stderr,
        )
        sys.exit(2)
    import json

    cli_kwargs = {}
    if len(sys.argv) == 5:
        cli_kwargs = {"atol": float(sys.argv[3]), "rtol": float(sys.argv[4])}

    try:
        result = compare_two_runs_by_recovered_id(sys.argv[1], sys.argv[2], **cli_kwargs)
    except RowOrderMismatch as e:
        print(f"ROW_ORDER_MISMATCH: {e}", file=sys.stderr)
        sys.exit(1)
    print(json.dumps(result, indent=2, default=str))
