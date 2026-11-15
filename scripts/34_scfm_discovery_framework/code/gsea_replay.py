"""Replays a previously-run GSEA result from its FROZEN on-disk artifacts
alone, in a brand-new process -- per `review:20261012_220407` F2's
instruction for "a command that replays from frozen artifacts in a new
process, never calling the synthetic ranking generator, comparing all
scientific result fields".

    python code/gsea_replay.py --from <original/gsea_out_dir> --out <new/empty/dir> \
        [--cache-dir <cache_dir_override>] [--allow-identity-drift]

Reads `RANKING_SNAPSHOT.json` and `GSEA_RESULT.json` from `--from` (NOT the
example script or any ranking-generating function), reconstructs the exact
`GeneRanking` and the exact run `params` that were recorded, re-runs
`LocalPrerankGSEAProvider` in `--out`, and compares every scientific field
(status, es, nes, p_value, fdr, fwer_p_value, input_set_size,
intersection_size, leading_edge_genes, excluded_gene_sets) between the
original and the replay. Exits 0 and prints "REPLAY MATCH" only if every
compared field is identical AND the replay is bound to the SAME inputs the
original run used; otherwise exits non-zero with the first mismatch.

Per `review:20261013_114111` R1: a numerically-identical result is
NOT sufficient evidence of a replay of the same analysis. Before ever
calling the backend, this module strictly verifies:

  1. The reconstructed ranking's semantic hash and the snapshot FILE's own
     hash both match what the original `GSEA_RESULT.json` itself recorded
     (`read_ranking_snapshot`'s internal self-check only proves the snapshot
     file is internally consistent -- it does NOT prove it is the SAME
     snapshot the original result was computed from; a legitimately
     re-hashed but DIFFERENT snapshot would otherwise pass unnoticed).
  2. The gene-set cache at `--cache-dir` (or the one derived from the
     original result) has the SAME collection name, content hash, species,
     namespace, and source version the original result recorded -- a cache
     directory may be copied/relocated for replay, but its semantic content
     must be the one the original run actually used, never a different
     collection/version that happens to produce the same numbers for one
     term.
  3. `original["params"]`'s required numeric keys are present with the
     correct type (enforced by `gsea_runner.validate_gsea_result`, which
     `read_gsea_result` already calls) -- a missing key is never silently
     defaulted here.

By default, a source-code-revision or dependency-version mismatch between
the original run and the CURRENT environment is also treated as a replay
failure (never silently tolerated as if numbers alone were proof of
identity). `--allow-identity-drift` is an explicit, differently-labeled
CROSS-VERSION COMPARISON mode ("REPLAY COMPARISON", never "REPLAY MATCH")
for the separate, deliberate case of checking whether two different
environments produce the same numbers -- this label now reflects what the
CALLER asked for, not merely whether a drift happened to be observed (see
R4 below).

Per `review:20261013_135634` R3: this strict replay supports ONLY
an original result with `status="complete"`. A probe found that for
`status="no_overlap"` (and any other non-"complete" status), the cache
identity check above was skipped entirely, so a cache swapped for a
different collection/version still produced a "REPLAY MATCH" as long as the
(zero-overlap) numeric outcome agreed. Replaying a "failed", "no_overlap",
or "not_run" original is refused outright, before the ranking snapshot is
even read or the backend is ever called -- the original record is left
exactly as it was (never rewritten, never re-attempted as if it were
"complete"). A future round that wants to support replaying those statuses
must bind the SAME full identity/parameter contract to them, not merely
delete this restriction.

Per `review:20261013_135634` R4: the printed label now reflects
which mode the CALLER explicitly asked for (`--allow-identity-drift` or
not), never a value inferred from whether a drift happened to be observed.
Passing `--allow-identity-drift` against two identical environments still
prints "REPLAY COMPARISON" -- the flag means "this is a cross-version
comparison run", not "a drift was detected".
"""
from __future__ import annotations

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import gsea_runner as R  # noqa: E402
from gene_set_cache import load_cached_gene_sets  # noqa: E402

_NUMERIC_PARAM_KEYS = ("min_set_size", "max_set_size", "random_seed",
                       "permutation_num", "threads", "weight")

_COMPARED_ROW_FIELDS = ("term", "es", "nes", "p_value", "fdr", "fwer_p_value",
                        "input_set_size", "intersection_size", "leading_edge_genes")

#: Per R1: the cache's full semantic identity (collection/content/species/
#: namespace/version), not just its collection name and content hash.
_CACHE_IDENTITY_FIELDS = ("gene_set_collection_name", "gene_set_content_sha256",
                          "gene_set_species", "gene_set_namespace",
                          "gene_set_source_version")


def _check_ranking_identity_or_raise(ranking, snapshot_path: str, original: dict) -> None:
    """Extracted so it can be monkeypatched and mutation-killed independently
    of the rest of `replay()` -- see
    test_mutation_kill_replay_ranking_identity_check in
    test_scfm_gsea_runner.py.
    """
    actual_ranking_sha = R._ranking_sha256(ranking)
    if actual_ranking_sha != original.get("ranking_sha256"):
        raise SystemExit(
            "replay refused: the ranking reconstructed from "
            f"{snapshot_path} has semantic hash {actual_ranking_sha!r}, which "
            "does not match the original result's recorded ranking_sha256 "
            f"{original.get('ranking_sha256')!r} -- this is not a replay of "
            "the same analysis, even if the numeric scores happen to be "
            "identical")
    actual_snapshot_file_sha = R.sha256_file(snapshot_path)
    if actual_snapshot_file_sha != original.get("ranking_snapshot_sha256"):
        raise SystemExit(
            f"replay refused: {snapshot_path} has file hash "
            f"{actual_snapshot_file_sha!r}, which does not match the "
            "original result's recorded ranking_snapshot_sha256 "
            f"{original.get('ranking_snapshot_sha256')!r}")


def _check_cache_identity_or_raise(cache_manifest: dict, original: dict, cache_dir: str) -> None:
    """Extracted so it can be monkeypatched and mutation-killed independently
    -- see test_mutation_kill_replay_cache_identity_check.
    """
    actual = {
        "gene_set_collection_name": cache_manifest["collection_name"],
        "gene_set_content_sha256": cache_manifest["cached_gene_sets_sha256"],
        "gene_set_species": cache_manifest["species"],
        "gene_set_namespace": cache_manifest["gene_id_namespace"],
        "gene_set_source_version": cache_manifest["source_version"],
    }
    for key in _CACHE_IDENTITY_FIELDS:
        if original.get(key) != actual[key]:
            raise SystemExit(
                f"replay refused: --cache-dir {cache_dir!r} has {key}="
                f"{actual[key]!r}, which does not match the original "
                f"result's recorded {key}={original.get(key)!r} -- a cache "
                "directory may be copied/relocated for replay, but its "
                "collection/version/species/namespace/content must be the "
                "SAME one the original run used")


#: Per R3: this strict replay only supports a "complete" original. A
#: "failed"/"no_overlap"/"not_run" original is refused before the ranking
#: snapshot is read or the backend is called -- never silently re-attempted
#: as if it had been "complete".
_REPLAYABLE_STATUSES = ("complete",)


def replay(from_dir: str, out_dir: str, cache_dir: str | None = None,
           *, allow_identity_drift: bool = False) -> dict:
    result_path = os.path.join(from_dir, "GSEA_RESULT.json")
    snapshot_path = os.path.join(from_dir, "RANKING_SNAPSHOT.json")
    original = R.read_gsea_result(result_path)

    if original["status"] not in _REPLAYABLE_STATUSES:
        raise SystemExit(
            f"replay refused: original result at {result_path!r} has "
            f"status={original['status']!r}, and this strict replay only "
            f"supports {_REPLAYABLE_STATUSES!r} -- the original record is "
            "left exactly as it was (not reread as a snapshot, not rerun "
            "through the backend); a 'failed'/'no_overlap'/'not_run' result "
            "has nothing a numeric replay can meaningfully re-verify by "
            "itself")

    ranking = R.read_ranking_snapshot(snapshot_path)

    _check_ranking_identity_or_raise(ranking, snapshot_path, original)

    if cache_dir is None:
        if not original.get("gene_set_cache_path"):
            raise SystemExit(
                "original result has no gene_set_cache_path recorded (likely "
                "status != 'complete') -- pass --cache-dir explicitly")
        cache_dir = os.path.dirname(original["gene_set_cache_path"])

    _, cache_manifest = load_cached_gene_sets(cache_dir)
    _check_cache_identity_or_raise(cache_manifest, original, cache_dir)

    # Per R1 point 3: original["params"]'s required numeric keys were already
    # enforced present-with-correct-type by validate_gsea_result inside
    # read_gsea_result above (for a 'complete' original) -- a missing key is
    # a read-time rejection, never a replay-time silent default.
    run_params = {k: v for k, v in original["params"].items() if k in _NUMERIC_PARAM_KEYS}
    provider = R.LocalPrerankGSEAProvider()
    replayed = provider.run(ranking, cache_dir, out_dir=out_dir, **run_params)
    replayed_dict = replayed.to_dict()

    mismatches = []
    identity_drifts = []
    if replayed_dict["status"] != original["status"]:
        mismatches.append(f"status: original={original['status']!r} "
                          f"replay={replayed_dict['status']!r}")

    if replayed_dict.get("source_code_revision_or_hash") != \
            original.get("source_code_revision_or_hash"):
        identity_drifts.append(
            "source_code_revision_or_hash: original="
            f"{original.get('source_code_revision_or_hash')!r} replay="
            f"{replayed_dict.get('source_code_revision_or_hash')!r}")
    if replayed_dict.get("dependency_versions") != original.get("dependency_versions"):
        identity_drifts.append(
            f"dependency_versions: original={original.get('dependency_versions')!r} "
            f"replay={replayed_dict.get('dependency_versions')!r}")
    # Per R1: a source/dependency drift fails a STRICT replay by default --
    # never silently accepted just because the numbers downstream happen to
    # agree. --allow-identity-drift is the separate, explicitly-labeled
    # cross-version comparison mode.
    if identity_drifts and not allow_identity_drift:
        mismatches.extend("STRICT IDENTITY: " + d for d in identity_drifts)

    orig_rows = {row["term"]: row for row in original.get("results_table", [])}
    replay_rows = {row["term"]: row for row in replayed_dict.get("results_table", [])}
    if set(orig_rows) != set(replay_rows):
        mismatches.append(f"result term sets differ: original={sorted(orig_rows)} "
                          f"replay={sorted(replay_rows)}")
    for term in sorted(set(orig_rows) & set(replay_rows)):
        o, r = orig_rows[term], replay_rows[term]
        for field in _COMPARED_ROW_FIELDS:
            if o.get(field) != r.get(field):
                mismatches.append(f"{term}.{field}: original={o.get(field)!r} "
                                  f"replay={r.get(field)!r}")

    orig_excluded = {(e["term"], e["reason"]) for e in original.get("excluded_gene_sets", [])}
    replay_excluded = {(e["term"], e["reason"])
                       for e in replayed_dict.get("excluded_gene_sets", [])}
    if orig_excluded != replay_excluded:
        mismatches.append(
            f"excluded_gene_sets differ: original={sorted(orig_excluded)} "
            f"replay={sorted(replay_excluded)}")

    return {"original": original, "replayed": replayed_dict, "mismatches": mismatches,
            "identity_drifts": identity_drifts, "allow_identity_drift": allow_identity_drift}


def _replay_label_and_extra(outcome: dict) -> tuple[str, str]:
    """Extracted so it can be independently unit-tested and mutation-killed
    -- see test_mutation_kill_replay_label_reflects_requested_mode in
    test_scfm_gsea_runner.py. Per R4: the label reflects which mode the
    CALLER explicitly asked for (`allow_identity_drift` or not), never
    whether a drift happened to be observed.
    """
    label = "REPLAY COMPARISON" if outcome["allow_identity_drift"] else "REPLAY MATCH"
    extra = ""
    if outcome["identity_drifts"]:
        extra = " (identity drift tolerated: " + "; ".join(outcome["identity_drifts"]) + ")"
    return label, extra


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--from", dest="from_dir", required=True,
                    help="directory containing the original GSEA_RESULT.json "
                         "and RANKING_SNAPSHOT.json")
    ap.add_argument("--out", required=True, help="new/empty directory for the replay run")
    ap.add_argument("--cache-dir", default=None,
                    help="override the gene-set cache directory (default: derived "
                         "from the original result's recorded cache path); must "
                         "point at the SAME collection/version the original run "
                         "used, a relocated copy of it is fine")
    ap.add_argument("--allow-identity-drift", action="store_true",
                    help="Explicit CROSS-VERSION COMPARISON mode: tolerate a "
                         "source-code-revision or dependency-version mismatch "
                         "between the original run and this environment, and "
                         "compare numeric results only. Reported as 'REPLAY "
                         "COMPARISON', never 'REPLAY MATCH', even when every "
                         "number agrees.")
    args = ap.parse_args()

    outcome = replay(args.from_dir, args.out, args.cache_dir,
                     allow_identity_drift=args.allow_identity_drift)
    if outcome["mismatches"]:
        print("REPLAY MISMATCH:")
        for m in outcome["mismatches"]:
            print(f"  - {m}")
        raise SystemExit(1)
    label, extra = _replay_label_and_extra(outcome)
    print(f"{label}: status={outcome['replayed']['status']!r}, "
          f"{len(outcome['replayed'].get('results_table', []))} term(s) compared, "
          f"all scientific fields identical.{extra}")


if __name__ == "__main__":
    main()
