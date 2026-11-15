"""Local preranked-GSEA backend for the scFM discovery framework.

Implements `enrichment_providers.GSEAProvider` via the installed `gseapy.prerank`
(version resolved and recorded at call time -- never assumed from online docs, per
`review:20261012_171916`'s own instruction). `gene_sets` is ALWAYS passed to
gseapy as an in-memory dict loaded from our own hash-verified local cache
(`gene_set_cache.load_cached_gene_sets`) -- never a path string or a library name,
which are the two call shapes gseapy could interpret as a request to fetch from
Enrichr. No network call is made anywhere in this module.

Distinguishes four execution states, never conflating any of them (DEVELOPMENT.md: a
failure must never be written as "no significant pathways"):
  - "not_run": never attempted (the caller's own concern, not produced here).
  - "failed": a structural problem (invalid ranking, species/namespace mismatch,
    cache hash mismatch, no gene set survives size filtering against the exact
    ranking intersection, an identity-consistency violation, or an unexpected
    backend error) stopped the run.
  - "no_overlap": the ranking shares zero genes with the ENTIRE gene-set collection
    -- distinct from "complete" with a non-significant per-term result.
  - "complete": ran end to end and produced a results table, whatever the p-values.

Per `review:20261012_220407` F1: the installed gseapy (confirmed via
`grep` of the real installed `gseapy/base.py`) silently uppercases gene IDs for
matching purposes when the ranking is not mostly-uppercase but the gene sets
sampled are (`Base.load_gmt`'s `_gene_toupper` branch) -- an implicit identity
rewrite this framework's gene-ID contract does not allow by default. Every gene
set is pre-filtered to its EXACT, case-sensitive intersection with the ranking
before being handed to gseapy, which eliminates CROSS-SET leakage (a gene set
being credited with IDs it never literally contained). Per
`review:20261013_114111` C1: this does NOT by itself resolve a case
COLLISION within the ranking's own distinct gene IDs (e.g. `ANCHOR` and
`Anchor` both present as separate ranking entries) -- gseapy's internal
uppercase lookup can still conflate two still-distinct IDs in that synthetic
boundary case. A post-hoc check additionally verifies gseapy's own reported
leading-edge genes never contain anything outside the exact, case-sensitive
intersection computed above, and fails closed (never a silently wrong
"complete") if they ever do -- this is the mechanism that actually catches a
case collision, not the pre-filter alone.

Per `review:20261013_114111` R1: a numerically-identical replay is
not evidence of a replay of the SAME analysis unless the ranking, the
gene-set cache's semantic identity, and the source/dependency identity are
all bound and compared against the original run's own recorded values --
see `gsea_replay.py`.
"""
from __future__ import annotations

import json
import math
import os
import platform
import sys
import time
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import gseapy  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from artifact_schema import (find_nonfinite_floats, resolve_source_revision,  # noqa: E402
                             sha256_file, sha256_text)
from enrichment_providers import GSEAProvider, GSEARunResult  # noqa: E402
from gene_set_cache import load_cached_gene_sets  # noqa: E402

RESULT_SCHEMA = "scfm-gsea-result/1"
RANKING_SNAPSHOT_SCHEMA = "scfm-gsea-ranking-snapshot/1"

#: This module's own source-identity binding (per F2) -- never just gsea_runner.py
#: alone, since gene_set_cache.py and enrichment_providers.py both affect the
#: numeric result (cache loading, ranking validation) just as directly.
#: `artifact_schema.py` is included per `review:20261013_114111` R1:
#: this module imports and calls it directly (find_nonfinite_floats,
#: resolve_source_revision, sha256_file/sha256_text) for every run, so it is
#: just as numerically-relevant as the other three and was missing from this
#: tuple despite being imported below.
SOURCE_MODULES = ("gsea_runner.py", "gene_set_cache.py", "enrichment_providers.py",
                  "artifact_schema.py")

_VALID_STATUSES = ("complete", "failed", "no_overlap", "not_run")


def _reject_nonfinite_constant(token: str):
    raise ValueError(f"non-finite JSON constant {token!r} is not accepted")


def _ranking_sha256(ranking) -> str:
    payload = json.dumps({
        "species": ranking.species,
        "gene_id_namespace": ranking.gene_id_namespace,
        "ranking_definition": ranking.ranking_definition,
        "genes": list(ranking.genes),
        "scores": list(ranking.scores),
    }, sort_keys=True)
    return sha256_text(payload)


def _dependency_versions() -> dict:
    return {
        "gseapy": gseapy.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "python": platform.python_version(),
    }


def write_ranking_snapshot(out_dir: str, ranking, ranking_sha256: str) -> tuple[str, str]:
    """Freezes the FULL ranking -- exact gene/score order, as actually sent to
    the backend -- to disk, per `review:20261012_220407` F2: a hash
    alone verifies content, it cannot RECOVER it. Called only AFTER
    `ranking.validate()` has already passed (per
    `review:20261013_114111` R2): a ranking with a non-finite score
    cannot be persisted here without relaxing `allow_nan=False`, which this
    framework never does (DEVELOPMENT.md: never smuggle NaN/Inf into a frozen
    artifact) -- an invalid ranking gets a clean `status="failed"` result
    with no snapshot and no leftover `.inprogress` file instead. Returns
    (path, file_sha256).
    """
    n_tied = sum(c for c in Counter(ranking.scores).values() if c > 1)
    payload = {
        "schema": RANKING_SNAPSHOT_SCHEMA,
        "species": ranking.species,
        "gene_id_namespace": ranking.gene_id_namespace,
        "ranking_definition": ranking.ranking_definition,
        "genes": list(ranking.genes),
        "scores": list(ranking.scores),
        "ranking_sha256": ranking_sha256,
        "n_genes": len(ranking.genes),
        "n_tied_score_values": n_tied,
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    path = os.path.join(out_dir, "RANKING_SNAPSHOT.json")
    tmp = path + ".inprogress"
    with open(tmp, "w") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True, allow_nan=False)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)
    return path, sha256_file(path)


def read_ranking_snapshot(path: str):
    """Reconstructs a `GeneRanking` from a frozen snapshot file -- used by
    `gsea_replay.py` so a replay never calls back into whatever generated the
    original ranking (per F2's "no replay calling the synthetic generator").
    """
    from enrichment_providers import GeneRanking
    if not os.path.exists(path):
        raise SystemExit(f"ranking snapshot not found: {path}")
    with open(path) as fh:
        payload = json.load(fh)
    if payload.get("schema") != RANKING_SNAPSHOT_SCHEMA:
        raise SystemExit(
            f"ranking snapshot schema {payload.get('schema')!r} != "
            f"expected {RANKING_SNAPSHOT_SCHEMA!r}")
    problems = find_nonfinite_floats(payload)
    if problems:
        raise SystemExit(f"ranking snapshot at {path} has non-finite value(s): {problems}")
    ranking = GeneRanking(
        species=payload["species"], gene_id_namespace=payload["gene_id_namespace"],
        ranking_definition=payload["ranking_definition"],
        genes=tuple(payload["genes"]), scores=tuple(payload["scores"]))
    actual_sha = _ranking_sha256(ranking)
    if actual_sha != payload.get("ranking_sha256"):
        raise SystemExit(
            f"ranking snapshot at {path} failed self-check: recomputed hash "
            f"{actual_sha} != recorded {payload.get('ranking_sha256')}")
    return ranking


def _exact_case_sensitive_filter(gene_sets: dict, ranking_genes) -> dict:
    """Pre-filters every gene set to its EXACT, case-sensitive intersection
    with the ranking's own gene IDs, per `review:20261012_220407`
    F1. Extracted into its own function specifically so it can be
    monkeypatched and mutation-killed (see
    test_mutation_kill_f1_exact_case_prefilter) -- without this filter,
    the installed gseapy's own implicit case-folding heuristic (confirmed
    in its real `base.py`) can silently report a different intersection
    than the one this framework's gene-ID contract promises.
    """
    ranking_gene_set = set(ranking_genes)
    return {term: sorted(set(genes) & ranking_gene_set)
            for term, genes in gene_sets.items()}


def _find_identity_mismatches(ranking, cache_manifest: dict) -> list[str]:
    """Extracted so it can be independently unit-tested and mutation-killed
    (see test_scfm_gsea_runner.py's test_mutation_kill_species_mismatch_check)
    -- never collapse species and namespace mismatches into one vague message,
    and never silently remap either side.
    """
    mismatches = []
    if cache_manifest["species"] != ranking.species:
        mismatches.append(
            f"species mismatch: ranking={ranking.species!r} vs "
            f"gene_set={cache_manifest['species']!r}")
    if cache_manifest["gene_id_namespace"] != ranking.gene_id_namespace:
        mismatches.append(
            f"gene_id_namespace mismatch: ranking={ranking.gene_id_namespace!r} "
            f"vs gene_set={cache_manifest['gene_id_namespace']!r}")
    return mismatches


#: Per `review:20261013_114111` R1: a 'complete' result's `params`
#: must carry every numeric key that could change the backend's output, with
#: the right type and a sane range -- not merely be "a non-empty dict" (a
#: probe found `{"unrelated": 1}` accepted by the previous check). A missing
#: key is REJECTED outright; `gsea_replay.py` must never fall back to a
#: default for a key the original record does not have.
_REQUIRED_NUMERIC_PARAMS = {
    "min_set_size": ((int,), lambda v: v >= 1),
    "max_set_size": ((int,), lambda v: v >= 1),
    "random_seed": ((int,), lambda v: True),
    "permutation_num": ((int,), lambda v: v >= 1),
    "threads": ((int,), lambda v: v >= 1),
    "weight": ((int, float), lambda v: math.isfinite(v)),
}


def _validate_complete_params(params) -> list[str]:
    if not isinstance(params, dict) or not params:
        return ["status='complete' requires a non-empty 'params' dict"]
    problems = []
    for key, (types, in_range) in _REQUIRED_NUMERIC_PARAMS.items():
        if key not in params:
            problems.append(f"status='complete' requires params[{key!r}]")
            continue
        v = params[key]
        if isinstance(v, bool) or not isinstance(v, types):
            problems.append(
                f"params[{key!r}]={v!r} must be of type {types}, got {type(v).__name__}")
        elif not in_range(v):
            problems.append(f"params[{key!r}]={v!r} is out of the accepted range")
    return problems


def validate_gsea_result(payload: dict) -> list[str]:
    """Shared by `_write_result_atomic` (write) and `read_gsea_result` (read),
    per `review:20261012_220407` F3: a result with only a `schema`
    and `status="complete"` and nothing else, an impossible `p_value=2.0`, or
    a non-finite value anywhere outside `results_table` must be rejected at
    BOTH boundaries -- not just the narrower checks the previous round had
    (schema, status enum, and a finiteness walk restricted to results_table).
    """
    problems = []
    if not isinstance(payload, dict):
        return [f"GSEA result must be a dict, got {type(payload).__name__}"]
    if payload.get("schema") != RESULT_SCHEMA:
        problems.append(f"schema {payload.get('schema')!r} != expected {RESULT_SCHEMA!r}")
    status = payload.get("status")
    if status not in _VALID_STATUSES:
        problems.append(f"status {status!r} not in {_VALID_STATUSES}")
    problems.extend(find_nonfinite_floats(payload))

    if status != "complete" and not payload.get("reason"):
        problems.append(f"status={status!r} requires a non-empty 'reason'")

    if status == "complete":
        for key in ("ranking_sha256", "backend", "backend_version",
                    "gene_set_collection_name", "gene_set_cache_path",
                    "gene_set_content_sha256", "gene_set_species",
                    "gene_set_namespace", "gene_set_source_version",
                    "ranking_snapshot_path", "ranking_snapshot_sha256",
                    "source_code_revision_or_hash"):
            v = payload.get(key)
            # Per review:20261013_114111 R1: a truthy NON-string
            # (e.g. the literal `True`) must not satisfy an identity field's
            # contract -- a probe found `source_code_revision_or_hash=True`
            # accepted under the previous `not payload.get(key)` check.
            if not isinstance(v, str) or not v:
                problems.append(
                    f"status='complete' requires a non-empty string {key!r}, got {v!r}")
        problems.extend(_validate_complete_params(payload.get("params")))
        dep = payload.get("dependency_versions")
        if not isinstance(dep, dict) or not dep:
            problems.append(
                "status='complete' requires a non-empty 'dependency_versions' dict")
        elif not all(isinstance(k, str) and isinstance(v, str) for k, v in dep.items()):
            problems.append(
                f"dependency_versions keys and values must all be strings, got {dep!r}")
        rows = payload.get("results_table")
        if not isinstance(rows, list) or not rows:
            problems.append("status='complete' requires a non-empty 'results_table' list")
        else:
            for i, row in enumerate(rows):
                problems.extend(_validate_result_row(row, i))
    return problems


def _validate_result_row(row: dict, i: int) -> list[str]:
    problems = []
    if not isinstance(row, dict):
        return [f"results_table[{i}] must be a dict, got {type(row).__name__}"]
    if not isinstance(row.get("term"), str) or not row["term"]:
        problems.append(f"results_table[{i}].term must be a non-empty string")
    for field in ("es", "nes"):
        v = row.get(field)
        if not isinstance(v, (int, float)) or isinstance(v, bool) or not math.isfinite(v):
            problems.append(f"results_table[{i}].{field} must be a finite number, got {v!r}")
    # ES is a GSEA running-sum statistic bounded in [-1, 1] by definition;
    # NES is normalized against a permutation null and is NOT bounded that
    # way (NES magnitudes routinely exceed 1, as this framework's own
    # synthetic example demonstrates) -- only ES gets a range check, with a
    # small epsilon for floating-point rounding, per the review's own note
    # not to add a bound intolerant of that rounding.
    es = row.get("es")
    if isinstance(es, (int, float)) and not isinstance(es, bool) and math.isfinite(es):
        if not (-1.0 - 1e-6 <= es <= 1.0 + 1e-6):
            problems.append(f"results_table[{i}].es={es!r} outside [-1, 1] (+/- rounding)")
    for field in ("p_value", "fdr", "fwer_p_value"):
        v = row.get(field)
        if not isinstance(v, (int, float)) or isinstance(v, bool) or not math.isfinite(v):
            problems.append(f"results_table[{i}].{field} must be a finite number, got {v!r}")
        elif not (-1e-9 <= v <= 1.0 + 1e-9):
            # A legitimate p-value/FDR of exactly 0 must NOT be rejected.
            problems.append(f"results_table[{i}].{field}={v!r} outside [0, 1] (+/- rounding)")
    input_size = row.get("input_set_size")
    inter_size = row.get("intersection_size")
    if not isinstance(input_size, int) or isinstance(input_size, bool) or input_size < 0:
        problems.append(f"results_table[{i}].input_set_size must be a non-negative int")
    if not isinstance(inter_size, int) or isinstance(inter_size, bool) or inter_size < 0:
        problems.append(f"results_table[{i}].intersection_size must be a non-negative int")
    elif isinstance(input_size, int) and inter_size > input_size:
        problems.append(
            f"results_table[{i}]: intersection_size={inter_size} > input_set_size={input_size}")
    leading_edge = row.get("leading_edge_genes")
    if not isinstance(leading_edge, list) or not all(isinstance(g, str) for g in leading_edge):
        problems.append(f"results_table[{i}].leading_edge_genes must be a list of strings")
    elif isinstance(inter_size, int) and len(leading_edge) > inter_size:
        problems.append(
            f"results_table[{i}]: {len(leading_edge)} leading_edge_genes > "
            f"intersection_size={inter_size}")
    return problems


def _write_result_atomic(path: str, result: GSEARunResult) -> None:
    payload = result.to_dict()
    payload["schema"] = RESULT_SCHEMA
    payload["generated_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    problems = validate_gsea_result(payload)
    if problems:
        raise ValueError(f"refusing to write an invalid GSEA result: {problems}")
    tmp = path + ".inprogress"
    with open(tmp, "w") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True, allow_nan=False)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def read_gsea_result(path: str) -> dict:
    """Strict re-read: rejects an illegal `status`, a missing `schema`, a
    non-finite value ANYWHERE in the payload, an impossible probability
    field, or a "complete" status missing its required identity/params --
    a hand-tampered or merely incomplete on-disk file is caught here, not
    only at write time.
    """
    if not os.path.exists(path):
        raise SystemExit(f"GSEA result not found: {path}")
    with open(path) as fh:
        try:
            payload = json.load(fh, parse_constant=_reject_nonfinite_constant)
        except ValueError as e:
            raise SystemExit(f"GSEA result at {path} failed to parse: {e}")
    problems = validate_gsea_result(payload)
    if problems:
        raise SystemExit(f"GSEA result at {path} failed validation: {problems}")
    return payload


class LocalPrerankGSEAProvider(GSEAProvider):
    """Real `GSEAProvider` implementation. See module docstring for the
    execution-state vocabulary. Genuinely unexpected backend errors (a real
    `gseapy` exception unrelated to the structural checks below) are recorded
    as a `failed` result for audit purposes and then RE-RAISED -- DEVELOPMENT.md's
    standing rule against swallowing failures into a misleading "complete",
    and per F5, never silently swallowed into no persisted record at all.
    """

    def run(self, ranking, gene_set_cache_dir, *, out_dir, min_set_size=10,
            max_set_size=500, random_seed=0, permutation_num=100, threads=1,
            weight=1.0) -> GSEARunResult:
        if os.path.exists(out_dir) and os.listdir(out_dir):
            raise SystemExit(f"{out_dir} must be new or empty")
        os.makedirs(out_dir, exist_ok=True)

        r_sha = _ranking_sha256(ranking)
        source_revision = resolve_source_revision(HERE, SOURCE_MODULES)
        dep_versions = _dependency_versions()
        n_tied = sum(c for c in Counter(ranking.scores).values() if c > 1)

        params = {
            "min_set_size": min_set_size, "max_set_size": max_set_size,
            "random_seed": random_seed, "permutation_num": permutation_num,
            "threads": threads, "weight": weight,
            "n_tied_score_values": n_tied,
            "ties_policy": "ties are broken by whatever relative order the tied "
                           "genes appear in in the FROZEN ranking_snapshot's "
                           "gene/score arrays, via gseapy's own internal (not "
                           "independently re-implemented) numpy argsort over "
                           "that exact array -- deterministic for a fixed input "
                           "order, NOT guaranteed invariant if that order is "
                           "permuted. This is a computational convention "
                           "recorded for reproducibility, not a biological "
                           "priority among tied genes.",
            "ranking_direction": "descending (highest score = most enriched "
                                 "at the top), per gseapy.prerank's own "
                                 "default ascending=False",
            "p_value_definition": "gseapy's own 'NOM p-val' (permutation-based "
                                  "nominal p-value, not independently recomputed)",
            "fdr_definition": "gseapy's own 'FDR q-val' (not independently recomputed)",
        }

        def finish(result: GSEARunResult) -> GSEARunResult:
            _write_result_atomic(os.path.join(out_dir, "GSEA_RESULT.json"), result)
            return result

        # Per review:20261013_114111 R2: validate the ranking BEFORE
        # ever writing a snapshot. A ranking with a non-finite score cannot be
        # persisted under allow_nan=False -- for an invalid ranking there is
        # no legitimate snapshot to write at all, so this branch returns a
        # clean 'failed' record (no ranking_snapshot_path/sha256, no
        # leftover '.inprogress' file) instead of letting write_ranking_snapshot
        # raise mid-write.
        problems = ranking.validate()
        if problems:
            return finish(GSEARunResult(
                status="failed",
                reason=f"invalid ranking, no snapshot was written: {problems}",
                ranking_sha256=r_sha, backend="gseapy", backend_version=gseapy.__version__,
                params=params, source_code_revision_or_hash=source_revision,
                dependency_versions=dep_versions))

        snapshot_path, snapshot_sha256 = write_ranking_snapshot(out_dir, ranking, r_sha)
        common_always = dict(
            ranking_sha256=r_sha, backend="gseapy", backend_version=gseapy.__version__,
            params=params, ranking_snapshot_path=snapshot_path,
            ranking_snapshot_sha256=snapshot_sha256,
            source_code_revision_or_hash=source_revision,
            dependency_versions=dep_versions)

        gene_sets, cache_manifest = load_cached_gene_sets(gene_set_cache_dir)

        mismatches = _find_identity_mismatches(ranking, cache_manifest)
        common = dict(
            common_always,
            gene_set_collection_name=cache_manifest["collection_name"],
            gene_set_cache_path=cache_manifest["cached_gene_sets_path"],
            gene_set_content_sha256=cache_manifest["cached_gene_sets_sha256"],
            # Per review:20261013_114111 R1: the cache's full semantic
            # identity, not just its collection name/content hash -- so a
            # strict replay can refuse a cache swapped for a different
            # species/namespace/version collection.
            gene_set_species=cache_manifest["species"],
            gene_set_namespace=cache_manifest["gene_id_namespace"],
            gene_set_source_version=cache_manifest["source_version"])
        if mismatches:
            return finish(GSEARunResult(status="failed", reason="; ".join(mismatches),
                                        **common))

        # Per F1: pre-filter every gene set to its EXACT, case-sensitive
        # intersection with the ranking before gseapy ever sees it -- this
        # eliminates CROSS-SET leakage (a gene set matching IDs it never
        # literally contained). It does NOT by itself resolve a case
        # COLLISION within the ranking's own distinct gene IDs (see C1 in
        # review:20261013_114111) -- the post-hoc leading-edge
        # check below is what fails closed on that remaining case.
        # input_set_size keeps the ORIGINAL (pre-filter) pathway size for
        # interpretability.
        exact_filtered = _exact_case_sensitive_filter(gene_sets, ranking.genes)
        input_sizes = {term: len(genes) for term, genes in gene_sets.items()}

        all_filtered_genes = {g for genes in exact_filtered.values() for g in genes}
        if not all_filtered_genes:
            return finish(GSEARunResult(
                status="no_overlap",
                reason="ranking shares zero genes (exact case-sensitive match) "
                       "with every gene set in the collection -- distinct from "
                       "a completed run with no significant terms",
                **common))

        # Per F5: replicate gseapy's own min/max-size *acceptance* condition
        # here, BEFORE calling it, using our exact intersection -- so "every
        # gene set was filtered out" becomes a structured, audited 'failed'
        # record with a per-term reason, never an unhandled LookupError and
        # never a silent empty 'complete'.
        included, excluded = {}, []
        for term, filtered_genes in exact_filtered.items():
            n = len(filtered_genes)
            if n == 0:
                reason = "zero_overlap_with_ranking"
            elif n < min_set_size:
                reason = "too_small"
            elif n > max_set_size:
                reason = "too_large"
            elif n >= len(ranking.genes):
                reason = "not_a_strict_subset_of_the_ranking"
            else:
                included[term] = filtered_genes
                continue
            excluded.append({"term": term, "input_set_size": input_sizes[term],
                             "intersection_size": n, "reason": reason})

        if not included:
            n_small = sum(1 for e in excluded if e["reason"] == "too_small")
            n_large = sum(1 for e in excluded if e["reason"] == "too_large")
            n_zero = sum(1 for e in excluded if e["reason"] == "zero_overlap_with_ranking")
            return finish(GSEARunResult(
                status="failed",
                reason=(f"no gene set passed size filtering against the exact "
                        f"ranking intersection (min_set_size={min_set_size}, "
                        f"max_set_size={max_set_size}): {n_zero} had zero "
                        f"overlap, {n_small} too_small, {n_large} too_large "
                        f"out of {len(exact_filtered)} gene set(s)"),
                excluded_gene_sets=tuple(excluded), **common))

        try:
            rnk = pd.Series(dict(zip(ranking.genes, ranking.scores)))
            pre = gseapy.prerank(
                rnk=rnk, gene_sets=included, min_size=min_set_size,
                max_size=max_set_size, permutation_num=permutation_num,
                threads=threads, weight=weight, seed=random_seed,
                outdir=None, no_plot=True, verbose=False)
        except Exception as e:
            finish(GSEARunResult(
                status="failed",
                reason=f"backend raised an unexpected exception: {e!r}",
                excluded_gene_sets=tuple(excluded), **common))
            raise

        rows = []
        identity_violations = []
        for _, row in pre.res2d.iterrows():
            term = row["Term"]
            lead_str = row["Lead_genes"] or ""
            leading_edge = [g for g in lead_str.split(";") if g]
            true_intersection = set(included.get(term, []))
            if not set(leading_edge) <= true_intersection:
                identity_violations.append(
                    f"{term}: leading edge {sorted(set(leading_edge) - true_intersection)} "
                    "not in the exact case-sensitive ranking intersection")
            rows.append({
                "term": term,
                "es": float(row["ES"]),
                "nes": float(row["NES"]),
                "p_value": float(row["NOM p-val"]),
                "fdr": float(row["FDR q-val"]),
                "fwer_p_value": float(row["FWER p-val"]),
                "input_set_size": input_sizes.get(term, len(included.get(term, []))),
                "intersection_size": len(true_intersection),
                "leading_edge_genes": leading_edge,
            })

        if identity_violations:
            # Defense in depth: the F1 pre-filter above should make this
            # unreachable, but if gseapy's own internal matching ever
            # disagrees with our exact-identity definition despite it, fail
            # closed rather than report a self-contradictory "complete".
            return finish(GSEARunResult(
                status="failed",
                reason="backend-reported leading edge violated exact gene "
                       f"identity: {identity_violations}",
                excluded_gene_sets=tuple(excluded), **common))

        return finish(GSEARunResult(
            status="complete", excluded_gene_sets=tuple(excluded),
            results_table=tuple(rows), **common))
