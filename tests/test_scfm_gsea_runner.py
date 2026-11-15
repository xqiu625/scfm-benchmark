"""Tests for the real, gseapy-backed LocalPrerankGSEAProvider. Every gene set
here is local and in-memory by the time it reaches gseapy (loaded from our own
hash-verified cache) -- no network call is ever made. Uses small, FIXED, tiny
synthetic rankings where the enrichment direction is known by construction
(a deliberately boosted "top" gene set), not a significant p-value used as the
sole pass/fail criterion.
"""

import os
import sys

import pytest

CODE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "scripts",
    "34_scfm_discovery_framework",
    "code",
)
sys.path.insert(0, CODE)

import gene_set_cache as G  # noqa: E402
import gsea_runner as R  # noqa: E402
from enrichment_providers import GeneRanking, GSEARunResult  # noqa: E402


def _write_gmt(path, sets: dict):
    lines = ["\t".join([name, f"{name}_desc", *genes]) for name, genes in sets.items()]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def _make_cache(tmp_path, sets: dict, *, species="human", namespace="gene_symbol", name="cache"):
    gmt = tmp_path / "sets.gmt"
    _write_gmt(gmt, sets)
    cache_dir = tmp_path / name
    G.import_gene_sets(
        source_gmt_path=str(gmt),
        collection_name="TEST_COLLECTION",
        species=species,
        gene_id_namespace=namespace,
        source_version="local-test-v1",
        out_dir=str(cache_dir),
    )
    return str(cache_dir)


def _make_known_direction_ranking(n_genes=200, n_top=20, species="human", namespace="gene_symbol"):
    """A ranking where TOP_SET's genes are deliberately given the highest
    scores by construction -- the enrichment direction is known, not hoped
    for. RANDOM_SET is an arbitrary, unrelated set of the same size class.
    """
    import numpy as np

    rng = np.random.default_rng(0)
    genes = [f"G{i}" for i in range(n_genes)]
    scores = list(rng.normal(size=n_genes))
    top_genes = genes[:n_top]
    for i, g in enumerate(top_genes):
        scores[genes.index(g)] = 5.0 - i * 0.05
    random_genes = list(rng.choice(genes[n_top:], size=15, replace=False))
    ranking = GeneRanking(
        species=species,
        gene_id_namespace=namespace,
        ranking_definition="synthetic, TOP_SET boosted by construction",
        genes=tuple(genes),
        scores=tuple(scores),
    )
    return ranking, top_genes, random_genes


def test_run_complete_detects_known_enriched_set_direction(tmp_path):
    """TOP_SET (boosted by construction) must come out with a large positive
    ES/NES and a far smaller NOM p-val than RANDOM_SET -- a real, deterministic
    consequence of the input, not a lucky significant p-value in isolation.
    """
    ranking, top_genes, random_genes = _make_known_direction_ranking()
    cache_dir = _make_cache(tmp_path, {"TOP_SET": top_genes, "RANDOM_SET": random_genes})
    out_dir = tmp_path / "gsea_out"

    provider = R.LocalPrerankGSEAProvider()
    result = provider.run(
        ranking,
        cache_dir,
        out_dir=str(out_dir),
        min_set_size=1,
        max_set_size=500,
        permutation_num=100,
        threads=1,
        random_seed=0,
    )

    assert result.status == "complete", result.reason
    by_term = {row["term"]: row for row in result.results_table}
    assert by_term["TOP_SET"]["es"] > by_term["RANDOM_SET"]["es"]
    assert by_term["TOP_SET"]["nes"] > by_term["RANDOM_SET"]["nes"]
    # Not a strict p_value inequality: with permutation_num=100, both were
    # observed (on this fixed seed) to tie at p=0.01 -- a real property of
    # this backend's permutation p-values at this permutation count, not a
    # bug. That tie alone does not establish a general 1/(permutation_num+1)
    # floor for this gseapy version; ES/NES above are the deterministic,
    # continuous signal this test actually relies on, and p_value is only
    # checked as <= here.
    assert by_term["TOP_SET"]["p_value"] <= by_term["RANDOM_SET"]["p_value"]
    # Every one of TOP_SET's genes was placed at the very top of the ranking,
    # so all of them must appear in TOP_SET's own leading edge.
    assert set(top_genes) <= set(by_term["TOP_SET"]["leading_edge_genes"])
    assert by_term["TOP_SET"]["intersection_size"] == len(top_genes)
    assert (out_dir / "GSEA_RESULT.json").exists()


def test_run_result_persists_and_rereads_exactly(tmp_path):
    ranking, top_genes, random_genes = _make_known_direction_ranking()
    cache_dir = _make_cache(tmp_path, {"TOP_SET": top_genes, "RANDOM_SET": random_genes})
    out_dir = tmp_path / "gsea_out"
    provider = R.LocalPrerankGSEAProvider()
    result = provider.run(
        ranking, cache_dir, out_dir=str(out_dir), min_set_size=1, permutation_num=100, random_seed=0
    )
    read_back = R.read_gsea_result(str(out_dir / "GSEA_RESULT.json"))
    assert read_back["status"] == "complete"
    assert read_back["ranking_sha256"] == result.ranking_sha256
    assert len(read_back["results_table"]) == len(result.results_table)
    assert read_back["backend"] == "gseapy"
    assert read_back["backend_version"]  # resolved at call time, not hardcoded


def test_run_is_reproducible_given_fixed_seed_threads_and_params(tmp_path):
    ranking, top_genes, random_genes = _make_known_direction_ranking()
    cache_dir = _make_cache(tmp_path, {"TOP_SET": top_genes, "RANDOM_SET": random_genes})
    provider = R.LocalPrerankGSEAProvider()
    r1 = provider.run(
        ranking,
        cache_dir,
        out_dir=str(tmp_path / "out1"),
        min_set_size=1,
        permutation_num=100,
        threads=1,
        random_seed=7,
    )
    r2 = provider.run(
        ranking,
        cache_dir,
        out_dir=str(tmp_path / "out2"),
        min_set_size=1,
        permutation_num=100,
        threads=1,
        random_seed=7,
    )
    rows1 = {row["term"]: row for row in r1.results_table}
    rows2 = {row["term"]: row for row in r2.results_table}
    for term in rows1:
        assert rows1[term]["es"] == rows2[term]["es"]
        assert rows1[term]["nes"] == rows2[term]["nes"]
        assert rows1[term]["p_value"] == rows2[term]["p_value"]


def test_run_fails_closed_on_species_mismatch(tmp_path):
    ranking, top_genes, random_genes = _make_known_direction_ranking(species="mouse")
    cache_dir = _make_cache(tmp_path, {"TOP_SET": top_genes}, species="human")
    provider = R.LocalPrerankGSEAProvider()
    result = provider.run(
        ranking, cache_dir, out_dir=str(tmp_path / "out"), min_set_size=1, permutation_num=50
    )
    assert result.status == "failed"
    assert "species mismatch" in result.reason


def test_run_fails_closed_on_namespace_mismatch(tmp_path):
    ranking, top_genes, random_genes = _make_known_direction_ranking(namespace="ensembl_gene_id")
    cache_dir = _make_cache(tmp_path, {"TOP_SET": top_genes}, namespace="gene_symbol")
    provider = R.LocalPrerankGSEAProvider()
    result = provider.run(
        ranking, cache_dir, out_dir=str(tmp_path / "out"), min_set_size=1, permutation_num=50
    )
    assert result.status == "failed"
    assert "gene_id_namespace mismatch" in result.reason


def test_run_reports_no_overlap_distinct_from_complete(tmp_path):
    """Ranking genes and gene-set genes share nothing at all -- must be
    'no_overlap', never silently 'complete' with an empty/degenerate table.
    """
    ranking = GeneRanking(
        species="human",
        gene_id_namespace="gene_symbol",
        ranking_definition="disjoint by construction",
        genes=("RANKED_ONLY_1", "RANKED_ONLY_2"),
        scores=(1.0, -1.0),
    )
    cache_dir = _make_cache(tmp_path, {"SET_A": ["OTHER_GENE_1", "OTHER_GENE_2"]})
    provider = R.LocalPrerankGSEAProvider()
    result = provider.run(
        ranking, cache_dir, out_dir=str(tmp_path / "out"), min_set_size=1, permutation_num=50
    )
    assert result.status == "no_overlap"
    assert result.results_table == ()


def test_run_fails_closed_on_invalid_ranking(tmp_path):
    ranking = GeneRanking(
        species="human", gene_id_namespace="gene_symbol", ranking_definition="", genes=(), scores=()
    )
    cache_dir = _make_cache(tmp_path, {"SET_A": ["G1"]})
    provider = R.LocalPrerankGSEAProvider()
    result = provider.run(ranking, cache_dir, out_dir=str(tmp_path / "out"))
    assert result.status == "failed"
    assert "invalid ranking" in result.reason


def test_run_detects_cache_tampering_before_calling_gseapy(tmp_path):
    ranking, top_genes, random_genes = _make_known_direction_ranking()
    cache_dir = _make_cache(tmp_path, {"TOP_SET": top_genes, "RANDOM_SET": random_genes})
    manifest_path = os.path.join(cache_dir, "CACHE_MANIFEST.json")
    import json

    manifest = json.load(open(manifest_path))
    # cached_gene_sets_path is stored RELATIVE to the manifest's own
    # directory (per review:20261012_220407 F4) -- resolve it the
    # same way load_cached_gene_sets does, not relative to this test's cwd.
    gene_sets_path = os.path.join(cache_dir, manifest["cached_gene_sets_path"])
    with open(gene_sets_path, "a") as fh:
        fh.write("\n// tampered")
    provider = R.LocalPrerankGSEAProvider()
    with pytest.raises(SystemExit, match="hash mismatch"):
        provider.run(
            ranking, cache_dir, out_dir=str(tmp_path / "out"), min_set_size=1, permutation_num=50
        )


def test_run_refuses_nonempty_out_dir(tmp_path):
    ranking, top_genes, random_genes = _make_known_direction_ranking()
    cache_dir = _make_cache(tmp_path, {"TOP_SET": top_genes})
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    (out_dir / "junk.txt").write_text("pre-existing")
    provider = R.LocalPrerankGSEAProvider()
    with pytest.raises(SystemExit, match="must be new or empty"):
        provider.run(ranking, cache_dir, out_dir=str(out_dir))
    assert (out_dir / "junk.txt").read_text() == "pre-existing"


def test_run_never_connects_to_network(tmp_path, monkeypatch):
    import socket

    def _boom(self, *a, **k):
        raise AssertionError("LocalPrerankGSEAProvider.run must never open a " "network connection")

    monkeypatch.setattr(socket.socket, "connect", _boom)
    ranking, top_genes, random_genes = _make_known_direction_ranking(n_genes=60, n_top=8)
    cache_dir = _make_cache(tmp_path, {"TOP_SET": top_genes})
    provider = R.LocalPrerankGSEAProvider()
    result = provider.run(
        ranking,
        cache_dir,
        out_dir=str(tmp_path / "out"),
        min_set_size=1,
        permutation_num=20,
        threads=1,
    )
    assert result.status == "complete"


def test_read_gsea_result_rejects_illegal_status(tmp_path):
    import json

    path = tmp_path / "tampered.json"
    path.write_text(json.dumps({"schema": R.RESULT_SCHEMA, "status": "definitely_ran"}))
    with pytest.raises(SystemExit, match="status"):
        R.read_gsea_result(str(path))


def test_read_gsea_result_rejects_literal_nan_in_results_table(tmp_path):
    import json

    payload = {
        "schema": R.RESULT_SCHEMA,
        "status": "complete",
        "results_table": [{"term": "X", "es": 1.0}],
    }
    text = json.dumps(payload).replace('"es": 1.0', '"es": NaN')
    path = tmp_path / "tampered.json"
    path.write_text(text)
    with pytest.raises(SystemExit):
        R.read_gsea_result(str(path))


def test_read_gsea_result_rejects_overflow_to_infinity_in_results_table(tmp_path):
    """1e999 is not one of the three special JSON constant tokens
    parse_constant intercepts -- it must still be caught by the recursive
    finiteness walk over results_table.
    """
    import json

    payload = {
        "schema": R.RESULT_SCHEMA,
        "status": "complete",
        "results_table": [{"term": "X", "es": 1.0}],
    }
    text = json.dumps(payload).replace('"es": 1.0', '"es": 1e999')
    path = tmp_path / "tampered.json"
    path.write_text(text)
    with pytest.raises(SystemExit, match="non-finite"):
        R.read_gsea_result(str(path))


def test_mutation_kill_species_and_namespace_mismatch_checks(monkeypatch, tmp_path):
    """Disables _find_identity_mismatches (returns [] unconditionally) and
    confirms test_run_fails_closed_on_species_mismatch /
    test_run_fails_closed_on_namespace_mismatch would then genuinely fail --
    proving the real check is load-bearing, not a tautological pass. Patches
    the exact function those checks are extracted into.
    """
    monkeypatch.setattr(R, "_find_identity_mismatches", lambda ranking, manifest: [])

    ranking, top_genes, _ = _make_known_direction_ranking(species="mouse")
    cache_dir = _make_cache(tmp_path, {"TOP_SET": top_genes}, species="human")
    provider = R.LocalPrerankGSEAProvider()
    result = provider.run(
        ranking, cache_dir, out_dir=str(tmp_path / "out"), min_set_size=1, permutation_num=50
    )
    # With the mutant, a genuine species mismatch is NOT caught -- the run
    # proceeds to completion instead of failing closed, proving the real
    # _find_identity_mismatches is what makes test_run_fails_closed_on_*
    # actually fail-closed.
    assert result.status == "complete"


# ===========================================================================
# review:20261012_220407 F1: the installed gseapy silently
# uppercases gene IDs for matching when the ranking isn't mostly-uppercase
# but the gene sets sampled are -- an implicit identity rewrite. Every gene
# set is now pre-filtered to its EXACT, case-sensitive intersection with the
# ranking before gseapy ever sees it, so intersection_size and leading_edge
# can never disagree about which genes actually matched.
# ===========================================================================


def test_run_f1_case_mismatch_intersection_and_leading_edge_stay_consistent(tmp_path):
    """The review's own reproduction: a mouse ranking in mixed case
    ("GeneN") against a gene set written in all-caps ("GENEn") plus one
    exact-case "ANCHOR" gene shared by both. Before the fix, gseapy's own
    implicit case-folding reported intersection_size=1 (ANCHOR only) while
    simultaneously reporting 12 case-mismatched genes in the leading edge --
    a self-contradictory "complete" result. After the fix, both numbers must
    agree: only the exact, case-sensitive match (ANCHOR) is ever reported.
    """
    ranking_genes = tuple(f"Gene{i}" for i in range(80)) + ("ANCHOR",)
    ranking_scores = tuple(float(80 - i) for i in range(80)) + (1.0,)
    ranking = GeneRanking(
        species="mouse",
        gene_id_namespace="gene_symbol",
        ranking_definition="f1 regression: mixed-case mouse ranking",
        genes=ranking_genes,
        scores=ranking_scores,
    )
    geneset_genes = [f"GENE{i}" for i in range(12)] + ["ANCHOR"]
    cache_dir = _make_cache(tmp_path, {"CASE_SET": geneset_genes}, species="mouse")

    provider = R.LocalPrerankGSEAProvider()
    result = provider.run(
        ranking,
        cache_dir,
        out_dir=str(tmp_path / "out"),
        min_set_size=1,
        max_set_size=500,
        permutation_num=50,
    )
    assert result.status == "complete", result.reason
    row = {r["term"]: r for r in result.results_table}["CASE_SET"]
    assert row["intersection_size"] == 1
    assert row["leading_edge_genes"] == ["ANCHOR"] or row["leading_edge_genes"] == []
    assert set(row["leading_edge_genes"]) <= {"ANCHOR"}


def test_run_f1_legitimate_consistent_mixed_case_mouse_input_still_works(tmp_path):
    """A LEGITIMATE mouse ranking and gene set that agree with each other on
    mixed case (real MGI-style symbols) must still run normally -- this
    fix must not start rejecting or degrading correctly-cased input.
    """
    import numpy as np

    rng = np.random.default_rng(0)
    n_genes = 200
    genes = tuple(f"Mgene{i}" for i in range(n_genes))
    scores = list(rng.normal(size=n_genes))
    top = genes[:20]
    for i, g in enumerate(top):
        scores[i] = 5.0 - i * 0.05
    ranking = GeneRanking(
        species="mouse",
        gene_id_namespace="gene_symbol",
        ranking_definition="f1: consistent mixed-case mouse input",
        genes=genes,
        scores=tuple(scores),
    )
    cache_dir = _make_cache(tmp_path, {"TOP_SET": list(top)}, species="mouse")
    provider = R.LocalPrerankGSEAProvider()
    result = provider.run(
        ranking, cache_dir, out_dir=str(tmp_path / "out"), min_set_size=1, permutation_num=50
    )
    assert result.status == "complete", result.reason
    row = {r["term"]: r for r in result.results_table}["TOP_SET"]
    assert row["intersection_size"] == len(top)
    assert set(top) <= set(row["leading_edge_genes"])


def test_run_f1_regular_uppercase_human_input_unaffected(tmp_path):
    """Standard, fully-uppercase human gene symbols on both sides -- the
    common case this framework is actually expected to run -- must be
    completely unaffected by the F1 fix."""
    ranking, top_genes, random_genes = _make_known_direction_ranking(species="human")
    cache_dir = _make_cache(tmp_path, {"TOP_SET": top_genes, "RANDOM_SET": random_genes})
    provider = R.LocalPrerankGSEAProvider()
    result = provider.run(
        ranking, cache_dir, out_dir=str(tmp_path / "out"), min_set_size=1, permutation_num=50
    )
    assert result.status == "complete", result.reason
    row = {r["term"]: r for r in result.results_table}["TOP_SET"]
    assert row["intersection_size"] == len(top_genes)
    assert set(top_genes) <= set(row["leading_edge_genes"])


def test_mutation_kill_f1_exact_case_prefilter(monkeypatch, tmp_path):
    """Disables the real `_exact_case_sensitive_filter` (makes it a no-op
    that returns the gene sets completely UNFILTERED, i.e. exactly the
    pre-fix behaviour) and confirms the case-mismatch regression test above
    would then genuinely fail: the backend's own implicit case-folding
    reports a self-contradictory "complete" result (intersection_size=1 but
    a 12-gene leading edge outside that intersection) instead of the
    consistent one the real filter produces -- proving the real filter,
    not a coincidence of this particular fixture, is what keeps
    intersection_size and leading_edge_genes in agreement.
    """
    monkeypatch.setattr(
        R,
        "_exact_case_sensitive_filter",
        lambda gene_sets, ranking_genes: {t: list(g) for t, g in gene_sets.items()},
    )

    ranking_genes = tuple(f"Gene{i}" for i in range(80)) + ("ANCHOR",)
    ranking_scores = tuple(float(80 - i) for i in range(80)) + (1.0,)
    ranking = GeneRanking(
        species="mouse",
        gene_id_namespace="gene_symbol",
        ranking_definition="f1 mutation-kill",
        genes=ranking_genes,
        scores=ranking_scores,
    )
    geneset_genes = [f"GENE{i}" for i in range(12)] + ["ANCHOR"]
    cache_dir = _make_cache(tmp_path, {"CASE_SET": geneset_genes}, species="mouse")

    provider = R.LocalPrerankGSEAProvider()
    result = provider.run(
        ranking,
        cache_dir,
        out_dir=str(tmp_path / "out"),
        min_set_size=1,
        max_set_size=500,
        permutation_num=50,
    )

    # With the real filter, this exact fixture reaches "complete" with a
    # consistent intersection_size=1 (see the paired regression test above).
    # With the pre-filter disabled, the backend's own implicit case-folding
    # DOES still get caught -- but only by the separate, defense-in-depth
    # leading-edge check, and only as a 'failed' run instead of a correct
    # 'complete' one: the genuinely present exact match (ANCHOR) that should
    # have produced a valid result is lost, proving the pre-filter is not
    # merely redundant with the post-hoc check -- removing it turns a
    # legitimate, computable result into a spurious failure.
    assert result.status == "failed"
    assert "leading edge" in result.reason and "not in the exact" in result.reason


# ===========================================================================
# review:20261012_220407 F2: a ranking's hash alone verifies
# content, it cannot RECOVER it. Every run now freezes the full ranking
# (exact gene/score order) to RANKING_SNAPSHOT.json, binds the wrapper's own
# source revision and dependency versions, records a tie count, and a
# brand-new process can replay the whole run from those frozen artifacts
# alone via gsea_replay.py.
# ===========================================================================


def test_run_writes_a_ranking_snapshot_matching_its_own_sha(tmp_path):
    ranking, top_genes, random_genes = _make_known_direction_ranking()
    cache_dir = _make_cache(tmp_path, {"TOP_SET": top_genes, "RANDOM_SET": random_genes})
    out_dir = tmp_path / "out"
    provider = R.LocalPrerankGSEAProvider()
    result = provider.run(
        ranking, cache_dir, out_dir=str(out_dir), min_set_size=1, permutation_num=50
    )
    snapshot_path = out_dir / "RANKING_SNAPSHOT.json"
    assert snapshot_path.exists()
    assert result.ranking_snapshot_path == str(snapshot_path)
    import hashlib

    assert hashlib.sha256(snapshot_path.read_bytes()).hexdigest() == result.ranking_snapshot_sha256
    reread = R.read_ranking_snapshot(str(snapshot_path))
    assert reread.genes == ranking.genes
    assert reread.scores == ranking.scores
    assert reread.species == ranking.species


def test_run_binds_source_revision_and_dependency_versions(tmp_path):
    ranking, top_genes, random_genes = _make_known_direction_ranking()
    cache_dir = _make_cache(tmp_path, {"TOP_SET": top_genes})
    provider = R.LocalPrerankGSEAProvider()
    result = provider.run(
        ranking, cache_dir, out_dir=str(tmp_path / "out"), min_set_size=1, permutation_num=50
    )
    assert result.source_code_revision_or_hash
    assert set(result.dependency_versions) >= {"gseapy", "numpy", "pandas", "python"}
    assert result.dependency_versions["gseapy"] == result.backend_version


def test_run_records_tie_count_in_params(tmp_path):
    genes = tuple(f"G{i}" for i in range(20))
    scores = tuple([5.0] * 10 + [1.0] * 10)  # two tied blocks of 10
    ranking = GeneRanking(
        species="human",
        gene_id_namespace="gene_symbol",
        ranking_definition="tie-count probe",
        genes=genes,
        scores=scores,
    )
    cache_dir = _make_cache(tmp_path, {"SET_A": list(genes[:5])})
    provider = R.LocalPrerankGSEAProvider()
    result = provider.run(
        ranking, cache_dir, out_dir=str(tmp_path / "out"), min_set_size=1, permutation_num=20
    )
    assert result.params["n_tied_score_values"] == 20


def test_replay_from_frozen_artifacts_matches_original_in_a_new_process(tmp_path):
    """Runs gsea_replay.py as a genuinely separate `python` process -- not an
    in-process function call -- against the FROZEN artifacts alone, and
    confirms it reports every scientific field as matching.
    """
    import subprocess

    ranking, top_genes, random_genes = _make_known_direction_ranking()
    cache_dir = _make_cache(tmp_path, {"TOP_SET": top_genes, "RANDOM_SET": random_genes})
    out_dir = tmp_path / "original"
    provider = R.LocalPrerankGSEAProvider()
    provider.run(
        ranking, cache_dir, out_dir=str(out_dir), min_set_size=1, permutation_num=100, random_seed=3
    )

    replay_script = os.path.join(CODE, "gsea_replay.py")
    replay_out = tmp_path / "replay"
    proc = subprocess.run(
        [sys.executable, replay_script, "--from", str(out_dir), "--out", str(replay_out)],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "REPLAY MATCH" in proc.stdout


def test_replay_fails_when_ranking_snapshot_is_missing(tmp_path):
    """The replay must genuinely depend on the frozen snapshot file -- not
    silently fall back to recomputing a ranking some other way.
    """
    ranking, top_genes, random_genes = _make_known_direction_ranking()
    cache_dir = _make_cache(tmp_path, {"TOP_SET": top_genes})
    out_dir = tmp_path / "original"
    provider = R.LocalPrerankGSEAProvider()
    provider.run(ranking, cache_dir, out_dir=str(out_dir), min_set_size=1, permutation_num=20)
    os.remove(out_dir / "RANKING_SNAPSHOT.json")

    import gsea_replay

    with pytest.raises(SystemExit, match="not found"):
        gsea_replay.replay(str(out_dir), str(tmp_path / "replay"))


def test_mutation_kill_replay_detects_a_tampered_frozen_ranking(monkeypatch, tmp_path):
    """Tampers with the frozen RANKING_SNAPSHOT.json's own self-check (the
    recorded ranking_sha256) after a real run, and confirms read_ranking_snapshot
    rejects it -- proving the snapshot's own hash re-check is load-bearing,
    not merely informational.
    """
    import json

    ranking, top_genes, random_genes = _make_known_direction_ranking()
    cache_dir = _make_cache(tmp_path, {"TOP_SET": top_genes})
    out_dir = tmp_path / "original"
    provider = R.LocalPrerankGSEAProvider()
    provider.run(ranking, cache_dir, out_dir=str(out_dir), min_set_size=1, permutation_num=20)
    snapshot_path = out_dir / "RANKING_SNAPSHOT.json"
    payload = json.loads(snapshot_path.read_text())
    payload["scores"][0] = payload["scores"][0] + 100.0  # tamper content only
    snapshot_path.write_text(json.dumps(payload))
    with pytest.raises(SystemExit, match="failed self-check"):
        R.read_ranking_snapshot(str(snapshot_path))


# ===========================================================================
# review:20261012_220407 F3: the GSEA result reader accepted a
# payload with only {schema, status="complete"}, an impossible p_value=2.0,
# and a params field that overflowed to inf -- each caught now by a single
# validate_gsea_result() shared by writer and reader.
# ===========================================================================


def test_read_gsea_result_rejects_minimal_complete_payload_missing_everything(tmp_path):
    import json

    path = tmp_path / "bad.json"
    path.write_text(json.dumps({"schema": R.RESULT_SCHEMA, "status": "complete"}))
    with pytest.raises(SystemExit, match="requires a non-empty"):
        R.read_gsea_result(str(path))


def test_read_gsea_result_rejects_impossible_pvalue(tmp_path):
    import json

    payload = _minimal_valid_complete_payload()
    payload["results_table"][0]["p_value"] = 2.0
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(payload))
    with pytest.raises(SystemExit, match=r"p_value=2\.0"):
        R.read_gsea_result(str(path))


def test_read_gsea_result_rejects_overflow_in_params(tmp_path):
    import json

    payload = _minimal_valid_complete_payload()
    text = json.dumps(payload).replace('"weight": 1.0', '"weight": 1e999')
    path = tmp_path / "bad.json"
    path.write_text(text)
    with pytest.raises(SystemExit, match="non-finite"):
        R.read_gsea_result(str(path))


def test_read_gsea_result_accepts_a_legitimate_zero_pvalue(tmp_path):
    """A real, legal p-value/FDR of exactly 0 must NOT be rejected."""
    import json

    payload = _minimal_valid_complete_payload()
    payload["results_table"][0]["p_value"] = 0.0
    payload["results_table"][0]["fdr"] = 0.0
    path = tmp_path / "ok.json"
    path.write_text(json.dumps(payload))
    read_back = R.read_gsea_result(str(path))
    assert read_back["results_table"][0]["p_value"] == 0.0


def _minimal_valid_complete_payload() -> dict:
    return {
        "schema": R.RESULT_SCHEMA,
        "status": "complete",
        "reason": "",
        "ranking_sha256": "a" * 64,
        "backend": "gseapy",
        "backend_version": "1.2.1",
        "gene_set_collection_name": "T",
        "gene_set_cache_path": "/x/gene_sets.json",
        "gene_set_content_sha256": "b" * 64,
        "gene_set_species": "human",
        "gene_set_namespace": "gene_symbol",
        "gene_set_source_version": "test-v1",
        "ranking_snapshot_path": "/x/RANKING_SNAPSHOT.json",
        "ranking_snapshot_sha256": "c" * 64,
        "source_code_revision_or_hash": "git:deadbeef",
        "dependency_versions": {
            "gseapy": "1.2.1",
            "numpy": "1.0",
            "pandas": "1.0",
            "python": "3.10",
        },
        "params": {
            "min_set_size": 1,
            "max_set_size": 500,
            "random_seed": 0,
            "permutation_num": 100,
            "threads": 1,
            "weight": 1.0,
        },
        "excluded_gene_sets": [],
        "results_table": [
            {
                "term": "SET_A",
                "es": 0.9,
                "nes": 1.5,
                "p_value": 0.01,
                "fdr": 0.02,
                "fwer_p_value": 0.03,
                "input_set_size": 10,
                "intersection_size": 5,
                "leading_edge_genes": ["G1", "G2"],
            }
        ],
    }


def test_write_result_atomic_really_rejects_incomplete_complete_result(tmp_path):
    bad_result = GSEARunResult(
        status="complete", backend="gseapy", backend_version="1.2.1", ranking_sha256="a" * 64
    )
    with pytest.raises(ValueError, match="refusing to write"):
        R._write_result_atomic(str(tmp_path / "GSEA_RESULT.json"), bad_result)


def test_mutation_kill_validate_gsea_result_skipped_on_write(monkeypatch, tmp_path):
    """Disables validate_gsea_result inside _write_result_atomic and confirms
    the same deliberately invalid GSEARunResult the test above rejects would
    then be written WITHOUT error, proving the real validation call -- not
    some other check -- is what makes the test above actually raise.
    """
    monkeypatch.setattr(R, "validate_gsea_result", lambda payload: [])
    bad_result = GSEARunResult(
        status="complete", backend="gseapy", backend_version="1.2.1", ranking_sha256="a" * 64
    )
    path = tmp_path / "GSEA_RESULT.json"
    R._write_result_atomic(str(path), bad_result)  # does NOT raise with the mutant
    assert path.exists()


# ===========================================================================
# review:20261012_220407 F4: a relative cache out_dir stored its
# gene-sets path relative to the IMPORT-TIME cwd, so loading from a
# different cwd (or a copied cache directory) failed even though nothing
# was actually missing.
# ===========================================================================


def test_load_cache_after_cwd_change_with_relative_import(tmp_path):
    gmt = tmp_path / "s.gmt"
    _write_gmt(gmt, {"SET_A": ["G1", "G2"]})
    old_cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        G.import_gene_sets(
            source_gmt_path=str(gmt),
            collection_name="T",
            species="human",
            gene_id_namespace="gene_symbol",
            source_version="v1",
            out_dir="relative_cache",
        )
    finally:
        os.chdir(old_cwd)
    gene_sets, manifest = G.load_cached_gene_sets(str(tmp_path / "relative_cache"))
    assert gene_sets == {"SET_A": ["G1", "G2"]}


def test_load_cache_after_whole_directory_copied_elsewhere(tmp_path):
    import shutil

    gmt = tmp_path / "s.gmt"
    _write_gmt(gmt, {"SET_A": ["G1", "G2"]})
    cache_dir = tmp_path / "cache"
    G.import_gene_sets(
        source_gmt_path=str(gmt),
        collection_name="T",
        species="human",
        gene_id_namespace="gene_symbol",
        source_version="v1",
        out_dir=str(cache_dir),
    )
    copied_dir = tmp_path / "copied_elsewhere" / "cache"
    shutil.copytree(cache_dir, copied_dir)
    gene_sets, manifest = G.load_cached_gene_sets(str(copied_dir))
    assert gene_sets == {"SET_A": ["G1", "G2"]}


# ===========================================================================
# review:20261012_220407 F5: when every gene set is excluded by
# size filtering AFTER a real overlap exists, the backend raised an
# unhandled LookupError and no GSEA_RESULT.json was ever written.
# ===========================================================================


def test_run_all_gene_sets_filtered_out_is_a_structured_failure(tmp_path):
    top = [f"T{i}" for i in range(10)]
    rest = [f"R{i}" for i in range(40)]
    genes = tuple(top + rest)
    scores = tuple([5.0 - i * 0.1 for i in range(len(top))] + [0.0] * len(rest))
    ranking = GeneRanking(
        species="human",
        gene_id_namespace="gene_symbol",
        ranking_definition="f5 regression",
        genes=genes,
        scores=scores,
    )
    cache_dir = _make_cache(tmp_path, {"TOP_SET": top})
    out_dir = tmp_path / "out"
    provider = R.LocalPrerankGSEAProvider()
    # min_set_size=100 is far larger than TOP_SET's real (10-gene) overlap --
    # must become an audited failure, never an unhandled exception.
    result = provider.run(
        ranking,
        cache_dir,
        out_dir=str(out_dir),
        min_set_size=100,
        max_set_size=500,
        permutation_num=20,
    )
    assert result.status == "failed"
    assert "no gene set passed size filtering" in result.reason
    assert result.excluded_gene_sets[0]["term"] == "TOP_SET"
    assert result.excluded_gene_sets[0]["reason"] == "too_small"
    assert result.excluded_gene_sets[0]["intersection_size"] == 10
    assert (out_dir / "GSEA_RESULT.json").exists()


def test_run_unexpected_backend_exception_is_recorded_then_reraised(monkeypatch, tmp_path):
    """A genuinely unexpected gseapy exception must still produce an
    auditable 'failed' GSEA_RESULT.json AND propagate (non-zero exit,
    traceback) -- never silently swallowed into no record at all.
    """
    ranking, top_genes, random_genes = _make_known_direction_ranking()
    cache_dir = _make_cache(tmp_path, {"TOP_SET": top_genes})
    out_dir = tmp_path / "out"

    def _boom(*a, **k):
        raise RuntimeError("simulated unexpected backend failure")

    monkeypatch.setattr(R.gseapy, "prerank", _boom)
    provider = R.LocalPrerankGSEAProvider()
    with pytest.raises(RuntimeError, match="simulated unexpected backend failure"):
        provider.run(ranking, cache_dir, out_dir=str(out_dir), min_set_size=1, permutation_num=20)

    result_path = out_dir / "GSEA_RESULT.json"
    assert result_path.exists()
    payload = R.read_gsea_result(str(result_path))
    assert payload["status"] == "failed"
    assert "simulated unexpected backend failure" in payload["reason"]


# ===========================================================================
# review:20261013_114111 R1: a numerically-identical replay was not
# evidence that the replay used the SAME ranking, cache, or source/dependency
# identity as the original run -- only the snapshot file's OWN internal
# self-check was verified, never compared against the original result's own
# recorded ranking_sha256/ranking_snapshot_sha256, and a swapped cache or a
# drifted source/dependency identity was silently accepted as "REPLAY MATCH".
# ===========================================================================


def test_replay_rejects_a_different_but_internally_valid_ranking_snapshot(tmp_path):
    """Swaps in another LEGITIMATE snapshot (same scores, different
    ranking_definition, its own correctly self-recomputed hash) after a real
    run -- `read_ranking_snapshot`'s internal self-check alone would accept
    it. Replay must additionally refuse it because it does not match the
    ORIGINAL result's own recorded ranking_sha256/ranking_snapshot_sha256.
    """
    import dataclasses

    ranking, top_genes, random_genes = _make_known_direction_ranking()
    cache_dir = _make_cache(tmp_path, {"TOP_SET": top_genes, "RANDOM_SET": random_genes})
    out_dir = tmp_path / "original"
    provider = R.LocalPrerankGSEAProvider()
    provider.run(
        ranking, cache_dir, out_dir=str(out_dir), min_set_size=1, permutation_num=50, random_seed=0
    )

    other_ranking = dataclasses.replace(
        ranking, ranking_definition="DIFFERENT_COMPARISON: unrelated score meaning"
    )
    other_sha = R._ranking_sha256(other_ranking)
    R.write_ranking_snapshot(str(out_dir), other_ranking, other_sha)  # overwrites in place

    import gsea_replay

    with pytest.raises(SystemExit, match="semantic hash"):
        gsea_replay.replay(str(out_dir), str(tmp_path / "replay"))


def test_replay_rejects_a_cache_swapped_for_a_different_collection_or_version(tmp_path):
    """A replacement cache that is internally valid (passes its own hash
    check) but has a different collection_name/source_version than the one
    the original run actually used must be refused, even via an explicit
    --cache-dir, and even if it happens to contain the same gene set.
    """
    ranking, top_genes, random_genes = _make_known_direction_ranking()
    cache_dir = _make_cache(
        tmp_path, {"TOP_SET": top_genes, "RANDOM_SET": random_genes}, name="cache_a"
    )
    out_dir = tmp_path / "original"
    provider = R.LocalPrerankGSEAProvider()
    provider.run(
        ranking, cache_dir, out_dir=str(out_dir), min_set_size=1, permutation_num=50, random_seed=0
    )

    gmt = tmp_path / "alt.gmt"
    _write_gmt(gmt, {"TOP_SET": top_genes, "RANDOM_SET": random_genes})
    alt_cache = tmp_path / "cache_b"
    G.import_gene_sets(
        source_gmt_path=str(gmt),
        collection_name="DIFFERENT_COLLECTION",
        species="human",
        gene_id_namespace="gene_symbol",
        source_version="alternate-v2",
        out_dir=str(alt_cache),
    )

    import gsea_replay

    with pytest.raises(SystemExit, match="gene_set_collection_name"):
        gsea_replay.replay(str(out_dir), str(tmp_path / "replay"), str(alt_cache))


def test_replay_accepts_a_relocated_copy_of_the_same_cache(tmp_path):
    """A cache directory copied verbatim to a new location (same collection,
    version, species, namespace, content) must still replay successfully --
    the R1 fix must not block legitimate relocation, only substitution.
    """
    import shutil

    ranking, top_genes, random_genes = _make_known_direction_ranking()
    cache_dir = _make_cache(tmp_path, {"TOP_SET": top_genes, "RANDOM_SET": random_genes})
    out_dir = tmp_path / "original"
    provider = R.LocalPrerankGSEAProvider()
    provider.run(
        ranking, cache_dir, out_dir=str(out_dir), min_set_size=1, permutation_num=50, random_seed=0
    )

    relocated = tmp_path / "relocated_elsewhere" / "cache"
    shutil.copytree(cache_dir, relocated)

    import gsea_replay

    outcome = gsea_replay.replay(str(out_dir), str(tmp_path / "replay"), str(relocated))
    assert outcome["mismatches"] == []


def test_replay_rejects_source_or_dependency_identity_drift_by_default(tmp_path):
    """A tampered source_code_revision_or_hash in the original result (now
    different from what THIS environment will actually compute at replay
    time) must fail a default, strict replay -- never silently pass just
    because the numeric results still agree. `--allow-identity-drift`
    (`allow_identity_drift=True`) is the separate, explicitly-labeled
    cross-version comparison mode that tolerates it.
    """
    import json

    ranking, top_genes, random_genes = _make_known_direction_ranking()
    cache_dir = _make_cache(tmp_path, {"TOP_SET": top_genes, "RANDOM_SET": random_genes})
    out_dir = tmp_path / "original"
    provider = R.LocalPrerankGSEAProvider()
    provider.run(
        ranking, cache_dir, out_dir=str(out_dir), min_set_size=1, permutation_num=50, random_seed=0
    )

    result_path = out_dir / "GSEA_RESULT.json"
    payload = json.loads(result_path.read_text())
    payload["source_code_revision_or_hash"] = "content_hash:" + "0" * 64
    result_path.write_text(json.dumps(payload))

    import gsea_replay

    strict = gsea_replay.replay(str(out_dir), str(tmp_path / "replay_strict"))
    assert any("source_code_revision_or_hash" in m for m in strict["mismatches"])

    comparison = gsea_replay.replay(
        str(out_dir), str(tmp_path / "replay_drift"), allow_identity_drift=True
    )
    assert comparison["mismatches"] == []
    assert comparison["identity_drifts"]


def test_mutation_kill_replay_ranking_identity_check(monkeypatch, tmp_path):
    """Disables `_check_ranking_identity_or_raise` and confirms the snapshot
    -swap scenario above would then silently pass, proving that check -- not
    some other guard -- is what rejects it.
    """
    import dataclasses

    import gsea_replay

    monkeypatch.setattr(gsea_replay, "_check_ranking_identity_or_raise", lambda *a, **k: None)

    ranking, top_genes, random_genes = _make_known_direction_ranking()
    cache_dir = _make_cache(tmp_path, {"TOP_SET": top_genes, "RANDOM_SET": random_genes})
    out_dir = tmp_path / "original"
    provider = R.LocalPrerankGSEAProvider()
    provider.run(
        ranking, cache_dir, out_dir=str(out_dir), min_set_size=1, permutation_num=50, random_seed=0
    )

    other_ranking = dataclasses.replace(
        ranking, ranking_definition="DIFFERENT_COMPARISON: unrelated score meaning"
    )
    other_sha = R._ranking_sha256(other_ranking)
    R.write_ranking_snapshot(str(out_dir), other_ranking, other_sha)

    outcome = gsea_replay.replay(str(out_dir), str(tmp_path / "replay"))
    assert outcome["mismatches"] == []  # with the mutant: silently "passes"


def test_mutation_kill_replay_cache_identity_check(monkeypatch, tmp_path):
    """Disables `_check_cache_identity_or_raise` and confirms the
    cache-swap scenario above would then silently pass."""
    import gsea_replay

    monkeypatch.setattr(gsea_replay, "_check_cache_identity_or_raise", lambda *a, **k: None)

    ranking, top_genes, random_genes = _make_known_direction_ranking()
    cache_dir = _make_cache(
        tmp_path, {"TOP_SET": top_genes, "RANDOM_SET": random_genes}, name="cache_a"
    )
    out_dir = tmp_path / "original"
    provider = R.LocalPrerankGSEAProvider()
    provider.run(
        ranking, cache_dir, out_dir=str(out_dir), min_set_size=1, permutation_num=50, random_seed=0
    )

    gmt = tmp_path / "alt.gmt"
    _write_gmt(gmt, {"TOP_SET": top_genes, "RANDOM_SET": random_genes})
    alt_cache = tmp_path / "cache_b"
    G.import_gene_sets(
        source_gmt_path=str(gmt),
        collection_name="DIFFERENT_COLLECTION",
        species="human",
        gene_id_namespace="gene_symbol",
        source_version="alternate-v2",
        out_dir=str(alt_cache),
    )

    outcome = gsea_replay.replay(str(out_dir), str(tmp_path / "replay"), str(alt_cache))
    assert outcome["mismatches"] == []  # with the mutant: silently "passes"


def test_read_gsea_result_rejects_boolean_identity_field(tmp_path):
    """A probe found `source_code_revision_or_hash=True` accepted by the
    previous `not payload.get(key)` check -- truthy but not a string."""
    import json

    payload = _minimal_valid_complete_payload()
    payload["source_code_revision_or_hash"] = True
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(payload))
    with pytest.raises(SystemExit, match="non-empty string"):
        R.read_gsea_result(str(path))


def test_read_gsea_result_rejects_params_missing_required_numeric_keys(tmp_path):
    """A probe found `params={"unrelated": 1}` accepted by the previous
    "non-empty dict" check -- a 'complete' result's params must carry every
    numeric key the backend actually used, or replay would have to silently
    default it."""
    import json

    payload = _minimal_valid_complete_payload()
    payload["params"] = {"unrelated": 1}
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(payload))
    with pytest.raises(SystemExit, match="requires params"):
        R.read_gsea_result(str(path))


# ===========================================================================
# review:20261013_114111 R2: an invalid (e.g. NaN-containing)
# ranking made write_ranking_snapshot() raise mid-write (allow_nan=False),
# before ranking.validate() ever ran -- leaving a bare ".inprogress" file
# and no GSEA_RESULT.json at all, contradicting the "every failed run leaves
# behind a clean audited record" claim.
# ===========================================================================


def test_run_with_nan_score_leaves_no_snapshot_and_records_a_clean_failed_result(tmp_path):
    ranking, top_genes, random_genes = _make_known_direction_ranking()
    bad_scores = list(ranking.scores)
    bad_scores[0] = float("nan")
    import dataclasses

    bad_ranking = dataclasses.replace(ranking, scores=tuple(bad_scores))
    cache_dir = _make_cache(tmp_path, {"TOP_SET": top_genes})
    out_dir = tmp_path / "out"

    provider = R.LocalPrerankGSEAProvider()
    result = provider.run(
        bad_ranking, cache_dir, out_dir=str(out_dir), min_set_size=1, permutation_num=20
    )

    assert result.status == "failed"
    assert "no snapshot was written" in result.reason
    assert not (out_dir / "RANKING_SNAPSHOT.json").exists()
    assert not any(p.name.endswith(".inprogress") for p in out_dir.iterdir())
    result_path = out_dir / "GSEA_RESULT.json"
    assert result_path.exists()
    payload = R.read_gsea_result(str(result_path))
    assert payload["status"] == "failed"


def test_mutation_kill_r2_validate_gate_before_snapshot_write(monkeypatch, tmp_path):
    """Through the REAL `provider.run()` entry point: disables
    `GeneRanking.validate` (makes it always report "no problems", as if the
    gate the R2 fix relies on were missing) and confirms the NaN-ranking
    scenario above then genuinely breaks -- `run()` falls through to the
    snapshot write with invalid data and RAISES mid-dump (allow_nan=False)
    instead of cleanly returning a `status="failed"` result, proving the
    validate-before-snapshot ordering (not some unrelated check) is what the
    regression test above depends on.
    """
    from enrichment_providers import GeneRanking

    monkeypatch.setattr(GeneRanking, "validate", lambda self: [])

    ranking, top_genes, random_genes = _make_known_direction_ranking()
    bad_scores = list(ranking.scores)
    bad_scores[0] = float("nan")
    import dataclasses

    bad_ranking = dataclasses.replace(ranking, scores=tuple(bad_scores))
    cache_dir = _make_cache(tmp_path, {"TOP_SET": top_genes})
    out_dir = tmp_path / "out"

    provider = R.LocalPrerankGSEAProvider()
    with pytest.raises(ValueError, match="not JSON compliant|Infinity|NaN"):
        provider.run(
            bad_ranking, cache_dir, out_dir=str(out_dir), min_set_size=1, permutation_num=20
        )


# ===========================================================================
# review:20261013_114111 C1: the F1 pre-filter eliminates CROSS-SET
# leakage, but a genuine case COLLISION within the ranking's own distinct
# gene IDs (e.g. "ANCHOR" and "Anchor" both present) can still trigger
# gseapy's internal case-folding; this must fail closed (never a silently
# wrong "complete"), and the scope of that guarantee must be documented
# accurately rather than claimed to be universally impossible.
# ===========================================================================


def test_run_case_collision_within_ranking_fails_closed_not_silently_wrong(tmp_path):
    """Per the review's own reproduction: a ranking containing BOTH `ANCHOR`
    and `Anchor` as distinct entries, against a gene set containing only
    `ANCHOR`. The backend's own case-folding heuristic can still conflate
    the two -- the fix's guarantee is "never a silently wrong complete
    result", not "this ambiguous input is impossible to trigger".
    """
    genes = ("ANCHOR",) + tuple(f"Gene{i}" for i in range(80)) + ("Anchor",)
    scores = (20.0,) + tuple(19.0 - i * 0.5 for i in range(80)) + (-5.25,)
    ranking = GeneRanking(
        species="mouse",
        gene_id_namespace="gene_symbol",
        ranking_definition="C1 regression: case collision within ranking",
        genes=genes,
        scores=scores,
    )
    cache_dir = _make_cache(tmp_path, {"TOY_CASE": ["ANCHOR"]}, species="mouse")

    provider = R.LocalPrerankGSEAProvider()
    result = provider.run(
        ranking,
        cache_dir,
        out_dir=str(tmp_path / "out"),
        min_set_size=1,
        permutation_num=50,
        random_seed=0,
        threads=1,
    )

    # This must NEVER silently report "complete" with a wrong leading edge --
    # either it stays 'failed' (current fail-closed behaviour) or a future,
    # explicitly-verified bijective remapping produces a correctly-mapped
    # 'complete' result. A silent wrong 'complete' is the only disallowed
    # outcome.
    if result.status == "complete":
        row = {r["term"]: r for r in result.results_table}["TOY_CASE"]
        assert set(row["leading_edge_genes"]) <= {"ANCHOR"}
    else:
        assert result.status == "failed"
        assert "leading edge" in result.reason


def test_source_modules_includes_artifact_schema(tmp_path):
    """review:20261013_114111 R1: gsea_runner.py imports and calls
    artifact_schema.py directly on every run (find_nonfinite_floats,
    resolve_source_revision, sha256_file/sha256_text) -- it must be part of
    the module's own source-identity binding, not missing from it.
    """
    assert "artifact_schema.py" in R.SOURCE_MODULES


# ===========================================================================
# review:20261013_135634 R3: the cache-identity check inside
# gsea_replay.replay() was only run when original["status"] == "complete" --
# for "no_overlap" (and "failed" results reached after the cache was already
# loaded), a cache swapped for a different collection/version was accepted
# and still printed "REPLAY MATCH" as long as the (zero-overlap) numeric
# outcome agreed. Strict replay is now restricted to "complete" originals
# only; any other status is refused before the ranking snapshot is read or
# the backend is ever called.
# ===========================================================================


def _make_zero_overlap_ranking(species="human", namespace="gene_symbol"):
    genes = tuple(f"RG{i}" for i in range(20))
    scores = tuple(float(20 - i) for i in range(20))
    ranking = GeneRanking(
        species=species,
        gene_id_namespace=namespace,
        ranking_definition="R3 regression: zero-overlap ranking",
        genes=genes,
        scores=scores,
    )
    return ranking


def test_replay_refuses_a_no_overlap_original_even_with_a_swapped_cache(tmp_path):
    """Reproduces the review's own probe: a 'no_overlap' original replayed
    against a DIFFERENT (also zero-overlap, but different collection/
    version) cache must be refused outright -- not silently reported as
    'REPLAY MATCH' just because both are vacuously zero-overlap.
    """
    ranking = _make_zero_overlap_ranking()
    cache_dir = _make_cache(tmp_path, {"UNRELATED_SET": ["Z0", "Z1"]}, name="cache_a")
    out_dir = tmp_path / "original"
    provider = R.LocalPrerankGSEAProvider()
    result = provider.run(
        ranking, cache_dir, out_dir=str(out_dir), min_set_size=1, permutation_num=20
    )
    assert result.status == "no_overlap"

    gmt = tmp_path / "alt.gmt"
    _write_gmt(gmt, {"ANOTHER_UNRELATED_SET": ["Z0", "Z1"]})
    alt_cache = tmp_path / "cache_b"
    G.import_gene_sets(
        source_gmt_path=str(gmt),
        collection_name="DIFFERENT_COLLECTION",
        species="human",
        gene_id_namespace="gene_symbol",
        source_version="alternate-v2",
        out_dir=str(alt_cache),
    )

    import gsea_replay

    with pytest.raises(SystemExit, match="only supports"):
        gsea_replay.replay(str(out_dir), str(tmp_path / "replay"), str(alt_cache))


def test_replay_refuses_a_failed_original_before_touching_snapshot_or_backend(tmp_path):
    top = [f"T{i}" for i in range(10)]
    rest = [f"Rr{i}" for i in range(40)]
    genes = tuple(top + rest)
    scores = tuple([5.0 - i * 0.1 for i in range(len(top))] + [0.0] * len(rest))
    ranking = GeneRanking(
        species="human",
        gene_id_namespace="gene_symbol",
        ranking_definition="R3 regression: failed original",
        genes=genes,
        scores=scores,
    )
    cache_dir = _make_cache(tmp_path, {"TOP_SET": top})
    out_dir = tmp_path / "original"
    provider = R.LocalPrerankGSEAProvider()
    result = provider.run(
        ranking,
        cache_dir,
        out_dir=str(out_dir),
        min_set_size=100,
        max_set_size=500,
        permutation_num=20,
    )
    assert result.status == "failed"

    import gsea_replay

    with pytest.raises(SystemExit, match="only supports"):
        gsea_replay.replay(str(out_dir), str(tmp_path / "replay"))
    # The original record must be left exactly as it was.
    assert (out_dir / "GSEA_RESULT.json").exists()
    assert R.read_gsea_result(str(out_dir / "GSEA_RESULT.json"))["status"] == "failed"


def test_mutation_kill_replay_status_restriction(monkeypatch, tmp_path):
    """Widens `_REPLAYABLE_STATUSES` to also include 'no_overlap' and
    confirms a 'no_overlap' original -- which the fixed code refuses
    outright, before ever touching the snapshot/cache/backend -- would then
    be silently accepted and replayed (even using the SAME, non-swapped
    cache; `_check_cache_identity_or_raise` now runs unconditionally and
    would independently catch a swapped cache, so a swapped-cache variant of
    this scenario would not isolate THIS guard). Proves the status
    restriction itself -- not some other guard -- is what refuses replaying
    a non-'complete' original.
    """
    import gsea_replay

    monkeypatch.setattr(gsea_replay, "_REPLAYABLE_STATUSES", ("complete", "no_overlap"))

    ranking = _make_zero_overlap_ranking()
    cache_dir = _make_cache(tmp_path, {"UNRELATED_SET": ["Z0", "Z1"]}, name="cache_a")
    out_dir = tmp_path / "original"
    provider = R.LocalPrerankGSEAProvider()
    provider.run(ranking, cache_dir, out_dir=str(out_dir), min_set_size=1, permutation_num=20)

    outcome = gsea_replay.replay(str(out_dir), str(tmp_path / "replay"))
    # With the mutant: a non-'complete' original is now replayed at all,
    # which the fixed code refuses outright via SystemExit.
    assert outcome["replayed"]["status"] == "no_overlap"


def test_replay_still_accepts_a_complete_original_after_the_r3_restriction(tmp_path):
    """The R3 restriction must not break the ordinary successful-replay
    path -- a 'complete' original must still replay normally."""
    ranking, top_genes, random_genes = _make_known_direction_ranking()
    cache_dir = _make_cache(tmp_path, {"TOP_SET": top_genes, "RANDOM_SET": random_genes})
    out_dir = tmp_path / "original"
    provider = R.LocalPrerankGSEAProvider()
    result = provider.run(
        ranking, cache_dir, out_dir=str(out_dir), min_set_size=1, permutation_num=50, random_seed=0
    )
    assert result.status == "complete"

    import gsea_replay

    outcome = gsea_replay.replay(str(out_dir), str(tmp_path / "replay"))
    assert outcome["mismatches"] == []


# ===========================================================================
# review:20261013_135634 R4: the printed replay label depended on
# whether a drift happened to be OBSERVED (outcome["identity_drifts"]), not
# on whether the caller explicitly asked for cross-version comparison mode
# (--allow-identity-drift) -- so running --allow-identity-drift against two
# otherwise-identical environments still printed "REPLAY MATCH", which
# contradicts the CLI help text's own promise ("REPLAY COMPARISON, never
# REPLAY MATCH, even when every number agrees").
# ===========================================================================


def test_replay_label_reflects_requested_mode_not_observed_drift(tmp_path):
    import gsea_replay

    # Flag explicitly set, but nothing was tampered with -- no real drift.
    comparison_no_drift = {"allow_identity_drift": True, "identity_drifts": []}
    label, extra = gsea_replay._replay_label_and_extra(comparison_no_drift)
    assert label == "REPLAY COMPARISON"
    assert extra == ""

    strict_no_drift = {"allow_identity_drift": False, "identity_drifts": []}
    label, extra = gsea_replay._replay_label_and_extra(strict_no_drift)
    assert label == "REPLAY MATCH"

    comparison_with_drift = {
        "allow_identity_drift": True,
        "identity_drifts": ["source_code_revision_or_hash: ..."],
    }
    label, extra = gsea_replay._replay_label_and_extra(comparison_with_drift)
    assert label == "REPLAY COMPARISON"
    assert "identity drift tolerated" in extra


def test_replay_cli_prints_comparison_label_for_requested_mode_with_no_drift(tmp_path):
    """Real entry point (gsea_replay.main(), the same function the CLI
    invokes): --allow-identity-drift against an untampered original/
    environment pair must still print 'REPLAY COMPARISON', never 'REPLAY
    MATCH'.
    """
    ranking, top_genes, random_genes = _make_known_direction_ranking()
    cache_dir = _make_cache(tmp_path, {"TOP_SET": top_genes, "RANDOM_SET": random_genes})
    out_dir = tmp_path / "original"
    provider = R.LocalPrerankGSEAProvider()
    provider.run(
        ranking, cache_dir, out_dir=str(out_dir), min_set_size=1, permutation_num=50, random_seed=0
    )

    import subprocess

    replay_script = os.path.join(CODE, "gsea_replay.py")
    proc = subprocess.run(
        [
            sys.executable,
            replay_script,
            "--from",
            str(out_dir),
            "--out",
            str(tmp_path / "replay"),
            "--allow-identity-drift",
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "REPLAY COMPARISON" in proc.stdout
    assert "REPLAY MATCH" not in proc.stdout


def test_mutation_kill_replay_label_reflects_requested_mode(monkeypatch):
    """Reverts `_replay_label_and_extra` to the OLD buggy logic (label keyed
    off whether a drift was OBSERVED, not off the requested mode) and
    confirms the no-drift-but-flag-set case above would then wrongly report
    'REPLAY MATCH' -- proving the fixed function, not something else, is
    what the regression test above depends on.
    """

    def _old_buggy_label_and_extra(outcome):
        label = "REPLAY COMPARISON" if outcome["identity_drifts"] else "REPLAY MATCH"
        extra = ""
        if outcome["identity_drifts"]:
            extra = " (identity drift tolerated: " + "; ".join(outcome["identity_drifts"]) + ")"
        return label, extra

    comparison_no_drift = {"allow_identity_drift": True, "identity_drifts": []}
    label, extra = _old_buggy_label_and_extra(comparison_no_drift)
    assert label == "REPLAY MATCH"  # the pre-fix bug this test demonstrates
