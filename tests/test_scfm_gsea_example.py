"""Runs the full offline GSEA example end to end, so it stays exercised by the
regular test suite. Per `review:20261012_171916` item 6: a local gene
set with a known direction, a complete ranking, a real backend call, persisted
results, and a manifest re-read -- not a mocked or synthetic-only unit test.
"""

import importlib
import json
import os
import subprocess
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLES_DIR = os.path.join(REPO_ROOT, "scripts", "34_scfm_discovery_framework", "examples")
CODE_DIR = os.path.join(REPO_ROOT, "scripts", "34_scfm_discovery_framework", "code")
EXAMPLE = os.path.join(EXAMPLES_DIR, "run_gsea_example.py")


def test_run_gsea_example_end_to_end(tmp_path):
    out_dir = tmp_path / "gsea_example_run"
    result = subprocess.run(
        [sys.executable, EXAMPLE, "--out-dir", str(out_dir)], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "OK: read_gsea_result" in result.stdout

    assert (out_dir / "toy_gene_sets.gmt").exists()
    assert (out_dir / "cache" / "CACHE_MANIFEST.json").exists()
    result_path = out_dir / "gsea" / "GSEA_RESULT.json"
    assert result_path.exists()
    payload = json.loads(result_path.read_text())
    assert payload["status"] == "complete"
    by_term = {row["term"]: row for row in payload["results_table"]}
    assert by_term["TOY_ENRICHED_SET"]["es"] > by_term["TOY_UNRELATED_SET"]["es"]


def test_run_gsea_example_refuses_nonempty_out_dir(tmp_path):
    out_dir = tmp_path / "gsea_example_run"
    out_dir.mkdir()
    (out_dir / "junk.txt").write_text("pre-existing")
    result = subprocess.run(
        [sys.executable, EXAMPLE, "--out-dir", str(out_dir)], capture_output=True, text=True
    )
    assert result.returncode != 0
    assert (out_dir / "junk.txt").read_text() == "pre-existing"


def _run_main_with_argv(module, out_dir):
    old_argv = sys.argv
    sys.argv = ["run_gsea_example.py", "--out-dir", out_dir]
    try:
        module.main()
    finally:
        sys.argv = old_argv


def test_run_gsea_example_detects_a_real_sign_flip_regression(tmp_path, monkeypatch):
    """Injects a sign flip into the ranking scores right before the real
    gseapy call and confirms the example's OWN internal consistency check
    (step 5) catches it -- proving the example actually verifies its
    known-direction claim, not merely "runs without crashing".
    """
    sys.path.insert(0, CODE_DIR)
    sys.path.insert(0, EXAMPLES_DIR)
    run_gsea_example = importlib.import_module("run_gsea_example")

    real_build_ranking = run_gsea_example.build_known_direction_ranking

    def _flipped(enriched_genes, unrelated_genes):
        ranking = real_build_ranking(enriched_genes, unrelated_genes)
        from enrichment_providers import GeneRanking

        # Flip every score's sign -- TOY_ENRICHED_SET's genes are now at the
        # BOTTOM of the ranking, not the top.
        return GeneRanking(
            species=ranking.species,
            gene_id_namespace=ranking.gene_id_namespace,
            ranking_definition=ranking.ranking_definition,
            genes=ranking.genes,
            scores=tuple(-s for s in ranking.scores),
        )

    monkeypatch.setattr(run_gsea_example, "build_known_direction_ranking", _flipped)
    out_dir = tmp_path / "gsea_example_flipped"
    with pytest.raises(SystemExit, match="MISMATCH"):
        _run_main_with_argv(run_gsea_example, str(out_dir))
