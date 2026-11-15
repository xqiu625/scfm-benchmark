"""Enrichment provider contract tests. Confirms GSEA (method) and GO/KEGG
(annotation source) stay separate classes, and that neither makes a live call
this round -- every entry point must raise NotImplementedError, not silently
return an empty/placeholder result that could be mistaken for 'ran, found
nothing'.
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

import enrichment_providers as E  # noqa: E402


def test_gsea_and_annotation_are_separate_classes():
    assert E.GSEAProvider is not E.AnnotationProvider
    assert not issubclass(E.GSEAProvider, E.AnnotationProvider)
    assert not issubclass(E.AnnotationProvider, E.GSEAProvider)


def test_gene_ranking_validate_catches_length_mismatch():
    ranking = E.GeneRanking(
        species="human",
        gene_id_namespace="ensembl_gene_id",
        ranking_definition="signed rho",
        genes=("A", "B", "C"),
        scores=(1.0, 2.0),
    )
    problems = ranking.validate()
    assert any("length mismatch" in p for p in problems)


def test_gene_ranking_validate_catches_duplicate_genes():
    ranking = E.GeneRanking(
        species="human",
        gene_id_namespace="ensembl_gene_id",
        ranking_definition="signed rho",
        genes=("A", "B", "A"),
        scores=(1.0, 2.0, 3.0),
    )
    problems = ranking.validate()
    assert any("duplicate gene IDs" in p for p in problems)


def test_gene_ranking_validate_catches_empty_definition():
    ranking = E.GeneRanking(
        species="human",
        gene_id_namespace="ensembl_gene_id",
        ranking_definition="",
        genes=("A",),
        scores=(1.0,),
    )
    problems = ranking.validate()
    assert any("ranking_definition is empty" in p for p in problems)


def test_valid_gene_ranking_has_no_problems():
    ranking = E.GeneRanking(
        species="human",
        gene_id_namespace="ensembl_gene_id",
        ranking_definition="signed rho",
        genes=("A", "B"),
        scores=(1.0, -2.0),
    )
    assert ranking.validate() == []


def test_gsea_provider_run_is_not_implemented_and_makes_no_network_call():
    """The base GSEAProvider stays a contract stub -- the real implementation
    is gsea_runner.LocalPrerankGSEAProvider (see test_scfm_gsea_runner.py).
    """
    provider = E.GSEAProvider()
    ranking = E.GeneRanking(
        species="human",
        gene_id_namespace="ensembl_gene_id",
        ranking_definition="signed rho",
        genes=("A", "B"),
        scores=(1.0, -2.0),
    )
    with pytest.raises(NotImplementedError, match="no gseapy call, no network access"):
        provider.run(ranking, "/tmp/does-not-matter", out_dir="/tmp/also-does-not-matter")


def test_gene_set_cache_entry_records_local_import_not_a_fabricated_retrieval():
    """review:20261012_171916 next-increment instruction: an offline
    import must record its OWN local import time under an explicit
    source_type, never a bare 'retrieved_utc' that implies a live fetch.
    """
    entry = E.GeneSetCacheEntry(
        collection_name="TEST_COLLECTION",
        species="human",
        gene_id_namespace="gene_symbol",
        source_release_or_version="v1-local-test",
        source_type="local_file",
        imported_or_retrieved_utc="2026-01-01T00:00:00Z",
        content_sha256="x" * 64,
        cache_path="/tmp/does-not-matter/gene_sets.json",
    )
    assert entry.source_type == "local_file"


def test_gene_ranking_requires_species_explicitly():
    with pytest.raises(TypeError):
        E.GeneRanking(
            gene_id_namespace="ensembl_gene_id",  # species omitted
            ranking_definition="signed rho",
            genes=("A",),
            scores=(1.0,),
        )


def test_gene_ranking_validate_catches_empty_species():
    ranking = E.GeneRanking(
        species="",
        gene_id_namespace="ensembl_gene_id",
        ranking_definition="signed rho",
        genes=("A",),
        scores=(1.0,),
    )
    problems = ranking.validate()
    assert any("species is empty" in p for p in problems)


def test_gene_ranking_validate_catches_empty_gene_id():
    ranking = E.GeneRanking(
        species="human",
        gene_id_namespace="ensembl_gene_id",
        ranking_definition="signed rho",
        genes=("A", ""),
        scores=(1.0, 2.0),
    )
    problems = ranking.validate()
    assert any("empty/non-string" in p for p in problems)


def test_annotation_provider_fetch_is_not_implemented_and_makes_no_network_call():
    provider = E.AnnotationProvider()
    with pytest.raises(NotImplementedError, match="no QuickGO/KEGG/MSigDB request"):
        provider.fetch("NCBITaxon:9606", cache_dir="/tmp")


def test_kegg_rate_limit_is_recorded_not_assumed():
    assert E.AnnotationProvider.KEGG_RATE_LIMIT_PER_SECOND == 3


# --- R4 (review:20261012_153316): the validator must not silently
# accept NaN scores or a fully-empty ranking (both previously returned []). ---


def test_gene_ranking_validate_catches_nan_score():
    ranking = E.GeneRanking(
        species="human",
        gene_id_namespace="ensembl_gene_id",
        ranking_definition="signed rho",
        genes=("A", "B"),
        scores=(1.0, float("nan")),
    )
    problems = ranking.validate()
    assert any("not finite" in p for p in problems)


def test_gene_ranking_validate_catches_inf_score():
    ranking = E.GeneRanking(
        species="human",
        gene_id_namespace="ensembl_gene_id",
        ranking_definition="signed rho",
        genes=("A", "B"),
        scores=(1.0, float("inf")),
    )
    problems = ranking.validate()
    assert any("not finite" in p for p in problems)


def test_gene_ranking_validate_catches_fully_empty_ranking():
    """The exact case review:20261012_153316 R4 reproduced as
    returning [] -- a non-empty ranking_definition with empty genes/scores
    and an empty namespace.
    """
    ranking = E.GeneRanking(
        species="human",
        gene_id_namespace="",
        ranking_definition="some ranking",
        genes=(),
        scores=(),
    )
    problems = ranking.validate()
    assert any("genes is empty" in p for p in problems)
    assert any("gene_id_namespace is empty" in p for p in problems)


def test_gene_ranking_validate_accepts_one_sided_ranking():
    """A ranking entirely on one side (all positive, say) is NOT itself a
    defect -- must not be flagged.
    """
    ranking = E.GeneRanking(
        species="human",
        gene_id_namespace="ensembl_gene_id",
        ranking_definition="signed rho",
        genes=("A", "B", "C"),
        scores=(1.0, 2.0, 3.0),
    )
    assert ranking.validate() == []


def test_mutation_kill_nan_check(monkeypatch):
    """Disables the finite-score check and confirms the NaN test above would
    then genuinely fail to detect the fault.
    """
    import math as math_mod

    real_isfinite = math_mod.isfinite
    monkeypatch.setattr(math_mod, "isfinite", lambda x: True)

    ranking = E.GeneRanking(
        species="human",
        gene_id_namespace="ensembl_gene_id",
        ranking_definition="signed rho",
        genes=("A", "B"),
        scores=(1.0, float("nan")),
    )
    problems = ranking.validate()
    assert not any(
        "not finite" in p for p in problems
    ), "mutant unexpectedly still caught the NaN -- mutation did not take effect"
    monkeypatch.setattr(math_mod, "isfinite", real_isfinite)
    # Restored: the real check must catch it again.
    assert any("not finite" in p for p in ranking.validate())
