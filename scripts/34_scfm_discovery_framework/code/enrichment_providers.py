"""Enrichment/annotation provider contracts for the scFM discovery framework.

Per `review:20261012_145004` item C.4 and `review:20261012_133205`
§3: GSEA is an analysis METHOD; GO and KEGG are annotation/pathway SOURCES. They are
kept as two different kinds of object on purpose, never one generic "enrichment API"
class. `AnnotationProvider` (GO/KEGG/MSigDB live fetch) is still a contract-only stub
this round -- every call raises `NotImplementedError` naming exactly what is missing,
so a caller cannot mistake "not implemented yet" for "ran and found nothing" (DEVELOPMENT.md:
"API 失败必须标为失败/部分完成，不能写成'无显著通路'"). `GSEAProvider` now has a REAL
local implementation (`gsea_runner.LocalPrerankGSEAProvider`, per
`review:20261012_171916`'s next-increment instruction) -- see that module.
This file stays the abstract contract layer.

This corrects, rather than reuses, `scripts/24_option_c_canonical_input/code/
gsea_per_fm_dim.py`'s pattern of calling Enrichr live by library name with no cached
gene-set file, no content hash, and no offline replay (noted as a real gap in
MODEL_SELECTION_v1.md's "reusable components" section) -- the ranking/prerank
mechanic there is fine; the live, unversioned fetch is not reused here.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field


@dataclass(frozen=True)
class GeneRanking:
    """A complete, signed, unique-gene-ID ranking -- the required GSEA input.

    Per `review:20261012_133205` §3: GSEA needs the full ranking, not a
    pre-filtered top-N list; duplicate IDs and the ranking's direction/definition
    must be recorded, not implied by column order. A ranking is NOT required to
    contain both positive and negative scores -- a real, valid ranking may land
    entirely on one side; that is not itself a defect (per
    `review:20261012_153316` R4's own instruction not to add that
    requirement). `species` is required (per
    `review:20261012_171916`'s next-increment instruction to check
    species consistency between a ranking and a gene set) -- never defaulted,
    since a silently-assumed species is exactly the kind of mistake this
    framework's other contracts (e.g. `model_registry.py`'s pretraining-vs-
    acceptance split) are built to avoid repeating.
    """
    species: str
    gene_id_namespace: str
    ranking_definition: str          # e.g. "signed Spearman rho vs. state score S"
    genes: tuple                     # gene IDs, same order as `scores`
    scores: tuple                    # same length as `genes`

    def validate(self) -> list[str]:
        problems = []
        if not self.species:
            problems.append("species is empty -- state it explicitly, never assumed")
        if not self.ranking_definition:
            problems.append("ranking_definition is empty -- state what the score IS")
        if not self.gene_id_namespace:
            problems.append("gene_id_namespace is empty -- state what kind of gene "
                             "ID these are")
        if not self.genes:
            problems.append("genes is empty -- no rankable input; report this as "
                             "not-runnable, never as a completed run with no "
                             "significant pathways")
        if len(self.genes) != len(self.scores):
            problems.append(
                f"genes ({len(self.genes)}) and scores ({len(self.scores)}) length mismatch")
        if len(set(self.genes)) != len(self.genes):
            dupes = [g for g in self.genes if self.genes.count(g) > 1]
            problems.append(f"duplicate gene IDs present: {sorted(set(dupes))}")
        empty_ids = [g for g in self.genes if not isinstance(g, str) or not g.strip()]
        if empty_ids:
            problems.append(f"{len(empty_ids)} gene ID(s) are empty/non-string")
        non_finite = [g for g, s in zip(self.genes, self.scores)
                      if not isinstance(s, (int, float)) or isinstance(s, bool)
                      or not math.isfinite(s)]
        if non_finite:
            problems.append(
                f"{len(non_finite)} score(s) are not finite numbers (NaN/Inf/"
                f"wrong type): genes {non_finite[:5]}"
                + (" (truncated)" if len(non_finite) > 5 else ""))
        return problems


@dataclass(frozen=True)
class GeneSetCacheEntry:
    """A single cached, versioned gene-set file -- what a GSEA run is actually
    allowed to read. Never a live library name resolved at run time.

    `source_type`/`imported_or_retrieved_utc` replace a bare `retrieved_utc`
    (per `review:20261012_171916`'s next-increment instruction): an
    offline import from a local file must record the LOCAL import time under
    `source_type="local_file"`, never fabricate a live retrieval timestamp
    that implies a network fetch happened.
    """
    collection_name: str              # e.g. "KEGG_2021_Human", "GO:BP"
    species: str
    gene_id_namespace: str
    source_release_or_version: str
    source_type: str                  # "local_file" -- the only type implemented this round
    imported_or_retrieved_utc: str
    content_sha256: str
    cache_path: str


@dataclass(frozen=True)
class GSEARunResult:
    #: "complete": ran, produced a results table (possibly with no
    #: significant terms -- that is a real, distinct outcome, not a failure).
    #: "failed": a structural problem (invalid ranking, species/namespace
    #: mismatch, cache hash mismatch, backend error) stopped the run before
    #: or during execution. "no_overlap": the ranking shares zero genes with
    #: every gene set in the collection -- distinct from "complete with no
    #: significant terms". "not_run": never attempted.
    status: str
    reason: str = ""                  # required when status != "complete"
    ranking_sha256: str = ""
    gene_set_collection_name: str = ""
    gene_set_cache_path: str = ""
    gene_set_content_sha256: str = ""
    backend: str = "unknown"
    backend_version: str = "unknown"
    #: seed, threads, permutation_num, min/max set size, weight, etc. -- every
    #: parameter that could change the numeric result, recorded explicitly.
    params: dict = field(default_factory=dict)
    #: Each row: term, es, nes, p_value, fdr, fwer_p_value, input_set_size,
    #: intersection_size, leading_edge_genes. Empty when status != "complete".
    results_table: tuple = ()
    #: Per review:20261012_220407 F2: the frozen full ranking's own
    #: path/hash (the exact gene/score order actually sent to the backend),
    #: the wrapper's own source-code revision/hash, and the numeric-result-
    #: relevant dependency versions (gseapy/numpy/pandas/python) -- all
    #: bound into THIS result, not left implicit in the surrounding code.
    ranking_snapshot_path: str = ""
    ranking_snapshot_sha256: str = ""
    source_code_revision_or_hash: str = ""
    dependency_versions: dict = field(default_factory=dict)
    #: Per `review:20261013_114111` R1: the cache's own semantic
    #: identity (not just its collection name/content hash) -- so a strict
    #: replay can refuse a cache that was swapped for a different
    #: species/namespace/version collection even if its content happens to
    #: produce byte-identical numeric results for the one term compared.
    gene_set_species: str = ""
    gene_set_namespace: str = ""
    gene_set_source_version: str = ""
    #: Per F5: every gene set that was excluded before the backend was ever
    #: called, with its own reason (too_small/too_large relative to the
    #: EXACT, case-sensitive ranking/gene-set intersection) -- kept distinct
    #: from "complete with no significant terms" and from a silently empty
    #: results_table.
    excluded_gene_sets: tuple = ()

    def to_dict(self) -> dict:
        d = dict(self.__dict__)
        d["results_table"] = [dict(row) for row in self.results_table]
        d["excluded_gene_sets"] = [dict(row) for row in self.excluded_gene_sets]
        return d


class GSEAProvider:
    """Method, not a data source. Must run against a `GeneRanking` and a cached
    gene-set collection ONLY -- never a live library name. This base class
    stays a contract stub; see `gsea_runner.LocalPrerankGSEAProvider` for the
    real, gseapy-backed implementation.
    """

    def run(self, ranking: GeneRanking, gene_set_cache_dir: str,
            *, out_dir: str, min_set_size: int = 10, max_set_size: int = 500,
            random_seed: int = 0, permutation_num: int = 100,
            threads: int = 1, weight: float = 1.0) -> GSEARunResult:
        raise NotImplementedError(
            "GSEAProvider.run is a contract stub -- no gseapy call, no "
            "network access, no gene-set fetch is made here. Use "
            "gsea_runner.LocalPrerankGSEAProvider for the real local "
            "implementation.")


@dataclass(frozen=True)
class AnnotationFetchResult:
    status: str                       # "complete" | "failed" | "rate_limited" | "not_run"
    species_taxonomy_id: str
    retrieved_utc: str
    content_sha256: str
    cache_path: str
    reason: str = ""                  # required when status != "complete"


class AnnotationProvider:
    """GO/KEGG/MSigDB annotation SOURCE -- fetch-and-cache only, never computes
    enrichment itself (that is GSEAProvider's/an ORA implementation's job, kept
    separate per `review:20261012_133205` §3).
    """
    #: KEGG's own terms of use cap academic REST use at 3 requests/second --
    #: recorded here so any future real implementation has to look at it, not
    #: discover it by getting rate-limited in production.
    KEGG_RATE_LIMIT_PER_SECOND = 3

    def fetch(self, species_taxonomy_id: str, *, cache_dir: str) -> AnnotationFetchResult:
        raise NotImplementedError(
            "AnnotationProvider.fetch is a contract stub this round -- no QuickGO/"
            "KEGG/MSigDB request is made here. A real implementation must cache the "
            "response with retrieval date + content hash + species, rate-limit "
            "itself, and report rate_limited/failed explicitly rather than an empty "
            "'no pathways found' result.")
