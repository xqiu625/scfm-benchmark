"""Runnable, fully offline, end-to-end example of the local gene-set cache +
GSEA chain, per `review:20261012_171916`'s next-increment instruction
(item 6): a local gene set with a KNOWN enrichment direction, a complete gene
ranking, a REAL backend call (no mocking), persisted results, and a manifest
re-read.

    python scripts/34_scfm_discovery_framework/examples/run_gsea_example.py \
        --out-dir <new/empty/directory>

Builds a tiny synthetic `.gmt` with two gene sets -- `TOY_ENRICHED_SET`
(deliberately given the highest scores in the ranking, BY CONSTRUCTION, not by
chance) and `TOY_UNRELATED_SET` (an arbitrary, unrelated set) -- named `TOY_*`
specifically so nobody mistakes this for real GO/KEGG/MSigDB pathway data.
Imports it into a local versioned cache (recording a LOCAL import time, never
a fabricated live "retrieval" timestamp), builds the full ranking, runs the
real gseapy-backed `LocalPrerankGSEAProvider`, persists `GSEA_RESULT.json`,
and reads it back -- asserting `TOY_ENRICHED_SET`'s ES/NES exceed
`TOY_UNRELATED_SET`'s. Never makes a network call, never reads real
expression data. Everything this script writes goes under `--out-dir`; it
refuses a non-empty one, same convention as the rest of this framework.
"""
from __future__ import annotations

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.join(os.path.dirname(HERE), "code")
sys.path.insert(0, CODE_DIR)


def build_synthetic_gmt(path: str) -> tuple[list[str], list[str]]:
    """Writes a tiny synthetic .gmt. Returns (enriched_genes, unrelated_genes)."""
    enriched_genes = [f"TOYGENE{i}" for i in range(15)]
    unrelated_genes = [f"TOYGENE{i}" for i in range(80, 90)]
    lines = [
        "\t".join(["TOY_ENRICHED_SET", "synthetic, boosted by construction", *enriched_genes]),
        "\t".join(["TOY_UNRELATED_SET", "synthetic, arbitrary", *unrelated_genes]),
    ]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return enriched_genes, unrelated_genes


def build_known_direction_ranking(enriched_genes: list[str], unrelated_genes: list[str]):
    import numpy as np
    from enrichment_providers import GeneRanking

    rng = np.random.default_rng(0)
    n_genes = 150
    all_genes = [f"TOYGENE{i}" for i in range(n_genes)]
    scores = {g: float(rng.normal()) for g in all_genes}
    # Deliberately place TOY_ENRICHED_SET's genes at the very top of the
    # ranking -- the enrichment direction is known by construction, not
    # hoped for from a random draw.
    for i, g in enumerate(enriched_genes):
        scores[g] = 5.0 - i * 0.05
    genes = tuple(all_genes)
    score_values = tuple(scores[g] for g in genes)
    return GeneRanking(
        species="human", gene_id_namespace="gene_symbol",
        ranking_definition="synthetic ranking; TOY_ENRICHED_SET boosted by "
                           "construction to validate the ES direction, not a "
                           "real biological measurement",
        genes=genes, scores=score_values)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    if os.path.exists(args.out_dir) and os.listdir(args.out_dir):
        raise SystemExit(f"{args.out_dir} must be new or empty")
    os.makedirs(args.out_dir, exist_ok=True)

    import gene_set_cache as G
    import gsea_runner as R

    print("## Step 1/6 -- build a tiny synthetic .gmt (clearly named TOY_*, "
          "never real GO/KEGG/MSigDB data)")
    gmt_path = os.path.join(args.out_dir, "toy_gene_sets.gmt")
    enriched_genes, unrelated_genes = build_synthetic_gmt(gmt_path)
    print(f"wrote: {gmt_path}")

    print("\n## Step 2/6 -- import into a local, versioned, hash-verified cache "
          "(no network, no live retrieval timestamp fabricated)")
    cache_dir = os.path.join(args.out_dir, "cache")
    manifest = G.import_gene_sets(
        source_gmt_path=gmt_path, collection_name="TOY_COLLECTION", species="human",
        gene_id_namespace="gene_symbol", source_version="example-v1",
        out_dir=cache_dir)
    print(f"cache manifest: imported_utc={manifest['imported_utc']!r} "
          f"source_type={manifest['source_type']!r} "
          f"source_content_sha256={manifest['source_content_sha256'][:16]}...")

    print("\n## Step 3/6 -- build the full gene ranking (TOY_ENRICHED_SET boosted "
          "by construction, not by chance)")
    ranking = build_known_direction_ranking(enriched_genes, unrelated_genes)
    print(f"ranking: {len(ranking.genes)} genes, species={ranking.species!r}, "
          f"namespace={ranking.gene_id_namespace!r}")

    print("\n## Step 4/6 -- run the REAL gseapy-backed LocalPrerankGSEAProvider "
          "(no mocking, no network)")
    gsea_out = os.path.join(args.out_dir, "gsea")
    provider = R.LocalPrerankGSEAProvider()
    result = provider.run(
        ranking, cache_dir, out_dir=gsea_out,
        min_set_size=1, max_set_size=500, permutation_num=500, threads=1,
        random_seed=0)
    if result.status != "complete":
        raise SystemExit(f"GSEA run did not complete: status={result.status} "
                         f"reason={result.reason}")
    by_term = {row["term"]: row for row in result.results_table}
    print(f"TOY_ENRICHED_SET:  ES={by_term['TOY_ENRICHED_SET']['es']:.3f}  "
          f"NES={by_term['TOY_ENRICHED_SET']['nes']:.3f}  "
          f"p={by_term['TOY_ENRICHED_SET']['p_value']:.4f}  "
          f"FDR={by_term['TOY_ENRICHED_SET']['fdr']:.4f}")
    print(f"TOY_UNRELATED_SET: ES={by_term['TOY_UNRELATED_SET']['es']:.3f}  "
          f"NES={by_term['TOY_UNRELATED_SET']['nes']:.3f}  "
          f"p={by_term['TOY_UNRELATED_SET']['p_value']:.4f}  "
          f"FDR={by_term['TOY_UNRELATED_SET']['fdr']:.4f}")

    print("\n## Step 5/6 -- read the persisted result back and verify the "
          "known enrichment direction")
    read_back = R.read_gsea_result(os.path.join(gsea_out, "GSEA_RESULT.json"))
    read_by_term = {row["term"]: row for row in read_back["results_table"]}
    enriched_row = read_by_term["TOY_ENRICHED_SET"]
    unrelated_row = read_by_term["TOY_UNRELATED_SET"]
    if not (enriched_row["es"] > unrelated_row["es"]):
        raise SystemExit(
            "MISMATCH: TOY_ENRICHED_SET's ES is not greater than TOY_UNRELATED_SET's "
            "-- the known-direction fixture did not reproduce; this would be a real bug")
    if not (enriched_row["nes"] > unrelated_row["nes"]):
        raise SystemExit(
            "MISMATCH: TOY_ENRICHED_SET's NES is not greater than TOY_UNRELATED_SET's "
            "-- the known-direction fixture did not reproduce; this would be a real bug")
    if set(enriched_genes) - set(enriched_row["leading_edge_genes"]):
        raise SystemExit(
            "MISMATCH: not all of TOY_ENRICHED_SET's own genes appear in its own "
            "leading edge, despite being placed at the very top of the ranking")
    print("OK: read_gsea_result's object confirms TOY_ENRICHED_SET's ES AND NES both "
          "exceed TOY_UNRELATED_SET's, and its genes are all in its own leading edge.")

    print("\n## Step 6/6 -- replay the run from its FROZEN artifacts alone, in a new "
          "process, and confirm every scientific field matches exactly")
    import subprocess
    replay_out = os.path.join(args.out_dir, "gsea_replay")
    replay_script = os.path.join(CODE_DIR, "gsea_replay.py")
    replay_result = subprocess.run(
        [sys.executable, replay_script, "--from", gsea_out, "--out", replay_out],
        capture_output=True, text=True)
    if replay_result.returncode != 0:
        raise SystemExit(
            f"MISMATCH: replay from frozen artifacts did not reproduce the original "
            f"result:\n{replay_result.stdout}{replay_result.stderr}")
    print(replay_result.stdout.strip())
    print(f"\nAll artifacts are under: {os.path.abspath(args.out_dir)}")


if __name__ == "__main__":
    main()
