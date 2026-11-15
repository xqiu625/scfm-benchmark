#!/usr/bin/env python3
"""GSEA on per-FM "fatal direction" gene rankings (component B of biological interp).

Input: top_genes_<fm>_dim<d>_rank<r>.csv from biological_interpretation_option_c.py
       (each contains per-gene Spearman correlation with one of the FM's top
       LR-weighted dimensions on the mortality binary task)

Per FM × top dim, run GSEA against:
  - KEGG_2021_Human
  - Reactome_2022
  - GO_Biological_Process_2023

Output: gsea_<fm>_dim<d>_rank<r>_<library>.csv

Uses gseapy.prerank — feeds the SIGNED Spearman rho per gene as the rank metric,
so positively-correlated genes drive enrichment for the "high-fatal-axis" pathway
direction and negatively-correlated genes drive the opposite.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("gsea")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-dir", required=True, type=Path,
                   help="Dir containing top_genes_<fm>_dim<d>_rank<r>.csv files")
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--libraries", nargs="+",
                   default=["KEGG_2021_Human", "Reactome_2022", "GO_Biological_Process_2023"],
                   help="Enrichr library names")
    p.add_argument("--overwrite", action="store_true",
                   help="Recompute even if the output CSV already exists. Required "
                        "when the input ranking changed but the output directory is "
                        "reused, otherwise every term is silently served from cache.")
    p.add_argument("--top-n", type=int, default=0,
                   help="Top genes per dim to use as the pre-ranked gene list")
    p.add_argument("--min-set-size", type=int, default=10)
    p.add_argument("--max-set-size", type=int, default=500)
    args = p.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    try:
        import gseapy as gp
    except ImportError:
        log.error("gseapy not installed in this env. Install via:")
        log.error("  uv pip install --python <env>/bin/python gseapy")
        sys.exit(2)

    gene_files = sorted(args.input_dir.glob("top_genes_*_dim*_rank*.csv"))
    if not gene_files:
        log.error("No top_genes files in %s", args.input_dir); sys.exit(1)
    log.info("Found %d top-gene files", len(gene_files))

    summary_rows = []

    for f in gene_files:
        # Parse filename: top_genes_<fm>_dim<d>_rank<r>.csv
        stem = f.stem  # top_genes_<fm>_dim<d>_rank<r>
        parts = stem.split("_")
        try:
            mk = parts[2]
            dim_idx = parts[3].replace("dim", "")
            rank = parts[4].replace("rank", "")
        except IndexError:
            log.warning("Cannot parse %s; skipping", f.name); continue

        df = pd.read_csv(f)
        if "gene" not in df.columns or "spearman_rho" not in df.columns:
            log.warning("%s missing required columns", f.name); continue

        # Use signed rho as rank metric (positive = aligned with dim, negative = anti-aligned)
        rnk = pd.DataFrame({"gene": df["gene"], "rank_metric": df["spearman_rho"]}) \
            .sort_values("rank_metric", ascending=False) \
            .reset_index(drop=True)
        # GSEA prerank requires the FULL signed ranking: the enrichment score walks the
        # whole list and the permutation null is drawn from it. Truncating to the top-N
        # (as this script did until 2026-08-06) discards the entire negatively-correlated
        # tail and computes NES/FDR against a few-hundred-gene universe, which makes both
        # statistics invalid. --top-n is retained for debugging only; 0 = no truncation.
        if args.top_n and args.top_n > 0:
            log.warning("TRUNCATING ranking to top %d - NES/FDR will NOT be valid", args.top_n)
            rnk = rnk.head(args.top_n).reset_index(drop=True)

        log.info("[%s dim%s rank%s] n=%d genes; rank_metric range=[%.3f, %.3f]",
                 mk, dim_idx, rank, len(rnk), rnk["rank_metric"].min(), rnk["rank_metric"].max())

        for lib in args.libraries:
            outname = f"gsea_{mk}_dim{dim_idx}_rank{rank}_{lib.replace('/', '-')}.csv"
            outpath = args.output_dir / outname
            if outpath.exists() and not args.overwrite:
                log.warning("  [%s] SKIPPED (file exists): %s -- pass --overwrite to "
                            "recompute. If the existing file came from a different "
                            "ranking, this silently keeps the stale result.",
                            lib, outname)
                continue
            try:
                pre = gp.prerank(
                    rnk=rnk,
                    gene_sets=lib,
                    outdir=None,
                    min_size=args.min_set_size,
                    max_size=args.max_set_size,
                    permutation_num=200,
                    seed=42,
                    verbose=False,
                )
                res = pre.res2d
                res.to_csv(outpath, index=False)
                top_term = res.iloc[0]["Term"] if len(res) else "(no enrichment)"
                top_nes = float(res.iloc[0]["NES"]) if len(res) else float("nan")
                log.info("  [%s] %d terms; top: %s (NES=%.2f)", lib, len(res), top_term, top_nes)
                summary_rows.append(dict(model=mk, dim=dim_idx, rank=rank, library=lib,
                                         n_terms=len(res), top_term=top_term, top_nes=top_nes,
                                         output_csv=outname))
            except Exception as e:
                log.exception("  [%s] failed: %s", lib, e)

    pd.DataFrame(summary_rows).to_csv(args.output_dir / "gsea_summary.csv", index=False)
    log.info("Wrote gsea_summary.csv with %d rows", len(summary_rows))


if __name__ == "__main__":
    main()
