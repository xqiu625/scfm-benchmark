#!/usr/bin/env python3
"""Option-C scGPT cell embedding on the canonical COVID input.

Generates the scGPT Option-C embedding that was deferred on 2026-05-01 (scGPT dropped
from Option C; entries shown as "—"). Same canonical raw-counts input as every other FM,
so the result is directly comparable. Run INSIDE the scGPT container:

    singularity exec /…/07_CONTAINERS/scgpt_0.1.7.sif python embed_scgpt_option_c.py \
        --input  /…/01_DATA/raw/human_platelet_covid_severity_raw_counts.h5ad \
        --model  /…/models/scGPT_human/scGPT_human \
        --output /…/02_EMBEDDINGS/option_c/scgpt_covid_canonical.h5ad \
        --batch-size 32 [--n-cells 1000]   # --n-cells for a quick smoke test

Uses scGPT's high-level `embed_data` (handles gene-vocab match + binning + forward) rather
than re-implementing the legacy manual pipeline. Preserves the canonical obs columns the
downstream Option-C benchmark needs (Cell_ID / Disease / covid_severity / Sample_ID / Data_NO)
and writes the cell embedding to obsm['X_scgpt'].
"""
import argparse
import sys

import numpy as np
import scanpy as sc

OBS_KEEP = ["Cell_ID", "Disease", "covid_severity", "Sample_ID", "Data_NO"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--model", required=True, help="scGPT model dir (vocab.json/best_model.pt/args.json)")
    ap.add_argument("--output", required=True)
    ap.add_argument("--batch-size", type=int, default=32)
    ap.add_argument("--n-cells", type=int, default=0, help="if >0, subsample this many cells (smoke test)")
    ap.add_argument("--gene-col", default="index", help="adata.var column with gene symbols (default: index)")
    args = ap.parse_args()

    import torch
    from scgpt.tasks import embed_data

    print(f"[scgpt-optc] torch {torch.__version__} cuda={torch.cuda.is_available()}", flush=True)
    adata = sc.read_h5ad(args.input)
    print(f"[scgpt-optc] loaded {adata.shape} from {args.input}", flush=True)

    if args.n_cells and args.n_cells < adata.n_obs:
        rng = np.random.default_rng(0)
        idx = np.sort(rng.choice(adata.n_obs, args.n_cells, replace=False))
        adata = adata[idx].copy()
        print(f"[scgpt-optc] SMOKE TEST subset -> {adata.shape}", flush=True)

    # Cast counts to float32 — embed_data builds the expression tensor from X's dtype and runs
    # under torch.cuda.amp.autocast; a float64 input reaches a Linear that autocast won't
    # downcast, giving "mat1 and mat2 must have the same dtype".
    adata.X = adata.X.astype(np.float32)

    # The container's embed_data asserts cell_type_key in adata.obs (used only for metadata,
    # not the embedding). All cells here are platelets — add a constant column to satisfy it.
    adata.obs["cell_type"] = "platelet"

    # gene_col='index' tokenizes against the model vocab on the var index (gene symbols);
    # return_new_adata=False adds obsm['X_scGPT'] to the (vocab-gene-subset) adata, keeping
    # all cells + obs so the canonical obs (Cell_ID/Disease/covid_severity/Sample_ID/Data_NO)
    # carry through to the downstream benchmark.
    result = embed_data(
        adata,
        args.model,
        gene_col=args.gene_col,
        batch_size=args.batch_size,
        return_new_adata=False,
    )
    X = np.asarray(result.obsm["X_scGPT"], dtype="float32")
    print(f"[scgpt-optc] embedding {X.shape}  nan={np.isnan(X).any()}  "
          f"global_std={X.std():.4f}  obs_cols={list(result.obs.columns)}", flush=True)

    result.obsm["X_scgpt"] = X          # standardized lowercase key for the benchmark
    if "X_scGPT" in result.obsm:
        del result.obsm["X_scGPT"]
    result.write_h5ad(args.output)
    print(f"[scgpt-optc] wrote {args.output} (obsm['X_scgpt'] = {X.shape}, cells={result.n_obs})", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
