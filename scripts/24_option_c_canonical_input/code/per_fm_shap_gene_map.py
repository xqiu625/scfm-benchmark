#!/usr/bin/env python3
"""Companion to per_fm_shap_attribution.py — gene-mapping only.

Reads existing shap_dims_<fm>.csv (per-FM top embedding dims by mean |SHAP|)
and produces shap_genes_<fm>.csv via Spearman correlation between each top
embedding dim and raw gene expression.

Skips the expensive 5-fold SHAP recomputation, which already succeeded.
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import spearmanr

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("shap_gene_map")

OBSM_KEYS = {
    "state": ["X_state", "X_emb", "X_STATE"],
    "uce": ["X_uce", "X_uce_4layer"],
    "geneformer": ["X_geneformer", "X_Geneformer"],
    "transcriptformer": ["X_transcriptformer", "X_TranscriptFormer"],
    "scprint2": ["X_scprint2", "X_scprint", "X_scPRINT", "X_scPRINT2", "scprint_emb"],
    "stack": ["X_stack", "X_STACK"],
    "scgpt": ["X_scgpt", "X_scGPT"],
}

TOP_GENES = 30
SPEARMAN_SUBSAMPLE = 10000
SEED = 42


def find_embedding(adata, model_key: str):
    for k in OBSM_KEYS.get(model_key, []):
        if k in adata.obsm:
            return np.asarray(adata.obsm[k]), k
    raise RuntimeError(f"No known obsm key for {model_key}; available: {list(adata.obsm.keys())}")


def jaccard_top1(genes_dfs: dict[str, pd.DataFrame], top_n=TOP_GENES):
    fms = sorted(genes_dfs.keys())
    sets = {fm: set(df.query("dim_rank == 1")["gene"].head(top_n)) for fm, df in genes_dfs.items()}
    rows = []
    for i, a in enumerate(fms):
        for b in fms[i:]:
            inter = len(sets[a] & sets[b]); union = len(sets[a] | sets[b])
            rows.append({"fm_a": a, "fm_b": b, "jaccard": inter / max(union, 1),
                         "n_intersect": inter, "n_union": union})
    return pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--shap-dir", required=True, type=Path,
                   help="Dir containing shap_dims_<fm>.csv from previous run")
    p.add_argument("--embedding-dir", required=True, type=Path)
    p.add_argument("--raw-input", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--fms", nargs="+",
                   default=["state", "uce", "transcriptformer", "scprint2", "stack", "geneformer",
                            "scgpt"])
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    log.info("Loading raw input: %s", args.raw_input)
    raw_adata = sc.read_h5ad(args.raw_input)
    log.info("  raw shape: %d × %d", raw_adata.n_obs, raw_adata.n_vars)
    gene_names = raw_adata.var_names.astype(str).to_numpy()

    rng = np.random.default_rng(SEED)
    genes_dfs = {}
    summary_rows = []

    for fm in args.fms:
        dims_csv = args.shap_dir / f"shap_dims_{fm}.csv"
        if not dims_csv.exists():
            log.warning("Missing %s; skip", dims_csv); continue
        df_dims = pd.read_csv(dims_csv)
        dim_order = df_dims["dim_index"].astype(int).to_numpy()
        log.info("%s: top dims = %s", fm, dim_order.tolist())

        cands = list(args.embedding_dir.glob(f"{fm}_*canonical*.h5ad"))
        if not cands:
            log.warning("No h5ad for %s; skip", fm); continue
        a = sc.read_h5ad(cands[0])
        X, obsm_key = find_embedding(a, fm)
        log.info("  obsm=%s shape=%s", obsm_key, X.shape)

        if "Cell_ID" not in a.obs.columns or "Cell_ID" not in raw_adata.obs.columns:
            log.warning("  no Cell_ID; skip"); continue

        a_cid = a.obs["Cell_ID"].astype(str).values
        r_cid = raw_adata.obs["Cell_ID"].astype(str).values
        a2pos = {c: i for i, c in enumerate(a_cid)}
        keep_pairs = [(a2pos[c], i) for i, c in enumerate(r_cid) if c in a2pos]
        a_idx = np.array([p[0] for p in keep_pairs])
        r_idx = np.array([p[1] for p in keep_pairs])
        log.info("  Cell_ID-aligned: %d", len(a_idx))

        if len(a_idx) > SPEARMAN_SUBSAMPLE:
            sel = rng.choice(len(a_idx), SPEARMAN_SUBSAMPLE, replace=False)
            a_idx, r_idx = a_idx[sel], r_idx[sel]

        raw_X = raw_adata.X[r_idx]
        if sp.issparse(raw_X):
            raw_X = raw_X.toarray()

        rows = []
        for rank_i, d in enumerate(dim_order):
            emb_col = X[a_idx, int(d)]
            rho = np.zeros(raw_X.shape[1], dtype=np.float32)
            for g in range(raw_X.shape[1]):
                col = raw_X[:, g]
                if col.std() == 0:
                    continue
                try:
                    rho[g] = spearmanr(emb_col, col).correlation
                except Exception:
                    rho[g] = np.nan
            rho_abs = np.abs(np.nan_to_num(rho))
            top = np.argsort(-rho_abs)[:TOP_GENES]
            for k, gi in enumerate(top):
                rows.append({
                    "fm": fm, "dim_rank": rank_i + 1, "dim_index": int(d),
                    "gene_rank": k + 1, "gene": gene_names[gi],
                    "spearman_rho": float(rho[gi]), "abs_rho": float(rho_abs[gi]),
                })
            log.info("  dim_rank=%d dim=%d top gene = %s (rho=%.3f)",
                     rank_i + 1, int(d), gene_names[top[0]], float(rho[top[0]]))

        df_genes = pd.DataFrame(rows)
        df_genes.to_csv(args.output_dir / f"shap_genes_{fm}.csv", index=False)
        genes_dfs[fm] = df_genes

        top_gene = df_genes.query("dim_rank == 1").iloc[0]
        summary_rows.append({
            "fm": fm,
            "top_dim": int(df_dims.iloc[0]["dim_index"]),
            "top_dim_mean_abs_shap": float(df_dims.iloc[0]["mean_abs_shap"]),
            "top_gene": top_gene["gene"], "top_gene_rho": float(top_gene["spearman_rho"]),
        })

    if summary_rows:
        pd.DataFrame(summary_rows).to_csv(args.output_dir / "shap_summary.csv", index=False)
    if len(genes_dfs) >= 2:
        jaccard_top1(genes_dfs).to_csv(args.output_dir / "shap_jaccard_overlap.csv", index=False)

    log.info("Done. %d FMs processed.", len(genes_dfs))


if __name__ == "__main__":
    main()
