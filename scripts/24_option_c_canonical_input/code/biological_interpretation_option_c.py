#!/usr/bin/env python3
"""Biological interpretation of Option C 5-FM embeddings (mortality binary).

Three components per FM:
  D. UMAP coordinates of the FM embedding (saved as CSV: cell × UMAP1, UMAP2,
     plus mortality, cohort, severity labels for downstream coloring)
  A. Per-FM "top mortality dimensions" — train an LR on the full embedding,
     extract top-K dimensions by absolute coefficient. For each top dim,
     compute Spearman correlation between every input gene's expression and
     that dim's value. Top 100 genes by |spearman| per top dim become the
     candidate biomarker set.
  B (deferred): GSEA on the per-dim gene rankings (separate script).

Outputs (in --output-dir):
  umap_<fm>.csv                cell × (UMAP1, UMAP2, mortality, cohort, severity)
  top_dims_<fm>.csv            dim_idx, lr_coef, |coef|
  top_genes_<fm>_dim<d>.csv    gene, spearman_rho, |rho|, p_value (full ranking)
  summary_<fm>.json            run metadata + LR train AUC
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("bio_interp")

OBSM_KEYS = {
    "state": ["X_state", "X_emb", "X_STATE"],
    "uce": ["X_uce", "X_uce_4layer"],
    "geneformer": ["X_geneformer", "X_Geneformer"],
    "transcriptformer": ["X_transcriptformer", "X_TranscriptFormer"],
    "scprint2": ["X_scprint2", "X_scprint", "X_scPRINT", "X_scPRINT2", "scprint_emb"],
    "stack": ["X_stack", "X_STACK"],
    "scgpt": ["X_scgpt", "X_scGPT"],
}

TOP_K_DIMS = 5
TOP_N_GENES = 0   # 0 = write the full ranked list (required for valid GSEA)
SPEARMAN_SUBSAMPLE = 10000  # cells used for spearman (full 47K is slow; subsample)


def find_embedding(adata, model_key: str):
    """Find an FM embedding by known obsm keys. Refuses the previous unsafe
    fallback (first 2D obsm) -- that could silently grab UMAP coords or other
    auxiliary arrays. Raises a clear error if no known key is present.
    Mirrors the policy in primary_benchmark_option_c.py."""
    keys = OBSM_KEYS.get(model_key, [])
    for k in keys:
        if k in adata.obsm:
            return np.asarray(adata.obsm[k]), k
    raise RuntimeError(
        f"No known embedding key for model_key={model_key!r}. "
        f"Tried {keys}; available obsm: {list(adata.obsm.keys())}. "
        f"If a new FM is added, register its obsm key(s) in OBSM_KEYS."
    )


def build_mortality(severity: np.ndarray):
    sev = np.asarray([str(s).lower() for s in severity])
    fatal = sev == "fatal"
    survived = np.isin(sev, ["mild", "moderate", "recovered", "severe"])
    valid = fatal | survived
    y = np.zeros(len(sev), dtype=np.int8)
    y[fatal] = 1
    return y, valid


def per_fm(fm_h5ad: Path, raw_adata, output_dir: Path):
    mk = fm_h5ad.stem.split("_")[0].lower()
    log.info("=" * 60); log.info("Model: %s", mk); log.info("=" * 60)

    a = sc.read_h5ad(fm_h5ad)
    X, obsm_key = find_embedding(a, mk)
    log.info("  obsm=%s shape=%s", obsm_key, X.shape)

    if "covid_severity" not in a.obs.columns or "Sample_ID" not in a.obs.columns:
        log.warning("  missing severity/Sample_ID; skip"); return None
    if "Cell_ID" not in a.obs.columns:
        log.warning("  missing Cell_ID; skip"); return None

    severity = a.obs["covid_severity"].astype(str).values
    cohort = a.obs["Data_NO"].astype(str).values if "Data_NO" in a.obs.columns else None
    cell_ids = a.obs["Cell_ID"].astype(str).values
    y_full, valid = build_mortality(severity)
    X_v = X[valid]
    y_v = y_full[valid]
    cell_ids_v = cell_ids[valid]
    sev_v = severity[valid]
    cohort_v = cohort[valid] if cohort is not None else None
    log.info("  CV input: n=%d, fatal=%d, survived=%d", len(y_v), int(y_v.sum()), int((y_v == 0).sum()))

    # --- D. UMAP ---
    log.info("[D] Computing UMAP")
    a_emb = sc.AnnData(X=X_v)
    a_emb.obsm["X_emb"] = X_v
    sc.pp.neighbors(a_emb, n_neighbors=30, use_rep="X_emb", random_state=42)
    sc.tl.umap(a_emb, random_state=42)
    umap_coords = a_emb.obsm["X_umap"]
    df_umap = pd.DataFrame({
        "Cell_ID": cell_ids_v,
        "UMAP1": umap_coords[:, 0],
        "UMAP2": umap_coords[:, 1],
        "mortality": np.where(y_v == 1, "fatal", "survived"),
        "severity": sev_v,
        "cohort": cohort_v if cohort_v is not None else "unknown",
    })
    df_umap.to_csv(output_dir / f"umap_{mk}.csv", index=False)
    log.info("  Wrote umap_%s.csv (n=%d)", mk, len(df_umap))

    # --- A. Top mortality dimensions ---
    log.info("[A] Training full-data LR for interpretation")
    scaler = StandardScaler().fit(X_v)
    X_std = scaler.transform(X_v)
    lr = LogisticRegression(penalty="l2", C=1.0, solver="lbfgs",
                            max_iter=2000, class_weight="balanced", n_jobs=4)
    lr.fit(X_std, y_v)
    train_auc = roc_auc_score(y_v, lr.predict_proba(X_std)[:, 1])
    log.info("  LR train AUC (no CV; for interpretation only): %.4f", train_auc)

    coef = lr.coef_.ravel()
    abs_coef = np.abs(coef)
    top_idx = np.argsort(-abs_coef)[:TOP_K_DIMS]
    df_top = pd.DataFrame({
        "rank": np.arange(1, TOP_K_DIMS + 1),
        "dim_idx": top_idx,
        "lr_coef": coef[top_idx],
        "abs_coef": abs_coef[top_idx],
    })
    df_top.to_csv(output_dir / f"top_dims_{mk}.csv", index=False)
    log.info("  Top dims: %s", df_top["dim_idx"].tolist())

    # --- A.cont: gene-loading per top dim via Spearman ---
    # Need gene expression aligned to the same cells. Use raw_adata.X with Cell_ID alignment.
    log.info("[A] Computing per-dim gene loadings (Spearman, subsample n=%d)", SPEARMAN_SUBSAMPLE)
    # Align FM cells to raw_adata via Cell_ID
    raw_cellids = raw_adata.obs["Cell_ID"].astype(str).values if "Cell_ID" in raw_adata.obs.columns else raw_adata.obs.index.values
    raw_cid_to_pos = {c: i for i, c in enumerate(raw_cellids)}
    fm_to_raw_pos = []
    keep_emb_pos = []
    for emb_pos, cid in enumerate(cell_ids_v):
        rp = raw_cid_to_pos.get(cid)
        if rp is not None:
            fm_to_raw_pos.append(rp)
            keep_emb_pos.append(emb_pos)
    log.info("  Aligned %d cells to raw expression", len(fm_to_raw_pos))
    if len(fm_to_raw_pos) < 1000:
        log.warning("  Too few aligned cells; skipping gene-loading")
        return dict(model=mk, n_cells=len(y_v), train_auc=train_auc, n_aligned=len(fm_to_raw_pos))

    # Subsample for spearman speed
    rng = np.random.RandomState(42)
    if len(keep_emb_pos) > SPEARMAN_SUBSAMPLE:
        sub = rng.choice(len(keep_emb_pos), size=SPEARMAN_SUBSAMPLE, replace=False)
        emb_sub = [keep_emb_pos[i] for i in sub]
        raw_sub = [fm_to_raw_pos[i] for i in sub]
    else:
        emb_sub = keep_emb_pos
        raw_sub = fm_to_raw_pos
    X_sub = X_v[emb_sub]              # (n_sub, n_dim)
    Xexpr_sub = raw_adata.X[raw_sub]  # (n_sub, n_genes), sparse
    if hasattr(Xexpr_sub, "toarray"):
        Xexpr_sub = Xexpr_sub.toarray()
    gene_names = np.asarray(raw_adata.var_names)
    n_genes = Xexpr_sub.shape[1]
    log.info("  Sub n=%d cells, n_genes=%d", X_sub.shape[0], n_genes)

    for d_rank, dim_idx in enumerate(top_idx, start=1):
        log.info("  [dim rank=%d, idx=%d] computing %d Spearman correlations", d_rank, dim_idx, n_genes)
        dim_vals = X_sub[:, dim_idx]
        rhos = np.zeros(n_genes, dtype=np.float32)
        pvals = np.ones(n_genes, dtype=np.float32)
        for g in range(n_genes):
            col = Xexpr_sub[:, g]
            if col.std() == 0:
                continue
            r, p = spearmanr(col, dim_vals)
            rhos[g] = r if not np.isnan(r) else 0.0
            pvals[g] = p if not np.isnan(p) else 1.0
        df_g = pd.DataFrame({
            "gene": gene_names,
            "spearman_rho": rhos,
            "abs_rho": np.abs(rhos),
            "p_value": pvals,
        }).sort_values("abs_rho", ascending=False)
        # Write the FULL ranked gene list. gsea_per_fm_dim.py consumes this file as the
        # prerank input, and GSEA needs the whole signed ranking: the enrichment score
        # walks the entire list and the permutation null is drawn from it. Truncating
        # here to the top 200 by |rho| (the behaviour until 2026-08-06) silently made
        # every downstream NES and FDR invalid, and removing the truncation inside the
        # GSEA script alone did not help because the truncation was already baked into
        # its input. Consumers that want a short display list should take .head(N)
        # themselves; SHAP/Jaccard top-30 lists are unaffected either way.
        if TOP_N_GENES and TOP_N_GENES > 0:
            log.warning("TRUNCATING gene list to top %d - GSEA NES/FDR will NOT be valid",
                        TOP_N_GENES)
            df_g = df_g.head(TOP_N_GENES)
        df_g.to_csv(output_dir / f"top_genes_{mk}_dim{dim_idx}_rank{d_rank}.csv", index=False)
        log.info("    Top 5 genes by |rho|: %s", df_g["gene"].head(5).tolist())

    # Summary metadata
    summary = dict(
        model=mk,
        obsm_key=obsm_key,
        n_cells_total=int(len(y_v)),
        n_aligned_to_raw=int(len(fm_to_raw_pos)),
        n_subsampled_for_spearman=int(X_sub.shape[0]),
        embedding_dim=int(X_v.shape[1]),
        n_genes=int(n_genes),
        train_auc=float(train_auc),
        top_dims=top_idx.tolist(),
        top_dim_coefs=coef[top_idx].tolist(),
    )
    with open(output_dir / f"summary_{mk}.json", "w") as f:
        json.dump(summary, f, indent=2)
    return summary


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--embedding-dir", required=True, type=Path)
    p.add_argument("--raw-input", required=True, type=Path,
                   help="Canonical raw input AnnData with gene expression (for per-dim gene loadings)")
    p.add_argument("--output-dir", required=True, type=Path)
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    log.info("Loading raw input: %s", args.raw_input)
    raw = sc.read_h5ad(args.raw_input)
    log.info("  raw: %d cells, %d genes", raw.n_obs, raw.n_vars)

    fm_files = sorted(args.embedding_dir.glob("*_covid_canonical.h5ad"))
    if not fm_files:
        log.error("No FMs in %s", args.embedding_dir); sys.exit(1)

    summaries = []
    for fm in fm_files:
        try:
            s = per_fm(fm, raw, args.output_dir)
            if s is not None:
                summaries.append(s)
        except Exception as e:
            log.exception("Failed %s: %s", fm.name, e)

    pd.DataFrame(summaries).to_csv(args.output_dir / "summaries_all.csv", index=False)
    log.info("Done. Wrote %d FM summaries.", len(summaries))


if __name__ == "__main__":
    main()
