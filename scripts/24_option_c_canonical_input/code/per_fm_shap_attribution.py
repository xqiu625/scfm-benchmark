#!/usr/bin/env python3
"""Per-FM SHAP attribution on the in-hospital mortality LR classifier.

For each Option C FM:
  1. Load the FM's COVID embedding (02_EMBEDDINGS/option_c/<fm>_covid_canonical.h5ad)
  2. Train a logistic-regression mortality classifier on out-of-fold splits
     (StratifiedGroupKFold on Sample_ID, mirroring tab:performance)
  3. Run shap.LinearExplainer on the LR classifier — yields per-cell, per-dim
     attribution
  4. Aggregate to mean |SHAP| per embedding dimension; pick top dimensions
  5. Map dimensions back to genes via Spearman correlation between each top
     embedding dimension and each raw-count gene; take top-30 genes per FM

Outputs (per FM):
  - shap_dims_<fm>.csv    : embedding-dimension importance ranking
  - shap_genes_<fm>.csv   : top-30 gene attribution
Plus aggregate:
  - shap_jaccard_overlap.csv : pairwise Jaccard of top-30 gene sets across FMs
  - shap_summary.csv        : per-FM | top dim | top gene | mean |SHAP|

Notes:
  - LinearExplainer is exact for LR; no Monte Carlo sampling needed.
  - Spearman computed on a 10K-cell subsample for speed (full 47K is slow).
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
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("shap_attr")

OBSM_KEYS = {
    "state": ["X_state", "X_emb", "X_STATE"],
    "uce": ["X_uce", "X_uce_4layer"],
    "geneformer": ["X_geneformer", "X_Geneformer"],
    "transcriptformer": ["X_transcriptformer", "X_TranscriptFormer"],
    "scprint2": ["X_scprint2", "X_scprint", "X_scPRINT", "X_scPRINT2", "scprint_emb"],
    "stack": ["X_stack", "X_STACK"],
    "scgpt": ["X_scgpt", "X_scGPT"],
}

TOP_DIMS = 10
TOP_GENES = 30
SPEARMAN_SUBSAMPLE = 10000
SEED = 42


def find_embedding(adata, model_key: str):
    for k in OBSM_KEYS.get(model_key, []):
        if k in adata.obsm:
            return np.asarray(adata.obsm[k]), k
    raise RuntimeError(f"No known obsm key for {model_key}; available: {list(adata.obsm.keys())}")


def build_mortality(severity: np.ndarray):
    sev = np.asarray([str(s).lower().strip() for s in severity])
    fatal = sev == "fatal"
    survived = np.isin(sev, ["mild", "moderate", "recovered", "severe"])
    valid = fatal | survived
    y = np.zeros(len(sev), dtype=np.int8)
    y[fatal] = 1
    return y, valid


def per_fm(fm_h5ad: Path, raw_adata, output_dir: Path, top_dims=TOP_DIMS, top_genes=TOP_GENES):
    mk = fm_h5ad.stem.split("_")[0].lower()
    log.info("=" * 60); log.info("Model: %s", mk); log.info("=" * 60)

    a = sc.read_h5ad(fm_h5ad)
    X, obsm_key = find_embedding(a, mk)
    log.info("  obsm=%s shape=%s", obsm_key, X.shape)

    if "covid_severity" not in a.obs.columns or "Sample_ID" not in a.obs.columns:
        log.warning("  missing severity/Sample_ID; skip"); return None

    y_full, valid = build_mortality(a.obs["covid_severity"].astype(str).values)
    X_v = X[valid]
    y_v = y_full[valid]
    groups = a.obs.loc[valid, "Sample_ID"].astype(str).values
    log.info("  valid cells=%d (fatal=%d, survived=%d, groups=%d)",
             len(y_v), int(y_v.sum()), int((1 - y_v).sum()), len(np.unique(groups)))

    # Train OOF LR; collect all per-cell SHAP attributions
    try:
        import shap
    except ImportError:
        log.error("shap not installed. uv pip install --python <env>/bin/python shap"); raise

    skf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=SEED)
    oof_shap = np.zeros_like(X_v, dtype=np.float32)
    oof_pred = np.zeros(len(y_v), dtype=np.float32)

    for fold, (tr, va) in enumerate(skf.split(np.zeros(len(y_v)), y_v, groups=groups)):
        scaler = StandardScaler().fit(X_v[tr])
        Xtr = scaler.transform(X_v[tr])
        Xva = scaler.transform(X_v[va])
        clf = LogisticRegression(max_iter=5000, C=1.0, solver="liblinear",
                                 class_weight="balanced", random_state=SEED).fit(Xtr, y_v[tr])
        oof_pred[va] = clf.predict_proba(Xva)[:, 1]

        explainer = shap.LinearExplainer(clf, Xtr, feature_perturbation="interventional")
        sv = explainer.shap_values(Xva)
        # Some shap versions return list for binary; take positive class
        if isinstance(sv, list):
            sv = sv[1]
        oof_shap[va] = sv.astype(np.float32)
        log.info("  fold %d: train=%d test=%d", fold + 1, len(tr), len(va))

    # Aggregate to per-dim mean |SHAP|
    dim_importance = np.abs(oof_shap).mean(axis=0)
    dim_order = np.argsort(-dim_importance)[:top_dims]
    df_dims = pd.DataFrame({
        "rank": np.arange(1, len(dim_order) + 1),
        "dim_index": dim_order,
        "mean_abs_shap": dim_importance[dim_order],
    })
    df_dims.to_csv(output_dir / f"shap_dims_{mk}.csv", index=False)
    log.info("  Top dim: %d (mean|SHAP|=%.4f)", int(dim_order[0]), float(dim_importance[dim_order[0]]))

    # Map top dims back to genes via Spearman
    if raw_adata is None:
        log.warning("  no raw adata supplied; skipping gene mapping"); return df_dims, None

    # Align cells via Cell_ID. There is no safe fallback: if Cell_ID is missing
    # on either side we cannot pair embedding rows with raw-count rows, and the
    # previous index-intersection branch left a_idx/r_idx undefined (NameError).
    if "Cell_ID" not in a.obs.columns or "Cell_ID" not in raw_adata.obs.columns:
        raise RuntimeError(
            f"Cannot map SHAP dims to genes for {mk}: 'Cell_ID' missing in obs "
            f"(embedding has it: {'Cell_ID' in a.obs.columns}; raw has it: "
            f"{'Cell_ID' in raw_adata.obs.columns}). Cell_ID-based alignment is "
            f"required; re-run with inputs that carry Cell_ID."
        )
    a_cid = a.obs["Cell_ID"].astype(str).values
    r_cid = raw_adata.obs["Cell_ID"].astype(str).values
    a2pos = {c: i for i, c in enumerate(a_cid)}
    keep_pairs = [(a2pos[c], i) for i, c in enumerate(r_cid) if c in a2pos]
    a_idx = np.array([p[0] for p in keep_pairs])
    r_idx = np.array([p[1] for p in keep_pairs])
    if len(a_idx) == 0:
        raise RuntimeError(
            f"Cannot map SHAP dims to genes for {mk}: zero cells shared between "
            f"the embedding and the raw input after Cell_ID alignment."
        )

    # Subsample for Spearman
    rng = np.random.default_rng(SEED)
    if len(a_idx) > SPEARMAN_SUBSAMPLE:
        sel = rng.choice(len(a_idx), SPEARMAN_SUBSAMPLE, replace=False)
        a_idx, r_idx = a_idx[sel], r_idx[sel]
    log.info("  Cell_ID-aligned cells for Spearman: %d", len(a_idx))

    raw_X = raw_adata.X[r_idx]
    if sp.issparse(raw_X):
        raw_X = raw_X.toarray()
    gene_names = raw_adata.var_names.astype(str).to_numpy()

    rows = []
    for rank_i, d in enumerate(dim_order):
        # a_idx indexes into the FULL adata (a), so use X (not X_v which is valid-filtered).
        emb_col = X[a_idx, d]
        # Vectorize: spearman per gene against this dim
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
        top = np.argsort(-rho_abs)[:top_genes]
        for k, gi in enumerate(top):
            rows.append({
                "fm": mk, "dim_rank": rank_i + 1, "dim_index": int(d),
                "gene_rank": k + 1, "gene": gene_names[gi],
                "spearman_rho": float(rho[gi]), "abs_rho": float(rho_abs[gi]),
            })

    df_genes = pd.DataFrame(rows)
    df_genes.to_csv(output_dir / f"shap_genes_{mk}.csv", index=False)
    return df_dims, df_genes


def jaccard_top1(genes_dfs: dict[str, pd.DataFrame], top_n=30):
    fms = sorted(genes_dfs.keys())
    sets = {fm: set(df.query("dim_rank == 1")["gene"].head(top_n)) for fm, df in genes_dfs.items()}
    rows = []
    for i, a in enumerate(fms):
        for b in fms[i:]:
            inter = len(sets[a] & sets[b]); union = len(sets[a] | sets[b])
            rows.append({"fm_a": a, "fm_b": b, "jaccard": inter / max(union, 1), "n_intersect": inter, "n_union": union})
    return pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--embedding-dir", required=True, type=Path)
    p.add_argument("--raw-input", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--fms", nargs="+",
                   default=["state", "uce", "transcriptformer", "scprint2", "stack", "geneformer", "scgpt"])
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    log.info("Loading raw input: %s", args.raw_input)
    raw_adata = sc.read_h5ad(args.raw_input)
    log.info("  raw shape: %d × %d", raw_adata.n_obs, raw_adata.n_vars)

    genes_dfs = {}
    summary_rows = []
    for fm in args.fms:
        candidates = list(args.embedding_dir.glob(f"{fm}_*canonical*.h5ad"))
        if not candidates:
            log.warning("No h5ad for %s in %s — skip", fm, args.embedding_dir); continue
        try:
            df_dims, df_genes = per_fm(candidates[0], raw_adata, args.output_dir)
        except Exception as e:
            log.error("Failed %s: %s", fm, e); continue
        if df_genes is not None:
            genes_dfs[fm] = df_genes
            top_gene = df_genes.query("dim_rank == 1").iloc[0]
            summary_rows.append({
                "fm": fm, "top_dim": int(df_dims.iloc[0]["dim_index"]),
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
