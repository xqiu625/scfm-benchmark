#!/usr/bin/env python3
"""Platelet activation marker analysis on FM embeddings (Option C, all 5 FMs + scPRINT-2).

For each FM's embedding:
  1. Compute platelet-activation module score per cell (sc.tl.score_genes on raw counts)
  2. Spearman correlation between every embedding dimension and the module score
  3. Rank dimensions by |correlation|; identify the top "activation axis" per FM
  4. Compare top-axis values between fatal and survived cells (Wilcoxon)
  5. Output per-FM CSV: {model, dim, spearman, abs_spearman, fatal_vs_survived_pval}
  6. Output a summary CSV: top dim per FM + module-score AUROC for fatal-vs-survived

Reference markers (literature consensus, COVID/sepsis platelets):
  Core activation:  PF4, PPBP, SELP, ITGA2B, ITGB3, GP1BA, GP1BB, GP9
  Hyper-reactivity: CCL5, CXCL4, THBS1, SERPINE1, F5
  Mitochondrial:    MT-CO1, MT-CO2, MT-CO3
  Immune:           IFITM1, IFITM2, IFITM3, HLA-A, HLA-B, HLA-C

Usage:
    python score_activation_markers.py \
        --embedding-dir /bigdata/.../02_EMBEDDINGS/option_c \
        --raw-input /bigdata/.../01_DATA/raw/human_platelet_covid_severity_raw_counts.h5ad \
        --label-col covid_severity \
        --output-dir /bigdata/.../05_RESULTS/option_c_activation_markers
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.metrics import roc_auc_score

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("activation_markers")


PLATELET_ACTIVATION_GENES = {
    "core_activation": ["PF4", "PPBP", "SELP", "ITGA2B", "ITGB3", "GP1BA", "GP1BB", "GP9"],
    "hyper_reactivity": ["CCL5", "CXCL4", "THBS1", "SERPINE1", "F5"],
    "mitochondrial": ["MT-CO1", "MT-CO2", "MT-CO3"],
    "interferon": ["IFITM1", "IFITM2", "IFITM3"],
    "hla": ["HLA-A", "HLA-B", "HLA-C"],
}


# Each FM stores its embedding under a different obsm key
OBSM_KEYS = {
    "state": ["X_state", "X_emb", "X_STATE"],
    "uce": ["X_uce", "X_uce_4layer"],
    "geneformer": ["X_geneformer", "X_Geneformer"],
    "transcriptformer": ["X_transcriptformer", "X_TranscriptFormer"],
    "scprint2": ["X_scprint2", "X_scprint", "X_scPRINT", "X_scPRINT2", "scprint_emb"],
    "stack": ["X_stack", "X_STACK"],
    "scgpt": ["X_scgpt", "X_scGPT"],
}


def find_embedding(adata, model_key: str) -> tuple[np.ndarray, str]:
    """Find an FM embedding by known obsm keys. Refuses the previous unsafe
    fallback (any 2D obsm with >1 column) -- that could silently grab UMAP
    coords or other auxiliary arrays. Raises a clear error if no known key is
    present. Mirrors the policy in primary_benchmark_option_c.py."""
    keys = OBSM_KEYS.get(model_key, [])
    for k in keys:
        if k in adata.obsm:
            X = np.asarray(adata.obsm[k])
            return X, k
    raise RuntimeError(
        f"No known embedding key for model_key={model_key!r}. "
        f"Tried {keys}; available obsm: {list(adata.obsm.keys())}. "
        f"If a new FM is added, register its obsm key(s) in OBSM_KEYS."
    )


def compute_module_score(adata_raw, gene_set_name: str, gene_list: list[str]) -> np.ndarray:
    """Run scanpy score_genes for one gene set, return a (n_obs,) score array."""
    valid = [g for g in gene_list if g in adata_raw.var_names]
    log.info("  %s: %d/%d markers in var_names", gene_set_name, len(valid), len(gene_list))
    if not valid:
        return np.zeros(adata_raw.n_obs, dtype=np.float32)
    sc.tl.score_genes(adata_raw, gene_list=valid, score_name=f"score_{gene_set_name}", use_raw=False)
    return np.asarray(adata_raw.obs[f"score_{gene_set_name}"].values, dtype=np.float32)


def per_fm_analysis(
    fm_h5ad_path: Path,
    raw_adata,
    module_scores: dict[str, np.ndarray],
    label_col: str,
    fatal_tokens: set[str],
    output_dir: Path,
) -> dict:
    model_key = fm_h5ad_path.stem.split("_")[0].lower()
    log.info("=== Model: %s (%s) ===", model_key, fm_h5ad_path.name)

    adata_emb = sc.read_h5ad(fm_h5ad_path)
    X, obsm_key = find_embedding(adata_emb, model_key)
    log.info("  obsm key=%s  shape=%s", obsm_key, X.shape)

    # Align to raw_adata if shapes differ. Try in order:
    #   1. obs.index intersection (TF/GF — preserved Cell_IDs)
    #   2. Cell_ID column intersection (scPRINT-2 — index is UUIDs but obs has Cell_ID)
    #   3. fallback: use FM file's own obs and recompute module scores from its .X
    if X.shape[0] != raw_adata.n_obs:
        log.warning("  embedding has %d cells, raw has %d — aligning",
                    X.shape[0], raw_adata.n_obs)

        common_idx = adata_emb.obs.index.intersection(raw_adata.obs.index)
        overlap_frac = len(common_idx) / max(X.shape[0], 1)
        log.info("  index overlap: %d cells (%.1f%% of FM cells)",
                 len(common_idx), 100 * overlap_frac)

        if overlap_frac >= 0.5:
            # Path 1: standard alignment via shared cell index
            emb_subset = adata_emb[common_idx].obsm[obsm_key]
            X = np.asarray(emb_subset)
            adata_subset = raw_adata[common_idx]
            scores_subset = {
                name: arr[adata_subset.obs.index.map(raw_adata.obs.index.get_loc).to_numpy()]
                for name, arr in module_scores.items()
            }
        elif "Cell_ID" in adata_emb.obs.columns and "Cell_ID" in raw_adata.obs.columns:
            # Path 2: align via Cell_ID column (scPRINT-2 case — index is UUIDs)
            log.info("  Falling back to Cell_ID column alignment")
            emb_cellids = adata_emb.obs["Cell_ID"].astype(str).values
            raw_cellids = raw_adata.obs["Cell_ID"].astype(str).values
            raw_id_to_pos = {cid: i for i, cid in enumerate(raw_cellids)}
            keep_emb_pos = [i for i, cid in enumerate(emb_cellids) if cid in raw_id_to_pos]
            keep_raw_pos = [raw_id_to_pos[emb_cellids[i]] for i in keep_emb_pos]
            log.info("  Cell_ID overlap: %d / %d FM cells matched", len(keep_emb_pos), len(emb_cellids))
            if len(keep_emb_pos) < 0.5 * X.shape[0]:
                raise RuntimeError(f"Cell_ID overlap too low ({len(keep_emb_pos)}/{X.shape[0]})")
            X = X[keep_emb_pos]
            adata_subset = raw_adata[keep_raw_pos]
            scores_subset = {name: arr[keep_raw_pos] for name, arr in module_scores.items()}
        else:
            # Path 3: no alignment possible — use FM file's own obs + recompute
            log.warning("  No alignment possible — using FM file's own obs + recomputing module scores from its .X")
            adata_subset = adata_emb
            scores_subset = {
                name: compute_module_score(adata_emb, name, gene_list)
                for name, gene_list in PLATELET_ACTIVATION_GENES.items()
            }
    else:
        adata_subset = raw_adata
        scores_subset = module_scores

    # Build a combined "any platelet activation" score (mean of all sets)
    combined = np.mean(np.column_stack(list(scores_subset.values())), axis=1)

    # Spearman per dim
    rows = []
    n_dim = X.shape[1]
    for d in range(n_dim):
        rho, _ = spearmanr(X[:, d], combined)
        rows.append({"model": model_key, "dim": d, "spearman_combined": rho,
                     "abs_spearman_combined": abs(rho)})
    df = pd.DataFrame(rows).sort_values("abs_spearman_combined", ascending=False)

    # Fatal vs survived test on the top axis
    sev = adata_subset.obs[label_col].astype(str).str.lower().values
    fatal_mask = np.isin(sev, list(fatal_tokens))
    surv_mask = np.isin(sev, ["mild", "moderate", "recovered", "healthy"])

    top_dim = int(df.iloc[0]["dim"])
    top_rho = float(df.iloc[0]["spearman_combined"])

    if fatal_mask.sum() > 0 and surv_mask.sum() > 0:
        fatal_vals = X[fatal_mask, top_dim]
        surv_vals = X[surv_mask, top_dim]
        u_stat, pval = mannwhitneyu(fatal_vals, surv_vals, alternative="two-sided")
        # AUROC of top-dim alone for fatal-vs-survived
        valid = fatal_mask | surv_mask
        y = fatal_mask[valid].astype(int)
        try:
            auroc_top_dim = roc_auc_score(y, X[valid, top_dim])
            auroc_module = roc_auc_score(y, combined[valid])
        except ValueError:
            auroc_top_dim = float("nan")
            auroc_module = float("nan")
    else:
        pval, auroc_top_dim, auroc_module = float("nan"), float("nan"), float("nan")

    df.to_csv(output_dir / f"{model_key}_dim_correlations.csv", index=False)
    log.info("  Top dim=%d, |rho|=%.3f, fatal-vs-survived MWU p=%.2e", top_dim, abs(top_rho), pval)
    log.info("  Top-dim AUROC=%.3f  Module-score AUROC=%.3f", auroc_top_dim, auroc_module)

    return {
        "model": model_key,
        "obsm_key": obsm_key,
        "n_dim": n_dim,
        "top_dim": top_dim,
        "top_dim_spearman": top_rho,
        "top_dim_fatal_vs_survived_pval": pval,
        "top_dim_auroc_fatal_vs_survived": auroc_top_dim,
        "module_score_auroc_fatal_vs_survived": auroc_module,
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--embedding-dir", required=True,
                   help="Dir with {state,uce,geneformer,transcriptformer,scprint2}_covid_canonical.h5ad")
    p.add_argument("--raw-input", required=True,
                   help="Raw AnnData with gene symbols in var_names")
    p.add_argument("--label-col", default="covid_severity")
    p.add_argument("--fatal-tokens", default="fatal", help="Comma-separated")
    p.add_argument("--output-dir", required=True)
    args = p.parse_args()

    emb_dir = Path(args.embedding_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fatal_tokens = {t.strip().lower() for t in args.fatal_tokens.split(",") if t.strip()}

    log.info("Reading raw AnnData: %s", args.raw_input)
    raw = sc.read_h5ad(args.raw_input)
    log.info("Raw shape: %d x %d", raw.n_obs, raw.n_vars)

    # Compute all module scores once on the raw expression
    log.info("Computing platelet activation module scores...")
    module_scores = {}
    for name, genes in PLATELET_ACTIVATION_GENES.items():
        module_scores[name] = compute_module_score(raw, name, genes)

    # Save module scores as a separate CSV (one row per cell)
    module_df = pd.DataFrame({k: v for k, v in module_scores.items()})
    module_df["combined"] = np.mean(module_df.values, axis=1)
    module_df.index = raw.obs.index
    module_df.to_csv(out_dir / "platelet_activation_module_scores.csv")
    log.info("Saved module scores to %s", out_dir / "platelet_activation_module_scores.csv")

    # Iterate over FM h5ads in embedding-dir
    summaries = []
    for fm_file in sorted(emb_dir.glob("*_covid_canonical.h5ad")):
        try:
            summaries.append(per_fm_analysis(fm_file, raw, module_scores,
                                             args.label_col, fatal_tokens, out_dir))
        except Exception as e:
            log.error("Failed on %s: %s", fm_file.name, e)
            import traceback
            traceback.print_exc()

    if summaries:
        summary_df = pd.DataFrame(summaries).sort_values(
            "top_dim_auroc_fatal_vs_survived", ascending=False
        )
        summary_df.to_csv(out_dir / "summary_per_fm.csv", index=False)
        log.info("\n=== SUMMARY ===\n%s", summary_df.to_string(index=False))

    return 0


if __name__ == "__main__":
    sys.exit(main())
