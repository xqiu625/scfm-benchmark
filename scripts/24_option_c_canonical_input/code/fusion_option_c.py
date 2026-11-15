#!/usr/bin/env python3
"""Multi-FM late-fusion analysis on Option C embeddings (mortality binary).

Strategies evaluated (out-of-fold logistic regression per FM, then combine):
  1. Late average      — mean of per-FM predicted probabilities
  2. Late rank         — mean of rank-transformed per-FM probabilities
  3. Late max          — max probability across FMs
  4. Late stacking     — LR meta-classifier on per-FM probability stack
                         (5-fold inner CV to avoid leakage)
  5. Late weighted     — weighted average using per-FM CV AUC as weights
  6. Early concat+PCA  — concatenate all 5 embeddings, PCA-256, LR
  7. Single-FM bests   — each FM alone (for reference)

All strategies evaluated on the SAME Cell_ID-aligned cell set
(intersection across all FMs) under 5-fold StratifiedGroupKFold(Sample_ID).

Output:
  fusion_summary.csv         strategy, AUC, BalAcc, n_cells
  per_fm_oof.csv             per-FM OOF probabilities (for downstream analysis)
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import rankdata
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("fusion")

OBSM_KEYS = {
    "state": ["X_state", "X_emb", "X_STATE"],
    "uce": ["X_uce", "X_uce_4layer"],
    "geneformer": ["X_geneformer", "X_Geneformer"],
    "transcriptformer": ["X_transcriptformer", "X_TranscriptFormer"],
    "scprint2": ["X_scprint2", "X_scprint", "X_scPRINT", "X_scPRINT2", "scprint_emb"],
    "stack": ["X_stack", "X_STACK"],
    "scgpt": ["X_scgpt", "X_scGPT"],
}
N_FOLDS = 5


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


def lr_oof(X, y, groups):
    """5-fold OOF predictions; standardize inside fold."""
    skf = StratifiedGroupKFold(n_splits=N_FOLDS, shuffle=True, random_state=42)
    oof = np.full(len(y), np.nan)
    for tr, te in skf.split(X, y, groups):
        scaler = StandardScaler().fit(X[tr])
        clf = LogisticRegression(penalty="l2", C=1.0, solver="lbfgs",
                                 max_iter=2000, class_weight="balanced", n_jobs=4)
        clf.fit(scaler.transform(X[tr]), y[tr])
        oof[te] = clf.predict_proba(scaler.transform(X[te]))[:, 1]
    return oof


def stacking_nested(per_fm_X: dict, fm_keys, y, groups):
    """Late stacking with TRUE nested cross-validation.

    `stacking_oof` below trains the meta-learner on base-model out-of-fold
    probabilities that were themselves produced using every cell in the dataset.
    A cell sitting in the meta-learner's test fold therefore contributed labels to
    the base models that generated its own meta-features, which leaks and inflates
    the stacking result.

    Here each outer fold is self-contained: base models are refit on the outer
    training block only, meta-features for the outer test block come from those
    refits, and the meta-learner is trained on an inner-CV OOF matrix built
    strictly inside the outer training block.
    """
    outer = StratifiedGroupKFold(n_splits=N_FOLDS, shuffle=True, random_state=42)
    ref = per_fm_X[fm_keys[0]]
    oof = np.full(len(y), np.nan)

    def _fit(Xtr, ytr):
        sc = StandardScaler().fit(Xtr)
        clf = LogisticRegression(penalty="l2", C=1.0, solver="lbfgs",
                                 max_iter=2000, class_weight="balanced", n_jobs=4)
        clf.fit(sc.transform(Xtr), ytr)
        return sc, clf

    for fold, (tr, te) in enumerate(outer.split(ref, y, groups), 1):
        y_tr, g_tr = y[tr], groups[tr]
        meta_tr = np.zeros((len(tr), len(fm_keys)))
        meta_te = np.zeros((len(te), len(fm_keys)))
        inner = StratifiedGroupKFold(n_splits=N_FOLDS, shuffle=True, random_state=43)
        for j, k in enumerate(fm_keys):
            Xk = per_fm_X[k]
            Xk_tr = Xk[tr]
            for itr, ite in inner.split(Xk_tr, y_tr, g_tr):
                sc, clf = _fit(Xk_tr[itr], y_tr[itr])
                meta_tr[ite, j] = clf.predict_proba(sc.transform(Xk_tr[ite]))[:, 1]
            sc, clf = _fit(Xk_tr, y_tr)
            meta_te[:, j] = clf.predict_proba(sc.transform(Xk[te]))[:, 1]
        meta = LogisticRegression(penalty="l2", C=1.0, solver="lbfgs",
                                  max_iter=2000, class_weight="balanced")
        meta.fit(meta_tr, y_tr)
        oof[te] = meta.predict_proba(meta_te)[:, 1]
        log.info("  [nested] outer fold %d/%d done", fold, N_FOLDS)
    return oof


def stacking_oof(per_fm_oof: np.ndarray, y, groups):
    """Train an LR meta-classifier on per-FM OOF probabilities; 5-fold inner CV.
    per_fm_oof: (n_cells, n_fms) — out-of-fold probabilities from each FM.
    Returns OOF predictions of the meta-classifier."""
    skf = StratifiedGroupKFold(n_splits=N_FOLDS, shuffle=True, random_state=42)
    oof = np.full(len(y), np.nan)
    for tr, te in skf.split(per_fm_oof, y, groups):
        meta = LogisticRegression(penalty="l2", C=1.0, solver="lbfgs",
                                  max_iter=2000, class_weight="balanced")
        meta.fit(per_fm_oof[tr], y[tr])
        oof[te] = meta.predict_proba(per_fm_oof[te])[:, 1]
    return oof


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--embedding-dir", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--models", nargs="+", default=None,
                   help="Subset of FMs to include in the fusion, by model_key "
                        "(e.g. --models state uce geneformer scprint2 scgpt for the "
                        "Option C primary 5-FM set). Default: all "
                        "*_covid_canonical.h5ad in --embedding-dir.")
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    fm_files = sorted(args.embedding_dir.glob("*_covid_canonical.h5ad"))
    if args.models:
        wanted = {m.lower() for m in args.models}
        fm_files = [f for f in fm_files if f.stem.split("_")[0].lower() in wanted]
    if not fm_files:
        log.error("No *_covid_canonical.h5ad in %s", args.embedding_dir); sys.exit(1)

    # --- Load each FM's embedding + Cell_ID + severity ---
    fm_data = {}  # mk -> {X, cellids, y, groups}
    for fm in fm_files:
        mk = fm.stem.split("_")[0].lower()
        log.info("Loading %s...", mk)
        a = sc.read_h5ad(fm)
        X, key = find_embedding(a, mk)
        if "Cell_ID" not in a.obs.columns:
            log.warning("  %s missing Cell_ID; skipping", mk); continue
        if "covid_severity" not in a.obs.columns or "Sample_ID" not in a.obs.columns:
            log.warning("  %s missing severity/Sample_ID; skipping", mk); continue
        cellids = a.obs["Cell_ID"].astype(str).values
        y_full, valid = build_mortality(a.obs["covid_severity"].astype(str).values)
        groups = a.obs["Sample_ID"].astype(str).values
        fm_data[mk] = dict(
            X=X[valid], cellids=cellids[valid], y=y_full[valid], groups=groups[valid],
            n_dim=X.shape[1], obsm_key=key,
        )
        log.info("  %s: %d valid cells, %d dim", mk, valid.sum(), X.shape[1])
        del a

    if len(fm_data) < 2:
        log.error("Need at least 2 FMs"); sys.exit(1)

    # --- Cell_ID intersection across all FMs ---
    common = set(next(iter(fm_data.values()))["cellids"])
    for d in fm_data.values():
        common &= set(d["cellids"])
    common_list = sorted(common)
    common_set = set(common_list)
    log.info("Common Cell_IDs across %d FMs: %d", len(fm_data), len(common_list))

    # Re-index each FM to the common cell set
    aligned = {}
    for mk, d in fm_data.items():
        cid_to_pos = {c: i for i, c in enumerate(d["cellids"])}
        keep = [cid_to_pos[c] for c in common_list]
        aligned[mk] = dict(
            X=d["X"][keep], y=d["y"][keep], groups=d["groups"][keep], n_dim=d["n_dim"],
        )
    # Reference y/groups (all FMs should agree on these)
    y_ref = aligned[next(iter(aligned))]["y"]
    g_ref = aligned[next(iter(aligned))]["groups"]
    for mk, d in aligned.items():
        if not np.array_equal(d["y"], y_ref):
            log.warning("%s: y mismatch on common cells", mk)

    # --- Per-FM OOF predictions ---
    per_fm_oof = {}
    per_fm_auc = {}
    for mk, d in aligned.items():
        log.info("[OOF] %s (n=%d, dim=%d)", mk, len(d["y"]), d["n_dim"])
        oof = lr_oof(d["X"], d["y"], d["groups"])
        per_fm_oof[mk] = oof
        per_fm_auc[mk] = roc_auc_score(y_ref, oof)
        log.info("  AUC: %.4f", per_fm_auc[mk])

    # Save per-FM OOF for downstream analysis
    df_oof = pd.DataFrame({"y": y_ref, **per_fm_oof})
    df_oof.to_csv(args.output_dir / "per_fm_oof.csv", index=False)

    fm_keys = sorted(per_fm_oof.keys())
    proba_stack = np.column_stack([per_fm_oof[k] for k in fm_keys])

    # --- Fusion strategies ---
    rows = []

    # 1. Late average
    avg = proba_stack.mean(axis=1)
    rows.append(dict(strategy="late_average", auc=roc_auc_score(y_ref, avg),
                     balanced_accuracy=balanced_accuracy_score(y_ref, (avg >= 0.5).astype(int)),
                     n_cells=len(y_ref), n_fms=len(fm_keys)))

    # 2. Late rank
    ranks = np.column_stack([rankdata(per_fm_oof[k], method="average") for k in fm_keys])
    rank_avg = ranks.mean(axis=1)
    rows.append(dict(strategy="late_rank", auc=roc_auc_score(y_ref, rank_avg),
                     balanced_accuracy=balanced_accuracy_score(y_ref, (rank_avg >= np.median(rank_avg)).astype(int)),
                     n_cells=len(y_ref), n_fms=len(fm_keys)))

    # 3. Late max
    mx = proba_stack.max(axis=1)
    rows.append(dict(strategy="late_max", auc=roc_auc_score(y_ref, mx),
                     balanced_accuracy=balanced_accuracy_score(y_ref, (mx >= 0.5).astype(int)),
                     n_cells=len(y_ref), n_fms=len(fm_keys)))

    # 4. Late stacking (meta-LR with inner 5-fold CV)
    log.info("[Stacking] training meta-LR on per-FM probabilities")
    stacked = stacking_oof(proba_stack, y_ref, g_ref)
    log.info("[stacking] nested-CV variant (leakage-free); this is the reportable number")
    stacked_nested = stacking_nested({k: aligned[k]["X"] for k in fm_keys},
                                     fm_keys, y_ref, g_ref)
    rows.append(dict(strategy="late_stacking_lr_nested",
                     auc=roc_auc_score(y_ref, stacked_nested),
                     balanced_accuracy=balanced_accuracy_score(
                         y_ref, (stacked_nested > 0.5).astype(int)),
                     n_cells=len(y_ref), n_fms=len(fm_keys)))
    rows.append(dict(strategy="late_stacking_lr", auc=roc_auc_score(y_ref, stacked),
                     balanced_accuracy=balanced_accuracy_score(y_ref, (stacked >= 0.5).astype(int)),
                     n_cells=len(y_ref), n_fms=len(fm_keys)))

    # 5. Late weighted (per-FM AUC as weight)
    weights = np.array([per_fm_auc[k] for k in fm_keys])
    weights = weights / weights.sum()
    weighted = (proba_stack * weights[None, :]).sum(axis=1)
    rows.append(dict(strategy="late_weighted_by_auc", auc=roc_auc_score(y_ref, weighted),
                     balanced_accuracy=balanced_accuracy_score(y_ref, (weighted >= 0.5).astype(int)),
                     n_cells=len(y_ref), n_fms=len(fm_keys)))

    # 6. Early concat + PCA-256 + LR (5-fold OOF)
    log.info("[Early] concatenating + PCA-256")
    X_concat = np.column_stack([aligned[k]["X"] for k in fm_keys])
    log.info("  Concat dim: %d", X_concat.shape[1])
    skf = StratifiedGroupKFold(n_splits=N_FOLDS, shuffle=True, random_state=42)
    early_oof = np.full(len(y_ref), np.nan)
    for tr, te in skf.split(X_concat, y_ref, g_ref):
        scaler = StandardScaler().fit(X_concat[tr])
        Xs_tr = scaler.transform(X_concat[tr])
        Xs_te = scaler.transform(X_concat[te])
        n_comp = min(256, Xs_tr.shape[1], Xs_tr.shape[0] - 1)
        pca = PCA(n_components=n_comp, random_state=42)
        Z_tr = pca.fit_transform(Xs_tr)
        Z_te = pca.transform(Xs_te)
        clf = LogisticRegression(penalty="l2", C=1.0, solver="lbfgs",
                                 max_iter=2000, class_weight="balanced", n_jobs=4)
        clf.fit(Z_tr, y_ref[tr])
        early_oof[te] = clf.predict_proba(Z_te)[:, 1]
    rows.append(dict(strategy="early_concat_pca256", auc=roc_auc_score(y_ref, early_oof),
                     balanced_accuracy=balanced_accuracy_score(y_ref, (early_oof >= 0.5).astype(int)),
                     n_cells=len(y_ref), n_fms=len(fm_keys)))

    # 7. Single-FM references (for the table)
    for mk, oof in per_fm_oof.items():
        rows.append(dict(strategy=f"single_{mk}", auc=per_fm_auc[mk],
                         balanced_accuracy=balanced_accuracy_score(y_ref, (oof >= 0.5).astype(int)),
                         n_cells=len(y_ref), n_fms=1))

    df = pd.DataFrame(rows).sort_values("auc", ascending=False)
    out_csv = args.output_dir / "fusion_summary.csv"
    df.to_csv(out_csv, index=False)
    log.info("Wrote %s", out_csv)
    print("\n=== FUSION SUMMARY ===")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
