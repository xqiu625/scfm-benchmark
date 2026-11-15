#!/usr/bin/env python3
"""Option C 5-FM embedding ablation analysis (in-hospital mortality binary).

Four ablation types per FM:
  1. PCA dimension reduction: 64, 128, 256, 512
  2. Top-PC removal: remove top 1, 2, 5, 10 PCs before classification
  3. Normalization: raw, L2-norm, z-score (StandardScaler)
  4. Random dimension selection: 64, 128, 256 random dims (no PCA)

Each ablation is evaluated under 5-fold StratifiedGroupKFold (Sample_ID groups,
mortality binary task), with PCA fit inside each fold to avoid leakage.

Output:
  ablation_summary.csv     # model × ablation_type × param → AUC, BalAcc
  per_fold_<model>.csv     # fold-level numbers per ablation
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from joblib import Parallel, delayed
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler, normalize

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("ablation")

OBSM_KEYS = {
    "state": ["X_state", "X_emb", "X_STATE"],
    "uce": ["X_uce", "X_uce_4layer"],
    "geneformer": ["X_geneformer", "X_Geneformer"],
    "transcriptformer": ["X_transcriptformer", "X_TranscriptFormer"],
    "scprint2": ["X_scprint2", "X_scprint", "X_scPRINT", "X_scPRINT2", "scprint_emb"],
    "stack": ["X_stack", "X_STACK"],
    "scgpt": ["X_scgpt", "X_scGPT"],
}

PCA_DIMS = [64, 128, 256, 512]
REMOVE_TOP_PC = [1, 2, 5, 10]
NORMALIZATIONS = ["raw", "l2", "zscore"]
RANDOM_DIMS = [64, 128, 256]
N_FOLDS = 5
N_JOBS = 8


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


def fit_predict(X_tr, y_tr, X_te):
    clf = LogisticRegression(
        penalty="l2", C=1.0, solver="lbfgs",
        max_iter=2000, class_weight="balanced",
    )
    clf.fit(X_tr, y_tr)
    return clf.predict_proba(X_te)[:, 1], clf.predict(X_te)


def cv_score(X, y, groups, transform_fn, fold_seed=42):
    """5-fold StratifiedGroupKFold; transform_fn applied INSIDE each fold (no leakage)."""
    skf = StratifiedGroupKFold(n_splits=N_FOLDS, shuffle=True, random_state=fold_seed)
    aucs, baccs = [], []
    for fold, (tr, te) in enumerate(skf.split(X, y, groups)):
        X_tr_raw, X_te_raw = X[tr], X[te]
        X_tr, X_te = transform_fn(X_tr_raw, X_te_raw)
        proba, pred = fit_predict(X_tr, y[tr], X_te)
        aucs.append(roc_auc_score(y[te], proba))
        baccs.append(balanced_accuracy_score(y[te], pred))
    return float(np.mean(aucs)), float(np.std(aucs)), float(np.mean(baccs))


def transform_pca(d):
    def _f(X_tr, X_te):
        n_comp = min(d, X_tr.shape[1], X_tr.shape[0] - 1)
        pca = PCA(n_components=n_comp, random_state=42)
        return pca.fit_transform(X_tr), pca.transform(X_te)
    return _f


def transform_remove_top_pc(k):
    def _f(X_tr, X_te):
        scaler = StandardScaler().fit(X_tr)
        Xt_tr = scaler.transform(X_tr)
        Xt_te = scaler.transform(X_te)
        n_comp = min(k + 50, X_tr.shape[1], X_tr.shape[0] - 1)
        pca = PCA(n_components=n_comp, random_state=42)
        Z_tr = pca.fit_transform(Xt_tr)
        Z_te = pca.transform(Xt_te)
        # Drop top-k components (column indices 0..k-1)
        Z_tr_drop = Z_tr[:, k:]
        Z_te_drop = Z_te[:, k:]
        return Z_tr_drop, Z_te_drop
    return _f


def transform_normalize(kind):
    def _f(X_tr, X_te):
        if kind == "raw":
            return X_tr, X_te
        elif kind == "l2":
            return normalize(X_tr, norm="l2"), normalize(X_te, norm="l2")
        elif kind == "zscore":
            scaler = StandardScaler().fit(X_tr)
            return scaler.transform(X_tr), scaler.transform(X_te)
        raise ValueError(kind)
    return _f


def transform_random_dims(d, seed):
    def _f(X_tr, X_te):
        rng = np.random.RandomState(seed)
        idx = rng.choice(X_tr.shape[1], size=min(d, X_tr.shape[1]), replace=False)
        return X_tr[:, idx], X_te[:, idx]
    return _f


def per_fm(fm_h5ad: Path, output_dir: Path) -> list[dict]:
    model_key = fm_h5ad.stem.split("_")[0].lower()
    log.info("=== %s (%s) ===", model_key, fm_h5ad.name)

    adata = sc.read_h5ad(fm_h5ad)
    X, obsm_key = find_embedding(adata, model_key)
    log.info("  obsm=%s shape=%s", obsm_key, X.shape)

    if "covid_severity" not in adata.obs.columns:
        log.warning("  no covid_severity; skipping")
        return []
    if "Sample_ID" not in adata.obs.columns:
        log.warning("  no Sample_ID; skipping")
        return []

    y_full, valid = build_mortality(adata.obs["covid_severity"].astype(str).values)
    groups = adata.obs["Sample_ID"].astype(str).values
    X_v, y_v, g_v = X[valid], y_full[valid], groups[valid]
    log.info("  CV input: n=%d, fatal=%d, survived=%d, groups=%d",
             len(y_v), int(y_v.sum()), int((y_v == 0).sum()), len(np.unique(g_v)))

    rows = []

    # 1. Baseline (no transform)
    log.info("  [baseline] no transform")
    auc, std, bacc = cv_score(X_v, y_v, g_v, lambda a, b: (a, b))
    rows.append(dict(model=model_key, ablation="baseline", param="none",
                     auc_mean=auc, auc_std=std, balanced_accuracy=bacc))

    # 2. PCA dimensionality
    for d in PCA_DIMS:
        if d >= X_v.shape[1]:
            continue
        log.info("  [pca_dims] d=%d", d)
        auc, std, bacc = cv_score(X_v, y_v, g_v, transform_pca(d))
        rows.append(dict(model=model_key, ablation="pca_dims", param=str(d),
                         auc_mean=auc, auc_std=std, balanced_accuracy=bacc))

    # 3. Remove top PCs
    for k in REMOVE_TOP_PC:
        log.info("  [remove_top_pc] k=%d", k)
        auc, std, bacc = cv_score(X_v, y_v, g_v, transform_remove_top_pc(k))
        rows.append(dict(model=model_key, ablation="remove_top_pc", param=str(k),
                         auc_mean=auc, auc_std=std, balanced_accuracy=bacc))

    # 4. Normalization
    for kind in NORMALIZATIONS:
        log.info("  [normalize] kind=%s", kind)
        auc, std, bacc = cv_score(X_v, y_v, g_v, transform_normalize(kind))
        rows.append(dict(model=model_key, ablation="normalize", param=kind,
                         auc_mean=auc, auc_std=std, balanced_accuracy=bacc))

    # 5. Random dim selection
    for d in RANDOM_DIMS:
        if d >= X_v.shape[1]:
            continue
        log.info("  [random_dims] d=%d", d)
        auc, std, bacc = cv_score(X_v, y_v, g_v, transform_random_dims(d, seed=42))
        rows.append(dict(model=model_key, ablation="random_dims", param=str(d),
                         auc_mean=auc, auc_std=std, balanced_accuracy=bacc))

    df = pd.DataFrame(rows)
    df.to_csv(output_dir / f"per_fm_{model_key}.csv", index=False)
    log.info("  Wrote per_fm_%s.csv (%d rows)", model_key, len(df))
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--embedding-dir", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--models", nargs="+", default=None,
                   help="Subset of FMs to run, by model_key (e.g. --models scgpt). "
                        "Default: all *_covid_canonical.h5ad in --embedding-dir. "
                        "On a subset run, new rows are merged into an existing "
                        "ablation_summary.csv (recomputed models replaced) instead "
                        "of overwriting it with subset-only rows.")
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    fm_files = sorted(args.embedding_dir.glob("*_covid_canonical.h5ad"))
    if args.models:
        wanted = {m.lower() for m in args.models}
        fm_files = [f for f in fm_files if f.stem.split("_")[0].lower() in wanted]
    if not fm_files:
        log.error("No *_covid_canonical.h5ad in %s", args.embedding_dir)
        sys.exit(1)
    log.info("Found %d FM files: %s", len(fm_files), [f.name for f in fm_files])

    all_rows = []
    for fm in fm_files:
        try:
            all_rows.extend(per_fm(fm, args.output_dir))
        except Exception as e:
            log.exception("Failed %s: %s", fm.name, e)

    if not all_rows:
        log.error("No FMs produced output")
        sys.exit(1)

    summary = pd.DataFrame(all_rows)
    out_csv = args.output_dir / "ablation_summary.csv"
    if args.models and out_csv.exists():
        existing = pd.read_csv(out_csv)
        recomputed = set(summary["model"].unique())
        kept = existing[~existing["model"].isin(recomputed)]
        summary = pd.concat([kept, summary], ignore_index=True)
        log.info("Merged into existing ablation_summary.csv: replaced models %s, "
                 "kept %d other models", sorted(recomputed), kept["model"].nunique())
    elif args.models:
        log.warning("No existing ablation_summary.csv; summary will contain only "
                    "the models computed in this run")
    summary.to_csv(out_csv, index=False)
    log.info("Wrote %s (%d rows)", out_csv, len(summary))
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
