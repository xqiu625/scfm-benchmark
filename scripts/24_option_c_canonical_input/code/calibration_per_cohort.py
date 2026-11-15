#!/usr/bin/env python3
"""Per-FM calibration + per-cohort performance for the in-hospital mortality task.

Addresses CZI Virtual Cells Workshop benchmarking recommendations
(Fahsbender et al. 2025, arXiv 2507.10502v2):
  - Multi-faceted evaluation beyond AUC (Brier, ECE, calibration plots)
  - Per-cohort performance breakdown (prevents over-claiming generalization)
  - Reproducible patient-level splits

For each Option-C FM h5ad:
  1. Build mortality binary labels (fatal=1; mild/moderate/recovered/severe=0;
     healthy/unknown excluded)
  2. 5-fold StratifiedGroupKFold(Sample_ID) — same as primary benchmark
  3. Logistic regression with class_weight='balanced'
  4. Out-of-fold predicted probabilities
  5. Compute: AUC, Brier score, ECE (10 bins)
  6. Calibration plot per FM (saved as PNG)
  7. Per-cohort AUC breakdown (cohort = `Data_NO` or `study_id` column)

Output:
  calibration_summary.csv         # per-FM AUC, Brier, ECE
  per_cohort_auc.csv              # FM × cohort AUC
  calibration_plot_<fm>.png       # per-FM reliability diagram

Usage:
    python calibration_per_cohort.py \
        --embedding-dir /bigdata/.../02_EMBEDDINGS/option_c \
        --label-col covid_severity \
        --group-col Sample_ID \
        --cohort-col Data_NO \
        --output-dir /bigdata/.../05_RESULTS/option_c_calibration
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.calibration import calibration_curve
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("calibration")


OBSM_KEYS = {
    "state": ["X_state", "X_emb", "X_STATE"],
    "uce": ["X_uce", "X_uce_4layer"],
    "geneformer": ["X_geneformer", "X_Geneformer"],
    "transcriptformer": ["X_transcriptformer", "X_TranscriptFormer"],
    "scprint2": ["X_scprint2", "X_scprint", "X_scPRINT", "X_scPRINT2", "scprint_emb"],
    "stack": ["X_stack", "X_STACK"],
    "scgpt": ["X_scgpt", "X_scGPT"],
}


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


def expected_calibration_error(y_true, y_prob, n_bins: int = 10) -> float:
    """ECE with equal-width bins."""
    bin_edges = np.linspace(0, 1, n_bins + 1)
    ece = 0.0
    n = len(y_true)
    for i in range(n_bins):
        mask = (y_prob >= bin_edges[i]) & (y_prob < bin_edges[i + 1])
        if i == n_bins - 1:
            mask |= y_prob == 1.0
        if mask.sum() == 0:
            continue
        bin_acc = y_true[mask].mean()
        bin_conf = y_prob[mask].mean()
        ece += (mask.sum() / n) * abs(bin_acc - bin_conf)
    return float(ece)


def build_mortality_labels(severity: np.ndarray):
    """fatal=1; mild/moderate/recovered/severe=0; healthy/unknown excluded."""
    sev = np.asarray([str(s).lower() for s in severity])
    fatal_mask = sev == "fatal"
    survived_tokens = {"mild", "moderate", "recovered", "severe"}
    survived_mask = np.isin(sev, list(survived_tokens))
    valid = fatal_mask | survived_mask
    y = np.zeros(len(sev), dtype=np.int8)
    y[fatal_mask] = 1
    return y, valid


def per_fm_calibration(
    fm_h5ad: Path,
    label_col: str,
    group_col: str,
    cohort_col: str,
    output_dir: Path,
):
    model_key = fm_h5ad.stem.split("_")[0].lower()
    log.info("=== %s (%s) ===", model_key, fm_h5ad.name)

    adata = sc.read_h5ad(fm_h5ad)
    X, obsm_key = find_embedding(adata, model_key)
    log.info("  obsm=%s shape=%s", obsm_key, X.shape)

    if label_col not in adata.obs.columns:
        raise ValueError(f"{label_col} not in obs")
    if group_col not in adata.obs.columns:
        raise ValueError(f"{group_col} not in obs")

    severity = adata.obs[label_col].astype(str).values
    groups = adata.obs[group_col].astype(str).values
    cohorts = adata.obs[cohort_col].astype(str).values if cohort_col in adata.obs.columns else None

    y, valid = build_mortality_labels(severity)
    X_v = X[valid]
    y_v = y[valid]
    g_v = groups[valid]
    c_v = cohorts[valid] if cohorts is not None else None
    log.info("  Valid cells: %d (fatal=%d, survived=%d)", valid.sum(), y_v.sum(), (1 - y_v).sum())

    skf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)
    oof_proba = np.full(len(y_v), np.nan, dtype=np.float64)

    for fold, (tr, te) in enumerate(skf.split(X_v, y_v, g_v)):
        scaler = StandardScaler()
        Xtr = scaler.fit_transform(X_v[tr])
        Xte = scaler.transform(X_v[te])
        clf = LogisticRegression(max_iter=2000, n_jobs=-1, class_weight="balanced",
                                 random_state=42)
        clf.fit(Xtr, y_v[tr])
        oof_proba[te] = clf.predict_proba(Xte)[:, 1]
        log.info("  Fold %d/5: train=%d, test=%d, AUC=%.3f",
                 fold + 1, len(tr), len(te), roc_auc_score(y_v[te], oof_proba[te]))

    # Aggregate metrics
    auc = float(roc_auc_score(y_v, oof_proba))
    brier = float(brier_score_loss(y_v, oof_proba))
    ece = expected_calibration_error(y_v, oof_proba, n_bins=10)
    log.info("  Overall: AUC=%.4f Brier=%.4f ECE=%.4f", auc, brier, ece)

    # Calibration plot
    fig, ax = plt.subplots(figsize=(5, 5))
    frac_pos, mean_pred = calibration_curve(y_v, oof_proba, n_bins=10, strategy="uniform")
    ax.plot([0, 1], [0, 1], "k--", alpha=0.5, label="Perfect")
    ax.plot(mean_pred, frac_pos, "o-", label=f"{model_key} (Brier={brier:.3f}, ECE={ece:.3f})")
    ax.set_xlabel("Mean predicted probability (fatal)")
    ax.set_ylabel("Fraction observed fatal")
    ax.set_title(f"Calibration: {model_key} (mortality binary)")
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    plot_path = output_dir / f"calibration_plot_{model_key}.png"
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)
    log.info("  Saved: %s", plot_path)

    # Per-cohort AUC
    cohort_aucs = {}
    if c_v is not None:
        for cohort in np.unique(c_v):
            mask = c_v == cohort
            if mask.sum() < 30 or len(np.unique(y_v[mask])) < 2:
                continue
            try:
                cohort_aucs[str(cohort)] = float(roc_auc_score(y_v[mask], oof_proba[mask]))
            except ValueError:
                pass
        log.info("  Per-cohort AUC: %s", cohort_aucs)

    return {
        "model": model_key,
        "obsm_key": obsm_key,
        "n_valid_cells": int(valid.sum()),
        "n_fatal": int(y_v.sum()),
        "auc": auc,
        "brier_score": brier,
        "ece_10bins": ece,
        "cohort_aucs": cohort_aucs,
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--embedding-dir", required=True)
    p.add_argument("--label-col", default="covid_severity")
    p.add_argument("--group-col", default="Sample_ID")
    p.add_argument("--cohort-col", default="Data_NO")
    p.add_argument("--output-dir", required=True)
    args = p.parse_args()

    emb_dir = Path(args.embedding_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    summaries = []
    for fm_file in sorted(emb_dir.glob("*_covid_canonical.h5ad")):
        try:
            summaries.append(per_fm_calibration(
                fm_file, args.label_col, args.group_col, args.cohort_col, out_dir
            ))
        except Exception as e:
            log.error("Failed %s: %s", fm_file.name, e)
            import traceback; traceback.print_exc()

    if summaries:
        # Top-level summary CSV
        summary_df = pd.DataFrame([
            {k: v for k, v in s.items() if k != "cohort_aucs"} for s in summaries
        ]).sort_values("auc", ascending=False)
        summary_df.to_csv(out_dir / "calibration_summary.csv", index=False)
        log.info("\n=== SUMMARY ===\n%s", summary_df.to_string(index=False))

        # Per-cohort CSV
        rows = []
        for s in summaries:
            for cohort, auc in s["cohort_aucs"].items():
                rows.append({"model": s["model"], "cohort": cohort, "auc": auc})
        if rows:
            cohort_df = pd.DataFrame(rows).pivot(index="cohort", columns="model", values="auc")
            cohort_df.to_csv(out_dir / "per_cohort_auc.csv")
            log.info("Per-cohort AUC matrix saved")

    return 0


if __name__ == "__main__":
    sys.exit(main())
