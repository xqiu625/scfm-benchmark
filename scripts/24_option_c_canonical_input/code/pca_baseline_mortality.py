#!/usr/bin/env python3
"""PCA baseline on raw gene expression for in-hospital mortality binary task.

Fits PCA dim reduction (50, 100, 200, 500 components) on the canonical raw
input matrix, then logistic regression under 5-fold patient-level CV.
Compares to the 5 Option C foundation models on the same task.
"""
from __future__ import annotations

import argparse, logging, sys
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("pca_baseline")

PCA_DIMS = [50, 100, 200, 500]
N_FOLDS = 5


def build_mortality(severity):
    sev = np.asarray([str(s).lower() for s in severity])
    fatal = sev == "fatal"
    survived = np.isin(sev, ["mild", "moderate", "recovered", "severe"])
    valid = fatal | survived
    y = np.zeros(len(sev), dtype=np.int8)
    y[fatal] = 1
    return y, valid


def cv_lr_with_pca(X, y, groups, n_components=None):
    skf = StratifiedGroupKFold(n_splits=N_FOLDS, shuffle=True, random_state=42)
    aucs, baccs = [], []
    for tr, te in skf.split(X, y, groups):
        scaler = StandardScaler(with_mean=False).fit(X[tr])  # with_mean=False for sparse
        Xs_tr = scaler.transform(X[tr])
        Xs_te = scaler.transform(X[te])
        if n_components is not None:
            pca = PCA(n_components=n_components, random_state=42)
            Z_tr = pca.fit_transform(Xs_tr.toarray() if hasattr(Xs_tr, "toarray") else Xs_tr)
            Z_te = pca.transform(Xs_te.toarray() if hasattr(Xs_te, "toarray") else Xs_te)
        else:
            Z_tr, Z_te = Xs_tr, Xs_te
        clf = LogisticRegression(penalty="l2", C=1.0, solver="lbfgs", max_iter=2000,
                                 class_weight="balanced", n_jobs=4)
        clf.fit(Z_tr, y[tr])
        p = clf.predict_proba(Z_te)[:, 1]
        aucs.append(roc_auc_score(y[te], p))
        baccs.append(balanced_accuracy_score(y[te], (p >= 0.5).astype(int)))
    return float(np.mean(aucs)), float(np.std(aucs)), float(np.mean(baccs))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    log.info("Loading: %s", args.input)
    a = sc.read_h5ad(args.input)
    log.info("  %d cells × %d genes", a.n_obs, a.n_vars)

    y_full, valid = build_mortality(a.obs["covid_severity"].astype(str).values)
    groups = a.obs["Sample_ID"].astype(str).values
    X_v, y_v, g_v = a.X[valid], y_full[valid], groups[valid]
    log.info("  Valid: n=%d, fatal=%d, survived=%d", len(y_v), int(y_v.sum()), int((y_v == 0).sum()))

    rows = []
    for d in PCA_DIMS:
        log.info("[PCA-%d]", d)
        auc, std, bacc = cv_lr_with_pca(X_v, y_v, g_v, n_components=d)
        rows.append(dict(method=f"PCA-{d} + LR", n_components=d, auc_mean=auc, auc_std=std, balanced_accuracy=bacc))
        log.info("  AUC=%.4f (+/- %.4f), BalAcc=%.4f", auc, std, bacc)

    df = pd.DataFrame(rows)
    df.to_csv(args.output, index=False)
    log.info("Wrote %s", args.output)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
