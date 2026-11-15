#!/usr/bin/env python3
"""Validate the geometric mechanism of cross-disease mortality transfer.

For each Option-C foundation-model embedding, train a within-disease binary
mortality classifier (logistic regression, fatal vs. survived) twice:
  - On COVID cells only -> w_covid_fatal (LR coefficient vector)
  - On sepsis cells only -> w_sepsis_fatal

Then compute cos(w_covid_fatal, w_sepsis_fatal). The hypothesis (Discussion
section "Cross-Disease Generalization Reveals Embedding-Level Differences"):

  - Bidirectional generalizers (STATE, STACK): cos > 0.3
  - Forward-only generalizers (UCE, TF, scPRINT-2): cos around 0
  - Anti-correlated transfer (legacy scGPT, if included): cos < 0

This is the geometric "smoking gun" for why some FMs transfer mortality
across diseases and others do not. Confirms that cross-disease asymmetry
is fixed by pre-training, not by the downstream classifier.

Inputs:
  - Each FM h5ad must contain the embedding in obsm and metadata
    (covid_severity, Disease, Sample_ID) in obs
  - Sepsis cells: Disease in {sepsis, Sepsis} (we accept both); fatal vs severe
  - COVID cells: Disease in {COVID, COVID-19, covid_19, COVID19}; fatal vs
    {mild, moderate, recovered, severe}

Output:
  cosine_fatal_direction.csv:
    model, n_covid, n_sepsis, cos_similarity, dot_product, w_covid_norm, w_sepsis_norm
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("cosine_validation")


OBSM_KEYS = {
    "state": ["X_state", "X_emb", "X_STATE"],
    "uce": ["X_uce", "X_uce_4layer"],
    "geneformer": ["X_geneformer", "X_Geneformer"],
    "transcriptformer": ["X_transcriptformer", "X_TranscriptFormer"],
    "scprint2": ["X_scprint2", "X_scprint", "X_scPRINT", "X_scPRINT2", "scprint_emb"],
    "stack": ["X_stack", "X_STACK"],
    "scgpt": ["X_scgpt", "X_scGPT"],
}

COVID_TOKENS = {"covid", "covid-19", "covid_19", "covid19", "sars-cov-2"}
SEPSIS_TOKENS = {"sepsis"}


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


def disease_mask(disease_col: np.ndarray, tokens: set[str]) -> np.ndarray:
    norm = np.array([str(d).strip().lower() for d in disease_col])
    return np.isin(norm, list(tokens))


def fit_lr_direction(X: np.ndarray, y: np.ndarray, scaler: StandardScaler | None = None) -> np.ndarray:
    """Fit LR(class_weight=balanced) and return the (D,) coefficient vector.

    If scaler is None, fit a fresh per-disease StandardScaler (legacy behavior;
    coefficients live in per-disease standardized space, NOT comparable across
    diseases). If scaler is provided, transform X with it first (coefficients
    live in the SHARED standardized space — the correct setting for cosine
    comparison across diseases)."""
    if scaler is None:
        scaler = StandardScaler(with_mean=True, with_std=True)
        Xs = scaler.fit_transform(X)
    else:
        Xs = scaler.transform(X)
    clf = LogisticRegression(
        penalty="l2",
        C=1.0,
        solver="lbfgs",
        max_iter=2000,
        class_weight="balanced",
    )
    clf.fit(Xs, y)
    return clf.coef_.ravel().astype(np.float64)


def cosine(u: np.ndarray, v: np.ndarray) -> float:
    nu = np.linalg.norm(u)
    nv = np.linalg.norm(v)
    if nu == 0 or nv == 0:
        return float("nan")
    return float(np.dot(u, v) / (nu * nv))


def per_fm(
    fm_h5ad: Path,
    label_col: str,
    disease_col: str,
    covid_survived: set[str],
    sepsis_survived: set[str],
) -> dict | None:
    model_key = fm_h5ad.stem.split("_")[0].lower()
    log.info("=== %s (%s) ===", model_key, fm_h5ad.name)

    adata = sc.read_h5ad(fm_h5ad)
    X, obsm_key = find_embedding(adata, model_key)
    log.info("  obsm=%s shape=%s", obsm_key, X.shape)

    if label_col not in adata.obs.columns:
        log.warning("  %s missing %s in obs; skipping", model_key, label_col)
        return None
    if disease_col not in adata.obs.columns:
        log.warning("  %s missing %s in obs; skipping", model_key, disease_col)
        return None

    severity = np.asarray([str(s).lower() for s in adata.obs[label_col].values])
    disease = adata.obs[disease_col].values

    covid_mask = disease_mask(disease, COVID_TOKENS)
    sepsis_mask = disease_mask(disease, SEPSIS_TOKENS)
    log.info("  covid cells=%d, sepsis cells=%d", covid_mask.sum(), sepsis_mask.sum())

    fatal = severity == "fatal"

    # COVID: fatal vs survived (mild/moderate/recovered/severe)
    covid_valid = covid_mask & (fatal | np.isin(severity, list(covid_survived)))
    yc = np.zeros(covid_valid.sum(), dtype=np.int8)
    yc[fatal[covid_valid]] = 1

    # Sepsis: fatal vs ICU-survived (severe)
    sepsis_valid = sepsis_mask & (fatal | np.isin(severity, list(sepsis_survived)))
    ys = np.zeros(sepsis_valid.sum(), dtype=np.int8)
    ys[fatal[sepsis_valid]] = 1

    log.info("  COVID train: n=%d, fatal=%d", len(yc), int(yc.sum()))
    log.info("  Sepsis train: n=%d, fatal=%d", len(ys), int(ys.sum()))

    if yc.sum() < 5 or (len(yc) - yc.sum()) < 5:
        log.warning("  COVID class imbalance too extreme; skipping")
        return None
    if ys.sum() < 5 or (len(ys) - ys.sum()) < 5:
        log.warning("  Sepsis class imbalance too extreme; skipping")
        return None

    Xc = X[covid_valid]
    Xs = X[sepsis_valid]

    # SHARED scaler: fit on the union so both LR coefficient vectors live in
    # the same standardized basis (required for a meaningful cosine comparison
    # across diseases).
    Xunion = np.vstack([Xc, Xs])
    shared_scaler = StandardScaler(with_mean=True, with_std=True)
    shared_scaler.fit(Xunion)

    log.info("  Fitting COVID LR (shared scaler)...")
    w_covid_shared = fit_lr_direction(Xc, yc, scaler=shared_scaler)
    log.info("  Fitting Sepsis LR (shared scaler)...")
    w_sepsis_shared = fit_lr_direction(Xs, ys, scaler=shared_scaler)

    cos_shared = cosine(w_covid_shared, w_sepsis_shared)
    dot_shared = float(np.dot(w_covid_shared, w_sepsis_shared))
    log.info("  cos_shared(w_covid_fatal, w_sepsis_fatal) = %.4f", cos_shared)

    # Legacy per-disease scaler (for comparison; not the right metric)
    log.info("  Fitting COVID LR (per-disease scaler, legacy)...")
    w_covid_local = fit_lr_direction(Xc, yc, scaler=None)
    log.info("  Fitting Sepsis LR (per-disease scaler, legacy)...")
    w_sepsis_local = fit_lr_direction(Xs, ys, scaler=None)
    cos_local = cosine(w_covid_local, w_sepsis_local)
    log.info("  cos_local(per-disease scalers, legacy) = %.4f", cos_local)

    return dict(
        model=model_key,
        n_covid=int(len(yc)),
        n_covid_fatal=int(yc.sum()),
        n_sepsis=int(len(ys)),
        n_sepsis_fatal=int(ys.sum()),
        embedding_dim=int(X.shape[1]),
        cos_shared_scaler=cos_shared,
        dot_shared_scaler=dot_shared,
        w_covid_shared_norm=float(np.linalg.norm(w_covid_shared)),
        w_sepsis_shared_norm=float(np.linalg.norm(w_sepsis_shared)),
        cos_local_scalers_legacy=cos_local,
        obsm_key=obsm_key,
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--embedding-dir", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--label-col", default="covid_severity")
    p.add_argument("--disease-col", default="Disease")
    p.add_argument("--covid-survived", default="mild,moderate,recovered,severe")
    p.add_argument("--sepsis-survived", default="severe")
    args = p.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    covid_survived = {t.strip().lower() for t in args.covid_survived.split(",")}
    sepsis_survived = {t.strip().lower() for t in args.sepsis_survived.split(",")}

    fm_files = sorted(args.embedding_dir.glob("*_covid_canonical.h5ad"))
    if not fm_files:
        log.error("No *_covid_canonical.h5ad files in %s", args.embedding_dir)
        sys.exit(1)

    log.info("Found %d FM files: %s", len(fm_files), [f.name for f in fm_files])

    rows = []
    for fm in fm_files:
        try:
            row = per_fm(fm, args.label_col, args.disease_col, covid_survived, sepsis_survived)
            if row is not None:
                rows.append(row)
        except Exception as e:
            log.exception("Failed %s: %s", fm.name, e)

    if not rows:
        log.error("No FMs produced valid output")
        sys.exit(1)

    df = pd.DataFrame(rows).sort_values("cos_shared_scaler", ascending=False)
    out_csv = args.output_dir / "cosine_fatal_direction.csv"
    df.to_csv(out_csv, index=False)
    log.info("Wrote %s", out_csv)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
