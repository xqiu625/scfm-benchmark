#!/usr/bin/env python3
"""Embedding-quality scIB-style metrics on Option C 5-FM (COVID).

Per FM, compute the embedding-quality metrics that scIB and the CZI Virtual
Cells Workshop recommend for assessing whether an embedding (a) preserves
biological signal and (b) is invariant to batch:

  Bio conservation:
    ARI(leiden, severity)      — overlap of unsupervised leiden clusters with severity
    NMI(leiden, severity)      — same, mutual-information based
    Bio ASW (severity)         — silhouette width using severity as labels (higher = better)
    Bio ASW (Disease)          — silhouette width using Disease as labels (5-class)

  Batch correction:
    Batch ASW (Data_NO, scIB-correct) — per-cell silhouette with cohort as the label,
                                  absolute value, computed WITHIN each severity stratum
                                  and averaged across strata (scIB's actual
                                  silhouette_batch definition: Luecken et al. 2022,
                                  Nat Methods -- "the absolute silhouette width is
                                  computed per cell type [here: severity stratum] ...
                                  and averaged across cell-type labels"), then reported
                                  as 1 - that mean so higher = better mixing.
    Leiden / Cohort overlap    — fraction of variance in leiden clusters NOT explained
                                  by cohort (higher = batch-invariant clusters)

Output:
  embedding_quality.csv  model × metric → value
  Suitable for the manuscript radar plot figure.

Notes:
- Silhouette is computed on a 5000-cell subsample for speed on 47K cells.
- Leiden resolution = 1.0 (default).
- ARI/NMI computed against severity bin (3-class: control, mild/moderate, severe/fatal/recovered).
- Fixed 2026-10-08: the previous `batch_asw`/`batch_asw_neg` computed ONE global
  silhouette_score(X, cohort) over the whole subsample and reported `1 - abs(mean)`.
  That is not scIB batch ASW -- it conflates genuine biological separation between
  cohorts (e.g. a cohort that is disease-enriched) with technical batch separation,
  because it never conditions on a biological label. The corrected `batch_asw_neg`
  below stratifies by the severity label before computing silhouette, per cell,
  per the scIB convention; the old unstratified value is kept as
  `batch_asw_global_unstratified_DEPRECATED` for before/after comparison only --
  do not read it as a batch-correction metric.
- Moved 2026-10-09: `scib_batch_asw()` itself now lives in
  `scripts/00_shared/code/batch_asw.py` so every caller (this script and the
  new stage-33 diagnostic CLI) shares one tested implementation instead of
  each keeping its own copy. Imported below via the same sys.path pattern
  `primary_benchmark_option_c.py` already uses for `embedding_sanity`.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    silhouette_score,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "00_shared/code"))
from batch_asw import scib_batch_asw  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("embedding_quality")

OBSM_KEYS = {
    "state": ["X_state", "X_emb", "X_STATE"],
    "uce": ["X_uce", "X_uce_4layer"],
    "geneformer": ["X_geneformer", "X_Geneformer"],
    "transcriptformer": ["X_transcriptformer", "X_TranscriptFormer"],
    "scprint2": ["X_scprint2", "X_scprint", "X_scPRINT", "X_scPRINT2", "scprint_emb"],
    "stack": ["X_stack", "X_STACK"],
    "scgpt": ["X_scgpt", "X_scGPT"],
}
SUBSAMPLE_N = 5000  # for silhouette + ARI/NMI to keep runtime tractable
LEIDEN_RES = 1.0
RNG_SEED = 42


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


def severity_to_3class(severity: np.ndarray) -> np.ndarray:
    """control (healthy/recovered) / mild (mild/moderate) / severe (severe/fatal)."""
    sev = np.asarray([str(s).lower() for s in severity])
    out = np.full(len(sev), "unknown", dtype=object)
    out[np.isin(sev, ["healthy", "recovered"])] = "control"
    out[np.isin(sev, ["mild", "moderate"])] = "mild"
    out[np.isin(sev, ["severe", "fatal"])] = "severe"
    return out


def compute_leiden(X: np.ndarray, n_neighbors: int = 30) -> np.ndarray:
    """Build kNN graph and run leiden. Uses sc.pp.neighbors + sc.tl.leiden."""
    a = sc.AnnData(X=X)
    a.obsm["X_emb"] = X
    sc.pp.neighbors(a, n_neighbors=n_neighbors, use_rep="X_emb", random_state=RNG_SEED)
    sc.tl.leiden(a, resolution=LEIDEN_RES, random_state=RNG_SEED, flavor="igraph",
                 n_iterations=2, directed=False)
    return np.asarray(a.obs["leiden"].values).astype(str)


def per_fm(fm_h5ad: Path) -> dict:
    mk = fm_h5ad.stem.split("_")[0].lower()
    log.info("=" * 60)
    log.info("Model: %s", mk)
    log.info("=" * 60)
    a = sc.read_h5ad(fm_h5ad)
    X, key = find_embedding(a, mk)
    log.info("  obsm=%s shape=%s", key, X.shape)

    severity_full = a.obs["covid_severity"].astype(str).values if "covid_severity" in a.obs.columns else None
    disease_full = a.obs["Disease"].astype(str).values if "Disease" in a.obs.columns else None
    cohort_full = a.obs["Data_NO"].astype(str).values if "Data_NO" in a.obs.columns else None

    if severity_full is None or cohort_full is None:
        log.warning("  missing severity or cohort col; skipping")
        return None

    sev3 = severity_to_3class(severity_full)
    valid = sev3 != "unknown"
    X_v = X[valid]
    sev_v = sev3[valid]
    cohort_v = cohort_full[valid]
    disease_v = disease_full[valid] if disease_full is not None else None

    log.info("  Valid cells: %d  | severity classes: %s  | cohorts: %d",
             len(sev_v), set(sev_v), len(set(cohort_v)))

    # Leiden clustering on full embedding (faster on dense)
    log.info("  Computing leiden...")
    leiden = compute_leiden(X_v)
    n_clusters = len(set(leiden))
    log.info("  Leiden clusters: %d", n_clusters)

    # ARI / NMI vs severity (full set)
    ari_sev = adjusted_rand_score(sev_v, leiden)
    nmi_sev = normalized_mutual_info_score(sev_v, leiden)
    log.info("  ARI(leiden, severity)=%.4f  NMI=%.4f", ari_sev, nmi_sev)

    # ARI / NMI vs disease (if available)
    ari_dis = nmi_dis = float("nan")
    if disease_v is not None:
        ari_dis = adjusted_rand_score(disease_v, leiden)
        nmi_dis = normalized_mutual_info_score(disease_v, leiden)
        log.info("  ARI(leiden, Disease)=%.4f  NMI=%.4f", ari_dis, nmi_dis)

    # Silhouette: subsample to SUBSAMPLE_N cells (silhouette is O(n²) — 47K too slow)
    rng = np.random.RandomState(RNG_SEED)
    if len(X_v) > SUBSAMPLE_N:
        sub = rng.choice(len(X_v), size=SUBSAMPLE_N, replace=False)
    else:
        sub = np.arange(len(X_v))
    X_sub = X_v[sub]
    sev_sub = sev_v[sub]
    cohort_sub = cohort_v[sub]

    # Bio ASW (severity): silhouette using severity labels
    asw_bio_sev = silhouette_score(X_sub, sev_sub, metric="euclidean", random_state=RNG_SEED) \
        if len(set(sev_sub)) > 1 else float("nan")
    log.info("  Bio ASW (severity, n=%d) = %.4f", len(sub), asw_bio_sev)

    # Bio ASW (Disease)
    asw_bio_dis = float("nan")
    if disease_v is not None:
        dis_sub = disease_v[sub]
        if len(set(dis_sub)) > 1:
            asw_bio_dis = silhouette_score(X_sub, dis_sub, metric="euclidean", random_state=RNG_SEED)
            log.info("  Bio ASW (Disease, n=%d) = %.4f", len(sub), asw_bio_dis)

    # Batch ASW (Data_NO), scIB-correct: per-cell |silhouette| with cohort as the
    # label, computed WITHIN each severity stratum, averaged across strata.
    # The legacy single global silhouette_score(X, cohort) is kept alongside for
    # before/after comparison only -- see module docstring 2026-10-08 fix note.
    asw_batch_legacy = silhouette_score(X_sub, cohort_sub, metric="euclidean", random_state=RNG_SEED) \
        if len(set(cohort_sub)) > 1 else float("nan")
    asw_batch_legacy_neg = 1 - abs(asw_batch_legacy) if not np.isnan(asw_batch_legacy) else float("nan")

    scib = scib_batch_asw(X_sub, cohort_sub, sev_sub)
    log.info("  Batch ASW legacy (global, cohort, n=%d) = %.4f → 1-|.| = %.4f "
             "[DEPRECATED, not scIB-correct, kept for comparison]",
             len(sub), asw_batch_legacy, asw_batch_legacy_neg)
    log.info("  Batch ASW scIB (per-severity-stratum, cohort) = %.4f → 1-|.| = %.4f "
             "| per-stratum=%s | skipped=%s",
             scib["batch_asw_scib"], scib["batch_asw_neg_scib"],
             scib["asw_batch_per_group"], scib["skipped_groups"])

    # Leiden/cohort overlap: 1 - NMI(leiden, cohort) (higher = clusters NOT explained by batch)
    nmi_cohort = normalized_mutual_info_score(cohort_v, leiden)
    leiden_cohort_invariance = 1 - nmi_cohort
    log.info("  NMI(leiden, cohort)=%.4f → invariance=%.4f", nmi_cohort, leiden_cohort_invariance)

    return dict(
        model=mk,
        obsm_key=key,
        n_dim=X.shape[1],
        n_valid=len(sev_v),
        n_leiden=n_clusters,
        ari_severity=ari_sev,
        nmi_severity=nmi_sev,
        ari_disease=ari_dis,
        nmi_disease=nmi_dis,
        bio_asw_severity=asw_bio_sev,
        bio_asw_disease=asw_bio_dis,
        batch_asw_neg=scib["batch_asw_neg_scib"],
        batch_asw_scib_raw=scib["batch_asw_scib"],
        batch_asw_skipped_strata=str(scib["skipped_groups"]),
        batch_asw_global_unstratified_DEPRECATED=asw_batch_legacy_neg,
        leiden_cohort_invariance=leiden_cohort_invariance,
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--embedding-dir", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    fm_files = sorted(args.embedding_dir.glob("*_covid_canonical.h5ad"))
    if not fm_files:
        log.error("No FMs"); sys.exit(1)

    rows = []
    for fm in fm_files:
        try:
            row = per_fm(fm)
            if row is not None:
                rows.append(row)
        except Exception as e:
            log.exception("Failed %s: %s", fm.name, e)

    df = pd.DataFrame(rows)
    out = args.output_dir / "embedding_quality.csv"
    df.to_csv(out, index=False)
    log.info("Wrote %s", out)
    print("\n=== EMBEDDING QUALITY ===")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
