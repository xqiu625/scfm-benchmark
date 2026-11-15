#!/usr/bin/env python3
"""Synthetic fixture for the stage-33 batch-effect diagnostic CLI
(`run_batch_diagnostics.py`), requested so the CLI can be validated
end-to-end (input -> computation -> output -> reload) without needing the
still-pending Gate-2 retained input or 5-FM embeddings.

Produces a real .h5ad (not an in-memory object the CLI special-cases) with:
  - 2 lineages (Platelet, Monocyte_Myeloid), unequal sizes.
  - 4 technical batches ("cohort"), standing in for the real cohorts.
  - Repeated patient draws: some patients contribute 2 samples (cells
    tagged with the same patient_id, different Sample_id), matching this
    project's actual CV-dependency concern (a patient's draws must stay
    together across folds).
  - A KNOWN, deterministic batch shift (added in expression/PCA space, per
    lineage, same construction already validated in tests/test_batch_asw.py)
    so a correct diagnostic run must detect it before correction and show
    reduced batch separation after Harmony.
  - A KNOWN biological signal independent of batch (outcome-correlated),
    so bio-conservation metrics have something real to detect and so
    correction must not destroy it.
  - A synthetic FM-like embedding (`obsm["X_synthFM"]`), separate from raw
    expression/PCA, to exercise the diagnostic CLI's optional-embedding path
    without needing a real foundation model.

Usage:
    python make_synthetic_batch_fixture.py --out <path.h5ad> --seed 0
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import anndata
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("make_synthetic_batch_fixture")

LINEAGES = ["Platelet", "Monocyte_Myeloid"]
COHORTS = ["SynCohortA", "SynCohortB", "SynCohortC", "SynCohortD"]
N_GENES = 150
N_PATIENTS_PER_COHORT = 20
CELLS_PER_PATIENT_RANGE = (15, 40)
BATCH_SHIFT_MAGNITUDE = 6.0       # mean log-count shift injected per cohort
BIO_SIGNAL_MAGNITUDE = 4.0        # mean log-count shift injected per outcome


def build(seed: int) -> anndata.AnnData:
    rng = np.random.RandomState(seed)

    # Outcome is a PATIENT-level clinical fact -- assign it ONCE per (cohort,
    # patient), outside and before the lineage loop, and reuse it for every
    # lineage/draw that patient contributes cells to. A 2026-10-09 review
    # found the previous version re-derived outcome inside the lineage loop
    # with `lineage_idx` folded into the seed, so the SAME patient_id could
    # get "fatal" in one lineage and "survived" in another -- biologically
    # incoherent and a real bug, not a stylistic nit. NOT Python's built-in
    # hash(), which is per-process randomized for strings (PYTHONHASHSEED)
    # and was found, in an earlier version of this fixture, to assign every
    # "fatal" patient to a single cohort.
    outcome_by_patient = {}
    n_draws_by_patient = {}
    for cohort_idx, cohort in enumerate(COHORTS):
        for p in range(N_PATIENTS_PER_COHORT):
            patient_id = f"{cohort}_P{p:03d}"
            patient_rng = np.random.RandomState(seed * 1_000_003 + cohort_idx * 101 + p)
            outcome_by_patient[patient_id] = "fatal" if patient_rng.rand() < 0.2 else "survived"
            n_draws_by_patient[patient_id] = 2 if p % 7 == 0 else 1  # ~15% repeated draw

    rows = []
    for lineage_idx, lineage in enumerate(LINEAGES):
        for cohort_idx, cohort in enumerate(COHORTS):
            n_patients = N_PATIENTS_PER_COHORT
            for p in range(n_patients):
                patient_id = f"{cohort}_P{p:03d}"
                outcome = outcome_by_patient[patient_id]
                n_draws = n_draws_by_patient[patient_id]
                for draw in range(n_draws):
                    sample_id = f"{patient_id}_D{draw}"
                    n_cells = rng.randint(*CELLS_PER_PATIENT_RANGE)
                    for c in range(n_cells):
                        rows.append(dict(
                            lineage=lineage, cohort=cohort, cohort_idx=cohort_idx,
                            patient_id=patient_id, sample_id=sample_id,
                            outcome=outcome,
                        ))

    obs = pd.DataFrame(rows)
    obs.index = [f"cell_{i:06d}" for i in range(len(obs))]
    n_cells_total = len(obs)
    log.info("Built obs: %d cells, %d patients, %d samples, lineages=%s, cohorts=%s",
              n_cells_total, obs["patient_id"].nunique(), obs["sample_id"].nunique(),
              obs["lineage"].unique().tolist(), obs["cohort"].unique().tolist())

    # Base expression: Poisson counts, log1p'd implicitly via the shift construction below.
    base = rng.poisson(lam=3.0, size=(n_cells_total, N_GENES)).astype(np.float32)

    # Lineage signal: first 1/3 of genes carry a lineage-specific mean shift
    # (so PCA/leiden on the FULL object would trivially separate lineages --
    # the diagnostic CLI must analyze lineages SEPARATELY, not pooled).
    lineage_genes = slice(0, N_GENES // 3)
    for lineage_idx, lineage in enumerate(LINEAGES):
        mask = (obs["lineage"] == lineage).values
        shift = lineage_idx * 4.0
        base[np.ix_(mask, range(*lineage_genes.indices(N_GENES)))] += shift

    # KNOWN batch shift: middle 1/3 of genes carry a cohort-specific mean
    # shift, independent of outcome -- this is the technical effect a
    # correct batch-ASW/Harmony pipeline must detect and can correct.
    batch_genes = range(N_GENES // 3, 2 * N_GENES // 3)
    for cohort_idx, cohort in enumerate(COHORTS):
        mask = (obs["cohort"] == cohort).values
        shift = cohort_idx * BATCH_SHIFT_MAGNITUDE
        base[np.ix_(mask, list(batch_genes))] += shift

    # KNOWN biological signal: last 1/3 of genes carry an outcome-specific
    # mean shift, independent of cohort -- a correction must NOT erase this.
    bio_genes = range(2 * N_GENES // 3, N_GENES)
    fatal_mask = (obs["outcome"] == "fatal").values
    base[np.ix_(fatal_mask, list(bio_genes))] += BIO_SIGNAL_MAGNITUDE

    base = np.clip(base, 0, None)
    counts = np.round(base).astype(np.int32)

    var = pd.DataFrame(index=[f"GENE{i:04d}" for i in range(N_GENES)])

    adata = anndata.AnnData(X=counts.astype(np.float32), obs=obs, var=var)
    adata.layers["counts"] = counts.copy()

    # A synthetic "FM-like" embedding: a random but DETERMINISTIC projection
    # of the same ground-truth batch+bio signal into a 32-D space, so the
    # diagnostic CLI's optional-embedding path has something non-trivial to
    # analyze without depending on any real foundation model.
    proj = rng.normal(size=(N_GENES, 32)).astype(np.float32)
    synth_fm = np.log1p(counts) @ proj
    adata.obsm["X_synthFM"] = synth_fm
    adata.uns["population_stage"] = "synthetic_fixture"
    adata.uns["population_stage_id"] = "synthetic_fixture"  # machine-checkable; see
    # scripts/00_shared/code/population_stage.py -- run_batch_diagnostics.py's --stage-id
    # is matched exactly against this and hard-fails on any mismatch (2026-10-09 review).
    adata.uns["expression_source"] = "counts"

    adata.obs["lineage"] = adata.obs["lineage"].astype(str)
    adata.obs["cohort"] = adata.obs["cohort"].astype(str)
    adata.obs["patient_id"] = adata.obs["patient_id"].astype(str)
    adata.obs["sample_id"] = adata.obs["sample_id"].astype(str)
    adata.obs["outcome"] = adata.obs["outcome"].astype(str)
    adata.obs = adata.obs.drop(columns=["cohort_idx"])
    return adata


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)

    adata = build(args.seed)

    n_outcomes_per_patient = adata.obs.groupby("patient_id")["outcome"].nunique()
    inconsistent = n_outcomes_per_patient[n_outcomes_per_patient > 1]
    if len(inconsistent):
        raise ValueError(
            f"{len(inconsistent)} patient_id(s) have more than one distinct outcome value "
            f"across their rows -- outcome must be assigned once per patient: "
            f"{inconsistent.index.tolist()[:10]}"
        )
    log.info("Verified: every patient_id has exactly one outcome value across all lineages/draws")

    adata2 = build(args.seed)
    pd.testing.assert_frame_equal(adata.obs, adata2.obs)
    assert np.array_equal(adata.X, adata2.X), "build() is not deterministic given the same seed"
    log.info("Determinism verified: build(seed=%d) called twice gives byte-identical obs and X", args.seed)

    adata.write_h5ad(args.out)
    log.info("Wrote %s: %s", args.out, adata)

    reloaded = anndata.read_h5ad(args.out)
    assert reloaded.shape == adata.shape, "write/reload shape mismatch"
    assert list(reloaded.obs.columns) == list(adata.obs.columns), "write/reload column mismatch"
    assert "X_synthFM" in reloaded.obsm, "write/reload lost the synthetic embedding"
    log.info("Write/reload verified: shape=%s, obs columns=%s, obsm keys=%s",
              reloaded.shape, list(reloaded.obs.columns), list(reloaded.obsm.keys()))


if __name__ == "__main__":
    main()
