#!/usr/bin/env python3
"""Embedding sanity-check module — fail-loud detector for degenerate FM embeddings.

This module exists because of a near-miss with TF-Metazoa where a gene-namespace
tokenization mismatch caused inference to produce a constant-vector output (every
cell mapped to nearly the same point in 2048-D space). The h5ad file existed, the
shape was correct, the LR classifier ran without error, and AUC reported as
0.500 +/- 0.000 across every endpoint — empirically indistinguishable from a real
"tissue-blind" biological negative result. We caught it only after a per-cell-pair
cosine-similarity diagnostic.

This module provides a 1-line check that benchmark scripts call right after
loading an embedding. It raises EmbeddingSanityError if the embedding is
degenerate; otherwise it logs a brief health summary and returns silently.

Usage:
    from scripts_00_shared.code.embedding_sanity import check_embedding_or_raise
    X = adata.obsm["X_state"]
    check_embedding_or_raise(X, fm_label="STATE")  # raises on degeneracy

Thresholds are calibrated from the comparison run we did in May 2026:
  - UCE (healthy):       per-dim std mean = 0.0137, cell-pair cos sim mean = 0.758
  - TF-Metazoa (broken): per-dim std mean = 3.4e-5, cell-pair cos sim mean = 0.985
We set thresholds 100x looser than UCE to flag only true collapse.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Iterable

import numpy as np

log = logging.getLogger(__name__)


class EmbeddingSanityError(RuntimeError):
    """Raised when an FM embedding looks degenerate (constant, collapsed, or wrong shape)."""


@dataclass
class EmbeddingHealth:
    n_cells: int
    n_dims: int
    global_std: float
    per_cell_magnitude_mean: float
    per_cell_magnitude_std: float
    per_dim_std_mean: float
    per_dim_std_min: float
    n_zero_var_dims: int
    cell_pair_cos_sim_mean: float
    cell_pair_cos_sim_min: float
    cell_pair_cos_sim_max: float
    has_nan: bool
    has_inf: bool

    def summary(self) -> str:
        return (
            f"n={self.n_cells}x{self.n_dims}, "
            f"global_std={self.global_std:.4f}, "
            f"per-cell mag std={self.per_cell_magnitude_std:.4f}, "
            f"per-dim std mean={self.per_dim_std_mean:.4f} "
            f"(zero-var dims: {self.n_zero_var_dims}/{self.n_dims}), "
            f"cell-pair cos sim mean={self.cell_pair_cos_sim_mean:.4f} "
            f"[min={self.cell_pair_cos_sim_min:.3f}, max={self.cell_pair_cos_sim_max:.3f}]"
        )


def compute_health(
    X: np.ndarray,
    n_sample_pairs: int = 200,
    seed: int = 0,
) -> EmbeddingHealth:
    """Compute a battery of health metrics on an embedding without raising."""
    X = np.asarray(X)
    if X.ndim != 2:
        raise EmbeddingSanityError(
            f"Embedding has ndim={X.ndim}, expected 2 (cells x dims)"
        )
    n_cells, n_dims = X.shape

    has_nan = bool(np.isnan(X).any())
    has_inf = bool(np.isinf(X).any())
    if has_nan or has_inf:
        return EmbeddingHealth(
            n_cells=n_cells, n_dims=n_dims,
            global_std=float("nan"),
            per_cell_magnitude_mean=float("nan"),
            per_cell_magnitude_std=float("nan"),
            per_dim_std_mean=float("nan"),
            per_dim_std_min=float("nan"),
            n_zero_var_dims=-1,
            cell_pair_cos_sim_mean=float("nan"),
            cell_pair_cos_sim_min=float("nan"),
            cell_pair_cos_sim_max=float("nan"),
            has_nan=has_nan, has_inf=has_inf,
        )

    # Magnitude statistics
    mags = np.linalg.norm(X, axis=1)
    per_cell_mag_mean = float(mags.mean())
    per_cell_mag_std = float(mags.std())

    # Per-dimension variance
    per_dim_std = X.std(axis=0)
    per_dim_std_mean = float(per_dim_std.mean())
    per_dim_std_min = float(per_dim_std.min())
    n_zero_var_dims = int((per_dim_std < 1e-9).sum())

    # Cell-pair cosine similarity on a sample
    rng = np.random.default_rng(seed)
    n_sample = min(n_sample_pairs, n_cells)
    if n_sample >= 2:
        idx = rng.choice(n_cells, n_sample, replace=False)
        sub = X[idx]
        norms = np.linalg.norm(sub, axis=1, keepdims=True)
        norms = np.where(norms < 1e-12, 1e-12, norms)
        sub_norm = sub / norms
        cos_sim = sub_norm @ sub_norm.T
        np.fill_diagonal(cos_sim, np.nan)
        cs_mean = float(np.nanmean(cos_sim))
        cs_min = float(np.nanmin(cos_sim))
        cs_max = float(np.nanmax(cos_sim))
    else:
        cs_mean = cs_min = cs_max = float("nan")

    return EmbeddingHealth(
        n_cells=n_cells, n_dims=n_dims,
        global_std=float(X.std()),
        per_cell_magnitude_mean=per_cell_mag_mean,
        per_cell_magnitude_std=per_cell_mag_std,
        per_dim_std_mean=per_dim_std_mean,
        per_dim_std_min=per_dim_std_min,
        n_zero_var_dims=n_zero_var_dims,
        cell_pair_cos_sim_mean=cs_mean,
        cell_pair_cos_sim_min=cs_min,
        cell_pair_cos_sim_max=cs_max,
        has_nan=has_nan, has_inf=has_inf,
    )


def check_embedding_or_raise(
    X: np.ndarray,
    fm_label: str = "<unknown>",
    *,
    min_cells: int = 100,
    min_dims: int = 16,
    min_per_cell_mag_std: float = 1e-3,
    min_per_dim_std_mean: float = 1e-4,
    max_zero_var_frac: float = 0.5,
    max_cell_pair_cos_mean: float = 0.95,
    cos_sim_discriminative_std: float = 5e-3,
    n_sample_pairs: int = 200,
    seed: int = 0,
) -> EmbeddingHealth:
    """Check an embedding for degeneracy. Raise EmbeddingSanityError on failure.

    Failure criteria (any one of these triggers EmbeddingSanityError):
      - per_dim_std_mean < min_per_dim_std_mean (default 1e-4)
        TF-Metazoa broken had 3.4e-5; UCE healthy had 0.014
      - cell_pair_cos_sim_mean > max_cell_pair_cos_mean (default 0.95)
        AND per_dim_std_mean < cos_sim_discriminative_std (default 5e-3) — i.e. high
        cosine similarity ONLY hard-fails when it co-occurs with collapsed per-dim
        variance (true direction-and-spread collapse). High cosine similarity *with*
        healthy per-dim variance is demoted to a WARNING, not a failure: some FMs are
        legitimately high-similarity yet discriminative (e.g. Geneformer on platelets —
        rank-value tokenization is dominated by a shared ribosomal signature, giving
        cos sim ~0.95 while still yielding real mortality AUC ~0.70–0.74). TF-Metazoa's
        true collapse (cos 0.985, per-dim std 3.4e-5) is still caught — both by this
        combined gate and independently by the per_dim_std_mean criterion above.
      - n_zero_var_dims / n_dims > max_zero_var_frac (default 0.5)
      - n_cells < min_cells, n_dims < min_dims, NaN/Inf present

    Per-cell magnitude std is REPORTED but is NOT a failure criterion on its own:
    L2-normalized embeddings (every cell on the unit sphere) have magnitude std = 0
    yet are healthy. We log an INFO if magnitude std is suspiciously low while
    direction-based checks pass (could indicate L2-normalization or an exotic
    output mode worth investigating).

    On success, logs a one-line summary and returns the EmbeddingHealth dataclass.
    On failure, raises EmbeddingSanityError with a multi-line diagnostic message.
    """
    health = compute_health(X, n_sample_pairs=n_sample_pairs, seed=seed)
    fails = []

    if health.n_cells < min_cells:
        fails.append(f"too few cells: {health.n_cells} < {min_cells}")
    if health.n_dims < min_dims:
        fails.append(f"too few dims: {health.n_dims} < {min_dims}")
    if health.has_nan:
        fails.append("contains NaN values")
    if health.has_inf:
        fails.append("contains Inf values")
    if not health.has_nan and not health.has_inf:
        # Per-dim std and cell-pair cosine are the primary failure criteria.
        # Per-cell magnitude std alone is NOT a failure — L2-normalized embeddings
        # have magnitude std = 0 by construction yet are perfectly healthy. We only
        # treat low magnitude std as confirming a failure that the direction-based
        # checks already flagged (i.e., when cells collapse to the same point in
        # both magnitude and direction).
        if health.per_dim_std_mean < min_per_dim_std_mean:
            fails.append(
                f"per-dim std mean {health.per_dim_std_mean:.6f} "
                f"< threshold {min_per_dim_std_mean} — most dimensions carry no variance"
            )
        zero_var_frac = health.n_zero_var_dims / health.n_dims
        if zero_var_frac > max_zero_var_frac:
            fails.append(
                f"{health.n_zero_var_dims}/{health.n_dims} ({100*zero_var_frac:.0f}%) "
                f"dims have zero variance — embedding mostly collapsed"
            )
        if health.cell_pair_cos_sim_mean > max_cell_pair_cos_mean:
            # High cosine similarity is only a *collapse* if the embedding ALSO lacks
            # per-dim variance. A high-similarity-but-discriminative embedding (real
            # per-dim spread) is legit biology, not a tokenization bug — demote to WARN.
            if health.per_dim_std_mean < cos_sim_discriminative_std:
                fails.append(
                    f"cell-pair cosine similarity mean {health.cell_pair_cos_sim_mean:.4f} "
                    f"> threshold {max_cell_pair_cos_mean} AND per-dim std mean "
                    f"{health.per_dim_std_mean:.6f} < {cos_sim_discriminative_std} — "
                    f"cells map to nearly the same vector with no discriminative spread"
                )
            else:
                log.warning(
                    "[%s] high cell-pair cosine similarity mean=%.4f > %g, but per-dim "
                    "std mean=%.4f >= %g (discriminative spread present) — NOT failing. "
                    "Expected for high-baseline-similarity FMs (e.g. Geneformer's "
                    "ribosomal-dominated platelet signature). Verify downstream AUC != 0.5.",
                    fm_label, health.cell_pair_cos_sim_mean, max_cell_pair_cos_mean,
                    health.per_dim_std_mean, cos_sim_discriminative_std,
                )
        # Magnitude std is reported but not a failure criterion on its own.
        # Log a warning if it's suspiciously low alongside passing direction checks
        # (suggests something exotic — could be L2-normalized, could be a bug).
        if (health.per_cell_magnitude_std < min_per_cell_mag_std
                and not fails  # only warn if direction checks passed
                and health.per_cell_magnitude_mean > 0):
            log.info(
                "[%s] per-cell magnitude std=%.6f < %g — embedding may be L2-normalized "
                "or otherwise constant-magnitude. Direction-based checks pass; not failing.",
                fm_label, health.per_cell_magnitude_std, min_per_cell_mag_std,
            )

    if fails:
        msg = (
            f"\n========================= EMBEDDING SANITY CHECK FAILED =========================\n"
            f"FM: {fm_label}\n"
            f"Health summary: {health.summary()}\n"
            f"Failures:\n  " + "\n  ".join(f"- {f}" for f in fails) + "\n"
            f"\n"
            f"This was the same failure mode as the TF-Metazoa tokenization-mismatch bug\n"
            f"(May 2026): mouse gene symbols passed to a tokenizer expecting ENSMUSG IDs\n"
            f"produced a constant baseline embedding. AUC = 0.500 +/- 0.000 was reported\n"
            f"across every endpoint, indistinguishable from a real biological negative\n"
            f"result. The embedding is degenerate, not the biology.\n"
            f"\n"
            f"Likely causes:\n"
            f"  - Gene namespace mismatch (e.g., symbols where Ensembl IDs are expected)\n"
            f"  - Wrong inference mode / output layer (e.g., shared 'reference' embedding\n"
            f"    instead of per-cell encoding)\n"
            f"  - Numerical underflow (fp16 + small magnitudes)\n"
            f"  - Empty / all-zero input cells\n"
            f"================================================================================="
        )
        raise EmbeddingSanityError(msg)

    log.info("[%s] embedding health OK: %s", fm_label, health.summary())
    return health


def check_embedding_or_warn(
    X: np.ndarray, fm_label: str = "<unknown>", **kwargs
) -> tuple[EmbeddingHealth, list[str]]:
    """Non-raising variant — returns (health, list_of_failures). Use when you
    want to log the failure but skip the FM rather than abort the whole batch.
    """
    try:
        h = check_embedding_or_raise(X, fm_label=fm_label, **kwargs)
        return h, []
    except EmbeddingSanityError as e:
        msg = str(e).split("\n")
        # Extract just the bullet failures
        fails = [m.strip("- ").strip() for m in msg if m.strip().startswith("- ")]
        log.error("[%s] embedding health FAILED: %s", fm_label, "; ".join(fails))
        return compute_health(X), fails


__all__ = [
    "EmbeddingHealth",
    "EmbeddingSanityError",
    "check_embedding_or_raise",
    "check_embedding_or_warn",
    "compute_health",
]
