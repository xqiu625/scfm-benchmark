#!/usr/bin/env python3
"""Shared scIB-correct batch ASW -- moved here 2026-10-09 so every caller in
this repo uses one tested implementation instead of each maintaining (or
forgetting to maintain) its own copy.

History: `scripts/24_option_c_canonical_input/code/embedding_quality_option_c.py`
originally computed "Batch ASW" as one GLOBAL, unstratified
`silhouette_score(X, cohort)` -- fixed 2026-10-08 to the scIB-correct,
biologically-stratified form below (see docs/decisions.md, 2026-10-08/09
entry). A second, independent, never-fixed copy of the same wrong pattern
was found in `scripts/12_umap_visualization/code/selected_code/
batch_effect_analysis.py::compute_batch_asw` and marked RETIRED, not deleted
(see that module's docstring) -- it now raises if called, pointing here,
rather than silently returning the old misleading number
(2026-10-09 review: "Ensure the retired legacy entry point cannot silently
emit the old misleading metric if invoked").

Usage:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "00_shared/code"))
    from batch_asw import scib_batch_asw
"""
from __future__ import annotations

import numpy as np
from sklearn.metrics import silhouette_samples


def scib_batch_asw(X: np.ndarray, batch_labels: np.ndarray, bio_labels: np.ndarray,
                    metric: str = "euclidean", min_group_size: int = 10) -> dict:
    """scIB-correct batch ASW: absolute per-cell silhouette (batch as label),
    computed within each biological stratum separately, then averaged across
    strata (unweighted mean of per-stratum means, per Luecken et al. 2022,
    https://scib.readthedocs.io/en/latest/api/scib.metrics.silhouette_batch.html).

    A stratum is SKIPPED, and the reason recorded, rather than silently
    omitted or raising, in three distinct cases -- callers must not conflate
    them, since only the first two are a hard mathematical requirement
    (sklearn.metrics.silhouette_samples needs 2 <= n_labels <= n_samples - 1;
    https://scikit-learn.org/stable/modules/generated/sklearn.metrics.silhouette_samples.html):

      1. `n_batches < 2`            -- mathematically undefined (no contrast).
      2. `n_batches == n` (every cell its own batch) -- also mathematically
         undefined (every "cluster" is a singleton); this is the n_labels ==
         n_samples boundary sklearn itself rejects, distinct from case 1.
      3. `n < min_group_size`       -- a PROJECT POLICY floor, not part of
         scIB's or sklearn's definition, kept separate so a reader can tell
         "undefined" apart from "defined but judged too small to trust."

    Validated before any computation (raises, does not warn-and-continue):
      - `X` must be a finite, 2-D array (`n_samples`, `n_features`).
      - `batch_labels` and `bio_labels` must be 1-D and the same length as X.
      - Neither label array may contain a missing value (NaN / None / pandas
        NA) -- a caller must decide how to handle unlabeled cells (drop them,
        or bin them into an explicit "unknown" stratum) BEFORE calling this
        function; silently including them here would pick a policy for the
        caller without saying so.

    Returns dict with per-stratum scores, the aggregate, and which strata were
    usable (tagged by which of the three reasons above), so a caller can see
    exactly what was averaged and why anything was left out.
    """
    X = np.asarray(X)
    batch_labels = np.asarray(batch_labels)
    bio_labels = np.asarray(bio_labels)

    if X.ndim != 2:
        raise ValueError(f"scib_batch_asw: X must be 2-D (n_samples, n_features), got shape {X.shape}")
    if batch_labels.ndim != 1:
        raise ValueError(f"scib_batch_asw: batch_labels must be 1-D, got shape {batch_labels.shape}")
    if bio_labels.ndim != 1:
        raise ValueError(f"scib_batch_asw: bio_labels must be 1-D, got shape {bio_labels.shape}")
    if not np.all(np.isfinite(X)):
        raise ValueError(
            f"scib_batch_asw: X contains non-finite values "
            f"(nan={np.isnan(X).sum()}, inf={np.isinf(X).sum()}); refusing to "
            f"silently hand this to sklearn, which does not validate it."
        )
    if len(batch_labels) != len(X) or len(bio_labels) != len(X):
        raise ValueError(
            f"scib_batch_asw: shape mismatch -- X has {len(X)} rows, "
            f"batch_labels has {len(batch_labels)}, bio_labels has {len(bio_labels)}."
        )
    n_missing_batch = sum(1 for v in batch_labels if v is None or (isinstance(v, float) and np.isnan(v)))
    n_missing_bio = sum(1 for v in bio_labels if v is None or (isinstance(v, float) and np.isnan(v)))
    if n_missing_batch or n_missing_bio:
        raise ValueError(
            f"scib_batch_asw: batch_labels has {n_missing_batch} missing value(s), "
            f"bio_labels has {n_missing_bio} -- filter or explicitly bin these cells "
            f"into a labeled stratum before calling this function; it does not decide "
            f"that policy for you."
        )

    per_stratum = {}
    skipped = {}
    for grp in sorted(set(bio_labels)):
        mask = bio_labels == grp
        n = int(mask.sum())
        batch_in_grp = batch_labels[mask]
        n_batches = len(set(batch_in_grp))
        if n_batches < 2:
            skipped[grp] = f"UNDEFINED: n_batches={n_batches} < 2 (no batch contrast)"
            continue
        if n_batches == n:
            skipped[grp] = (
                f"UNDEFINED: n_batches={n_batches} == n_cells={n} "
                f"(every cell is its own batch; silhouette needs n_labels <= n_samples-1)"
            )
            continue
        if n < min_group_size:
            skipped[grp] = f"POLICY: n={n} < min_group_size={min_group_size}"
            continue
        s = silhouette_samples(X[mask], batch_in_grp, metric=metric)
        per_stratum[grp] = float(np.mean(np.abs(s)))

    if not per_stratum:
        return dict(asw_batch_per_group=per_stratum, skipped_groups=skipped,
                    batch_asw_scib=float("nan"), batch_asw_neg_scib=float("nan"))

    mean_abs_silhouette = float(np.mean(list(per_stratum.values())))
    return dict(
        asw_batch_per_group=per_stratum,
        skipped_groups=skipped,
        batch_asw_scib=mean_abs_silhouette,
        batch_asw_neg_scib=1.0 - mean_abs_silhouette,
    )
