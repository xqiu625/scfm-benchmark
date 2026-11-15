"""Tests for scib_batch_asw(), moved 2026-10-09 to
scripts/00_shared/code/batch_asw.py so every caller shares one
implementation (embedding_quality_option_c.py imports it from there now;
see that module's docstring and docs/decisions.md).

Pure-math unit tests on synthetic data -- no HPCC embeddings needed. These are
the cases a 2026-10-09 review of the 2026-10-08 batch-ASW fix asked to see
committed and reproducible, not just reported as a one-off smoke-test log:
bio-driven separation without batch separation, batch separation within
biology, a bio-confounded batch label (every stratum must be skipped, not
scored), the n_batches==n_cells sklearn boundary, the min_group_size policy
floor, and input-validation failures.
"""

import math
import sys
from pathlib import Path

import numpy as np
import pytest
from sklearn.metrics import silhouette_samples

sys.path.insert(
    0,
    str(Path(__file__).parent.parent / "scripts/00_shared/code"),
)
from batch_asw import scib_batch_asw


def _two_bio_two_batch(n_per_cell=400, batch_shift=0.0, bio_shift=5.0, seed=0):
    """3 bio strata x 2 batches, each stratum balanced across batches."""
    rng = np.random.RandomState(seed)
    bio = np.repeat(["control", "mild", "severe"], n_per_cell)
    batch = np.tile(["A", "B"], (n_per_cell * 3) // 2)
    bio_offset = np.where(bio == "control", 0, np.where(bio == "mild", bio_shift, 2 * bio_shift))
    batch_offset = np.where(batch == "A", 0.0, batch_shift)
    X = rng.normal(size=(n_per_cell * 3, 8))
    X = X * 0.5 + bio_offset[:, None] + batch_offset[:, None]
    return X, batch, bio


def test_bio_separated_batch_mixed_scores_near_perfect_mixing():
    X, batch, bio = _two_bio_two_batch(batch_shift=0.0, bio_shift=8.0)
    r = scib_batch_asw(X, batch, bio)
    assert r["skipped_groups"] == {}
    assert r["batch_asw_neg_scib"] > 0.8, r


def test_batch_separated_within_bio_scores_near_zero_mixing():
    X, batch, bio = _two_bio_two_batch(batch_shift=40.0, bio_shift=8.0)
    r = scib_batch_asw(X, batch, bio)
    assert r["skipped_groups"] == {}
    assert r["batch_asw_neg_scib"] < 0.2, r


def test_batch_label_identical_to_bio_label_all_strata_skipped_as_undefined():
    """The exact failure mode the 2026-10-08 fix targets: if batch == bio,
    the OLD global silhouette_score(X, batch) would read pure biology as a
    severe batch effect. The corrected function must decline to score every
    stratum (n_batches=1 within each) rather than invent a number."""
    rng = np.random.RandomState(1)
    n = 300
    bio = np.repeat(["control", "mild", "severe"], n // 3)
    X = rng.normal(size=(n, 8)) + np.where(
        bio[:, None] == "control", 0, np.where(bio[:, None] == "mild", 5, 10)
    )
    batch = bio.copy()
    r = scib_batch_asw(X, batch, bio)
    assert r["asw_batch_per_group"] == {}
    assert all(v.startswith("UNDEFINED: n_batches=1") for v in r["skipped_groups"].values())
    assert math.isnan(r["batch_asw_neg_scib"])


def test_every_cell_its_own_batch_is_skipped_not_crashed():
    """sklearn silhouette requires n_labels <= n_samples - 1; if every cell in
    a stratum has a unique batch label, n_labels == n_samples and the call
    must be skipped with an explicit reason, not raise or silently score."""
    rng = np.random.RandomState(2)
    n = 50
    X = rng.normal(size=(n, 4))
    bio = np.array(["only_group"] * n)
    batch = np.array([f"batch_{i}" for i in range(n)])  # n unique batches == n cells
    r = scib_batch_asw(X, batch, bio, min_group_size=1)
    assert r["asw_batch_per_group"] == {}
    assert "n_batches=50 == n_cells=50" in r["skipped_groups"]["only_group"]
    assert math.isnan(r["batch_asw_neg_scib"])


def test_small_stratum_skipped_by_policy_not_by_math():
    """A stratum with 2 batches and >=2 cells is mathematically scoreable;
    min_group_size is a separate, labeled policy floor."""
    rng = np.random.RandomState(3)
    n = 4
    X = rng.normal(size=(n, 4))
    bio = np.array(["tiny"] * n)
    batch = np.array(["A", "A", "B", "B"])
    r = scib_batch_asw(X, batch, bio, min_group_size=10)
    assert r["skipped_groups"]["tiny"].startswith("POLICY:")
    # same data, policy floor lowered -- now scoreable, not skipped
    r2 = scib_batch_asw(X, batch, bio, min_group_size=2)
    assert "tiny" in r2["asw_batch_per_group"]


def test_unequal_stratum_sizes_each_stratum_weighted_equally_not_by_cell_count():
    """scIB's convention averages per-stratum means with equal weight per
    stratum, not pooled over cells -- a large well-mixed stratum must not
    drown out a small poorly-mixed one."""
    rng = np.random.RandomState(4)
    n_big, n_small = 2000, 20
    bio = np.array(["big"] * n_big + ["small"] * n_small)
    batch_big = np.tile(["A", "B"], n_big // 2)  # perfectly mixed
    batch_small = np.array(
        ["A"] * (n_small // 2) + ["B"] * (n_small // 2)
    )  # separated by construction below
    X_big = rng.normal(size=(n_big, 4))  # batch has no effect -> well mixed
    X_small = rng.normal(size=(n_small, 4)) * 0.1
    X_small += np.where(batch_small == "A", 0.0, 20.0)[:, None]  # batch fully separates "small"
    X = np.vstack([X_big, X_small])
    batch = np.concatenate([batch_big, batch_small])
    r = scib_batch_asw(X, batch, bio, min_group_size=10)
    assert set(r["asw_batch_per_group"]) == {"big", "small"}
    assert r["asw_batch_per_group"]["big"] < 0.1
    assert r["asw_batch_per_group"]["small"] > 0.9
    # equal-weight average of a near-0 and a near-1 stratum must land near 0.5,
    # NOT near 0 (which a cell-count-weighted average would give, since the
    # "big" stratum has 100x the cells).
    assert 0.35 < r["batch_asw_scib"] < 0.65, r


def test_per_stratum_score_matches_abs_silhouette_exactly_not_the_signed_mean():
    """Direct, exact check against the reference computation
    mean(|silhouette_samples(...)|), on a fixture engineered to contain a
    real mix of positive and negative per-cell silhouette values (a few
    points mislabeled into the wrong, well-separated blob) -- so this test
    only passes if np.abs() is actually applied before averaging, not after,
    and not skipped. This is the exact mutation this test is meant to kill
    (removing `np.abs` in scib_batch_asw's per-stratum reduction): dropping
    abs() changes the per-stratum value from mean(|s|) to mean(s), and for
    this fixture those two differ by construction (verified below)."""
    rng = np.random.RandomState(7)
    n_blob = 50
    blob1 = rng.normal(0, 0.3, size=(n_blob, 4))  # near x=0
    blob2 = rng.normal(0, 0.3, size=(n_blob, 4)) + 100.0  # near x=100, all dims shifted
    X = np.vstack([blob1, blob2])
    batch = np.array(["A"] * n_blob + ["B"] * n_blob)
    # Mislabel 10 of blob1's physically-A points as "B": they sit deep inside
    # the A blob but are scored against the far B blob as their "own" label.
    batch[:10] = "B"
    bio = np.array(["g"] * len(X))

    raw = silhouette_samples(X, batch)
    assert (raw < -0.5).sum() >= 8, (
        f"fixture assumption failed: expected most of the 10 mislabeled "
        f"points to have strongly negative silhouette, got {sorted(raw[:10])}"
    )
    assert (
        raw.mean() > 0
    ), "fixture assumption failed: expected the majority-correct labeling to keep the signed mean positive"
    reference_abs_mean = float(np.mean(np.abs(raw)))
    assert reference_abs_mean > raw.mean() + 0.05, (
        "fixture assumption failed: expected abs-mean to exceed signed-mean "
        "by a clear margin given the negative-silhouette minority"
    )

    r = scib_batch_asw(X, batch, bio, min_group_size=1)
    assert r["asw_batch_per_group"]["g"] == pytest.approx(reference_abs_mean, abs=1e-9), (
        f"scib_batch_asw's per-stratum value {r['asw_batch_per_group']['g']} does not match "
        f"the reference mean(|silhouette|) {reference_abs_mean} -- abs() was dropped, "
        f"applied in the wrong place, or something else changed the reduction."
    )
    assert 0.0 <= r["batch_asw_scib"] <= 1.0, r
    assert 0.0 <= r["batch_asw_neg_scib"] <= 1.0, r


def test_non_finite_embedding_raises_rather_than_silently_scoring():
    X = np.array([[1.0, 2.0], [np.nan, 4.0], [5.0, 6.0], [7.0, 8.0]])
    batch = np.array(["A", "A", "B", "B"])
    bio = np.array(["g", "g", "g", "g"])
    with pytest.raises(ValueError, match="non-finite"):
        scib_batch_asw(X, batch, bio)


def test_shape_mismatch_raises():
    X = np.zeros((10, 4))
    batch = np.array(["A"] * 5)  # wrong length
    bio = np.array(["g"] * 10)
    with pytest.raises(ValueError, match="shape mismatch"):
        scib_batch_asw(X, batch, bio)


def test_non_2d_embedding_raises():
    X_1d = np.zeros(10)
    batch = np.array(["A"] * 10)
    bio = np.array(["g"] * 10)
    with pytest.raises(ValueError, match="must be 2-D"):
        scib_batch_asw(X_1d, batch, bio)


def test_non_1d_label_arrays_raise():
    X = np.zeros((10, 4))
    bio = np.array(["g"] * 10)
    batch_2d = np.zeros((10, 2))
    with pytest.raises(ValueError, match="must be 1-D"):
        scib_batch_asw(X, batch_2d, bio)


def test_missing_batch_label_raises_rather_than_silently_including_it():
    X = np.zeros((10, 4))
    batch = np.array(["A", "B"] * 4 + [None, "A"], dtype=object)
    bio = np.array(["g"] * 10)
    with pytest.raises(ValueError, match="missing value"):
        scib_batch_asw(X, batch, bio)


def test_missing_bio_label_raises_rather_than_silently_including_it():
    X = np.zeros((10, 4))
    batch = np.array(["A", "B"] * 5)
    bio = np.array(["g"] * 9 + [None], dtype=object)
    with pytest.raises(ValueError, match="missing value"):
        scib_batch_asw(X, batch, bio)
