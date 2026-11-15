"""Tests for calibration_per_cohort.py — ECE math + mortality label construction.

Pure-math unit tests that don't require Option C embeddings to exist.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(
    0,
    str(Path(__file__).parent.parent / "scripts/24_option_c_canonical_input/code"),
)
from calibration_per_cohort import (
    build_mortality_labels,
    expected_calibration_error,
    find_embedding,
)


class TestExpectedCalibrationError:
    """ECE = expected gap between predicted confidence and observed accuracy per bin."""

    def test_perfect_calibration_gives_low_ece(self):
        """Predicted probabilities matching observed frequencies should yield ECE near 0."""
        rng = np.random.default_rng(42)
        n = 5000
        # Generate well-calibrated probabilities: y_prob ~ U[0,1], y_true ~ Bernoulli(y_prob)
        y_prob = rng.uniform(0, 1, n)
        y_true = (rng.uniform(0, 1, n) < y_prob).astype(int)

        ece = expected_calibration_error(y_true, y_prob, n_bins=10)
        assert 0.0 <= ece <= 0.05  # Well-calibrated → small ECE

    def test_overconfident_predictions_give_high_ece(self):
        """If model says 0.9 confidence but only 0.5 are actually positive, ECE should be high."""
        n = 1000
        y_prob = np.full(n, 0.9)
        y_true = np.zeros(n, dtype=int)
        y_true[: n // 2] = 1  # 50% positive but model says 90%
        ece = expected_calibration_error(y_true, y_prob, n_bins=10)
        # Gap is |0.5 - 0.9| = 0.4
        assert ece > 0.3

    def test_perfect_separator_gives_zero_ece(self):
        """Hard predictions (0 or 1) that match labels exactly → ECE near 0."""
        y_true = np.array([1] * 50 + [0] * 50)
        y_prob = y_true.astype(float)  # exactly 0 or 1
        ece = expected_calibration_error(y_true, y_prob, n_bins=10)
        assert ece < 1e-6

    def test_ece_in_zero_one_range(self):
        rng = np.random.default_rng(7)
        for _ in range(5):
            n = rng.integers(50, 500)
            y_true = rng.integers(0, 2, n)
            y_prob = rng.uniform(0, 1, n)
            ece = expected_calibration_error(y_true, y_prob, n_bins=10)
            assert 0.0 <= ece <= 1.0

    def test_more_bins_produces_finer_estimate(self):
        """More bins → at least as sensitive to local miscalibration."""
        rng = np.random.default_rng(0)
        n = 2000
        y_prob = rng.uniform(0, 1, n)
        y_true = (rng.uniform(0, 1, n) < y_prob).astype(int)

        ece_5 = expected_calibration_error(y_true, y_prob, n_bins=5)
        ece_20 = expected_calibration_error(y_true, y_prob, n_bins=20)
        # Both should be small, but more bins shouldn't blow up
        assert 0 <= ece_5 < 0.1
        assert 0 <= ece_20 < 0.15


class TestBuildMortalityLabels:
    """build_mortality_labels() encodes the 2026-05-01 primary task."""

    def test_fatal_one_survived_zero(self):
        sev = np.array(["fatal", "fatal", "mild", "moderate", "recovered", "severe"])
        y, valid = build_mortality_labels(sev)
        assert (y[:2] == 1).all()
        assert (y[2:] == 0).all()
        assert valid.all()

    def test_healthy_excluded(self):
        sev = np.array(["fatal", "healthy", "severe"])
        y, valid = build_mortality_labels(sev)
        assert valid.tolist() == [True, False, True]

    def test_unknown_excluded(self):
        sev = np.array(["fatal", "unknown", "Unknown", "UNKNOWN"])
        y, valid = build_mortality_labels(sev)
        # Function lowercases first → only "fatal" survives
        assert valid.tolist() == [True, False, False, False]

    def test_severe_counts_as_survived(self):
        """COVID-severe and sepsis-severe both count as ICU-survived (label=0)."""
        sev = np.array(["severe"] * 10)
        y, valid = build_mortality_labels(sev)
        assert valid.all()
        assert (y == 0).all()


class TestFindEmbedding:
    """find_embedding() locates a model's embedding via OBSM_KEYS or fallback heuristic."""

    def test_finds_canonical_key(self):
        import anndata

        n = 50
        d = 2058
        a = anndata.AnnData(
            X=np.zeros((n, 100)),
            obsm={"X_state": np.random.randn(n, d).astype(np.float32)},
        )
        X, key = find_embedding(a, "state")
        assert key == "X_state"
        assert X.shape == (n, d)

    def test_falls_back_to_alternate_key(self):
        import anndata

        n = 50
        d = 1280
        a = anndata.AnnData(
            X=np.zeros((n, 100)),
            obsm={"X_uce_4layer": np.random.randn(n, d).astype(np.float32)},
        )
        X, key = find_embedding(a, "uce")
        assert key == "X_uce_4layer"

    def test_raises_when_unknown_key_only(self):
        """Unfamiliar obsm key is NOT silently picked up (fallback was removed
        per the 2026-05-08 obsm-safety review: silent fallback could grab UMAP
        coords or other auxiliary 2D arrays). find_embedding must raise."""
        import anndata

        n = 50
        a = anndata.AnnData(
            X=np.zeros((n, 100)),
            obsm={"unfamiliar_emb_name": np.random.randn(n, 256).astype(np.float32)},
        )
        with pytest.raises(RuntimeError, match="No known embedding key"):
            find_embedding(a, "state")

    def test_raises_when_no_valid_embedding(self):
        import anndata

        n = 50
        a = anndata.AnnData(X=np.zeros((n, 100)))
        with pytest.raises(RuntimeError, match="No known embedding key"):
            find_embedding(a, "state")
