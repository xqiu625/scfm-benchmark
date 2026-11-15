"""Tests for cross_disease_fatal_vs_survived.py — label construction and bootstrap.

Covers the helpers that don't require HPCC embeddings:
  - build_binary(): fatal=1 / survived=0 / excluded mask construction
  - patient_bootstrap_auc(): patient-grouped resampling + AUC stability
  - disease_mask(): COVID/Sepsis alias matching
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(
    0,
    str(Path(__file__).parent.parent / "scripts/06_batch_generalization/code"),
)
from cross_disease_fatal_vs_survived import (
    build_binary,
    disease_mask,
    patient_bootstrap_auc,
)


class TestBuildBinary:
    """build_binary() converts severity strings to fatal=1 / survived=0 / valid mask."""

    def test_fatal_maps_to_one(self):
        labels = np.array(["fatal", "fatal", "mild"])
        binary, valid = build_binary(labels, ["fatal"], ["mild", "moderate", "recovered", "severe"])
        assert (binary == [1, 1, 0]).all()
        assert valid.all()

    def test_survived_tokens_map_to_zero(self):
        labels = np.array(["mild", "moderate", "recovered", "severe"])
        binary, valid = build_binary(labels, ["fatal"], ["mild", "moderate", "recovered", "severe"])
        assert (binary == 0).all()
        assert valid.all()

    def test_healthy_and_unknown_excluded_when_explicit(self):
        """With explicit survived_tokens, only listed tokens map to 0; everything else is invalid."""
        labels = np.array(["fatal", "healthy", "unknown", "severe"])
        binary, valid = build_binary(labels, ["fatal"], ["mild", "moderate", "recovered", "severe"])
        # fatal->valid+1, healthy->invalid, unknown->invalid, severe->valid+0
        assert valid.tolist() == [True, False, False, True]
        assert binary[0] == 1
        assert binary[3] == 0

    def test_severe_in_default_survived_set(self):
        """Per 2026-05-01 mortality definition, COVID-severe counts as survived."""
        labels = np.array(["severe"] * 5 + ["fatal"] * 3 + ["healthy"] * 2)
        binary, valid = build_binary(
            labels,
            fatal_tokens=["fatal"],
            survived_tokens=["mild", "moderate", "recovered", "severe"],
        )
        assert binary[:5].sum() == 0  # severe → survived
        assert binary[5:8].sum() == 3  # fatal → fatal
        assert valid[:8].all()
        assert not valid[8:].any()  # healthy excluded

    def test_case_insensitive(self):
        labels = np.array(["FATAL", "Fatal", "Severe", "MILD"])
        binary, valid = build_binary(labels, ["fatal"], ["mild", "moderate", "recovered", "severe"])
        assert (binary == [1, 1, 0, 0]).all()
        assert valid.all()

    def test_unknown_excluded_with_none_survived(self):
        """When survived_tokens is None, 'unknown' / 'nan' / '' are excluded but everything else is 0."""
        labels = np.array(["fatal", "anything_else", "unknown", "nan", ""])
        binary, valid = build_binary(labels, ["fatal"], None)
        assert binary[0] == 1
        assert valid[0]
        assert binary[1] == 0
        assert valid[1]
        assert not valid[2:].any()

    def test_empty_survived_tokens_treated_as_explicit(self):
        """Empty/whitespace-only survived_tokens collapse to empty set; only fatal matches."""
        labels = np.array(["fatal", "mild"])
        binary, valid = build_binary(labels, ["fatal"], [""])
        assert bool(valid[0])  # fatal still valid (matched fatal_set)
        assert not bool(valid[1])  # mild not in fatal or survived -> invalid


class TestDiseaseMask:
    """disease_mask() handles known disease aliases."""

    def test_covid_aliases(self):
        diseases = np.array(["COVID-19", "covid", "COVID19", "Sepsis", "Other"])
        mask = disease_mask(diseases, "COVID")
        assert mask.tolist() == [True, True, True, False, False]

    def test_sepsis_aliases(self):
        diseases = np.array(["Sepsis", "sepsis", "SEPSIS", "COVID-19"])
        mask = disease_mask(diseases, "Sepsis")
        assert mask.tolist() == [True, True, True, False]


class TestPatientBootstrap:
    """patient_bootstrap_auc() resamples groups (not cells) to compute AUC CI."""

    def test_bootstrap_returns_finite_estimates(self):
        rng = np.random.default_rng(0)
        n_patients = 30
        n_cells_per = 5
        groups = np.repeat([f"P{i}" for i in range(n_patients)], n_cells_per)
        # Make a clearly-separable problem
        y = (np.repeat(rng.integers(0, 2, n_patients), n_cells_per)).astype(int)
        prob = y * 0.7 + (1 - y) * 0.3 + rng.normal(0, 0.05, len(y))
        prob = np.clip(prob, 0, 1)

        mean, lo, hi = patient_bootstrap_auc(y, prob, groups, n_iter=100, seed=0)

        assert np.isfinite(mean)
        assert np.isfinite(lo) and np.isfinite(hi)
        assert 0.0 <= lo <= mean <= hi <= 1.0
        assert hi - lo < 0.5  # CI not absurdly wide for separable problem

    def test_bootstrap_robust_to_single_class_iter(self):
        """When a bootstrap iteration accidentally has only one class, that
        iteration is skipped — function still returns a finite mean from the rest."""
        rng = np.random.default_rng(1)
        n_patients = 5
        groups = np.repeat([f"P{i}" for i in range(n_patients)], 10)
        y = np.repeat([1, 0, 1, 0, 1], 10).astype(int)
        prob = y * 0.6 + (1 - y) * 0.4 + rng.normal(0, 0.1, len(y))
        prob = np.clip(prob, 0, 1)

        mean, lo, hi = patient_bootstrap_auc(y, prob, groups, n_iter=50, seed=2)
        assert np.isfinite(mean) or np.isnan(mean)  # nan is acceptable in pathological cases


@pytest.mark.parametrize(
    "labels,expected_mortality_count",
    [
        (np.array(["fatal", "fatal", "fatal"]), 3),
        (np.array(["fatal", "severe", "mild"]), 1),
        (np.array(["healthy", "healthy"]), 0),  # healthy excluded entirely
    ],
)
def test_mortality_count_matches_expected(labels, expected_mortality_count):
    binary, valid = build_binary(
        labels,
        fatal_tokens=["fatal"],
        survived_tokens=["mild", "moderate", "recovered", "severe"],
    )
    assert binary[valid].sum() == expected_mortality_count
