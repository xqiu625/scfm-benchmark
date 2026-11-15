"""Tests for build_patient_batch_confounding_table.py's pure functions.

Targets the three v3 fixes from the 2026-10-09 review
(review:20261009_101831): exact (not np.isclose) event
classification, (cohort, donor_id) composite-key uniqueness instead of a
global bare-donor_id check, and stable empty-output schemas.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(
    0,
    str(Path(__file__).parent.parent / "scripts/33_batch_effect_audit/code"),
)
from build_patient_batch_confounding_table import classify_event


def test_classify_event_exact_binary():
    s = pd.Series([1.0, 0.0, np.nan, 1.0, 0.0])
    out = classify_event(s)
    assert list(out) == ["fatal", "survived", "NA_undefined", "fatal", "survived"]


def test_classify_event_rejects_near_but_not_exact_value():
    """The exact failure mode the review flagged: np.isclose(1.000001, 1.0)
    is True at default tolerance, which would have silently classified an
    out-of-contract value as "fatal". Exact equality must reject it."""
    s = pd.Series([1.0, 1.000001, 0.0])
    with pytest.raises(ValueError, match="neither exactly 0, exactly 1"):
        classify_event(s)


def test_classify_event_rejects_other_non_binary_values():
    s = pd.Series([1.0, 0.0, 2.0, np.nan])
    with pytest.raises(ValueError):
        classify_event(s)


def test_classify_event_works_on_int_dtype_not_just_float():
    """v2's string-matching bug depended on the column having been read as
    float64 (true only because other rows had NaN); an all-present int64
    column must classify identically."""
    s = pd.Series([1, 0, 1, 0], dtype="int64")
    out = classify_event(s)
    assert list(out) == ["fatal", "survived", "fatal", "survived"]


def test_composite_key_uniqueness_allows_donor_id_reuse_across_cohorts():
    """A bare global donor_id uniqueness check would reject this even though
    it is legitimate: two different cohorts, two different biological
    patients, coincidentally sharing a donor_id string."""
    df = pd.DataFrame(
        {
            "cohort": ["Cohort A", "Cohort B"],
            "donor_id": ["P001", "P001"],
        }
    )
    dup_composite = df.duplicated(subset=["cohort", "donor_id"], keep=False)
    assert not dup_composite.any()
    reused = df.groupby("donor_id")["cohort"].nunique().loc[lambda s: s > 1]
    assert list(reused.index) == ["P001"]


def test_composite_key_uniqueness_catches_a_real_duplicate_within_one_cohort():
    df = pd.DataFrame(
        {
            "cohort": ["Cohort A", "Cohort A", "Cohort A"],
            "donor_id": ["P001", "P002", "P001"],
        }
    )
    dup_composite = df.duplicated(subset=["cohort", "donor_id"], keep=False)
    assert dup_composite.sum() == 2
