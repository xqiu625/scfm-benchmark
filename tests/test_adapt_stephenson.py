"""Tests calling adapt_stephenson() itself, at the actual production
boundary -- a 2026-10-09 review noted the previous test round added
collision/union-find/event-classification tests but none that actually
invoked this adapter, so a regression in it could pass every other test.

Builds a minimal T3/T4 fixture under pytest's tmp_path, matching the real
column schemas exactly (stephenson_patient_outcomes.csv,
stephenson_pool_composition.csv, stephenson_sample_severity.csv).
"""

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(
    0,
    str(Path(__file__).parent.parent / "scripts/33_batch_effect_audit/code"),
)
from build_batch_metadata_inventory import T3_DIR, T4_DIR, adapt_stephenson


def _write_fixture(root: Path, patients: list[dict], pools: list[dict], severities: list[dict]):
    """pd.DataFrame([]).to_csv(...) with an empty list writes a file with no
    header at all, which pd.read_csv then refuses with EmptyDataError -- not
    a production bug, just a fixture-construction detail. Give an empty
    pools/severities list the right header-only schema instead."""
    t3 = root / T3_DIR
    t4 = root / T4_DIR
    t3.mkdir(parents=True, exist_ok=True)
    t4.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(patients).to_csv(t3 / "stephenson_patient_outcomes.csv", index=False)
    pools_df = (
        pd.DataFrame(pools)
        if pools
        else pd.DataFrame(columns=["pool", "n_donors", "donors", "distinct_outcomes", "outcomes"])
    )
    pools_df.to_csv(t3 / "stephenson_pool_composition.csv", index=False)
    sev_df = (
        pd.DataFrame(severities) if severities else pd.DataFrame(columns=["sample_id", "n_cells"])
    )
    sev_df.to_csv(t4 / "stephenson_sample_severity.csv", index=False)


def _base_patient(**overrides):
    row = dict(
        patient_id="P1",
        n_samples=1,
        sample_ids="P1",
        pools="",
        Outcome="Home",
        outcome_unique_within_patient=True,
        Days_from_onset=1,
        Collection_Days="D0",
        Resample="Initial",
        Site="Sanger",
        n_cells=100,
    )
    row.update(overrides)
    return row


def test_unpooled_donor_gets_unknown_status_not_a_literal_nan_string(tmp_path):
    """The exact bug a 2026-10-09 review found: a blank `pools` CSV field
    must not become the literal string "nan" in pool_id, and must not be
    read as confirmed_single_donor."""
    _write_fixture(
        tmp_path,
        patients=[_base_patient(patient_id="AP1", sample_ids="AP1", pools="", Site="Sanger")],
        pools=[],
        severities=[dict(sample_id="AP1", n_cells=100)],
    )
    out = adapt_stephenson(tmp_path)
    assert len(out) == 1
    row = out.iloc[0]
    assert pd.isna(row["pool_id"])  # genuinely missing, not the 3-char string "nan"
    assert type(row["pool_id"]) is not str
    assert row["pooling_status"] == "unknown"
    # is_multi_donor_capture is nullable: "unknown" pooling must be pd.NA, never False
    # (False would silently read as "confirmed independent" -- a 2026-10-09 review finding).
    assert pd.isna(row["is_multi_donor_capture"])


def test_cambridge_pooled_donor_gets_confirmed_pooled_status(tmp_path):
    _write_fixture(
        tmp_path,
        patients=[
            _base_patient(
                patient_id="CV0001", sample_ids="BGCV01_CV0001", pools="BGCV01", Site="Cambridge"
            )
        ],
        pools=[
            dict(pool="BGCV01", n_donors=1, donors="CV0001", distinct_outcomes=1, outcomes="Home")
        ],
        severities=[dict(sample_id="BGCV01_CV0001", n_cells=500)],
    )
    out = adapt_stephenson(tmp_path)
    assert len(out) == 1
    row = out.iloc[0]
    assert row["pool_id"] == "BGCV01"
    assert row["pooling_status"] == "confirmed_pooled"
    assert row["is_multi_donor_capture"] == True  # noqa: E712
    assert row["n_cells_capture"] == 500


def test_non_cambridge_site_with_a_pool_record_is_not_confirmed_pooled(tmp_path):
    """A pool record exists but the site has no primary-source confirmation
    (only Cambridge does) -- must not be upgraded to confirmed_pooled."""
    _write_fixture(
        tmp_path,
        patients=[_base_patient(patient_id="N1", sample_ids="X1", pools="POOLX", Site="Ncl")],
        pools=[dict(pool="POOLX", n_donors=1, donors="N1", distinct_outcomes=1, outcomes="Home")],
        severities=[dict(sample_id="X1", n_cells=200)],
    )
    out = adapt_stephenson(tmp_path)
    row = out.iloc[0]
    assert row["pooling_status"] == "unknown"
    assert pd.isna(row["is_multi_donor_capture"])  # NA, not False -- see fixture above


def test_sample_to_pool_positional_mismatch_against_pool_composition_raises(tmp_path):
    """sample_ids and pools have equal length (1 and 1) but the claimed pool
    does not actually list this patient as a member -- must raise, not
    silently trust positional alignment."""
    _write_fixture(
        tmp_path,
        patients=[
            _base_patient(
                patient_id="CV0001", sample_ids="BGCV01_CV0001", pools="BGCV01", Site="Cambridge"
            )
        ],
        pools=[
            dict(
                pool="BGCV01",
                n_donors=1,
                donors="CV9999",  # wrong donor!
                distinct_outcomes=1,
                outcomes="Home",
            )
        ],
        severities=[dict(sample_id="BGCV01_CV0001", n_cells=500)],
    )
    with pytest.raises(ValueError, match="does not contain"):
        adapt_stephenson(tmp_path)


def test_duplicate_t4_sample_id_raises_before_building_count_lookup(tmp_path):
    _write_fixture(
        tmp_path,
        patients=[_base_patient(patient_id="P1", sample_ids="S1", pools="")],
        pools=[],
        severities=[dict(sample_id="S1", n_cells=100), dict(sample_id="S1", n_cells=999)],
    )
    with pytest.raises(ValueError, match="duplicate sample_id"):
        adapt_stephenson(tmp_path)


def test_repeated_draw_patient_produces_one_row_per_sample_with_correct_counts(tmp_path):
    """The exact review example: CV0068 with two samples, 2,701 and 2,250
    cells, must appear as TWO rows, not one row with the 4,951 patient
    aggregate."""
    _write_fixture(
        tmp_path,
        patients=[
            _base_patient(
                patient_id="CV0068",
                n_samples=2,
                sample_ids="BGCV02_CV0068;BGCV12_CV0068",
                pools="BGCV02;BGCV12",
                Site="Cambridge",
                n_cells=4951,
            )
        ],
        pools=[
            dict(pool="BGCV02", n_donors=1, donors="CV0068", distinct_outcomes=1, outcomes="Home"),
            dict(pool="BGCV12", n_donors=1, donors="CV0068", distinct_outcomes=1, outcomes="Home"),
        ],
        severities=[
            dict(sample_id="BGCV02_CV0068", n_cells=2701),
            dict(sample_id="BGCV12_CV0068", n_cells=2250),
        ],
    )
    out = adapt_stephenson(tmp_path)
    assert len(out) == 2
    counts = sorted(out["n_cells_capture"].tolist())
    assert counts == [2250, 2701]
    assert 4951 not in counts


def test_swapped_pool_order_is_not_caught_known_limitation(tmp_path):
    """Documents a real, unresolved validation gap (2026-10-09 review): when
    a patient is a genuine member of TWO pools, swapping the order of the
    `pools` column relative to `sample_ids` is NOT caught, because the
    check is set-membership ("is pid a member of this pool"), not a
    sample-specific pairing proof -- no local evidence source ties a
    specific sample_id to a specific pool more directly than this same
    (sample_ids, pools) column pair. This test exists to keep that gap
    visible and regression-tested as a known limitation, not to claim it
    is fixed."""
    _write_fixture(
        tmp_path,
        patients=[
            _base_patient(
                patient_id="CV0068",
                n_samples=2,
                sample_ids="BGCV02_CV0068;BGCV12_CV0068",
                pools="BGCV12;BGCV02",  # SWAPPED relative to sample_ids
                Site="Cambridge",
                n_cells=4951,
            )
        ],
        pools=[
            dict(pool="BGCV02", n_donors=1, donors="CV0068", distinct_outcomes=1, outcomes="Home"),
            dict(pool="BGCV12", n_donors=1, donors="CV0068", distinct_outcomes=1, outcomes="Home"),
        ],
        severities=[
            dict(sample_id="BGCV02_CV0068", n_cells=2701),
            dict(sample_id="BGCV12_CV0068", n_cells=2250),
        ],
    )
    # Does NOT raise, even though pools/sample_ids are swapped -- because CV0068
    # is a genuine member of both pools, membership-only validation cannot tell.
    out = adapt_stephenson(tmp_path)
    assert len(out) == 2
    row_for_bgcv02_sample = out[out["capture_id"] == "BGCV02_CV0068"].iloc[0]
    # The swap means this row's pool_id ends up "BGCV12", not "BGCV02" --
    # exactly the undetected mis-pairing this test documents.
    assert row_for_bgcv02_sample["pool_id"] == "BGCV12"


def test_whitespace_and_blank_tokens_in_semicolon_fields_are_dropped(tmp_path):
    _write_fixture(
        tmp_path,
        patients=[_base_patient(patient_id="P1", sample_ids=" S1 ; ;S2", pools="")],
        pools=[],
        severities=[dict(sample_id="S1", n_cells=10), dict(sample_id="S2", n_cells=20)],
    )
    out = adapt_stephenson(tmp_path)
    assert sorted(out["capture_id"].tolist()) == ["S1", "S2"]
