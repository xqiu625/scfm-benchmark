"""Tests for build_batch_metadata_inventory.py's and
resolve_stephenson_grouping.py's pure functions.

Targets the 2026-10-09 review's P1 finding: the v2 alias-normalization dict
comprehensions silently dropped one of several IDs sharing a normalized key.
v3's build_alias_crosswalk() must flag that as AMBIGUOUS, not drop anything.
Also covers the Stephenson union-find component builder's determinism and
cross-pool merging.
"""

import sys
from pathlib import Path

sys.path.insert(
    0,
    str(Path(__file__).parent.parent / "scripts/33_batch_effect_audit/code"),
)
from build_batch_metadata_inventory import _normalize_liu, build_alias_crosswalk
from resolve_stephenson_grouping import union_find_components


def test_crosswalk_unambiguous_normalized_match():
    canon = {"HGR0000042"}
    inv = {"HGR0000042_T0"}
    cw = build_alias_crosswalk("Liu 2021", canon, inv)
    assert len(cw) == 1
    row = cw.iloc[0]
    assert row["status"] == "normalized_match"
    assert row["cardinality"] == "1:1"
    assert row["canonical_donor_ids"] == "HGR0000042"
    assert row["capture_file_donor_ids"] == "HGR0000042_T0"


def test_crosswalk_includes_exact_matches_too():
    """Per the review: an ID that already matches exactly must still appear
    in the crosswalk (not only the previously-unmatched subset), so an
    EXTRA alias of an already-matched patient isn't invisible."""
    canon = {"HGR0000042"}
    inv = {"HGR0000042"}
    cw = build_alias_crosswalk("Liu 2021", canon, inv)
    assert len(cw) == 1
    assert cw.iloc[0]["status"] == "exact_match"


def test_crosswalk_flags_collision_instead_of_silently_picking_one():
    """The exact hazard the review found: two DIFFERENT canonical donors
    whose IDs both strip down to the same normalized key (e.g. two patients
    each resampled, '_T0' suffix). v2's dict comprehension would have kept
    only one, depending on set iteration order; v3 must flag it as
    AMBIGUOUS_NOT_RESOLVED and must not merge them."""
    # Construct a genuine normalized-key collision directly:
    assert _normalize_liu("HGR0000042_T0") == _normalize_liu("HGR0000042_T1") == "HGR0000042"
    canon_colliding = {"HGR0000042"}
    inv_colliding = {"HGR0000042_T0", "HGR0000042_T1"}  # two DISTINCT capture-file IDs, same key
    cw = build_alias_crosswalk("Liu 2021", canon_colliding, inv_colliding)
    assert len(cw) == 1
    row = cw.iloc[0]
    assert row["status"] == "AMBIGUOUS_NOT_RESOLVED"
    assert row["cardinality"] == "1:2"
    assert set(row["capture_file_donor_ids"].split(";")) == {"HGR0000042_T0", "HGR0000042_T1"}


def test_crosswalk_two_different_canonical_patients_colliding_is_also_ambiguous():
    """The other direction of the same hazard: two different BIOLOGICAL
    patients whose canonical IDs normalize to the same key must not be
    silently treated as the same person."""
    canon_colliding = {
        "HGR0000042_T0",
        "HGR0000042_T1",
    }  # pretend these are 2 distinct canonical IDs
    inv_colliding = {"HGR0000042"}
    cw = build_alias_crosswalk("Liu 2021", canon_colliding, inv_colliding)
    assert len(cw) == 1
    assert cw.iloc[0]["status"] == "AMBIGUOUS_NOT_RESOLVED"
    assert cw.iloc[0]["cardinality"] == "2:1"


def test_crosswalk_unmatched_canonical_only_gets_its_own_explicit_row():
    """v4: an ID with no counterpart on the other side at all still gets a
    row (CANONICAL_ONLY_UNMATCHED / CAPTURE_FILE_ONLY_UNMATCHED), not a
    silent drop -- a 2026-10-09 review wanted canonical-only and
    source-only records explicit in the crosswalk itself, not only in the
    caller's separate set-difference summary."""
    cw = build_alias_crosswalk("Liu 2021", {"HGR0000042"}, set())
    assert len(cw) == 1
    assert cw.iloc[0]["status"] == "CANONICAL_ONLY_UNMATCHED"
    assert cw.iloc[0]["cardinality"] == "1:0"


def test_crosswalk_unmatched_capture_file_only_gets_its_own_explicit_row():
    cw = build_alias_crosswalk("Liu 2021", set(), {"HGR0000042"})
    assert len(cw) == 1
    assert cw.iloc[0]["status"] == "CAPTURE_FILE_ONLY_UNMATCHED"
    assert cw.iloc[0]["cardinality"] == "0:1"


def test_crosswalk_unknown_cohort_returns_empty_fixed_schema():
    cw = build_alias_crosswalk("Not A Real Cohort", {"A"}, {"B"})
    assert len(cw) == 0
    assert list(cw.columns)


def test_union_find_merges_transitively_across_pools():
    """CV0902 appears in both BGCV01 and BGCV02 in the real data -- a donor
    shared between two pools must merge BOTH pools into one component."""
    edges = [
        ("A", "B"),
        ("A", "C"),
        ("A", "D"),  # pool 1: A,B,C,D
        ("E", "B"),
        ("E", "F"),
        ("E", "G"),  # pool 2: B,E,F,G -- B bridges pool 1 and 2
    ]
    nodes = {"A", "B", "C", "D", "E", "F", "G"}
    components = union_find_components(edges, nodes)
    assert len(set(components.values())) == 1  # everything merges into one component


def test_union_find_keeps_disjoint_pools_separate():
    edges = [("A", "B"), ("A", "C")]  # one pool
    nodes = {"A", "B", "C", "D"}  # D is unconnected
    components = union_find_components(edges, nodes)
    assert components["A"] == components["B"] == components["C"]
    assert components["D"] != components["A"]


def test_union_find_is_deterministic_regardless_of_edge_order():
    edges_order_1 = [("B", "A"), ("C", "A")]
    edges_order_2 = [("C", "A"), ("B", "A")]
    nodes = {"A", "B", "C"}
    c1 = union_find_components(edges_order_1, nodes)
    c2 = union_find_components(edges_order_2, nodes)
    assert c1 == c2
    assert c1["A"] == "A"  # lexicographically smallest wins as root
