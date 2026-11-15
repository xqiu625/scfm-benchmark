"""Fixture-based tests for `audit_stephenson_source_identity.py`'s real
parsing/classification/join functions (not mocks). Covers the 5 defects
`review:20261010_202340` found in v1 (F1-F5) plus the v1 mutation
scenarios that still apply.
"""

import csv
import gzip
import json
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

CODE_DIR = Path(__file__).parent.parent / "scripts/33_batch_effect_audit/code"
sys.path.insert(0, str(CODE_DIR))

from audit_stephenson_source_identity import (  # noqa: E402
    build_dependency_components,
    build_identity_conflicts,
    classify_multi_pool_donor,
    classify_pool_evidence,
    cross_source_cell_join,
    diff_component_mappings,
    exclude_donors_from_pools,
    parse_annotation_cellids,
    parse_bcr_pool_membership,
    parse_canonical_donor_pools,
    parse_gex_cambridge,
    parse_tcr_cambridge_and_sanger,
)

SDRF_HEADER = [
    "Source Name",
    "Characteristics[organism]",
    "Characteristics[individual]",
    "Characteristics[sequence pool]",
]


def _write_sdrf(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(SDRF_HEADER)
        for donor, pool in rows:
            w.writerow([f"{pool}_{donor}", "Homo sapiens", donor, pool])


def _write_bcr_gz(path, rows):
    with gzip.open(path, "wt", newline="") as f:
        w = csv.writer(f)
        w.writerow(["", "barcode", "is_cell", "patient_id"])
        for i, (barcode, patient_id) in enumerate(rows):
            w.writerow([i, barcode, "True", patient_id])


def _write_tcr(path, rows):
    fieldnames = ["CellID", "Sample", "patient_id", "Site", "Resample", "Collection_Day"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for row in rows:
            row.setdefault("CellID", f"{row['Sample']}_{row['patient_id']}cell-1")
            w.writerow(row)


# ---- F1: multi-pool donor classification -----------------------------------


def test_classify_multi_pool_donor_same_day_is_not_a_resample():
    entries = [("BGCV01", "Initial", "D0"), ("BGCV02", "Initial", "D0")]
    assert classify_multi_pool_donor(entries) == "same_day_multi_pool"


def test_classify_multi_pool_donor_d0_d28_is_a_resample():
    entries = [("BGCV02", "Initial", "D0"), ("BGCV12", "Resample", "D28")]
    assert classify_multi_pool_donor(entries) == "d0_d28_style_resample"


def test_classify_multi_pool_donor_inconsistent_combo_is_named_not_bucketed():
    """Two pools, same day, but DIFFERENT flags -- does not fit either
    clean pattern; must not be silently called a resample."""
    entries = [("BGCV01", "Initial", "D0"), ("BGCV02", "Resample", "D0")]
    assert classify_multi_pool_donor(entries) == "inconsistent_day_flag_combination"


def test_classify_multi_pool_donor_single_pool_is_not_flagged():
    assert classify_multi_pool_donor([("BGCV01", "Initial", "D0")]) == "single_pool"


def test_classify_multi_pool_donor_inverted_day_flag_pairing_is_not_a_resample():
    """F5a: two pools, flags {Initial,Resample}, two distinct days -- but
    paired backwards (Initial is the LATER day). The pre-fix version only
    checked the two SETS and would have called this a resample too."""
    entries = [("BGCV01", "Initial", "D28"), ("BGCV02", "Resample", "D0")]
    assert classify_multi_pool_donor(entries) == "inconsistent_day_flag_combination"


def test_classify_multi_pool_donor_handles_non_d_prefixed_or_multi_digit_days():
    """Day-string sort must be numeric, not lexical ('D14' lexically
    precedes 'D7')."""
    entries = [("BGCV01", "Initial", "D7"), ("BGCV02", "Resample", "D14")]
    assert classify_multi_pool_donor(entries) == "d0_d28_style_resample"


# ---- F2: declared vs observed conflict semantics ----------------------------


def test_two_declared_sources_disagreeing_is_a_real_conflict():
    result = classify_pool_evidence(
        "P1",
        {
            "sdrf": {"D1", "D2"},
            "local_pool_composition": {"D1", "D3"},
            "gex_processed": None,
            "bcr": None,
            "tcr_cambridge": None,
        },
    )
    assert result["conflict_type"] == "declared_membership_disagreement"


def test_observed_subset_of_declared_is_coverage_gap_not_conflict():
    """The exact static counterexample automated review named: SDRF declares {D1,D2},
    TCR observes only {D1} -- this must NOT be reported as a conflict."""
    result = classify_pool_evidence(
        "P1",
        {
            "sdrf": {"D1", "D2"},
            "local_pool_composition": None,
            "gex_processed": None,
            "bcr": None,
            "tcr_cambridge": {"D1"},
        },
    )
    assert result["conflict_type"] == "none"
    assert "subset_of_sdrf" in result["observed_coverage"]["tcr_cambridge"]
    assert "missing=['D2']" in result["observed_coverage"]["tcr_cambridge"]


def test_observed_member_not_in_any_declared_is_a_real_contradiction():
    result = classify_pool_evidence(
        "P1",
        {
            "sdrf": {"D1", "D2"},
            "local_pool_composition": None,
            "gex_processed": None,
            "bcr": {"D1", "D9"},
            "tcr_cambridge": None,  # D9 is nowhere in declared
        },
    )
    assert result["conflict_type"] == "observed_contradicts_declared"
    assert result["observed_contradictions"]["bcr"] == ["D9"]


def test_pool_with_no_declared_reference_is_reported_not_silently_passed():
    result = classify_pool_evidence(
        "P1",
        {
            "sdrf": None,
            "local_pool_composition": None,
            "gex_processed": None,
            "bcr": {"D1"},
            "tcr_cambridge": None,
        },
    )
    assert result["conflict_type"] == "insufficient_declared_reference"


def test_build_identity_conflicts_end_to_end_matches_real_bgcv06_pattern():
    """Reproduces the real BGCV06 shape: SDRF and local disagree (declared
    conflict); TCR's subset doesn't add a second conflict type."""
    sources = {
        "sdrf": {"BGCV06": {"CV0037", "CV0178", "CV0234", "CV0326"}},
        "local_pool_composition": {"BGCV06": {"CV0037", "CV0178", "CV0201", "CV0234"}},
        "gex_processed": {"BGCV06": {"CV0037", "CV0178", "CV0201", "CV0234"}},
        "bcr": {"BGCV06": {"CV0037", "CV0178", "CV0201", "CV0234"}},
        "tcr_cambridge": {"BGCV06": {"CV0037", "CV0178", "CV0234"}},  # subset, no CV0201/CV0326
    }
    rows = build_identity_conflicts(sources)
    row = next(r for r in rows if r["pool"] == "BGCV06")
    assert row["conflict_type"] == "declared_membership_disagreement"
    coverage = json.loads(row["observed_coverage"])
    assert "subset_of" in coverage["tcr_cambridge"]


# ---- F3: GEX / annotation / BCR / TCR per-cell join -------------------------


def _write_minimal_gex_h5ad(path, cells):
    """cells: list of (cell_id, site, sample_id_label). Builds the LEGACY
    AnnData 0.1.0 encoding (__categories sibling group) this real deposit
    uses, via raw h5py writes -- not via anndata.write_h5ad (which would
    produce the MODERN encoding this file does not actually use)."""
    with h5py.File(path, "w") as f:
        obs = f.create_group("obs")
        obs.attrs["_index"] = "covid_index"
        ids = np.array([c[0] for c in cells], dtype="S")
        obs.create_dataset("covid_index", data=ids)
        cats_grp = obs.create_group("__categories")
        for field_idx, field in enumerate(["Site", "sample_id"]):
            values = [c[1 + field_idx] for c in cells]
            uniq = sorted(set(values))
            code_of = {v: i for i, v in enumerate(uniq)}
            codes = np.array([code_of[v] for v in values], dtype="i2")
            obs.create_dataset(field, data=codes)
            cats_grp.create_dataset(field, data=np.array(uniq, dtype="S"))


def test_parse_gex_cambridge_reads_legacy_categorical_encoding(tmp_path):
    h5ad_path = tmp_path / "fixture_gex.h5ad"
    _write_minimal_gex_h5ad(
        h5ad_path,
        [
            ("BGCV06_AAA-1", "Cambridge", "BGCV06_CV0201"),
            ("BGCV06_BBB-1", "Cambridge", "BGCV06_CV0037"),
            ("MH123_CCC-1", "Ncl", "MH123"),  # non-Cambridge, no "_" donor split -> excluded
        ],
    )
    pools, cell_counts, cell_donor, conflicts, coverage = parse_gex_cambridge(h5ad_path)
    assert pools == {"BGCV06": {"CV0201", "CV0037"}}
    assert cell_counts == {"BGCV06": 2}
    assert cell_donor == {"BGCV06_AAA-1": "CV0201", "BGCV06_BBB-1": "CV0037"}
    assert conflicts == []
    assert coverage == {
        "n_missing_site_or_cell_id": 0,
        "n_missing_sample_id_label": 0,
        "n_cell_id_conflicts": 0,
    }


# ---- F2a: legacy -1 missing-code sentinel must not become the last category --


def test_legacy_categorical_missing_code_becomes_none_not_last_category(tmp_path):
    """The real bug: `cats[codes]` with `codes == -1` silently indexed the
    LAST category. A missing Site/sample_id code must decode to None, and
    a Cambridge-site cell with a missing sample_id label must be excluded
    (not crash, not get the wrong donor)."""
    h5ad_path = tmp_path / "fixture_missing_code.h5ad"
    with h5py.File(h5ad_path, "w") as f:
        obs = f.create_group("obs")
        obs.attrs["_index"] = "covid_index"
        obs.create_dataset("covid_index", data=np.array(["c1", "c2"], dtype="S"))
        cats_grp = obs.create_group("__categories")
        # Site: both real cells are Cambridge (code 0); "Zzz_not_cambridge"
        # deliberately sorts LAST so a -1-as-cats[-1] bug would misroute.
        obs.create_dataset("Site", data=np.array([0, 0], dtype="i2"))
        cats_grp.create_dataset(
            "Site", data=np.array(["Cambridge", "Zzz_not_cambridge"], dtype="S")
        )
        # sample_id: c1 has a REAL code (0); c2 has the missing sentinel (-1).
        obs.create_dataset("sample_id", data=np.array([0, -1], dtype="i2"))
        cats_grp.create_dataset(
            "sample_id", data=np.array(["BGCV06_CV0201", "Zzz_last_label"], dtype="S")
        )

    pools, cell_counts, cell_donor, conflicts, coverage = parse_gex_cambridge(h5ad_path)
    assert pools == {
        "BGCV06": {"CV0201"}
    }  # c2's missing label must not contribute "Zzz_last_label"
    assert cell_donor == {"c1": "CV0201"}
    assert "c2" not in cell_donor
    assert coverage["n_missing_sample_id_label"] == 1


def test_legacy_categorical_out_of_range_code_raises_not_wraps(tmp_path):
    h5ad_path = tmp_path / "fixture_bad_code.h5ad"
    with h5py.File(h5ad_path, "w") as f:
        obs = f.create_group("obs")
        obs.attrs["_index"] = "covid_index"
        obs.create_dataset("covid_index", data=np.array(["c1"], dtype="S"))
        cats_grp = obs.create_group("__categories")
        obs.create_dataset("Site", data=np.array([5], dtype="i2"))  # out of range, not -1
        cats_grp.create_dataset("Site", data=np.array(["Cambridge"], dtype="S"))
        obs.create_dataset("sample_id", data=np.array([0], dtype="i2"))
        cats_grp.create_dataset("sample_id", data=np.array(["BGCV06_CV0201"], dtype="S"))

    with pytest.raises(ValueError):
        parse_gex_cambridge(h5ad_path)


# ---- F2c: a repeated CellID must not silently overwrite a prior donor ------


def test_parse_gex_cambridge_duplicate_cell_id_different_donor_is_a_conflict(tmp_path):
    h5ad_path = tmp_path / "fixture_dup.h5ad"
    _write_minimal_gex_h5ad(
        h5ad_path,
        [
            ("BGCV06_AAA-1", "Cambridge", "BGCV06_CV0201"),
            ("BGCV06_AAA-1", "Cambridge", "BGCV06_CV0037"),  # same CellID, different donor
        ],
    )
    pools, cell_counts, cell_donor, conflicts, coverage = parse_gex_cambridge(h5ad_path)
    assert conflicts == [("BGCV06_AAA-1", {"CV0201", "CV0037"})]
    assert "BGCV06_AAA-1" not in cell_donor  # excluded, not silently picked
    assert coverage["n_cell_id_conflicts"] == 1


# ---- F2b: GEX is OBSERVED, not DECLARED -- a GEX subset must not become --
# ---- a declared_membership_disagreement (the same coverage-gap rule BCR/TCR get)


def test_gex_subset_of_declared_roster_is_not_a_declared_conflict():
    result = classify_pool_evidence(
        "P1",
        {
            "sdrf": {"D1", "D2"},
            "local_pool_composition": {"D1", "D2"},
            "gex_processed": {"D1"},  # GEX only observed D1 -- a coverage gap, not a conflict
            "bcr": None,
            "tcr_cambridge": None,
        },
    )
    assert result["conflict_type"] == "none"
    assert "subset_of" in result["observed_coverage"]["gex_processed"]


def test_gex_processed_is_in_observed_not_declared_sources():
    import audit_stephenson_source_identity as mod

    assert "gex_processed" not in mod.DECLARED_SOURCES
    assert "gex_processed" in mod.OBSERVED_SOURCES


def test_parse_annotation_cellids_streams_first_column_only(tmp_path):
    csv_path = tmp_path / "annotation.csv"
    csv_path.write_text(",full_clustering\nBGCV06_AAA-1,MAIT\nBGCV06_BBB-1,B_memory\n")
    ids = parse_annotation_cellids(csv_path)
    assert ids == {"BGCV06_AAA-1", "BGCV06_BBB-1"}


def test_cross_source_cell_join_detects_real_per_cell_agreement_and_disagreement():
    gex_cell_donor = {"c1": "CV0201", "c2": "CV0037", "c3": "CV0178"}
    bcr_cell_donor = {"c1": "CV0201", "c2": "CV0999"}  # c2 disagrees with GEX
    tcr_cell_donor = {"c3": "CV0178"}
    annotation_ids = {"c1", "c2"}  # c3 deliberately absent from annotation
    result = cross_source_cell_join(gex_cell_donor, bcr_cell_donor, tcr_cell_donor, annotation_ids)
    assert result["bcr"]["n_shared_cell_ids"] == 2
    assert result["bcr"]["n_agree"] == 1
    assert result["bcr"]["n_disagree"] == 1
    assert result["bcr"]["disagree_examples"][0]["cell_id"] == "c2"
    assert result["tcr_cambridge"]["n_agree"] == 1
    assert result["annotation_coverage"]["n_gex_cellids_found_in_annotation"] == 2  # c1, c2
    assert result["annotation_coverage"]["n_gex_cellids_missing_from_annotation"] == 1  # c3


def test_cross_source_cell_join_combined_overlap_is_not_a_naive_sum():
    """F1: the prior decision sheet added n_bcr_shared + n_tcr_shared as
    if the two intersections with GEX were disjoint. Here c2 is shared
    with BOTH BCR and TCR, so the naive sum (2+2=4) would overcount the
    real distinct-cell union (3)."""
    gex_cell_donor = {"c1": "D1", "c2": "D2", "c3": "D3"}
    bcr_cell_donor = {"c1": "D1", "c2": "D2"}
    tcr_cell_donor = {"c2": "D2", "c3": "D3"}
    result = cross_source_cell_join(gex_cell_donor, bcr_cell_donor, tcr_cell_donor, set())
    combined = result["combined"]
    assert combined["n_bcr_shared_with_gex"] == 2
    assert combined["n_tcr_shared_with_gex"] == 2
    assert combined["n_overlap_between_bcr_and_tcr_shared_sets"] == 1  # c2
    assert combined["n_union_shared_distinct_cells"] == 3  # c1, c2, c3 -- NOT 2+2=4
    assert combined["n_union_agree_with_every_shared_observed_source"] == 3
    assert combined["n_union_disagree_with_at_least_one_shared_observed_source"] == 0


def test_cross_source_cell_join_combined_flags_disagreement_with_either_source():
    gex_cell_donor = {"c1": "D1"}
    bcr_cell_donor = {"c1": "D1"}
    tcr_cell_donor = {"c1": "D_WRONG"}  # disagrees with GEX
    result = cross_source_cell_join(gex_cell_donor, bcr_cell_donor, tcr_cell_donor, set())
    combined = result["combined"]
    assert combined["n_union_shared_distinct_cells"] == 1
    assert combined["n_union_agree_with_every_shared_observed_source"] == 0
    assert combined["n_union_disagree_with_at_least_one_shared_observed_source"] == 1


# ---- TCR per-cell conflicts (extends the v1 BCR-only check) -----------------


def test_tcr_detects_per_cell_identity_conflict_via_cellid():
    tcr_path_rows = [
        {
            "CellID": "BGCV01_AAA-1",
            "Sample": "BGCV01",
            "patient_id": "D1",
            "Site": "Cambridge",
            "Resample": "Initial",
            "Collection_Day": "D0",
        },
        {
            "CellID": "BGCV01_AAA-1",
            "Sample": "BGCV01",
            "patient_id": "D2",
            "Site": "Cambridge",
            "Resample": "Initial",
            "Collection_Day": "D0",
        },
    ]
    import tempfile

    with tempfile.TemporaryDirectory() as d:
        tcr_path = Path(d) / "tcr.tsv"
        _write_tcr(tcr_path, tcr_path_rows)
        _pools, _resample, _sanger, cell_conflicts, _cell_donor = parse_tcr_cambridge_and_sanger(
            tcr_path
        )
    assert len(cell_conflicts) == 1
    cell_id, ids = cell_conflicts[0]
    assert cell_id == "BGCV01_AAA-1"
    assert ids == {"D1", "D2"}


# ---- carried-over v1 coverage: bare-barcode merge, canonical "nan" sentinel --


def test_bcr_bare_barcode_reused_across_pools_is_not_merged(tmp_path):
    bcr_path = tmp_path / "fixture_bcr.csv.gz"
    _write_bcr_gz(
        bcr_path,
        [
            ("P1_AAAACCC-1", "D1"),
            ("P1_AAAACCC-1", "D1"),
            ("P2_AAAACCC-1", "D9"),
        ],
    )
    pool_donors, cell_conflicts, _cell_donor = parse_bcr_pool_membership(bcr_path)
    assert pool_donors == {"P1": {"D1"}, "P2": {"D9"}}
    assert cell_conflicts == []


def test_canonical_notes_nan_sentinel_is_not_an_answer(tmp_path):
    csv_path = tmp_path / "canonical_donor_labels.csv"
    fieldnames = ["cohort", "donor_id", "notes"]
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerow({"cohort": "Stephenson 2021", "donor_id": "AP1", "notes": "pools: nan"})
        w.writerow({"cohort": "Stephenson 2021", "donor_id": "CV0037", "notes": "pools: BGCV06"})
    donor_pools = parse_canonical_donor_pools(csv_path)
    assert "AP1" not in donor_pools
    assert donor_pools["CV0037"] == {"BGCV06"}


# ---- item 3: two-mapping dependency component diff -------------------------


def test_build_dependency_components_transitive_closure():
    pools = {"P1": {"D1", "D2"}, "P2": {"D2", "D3"}, "P3": {"D9"}}
    comps = build_dependency_components(pools)
    members = sorted(frozenset(v) for v in comps.values())
    assert frozenset({"D1", "D2", "D3"}) in members  # P1+P2 merged via shared D2
    assert frozenset({"D9"}) in members


def test_exclude_donors_from_pools_removes_from_every_pool():
    pools = {"P1": {"D1", "D2"}, "P2": {"D2", "D3"}}
    result = exclude_donors_from_pools(pools, {"D2"})
    assert result == {"P1": {"D1"}, "P2": {"D3"}}


def test_excluding_disputed_donors_does_not_change_an_unrelated_donors_component():
    """Decision sheet option 3, actually tested: removing the disputed
    donors from the pool graph must not perturb a donor in an unrelated
    pool that never shared a pool with them."""
    pools = {
        "BGCV06": {"D_shared6", "D326"},
        "BGCV08": {"D_shared8", "D201"},
        "BGCV13": {"D_shared8", "D_shared13"},
        "BGCV99_unrelated": {"D_other1", "D_other2"},
    }
    excluded = exclude_donors_from_pools(pools, {"D201", "D326"})
    result = diff_component_mappings(pools, excluded, "with_disputed", "excluding_disputed")
    changed = {r["donor_id"] for r in result["donor_rows"] if r["component_changed"]}
    assert "D_other1" not in changed and "D_other2" not in changed
    assert {"D201", "D326"} <= changed  # they trivially change (removed entirely)


def test_diff_component_mappings_flags_only_donors_whose_component_actually_changes():
    """Mirrors the real BGCV06/BGCV13 shape in miniature: swapping D201/
    D326 between BGCV06 and BGCV13 changes which cluster BGCV06 joins (it
    has no OTHER shared donor with the BGCV08/13 cluster), so D_shared6
    and D_shared8/D_shared13 legitimately change too -- but a donor in a
    completely unrelated pool must NOT be flagged."""
    pools_a = {  # "sdrf"-style
        "BGCV06": {"D_shared6", "D326"},
        "BGCV08": {"D_shared8", "D201"},
        "BGCV13": {"D_shared8", "D_shared13"},
        "BGCV99_unrelated": {"D_other1", "D_other2"},
    }
    pools_b = {  # "local"-style: D201/D326 swapped
        "BGCV06": {"D_shared6", "D201"},
        "BGCV08": {"D_shared8", "D201"},
        "BGCV13": {"D_shared8", "D326", "D_shared13"},
        "BGCV99_unrelated": {"D_other1", "D_other2"},
    }
    result = diff_component_mappings(pools_a, pools_b, "a", "b")
    changed = {r["donor_id"] for r in result["donor_rows"] if r["component_changed"]}
    assert {"D201", "D326", "D_shared6"} <= changed
    assert "D_other1" not in changed and "D_other2" not in changed


# ---- F5: no-clobber --out-dir, via the real CLI entry point ----------------


def test_cli_refuses_to_overwrite_a_nonempty_out_dir(tmp_path):
    sdrf_path = tmp_path / "f.sdrf.txt"
    _write_sdrf(sdrf_path, [("D1", "P1")])
    local_csv = tmp_path / "local.csv"
    local_csv.write_text("pool,n_donors,donors\nP1,1,D1\n")
    canonical_csv = tmp_path / "canonical.csv"
    canonical_csv.write_text("cohort,donor_id,notes\n")
    gex_path = tmp_path / "gex.h5ad"
    _write_minimal_gex_h5ad(gex_path, [("P1_AAA-1", "Cambridge", "P1_D1")])
    annotation_csv = tmp_path / "annotation.csv"
    annotation_csv.write_text(",full_clustering\n")
    bcr_path = tmp_path / "bcr.csv.gz"
    _write_bcr_gz(bcr_path, [])
    tcr_path = tmp_path / "tcr.tsv"
    _write_tcr(tcr_path, [])

    out_dir = tmp_path / "out"
    out_dir.mkdir()
    (out_dir / "preexisting.txt").write_text("already here")

    result = subprocess.run(
        [
            sys.executable,
            str(CODE_DIR / "audit_stephenson_source_identity.py"),
            "--sdrf",
            str(sdrf_path),
            "--local-pool-composition",
            str(local_csv),
            "--canonical-donor-labels",
            str(canonical_csv),
            "--gex-h5ad",
            str(gex_path),
            "--annotation-csv",
            str(annotation_csv),
            "--bcr-gz",
            str(bcr_path),
            "--tcr-tsv",
            str(tcr_path),
            "--out-dir",
            str(out_dir),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "already exists" in (result.stdout + result.stderr)


def test_cli_allows_reusing_a_preexisting_but_empty_out_dir(tmp_path):
    """F5b: the exclusive-mkdir fix must not regress the legitimate case
    of an empty --out-dir that already exists (e.g. created by a prior
    aborted run, or by the caller itself)."""
    sdrf_path = tmp_path / "f.sdrf.txt"
    _write_sdrf(sdrf_path, [("D1", "P1")])
    local_csv = tmp_path / "local.csv"
    local_csv.write_text("pool,n_donors,donors\nP1,1,D1\n")
    canonical_csv = tmp_path / "canonical.csv"
    canonical_csv.write_text("cohort,donor_id,notes\n")
    gex_path = tmp_path / "gex.h5ad"
    _write_minimal_gex_h5ad(gex_path, [("P1_AAA-1", "Cambridge", "P1_D1")])
    annotation_csv = tmp_path / "annotation.csv"
    annotation_csv.write_text(",full_clustering\n")
    bcr_path = tmp_path / "bcr.csv.gz"
    _write_bcr_gz(bcr_path, [])
    tcr_path = tmp_path / "tcr.tsv"
    _write_tcr(tcr_path, [])

    out_dir = tmp_path / "out_empty"
    out_dir.mkdir()  # pre-exists, but empty

    result = subprocess.run(
        [
            sys.executable,
            str(CODE_DIR / "audit_stephenson_source_identity.py"),
            "--sdrf",
            str(sdrf_path),
            "--local-pool-composition",
            str(local_csv),
            "--canonical-donor-labels",
            str(canonical_csv),
            "--gex-h5ad",
            str(gex_path),
            "--annotation-csv",
            str(annotation_csv),
            "--bcr-gz",
            str(bcr_path),
            "--tcr-tsv",
            str(tcr_path),
            "--out-dir",
            str(out_dir),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert (out_dir / "summary.json").exists()


def test_claim_ownership_is_atomic_for_two_truly_simultaneous_starts(tmp_path):
    """R1 (review:20261011_124253): the real guarantee the fix
    relies on. `claim_ownership`'s `os.O_CREAT | os.O_EXCL` is an
    OS-atomic primitive -- calling it twice on the SAME fresh lock path
    (modeling two processes that both passed the outer non-empty check
    before either wrote anything) must let exactly one through."""
    sys.path.insert(0, str(Path(__file__).parent.parent / "scripts/00_shared/code"))
    from run_ownership import OwnershipError, claim_ownership  # noqa: E402

    lock_path = tmp_path / "race" / ".run_owner.json"
    lock_path.parent.mkdir()
    a_owner = claim_ownership(lock_path, command="process-A")
    with pytest.raises(OwnershipError):
        claim_ownership(lock_path, command="process-B")
    # A's own claim is untouched by B's failed attempt.
    assert json.loads(lock_path.read_text()) == a_owner


def test_cli_refuses_a_dir_already_claimed_by_another_run(tmp_path):
    """CLI-level companion to the direct atomicity test above: "A" claims
    ownership of the out-dir first (modeling a real concurrent process
    that has claimed but not yet written), then "B" (the real CLI) must
    fail -- via the outer non-empty-dir check in THIS exact scenario,
    since A's lock file already makes the dir non-empty by the time B's
    CLI starts -- and A's claimed lock file must be untouched by B's
    attempt. The underlying race for two TRULY simultaneous starts (both
    passing the outer check before either writes) is what the direct
    claim_ownership test above proves atomically."""
    sdrf_path = tmp_path / "f.sdrf.txt"
    _write_sdrf(sdrf_path, [("D1", "P1")])
    local_csv = tmp_path / "local.csv"
    local_csv.write_text("pool,n_donors,donors\nP1,1,D1\n")
    canonical_csv = tmp_path / "canonical.csv"
    canonical_csv.write_text("cohort,donor_id,notes\n")
    gex_path = tmp_path / "gex.h5ad"
    _write_minimal_gex_h5ad(gex_path, [("P1_AAA-1", "Cambridge", "P1_D1")])
    annotation_csv = tmp_path / "annotation.csv"
    annotation_csv.write_text(",full_clustering\n")
    bcr_path = tmp_path / "bcr.csv.gz"
    _write_bcr_gz(bcr_path, [])
    tcr_path = tmp_path / "tcr.tsv"
    _write_tcr(tcr_path, [])

    out_dir = tmp_path / "out_race"
    out_dir.mkdir()
    sys.path.insert(0, str(Path(__file__).parent.parent / "scripts/00_shared/code"))
    from run_ownership import claim_ownership  # noqa: E402

    claim_ownership(out_dir / ".run_owner.json", command="process-A (simulated concurrent run)")
    lock_path = out_dir / ".run_owner.json"
    before = lock_path.read_text()

    result = subprocess.run(
        [
            sys.executable,
            str(CODE_DIR / "audit_stephenson_source_identity.py"),
            "--sdrf",
            str(sdrf_path),
            "--local-pool-composition",
            str(local_csv),
            "--canonical-donor-labels",
            str(canonical_csv),
            "--gex-h5ad",
            str(gex_path),
            "--annotation-csv",
            str(annotation_csv),
            "--bcr-gz",
            str(bcr_path),
            "--tcr-tsv",
            str(tcr_path),
            "--out-dir",
            str(out_dir),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "already exists and is non-empty" in (result.stdout + result.stderr)
    after = lock_path.read_text()
    assert after == before  # B's failed attempt did not touch A's claimed lock
    assert not (out_dir / "summary.json").exists()  # B never got to write anything
