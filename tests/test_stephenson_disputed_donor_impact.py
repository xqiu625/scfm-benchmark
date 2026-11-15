"""Fixture-based tests for `stephenson_disputed_donor_impact.py`, the
minimal real per-donor/per-pool impact table review:20261011_113512
§4 item 3 asked for: does GEX's independent `patient_id` field agree with
the sample_id-derived donor, which pool/day/Resample combinations are
actually observed per source, and do the disputed donors' cells actually
participate in Job 1's output (0 cells recorded explicitly when absent).
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
TESTS_DIR = Path(__file__).parent
sys.path.insert(0, str(CODE_DIR))
sys.path.insert(0, str(TESTS_DIR))

from stephenson_disputed_donor_impact import (  # noqa: E402
    bcr_donor_detail,
    gex_donor_pool_detail,
    job1_gex_cross_check,
    job1_output_detail,
    tcr_donor_detail,
)
from test_prepare_provisional_diagnostic_input import _build_fixture, _run_cli  # noqa: E402


def _write_minimal_job1_h5ad(path, rows):
    """rows: dicts with keys cell_id, cohort, donor_id, lineage, outcome.
    Plain string datasets (not the categorical Group encoding) --
    validate_job1_output.read_categorical_or_plain's non-Group branch
    (arr.dtype.kind in ('S','O')) handles this directly, and this is a
    real, valid (if minimal) modern-AnnData-style obs layout."""
    with h5py.File(path, "w") as f:
        obs = f.create_group("obs")
        obs.attrs["_index"] = "cell_id"
        obs.create_dataset("cell_id", data=np.array([r["cell_id"] for r in rows], dtype="S"))
        for field in ["cohort", "donor_id", "lineage", "outcome"]:
            obs.create_dataset(field, data=np.array([r[field] for r in rows], dtype="S"))


def _write_gex_fixture(path, cells):
    """cells: list of dicts with keys Site, sample_id, patient_id,
    Resample, Collection_Day, cell_id. Raw h5py writes reproducing the
    REAL deposited file's legacy AnnData 0.1.0 encoding (confirmed on
    HPCC 2026-10-11) -- not anndata.write_h5ad (modern encoding)."""
    fields = ["Site", "sample_id", "patient_id", "Resample", "Collection_Day"]
    with h5py.File(path, "w") as f:
        obs = f.create_group("obs")
        obs.attrs["_index"] = "covid_index"
        obs.create_dataset("covid_index", data=np.array([c["cell_id"] for c in cells], dtype="S"))
        cats_grp = obs.create_group("__categories")
        for field in fields:
            values = [c[field] for c in cells]
            uniq = sorted(set(values))
            code_of = {v: i for i, v in enumerate(uniq)}
            codes = np.array([code_of[v] for v in values], dtype="i2")
            obs.create_dataset(field, data=codes)
            cats_grp.create_dataset(field, data=np.array(uniq, dtype="S"))


def test_gex_donor_pool_detail_reports_per_pool_resample_day_combinations(tmp_path):
    """Mirrors the real CV0201 shape found on HPCC: Initial/D0 in one
    pool, Resample/D28 in another -- a genuine resample TCR never
    observed (coverage gap), now surfaced directly from GEX."""
    h5ad_path = tmp_path / "gex.h5ad"
    _write_gex_fixture(
        h5ad_path,
        [
            {
                "cell_id": "BGCV08_AAA-1",
                "Site": "Cambridge",
                "sample_id": "BGCV08_CV0201",
                "patient_id": "CV0201",
                "Resample": "Initial",
                "Collection_Day": "D0",
            },
            {
                "cell_id": "BGCV06_BBB-1",
                "Site": "Cambridge",
                "sample_id": "BGCV06_CV0201",
                "patient_id": "CV0201",
                "Resample": "Resample",
                "Collection_Day": "D28",
            },
            {
                "cell_id": "BGCV08_CCC-1",
                "Site": "Cambridge",
                "sample_id": "BGCV08_CV0999",
                "patient_id": "CV0999",
                "Resample": "Initial",
                "Collection_Day": "D0",
            },
        ],
    )
    result = gex_donor_pool_detail(h5ad_path, ["CV0201"])
    detail = result["per_donor"]["CV0201"]
    assert detail["n_cells_total"] == 2
    assert detail["n_cells_by_pool"] == {"BGCV08": 1, "BGCV06": 1}
    combos = {
        (c["pool"], c["resample_flag"], c["collection_day"])
        for c in detail["pool_resample_day_combinations"]
    }
    assert combos == {("BGCV08", "Initial", "D0"), ("BGCV06", "Resample", "D28")}
    # whole-pool count includes the OTHER donor too -- the F5 "whole pool
    # vs disputed donor's share" distinction.
    assert result["whole_pool_cell_counts_all_donors"]["BGCV08"] == 2
    assert result["patient_id_vs_sample_id_derived_donor_check"]["n_mismatches"] == 0


def test_gex_donor_pool_detail_detects_a_real_patient_id_sample_id_mismatch(tmp_path):
    """Mutation check: if sample_id's embedded donor and the independent
    patient_id field ever disagreed, this must be caught, not averaged
    away or silently trusted from one field alone."""
    h5ad_path = tmp_path / "gex_mismatch.h5ad"
    _write_gex_fixture(
        h5ad_path,
        [
            {
                "cell_id": "BGCV06_AAA-1",
                "Site": "Cambridge",
                "sample_id": "BGCV06_CV0201",
                "patient_id": "CV0999",
                "Resample": "Initial",
                "Collection_Day": "D0",
            },  # deliberately disagrees
        ],
    )
    result = gex_donor_pool_detail(h5ad_path, ["CV0201", "CV0999"])
    check = result["patient_id_vs_sample_id_derived_donor_check"]
    assert check["n_mismatches"] == 1
    assert check["mismatch_examples"][0]["sample_id_derived_donor"] == "CV0201"
    assert check["mismatch_examples"][0]["patient_id_field"] == "CV0999"


def _write_tcr_fixture(path, rows):
    fieldnames = ["CellID", "Sample", "patient_id", "Site", "Resample", "Collection_Day"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for row in rows:
            row.setdefault("CellID", f"{row['Sample']}_{row['patient_id']}cell-1")
            w.writerow(row)


def test_tcr_donor_detail_distinguishes_which_pool_is_covered(tmp_path):
    """Mirrors the real finding: TCR observes CV0201 in BGCV08 but has
    NO rows for it in BGCV06 -- a per-pool coverage gap, not a blanket
    presence/absence boolean."""
    tcr_path = tmp_path / "tcr.tsv"
    _write_tcr_fixture(
        tcr_path,
        [
            {
                "Sample": "BGCV08",
                "patient_id": "CV0201",
                "Site": "Cambridge",
                "Resample": "Initial",
                "Collection_Day": "D0",
            },
            {
                "Sample": "BGCV13",
                "patient_id": "CV0999",
                "Site": "Cambridge",
                "Resample": "Initial",
                "Collection_Day": "D0",
            },
        ],
    )
    result = tcr_donor_detail(tcr_path, ["CV0201", "CV0326"])
    assert len(result["CV0201"]) == 1
    assert result["CV0201"][0]["sample"] == "BGCV08"
    assert result["CV0201"][0]["n_unique_cell_ids"] == 1
    assert result["CV0326"] == []  # explicitly empty, not absent from the dict


def test_tcr_donor_detail_distinguishes_rows_from_unique_cells(tmp_path):
    """Two rows sharing the exact same combo but different CellIDs are 2
    unique cells; two rows sharing the same CellID (two chains of the
    same cell) are 1."""
    tcr_path = tmp_path / "tcr2.tsv"
    _write_tcr_fixture(
        tcr_path,
        [
            {
                "CellID": "BGCV08_AAA-1",
                "Sample": "BGCV08",
                "patient_id": "CV0201",
                "Site": "Cambridge",
                "Resample": "Initial",
                "Collection_Day": "D0",
            },
            {
                "CellID": "BGCV08_AAA-1",
                "Sample": "BGCV08",
                "patient_id": "CV0201",
                "Site": "Cambridge",
                "Resample": "Initial",
                "Collection_Day": "D0",
            },  # 2nd chain, same cell
            {
                "CellID": "BGCV08_BBB-1",
                "Sample": "BGCV08",
                "patient_id": "CV0201",
                "Site": "Cambridge",
                "Resample": "Initial",
                "Collection_Day": "D0",
            },  # different cell
        ],
    )
    result = tcr_donor_detail(tcr_path, ["CV0201"])
    assert len(result["CV0201"]) == 1
    assert result["CV0201"][0]["n_rows"] == 3
    assert result["CV0201"][0]["n_unique_cell_ids"] == 2


def _write_bcr_fixture(path, rows):
    with gzip.open(path, "wt", newline="") as f:
        w = csv.writer(f)
        w.writerow(["", "barcode", "is_cell", "patient_id"])
        for i, (barcode, patient_id) in enumerate(rows):
            w.writerow([i, barcode, "True", patient_id])


def test_bcr_donor_detail_separates_contig_rows_from_unique_cell_ids(tmp_path):
    """R3a: a real barcode with two chain rows (heavy+light) must be
    counted as 1 cell, not 2 -- n_contig_rows and n_unique_cell_ids must
    diverge in exactly this case."""
    bcr_path = tmp_path / "bcr.csv.gz"
    _write_bcr_fixture(
        bcr_path,
        [
            ("BGCV06_AAA-1", "CV0201"),
            ("BGCV06_AAA-1", "CV0201"),  # same cell, 2 chains
            ("BGCV08_BBB-1", "CV0201"),
            ("BGCV13_CCC-1", "CV0999"),
        ],
    )
    result = bcr_donor_detail(bcr_path, ["CV0201", "CV0326"])
    assert result["CV0201"]["BGCV06"] == {"n_contig_rows": 2, "n_unique_cell_ids": 1}
    assert result["CV0201"]["BGCV08"] == {"n_contig_rows": 1, "n_unique_cell_ids": 1}
    assert result["CV0326"] == {}


def test_job1_output_detail_records_zero_cells_explicitly_when_donor_absent(tmp_path):
    input_root = tmp_path / "input"
    out_h5ad = tmp_path / "out" / "prepared.h5ad"
    _build_fixture(input_root)
    result = _run_cli(
        [
            "--input-root",
            str(input_root),
            "--out",
            str(out_h5ad),
            "--n-per-donor",
            "3",
            "--seed",
            "0",
        ]
    )
    assert result.returncode == 0, result.stdout + result.stderr

    detail = job1_output_detail(out_h5ad, ["NOT_A_REAL_DONOR"], cohort="CohortA")
    assert detail["NOT_A_REAL_DONOR"] == {"n_cells_in_job1_output": 0}


def test_job1_output_detail_reports_real_selected_cells(tmp_path):
    input_root = tmp_path / "input"
    out_h5ad = tmp_path / "out" / "prepared.h5ad"
    _build_fixture(input_root)
    result = _run_cli(
        [
            "--input-root",
            str(input_root),
            "--out",
            str(out_h5ad),
            "--n-per-donor",
            "3",
            "--seed",
            "0",
        ]
    )
    assert result.returncode == 0, result.stdout + result.stderr

    detail = job1_output_detail(out_h5ad, ["D00"], cohort="CohortA")
    entry = detail["D00"]
    assert entry["n_cells_in_job1_output"] > 0
    assert set(entry["lineages"]) <= {"Platelet", "Monocyte_Myeloid"}
    assert entry["outcome_eligibility"]["column_loaded_and_inspected"] is True
    assert len(entry["example_cell_ids"]) <= 5


def test_job1_output_detail_outcome_eligibility_excludes_unknown_not_just_null(tmp_path):
    """R3b, the real defect: the prior version treated ANY non-null value
    (including the literal string "unknown", which prepare_provisional_
    diagnostic_input.py's own validate_event_column writes for missing
    events) as "known". The real run_batch_diagnostics.py bio_known_mask
    excludes --exclude-outcome-values (default ["unknown"]), not merely
    nulls -- this must match that definition."""
    h5ad_path = tmp_path / "job1_outcome_fixture.h5ad"
    _write_minimal_job1_h5ad(
        h5ad_path,
        [
            {
                "cell_id": "BGCV06_AAA-1",
                "cohort": "Stephenson 2021",
                "donor_id": "CV0201",
                "lineage": "Platelet",
                "outcome": "fatal",
            },
            {
                "cell_id": "BGCV06_BBB-1",
                "cohort": "Stephenson 2021",
                "donor_id": "CV0201",
                "lineage": "Platelet",
                "outcome": "unknown",
            },  # non-null, but NOT known
            {
                "cell_id": "BGCV08_CCC-1",
                "cohort": "Stephenson 2021",
                "donor_id": "CV0201",
                "lineage": "Monocyte_Myeloid",
                "outcome": "survived",
            },
        ],
    )
    entry = job1_output_detail(h5ad_path, ["CV0201"], cohort="Stephenson 2021")["CV0201"]
    assert entry["n_cells_in_job1_output"] == 3
    elig = entry["outcome_eligibility"]
    assert elig["n_known_per_bio_known_mask_definition"] == 2  # "unknown" excluded
    assert elig["n_excluded_or_null"] == 1


def test_job1_gex_cross_check_reports_matched_unmatched_and_donor_mismatches(tmp_path):
    h5ad_path = tmp_path / "job1_join_fixture.h5ad"
    _write_minimal_job1_h5ad(
        h5ad_path,
        [
            {
                "cell_id": "BGCV06_AAA-1",
                "cohort": "Stephenson 2021",
                "donor_id": "CV0201",
                "lineage": "Platelet",
                "outcome": "fatal",
            },
            {
                "cell_id": "BGCV08_BBB-1",
                "cohort": "Stephenson 2021",
                "donor_id": "CV0201",
                "lineage": "Monocyte_Myeloid",
                "outcome": "survived",
            },
            {
                "cell_id": "BGCV06_CCC-1",
                "cohort": "Stephenson 2021",
                "donor_id": "CV0201",
                "lineage": "Platelet",
                "outcome": "unknown",
            },  # not in gex_cell_donor -> unmatched
        ],
    )
    gex_cell_donor = {
        "BGCV06_AAA-1": "CV0201",
        "BGCV08_BBB-1": "CV0999",  # deliberately disagrees with Job 1's own label -> mismatch
    }
    result = job1_gex_cross_check(h5ad_path, gex_cell_donor, ["CV0201"])
    entry = result["CV0201"]
    assert entry["n_cells_total"] == 3
    assert entry["n_matched_to_gex_cell_id"] == 2
    assert entry["n_unmatched_to_gex_cell_id"] == 1
    assert entry["n_donor_label_mismatches"] == 1
    assert entry["donor_label_mismatch_examples"][0]["cell_id"] == "BGCV08_BBB-1"
    assert entry["n_cells_by_pool_and_lineage"] == {
        "BGCV06::Platelet": 2,
        "BGCV08::Monocyte_Myeloid": 1,
    }


def _minimal_cli_inputs(tmp_path):
    gex_path = tmp_path / "gex.h5ad"
    _write_gex_fixture(
        gex_path,
        [
            {
                "cell_id": "BGCV06_AAA-1",
                "Site": "Cambridge",
                "sample_id": "BGCV06_CV0201",
                "patient_id": "CV0201",
                "Resample": "Initial",
                "Collection_Day": "D0",
            },
        ],
    )
    tcr_path = tmp_path / "tcr.tsv"
    _write_tcr_fixture(tcr_path, [])
    bcr_path = tmp_path / "bcr.csv.gz"
    _write_bcr_fixture(bcr_path, [])
    job1_path = tmp_path / "job1.h5ad"
    _write_minimal_job1_h5ad(
        job1_path,
        [
            {
                "cell_id": "BGCV06_AAA-1",
                "cohort": "Stephenson 2021",
                "donor_id": "CV0201",
                "lineage": "Platelet",
                "outcome": "fatal",
            },
        ],
    )
    return gex_path, job1_path, tcr_path, bcr_path


def test_cli_writes_source_manifest_with_full_provenance(tmp_path):
    """R5: the prior version had no run provenance at all."""
    gex_path, job1_path, tcr_path, bcr_path = _minimal_cli_inputs(tmp_path)
    out_dir = tmp_path / "out"
    result = subprocess.run(
        [
            sys.executable,
            str(CODE_DIR / "stephenson_disputed_donor_impact.py"),
            "--gex-h5ad",
            str(gex_path),
            "--job1-h5ad",
            str(job1_path),
            "--tcr-tsv",
            str(tcr_path),
            "--bcr-gz",
            str(bcr_path),
            "--donors",
            "CV0201",
            "--out-dir",
            str(out_dir),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    manifest = json.loads((out_dir / "source_manifest.json").read_text())
    assert manifest["script_sha256"]
    assert len(manifest["dependency_modules"]) == 2
    assert manifest["run_ownership"]["pid"]
    assert (out_dir / ".run_owner.json").exists()
    assert (out_dir / "disputed_donor_impact.json").exists()


def test_claim_ownership_is_atomic_for_two_truly_simultaneous_starts(tmp_path):
    """R1 (review:20261011_124253): the real guarantee the fix
    relies on. Calling `claim_ownership` twice on the SAME fresh lock path
    (modeling two processes that both passed the outer non-empty check
    before either wrote anything) must let exactly one through."""
    sys.path.insert(0, str(Path(__file__).parent.parent / "scripts/00_shared/code"))
    from run_ownership import OwnershipError, claim_ownership  # noqa: E402

    lock_path = tmp_path / "race" / ".run_owner.json"
    lock_path.parent.mkdir()
    a_owner = claim_ownership(lock_path, command="process-A")
    with pytest.raises(OwnershipError):
        claim_ownership(lock_path, command="process-B")
    assert json.loads(lock_path.read_text()) == a_owner


def test_cli_refuses_a_dir_already_claimed_by_another_run(tmp_path):
    """CLI-level companion: "A" claims ownership of the out-dir first
    (modeling a real concurrent process that has claimed but not yet
    written), then "B" (the real CLI) must fail -- via the outer
    non-empty-dir check in this exact scenario, since A's lock file
    already makes the dir non-empty -- and A's lock file must be
    untouched by B's attempt."""
    gex_path, job1_path, tcr_path, bcr_path = _minimal_cli_inputs(tmp_path)
    out_dir = tmp_path / "out_race"
    out_dir.mkdir()
    sys.path.insert(0, str(Path(__file__).parent.parent / "scripts/00_shared/code"))
    from run_ownership import claim_ownership  # noqa: E402

    claim_ownership(out_dir / ".run_owner.json", command="process-A (simulated concurrent run)")
    before = (out_dir / ".run_owner.json").read_text()

    result = subprocess.run(
        [
            sys.executable,
            str(CODE_DIR / "stephenson_disputed_donor_impact.py"),
            "--gex-h5ad",
            str(gex_path),
            "--job1-h5ad",
            str(job1_path),
            "--tcr-tsv",
            str(tcr_path),
            "--bcr-gz",
            str(bcr_path),
            "--donors",
            "CV0201",
            "--out-dir",
            str(out_dir),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "already exists and is non-empty" in (result.stdout + result.stderr)
    after = (out_dir / ".run_owner.json").read_text()
    assert after == before
    assert not (out_dir / "disputed_donor_impact.json").exists()
