#!/usr/bin/env python3
"""Minimal, real, per-donor/per-pool impact table for the disputed
Stephenson Cambridge donors (default: CV0201/CV0326, BGCV06/BGCV13, with
BGCV08 kept as associated evidence -- review:20261011_113512 §4
item 3). Bounded h5py reads only; never materializes GEX `X`/`layers` or
Job 1's `X`/`layers`. Does not infer any missing value and does not
change any real label -- only reports, per (donor, pool), what each real
source directly observed.

Separates two things the prior round's decision sheet conflated:
  - the WHOLE disputed pool's cell count (every donor in it), vs.
  - the disputed donor's OWN cell count within that pool.

Reads `obs['patient_id']` from the real GEX h5ad as an INDEPENDENT field
from `obs['sample_id']` (confirmed on HPCC 2026-10-11: these two fields
agree on all 130,644 real Cambridge cells, 0 mismatches) -- this answers
"does GEX's own patient_id agree with the sample_id-derived donor" with a
real, computed comparison rather than an assumption.

v2 (this file) responds to review:20261011_124253 (R1/R3/R5),
which reviewed v1 and found:

  R1: `--out-dir`'s check-then-mkdir was not atomic against two fresh
      concurrent starts. Now uses the project's existing
      `claim_ownership()` (same primitive `run_batch_diagnostics.py`
      already uses), on top of the pre-existing non-empty-dir check.
  R3a: `bcr_donor_detail` counted raw CSV ROWS per donor/pool and reported
      them as if they were cell counts -- a BCR contig file has multiple
      rows per cell (one per chain). Now reports `n_contig_rows` AND
      `n_unique_cell_ids` separately, never conflating the two.
  R3b: the docstring claimed "never reads outcome values"; the code
      actually read and decoded the full `outcome` column, just never
      printed individual values. The claim was corrected to describe what
      actually happens: the column IS loaded and inspected to compute an
      AGGREGATE known/unknown eligibility count, using the SAME
      `exclude_outcome_values` definition `run_batch_diagnostics.py`'s
      `bio_known_mask` uses (a value is "known" only if it is not in
      `--exclude-outcome-values`, default `["unknown"]` -- NOT merely
      non-null, since the string `"unknown"` is itself non-null).
      Individual outcome values are still never printed, returned, or
      used to select/filter any cell or change any label.
  R3c: added `job1_gex_cross_check` -- a real full-CellID join between
      Job 1's actually-selected cells for the disputed donors and the
      real GEX's own cell_donor map, reporting matched/unmatched counts,
      donor-label mismatches, and per (pool, lineage) cell counts -- the
      previous version only reported pool PREFIXES of selected cell_ids,
      never actually joined them against GEX.
  R5: added a `provenance` block (this script's own sha256, git revision
      and dirty status, the sha256 of the two modules it imports from,
      full argv, and the run-ownership record) -- the prior version had
      no provenance at all.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path

import h5py

CODE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE_DIR))
from audit_stephenson_source_identity import (  # noqa: E402
    _read_legacy_categorical, parse_gex_cambridge, write_source_manifest,
)
from validate_job1_output import read_categorical_or_plain  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "00_shared/code"))
from run_ownership import claim_ownership  # noqa: E402

DEFAULT_EXCLUDE_OUTCOME_VALUES = ["unknown"]

DEFAULT_DONORS = ["CV0201", "CV0326"]
DEFAULT_POOLS = ["BGCV06", "BGCV08", "BGCV13"]


def gex_donor_pool_detail(h5ad_path: Path, donors: list[str]) -> dict:
    """Reads Site, sample_id, patient_id, Resample, Collection_Day, and the
    obs index from the real GEX h5ad (legacy encoding) -- bounded to
    Cambridge-site cells only, never `X`/`layers`. Checks patient_id
    against the sample_id-derived donor as an independent cross-check,
    over ALL Cambridge cells (not just the disputed ones), and reports
    per-(donor, pool) cell counts and observed (Resample, Collection_Day)
    combinations for the requested donors specifically."""
    with h5py.File(h5ad_path, "r") as f:
        site = _read_legacy_categorical(f, "Site")
        sample_id = _read_legacy_categorical(f, "sample_id")
        patient_id = _read_legacy_categorical(f, "patient_id")
        resample = _read_legacy_categorical(f, "Resample")
        day = _read_legacy_categorical(f, "Collection_Day")
        index_key = f["obs"].attrs.get("_index", "covid_index")
        index_key = index_key.decode() if isinstance(index_key, bytes) else index_key
        cell_id = _read_legacy_categorical(f, index_key)

    cambridge = site == "Cambridge"
    sid_c = sample_id[cambridge]
    pid_c = patient_id[cambridge]

    n_compared = 0
    n_mismatch = 0
    mismatch_examples = []
    for sid, pid in zip(sid_c, pid_c):
        if sid is None or pid is None or "_" not in sid:
            continue
        derived = sid.split("_", 1)[1]
        n_compared += 1
        if derived != pid:
            n_mismatch += 1
            if len(mismatch_examples) < 10:
                mismatch_examples.append({"sample_id": sid, "sample_id_derived_donor": derived,
                                           "patient_id_field": pid})

    pool_whole_counts = defaultdict(int)
    for sid in sid_c:
        if sid and "_" in sid:
            pool_whole_counts[sid.split("_", 1)[0]] += 1

    per_donor = {}
    for donor in donors:
        mask = cambridge & (patient_id == donor)
        pools_present = sample_id[mask]
        rs = resample[mask]
        dy = day[mask]
        combos = defaultdict(int)
        for pool_label, r, d in zip(pools_present, rs, dy):
            pool = pool_label.split("_", 1)[0] if pool_label and "_" in pool_label else "UNKNOWN"
            combos[(pool, r, d)] += 1
        per_pool = defaultdict(int)
        for pool_label in pools_present:
            if pool_label and "_" in pool_label:
                per_pool[pool_label.split("_", 1)[0]] += 1
        per_donor[donor] = {
            "n_cells_total": int(mask.sum()),
            "n_cells_by_pool": dict(per_pool),
            "pool_resample_day_combinations": [
                {"pool": k[0], "resample_flag": k[1], "collection_day": k[2], "n_cells": v}
                for k, v in sorted(combos.items())
            ],
        }

    return {
        "patient_id_vs_sample_id_derived_donor_check": {
            "n_cambridge_cells_compared": n_compared,
            "n_mismatches": n_mismatch,
            "mismatch_examples": mismatch_examples,
        },
        "whole_pool_cell_counts_all_donors": dict(pool_whole_counts),
        "per_donor": per_donor,
    }


def tcr_donor_detail(tcr_tsv_path: Path, donors: list[str]) -> dict:
    """Streams the TCR TSV once, restricted to the requested donors --
    reports every (Site, Sample, Resample, Collection_Day) combination
    actually observed, rather than only a pool-level presence/absence
    boolean (which would hide that TCR observes a donor in ONE pool but
    not another). `n_rows` is honestly labeled as rows (a TCR file can
    also have multiple chain rows per cell); `n_unique_cell_ids` is the
    deduplicated-by-CellID count for the same combination."""
    donors_set = set(donors)
    per_donor_rows: dict[str, dict] = {d: defaultdict(int) for d in donors}
    per_donor_cellids: dict[str, dict] = {d: defaultdict(set) for d in donors}
    with open(tcr_tsv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            pid = row.get("patient_id")
            if pid not in donors_set:
                continue
            key = (row.get("Site"), row.get("Sample"), row.get("Resample"), row.get("Collection_Day"))
            per_donor_rows[pid][key] += 1
            cell_id = row.get("CellID")
            if cell_id:
                per_donor_cellids[pid][key].add(cell_id)
    return {
        donor: [
            {"site": k[0], "sample": k[1], "resample_flag": k[2], "collection_day": k[3],
             "n_rows": v, "n_unique_cell_ids": len(per_donor_cellids[donor][k])}
            for k, v in sorted(per_donor_rows[donor].items())
        ]
        for donor in donors
    }


def bcr_donor_detail(bcr_gz_path: Path, donors: list[str]) -> dict:
    """Streams the BCR gzip CSV once, restricted to the requested donors
    -- reports every pool (derived from the barcode prefix) actually
    observed, by row NAME lookup (never a hard-coded column position).

    R3a fix (review:20261011_124253): a BCR contig annotation
    file has MULTIPLE ROWS per real cell (one per chain -- heavy+light,
    etc.), so counting rows and calling the result "cells" overstates the
    real cell count. Now reports `n_contig_rows` (the raw row count) and
    `n_unique_cell_ids` (the real cell count, deduplicated by full
    barcode) separately -- never conflated into one number."""
    donors_set = set(donors)
    per_donor_rows: dict[str, dict] = {d: defaultdict(int) for d in donors}
    per_donor_cellids: dict[str, dict] = {d: defaultdict(set) for d in donors}
    with gzip.open(bcr_gz_path, "rt", newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        barcode_i = header.index("barcode")
        patient_i = header.index("patient_id")
        for row in reader:
            pid = row[patient_i]
            if pid not in donors_set:
                continue
            barcode = row[barcode_i]
            pool = barcode.split("_", 1)[0] if "_" in barcode else "UNKNOWN"
            per_donor_rows[pid][pool] += 1
            per_donor_cellids[pid][pool].add(barcode)
    return {
        donor: {
            pool: {
                "n_contig_rows": per_donor_rows[donor][pool],
                "n_unique_cell_ids": len(per_donor_cellids[donor][pool]),
            }
            for pool in per_donor_rows[donor]
        }
        for donor in donors
    }


def _read_job1_obs_arrays(job1_h5ad_path: Path):
    with h5py.File(job1_h5ad_path, "r") as f:
        obs = f["obs"]
        cohort_arr = read_categorical_or_plain(obs["cohort"])
        donor_arr = read_categorical_or_plain(obs["donor_id"])
        lineage_arr = read_categorical_or_plain(obs["lineage"])
        outcome_arr = read_categorical_or_plain(obs["outcome"]) if "outcome" in obs else None
        index_key = obs.attrs.get("_index", "cell_id")
        index_key = index_key.decode() if isinstance(index_key, bytes) else index_key
        cell_id_arr = read_categorical_or_plain(obs[index_key])
    return cohort_arr, donor_arr, lineage_arr, outcome_arr, cell_id_arr


def job1_output_detail(job1_h5ad_path: Path, donors: list[str], cohort: str = "Stephenson 2021",
                        exclude_outcome_values=DEFAULT_EXCLUDE_OUTCOME_VALUES) -> dict:
    """Reads Job 1's actual output obs (modern AnnData encoding, via
    `validate_job1_output.read_categorical_or_plain` -- bounded, never
    `X`/`layers`) and reports, per requested donor, whether any of its
    cells were actually selected into Job 1's output, their lineage
    composition, and their real cell_id values.

    R3b fix (review:20261011_124253): the prior docstring claimed
    this function never reads outcome VALUES, but it did read and decode
    the full `outcome` column -- it just never printed a per-cell value.
    The corrected claim: the column IS loaded and inspected, to compute an
    AGGREGATE known/unknown eligibility count using the SAME definition
    `run_batch_diagnostics.py`'s `bio_known_mask` uses (a value is "known"
    only if it is NOT in `exclude_outcome_values`, default `["unknown"]`
    -- not merely non-null, since the string `"unknown"` is itself
    non-null). No individual outcome value is ever printed, returned, or
    used to select/filter a cell or change any label."""
    exclude_set = {str(v) for v in exclude_outcome_values}
    cohort_arr, donor_arr, lineage_arr, outcome_arr, cell_id_arr = _read_job1_obs_arrays(job1_h5ad_path)

    cohort_mask = cohort_arr == cohort
    result = {}
    for donor in donors:
        mask = cohort_mask & (donor_arr == donor)
        n_cells = int(mask.sum())
        entry = {"n_cells_in_job1_output": n_cells}
        if n_cells:
            entry["lineages"] = sorted(set(lineage_arr[mask].tolist()))
            entry["example_cell_ids"] = cell_id_arr[mask][:5].tolist()
            entry["pool_prefixes_of_selected_cells"] = sorted({
                cid.split("_", 1)[0] for cid in cell_id_arr[mask].tolist() if cid and "_" in cid
            })
            if outcome_arr is not None:
                values = outcome_arr[mask]
                n_known = sum(1 for v in values if v is not None and str(v) not in exclude_set)
                entry["outcome_eligibility"] = {
                    "column_loaded_and_inspected": True,
                    "exclude_outcome_values_used": sorted(exclude_set),
                    "n_known_per_bio_known_mask_definition": int(n_known),
                    "n_excluded_or_null": int(n_cells - n_known),
                    "note": (
                        "Aggregate counts only -- individual outcome values are never "
                        "printed, returned, or used to select/filter any cell. 'known' "
                        "matches run_batch_diagnostics.py's bio_known_mask definition "
                        "(not in exclude_outcome_values), not merely non-null."
                    ),
                }
            else:
                entry["outcome_eligibility"] = "NOT_AVAILABLE (no outcome column)"
        result[donor] = entry
    return result


def job1_gex_cross_check(job1_h5ad_path: Path, gex_cell_donor: dict[str, str], donors: list[str],
                          cohort: str = "Stephenson 2021") -> dict:
    """R3c (review:20261011_124253): the prior version only
    reported POOL PREFIXES of Job 1's selected cell_ids for the disputed
    donors -- it never actually joined them against GEX's own per-cell
    donor map. This does the real full-CellID join: for every cell Job 1
    actually selected for `donors`, checks whether that exact CellID is
    in `gex_cell_donor` (it should always be, since Job 1 is itself built
    from this GEX file) and whether the donor label agrees, and reports
    per-(pool, lineage) cell counts -- not just a set of pool prefixes."""
    cohort_arr, donor_arr, lineage_arr, _outcome_arr, cell_id_arr = _read_job1_obs_arrays(job1_h5ad_path)
    cohort_mask = cohort_arr == cohort

    result = {}
    for donor in donors:
        mask = cohort_mask & (donor_arr == donor)
        cell_ids = cell_id_arr[mask].tolist()
        lineages = lineage_arr[mask].tolist()

        n_matched = 0
        n_unmatched = 0
        mismatches = []
        per_pool_lineage: dict[tuple, int] = defaultdict(int)
        for cid, lineage in zip(cell_ids, lineages):
            pool = cid.split("_", 1)[0] if cid and "_" in cid else "UNKNOWN"
            per_pool_lineage[(pool, lineage)] += 1
            gex_donor = gex_cell_donor.get(cid)
            if gex_donor is None:
                n_unmatched += 1
                continue
            n_matched += 1
            if gex_donor != donor:
                mismatches.append({"cell_id": cid, "job1_donor_id": donor, "gex_donor_id": gex_donor})

        result[donor] = {
            "n_cells_total": len(cell_ids),
            "n_matched_to_gex_cell_id": n_matched,
            "n_unmatched_to_gex_cell_id": n_unmatched,
            "n_donor_label_mismatches": len(mismatches),
            "donor_label_mismatch_examples": mismatches[:10],
            "n_cells_by_pool_and_lineage": {
                f"{pool}::{lineage}": n for (pool, lineage), n in sorted(per_pool_lineage.items())
            },
        }
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gex-h5ad", required=True, type=Path)
    p.add_argument("--job1-h5ad", required=True, type=Path)
    p.add_argument("--tcr-tsv", required=True, type=Path)
    p.add_argument("--bcr-gz", required=True, type=Path)
    p.add_argument("--donors", default=",".join(DEFAULT_DONORS))
    p.add_argument("--exclude-outcome-values", nargs="*", default=DEFAULT_EXCLUDE_OUTCOME_VALUES,
                    help="Must match run_batch_diagnostics.py's own flag/default for the "
                         "outcome-eligibility counts to mean the same thing.")
    p.add_argument("--out-dir", required=True, type=Path)
    args = p.parse_args()
    donors = [d.strip() for d in args.donors.split(",") if d.strip()]

    # R1 fix (review:20261011_124253): same pattern as
    # audit_stephenson_source_identity.py / run_batch_diagnostics.py --
    # the non-empty check catches leftover files from a completed/crashed
    # run; claim_ownership() closes the race for two fresh concurrent
    # starts via an OS-atomic exclusive create.
    if args.out_dir.exists() and any(args.out_dir.iterdir()):
        raise FileExistsError(
            f"{args.out_dir} already exists and is non-empty -- refusing to overwrite a "
            f"prior result; use a fresh, newly timestamped --out-dir."
        )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    ownership = claim_ownership(args.out_dir / ".run_owner.json", command=" ".join(sys.argv))

    write_source_manifest(
        args.out_dir,
        input_paths={
            "gex_h5ad": args.gex_h5ad, "job1_h5ad": args.job1_h5ad,
            "tcr_tsv": args.tcr_tsv, "bcr_gz": args.bcr_gz,
        },
        argv=sys.argv,
        script_path=Path(__file__),
        dependency_module_paths=(
            CODE_DIR / "audit_stephenson_source_identity.py",
            CODE_DIR / "validate_job1_output.py",
        ),
        run_ownership=ownership,
    )

    # gex_cell_donor is reused (not re-derived) for the Job1<->GEX join below
    # -- the same validated, duplicate-checked map parse_gex_cambridge builds
    # for the main audit tool, rather than a second, possibly-divergent read.
    _gex_pools, _gex_pool_cell_counts, gex_cell_donor, _gex_conflicts, _gex_coverage = (
        parse_gex_cambridge(args.gex_h5ad)
    )

    report = {
        "donors": donors,
        "gex": gex_donor_pool_detail(args.gex_h5ad, donors),
        "tcr": tcr_donor_detail(args.tcr_tsv, donors),
        "bcr": bcr_donor_detail(args.bcr_gz, donors),
        "job1_output": job1_output_detail(args.job1_h5ad, donors,
                                           exclude_outcome_values=args.exclude_outcome_values),
        "job1_gex_cross_check": job1_gex_cross_check(args.job1_h5ad, gex_cell_donor, donors),
    }
    (args.out_dir / "disputed_donor_impact.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
