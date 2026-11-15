"""Runnable, fully offline, end-to-end example of the raw-count-sum
pseudobulk aggregation increment, per `review:20261014_100436`'s
closing instruction: a small synthetic HUMAN dataset and a small synthetic
MOUSE dataset through the SAME interface
(`pseudobulk_providers.PseudobulkAggregationProvider`), a real aggregation,
persisted artifacts, a manifest+npz read-back with the F1 identity cross-
checks exercised, AND a genuinely invalid input producing a readable
'failed' manifest (never a crash, never an orphan npz) -- all through the
REAL production entry points, never a test-only helper.

    python scripts/34_scfm_discovery_framework/examples/run_pseudobulk_example.py \
        --out-dir <new/empty/directory>

Every input here is clearly synthetic (`TOYGENE*`/`TOYCELL*` IDs, small
hand-checkable integer counts) -- this is NOT a run on real data, and this
script computes NO DE statistic, p-value, or `GeneRanking` -- aggregation
only, per this increment's explicit scope.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.join(os.path.dirname(HERE), "code")
sys.path.insert(0, CODE_DIR)


def build_synthetic_h5ad(path: str, *, species: str, rawslot: bool) -> dict:
    """8 cells, 2 donors (D0/D1), 4 samples (S0..S3) -- each donor
    contributes TWO samples (a paired visit structure), 3 genes, small
    hand-checkable integer counts. `rawslot=True` (used for the mouse
    case) puts the counts in `.raw.X` instead of `.X`, exercising the same
    gene-axis-via-`raw` path `embedding_providers.py`'s F2 fix covers.
    """
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse

    cell_ids = [f"TOYCELL{species.upper()}{i}" for i in range(8)]
    sample_ids = ["S0", "S0", "S1", "S1", "S2", "S2", "S3", "S3"]
    donor_ids = ["D0", "D0", "D0", "D0", "D1", "D1", "D1", "D1"]
    conditions = ["healthy", "healthy", "disease", "disease",
                 "healthy", "healthy", "disease", "disease"]
    visits = ["v1", "v1", "v2", "v2", "v1", "v1", "v2", "v2"]
    X = np.array([
        [1, 4, 0], [3, 2, 1],   # S0 sum = [4, 6, 1], lib=11
        [2, 0, 1], [0, 1, 2],   # S1 sum = [2, 1, 3], lib=6
        [1, 1, 1], [1, 1, 1],   # S2 sum = [2, 2, 2], lib=6
        [0, 0, 5], [5, 0, 0],   # S3 sum = [5, 0, 5], lib=10
    ], dtype="float32")
    obs = pd.DataFrame({
        "cell_barcode": cell_ids, "donor_id": donor_ids, "sample_id": sample_ids,
        "condition": conditions, "visit": visits,
    })
    var = pd.DataFrame(index=[f"TOYGENE{i}" for i in range(3)])
    a = ad.AnnData(X=sparse.csr_matrix(X) if rawslot else X.copy(), obs=obs, var=var)
    if rawslot:
        a.raw = a.copy()
        a = a[:, :1].copy()  # .var now disagrees with .raw.var -- exercises F2's gene-axis-via-raw path
    a.write_h5ad(path)
    expected = {"S0": [4.0, 6.0, 1.0], "S1": [2.0, 1.0, 3.0],
               "S2": [2.0, 2.0, 2.0], "S3": [5.0, 0.0, 5.0]}
    return {"cell_ids": cell_ids, "expected_sample_sums": expected}


def build_dataset_config(dataset_id: str, species: str, matrix_path: str, rawslot: bool) -> dict:
    from dataset_schema import NO_BATCH_SENTINEL
    return {
        "dataset_id": dataset_id, "species": species,
        "tissue_or_target_cell_type": "synthetic_toy_tissue",
        "subject_id_column": "donor_id", "sample_id_column": "sample_id",
        "biological_replicate_unit": "donor_id", "technical_batch_column": NO_BATCH_SENTINEL,
        "expression_matrix_path": matrix_path, "counts_semantics": "raw_counts",
        "counts_slot": "raw.X" if rawslot else "X", "gene_id_namespace": "gene_symbol",
        "cell_id_column": "cell_barcode", "qc_status": "provisional",
        "qc_version_ref": "none (synthetic example)",
    }


def run_one_species(species: str, out_dir: str, rawslot: bool) -> None:
    from dataset_schema import validate_dataset_config
    from pseudobulk_providers import (PseudobulkAggregationProvider, read_pseudobulk_result,
                                      load_pseudobulk_counts)
    import cli as C

    print(f"\n### species={species!r} (counts_slot={'raw.X' if rawslot else 'X'}) ###")
    print(f"## Step 1/5 -- build a tiny synthetic {species} h5ad (TOYGENE*/TOYCELL* IDs, "
          "hand-checkable small integer counts -- never real data)")
    h5ad_path = os.path.join(out_dir, f"{species}.h5ad")
    fixture = build_synthetic_h5ad(h5ad_path, species=species, rawslot=rawslot)
    config = build_dataset_config(f"toy_{species}_pilot", species,
                                  os.path.abspath(h5ad_path), rawslot)
    config_path = os.path.join(out_dir, f"{species}_dataset_config.json")
    with open(config_path, "w") as fh:
        json.dump(config, fh, indent=2)
    print(f"wrote: {h5ad_path}, {config_path}")

    print("\n## Step 2/5 -- validate via the SAME dataset_schema/cli checks the rest "
          "of this framework uses")
    result = validate_dataset_config(config)
    if not result.ok:
        raise SystemExit(f"dataset config invalid: {[str(i) for i in result.issues]}")
    structural_problems = C.validate_matrix_structure(h5ad_path, config)
    if structural_problems:
        raise SystemExit(f"matrix structure invalid: {structural_problems}")
    print("OK: dataset_schema.validate_dataset_config and cli.validate_matrix_structure both pass.")

    print("\n## Step 3/5 -- run the REAL PseudobulkAggregationProvider (sample x cell-"
          "type raw count SUM, never within-donor) through the production entry point")
    pb_out = os.path.join(out_dir, f"{species}_pseudobulk")
    provider = PseudobulkAggregationProvider()
    pb_result = provider.aggregate(
        config, config_path, out_dir=pb_out, cell_ids_subset=tuple(fixture["cell_ids"]),
        cell_type_or_state_label="toy_cell_type", condition_column="condition",
        visit_or_state_column="visit")
    if pb_result.status != "complete":
        raise SystemExit(f"aggregation did not complete: status={pb_result.status} "
                         f"reason={pb_result.reason}")
    print(f"aggregation: {pb_result.n_pseudobulk_samples} pseudobulk sample(s), "
          f"{pb_result.n_distinct_donors} distinct donor(s) (reported SEPARATELY -- "
          "D0/D1 each contribute TWO samples here), "
          f"total_counts_aggregated={pb_result.total_counts_aggregated}")

    print("\n## Step 4/5 -- read the persisted manifest + npz back; the F1 content-"
          "binding checks (per-row library size, per-gene totals, identity_sha256) "
          "all re-verify against the actual loaded array, not merely the JSON's own "
          "self-consistency")
    manifest = read_pseudobulk_result(os.path.join(pb_out, "PSEUDOBULK_RESULT.json"))
    counts = load_pseudobulk_counts(manifest)
    row_by_sid = {r["sample_id"]: i for i, r in enumerate(manifest["sample_rows"])}
    for sid, expected_row in fixture["expected_sample_sums"].items():
        actual_row = counts[row_by_sid[sid]].tolist()
        if actual_row != expected_row:
            raise SystemExit(f"MISMATCH: sample {sid} expected {expected_row}, got {actual_row}")
    print(f"OK: read_pseudobulk_result + load_pseudobulk_counts round-trip confirms every "
          f"hand-computed sample sum, hash-checked against "
          f"aggregation_sha256={manifest['aggregation_sha256'][:16]}..., "
          f"identity_sha256={manifest['identity_sha256'][:16]}...")

    print("\n## Step 5/5 -- a genuinely invalid input (cell_type_or_state_label=True, "
          "not a string) through the SAME production entry point must produce a "
          "READABLE 'failed' manifest -- never a crash, never an orphan npz")
    bad_out = os.path.join(out_dir, f"{species}_bad_label")
    bad_result = provider.aggregate(
        config, config_path, out_dir=bad_out, cell_ids_subset=tuple(fixture["cell_ids"]),
        cell_type_or_state_label=True)
    if bad_result.status != "failed":
        raise SystemExit(f"expected status='failed' for a non-string label, got "
                         f"{bad_result.status!r}")
    if os.path.exists(os.path.join(bad_out, "pseudobulk_counts.npz")):
        raise SystemExit("an orphan npz was written for a failed aggregation -- F2 regressed")
    failed_manifest = read_pseudobulk_result(os.path.join(bad_out, "PSEUDOBULK_RESULT.json"))
    print(f"OK: status='failed', reason={failed_manifest['reason'][:70]!r}..., "
          "no npz file was written, and the failed manifest itself reads back cleanly.")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    if os.path.exists(args.out_dir) and os.listdir(args.out_dir):
        raise SystemExit(f"{args.out_dir} must be new or empty")
    os.makedirs(args.out_dir, exist_ok=True)

    run_one_species("human", args.out_dir, rawslot=False)
    run_one_species("mouse", args.out_dir, rawslot=True)

    print(f"\nAll artifacts are under: {os.path.abspath(args.out_dir)}")
    print("NOTE: both runs are fully synthetic (TOYGENE*/TOYCELL* IDs, hand-checkable "
          "small integer counts). This increment computes NO DE statistic, p-value, or "
          "GeneRanking -- raw-count-sum aggregation only, per this round's explicit scope.")


if __name__ == "__main__":
    main()
