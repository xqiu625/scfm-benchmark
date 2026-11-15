"""Runnable, fully offline, end-to-end example of the pseudobulk -> DE
design-matrix PREPARATION increment (`de_design_providers.py`), per
`review:20261014_203915`'s closing instruction B: a small
synthetic HUMAN independent-donors-unpaired case and a small synthetic
MOUSE fully-paired case, through the SAME interface
(`DEDesignPreparationProvider`), BOTH via the real
`aggregate -> read_pseudobulk_result/load_pseudobulk_counts -> prepare ->
write_prepared_de_design -> read_prepared_de_design/load_prepared_de_counts`
chain -- never a hand-built pseudobulk manifest, and never the OLD
`run_pseudobulk_example.py` (which stops at aggregation and never builds a
design).

    python scripts/34_scfm_discovery_framework/examples/run_de_design_example.py \
        --out-dir <new/empty/directory>

The human case deliberately picks `sample_id`s that sort ALPHABETICALLY
in the OPPOSITE order from their donors (`PseudobulkAggregationProvider`
always orders its own `sample_rows` by `sample_id`, per
`pseudobulk_providers.py`'s `for sid in sorted(by_sample)`), so the
aggregation's OWN `sample_rows` order is the exact reverse of the design's
donor-sorted row order -- confirming `prepare` aligns by `sample_id`
LOOKUP, not by position, through the real producer (never by tampering an
already-frozen manifest after the fact).

Every input here is clearly synthetic (`TOYGENE*`/`TOYCELL*`/`S*`/`D*`
IDs, small hand-checkable integer counts) -- this is NOT a run on real
data. `status='prepared'` here means "a DE input/design was prepared",
nothing more: `ExpressionDERankingProvider` still has no real statistical
backend, and this script computes NO LFC, p-value, or `GeneRanking`.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.join(os.path.dirname(HERE), "code")
sys.path.insert(0, CODE_DIR)


def _dataset_config(dataset_id: str, species: str, matrix_path: str, rawslot: bool) -> dict:
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


def _build_h5ad(path: str, *, species: str, X, sample_ids, donor_ids, conditions, rawslot: bool):
    import anndata as ad
    import numpy as np
    import pandas as pd

    cell_ids = [f"TOYCELL{species.upper()}{i}" for i in range(len(sample_ids))]
    obs = pd.DataFrame({
        "cell_barcode": cell_ids, "donor_id": donor_ids, "sample_id": sample_ids,
        "condition": conditions,
    })
    var = pd.DataFrame(index=[f"TOYGENE{i}" for i in range(X.shape[1])])
    a = ad.AnnData(X=np.asarray(X, dtype="float32"), obs=obs, var=var)
    if rawslot:
        a.raw = a.copy()
    a.write_h5ad(path)
    return cell_ids


def _aggregate_and_load(out_dir: str, *, dataset_id: str, species: str, X, sample_ids, donor_ids,
                        conditions, rawslot: bool):
    from dataset_schema import validate_dataset_config
    from pseudobulk_providers import (PseudobulkAggregationProvider, read_pseudobulk_result,
                                      load_pseudobulk_counts)
    import cli as C

    h5ad_path = os.path.join(out_dir, f"{species}.h5ad")
    cell_ids = _build_h5ad(h5ad_path, species=species, X=X, sample_ids=sample_ids,
                           donor_ids=donor_ids, conditions=conditions, rawslot=rawslot)
    config = _dataset_config(dataset_id, species, os.path.abspath(h5ad_path), rawslot)
    config_path = os.path.join(out_dir, f"{species}_dataset_config.json")
    with open(config_path, "w") as fh:
        json.dump(config, fh, indent=2)

    cfg_result = validate_dataset_config(config)
    if not cfg_result.ok:
        raise SystemExit(f"dataset config invalid: {[str(i) for i in cfg_result.issues]}")
    structural_problems = C.validate_matrix_structure(h5ad_path, config)
    if structural_problems:
        raise SystemExit(f"matrix structure invalid: {structural_problems}")

    pb_out = os.path.join(out_dir, f"{species}_pseudobulk")
    provider = PseudobulkAggregationProvider()
    pb_result = provider.aggregate(
        config, config_path, out_dir=pb_out, cell_ids_subset=tuple(cell_ids),
        cell_type_or_state_label="toy_cell_type", condition_column="condition")
    if pb_result.status != "complete":
        raise SystemExit(f"aggregation did not complete: status={pb_result.status} "
                         f"reason={pb_result.reason}")

    pb_path = os.path.join(pb_out, "PSEUDOBULK_RESULT.json")
    manifest = read_pseudobulk_result(pb_path)
    load_pseudobulk_counts(manifest)  # real read-back check, same as run_pseudobulk_example.py
    return pb_path, manifest


def run_human_unpaired(out_dir: str) -> None:
    import numpy as np
    from gene_ranking_providers import ComparisonSpec, SampleStructure
    from de_design_providers import (DEDesignPreparationProvider, read_prepared_de_design,
                                     load_prepared_de_counts)

    print("\n### species='human', design_type='independent_donors_unpaired' ###")
    print("## Step 1/4 -- aggregate a tiny synthetic human dataset; sample_ids (Z0/Y1/X2/W3) "
          "sort ALPHABETICALLY in the opposite order from their donors (D0/D1/D2/D3), so the "
          "pseudobulk's OWN sample_rows order is the exact reverse of the design's "
          "donor-sorted row order")
    X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]], dtype="float32")
    pb_path, manifest = _aggregate_and_load(
        out_dir, dataset_id="toy_human_monocyte_pilot", species="human", X=X,
        sample_ids=["Z0", "Y1", "X2", "W3"], donor_ids=["D0", "D1", "D2", "D3"],
        conditions=["healthy", "healthy", "disease", "disease"], rawslot=False)
    actual_order = [r["sample_id"] for r in manifest["sample_rows"]]
    if actual_order != ["W3", "X2", "Y1", "Z0"]:
        raise SystemExit(f"fixture bug: expected pseudobulk's own order ['W3','X2','Y1','Z0'], "
                         f"got {actual_order}")
    print(f"aggregation OK: sample_rows order is {actual_order} (the REVERSE of the "
          "design's donor-sorted order)")

    print("\n## Step 2/4 -- prepare the DE design through the REAL production entry point "
          "(DEDesignPreparationProvider.prepare, which itself calls read_pseudobulk_result + "
          "load_pseudobulk_counts -- per codex_output/20261014_203915.md F1, there is no "
          "caller-supplied counts parameter any more)")
    comparison = ComparisonSpec(
        group_column="condition", reference_level="healthy", comparison_level="disease",
        direction_definition="disease vs. healthy, positive = higher in disease")
    sample_structure = SampleStructure(
        donor_or_animal_id_column="donor_id", sample_id_column="sample_id",
        repeated_sampling_or_pairing="none_single_sample_per_donor", group_column="condition",
        design_type="independent_donors_unpaired")
    prep_out = os.path.join(out_dir, "human_de_design")
    provider = DEDesignPreparationProvider()
    result = provider.prepare(pseudobulk_result_path=pb_path, comparison=comparison,
                              sample_structure=sample_structure, out_dir=prep_out)
    if result.status != "prepared":
        raise SystemExit(f"expected status='prepared', got {result.status!r}: {result.reason}")

    expected_sample_ids = ("Z0", "Y1", "X2", "W3")  # re-sorted to donor order, NOT pseudobulk's own
    expected_matrix = ((1.0, 0.0), (1.0, 0.0), (1.0, 1.0), (1.0, 1.0))
    if (result.sample_ids, result.design_matrix, result.rank, result.residual_degrees_of_freedom,
       result.contrast_vector) != (expected_sample_ids, expected_matrix, 2, 2, (0.0, 1.0)):
        raise SystemExit(
            f"MISMATCH against hand-computed expectation: sample_ids={result.sample_ids}, "
            f"design_matrix={result.design_matrix}, rank={result.rank}, "
            f"residual_df={result.residual_degrees_of_freedom}, "
            f"contrast_vector={result.contrast_vector}")
    print(f"OK: realigned to donor-sorted order {result.sample_ids}, design_matrix="
          f"{result.design_matrix}, rank={result.rank}, "
          f"residual_degrees_of_freedom={result.residual_degrees_of_freedom}, "
          f"contrast={result.contrast_name}={result.contrast_vector} "
          f"({result.contrast_direction})")

    print("\n## Step 3/4 -- write/read/load round trip through the REAL production entry "
          "points; the prepared counts submatrix must match the REORDERED hand-computed rows")
    reloaded = read_prepared_de_design(os.path.join(prep_out, "PREPARED_DE_DESIGN.json"))
    prepared_counts = load_prepared_de_counts(reloaded)
    expected_counts = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0], [7.0, 8.0]]  # Z0,Y1,X2,W3, donor-sorted
    if prepared_counts.tolist() != expected_counts:
        raise SystemExit(f"MISMATCH: expected prepared counts {expected_counts}, "
                         f"got {prepared_counts.tolist()}")
    print(f"OK: round trip confirms prepared_de_counts == {expected_counts}, species="
          f"{reloaded['dataset_identity']['species']!r}, "
          f"gene_id_namespace={reloaded['dataset_identity']['gene_id_namespace']!r}")

    print("\n## Step 4/4 -- a declared column name that does NOT match what the pseudobulk "
          "aggregation actually used (per F2) must be refused, never silently accepted just "
          "because condition_value happens to already be populated")
    bad_comparison = ComparisonSpec(
        group_column="not_the_real_column", reference_level="healthy", comparison_level="disease",
        direction_definition="disease vs. healthy, positive = higher in disease")
    bad_out = os.path.join(out_dir, "human_de_design_bad_column")
    bad_result = provider.prepare(pseudobulk_result_path=pb_path, comparison=bad_comparison,
                                  sample_structure=sample_structure, out_dir=bad_out)
    if bad_result.status != "failed":
        raise SystemExit(f"expected status='failed' for a mismatched group_column, got "
                         f"{bad_result.status!r}")
    if os.path.exists(os.path.join(bad_out, "prepared_de_counts.npz")):
        raise SystemExit("an orphan counts npz was written for a failed preparation")
    print(f"OK: status='failed', reason={bad_result.reason[:90]!r}..., no orphan npz written.")


def run_mouse_paired(out_dir: str) -> None:
    import numpy as np
    from gene_ranking_providers import ComparisonSpec, SampleStructure
    from de_design_providers import (DEDesignPreparationProvider, read_prepared_de_design,
                                     load_prepared_de_counts)

    print("\n### species='mouse', design_type='paired_or_repeated_within_donor' ###")
    print("## Step 1/3 -- aggregate a tiny synthetic mouse dataset (raw.X slot, exercising "
          "the gene-axis-via-raw path too), 2 donors each FULLY paired (one healthy + one "
          "disease sample)")
    X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]], dtype="float32")
    pb_path, manifest = _aggregate_and_load(
        out_dir, dataset_id="toy_mouse_monocyte_pilot", species="mouse", X=X,
        sample_ids=["S0", "S1", "S2", "S3"], donor_ids=["D0", "D0", "D1", "D1"],
        conditions=["healthy", "disease", "healthy", "disease"], rawslot=True)

    print("\n## Step 2/3 -- prepare the FULLY-PAIRED DE design through the SAME interface "
          "as the human unpaired case above")
    comparison = ComparisonSpec(
        group_column="condition", reference_level="healthy", comparison_level="disease",
        direction_definition="disease vs. healthy, positive = higher in disease")
    sample_structure = SampleStructure(
        donor_or_animal_id_column="donor_id", sample_id_column="sample_id",
        repeated_sampling_or_pairing="paired:pre_post", group_column="condition",
        design_type="paired_or_repeated_within_donor")
    prep_out = os.path.join(out_dir, "mouse_de_design")
    provider = DEDesignPreparationProvider()
    result = provider.prepare(pseudobulk_result_path=pb_path, comparison=comparison,
                              sample_structure=sample_structure, out_dir=prep_out)
    if result.status != "prepared":
        raise SystemExit(f"expected status='prepared', got {result.status!r}: {result.reason}")

    expected_matrix = ((1.0, 0.0, 0.0), (1.0, 0.0, 1.0), (1.0, 1.0, 0.0), (1.0, 1.0, 1.0))
    if (result.design_matrix, result.rank, result.residual_degrees_of_freedom,
       result.contrast_vector) != (expected_matrix, 3, 1, (0.0, 0.0, 1.0)):
        raise SystemExit(
            f"MISMATCH against hand-computed expectation: design_matrix={result.design_matrix}, "
            f"rank={result.rank}, residual_df={result.residual_degrees_of_freedom}, "
            f"contrast_vector={result.contrast_vector}")
    print(f"OK (full rank, residual_degrees_of_freedom={result.residual_degrees_of_freedom} -- "
          "exactly at the >=1 boundary): design_matrix="
          f"{result.design_matrix}, contrast={result.contrast_name}={result.contrast_vector}")

    print("\n## Step 3/3 -- write/read/load round trip")
    reloaded = read_prepared_de_design(os.path.join(prep_out, "PREPARED_DE_DESIGN.json"))
    prepared_counts = load_prepared_de_counts(reloaded)
    expected_counts = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0], [7.0, 8.0]]
    if prepared_counts.tolist() != expected_counts:
        raise SystemExit(f"MISMATCH: expected prepared counts {expected_counts}, "
                         f"got {prepared_counts.tolist()}")
    print(f"OK: round trip confirms prepared_de_counts == {expected_counts}, species="
          f"{reloaded['dataset_identity']['species']!r} (SAME interface as the human case, "
          "a different species through the gene-axis-via-raw path).")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    if os.path.exists(args.out_dir) and os.listdir(args.out_dir):
        raise SystemExit(f"{args.out_dir} must be new or empty")
    os.makedirs(args.out_dir, exist_ok=True)

    run_human_unpaired(args.out_dir)
    run_mouse_paired(args.out_dir)

    print(f"\nAll artifacts are under: {os.path.abspath(args.out_dir)}")
    print("NOTE: both runs are fully synthetic. `status='prepared'` means a DE input/design "
          "was prepared, NOTHING more -- ExpressionDERankingProvider still has no real "
          "statistical backend; this script computes NO LFC, p-value, or GeneRanking.")


if __name__ == "__main__":
    main()
