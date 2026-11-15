"""Runnable, fully offline, end-to-end example of the PCA embedding CONTROL
baseline, per `review:20261013_135634` section 4 item 5: a small
synthetic HUMAN dataset and a small synthetic MOUSE dataset through the
SAME interface (`embedding_providers.PCAEmbeddingProvider`), validated via
the real `dataset_schema`/`cli` input checks, a REAL `sklearn.decomposition
.PCA` call (no mocking), persisted artifacts, and a manifest re-read.

    python scripts/34_scfm_discovery_framework/examples/run_pca_example.py \
        --out-dir <new/empty/directory>

Every input here is clearly synthetic (`TOYGENE*`/`TOYCELL*` IDs, random
counts) -- this is NOT a run on real study data or any other real
dataset, and this script makes no claim that PCA "discovers" anything
biological; it only demonstrates the input -> PCA -> readable-artifact
path working identically for two species through one interface. Never
makes a network call. Everything this script writes goes under
`--out-dir`; it refuses a non-empty one, same convention as the rest of
this framework.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.join(os.path.dirname(HERE), "code")
sys.path.insert(0, CODE_DIR)


def build_synthetic_h5ad(path: str, *, species: str, n_cells: int = 24, n_genes: int = 30,
                         seed: int) -> None:
    import anndata as ad
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(seed)
    cell_ids = [f"TOYCELL{species.upper()}{i}" for i in range(n_cells)]
    donor_ids = [f"D{i % 4}" for i in range(n_cells)]
    sample_ids = [f"S{i}" for i in range(n_cells)]
    chip_ids = ["CHIP1" if i % 2 == 0 else "CHIP2" for i in range(n_cells)]
    X = rng.integers(0, 8, size=(n_cells, n_genes)).astype("float32")
    obs = pd.DataFrame({
        "cell_barcode": cell_ids, "donor_id": donor_ids,
        "sample_id": sample_ids, "chip_id": chip_ids,
    })
    var = pd.DataFrame(index=[f"TOYGENE{i}" for i in range(n_genes)])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.write_h5ad(path)


def build_dataset_config(dataset_id: str, species: str, matrix_path: str) -> dict:
    return {
        "dataset_id": dataset_id, "species": species,
        "tissue_or_target_cell_type": "synthetic_toy_tissue",
        "subject_id_column": "donor_id", "sample_id_column": "sample_id",
        "biological_replicate_unit": "donor_id", "technical_batch_column": "chip_id",
        "expression_matrix_path": matrix_path, "counts_semantics": "raw_counts",
        "counts_slot": "X", "gene_id_namespace": "gene_symbol",
        "cell_id_column": "cell_barcode", "qc_status": "provisional",
        "qc_version_ref": "none (synthetic example)",
    }


def run_one_species(species: str, out_dir: str, seed: int) -> None:
    from dataset_schema import validate_dataset_config
    from embedding_providers import PCAEmbeddingProvider, read_embedding_result, load_embedding_matrices
    import cli as C

    print(f"\n### species={species!r} ###")
    print(f"## Step 1/4 -- build a tiny synthetic {species} h5ad (TOYGENE*/TOYCELL* IDs, "
          "random counts -- never real data)")
    h5ad_path = os.path.join(out_dir, f"{species}.h5ad")
    build_synthetic_h5ad(h5ad_path, species=species, seed=seed)
    # Per cli.resolve_matrix_path's own contract: a RELATIVE
    # expression_matrix_path resolves against the dataset config FILE's own
    # directory, not the process cwd -- store it absolute here to avoid
    # re-deriving that directory twice (out_dir is already baked into
    # h5ad_path above).
    config = build_dataset_config(f"toy_{species}_pilot", species, os.path.abspath(h5ad_path))
    config_path = os.path.join(out_dir, f"{species}_dataset_config.json")
    with open(config_path, "w") as fh:
        json.dump(config, fh, indent=2)
    print(f"wrote: {h5ad_path}, {config_path}")

    print("\n## Step 2/4 -- validate via the SAME dataset_schema/cli checks the rest "
          "of this framework uses (never a bare file path)")
    result = validate_dataset_config(config)
    if not result.ok:
        raise SystemExit(f"dataset config invalid: {[str(i) for i in result.issues]}")
    structural_problems = C.validate_matrix_structure(h5ad_path, config)
    if structural_problems:
        raise SystemExit(f"matrix structure invalid: {structural_problems}")
    print("OK: dataset_schema.validate_dataset_config and cli.validate_matrix_structure "
          "both pass.")

    print("\n## Step 3/4 -- run the REAL sklearn-backed PCAEmbeddingProvider "
          "(no mocking, no network) -- an explicit CONTROL, not a scFM")
    pca_out = os.path.join(out_dir, f"{species}_pca")
    provider = PCAEmbeddingProvider()
    pca_result = provider.fit_transform(
        config, config_path, out_dir=pca_out, n_components=3, log1p=True,
        random_seed=0, svd_solver="full")
    if pca_result.status != "complete":
        raise SystemExit(f"PCA run did not complete: status={pca_result.status} "
                         f"reason={pca_result.reason}")
    print(f"embedding: {len(pca_result.cell_ids)} cells x {pca_result.n_components} "
          f"components, {len(pca_result.feature_ids)} genes used, "
          f"explained_variance_ratio={[round(v, 3) for v in pca_result.explained_variance_ratio]}")

    print("\n## Step 4/4 -- read the persisted manifest + npz back and verify shapes")
    manifest = read_embedding_result(os.path.join(pca_out, "EMBEDDING_RESULT.json"))
    arrays = load_embedding_matrices(manifest)
    expected_shape = (len(manifest["cell_ids"]), manifest["n_components"])
    if arrays["embedding"].shape != expected_shape:
        raise SystemExit(
            f"MISMATCH: read-back embedding shape {arrays['embedding'].shape} != "
            f"expected {expected_shape}")
    print(f"OK: read_embedding_result + load_embedding_matrices round-trip confirms "
          f"embedding.shape == {expected_shape}, hash-checked against "
          f"embedding_sha256={manifest['embedding_sha256'][:16]}...")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    if os.path.exists(args.out_dir) and os.listdir(args.out_dir):
        raise SystemExit(f"{args.out_dir} must be new or empty")
    os.makedirs(args.out_dir, exist_ok=True)

    run_one_species("human", args.out_dir, seed=0)
    run_one_species("mouse", args.out_dir, seed=1)

    print(f"\nAll artifacts are under: {os.path.abspath(args.out_dir)}")
    print("NOTE: both runs are fully synthetic (TOYGENE*/TOYCELL* IDs, random counts). "
          "This does NOT claim to have run real study data or any real dataset, "
          "and PCA is a baseline CONTROL here, never evidence of scFM superiority or of "
          "any biological finding.")


if __name__ == "__main__":
    main()
