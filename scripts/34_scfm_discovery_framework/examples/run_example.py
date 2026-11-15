"""Runnable, end-to-end example of the scFM discovery framework's config/
validate/plan contracts, per `review:20261012_163250`'s "交付缺口"
item: a real config template, a real small fixture, and the full
`validate --check-matrix-structure -> plan --check-matrix-structure ->
read_manifest` chain with an actual object-equality check, not just "a test
fixture exists somewhere."

    python scripts/34_scfm_discovery_framework/examples/run_example.py \
        --out-dir <new/empty/directory>

Builds a tiny SYNTHETIC human monocyte fixture (6 cells): donor `D1`
contributes TWO samples (`S1`, `S2`) -- a genuine biological-replicate example,
never a cell or a sample standing in for an independent subject -- and donor
`D2` contributes one sample (`S3`). Never touches real data, never submits a
job, never makes a network call. Everything this script writes goes under
`--out-dir`; it refuses a non-empty one, same convention as `cli.py plan`.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CODE_DIR = os.path.join(os.path.dirname(HERE), "code")
CLI = os.path.join(CODE_DIR, "cli.py")
sys.path.insert(0, CODE_DIR)


def build_fixture(data_dir: str) -> str:
    """Writes the tiny synthetic h5ad and its matching dataset config (with a
    RELATIVE expression_matrix_path, to exercise config-dir-relative
    resolution, not just absolute paths) into `data_dir`. Returns the config
    path.
    """
    import anndata as ad
    import numpy as np
    import pandas as pd

    cell_barcodes = [f"CELL{i}" for i in range(6)]
    donor_ids = ["D1", "D1", "D1", "D1", "D2", "D2"]
    sample_ids = ["S1", "S1", "S2", "S2", "S3", "S3"]
    chip_ids = ["CHIP1", "CHIP1", "CHIP2", "CHIP2", "CHIP1", "CHIP2"]
    n_vars = 8

    obs = pd.DataFrame({
        "cell_barcode": cell_barcodes, "donor_id": donor_ids,
        "sample_id": sample_ids, "chip_id": chip_ids,
    })
    var = pd.DataFrame(index=[f"ENSG{i:08d}" for i in range(n_vars)])
    X = np.random.default_rng(0).integers(0, 5, size=(6, n_vars)).astype("float32")
    a = ad.AnnData(X=X, obs=obs, var=var)
    h5ad_name = "tiny_monocyte_example.h5ad"
    a.write_h5ad(os.path.join(data_dir, h5ad_name))

    config = {
        "dataset_id": "example_monocyte_pilot",
        "species": "human",
        "tissue_or_target_cell_type": "monocyte",
        "subject_id_column": "donor_id",
        "sample_id_column": "sample_id",
        "biological_replicate_unit": "donor_id",
        "technical_batch_column": "chip_id",
        # Relative to THIS config file's own directory, not the caller's cwd
        # (review:20261012_163250 finding 2) -- resolved by
        # cli.py's resolve_matrix_path, never by the shell's cwd.
        "expression_matrix_path": h5ad_name,
        "counts_semantics": "raw_counts",
        "counts_slot": "X",
        "gene_id_namespace": "ensembl_gene_id",
        "cell_id_column": "cell_barcode",
        "qc_status": "provisional",
        "qc_version_ref": "example fixture, not a real QC pass",
    }
    config_path = os.path.join(data_dir, "dataset_config.json")
    with open(config_path, "w") as fh:
        json.dump(config, fh, indent=2, sort_keys=True)
    return config_path


def run(cmd: list[str], cwd: str | None = None) -> subprocess.CompletedProcess:
    print(f"$ {' '.join(cmd)}" + (f"   (cwd={cwd})" if cwd else ""))
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    if os.path.exists(args.out_dir) and os.listdir(args.out_dir):
        raise SystemExit(f"{args.out_dir} must be new or empty")
    data_dir = os.path.join(args.out_dir, "data")
    plan_dir = os.path.join(args.out_dir, "plan")
    os.makedirs(data_dir, exist_ok=True)

    print("## Step 1/4 -- build the tiny synthetic fixture + config template")
    config_path = build_fixture(data_dir)
    print(f"wrote: {config_path}")
    print(f"wrote: {os.path.join(data_dir, 'tiny_monocyte_example.h5ad')}\n")

    print("## Step 2/4 -- validate --check-matrix-structure "
          "(run from the config's OWN directory, to exercise relative-path "
          "resolution, not the repo root)")
    r1 = run([sys.executable, CLI, "validate", "--dataset-config",
             "dataset_config.json", "--check-matrix-structure"], cwd=data_dir)
    if r1.returncode != 0:
        raise SystemExit("validate failed -- example is not self-consistent")

    print("## Step 3/4 -- plan --check-matrix-structure")
    r2 = run([sys.executable, CLI, "plan", "--dataset-config",
             "dataset_config.json", "--out", os.path.abspath(plan_dir),
             "--check-matrix-structure"], cwd=data_dir)
    if r2.returncode != 0:
        raise SystemExit("plan failed -- example is not self-consistent")

    print("## Step 4/4 -- read_manifest full-object round-trip check")
    import artifact_schema as A
    plan_path = os.path.join(plan_dir, "PLAN.json")
    with open(plan_path) as fh:
        written_raw = json.load(fh)
    read_back = A.read_manifest(plan_path)
    if read_back != written_raw:
        raise SystemExit(
            "MISMATCH: read_manifest's object does not equal the raw written JSON "
            f"-- written={written_raw}\nread_back={read_back}")
    print(f"OK: read_manifest({plan_path!r}) == the exact JSON written to disk")
    print(f"stage_status={read_back['stage_status']!r}  "
          f"model_identity_or_pending={read_back['model_identity_or_pending']!r}")
    print(f"\nAll artifacts are under: {os.path.abspath(args.out_dir)}")


if __name__ == "__main__":
    main()
