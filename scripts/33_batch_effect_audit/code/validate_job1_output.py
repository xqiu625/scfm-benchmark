#!/usr/bin/env python3
"""Bounded, read-only validation of `prepare_provisional_diagnostic_input.py`'s
output h5ad, reading the on-disk HDF5 structure directly via `h5py` --
NOT `anndata.read_h5ad(..., backed="r")`.

Why not backed mode: a 2026-10-10 review (review:20261010_191146,
finding R1) read the installed anndata source on HPCC directly and found
`backed="r"` already materializes every sparse `layers` entry in full.

v3 (this file) responds to review:20261011_113512 finding F4,
which reviewed v2 and found the structural checks still had gaps:

  - `_n_obs_from_metadata_only` trusted `X.attrs['shape']` without cross-
    checking it against the obs index dataset's own independent length.
    Now computes both (when both are available) and raises if they
    disagree, instead of silently preferring one.
  - `_check_csr_structure` checked `indptr` length/endpoints but never the
    group's own declared `shape` attribute against `(n_obs, n_vars)`, nor
    `indptr`/`indices` dtype (must be an integer kind) or rank (must be
    1-D). A corrupted `shape` attribute, a float-typed `indptr`, or a 2-D
    `indices` array could all pass the v2 checks. Now checked explicitly.
  - `--out`'s check-then-write had a TOCTOU window between `Path.exists()`
    and `write_text()`. Now opens with `os.O_CREAT | os.O_EXCL` (atomic
    exclusive create), which fails if anything creates the file in between.

v2 (this file) responds to review:20261010_202340 finding F4,
which reviewed v1 and found the claimed bounds/rejection guarantees were
not actually enforced:

  - the `n_obs` ceiling only checked `X.attrs['shape']`; if that attribute
    were absent it silently proceeded. Now also falls back to the obs
    index dataset's own `.shape` (metadata-only, no data read) and refuses
    if NEITHER is available.
  - `indptr` was read and checked for monotonicity only; `len(indptr) ==
    n_obs+1`, `indptr[0] == 0`, and `indptr[-1] == len(data) == len(indices)`
    were never checked -- an HDF5 slice past an array's real end can
    silently truncate rather than error, so a slice's returned length is
    not itself proof the array was as long as assumed. Now checked via
    `.shape` (cheap) before any `data`/`indices` slicing.
  - `--cells-per-group 0` produced an empty sample whose (vacuously true)
    numeric checks still reported `True`. Now rejected outright.
  - no per-row or cumulative nnz budget was enforced before slicing --
    corrupted `indptr` metadata could not be caught structurally. Now
    `MAX_NNZ_PER_SAMPLED_ROW` / `MAX_TOTAL_SAMPLED_NNZ` are enforced BEFORE
    any `data`/`indices` read for the rows in question.
  - `stage_id`/`expression_source` were read and recorded but never
    compared to Job 1's actual documented contract. Now compared, with a
    hard-coded expected-value default overridable via `--expected-*`.
  - the manifest was only checked for row count and `cell_id` order;
    `cohort`/`donor_id`/`harmonized_lineage` were never compared row-by-row
    against `obs`. Now they are.
  - "layer-copy" checks (`X` vs `counts_decontx`, `counts` vs
    `counts_decontx`) used `np.allclose` (a tolerance). The preparer writes
    these as literal duplicates of the same in-memory array, so they are
    now compared with exact `np.array_equal` -- NOT a tolerance.
  - the sampled-read byte accounting was named `..._nnz_entries` but
    reported nnz counts, not bytes. Renamed `nnz_entries_read_per_array`;
    an actual byte estimate (`dtype.itemsize`-derived) is reported
    separately.
  - `--out` could silently overwrite a prior report. Now refused.

This tool never opens `layers`/`X` wholesale. `obs`/`var`/`uns` are read
directly (cheap: bounded by `n_obs`) and reported as **metadata checks**,
run over ALL cells. Sparse layers are read ONLY as CSR row-slices for a
small, fixed, deterministic, outcome-blind sample of cells, reported
separately as **sampled numeric checks** -- never merged into one verdict.

Explicitly NOT verified: full numeric equality between every cell's
`counts_raw`/`counts_decontx` and the source Matrix Market files.
"""
from __future__ import annotations

import argparse
import errno
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "00_shared/code"))
from stable_seed import derive_seed  # noqa: E402

MAX_OBS_ROWS_FOR_METADATA_READ = 200_000
MAX_NNZ_PER_SAMPLED_ROW = 50_000
MAX_TOTAL_SAMPLED_NNZ = 2_000_000

REQUIRED_OBS_COLUMNS = [
    "cohort", "donor_id", "harmonized_lineage", "lineage", "patient_id", "outcome",
]
REQUIRED_UNS_KEYS = ["population_stage_id", "expression_source"]
REQUIRED_LAYERS = ["counts_raw", "counts_decontx", "counts"]
EXPECTED_STAGE_ID_DEFAULT = "pre_gate2_decontx_provisional"
EXPECTED_EXPRESSION_SOURCE_DEFAULT = "decontx_counts"


class ValidationError(ValueError):
    pass


def _decode_bytes(value):
    return value.decode() if isinstance(value, bytes) else value


def _file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(1 << 20)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def _script_git_revision() -> str:
    try:
        out = subprocess.run(
            ["git", "-C", str(Path(__file__).resolve().parent), "rev-parse", "HEAD"],
            capture_output=True, text=True, timeout=10,
        )
        return out.stdout.strip() if out.returncode == 0 else "UNVERIFIED (git rev-parse failed)"
    except Exception as exc:  # pragma: no cover
        return f"UNVERIFIED ({exc})"


def read_categorical_or_plain(node: "h5py.Group | h5py.Dataset") -> np.ndarray:
    if isinstance(node, h5py.Group):
        if "categories" not in node or "codes" not in node:
            raise ValidationError(
                f"{node.name}: a Group without both 'categories' and 'codes' is not "
                f"the modern AnnData categorical encoding this tool supports -- "
                f"refusing to guess at a legacy encoding."
            )
        enc = node.attrs.get("encoding-type")
        if enc not in (None, "categorical"):
            raise ValidationError(f"{node.name}: unexpected encoding-type {enc!r}")
        cats = np.array([_decode_bytes(c) for c in node["categories"][:]])
        codes = node["codes"][:]
        out = np.empty(len(codes), dtype=object)
        valid = codes >= 0
        out[valid] = cats[codes[valid]]
        out[~valid] = None
        return out
    arr = node[:]
    if arr.dtype.kind in ("S", "O"):
        return np.array([_decode_bytes(v) for v in arr])
    return arr


def _obs_index_key(f: h5py.File) -> str:
    return _decode_bytes(f["obs"].attrs.get("_index", "cell_id"))


def _n_obs_from_obs_index_shape_only(f: h5py.File) -> int:
    """Metadata-only (shape, never contents) length of the obs index."""
    index_key = _obs_index_key(f)
    index_node = f["obs"][index_key]
    if isinstance(index_node, h5py.Group):
        return int(index_node["codes"].shape[0])
    return int(index_node.shape[0])


def _n_obs_from_metadata_only(f: h5py.File) -> int:
    """Both branches read only dataset/attribute SHAPE, never contents.
    When X's own declared shape IS available, it is cross-checked against
    the obs index's independent length rather than trusted alone -- a
    corrupted `X.attrs['shape']` that still passed v2's checks would
    otherwise go undetected (F4, 20261011_113512.md)."""
    x_shape = f["X"].attrs.get("shape")
    obs_index_n = _n_obs_from_obs_index_shape_only(f)
    if x_shape is not None:
        x_n_obs = int(x_shape[0])
        if x_n_obs != obs_index_n:
            raise ValidationError(
                f"X.attrs['shape'][0]={x_n_obs} disagrees with the obs index "
                f"dataset's own length={obs_index_n} -- refusing to pick one."
            )
        return x_n_obs
    return obs_index_n


def read_obs_index(f: h5py.File) -> np.ndarray:
    return read_categorical_or_plain(f["obs"][_obs_index_key(f)])


def read_obs_columns(f: h5py.File, columns) -> pd.DataFrame:
    obs = f["obs"]
    data = {}
    for col in columns:
        if col not in obs:
            raise ValidationError(f"obs/{col} not present")
        data[col] = read_categorical_or_plain(obs[col])
    return pd.DataFrame(data, index=read_obs_index(f))


def read_var_names(f: h5py.File) -> np.ndarray:
    var = f["var"]
    index_key = _decode_bytes(var.attrs.get("_index", "_index"))
    return read_categorical_or_plain(var[index_key])


def read_uns_scalar(f: h5py.File, key: str) -> str:
    node = f["uns"][key]
    value = node[()] if isinstance(node, h5py.Dataset) else None
    return _decode_bytes(value)


def _check_csr_structure(group: h5py.Group, n_obs: int, n_vars: int) -> np.ndarray:
    """Structural pre-check using ONLY `.shape`/`.dtype`/cheap scalar reads
    -- never slices `data`/`indices` by content length. An HDF5 slice past
    an array's real end can silently truncate rather than raise, so the
    *returned length* of a later slice is not itself proof the array was
    as long as assumed; this checks the real on-disk shapes/dtypes first.

    v3 additions (20261011_113512.md F4): the group's own declared `shape`
    attribute is checked against `(n_obs, n_vars)` -- a corrupted column
    count would otherwise pass unnoticed -- and `indptr`/`indices` are
    checked for integer dtype and 1-D rank before any content is read."""
    enc = group.attrs.get("encoding-type")
    if enc != "csr_matrix":
        raise ValidationError(f"{group.name}: expected encoding-type 'csr_matrix', got {enc!r}")

    declared_shape = tuple(int(v) for v in group.attrs.get("shape", ()))
    if declared_shape != (n_obs, n_vars):
        raise ValidationError(
            f"{group.name}: attrs['shape']={declared_shape} != expected (n_obs, n_vars)="
            f"({n_obs}, {n_vars})"
        )

    for name in ("indptr", "indices"):
        ds = group[name]
        if ds.ndim != 1:
            raise ValidationError(f"{group.name}/{name}: expected rank 1, got shape {ds.shape}")
        if ds.dtype.kind not in ("i", "u"):
            raise ValidationError(f"{group.name}/{name}: expected an integer dtype, got {ds.dtype}")

    indptr_len = group["indptr"].shape[0]
    if indptr_len != n_obs + 1:
        raise ValidationError(
            f"{group.name}: indptr length {indptr_len} != n_obs+1 ({n_obs + 1})"
        )
    indptr = group["indptr"][:]  # length n_obs+1: a small int array, not the data array
    if indptr[0] != 0:
        raise ValidationError(f"{group.name}: indptr[0] = {indptr[0]} != 0")
    if not np.all(np.diff(indptr) >= 0):
        raise ValidationError(f"{group.name}: indptr is not monotonically non-decreasing")
    data_len = group["data"].shape[0]
    indices_len = group["indices"].shape[0]
    if not (indptr[-1] == data_len == indices_len):
        raise ValidationError(
            f"{group.name}: indptr[-1]={indptr[-1]} but data has {data_len} entries and "
            f"indices has {indices_len} -- these must all agree."
        )
    return indptr


def metadata_checks(h5ad_path: Path, manifest_path: Path, genes_csv: Path | None,
                     expected_stage_id: str = EXPECTED_STAGE_ID_DEFAULT,
                     expected_expression_source: str = EXPECTED_EXPRESSION_SOURCE_DEFAULT) -> dict:
    report: dict = {"scope": "metadata (all cells, obs/var/uns only -- no layer/X read)"}
    with h5py.File(h5ad_path, "r") as f:
        n_obs = _n_obs_from_metadata_only(f)
        if n_obs > MAX_OBS_ROWS_FOR_METADATA_READ:
            raise ValidationError(
                f"n_obs={n_obs} exceeds MAX_OBS_ROWS_FOR_METADATA_READ="
                f"{MAX_OBS_ROWS_FOR_METADATA_READ} -- login-node safety ceiling; "
                f"re-run on a compute node with a raised limit if genuinely intended."
            )
        var_names = read_var_names(f)
        n_vars = len(var_names)
        report["n_genes"] = int(n_vars)

        for layer_name in ["X"] + REQUIRED_LAYERS:
            grp = f["X"] if layer_name == "X" else f["layers"].get(layer_name)
            if grp is None:
                raise ValidationError(f"layers/{layer_name} missing")
            _check_csr_structure(grp, n_obs, n_vars)

        obs_names = read_obs_index(f)
        n_obs = len(obs_names)
        report["n_obs"] = int(n_obs)
        report["cell_id_unique"] = bool(len(set(obs_names)) == n_obs)

        obs_df = read_obs_columns(f, REQUIRED_OBS_COLUMNS)
        obs_df.index = obs_names

        uns = {}
        for key in REQUIRED_UNS_KEYS:
            if key not in f["uns"]:
                raise ValidationError(f"uns/{key} missing")
            uns[key] = read_uns_scalar(f, key)
        report["uns"] = uns
        report["stage_id_matches_expected"] = bool(uns["population_stage_id"] == expected_stage_id)
        report["expression_source_matches_expected"] = bool(
            uns["expression_source"] == expected_expression_source
        )
        report["layers_present"] = list(f["layers"].keys())

    manifest = pd.read_csv(manifest_path)
    report["manifest_rows"] = int(len(manifest))
    report["manifest_rows_match_n_obs"] = bool(len(manifest) == n_obs)

    order_matches = bool(
        len(manifest) == n_obs
        and (manifest["cell_id"].astype(str).to_numpy() == np.asarray(obs_names, dtype=str)).all()
    )
    report["cell_id_order_matches_selection_manifest"] = order_matches

    manifest_fields_match = "UNVERIFIED (cell_id order mismatch; row-wise field comparison skipped)"
    if order_matches:
        checks = []
        for field in ["cohort", "donor_id", "harmonized_lineage"]:
            checks.append((manifest[field].astype(str).to_numpy()
                            == obs_df[field if field != "harmonized_lineage" else "lineage"]
                            .astype(str).to_numpy()).all())
        manifest_fields_match = bool(all(checks))
    report["manifest_cohort_donor_lineage_match_obs_independent"] = manifest_fields_match

    patient_id_ok = bool(
        (obs_df["patient_id"].astype(str)
         == (obs_df["cohort"].astype(str) + "::" + obs_df["donor_id"].astype(str))).all()
    )
    report["patient_id_cohort_qualified_independent"] = patient_id_ok

    n_platelet = int((obs_df["lineage"] == "Platelet").sum())
    n_monocyte = int((obs_df["lineage"] == "Monocyte_Myeloid").sum())
    report["n_platelet"] = n_platelet
    report["n_monocyte"] = n_monocyte
    report["per_lineage_counts_reconciled_independent"] = bool(n_platelet + n_monocyte == n_obs)

    if genes_csv is not None:
        genes = pd.read_csv(genes_csv)
        gene_col = "gene" if "gene" in genes.columns else genes.columns[0]
        expected = genes[gene_col].astype(str).tolist()
        report["var_order_matches_source_genes_csv"] = bool(
            len(expected) == len(var_names) and list(map(str, var_names)) == expected
        )
    else:
        report["var_order_matches_source_genes_csv"] = "UNVERIFIED (no --genes-csv given)"

    report["_obs_for_sampling"] = obs_df
    return report


def _select_sample_rows(obs_df: pd.DataFrame, cells_per_group: int, seed: int) -> dict:
    if cells_per_group <= 0:
        raise ValidationError(
            f"--cells-per-group must be >= 1, got {cells_per_group} -- an empty sample "
            f"would make every sampled-numeric check vacuously True."
        )
    chosen = {}
    positions = np.arange(len(obs_df))
    for (cohort, lineage), _group in obs_df.groupby(["cohort", "lineage"]):
        idx = positions[(obs_df["cohort"].to_numpy() == cohort) & (obs_df["lineage"].to_numpy() == lineage)]
        rng = np.random.RandomState(derive_seed(seed, "job1_validate_sample", cohort, lineage))
        n_take = min(len(idx), cells_per_group)
        chosen[(cohort, lineage)] = sorted(rng.choice(idx, size=n_take, replace=False).tolist())
    return chosen


def _read_csr_rows(group: h5py.Group, row_positions, indptr: np.ndarray) -> dict:
    """`indptr` must already be the real, structurally-checked array from
    `_check_csr_structure`. Enforces a per-row and cumulative nnz budget
    BEFORE any `data`/`indices` slice is requested -- corrupted/oversized
    metadata is refused rather than silently read."""
    out = {}
    total = 0
    for pos in row_positions:
        lo, hi = int(indptr[pos]), int(indptr[pos + 1])
        n = hi - lo
        if n > MAX_NNZ_PER_SAMPLED_ROW:
            raise ValidationError(
                f"{group.name} row {pos}: {n} nnz exceeds MAX_NNZ_PER_SAMPLED_ROW="
                f"{MAX_NNZ_PER_SAMPLED_ROW} -- refusing to read before slicing."
            )
        total += n
        if total > MAX_TOTAL_SAMPLED_NNZ:
            raise ValidationError(
                f"{group.name}: cumulative sampled nnz {total} exceeds "
                f"MAX_TOTAL_SAMPLED_NNZ={MAX_TOTAL_SAMPLED_NNZ} -- refusing further reads."
            )
        out[pos] = (group["indices"][lo:hi], group["data"][lo:hi], n)
    return out


def sampled_numeric_checks(h5ad_path: Path, obs_df: pd.DataFrame, cells_per_group: int, seed: int) -> dict:
    report: dict = {
        "scope": "sampled numeric (fixed, deterministic, outcome-blind subset of cells only)",
        "cells_per_group": cells_per_group,
        "seed": seed,
    }
    sample = _select_sample_rows(obs_df, cells_per_group, seed)
    all_positions = sorted({p for positions in sample.values() for p in positions})
    sampled_cell_ids = {int(p): str(obs_df.index[p]) for p in all_positions}
    report["sampled_cell_ids_by_group"] = {
        f"{cohort}::{lineage}": [sampled_cell_ids[p] for p in positions]
        for (cohort, lineage), positions in sample.items()
    }

    group_entries = {}
    with h5py.File(h5ad_path, "r") as f:
        n_obs = _n_obs_from_metadata_only(f)
        n_genes_total = int(len(read_var_names(f)))
        for name in ["X"] + REQUIRED_LAYERS:
            grp = f["X"] if name == "X" else f["layers"][name]
            indptr = _check_csr_structure(grp, n_obs, n_genes_total)
            group_entries[name] = _read_csr_rows(grp, all_positions, indptr)

    nnz_read = {name: sum(n for _i, _v, n in rows.values()) for name, rows in group_entries.items()}
    report["nnz_entries_read_per_array"] = nnz_read
    report["approx_bytes_read_per_array"] = {
        name: int(sum(
            n * (group_entries[name][p][0].dtype.itemsize + group_entries[name][p][1].dtype.itemsize)
            for p, (_i, _v, n) in group_entries[name].items()
        ))
        for name in group_entries
    }

    indices_in_bounds = True
    values_finite_nonneg = True
    x_equals_counts_decontx = True
    counts_equals_counts_decontx = True
    for pos in all_positions:
        for name in ["X", "counts_raw", "counts_decontx", "counts"]:
            idx_arr, val_arr, _n = group_entries[name][pos]
            if idx_arr.size and (idx_arr.min() < 0 or idx_arr.max() >= n_genes_total):
                indices_in_bounds = False
            if not np.all(np.isfinite(val_arr)) or np.any(val_arr < 0):
                values_finite_nonneg = False
        x_idx, x_val, _ = group_entries["X"][pos]
        d_idx, d_val, _ = group_entries["counts_decontx"][pos]
        if not (np.array_equal(x_idx, d_idx) and np.array_equal(x_val, d_val)):
            x_equals_counts_decontx = False
        c_idx, c_val, _ = group_entries["counts"][pos]
        if not (np.array_equal(c_idx, d_idx) and np.array_equal(c_val, d_val)):
            counts_equals_counts_decontx = False

    report["sampled_indices_within_gene_bounds"] = bool(indices_in_bounds)
    report["sampled_values_finite_and_nonnegative"] = bool(values_finite_nonneg)
    report["sampled_X_equals_counts_decontx_exact"] = bool(x_equals_counts_decontx)
    report["sampled_counts_layer_equals_counts_decontx_exact"] = bool(counts_equals_counts_decontx)
    report["note_storage_dtype"] = (
        "Values are stored as float32; exact np.array_equal is used for the layer-copy "
        "checks above (not a tolerance) because the producer writes these as literal "
        "duplicates of the same in-memory array. 'non-negative' is checked exactly. "
        "See source_mtx_full_numeric_equality below for the boundary this does not cross."
    )
    report["source_mtx_full_numeric_equality"] = (
        "NOT_VERIFIED -- requires a full streamed re-read of raw/counts.mtx and "
        "decontx/counts.mtx against every selected column; out of scope for this tool."
    )
    return report


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", required=True, type=Path)
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--genes-csv", type=Path, default=None)
    p.add_argument("--cells-per-group", type=int, default=2)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--expected-stage-id", default=EXPECTED_STAGE_ID_DEFAULT)
    p.add_argument("--expected-expression-source", default=EXPECTED_EXPRESSION_SOURCE_DEFAULT)
    p.add_argument("--out", type=Path, default=None)
    args = p.parse_args()

    out_fd = None
    if args.out is not None:
        try:
            out_fd = os.open(str(args.out), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except OSError as exc:
            if exc.errno == errno.EEXIST:
                raise FileExistsError(
                    f"{args.out} already exists -- refusing to overwrite a prior report; "
                    f"use a fresh, newly timestamped --out path."
                ) from exc
            raise

    meta = metadata_checks(args.h5ad, args.manifest, args.genes_csv,
                            args.expected_stage_id, args.expected_expression_source)
    obs_df = meta.pop("_obs_for_sampling")
    numeric = sampled_numeric_checks(args.h5ad, obs_df, args.cells_per_group, args.seed)

    provenance = {
        "script": str(Path(__file__).resolve()),
        "script_git_revision": _script_git_revision(),
        "argv": sys.argv,
        "h5ad_sha256": _file_sha256(args.h5ad),
        "manifest_sha256": _file_sha256(args.manifest),
    }
    report = {"provenance": provenance, "metadata_checks": meta, "sampled_numeric_checks": numeric}
    text = json.dumps(report, indent=2, default=str)
    if out_fd is not None:
        with os.fdopen(out_fd, "w") as fh:
            fh.write(text)
        print(f"Wrote report to {args.out}")
    else:
        print(text)

    boolean_check_keys = [
        "cell_id_unique", "manifest_rows_match_n_obs",
        "cell_id_order_matches_selection_manifest",
        "manifest_cohort_donor_lineage_match_obs_independent",
        "patient_id_cohort_qualified_independent",
        "per_lineage_counts_reconciled_independent",
        "var_order_matches_source_genes_csv",
        "stage_id_matches_expected", "expression_source_matches_expected",
    ]
    numeric_check_keys = [
        "sampled_indices_within_gene_bounds", "sampled_values_finite_and_nonnegative",
        "sampled_X_equals_counts_decontx_exact", "sampled_counts_layer_equals_counts_decontx_exact",
    ]
    failed = [k for k in boolean_check_keys if meta.get(k) is False]
    failed += [k for k in numeric_check_keys if numeric.get(k) is False]
    if failed:
        print(f"FAILED checks: {failed}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
