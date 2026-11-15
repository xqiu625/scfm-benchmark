"""Identity-recovery adapter for the Stephenson monocyte scPRINT-2 small-v2 pilot
(job 30044888, FAILED on `stephenson_monocyte_scprint2_pilot_validate.py`'s hard
identity check -- see `implementation-note:20261017_164353` and
`review:20261017_213117`).

## What this fixes and what it does not

The installed `scdataloader==2.1.2` `Preprocessor.__call__` (`preprocess.py:424`)
unconditionally overwrites `adata.obs_names` with fresh UUIDs before the model ever
sees the data. This destroys the ORIGINAL `obs_names` label but never touches any
OTHER obs column -- in particular, `stephenson_monocyte_scprint2_pilot_prep.py`
writes `obs["cell_id_stripped"]` (the barcode-donor suffix with the `steph_` prefix
already removed) onto the same dataframe, and that column survives the entire
`scprint2 embed` round trip untouched: read `scprint2/tasks/cell_emb.py` end to end
(installed copy, 2026-10-17) and the ORIGINAL `adata.obs` (containing
`cell_id_stripped`) is never reordered or rebuilt -- only `pred_adata.obs` (the
per-shard round-tripped frame) is index-reassigned and then concatenated onto it as
EXTRA COLUMNS via `pd.concat([adata.obs, pred_adata.obs], axis=1)`. So recovering
`obs_names` from `cell_id_stripped` is a pure relabeling of rows that were never
reordered -- it does not depend on, and cannot by itself verify, whether the
EMBEDDING MATRIX (`obsm["scprint_emb"]`, built from the round-tripped `pred_adata`)
is row-aligned with those now-correctly-labeled rows. That is a separate question,
addressed by `stephenson_monocyte_scprint2_pilot_row_order_smoke_test.py`
(a real-pipeline permutation-invariance test, not run this round -- it needs a new
SLURM allocation that this round does not submit). This adapter and its tests
therefore CANNOT and do NOT claim to detect a fault that swaps two same-sample rows'
embedding values while leaving every obs column (including `cell_id_stripped`)
exactly in place -- see the dedicated test for that scope limit, asserted directly
rather than left implicit.

This script writes to a NEW path and never modifies the frozen FAILED artifacts.
"""
import hashlib
import json
import os
import shutil
import sys
import time

import h5py
import numpy as np
import pandas as pd

CELL_ID_PREFIX = "steph_"
STRIPPED_COL = "cell_id_stripped"
STRING_ARRAY_ENCODING_ATTRS = {"encoding-type": "string-array", "encoding-version": "0.2.0"}


class RecoveryRefused(Exception):
    """Raised whenever the recovered identity cannot be trusted. Never caught and
    silently downgraded -- every refusal here is a hard failure for the caller."""


def recover_cell_ids(obs: pd.DataFrame, prefix: str = CELL_ID_PREFIX, stripped_col: str = STRIPPED_COL) -> pd.Index:
    """Pure function: `obs` -> recovered `pd.Index`, same length and row order as
    `obs`. Refuses (does not silently drop/dedupe/default) on any of: missing
    column, null values, non-string values, empty-string values, or duplicate
    recovered ids."""
    if stripped_col not in obs.columns:
        raise RecoveryRefused(f"obs has no {stripped_col!r} column -- nothing to recover from")

    col = obs[stripped_col]
    if len(col) == 0:
        raise RecoveryRefused("obs has 0 rows -- nothing to recover")

    n_null = int(col.isna().sum())
    if n_null:
        raise RecoveryRefused(f"{n_null} null value(s) in obs[{stripped_col!r}] -- refusing a partial recovery")

    non_str = [v for v in col if not isinstance(v, str)]
    if non_str:
        raise RecoveryRefused(
            f"{len(non_str)} non-string value(s) in obs[{stripped_col!r}] "
            f"(e.g. {non_str[:3]!r}) -- refusing"
        )

    n_empty = int((col == "").sum())
    if n_empty:
        raise RecoveryRefused(f"{n_empty} empty-string value(s) in obs[{stripped_col!r}] -- refusing")

    recovered = pd.Index([prefix + v for v in col])
    dup_mask = recovered.duplicated()
    if dup_mask.any():
        raise RecoveryRefused(
            f"{int(dup_mask.sum())} duplicate recovered id(s) "
            f"(e.g. {list(recovered[dup_mask][:5])}) -- refusing"
        )
    return recovered


def check_against_expected_manifest(
    recovered_ids: pd.Index,
    embedded_obs: pd.DataFrame,
    expected_obs: pd.DataFrame,
    required_columns=("patient_id", "Sample_ID", "gate1_flag_would_delete", "decontx_estimated_contamination"),
) -> dict:
    """Refuses unless the recovered id set is EXACTLY the expected manifest's id
    set (no missing, no extra), `expected_obs`/`recovered_ids` carry no duplicate
    ids, and every column in `required_columns` is PRESENT on both sides and
    agrees exactly, keyed by id, between the embedded output (relabeled) and the
    frozen prepared-input manifest. These columns are required, not optional --
    a column missing from either side is a hard refusal, never a silently
    skipped cross-check (see review:20261017_223541 R2: the prior
    version's `skipped` branch let a missing `patient_id` through as
    `status=recovered`). Returns a report dict on success; raises
    RecoveryRefused on any mismatch."""
    if expected_obs.index.duplicated().any():
        n_dup = int(expected_obs.index.duplicated().sum())
        raise RecoveryRefused(
            f"expected (prepared-input) manifest has {n_dup} duplicate id(s) -- "
            f"refusing rather than risk a pandas .loc broadcast across them"
        )

    expected_ids = expected_obs.index
    recovered_set, expected_set = set(recovered_ids), set(expected_ids)
    missing = expected_set - recovered_set
    extra = recovered_set - expected_set
    if missing or extra:
        raise RecoveryRefused(
            f"recovered id set does not exactly match the expected manifest: "
            f"{len(missing)} missing (e.g. {sorted(missing)[:3]}), "
            f"{len(extra)} extra (e.g. {sorted(extra)[:3]})"
        )

    relabeled_obs = embedded_obs.copy()
    relabeled_obs.index = recovered_ids

    missing_required = [
        col for col in required_columns
        if col not in relabeled_obs.columns or col not in expected_obs.columns
    ]
    if missing_required:
        raise RecoveryRefused(
            f"required column(s) {missing_required} not present on both the "
            f"recovered output and the expected manifest -- refusing, not skipping"
        )

    report = {"n_cells": len(recovered_ids), "cross_check": {}}
    for col in required_columns:
        actual = relabeled_obs.loc[expected_ids, col]
        expected = expected_obs.loc[expected_ids, col]
        if actual.dtype == bool or expected.dtype == bool:
            n_mismatch = int((actual.to_numpy() != expected.to_numpy()).sum())
        elif pd.api.types.is_numeric_dtype(actual) and pd.api.types.is_numeric_dtype(expected):
            n_mismatch = int((~np.isclose(actual.to_numpy(dtype=float), expected.to_numpy(dtype=float), atol=1e-6, equal_nan=False)).sum())
        else:
            n_mismatch = int((actual.astype(str).to_numpy() != expected.astype(str).to_numpy()).sum())
        report["cross_check"][col] = {"n_common_ids": len(expected_ids), "n_mismatching": n_mismatch}
        if n_mismatch:
            raise RecoveryRefused(
                f"obs[{col!r}] diverges from the frozen prepared-input manifest for "
                f"{n_mismatch} recovered id(s) -- refusing"
            )
    return report


def _read_obs_column_h5(group, key):
    x = group[key]
    if isinstance(x, h5py.Dataset):
        return x.asstr()[:] if x.dtype.kind in "OSU" else x[:]
    cats = x["categories"].asstr()[:]
    codes = x["codes"][:]
    n_cats = len(cats)
    out_of_range = [int(c) for c in codes if c >= n_cats]
    if out_of_range:
        raise RecoveryRefused(
            f"categorical column {key!r} has {len(out_of_range)} code(s) >= "
            f"{n_cats} categories (e.g. {out_of_range[:3]}) -- corrupt categorical, refusing"
        )
    return np.array([cats[c] if c >= 0 else None for c in codes], dtype=object)


def _read_obs_dataframe_h5(h5ad_path, columns):
    """Reads ONLY the named obs columns (plus the index) directly via h5py --
    never opens the file through anndata, never touches X/obsm/varm/layers. This
    is the login-node-safe path: a 5-6 GB h5ad's obs group is a few MB at most,
    independent of how large its expression matrix is."""
    with h5py.File(h5ad_path, "r") as f:
        obs_group = f["obs"]
        index_name = obs_group.attrs["_index"]
        ids = _read_obs_column_h5(obs_group, index_name)
        data = {}
        for col in columns:
            if col in obs_group:
                data[col] = _read_obs_column_h5(obs_group, col)
    return pd.DataFrame(data, index=list(ids))


REQUIRED_COLUMNS = ("patient_id", "Sample_ID", "gate1_flag_would_delete", "decontx_estimated_contamination")


def _sha256_file(path, chunk_size=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def _refuse_if_path_conflict(out_h5ad_path, *other_paths):
    """Cheap, EARLY, advisory check (not the race-proof guard -- see
    _exclusive_create_and_copy/_exclusive_write_json below for that): refuses if
    `out_h5ad_path` already exists, or resolves (realpath) to the same file as
    any input path, before this function does any of the expensive obs reading
    or validation work. Per review:20261017_223541 R1."""
    out_real = os.path.realpath(out_h5ad_path)
    for p in other_paths:
        if os.path.realpath(p) == out_real:
            raise RecoveryRefused(
                f"output path {out_h5ad_path!r} resolves to the same file as input path {p!r} -- refusing"
            )
    if os.path.exists(out_h5ad_path):
        raise RecoveryRefused(
            f"refusing to write a recovered output to an EXISTING path {out_h5ad_path!r} -- "
            f"a rerun must use a new, not-yet-existing output path"
        )


def _refuse_if_manifest_exists(manifest_path):
    """Early, advisory check (not the race-proof guard -- _exclusive_write_json
    is): avoids doing the real h5ad copy only to fail afterward on an
    already-taken manifest sidecar path."""
    if os.path.exists(manifest_path):
        raise RecoveryRefused(
            f"refusing: manifest sidecar {manifest_path!r} already exists -- "
            f"a rerun must use a new, not-yet-existing output path"
        )


def _exclusive_create_and_copy(src_path, dst_path):
    """Atomically claims `dst_path` (os.O_CREAT | os.O_EXCL) and only then copies
    `src_path`'s bytes into it -- closes the TOCTOU race
    review:20261018_000122's independent probe demonstrated: an earlier
    os.path.exists() check followed by shutil.copy2() lets a writer that creates
    dst_path IN BETWEEN those two steps get silently overwritten, because
    copy2() opens its destination non-exclusively. Raises RecoveryRefused
    (never a raw OSError) if dst_path already exists at the moment of the
    exclusive create -- including a destination created after any earlier
    check in this module. On any failure after claiming the path, removes our
    own (and only our own -- the exclusive create guarantees no one else could
    be writing here) partial file rather than leaving a half-written claim
    behind for a later, unrelated attempt to misread as "already exists"."""
    try:
        fd = os.open(dst_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
    except FileExistsError:
        raise RecoveryRefused(
            f"refusing to write to {dst_path!r}: it was created (by this call or another "
            f"writer) between the earlier existence check and this exclusive-create -- "
            f"a rerun must use a new, not-yet-existing output path"
        )
    try:
        with os.fdopen(fd, "wb") as dst, open(src_path, "rb") as src:
            shutil.copyfileobj(src, dst, length=1 << 20)
        shutil.copystat(src_path, dst_path)
    except Exception:
        try:
            os.remove(dst_path)
        except OSError:
            pass
        raise


def _exclusive_write_json(path, payload):
    """Same atomic-claim discipline as _exclusive_create_and_copy, for the
    manifest sidecar: review:20261018_000122's independent probe
    showed the prior `open(path, "w")` silently overwrote a pre-existing
    sidecar (a sentinel placed before the call was gone afterward). Raises
    RecoveryRefused, never leaves a partial file, on any conflict or failure."""
    try:
        fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
    except FileExistsError:
        raise RecoveryRefused(
            f"refusing to write manifest {path!r}: it already exists -- "
            f"a rerun must use a new, not-yet-existing output path"
        )
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(payload, f, indent=2, sort_keys=True)
            f.write("\n")
    except Exception:
        try:
            os.remove(path)
        except OSError:
            pass
        raise


def apply_recovery_to_h5ad(embedded_h5ad_path: str, expected_input_h5ad_path: str, out_h5ad_path: str) -> dict:
    """Surgical, X/obsm-untouched identity fix: copies `embedded_h5ad_path` to
    `out_h5ad_path` at the file level (never loads the expression matrix or the
    embedding array into memory), then patches ONLY `obs/_index` (in place, same
    length/order) and adds `obs/original_model_uuid` to the copy. The original
    FAILED artifact at `embedded_h5ad_path` is opened read-only and never modified.
    Every obs read in this function goes through `_read_obs_dataframe_h5`
    (h5py-only) -- login-node-safe even for multi-GB inputs, unlike an
    `anndata.read_h5ad(..., backed="r")` call, which still showed up as an
    OOM-kill (exit 137) on this cluster's login node when tried against this
    exact pair of files (2026-10-17) -- that attempt is not repeated here.
    Refuses up front (before touching `out_h5ad_path` at all) if that path
    already exists or aliases an input path -- never silently overwrites an
    existing result (review:20261017_223541 R1). Also refuses
    (raises RecoveryRefused) if identity validation fails; in that case nothing
    is written to `out_h5ad_path` at all (the path-conflict/validation checks
    all happen before `shutil.copy2`)."""
    _refuse_if_path_conflict(out_h5ad_path, embedded_h5ad_path, expected_input_h5ad_path)

    embedded_hash_before = _sha256_file(embedded_h5ad_path)

    expected_obs = _read_obs_dataframe_h5(expected_input_h5ad_path, REQUIRED_COLUMNS)

    embedded_obs_full = _read_obs_dataframe_h5(
        embedded_h5ad_path, (STRIPPED_COL,) + REQUIRED_COLUMNS
    )
    original_ids = list(embedded_obs_full.index)

    if STRIPPED_COL not in embedded_obs_full.columns:
        raise RecoveryRefused(f"{STRIPPED_COL!r} not present in {embedded_h5ad_path}")

    recovered_ids = recover_cell_ids(embedded_obs_full)

    cross_check_report = check_against_expected_manifest(recovered_ids, embedded_obs_full, expected_obs)

    manifest_path = out_h5ad_path + ".manifest.json"
    _refuse_if_manifest_exists(manifest_path)  # early, advisory -- avoids a wasted copy; see that function's docstring

    _exclusive_create_and_copy(embedded_h5ad_path, out_h5ad_path)
    with h5py.File(out_h5ad_path, "r+") as f:
        obs_group = f["obs"]
        index_name = obs_group.attrs["_index"]
        del obs_group[index_name]
        idx_ds = obs_group.create_dataset(
            index_name, data=np.array(list(recovered_ids), dtype=object), dtype=h5py.string_dtype()
        )
        idx_ds.attrs.update(STRING_ARRAY_ENCODING_ATTRS)
        if "original_model_uuid" in obs_group:
            del obs_group["original_model_uuid"]
        uuid_ds = obs_group.create_dataset(
            "original_model_uuid", data=np.array(list(original_ids), dtype=object), dtype=h5py.string_dtype()
        )
        uuid_ds.attrs.update(STRING_ARRAY_ENCODING_ATTRS)
        attrs_cols = obs_group.attrs.get("column-order")
        if attrs_cols is not None and "original_model_uuid" not in list(attrs_cols):
            obs_group.attrs["column-order"] = list(attrs_cols) + ["original_model_uuid"]

    embedded_hash_after = _sha256_file(embedded_h5ad_path)
    if embedded_hash_after != embedded_hash_before:
        raise RecoveryRefused(
            f"source artifact {embedded_h5ad_path!r} changed during recovery "
            f"({embedded_hash_before} -> {embedded_hash_after}) -- this should be impossible "
            f"and is reported rather than silently accepted"
        )

    manifest = {
        "schema": "stephenson_monocyte_scprint2_pilot_recover_cell_ids_manifest_v1",
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "code_sha256": _sha256_file(__file__),
        "embedded_h5ad_path": embedded_h5ad_path,
        "embedded_h5ad_sha256": embedded_hash_before,
        "expected_input_h5ad_path": expected_input_h5ad_path,
        "expected_input_h5ad_sha256": _sha256_file(expected_input_h5ad_path),
        "out_h5ad_path": out_h5ad_path,
        "out_h5ad_sha256": _sha256_file(out_h5ad_path),
        "n_cells": len(recovered_ids),
        "cell_id_prefix": CELL_ID_PREFIX,
        "mapping_rule": f"obs_names := {CELL_ID_PREFIX!r} + obs[{STRIPPED_COL!r}]; original uuid kept in obs['original_model_uuid']",
        "required_columns_cross_checked": list(REQUIRED_COLUMNS),
        "cross_check": cross_check_report["cross_check"],
    }
    _exclusive_write_json(manifest_path, manifest)

    return {
        "status": "recovered",
        "n_cells": len(recovered_ids),
        "out_h5ad_path": out_h5ad_path,
        "manifest_path": manifest_path,
        "cross_check": cross_check_report["cross_check"],
    }


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(
            "usage: stephenson_monocyte_scprint2_pilot_recover_cell_ids.py "
            "<embedded.h5ad> <expected_input.h5ad> <out_recovered.h5ad>",
            file=sys.stderr,
        )
        sys.exit(2)
    try:
        result = apply_recovery_to_h5ad(sys.argv[1], sys.argv[2], sys.argv[3])
    except RecoveryRefused as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        sys.exit(1)
    print(json.dumps(result, indent=2, default=str))
