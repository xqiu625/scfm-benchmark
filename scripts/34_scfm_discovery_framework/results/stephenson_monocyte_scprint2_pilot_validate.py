"""Step 3 of 3 in the Stephenson monocyte scPRINT-2 technical-pilot execution packet
(`STEPHENSON_MONOCYTE_SCPRINT2_PILOT_EXECUTION_PACKET_v4.md`): checks the ACTUAL
`scprint2 embed` output, not just that a file with some name exists and has roughly
the right cell count.

v2 -> v3 fix (per review:20261016_152556, "剩余阻断 1"): the v2 version
only took `expected_n_obs` -- an integer -- so an output with the right COUNT but the
WRONG cells (unrelated ids, duplicated ids, misrouted QC annotations) would pass. This
version takes `expected_obs`, the frozen prepared-input object's real `.obs` (the same
object `scprint2 embed` was given), and checks coverage/duplicates/extra-ids/order by
cell ID, never by position or count alone. It also fixes two further gaps the same
review found: `find_embedding_keys()`'s prefix match previously let an `X_scprint_umap`
-only file (no primary embedding) pass -- the installed producer
(`scprint2/tasks/cell_emb.py:214,227-251`, read directly 2026-10-16) writes the UMAP to
`X_scprint_umap` and the actual primary cell embedding to `scprint_emb`; this version
gates PASS/FAIL on `scprint_emb` specifically. And the embedding-value check previously
skipped any non-float/complex dtype (so a string or all-zero-column embedding could
pass); this version requires a real-valued, non-empty, 2-D array.

v3 -> v4 fix (per review:20261016_202903, "R2"+"R3"): v3's
`_qc_column_report` set `recovered_from_expected_manifest=True`/`n_recovered` in the
JSON report when a required column was missing but every cell id matched -- but
`main()` never actually wrote a recovered column back into any output or sidecar file,
so a real missing-column case produced a technical PASS whose "recovered" claim was
pure fiction (the review's own `persisted_missing_qc` probe reopened the output file
and found both columns still absent). That fake-success branch is removed: a required
column missing from the real output is now ALWAYS a hard failure, independent of
whether the cell-id manifest otherwise matches -- "recoverable in principle" is not the
same as "recovered", and this validator does not write or persist anything, so it must
not claim either. The same review's independent probes (`swapped_patient_id`,
`swapped_Sample_ID`, `swapped_visit`) also showed the old `DEFAULT_QC_OBS_COLUMNS`
checked only the 2 QC columns, so a present-but-misrouted donor/sample identity column
produced `hard_failures=[]` even with a correct primary embedding and correct QC
values. `stephenson_monocyte_raw.h5ad`'s real obs columns were read directly via a
bounded h5py read on 2026-10-16 (`['Data_NO', 'Sample_ID', '_index', 'disease',
'monocyte_source_provenance', 'orig_celltype', 'patient_id', 'source_study', 'y']`) --
there is no `visit` column on this object (the review's own `visit` probe field name
is a synthetic stand-in, not a claim that the real object has it); `patient_id` and
`Sample_ID` are real and are now checked by id exactly like the QC columns, by the
same (now simplified) required-column logic.

Scope: technical validation only (traceability, identity, numeric finiteness) -- NOT a
scientific-quality judgment. "Looks like good clusters" is explicitly never a criterion
checked here. Missing or extra cell ids are NOT treated as an acceptable partial
success under this round's authorized scope -- Cindy has not authorized any cell
deletion or partial-coverage acceptance for this pilot.
"""
import json
import os
import sys

import numpy as np
import pandas as pd

PRIMARY_EMBEDDING_KEY = "scprint_emb"
AUXILIARY_EMBEDDING_KEY_PREFIXES = ("X_scprint", "scprint_emb")

# QC annotations (prep-step-attached) AND donor/sample identity columns (already
# present on the raw object, carried through prep unchanged) -- both are required by
# id, with no recovery branch: see the v3 -> v4 fix note above.
DEFAULT_REQUIRED_OBS_COLUMNS = (
    "gate1_flag_would_delete",
    "decontx_estimated_contamination",
    "patient_id",
    "Sample_ID",
)

_DIAGNOSTICS_SAMPLE_SIZE = 50


def find_embedding_keys(obsm_keys):
    """Informational only -- any obsm key that LOOKS like an scPRINT-produced key
    (prefix match, e.g. `X_scprint_umap`, `scprint_emb`). Never used to decide
    pass/fail: `X_scprint_umap` is a derived UMAP, not the primary embedding, and a
    file that has only a UMAP must not pass as if it had the real embedding -- see
    PRIMARY_EMBEDDING_KEY, checked directly in `validate_embedding_output`."""
    return [k for k in obsm_keys if any(k.startswith(p) for p in AUXILIARY_EMBEDDING_KEY_PREFIXES)]


def _primary_embedding_failures(arr, n_obs_expected_rows):
    """Returns (failures, n_non_finite_or_None). Requires a real-valued (int/float),
    non-empty, 2-D array with one row per cell -- a string array, a zero-column
    array, or a row-count mismatch are all hard failures, not silently skipped."""
    failures = []
    arr = np.asarray(arr)
    if arr.ndim != 2 or arr.shape[1] == 0:
        failures.append(
            f"obsm[{PRIMARY_EMBEDDING_KEY!r}] has invalid shape {arr.shape} "
            f"(expected 2-D with at least 1 column)"
        )
        return failures, None
    if arr.shape[0] != n_obs_expected_rows:
        failures.append(
            f"obsm[{PRIMARY_EMBEDDING_KEY!r}] has {arr.shape[0]} rows, which does not "
            f"match {n_obs_expected_rows} cells in obs"
        )
    if arr.dtype.kind not in "fiu":
        failures.append(
            f"obsm[{PRIMARY_EMBEDDING_KEY!r}] has non-numeric dtype {arr.dtype} "
            f"(expected a real-valued embedding)"
        )
        return failures, None
    n_non_finite = int((~np.isfinite(arr.astype(np.float64))).sum())
    if n_non_finite > 0:
        failures.append(f"obsm[{PRIMARY_EMBEDDING_KEY!r}] contains {n_non_finite} non-finite value(s)")
    return failures, n_non_finite


def _values_diverge_by_id(actual_vals, expected_vals):
    """Compares two same-length, id-aligned Series. Boolean columns compare exactly;
    numeric columns compare with a small float tolerance; anything else (e.g. the
    string-valued `patient_id`/`Sample_ID` identity columns) compares as exact string
    equality -- there is no meaningful "close enough" for an identity label."""
    if actual_vals.dtype == bool or expected_vals.dtype == bool:
        return int((actual_vals.to_numpy() != expected_vals.to_numpy()).sum())
    if pd.api.types.is_numeric_dtype(actual_vals) and pd.api.types.is_numeric_dtype(expected_vals):
        return int(
            (~np.isclose(
                actual_vals.to_numpy(dtype=float),
                expected_vals.to_numpy(dtype=float),
                atol=1e-6,
                equal_nan=False,
            )).sum()
        )
    return int((actual_vals.astype(str).to_numpy() != expected_vals.astype(str).to_numpy()).sum())


def _required_obs_column_report(embedded_obs, expected_obs, common_ids, col, hard_failures):
    """Checks one required obs column (a QC annotation or a donor/sample identity
    column) by cell id against the frozen prepared-input manifest.

    Fail-closed, no recovery branch (per review:20261016_202903 "R2"):
    a missing required column is ALWAYS a hard failure. The prior version's "recovered
    from the expected manifest" claim described what COULD be reconstructed, not
    anything this function or its caller actually writes anywhere -- that review
    proved a real missing-column case opened a technical PASS whose recovered columns
    did not exist in the reopened output. If a future round adds an actual persisted
    recovery (a sidecar file, written and read back), that is a different, additional
    code path with its own round-trip test -- not a report-only boolean.
    """
    present = col in embedded_obs.columns
    entry = {"present_in_output": present}

    if col not in expected_obs.columns:
        entry["note"] = "column not present in the expected (prepared-input) manifest either"
        return entry

    if not present:
        hard_failures.append(
            f"required obs column {col!r} is missing from the output -- fail closed "
            f"(no recovery is claimed or performed by this validator)"
        )
        return entry

    if common_ids:
        actual_vals = embedded_obs.loc[common_ids, col]
        expected_vals = expected_obs.loc[common_ids, col]
        n_diverging = _values_diverge_by_id(actual_vals, expected_vals)
        entry["n_values_diverging_from_expected_manifest"] = n_diverging
        if n_diverging:
            hard_failures.append(
                f"obs[{col!r}] diverges from the frozen prepared-input manifest for "
                f"{n_diverging} cell id(s) -- misrouted annotation (checked by id, not position)"
            )
    return entry


def validate_embedding_output(embedded_adata, expected_obs, required_obs_columns=DEFAULT_REQUIRED_OBS_COLUMNS):
    """Pure function: `embedded_adata` is the real `scprint2 embed` output (an opened
    AnnData); `expected_obs` is the frozen prepared-input object's real `.obs` (same
    object given to `scprint2 embed`, indexed by cell id) -- NOT merely an expected
    count. Returns (report, hard_failures); `hard_failures` non-empty means the output
    fails a HARD criterion and must not be treated as a technical PASS.

    Every identity/annotation check below is keyed by cell id, never by row position
    or count alone: a same-count output made of unrelated or duplicated ids, or one
    whose QC annotations silently swapped rows, is rejected exactly like a
    count mismatch is.
    """
    hard_failures = []

    actual_ids = list(embedded_adata.obs_names)
    actual_index = pd.Index(actual_ids)
    dup_mask = actual_index.duplicated()
    n_duplicates = int(dup_mask.sum())
    if n_duplicates:
        hard_failures.append(
            f"{n_duplicates} duplicate cell id(s) in the output "
            f"(e.g. {list(actual_index[dup_mask][:5])})"
        )

    expected_ids = list(expected_obs.index)
    expected_set, actual_set = set(expected_ids), set(actual_ids)
    missing_ids = sorted(expected_set - actual_set)
    extra_ids = sorted(actual_set - expected_set)

    if extra_ids:
        hard_failures.append(
            f"{len(extra_ids)} cell id(s) in the output are not in the expected (prepared-input) "
            f"manifest -- unrelated or corrupted ids (e.g. {extra_ids[:5]})"
        )
    if missing_ids:
        hard_failures.append(
            f"{len(missing_ids)} expected cell id(s) are missing from the output -- this round's "
            f"authorized scope does not accept a partial embedding as a technical PASS "
            f"(e.g. {missing_ids[:5]})"
        )

    ids_match_as_sets = not missing_ids and not extra_ids
    order_matches_expected = ids_match_as_sets and actual_ids == expected_ids

    n_obs_actual = int(embedded_adata.n_obs)
    if n_obs_actual == 0:
        hard_failures.append("output has 0 cells")

    all_recognized_embedding_keys = find_embedding_keys(list(embedded_adata.obsm.keys()))
    non_finite_by_key = {}
    if PRIMARY_EMBEDDING_KEY not in embedded_adata.obsm:
        hard_failures.append(
            f"primary embedding key {PRIMARY_EMBEDDING_KEY!r} not found in obsm -- obsm keys "
            f"were: {list(embedded_adata.obsm.keys())}"
        )
    else:
        failures, n_non_finite = _primary_embedding_failures(
            embedded_adata.obsm[PRIMARY_EMBEDDING_KEY], n_obs_actual
        )
        hard_failures.extend(failures)
        if n_non_finite is not None:
            non_finite_by_key[PRIMARY_EMBEDDING_KEY] = n_non_finite

    # Required obs columns (QC annotations + donor/sample identity), recomputed by
    # cell id (never by position) against the frozen prepared-input manifest -- a
    # missing one is always a hard failure (no recovery claim, see the function
    # docstring); a present-but-diverging one is a hard failure too (misrouted
    # annotation). Skipped when the output has duplicate ids: a duplicate index makes
    # `.loc` by id fan out (pandas cross-matches every occurrence), so a by-id
    # comparison would be meaningless -- the duplicate itself is already a hard
    # failure above.
    common_ids = [i for i in actual_ids if i in expected_set] if n_duplicates == 0 else []
    required_columns_report = {
        col: _required_obs_column_report(embedded_adata.obs, expected_obs, common_ids, col, hard_failures)
        for col in required_obs_columns
    }

    report = {
        "primary_embedding_key": PRIMARY_EMBEDDING_KEY,
        "embedding_keys_found_all": all_recognized_embedding_keys,
        "n_obs_expected": len(expected_ids),
        "n_obs_actual": n_obs_actual,
        "n_obs_difference": len(expected_ids) - n_obs_actual,
        "n_missing_cell_ids": len(missing_ids),
        "n_extra_cell_ids": len(extra_ids),
        "n_duplicate_cell_ids_in_output": n_duplicates,
        "missing_cell_ids_sample": missing_ids[:_DIAGNOSTICS_SAMPLE_SIZE],
        "extra_cell_ids_sample": extra_ids[:_DIAGNOSTICS_SAMPLE_SIZE],
        "ids_match_as_sets": ids_match_as_sets,
        "order_matches_expected": order_matches_expected,
        "non_finite_values_by_embedding_key": non_finite_by_key,
        "required_obs_columns": required_columns_report,
    }
    return report, hard_failures


def main(embedded_h5ad_path, expected_input_h5ad_path, out_json, required_obs_columns=DEFAULT_REQUIRED_OBS_COLUMNS):
    import anndata as ad

    if not os.path.exists(embedded_h5ad_path):
        print(f"FAIL: embedding output not found: {embedded_h5ad_path}", file=sys.stderr)
        return 1
    if not os.path.exists(expected_input_h5ad_path):
        print(f"FAIL: expected (prepared-input) manifest not found: {expected_input_h5ad_path}", file=sys.stderr)
        return 1

    # backed="r": obs/obsm are loaded eagerly (small -- what we need), X stays
    # lazy/unread -- this script never needs gene expression values, and the
    # prepared/embedded objects can be large.
    embedded = ad.read_h5ad(embedded_h5ad_path, backed="r")
    expected_obs = ad.read_h5ad(expected_input_h5ad_path, backed="r").obs

    report, hard_failures = validate_embedding_output(embedded, expected_obs, required_obs_columns)
    report["embedded_h5ad_path"] = embedded_h5ad_path
    report["expected_input_h5ad_path"] = expected_input_h5ad_path
    report["hard_failures"] = hard_failures

    with open(out_json, "w") as f:
        json.dump(report, f, indent=2, sort_keys=True, default=str)
        f.write("\n")

    print(json.dumps(report, indent=2, sort_keys=True, default=str))
    if hard_failures:
        print(f"FAIL: {len(hard_failures)} hard failure(s): {hard_failures}", file=sys.stderr)
        return 1
    print("PASS (technical criteria only -- not a scientific-quality judgment)")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(
            "usage: stephenson_monocyte_scprint2_pilot_validate.py "
            "<embedded.h5ad> <expected_input.h5ad> <out.json>",
            file=sys.stderr,
        )
        sys.exit(2)
    sys.exit(main(sys.argv[1], sys.argv[2], sys.argv[3]))
