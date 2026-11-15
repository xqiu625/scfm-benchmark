"""Input preparation for the Stephenson monocyte scPRINT-2 TECHNICAL pilot authorized by
Cindy's 2026-10-16 "允许" decision (docs/decisions.md, "Stephenson monocytes: provisional
raw counts accepted for a technical pilot"). This is step 1 of 2 in the execution packet
(`STEPHENSON_MONOCYTE_SCPRINT2_PILOT_EXECUTION_PACKET_v1.md`); step 2 is the `scprint2 embed`
CLI call documented there, not performed by this script.

Scope, exactly as authorized -- this module does each of these and nothing more:
  - Subsets the var (gene/feature) axis to `feature_types == "Gene Expression"` only,
    dropping the 192 "Antibody Capture" (CITE-seq) features mixed into
    `stephenson_monocyte_raw.h5ad`'s var (confirmed via a bounded h5py read, 2026-10-16;
    `build_stephenson_decontx_input.py` already does the identical GEX-only subset for the
    full-cohort DecontX input, for the same documented reason -- a scFM operates on gene
    expression, not antibody-capture features). Verified by direct comparison against
    `human_platelet_covid_severity_raw_counts.h5ad`, the ALREADY-WORKING scPRINT-2 input for
    a different cohort: that object's var carries ONLY `_index` (no `feature_types` at all),
    i.e. it was ALREADY GEX-only before scPRINT-2 ever saw it -- this script performs the
    equivalent step for Stephenson's still-mixed object.
  - Attaches Gate 1 QC flags and historical DecontX contamination scores as NEW obs columns,
    joined by the already-verified (100%, 0 mismatches, `stephenson_qc_source_crosscheck.py`)
    `cell_id` key -- ANNOTATION ONLY. Never filters a cell by these flags and never replaces
    X with decontaminated counts -- both explicitly excluded by the owner's approval scope.
  - Does NOT touch the obs (cell) axis at all: same 128,503 cells, same row order, same
    donor/sample/visit/cell identities as the source file -- verified by
    `validate_prepared_adata()` below.
  - Does NOT touch X's VALUES at all -- same raw integer counts, only a column (var) subset.

NOT run against the real 5.7 GB object this round: DEVELOPMENT.md's standing rule is "do not run
the real 5 GB input on the login node" (OOM risk), and no SLURM job is authorized this round.
`build_prepared_adata()` is implemented and tested (`tests/test_stephenson_monocyte_scprint2_
pilot_prep.py`) only against a small, in-memory synthetic fixture. `main()` wires real file
paths and would need to run inside a SLURM allocation (same resource class already proven
sufficient for reading this exact file in `stephenson_sample_subtype_cell_counts.py` is NOT
evidence here, since that script only ever read 2 small obs columns via bounded h5py slicing
-- this script's real run needs to materialize X, a genuinely different, larger operation).
"""
import csv
import hashlib
import json
import os

import numpy as np
import pandas as pd

RAW_H5AD = "/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c_monocyte/stephenson_monocyte_raw.h5ad"
GATE1_FLAGS_CSV = (
    "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/gate1_qc/stephenson/"
    "qc_flags_per_cell.csv"
)
DECONTX_SCORES_CSV = (
    "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/stephenson_decontx_input/"
    "stephenson_decontx_contamination_scores.csv"
)
OUT_H5AD = (
    "/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c_monocyte/"
    "stephenson_monocyte_scprint2_pilot_input_v1.h5ad"
)
CELL_ID_PREFIX = "steph_"


def _strip_prefix(obs_names, prefix=CELL_ID_PREFIX):
    out = []
    for name in obs_names:
        if not str(name).startswith(prefix):
            raise SystemExit(f"refusing: obs name {name!r} does not start with {prefix!r}")
        out.append(str(name)[len(prefix):])
    return np.array(out)


def _load_gate1_flags(path):
    rows = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            rows.append({
                "cell_id": row["cell_id"],
                "flag_would_delete": row["flag_would_delete"].strip().lower() == "true",
                "n_flags": int(row["n_flags"]),
            })
    df = pd.DataFrame(rows).set_index("cell_id")
    return df


def _load_decontx_scores(path):
    rows = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            rows.append({"cell_id": row["cell_id"], "decontX_contamination": float(row["decontX_contamination"])})
    df = pd.DataFrame(rows).set_index("cell_id")
    return df["decontX_contamination"]


def build_prepared_adata(adata, gate1_flags, decontx_scores):
    """Pure transform: `adata` (an in-memory AnnData) -> a NEW AnnData. Never mutates
    `adata` in place. `gate1_flags` is a DataFrame indexed by cell_id with columns
    flag_would_delete/n_flags; `decontx_scores` is a Series indexed by cell_id.

    Returns (prepared_adata, report). `report` states exactly what changed, including
    any cell missing a QC annotation (reported as pandas NA, never silently 0/False).
    """
    if "feature_types" not in adata.var.columns:
        raise SystemExit(
            "refusing: adata.var has no 'feature_types' column -- cannot confirm which "
            "features are Gene Expression vs Antibody Capture; this object may already be "
            "GEX-only (check directly) but this function refuses to assume that silently"
        )
    gex_mask = (adata.var["feature_types"].astype(str) == "Gene Expression").to_numpy()
    n_dropped = int((~gex_mask).sum())

    out = adata[:, gex_mask].copy()

    cell_ids = _strip_prefix(out.obs_names.to_numpy())
    flags = gate1_flags.reindex(cell_ids)
    scores = decontx_scores.reindex(cell_ids)

    out.obs = out.obs.copy()
    out.obs["cell_id_stripped"] = cell_ids
    out.obs["gate1_flag_would_delete"] = flags["flag_would_delete"].to_numpy()
    out.obs["gate1_n_flags"] = flags["n_flags"].to_numpy()
    out.obs["decontx_estimated_contamination"] = scores.to_numpy()

    report = {
        "n_cells_in": int(adata.n_obs),
        "n_cells_out": int(out.n_obs),
        "n_vars_in": int(adata.n_vars),
        "n_vars_out": int(out.n_vars),
        "n_vars_dropped_non_gene_expression": n_dropped,
        "n_cells_missing_gate1_flag": int(flags["flag_would_delete"].isna().sum()),
        "n_cells_missing_decontx_score": int(scores.isna().sum()),
        "obs_row_order_preserved": bool(list(out.obs_names) == list(adata.obs_names)),
    }
    return out, report


def _require(condition, message):
    """Raises AssertionError on a false condition, exactly like `assert cond, message`
    -- but NOT stripped under `python -O`/`PYTHONOPTIMIZE`, unlike a bare `assert`.

    Per review:20261016_152556 ('两个小问题' #1): the prior version of
    `validate_prepared_adata` used bare `assert` statements for every one of its
    acceptance checks, which `python -O` (or `PYTHONOPTIMIZE=1`) silently removes --
    so a validator meant to be the last line of defense before writing output could
    be made a no-op by the interpreter flag alone, with no change to this file.
    """
    if not condition:
        raise AssertionError(message)


def validate_prepared_adata(original_adata, prepared_adata, report):
    """Identity/numeric acceptance checks -- NOT a scientific-quality judgment. Raises
    AssertionError (never silently passes) if any of these fail. Survives
    `python -O`/`PYTHONOPTIMIZE` -- see `_require` above.

    R4 fix (per review:20261016_141242): the prior version trusted
    `report["obs_row_order_preserved"]`/`report["n_cells_missing_*"]` -- precomputed
    booleans/counts the CALLER supplied -- without re-deriving them from the actual
    `original_adata`/`prepared_adata` objects passed in. An independent probe proved
    this let a reversed-cell-order object and a corrupted-counts object both pass as
    long as the (stale) `report` said everything was fine. Every check below is now
    recomputed directly from the two real objects; `report` is accepted for the
    caller's convenience but never trusted.
    """
    import scipy.sparse as sp

    # 1. Row identity/order, recomputed directly -- not read from `report`.
    _require(prepared_adata.n_obs == original_adata.n_obs, "cell count changed")
    _require(
        list(prepared_adata.obs_names) == list(original_adata.obs_names),
        "cell identity/order differs between the original and prepared object",
    )

    # 2. GEX var axis: non-empty, and an EXACT, in-order match to the original's
    # Gene-Expression-only var names -- not merely "some genes remain".
    if "feature_types" not in original_adata.var.columns:
        raise SystemExit("refusing: original_adata.var has no 'feature_types' column")
    orig_gex_mask = (original_adata.var["feature_types"].astype(str) == "Gene Expression").to_numpy()
    orig_gex_names = list(np.asarray(original_adata.var_names)[orig_gex_mask])
    _require(len(orig_gex_names) > 0, "the original object has zero Gene Expression features")
    _require(prepared_adata.n_vars > 0, "prepared object has zero genes")
    _require(
        list(prepared_adata.var_names) == orig_gex_names,
        "prepared var axis does not exactly match the original's Gene-Expression-only "
        "var names, in the same order",
    )

    # 3. Original GEX counts preserved EXACTLY (not merely finite) -- compare against
    # the GEX-subset slice of the ORIGINAL X directly, without densifying the full
    # matrix (sparse difference stays bounded by nnz).
    orig_gex_X = original_adata[:, orig_gex_mask].X
    prep_X = prepared_adata.X
    if sp.issparse(orig_gex_X) or sp.issparse(prep_X):
        diff = sp.csr_matrix(orig_gex_X) - sp.csr_matrix(prep_X)
        _require(diff.nnz == 0 or not np.any(diff.data), "X values changed after the var subset")
    else:
        _require(
            np.array_equal(np.asarray(orig_gex_X), np.asarray(prep_X)),
            "X values changed after the var subset",
        )

    # 4. QC annotations: no missing, finite, and within their known-valid ranges --
    # recomputed from `prepared_adata.obs` directly, never from `report`'s counts.
    _require(
        prepared_adata.obs["gate1_flag_would_delete"].isna().sum() == 0,
        "at least one cell is missing a Gate 1 flag",
    )
    _require(
        prepared_adata.obs["decontx_estimated_contamination"].isna().sum() == 0,
        "at least one cell is missing a DecontX score",
    )
    n_flags = prepared_adata.obs["gate1_n_flags"].to_numpy(dtype=float)
    _require(
        np.all(np.isfinite(n_flags)) and np.all(n_flags >= 0),
        "gate1_n_flags contains a non-finite or negative value",
    )
    contam = prepared_adata.obs["decontx_estimated_contamination"].to_numpy(dtype=float)
    _require(np.all(np.isfinite(contam)), "decontx_estimated_contamination contains a non-finite value")
    _require(
        np.all((contam >= 0) & (contam <= 1)),
        "decontx_estimated_contamination outside its valid [0, 1] range",
    )

    return True


def _require_output_does_not_already_exist(out_h5ad, out_summary=None):
    """R5 fix (per review:20261016_141242): the prior version only
    checked `out_h5ad` -- an independent fixture proved that if a `_summary.json`
    already existed but the h5ad did not, `main()` would overwrite the summary.
    Both target paths are now checked before anything is written."""
    existing = [p for p in (out_h5ad, out_summary) if p and os.path.exists(p)]
    if existing:
        raise SystemExit(f"refusing to overwrite existing prepared-pilot-input output(s): {existing}")


def _sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(h5ad_path=RAW_H5AD, gate1_flags_path=GATE1_FLAGS_CSV, decontx_scores_path=DECONTX_SCORES_CSV,
         out_h5ad=OUT_H5AD):
    """NOT exercised against real data this round -- see module docstring. Requires
    `anndata` and enough memory to materialize X (a SLURM allocation, not the login node)."""
    import anndata as ad

    summary_path = out_h5ad.replace(".h5ad", "_summary.json")
    _require_output_does_not_already_exist(out_h5ad, summary_path)

    adata = ad.read_h5ad(h5ad_path)
    gate1_flags = _load_gate1_flags(gate1_flags_path)
    decontx_scores = _load_decontx_scores(decontx_scores_path)

    prepared, report = build_prepared_adata(adata, gate1_flags, decontx_scores)
    validate_prepared_adata(adata, prepared, report)

    prepared.write_h5ad(out_h5ad)

    report["source_files"] = {
        "h5ad_path": h5ad_path,
        "h5ad_sha256": _sha256_file(h5ad_path),
        "gate1_flags_csv": gate1_flags_path,
        "gate1_flags_csv_sha256": _sha256_file(gate1_flags_path),
        "decontx_scores_csv": decontx_scores_path,
        "decontx_scores_csv_sha256": _sha256_file(decontx_scores_path),
    }
    report["output_file"] = out_h5ad
    report["output_file_sha256"] = _sha256_file(out_h5ad)

    with open(summary_path, "w") as f:
        json.dump(report, f, indent=2, sort_keys=True)
        f.write("\n")

    print(json.dumps(report, indent=2, sort_keys=True))
    return report


if __name__ == "__main__":
    main()
