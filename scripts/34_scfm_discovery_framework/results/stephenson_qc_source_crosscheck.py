"""Per review:20261016_120026 item 2 and review:20261016_130032
R3: a REPRODUCIBLE script (not prose in a handoff) that performs the exact cell-identity
cross-check between the monocyte pilot's raw object and Stephenson's historical Gate 1
QC / DecontX artifacts. Running this script is how the claims in
`STEPHENSON_QC_SOURCE_CROSSCHECK_v1.md` and `stephenson_monocyte_pilot_design_readiness_
table.py`'s `QC_SOURCE_EVIDENCE` can be independently re-derived, rather than re-read as text.

Bounded, read-only: reads 3 small `obs` columns from the 5.7 GB raw h5ad (never `X`), and 2
small CSVs (132 MB, 44 MB) -- never the 12 G/25 G count matrices. Never run on data this large
without checking it is actually bounded first (DEVELOPMENT.md: do not run the real multi-GB input
unboundedly on the login node -- this script's own reads stay under ~200 MB total).

R6 fix (per review:20261016_141242): this NEW script reintroduced the exact
categorical-code defect fixed elsewhere in an earlier round (`stephenson_sample_subtype_cell_
counts.py`'s `_validate_codes`) -- `categories[codes]` indexing without validating for `-1`
(anndata's "missing" sentinel), which Python's negative indexing silently resolves to the LAST
category. `read_monocyte_obs` now validates both code arrays before indexing. Also: duplicate
`cell_id` rows in a CSV used to silently overwrite each other in a plain `dict` comprehension --
`_load_csv_rows_unique` now refuses on any duplicate key. And `flag_would_delete` text now must
be exactly `"true"`/`"false"` (case-insensitive) -- any other value refuses rather than being
read as a silent `False`.
"""
import csv
import hashlib
import json
import os

import h5py
import numpy as np


def _validate_codes(codes, categories, field_name):
    """Refuses (never silently remaps) any code outside [0, len(categories)) -- in
    particular -1 ("missing", anndata's on-disk sentinel), which Python negative
    indexing would otherwise resolve to the LAST category."""
    codes = np.asarray(codes)
    if codes.ndim != 1 or not np.issubdtype(codes.dtype, np.integer):
        raise SystemExit(f"refusing: {field_name} codes array is not a 1-D integer array")
    bad = (codes < 0) | (codes >= len(categories))
    n_bad = int(np.sum(bad))
    if n_bad:
        raise SystemExit(
            f"refusing: {field_name} has {n_bad} missing/out-of-range code(s) "
            f"(valid range 0..{len(categories) - 1})"
        )
    return codes

RAW_H5AD = "/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c_monocyte/stephenson_monocyte_raw.h5ad"
GATE1_FLAGS_CSV = (
    "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/gate1_qc/stephenson/"
    "qc_flags_per_cell.csv"
)
DECONTX_SCORES_CSV = (
    "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/stephenson_decontx_input/"
    "stephenson_decontx_contamination_scores.csv"
)
CELL_ID_PREFIX = "steph_"

HERE = os.path.dirname(__file__)
OUT_JSON = os.path.join(HERE, "stephenson_qc_source_crosscheck_v1.json")


def _sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_monocyte_obs(h5ad_path=RAW_H5AD):
    """Bounded read: obs/_index, obs/Sample_ID, obs/orig_celltype only -- never X."""
    with h5py.File(h5ad_path, "r") as f:
        idx = f["obs"]["_index"][:]
        idx = np.array([x.decode() if isinstance(x, bytes) else x for x in idx])

        sg = f["obs"]["Sample_ID"]
        sample_cats = [c.decode() if isinstance(c, bytes) else c for c in sg["categories"][:]]
        sample_codes = _validate_codes(sg["codes"][:], sample_cats, "Sample_ID")
        sample_of_cell = np.array(sample_cats)[sample_codes]

        tg = f["obs"]["orig_celltype"]
        subtype_cats = [c.decode() if isinstance(c, bytes) else c for c in tg["categories"][:]]
        subtype_codes = _validate_codes(tg["codes"][:], subtype_cats, "orig_celltype")
        subtype_of_cell = np.array(subtype_cats)[subtype_codes]

    if not all(i.startswith(CELL_ID_PREFIX) for i in idx):
        raise SystemExit(f"refusing: not every obs_name starts with {CELL_ID_PREFIX!r}")
    stripped = np.array([i[len(CELL_ID_PREFIX):] for i in idx])
    return idx, stripped, sample_of_cell, subtype_of_cell


def _load_csv_rows(path, fieldnames):
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        return [{k: row[k] for k in fieldnames} for row in reader]


def _load_csv_rows_unique_by_id(path, fieldnames, id_field="cell_id"):
    """Like `_load_csv_rows`, but refuses (never silently overwrites) if `id_field`
    is not unique -- a plain `{r[id_field]: r for r in rows}` dict comprehension
    would otherwise keep only the LAST duplicate row with no warning."""
    rows = _load_csv_rows(path, fieldnames)
    by_id = {}
    duplicates = set()
    for r in rows:
        key = r[id_field]
        if key in by_id:
            duplicates.add(key)
        by_id[key] = r
    if duplicates:
        raise SystemExit(
            f"refusing: {path} has {len(duplicates)} duplicate {id_field} value(s), "
            f"e.g. {sorted(duplicates)[:5]!r} -- a dict keyed by {id_field} would "
            "silently keep only the last row"
        )
    return by_id


def _parse_strict_bool(text, field_name):
    v = text.strip().lower()
    if v == "true":
        return True
    if v == "false":
        return False
    raise SystemExit(f"refusing: {field_name} has an unrecognized boolean value {text!r}")


def crosscheck(h5ad_path=RAW_H5AD, gate1_flags_path=GATE1_FLAGS_CSV, decontx_scores_path=DECONTX_SCORES_CSV):
    """Returns a dict report. Pure computation given the 3 input paths -- no writes."""
    idx, stripped, sample_of_cell, subtype_of_cell = read_monocyte_obs(h5ad_path)

    gate1_by_id = _load_csv_rows_unique_by_id(
        gate1_flags_path, ["cell_id", "channel_id", "n_flags", "flag_would_delete"])
    decontx_by_id = _load_csv_rows_unique_by_id(
        decontx_scores_path, ["cell_id", "channel_id", "decontX_contamination"])

    n_total = len(stripped)
    n_gate1_missing = 0
    n_gate1_channel_mismatch = 0
    n_decontx_missing = 0
    n_decontx_channel_mismatch = 0

    per_subtype = {}
    for cid, sample, subtype in zip(stripped, sample_of_cell, subtype_of_cell):
        bucket = per_subtype.setdefault(subtype, {
            "n_cells": 0, "n_gate1_flag_would_delete_true": 0,
            "decontx_scores": [],
        })
        bucket["n_cells"] += 1

        g = gate1_by_id.get(cid)
        if g is None:
            n_gate1_missing += 1
        else:
            if g["channel_id"] != sample:
                n_gate1_channel_mismatch += 1
            if _parse_strict_bool(g["flag_would_delete"], "flag_would_delete"):
                bucket["n_gate1_flag_would_delete_true"] += 1

        d = decontx_by_id.get(cid)
        if d is None:
            n_decontx_missing += 1
        else:
            if d["channel_id"] != sample:
                n_decontx_channel_mismatch += 1
            bucket["decontx_scores"].append(float(d["decontX_contamination"]))

    per_subtype_summary = {}
    for subtype, bucket in sorted(per_subtype.items()):
        scores = bucket["decontx_scores"]
        per_subtype_summary[subtype] = {
            "n_cells": bucket["n_cells"],
            "n_gate1_flag_would_delete_true": bucket["n_gate1_flag_would_delete_true"],
            "gate1_flag_would_delete_rate": bucket["n_gate1_flag_would_delete_true"] / bucket["n_cells"],
            "decontx_mean_estimated_contamination": (sum(scores) / len(scores)) if scores else None,
            "n_decontx_scores_available": len(scores),
        }

    return {
        "n_cells_total": n_total,
        "n_gate1_flag_missing": n_gate1_missing,
        "n_gate1_channel_mismatch": n_gate1_channel_mismatch,
        "n_decontx_score_missing": n_decontx_missing,
        "n_decontx_channel_mismatch": n_decontx_channel_mismatch,
        "per_subtype": per_subtype_summary,
    }


def _require_output_does_not_already_exist(out_json):
    if os.path.exists(out_json):
        raise SystemExit(f"refusing to overwrite an existing cross-check output: {out_json}")


def main(h5ad_path=RAW_H5AD, gate1_flags_path=GATE1_FLAGS_CSV, decontx_scores_path=DECONTX_SCORES_CSV,
         out_json=OUT_JSON):
    _require_output_does_not_already_exist(out_json)

    report = crosscheck(h5ad_path, gate1_flags_path, decontx_scores_path)
    report["source_files"] = {
        "h5ad_path": h5ad_path,
        "gate1_flags_csv": gate1_flags_path,
        "gate1_flags_csv_sha256": _sha256_file(gate1_flags_path),
        "decontx_scores_csv": decontx_scores_path,
        "decontx_scores_csv_sha256": _sha256_file(decontx_scores_path),
    }
    with open(out_json, "w") as f:
        json.dump(report, f, indent=2, sort_keys=True)
        f.write("\n")

    print(json.dumps(report, indent=2, sort_keys=True))
    return report


if __name__ == "__main__":
    main()
