"""Per review:20261015_230801 increment 2: a BOUNDED, metadata-only
read of the frozen per-cohort raw object's cell-level `obs` annotations -- never
`X`/`layers`/`obsm`, never a full-array load. Produces a frozen sample x subtype
cell-count table so the longitudinal feasibility note can answer "does each
multi-visit donor's repeated visit actually have comparable monocyte subtype
coverage" with real numbers instead of leaving it an assumed yes.

Source object (read-only, not modified, not copied):
`/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c_monocyte/
stephenson_monocyte_raw.h5ad` -- the SAME file `monocyte_pilot_feasibility_v3.tsv`'s
Stephenson row already names as having `orig_celltype` "in per-cohort raw file
(not yet in the merged monocyte_canonical_raw.h5ad)". That row's own
`qc_input_admissibility_status` caveat ("provisional... a feasibility read, not
an admitted input") applies here unchanged: this table is frozen cell-count
EVIDENCE for the feasibility question, never a QC-admitted analysis input, and
never silently promoted into one.

Only two `obs` columns are read (`Sample_ID`, `orig_celltype`), both stored as
anndata categorical groups (`categories` + small `int8` `codes` arrays, ~128KB
each for this object's 128,503 cells) -- h5py reads exactly the bytes sliced,
so this never touches the 6+ GB `X` matrix. Confirmed by inspection
(`h5py.File(path, "r")["obs"].keys()` run against this exact file, 2026-10-16):
top-level keys are `X, layers, obs, obsm, obsp, uns, var, varm, varp`; this
script opens none of them except `obs/Sample_ID` and `obs/orig_celltype`.

R1 fix (review:20261016_102834) -- the counting loop below used
to index `sample_cats[s_code]`/`subtype_cats[t_code]` directly. Anndata's own
on-disk categorical format reserves code `-1` for "missing" (see
https://anndata.readthedocs.io/en/stable/fileformat-prose.html#categorical-arrays);
Python's negative indexing silently maps `-1` to the LAST category rather than
raising, so a missing code would have been counted as a real observation of
whichever category happened to sort last. Reproduced on an isolated two-cell
fixture (`sample_codes=[0, -1]`, `subtype_codes=[-1, 0]`): the unfixed loop
returned `s0/B=1, s1/A=1` with no error for either cell. `_validate_codes` now
rejects (refuses, does not silently drop) any code outside `[0, len(categories))`
before the counting loop runs. Confirmed on the real object (2026-10-16): both
`Sample_ID` codes (range 0-101) and `orig_celltype` codes (range 0-4) are fully
in-range for all 128,503 cells -- this fix changes no real count already
committed, only how a future missing/out-of-range code would be handled.
"""
import csv
import hashlib
import json
import os

import h5py
import numpy as np

H5AD_PATH = "/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c_monocyte/stephenson_monocyte_raw.h5ad"

HERE = os.path.dirname(__file__)
OUT_TSV = os.path.join(HERE, "stephenson_sample_subtype_cell_counts_v1.tsv")
OUT_SUMMARY = os.path.join(HERE, "stephenson_sample_subtype_cell_counts_v1_summary.json")


def _sha256_file(path, limit_bytes=None):
    """Hashes the whole file by default. `limit_bytes` exists only so a
    mutation/smoke test can hash a small fixture quickly -- never used
    against the real 6 GB object (hashing it fully is a one-time, bounded,
    sequential read, not a parse of its structured contents)."""
    h = hashlib.sha256()
    read = 0
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
            read += len(chunk)
            if limit_bytes is not None and read >= limit_bytes:
                break
    return h.hexdigest()


def _validate_codes(codes, categories, field_name):
    """Refuses (never silently drops or remaps) any code that is not a
    well-formed index into `categories`. Anndata's on-disk format uses -1 for
    a missing value; Python's `categories[-1]` would resolve that to the LAST
    category instead of raising, which is exactly the R1 defect this guards
    against (review:20261016_102834)."""
    codes = np.asarray(codes)
    if codes.ndim != 1:
        raise SystemExit(f"refusing to count: obs/{field_name} codes array is not 1-D (shape={codes.shape})")
    if not np.issubdtype(codes.dtype, np.integer):
        raise SystemExit(
            f"refusing to count: obs/{field_name} codes array is not an integer dtype (dtype={codes.dtype})"
        )
    n_missing = int(np.sum(codes == -1))
    out_of_range = (codes < -1) | (codes >= len(categories))
    n_out_of_range = int(np.sum(out_of_range))
    if n_missing or n_out_of_range:
        raise SystemExit(
            f"refusing to count: obs/{field_name} has {n_missing} missing (-1) code(s) and "
            f"{n_out_of_range} out-of-range code(s) (valid range 0..{len(categories) - 1}); "
            "missing/invalid codes must not be mapped to any category"
        )
    return codes


def read_sample_subtype_counts(h5ad_path=H5AD_PATH):
    """Returns (counts: dict[sample_id][subtype] -> int, n_cells_total,
    sample_ids: list, subtypes: list). Reads ONLY obs/Sample_ID and
    obs/orig_celltype -- two categorical code arrays, never X/layers/obsm."""
    with h5py.File(h5ad_path, "r") as f:
        sample_group = f["obs"]["Sample_ID"]
        subtype_group = f["obs"]["orig_celltype"]
        sample_cats = [c.decode() if isinstance(c, bytes) else c for c in sample_group["categories"][:]]
        subtype_cats = [c.decode() if isinstance(c, bytes) else c for c in subtype_group["categories"][:]]
        sample_codes = sample_group["codes"][:]
        subtype_codes = subtype_group["codes"][:]

    if len(sample_codes) != len(subtype_codes):
        raise SystemExit(
            f"refusing to count: obs/Sample_ID has {len(sample_codes)} cells but "
            f"obs/orig_celltype has {len(subtype_codes)} -- not the same length"
        )

    sample_codes = _validate_codes(sample_codes, sample_cats, "Sample_ID")
    subtype_codes = _validate_codes(subtype_codes, subtype_cats, "orig_celltype")

    counts = {sid: {st: 0 for st in subtype_cats} for sid in sample_cats}
    for s_code, t_code in zip(sample_codes, subtype_codes):
        counts[sample_cats[s_code]][subtype_cats[t_code]] += 1

    return counts, len(sample_codes), sample_cats, subtype_cats


def _require_output_does_not_already_exist(out_tsv, out_summary):
    existing = [p for p in (out_tsv, out_summary) if os.path.exists(p)]
    if existing:
        raise SystemExit(f"refusing to overwrite an existing cell-count output: {existing}")


def main(h5ad_path=H5AD_PATH, out_tsv=OUT_TSV, out_summary=OUT_SUMMARY):
    _require_output_does_not_already_exist(out_tsv, out_summary)

    counts, n_cells_total, sample_cats, subtype_cats = read_sample_subtype_counts(h5ad_path)

    fieldnames = ["sample_id", "n_cells_this_sample_total"] + subtype_cats
    with open(out_tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for sid in sorted(sample_cats):
            row = {"sample_id": sid, "n_cells_this_sample_total": sum(counts[sid].values())}
            row.update(counts[sid])
            w.writerow(row)

    summary = {
        "source_h5ad": h5ad_path,
        "source_h5ad_sha256": _sha256_file(h5ad_path),
        "source_h5ad_bytes_read_for_hash": "full file, one sequential pass",
        "fields_read_from_obs": ["Sample_ID", "orig_celltype"],
        "fields_NOT_read": ["X", "layers", "obsm", "obsp", "var", "varm", "varp", "uns",
                             "Data_NO", "disease", "patient_id", "source_study", "y", "_index"],
        "n_cells_total_this_object": n_cells_total,
        "n_samples": len(sample_cats),
        "subtypes": subtype_cats,
        "provenance_caveat": (
            "This is the per-cohort RAW object's own obs['orig_celltype'] (curated_atlas_label, "
            "the authors' published full_clustering annotation) -- not a scFM-discovered state, "
            "and not QC-filtered for this project (qc_status: provisional, same caveat as "
            "monocyte_pilot_feasibility_v3.tsv's Stephenson row). A sample's count for a given "
            "subtype in THIS object is a real, confirmed number for THIS specific file/version -- "
            "it is not the same claim as 'this is the final QC-admitted cell count for this "
            "analysis', and a sample absent from this object entirely (none observed here) would "
            "be a mapping gap, reported as such, never silently folded into a subtype's zero count."
        ),
    }
    with open(out_summary, "w") as f:
        json.dump(summary, f, indent=2, sort_keys=True)
        f.write("\n")

    print(json.dumps(summary, indent=2, sort_keys=True))
    return summary


if __name__ == "__main__":
    main()
