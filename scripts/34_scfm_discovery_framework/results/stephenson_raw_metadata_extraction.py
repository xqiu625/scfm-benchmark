"""Per `review:20261014_221155` R2: read-only, obs-only extraction
(h5py, never loads X) of the REAL per-sampling-day clinical fields from
Stephenson's own raw object (`covid_portal_210320_with_raw.h5ad`), joined
by the EXACT CellID crosswalk `build_monocyte_canonical.py` itself uses
(`steph_<covid_index>`) to the monocyte-annotated cells already in
`monocyte_canonical_raw.h5ad`/`stephenson_monocyte_raw.h5ad`.

Never loads X on either file. Reports join uniqueness (unmatched/one-to-
many/conflicts), the real donor->visit->sample structure, and a
cross-tab of the ORIGINAL per-sampling-day severity field against the
monocyte-canonical object's own collapsed `covid_severity` (survived/
fatal) -- never silently prefers one over the other.
"""
import json
import h5py
import numpy as np
import pandas as pd

RAW_H5AD = ("/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/monocyte_source/"
           "stephenson/covid_portal_210320_with_raw.h5ad")
MONO_H5AD = "/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c_monocyte/stephenson_monocyte_raw.h5ad"
MAP = ("/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/"
      "31_donor_level_revision/results/donor_meta_with_patient_id.csv")
MONO_LABELS_STEPHENSON = ["C1_CD16_mono", "CD14_mono", "CD16_mono", "CD83_CD14_mono", "Mono_prolif"]


def _decode_codes(cats: np.ndarray, codes: np.ndarray) -> np.ndarray:
    """Shared, safe code->category decode for BOTH categorical encodings
    below. Per `review:20261014_231259` R1 (metadata-script item):
    a bare `cats[codes]` silently wraps a missing code (`-1`, pandas/
    anndata's own missing-value sentinel) around to the LAST real
    category via numpy negative indexing -- confirmed by the review's own
    synthetic reproduction (`categories=['Mild','Severe']`,
    `codes=[0,-1,1]` decoded to `['Mild','Severe','Severe']`, silently
    turning a MISSING value into a claimed 'Severe'). Any code outside
    `[-1, len(cats))` is rejected outright (never silently clamped/wrapped)
    -- a value this far out of range indicates a genuinely malformed file,
    not an ordinary missing-value case.
    """
    codes = np.asarray(codes)
    out_of_range = (codes < -1) | (codes >= len(cats))
    if out_of_range.any():
        bad = sorted(set(codes[out_of_range].tolist()))
        raise ValueError(
            f"categorical codes out of range [-1, {len(cats)}): {bad!r} -- refusing to "
            "decode rather than silently wrapping/clamping")
    out = np.empty(codes.shape, dtype=object)
    valid = codes >= 0
    out[valid] = cats[codes[valid]]
    out[~valid] = "nan"
    return out.astype(str)


def rc(g):
    """Robust obs-column decoder: handles the NEWER anndata encoding (a
    group with 'categories'/'codes' sub-datasets), the LEGACY pandas/h5ad
    categorical encoding (a plain int dataset with an `attrs['categories']`
    HDF5 object reference, codes==-1 meaning missing), and a plain
    non-categorical dataset. `covid_portal_210320_with_raw.h5ad` uses the
    LEGACY encoding for several columns -- confirmed directly this round
    (`sample_id`'s own `attrs` is `{'categories': <HDF5 object reference>}`,
    not a sibling group). Both categorical branches share `_decode_codes`'s
    missing/out-of-range handling -- the first version of this function
    only added that handling to the LEGACY branch, leaving the newer
    group-style branch vulnerable to the same wraparound bug.
    """
    if isinstance(g, h5py.Group) and "categories" in g:
        return _decode_codes(g["categories"][:].astype(str), g["codes"][:])
    if "categories" in g.attrs:
        cats = g.file[g.attrs["categories"]][:].astype(str)
        return _decode_codes(cats, g[:])
    return g[:].astype(str)


def _missing_code_counts(h, cols: list) -> dict:
    """Per `review:20261014_231259` R1: before trusting any
    earlier run's real statistics, check whether any of the ACTUAL columns
    this script reads hit the missing/out-of-range-code case at all --
    distinguishes a latent bug in the decoder from an actual impact on
    this round's real findings.
    """
    report = {}
    for c in cols:
        g = h["obs"][c]
        if isinstance(g, h5py.Group) and "categories" in g:
            codes = g["codes"][:]
            n_cats = len(g["categories"])
        elif "categories" in g.attrs:
            codes = g[:]
            n_cats = len(g.file[g.attrs["categories"]])
        else:
            continue
        report[c] = {
            "n_missing_(-1)": int((codes == -1).sum()),
            "n_out_of_range": int(((codes < -1) | (codes >= n_cats)).sum()),
        }
    return report


def main():
    meta = pd.read_csv(MAP).set_index("sample_id")
    target = meta[meta["Data_NO"] == "302004data17"]
    target_sample_ids = set(target.index)
    print(f"target samples (Data_NO=302004data17): {len(target_sample_ids)}")

    print("\n--- Step 1: obs-only read of the REAL raw Stephenson h5ad (never X) ---")
    print("(filters to the target mask BEFORE building any full-atlas DataFrame, to keep "
         "peak memory low on the login node -- a full 647k-row x 13-col DataFrame previously "
         "triggered this session's login-node memory cap; see claude-code_output for the trace)")
    other_cols = ["covid_index", "patient_id", "Resample", "Collection_Day",
                 "Status_on_day_collection", "Status_on_day_collection_summary",
                 "Days_from_onset", "Site", "Worst_Clinical_Status", "Outcome", "Status"]
    all_decoded_cols = ["sample_id", "full_clustering"] + other_cols
    with h5py.File(RAW_H5AD, "r") as h:
        n_total = h["obs"]["sample_id"].shape[0]
        print("\nmissing/out-of-range categorical code counts, ACTUAL columns read this run "
             "(per codex_output/20261014_231259.md R1 -- distinguishes a latent decoder bug "
             "from a real impact on this round's own findings):")
        missing_report = _missing_code_counts(h, all_decoded_cols)
        for c, counts in missing_report.items():
            print(f"  {c}: {counts}")
        any_missing_or_oor = any(v["n_missing_(-1)"] or v["n_out_of_range"] for v in missing_report.values())
        print(f"any column hit a missing/out-of-range code: {any_missing_or_oor}")

        full_clustering = rc(h["obs"]["full_clustering"])
        sample_id = rc(h["obs"]["sample_id"])
        is_mono = np.isin(full_clustering, MONO_LABELS_STEPHENSON)
        in_cohort = np.isin(sample_id, np.array(list(target_sample_ids)))
        mask = is_mono & in_cohort
        print(f"raw obs rows (total atlas cells): {n_total}")
        print(f"monocyte AND target-cohort cells (fresh, from raw file): {int(mask.sum())}")
        data = {"sample_id": sample_id[mask], "full_clustering": full_clustering[mask]}
        for c in other_cols:
            data[c] = rc(h["obs"][c])[mask]
    sel = pd.DataFrame(data)
    sel["steph_cellid"] = "steph_" + sel["covid_index"]

    print("\n--- Step 2: join-uniqueness check against the EXISTING monocyte canonical object ---")
    with h5py.File(MONO_H5AD, "r") as h:
        mono_index = rc(h["obs"]["_index"]) if "_index" in h["obs"] else rc(h["obs"]["covid_index"])
        mono_sample_id = rc(h["obs"]["Sample_ID"])
        mono_patient_id = rc(h["obs"]["patient_id"])
    existing = pd.DataFrame({"steph_cellid": mono_index, "existing_sample_id": mono_sample_id,
                            "existing_patient_id": mono_patient_id})
    print(f"existing monocyte-canonical Stephenson cells: {len(existing)}")

    dup_fresh = sel["steph_cellid"].duplicated().sum()
    dup_existing = existing["steph_cellid"].duplicated().sum()
    print(f"duplicate steph_cellid within fresh extraction: {dup_fresh}")
    print(f"duplicate steph_cellid within existing canonical object: {dup_existing}")

    merged = existing.merge(sel, on="steph_cellid", how="left", indicator=True)
    unmatched = merged[merged["_merge"] == "left_only"]
    print(f"existing cells with NO match in the fresh raw-file extraction: {len(unmatched)} "
         f"of {len(existing)}")
    if len(unmatched):
        print("  sample of unmatched steph_cellid:", unmatched["steph_cellid"].head(5).tolist())
    sample_id_conflicts = merged[(merged["_merge"] == "both") &
                                 (merged["existing_sample_id"] != merged["sample_id"])]
    print(f"existing vs. fresh sample_id CONFLICTS (same cell, different sample_id): "
         f"{len(sample_id_conflicts)}")
    patient_id_conflicts = merged[(merged["_merge"] == "both") &
                                  (merged["existing_patient_id"] != merged["patient_id"])]
    print(f"existing vs. fresh patient_id CONFLICTS: {len(patient_id_conflicts)}")

    # Per `review:20261014_231259` R1: these must become hard failures before this
    # script's output is trusted as a formal metadata snapshot -- printing and continuing to
    # `drop_duplicates()`/write the output anyway would silently paper over exactly the
    # problems this join-diagnostic step exists to catch.
    join_problems = []
    if dup_fresh:
        join_problems.append(f"{dup_fresh} duplicate steph_cellid value(s) within the fresh extraction")
    if dup_existing:
        join_problems.append(f"{dup_existing} duplicate steph_cellid value(s) within the existing canonical object")
    if len(unmatched):
        join_problems.append(f"{len(unmatched)} existing cell(s) with no match in the fresh extraction")
    if len(sample_id_conflicts):
        join_problems.append(f"{len(sample_id_conflicts)} existing-vs-fresh sample_id conflict(s)")
    if len(patient_id_conflicts):
        join_problems.append(f"{len(patient_id_conflicts)} existing-vs-fresh patient_id conflict(s)")
    if join_problems:
        raise SystemExit(f"refusing to produce a metadata snapshot: {join_problems}")

    print("\n--- Step 3: real donor -> visit -> sample structure (sample-level, deduplicated) ---")
    sample_level = sel.drop_duplicates(subset="sample_id")[
        ["sample_id", "patient_id", "Resample", "Collection_Day", "Days_from_onset", "Site",
         "Status_on_day_collection", "Status_on_day_collection_summary", "Worst_Clinical_Status"]
    ].sort_values(["patient_id", "Collection_Day"])
    print(f"distinct samples (monocyte-contributing, this cohort): {len(sample_level)}")
    print(f"distinct patients: {sample_level['patient_id'].nunique()}")
    visits_per_patient = sample_level.groupby("patient_id").size()
    print(f"patients with >1 sample (repeated visit): {int((visits_per_patient > 1).sum())} "
         f"of {len(visits_per_patient)}")
    multi = sample_level[sample_level["patient_id"].isin(
        visits_per_patient[visits_per_patient > 1].index)]
    print("\nmulti-visit patients' real sample-level rows:")
    print(multi.to_string(index=False))

    # Per `review:20261014_231259` R2: a patient_id with exactly ONE record is NOT
    # automatically a clean single-visit case -- if that one record is itself tagged
    # 'Resample' rather than 'Initial', its real Initial visit may simply be missing from
    # THIS monocyte-restricted subset (e.g. it had no monocyte cells meeting the mask), or
    # the same biological donor may be recorded under a DIFFERENT patient_id string for a
    # different visit (e.g. a 'v2'-suffixed id). This table alone cannot adjudicate that --
    # it only reports the fact plainly, as a real evidence gap, never silently merging by
    # suffix or assuming independence.
    single_record = sample_level[sample_level["patient_id"].isin(
        visits_per_patient[visits_per_patient == 1].index)]
    resample_only = single_record[single_record["Resample"] != "Initial"]
    print(f"\nsingle-record patient_ids whose ONLY record is tagged Resample (not Initial) -- "
         f"identity status unresolved, see evidence gap note: {len(resample_only)}")
    if len(resample_only):
        print(resample_only.to_string(index=False))

    # Internal consistency: per `review:20261014_231259` R1, this must cover EVERY
    # kept clinical field (not just 3 of the 9) -- a narrower check that reports "0
    # inconsistent" is not evidence the OTHER fields are consistent too.
    consistency_cols = ["patient_id", "Status_on_day_collection_summary", "Status_on_day_collection",
                        "Collection_Day", "Resample", "Days_from_onset", "Site",
                        "Worst_Clinical_Status", "Outcome", "Status"]
    per_sample_nunique = sel.groupby("sample_id")[consistency_cols].nunique()
    inconsistent_samples = per_sample_nunique[(per_sample_nunique > 1).any(axis=1)]
    print(f"\nsamples with internally INCONSISTENT clinical fields across their own cells "
         f"(checked: {consistency_cols}): {len(inconsistent_samples)}")
    if len(inconsistent_samples):
        raise SystemExit(
            f"refusing to produce a metadata snapshot: {len(inconsistent_samples)} sample(s) "
            f"have internally inconsistent clinical fields: {inconsistent_samples.index.tolist()}")

    print("\n--- Step 4: cross-tab, ORIGINAL per-sampling-day severity vs. the monocyte-canonical "
         "object's own collapsed covid_severity (survived/fatal) -- patient level ---")
    # covid_severity in the EXISTING merged canonical object is np.where(y==1, 'fatal', 'survived').
    patient_level = sel.drop_duplicates(subset="patient_id")[
        ["patient_id", "Status_on_day_collection_summary"]].merge(
        meta.reset_index()[["patient_id", "y"]].drop_duplicates(subset="patient_id"),
        on="patient_id", how="left")
    patient_level["covid_severity_collapsed"] = np.where(
        patient_level["y"] == 1, "fatal", "survived")
    tab = patient_level.groupby(
        ["Status_on_day_collection_summary", "covid_severity_collapsed"]).size().unstack(fill_value=0)
    print(tab.to_string())

    print("\n--- Step 5: known identity caveat (stage33 decision sheet v3) -- CV0201/CV0326 ---")
    for cv in ("CV0201", "CV0326"):
        present = cv in set(sel["patient_id"])
        print(f"  {cv} present among this cohort's monocyte-contributing patients: {present}")

    print("\n--- Step 6: Days_from_onset -- numeric vs. non-numeric values (never blanket-cast) ---")
    numeric_mask = pd.to_numeric(sel["Days_from_onset"], errors="coerce").notna()
    print(f"Days_from_onset numeric: {int(numeric_mask.sum())}, "
         f"non-numeric (e.g. 'Healthy'): {int((~numeric_mask).sum())}")
    print("non-numeric distinct values:", sorted(sel.loc[~numeric_mask, "Days_from_onset"].unique()))

    out = {
        "missing_code_report": missing_report,
        "n_fresh_monocyte_target_cells": int(mask.sum()),
        "n_existing_canonical_cells": len(existing),
        "n_unmatched_existing_cells": len(unmatched),
        "n_sample_id_conflicts": len(sample_id_conflicts),
        "n_patient_id_conflicts": len(patient_id_conflicts),
        "n_resample_only_single_record_patients": len(resample_only),
        "resample_only_patient_ids": resample_only["patient_id"].tolist(),
        "n_distinct_samples": len(sample_level),
        "n_distinct_patients": int(sample_level["patient_id"].nunique()),
        "n_multi_visit_patients": int((visits_per_patient > 1).sum()),
        "n_internally_inconsistent_samples": len(inconsistent_samples),
    }
    with open("stephenson_raw_metadata_extraction_summary.json", "w") as fh:
        json.dump(out, fh, indent=2)
    sample_level.to_csv("stephenson_sample_level_real_fields.csv", index=False)
    sel.to_csv("stephenson_cell_level_real_fields_monocyte_target.csv", index=False)
    print("\nwrote: stephenson_raw_metadata_extraction_summary.json, "
         "stephenson_sample_level_real_fields.csv, "
         "stephenson_cell_level_real_fields_monocyte_target.csv")


if __name__ == "__main__":
    main()
