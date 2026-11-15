#!/usr/bin/env python3
"""
Multi-cell-type DecontX rerun (response to decontx_audit_review.md point B,
endorsed by reprocessing_protocol.md Stage 4: "DecontX on the full
multi-cell-type object -- never on a platelet-only matrix"). Step 1: build
the column-index keep-list for Ren's full raw matrix (GSE158055_covid19_
counts.mtx.gz, genes x 1,462,702 cells), restricted to non-platelet lineage
cells (CD4/CD8/B/Mono/NK/Macro/DC/Mega) for the platelet-cohort's Ren-sourced
patients, stratified per patient.

Column-index mapping verified empirically (2026-09-14): row order in
GSE158055_cell_annotation.csv.gz IS the 1-based column index into the mtx
(cross-checked against 20,000 already-trusted pairs in the existing
ren_target_monocyte_cols.tsv with zero mismatches) -- reusing the exact same
convention build_monocyte_canonical.py's build_ren() relies on.
"""
import csv
import gzip
import os

import numpy as np
import pandas as pd

RAW = "/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/monocyte_source"
MAPPING_CSV = ("/bigdata/godziklab/shared/Xinru/302006/302006_updategit/"
               "scripts/31_donor_level_revision/results/donor_meta_with_patient_id.csv")
OUT_DIR = "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/bib_revision_tier2/contamination_audit/decontx_multicelltype"
os.makedirs(OUT_DIR, exist_ok=True)

LINEAGE_BUCKETS = ["CD4", "CD8", "B", "Mono", "NK", "Macro", "DC", "Mega"]
CELLS_PER_BUCKET_PER_PATIENT = 6
RNG_SEED = 20260914

meta = pd.read_csv(MAPPING_CSV).set_index("sample_id")
target = meta[meta["Data_NO"] == "302004data15"]
target_sample_ids = set(target.index.astype(str))
print(f"Target Ren sample_ids (platelet cohort): {len(target_sample_ids)}", flush=True)

rng = np.random.default_rng(RNG_SEED)
buckets = {(sid, lt): [] for sid in target_sample_ids for lt in LINEAGE_BUCKETS}

print("Streaming annotation CSV (1.46M rows)...", flush=True)
with gzip.open(f"{RAW}/ren/GSE158055_cell_annotation.csv.gz", "rt") as f:
    r = csv.reader(f)
    header = next(r)
    assert header == ["cellName", "sampleID", "celltype", "majorType"], header
    for row_idx, row in enumerate(r, start=1):
        cell_name, sample_id, celltype, major_type = row
        if sample_id in target_sample_ids and major_type in LINEAGE_BUCKETS:
            buckets[(sample_id, major_type)].append((row_idx, cell_name))

keep = []
for (sid, lt), items in buckets.items():
    if not items:
        continue
    take = min(len(items), CELLS_PER_BUCKET_PER_PATIENT)
    idx = rng.choice(len(items), size=take, replace=False)
    for i in idx:
        col, cell_name = items[i]
        keep.append((col, cell_name, sid, lt))

keep_df = pd.DataFrame(keep, columns=["col", "cellName", "sampleID", "majorType"])
keep_df = keep_df.sort_values("col").drop_duplicates("col")
print(f"Total non-platelet-lineage cells selected: {len(keep_df)} across "
      f"{keep_df['sampleID'].nunique()} Ren patients", flush=True)
print(keep_df.groupby("majorType").size(), flush=True)

keep_df.to_csv(f"{OUT_DIR}/ren_multilineage_cols.tsv", sep="\t", index=False,
               columns=["col", "cellName"])
keep_df.to_csv(f"{OUT_DIR}/ren_multilineage_meta.csv", index=False)
print("SAVED ren_multilineage_cols.tsv, ren_multilineage_meta.csv", flush=True)
