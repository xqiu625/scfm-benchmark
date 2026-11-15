#!/usr/bin/env python3
"""Prepares a BOUNDED, patient-balanced h5ad for `run_batch_diagnostics.py`
from the existing real (pre-Gate-2, PROVISIONAL) assembled benchmark input
-- never reads the full 22GB/50GB counts.mtx into memory; streams through it
column-chunk by column-chunk via the already-built, already-used
`read_mtx_columns()` (scripts/32_rebuild_upfront_qc/code/
assemble_platelet_monocyte_benchmark_input.py), selecting only the wanted
cell columns while reading.

This is explicitly a PROVISIONAL input, not a final-benchmark artifact: the
source matrices (`05_RESULTS/rebuild_upfront_qc/platelet_monocyte_benchmark_
input_v2/{raw,decontx}/`) predate Gate 2 (doublet scoring is not yet
resolved; see docs/current-status.md) -- every output h5ad this script
writes gets `uns["population_stage"] = "pre-Gate-2 PROVISIONAL, not final"`,
and `run_batch_diagnostics.py` is passed this stage label explicitly, so no
downstream output can be mistaken for final-benchmark evidence.

v2 fixes defects a 2026-10-09 review found in v1
(review:20261009_133450), confirmed independently before fixing:

- Seeding used Python's built-in `hash()` (per-process-randomized for
  strings) -- replaced with `stable_seed.derive_seed()`.
- The patient-balanced cap was applied ACROSS both lineages combined
  (grouping by bare `donor_id` over the already lineage-filtered rows), so
  a donor with cells in both lineages could lose most/all of a rarer
  lineage's cells to the more abundant one. Fixed: sampling is now per
  (cohort, donor_id, lineage) -- the cap applies within each lineage
  separately, and the groupby key includes `cohort` so a donor_id string
  reused (even coincidentally) across two different cohorts cannot merge
  two different biological patients' cells into one quota.
- `raw/cells.csv` and `decontx/cells.csv` were assumed, not verified, to
  have identical row ordering before applying the SAME `idx_in_full_file`
  positions to both matrices. Fixed: both are read and compared directly.
- `event` was mapped with `.map({1.0:"fatal", 0.0:"survived"}).fillna("unknown")`
  with no check that only {0.0, 1.0, NaN} are actually present -- an
  unexpected value would have silently become "unknown" instead of being
  caught. Fixed: validated with the same exact-equality discipline as
  `build_patient_batch_confounding_table.classify_event`.
- No no-clobber check on `--out`. Fixed: refuses if the output file (or a
  sibling `<out>.selection_manifest.csv`) already exists.
- No persisted, reviewable selection manifest -- a reader had to trust the
  printed per-cohort/lineage counts rather than see the exact cell_ids
  chosen. Fixed: writes `<out>.selection_manifest.csv` (cell_id, cohort,
  donor_id, lineage, orig_row_index) alongside the h5ad.
- `adata.X`/layers had no declared expression-source contract (silently
  "whatever was assigned to X", here decontx). Fixed: `uns["expression_source"]`
  states it explicitly, and a `"counts"` layer alias is added so
  `run_batch_diagnostics.compute_pca`'s layer lookup does not depend on an
  implicit X-already-equals-the-right-thing fallback.

v3 fixes defects a 2026-10-09 review found (review:20261009_230728),
confirmed independently before fixing:

- `meta["patient_id"] = meta["donor_id"]` dropped the cohort qualifier this
  same script's own sampler uses -- a donor_id string reused (even
  coincidentally) across two cohorts would be treated as ONE patient
  downstream. Fixed: `patient_id` is now `f"{cohort}::{donor_id}"`.
- A `--dry-run`'s selection manifest was computed then immediately deleted --
  no durable, reviewable artifact existed to bind a real extraction to. Fixed:
  the manifest is always kept, with a `.status.json` sidecar distinguishing
  "preflight_only" from "extracted"; a new `--selection-manifest <path>`
  lets a real run consume and re-validate a previously-reviewed one instead
  of only recomputing deterministically.
- No atomic ownership of `--out` before the streamed matrix reads begin --
  only a check-then-act existence check, the same race class found in
  `run_batch_diagnostics.py`. Fixed: `claim_ownership()` on a non-dry-run.
- `uns["population_stage"]`'s free text could drift from a caller's
  documented wording with no hard check anywhere. Fixed: also writes the
  fixed, machine-checkable `uns["population_stage_id"]` (see
  `scripts/00_shared/code/population_stage.py`).

v4 fixes a defect a 2026-10-10 review found (review:20261010_120725),
confirmed independently before fixing:

- The non-dry-run path wrote the status sidecar's `status="extracted"`
  BEFORE either `counts.mtx` had been streamed or the h5ad written/reloaded
  -- a process killed or raising during that work left an "extracted"
  sidecar pointing at an h5ad that was never written, or only partially
  written. Fixed: `status="selected"` is published first (true at that
  point -- nothing downstream has happened yet); `status="extracted"` is
  published only after the streamed reads, the h5ad write, the write/reload
  shape+order check, AND the patient-id/lineage-count reconciliation checks
  below all succeed. Any exception raised in between is caught, published
  as `status="failed"` with the error recorded, and re-raised (the process
  still exits non-zero). Status sidecar writes are now atomic
  (same-directory temp file, fsynced, `os.replace`'d -- mirrors
  `scripts/32_rebuild_upfront_qc/code/gate2_checkpoint.py`'s existing
  convention), so a reader never observes a partially-written sidecar.
- The provisional job spec had pointed Job 1's acceptance criteria at
  `manifest.json` fields (`code_revision_verified`, etc.) that only
  `run_batch_diagnostics.py` (Job 2) writes -- this script never writes a
  `manifest.json` at all. Fixed by giving the `"extracted"` status payload
  its own completion record instead: `output_sha256`, `n_cells`, `n_genes`,
  `stage_id`, `expression_source`, and explicit boolean fields for the
  patient-id-cohort-qualified and per-lineage-count-reconciliation checks --
  each of which is RAISED on, not merely recorded, if it fails (see the
  `raise ValueError` calls immediately before `final_status` is built).

Sampling rule: deterministic (`stable_seed.derive_seed`, not Python's
`hash()`), capped at `--n-per-donor` cells per (cohort, donor_id) WITHIN
EACH requested lineage separately. This caps estimated matrix size to
roughly (n_distinct (cohort,donor_id,lineage) combinations * n_per_donor)
rows at most, computed up front from `cells.csv` alone and logged/printed
before the (streamed) matrix read begins, so a reviewer knows the real
output size before approving the job. `--dry-run` stops right there.

Usage:
    python prepare_provisional_diagnostic_input.py \\
        --input-root /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/platelet_monocyte_benchmark_input_v2 \\
        --out <new path.h5ad> --n-per-donor 20 --seed 0
"""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import sys
import uuid
from pathlib import Path

import anndata
import numpy as np
import pandas as pd

sys.path.insert(
    0, str(Path(__file__).resolve().parents[2] / "32_rebuild_upfront_qc/code")
)
from assemble_platelet_monocyte_benchmark_input import read_mtx_columns  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "00_shared/code"))
from stable_seed import derive_seed  # noqa: E402
from run_ownership import claim_ownership  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("prepare_provisional_diagnostic_input")

STAGE_ID = "pre_gate2_decontx_provisional"
SELECTION_MANIFEST_COLS = ["cell_id", "cohort", "donor_id", "harmonized_lineage", "orig_row_index"]


def _file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(1 << 20)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def _atomic_write_json(path: Path, payload: dict) -> None:
    """Same-directory temp file, flushed and fsynced, then `os.replace`d, so a
    reader never observes a partially-written or truncated status sidecar --
    mirrors this project's existing atomic-write convention
    (scripts/32_rebuild_upfront_qc/code/gate2_checkpoint.py). The unique
    pid+uuid4 suffix means a failed write can only ever delete its own temp."""
    tmp = path.with_name(f".{path.name}.{os.getpid()}.{uuid.uuid4().hex}.tmp")
    with open(tmp, "w") as f:
        json.dump(payload, f, indent=2)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def validate_event_column(event: pd.Series) -> pd.Series:
    """Exact-equality 0/1/NaN validation -- same discipline as
    build_patient_batch_confounding_table.classify_event. Raises on any
    other value rather than silently mapping it to "unknown"."""
    numeric = pd.to_numeric(event, errors="raise")
    is_one, is_zero, is_na = numeric == 1.0, numeric == 0.0, numeric.isna()
    unexpected = ~(is_one | is_zero | is_na)
    if unexpected.any():
        raise ValueError(
            f"event column has {int(unexpected.sum())} value(s) that are neither exactly 0, "
            f"exactly 1, nor NaN: {sorted(numeric[unexpected].unique())}"
        )
    return pd.Series(np.where(is_one, "fatal", np.where(is_zero, "survived", "unknown")),
                      index=event.index)


def patient_balanced_index_per_lineage(cells: pd.DataFrame, n_per_donor: int,
                                        seed: int) -> np.ndarray:
    """Deterministic (stable_seed.derive_seed, not hash()). Groups by
    (cohort, donor_id, harmonized_lineage) -- NOT bare donor_id -- so the
    cap applies within each lineage separately (a donor with cells in two
    lineages gets up to n_per_donor of EACH, not n_per_donor total split
    unpredictably between them) and a donor_id string accidentally reused
    across two cohorts cannot merge their quotas. Operates on row POSITIONS
    (0-indexed into `cells`), since that is what read_mtx_columns needs."""
    chosen = []
    groups = cells.groupby(["cohort", "donor_id", "harmonized_lineage"]).groups
    for key, idx in groups.items():
        idx = np.asarray(idx)
        n_take = min(len(idx), n_per_donor)
        rng = np.random.RandomState(derive_seed(seed, "cohort_donor_lineage", *key))
        chosen.extend(rng.choice(idx, size=n_take, replace=False))
    return np.array(sorted(chosen))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-root", required=True, type=Path,
                    help="Directory containing raw/{cells.csv,counts.mtx,genes.csv} "
                         "and decontx/{cells.csv,counts.mtx} (genes.csv shared).")
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--n-per-donor", type=int, default=20)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--lineages", nargs="*", default=["Platelet", "Monocyte_Myeloid"])
    p.add_argument("--dry-run", action="store_true",
                    help="Compute and report the patient-balanced selection and estimated "
                         "output size from cells.csv only; do not stream the counts.mtx "
                         "files at all. Safe to run on a login node (cells.csv is small); "
                         "the real run (streaming 22GB/50GB counts.mtx) is not. The selection "
                         "manifest IS written (status 'preflight_only', not deleted -- a "
                         "2026-10-09 review asked for a durable, reviewable artifact a later "
                         "real run can be bound to via --selection-manifest).")
    p.add_argument("--selection-manifest", type=Path, default=None,
                    help="Consume a previously-written <out>.selection_manifest.csv (e.g. "
                         "from a reviewed --dry-run) instead of recomputing the sample. "
                         "Re-validates every cell_id/orig_row_index pair against the CURRENT "
                         "raw/cells.csv before use and refuses on any drift -- this binds a "
                         "real extraction to a specific reviewed list, not merely to "
                         "'the same seed would reproduce this'.")
    args = p.parse_args()

    manifest_path = args.out.with_suffix(args.out.suffix + ".selection_manifest.csv")
    status_path = args.out.with_suffix(args.out.suffix + ".selection_manifest.status.json")
    if args.out.exists() or manifest_path.exists() or status_path.exists():
        raise FileExistsError(
            f"{args.out}, its selection manifest, or its status sidecar already exists -- "
            f"refusing to overwrite (this project's standing rule: never overwrite a result; "
            f"use a fresh path)."
        )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    if not args.dry_run:
        ownership = claim_ownership(
            args.out.with_suffix(args.out.suffix + ".run_owner.json"),
            command=" ".join(sys.argv),
        )
        log.info("Claimed run ownership: %s", ownership)

    raw_cells_path = args.input_root / "raw" / "cells.csv"
    decontx_cells_path = args.input_root / "decontx" / "cells.csv"
    genes_path = args.input_root / "genes.csv"  # shared between raw/ and decontx/, lives at input_root top level
    raw_mtx_path = args.input_root / "raw" / "counts.mtx"
    decontx_mtx_path = args.input_root / "decontx" / "counts.mtx"
    for req in [raw_cells_path, decontx_cells_path, genes_path, raw_mtx_path, decontx_mtx_path]:
        if not req.exists():
            raise FileNotFoundError(f"Required input not found: {req}")

    # Identity check done FIRST, with single-column reads of BOTH files, so the
    # (larger, 11-column) full raw/cells.csv load below does not have to coexist
    # in memory with a second full decontx/cells.csv load -- a real 2026-10-09
    # login-node SIGKILL (exit 137) was observed with the heavier read pattern,
    # consistent with this project's own documented login-node memory ceiling
    # (DEVELOPMENT.md: "Do not run the real 5GB input on the login node"). This is a
    # memory-pattern fix, not a logic change -- the check itself is unchanged
    # and still covered by tests/test_prepare_provisional_diagnostic_input.py.
    log.info("Verifying decontx/cells.csv has the IDENTICAL cell_id ordering as raw/cells.csv "
              "(cell_id-only reads of both files first, before loading the full raw/cells.csv)")
    raw_ids = pd.read_csv(raw_cells_path, usecols=["cell_id"])["cell_id"].to_numpy()
    decontx_ids = pd.read_csv(decontx_cells_path, usecols=["cell_id"])["cell_id"].to_numpy()
    if len(raw_ids) != len(decontx_ids) or not (raw_ids == decontx_ids).all():
        raise ValueError(
            f"raw/cells.csv and decontx/cells.csv do not have identical cell_id ordering "
            f"({len(raw_ids)} vs {len(decontx_ids)} rows) -- refusing to apply the same "
            f"column-position selection to both counts.mtx files without this guarantee."
        )
    log.info("Confirmed: %d rows, identical cell_id ordering in raw/ and decontx/ cells.csv", len(raw_ids))
    del raw_ids, decontx_ids

    log.info("Reading %s (full columns, cheap -- the only full file read before --dry-run returns)",
              raw_cells_path)
    cells = pd.read_csv(raw_cells_path, low_memory=False)
    for required_col in ("harmonized_lineage", "donor_id", "cohort", "cell_id", "event"):
        if required_col not in cells.columns:
            raise ValueError(f"Expected column '{required_col}' in {raw_cells_path}, "
                              f"got {list(cells.columns)}")

    validate_event_column(cells["event"])  # raises on any unexpected value; result not needed yet

    if args.selection_manifest is not None:
        # 2026-10-09 review: a dry-run's selection was computed then immediately
        # deleted, so there was no way to bind a real extraction to a specific
        # REVIEWED list -- only to "the same seed would deterministically reproduce
        # this", which does not protect against the source data itself having
        # drifted between review and extraction. This path consumes a frozen
        # manifest and re-validates every row against the CURRENT cells.csv before
        # trusting it.
        log.info("Consuming frozen selection manifest %s instead of recomputing", args.selection_manifest)
        frozen = pd.read_csv(args.selection_manifest)
        missing_cols = [c for c in SELECTION_MANIFEST_COLS if c not in frozen.columns]
        if missing_cols:
            raise ValueError(f"--selection-manifest {args.selection_manifest} missing "
                              f"column(s) {missing_cols}")
        idx_in_full_file = frozen["orig_row_index"].to_numpy()
        if len(idx_in_full_file) == 0 or idx_in_full_file.max() >= len(cells):
            raise ValueError(f"--selection-manifest {args.selection_manifest} has an "
                              f"orig_row_index out of range for the current {raw_cells_path} "
                              f"({len(cells)} rows) -- refusing to extract against a selection "
                              f"that cannot even index into the current source data.")
        current = cells.iloc[idx_in_full_file]
        drifted = (
            (current["cell_id"].to_numpy() != frozen["cell_id"].to_numpy())
            | (current["cohort"].to_numpy() != frozen["cohort"].to_numpy())
            | (current["donor_id"].to_numpy() != frozen["donor_id"].to_numpy())
            | (current["harmonized_lineage"].to_numpy() != frozen["harmonized_lineage"].to_numpy())
        )
        if drifted.any():
            raise ValueError(
                f"{int(drifted.sum())} row(s) of --selection-manifest {args.selection_manifest} "
                f"no longer match {raw_cells_path} at their recorded orig_row_index -- the "
                f"source data has drifted since this selection was produced; refusing to "
                f"extract against a stale selection. Re-run --dry-run to produce a fresh one."
            )
        selected = current.copy()
        log.info("Consumed and re-validated %d rows from %s against the current cells.csv",
                  len(selected), args.selection_manifest)
    else:
        lineage_mask = cells["harmonized_lineage"].isin(args.lineages)
        log.info("%d / %d cells match requested lineages %s", int(lineage_mask.sum()),
                  len(cells), args.lineages)
        cells_sub = cells[lineage_mask].reset_index(drop=True)
        orig_positions = np.where(lineage_mask.to_numpy())[0]

        idx_within_sub = patient_balanced_index_per_lineage(cells_sub, args.n_per_donor, args.seed)
        idx_in_full_file = orig_positions[idx_within_sub]
        selected = cells.iloc[idx_in_full_file].copy()
    per_cohort_lineage = selected.groupby(["cohort", "harmonized_lineage"]).size()
    n_monocyte = int(per_cohort_lineage.loc[(slice(None), "Monocyte_Myeloid")].sum()) \
        if "Monocyte_Myeloid" in per_cohort_lineage.index.get_level_values(1) else 0
    n_platelet = int(per_cohort_lineage.loc[(slice(None), "Platelet")].sum()) \
        if "Platelet" in per_cohort_lineage.index.get_level_values(1) else 0
    log.info("Selection (%s): %d distinct (cohort,donor,lineage) groups -> %d cells total "
              "(%d Platelet + %d Monocyte_Myeloid; estimated output matrix size, computed "
              "BEFORE the streamed matrix read)",
              f"consumed from {args.selection_manifest}" if args.selection_manifest is not None
              else f"<={args.n_per_donor} cells/donor/lineage, seed={args.seed}",
              len(selected.groupby(["cohort", "donor_id", "harmonized_lineage"])),
              len(idx_in_full_file), n_platelet, n_monocyte)
    print("\n=== ESTIMATED OUTPUT SIZE (cells.csv only -- no counts.mtx read) ===")
    print(f"{len(idx_in_full_file)} cells ({n_platelet} Platelet + {n_monocyte} Monocyte_Myeloid) "
          f"x (n_genes from {genes_path})")
    print(per_cohort_lineage.to_string())

    selection_manifest = selected[["cell_id", "cohort", "donor_id", "harmonized_lineage"]].copy()
    selection_manifest["orig_row_index"] = idx_in_full_file
    selection_manifest.to_csv(manifest_path, index=False)
    # 2026-10-09 review: a dry-run previously computed this manifest then deleted it,
    # leaving no durable, reviewable selection artifact -- only a claim that the same
    # seed would reproduce it. Fixed: the manifest is always kept; a sidecar status
    # file distinguishes "preflight_only" (no matrices extracted yet -- this is what a
    # --dry-run produces) from "extracted" (a real run consumed/produced it), so a
    # reader cannot mistake one for the other.
    base_status = dict(
        n_rows=len(selection_manifest),
        source_cells_csv_sha256=_file_sha256(raw_cells_path),
        produced_by_command=" ".join(sys.argv),
    )

    if args.dry_run:
        _atomic_write_json(status_path, dict(base_status, status="preflight_only"))
        log.info("Wrote selection manifest %s (%d rows) and status preflight_only (%s)",
                  manifest_path, len(selection_manifest), status_path)
        log.info("--dry-run: stopping here. The counts.mtx streaming read below (%s, %s) "
                  "was NOT performed. The selection manifest above is the exact, reviewable "
                  "selection -- pass it to a real run via --selection-manifest to bind that "
                  "run to this specific reviewed list (re-validated against cells.csv at "
                  "that time), or omit --selection-manifest to let a real run recompute "
                  "deterministically from the same --seed instead.")
        return

    # 2026-10-10 review: this previously published status="extracted" HERE, before
    # either counts.mtx had been streamed or the h5ad written/reloaded -- a process
    # killed or raising during the work below left an "extracted" sidecar pointing
    # at an h5ad that was never written, or only partially written. Fixed:
    # "selected" is published first (true at this point -- nothing downstream has
    # happened yet); "extracted" is published only after streaming, the h5ad
    # write/reload check, and the reconciliation checks below all succeed. Any
    # exception in between is caught, published as "failed" with the error
    # recorded, and re-raised (the process still exits non-zero).
    _atomic_write_json(status_path, dict(base_status, status="selected"))
    log.info("Wrote selection manifest %s (%d rows); status sidecar %s: 'selected' "
              "(extraction not yet attempted)", manifest_path, len(selection_manifest), status_path)

    try:
        log.info("Streaming %s (this is the bounded, column-filtered read -- the full file "
                  "is never materialized)", raw_mtx_path)
        raw = read_mtx_columns(str(raw_mtx_path), idx_in_full_file)
        log.info("Streaming %s", decontx_mtx_path)
        decontx = read_mtx_columns(str(decontx_mtx_path), idx_in_full_file)

        genes = pd.read_csv(genes_path)
        gene_col = "gene" if "gene" in genes.columns else genes.columns[0]
        var = pd.DataFrame(index=genes[gene_col].astype(str).tolist())

        meta = selected.reset_index(drop=True)
        meta["lineage"] = meta["harmonized_lineage"]
        # 2026-10-09 review: this was `meta["donor_id"]` -- NOT cohort-qualified, so a
        # donor_id string coincidentally reused across two different cohorts would be
        # treated as ONE patient by run_batch_diagnostics.py's patient-balanced sampler,
        # which groups by this column globally (not cohort-qualified). The sampler used
        # to build THIS selection already qualifies by (cohort, donor_id, lineage) --
        # this line was the one place that qualifier was dropped before persisting.
        meta["patient_id"] = meta["cohort"].astype(str) + "::" + meta["donor_id"].astype(str)
        meta["outcome"] = validate_event_column(meta["event"]).values
        meta.index = meta["cell_id"].astype(str)
        if not meta.index.is_unique:
            raise ValueError("cell_id is not unique in the patient-balanced subsample -- "
                              "refusing to build an h5ad with ambiguous cell identity.")

        adata = anndata.AnnData(X=decontx.T.tocsr().astype(np.float32), obs=meta, var=var)
        adata.layers["counts_raw"] = raw.T.tocsr().astype(np.float32)
        adata.layers["counts_decontx"] = decontx.T.tocsr().astype(np.float32)
        adata.layers["counts"] = adata.layers["counts_decontx"]  # explicit alias: this IS what X is
        adata.uns["population_stage"] = "pre-Gate-2 PROVISIONAL, not final (see docs/current-status.md)"
        adata.uns["population_stage_id"] = STAGE_ID  # machine-checkable; see 00_shared/code/population_stage.py
        adata.uns["expression_source"] = "decontx_counts"
        adata.uns["source_root"] = str(args.input_root)
        adata.uns["sampling_rule"] = (
            f"consumed from {args.selection_manifest}" if args.selection_manifest is not None
            else f"<= {args.n_per_donor} cells/(cohort,donor,lineage), seed={args.seed}, "
                 f"stable_seed.derive_seed"
        )
        adata.uns["selection_manifest_path"] = str(manifest_path)

        adata.write_h5ad(args.out)
        reloaded = anndata.read_h5ad(args.out)
        if reloaded.shape != adata.shape:
            raise ValueError(f"write/reload shape mismatch: {reloaded.shape} != {adata.shape}")
        if not (reloaded.obs_names == adata.obs_names).all():
            raise ValueError("write/reload obs_names (cell_id) order mismatch")

        # Job 1's own completion record -- distinct from run_batch_diagnostics.py's
        # (Job 2's) manifest.json, which this script never writes. A 2026-10-10
        # review found an earlier handoff had pointed Job 1's acceptance criteria at
        # fields only Job 2 produces. These two checks are RAISED on, not merely
        # recorded, if violated -- a recorded-but-unenforced boolean is not a guard.
        if not all(pid == f"{c}::{d}" for pid, c, d in
                   zip(meta["patient_id"], meta["cohort"], meta["donor_id"])):
            raise ValueError("patient_id does not match the expected cohort::donor_id form")
        n_platelet_final = int((meta["lineage"] == "Platelet").sum())
        n_monocyte_final = int((meta["lineage"] == "Monocyte_Myeloid").sum())
        if n_platelet_final + n_monocyte_final != len(meta):
            raise ValueError("per-lineage counts do not reconcile with the total selection")

        final_status = dict(
            base_status, status="extracted",
            output_h5ad=str(args.out), output_sha256=_file_sha256(args.out),
            n_cells=int(adata.shape[0]), n_genes=int(adata.shape[1]),
            stage_id=STAGE_ID, expression_source="decontx_counts",
            cell_id_order_matches_selection=bool((reloaded.obs_names == adata.obs_names).all()),
            patient_id_cohort_qualified=True,
            per_lineage_counts_reconciled=True,
        )
        _atomic_write_json(status_path, final_status)
    except Exception as exc:
        _atomic_write_json(status_path, dict(base_status, status="failed", error=str(exc)))
        log.error("Extraction failed after selection; status sidecar %s published as "
                  "'failed', not 'extracted': %s", status_path, exc)
        raise

    log.info("Wrote and verified %s: %s", args.out, adata)
    log.info("Status sidecar %s: 'extracted' (published only after all checks above succeeded)",
              status_path)
    print(f"\n=== PROVISIONAL DIAGNOSTIC INPUT: {adata.shape[0]} cells x {adata.shape[1]} genes ===")
    print(meta.groupby(["cohort", "lineage"]).size().to_string())


if __name__ == "__main__":
    main()
