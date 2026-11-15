#!/usr/bin/env python3
"""Stephenson (E-MTAB-10026) pool-membership crosswalk: SDRF (deposited
declaration) vs this project's existing local table
(`stephenson_pool_composition.csv` / `canonical_donor_labels.csv` "notes")
-- the two DECLARED sources -- vs three OBSERVED, per-cell sources: the
REAL processed GEX h5ad's own `obs['sample_id']`, BCR, and TCR -- for the
Cambridge-site pools, plus a separate Sanger-site AP-patient observed-
membership table.

v3 (this file) responds to review:20261011_113512, which reviewed
v2 and found further defects, fixed here:

  F2a: `_read_legacy_categorical` indexed `cats[codes]` directly; a legacy
      missing-category code (-1) silently returned the LAST category
      instead of a real missing value. Now returns `None` for `-1` and
      raises on any other out-of-range code.
  F2b: `gex_processed` was in `DECLARED_SOURCES`, treating the real GEX
      h5ad's per-cell label as an independent declaration of full pool
      membership rather than one more OBSERVED, processed label that can
      legitimately be a subset of a declared roster (QC/filtering loss).
      Moved to `OBSERVED_SOURCES`.
  F2c/F5: `parse_gex_cambridge` silently overwrote a repeated CellID's
      donor label; now collects candidate donors per CellID and reports
      a real per-cell conflict (mirroring BCR/TCR) instead of overwriting.
  F5a: `classify_multi_pool_donor` treated "two flags, two days" as
      sufficient for `d0_d28_style_resample` without checking the flags
      and days are actually PAIRED the expected way (`Initial`+`D0` with
      `Resample`+`D28`); an inverted pairing (`Initial`+`D28` with
      `Resample`+`D0`) would have been misclassified as a resample. Now
      checks the actual (flag, day) pairs, not just the two sets.
  F1: `cross_source_cell_join` is extended with a `combined` block that
      reports the ACTUAL union/overlap between the GEX-BCR and GEX-TCR
      shared-cell sets, rather than letting a caller add their sizes
      (21,481 + 74,761) as if they were disjoint and calling the sum a
      single "96,242 distinct cells" figure -- which the review
      caught as unsupported (the decision sheet never checked whether the
      two intersections overlap).
  F5b: `--out-dir`'s check-then-create had a TOCTOU window between
      `Path.exists()` and `mkdir()`. Now attempts an exclusive `mkdir`
      first and only inspects contents if that fails because the
      directory already existed.

v2 (this file) responds to review:20261010_202340, which reviewed
v1 (commits `400126c`/`ed3c634`) and found 5 real defects, fixed here:

  F1: a donor seen in >1 Cambridge pool was reported as a "resample" even
      when both observations were the SAME day/flag (`CV0902`:
      `BGCV01`+`BGCV02`, both `Initial`/`D0`) -- `classify_multi_pool_donor`
      now reports a category per donor instead of a blanket "resample".
  F2: ANY two-source disagreement was reported as `conflict=True`, even
      when one source (BCR/TCR) had simply not observed every declared
      member (coverage gap, not a disagreement). `classify_pool_evidence`
      now distinguishes `declared_membership_disagreement` (two DECLARED
      sources -- SDRF, the local table, or the real GEX's own per-cell
      `sample_id` -- give different full rosters) from
      `observed_subset_of_declared` (an OBSERVED source, BCR/TCR, is a
      subset of some declared roster: a coverage gap) from
      `observed_contradicts_declared` (an OBSERVED source has a member no
      DECLARED source agrees with).
  F3: no actual GEX/annotation/BCR/TCR cell-level join existed -- added
      `parse_gex_cambridge` (the real processed h5ad, legacy `__categories`
      AnnData encoding, read via h5py only -- never materializes `X`/
      `layers`) and `cross_source_cell_join`, which joins on the FULL
      CellID (e.g. `BGCV06_AAACCTGAGCCTATGT-1`) across GEX/annotation/BCR/
      TCR -- never a bare, cross-library 10x barcode.
  F4/F5: addressed in `validate_job1_output.py` and this file's new
      no-clobber `--out-dir` guard + `write_source_manifest` (sha256 of
      every input, this script's own git revision, full argv) --
      see the companion handoff for exact diffs.

Streams BCR (gzip CSV, ~13MB) and TCR (247MB TSV, 547,631 rows) row-by-row
by column NAME, never loading either into memory whole. Never reads
Outcome/Status columns from any source for any identity decision, per this
project's outcome-data boundary (DEVELOPMENT.md).
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "00_shared/code"))
from run_ownership import claim_ownership  # noqa: E402

# F2 correction (review:20261011_113512): `gex_processed` moved
# OUT of DECLARED_SOURCES. The real processed GEX h5ad's `sample_id` is a
# per-cell OBSERVED, processed label, not an independent declaration of
# full pool membership -- if upstream QC/filtering dropped some donor's
# cells before this file was built, GEX's observed roster would be a
# genuine (and expected) SUBSET of the declared roster, not a conflict.
# Treating it as declared would misclassify that coverage gap as a
# `declared_membership_disagreement`, exactly as BCR/TCR's partial
# coverage would be if they were declared sources.
DECLARED_SOURCES = ("sdrf", "local_pool_composition")
OBSERVED_SOURCES = ("bcr", "tcr_cambridge", "gex_processed")


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(1 << 20)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def _git_revision_of(dir_path: Path) -> str:
    try:
        out = subprocess.run(
            ["git", "-C", str(dir_path), "rev-parse", "HEAD"],
            capture_output=True, text=True, timeout=10,
        )
        return out.stdout.strip() if out.returncode == 0 else "UNVERIFIED (git rev-parse failed)"
    except Exception as exc:  # pragma: no cover -- environment-dependent
        return f"UNVERIFIED ({exc})"


def _git_dirty_status_of(dir_path: Path) -> str:
    """R5 (review:20261011_124253): a bare HEAD hash does not prove
    the file actually executed matched that commit -- a hand-patched
    checkout can still report a clean-looking HEAD. Records
    `git status --porcelain` too, so a dirty working tree is visible in the
    manifest itself rather than requiring a separate SSH check."""
    try:
        out = subprocess.run(
            ["git", "-C", str(dir_path), "status", "--porcelain"],
            capture_output=True, text=True, timeout=10,
        )
        if out.returncode != 0:
            return "UNVERIFIED (git status failed)"
        return "clean" if not out.stdout.strip() else f"DIRTY:\n{out.stdout.strip()}"
    except Exception as exc:  # pragma: no cover -- environment-dependent
        return f"UNVERIFIED ({exc})"


def _script_git_revision() -> str:
    return _git_revision_of(Path(__file__).resolve().parent)


def write_source_manifest(out_dir: Path, input_paths: dict[str, Path], argv: list[str],
                           script_path: Path | None = None,
                           dependency_module_paths: tuple[Path, ...] = (),
                           run_ownership: dict | None = None) -> Path:
    """R5 fix: previously recorded only `git rev-parse HEAD` for THIS
    module, which does not prove the file that actually ran matched that
    commit (a hand-patched checkout can still report a clean HEAD) and said
    nothing about modules this run imports from. Now also records this
    run's own script sha256, the working tree's dirty state, and the
    sha256 of every `dependency_module_paths` entry (e.g. a caller script
    that imports functions from this module records THIS file's sha256
    too, not just its own)."""
    script_path = Path(script_path).resolve() if script_path is not None else Path(__file__).resolve()
    manifest = {
        "script": str(script_path),
        "script_sha256": file_sha256(script_path),
        "script_git_revision": _git_revision_of(script_path.parent),
        "script_git_dirty_status": _git_dirty_status_of(script_path.parent),
        "argv": argv,
        "dependency_modules": {
            str(p): {"sha256": file_sha256(p)} for p in dependency_module_paths
        },
        "run_ownership": run_ownership,
        "inputs": {
            name: {"path": str(p), "sha256": file_sha256(p)}
            for name, p in input_paths.items()
        },
    }
    manifest_path = out_dir / "source_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    return manifest_path


def parse_sdrf_pool_membership(sdrf_path: Path) -> dict[str, set[str]]:
    """{pool_label: {donor_id, ...}} from the SDRF's own
    `Characteristics[individual]`/`Characteristics[sequence pool]` -- the
    officially deposited declaration. Read by column NAME."""
    pool_donors: dict[str, set[str]] = defaultdict(set)
    with open(sdrf_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"Characteristics[individual]", "Characteristics[sequence pool]"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{sdrf_path}: missing required column(s) {missing}")
        for row in reader:
            donor = row["Characteristics[individual]"].strip()
            pool = row["Characteristics[sequence pool]"].strip()
            if donor and pool:
                pool_donors[pool].add(donor)
    return dict(pool_donors)


def parse_local_pool_composition(csv_path: Path) -> dict[str, set[str]]:
    pool_donors: dict[str, set[str]] = {}
    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            pool_donors[row["pool"].strip()] = {
                d.strip() for d in row["donors"].split(";") if d.strip()
            }
    return pool_donors


def parse_canonical_donor_pools(csv_path: Path) -> dict[str, set[str]]:
    """'pools: nan' is this project's OWN sentinel for "never looked up"
    (Sanger/Newcastle donors) -- excluded, not compared as an empty claim."""
    donor_pools: dict[str, set[str]] = {}
    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("cohort") != "Stephenson 2021":
                continue
            notes = row.get("notes", "")
            if not notes.startswith("pools:"):
                continue
            pools = {p.strip() for p in notes[len("pools:"):].split(";") if p.strip()}
            if pools == {"nan"}:
                continue
            donor_pools[row["donor_id"].strip()] = pools
    return donor_pools


def _read_legacy_categorical(f: h5py.File, key: str, group: str = "obs") -> np.ndarray:
    """E-MTAB-10026's deposited h5ad uses the LEGACY AnnData 0.1.0 encoding
    (`gate2_source_crosswalk.py`'s documented concern): a plain int
    codes dataset at `<group>/<key>` with labels in the SIBLING group
    `<group>/__categories/<key>` -- NOT the modern Group-with-categories+
    codes encoding `validate_job1_output.py` reads. Decoding this as the
    modern encoding would silently yield raw integers, not labels.

    Fixes F2 (review:20261011_113512): `codes == -1` is the
    legacy pandas/AnnData sentinel for a MISSING category, not a real
    index. The previous `cats[codes]` indexed `cats[-1]` for any such row
    -- silently returning the LAST category as if it were real. Missing
    codes now become `None`; any other out-of-range code (neither a valid
    index into `cats` nor exactly -1) raises rather than wrapping or
    guessing."""
    grp = f[group]
    if "__categories" in grp and key in grp["__categories"]:
        codes = grp[key][:]
        cats = grp["__categories"][key][:]
        cats = np.array([c.decode() if isinstance(c, bytes) else c for c in cats])
        out = np.empty(len(codes), dtype=object)
        missing = codes == -1
        valid = (~missing) & (codes >= 0) & (codes < len(cats))
        bad = ~(missing | valid)
        if np.any(bad):
            bad_values = sorted(set(int(c) for c in codes[bad]))[:10]
            raise ValueError(
                f"{group}/{key}: {int(bad.sum())} code(s) out of range [0, {len(cats)}) "
                f"and != -1 (the legacy missing sentinel); examples: {bad_values}"
            )
        out[valid] = cats[codes[valid]]
        out[missing] = None
        return out
    arr = grp[key][:]
    if arr.dtype.kind in ("S", "O"):
        return np.array([v.decode() if isinstance(v, bytes) else v for v in arr])
    return arr


def parse_gex_cambridge(h5ad_path: Path):
    """Reads ONLY `obs['sample_id']`, `obs['Site']`, and the obs index
    (`covid_index`) via h5py -- never `X`/`layers` (the real counts, tens
    of GB). `sample_id` for Cambridge cells is the full `<pool>_<donor>`
    label (e.g. `BGCV06_CV0201`) this project's own GEX-processing
    pipeline assigned to that real, already-sequenced cell.

    F1 correction (review:20261011_113512): this is NOT
    independent biological ground truth -- it is one more OBSERVED,
    processed label, which may itself derive from the same donor-
    demultiplexing process BCR/TCR inherit. Cross-file agreement with
    BCR/TCR is internal-consistency evidence, not proof of independent
    identity. `gex_processed` is accordingly in `OBSERVED_SOURCES`, not
    `DECLARED_SOURCES` (see module-level comment).

    Returns (gex_pool_donors, gex_pool_cell_counts, gex_cell_donor,
    cell_id_conflicts, coverage) where gex_cell_donor is {full_CellID:
    donor} for Cambridge cells with a SINGLE consistent donor label
    (bounded: tens of thousands of entries, not the full 647k-cell file).
    F2 fix: a CellID seen more than once is no longer silently overwritten
    -- if its donor label disagrees across occurrences it is reported in
    `cell_id_conflicts` and excluded from `gex_cell_donor`, mirroring how
    `parse_bcr_pool_membership`/`parse_tcr_cambridge_and_sanger` already
    handle a repeated barcode/CellID. `coverage` reports how many rows had
    a missing (legacy-encoding `-1`) Site/sample_id/index value."""
    with h5py.File(h5ad_path, "r") as f:
        site = _read_legacy_categorical(f, "Site")
        sample_id = _read_legacy_categorical(f, "sample_id")
        index_key = f["obs"].attrs.get("_index", "_index")
        index_key = index_key.decode() if isinstance(index_key, bytes) else index_key
        cell_ids = _read_legacy_categorical(f, index_key)

    cambridge_mask = site == "Cambridge"  # None (missing Site) compares False, so excluded
    labels = sample_id[cambridge_mask]
    ids = cell_ids[cambridge_mask]

    pool_donors: dict[str, set[str]] = defaultdict(set)
    pool_cell_counts: dict[str, int] = defaultdict(int)
    cell_candidate_donors: dict[str, set[str]] = defaultdict(set)
    n_missing_cell_id = 0
    n_missing_label = 0
    for cell_id, label in zip(ids, labels):
        if cell_id is None:
            n_missing_cell_id += 1
            continue
        if label is None:
            n_missing_label += 1
            continue
        if "_" not in label:
            continue  # not a <pool>_<donor> label; skip rather than guess
        pool, donor = label.split("_", 1)
        pool_donors[pool].add(donor)
        pool_cell_counts[pool] += 1
        cell_candidate_donors[cell_id].add(donor)

    cell_id_conflicts = [(cid, ds) for cid, ds in cell_candidate_donors.items() if len(ds) > 1]
    cell_donor = {cid: next(iter(ds)) for cid, ds in cell_candidate_donors.items() if len(ds) == 1}
    coverage = {
        "n_missing_site_or_cell_id": n_missing_cell_id,
        "n_missing_sample_id_label": n_missing_label,
        "n_cell_id_conflicts": len(cell_id_conflicts),
    }
    return dict(pool_donors), dict(pool_cell_counts), cell_donor, cell_id_conflicts, coverage


def parse_annotation_cellids(csv_path: Path) -> set[str]:
    """`annotation_200112.csv`: index column is the full CellID
    (`<pool>_<barcode>-1`), second column `full_clustering` (cell type) --
    no donor field (confirmed). Streamed row-by-row; only the CellID is
    kept (a set of strings, bounded by the file's real row count)."""
    cell_ids: set[str] = set()
    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        if len(header) < 1:
            raise ValueError(f"{csv_path}: empty header")
        for row in reader:
            if row:
                cell_ids.add(row[0])
    return cell_ids


def _bare_barcode_pool_prefix(full_barcode: str) -> tuple[str, str] | None:
    if "_" not in full_barcode:
        return None
    pool, bare = full_barcode.split("_", 1)
    return pool, bare


def parse_bcr_pool_membership(bcr_gz_path: Path):
    """Returns (pool_donors, per_cell_conflicts, cell_donor) -- cell_donor
    is {full_barcode: patient_id} for barcodes with a SINGLE consistent
    patient_id across all their chain rows (conflicting ones are excluded
    here and reported only in per_cell_conflicts, never silently picked)."""
    pool_donors: dict[str, set[str]] = defaultdict(set)
    cell_patient_ids: dict[str, set[str]] = defaultdict(set)
    with gzip.open(bcr_gz_path, "rt", newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        try:
            barcode_i = header.index("barcode")
            patient_i = header.index("patient_id")
        except ValueError as exc:
            raise ValueError(f"{bcr_gz_path}: expected 'barcode' and 'patient_id' columns") from exc
        for row in reader:
            full_barcode = row[barcode_i]
            patient_id = row[patient_i]
            if not patient_id:
                continue
            split = _bare_barcode_pool_prefix(full_barcode)
            if split is not None:
                pool, _bare = split
                pool_donors[pool].add(patient_id)
            cell_patient_ids[full_barcode].add(patient_id)

    conflicts = [(b, ids) for b, ids in cell_patient_ids.items() if len(ids) > 1]
    cell_donor = {b: next(iter(ids)) for b, ids in cell_patient_ids.items() if len(ids) == 1}
    return dict(pool_donors), conflicts, cell_donor


def parse_tcr_cambridge_and_sanger(tcr_tsv_path: Path):
    """Returns (cambridge_pool_donors, cambridge_resample_entries,
    sanger_membership, cambridge_cell_conflicts, cambridge_cell_donor).
    `CellID` (the real full per-cell identifier, distinct from `Sample`) is
    now read and used for a genuine per-cell donor-consistency check,
    mirroring the BCR check -- F3's "TCR 同一细胞的 donor 冲突" request."""
    cambridge_pool_donors: dict[str, set[str]] = defaultdict(set)
    cambridge_resample: dict[str, list[tuple]] = defaultdict(list)
    sanger_membership: dict[str, set[str]] = defaultdict(set)
    seen_resample_rows: set[tuple] = set()
    cell_patient_ids: dict[str, set[str]] = defaultdict(set)

    with open(tcr_tsv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"CellID", "Sample", "patient_id", "Site", "Resample", "Collection_Day"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{tcr_tsv_path}: missing required column(s) {missing}")
        for row in reader:
            site = row["Site"]
            sample = row["Sample"]
            patient_id = row["patient_id"]
            cell_id = row["CellID"]
            if not sample or not patient_id:
                continue
            if site == "Cambridge":
                cambridge_pool_donors[sample].add(patient_id)
                if cell_id:
                    cell_patient_ids[cell_id].add(patient_id)
                key = (patient_id, sample, row["Resample"], row["Collection_Day"])
                if key not in seen_resample_rows:
                    seen_resample_rows.add(key)
                    cambridge_resample[patient_id].append((sample, row["Resample"], row["Collection_Day"]))
            elif site == "Sanger":
                sanger_membership[sample].add(patient_id)

    cell_conflicts = [(c, ids) for c, ids in cell_patient_ids.items() if len(ids) > 1]
    cell_donor = {c: next(iter(ids)) for c, ids in cell_patient_ids.items() if len(ids) == 1}
    return (dict(cambridge_pool_donors), dict(cambridge_resample), dict(sanger_membership),
            cell_conflicts, cell_donor)


def _day_sort_key(day: str):
    """`Collection_Day` values look like 'D0'/'D28'; sort numerically on
    the trailing digits when possible so 'D7' sorts before 'D14' (a plain
    lexical sort would not). Falls back to the raw string, grouped after
    any numerically-parseable value, rather than guessing."""
    digits = "".join(ch for ch in day if ch.isdigit())
    if digits:
        return (0, int(digits))
    return (1, day)


def classify_multi_pool_donor(entries: list[tuple]) -> str:
    """`entries`: [(sample_pool, resample_flag, collection_day), ...] for
    ONE donor seen in >1 distinct Sample pool (TCR evidence). Fixes F1:
    previously every such donor was called a "resample" even when both
    observations shared the same day AND flag (`CV0902`:
    `Initial`/`D0` + `Initial`/`D0` in two different pools).

    F5a fix (review:20261011_113512): the previous version
    accepted ANY combination with flags == {"Initial", "Resample"} and 2
    distinct days as a "D0/D28-style resample", without checking the
    flags and days are actually PAIRED as `Initial`+earlier-day with
    `Resample`+later-day. An inverted pairing (`Initial`/`D28` with
    `Resample`/`D0`) is a different, suspicious pattern and must not be
    silently called the same resample shape."""
    distinct_pools = {e[0] for e in entries}
    if len(distinct_pools) < 2:
        return "single_pool"
    flags_days = {(e[1], e[2]) for e in entries}
    if len(flags_days) == 1:
        return "same_day_multi_pool"  # e.g. CV0902: two pools, both Initial/D0
    if len(distinct_pools) == 2 and len(flags_days) == 2:
        flags = {p[0] for p in flags_days}
        days = {p[1] for p in flags_days}
        if flags == {"Initial", "Resample"} and len(days) == 2:
            initial_day = next(p[1] for p in flags_days if p[0] == "Initial")
            resample_day = next(p[1] for p in flags_days if p[0] == "Resample")
            if _day_sort_key(initial_day) < _day_sort_key(resample_day):
                return "d0_d28_style_resample"
    return "inconsistent_day_flag_combination"  # needs explicit review, not silently bucketed


def classify_pool_evidence(pool: str, per_source: dict[str, set[str] | None]) -> dict:
    """Fixes F2. Returns a dict with `conflict_type` (one of 'none',
    'declared_membership_disagreement', 'observed_contradicts_declared',
    'insufficient_declared_reference') plus per-observed-source coverage/
    contradiction detail. A DECLARED source (sdrf, local_pool_composition)
    disagreeing with another DECLARED source is the only thing called a
    membership disagreement. An OBSERVED source (bcr, tcr_cambridge,
    gex_processed) that is merely a SUBSET of some declared roster is a
    coverage gap, not a conflict; one with a member no declared source
    has is a real contradiction, reported separately from a subset gap."""
    declared = {k: per_source[k] for k in DECLARED_SOURCES if per_source.get(k) is not None}
    observed = {k: per_source[k] for k in OBSERVED_SOURCES if per_source.get(k) is not None}

    declared_sets = {frozenset(v) for v in declared.values()}
    declared_disagree = len(declared_sets) > 1
    result = {
        "pool": pool,
        "conflict_type": "declared_membership_disagreement" if declared_disagree else "none",
        "declared_sources_observed": sorted(declared),
    }

    reference_sets = declared if declared else {}
    union_of_declared = set().union(*reference_sets.values()) if reference_sets else set()

    coverage = {}
    contradictions = {}
    for obs_name, obs_set in observed.items():
        if not reference_sets:
            coverage[obs_name] = "no_declared_reference"
            continue
        subset_of_any = [name for name, ref in reference_sets.items() if obs_set <= ref]
        extra = obs_set - union_of_declared
        if subset_of_any:
            best = max(subset_of_any, key=lambda n: len(reference_sets[n] - obs_set) * -1)
            missing = sorted(reference_sets[best] - obs_set)
            coverage[obs_name] = "exact_match" if not missing else f"subset_of_{best};missing={missing}"
        else:
            coverage[obs_name] = "not_a_subset_of_any_declared_source"
        if extra:
            contradictions[obs_name] = sorted(extra)

    result["observed_coverage"] = coverage
    result["observed_contradictions"] = contradictions
    if contradictions and result["conflict_type"] == "none":
        result["conflict_type"] = "observed_contradicts_declared"
    if not reference_sets and observed:
        result["conflict_type"] = "insufficient_declared_reference"
    return result


def build_identity_conflicts(sources: dict[str, dict[str, set[str]]]) -> list[dict]:
    all_pools = sorted({pool for d in sources.values() for pool in d})
    rows = []
    for pool in all_pools:
        per_source = {name: d.get(pool) for name, d in sources.items()}
        classification = classify_pool_evidence(pool, per_source)
        row = {
            "pool": pool,
            "conflict_type": classification["conflict_type"],
            "conflict": classification["conflict_type"] != "none",
        }
        for name in sources:
            donors = per_source[name]
            row[f"donors_{name}"] = ";".join(sorted(donors)) if donors is not None else f"not_observed_in_{name}"
        row["observed_coverage"] = json.dumps(classification["observed_coverage"])
        row["observed_contradictions"] = json.dumps(classification["observed_contradictions"])
        rows.append(row)
    return rows


def build_dependency_components(pool_donors: dict[str, set[str]]) -> dict[str, set[str]]:
    """Union-find over the pool-membership graph: two donors are in the
    same component iff they share at least one pool, transitively. Returns
    {representative_donor: {member, ...}}. Pure stdlib, no canonical
    dataset is modified -- this is a read-only comparison tool."""
    parent: dict[str, str] = {}

    def find(x: str) -> str:
        root = x
        while parent.get(root, root) != root:
            root = parent[root]
        while parent.get(x, x) != root:
            parent[x], x = root, parent.get(x, x)
        return root

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for donors in pool_donors.values():
        donor_list = sorted(donors)
        for d in donor_list:
            parent.setdefault(d, d)
        for i in range(1, len(donor_list)):
            union(donor_list[0], donor_list[i])

    groups: dict[str, set[str]] = defaultdict(set)
    for d in parent:
        groups[find(d)].add(d)
    return dict(groups)


def exclude_donors_from_pools(pool_donors: dict[str, set[str]], donors: set[str]) -> dict[str, set[str]]:
    """Returns a new pool-membership mapping with `donors` removed from
    EVERY pool's member set. Used to actually test decision sheet §2's
    option 3 ("exclude the disputed donors from any pool-dependent
    grouping decision") via `diff_component_mappings`, rather than
    asserting its consequence without computing it."""
    return {pool: members - donors for pool, members in pool_donors.items()}


def diff_component_mappings(pools_a: dict[str, set[str]], pools_b: dict[str, set[str]],
                             label_a: str, label_b: str) -> dict:
    """Compares the connected-component structure built from two different
    pool-membership rosters for the SAME 15 Cambridge pools (e.g. SDRF's
    declared membership vs. the local/BCR-matching membership) -- member-
    by-member and edge-by-edge, not just component counts/sizes."""
    comps_a = build_dependency_components(pools_a)
    comps_b = build_dependency_components(pools_b)
    all_donors = sorted(set().union(*pools_a.values()) | set().union(*pools_b.values()))

    def component_of(donor: str, comps: dict[str, set[str]]) -> frozenset:
        for members in comps.values():
            if donor in members:
                return frozenset(members)
        return frozenset({donor})

    donor_rows = []
    for donor in all_donors:
        comp_a = component_of(donor, comps_a)
        comp_b = component_of(donor, comps_b)
        donor_rows.append({
            "donor_id": donor,
            f"component_size_{label_a}": len(comp_a),
            f"component_size_{label_b}": len(comp_b),
            "component_changed": comp_a != comp_b,
            f"component_members_{label_a}": ";".join(sorted(comp_a)),
            f"component_members_{label_b}": ";".join(sorted(comp_b)),
        })
    return {
        "donor_rows": donor_rows,
        f"component_sizes_{label_a}": sorted((len(c) for c in comps_a.values()), reverse=True),
        f"component_sizes_{label_b}": sorted((len(c) for c in comps_b.values()), reverse=True),
        "n_donors_with_changed_component": sum(1 for r in donor_rows if r["component_changed"]),
    }


def cross_source_cell_join(gex_cell_donor: dict[str, str], bcr_cell_donor: dict[str, str],
                            tcr_cell_donor: dict[str, str], annotation_cellids: set[str]) -> dict:
    """Full-CellID join (never a bare barcode) between the real processed
    GEX and each of BCR/TCR/annotation, restricted to Cambridge. Reports
    overlap and per-cell donor agreement -- the real identity evidence
    F3 asked for, not a pool-level set comparison.

    F1 fix (review:20261011_113512): `bcr`'s and
    `tcr_cambridge`'s `n_shared_cell_ids` are two DIFFERENT intersections
    with GEX; a cell can be in both (BCR+TCR both exist for the same
    chain-rich cell) or in neither. Adding the two sizes and calling the
    sum a count of distinct cells (what the prior decision sheet did --
    21,481 + 74,761 = "96,242") is unsupported without knowing their
    overlap. The new `combined` block computes the actual union and
    overlap instead of leaving a caller to add sizes."""
    gex_ids = set(gex_cell_donor)

    def _compare(other_name: str, other_donor: dict[str, str]) -> dict:
        shared = gex_ids & set(other_donor)
        agree = sum(1 for cid in shared if gex_cell_donor[cid] == other_donor[cid])
        disagree_examples = [
            {"cell_id": cid, "gex_donor": gex_cell_donor[cid], f"{other_name}_donor": other_donor[cid]}
            for cid in shared if gex_cell_donor[cid] != other_donor[cid]
        ][:20]
        return {
            "n_gex_cambridge_cells": len(gex_ids),
            f"n_{other_name}_cells": len(other_donor),
            "n_shared_cell_ids": len(shared),
            "n_agree": agree,
            "n_disagree": len(shared) - agree,
            "disagree_examples": disagree_examples,
        }

    annotation_overlap = {
        "n_gex_cambridge_cells": len(gex_ids),
        "n_annotation_cellids": len(annotation_cellids),
        "n_gex_cellids_found_in_annotation": len(gex_ids & annotation_cellids),
        "n_gex_cellids_missing_from_annotation": len(gex_ids - annotation_cellids),
    }

    bcr_shared = gex_ids & set(bcr_cell_donor)
    tcr_shared = gex_ids & set(tcr_cell_donor)
    overlap_bcr_and_tcr_shared = bcr_shared & tcr_shared
    union_shared = bcr_shared | tcr_shared

    def _agrees_with_every_shared_observed_source(cid: str) -> bool:
        checks = []
        if cid in bcr_cell_donor:
            checks.append(gex_cell_donor[cid] == bcr_cell_donor[cid])
        if cid in tcr_cell_donor:
            checks.append(gex_cell_donor[cid] == tcr_cell_donor[cid])
        return all(checks)  # non-empty by construction (cid is in union_shared)

    n_agree_all = sum(1 for cid in union_shared if _agrees_with_every_shared_observed_source(cid))
    combined = {
        "n_bcr_shared_with_gex": len(bcr_shared),
        "n_tcr_shared_with_gex": len(tcr_shared),
        "n_overlap_between_bcr_and_tcr_shared_sets": len(overlap_bcr_and_tcr_shared),
        "n_union_shared_distinct_cells": len(union_shared),
        "n_union_agree_with_every_shared_observed_source": n_agree_all,
        "n_union_disagree_with_at_least_one_shared_observed_source": len(union_shared) - n_agree_all,
        "note": (
            "n_union_shared_distinct_cells is NOT n_bcr_shared_with_gex + "
            "n_tcr_shared_with_gex -- those two sets can overlap; see "
            "n_overlap_between_bcr_and_tcr_shared_sets for the real intersection size."
        ),
    }
    return {
        "bcr": _compare("bcr", bcr_cell_donor),
        "tcr_cambridge": _compare("tcr_cambridge", tcr_cell_donor),
        "annotation_coverage": annotation_overlap,
        "combined": combined,
    }


def main():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--sdrf", required=True, type=Path)
    p.add_argument("--local-pool-composition", required=True, type=Path)
    p.add_argument("--canonical-donor-labels", required=True, type=Path)
    p.add_argument("--gex-h5ad", required=True, type=Path)
    p.add_argument("--annotation-csv", required=True, type=Path)
    p.add_argument("--bcr-gz", required=True, type=Path)
    p.add_argument("--tcr-tsv", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    args = p.parse_args()

    # F5b/R1 fix (review:20261011_124253): the non-empty check alone
    # catches leftover files from a COMPLETED or crashed run, but is check-
    # then-act and does not prevent two FRESH concurrent starts from both
    # passing it before either writes -- the same gap `run_ownership.py`'s
    # docstring documents and `run_batch_diagnostics.py` already closes with
    # `claim_ownership()`'s OS-atomic exclusive create. Reused here rather
    # than reimplemented, for the same guarantee.
    if args.out_dir.exists() and any(args.out_dir.iterdir()):
        raise FileExistsError(
            f"{args.out_dir} already exists and is non-empty -- refusing to overwrite a "
            f"prior result (this project's standing rule); use a fresh, newly timestamped "
            f"--out-dir."
        )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    ownership = claim_ownership(args.out_dir / ".run_owner.json", command=" ".join(sys.argv))

    input_paths = {
        "sdrf": args.sdrf, "local_pool_composition": args.local_pool_composition,
        "canonical_donor_labels": args.canonical_donor_labels, "gex_h5ad": args.gex_h5ad,
        "annotation_csv": args.annotation_csv, "bcr_gz": args.bcr_gz, "tcr_tsv": args.tcr_tsv,
    }
    write_source_manifest(args.out_dir, input_paths, sys.argv, run_ownership=ownership)

    sdrf_pools = parse_sdrf_pool_membership(args.sdrf)
    local_pools = parse_local_pool_composition(args.local_pool_composition)
    canonical_donor_pools = parse_canonical_donor_pools(args.canonical_donor_labels)
    (gex_pools, gex_pool_cell_counts, gex_cell_donor,
     gex_cell_id_conflicts, gex_coverage) = parse_gex_cambridge(args.gex_h5ad)
    annotation_cellids = parse_annotation_cellids(args.annotation_csv)
    bcr_pools, bcr_cell_conflicts, bcr_cell_donor = parse_bcr_pool_membership(args.bcr_gz)
    (tcr_cambridge_pools, tcr_resample, sanger_membership,
     tcr_cell_conflicts, tcr_cell_donor) = parse_tcr_cambridge_and_sanger(args.tcr_tsv)

    sources = {
        "sdrf": sdrf_pools, "local_pool_composition": local_pools,
        "gex_processed": gex_pools, "bcr": bcr_pools, "tcr_cambridge": tcr_cambridge_pools,
    }
    conflict_rows = build_identity_conflicts(sources)
    with open(args.out_dir / "identity_conflicts.csv", "w", newline="") as f:
        fieldnames = (["pool", "conflict_type", "conflict"] + [f"donors_{n}" for n in sources]
                      + ["observed_coverage", "observed_contradictions"])
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(conflict_rows)

    with open(args.out_dir / "sanger_observed_membership.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample", "observed_ap_patients", "n_observed"])
        for sample in sorted(sanger_membership):
            aps = sorted(sanger_membership[sample])
            w.writerow([sample, ";".join(aps), len(aps)])

    with open(args.out_dir / "multi_pool_donor_classification.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["patient_id", "category", "sample_pool", "resample_flag", "collection_day"])
        for patient_id, entries in sorted(tcr_resample.items()):
            category = classify_multi_pool_donor(entries)
            if category == "single_pool":
                continue
            for sample, resample_flag, day in sorted(entries):
                w.writerow([patient_id, category, sample, resample_flag, day])

    with open(args.out_dir / "bcr_per_cell_identity_conflicts.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["full_barcode", "conflicting_patient_ids"])
        for barcode, ids in sorted(bcr_cell_conflicts):
            w.writerow([barcode, ";".join(sorted(ids))])

    with open(args.out_dir / "tcr_per_cell_identity_conflicts.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cell_id", "conflicting_patient_ids"])
        for cell_id, ids in sorted(tcr_cell_conflicts):
            w.writerow([cell_id, ";".join(sorted(ids))])

    with open(args.out_dir / "gex_per_cell_identity_conflicts.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cell_id", "conflicting_donor_labels"])
        for cell_id, ids in sorted(gex_cell_id_conflicts):
            w.writerow([cell_id, ";".join(sorted(ids))])

    donor_level_mismatches = []
    for donor, local_note_pools in sorted(canonical_donor_pools.items()):
        sdrf_pools_for_donor = {pool for pool, donors in sdrf_pools.items() if donor in donors}
        if sdrf_pools_for_donor != local_note_pools:
            donor_level_mismatches.append({
                "donor_id": donor,
                "canonical_labels_notes_pools": ";".join(sorted(local_note_pools)),
                "sdrf_observed_pools": ";".join(sorted(sdrf_pools_for_donor)),
            })
    with open(args.out_dir / "donor_level_sdrf_vs_canonical_notes_mismatches.csv", "w", newline="") as f:
        fieldnames = ["donor_id", "canonical_labels_notes_pools", "sdrf_observed_pools"]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(donor_level_mismatches)

    cell_join = cross_source_cell_join(gex_cell_donor, bcr_cell_donor, tcr_cell_donor, annotation_cellids)
    (args.out_dir / "gex_cell_level_crosswalk.json").write_text(json.dumps(cell_join, indent=2))

    bgcv_sdrf = {p: d for p, d in sdrf_pools.items() if p.startswith("BGCV")}
    bgcv_local = {p: d for p, d in local_pools.items() if p.startswith("BGCV")}
    component_diff = diff_component_mappings(bgcv_sdrf, bgcv_local, "sdrf_declared", "local_bcr_gex_matching")
    with open(args.out_dir / "component_diff_sdrf_vs_local.csv", "w", newline="") as f:
        fieldnames = list(component_diff["donor_rows"][0].keys()) if component_diff["donor_rows"] else []
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(component_diff["donor_rows"])
    (args.out_dir / "component_diff_summary.json").write_text(json.dumps({
        k: v for k, v in component_diff.items() if k != "donor_rows"
    }, indent=2))

    conflicting_pools = [r["pool"] for r in conflict_rows if r["conflict"]]
    declared_disagreement_pools = [
        r["pool"] for r in conflict_rows if r["conflict_type"] == "declared_membership_disagreement"
    ]
    insufficient_reference_pools = [
        r["pool"] for r in conflict_rows if r["conflict_type"] == "insufficient_declared_reference"
    ]

    # F1 (decision sheet §2 option 3, now actually tested rather than
    # asserted "zero consequence" without evidence): what if the disputed
    # donors are removed entirely from the pool-dependency graph instead
    # of being assigned to either side? disputed_donors is derived from
    # the real symmetric difference between SDRF and the local roster for
    # each declared_membership_disagreement BGCV pool, not hard-coded.
    disputed_donors = set()
    for pool in declared_disagreement_pools:
        if pool.startswith("BGCV"):
            disputed_donors |= (bgcv_sdrf.get(pool, set()) ^ bgcv_local.get(pool, set()))
    bgcv_local_excluding_disputed = exclude_donors_from_pools(bgcv_local, disputed_donors)
    option3_diff = diff_component_mappings(
        bgcv_local, bgcv_local_excluding_disputed,
        "local_bcr_gex_matching", "local_excluding_disputed_donors",
    )
    with open(args.out_dir / "component_diff_excluding_disputed_donors.csv", "w", newline="") as f:
        fieldnames = list(option3_diff["donor_rows"][0].keys()) if option3_diff["donor_rows"] else []
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(option3_diff["donor_rows"])
    (args.out_dir / "component_diff_excluding_disputed_donors_summary.json").write_text(json.dumps({
        "disputed_donors_excluded": sorted(disputed_donors),
        **{k: v for k, v in option3_diff.items() if k != "donor_rows"},
    }, indent=2))
    affected_cells = {
        pool: {
            "gex_cells": gex_pool_cell_counts.get(pool, 0),
            "gex_donors": sorted(gex_pools.get(pool, set())),
        }
        for pool in declared_disagreement_pools
    }
    summary = {
        "n_pools_total": len(conflict_rows),
        "n_pools_in_conflict": len(conflicting_pools),
        "conflicting_pools": conflicting_pools,
        "note": (
            "'conflicting_pools' includes BOTH real declared_membership_disagreement "
            "pools AND insufficient_declared_reference pools (a label only observed "
            "sources use, like the Sanger S-sample names, which is a namespace gap, "
            "not a disagreement). See declared_membership_disagreement_pools for the "
            "real conflicts only."
        ),
        "declared_membership_disagreement_pools": declared_disagreement_pools,
        "insufficient_declared_reference_pools": insufficient_reference_pools,
        "conflict_type_counts": {
            ct: sum(1 for r in conflict_rows if r["conflict_type"] == ct)
            for ct in sorted({r["conflict_type"] for r in conflict_rows})
        },
        "affected_cells_in_conflicting_pools": affected_cells,
        "n_bcr_per_cell_conflicts": len(bcr_cell_conflicts),
        "n_tcr_per_cell_conflicts": len(tcr_cell_conflicts),
        "n_donor_level_mismatches": len(donor_level_mismatches),
        "n_sanger_samples_observed": len(sanger_membership),
        "n_cambridge_pools_observed_in_tcr": len(tcr_cambridge_pools),
        "n_cambridge_pools_observed_in_gex": len(gex_pools),
        "gex_coverage": gex_coverage,
        "n_gex_per_cell_conflicts": len(gex_cell_id_conflicts),
        "disputed_donors_excluded_from_grouping_test": {
            "disputed_donors": sorted(disputed_donors),
            "n_donors_with_changed_component": option3_diff["n_donors_with_changed_component"],
            "note": (
                "n_donors_with_changed_component here always counts the disputed donors "
                "themselves (removing them from every pool trivially changes their own "
                "component to a singleton); what matters is whether any OTHER donor_id is "
                "also in component_diff_excluding_disputed_donors.csv with "
                "component_changed=True -- see that file."
            ),
        },
        "cell_level_crosswalk": cell_join,
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
