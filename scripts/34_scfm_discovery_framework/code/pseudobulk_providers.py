"""Raw-count-sum pseudobulk aggregation for the scFM discovery framework.

Per `review:20261013_222450` B: the first real increment toward
`ExpressionDERankingProvider`'s intended backend (donor/animal-aware
pseudobulk aggregation + a shrunken LFC-style estimator, see
`gene_ranking_providers.py`'s module docstring) -- THIS module implements
ONLY the aggregation step: sample x (caller-selected) cell-type/state raw
count SUMS, with explicit sample/donor/condition/visit identity bound and
validated. It does NOT compute any DE statistic, p-value, log-fold-change,
or `GeneRanking` -- `PseudobulkAggregationProvider.aggregate` never returns
anything that could be mistaken for a DE result, and `ExpressionDERankingProvider`
is untouched by this module. A `status='complete'` result here means
"aggregation complete", nothing about differential expression.

Per Bioconductor OSCA's own convention (Creating pseudo-bulk samples,
https://bioconductor.org/books/release/OSCA.multisample/multi-sample-comparisons.html#creating-pseudo-bulk-samples):
pseudobulk aggregation sums raw counts within (sample, cell label) -- never
within (donor, cell label), since one donor may contribute more than one
real sample (different visits/timepoints), and merging those is a decision
this module refuses to make implicitly. The aggregation unit here is
always `dataset_config["sample_id_column"]`'s value, the dataset's own
declared sample identity column -- grouping by donor/animal ID instead is
not offered as an option, so repeat visits from the same donor never
silently collapse into one row.

Scope this round: RAW, non-negative, integer-valued counts only
(`dataset_config["counts_semantics"] == "raw_counts"`, enforced, never
coerced via rounding or a changed `counts_semantics` declaration); the
caller must explicitly select which cells go into this aggregation --
there is no "all cells" default, matching this framework's existing
refusal to assume an unstated cell-type/state selection. Metadata that
conflicts across the cells of one sample_id (a different donor, condition,
or visit/state value among them) is refused outright, never resolved by
taking the first/last value seen. A sample_id with zero selected cells is
simply absent from the output -- never fabricated as an observed all-zero
row. Opens the input matrix read-only throughout, same as
`embedding_providers.py`; reuses that module's own gene-axis validation
(F2) and dense counts reader rather than re-implementing them.

Per `review:20261014_100436` (review of `32b22fa`): two gaps.

  - F1: the manifest's `sample_rows`/`gene_ids` carried no binding to the
    npz array's actual row/column CONTENT -- reversing `sample_rows` or
    `gene_ids`, or editing `n_distinct_donors`/a row's `library_size`,
    passed both the writer's and the reader's validation unchanged,
    because the npz file-hash check only proves the npz BYTES are
    unmodified, never that the JSON's LABELS for those bytes are correct.
    Fixed with two independent mechanisms: (1) `identity_sha256`, a hash
    over a canonical snapshot of `gene_ids`/`gene_total_counts`/
    `sample_rows`/the derived counts, recomputed and compared on every
    read -- catches any one of these fields changing without every other
    dependent field being rewritten consistently; (2) `load_pseudobulk_
    counts` now recomputes each row's library size and each gene's total
    count DIRECTLY from the loaded npz array and compares them against the
    manifest's claims -- this binds the JSON labels to the array's actual
    content, not just to the JSON's own internal self-consistency.
    `dataset_identity` also gained `selected_raw_counts_sha256`/
    `selected_cell_metadata_sha256`/`gene_ids_sha256` -- content hashes of
    the EXACT selected-cell raw-count block and per-cell sample/donor/
    condition/visit values actually read this run (never re-reads the
    whole source file to get these; they are hashes of what `aggregate`
    already held in memory). These record what was read; verifying them
    against an independent re-read of the ORIGINAL source file is a
    separate, not-yet-implemented capability -- this round states that
    distinction explicitly rather than implying a replay guarantee it does
    not have.

  - F2: type/null validation ran too late or not at all in three ways: (1)
    `cell_type_or_state_label=True` passed the entry TRUTHINESS check
    (`not True` is `False`), so aggregation ran to completion and the npz
    was written BEFORE the writer's own type check finally rejected it,
    leaving an orphan npz with no manifest; (2) a sample whose selected
    cells all had a missing `condition`/`visit` value (the column was
    explicitly requested, not merely unprovided) silently wrote
    `condition_value: null` with no diagnostic; (3) a sample mixing a
    missing value with a real string value crashed `_sample_metadata_
    conflict_problem`'s `sorted()` call with an uncaught `TypeError`
    (`None` is not orderable against `str`), leaving an empty output
    directory. All three are fixed by validating `cell_type_or_state_
    label`/`cell_ids_subset`/`condition_column`/`visit_or_state_column`'s
    TYPES (not just truthiness) at entry, and by refusing outright --
    returning a structured `failed` record, never raising -- any selected
    cell whose EXPLICITLY-requested condition/visit value is missing or
    non-string, before any grouping or `sorted()` call runs. `_write_
    pseudobulk_result` also now validates the complete record's
    non-npz-dependent fields BEFORE writing the npz file, as a second,
    independent guard against ever leaving an orphan npz behind.

Per `review:20261014_114357` (review of `0a40d4f`): F1's own fix
was itself insufficient against an EQUAL-MARGINALS swap -- a fixture with
two rows whose library sizes AND per-gene totals are identical
([[1, 9], [9, 1]], both summing to 10) can have `sample_rows` reversed, or
`gene_ids` reversed, or a `condition_value` rewritten, while every
recomputed row/column TOTAL (and a correspondingly self-consistently
rehashed `identity_sha256`) still agrees -- a total is a many-to-one
projection of content, never a content identity. Fixed by storing
`sample_ids`/`gene_ids`/a normalized per-row `sample_identity_json`
snapshot (donor/condition/visit/cell_ids) INSIDE the npz itself (plain
numpy string arrays, no pickle), covered by the same whole-file
`aggregation_sha256`; `load_pseudobulk_counts` now compares the manifest's
claimed order/metadata DIRECTLY against this independent, hash-protected
in-npz copy, not only against totals recomputed from `counts`. Separately,
`n_pseudobulk_samples`/`n_distinct_donors` in `_validate_pseudobulk_
result_core` had the same `isinstance(x, int) and ... != ...` shape as the
review's R2 (a non-int value skips the consistency check with no `else`
branch reporting the type problem) -- `read_pseudobulk_result` alone (no
npz access) would have accepted e.g. `n_pseudobulk_samples="two"`
unchanged; both fields are now type-checked unconditionally first.
"""
from __future__ import annotations

import json
import math
import os
import platform
import sys
import time
from dataclasses import dataclass, field, replace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from artifact_schema import (array_content_sha256, find_nonfinite_floats,  # noqa: E402
                             resolve_source_revision, safe_json_identity,
                             sha256_file, sha256_text)
from dataset_schema import validate_dataset_config  # noqa: E402
import cli as _cli  # noqa: E402
import embedding_providers as _emb  # noqa: E402

PSEUDOBULK_RESULT_SCHEMA = "scfm-pseudobulk-result/1"

#: Per `review:20261013_135634`'s convention: this module's own
#: source-identity binding, plus the modules it actually reuses.
SOURCE_MODULES = ("pseudobulk_providers.py", "embedding_providers.py",
                  "dataset_schema.py", "cli.py", "artifact_schema.py")

_VALID_STATUSES = ("complete", "failed")


def _reject_nonfinite_constant(token: str):
    raise ValueError(f"non-finite JSON constant {token!r} is not accepted")


def _dependency_versions() -> dict:
    import h5py
    return {"numpy": np.__version__, "h5py": h5py.__version__,
           "python": platform.python_version()}


def _raw_counts_validity_problem(raw: np.ndarray) -> str | None:
    """Per `review:20261013_222450` B.1: refuses a non-negative-
    OR-non-integer counts block outright -- never rounds or coerces it to
    fit. Factored out so a mutation that disables this check can be
    demonstrated directly (see
    `test_mutation_kill_pseudobulk_raw_counts_validity_guard`).
    """
    if np.any(raw < 0):
        return ("the selected counts block contains negative value(s) -- raw "
               "counts must be non-negative")
    if not np.all(raw == np.round(raw)):
        return ("the selected counts block contains non-integer value(s) -- "
               "this round's aggregation requires verifiable raw integer "
               "counts; it never rounds to fit")
    return None


def _sample_metadata_conflict_problem(sample_id: str, donors_here: set,
                                      donor_col: str, conditions_here: set | None,
                                      condition_column: str | None,
                                      visits_here: set | None,
                                      visit_or_state_column: str | None) -> str | None:
    """Per `review:20261013_222450` B.2: a sample_id whose selected
    cells disagree on donor/condition/visit must be refused outright, never
    resolved by taking the first/last value seen. Factored out so a
    mutation that disables this check can be demonstrated directly (see
    `test_mutation_kill_pseudobulk_sample_metadata_conflict_guard`).
    """
    # Per `review:20261014_100436` F2: `key=str` -- by the time this
    # function runs, `aggregate` has already refused any missing/non-string
    # condition/visit value among the SELECTED cells of a column the caller
    # explicitly requested, so a `None`/`str` mix should not reach here. This
    # `key=str` is still a defense-in-depth measure for a direct caller of
    # this function (e.g. a test) that bypasses that upstream guard --
    # `sorted()` alone raises `TypeError` comparing `None` to `str`.
    if len(donors_here) > 1:
        return (f"sample_id={sample_id!r} has conflicting {donor_col!r} values "
               f"across its selected cells: {sorted(donors_here, key=str)} -- "
               "refusing to guess which is correct")
    if conditions_here is not None and len(conditions_here) > 1:
        return (f"sample_id={sample_id!r} has conflicting {condition_column!r} "
               f"values across its selected cells: {sorted(conditions_here, key=str)}")
    if visits_here is not None and len(visits_here) > 1:
        return (f"sample_id={sample_id!r} has conflicting {visit_or_state_column!r} "
               f"values across its selected cells: {sorted(visits_here, key=str)}")
    return None


def _counts_conservation_problem(total_selected: float, total_aggregated: float) -> str | None:
    """Per `review:20261013_222450` B.4: the sum of ALL selected
    cells' counts must equal the sum of the aggregated output -- a defensive
    check that this module's own grouping/summing never silently drops or
    double-counts a cell. Factored out so a mutation that disables this
    check can be demonstrated directly (see
    `test_mutation_kill_pseudobulk_counts_conservation_guard`).
    """
    if abs(total_selected - total_aggregated) > 1e-6:
        return (f"counts conservation check failed: total selected-cell counts "
               f"{total_selected} != total aggregated counts {total_aggregated} -- "
               "refusing to persist a result that does not account for all "
               "selected counts")
    return None


def _pseudobulk_identity_snapshot(gene_ids: tuple, gene_total_counts: tuple,
                                  sample_rows: tuple, n_pseudobulk_samples: int,
                                  n_distinct_donors: int, total_counts_selected_cells: float,
                                  total_counts_aggregated: float) -> dict:
    """Per `review:20261014_100436` F1: a canonical snapshot of
    every field whose correctness depends on `gene_ids`/`sample_rows`
    staying mutually consistent -- reversing `sample_rows`, reversing
    `gene_ids`, or editing `n_distinct_donors`/a row's `library_size`/
    `condition_value` WITHOUT updating every other dependent field changes
    this snapshot's hash. This catches a SLOPPY or PARTIAL edit (one field
    changed, the hash not recomputed); it does not defend against a fully
    self-consistent forgery, the same limitation `_input_mapping_sha256`
    in `gene_ranking_providers.py` already accepts -- a SHA-256 recompute-
    and-compare is a tamper/corruption detector, not a cryptographic
    signature.
    """
    return {
        "gene_ids": list(gene_ids),
        "gene_total_counts": list(gene_total_counts),
        "sample_rows": [
            {k: row.get(k) for k in ("sample_id", "donor_or_animal_id", "condition_value",
                                     "visit_or_state_value", "n_cells", "library_size", "cell_ids")}
            for row in sample_rows
        ],
        "n_pseudobulk_samples": n_pseudobulk_samples,
        "n_distinct_donors": n_distinct_donors,
        "total_counts_selected_cells": total_counts_selected_cells,
        "total_counts_aggregated": total_counts_aggregated,
    }


def _pseudobulk_identity_sha256(*args, **kwargs) -> str:
    return sha256_text(json.dumps(_pseudobulk_identity_snapshot(*args, **kwargs), sort_keys=True))


def _npz_sample_identity_snapshot(sample_rows) -> list:
    """Per `review:20261014_114357` R1: the per-row identity this
    npz's rows actually carry, in ARRAY ROW ORDER -- stored INSIDE the
    hash-protected npz itself (see `_write_pseudobulk_result`), not only in
    the JSON manifest. `identity_sha256` and `load_pseudobulk_counts`'s
    library-size/gene-total recompute (F1, prior round) both proved
    insufficient against an EQUAL-MARGINALS swap: a fixture with two rows
    whose library sizes and per-gene totals are identical ([[1, 9], [9, 1]])
    can have `sample_rows` reversed, or a condition_value rewritten, with
    every recomputed TOTAL still matching -- totals are a many-to-one
    projection of content, never a content identity. This snapshot is
    compared directly against the npz's OWN stored copy (not merely the
    manifest's self-consistently-rehashed claim) in `load_pseudobulk_
    counts`, so a manifest edit that does not also rewrite the npz bytes is
    caught regardless of whether row/column totals happen to agree.
    """
    return [
        {"sample_id": row.get("sample_id"), "donor_or_animal_id": row.get("donor_or_animal_id"),
         "condition_value": row.get("condition_value"),
         "visit_or_state_value": row.get("visit_or_state_value"),
         "cell_ids": list(row.get("cell_ids") or ())}
        for row in sample_rows
    ]


@dataclass(frozen=True)
class PseudobulkResult:
    status: str  # "complete" | "failed" -- aggregation status ONLY, never a DE status
    reason: str = ""
    aggregation_path: str = ""
    aggregation_sha256: str = ""
    #: What predefined cell type/state `cell_ids_subset` was selected for --
    #: stated explicitly by the caller, never inferred.
    cell_type_or_state_label: str = ""
    gene_ids: tuple = ()
    #: Per `review:20261014_100436` F1: the per-gene sum ACROSS all
    #: pseudobulk samples, same order as `gene_ids` -- `load_pseudobulk_
    #: counts` recomputes this directly from the loaded npz array's columns
    #: and compares it against this field, binding `gene_ids`' order to the
    #: array's actual column content (a hash of the npz FILE alone cannot
    #: catch `gene_ids` being reordered without the array itself changing).
    gene_total_counts: tuple = ()
    #: One dict per pseudobulk row: sample_id, donor_or_animal_id,
    #: condition_value, visit_or_state_value, n_cells, library_size,
    #: cell_ids (the ACTUAL cells summed into this row).
    sample_rows: tuple = ()
    #: Reported SEPARATELY, per `review:20261013_222450` B.2 -- a
    #: pseudobulk SAMPLE count must never be read as an independent-donor
    #: count (one donor may contribute more than one sample/visit).
    n_pseudobulk_samples: int = 0
    n_distinct_donors: int = 0
    #: Conservation check: the sum of ALL selected cells' counts (every
    #: gene, every selected cell) vs. the sum of the aggregated output --
    #: these must agree exactly; see `validate_pseudobulk_result`.
    total_counts_selected_cells: float = 0.0
    total_counts_aggregated: float = 0.0
    #: Per `review:20261014_100436` F1: a hash over a canonical
    #: snapshot of `gene_ids`/`gene_total_counts`/`sample_rows` and the
    #: derived counts -- see `_pseudobulk_identity_snapshot`. Recomputed and
    #: compared in `validate_pseudobulk_result`, so a hand-edited field
    #: (e.g. a reordered `sample_rows`, or `n_distinct_donors` changed
    #: without the donor labels actually changing) is caught on read.
    identity_sha256: str = ""
    dataset_identity: dict = field(default_factory=dict)
    backend: str = "raw_count_sum_pseudobulk"
    backend_version: str = "1"
    params: dict = field(default_factory=dict)
    source_code_revision_or_hash: str = ""
    dependency_versions: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        d = dict(self.__dict__)
        d["gene_ids"] = list(self.gene_ids)
        d["gene_total_counts"] = list(self.gene_total_counts)
        d["sample_rows"] = list(self.sample_rows)
        return d


#: Per `review:20261014_100436` F1: content-identity keys that
#: must be present in a `'complete'` result's `dataset_identity` -- bind
#: the EXACT selected raw-count block, per-cell metadata, and gene axis
#: actually read this run; a config hash alone cannot stand in for a
#: content hash of the expression matrix itself.
_DATASET_IDENTITY_CONTENT_KEYS = (
    "selected_raw_counts_sha256", "selected_cell_metadata_sha256", "gene_ids_sha256")


def _validate_pseudobulk_result_core(payload: dict) -> list[str]:
    """Everything `validate_pseudobulk_result` checks EXCEPT `aggregation_
    path`/`aggregation_sha256` -- those two are only known AFTER the npz
    file is written, so `_write_pseudobulk_result` runs this core check
    BEFORE writing the npz (per `review:20261014_100436` F2: a
    structurally-invalid result must never cause an npz file to be written
    and then orphaned by a validation failure on the JSON manifest).
    """
    if not isinstance(payload, dict):
        return [f"result must be a dict, got {type(payload).__name__}"]
    problems = []
    if payload.get("schema") != PSEUDOBULK_RESULT_SCHEMA:
        problems.append(f"schema {payload.get('schema')!r} != expected {PSEUDOBULK_RESULT_SCHEMA!r}")
    status = payload.get("status")
    if status not in _VALID_STATUSES:
        problems.append(f"status {status!r} not in {_VALID_STATUSES}")
    problems.extend(find_nonfinite_floats(payload))
    if status != "complete" and not payload.get("reason"):
        problems.append(f"status={status!r} requires a non-empty 'reason'")

    if status == "complete":
        for key in ("cell_type_or_state_label", "backend", "backend_version",
                    "source_code_revision_or_hash", "identity_sha256"):
            v = payload.get(key)
            if not isinstance(v, str) or not v:
                problems.append(
                    f"status='complete' requires a non-empty string {key!r}, got {v!r}")

        gene_ids = payload.get("gene_ids")
        gene_ids_ok = isinstance(gene_ids, list) and gene_ids \
            and all(isinstance(g, str) and g for g in gene_ids)
        if not gene_ids_ok:
            problems.append(
                "status='complete' requires a non-empty list of non-empty string gene_ids")
        elif len(set(gene_ids)) != len(gene_ids):
            problems.append("gene_ids contains duplicate value(s)")
            gene_ids_ok = False

        gene_total_counts = payload.get("gene_total_counts")
        gtc_ok = isinstance(gene_total_counts, list) and gene_ids_ok \
            and len(gene_total_counts) == len(gene_ids)
        if not isinstance(gene_total_counts, list) or not gene_total_counts:
            problems.append("status='complete' requires a non-empty 'gene_total_counts' list")
        elif gene_ids_ok and len(gene_total_counts) != len(gene_ids):
            problems.append(
                f"gene_total_counts has {len(gene_total_counts)} entries, expected "
                f"len(gene_ids)={len(gene_ids)}")
        if isinstance(gene_total_counts, list):
            for i, v in enumerate(gene_total_counts):
                if not isinstance(v, (int, float)) or isinstance(v, bool) \
                        or not math.isfinite(v) or v < 0:
                    problems.append(
                        f"gene_total_counts[{i}]={v!r} must be a non-negative finite number")
                    gtc_ok = False

        rows = payload.get("sample_rows")
        rows_ok = False
        all_cell_ids: list = []
        if not isinstance(rows, list) or not rows:
            problems.append("status='complete' requires a non-empty 'sample_rows' list")
        else:
            rows_ok = True
            seen_samples: set = set()
            for i, row in enumerate(rows):
                if not isinstance(row, dict):
                    problems.append(f"sample_rows[{i}] must be a dict")
                    rows_ok = False
                    continue
                for key in ("sample_id", "donor_or_animal_id"):
                    if not isinstance(row.get(key), str) or not row.get(key):
                        problems.append(f"sample_rows[{i}][{key!r}] must be a non-empty string")
                        rows_ok = False
                sid = row.get("sample_id")
                if sid in seen_samples:
                    problems.append(f"sample_rows has duplicate sample_id {sid!r}")
                    rows_ok = False
                elif isinstance(sid, str):
                    seen_samples.add(sid)
                n_cells = row.get("n_cells")
                if not isinstance(n_cells, int) or isinstance(n_cells, bool) or n_cells < 1:
                    problems.append(f"sample_rows[{i}]['n_cells'] must be a positive int")
                    rows_ok = False
                cell_ids = row.get("cell_ids")
                if not isinstance(cell_ids, list) or not cell_ids \
                        or not all(isinstance(c, str) and c for c in cell_ids):
                    problems.append(
                        f"sample_rows[{i}]['cell_ids'] must be a non-empty list of "
                        "non-empty strings")
                    rows_ok = False
                else:
                    all_cell_ids.extend(cell_ids)
                    if isinstance(n_cells, int) and not isinstance(n_cells, bool) \
                            and len(cell_ids) != n_cells:
                        problems.append(
                            f"sample_rows[{i}]: len(cell_ids)={len(cell_ids)} != "
                            f"n_cells={n_cells}")
                        rows_ok = False
                lib = row.get("library_size")
                if not isinstance(lib, (int, float)) or isinstance(lib, bool) \
                        or not math.isfinite(lib) or lib < 0:
                    problems.append(
                        f"sample_rows[{i}]['library_size'] must be a non-negative "
                        "finite number")
                    rows_ok = False
            # Per `review:20261014_114357` R2: check the TYPE
            # first, unconditionally -- a non-int `n_pseudobulk_samples`
            # (e.g. `"two"`/`None`) must never silently SKIP the
            # consistency check below (an `isinstance` guard with no `else`
            # branch is a type check that only fires for already-valid
            # types, never for the illegal ones it exists to catch).
            n_declared = payload.get("n_pseudobulk_samples")
            if not isinstance(n_declared, int) or isinstance(n_declared, bool):
                problems.append(f"n_pseudobulk_samples must be an int, got {n_declared!r}")
            elif n_declared != len(rows):
                problems.append(
                    f"n_pseudobulk_samples={n_declared} does not match "
                    f"len(sample_rows)={len(rows)}")
            # Per `review:20261014_100436` F1's "核对细胞分组唯一性":
            # a cell counted into two different pseudobulk rows would mean
            # the same cell's counts were summed twice into different
            # samples -- never valid for a partition of selected cells.
            if len(set(all_cell_ids)) != len(all_cell_ids):
                problems.append(
                    "the same cell_id appears in more than one sample_rows entry -- "
                    "cell grouping must be a partition, not an overlap")

            donors_declared = {row.get("donor_or_animal_id") for row in rows
                               if isinstance(row, dict)}
            # Per `review:20261014_114357` R2: same unconditional-
            # type-check-first pattern as `n_pseudobulk_samples` above --
            # `load_pseudobulk_counts` already rejects a non-int
            # `n_distinct_donors` independently, but `read_pseudobulk_
            # result` alone (no npz access) must reject it too, so the two
            # entry points do not disagree about the same field.
            n_distinct_donors = payload.get("n_distinct_donors")
            if not isinstance(n_distinct_donors, int) or isinstance(n_distinct_donors, bool):
                problems.append(f"n_distinct_donors must be an int, got {n_distinct_donors!r}")
            elif n_distinct_donors != len(donors_declared):
                problems.append(
                    f"n_distinct_donors={n_distinct_donors} does not match the "
                    f"{len(donors_declared)} distinct donor_or_animal_id value(s) actually "
                    "present in sample_rows")

        total_sel = payload.get("total_counts_selected_cells")
        total_agg = payload.get("total_counts_aggregated")
        both_numeric = (isinstance(total_sel, (int, float)) and not isinstance(total_sel, bool)
                       and isinstance(total_agg, (int, float)) and not isinstance(total_agg, bool))
        if not both_numeric:
            problems.append(
                "status='complete' requires numeric 'total_counts_selected_cells' and "
                "'total_counts_aggregated'")
        elif abs(total_sel - total_agg) > 1e-6:
            problems.append(
                f"total_counts_aggregated={total_agg} does not match "
                f"total_counts_selected_cells={total_sel} -- counts conservation violated")

        ds = payload.get("dataset_identity")
        if not isinstance(ds, dict) or not ds:
            problems.append("status='complete' requires a non-empty 'dataset_identity' dict")
        else:
            for key in _DATASET_IDENTITY_CONTENT_KEYS:
                v = ds.get(key)
                if not isinstance(v, str) or not v:
                    problems.append(
                        f"dataset_identity[{key!r}] must be a non-empty string -- a "
                        "config hash alone is not a content hash of the actual "
                        "selected raw-count block/cell metadata/gene axis read")

        dep = payload.get("dependency_versions")
        if not isinstance(dep, dict) or not dep:
            problems.append("status='complete' requires a non-empty 'dependency_versions' dict")
        elif not all(isinstance(k, str) and isinstance(v, str) for k, v in dep.items()):
            problems.append(f"dependency_versions keys/values must all be strings, got {dep!r}")

        # Per `review:20261014_100436` F1: recompute `identity_sha256`
        # from the payload's OWN fields -- any single field above changing
        # without this hash being correspondingly recomputed is caught here,
        # independent of whether the npz file's bytes changed at all.
        if gene_ids_ok and gtc_ok and rows_ok \
                and isinstance(payload.get("n_pseudobulk_samples"), int) \
                and isinstance(payload.get("n_distinct_donors"), int) \
                and both_numeric:
            try:
                expected_identity = _pseudobulk_identity_sha256(
                    gene_ids, gene_total_counts, rows, payload["n_pseudobulk_samples"],
                    payload["n_distinct_donors"], total_sel, total_agg)
            except (TypeError, KeyError):
                problems.append("identity_sha256 could not be recomputed from the payload")
            else:
                if payload.get("identity_sha256") != expected_identity:
                    problems.append(
                        "identity_sha256 does not match a recomputed hash of "
                        "gene_ids/gene_total_counts/sample_rows/the derived counts "
                        "-- tampered or stale")
    return problems


def validate_pseudobulk_result(payload: dict) -> list[str]:
    """Shared by the writer and `read_pseudobulk_result` -- a tampered or
    merely incomplete on-disk file is rejected at BOTH boundaries, matching
    `embedding_providers.validate_embedding_result`'s rigor. Also re-checks
    the counts-conservation identity (`total_counts_aggregated` ==
    `total_counts_selected_cells`), so a hand-edited file cannot claim a
    conservation that was never actually true. See `_validate_pseudobulk_
    result_core` for every check except `aggregation_path`/`aggregation_
    sha256`, which this adds on top (those two are only meaningful once an
    npz file actually exists, so a PRE-write check cannot include them --
    see `_write_pseudobulk_result`).
    """
    problems = _validate_pseudobulk_result_core(payload)
    if not isinstance(payload, dict):
        return problems
    if payload.get("status") == "complete":
        for key in ("aggregation_path", "aggregation_sha256"):
            v = payload.get(key)
            if not isinstance(v, str) or not v:
                problems.append(
                    f"status='complete' requires a non-empty string {key!r}, got {v!r}")
    return problems


def _write_pseudobulk_result(out_dir: str, result: PseudobulkResult, counts=None) -> PseudobulkResult:
    """Per `review:20261014_100436` F2: runs `_validate_pseudobulk_
    result_core` -- every check EXCEPT `aggregation_path`/`aggregation_
    sha256`, which cannot be known until after the npz is written -- BEFORE
    writing any npz bytes. A structurally-invalid `PseudobulkResult` (e.g.
    a non-string `cell_type_or_state_label`) is refused at this pre-write
    point, never leaving an orphan npz with no corresponding manifest.
    """
    pre_payload = result.to_dict()
    pre_payload["schema"] = PSEUDOBULK_RESULT_SCHEMA
    pre_payload["generated_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    pre_problems = _validate_pseudobulk_result_core(pre_payload)
    if pre_problems:
        raise ValueError(f"refusing to write an invalid pseudobulk result: {pre_problems}")

    if result.status == "complete":
        npz_path = os.path.join(out_dir, "pseudobulk_counts.npz")
        tmp = npz_path + ".inprogress"
        # Per `review:20261014_114357` R1: `sample_ids`/`gene_ids`/
        # `sample_identity_json` are saved INSIDE the npz, in the same row
        # order as `counts` -- plain numpy string arrays, never pickled.
        # `aggregation_sha256` (a whole-file hash) already covers these
        # bytes, so `load_pseudobulk_counts` can compare the manifest's
        # claimed labels against this independent, hash-protected copy
        # directly, rather than only against row/column TOTALS recomputed
        # from `counts` (which an equal-marginals swap cannot distinguish).
        sample_ids_arr = np.array([row.get("sample_id") for row in result.sample_rows])
        gene_ids_arr = np.array(list(result.gene_ids))
        identity_json = json.dumps(
            _npz_sample_identity_snapshot(result.sample_rows), sort_keys=True)
        with open(tmp, "wb") as fh:
            np.savez(fh, counts=np.ascontiguousarray(counts, dtype=np.float64),
                    sample_ids=sample_ids_arr, gene_ids=gene_ids_arr,
                    sample_identity_json=np.array([identity_json]))
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, npz_path)
        result = replace(result, aggregation_path=npz_path, aggregation_sha256=sha256_file(npz_path))

    payload = result.to_dict()
    payload["schema"] = PSEUDOBULK_RESULT_SCHEMA
    payload["generated_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    problems = validate_pseudobulk_result(payload)
    if problems:
        raise ValueError(f"refusing to write an invalid pseudobulk result: {problems}")

    result_path = os.path.join(out_dir, "PSEUDOBULK_RESULT.json")
    tmp2 = result_path + ".inprogress"
    with open(tmp2, "w") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True, allow_nan=False)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp2, result_path)
    return result


def read_pseudobulk_result(path: str) -> dict:
    if not os.path.exists(path):
        raise SystemExit(f"pseudobulk result not found: {path}")
    with open(path) as fh:
        try:
            payload = json.load(fh, parse_constant=_reject_nonfinite_constant)
        except ValueError as e:
            raise SystemExit(f"pseudobulk result at {path} failed to parse: {e}")
    problems = validate_pseudobulk_result(payload)
    if problems:
        raise SystemExit(f"pseudobulk result at {path} failed validation: {problems}")
    return payload


def _npz_identity_binding_problems(npz_sample_ids: list, npz_gene_ids: list,
                                   npz_identity_json: str, sample_rows: list,
                                   gene_ids: list) -> list[str]:
    """Per `review:20261014_114357` R1: compares the manifest's
    claimed `sample_id` order, `gene_id` order, and per-row donor/condition/
    visit/cell_id identity DIRECTLY against the npz's OWN stored arrays --
    the direct manifest<->npz binding mechanism, independent of (and in
    addition to) any row/column TOTAL recompute, which an equal-marginals
    swap cannot distinguish at all. Factored out so a mutation that
    disables ONLY this comparison can be demonstrated directly (see
    `test_mutation_kill_load_pseudobulk_counts_npz_identity_binding`).
    """
    problems = []
    manifest_sample_ids = [row.get("sample_id") for row in sample_rows]
    if npz_sample_ids != manifest_sample_ids:
        problems.append(
            f"the npz's own stored sample_ids {npz_sample_ids!r} do not match the "
            f"manifest's sample_rows order {manifest_sample_ids!r} -- the manifest's "
            "row labels are not the ones actually bound to this npz's rows (an "
            "equal-row-total swap would not be caught by a totals-only check)")
    if npz_gene_ids != list(gene_ids):
        problems.append(
            f"the npz's own stored gene_ids {npz_gene_ids!r} do not match the "
            f"manifest's gene_ids {list(gene_ids)!r} -- the manifest's gene axis is "
            "not the one actually bound to this npz's columns (an equal-column-"
            "total swap would not be caught by a totals-only check)")
    expected_identity_json = json.dumps(_npz_sample_identity_snapshot(sample_rows), sort_keys=True)
    if npz_identity_json != expected_identity_json:
        problems.append(
            "the npz's own stored sample identity/metadata snapshot does not match "
            "the manifest's sample_rows (donor_or_animal_id/condition_value/"
            "visit_or_state_value/cell_ids) -- the manifest's per-row metadata is "
            "not the metadata actually bound to this npz's rows")
    return problems


def load_pseudobulk_counts(manifest: dict) -> np.ndarray:
    """Loads the npz artifact a 'complete' manifest references, re-checks
    its file hash against `aggregation_sha256`, and checks the array's
    shape against the manifest's own `sample_rows`/`gene_ids` -- the
    read-back contract this round's review explicitly asked for (mirrors
    `embedding_providers.load_embedding_matrices`).

    Per `review:20261014_100436` F1: `aggregation_sha256` proves
    the npz BYTES are unmodified, never that `sample_rows`/`gene_ids`'
    LABELS are correctly assigned to those bytes -- reversing `sample_rows`
    or `gene_ids` changes no byte of the npz file at all. This function
    therefore RECOMPUTES each row's library size and each gene's total
    count directly from the loaded array and compares them against the
    manifest's own claims, binding the JSON labels to the array's actual
    content rather than only to the JSON's internal self-consistency
    (`validate_pseudobulk_result`'s `identity_sha256` recompute already
    catches most single-field edits; this catches a swap that happens to
    leave `identity_sha256` stale too, and is the direct npz<->JSON
    binding the review asked for).

    Per `review:20261014_114357` R1: a row/column-total recompute
    is a MANY-TO-ONE projection of content -- a fixture with equal row and
    column totals ([[1, 9], [9, 1]]) can have `sample_rows` reversed, or a
    `condition_value` rewritten, with every recomputed total still agreeing.
    This function therefore ALSO compares the manifest's claimed
    `sample_id`/`gene_id` order and per-row donor/condition/visit/cell_id
    identity DIRECTLY against the npz's OWN stored `sample_ids`/`gene_ids`/
    `sample_identity_json` arrays (see `_write_pseudobulk_result`) -- an
    independent, hash-protected copy inside the same file, not merely the
    manifest's own (possibly self-consistently re-hashed) claim about
    itself.
    """
    if manifest.get("status") != "complete":
        raise SystemExit(
            f"load_pseudobulk_counts requires status='complete', got "
            f"{manifest.get('status')!r}")
    path = manifest["aggregation_path"]
    actual_sha = sha256_file(path)
    if actual_sha != manifest["aggregation_sha256"]:
        raise SystemExit(
            f"{path} has file hash {actual_sha!r}, which does not match the "
            f"manifest's recorded aggregation_sha256 {manifest['aggregation_sha256']!r}")
    with np.load(path, allow_pickle=False) as npz:
        counts = npz["counts"]
        missing_keys = [k for k in ("sample_ids", "gene_ids", "sample_identity_json")
                        if k not in npz.files]
        if missing_keys:
            raise SystemExit(
                f"{path} is missing {missing_keys} -- a 'complete' npz must carry "
                "its own sample/gene identity arrays alongside 'counts', bound by "
                "the same file hash, for direct manifest<->npz identity binding")
        npz_sample_ids = [str(x) for x in npz["sample_ids"]]
        npz_gene_ids = [str(x) for x in npz["gene_ids"]]
        npz_identity_json = str(npz["sample_identity_json"][0])
    sample_rows = manifest["sample_rows"]
    gene_ids = manifest["gene_ids"]
    n_rows = len(sample_rows)
    n_genes = len(gene_ids)
    if counts.shape != (n_rows, n_genes):
        raise SystemExit(
            f"counts array shape {counts.shape} != expected ({n_rows}, {n_genes})")
    if not np.all(np.isfinite(counts)):
        raise SystemExit("counts array contains non-finite value(s)")

    identity_problems = _npz_identity_binding_problems(
        npz_sample_ids, npz_gene_ids, npz_identity_json, sample_rows, gene_ids)
    if identity_problems:
        raise SystemExit("; ".join(identity_problems))

    recomputed_row_sums = counts.sum(axis=1)
    for i, row in enumerate(sample_rows):
        claimed = row.get("library_size")
        if not isinstance(claimed, (int, float)) or isinstance(claimed, bool) \
                or abs(float(recomputed_row_sums[i]) - claimed) > 1e-6:
            raise SystemExit(
                f"sample_rows[{i}] (sample_id={row.get('sample_id')!r}) claims "
                f"library_size={claimed!r}, but the loaded counts array's row {i} "
                f"actually sums to {float(recomputed_row_sums[i])!r} -- the manifest's "
                "row labels do not match the npz array's row content")

    gene_total_counts = manifest.get("gene_total_counts")
    if not isinstance(gene_total_counts, list) or len(gene_total_counts) != n_genes:
        raise SystemExit(
            f"manifest 'gene_total_counts' must be a list of length {n_genes}, "
            f"got {gene_total_counts!r}")
    recomputed_gene_sums = counts.sum(axis=0)
    for j, claimed in enumerate(gene_total_counts):
        if not isinstance(claimed, (int, float)) or isinstance(claimed, bool) \
                or abs(float(recomputed_gene_sums[j]) - claimed) > 1e-6:
            raise SystemExit(
                f"gene_total_counts[{j}] (gene_id={gene_ids[j]!r}) claims "
                f"{claimed!r}, but the loaded counts array's column {j} actually "
                f"sums to {float(recomputed_gene_sums[j])!r} -- the manifest's gene "
                "axis does not match the npz array's column content")

    recomputed_total = float(counts.sum())
    claimed_total = manifest.get("total_counts_aggregated")
    if not isinstance(claimed_total, (int, float)) or isinstance(claimed_total, bool) \
            or abs(recomputed_total - claimed_total) > 1e-6:
        raise SystemExit(
            f"manifest claims total_counts_aggregated={claimed_total!r}, but the "
            f"loaded counts array actually sums to {recomputed_total!r}")

    donors = {row.get("donor_or_animal_id") for row in sample_rows}
    claimed_donors = manifest.get("n_distinct_donors")
    if not isinstance(claimed_donors, int) or isinstance(claimed_donors, bool) \
            or claimed_donors != len(donors):
        raise SystemExit(
            f"manifest claims n_distinct_donors={claimed_donors!r}, but "
            f"sample_rows actually names {len(donors)} distinct donor_or_animal_id "
            "value(s)")
    return counts


class PseudobulkAggregationProvider:
    """Real, numpy-backed raw-count-sum pseudobulk aggregation. See module
    docstring for the scope this round explicitly supports and explicitly
    refuses. Reuses `embedding_providers`'s gene-axis validation (F2) and
    dense counts reader -- never re-implements them.
    """

    def aggregate(self, dataset_config: dict, dataset_config_path: str, *,
                 out_dir: str, cell_ids_subset: tuple, cell_type_or_state_label: str,
                 condition_column: str | None = None,
                 visit_or_state_column: str | None = None) -> PseudobulkResult:
        if os.path.exists(out_dir) and os.listdir(out_dir):
            raise SystemExit(f"{out_dir} must be new or empty")
        os.makedirs(out_dir, exist_ok=True)

        source_revision = resolve_source_revision(HERE, SOURCE_MODULES)
        common = dict(
            backend="raw_count_sum_pseudobulk", backend_version="1",
            source_code_revision_or_hash=source_revision,
            dependency_versions=_dependency_versions())

        def finish(result: PseudobulkResult, **arrays) -> PseudobulkResult:
            return _write_pseudobulk_result(out_dir, result, **arrays)

        def fail(reason: str) -> PseudobulkResult:
            return finish(PseudobulkResult(status="failed", reason=reason, **common))

        # Per `review:20261014_100436` F2: `isinstance(..., str)`,
        # not a bare truthiness check -- `cell_type_or_state_label=True` is
        # truthy but not a string, and a truthiness-only check let this
        # aggregation run to completion (writing an npz) before the
        # writer's own type check rejected it, leaving an orphan npz.
        if not isinstance(cell_type_or_state_label, str) or not cell_type_or_state_label:
            return fail(
                f"cell_type_or_state_label must be a non-empty string, got "
                f"{cell_type_or_state_label!r} -- this aggregation is scoped to a "
                "predefined cell type/state, never 'all cells' by default")
        if not cell_ids_subset:
            return fail(
                "cell_ids_subset is empty -- the caller must explicitly select "
                "which cells go into this aggregation; there is no 'all cells' "
                "default")
        if not all(isinstance(c, str) and c for c in cell_ids_subset):
            return fail("cell_ids_subset must contain only non-empty strings")
        if len(set(cell_ids_subset)) != len(cell_ids_subset):
            return fail("cell_ids_subset contains duplicate value(s)")
        if condition_column is not None and (not isinstance(condition_column, str) or not condition_column):
            return fail(
                f"condition_column must be None or a non-empty string, got "
                f"{condition_column!r}")
        if visit_or_state_column is not None \
                and (not isinstance(visit_or_state_column, str) or not visit_or_state_column):
            return fail(
                f"visit_or_state_column must be None or a non-empty string, got "
                f"{visit_or_state_column!r}")

        val = validate_dataset_config(dataset_config)
        if not val.ok:
            return fail(f"dataset_config failed validation: {[str(i) for i in val.issues]}")
        if dataset_config["counts_semantics"] != "raw_counts":
            return fail(
                f"counts_semantics={dataset_config['counts_semantics']!r} -- this "
                "aggregation requires 'raw_counts' explicitly; it never coerces "
                "log/normalized values by rounding, and never relabels "
                "counts_semantics to fit")

        matrix_path = _cli.resolve_matrix_path(dataset_config_path, dataset_config)
        if not matrix_path:
            return fail("dataset_config has no expression_matrix_path -- nothing to read")

        structural_problems = _cli.validate_matrix_structure(matrix_path, dataset_config)
        if structural_problems:
            return fail(f"matrix structure invalid: {structural_problems}")

        counts_slot = dataset_config["counts_slot"]
        cell_id_col = dataset_config["cell_id_column"]
        sample_id_col = dataset_config["sample_id_column"]
        donor_col = dataset_config["subject_id_column"]

        all_cell_ids = _cli.read_obs_column_values(matrix_path, cell_id_col)
        pos = {c: i for i, c in enumerate(all_cell_ids)}
        missing = [c for c in cell_ids_subset if c not in pos]
        if missing:
            return fail(
                f"cell_ids_subset contains {len(missing)} ID(s) not present in "
                f"{cell_id_col!r}: {missing[:5]}")
        row_index = [pos[c] for c in cell_ids_subset]

        try:
            sample_values = _cli.read_obs_column_values(matrix_path, sample_id_col)
            donor_values = _cli.read_obs_column_values(matrix_path, donor_col)
            condition_values = (_cli.read_obs_column_values(matrix_path, condition_column)
                                if condition_column else None)
            visit_values = (_cli.read_obs_column_values(matrix_path, visit_or_state_column)
                            if visit_or_state_column else None)
        except SystemExit as e:
            return fail(str(e))

        var_group_path = _cli.gene_axis_group_for_slot(counts_slot)
        try:
            gene_ids = _emb._read_var_index_values(matrix_path, var_group_path)
        except SystemExit as e:
            return fail(str(e))
        gene_axis_problems = _emb._validate_gene_axis_values(gene_ids, var_group_path)
        if gene_axis_problems:
            return fail("gene axis invalid: " + "; ".join(gene_axis_problems))

        try:
            raw = _emb._read_counts_slot_dense(matrix_path, counts_slot, row_index)
        except SystemExit as e:
            return fail(str(e))

        bad = np.where(~np.isfinite(raw))
        if bad[0].size:
            return fail(
                f"counts_slot={counts_slot!r} contains {bad[0].size} non-finite "
                "value(s) in the selected cells")
        raw_counts_problem = _raw_counts_validity_problem(raw)
        if raw_counts_problem is not None:
            return fail(raw_counts_problem)

        sel_sample = [sample_values[i] for i in row_index]
        sel_donor = [donor_values[i] for i in row_index]
        sel_condition = [condition_values[i] for i in row_index] if condition_values is not None else None
        sel_visit = [visit_values[i] for i in row_index] if visit_values is not None else None

        bad_sample = [i for i, s in enumerate(sel_sample) if not isinstance(s, str) or not s]
        if bad_sample:
            return fail(
                f"{len(bad_sample)} selected cell(s) have an empty/non-string "
                f"{sample_id_col!r} value")
        bad_donor = [i for i, d in enumerate(sel_donor) if not isinstance(d, str) or not d]
        if bad_donor:
            return fail(
                f"{len(bad_donor)} selected cell(s) have an empty/non-string "
                f"{donor_col!r} value")
        # Per `review:20261014_100436` F2: `condition_column=None`
        # ("the caller did not provide this optional column") and
        # `condition_column` explicitly provided but some selected cell's
        # VALUE is missing/non-string are two DIFFERENT situations -- the
        # latter is refused outright this round, never silently written as
        # `condition_value: null` and never silently filled/dropped. This
        # upfront check also removes the possibility of a `None`/`str` mix
        # ever reaching `_sample_metadata_conflict_problem`'s `sorted()`
        # call below, which previously raised an uncaught `TypeError`.
        if sel_condition is not None:
            bad_condition = [i for i, c in enumerate(sel_condition) if not isinstance(c, str) or not c]
            if bad_condition:
                return fail(
                    f"condition_column={condition_column!r} was explicitly provided, "
                    f"but {len(bad_condition)} selected cell(s) have a missing/"
                    "non-string value -- this round refuses rather than silently "
                    "filling or dropping cells; do not pass condition_column if its "
                    "value is not available for every selected cell")
        if sel_visit is not None:
            bad_visit = [i for i, v in enumerate(sel_visit) if not isinstance(v, str) or not v]
            if bad_visit:
                return fail(
                    f"visit_or_state_column={visit_or_state_column!r} was explicitly "
                    f"provided, but {len(bad_visit)} selected cell(s) have a missing/"
                    "non-string value -- this round refuses rather than silently "
                    "filling or dropping cells; do not pass visit_or_state_column if "
                    "its value is not available for every selected cell")

        by_sample: dict = {}
        for local_i, sid in enumerate(sel_sample):
            by_sample.setdefault(sid, []).append(local_i)

        sample_rows = []
        counts_rows = []
        for sid in sorted(by_sample):
            idxs = by_sample[sid]
            donors_here = {sel_donor[i] for i in idxs}
            conditions_here = {sel_condition[i] for i in idxs} if sel_condition is not None else None
            visits_here = {sel_visit[i] for i in idxs} if sel_visit is not None else None
            conflict_problem = _sample_metadata_conflict_problem(
                sid, donors_here, donor_col, conditions_here, condition_column,
                visits_here, visit_or_state_column)
            if conflict_problem is not None:
                return fail(conflict_problem)

            row_cell_ids = [cell_ids_subset[i] for i in idxs]
            row_counts = raw[idxs, :].sum(axis=0)
            counts_rows.append(row_counts)
            sample_rows.append({
                "sample_id": sid,
                "donor_or_animal_id": next(iter(donors_here)),
                "condition_value": next(iter(conditions_here)) if conditions_here else "",
                "visit_or_state_value": next(iter(visits_here)) if visits_here else "",
                "n_cells": len(idxs),
                "library_size": float(row_counts.sum()),
                "cell_ids": row_cell_ids,
            })

        counts_matrix = np.stack(counts_rows).astype(np.float64)
        total_counts_selected_cells = float(raw.sum())
        total_counts_aggregated = float(counts_matrix.sum())
        conservation_problem = _counts_conservation_problem(
            total_counts_selected_cells, total_counts_aggregated)
        if conservation_problem is not None:
            return fail(conservation_problem)

        n_distinct_donors = len({row["donor_or_animal_id"] for row in sample_rows})
        gene_total_counts = tuple(float(x) for x in counts_matrix.sum(axis=0))

        # Per `review:20261014_100436` F1: content hashes of the
        # EXACT selected raw-count block and per-cell metadata actually
        # read this run -- never re-reads the whole source file to get
        # these (they are hashes of what is ALREADY in memory). This binds
        # `dataset_identity` to the expression matrix's actual content,
        # not only to `dataset_config`'s own hash (a config hash cannot
        # catch the matrix FILE changing while the config stays the same).
        # These record what was read; they do not themselves prove that an
        # INDEPENDENT re-read of the original source would reproduce the
        # same content -- that would require a separate replay step this
        # round does not implement.
        selected_cell_metadata = [
            {"cell_id": cell_ids_subset[i], "sample_id": sel_sample[i],
             "donor_id": sel_donor[i],
             "condition_value": sel_condition[i] if sel_condition is not None else None,
             "visit_or_state_value": sel_visit[i] if sel_visit is not None else None}
            for i in range(len(cell_ids_subset))
        ]
        dataset_identity = dict(
            dataset_id=dataset_config["dataset_id"], species=dataset_config["species"],
            gene_id_namespace=dataset_config["gene_id_namespace"], counts_slot=counts_slot,
            counts_semantics=dataset_config["counts_semantics"],
            expression_matrix_path_resolved=matrix_path,
            dataset_config_sha256=sha256_text(json.dumps(dataset_config, sort_keys=True)),
            selected_raw_counts_sha256=array_content_sha256(raw),
            selected_cell_metadata_sha256=safe_json_identity(selected_cell_metadata) or "",
            gene_ids_sha256=safe_json_identity(list(gene_ids)) or "")

        params = dict(
            sample_id_column=sample_id_col, donor_or_animal_id_column=donor_col,
            condition_column=condition_column or "",
            visit_or_state_column=visit_or_state_column or "",
            n_cells_selected=len(cell_ids_subset))

        identity_sha256 = _pseudobulk_identity_sha256(
            tuple(gene_ids), gene_total_counts, tuple(sample_rows), len(sample_rows),
            n_distinct_donors, total_counts_selected_cells, total_counts_aggregated)

        result = PseudobulkResult(
            status="complete", cell_type_or_state_label=cell_type_or_state_label,
            gene_ids=tuple(gene_ids), gene_total_counts=gene_total_counts,
            sample_rows=tuple(sample_rows),
            n_pseudobulk_samples=len(sample_rows), n_distinct_donors=n_distinct_donors,
            total_counts_selected_cells=total_counts_selected_cells,
            total_counts_aggregated=total_counts_aggregated,
            identity_sha256=identity_sha256,
            dataset_identity=dataset_identity, params=params, **common)

        return finish(result, counts=counts_matrix)
