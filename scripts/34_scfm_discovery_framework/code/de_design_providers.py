"""Pseudobulk -> DE input/design-matrix PREPARATION ONLY for the scFM
discovery framework. Per `review:20261014_151747` step 2 (the "B"
increment, gated on the A increment -- `review:20261014_114357`
R1-R3 plus `review:20261014_151747` P2 -- having been fixed and
tested first).

This module performs NO statistical fitting. It consumes a VALIDATED
pseudobulk aggregation (`pseudobulk_providers.read_pseudobulk_result` +
`pseudobulk_providers.load_pseudobulk_counts` -- never a raw/unvalidated
manifest or npz), an EXPLICIT `gene_ranking_providers.ComparisonSpec` and
`gene_ranking_providers.SampleStructure` (the SAME interfaces
`ExpressionDERankingProvider` already uses -- never a parallel/duplicate
pair of types), and produces a `PreparedDEDesign`: an ordered sample/gene
identity, a real numeric design matrix, a named target contrast, and the
design's rank/residual-degrees-of-freedom -- all either genuinely usable
(`status='prepared'`) or explicitly refused (`status='failed'`), never a
partial/degenerate matrix silently accepted.

`status` here is `'prepared'`/`'failed'` ONLY -- deliberately NOT the
`'complete'`/`'not_run'`/`'failed'` triple `ExpressionDEResult` uses, so a
caller can never confuse "a design matrix was prepared" with "a DE
ranking was computed" (`ExpressionDERankingProvider.run` still always
returns `status in ('failed', 'not_run')` and never a real ranking -- see
`gene_ranking_providers.py`'s module docstring). This module's OWN output
is explicitly "DE input/design prepared" and nothing more; no LFC, no
p-value, no `GeneRanking`.

First version supports exactly two synthetic designs, matching
`gene_ranking_providers._VALID_DESIGN_TYPES` (reused directly, never
redefined):

  - `"independent_donors_unpaired"`: every donor contributes exactly ONE
    sample, to exactly ONE of the two compared groups. Design matrix:
    intercept + a single `condition[comparison_level]` dummy (treatment
    coding, `reference_level=0`).
  - `"paired_or_repeated_within_donor"`: every donor contributes exactly
    ONE sample to EACH of the two compared groups (fully paired; this
    round supports no partial pairing -- see `gene_ranking_providers.py`'s
    R2). Design matrix: intercept + donor fixed effects (one donor's
    reference level dropped) + the same condition dummy.

Explicitly REFUSED, never silently pooled or dropped: declared
`sample_structure.covariate_columns` (no covariate VALUES are taken this
round); more than one sample for the same donor within the same compared
group (covers extra visits/repeated sampling the two supported designs do
not model); a donor that straddles both groups under
`"independent_donors_unpaired"`; a donor missing either group under
`"paired_or_repeated_within_donor"`; a missing comparison level entirely;
a design matrix that is not full rank; a design with fewer than 1 residual
degree of freedom (full rank is NOT the same claim as adequate power --
see `_prepare_de_design_core`); a target contrast that is not estimable
from the design matrix (checked directly via its row-space projection,
never assumed from rank alone).

Counts/metadata alignment is by `sample_id` LOOKUP, never by position:
`DEDesignPreparationProvider.prepare` builds an explicit `sample_id ->
pseudobulk-row-index` map and selects/reorders the counts submatrix
through it, so a caller's `sample_rows` order and a design's row order
disagreeing is never silently masked by them happening to already match.

Persistence mirrors `pseudobulk_providers.py`'s own R1 fix (direct
in-npz identity binding, never total-only): the selected/reordered counts
submatrix is saved to its own npz, alongside its own copies of
`sample_ids`/`gene_ids`, covered by BOTH a whole-file hash
(`counts_file_sha256`) and a content-only hash
(`counts_content_sha256`, `artifact_schema.array_content_sha256`) -- so a
label swap that happens to leave file bytes self-consistent (the exact
equal-marginals class of bug `review:20261014_151747`'s Material
Passport attributes to this round's prior A increment) is caught the same
way here, proactively, rather than waiting for a follow-up review to name
it. See `load_prepared_de_counts`.

Per `review:20261014_203915` F1/F2 (a follow-up review of this
module's first version): two gaps were found and fixed here.

F1 -- `prepare()` previously took an independently-supplied `counts`
ndarray alongside the manifest `dict`, checking only their `status` and
`shape` agreement; a caller could pass counts whose rows/sign had been
altered (same shape, wrong content or wrong row order) and it was
accepted, still recorded under the real upstream `aggregation_sha256`.
`prepare()` now takes `pseudobulk_result_path` and calls
`pseudobulk_providers.read_pseudobulk_result` +
`pseudobulk_providers.load_pseudobulk_counts` ITSELF -- the counts array
a caller could tamper with independently of the manifest no longer
exists as an input at all; the only inputs are a path and the hash/
identity-bound loader's own output. `PreparedDEDesign.pseudobulk_identity`
also now records the manifest's OWN path and whole-file hash
(`pseudobulk_manifest_path`/`pseudobulk_manifest_sha256`), and
`PreparedDEDesign.dataset_identity` carries the upstream `species`/
`gene_id_namespace`/`counts_semantics` forward -- a downstream consumer
(GeneRanking/GSEA) is never left to guess the species from a file name.

F2 -- declared `comparison.group_column`/`sample_structure.group_column`/
`sample_structure.sample_id_column`/`sample_structure.
donor_or_animal_id_column` were never checked against the pseudobulk
aggregation's ACTUAL column names (`pseudobulk_manifest["params"]`); and
at read time, `rank`/`residual_degrees_of_freedom`/contrast estimability
were recomputed from `design_matrix`/`contrast_vector` directly, but
nothing rebuilt the EXPECTED `design_matrix`/`column_names`/`formula`/
`contrast_vector`/`contrast_direction` from the persisted
`sample_ids`/`donor_ids`/`group_values`/`comparison`/`sample_structure`
and compared them -- and the "estimability" check itself is VACUOUS
whenever the design is full rank (`c @ pinv(X) @ X == c` identically for
ANY `c` when `X` has full column rank), so a reversed or even a zero
contrast vector passed. `prepare()` now rejects a column-name mismatch
against the upstream `params` before building anything, and
`_validate_prepared_de_design_core` now reconstructs
`ComparisonSpec`/`SampleStructure` from the persisted payload and reruns
`_prepare_de_design_core` on the persisted identity triples
(`_expected_design_fields_from_identity`), rejecting ANY mismatch between
the recomputed and the stored semantic fields -- the same
recompute-and-compare principle this framework already applies to
hashes, applied here to the design's own biological meaning.
"""
from __future__ import annotations

import dataclasses
import json
import os
import platform
import sys
import time
from dataclasses import dataclass, field, replace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from artifact_schema import (array_content_sha256, find_nonfinite_floats,  # noqa: E402
                             resolve_source_revision, sha256_file)
from gene_ranking_providers import ComparisonSpec, SampleStructure  # noqa: E402
import pseudobulk_providers as _pb  # noqa: E402

SOURCE_MODULES = ("de_design_providers.py", "gene_ranking_providers.py",
                  "pseudobulk_providers.py", "artifact_schema.py")

#: Maps a declared config field to the pseudobulk aggregation's own
#: `params` key it must agree with -- per `review:20261014_203915`
#: F2: a declared column name that does not match what the aggregation
#: ACTUALLY used must be refused before any design is built, never
#: trusted as a caller's unverified label.
_COLUMN_NAME_BINDINGS = (
    ("sample_structure.sample_id_column", "sample_structure", "sample_id_column",
     "sample_id_column"),
    ("sample_structure.donor_or_animal_id_column", "sample_structure",
     "donor_or_animal_id_column", "donor_or_animal_id_column"),
    ("sample_structure.group_column", "sample_structure", "group_column", "condition_column"),
    ("comparison.group_column", "comparison", "group_column", "condition_column"),
)

_VALID_PREP_STATUSES = ("prepared", "failed")
PREPARED_DE_DESIGN_SCHEMA = "scfm-prepared-de-design/1"

#: Currently identical to `gene_ranking_providers._VALID_DESIGN_TYPES`, but
#: declared INDEPENDENTLY on purpose: `SampleStructure.validate()` accepting
#: a `design_type` means it is a design `ExpressionDERankingProvider`'s
#: estimability check understands, never a promise that THIS module's
#: design-matrix builder has an implementation for it. If a third design
#: type is ever added there before this module supports it, this tuple
#: must stay unchanged and the check below must still refuse it.
_SUPPORTED_DESIGN_TYPES = ("independent_donors_unpaired", "paired_or_repeated_within_donor")


def _source_revision() -> str:
    return resolve_source_revision(HERE, SOURCE_MODULES)


def _reject_nonfinite_constant(token: str):
    raise ValueError(f"non-finite JSON constant {token!r} is not accepted")


def _prepare_de_design_core(sample_rows: list, comparison: ComparisonSpec,
                            sample_structure: SampleStructure) -> tuple[list[str], dict]:
    """Pure, I/O-free: the actual design-building logic, factored out so
    every refusal path can be tested directly (and mutation-killed) without
    constructing a real pseudobulk artifact on disk. Returns
    `(problems, fields)` -- `fields` is only meaningful when `problems`
    is empty, and never contains partial/best-effort values otherwise.
    """
    problems = comparison.validate() + sample_structure.validate()
    if problems:
        return problems, {}

    if sample_structure.covariate_columns:
        return [
            f"sample_structure.covariate_columns={sample_structure.covariate_columns!r} is "
            "declared, but this round's DE design preparation takes no covariate VALUES "
            "to build into the design matrix -- refuse rather than silently ignore"], {}

    design_type = sample_structure.design_type
    if design_type not in _SUPPORTED_DESIGN_TYPES:
        return [
            f"design_type={design_type!r} not in {_SUPPORTED_DESIGN_TYPES!r} -- this round's "
            "DE design preparation supports only these two synthetic designs"], {}

    ref, cmp = comparison.reference_level, comparison.comparison_level
    selected = [r for r in sample_rows if r.get("condition_value") in (ref, cmp)]
    sids = [r.get("sample_id") for r in selected]
    if len(set(sids)) != len(sids):
        return [f"duplicate sample_id(s) among the comparison's selected sample_rows: {sids!r}"], {}

    by_group = {ref: [], cmp: []}
    for r in selected:
        by_group[r["condition_value"]].append(r)
    missing_groups = [lvl for lvl in (ref, cmp) if not by_group[lvl]]
    if missing_groups:
        return [
            f"comparison level(s) {missing_groups} have zero pseudobulk samples -- BOTH "
            f"{ref!r} and {cmp!r} must be present to prepare this comparison"], {}

    donor_group_counts: dict = {}
    for r in selected:
        key = (r.get("donor_or_animal_id"), r.get("condition_value"))
        donor_group_counts[key] = donor_group_counts.get(key, 0) + 1
    dup = {k: v for k, v in donor_group_counts.items() if v > 1}
    if dup:
        return [
            f"donor/group combination(s) {dup!r} have more than one pseudobulk sample -- "
            "repeated samples (e.g. extra visits) for the same donor within the same "
            "compared group are not supported this round, and are never silently pooled"], {}

    donor_levels: dict = {}
    for r in selected:
        donor_levels.setdefault(r["donor_or_animal_id"], set()).add(r["condition_value"])

    if design_type == "independent_donors_unpaired":
        straddling = {d: sorted(lv) for d, lv in donor_levels.items() if len(lv) > 1}
        if straddling:
            return [
                f"design_type='independent_donors_unpaired' but donor(s) {straddling!r} "
                "contribute samples to BOTH compared groups -- contradicts the declared "
                "design; a cell or repeated sample is never treated as an independent donor"], {}
    else:  # paired_or_repeated_within_donor
        unpaired = {d: sorted(lv) for d, lv in donor_levels.items() if lv != {ref, cmp}}
        if unpaired:
            return [
                f"design_type='paired_or_repeated_within_donor' but donor(s) {unpaired!r} do "
                f"not carry BOTH {ref!r} and {cmp!r} -- this round supports only the "
                "FULLY-paired case, never partial pairing"], {}

    # Deterministic row order: donor id, then group (reference before
    # comparison), then sample id -- NEVER the caller's `sample_rows` input
    # order, which carries no ordering guarantee of its own.
    selected_sorted = sorted(
        selected,
        key=lambda r: (r["donor_or_animal_id"], 0 if r["condition_value"] == ref else 1, r["sample_id"]))

    sample_ids = [r["sample_id"] for r in selected_sorted]
    donor_ids = [r["donor_or_animal_id"] for r in selected_sorted]
    group_values = [r["condition_value"] for r in selected_sorted]
    n_samples = len(selected_sorted)
    distinct_donors = sorted(set(donor_ids))
    n_distinct_donors = len(distinct_donors)

    condition_dummy = [1.0 if g == cmp else 0.0 for g in group_values]
    if design_type == "independent_donors_unpaired":
        column_names = ["intercept", f"condition[{cmp}]"]
        design_matrix = [[1.0, condition_dummy[i]] for i in range(n_samples)]
        formula = f"~ 1 + condition[{cmp}]"
        encoding_convention = (
            f"treatment coding: condition {ref!r}=0 (reference_level), {cmp!r}=1 "
            "(comparison_level); no donor term (independent_donors_unpaired)")
    else:
        donor_ref = distinct_donors[0]
        donor_cols = distinct_donors[1:]
        column_names = ["intercept"] + [f"donor[{d}]" for d in donor_cols] + [f"condition[{cmp}]"]
        design_matrix = []
        for i in range(n_samples):
            row = ([1.0] + [1.0 if donor_ids[i] == d else 0.0 for d in donor_cols]
                  + [condition_dummy[i]])
            design_matrix.append(row)
        formula = f"~ 1 + donor + condition[{cmp}]"
        encoding_convention = (
            f"treatment coding: condition {ref!r}=0 (reference_level), {cmp!r}=1 "
            f"(comparison_level); donor fixed effects, reference donor={donor_ref!r} "
            "(dropped column, absorbed into intercept)")

    X = np.array(design_matrix, dtype=float)
    rank = int(np.linalg.matrix_rank(X))
    residual_df = n_samples - rank
    if rank < len(column_names):
        return [
            f"design matrix is not full rank: recomputed rank={rank}, expected "
            f"{len(column_names)} columns {column_names!r}"], {}
    if residual_df < 1:
        return [
            f"residual_degrees_of_freedom={residual_df} (<1) with n_samples={n_samples}, "
            f"rank={rank} -- a full-rank design with no residual degrees of freedom cannot "
            "support real dispersion estimation; full rank is never itself treated as "
            "adequate power"], {}

    contrast_vector = [1.0 if name == f"condition[{cmp}]" else 0.0 for name in column_names]
    c = np.array(contrast_vector, dtype=float)
    pinv = np.linalg.pinv(X)
    projected = c @ pinv @ X
    if not np.allclose(projected, c, atol=1e-8):
        return [
            f"target contrast {contrast_vector!r} (column {column_names.index(f'condition[{cmp}]')!r}) "
            "is not estimable from this design matrix (outside its row space)"], {}

    fields = dict(
        sample_ids=tuple(sample_ids), donor_ids=tuple(donor_ids), group_values=tuple(group_values),
        design_matrix=tuple(tuple(row) for row in design_matrix), column_names=tuple(column_names),
        formula=formula, encoding_convention=encoding_convention,
        contrast_name=f"condition[{cmp}]", contrast_vector=tuple(contrast_vector),
        contrast_direction=f"{cmp!r} minus {ref!r} ({comparison.direction_definition})",
        rank=rank, residual_degrees_of_freedom=residual_df,
        n_samples=n_samples, n_distinct_donors=n_distinct_donors)
    return [], fields


#: Fields `_expected_design_fields_from_identity` compares, verbatim,
#: between what `_prepare_de_design_core` recomputes from the PERSISTED
#: `sample_ids`/`donor_ids`/`group_values`/`comparison`/`sample_structure`
#: and what is actually stored in a `PreparedDEDesign` payload. Listed
#: once so `_validate_prepared_de_design_core` and its tests refer to the
#: same authoritative set -- see `review:20261014_203915` F2.
_SEMANTIC_FIELDS_BOUND_TO_IDENTITY = (
    "sample_ids", "donor_ids", "group_values", "design_matrix", "column_names",
    "formula", "encoding_convention", "contrast_name", "contrast_vector",
    "contrast_direction", "rank", "residual_degrees_of_freedom", "n_samples",
    "n_distinct_donors")


def _expected_design_fields_from_identity(sample_ids: list, donor_ids: list, group_values: list,
                                          comparison: ComparisonSpec,
                                          sample_structure: SampleStructure) -> tuple[list[str], dict]:
    """Rebuilds `sample_rows`-shaped dicts directly from a PERSISTED
    payload's own `sample_ids`/`donor_ids`/`group_values` and reruns
    `_prepare_de_design_core` on them -- the same recompute-and-compare
    principle this framework already applies to hashes/rank/residual-df,
    applied here to the design's biological MEANING: which samples got
    which group/donor role, which column names and contrast direction
    those roles produce, and whether the declared `comparison`/
    `sample_structure` (also reconstructed from the SAME persisted
    payload, never from a separately-trusted copy) are what actually
    produced the stored `design_matrix`/`contrast_vector`. Per
    `review:20261014_203915` F2: a hash-of-design_matrix plus an
    estimability check (which is VACUOUS for any full-rank design) cannot
    catch a reversed/zero contrast, a swapped reference/comparison level,
    or reordered group/donor roles -- rebuilding the whole design from
    the claimed identity and comparing it field-by-field can.
    """
    rows = [{"sample_id": sid, "donor_or_animal_id": did, "condition_value": gv}
           for sid, did, gv in zip(sample_ids, donor_ids, group_values)]
    return _prepare_de_design_core(rows, comparison, sample_structure)


@dataclass(frozen=True)
class PreparedDEDesign:
    status: str  # "prepared" | "failed" -- see module docstring; NEVER "complete"/"not_run"
    reason: str = ""
    evidence_type: str = "de_design_preparation"
    comparison: dict = field(default_factory=dict)
    sample_structure: dict = field(default_factory=dict)
    #: The EXACT upstream pseudobulk artifact this design was built from --
    #: never just a config hash. Carries the manifest's OWN path/whole-file
    #: hash (`pseudobulk_manifest_path`/`pseudobulk_manifest_sha256`) in
    #: addition to the npz-side `aggregation_path`/`aggregation_sha256`/
    #: `identity_sha256`. See `pseudobulk_providers.PseudobulkResult`.
    pseudobulk_identity: dict = field(default_factory=dict)
    #: Carried forward from `pseudobulk_providers.PseudobulkResult.
    #: dataset_identity` (species/gene_id_namespace/counts_semantics/...)
    #: so a downstream GeneRanking/GSEA consumer never has to guess the
    #: species or gene namespace from a file name. Per
    #: `review:20261014_203915` F1.
    dataset_identity: dict = field(default_factory=dict)
    design_type: str = ""
    gene_ids: tuple = ()
    #: Ordered row identity ACTUALLY used in `design_matrix` -- a SUBSET of
    #: the source pseudobulk's `sample_rows`, selected and reordered by
    #: `sample_id` lookup, never by position. Empty on `status='failed'`.
    sample_ids: tuple = ()
    donor_ids: tuple = ()
    group_values: tuple = ()
    counts_path: str = ""
    counts_file_sha256: str = ""
    counts_content_sha256: str = ""
    design_matrix: tuple = ()
    design_matrix_sha256: str = ""
    column_names: tuple = ()
    formula: str = ""
    encoding_convention: str = ""
    contrast_name: str = ""
    contrast_vector: tuple = ()
    contrast_direction: str = ""
    rank: int = 0
    residual_degrees_of_freedom: int = 0
    n_samples: int = 0
    n_distinct_donors: int = 0
    backend: str = "design_matrix_preparation_only"
    backend_version: str = "1"
    source_code_revision_or_hash: str = ""
    dependency_versions: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        d = dict(self.__dict__)
        for k in ("sample_ids", "donor_ids", "group_values", "gene_ids", "column_names",
                 "contrast_vector"):
            d[k] = list(getattr(self, k))
        d["design_matrix"] = [list(row) for row in self.design_matrix]
        return d


class DEDesignPreparationProvider:
    """See module docstring. `prepare` takes a PATH to an already-written
    `PSEUDOBULK_RESULT.json` and loads it ITSELF, via
    `pseudobulk_providers.read_pseudobulk_result` +
    `pseudobulk_providers.load_pseudobulk_counts` -- per
    `review:20261014_203915` F1, there is no longer a bare
    `counts` ndarray parameter a caller could supply independently of the
    manifest (and so no way for a caller's counts array to disagree with
    the manifest's content/order without `load_pseudobulk_counts`'s own
    hash/identity-binding guards catching it first).
    """

    def prepare(self, *, pseudobulk_result_path: str,
               comparison: ComparisonSpec, sample_structure: SampleStructure,
               out_dir: str) -> PreparedDEDesign:
        if os.path.exists(out_dir) and os.listdir(out_dir):
            raise SystemExit(f"{out_dir} must be new or empty")
        os.makedirs(out_dir, exist_ok=True)

        common = dict(
            comparison=dataclasses.asdict(comparison),
            sample_structure=dataclasses.asdict(sample_structure),
            design_type=sample_structure.design_type,
            backend="design_matrix_preparation_only", backend_version="1",
            source_code_revision_or_hash=_source_revision(),
            dependency_versions={"python": platform.python_version(), "numpy": np.__version__})

        def finish(result: PreparedDEDesign, selected_counts=None) -> PreparedDEDesign:
            return write_prepared_de_design(out_dir, result, counts=selected_counts)

        def fail(reason: str, gene_ids: tuple = ()) -> PreparedDEDesign:
            return finish(PreparedDEDesign(status="failed", reason=reason, gene_ids=gene_ids, **common))

        # `read_pseudobulk_result` is the SAME validated reader every other
        # consumer of a pseudobulk artifact uses -- never a dict the caller
        # assembled or tampered with itself.
        pseudobulk_manifest = _pb.read_pseudobulk_result(pseudobulk_result_path)

        # Per `review:20261014_151747` step 2: a failed/incomplete
        # aggregation must never enter DE design preparation at all.
        if pseudobulk_manifest.get("status") != "complete":
            return fail(
                f"pseudobulk_manifest.status={pseudobulk_manifest.get('status')!r} != "
                "'complete' -- a failed/incomplete pseudobulk aggregation must never "
                "enter DE design preparation")

        gene_ids = tuple(pseudobulk_manifest.get("gene_ids") or ())

        # Per `review:20261014_203915` F2: a declared column name
        # that does not match what the pseudobulk aggregation ACTUALLY used
        # must be refused before any design is built -- `comparison`/
        # `sample_structure`'s column-name fields are never used to
        # extract data (the pseudobulk aggregation already extracted
        # `condition_value`/`donor_or_animal_id`/`sample_id` into
        # `sample_rows`), so a mismatched declaration would otherwise be
        # silently irrelevant rather than refused.
        params = pseudobulk_manifest.get("params") or {}
        mismatches = []
        for declared_name, obj_name, field_name, params_key in _COLUMN_NAME_BINDINGS:
            obj = comparison if obj_name == "comparison" else sample_structure
            declared_value = getattr(obj, field_name)
            actual_value = params.get(params_key)
            if declared_value != actual_value:
                mismatches.append(
                    f"{declared_name}={declared_value!r} != pseudobulk_manifest.params."
                    f"{params_key}={actual_value!r}")
        if mismatches:
            return fail(
                "declared column name(s) do not match what the pseudobulk aggregation "
                f"actually used: {mismatches}", gene_ids=gene_ids)

        # The hash/identity-bound loader this module's own F1 fix now
        # routes through -- counts can no longer disagree with the
        # manifest's content/order without this call itself refusing.
        counts = _pb.load_pseudobulk_counts(pseudobulk_manifest)

        sample_rows = pseudobulk_manifest.get("sample_rows") or []
        problems, fields = _prepare_de_design_core(sample_rows, comparison, sample_structure)
        if problems:
            return fail("; ".join(problems), gene_ids=gene_ids)

        # Per the module docstring: select/reorder the counts submatrix by
        # an explicit sample_id -> row-index LOOKUP, never by position --
        # `sample_rows`' order and the design's row order are never assumed
        # to already coincide just because they often will.
        sample_id_to_index = {row.get("sample_id"): i for i, row in enumerate(sample_rows)}
        idx = [sample_id_to_index[sid] for sid in fields["sample_ids"]]
        selected_counts = counts[idx, :]

        pseudobulk_identity = dict(
            pseudobulk_manifest_path=os.path.abspath(pseudobulk_result_path),
            pseudobulk_manifest_sha256=sha256_file(pseudobulk_result_path),
            aggregation_path=pseudobulk_manifest["aggregation_path"],
            aggregation_sha256=pseudobulk_manifest["aggregation_sha256"],
            identity_sha256=pseudobulk_manifest["identity_sha256"])
        dataset_identity = dict(pseudobulk_manifest.get("dataset_identity") or {})
        design_matrix_sha256 = array_content_sha256(np.array(fields["design_matrix"], dtype=float))

        result = PreparedDEDesign(
            status="prepared", reason="", gene_ids=gene_ids,
            pseudobulk_identity=pseudobulk_identity, dataset_identity=dataset_identity,
            design_matrix_sha256=design_matrix_sha256,
            **common, **fields)
        return finish(result, selected_counts=selected_counts)


def _validate_prepared_de_design_core(payload: dict) -> list[str]:
    """Checks every field EXCEPT `counts_path`/`counts_file_sha256`/
    `counts_content_sha256`, which cannot be known until after the counts
    npz is written -- mirrors `pseudobulk_providers._validate_pseudobulk_
    result_core` (pre-write check, before any npz bytes exist).
    """
    if not isinstance(payload, dict):
        return [f"result must be a dict, got {type(payload).__name__}"]
    problems = []
    if payload.get("schema") != PREPARED_DE_DESIGN_SCHEMA:
        problems.append(f"schema {payload.get('schema')!r} != expected {PREPARED_DE_DESIGN_SCHEMA!r}")
    if payload.get("evidence_type") != "de_design_preparation":
        problems.append(
            f"evidence_type {payload.get('evidence_type')!r} != expected 'de_design_preparation'")
    status = payload.get("status")
    if status not in _VALID_PREP_STATUSES:
        problems.append(f"status {status!r} not in {_VALID_PREP_STATUSES}")
    problems.extend(find_nonfinite_floats(payload))

    for key in ("comparison", "sample_structure"):
        if not isinstance(payload.get(key), dict) or not payload.get(key):
            problems.append(f"{key!r} must be a non-empty dict")

    if status != "prepared":
        if not payload.get("reason"):
            problems.append(f"status={status!r} requires a non-empty 'reason'")
        return problems

    design_type = payload.get("design_type")
    if design_type not in _SUPPORTED_DESIGN_TYPES:
        problems.append(f"design_type={design_type!r} not in {_SUPPORTED_DESIGN_TYPES!r}")

    sample_ids = payload.get("sample_ids")
    donor_ids = payload.get("donor_ids")
    group_values = payload.get("group_values")
    gene_ids = payload.get("gene_ids")
    for name, val in (("sample_ids", sample_ids), ("donor_ids", donor_ids),
                     ("group_values", group_values), ("gene_ids", gene_ids)):
        if not isinstance(val, list) or not val or not all(isinstance(x, str) and x for x in val):
            problems.append(f"{name} must be a non-empty list of non-empty strings, got {val!r}")

    if isinstance(sample_ids, list):
        if len(set(sample_ids)) != len(sample_ids):
            problems.append("sample_ids contains duplicate value(s)")
        for name, val in (("donor_ids", donor_ids), ("group_values", group_values)):
            if isinstance(val, list) and len(val) != len(sample_ids):
                problems.append(f"{name} length {len(val)} does not match len(sample_ids)={len(sample_ids)}")

    n_samples = payload.get("n_samples")
    if not isinstance(n_samples, int) or isinstance(n_samples, bool):
        problems.append(f"n_samples must be an int, got {n_samples!r}")
    elif isinstance(sample_ids, list) and n_samples != len(sample_ids):
        problems.append(f"n_samples={n_samples} does not match len(sample_ids)={len(sample_ids)}")

    n_distinct_donors = payload.get("n_distinct_donors")
    if not isinstance(n_distinct_donors, int) or isinstance(n_distinct_donors, bool):
        problems.append(f"n_distinct_donors must be an int, got {n_distinct_donors!r}")
    elif isinstance(donor_ids, list) and n_distinct_donors != len(set(donor_ids)):
        problems.append(
            f"n_distinct_donors={n_distinct_donors} does not match the distinct count "
            f"in donor_ids ({len(set(donor_ids))})")

    column_names = payload.get("column_names")
    design_matrix = payload.get("design_matrix")
    if not isinstance(column_names, list) or not column_names \
            or not all(isinstance(c, str) and c for c in column_names):
        problems.append(f"column_names must be a non-empty list of non-empty strings, got {column_names!r}")
    dm_rows_ok = isinstance(design_matrix, list) and bool(design_matrix) \
        and all(isinstance(row, list) for row in design_matrix)
    if not dm_rows_ok:
        problems.append("design_matrix must be a non-empty list of row-lists")
    else:
        if isinstance(sample_ids, list) and len(design_matrix) != len(sample_ids):
            problems.append(
                f"design_matrix has {len(design_matrix)} rows, expected "
                f"len(sample_ids)={len(sample_ids)}")
        if isinstance(column_names, list) and any(len(row) != len(column_names) for row in design_matrix):
            problems.append(
                f"design_matrix row length(s) do not all match len(column_names)={len(column_names)}")

    # Per `review:20261014_203915`'s minor item: a non-numeric
    # design_matrix element must become a clean validation problem, never
    # an uncaught `ValueError` from a bare `np.array(..., dtype=float)` --
    # every downstream numeric use of `design_matrix` is gated on this
    # flag rather than on `dm_rows_ok` alone.
    dm_numeric_ok = dm_rows_ok and all(
        isinstance(v, (int, float)) and not isinstance(v, bool)
        for row in design_matrix for v in row)
    if dm_rows_ok and not dm_numeric_ok:
        problems.append("design_matrix contains a non-numeric element")

    design_matrix_sha256 = payload.get("design_matrix_sha256")
    if dm_numeric_ok:
        expected_dm_hash = array_content_sha256(np.array(design_matrix, dtype=float))
        if design_matrix_sha256 != expected_dm_hash:
            problems.append(
                "design_matrix_sha256 does not match a recomputed hash of design_matrix "
                "-- tampered or stale")

    contrast_vector = payload.get("contrast_vector")
    cv_ok = isinstance(contrast_vector, list) and bool(contrast_vector) \
        and all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in contrast_vector)
    if not cv_ok:
        problems.append(f"contrast_vector must be a non-empty list of numbers, got {contrast_vector!r}")
    elif isinstance(column_names, list) and len(contrast_vector) != len(column_names):
        problems.append("contrast_vector length does not match len(column_names)")

    rank = payload.get("rank")
    residual_df = payload.get("residual_degrees_of_freedom")
    if not isinstance(rank, int) or isinstance(rank, bool):
        problems.append(f"rank must be an int, got {rank!r}")
    if not isinstance(residual_df, int) or isinstance(residual_df, bool):
        problems.append(f"residual_degrees_of_freedom must be an int, got {residual_df!r}")

    # Per `review:20261014_203915` F2: `_prepare_de_design_core`
    # is RERUN on the PERSISTED `sample_ids`/`donor_ids`/`group_values`
    # and the PERSISTED `comparison`/`sample_structure` (reconstructed
    # from the SAME payload, never a separately-trusted copy) and every
    # recomputed field is compared against the stored one for EXACT
    # equality -- not merely a design_matrix hash plus an "estimability"
    # check, which is VACUOUS whenever the design is full rank
    # (`c @ pinv(X) @ X == c` identically for ANY `c`, so a reversed or
    # zero contrast previously passed). This subsumes the rank/residual-df/
    # estimability recompute `_prepare_de_design_core` already performs as
    # a side effect of building the design in the first place.
    prerequisites_ok = (
        isinstance(sample_ids, list) and isinstance(donor_ids, list)
        and isinstance(group_values, list) and len(sample_ids) == len(donor_ids) == len(group_values)
        and isinstance(column_names, list) and dm_numeric_ok and cv_ok
        and len(contrast_vector) == len(column_names)
        and all(len(row) == len(column_names) for row in design_matrix)
        and isinstance(n_samples, int) and not isinstance(n_samples, bool)
        and isinstance(n_distinct_donors, int) and not isinstance(n_distinct_donors, bool))
    if prerequisites_ok:
        try:
            comparison_obj = ComparisonSpec(**payload["comparison"])
            sample_structure_obj = SampleStructure(**payload["sample_structure"])
        except TypeError as e:
            problems.append(
                f"comparison/sample_structure could not be parsed back into "
                f"ComparisonSpec/SampleStructure: {e}")
        else:
            recompute_problems, expected = _expected_design_fields_from_identity(
                sample_ids, donor_ids, group_values, comparison_obj, sample_structure_obj)
            if recompute_problems:
                problems.append(
                    "rebuilding the design from the persisted sample_ids/donor_ids/"
                    f"group_values/comparison/sample_structure was refused: {recompute_problems}")
            else:
                stored = dict(
                    sample_ids=tuple(sample_ids), donor_ids=tuple(donor_ids),
                    group_values=tuple(group_values),
                    design_matrix=tuple(tuple(row) for row in design_matrix),
                    column_names=tuple(column_names), formula=payload.get("formula"),
                    encoding_convention=payload.get("encoding_convention"),
                    contrast_name=payload.get("contrast_name"),
                    contrast_vector=tuple(contrast_vector),
                    contrast_direction=payload.get("contrast_direction"),
                    rank=rank, residual_degrees_of_freedom=residual_df,
                    n_samples=n_samples, n_distinct_donors=n_distinct_donors)
                for key in _SEMANTIC_FIELDS_BOUND_TO_IDENTITY:
                    exp_val, got_val = expected[key], stored[key]
                    if key in ("design_matrix", "contrast_vector"):
                        matches = np.allclose(
                            np.array(exp_val, dtype=float), np.array(got_val, dtype=float))
                    else:
                        matches = exp_val == got_val
                    if not matches:
                        problems.append(
                            f"{key} stored as {got_val!r} does not match the value "
                            f"recomputed from the persisted sample_ids/donor_ids/"
                            f"group_values/comparison/sample_structure ({exp_val!r}) -- "
                            "the declared comparison/config does not match the "
                            "persisted design")

    pb_identity = payload.get("pseudobulk_identity")
    if not isinstance(pb_identity, dict) or not all(
            isinstance(pb_identity.get(k), str) and pb_identity.get(k)
            for k in ("pseudobulk_manifest_path", "pseudobulk_manifest_sha256",
                     "aggregation_path", "aggregation_sha256", "identity_sha256")):
        problems.append(
            "pseudobulk_identity must be a dict with non-empty string "
            "'pseudobulk_manifest_path'/'pseudobulk_manifest_sha256'/'aggregation_path'/"
            "'aggregation_sha256'/'identity_sha256' -- the exact upstream pseudobulk "
            "artifact this design was built from")

    dataset_identity = payload.get("dataset_identity")
    if not isinstance(dataset_identity, dict) or not all(
            isinstance(dataset_identity.get(k), str) and dataset_identity.get(k)
            for k in ("species", "gene_id_namespace", "counts_semantics")):
        problems.append(
            "dataset_identity must be a dict with non-empty string 'species'/"
            "'gene_id_namespace'/'counts_semantics', carried forward from the "
            "upstream pseudobulk aggregation -- a downstream consumer must never "
            "guess these from a file name")

    return problems


def validate_prepared_de_design(payload: dict) -> list[str]:
    """Full check, including the fields `_validate_prepared_de_design_core`
    cannot check before the counts npz exists. See
    `load_prepared_de_counts` for the npz<->manifest identity binding
    itself (not re-done here, which has no file access).
    """
    problems = _validate_prepared_de_design_core(payload)
    if not isinstance(payload, dict) or payload.get("status") != "prepared":
        return problems
    for key in ("counts_path", "counts_file_sha256", "counts_content_sha256"):
        if not isinstance(payload.get(key), str) or not payload.get(key):
            problems.append(f"{key!r} must be a non-empty string, got {payload.get(key)!r}")
    return problems


def _nonneg_finite_integer_counts_problems(counts: np.ndarray) -> list[str]:
    """Per `review:20261014_221155` R1: the raw-count-sum contract
    (non-negative, finite, integer-VALUED -- never a dtype requirement,
    since the npz stores `float64`) is checked explicitly and cheaply,
    independent of (and before) the expensive upstream source-binding
    recompute below -- a defense-in-depth net, not a replacement for it.

    Per `review:20261014_231259` R1/P2: uses EXACT equality
    (`==`), never `np.allclose`'s default relative tolerance -- at a
    realistic raw-count scale (e.g. 100000), `np.allclose`'s default
    `rtol=1e-5` alone already tolerates an absolute difference of ~1,
    so a genuinely-wrong neighbouring integer (100001) OR a non-integer
    value (100000.25) would both be silently accepted. Counts identity is
    never an approximate-equality question.
    """
    problems = []
    if not np.all(np.isfinite(counts)):
        problems.append("counts array contains non-finite value(s)")
        return problems
    if np.any(counts < 0):
        problems.append("counts array contains negative value(s) -- raw counts are never negative")
    if not np.all(counts == np.round(counts)):
        problems.append("counts array contains non-integer-valued entries -- raw counts are always integer-valued")
    return problems


def _source_binding_problems(pseudobulk_identity: dict, dataset_identity: dict, comparison: dict,
                             sample_structure: dict, sample_ids: list, donor_ids: list,
                             group_values: list, gene_ids: list, counts: np.ndarray) -> list[str]:
    """Per `review:20261014_221155` R1: every claim a `PreparedDEDesign`
    makes about its upstream pseudobulk source -- `dataset_identity`, the
    declared column names, `sample_ids`/`donor_ids`/`group_values`, `gene_ids`,
    and the counts submatrix itself -- is RE-RESOLVED against the LIVE
    upstream pseudobulk manifest here, never trusted merely because the
    `PreparedDEDesign` payload is internally self-consistent. Factored out
    so both the public write path (`write_prepared_de_design`) and the
    public read/load path (`load_prepared_de_counts`) go through the SAME
    check -- species/namespace/column-name/donor-identity/manifest-hash
    tampering, and a same-shape reversed/negated counts array handed
    directly to the low-level writer, are all caught here, not only at the
    higher-level `DEDesignPreparationProvider.prepare()` entry point.
    """
    manifest_path = pseudobulk_identity.get("pseudobulk_manifest_path")
    try:
        manifest = _pb.read_pseudobulk_result(manifest_path)
    except SystemExit as e:
        return [f"upstream pseudobulk manifest at {manifest_path!r} could not be resolved: {e}"]

    problems = []
    actual_manifest_sha = sha256_file(manifest_path)
    if actual_manifest_sha != pseudobulk_identity.get("pseudobulk_manifest_sha256"):
        problems.append(
            "pseudobulk_manifest_sha256 does not match a freshly recomputed hash of the "
            "live upstream manifest file -- tampered or stale")
    for key in ("aggregation_path", "aggregation_sha256", "identity_sha256"):
        if manifest.get(key) != pseudobulk_identity.get(key):
            problems.append(
                f"pseudobulk_identity[{key!r}]={pseudobulk_identity.get(key)!r} does not "
                f"match the live upstream manifest's own {key!r}={manifest.get(key)!r}")

    live_dataset_identity = manifest.get("dataset_identity") or {}
    for key in ("species", "gene_id_namespace", "counts_semantics"):
        if dataset_identity.get(key) != live_dataset_identity.get(key):
            problems.append(
                f"dataset_identity[{key!r}]={dataset_identity.get(key)!r} does not match "
                f"the live upstream manifest's {key!r}={live_dataset_identity.get(key)!r}")

    params = manifest.get("params") or {}
    for declared_name, obj_name, field_name, params_key in _COLUMN_NAME_BINDINGS:
        obj = comparison if obj_name == "comparison" else sample_structure
        declared_value = obj.get(field_name)
        actual_value = params.get(params_key)
        if declared_value != actual_value:
            problems.append(
                f"{declared_name}={declared_value!r} != the live upstream manifest's "
                f"params.{params_key}={actual_value!r}")

    if list(gene_ids) != list(manifest.get("gene_ids") or ()):
        problems.append("gene_ids does not match the live upstream manifest's gene_ids order")

    sample_rows = manifest.get("sample_rows") or []
    row_by_sid = {r.get("sample_id"): r for r in sample_rows}
    for sid, did, gv in zip(sample_ids, donor_ids, group_values):
        row = row_by_sid.get(sid)
        if row is None:
            problems.append(f"sample_id {sid!r} not found in the live upstream manifest's sample_rows")
            continue
        if row.get("donor_or_animal_id") != did:
            problems.append(
                f"donor_ids entry for sample_id {sid!r} is {did!r}, but the live upstream "
                f"manifest records donor_or_animal_id={row.get('donor_or_animal_id')!r}")
        if row.get("condition_value") != gv:
            problems.append(
                f"group_values entry for sample_id {sid!r} is {gv!r}, but the live upstream "
                f"manifest records condition_value={row.get('condition_value')!r}")

    if problems:
        # The upstream identity/column/row claims are already wrong -- do not also attempt
        # a counts recompute against a manifest that has not been shown trustworthy.
        return problems

    try:
        live_counts = _pb.load_pseudobulk_counts(manifest)
    except SystemExit as e:
        return [f"upstream pseudobulk counts could not be loaded: {e}"]

    idx_by_sid = {r.get("sample_id"): i for i, r in enumerate(sample_rows)}
    idx = [idx_by_sid[sid] for sid in sample_ids]
    expected_counts = live_counts[idx, :]
    # Per `review:20261014_231259` R1/P2: EXACT equality, never `np.allclose` --
    # this compares the SAME upstream content, not two independently-computed floats that
    # could legitimately differ by a tiny rounding amount. See
    # `_nonneg_finite_integer_counts_problems` for why a relative tolerance is unsafe here.
    if not np.array_equal(np.asarray(counts, dtype=float), expected_counts):
        problems.append(
            "counts does not match the submatrix freshly recomputed, by sample_id, from the "
            "live upstream pseudobulk source -- tampered or stale")
    return problems


def write_prepared_de_design(out_dir: str, result: PreparedDEDesign,
                             counts: np.ndarray | None = None) -> PreparedDEDesign:
    """Mirrors `pseudobulk_providers._write_pseudobulk_result`: the
    pre-npz-bytes structural check runs BEFORE any npz is written, so an
    invalid `PreparedDEDesign` never leaves an orphan counts npz with no
    corresponding manifest. Per `review:20261014_221155` R1, this
    is ALSO the public entry point a caller could use directly (bypassing
    `DEDesignPreparationProvider.prepare()` entirely) to persist a
    `PreparedDEDesign` alongside an arbitrary `counts` array -- so a
    `status='prepared'` write now re-resolves the claimed upstream
    pseudobulk source and refuses if `counts` does not match what that
    live source actually contains for the claimed `sample_ids`, never
    trusting the caller's `counts` argument merely because it has the
    right shape.
    """
    pre_payload = result.to_dict()
    pre_payload["schema"] = PREPARED_DE_DESIGN_SCHEMA
    pre_payload["generated_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    pre_problems = _validate_prepared_de_design_core(pre_payload)
    if pre_problems:
        raise ValueError(f"refusing to write an invalid prepared DE design: {pre_problems}")

    if result.status == "prepared":
        if counts is None:
            raise ValueError("status='prepared' requires the selected counts submatrix (counts=...)")
        expected_shape = (len(result.sample_ids), len(result.gene_ids))
        if counts.shape != expected_shape:
            raise ValueError(f"counts shape {counts.shape} != expected {expected_shape}")
        numeric_problems = _nonneg_finite_integer_counts_problems(counts)
        if numeric_problems:
            raise ValueError(f"refusing to write an invalid prepared DE design: {numeric_problems}")
        binding_problems = _source_binding_problems(
            result.pseudobulk_identity, result.dataset_identity, result.comparison,
            result.sample_structure, list(result.sample_ids), list(result.donor_ids),
            list(result.group_values), list(result.gene_ids), counts)
        if binding_problems:
            raise ValueError(f"refusing to write an invalid prepared DE design: {binding_problems}")
        npz_path = os.path.join(out_dir, "prepared_de_counts.npz")
        tmp = npz_path + ".inprogress"
        sample_ids_arr = np.array(list(result.sample_ids))
        gene_ids_arr = np.array(list(result.gene_ids))
        with open(tmp, "wb") as fh:
            np.savez(fh, counts=np.ascontiguousarray(counts, dtype=np.float64),
                    sample_ids=sample_ids_arr, gene_ids=gene_ids_arr)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, npz_path)
        result = replace(result, counts_path=npz_path, counts_file_sha256=sha256_file(npz_path),
                         counts_content_sha256=array_content_sha256(counts))

    payload = result.to_dict()
    payload["schema"] = PREPARED_DE_DESIGN_SCHEMA
    payload["generated_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    problems = validate_prepared_de_design(payload)
    if problems:
        raise ValueError(f"refusing to write an invalid prepared DE design: {problems}")

    result_path = os.path.join(out_dir, "PREPARED_DE_DESIGN.json")
    tmp2 = result_path + ".inprogress"
    with open(tmp2, "w") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True, allow_nan=False)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp2, result_path)
    return result


def read_prepared_de_design(path: str) -> dict:
    if not os.path.exists(path):
        raise SystemExit(f"prepared DE design not found: {path}")
    with open(path) as fh:
        try:
            payload = json.load(fh, parse_constant=_reject_nonfinite_constant)
        except ValueError as e:
            raise SystemExit(f"prepared DE design at {path} failed to parse: {e}")
    problems = validate_prepared_de_design(payload)
    if problems:
        raise SystemExit(f"prepared DE design at {path} failed validation: {problems}")
    return payload


def _prepared_counts_identity_binding_problems(npz_sample_ids: list, npz_gene_ids: list,
                                               manifest_sample_ids: list,
                                               manifest_gene_ids: list) -> list[str]:
    """Factored out so this specific comparison can be mutation-killed
    directly. See module docstring: direct npz<->manifest identity binding,
    independent of the whole-file hash and of the array's own content hash.
    """
    problems = []
    if npz_sample_ids != manifest_sample_ids:
        problems.append(
            f"the npz's own stored sample_ids {npz_sample_ids!r} do not match the "
            f"manifest's sample_ids {manifest_sample_ids!r}")
    if npz_gene_ids != manifest_gene_ids:
        problems.append(
            f"the npz's own stored gene_ids {npz_gene_ids!r} do not match the "
            f"manifest's gene_ids {manifest_gene_ids!r}")
    return problems


def load_prepared_de_counts(manifest: dict) -> np.ndarray:
    """Loads the counts npz a `status='prepared'` manifest references.
    Three independent npz-level checks, mirroring `pseudobulk_providers.
    load_pseudobulk_counts`'s R1 fix: (1) the whole-file hash
    (`counts_file_sha256`), (2) the npz's OWN stored `sample_ids`/
    `gene_ids` arrays compared DIRECTLY against the manifest's claims (not
    merely a total/shape recompute, which an equal-marginals-style swap
    cannot distinguish), and (3) a content-only hash of the array itself
    (`counts_content_sha256`) recomputed from the loaded array. Per
    `review:20261014_221155` R1, a FOURTH check then re-resolves
    the claimed upstream pseudobulk source LIVE and recomputes the
    expected counts submatrix from it (`_source_binding_problems`) -- the
    first three checks alone only prove the on-disk artifact is internally
    self-consistent (its own JSON agrees with its own npz), never that the
    identity/column-name/counts claims it makes about its UPSTREAM source
    are actually still true.
    """
    if manifest.get("status") != "prepared":
        raise SystemExit(
            f"load_prepared_de_counts requires status='prepared', got "
            f"{manifest.get('status')!r}")
    path = manifest["counts_path"]
    actual_file_sha = sha256_file(path)
    if actual_file_sha != manifest["counts_file_sha256"]:
        raise SystemExit(
            f"{path} has file hash {actual_file_sha!r}, which does not match the "
            f"manifest's recorded counts_file_sha256 {manifest['counts_file_sha256']!r}")
    with np.load(path, allow_pickle=False) as npz:
        counts = npz["counts"]
        missing = [k for k in ("sample_ids", "gene_ids") if k not in npz.files]
        if missing:
            raise SystemExit(f"{path} is missing {missing}")
        npz_sample_ids = [str(x) for x in npz["sample_ids"]]
        npz_gene_ids = [str(x) for x in npz["gene_ids"]]

    identity_problems = _prepared_counts_identity_binding_problems(
        npz_sample_ids, npz_gene_ids, list(manifest["sample_ids"]), list(manifest["gene_ids"]))
    if identity_problems:
        raise SystemExit("; ".join(identity_problems))

    expected_shape = (len(manifest["sample_ids"]), len(manifest["gene_ids"]))
    if counts.shape != expected_shape:
        raise SystemExit(f"counts array shape {counts.shape} != expected {expected_shape}")
    numeric_problems = _nonneg_finite_integer_counts_problems(counts)
    if numeric_problems:
        raise SystemExit("; ".join(numeric_problems))

    recomputed_content_sha = array_content_sha256(counts)
    if recomputed_content_sha != manifest["counts_content_sha256"]:
        raise SystemExit(
            f"loaded counts array content hash {recomputed_content_sha!r} does not match "
            f"the manifest's recorded counts_content_sha256 "
            f"{manifest['counts_content_sha256']!r}")

    binding_problems = _source_binding_problems(
        manifest["pseudobulk_identity"], manifest["dataset_identity"], manifest["comparison"],
        manifest["sample_structure"], list(manifest["sample_ids"]), list(manifest["donor_ids"]),
        list(manifest["group_values"]), list(manifest["gene_ids"]), counts)
    if binding_problems:
        raise SystemExit("; ".join(binding_problems))
    return counts
