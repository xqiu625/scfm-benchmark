"""Gene-ranking evidence-stream contracts for the scFM discovery framework.

Per `review:20261013_135634` section 3: two SEPARATE evidence
streams, never conflated with each other and never conflated with the PCA
baseline (`embedding_providers.py`, a control, not a source of either kind
of evidence):

  - `ExpressionDERankingProvider` ("表达流"): under an EXPLICIT biological
    comparison, which genes' expression changed. The INTENDED estimator is
    donor/animal-aware pseudobulk aggregation plus a shrunken log-fold-
    change-style estimator -- NOT independent-cell Wilcoxon, NOT an
    arbitrary signed -log10(p) substitute, by default (the prior handoff's
    own overclaim this round corrects). A design where every donor
    contributes to only one compared group cannot support a donor-blocked
    term; `check_design_estimability` refuses that case rather than
    silently fitting a degenerate model or silently treating cells as
    independent replicates.

  - `ModelMaskingRankingProvider` ("模型流"): `Delta(g) = s_full -
    s_masked(g)` for a FIXED, bound scoring function. Masking gene `g`
    means masking a GENE INPUT to the model -- this interface has no
    cell-selection/QC-mask parameter at all, by construction. The score
    function, the masked-position reference/token semantics, the
    perturbable gene set, and the donor/animal aggregation rule must all be
    bound to the result, never re-decided per call or per gene.

Both streams produce the SAME `enrichment_providers.GeneRanking` object the
already-fixed/tested GSEA chain consumes, but each binds its OWN
`evidence_type` and comparison/score/sample-structure metadata alongside
it -- a downstream GSEA run can never mistake one stream's ranking for the
other's just because they share a `GeneRanking` shape, and the two streams'
enrichment results overlapping is never itself treated as independent
corroboration (they may share cells/donors and are not independent
measurements).

This round implements TYPES, VALIDATION, and SERIALIZATION ONLY. Neither
provider's `run()` performs a real statistical or model computation -- both
always return `status="failed"` (a structural/design problem was found) or
`status="not_run"` (inputs passed validation, but no real backend exists
yet) and NEVER a fabricated `status="complete"` result. No main model
checkpoint is pinned this round, and no donor/animal grouping, masking
score, or comparison has been confirmed by the owner -- see the decision
table in this round's implementation handoff.

Per `review:20261013_210131` F1: a donor-blocked `~ donor +
condition` estimability check must not be applied unconditionally to every
comparison. `SampleStructure.design_type` is now an EXPLICIT, code-
recognized field (one of `_VALID_DESIGN_TYPES`) that `check_design_
estimability` branches on -- never a free-text `repeated_sampling_or_
pairing`/`estimator=` string, which is descriptive only and never changes
which check runs. `"independent_donors_unpaired"` refuses a donor that
actually straddles both compared groups (that contradicts the declared
design); `"paired_or_repeated_within_donor"` refuses unless EVERY donor in
the comparison carries BOTH levels -- a single cross-group donor is no
longer enough to call a donor-blocked design estimable.
`ExpressionDERankingProvider.run` also refuses a `comparison.group_column`
that disagrees with `sample_structure.group_column`, and flags declared-
but-unchecked `covariate_columns` as explicitly UNVERIFIED in the reason,
rather than silently passing them.

Per `review:20261013_210131` F4: `ExpressionDEResult` and
`ModelMaskingResult` now record a content identity of their actual inputs
(`perturbable_genes` + its hash; `group_values`/`donor_values` hashes) --
two different inputs to the same call shape now serialize to two different
records, even while no real backend exists yet.

Per `review:20261013_222450` (review of `0fec319`): three gaps in
the above remained.

  - R1: `ExpressionDERankingProvider.run`'s entry point only checked that
    SOME sample carried a compared group value, never that BOTH
    `reference_level` and `comparison_level` were present; an empty-string
    or `None` donor/group value reached Python's `sorted()` uncaught (a
    `None`/`str` mix raises `TypeError`) instead of a structured `failed`
    record. `check_design_estimability` now runs `_validate_group_donor_
    values` first -- length match, every group/donor value a non-empty
    string, BOTH comparison levels actually present -- before building any
    donor grouping or calling `sorted()`.

  - R2: the `paired_or_repeated_within_donor` refusal message claimed a
    partially-paired `~ donor + condition` design "cannot be estimated",
    which is not generally true (a design matrix can remain full rank with
    some, not all, donors paired). This round only SUPPORTS the
    fully-paired case -- the message now says so explicitly
    ("UNSUPPORTED by this round's check", not "mathematically
    inestimable") and states that straddling/non-straddling donors are
    never silently dropped by the check itself.

  - R3: `ExpressionDEResult`/`ModelMaskingResult` carried only a length +
    hash summary, with no way to recover or independently re-check the
    actual ordered input from the record, no dedicated writer/reader, and
    an identity hash built with `json.dumps(..., default=str)` (silently
    stringifies an unrecognized object instead of refusing it). Both
    results now carry an explicit `input_mapping`/`perturbable_genes`
    round-trippable list plus a recomputed-and-compared content hash, and
    `write_expression_de_result`/`read_expression_de_result` and
    `write_model_masking_result`/`read_model_masking_result` validate the
    SAME way at both the write and the read boundary -- a tampered mapping,
    hash, or identity field is rejected on read, matching
    `embedding_providers.write/read_embedding_result`'s pattern. Identity
    hashing no longer uses `default=str`; an input that is not natively
    JSON-serializable is a distinct, explicit "identity unavailable" state,
    never a silently stringified guess.

Per `review:20261014_100436` F3: the R3 envelope's write/read
validators were still asymmetric -- an illegitimate "complete" record
(`ranking={}`, or `backend='not_implemented'`) was accepted, while a
genuinely `failed` record carrying a diagnostic `None` donor_value (the
REAL provider's own record of that exact invalid input) was REJECTED by
the same strict checks a `complete`/`not_run` record needs. Both are fixed
together: `_validate_evidence_result_common` now re-validates `ranking`
against `GeneRanking.validate()` itself (an empty dict is no longer
accepted as "non-null") and refuses the `backend='not_implemented'`/
`backend_version='n/a'` placeholders on a `status='complete'` record;
`validate_expression_de_result`/`validate_model_masking_result` now gate
their strict per-row/per-field VALUE checks (not the key-presence/hash
checks, which still apply unconditionally) to `status in ('complete',
'not_run')` -- a `'failed'` record's `input_mapping`/`masking_score` is a
diagnostic snapshot of the actual (possibly invalid) input, never a
record that claims to BE a validated one. `MaskingScoreSpec.validate()`
and `SampleStructure`/`ComparisonSpec`'s dict forms are now checked by
TYPE, not truthiness (`score_definition=True` is truthy but not a
string); nested `sample_structure.design_type` is re-checked against
`_VALID_DESIGN_TYPES` on every read, not only inside `SampleStructure
.validate()` at construction time.

Per `review:20261014_114357` (review of `0a40d4f`): two further
gaps, both the SAME root cause recurring in a new place.

  - R2: `validate_expression_de_result`'s `n_samples` check had the shape
    `isinstance(n_samples, int) and ... != ...` with no `else` branch, so
    `n_samples='two'`/`None` silently SKIPPED the consistency check
    instead of being rejected for its type -- the exact pattern F3 fixed
    for `backend`/`score_definition` elsewhere in this module, recurring
    here unaddressed. Now checked unconditionally first.

  - R3: that same `n_samples == len(input_mapping)` check ran
    UNCONDITIONALLY regardless of status, which made a genuinely `'failed'`
    length-mismatch record (`group_values`/`donor_values` of DIFFERENT
    lengths, so `_build_input_mapping` correctly returns `()` while
    `n_samples=len(group_values)` is non-zero) unwritable -- comparing
    `n_samples` to an intentionally-empty mapping's length always failed.
    `ModelMaskingResult`'s `perturbable_genes` had the matching gap: its
    per-entry value check ran unconditionally, so a genuinely `'failed'`
    record for `perturbable_genes=(None,)` could not be written either,
    even though `masking_score`'s OWN nested fields already got the F3
    gating treatment. Both fixed the same way F3 fixed `masking_score`/
    `input_mapping`'s per-row values: a NEW `input_lengths`/`input_mapping_
    buildable` pair of fields records the raw lengths and whether they
    agreed, ALWAYS (never dropped on failure); the length-consistency
    check now runs only when a mapping was actually buildable, and
    `perturbable_genes`' value check is gated by `require_valid_values`
    exactly like `masking_score`'s fields already were.

Per `review:20261014_151747` P2 (review of `3131146`): R3's
`input_lengths`/`input_mapping_buildable` fix closed the length-vs-mapping
cross-check, but nothing tied `input_mapping_buildable=False` to `status`
itself -- a genuinely `'failed'` length-mismatch record, written and read
back correctly, could have ONLY its `status` field rewritten to `'not_run'`
(or `'complete'`) on disk and still be ACCEPTED, since every other
cross-check (mapping emptiness, hash, `n_samples` vs. `input_lengths`)
still held. `validate_expression_de_result` now rejects
`input_mapping_buildable=False` combined with `status in ('not_run',
'complete')` directly: this module defines `'not_run'` as "inputs PASSED
validation, no backend ran yet" (see `ExpressionDERankingProvider.run`),
and an unbuildable input never passed validation.
"""
from __future__ import annotations

import dataclasses
import json
import math
import os
import platform
import sys
import time
from dataclasses import dataclass, field

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from artifact_schema import (find_nonfinite_floats, resolve_source_revision,  # noqa: E402
                             safe_json_identity, sha256_text)
from enrichment_providers import GeneRanking  # noqa: E402

SOURCE_MODULES = ("gene_ranking_providers.py", "enrichment_providers.py", "artifact_schema.py")

_VALID_RESULT_STATUSES = ("complete", "failed", "not_run")
EXPRESSION_DE_RESULT_SCHEMA = "scfm-expression-de-result/1"
MODEL_MASKING_RESULT_SCHEMA = "scfm-model-masking-result/1"


def _reject_nonfinite_constant(token: str):
    raise ValueError(f"non-finite JSON constant {token!r} is not accepted")


def _bounded_n_cpus() -> int:
    """Per `review:20261015_142903` R4: pydeseq2's own
    `DefaultInference` (`utils.get_num_processes`, verified directly
    against the real v0.5.4 source) defaults `n_cpus=None` to "the maximum
    number of available cpus" -- i.e. `joblib.Parallel(n_jobs=...)` would
    try to use every core the NODE has, not just the cores this project's
    SLURM allocation actually reserved. Setting `OMP_NUM_THREADS`/
    `OPENBLAS_NUM_THREADS` (this project's standing convention) bounds BLAS
    thread counts, never joblib's own process pool -- a SEPARATE knob.

    Reads the SLURM allocation's own CPU count (`SLURM_CPUS_PER_TASK`, then
    `SLURM_CPUS_ON_NODE`) and returns it; if NEITHER is set (not running
    under SLURM -- e.g. a local/interactive check), returns 1, the SAFE
    default, never `os.cpu_count()` (which would silently oversubscribe a
    shared, unallocated environment).
    """
    for var in ("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE"):
        value = os.environ.get(var)
        if value:
            try:
                n = int(value)
            except ValueError:
                continue
            if n > 0:
                return n
    return 1


#: Per `review:20261014_100436` F1/F3: moved to `artifact_schema.py`
#: so `pseudobulk_providers.py` can reuse the same strict-identity policy
#: without a backwards layering dependency on this module; kept under its
#: original private name here so every existing call site is unchanged.
_safe_json_identity = safe_json_identity


def _ranking_to_dict(ranking) -> dict | None:
    if ranking is None:
        return None
    return {
        "species": ranking.species, "gene_id_namespace": ranking.gene_id_namespace,
        "ranking_definition": ranking.ranking_definition,
        "genes": list(ranking.genes), "scores": list(ranking.scores),
    }


def _source_revision() -> str:
    return resolve_source_revision(HERE, SOURCE_MODULES)


# ===========================================================================
# Expression evidence stream
# ===========================================================================

#: Per `review:20261013_210131` F1: a fixed, code-recognized set --
#: never free text, never inferred from a caller's `estimator=` string. An
#: unknown/unstated value is refused rather than silently defaulting to
#: either branch.
_VALID_DESIGN_TYPES = ("independent_donors_unpaired", "paired_or_repeated_within_donor")


@dataclass(frozen=True)
class SampleStructure:
    """Donor/animal vs. sample vs. cell identity, EXPLICIT per comparison --
    this interface requires the caller to state it, never assumes "every
    cell is an independent biological replicate" as an implicit default.
    """
    donor_or_animal_id_column: str
    sample_id_column: str
    #: Free text, descriptive only -- e.g. "none_single_sample_per_donor",
    #: "paired:pre_post_treatment", "repeated:3_timepoints_per_donor". Never
    #: read by `check_design_estimability`; see `design_type` for the field
    #: that actually drives which check runs.
    repeated_sampling_or_pairing: str
    group_column: str
    #: One of `_VALID_DESIGN_TYPES`, EXPLICIT -- the code contract (not this
    #: text, and not a caller's `estimator=` string) decides which
    #: estimability check `check_design_estimability` applies.
    design_type: str
    covariate_columns: tuple = ()

    def validate(self) -> list[str]:
        problems = []
        for f in ("donor_or_animal_id_column", "sample_id_column",
                  "repeated_sampling_or_pairing", "group_column", "design_type"):
            if not getattr(self, f):
                problems.append(f"{f} is empty -- state it explicitly")
        if self.design_type and self.design_type not in _VALID_DESIGN_TYPES:
            problems.append(
                f"design_type={self.design_type!r} not in {_VALID_DESIGN_TYPES!r} -- "
                "an unrecognized design type cannot be checked for estimability")
        return problems


@dataclass(frozen=True)
class ComparisonSpec:
    """The EXACT biological comparison -- never implied by column order or
    left to a backend's default factor-level ordering.
    """
    group_column: str
    reference_level: str
    comparison_level: str
    #: e.g. "comparison_level vs. reference_level, positive = higher in
    #: comparison_level" -- states what a POSITIVE score means.
    direction_definition: str

    def validate(self) -> list[str]:
        problems = []
        if not self.group_column:
            problems.append("group_column is empty")
        if not self.reference_level:
            problems.append("reference_level is empty")
        if not self.comparison_level:
            problems.append("comparison_level is empty")
        if self.reference_level and self.comparison_level \
                and self.reference_level == self.comparison_level:
            problems.append("reference_level must differ from comparison_level")
        if not self.direction_definition:
            problems.append(
                "direction_definition is empty -- state what a positive score means")
        return problems


def _validate_group_donor_values(group_values: list, donor_values: list,
                                 comparison: ComparisonSpec) -> list[str]:
    """Structural validation of the ACTUAL group/donor LABEL VALUES, run
    before any donor grouping is built or `sorted()` is called -- per
    `review:20261013_222450` R1: a missing comparison group (only
    `reference_level` samples present, say, with zero `comparison_level`
    samples), an empty-string donor ID, or a non-string donor/group value
    must become a structured `failed` reason here, never a silent
    "not_run, design validated" (the prior round only checked that SOME
    compared-group sample existed, not that BOTH levels were present) or an
    uncaught `TypeError` from `sorted()` comparing `None` against `str`.
    """
    problems = []
    if len(group_values) != len(donor_values):
        problems.append(
            f"group_values ({len(group_values)}) and donor_values "
            f"({len(donor_values)}) length mismatch")
        return problems  # nothing further can be checked pairwise
    bad_group = [i for i, g in enumerate(group_values) if not isinstance(g, str) or not g]
    if bad_group:
        problems.append(
            f"group_values has {len(bad_group)} non-string/empty value(s) at "
            f"index/indices {bad_group[:5]}" + (" (truncated)" if len(bad_group) > 5 else ""))
    bad_donor = [i for i, d in enumerate(donor_values) if not isinstance(d, str) or not d]
    if bad_donor:
        problems.append(
            f"donor_values has {len(bad_donor)} non-string/empty value(s) at "
            f"index/indices {bad_donor[:5]}" + (" (truncated)" if len(bad_donor) > 5 else ""))
    if problems:
        # Type/empty-value problems must be resolved before any sorted()
        # level-membership check below, which a None/non-string value would
        # make raise instead of returning a structured reason.
        return problems
    present_levels = set(group_values)
    missing = [lvl for lvl in (comparison.reference_level, comparison.comparison_level)
              if lvl and lvl not in present_levels]
    if missing:
        problems.append(
            f"group_values is missing comparison level(s) {missing} -- BOTH "
            f"{comparison.reference_level!r} and {comparison.comparison_level!r} must "
            "be present to check this comparison, not just one of them")
    return problems


def check_design_estimability(group_values: list, donor_values: list,
                              comparison: ComparisonSpec, design_type: str) -> list[str]:
    """Branches on the EXPLICIT `design_type` (one of `_VALID_DESIGN_TYPES`)
    rather than mechanically applying a donor-blocked check to every
    comparison -- per `review:20261013_210131` F1, correcting the
    prior round's own gap (the check used to run unconditionally and
    treated even ONE cross-group donor as sufficient to call the whole
    design estimable). Per `review:20261013_222450` R1: runs
    `_validate_group_donor_values` FIRST -- length match, non-empty string
    group/donor values, BOTH comparison levels present -- before building
    any donor grouping or calling `sorted()`.

    `"paired_or_repeated_within_donor"`: refuses unless EVERY donor
    contributing to this comparison carries BOTH `reference_level` and
    `comparison_level`. Per `review:20261013_222450` R2: this is a
    SUPPORT-SCOPE limit of this round's check, not a mathematical claim --
    a `~ donor + condition` design matrix CAN remain full rank with some,
    not all, donors paired; this round simply does not recognize that case
    as a supported design, and does not silently drop the straddling or
    non-straddling donor(s) that make it partially paired.

    `"independent_donors_unpaired"`: refuses if ANY donor actually
    straddles both compared levels -- that contradicts the declared
    independent-donor design (checked against the REAL input, not just
    trusted from the declaration); such a donor must be either excluded or
    the design redeclared as `"paired_or_repeated_within_donor"`.
    """
    if design_type not in _VALID_DESIGN_TYPES:
        return [f"design_type={design_type!r} not in {_VALID_DESIGN_TYPES!r} -- "
               "an unrecognized design type cannot be checked for estimability"]
    value_problems = _validate_group_donor_values(group_values, donor_values, comparison)
    if value_problems:
        return value_problems
    by_donor: dict = {}
    for g, d in zip(group_values, donor_values):
        if g not in (comparison.reference_level, comparison.comparison_level):
            continue
        by_donor.setdefault(d, set()).add(g)
    donors_in_both = sorted(d for d, levels in by_donor.items() if len(levels) > 1)
    donors_in_one = sorted(d for d, levels in by_donor.items() if len(levels) == 1)

    if design_type == "paired_or_repeated_within_donor":
        if donors_in_one:
            return [
                f"design_type='paired_or_repeated_within_donor' declared, but "
                f"{len(donors_in_one)} donor(s) contribute to only ONE of "
                f"{comparison.reference_level!r}/{comparison.comparison_level!r}: "
                f"{donors_in_one[:5]}" + (" (truncated)" if len(donors_in_one) > 5 else "") +
                " -- partially-paired designs are UNSUPPORTED by this round's "
                "check (NOT necessarily mathematically inestimable: a "
                "'~ donor + condition' design can remain full rank with some, "
                "not all, donors paired -- this round only recognizes the "
                "fully-paired case). Declare design_type="
                "'independent_donors_unpaired' if these are genuinely "
                "independent samples, or restrict the comparison to the "
                "fully-paired donor subset explicitly; the partially-paired "
                "donor(s) listed above are NOT silently dropped by this check"]
        return []

    # design_type == "independent_donors_unpaired"
    if donors_in_both:
        return [
            f"design_type='independent_donors_unpaired' declared, but "
            f"{len(donors_in_both)} donor(s) actually contribute to BOTH "
            f"{comparison.reference_level!r} and {comparison.comparison_level!r}: "
            f"{donors_in_both[:5]}" + (" (truncated)" if len(donors_in_both) > 5 else "") +
            " -- this contradicts an independent-donor design; declare "
            "design_type='paired_or_repeated_within_donor' instead if these are "
            "genuinely repeated measures from the same biological unit, or "
            "exclude the straddling donor(s) explicitly"]
    return []


def _sample_label_identity(group_values: list, donor_values: list) -> tuple:
    """Content identity of the ACTUAL `group_values`/`donor_values` this run
    was called with -- per `review:20261013_210131` F4, so two
    different inputs to the same call shape serialize to two different
    `ExpressionDEResult` records, even while no real backend exists yet.
    Factored out as its own function so a mutation that stops binding this
    identity can be demonstrated directly (see
    `test_mutation_kill_expression_de_sample_label_identity`). Per
    `review:20261013_222450` R3: uses `_safe_json_identity`, never
    `default=str` -- this is a CONVENIENCE summary hash only; the
    authoritative, round-trippable record of the actual input is
    `_build_input_mapping`/`ExpressionDEResult.input_mapping` below.
    """
    return (
        len(group_values),
        _safe_json_identity(list(group_values)) or "",
        _safe_json_identity(list(donor_values)) or "",
    )


def _input_mapping_buildable(sample_ids: list | None, group_values: list, donor_values: list) -> bool:
    """Per `review:20261014_114357` R3: whether `_build_input_
    mapping` CAN build a real row-aligned mapping from these three inputs,
    i.e. whether their lengths actually agree -- shared by `_build_input_
    mapping` (so the two never disagree with each other) and `run()` (so
    this fact is recorded on the result as `input_mapping_buildable`,
    separately from `n_samples`/`input_lengths`). A caller passing inputs
    whose lengths disagree (e.g. 2 `group_values` but 1 `donor_values`) is
    a genuine, diagnosable input error -- NOT a 0-sample run -- and must be
    distinguishable from one on the persisted record.
    """
    return len(group_values) == len(donor_values) \
        and (sample_ids is None or len(sample_ids) == len(group_values))


def _build_input_mapping(sample_ids: list | None, group_values: list, donor_values: list) -> tuple:
    """The ACTUAL ordered (row identity, group_value, donor_value) triples
    this run was called with -- per `review:20261013_222450` R3: a
    length+hash summary alone cannot be read back into "what the input
    actually was", so a reader has nothing to independently re-check it
    against. `sample_ids=None` means the caller did not supply real sample
    identity; each row then carries only its `row_index` (never a
    FABRICATED sample ID), and `ExpressionDEResult.sample_identity_
    available=False` records that explicitly so a reader can tell "no
    sample identity given" apart from "sample identity given but lost".
    Returns `()` (no mapping built) when `group_values`/`donor_values`/
    `sample_ids` lengths disagree -- that mismatch is reported as its own
    `problems` entry elsewhere; a partial/misaligned mapping must not be
    built from it.
    """
    if not _input_mapping_buildable(sample_ids, group_values, donor_values):
        return ()
    rows = []
    for i in range(len(group_values)):
        row = {"row_index": i, "group_value": group_values[i], "donor_value": donor_values[i]}
        if sample_ids is not None:
            row["sample_id"] = sample_ids[i]
        rows.append(row)
    return tuple(rows)


def _input_mapping_sha256(input_mapping: tuple) -> str:
    return sha256_text(json.dumps(list(input_mapping), sort_keys=True))


@dataclass(frozen=True)
class ExpressionDEResult:
    status: str  # "complete" | "failed" | "not_run"
    reason: str = ""
    evidence_type: str = "expression_de"
    ranking: object = None  # GeneRanking, only when status == "complete"
    comparison: dict = field(default_factory=dict)
    sample_structure: dict = field(default_factory=dict)
    estimator: str = ""
    #: See `_sample_label_identity` -- a convenience summary hash of the
    #: actual group/donor label inputs, not just the comparison/sample_
    #: structure CONFIGURATION. Superseded as the AUTHORITATIVE record by
    #: `input_mapping` below (per R3); kept so existing callers/records are
    #: not broken.
    n_samples: int = 0
    group_values_sha256: str = ""
    donor_values_sha256: str = ""
    #: Per `review:20261014_114357` R3: the RAW lengths of the
    #: actual `group_values`/`donor_values`/`sample_ids` this run was
    #: called with, recorded REGARDLESS of whether those lengths agree --
    #: `n_samples` alone (just `len(group_values)`) cannot tell a reader
    #: whether `input_mapping` is empty because there were genuinely 0
    #: samples, or because the inputs disagreed in length and no row-
    #: aligned mapping could be built at all. `sample_ids` is `None` when
    #: the caller did not supply sample identity (not 0 -- see
    #: `sample_identity_available`).
    input_lengths: dict = field(default_factory=dict)
    #: Per `review:20261014_114357` R3: whether `input_mapping`
    #: below is a REAL row-aligned mapping (lengths agreed) or a FORCED
    #: empty diagnostic placeholder (lengths disagreed) -- see
    #: `_input_mapping_buildable`. When `False`, `input_mapping == ()` is
    #: not "0 samples ran cleanly"; it is "no row-aligned mapping could be
    #: built from this invalid input", and the record's `status` must be
    #: `'failed'`.
    input_mapping_buildable: bool = True
    #: Per `review:20261013_222450` R3: the ACTUAL ordered input
    #: this run was called with -- `row_index`/`group_value`/`donor_value`
    #: per row, plus `sample_id` per row when `sample_identity_available`.
    #: See `_build_input_mapping`.
    sample_identity_available: bool = False
    input_mapping: tuple = ()
    input_mapping_sha256: str = ""
    backend: str = "unknown"
    backend_version: str = "unknown"
    params: dict = field(default_factory=dict)
    source_code_revision_or_hash: str = ""
    dependency_versions: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        d = dict(self.__dict__)
        d["ranking"] = _ranking_to_dict(self.ranking)
        d["input_mapping"] = list(self.input_mapping)
        return d


class ExpressionDERankingProvider:
    """Expression evidence stream. See module docstring. Implements
    validation/estimability-checking/serialization only this round -- `run`
    never performs a real pseudobulk/shrunken-LFC computation.
    """

    #: The DEFAULT/INTENDED estimator direction. A caller naming a different
    #: `estimator=` must say so explicitly; this constant is never silently
    #: substituted for a real computation, and is not itself a computation.
    INTENDED_ESTIMATOR = (
        "donor/animal-aware pseudobulk aggregation + a shrunken "
        "log-fold-change-style estimator (NOT independent-cell Wilcoxon, NOT "
        "an arbitrary signed -log10(p) substitute) -- not yet implemented")

    def run(self, *, comparison: ComparisonSpec, sample_structure: SampleStructure,
            group_values: list, donor_values: list, sample_ids: list | None = None,
            estimator: str = INTENDED_ESTIMATOR, **params) -> ExpressionDEResult:
        problems = comparison.validate() + sample_structure.validate()
        if comparison.group_column and sample_structure.group_column \
                and comparison.group_column != sample_structure.group_column:
            problems.append(
                f"comparison.group_column={comparison.group_column!r} != "
                f"sample_structure.group_column={sample_structure.group_column!r} -- "
                "the compared group and the declared sample-structure group "
                "column must be the SAME column")
        if not estimator:
            problems.append(
                "estimator must be stated explicitly, even to name a "
                "placeholder -- never left blank")
        if sample_ids is not None:
            if len(sample_ids) != len(group_values):
                problems.append(
                    f"sample_ids ({len(sample_ids)}) and group_values "
                    f"({len(group_values)}) length mismatch")
            else:
                bad_sid = [i for i, s in enumerate(sample_ids)
                          if not isinstance(s, str) or not s]
                if bad_sid:
                    problems.append(
                        f"sample_ids has {len(bad_sid)} non-string/empty value(s) at "
                        f"index/indices {bad_sid[:5]}" + (" (truncated)" if len(bad_sid) > 5 else ""))
                elif len(set(sample_ids)) != len(sample_ids):
                    problems.append("sample_ids contains duplicate value(s)")
        if not problems:
            problems = check_design_estimability(
                group_values, donor_values, comparison, sample_structure.design_type)

        n_samples, group_values_sha256, donor_values_sha256 = _sample_label_identity(
            group_values, donor_values)
        input_mapping = _build_input_mapping(sample_ids, group_values, donor_values)
        # Per `review:20261014_114357` R3: `input_lengths`
        # preserves the RAW lengths even when they disagree (never
        # fabricated, dropped, or padded to agree); `input_mapping_
        # buildable` records SEPARATELY whether a row-aligned mapping
        # could be built from them -- so an empty `input_mapping` on a
        # genuinely `'failed'` length-mismatch record is never confused
        # with "0 samples" on a clean run.
        input_lengths = dict(
            group_values=len(group_values), donor_values=len(donor_values),
            sample_ids=(len(sample_ids) if sample_ids is not None else None))
        input_mapping_buildable = _input_mapping_buildable(sample_ids, group_values, donor_values)
        common = dict(
            comparison=dataclasses.asdict(comparison),
            sample_structure=dataclasses.asdict(sample_structure),
            estimator=estimator, n_samples=n_samples,
            group_values_sha256=group_values_sha256, donor_values_sha256=donor_values_sha256,
            input_lengths=input_lengths, input_mapping_buildable=input_mapping_buildable,
            sample_identity_available=sample_ids is not None,
            input_mapping=input_mapping, input_mapping_sha256=_input_mapping_sha256(input_mapping),
            backend="not_implemented", backend_version="n/a",
            params=dict(params), source_code_revision_or_hash=_source_revision(),
            dependency_versions={"python": platform.python_version()})

        if problems:
            return ExpressionDEResult(status="failed", reason="; ".join(problems), **common)

        reason = ("ExpressionDERankingProvider.run is a contract/validation "
                 "stub this round -- inputs passed structural and "
                 "design-estimability validation, but no real pseudobulk/"
                 "shrunken-LFC computation was performed. No fabricated "
                 "ranking is returned.")
        if sample_structure.covariate_columns:
            reason += (
                f" NOTE: sample_structure.covariate_columns="
                f"{sample_structure.covariate_columns!r} is DECLARED, but this "
                "round's contract takes no covariate VALUES to check estimability "
                "against -- the covariate-adjusted design is UNVERIFIED; only the "
                "group/donor structure above was checked.")

        return ExpressionDEResult(status="not_run", reason=reason, **common)

    def run_pydeseq2_from_prepared_design(self, *, prepared_design: dict, counts,
                                          gene_ranking_definition: str = "") -> ExpressionDEResult:
        """The REAL backend call path, per `review:20261015_104226`
        item 4. Consumes an already-validated `PreparedDEDesign` manifest
        (`de_design_providers.read_prepared_de_design`) and its hash-bound
        counts (`de_design_providers.load_prepared_de_counts`) -- never
        re-derives sample order, gene order, the design matrix, or the
        contrast; all four are taken VERBATIM from `prepared_design`. The
        caller passes the already-loaded manifest dict and counts array.
        `de_design_providers` is imported LAZILY inside this method (not at
        this module's top level) -- that module already imports FROM this
        one at ITS top level, so a top-level import back here would create
        a cycle; by the time this method actually RUNS, both modules'
        top-level code has already finished, so a deferred import is safe.

        Per `review:20261015_142903` R1: a caller-supplied
        `prepared_design` dict and `counts` array are NOT trusted merely
        because they have the right keys/shape -- this method re-runs the
        SAME semantic/source-binding validation every other consumer of a
        prepared design goes through (`de_design_providers.
        validate_prepared_de_design`, the identical check
        `read_prepared_de_design` already applies; and
        `de_design_providers._source_binding_problems`, the identical check
        `load_prepared_de_counts` already applies), rather than inventing a
        second, parallel validation path. This closes a probe that found a
        same-shape counts-content change, a fractional count value, a
        flipped contrast sign, and an all-zero design matrix all silently
        reaching the backend call unrejected -- each of those is now caught
        by the SAME validator/binding-check functions the write/load path
        already relies on, never a new bespoke check.

        Runs this stream's OWN, already-tested structural/design-
        estimability validation first (`self.run`, UNCHANGED) on the
        manifest's own `group_values`/`donor_ids`/`sample_ids` -- a
        structurally invalid design is refused by the SAME checks the
        contract-stub path already applies; this method adds no second,
        parallel validation path. Only once `self.run` returns
        `status='not_run'` (inputs passed validation, stub placeholder) does
        this method attempt the real `pydeseq2` call.

        `pydeseq2` is NOT installed in `302006_env_clean` as of this round
        (no isolated `302006_de_backend_v1` environment exists yet -- see
        `DE_BACKEND_EXECUTION_PACKET_v1.md`). On `ImportError`, this returns
        `status='not_run'` with an explicit reason naming the missing
        dependency -- this is the HONEST not_run reason for the REAL
        backend path, distinct from (and replacing) the earlier generic
        contract-stub placeholder reason. It never returns a fabricated
        `status='complete'` result merely because inputs validated.

        UNVERIFIED this round: the actual `DeseqDataSet`/`DeseqStats`/
        `lfc_shrink` call sequence below has not been exercised against a
        real `pydeseq2` install (none exists in any environment this round
        has access to) -- it is reviewed code, not executed code, until an
        isolated environment with `pydeseq2` is created and this path is
        re-run there.
        """
        if prepared_design.get("status") != "prepared":
            raise ValueError(
                "run_pydeseq2_from_prepared_design requires a status='prepared' "
                f"PreparedDEDesign manifest, got {prepared_design.get('status')!r}")

        import de_design_providers as DD

        comparison = ComparisonSpec(**prepared_design["comparison"])
        sample_structure_fields = dict(prepared_design["sample_structure"])
        # The JSON round-trip through `write_prepared_de_design`/
        # `read_prepared_de_design` turns this tuple field into a list;
        # restore the tuple type `SampleStructure`'s own field declares.
        sample_structure_fields["covariate_columns"] = tuple(
            sample_structure_fields.get("covariate_columns") or ())
        sample_structure = SampleStructure(**sample_structure_fields)
        sample_ids = list(prepared_design["sample_ids"])
        gene_ids = list(prepared_design["gene_ids"])

        stub = self.run(
            comparison=comparison, sample_structure=sample_structure,
            group_values=list(prepared_design["group_values"]),
            donor_values=list(prepared_design["donor_ids"]), sample_ids=sample_ids,
            estimator="pydeseq2 Wald test + lfc_shrink, consuming a PreparedDEDesign "
                     "verbatim (DeseqDataSet/DeseqStats/lfc_shrink) -- not yet executed "
                     "against a real pydeseq2 install this round")
        if stub.status != "not_run":
            # The SAME structural/design-estimability problem self.run already
            # detected from these group/donor values -- never attempt a real
            # backend call on an input run() itself refused.
            return stub

        # Per review:20261015_142903 R1: re-run the SAME semantic
        # recompute-and-compare `read_prepared_de_design` already applies --
        # rebuilds design_matrix/column_names/contrast_vector/contrast_name/
        # rank/residual_degrees_of_freedom from the payload's OWN sample_ids/
        # donor_ids/group_values/comparison/sample_structure and refuses on
        # any mismatch. This is what actually catches a flipped contrast
        # sign or an all-zero design matrix -- the estimability check alone
        # (already run above via `self.run`/`check_design_estimability`) is
        # VACUOUS for any full-rank design, per this module's own documented
        # limitation (de_design_providers.py's module docstring).
        manifest_problems = DD.validate_prepared_de_design(prepared_design)
        if manifest_problems:
            return dataclasses.replace(
                stub, status="failed",
                reason="prepared_design failed re-validation at the backend "
                      f"consumption entry point: {manifest_problems}")

        expected_shape = (len(sample_ids), len(gene_ids))
        counts_shape = tuple(getattr(counts, "shape", ()))
        if counts_shape != expected_shape:
            return dataclasses.replace(
                stub, status="failed",
                reason=f"counts shape {counts_shape} != expected {expected_shape} from "
                      "prepared_design's own sample_ids/gene_ids -- refusing to call "
                      "the backend on a mismatched counts array")

        # Per review:20261015_142903 R1: reject non-integer/negative/
        # non-finite counts BEFORE the int64 cast below -- `np.asarray(...,
        # dtype=np.int64)` TRUNCATES a fractional value (50.25 -> 50) rather
        # than raising, so this check must run first, not be implied by it.
        numeric_problems = DD._nonneg_finite_integer_counts_problems(counts)
        if numeric_problems:
            return dataclasses.replace(
                stub, status="failed",
                reason=f"counts failed validation: {numeric_problems}")

        # Per review:20261015_142903 R1: bind the passed-in `counts`
        # to the LIVE upstream pseudobulk source -- the SAME check
        # `load_prepared_de_counts` already performs. This is strictly
        # stronger than comparing against `prepared_design`'s own recorded
        # `counts_content_sha256`: a caller who tampered with `counts` AND
        # consistently updated that recorded hash to match would pass a
        # stored-hash comparison but still fails here, because the expected
        # counts submatrix is re-derived independently from the live
        # upstream manifest every time, never from anything `prepared_design`
        # itself claims about its own counts.
        binding_problems = DD._source_binding_problems(
            prepared_design["pseudobulk_identity"], prepared_design["dataset_identity"],
            prepared_design["comparison"], prepared_design["sample_structure"],
            sample_ids, list(prepared_design["donor_ids"]), list(prepared_design["group_values"]),
            gene_ids, counts)
        if binding_problems:
            return dataclasses.replace(
                stub, status="failed",
                reason=f"counts failed source-binding re-verification: {binding_problems}")

        try:
            import pydeseq2
            from pydeseq2.dds import DeseqDataSet
            from pydeseq2.ds import DeseqStats
        except ImportError as e:
            return dataclasses.replace(
                stub, backend="pydeseq2", backend_version="unavailable",
                reason=f"pydeseq2 is not importable in this interpreter ({e!r}) -- "
                      "no isolated DE-backend environment exists yet (see "
                      "DE_BACKEND_EXECUTION_PACKET_v1.md). This is the honest "
                      "not_run reason for the REAL backend call path.")

        import numpy as np
        import pandas as pd

        counts_df = pd.DataFrame(
            np.asarray(counts, dtype=np.int64), index=sample_ids, columns=gene_ids)
        metadata_df = pd.DataFrame(
            {"condition": list(prepared_design["group_values"]),
            "donor": list(prepared_design["donor_ids"])},
            index=sample_ids)
        # Per DE_BACKEND_EXECUTION_PACKET_v1.md section 5 (verified against
        # pydeseq2 v0.5.4's own source in a prior round, NOT re-verified by
        # this method): DeseqDataSet accepts the design as a DataFrame
        # directly, so PreparedDEDesign's own design_matrix/column_names
        # feed it with no re-derivation through a formula engine.
        design_df = pd.DataFrame(
            [list(row) for row in prepared_design["design_matrix"]],
            index=sample_ids, columns=list(prepared_design["column_names"]))
        # Per review:20261015_142903 R2: pydeseq2 v0.5.4's own
        # `DeseqStats` constructor branches on the Python TYPE of `contrast`
        # -- a numpy ndarray is treated as a NUMERIC contrast, while a plain
        # list (even one holding only floats) is treated as a 3-item
        # CATEGORICAL contrast spec and indexed at position 2, raising
        # IndexError for this project's 2-column design. Verified directly
        # against pydeseq2 v0.5.4's own source (DE_BACKEND_EXECUTION_PACKET_
        # v1.md section 5); a plain list here is a real crash, not a style
        # preference.
        contrast_vector = np.asarray(prepared_design["contrast_vector"], dtype=float)
        contrast_name = prepared_design["contrast_name"]

        # Per review:20261015_142903 R4: bound pydeseq2's OWN
        # joblib worker count to this job's actual SLURM allocation --
        # `OMP_NUM_THREADS`/`OPENBLAS_NUM_THREADS` alone do not constrain
        # joblib's separate process pool, which otherwise defaults to every
        # core on the node per pydeseq2's own `DefaultInference`.
        n_cpus = _bounded_n_cpus()

        try:
            dds = DeseqDataSet(
                counts=counts_df, metadata=metadata_df, design=design_df, n_cpus=n_cpus)
            dds.deseq2()
            stats = DeseqStats(dds, contrast=contrast_vector, n_cpus=n_cpus)
            stats.summary()
            raw_results = stats.results_df.copy()
            stats.lfc_shrink(coeff=contrast_name)
            shrunk_results = stats.results_df.copy()
        except Exception as e:
            # pydeseq2's own fit/convergence failure -- never upgraded to a
            # fabricated 'complete' result.
            return dataclasses.replace(
                stub, backend="pydeseq2",
                backend_version=getattr(pydeseq2, "__version__", "unknown"),
                status="failed", reason=f"pydeseq2 raised during fit/stats: {e!r}")

        dropped_before_fit = sorted(set(gene_ids) - set(raw_results.index))
        usable_genes, scores, nonfinite_genes = [], [], []
        for g in gene_ids:
            if g not in raw_results.index:
                continue
            raw_lfc = raw_results.loc[g, "log2FoldChange"]
            raw_stat = raw_results.loc[g, "stat"]
            shrunk_lfc = shrunk_results.loc[g, "log2FoldChange"]
            if all(isinstance(v, (int, float)) and math.isfinite(v)
                  for v in (raw_lfc, raw_stat, shrunk_lfc)):
                usable_genes.append(g)
                scores.append(float(shrunk_lfc))
            else:
                nonfinite_genes.append(g)
        genes_backend_dropped = dropped_before_fit + nonfinite_genes

        ranking = GeneRanking(
            species=prepared_design.get("dataset_identity", {}).get("species", "unknown"),
            gene_id_namespace=prepared_design.get("dataset_identity", {}).get(
                "gene_id_namespace", "unknown"),
            ranking_definition=(
                gene_ranking_definition or
                f"expression_de: pydeseq2 shrunken log2FoldChange, contrast={contrast_name}"),
            genes=tuple(usable_genes), scores=tuple(scores))

        params = dict(
            stub.params,
            raw_wald_stat={g: float(raw_results.loc[g, "stat"]) for g in usable_genes},
            raw_lfc={g: float(raw_results.loc[g, "log2FoldChange"]) for g in usable_genes},
            shrunken_lfc={g: float(shrunk_results.loc[g, "log2FoldChange"]) for g in usable_genes},
            shrinkage_method="pydeseq2.DeseqStats.lfc_shrink",
            genes_backend_dropped=genes_backend_dropped,
            contrast_name=contrast_name, contrast_vector=contrast_vector.tolist(),
            n_cpus_bound_to_slurm_allocation=n_cpus,
            prepared_design_counts_content_sha256=prepared_design.get(
                "counts_content_sha256", ""))

        return dataclasses.replace(
            stub, status="complete",
            reason="real pydeseq2 Wald test + lfc_shrink executed, consuming the "
                  "PreparedDEDesign's own sample/gene order, design matrix, and "
                  "contrast verbatim",
            ranking=ranking, backend="pydeseq2",
            backend_version=getattr(pydeseq2, "__version__", "unknown"),
            params=params,
            dependency_versions=dict(
                stub.dependency_versions,
                pydeseq2=getattr(pydeseq2, "__version__", "unknown")))


# ===========================================================================
# Model-masking evidence stream
# ===========================================================================

@dataclass(frozen=True)
class MaskingScoreSpec:
    """Binds `Delta(g) = s_full - s_masked(g)` to a FIXED scoring function --
    never re-fit per gene, never re-decided per call.
    """
    #: What s(.) IS -- e.g. "state-attribution direction score" or
    #: "disease-contrast direction score" (Cindy's decision, pending); never
    #: left to be inferred from the ranking direction after the fact.
    score_definition: str
    #: What "masking gene g" does to the model's INPUT -- e.g. "set gene g's
    #: expression value to 0 before tokenization" -- never a cell-selection
    #: operation; this spec has no field that could mean that.
    masking_operation: str
    #: What value/token a masked gene's input position takes (e.g. a
    #: model-specific [MASK] token id, or a literal 0 in expression space).
    reference_or_token_semantics: str
    #: e.g. "all genes in the model's input vocabulary present in this
    #: dataset" -- states exactly which genes are even eligible for g.
    perturbable_gene_set_definition: str
    #: How per-cell Delta(g) rolls up to a donor/animal-level score (e.g.
    #: "mean over cells within donor, then mean over donors" -- never
    #: implied by whatever the backend's default reduction happens to be).
    donor_or_animal_aggregation_rule: str

    def validate(self) -> list[str]:
        problems = []
        for f in ("score_definition", "masking_operation", "reference_or_token_semantics",
                  "perturbable_gene_set_definition", "donor_or_animal_aggregation_rule"):
            v = getattr(self, f)
            # Per `review:20261014_100436` F3: a truthiness check
            # alone passes `score_definition=True` (a non-empty bool is
            # truthy) -- this field must be the actual descriptive string,
            # never any other truthy value.
            if not isinstance(v, str) or not v:
                problems.append(f"{f} must be a non-empty string, got {v!r} -- state it explicitly")
        return problems


def _gene_set_identity(perturbable_genes) -> tuple:
    """Content identity of the ACTUAL `perturbable_genes` this run was
    called with -- per `review:20261013_210131` F4, so e.g.
    `perturbable_genes=("A",)` vs. `("B",)` serialize to two DIFFERENT
    `ModelMaskingResult` records, even while no real backend exists yet.
    Factored out as its own function so a mutation that stops binding this
    identity can be demonstrated directly (see
    `test_mutation_kill_model_masking_perturbable_genes_identity`). Per
    `review:20261013_222450` R3: uses `_safe_json_identity`, never
    `default=str` -- a gene "ID" that is not already a plain string is a
    structural validation failure (see `ModelMaskingRankingProvider.run`),
    never silently stringified into this hash.
    """
    genes = tuple(perturbable_genes)
    return genes, _safe_json_identity(list(genes)) or ""


@dataclass(frozen=True)
class ModelMaskingResult:
    status: str  # "complete" | "failed" | "not_run"
    reason: str = ""
    evidence_type: str = "model_masking"
    ranking: object = None  # GeneRanking, only when status == "complete"
    masking_score: dict = field(default_factory=dict)
    #: `model_registry.ModelRegistryEntry.to_dict()`-shaped dict, or the
    #: literal string "model_undecided" -- mirrors `artifact_schema
    #: .build_manifest`'s own contract so this cannot silently diverge from
    #: how the rest of the framework states "no model pinned yet".
    model_identity: object = "model_undecided"
    #: See `_gene_set_identity` -- the ACTUAL, ORDERED gene set this run was
    #: called with, bound alongside its content hash.
    perturbable_genes: tuple = ()
    perturbable_genes_sha256: str = ""
    backend: str = "unknown"
    backend_version: str = "unknown"
    params: dict = field(default_factory=dict)
    source_code_revision_or_hash: str = ""
    dependency_versions: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        d = dict(self.__dict__)
        d["ranking"] = _ranking_to_dict(self.ranking)
        d["perturbable_genes"] = list(self.perturbable_genes)
        return d


class ModelMaskingRankingProvider:
    """Model evidence stream. See module docstring. Implements validation/
    serialization only this round -- `run` never performs a real full/
    masked forward pass.
    """

    def run(self, *, masking_score: MaskingScoreSpec, model_identity_or_pending,
            perturbable_genes: tuple, **params) -> ModelMaskingResult:
        problems = masking_score.validate()
        if model_identity_or_pending != "model_undecided" \
                and not isinstance(model_identity_or_pending, dict):
            problems.append(
                "model_identity_or_pending must be the literal string "
                "'model_undecided' or a dict describing a resolved model "
                "registry entry")
        if not perturbable_genes:
            problems.append(
                "perturbable_genes is empty -- state which genes are actually "
                "perturbable under this masking operation")
        else:
            bad_gene = [i for i, g in enumerate(perturbable_genes)
                       if not isinstance(g, str) or not g]
            if bad_gene:
                problems.append(
                    f"perturbable_genes has {len(bad_gene)} non-string/empty value(s) "
                    f"at index/indices {bad_gene[:5]}" + (" (truncated)" if len(bad_gene) > 5 else ""))
            elif len(set(perturbable_genes)) != len(perturbable_genes):
                problems.append("perturbable_genes contains duplicate gene ID(s)")

        gene_set, gene_set_sha256 = _gene_set_identity(perturbable_genes)
        common = dict(
            masking_score=dataclasses.asdict(masking_score),
            model_identity=model_identity_or_pending,
            perturbable_genes=gene_set, perturbable_genes_sha256=gene_set_sha256,
            params=dict(params),
            backend="not_implemented", backend_version="n/a",
            source_code_revision_or_hash=_source_revision(),
            dependency_versions={"python": platform.python_version()})

        if problems:
            return ModelMaskingResult(status="failed", reason="; ".join(problems), **common)

        return ModelMaskingResult(
            status="not_run",
            reason="ModelMaskingRankingProvider.run is a contract/validation "
                   "stub this round -- inputs passed structural validation, "
                   "but no real full/masked forward pass was performed (no "
                   "main model checkpoint is pinned yet). No fabricated "
                   "Delta(g) ranking is returned.",
            **common)


def make_evidence_ranking(evidence_type: str, ranking: GeneRanking) -> GeneRanking:
    """Proves the "shared `GeneRanking` output" contract: whichever evidence
    stream eventually produces a real ranking, it is bound through this one
    choke point, which requires `ranking_definition` to actually NAME the
    evidence stream it came from -- a bare, evidence-agnostic
    `ranking_definition` string can never reach the GSEA chain through this
    function. Does not itself compute anything; a real backend calls this
    once it has a genuine `GeneRanking` to bind.

    Per `review:20261013_210131` F4: this is a LINT-LEVEL tripwire
    on `ranking_definition`'s text, not a structured source/provenance
    check -- `GeneRanking` (shared with the already-tested GSEA chain) is
    deliberately left unchanged this round, so it has no structured
    `evidence_type`/comparison/score-identity fields of its own. The real,
    structured input identity for each stream lives on `ExpressionDEResult`
    (`group_values_sha256`/`donor_values_sha256`) and `ModelMaskingResult`
    (`perturbable_genes`/`perturbable_genes_sha256`) -- this substring check
    must not be read as having verified those.
    """
    if evidence_type not in ("expression_de", "model_masking"):
        raise ValueError(f"evidence_type={evidence_type!r} not in "
                         "('expression_de', 'model_masking')")
    problems = ranking.validate()
    if problems:
        raise ValueError(f"cannot bind an invalid ranking: {problems}")
    if evidence_type not in ranking.ranking_definition:
        raise ValueError(
            f"ranking_definition {ranking.ranking_definition!r} does not name "
            f"its own evidence_type {evidence_type!r} -- a ranking must be "
            "traceable to which evidence stream produced it from the "
            "ranking_definition text alone, not only from the result "
            "wrapper's own evidence_type field")
    return ranking


# ===========================================================================
# Evidence envelope: strict write/read contract for both ranking results
# ===========================================================================
#
# Per `review:20261013_222450` R3/A.3: a unified validator shared
# by construction, write, and read-back -- a tampered on-disk file (mapping,
# hash, or identity field changed) is rejected on read, never silently
# accepted because only the writer checked. Mirrors
# `embedding_providers.validate_embedding_result`/`write_embedding_result`/
# `read_embedding_result`'s atomic-write + dual-validation pattern. A
# `not_run`/`failed` record can never carry a `ranking` that would make it
# LOOK like a completed computation.


def _ranking_payload_problems(ranking) -> list[str]:
    """Per `review:20261014_100436` F3: a prior round's reader
    accepted `status='complete', ranking={}` -- `ranking is not None` was
    the only check, and an EMPTY dict is not `None`. Reuses `GeneRanking
    .validate()` (the same check the GSEA chain already trusts) rather than
    reimplementing its field rules here; never re-validates GeneRanking
    itself.
    """
    if not isinstance(ranking, dict):
        return [f"ranking must be a dict, got {type(ranking).__name__}"]
    required = ("species", "gene_id_namespace", "ranking_definition", "genes", "scores")
    missing = [k for k in required if k not in ranking]
    if missing:
        return [f"ranking is missing required key(s): {missing}"]
    try:
        gr = GeneRanking(
            species=ranking["species"], gene_id_namespace=ranking["gene_id_namespace"],
            ranking_definition=ranking["ranking_definition"],
            genes=tuple(ranking["genes"]) if isinstance(ranking["genes"], list) else ranking["genes"],
            scores=tuple(ranking["scores"]) if isinstance(ranking["scores"], list) else ranking["scores"])
    except TypeError as e:
        return [f"ranking fields have the wrong type to construct a GeneRanking: {e}"]
    return gr.validate()


def _validate_evidence_result_common(payload: dict, schema: str,
                                     expected_evidence_type: str) -> list[str]:
    """Shared by `validate_expression_de_result`/`validate_model_masking_
    result` and both writer/reader pairs.
    """
    if not isinstance(payload, dict):
        return [f"result must be a dict, got {type(payload).__name__}"]
    problems = []
    if payload.get("schema") != schema:
        problems.append(f"schema {payload.get('schema')!r} != expected {schema!r}")
    if payload.get("evidence_type") != expected_evidence_type:
        problems.append(
            f"evidence_type {payload.get('evidence_type')!r} != expected "
            f"{expected_evidence_type!r}")
    status = payload.get("status")
    if status not in _VALID_RESULT_STATUSES:
        problems.append(f"status {status!r} not in {_VALID_RESULT_STATUSES}")
    problems.extend(find_nonfinite_floats(payload))
    if status != "complete":
        if not payload.get("reason"):
            problems.append(f"status={status!r} requires a non-empty 'reason'")
        if payload.get("ranking") is not None:
            problems.append(
                f"status={status!r} must not carry a non-null 'ranking' -- a "
                "result that is not status='complete' must never look like a "
                "computed one")
    else:
        if payload.get("ranking") is None:
            problems.append("status='complete' requires a non-null 'ranking'")
        else:
            problems.extend(_ranking_payload_problems(payload.get("ranking")))
        # Per `review:20261014_100436` F3: `backend='not_implemented'`/
        # `backend_version='n/a'` are this module's OWN placeholders for "no
        # real backend yet" (see `ExpressionDERankingProvider.run`/
        # `ModelMaskingRankingProvider.run`) -- a record claiming
        # status='complete' can never carry either placeholder, since no
        # backend this round can have actually produced a real computation.
        backend = payload.get("backend")
        if not isinstance(backend, str) or not backend or backend == "not_implemented":
            problems.append(
                "status='complete' requires a real 'backend' name, not the "
                "'not_implemented' placeholder -- no backend exists yet this "
                "round, so no 'complete' result can be genuine")
        backend_version = payload.get("backend_version")
        if not isinstance(backend_version, str) or not backend_version or backend_version == "n/a":
            problems.append(
                "status='complete' requires a real 'backend_version', not the "
                "'n/a' placeholder")
    return problems


def _require_nonempty_str_fields(d: dict, fields: tuple, label: str) -> list[str]:
    return [f"{label}[{f!r}] must be a non-empty string"
           for f in fields if not isinstance(d.get(f), str) or not d.get(f)]


def _unbuildable_status_problems(expected_buildable: bool, status, input_lengths: dict) -> list[str]:
    """Per `review:20261014_151747` P2: `expected_buildable=False`
    (recomputed from `input_lengths` -- a genuine length disagreement) ties
    directly to `status`. This module defines `'not_run'` as "inputs
    PASSED validation, no backend ran yet" (see `ExpressionDERankingProvider
    .run`) and `'complete'` as a real computation having run; an unbuildable
    input never passed validation, so it can only ever be `'failed'`.
    Without this check, a genuinely `'failed'` length-mismatch record could
    have ONLY its `status` field rewritten to `'not_run'`/`'complete'` on
    disk and still be accepted, since the mapping-emptiness/hash/`n_samples`
    cross-checks all still hold for an unchanged, merely status-relabelled
    record. Factored out so the guard itself can be mutation-killed
    (see `test_mutation_kill_validate_expression_de_result_unbuildable_status`).
    """
    if not expected_buildable and status in ("not_run", "complete"):
        return [
            f"input_lengths={input_lengths!r} cannot build a real "
            f"input_mapping (input_mapping_buildable=False), so status "
            f"must be 'failed', not {status!r} -- an unbuildable input "
            "never passed validation"]
    return []


def validate_expression_de_result(payload: dict) -> list[str]:
    """Per `review:20261013_222450` R3 and
    `review:20261014_100436` F3: checks field types/enums, the
    `reason` requirement, `ranking`/`status` consistency (shared via
    `_validate_evidence_result_common`), AND re-derives `input_mapping_
    sha256` from `input_mapping` itself -- a caller that tampers with
    either the mapping or the hash alone is caught here, matching the
    review's explicit "篡改映射、哈希或身份字段时 reader 必须拒绝" requirement.

    Per F3: a `status='failed'` record's `input_mapping` is a DIAGNOSTIC
    record of the ACTUAL (possibly invalid -- e.g. a `None` donor_value)
    input that caused the failure, never a validated mapping -- only the
    KEYS are required to be present for `status='failed'`; the per-row
    VALUE checks (group_value/donor_value/sample_id must be non-empty
    strings) apply only when `status` is `'complete'` or `'not_run'`,
    since both of those states are contracts that the raw values WERE
    already checked valid (`check_design_estimability` cannot return
    `'not_run'` with an invalid group/donor value -- see
    `_validate_group_donor_values`). This is the fix for the prior round's
    gap where a genuinely `failed` record (produced by the REAL provider
    path for an invalid `donor_value=None` input) could not be written at
    all, because the writer applied the SAME strict checks to every status.
    `comparison`/`sample_structure` are DECLARED CONFIGURATION, not raw
    per-row data -- their nested field/enum checks apply whenever `status`
    is `'complete'` or `'not_run'`, catching a tampered on-disk
    `sample_structure.design_type` even on a `'not_run'` record.

    Per `review:20261014_114357` R2/R3: two further gaps. (1)
    `n_samples`'s own TYPE was never checked -- `isinstance(n_samples, int)
    and ... != ...` had no `else` branch, so `n_samples='two'`/`None` just
    skipped the consistency check silently instead of being rejected; now
    checked unconditionally first. (2) that same consistency check
    (`n_samples == len(input_mapping)`) used to run UNCONDITIONALLY, which
    made a genuinely `'failed'` length-mismatch record (`group_values`/
    `donor_values` of different lengths, so `input_mapping` is correctly
    forced empty while `n_samples=len(group_values)` is non-zero)
    unwritable -- the same "illegitimate-accept vs. genuine-reject"
    asymmetry F3 fixed for `input_mapping`/`masking_score` values, but for
    a LENGTH relationship instead of a value. Fixed by recording
    `input_lengths` (the raw lengths, always) and `input_mapping_buildable`
    (whether they agreed, always) as their own fields, recomputing both
    from the payload and cross-checking before applying the length
    comparison -- which now only runs when a mapping was actually
    buildable; when it was not, only "is `input_mapping` empty" is
    required, never "does its length match n_samples".

    Per `review:20261014_151747` P2: additionally rejects
    `input_mapping_buildable=False` (recomputed from `input_lengths`)
    combined with `status in ('not_run', 'complete')` -- an unbuildable
    input can only ever be `status='failed'`; see the module docstring.
    """
    problems = _validate_evidence_result_common(
        payload, EXPRESSION_DE_RESULT_SCHEMA, "expression_de")
    if not isinstance(payload, dict):
        return problems
    status = payload.get("status")
    require_valid_values = status in ("complete", "not_run")

    mapping = payload.get("input_mapping")
    if not isinstance(mapping, list):
        problems.append("input_mapping must be a list")
        mapping = []
    seen_sample_ids: set = set()
    for i, row in enumerate(mapping):
        if not isinstance(row, dict):
            problems.append(f"input_mapping[{i}] must be a dict")
            continue
        if row.get("row_index") != i:
            problems.append(f"input_mapping[{i}]['row_index']={row.get('row_index')!r} != {i}")
        if "group_value" not in row or "donor_value" not in row:
            problems.append(
                f"input_mapping[{i}] must carry 'group_value' and 'donor_value' "
                "keys even when diagnosing an invalid input")
        elif require_valid_values:
            if not isinstance(row.get("group_value"), str) or not row.get("group_value"):
                problems.append(f"input_mapping[{i}]['group_value'] must be a non-empty string")
            if not isinstance(row.get("donor_value"), str) or not row.get("donor_value"):
                problems.append(f"input_mapping[{i}]['donor_value'] must be a non-empty string")
        if payload.get("sample_identity_available"):
            if "sample_id" not in row:
                problems.append(
                    f"input_mapping[{i}] must carry a 'sample_id' key when "
                    "sample_identity_available")
            elif require_valid_values:
                sid = row.get("sample_id")
                if not isinstance(sid, str) or not sid:
                    problems.append(f"input_mapping[{i}]['sample_id'] must be a non-empty string")
                elif sid in seen_sample_ids:
                    problems.append(f"input_mapping has duplicate sample_id {sid!r}")
                else:
                    seen_sample_ids.add(sid)
    if isinstance(payload.get("input_mapping"), list):
        expected_hash = _input_mapping_sha256(tuple(payload["input_mapping"]))
        if payload.get("input_mapping_sha256") != expected_hash:
            problems.append(
                "input_mapping_sha256 does not match a recomputed hash of "
                "input_mapping -- tampered or stale")

    # Per `review:20261014_114357` R2/R3: `n_samples`'s TYPE is
    # checked unconditionally, FIRST -- `n_samples='two'`/`None` must never
    # silently skip the consistency check below the way an `isinstance`
    # guard with no `else` branch did previously. `input_lengths` is
    # checked the same way, then `n_samples` is cross-checked against it
    # (ALWAYS true by construction in `run()`, regardless of status or
    # buildability -- `n_samples` is defined as `len(group_values)`).
    n_samples = payload.get("n_samples")
    if not isinstance(n_samples, int) or isinstance(n_samples, bool):
        problems.append(f"n_samples must be an int, got {n_samples!r}")

    il = payload.get("input_lengths")
    il_ok = (isinstance(il, dict)
            and isinstance(il.get("group_values"), int) and not isinstance(il.get("group_values"), bool)
            and isinstance(il.get("donor_values"), int) and not isinstance(il.get("donor_values"), bool)
            and (il.get("sample_ids") is None
                 or (isinstance(il.get("sample_ids"), int) and not isinstance(il.get("sample_ids"), bool))))
    if not il_ok:
        problems.append(
            "input_lengths must be a dict with int 'group_values'/'donor_values' and "
            f"an int-or-None 'sample_ids', got {il!r}")
    elif isinstance(n_samples, int) and not isinstance(n_samples, bool) \
            and n_samples != il["group_values"]:
        problems.append(
            f"n_samples={n_samples} does not match input_lengths['group_values']="
            f"{il['group_values']}")

    # Per `review:20261014_114357` R3: whether `input_mapping`
    # SHOULD be a real row-aligned mapping is recomputed from `input_
    # lengths` itself (never trusted from the payload's own `input_mapping_
    # buildable` claim alone) and compared against it -- a tampered
    # `input_mapping_buildable` is caught the same way a tampered hash is.
    # The length-vs-`n_samples` consistency check that used to run
    # unconditionally now applies ONLY when a mapping was actually
    # buildable; when it was NOT (the inputs' lengths genuinely disagreed),
    # `input_mapping` must be the forced-empty diagnostic placeholder, and
    # comparing its length to `n_samples` would wrongly demand 0 raw
    # samples for what is actually a length-mismatch failure.
    buildable = payload.get("input_mapping_buildable")
    if not isinstance(buildable, bool):
        problems.append(f"input_mapping_buildable must be a bool, got {buildable!r}")
    elif il_ok:
        expected_buildable = (il["group_values"] == il["donor_values"]
                              and (il["sample_ids"] is None or il["sample_ids"] == il["group_values"]))
        if buildable != expected_buildable:
            problems.append(
                f"input_mapping_buildable={buildable} does not match a recomputed "
                f"value of {expected_buildable} from input_lengths={il!r}")
        elif isinstance(payload.get("input_mapping"), list):
            mapping_list = payload["input_mapping"]
            if expected_buildable:
                if isinstance(n_samples, int) and not isinstance(n_samples, bool) \
                        and n_samples != len(mapping_list):
                    problems.append(
                        f"n_samples={n_samples} does not match "
                        f"len(input_mapping)={len(mapping_list)}")
            elif mapping_list:
                problems.append(
                    "input_mapping_buildable=False but input_mapping is non-empty -- "
                    "an unbuildable (length-mismatched) input must be recorded as an "
                    "empty diagnostic list, never a partial/fabricated mapping")
        problems.extend(_unbuildable_status_problems(expected_buildable, status, il))

    for key in ("comparison", "sample_structure"):
        if not isinstance(payload.get(key), dict) or not payload.get(key):
            problems.append(f"{key!r} must be a non-empty dict")

    if require_valid_values:
        comparison = payload.get("comparison")
        if isinstance(comparison, dict):
            problems.extend(_require_nonempty_str_fields(
                comparison, ("group_column", "reference_level", "comparison_level",
                            "direction_definition"), "comparison"))
        sample_structure = payload.get("sample_structure")
        if isinstance(sample_structure, dict):
            problems.extend(_require_nonempty_str_fields(
                sample_structure, ("donor_or_animal_id_column", "sample_id_column",
                                   "repeated_sampling_or_pairing", "group_column",
                                   "design_type"), "sample_structure"))
            design_type = sample_structure.get("design_type")
            if design_type and design_type not in _VALID_DESIGN_TYPES:
                problems.append(
                    f"sample_structure['design_type']={design_type!r} not in "
                    f"{_VALID_DESIGN_TYPES!r}")
    return problems


def write_expression_de_result(out_dir: str, result: ExpressionDEResult) -> ExpressionDEResult:
    """Writes the real input->provider->persisted-output boundary (not just
    an in-memory dataclass) -- `validate_expression_de_result` runs here AND
    in `read_expression_de_result`.
    """
    if os.path.exists(out_dir) and os.listdir(out_dir):
        raise SystemExit(f"{out_dir} must be new or empty")
    os.makedirs(out_dir, exist_ok=True)
    payload = result.to_dict()
    payload["schema"] = EXPRESSION_DE_RESULT_SCHEMA
    payload["generated_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    problems = validate_expression_de_result(payload)
    if problems:
        raise ValueError(f"refusing to write an invalid expression DE result: {problems}")
    path = os.path.join(out_dir, "EXPRESSION_DE_RESULT.json")
    tmp = path + ".inprogress"
    with open(tmp, "w") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True, allow_nan=False)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)
    return result


def read_expression_de_result(path: str) -> dict:
    if not os.path.exists(path):
        raise SystemExit(f"expression DE result not found: {path}")
    with open(path) as fh:
        try:
            payload = json.load(fh, parse_constant=_reject_nonfinite_constant)
        except ValueError as e:
            raise SystemExit(f"expression DE result at {path} failed to parse: {e}")
    problems = validate_expression_de_result(payload)
    if problems:
        raise SystemExit(f"expression DE result at {path} failed validation: {problems}")
    return payload


def validate_model_masking_result(payload: dict) -> list[str]:
    """See `validate_expression_de_result`; the masking-stream analogue,
    checking `perturbable_genes`/`perturbable_genes_sha256` round-trip
    consistency instead of `input_mapping`.

    Per `review:20261014_100436` F3: `masking_score` is DECLARED
    CONFIGURATION, same status as `comparison`/`sample_structure` above --
    its nested field checks (every field a non-empty string, e.g. rejecting
    `score_definition=True`) apply when `status` is `'complete'` or
    `'not_run'`. A genuinely `'failed'` record (the REAL provider path for
    an invalid `MaskingScoreSpec`, e.g. `score_definition=True`) must still
    be writable/readable with the raw invalid value preserved as
    diagnostic -- it is the record OF that failure, not a record that
    claims the configuration was valid.

    Per `review:20261014_114357` R3: `perturbable_genes`' own
    per-entry VALUE check (non-empty string, no duplicates) used to run
    UNCONDITIONALLY, so a genuinely `'failed'` record from the REAL
    provider path for `perturbable_genes=(None,)` could not be written --
    the exact same asymmetry F3 already fixed for `masking_score`, missed
    here. Now gated by `require_valid_values`, same as `masking_score`;
    the hash recompute-and-compare stays unconditional, since it is value-
    agnostic (`_safe_json_identity` hashes a `None` entry without error).
    """
    problems = _validate_evidence_result_common(
        payload, MODEL_MASKING_RESULT_SCHEMA, "model_masking")
    if not isinstance(payload, dict):
        return problems
    status = payload.get("status")
    require_valid_values = status in ("complete", "not_run")
    genes = payload.get("perturbable_genes")
    if not isinstance(genes, list):
        problems.append("perturbable_genes must be a list")
    else:
        # Per `review:20261014_114357` R3: the per-entry VALUE
        # check (non-empty string, no duplicates) is gated the same way
        # `masking_score`'s fields are -- a genuinely `'failed'` record
        # (the REAL provider path for e.g. `perturbable_genes=(None,)`)
        # must preserve that raw invalid value as diagnostic, never be
        # rejected by the writer for carrying it. The hash recompute-and-
        # compare stays UNCONDITIONAL: it is a tamper/staleness check on
        # whatever `perturbable_genes` actually contains, valid or not,
        # and `_safe_json_identity` hashes `None` entries without error.
        if require_valid_values:
            if not all(isinstance(g, str) and g for g in genes):
                problems.append("perturbable_genes entries must all be non-empty strings")
            elif len(set(genes)) != len(genes):
                problems.append("perturbable_genes contains duplicate value(s)")
        expected_hash = _safe_json_identity(genes) or ""
        if payload.get("perturbable_genes_sha256") != expected_hash:
            problems.append(
                "perturbable_genes_sha256 does not match a recomputed hash of "
                "perturbable_genes -- tampered or stale")
    mip = payload.get("model_identity")
    if mip != "model_undecided" and not isinstance(mip, dict):
        problems.append(
            "model_identity must be the literal string 'model_undecided' or a dict")
    masking_score = payload.get("masking_score")
    if not isinstance(masking_score, dict) or not masking_score:
        problems.append("'masking_score' must be a non-empty dict")
    elif require_valid_values:
        problems.extend(_require_nonempty_str_fields(
            masking_score, ("score_definition", "masking_operation",
                           "reference_or_token_semantics",
                           "perturbable_gene_set_definition",
                           "donor_or_animal_aggregation_rule"), "masking_score"))
    return problems


def write_model_masking_result(out_dir: str, result: ModelMaskingResult) -> ModelMaskingResult:
    if os.path.exists(out_dir) and os.listdir(out_dir):
        raise SystemExit(f"{out_dir} must be new or empty")
    os.makedirs(out_dir, exist_ok=True)
    payload = result.to_dict()
    payload["schema"] = MODEL_MASKING_RESULT_SCHEMA
    payload["generated_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    problems = validate_model_masking_result(payload)
    if problems:
        raise ValueError(f"refusing to write an invalid model masking result: {problems}")
    path = os.path.join(out_dir, "MODEL_MASKING_RESULT.json")
    tmp = path + ".inprogress"
    with open(tmp, "w") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True, allow_nan=False)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)
    return result


def read_model_masking_result(path: str) -> dict:
    if not os.path.exists(path):
        raise SystemExit(f"model masking result not found: {path}")
    with open(path) as fh:
        try:
            payload = json.load(fh, parse_constant=_reject_nonfinite_constant)
        except ValueError as e:
            raise SystemExit(f"model masking result at {path} failed to parse: {e}")
    problems = validate_model_masking_result(payload)
    if problems:
        raise SystemExit(f"model masking result at {path} failed validation: {problems}")
    return payload
