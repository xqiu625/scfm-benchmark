"""Dataset manifest/config schema for the scFM discovery framework.

Per `review:20261012_145004` item C.1 and `review:20261012_153316`
R1: species, tissue/target cell type, sample/subject/biological-replicate/
technical-batch identity, expression matrix location AND the specific slot counts
live in, counts semantics, gene ID type, cell ID column, and input QC/provenance.
Samples and subjects are never treated as equivalent -- a subject (patient/donor/
animal) may contribute more than one sample. Outcome/severity is accepted only as
optional analysis metadata; it is never a required field, and nothing here reads or
needs it to validate or plan.

No field is defaulted silently, and no required field is accepted as merely
"present and not None/''" -- every one of them is also type-checked (R1's reproduced
finding: `dataset_id=[]`, `subject_id_column=17`, `expression_matrix_path=False` all
previously passed because the empty-check alone does not catch wrong types).
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field

#: Known at this stage. An unlisted species is not an error by itself -- the
#: model registry decides compatibility (see model_registry.py) -- but an
#: unknown-to-us species in `species` must be spelled explicitly, not implied.
KNOWN_SPECIES = ("human", "mouse")

#: Explicit, closed vocabulary. "raw" and "normalized" are deliberately spelled
#: out rather than inferred from the file -- a config author must say which one
#: the referenced matrix actually holds.
KNOWN_COUNTS_SEMANTICS = ("raw_counts", "normalized_counts")

#: Explicit gene ID namespaces this framework can reason about. A config using
#: something else must still name it (not blank) so a downstream adapter can
#: decide it is unsupported, rather than silently mis-joining gene IDs.
KNOWN_GENE_ID_NAMESPACES = ("ensembl_gene_id", "gene_symbol", "entrez_id")

#: Explicit QC status vocabulary -- `analysis-ready` is never the default.
KNOWN_QC_STATUS = ("provisional", "analysis-ready")

#: Explicit sentinel for "this dataset genuinely has no technical-batch
#: variation to declare" -- distinct from an omitted/empty field. A config
#: author must say this on purpose; it is never inferred from a missing column.
NO_BATCH_SENTINEL = "none_single_batch"

#: Explicit literal for "this subject does not actually have repeat sampling" --
#: distinct from an omitted/empty field, mirroring NO_BATCH_SENTINEL.
NO_REPLICATE_SENTINEL = "none_single_sample_per_subject"

REQUIRED_FIELDS = (
    "dataset_id",
    "species",
    "tissue_or_target_cell_type",
    "subject_id_column",
    "sample_id_column",
    "biological_replicate_unit",
    "technical_batch_column",
    "expression_matrix_path",
    "counts_semantics",
    "counts_slot",
    "gene_id_namespace",
    "cell_id_column",
    "qc_status",
    "qc_version_ref",
)

#: Every required field must be a string, in addition to being non-empty. This
#: is enforced uniformly rather than per-field, since every field defined so
#: far is identifier/path/enum-like text -- a future non-string field must be
#: added to its own explicit check, not silently exempted here.
assert set(REQUIRED_FIELDS) == {
    "dataset_id", "species", "tissue_or_target_cell_type", "subject_id_column",
    "sample_id_column", "biological_replicate_unit", "technical_batch_column",
    "expression_matrix_path", "counts_semantics", "counts_slot",
    "gene_id_namespace", "cell_id_column", "qc_status", "qc_version_ref",
}

#: Accepted but never required. Present only as optional analysis metadata --
#: validate()/plan() must succeed with every one of these absent.
OPTIONAL_FIELDS = ("outcome_or_severity_column", "tissue_detail", "sex_column",
                   "collection_timepoint_column", "notes")


@dataclass(frozen=True)
class ValidationIssue:
    field: str
    reason: str

    def __str__(self) -> str:
        return f"{self.field}: {self.reason}"


@dataclass(frozen=True)
class ValidationResult:
    ok: bool
    issues: tuple[ValidationIssue, ...] = field(default_factory=tuple)

    def to_dict(self) -> dict:
        return {"ok": self.ok, "issues": [str(i) for i in self.issues]}


def is_valid_counts_slot(value: str) -> bool:
    """`"X"`, `"raw.X"`, or `"layers/<name>"` with a non-empty name -- states
    exactly where counts live, not merely that `counts_semantics` names what
    kind of value they are.
    """
    if value in ("X", "raw.X"):
        return True
    if value.startswith("layers/") and len(value) > len("layers/"):
        return True
    return False


def validate_dataset_config(config: dict) -> ValidationResult:
    """Structural validation only -- never opens the referenced expression
    matrix. Reading the matrix (even just its header/shape) is `cli.py`'s job,
    using the bounded h5py reads this project already established are
    necessary (`anndata(..., backed="r")` was found, in the batch-effect audit
    work, to still load every layer in full -- never used here or downstream).
    """
    issues: list[ValidationIssue] = []

    for f in REQUIRED_FIELDS:
        if f not in config:
            issues.append(ValidationIssue(f, "required field is missing"))
            continue
        value = config[f]
        if not isinstance(value, str):
            issues.append(ValidationIssue(
                f, f"required field must be a string, got {type(value).__name__} "
                   f"({value!r})"))
            continue
        if value == "":
            issues.append(ValidationIssue(f, "required field is present but empty"))

    def _is_present_string(key: str) -> bool:
        return key in config and isinstance(config[key], str) and config[key] != ""

    if _is_present_string("species"):
        if config["species"] not in KNOWN_SPECIES:
            issues.append(ValidationIssue(
                "species",
                f"{config['species']!r} is not in KNOWN_SPECIES {KNOWN_SPECIES} -- "
                "spelled explicitly is fine, but it will be reported as "
                "unsupported by every model registry entry until added there"))

    if _is_present_string("counts_semantics"):
        if config["counts_semantics"] not in KNOWN_COUNTS_SEMANTICS:
            issues.append(ValidationIssue(
                "counts_semantics",
                f"{config['counts_semantics']!r} not in {KNOWN_COUNTS_SEMANTICS} -- "
                "state exactly what the matrix holds, do not guess"))

    if _is_present_string("counts_slot"):
        if not is_valid_counts_slot(config["counts_slot"]):
            issues.append(ValidationIssue(
                "counts_slot",
                f"{config['counts_slot']!r} must be 'X', 'raw.X', or "
                "'layers/<name>' -- state exactly which slot counts live in"))

    if _is_present_string("gene_id_namespace"):
        if config["gene_id_namespace"] not in KNOWN_GENE_ID_NAMESPACES:
            issues.append(ValidationIssue(
                "gene_id_namespace",
                f"{config['gene_id_namespace']!r} not in {KNOWN_GENE_ID_NAMESPACES}"))

    if _is_present_string("qc_status"):
        if config["qc_status"] not in KNOWN_QC_STATUS:
            issues.append(ValidationIssue(
                "qc_status", f"{config['qc_status']!r} not in {KNOWN_QC_STATUS}"))
        elif config["qc_status"] == "analysis-ready" and not config.get("qc_version_ref"):
            issues.append(ValidationIssue(
                "qc_version_ref",
                "qc_status='analysis-ready' requires a non-empty qc_version_ref "
                "naming the specific QC pass/commit this status refers to"))

    if _is_present_string("subject_id_column") and _is_present_string("sample_id_column"):
        if config["subject_id_column"] == config["sample_id_column"]:
            issues.append(ValidationIssue(
                "subject_id_column",
                "subject_id_column must not equal sample_id_column -- a subject "
                "(patient/donor/animal) and a sample are not the same identity; "
                "one subject may contribute more than one sample"))

    unknown_extra = set(config) - set(REQUIRED_FIELDS) - set(OPTIONAL_FIELDS)
    for extra in sorted(unknown_extra):
        issues.append(ValidationIssue(
            extra, "not a recognized field (typo, or add it to REQUIRED_FIELDS/"
                   "OPTIONAL_FIELDS deliberately rather than letting it pass "
                   "silently)"))

    return ValidationResult(ok=not issues, issues=tuple(issues))


def load_dataset_config(path: str) -> dict:
    if not os.path.exists(path):
        raise SystemExit(f"dataset config not found: {path}")
    with open(path) as fh:
        return json.load(fh)
