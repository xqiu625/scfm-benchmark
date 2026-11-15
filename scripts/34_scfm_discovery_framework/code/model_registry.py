"""Model/checkpoint identity registry for the scFM discovery framework.

Per `review:20261012_145004` item B/C.1: model family/generation, official
repo + code revision, checkpoint file + revision/SHA-256, specific species-training
evidence, input gene-ID/preprocessing requirements, gene/cell output interfaces,
license, and resource evidence. Unknown fields are the literal string "unknown",
never left absent or guessed. "Pretraining includes species X" and "can technically
accept species X as input" are kept as two different fields -- never collapsed into
one boolean, per `review:20261012_133205` §4's own distinction.
"""
from __future__ import annotations

from dataclasses import dataclass, field

UNKNOWN = "unknown"

#: Confirmed (quoted from an official source) vs. not -- never merged. A
#: generation mismatch (quoting model N's doc while registering model N+1) is
#: exactly the mistake `review:20261012_145004` found and this schema
#: is built to make structurally harder to repeat silently.
SPECIES_EVIDENCE_STATES = (
    "checkpoint_metadata_confirmed",  # the SPECIFIC checkpoint file's own saved
                                       # hyperparameters/metadata name this species --
                                       # stronger than a README, since it is the
                                       # artifact itself, not prose about it
    "documented_confirmed",        # the model's OWN doc/paper names this species
    "documented_unconfirmed",      # secondary-sourced only; primary text not read
    "not_documented",              # checked, and the species is not named
    UNKNOWN,                       # not checked yet
)


@dataclass(frozen=True)
class ModelRegistryEntry:
    model_family: str
    model_generation: str                 # e.g. "scPRINT-2", not just "scPRINT"
    official_repo_url: str
    code_revision: str = UNKNOWN          # commit/tag; UNKNOWN is a valid, explicit value
    checkpoint_filename: str = UNKNOWN
    checkpoint_sha256: str = UNKNOWN
    checkpoint_sha256_source: str = UNKNOWN   # e.g. "verified: local file hashed and
                                               # matched against the official page" --
                                               # never silently assumed
    #: Per-species: "documented_confirmed" | "documented_unconfirmed" |
    #: "not_documented" | "unknown". A species absent from this dict is
    #: implicitly "unknown", not "not_documented" -- absence of a check is not
    #: a negative finding (DEVELOPMENT.md rule 3).
    species_pretraining_evidence: dict = field(default_factory=dict)
    #: Species the adapter can technically be POINTED at (e.g. a free-text
    #: --species flag accepts anything) -- kept separate from pretraining
    #: evidence above on purpose.
    species_input_acceptance: str = UNKNOWN
    input_gene_id_namespace_required: str = UNKNOWN
    input_preprocessing_notes: str = UNKNOWN
    #: e.g. ["cell_embedding"], ["cell_embedding", "gene_network"]. A causal
    #: claim is never implied by this list -- a gene-network/attention output
    #: is model-inferred structure, not evidence of real regulation.
    output_interfaces: tuple = ()
    license: str = UNKNOWN
    resource_notes: str = UNKNOWN
    #: Free-text pointer to where this entry's claims were checked (so a
    #: future edit can be traced back to its evidence, same convention as the
    #: Gate 2 manifests' `provenance` fields).
    evidence_source: str = UNKNOWN

    def to_dict(self) -> dict:
        d = dict(self.__dict__)
        d["output_interfaces"] = list(self.output_interfaces)
        return d


def validate_registry_entry(entry: ModelRegistryEntry) -> list[str]:
    """Returns a list of problems (empty = no structural problems found).

    This checks internal consistency, not biological truth: e.g. it will not
    accept a checkpoint_sha256_source that claims verification without an
    actual checkpoint_sha256 value, but it cannot tell you whether the model
    actually trained on the species it claims to.
    """
    problems = []
    if entry.checkpoint_sha256_source.startswith("verified") and entry.checkpoint_sha256 == UNKNOWN:
        problems.append(
            "checkpoint_sha256_source claims verification but checkpoint_sha256 is unknown")
    for species, state in entry.species_pretraining_evidence.items():
        if state not in SPECIES_EVIDENCE_STATES:
            problems.append(
                f"species_pretraining_evidence[{species!r}]={state!r} not in {SPECIES_EVIDENCE_STATES}")
    if not entry.output_interfaces:
        problems.append("output_interfaces is empty -- state at least one, even if "
                         "just 'cell_embedding', so downstream code cannot assume one")
    return problems


#: Two entries reflecting this round's independently-verified findings
#: (`scripts/34_scfm_discovery_framework/MODEL_SELECTION_v1.md`). Neither is
#: pinned as THE framework backbone -- that confirmation is the owner's,
#: still pending. These exist so `cli.py validate`/`plan` have real, not
#: placeholder, registry entries to exercise.
SCPRINT2_SMALL_V2 = ModelRegistryEntry(
    model_family="scPRINT",
    model_generation="scPRINT-2",
    official_repo_url="https://github.com/cantinilab/scPRINT-2",
    checkpoint_filename="small-v2.ckpt",
    checkpoint_sha256="2b586c144cf9a1f638b4b3e803554ebf6ce389f81c5f34179288a970addfd822",
    checkpoint_sha256_source=(
        "verified 2026-10-12: local cached file "
        "/bigdata/godziklab/shared/Xinru/.huggingface/scprint2_ckpts/small-v2.ckpt "
        "(889,706,878 bytes) hashed and matched exactly against the official "
        "huggingface.co/jkobject/scPRINT page's published SHA-256"),
    species_pretraining_evidence={
        # Upgraded 2026-10-12 per review:20261012_153316: the checkpoint
        # FILE's own pickled hyperparameters (not the README) contain an
        # "organisms" list. Independently re-verified this round by scanning the
        # checkpoint's archive/data.pkl bytes directly for NCBITaxon:* patterns
        # (no pickle deserialization, no REDUCE/GLOBAL/BUILD opcodes executed --
        # a pure byte/opcode scan): found exactly 16 unique NCBITaxon IDs,
        # including NCBITaxon:9606 (human, 1 occurrence) and NCBITaxon:10090
        # (mouse, 1 occurrence), matching the pickle disassembly's own
        # 16-element 'organisms' list exactly. This is checkpoint-metadata
        # evidence that these two species are IN the declared coverage list --
        # it is NOT a per-species training-cell count, train/val split, or any
        # downstream downstream efficacy claim, none of which are checked.
        "human": "checkpoint_metadata_confirmed",
        "mouse": "checkpoint_metadata_confirmed",
    },
    species_input_acceptance="ontology ID at inference, e.g. --species NCBITaxon:9606",
    input_gene_id_namespace_required=(
        "ensembl_gene_id OR gene_symbol (corrected 2026-10-12 per "
        "codex_output/20261012_153316.md -- the official README FAQ states "
        "var_names/var.index must be 'ENSEMBL_IDs or HUGO_SYMBOL', re-verified "
        "this round by re-fetching that FAQ entry directly; the first version "
        "of this registry wrongly required gene_symbol only)"),
    input_preprocessing_notes=(
        "raw counts accepted directly per README ('directly takes raw counts as "
        "input'); genes not recognized by the model's internal ontology are "
        "dropped at preprocessing -- observed in this project's own real mouse "
        "attempts (jobs 23758826/23760811/23761910) to drop tens of thousands of "
        "genes before failing or producing an unverified result"),
    output_interfaces=("cell_embedding", "gene_network"),
    license=UNKNOWN,
    resource_notes="ran on a single GPU node in this project's prior attempts; no "
                   "systematic resource ceiling established",
    evidence_source=(
        "scripts/34_scfm_discovery_framework/MODEL_SELECTION_v1.md (2026-10-12); "
        "species evidence upgraded and gene-ID requirement corrected per "
        "codex_output/20261012_153316.md, independently re-verified this round "
        "via a direct byte-level scan of the cached checkpoint's own pickle"),
)

UCE_4LAYER = ModelRegistryEntry(
    model_family="UCE",
    model_generation="UCE (4-layer checkpoint)",
    official_repo_url="https://github.com/snap-stanford/UCE",
    checkpoint_filename="4layer_model.torch",
    checkpoint_sha256="acb28f3f0a1d803e4a4ffe891b9bab38bf93c84762dc06b2452f0d515da91560",
    checkpoint_sha256_source=(
        "recorded, not verified: local file "
        "/bigdata/godziklab/shared/Xinru/302006/03_MODELS/uce/UCE/model_files/"
        "4layer_model.torch hashed this round, but the official UCE README gives "
        "no published checkpoint hash to compare against (it states files are "
        "'downloaded automatically', naming no SHA-256) -- unlike scPRINT-2's "
        "small-v2.ckpt above, this hash cannot be cross-checked against an "
        "official source yet"),
    species_pretraining_evidence={
        # review:20261012_153316 reports reading the UCE Nature paper
        # (s41586-026-10689-z) directly and finding its training-species
        # paragraph explicitly names human and mouse. My own two attempts to
        # fetch that URL this round were both blocked (an authentication
        # redirect to idp.nature.com) and a retry of the bioRxiv preprint was
        # rate-limited (HTTP 429) -- I could not independently reproduce this
        # specific read myself. Left at "documented_unconfirmed" from my own
        # verification standard rather than upgraded on attribution alone;
        # the secondary-sourced figure (human+mouse ~= 33.9M/36.2M training
        # cells) I DID independently obtain is consistent with the review's
        # claim, for what that is worth.
        "human": "documented_unconfirmed",
        "mouse": "documented_unconfirmed",
    },
    species_input_acceptance="free-text --species flag (no choices= restriction); "
                              "actual support depends on a per-species protein-"
                              "embedding asset file being present",
    input_gene_id_namespace_required="gene_symbol",
    input_preprocessing_notes=(
        "requires a per-species gene_symbol_to_embedding_ESM2.pt protein-embedding "
        "file; the mouse file (Mus_musculus.GRCm39...) is confirmed present in this "
        "project's local install, so a mouse run needs no new download"),
    output_interfaces=("cell_embedding",),
    license=UNKNOWN,
    resource_notes="200G RAM / 1 GPU / 10h ceiling used in this project's prior "
                   "mouse-full submit script (uce_mouse_full.slurm); not confirmed "
                   "sufficient, since every logged attempt failed before completion",
    evidence_source="scripts/34_scfm_discovery_framework/MODEL_SELECTION_v1.md (2026-10-12)",
)

REGISTRY = {
    "scprint2_small_v2": SCPRINT2_SMALL_V2,
    "uce_4layer": UCE_4LAYER,
}
