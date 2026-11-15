"""Model registry schema tests. Checks internal consistency only -- these tests
do not and cannot verify that a model actually trained on the species it claims;
they verify the registry cannot silently assert verification it never did.
"""

import os
import sys

CODE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "scripts",
    "34_scfm_discovery_framework",
    "code",
)
sys.path.insert(0, CODE)

import model_registry as R  # noqa: E402


def test_builtin_registry_entries_have_no_structural_problems():
    for key, entry in R.REGISTRY.items():
        problems = R.validate_registry_entry(entry)
        assert not problems, f"{key}: {problems}"


def test_scprint2_entry_species_evidence_is_checkpoint_metadata_confirmed():
    """Upgraded 2026-10-12 per review:20261012_153316: the checkpoint
    FILE's own pickled hyperparameters (re-verified this round via a direct
    byte scan, no deserialization) declare human+mouse among 16 organisms --
    stronger than a README/paper claim, since it is the specific checkpoint's
    own saved metadata, not prose about the model family.
    """
    entry = R.REGISTRY["scprint2_small_v2"]
    assert entry.checkpoint_sha256 != R.UNKNOWN
    assert "verified" in entry.checkpoint_sha256_source
    assert entry.species_pretraining_evidence["human"] == "checkpoint_metadata_confirmed"
    assert entry.species_pretraining_evidence["mouse"] == "checkpoint_metadata_confirmed"
    assert "checkpoint_metadata_confirmed" in R.SPECIES_EVIDENCE_STATES


def test_scprint2_entry_accepts_ensembl_or_hugo_gene_ids():
    """Corrected 2026-10-12 per review:20261012_153316: the official
    README FAQ allows Ensembl gene IDs OR HUGO symbols, not gene_symbol only.
    """
    entry = R.REGISTRY["scprint2_small_v2"]
    namespace_note = entry.input_gene_id_namespace_required.lower()
    assert "ensembl" in namespace_note and "hugo" in namespace_note


def test_uce_entry_hash_is_recorded_but_not_claimed_verified():
    entry = R.REGISTRY["uce_4layer"]
    assert entry.checkpoint_sha256 != R.UNKNOWN
    assert "not" in entry.checkpoint_sha256_source  # "recorded, not verified"
    assert "verified" not in entry.checkpoint_sha256_source.split("not")[0]


def test_claiming_verification_without_a_hash_is_rejected():
    bad = R.ModelRegistryEntry(
        model_family="X",
        model_generation="X-1",
        official_repo_url="https://example.invalid",
        checkpoint_sha256=R.UNKNOWN,
        checkpoint_sha256_source="verified somehow",
        output_interfaces=("cell_embedding",),
    )
    problems = R.validate_registry_entry(bad)
    assert any("claims verification" in p for p in problems)


def test_empty_output_interfaces_is_rejected():
    bad = R.ModelRegistryEntry(
        model_family="X",
        model_generation="X-1",
        official_repo_url="https://example.invalid",
        output_interfaces=(),
    )
    problems = R.validate_registry_entry(bad)
    assert any("output_interfaces is empty" in p for p in problems)


def test_invalid_species_evidence_state_is_rejected():
    bad = R.ModelRegistryEntry(
        model_family="X",
        model_generation="X-1",
        official_repo_url="https://example.invalid",
        species_pretraining_evidence={"human": "definitely yes"},
        output_interfaces=("cell_embedding",),
    )
    problems = R.validate_registry_entry(bad)
    assert any("not in" in p for p in problems)


def test_pretraining_evidence_and_input_acceptance_are_separate_fields():
    """review:20261012_133205 §4's own distinction: 'pretraining
    includes species X' and 'can technically accept species X as input' must
    never be collapsed into one field.
    """
    entry = R.REGISTRY["uce_4layer"]
    assert hasattr(entry, "species_pretraining_evidence")
    assert hasattr(entry, "species_input_acceptance")
    assert entry.species_pretraining_evidence is not entry.species_input_acceptance
