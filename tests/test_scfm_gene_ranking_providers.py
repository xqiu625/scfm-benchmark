"""Tests for the two gene-ranking evidence-stream contracts
(`gene_ranking_providers.py`), per `review:20261013_135634` section
3: type/validation/serialization only this round -- no real statistical or
model backend. Every `run()` call here must return `status in ('failed',
'not_run')`, never a fabricated `'complete'`.
"""

import os
import sys

import pytest

CODE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "scripts",
    "34_scfm_discovery_framework",
    "code",
)
sys.path.insert(0, CODE)

import gene_ranking_providers as GR  # noqa: E402
from enrichment_providers import GeneRanking  # noqa: E402

# ===========================================================================
# ExpressionDERankingProvider
# ===========================================================================


def _comparison(**overrides):
    base = dict(
        group_column="condition",
        reference_level="healthy",
        comparison_level="disease",
        direction_definition="disease vs. healthy, positive = higher in disease",
    )
    base.update(overrides)
    return GR.ComparisonSpec(**base)


def _sample_structure(**overrides):
    base = dict(
        donor_or_animal_id_column="donor_id",
        sample_id_column="sample_id",
        repeated_sampling_or_pairing="none_single_sample_per_donor",
        group_column="condition",
        design_type="independent_donors_unpaired",
    )
    base.update(overrides)
    return GR.SampleStructure(**base)


def test_expression_de_run_never_fabricates_a_complete_result(tmp_path):
    """Per `review:20261013_210131` F1: three independent healthy
    donors vs. three independent disease donors, under an EXPLICIT unpaired
    design, must pass the applicable structural check and land on the
    'inputs are fine, but no real backend exists yet' not_run path -- a
    donor-blocked term must never be forced onto this design.
    """
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    group_values = ["healthy", "healthy", "healthy", "disease", "disease", "disease"]
    donor_values = ["D1", "D2", "D3", "D4", "D5", "D6"]  # each donor only ONE level
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=group_values,
        donor_values=donor_values,
    )
    assert result.status == "not_run"
    assert result.ranking is None
    assert "ExpressionDERankingProvider.run is a contract" in result.reason
    assert result.estimator == GR.ExpressionDERankingProvider.INTENDED_ESTIMATOR
    assert result.n_samples == 6


def test_expression_de_run_accepts_a_genuinely_fully_paired_design():
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="paired_or_repeated_within_donor")
    # Each donor contributes BOTH compared levels -- an ESTIMABLE paired design.
    group_values = ["healthy", "disease", "healthy", "disease"]
    donor_values = ["D1", "D1", "D2", "D2"]
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=group_values,
        donor_values=donor_values,
    )
    assert result.status == "not_run", result.reason


def test_expression_de_run_rejects_a_partially_paired_design():
    """Per F1: one cross-group donor out of several is NOT sufficient to
    call a declared donor-blocked design estimable -- the prior round's own
    gap (`donors_in_both` nonempty was treated as a full pass).
    """
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="paired_or_repeated_within_donor")
    group_values = ["healthy", "disease", "healthy", "disease", "disease"]
    donor_values = ["D1", "D1", "D2", "D3", "D3"]  # only D1 is actually paired
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=group_values,
        donor_values=donor_values,
    )
    assert result.status == "failed"
    assert "contribute to only ONE" in result.reason


def test_expression_de_run_rejects_an_unpaired_design_with_a_straddling_donor():
    """Per F1: a declared independent-donor design must be checked against
    the REAL input, not just trusted -- a donor that actually appears in
    both compared groups contradicts the declaration.
    """
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    group_values = ["healthy", "disease", "healthy", "disease"]
    donor_values = ["D1", "D1", "D2", "D2"]  # both donors actually straddle
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=group_values,
        donor_values=donor_values,
    )
    assert result.status == "failed"
    assert "contribute to BOTH" in result.reason


def test_expression_de_rejects_an_unknown_design_type():
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="some_other_design")
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "disease"],
        donor_values=["D1", "D2"],
    )
    assert result.status == "failed"
    assert "design_type" in result.reason


def test_expression_de_rejects_inconsistent_group_columns():
    comparison = _comparison(group_column="condition")
    sample_structure = _sample_structure(
        group_column="different_column", design_type="independent_donors_unpaired"
    )
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "disease"],
        donor_values=["D1", "D2"],
    )
    assert result.status == "failed"
    assert "must be the SAME column" in result.reason


def test_expression_de_run_flags_declared_covariates_as_unverified():
    comparison = _comparison()
    sample_structure = _sample_structure(
        design_type="independent_donors_unpaired", covariate_columns=("age",)
    )
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "disease"],
        donor_values=["D1", "D2"],
    )
    assert result.status == "not_run"
    assert "UNVERIFIED" in result.reason
    assert "age" in result.reason


def test_expression_de_rejects_an_empty_comparison_spec_field():
    comparison = _comparison(direction_definition="")
    sample_structure = _sample_structure()
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "disease"],
        donor_values=["D1", "D2"],
    )
    assert result.status == "failed"
    assert "direction_definition is empty" in result.reason


def test_expression_de_rejects_reference_equal_to_comparison_level():
    comparison = _comparison(comparison_level="healthy")  # same as reference_level
    sample_structure = _sample_structure()
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy"],
        donor_values=["D1"],
    )
    assert result.status == "failed"
    assert "must differ" in result.reason


def test_expression_de_rejects_a_blank_estimator():
    comparison = _comparison()
    sample_structure = _sample_structure()
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "disease"],
        donor_values=["D1", "D2"],
        estimator="",
    )
    assert result.status == "failed"
    assert "estimator must be stated explicitly" in result.reason


# ===========================================================================
# R1 (review:20261013_222450) -- the REAL run() entry point must
# require BOTH compared groups present, and reject malformed donor/group
# values BEFORE check_design_estimability's internal sorted() can raise an
# uncaught TypeError on a None/str mix.
# ===========================================================================


def test_expression_de_rejects_a_comparison_missing_one_entire_group():
    """Per R1: the prior round's check only asked 'does ANY sample carry a
    compared group value', which `healthy`-only input passed -- it never
    asked whether BOTH `reference_level` and `comparison_level` were
    present.
    """
    comparison = _comparison()
    sample_structure = _sample_structure()
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "healthy", "healthy"],
        donor_values=["D1", "D2", "D3"],
    )
    assert result.status == "failed"
    assert "missing comparison level" in result.reason


def test_expression_de_rejects_an_empty_string_donor_id():
    comparison = _comparison()
    sample_structure = _sample_structure()
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "disease"],
        donor_values=["", "D2"],
    )
    assert result.status == "failed"
    assert "donor_values has 1 non-string/empty value" in result.reason


def test_expression_de_rejects_a_none_donor_id_without_crashing():
    """Per R1: a `None` donor mixed with string donors used to reach
    `sorted()` inside `check_design_estimability` uncaught, raising a bare
    `TypeError` instead of a structured `failed` record.
    """
    comparison = _comparison()
    sample_structure = _sample_structure()
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "disease"],
        donor_values=[None, "D2"],
    )
    assert result.status == "failed"
    assert "donor_values has 1 non-string/empty value" in result.reason


def test_expression_de_rejects_a_non_string_group_value():
    comparison = _comparison()
    sample_structure = _sample_structure()
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", 7],
        donor_values=["D1", "D2"],
    )
    assert result.status == "failed"
    assert "group_values has 1 non-string/empty value" in result.reason


def test_design_estimability_rejects_a_comparison_missing_one_entire_group():
    comparison = _comparison()
    problems = GR.check_design_estimability(
        ["healthy", "healthy"], ["D1", "D2"], comparison, "independent_donors_unpaired"
    )
    assert problems
    assert "missing comparison level" in problems[0]


def test_mutation_kill_expression_de_group_donor_value_validation(monkeypatch):
    """Disables `_validate_group_donor_values` (as if R1's entry-point check
    did not exist) and confirms the missing-comparison-group scenario above
    would then be silently accepted as `not_run` -- reproducing the
    original R1 defect through the real `run()` entry point.
    """
    monkeypatch.setattr(GR, "_validate_group_donor_values", lambda g, d, c: [])
    comparison = _comparison()
    sample_structure = _sample_structure()
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "healthy", "healthy"],
        donor_values=["D1", "D2", "D3"],
    )
    assert result.status == "not_run"  # with the mutant: a one-sided comparison is no longer caught


# ===========================================================================
# check_design_estimability -- review:20261013_135634 +
# review:20261013_210131 F1: a donor-blocked design must not be
# mechanically applied to every comparison, and the declared design_type
# (not free text, not an estimator= string) decides which check applies.
# ===========================================================================


def test_design_estimability_accepts_a_fully_paired_design():
    comparison = _comparison()
    group_values = ["healthy", "disease", "healthy", "disease"]
    donor_values = ["D1", "D1", "D2", "D2"]  # each donor contributes BOTH levels
    problems = GR.check_design_estimability(
        group_values, donor_values, comparison, "paired_or_repeated_within_donor"
    )
    assert problems == []


def test_design_estimability_rejects_a_design_fully_confounded_with_donor_when_paired_declared():
    comparison = _comparison()
    group_values = ["healthy", "healthy", "disease", "disease"]
    donor_values = ["D1", "D1", "D2", "D2"]  # D1 only healthy, D2 only disease
    problems = GR.check_design_estimability(
        group_values, donor_values, comparison, "paired_or_repeated_within_donor"
    )
    assert problems
    assert "contribute to only ONE" in problems[0]


def test_design_estimability_accepts_the_same_data_when_independent_declared():
    """Per F1: the SAME donor-confounded-by-construction data is a NORMAL,
    valid unpaired design when declared as such -- the design_type, not the
    data shape alone, decides which check applies.
    """
    comparison = _comparison()
    group_values = ["healthy", "healthy", "disease", "disease"]
    donor_values = ["D1", "D1", "D2", "D2"]
    problems = GR.check_design_estimability(
        group_values, donor_values, comparison, "independent_donors_unpaired"
    )
    assert problems == []


def test_design_estimability_rejects_partial_pairing_even_with_one_crossing_donor():
    """Per F1: previously, `donors_in_both` being non-empty (even for just
    ONE donor out of several) was treated as a full pass. One crossing
    donor is no longer sufficient to call the whole design estimable.
    """
    comparison = _comparison()
    group_values = ["healthy", "disease", "healthy", "disease", "disease"]
    donor_values = ["D1", "D1", "D2", "D3", "D3"]  # only D1 is paired
    problems = GR.check_design_estimability(
        group_values, donor_values, comparison, "paired_or_repeated_within_donor"
    )
    assert problems
    assert "contribute to only ONE" in problems[0]


def test_design_estimability_rejects_an_unknown_design_type():
    comparison = _comparison()
    problems = GR.check_design_estimability(
        ["healthy", "disease"], ["D1", "D2"], comparison, "something_else"
    )
    assert problems
    assert "not in" in problems[0]


def test_expression_de_run_rejects_an_unestimable_paired_design_end_to_end():
    """Real entry point: ExpressionDERankingProvider.run itself refuses an
    unestimable donor-confounded PAIRED design, not just the standalone
    checker.
    """
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="paired_or_repeated_within_donor")
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "healthy", "disease", "disease"],
        donor_values=["D1", "D1", "D2", "D2"],
    )
    assert result.status == "failed"
    assert "contribute to only ONE" in result.reason


# ===========================================================================
# R2 (review:20261013_222450) -- a partially-paired design must be
# refused as UNSUPPORTED by this round's check, not described as
# mathematically inestimable (which is not generally true).
# ===========================================================================


def test_partial_pairing_rejection_message_states_unsupported_not_inestimable():
    comparison = _comparison()
    group_values = ["healthy", "disease", "healthy", "disease", "disease"]
    donor_values = ["D1", "D1", "D2", "D3", "D3"]  # only D1 is paired
    problems = GR.check_design_estimability(
        group_values, donor_values, comparison, "paired_or_repeated_within_donor"
    )
    assert problems
    assert "UNSUPPORTED" in problems[0]
    assert "NOT necessarily mathematically inestimable" in problems[0]
    assert "NOT silently dropped" in problems[0]


def test_partial_pairing_can_still_be_full_rank_demonstrating_R2_is_a_scope_limit():
    """Per R2: this is the review's own counterexample -- D1/D2 fully
    paired, D3 healthy-only, D4 disease-only. The `~ donor + condition`
    design matrix (intercept, D2/D3/D4 dummies, condition indicator) is 6
    rows x 5 columns with full column rank, i.e. genuinely estimable --
    even though this round's check still refuses it as UNSUPPORTED (a
    scope limit, not a mathematical necessity).
    """
    import numpy as np

    donors = ["D1", "D1", "D2", "D2", "D3", "D4"]
    conditions = ["healthy", "disease", "healthy", "disease", "healthy", "disease"]
    donor_levels = ["D1", "D2", "D3", "D4"]
    rows = []
    for d, c in zip(donors, conditions):
        intercept = [1.0]
        donor_dummies = [1.0 if d == lvl else 0.0 for lvl in donor_levels[1:]]
        condition_indicator = [1.0 if c == "disease" else 0.0]
        rows.append(intercept + donor_dummies + condition_indicator)
    design_matrix = np.array(rows)
    assert design_matrix.shape == (6, 5)
    assert np.linalg.matrix_rank(design_matrix) == 5  # full column rank -- genuinely estimable

    comparison = _comparison()
    problems = GR.check_design_estimability(
        conditions, donors, comparison, "paired_or_repeated_within_donor"
    )
    assert problems  # still refused: this round only supports FULL pairing
    assert "UNSUPPORTED" in problems[0]


def test_mutation_kill_expression_de_estimability_gate(monkeypatch):
    """Through the REAL `provider.run()` entry point: disables
    `check_design_estimability` (as if the estimability gate did not exist)
    and confirms the unestimable-design scenario above would then be
    silently accepted as `not_run` instead of refused -- proving the real
    check, not some other validation, is what the test above depends on.
    """
    monkeypatch.setattr(GR, "check_design_estimability", lambda *a, **k: [])
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="paired_or_repeated_within_donor")
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "healthy", "disease", "disease"],
        donor_values=["D1", "D1", "D2", "D2"],
    )
    assert result.status == "not_run"  # with the mutant: the unestimable design is no longer caught


# ===========================================================================
# ModelMaskingRankingProvider
# ===========================================================================


def _masking_score(**overrides):
    base = dict(
        score_definition="state-attribution direction score (pending Cindy confirmation)",
        masking_operation="set gene g's expression value to 0 before tokenization",
        reference_or_token_semantics="literal 0 in expression space",
        perturbable_gene_set_definition="all genes in the model's input vocabulary "
        "present in this dataset",
        donor_or_animal_aggregation_rule="mean over cells within donor, then mean " "over donors",
    )
    base.update(overrides)
    return GR.MaskingScoreSpec(**base)


def test_model_masking_run_never_fabricates_a_complete_result():
    masking_score = _masking_score()
    provider = GR.ModelMaskingRankingProvider()
    result = provider.run(
        masking_score=masking_score,
        model_identity_or_pending="model_undecided",
        perturbable_genes=("G1", "G2", "G3"),
    )
    assert result.status == "not_run"
    assert result.ranking is None
    assert "ModelMaskingRankingProvider.run is a contract" in result.reason


def test_model_masking_has_no_cell_selection_parameter():
    """The interface itself has no parameter that could mean 'mask a cell' --
    masking is a gene-input operation only, per the review's own correction
    of the prior handoff's cell-level-mask framing.
    """
    import inspect

    sig = inspect.signature(GR.ModelMaskingRankingProvider.run)
    param_names = set(sig.parameters) - {"self", "params"}
    assert not any("cell" in p.lower() for p in param_names)


def test_model_masking_rejects_an_empty_masking_score_field():
    masking_score = _masking_score(masking_operation="")
    provider = GR.ModelMaskingRankingProvider()
    result = provider.run(
        masking_score=masking_score,
        model_identity_or_pending="model_undecided",
        perturbable_genes=("G1",),
    )
    assert result.status == "failed"
    assert "masking_operation must be a non-empty string" in result.reason


def test_model_masking_rejects_empty_perturbable_genes():
    masking_score = _masking_score()
    provider = GR.ModelMaskingRankingProvider()
    result = provider.run(
        masking_score=masking_score,
        model_identity_or_pending="model_undecided",
        perturbable_genes=(),
    )
    assert result.status == "failed"
    assert "perturbable_genes is empty" in result.reason


def test_model_masking_rejects_duplicate_perturbable_genes():
    masking_score = _masking_score()
    provider = GR.ModelMaskingRankingProvider()
    result = provider.run(
        masking_score=masking_score,
        model_identity_or_pending="model_undecided",
        perturbable_genes=("G1", "G1"),
    )
    assert result.status == "failed"
    assert "duplicate" in result.reason


def test_model_masking_rejects_a_non_string_perturbable_gene():
    """Per R3: `_gene_set_identity` no longer uses `json.dumps(...,
    default=str)`, so a non-string gene ID must be refused structurally
    here instead of being silently stringified into the identity hash.
    """
    masking_score = _masking_score()
    provider = GR.ModelMaskingRankingProvider()
    result = provider.run(
        masking_score=masking_score,
        model_identity_or_pending="model_undecided",
        perturbable_genes=("G1", 7),
    )
    assert result.status == "failed"
    assert "perturbable_genes has 1 non-string/empty value" in result.reason


def test_model_masking_rejects_a_model_identity_that_is_not_the_sentinel_or_a_dict():
    masking_score = _masking_score()
    provider = GR.ModelMaskingRankingProvider()
    result = provider.run(
        masking_score=masking_score,
        model_identity_or_pending="scprint2_small_v2",
        perturbable_genes=("G1",),
    )
    assert result.status == "failed"
    assert "model_undecided" in result.reason


# ===========================================================================
# Shared GeneRanking binding -- both streams must produce the SAME
# GeneRanking type, each traceable to its own evidence_type.
# ===========================================================================


def test_make_evidence_ranking_binds_expression_de_ranking():
    ranking = GeneRanking(
        species="human",
        gene_id_namespace="gene_symbol",
        ranking_definition="expression_de: disease vs. healthy LFC",
        genes=("G1", "G2"),
        scores=(1.0, -1.0),
    )
    bound = GR.make_evidence_ranking("expression_de", ranking)
    assert bound is ranking


def test_make_evidence_ranking_rejects_a_ranking_not_naming_its_own_evidence_type():
    ranking = GeneRanking(
        species="human",
        gene_id_namespace="gene_symbol",
        ranking_definition="some unrelated score",
        genes=("G1", "G2"),
        scores=(1.0, -1.0),
    )
    with pytest.raises(ValueError, match="does not name"):
        GR.make_evidence_ranking("expression_de", ranking)


def test_make_evidence_ranking_rejects_an_invalid_ranking():
    ranking = GeneRanking(
        species="human",
        gene_id_namespace="gene_symbol",
        ranking_definition="model_masking: Delta(g) score",
        genes=("G1", "G1"),
        scores=(1.0, -1.0),
    )  # duplicate gene ID
    with pytest.raises(ValueError, match="cannot bind an invalid ranking"):
        GR.make_evidence_ranking("model_masking", ranking)


def test_make_evidence_ranking_rejects_an_unknown_evidence_type():
    ranking = GeneRanking(
        species="human",
        gene_id_namespace="gene_symbol",
        ranking_definition="expression_de: x",
        genes=("G1",),
        scores=(1.0,),
    )
    with pytest.raises(ValueError, match="not in"):
        GR.make_evidence_ranking("something_else", ranking)


def test_expression_de_result_to_dict_is_json_serializable():
    import dataclasses
    import json

    comparison = _comparison()
    sample_structure = _sample_structure()
    result = GR.ExpressionDEResult(
        status="complete",
        evidence_type="expression_de",
        ranking=GeneRanking(
            species="human",
            gene_id_namespace="gene_symbol",
            ranking_definition="expression_de: x",
            genes=("G1",),
            scores=(1.0,),
        ),
        comparison=dataclasses.asdict(comparison),
        sample_structure=dataclasses.asdict(sample_structure),
        estimator="placeholder",
        backend="placeholder_backend",
        backend_version="0.0",
    )
    json.dumps(result.to_dict())  # must not raise


def test_model_masking_result_to_dict_is_json_serializable():
    import dataclasses
    import json

    masking_score = _masking_score()
    result = GR.ModelMaskingResult(
        status="complete",
        evidence_type="model_masking",
        ranking=GeneRanking(
            species="human",
            gene_id_namespace="gene_symbol",
            ranking_definition="model_masking: x",
            genes=("G1",),
            scores=(1.0,),
        ),
        masking_score=dataclasses.asdict(masking_score),
        model_identity="model_undecided",
        backend="placeholder_backend",
        backend_version="0.0",
    )
    json.dumps(result.to_dict())  # must not raise


# ===========================================================================
# F4 (review:20261013_210131) -- the actual inputs this run was
# called with must be traceable in the serialized record, not just the
# comparison/masking-score CONFIGURATION.
# ===========================================================================


def test_model_masking_result_distinguishes_different_perturbable_gene_sets():
    """Per F4: two calls differing only in `perturbable_genes` must produce
    two DIFFERENT serialized records -- the prior round's own gap (neither
    the gene set nor a hash of it was ever stored).
    """
    masking_score = _masking_score()
    provider = GR.ModelMaskingRankingProvider()
    r1 = provider.run(
        masking_score=masking_score,
        model_identity_or_pending="model_undecided",
        perturbable_genes=("A",),
    )
    r2 = provider.run(
        masking_score=masking_score,
        model_identity_or_pending="model_undecided",
        perturbable_genes=("B",),
    )
    assert r1.perturbable_genes == ("A",)
    assert r2.perturbable_genes == ("B",)
    assert r1.perturbable_genes_sha256 != r2.perturbable_genes_sha256
    assert r1.to_dict() != r2.to_dict()


def test_mutation_kill_model_masking_perturbable_genes_identity(monkeypatch):
    """Disables `_gene_set_identity` (as if nothing about `perturbable_genes`
    ever reached the record) and confirms the two different-gene-set calls
    above would then serialize to the SAME record -- reproducing the
    original F4 defect exactly.
    """
    monkeypatch.setattr(GR, "_gene_set_identity", lambda genes: ((), "constant"))
    masking_score = _masking_score()
    provider = GR.ModelMaskingRankingProvider()
    r1 = provider.run(
        masking_score=masking_score,
        model_identity_or_pending="model_undecided",
        perturbable_genes=("A",),
    )
    r2 = provider.run(
        masking_score=masking_score,
        model_identity_or_pending="model_undecided",
        perturbable_genes=("B",),
    )
    assert r1.to_dict() == r2.to_dict()  # with the mutant: the original F4 bug reproduces


def test_expression_de_result_distinguishes_different_group_value_inputs():
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = GR.ExpressionDERankingProvider()
    r1 = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "disease"],
        donor_values=["D1", "D2"],
    )
    r2 = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["disease", "healthy"],
        donor_values=["D1", "D2"],
    )
    assert r1.group_values_sha256 != r2.group_values_sha256
    assert r1.to_dict() != r2.to_dict()


def test_mutation_kill_expression_de_sample_label_identity(monkeypatch):
    """Disables `_sample_label_identity` and confirms its OWN two fields
    (`group_values_sha256`/`donor_values_sha256`) would then be identical
    for two different `group_values` inputs -- reproducing the original F4
    defect for this specific summary-hash mechanism. Per
    `review:20261013_222450` R3, `input_mapping` (built
    independently by `_build_input_mapping`, see the mutation-kill test
    below) is now the AUTHORITATIVE identity carrier and still distinguishes
    the two records even with this one mechanism alone disabled -- so this
    test checks `_sample_label_identity`'s own two fields specifically,
    not the whole `to_dict()`.
    """
    monkeypatch.setattr(GR, "_sample_label_identity", lambda g, d: (0, "constant", "constant"))
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = GR.ExpressionDERankingProvider()
    r1 = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "disease"],
        donor_values=["D1", "D2"],
    )
    r2 = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["disease", "healthy"],
        donor_values=["D1", "D2"],
    )
    assert r1.group_values_sha256 == r2.group_values_sha256 == "constant"
    assert r1.donor_values_sha256 == r2.donor_values_sha256 == "constant"


def test_mutation_kill_expression_de_input_mapping_identity(monkeypatch):
    """Disables `_build_input_mapping` (as if R3's authoritative,
    round-trippable identity carrier did not exist) and confirms two
    different `group_values` inputs would then serialize to the SAME
    record overall -- reproducing the original F4 defect end to end,
    now that `input_mapping` is what actually carries it.
    """
    monkeypatch.setattr(GR, "_build_input_mapping", lambda sample_ids, g, d: ())
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = GR.ExpressionDERankingProvider()
    r1 = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "disease"],
        donor_values=["D1", "D2"],
    )
    r2 = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["disease", "healthy"],
        donor_values=["D1", "D2"],
    )
    # with the mutant: input_mapping is empty for BOTH records -- it no
    # longer carries any identity at all, reproducing R3's exact gap for
    # this specific field (even though other fields remain real here)
    assert r1.input_mapping == r2.input_mapping == ()


# ===========================================================================
# R3 (review:20261013_222450) -- a unified, strict write/read
# contract: the ACTUAL ordered input must round-trip through the persisted
# record, and a tampered mapping/hash must be rejected on read.
# ===========================================================================


def test_expression_de_result_records_an_ordered_input_mapping_with_sample_identity(tmp_path):
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "healthy", "disease", "disease"],
        donor_values=["D1", "D2", "D3", "D4"],
        sample_ids=["S1", "S2", "S3", "S4"],
    )
    assert result.sample_identity_available is True
    assert len(result.input_mapping) == 4
    assert result.input_mapping[0] == {
        "row_index": 0,
        "group_value": "healthy",
        "donor_value": "D1",
        "sample_id": "S1",
    }
    assert result.input_mapping[3] == {
        "row_index": 3,
        "group_value": "disease",
        "donor_value": "D4",
        "sample_id": "S4",
    }

    out_dir = tmp_path / "de_result"
    GR.write_expression_de_result(str(out_dir), result)
    path = out_dir / "EXPRESSION_DE_RESULT.json"
    assert path.exists()
    payload = GR.read_expression_de_result(str(path))
    assert payload["input_mapping"] == list(result.input_mapping)
    assert payload["input_mapping_sha256"] == result.input_mapping_sha256
    assert payload["status"] == result.status


def test_expression_de_result_without_sample_ids_never_fabricates_a_sample_id():
    """Per R3: `sample_ids=None` must never be papered over with a made-up
    sample ID -- each row carries only its `row_index`, and
    `sample_identity_available=False` records the gap explicitly.
    """
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "disease"],
        donor_values=["D1", "D2"],
    )
    assert result.sample_identity_available is False
    assert all("sample_id" not in row for row in result.input_mapping)
    assert [row["row_index"] for row in result.input_mapping] == [0, 1]


def test_expression_de_run_rejects_duplicate_sample_ids():
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "disease"],
        donor_values=["D1", "D2"],
        sample_ids=["S1", "S1"],
    )
    assert result.status == "failed"
    assert "sample_ids contains duplicate" in result.reason


def test_read_expression_de_result_rejects_a_tampered_input_mapping(tmp_path):
    """The real reader, not just the validator function in isolation: a
    hand-edited on-disk file with a changed `group_value` (and therefore a
    now-stale `input_mapping_sha256`) must be refused.
    """
    import json

    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "disease"],
        donor_values=["D1", "D2"],
        sample_ids=["S1", "S2"],
    )
    out_dir = tmp_path / "de_result"
    GR.write_expression_de_result(str(out_dir), result)
    path = out_dir / "EXPRESSION_DE_RESULT.json"

    payload = json.loads(path.read_text())
    payload["input_mapping"][0]["group_value"] = "disease"  # tamper, hash now stale
    path.write_text(json.dumps(payload))
    with pytest.raises(SystemExit, match="failed validation"):
        GR.read_expression_de_result(str(path))


def test_read_expression_de_result_rejects_a_tampered_hash_alone(tmp_path):
    import json

    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "disease"],
        donor_values=["D1", "D2"],
        sample_ids=["S1", "S2"],
    )
    out_dir = tmp_path / "de_result"
    GR.write_expression_de_result(str(out_dir), result)
    path = out_dir / "EXPRESSION_DE_RESULT.json"

    payload = json.loads(path.read_text())
    payload["input_mapping_sha256"] = "0" * 64  # tamper the hash alone, mapping untouched
    path.write_text(json.dumps(payload))
    with pytest.raises(SystemExit, match="failed validation"):
        GR.read_expression_de_result(str(path))


def test_mutation_kill_validate_expression_de_result_hash_recompute(monkeypatch, tmp_path):
    """Proves the hash-recompute check in `validate_expression_de_result` is
    what actually catches tampering -- first shows the real reader rejects
    the tampered file, then disables ONLY the recompute (as if it always
    agreed with whatever hash is on disk) and confirms the SAME tampered
    file is then silently accepted, reproducing R3's exact gap.
    """
    import json

    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "disease"],
        donor_values=["D1", "D2"],
        sample_ids=["S1", "S2"],
    )
    out_dir = tmp_path / "de_result"
    GR.write_expression_de_result(str(out_dir), result)
    path = out_dir / "EXPRESSION_DE_RESULT.json"

    payload = json.loads(path.read_text())
    payload["input_mapping"][0]["group_value"] = "disease"
    constant_hash = "0" * 64
    payload["input_mapping_sha256"] = constant_hash
    path.write_text(json.dumps(payload))

    with pytest.raises(SystemExit, match="failed validation"):
        GR.read_expression_de_result(str(path))  # real guard: rejects

    monkeypatch.setattr(GR, "_input_mapping_sha256", lambda mapping: constant_hash)
    reread = GR.read_expression_de_result(str(path))  # with the mutant: silently accepted
    assert reread["input_mapping"][0]["group_value"] == "disease"


def test_write_expression_de_result_refuses_a_nonfinite_param_value(tmp_path):
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "disease"],
        donor_values=["D1", "D2"],
        some_future_param=float("nan"),
    )
    out_dir = tmp_path / "de_result"
    with pytest.raises(ValueError, match="refusing to write an invalid"):
        GR.write_expression_de_result(str(out_dir), result)


def test_model_masking_result_write_read_round_trips(tmp_path):
    masking_score = _masking_score()
    provider = GR.ModelMaskingRankingProvider()
    result = provider.run(
        masking_score=masking_score,
        model_identity_or_pending="model_undecided",
        perturbable_genes=("A", "B"),
    )
    out_dir = tmp_path / "masking_result"
    GR.write_model_masking_result(str(out_dir), result)
    path = out_dir / "MODEL_MASKING_RESULT.json"
    assert path.exists()
    payload = GR.read_model_masking_result(str(path))
    assert payload["perturbable_genes"] == ["A", "B"]
    assert payload["perturbable_genes_sha256"] == result.perturbable_genes_sha256


def test_read_model_masking_result_rejects_a_tampered_gene_set(tmp_path):
    import json

    masking_score = _masking_score()
    provider = GR.ModelMaskingRankingProvider()
    result = provider.run(
        masking_score=masking_score,
        model_identity_or_pending="model_undecided",
        perturbable_genes=("A", "B"),
    )
    out_dir = tmp_path / "masking_result"
    GR.write_model_masking_result(str(out_dir), result)
    path = out_dir / "MODEL_MASKING_RESULT.json"

    payload = json.loads(path.read_text())
    payload["perturbable_genes"] = ["A", "C"]  # tamper, hash now stale
    path.write_text(json.dumps(payload))
    with pytest.raises(SystemExit, match="failed validation"):
        GR.read_model_masking_result(str(path))


# ===========================================================================
# F3 (review:20261014_100436) -- an illegitimate 'complete'/
# 'not_run' record must be rejected, while a GENUINELY 'failed' record
# (the real provider's own diagnostic of an invalid input) must still be
# writable/readable with the raw invalid value preserved.
# ===========================================================================


def test_read_expression_de_result_rejects_an_empty_ranking_dict_claiming_complete(tmp_path):
    """Per F3: the prior round's reader only checked `ranking is not None`
    -- an EMPTY dict passed that check. `ranking={}` plus the stub's own
    `backend='not_implemented'` placeholder must both be refused.
    """
    import json

    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "disease"],
        donor_values=["D1", "D2"],
        sample_ids=["S1", "S2"],
    )
    assert result.status == "not_run"
    out_dir = tmp_path / "de_result"
    GR.write_expression_de_result(str(out_dir), result)
    path = out_dir / "EXPRESSION_DE_RESULT.json"

    payload = json.loads(path.read_text())
    payload["status"] = "complete"
    payload["ranking"] = {}
    path.write_text(json.dumps(payload))
    with pytest.raises(SystemExit, match="failed validation"):
        GR.read_expression_de_result(str(path))


def test_read_expression_de_result_rejects_a_complete_record_with_placeholder_backend(tmp_path):
    import json

    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "disease"],
        donor_values=["D1", "D2"],
        sample_ids=["S1", "S2"],
    )
    out_dir = tmp_path / "de_result"
    GR.write_expression_de_result(str(out_dir), result)
    path = out_dir / "EXPRESSION_DE_RESULT.json"

    payload = json.loads(path.read_text())
    payload["status"] = "complete"
    payload["ranking"] = {
        "species": "human",
        "gene_id_namespace": "gene_symbol",
        "ranking_definition": "expression_de: x",
        "genes": ["G1"],
        "scores": [1.0],
    }
    # backend is still the 'not_implemented' placeholder -- must be rejected
    # even though the ranking dict itself is now well-formed.
    path.write_text(json.dumps(payload))
    with pytest.raises(SystemExit, match="failed validation"):
        GR.read_expression_de_result(str(path))


def test_read_expression_de_result_rejects_a_wrong_n_samples(tmp_path):
    import json

    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "disease"],
        donor_values=["D1", "D2"],
        sample_ids=["S1", "S2"],
    )
    out_dir = tmp_path / "de_result"
    GR.write_expression_de_result(str(out_dir), result)
    path = out_dir / "EXPRESSION_DE_RESULT.json"

    payload = json.loads(path.read_text())
    payload["n_samples"] = 999
    path.write_text(json.dumps(payload))
    with pytest.raises(SystemExit, match="failed validation"):
        GR.read_expression_de_result(str(path))


def test_read_expression_de_result_rejects_a_tampered_design_type_on_a_not_run_record(tmp_path):
    """Per F3: the baseline record here has `status='not_run'`, not
    'complete' -- the tampered `design_type` must still be caught, because
    `sample_structure`/`comparison` are DECLARED CONFIGURATION, checked
    regardless of status (unlike `input_mapping`'s raw per-row values,
    which are a diagnostic area only for 'failed').
    """
    import json

    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "disease"],
        donor_values=["D1", "D2"],
        sample_ids=["S1", "S2"],
    )
    assert result.status == "not_run"
    out_dir = tmp_path / "de_result"
    GR.write_expression_de_result(str(out_dir), result)
    path = out_dir / "EXPRESSION_DE_RESULT.json"

    payload = json.loads(path.read_text())
    payload["sample_structure"]["design_type"] = "not_a_design"
    path.write_text(json.dumps(payload))
    with pytest.raises(SystemExit, match="failed validation"):
        GR.read_expression_de_result(str(path))


def test_expression_de_failed_record_with_none_donor_value_still_writes_and_reads_back(tmp_path):
    """Per F3: the REAL provider path for an invalid `donor_value=None`
    input must produce a 'failed' record that is STILL writable/readable
    -- the prior round's gap applied the SAME strict per-row value checks
    to a 'failed' record's `input_mapping`, which is a diagnostic snapshot
    of the exact (invalid) input, not a record claiming to BE valid.
    """
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "disease"],
        donor_values=[None, "D2"],
        sample_ids=["S1", "S2"],
    )
    assert result.status == "failed"
    out_dir = tmp_path / "de_result"
    GR.write_expression_de_result(str(out_dir), result)  # must not raise
    path = out_dir / "EXPRESSION_DE_RESULT.json"
    payload = GR.read_expression_de_result(str(path))  # must not raise
    assert payload["status"] == "failed"
    assert payload["input_mapping"][0]["donor_value"] is None  # diagnostic value preserved


def test_masking_score_spec_rejects_a_non_string_score_definition_type():
    masking_score = _masking_score(score_definition=True)
    problems = masking_score.validate()
    assert any("score_definition" in p for p in problems)


def test_model_masking_run_with_invalid_masking_score_returns_failed_not_not_run():
    """Per F3: the prior round's `MaskingScoreSpec.validate()` used a bare
    truthiness check, so `score_definition=True` (truthy, not a string)
    passed validation and the provider returned 'not_run' -- as if the
    configuration had been checked and found fine.
    """
    masking_score = _masking_score(score_definition=True)
    provider = GR.ModelMaskingRankingProvider()
    result = provider.run(
        masking_score=masking_score,
        model_identity_or_pending="model_undecided",
        perturbable_genes=("G0",),
    )
    assert result.status == "failed"
    assert "score_definition" in result.reason


def test_model_masking_failed_record_with_invalid_score_definition_still_writes_and_reads_back(
    tmp_path,
):
    masking_score = _masking_score(score_definition=True)
    provider = GR.ModelMaskingRankingProvider()
    result = provider.run(
        masking_score=masking_score,
        model_identity_or_pending="model_undecided",
        perturbable_genes=("G0",),
    )
    assert result.status == "failed"
    out_dir = tmp_path / "masking_result"
    GR.write_model_masking_result(str(out_dir), result)  # must not raise
    path = out_dir / "MODEL_MASKING_RESULT.json"
    payload = GR.read_model_masking_result(str(path))  # must not raise
    assert payload["status"] == "failed"
    assert payload["masking_score"]["score_definition"] is True  # diagnostic value preserved


# ===========================================================================
# R2/R3 (review:20261014_114357) -- `n_samples`'s own TYPE must be
# checked unconditionally (not silently skipped for a non-int value), and a
# GENUINELY `'failed'` length-mismatch/invalid-gene record (produced by the
# REAL provider path) must still be writable/readable, mirroring F3's fix
# for `input_mapping`/`masking_score` values but for the length/presence
# relationship between `n_samples` and `input_mapping` instead.
# ===========================================================================


def test_read_expression_de_result_rejects_a_non_int_n_samples(tmp_path):
    import json

    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "disease"],
        donor_values=["D1", "D2"],
        sample_ids=["S1", "S2"],
    )
    out_dir = tmp_path / "de_result"
    GR.write_expression_de_result(str(out_dir), result)
    path = out_dir / "EXPRESSION_DE_RESULT.json"

    for bad_value in ("two", None):
        payload = json.loads(path.read_text())
        payload["n_samples"] = bad_value
        tampered = out_dir / f"tampered_{bad_value}.json"
        tampered.write_text(json.dumps(payload))
        with pytest.raises(SystemExit, match="n_samples must be an int"):
            GR.read_expression_de_result(str(tampered))


def test_expression_de_length_mismatch_failed_record_still_writes_and_reads_back(tmp_path):
    """Per R3: `group_values`/`donor_values` of DIFFERENT lengths is a
    genuine input error the real provider correctly flags as 'failed' --
    `_build_input_mapping` returns `()` (no row-aligned mapping can be
    built), but `n_samples=len(group_values)` is still 2, not 0. The prior
    round's writer applied an unconditional `n_samples == len(input_
    mapping)` check to EVERY status, which rejected this genuinely failed
    record for the same reason F3 found with `donor_value=None`.
    """
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "disease"],
        donor_values=["D0"],
    )
    assert result.status == "failed"
    assert "length mismatch" in result.reason
    assert result.input_mapping == ()
    assert result.input_mapping_buildable is False
    assert result.input_lengths == {"group_values": 2, "donor_values": 1, "sample_ids": None}
    assert result.n_samples == 2  # the RAW group_values length, never zeroed to "match"

    out_dir = tmp_path / "de_result"
    GR.write_expression_de_result(str(out_dir), result)  # must not raise
    path = out_dir / "EXPRESSION_DE_RESULT.json"
    payload = GR.read_expression_de_result(str(path))  # must not raise
    assert payload["status"] == "failed"
    assert payload["input_mapping"] == []
    assert payload["input_mapping_buildable"] is False
    assert payload["n_samples"] == 2


def test_model_masking_invalid_gene_failed_record_still_writes_and_reads_back(tmp_path):
    """Per R3: `perturbable_genes=(None,)` is a genuine, diagnosable input
    error from the REAL provider path -- its per-entry value check used to
    run unconditionally, rejecting this genuinely 'failed' record the same
    way `masking_score`'s fields were rejected before F3 gated them.
    """
    masking_score = _masking_score()
    provider = GR.ModelMaskingRankingProvider()
    result = provider.run(
        masking_score=masking_score,
        model_identity_or_pending="model_undecided",
        perturbable_genes=(None,),
    )
    assert result.status == "failed"
    assert "perturbable_genes" in result.reason

    out_dir = tmp_path / "masking_result"
    GR.write_model_masking_result(str(out_dir), result)  # must not raise
    path = out_dir / "MODEL_MASKING_RESULT.json"
    payload = GR.read_model_masking_result(str(path))  # must not raise
    assert payload["status"] == "failed"
    assert payload["perturbable_genes"] == [None]  # diagnostic value preserved


def test_read_expression_de_result_rejects_a_tampered_input_lengths(tmp_path):
    import json

    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "disease"],
        donor_values=["D1", "D2"],
        sample_ids=["S1", "S2"],
    )
    out_dir = tmp_path / "de_result"
    GR.write_expression_de_result(str(out_dir), result)
    path = out_dir / "EXPRESSION_DE_RESULT.json"

    payload = json.loads(path.read_text())
    payload["input_lengths"]["group_values"] = 999  # n_samples now disagrees
    path.write_text(json.dumps(payload))
    with pytest.raises(SystemExit, match="failed validation"):
        GR.read_expression_de_result(str(path))


def test_read_expression_de_result_rejects_a_tampered_input_mapping_buildable(tmp_path):
    """A length-mismatch failed record claiming `input_mapping_buildable=
    True` (contradicting its own `input_lengths`) must be rejected -- the
    flag is recomputed from `input_lengths`, never trusted on its own.
    """
    import json

    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "disease"],
        donor_values=["D0"],
    )
    assert result.status == "failed"
    out_dir = tmp_path / "de_result"
    GR.write_expression_de_result(str(out_dir), result)
    path = out_dir / "EXPRESSION_DE_RESULT.json"

    payload = json.loads(path.read_text())
    payload["input_mapping_buildable"] = True  # contradicts input_lengths (2 != 1)
    path.write_text(json.dumps(payload))
    with pytest.raises(SystemExit, match="failed validation"):
        GR.read_expression_de_result(str(path))


def test_validate_expression_de_result_rejects_a_tampered_buildable_flag():
    """`input_mapping_buildable` is recomputed-and-compared against
    `input_lengths`, never trusted as given -- a tampered `True` on a
    length-mismatched record (2 group_values, 1 donor_value: `False` is
    the only correct value) must be rejected with a message naming the
    recomputed value, mirroring every other recompute-and-compare field
    in this module (`input_mapping_sha256`, `identity_sha256`, etc).
    """
    import dataclasses

    comparison = dataclasses.asdict(_comparison())
    sample_structure = dataclasses.asdict(
        _sample_structure(design_type="independent_donors_unpaired")
    )
    payload = {
        "schema": GR.EXPRESSION_DE_RESULT_SCHEMA,
        "evidence_type": "expression_de",
        "status": "failed",
        "reason": "length mismatch",
        "ranking": None,
        "comparison": comparison,
        "sample_structure": sample_structure,
        "input_mapping": [],
        "input_mapping_sha256": GR._input_mapping_sha256(()),
        "n_samples": 2,
        "input_lengths": {"group_values": 2, "donor_values": 1, "sample_ids": None},
        "input_mapping_buildable": True,  # tampered: 2 != 1, this should be False
        "backend": "not_implemented",
        "backend_version": "n/a",
    }
    problems = GR.validate_expression_de_result(payload)
    assert any("input_mapping_buildable" in p and "does not match" in p for p in problems)


def test_mutation_kill_validate_expression_de_result_ranking_shape_check(monkeypatch):
    """Proves `_ranking_payload_problems` is what actually catches
    `ranking={}` claiming 'complete' -- disabling it (as if any non-None
    ranking were accepted) lets the SAME, otherwise-fully-valid payload
    pass with zero problems.
    """
    import dataclasses

    comparison = dataclasses.asdict(_comparison())
    sample_structure = dataclasses.asdict(
        _sample_structure(design_type="independent_donors_unpaired")
    )
    payload = {
        "schema": GR.EXPRESSION_DE_RESULT_SCHEMA,
        "evidence_type": "expression_de",
        "status": "complete",
        "ranking": {},
        "comparison": comparison,
        "sample_structure": sample_structure,
        "input_mapping": [],
        "input_mapping_sha256": GR._input_mapping_sha256(()),
        "n_samples": 0,
        "input_lengths": {"group_values": 0, "donor_values": 0, "sample_ids": None},
        "input_mapping_buildable": True,
        "backend": "real_backend",
        "backend_version": "1",
    }
    assert GR.validate_expression_de_result(payload)  # real guard: rejects

    monkeypatch.setattr(GR, "_ranking_payload_problems", lambda ranking: [])
    assert GR.validate_expression_de_result(payload) == []  # with the mutant: silently accepted


# ===========================================================================
# P2 (`review:20261014_151747`): an unbuildable length-mismatch
# record's `status` must be tied to `input_mapping_buildable`, not just
# cross-checked against `input_mapping`'s own shape.
# ===========================================================================


def test_read_expression_de_result_rejects_a_failed_record_promoted_to_not_run(tmp_path):
    """Mirrors the review's own `probe_status.py`: a genuinely `'failed'`
    length-mismatch record (2 `group_values`, 1 `donor_value`) is written
    and read back correctly as `'failed'`; rewriting ONLY its `status` to
    `'not_run'` on disk (every other field -- the stale `reason`, the empty
    `input_mapping`, `input_lengths`, `input_mapping_buildable=False` --
    left untouched) must now be rejected, not silently accepted.
    """
    import json

    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = GR.ExpressionDERankingProvider()
    result = provider.run(
        comparison=comparison,
        sample_structure=sample_structure,
        group_values=["healthy", "disease"],
        donor_values=["D0"],
    )
    assert result.status == "failed"
    out_dir = tmp_path / "de_result"
    GR.write_expression_de_result(str(out_dir), result)
    path = out_dir / "EXPRESSION_DE_RESULT.json"

    payload = json.loads(path.read_text())
    assert payload["input_mapping_buildable"] is False
    payload["status"] = "not_run"
    path.write_text(json.dumps(payload))
    with pytest.raises(SystemExit, match="failed validation"):
        GR.read_expression_de_result(str(path))


def test_validate_expression_de_result_rejects_an_unbuildable_record_claiming_complete():
    """Direct validator test with an otherwise-FULLY-VALID 'complete'
    payload (real backend/version, a well-formed ranking) -- the only
    defect is `input_lengths`/`input_mapping_buildable=False` disagreeing
    with `status='complete'`. Isolates that this specific new check fires
    independently of the pre-existing 'complete'-specific checks (real
    ranking/backend), which this payload otherwise satisfies.
    """
    import dataclasses

    comparison = dataclasses.asdict(_comparison())
    sample_structure = dataclasses.asdict(
        _sample_structure(design_type="independent_donors_unpaired")
    )
    payload = {
        "schema": GR.EXPRESSION_DE_RESULT_SCHEMA,
        "evidence_type": "expression_de",
        "status": "complete",
        "ranking": {
            "species": "human",
            "gene_id_namespace": "gene_symbol",
            "ranking_definition": "expression_de: x",
            "genes": ["G1"],
            "scores": [1.0],
        },
        "comparison": comparison,
        "sample_structure": sample_structure,
        "input_mapping": [],
        "input_mapping_sha256": GR._input_mapping_sha256(()),
        "n_samples": 2,
        "input_lengths": {"group_values": 2, "donor_values": 1, "sample_ids": None},
        "input_mapping_buildable": False,
        "backend": "real_backend",
        "backend_version": "1",
    }
    problems = GR.validate_expression_de_result(payload)
    assert any("status must be 'failed'" in p for p in problems)


def test_mutation_kill_validate_expression_de_result_unbuildable_status(monkeypatch):
    """Proves `_unbuildable_status_problems` is what actually ties
    `input_mapping_buildable=False` to `status` -- disabling it lets the
    SAME tampered-to-`'not_run'` payload from the round-trip test above
    (reconstructed by hand here) pass with zero problems.
    """
    import dataclasses

    comparison = dataclasses.asdict(_comparison())
    sample_structure = dataclasses.asdict(
        _sample_structure(design_type="independent_donors_unpaired")
    )
    payload = {
        "schema": GR.EXPRESSION_DE_RESULT_SCHEMA,
        "evidence_type": "expression_de",
        "status": "not_run",
        "reason": "group_values (2) and donor_values (1) length mismatch",
        "ranking": None,
        "comparison": comparison,
        "sample_structure": sample_structure,
        "input_mapping": [],
        "input_mapping_sha256": GR._input_mapping_sha256(()),
        "n_samples": 2,
        "input_lengths": {"group_values": 2, "donor_values": 1, "sample_ids": None},
        "input_mapping_buildable": False,
        "backend": "not_implemented",
        "backend_version": "n/a",
    }
    assert GR.validate_expression_de_result(payload)  # real guard: rejects

    monkeypatch.setattr(GR, "_unbuildable_status_problems", lambda *a, **k: [])
    assert GR.validate_expression_de_result(payload) == []  # with the mutant: silently accepted


# ===========================================================================
# run_pydeseq2_from_prepared_design: the real backend call path.
# Per `review:20261015_104226` item 4. Reuses the EXACT PreparedDEDesign
# fixture machinery `test_scfm_de_design_providers.py` already uses -- a tiny,
# synthetic, local h5ad with hand-computable sums, no network, no real data.
# pydeseq2 is NOT installed in this environment, so only the
# validation-short-circuit and ImportError->not_run paths are exercised for
# real here; the full numeric path is written but skipped
# (`pytest.importorskip`) until an isolated env with pydeseq2 exists.
# ===========================================================================

import dataclasses as _dc  # noqa: E402

import de_design_providers as DD  # noqa: E402
import numpy as np  # noqa: E402
import pseudobulk_providers as PB  # noqa: E402
from dataset_schema import NO_BATCH_SENTINEL  # noqa: E402


def _rb_base_config(matrix_path, **overrides):
    config = {
        "dataset_id": "test_real_backend",
        "species": "human",
        "tissue_or_target_cell_type": "monocyte",
        "subject_id_column": "donor_id",
        "sample_id_column": "sample_id",
        "biological_replicate_unit": "donor_id",
        "technical_batch_column": NO_BATCH_SENTINEL,
        "expression_matrix_path": matrix_path,
        "counts_semantics": "raw_counts",
        "counts_slot": "X",
        "gene_id_namespace": "gene_symbol",
        "cell_id_column": "cell_barcode",
        "qc_status": "provisional",
        "qc_version_ref": "none (provisional)",
    }
    config.update(overrides)
    return config


#: Per review:20261015_142903 R3: the panel's gene names, in the
#: SAME order as `_RB_X`'s columns -- fixed here (not re-derived from X's
#: shape) so a reader can see which column is which gene without counting.
_RB_GENE_NAMES = (
    "GENE_EFFECT_UP",
    "GENE_EFFECT_DOWN",
    "GENE_NULL",
    "GENE_BG1",
    "GENE_BG2",
    "GENE_BG3",
    "GENE_BG4",
)


def _rb_make_h5ad(
    tmp_path, *, X, cell_ids, sample_ids, donor_ids, conditions, name="real_backend.h5ad"
):
    import anndata as ad
    import pandas as pd

    obs = {
        "cell_barcode": cell_ids,
        "donor_id": donor_ids,
        "sample_id": sample_ids,
        "condition": conditions,
    }
    obs_df = pd.DataFrame(obs, index=[f"idx{i}" for i in range(len(cell_ids))])
    var = pd.DataFrame(index=list(_RB_GENE_NAMES[: X.shape[1]]))
    a = ad.AnnData(X=np.asarray(X, dtype="float32"), obs=obs_df, var=var)
    path = tmp_path / name
    a.write_h5ad(path)
    return str(path)


#: Per review:20261015_142903 R3, wording corrected per
#: review:20261015_150231 item 4: the PRIOR 3-gene fixture (one
#: gene fixed at the SAME value, 80, in every single sample) distorted
#: median-ratio size-factor estimation -- the review's own independent
#: diagnostic found a nominally-null gene normalizing to log2FC approx -2.68,
#: not approx 0. The root cause was NOT "a constant gene can never be a
#: normalization anchor" in general -- with only 3 genes total, 2 of which
#: carried a large injected shift, the single constant gene ended up as
#: essentially the ONLY stable reference the geometric-mean/median-of-ratios
#: calculation had to work with, which skewed the resulting per-sample size
#: factor. This fixture instead has FOUR stable background
#: genes (GENE_BG1-4) whose per-sample values vary (never a repeated
#: constant) but whose GROUP MEANS match within ~2% -- a real, if synthetic,
#: normalization anchor -- plus one up-regulated, one down-regulated, and
#: one genuinely null (varying, group-mean-matched) gene of interest.
#: Independently verified (median-of-ratios against each gene's own
#: geometric mean across all 6 samples, the same normalization DESeq2 itself
#: uses -- NOT a pydeseq2 run, a hand-checkable arithmetic diagnostic):
#: resulting size factors are 0.98-1.02 for all 6 samples (no sample
#: dominates); normalized disease-vs-healthy log2FC is approx +3.29 for
#: GENE_EFFECT_UP, approx -3.28 for GENE_EFFECT_DOWN, approx -0.02 for
#: GENE_NULL, and |log2FC| < 0.02 for all four background genes. 6
#: independent donors, 1 sample each (so pseudobulk counts == the cell's own
#: raw X): D0-D2 healthy, D3-D5 disease.
_RB_X = np.array(
    [
        #  EFFECT_UP  EFFECT_DOWN  NULL   BG1    BG2    BG3    BG4
        [50, 520, 98, 300, 290, 310, 295],  # D0 healthy
        [55, 500, 103, 310, 305, 295, 308],  # D1 healthy
        [52, 510, 101, 295, 310, 300, 302],  # D2 healthy
        [520, 50, 95, 305, 295, 308, 300],  # D3 disease
        [500, 55, 104, 290, 300, 312, 296],  # D4 disease
        [510, 52, 99, 315, 308, 296, 304],  # D5 disease
    ],
    dtype="float32",
)
_RB_CELL_IDS = [f"CELL{i}" for i in range(6)]
_RB_SAMPLE = [f"S{i}" for i in range(6)]
_RB_DONOR = [f"D{i}" for i in range(6)]
_RB_CONDITION = ["healthy", "healthy", "healthy", "disease", "disease", "disease"]


def _rb_prepared_fixture(tmp_path):
    """Builds a REAL `status='prepared'` PreparedDEDesign + its real,
    hash-bound counts, via the actual production `PseudobulkAggregationProvider`
    / `DEDesignPreparationProvider` -- never a hand-built dict standing in for
    either producer's output.
    """
    h5ad = _rb_make_h5ad(
        tmp_path,
        X=_RB_X,
        cell_ids=_RB_CELL_IDS,
        sample_ids=_RB_SAMPLE,
        donor_ids=_RB_DONOR,
        conditions=_RB_CONDITION,
    )
    config = _rb_base_config(h5ad)
    config_path = str(tmp_path / "dataset_config.json")
    pb_provider = PB.PseudobulkAggregationProvider()
    pb_out = str(tmp_path / "pb_real_backend")
    pb_result = pb_provider.aggregate(
        config,
        config_path,
        out_dir=pb_out,
        cell_ids_subset=tuple(_RB_CELL_IDS),
        cell_type_or_state_label="monocyte",
        condition_column="condition",
    )
    assert pb_result.status == "complete", pb_result.reason
    pb_path = os.path.join(pb_out, "PSEUDOBULK_RESULT.json")

    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    de_provider = DD.DEDesignPreparationProvider()
    prep_out = str(tmp_path / "prep_real_backend")
    prepared = de_provider.prepare(
        pseudobulk_result_path=pb_path,
        comparison=comparison,
        sample_structure=sample_structure,
        out_dir=prep_out,
    )
    assert prepared.status == "prepared", prepared.reason

    manifest = DD.read_prepared_de_design(os.path.join(prep_out, "PREPARED_DE_DESIGN.json"))
    counts = DD.load_prepared_de_counts(manifest)
    return manifest, counts


def test_run_pydeseq2_from_prepared_design_requires_status_prepared():
    provider = GR.ExpressionDERankingProvider()
    with pytest.raises(ValueError, match="status='prepared'"):
        provider.run_pydeseq2_from_prepared_design(
            prepared_design={"status": "failed"}, counts=np.zeros((1, 1))
        )


def test_run_pydeseq2_from_prepared_design_returns_not_run_without_pydeseq2(tmp_path):
    """REAL execution against the REAL PreparedDEDesign/counts: no mock, no
    stand-in. This environment (`302006_env_clean`) genuinely has no
    `pydeseq2` install, so this is the actual honest not_run path, not a
    simulated one.
    """
    manifest, counts = _rb_prepared_fixture(tmp_path)
    provider = GR.ExpressionDERankingProvider()
    result = provider.run_pydeseq2_from_prepared_design(prepared_design=manifest, counts=counts)
    try:
        import pydeseq2  # noqa: F401

        pytest.skip(
            "pydeseq2 IS installed in this environment; the ImportError "
            "path this test checks does not apply here"
        )
    except ImportError:
        pass
    assert result.status == "not_run", result.reason
    assert "pydeseq2" in result.reason and "not importable" in result.reason
    assert result.backend == "pydeseq2"
    assert result.backend_version == "unavailable"


def test_run_pydeseq2_from_prepared_design_short_circuits_on_a_failed_stub(tmp_path):
    """A structurally invalid input (here: a comparison whose group_column
    disagrees with sample_structure's) must be refused by the SAME `self.run`
    validation the contract-stub path already applies, with pydeseq2 never
    even attempted -- proven by using a comparison/sample_structure pairing
    that `self.run` itself rejects BEFORE any pydeseq2 import is reached.
    """
    manifest, counts = _rb_prepared_fixture(tmp_path)
    tampered = dict(manifest)
    tampered["comparison"] = dict(manifest["comparison"])
    tampered["comparison"]["group_column"] = "not_the_real_column"
    provider = GR.ExpressionDERankingProvider()
    result = provider.run_pydeseq2_from_prepared_design(prepared_design=tampered, counts=counts)
    assert result.status == "failed"
    assert "group_column" in result.reason
    assert "pydeseq2" not in result.reason.lower()


def test_run_pydeseq2_from_prepared_design_rejects_a_counts_shape_mismatch(tmp_path):
    manifest, counts = _rb_prepared_fixture(tmp_path)
    wrong_shape_counts = np.hstack([counts, np.ones((counts.shape[0], 1))])
    provider = GR.ExpressionDERankingProvider()
    result = provider.run_pydeseq2_from_prepared_design(
        prepared_design=manifest, counts=wrong_shape_counts
    )
    assert result.status == "failed"
    assert "shape" in result.reason


# Per review:20261015_142903 R1: a probe found a same-shape counts
# content change, a fractional count value, a flipped contrast sign, and an
# all-zero design matrix all silently reaching the backend call unrejected
# (the "failed" result the probe saw was a SPY forcibly raising inside the
# constructor, never production code refusing the bad input). These 4 tests
# exercise the REAL production code path -- no spy -- and prove each is now
# refused before pydeseq2 is even imported.


def test_run_pydeseq2_from_prepared_design_rejects_same_shape_counts_content_tamper(tmp_path):
    manifest, counts = _rb_prepared_fixture(tmp_path)
    tampered_counts = counts.copy()
    tampered_counts[0, 0] += 1  # same shape, one value off by one
    provider = GR.ExpressionDERankingProvider()
    result = provider.run_pydeseq2_from_prepared_design(
        prepared_design=manifest, counts=tampered_counts
    )
    assert result.status == "failed"
    assert "source-binding" in result.reason
    assert "pydeseq2" not in result.reason.lower()


def test_run_pydeseq2_from_prepared_design_rejects_fractional_counts(tmp_path):
    manifest, counts = _rb_prepared_fixture(tmp_path)
    tampered_counts = counts.copy()
    tampered_counts[0, 0] += 0.25  # same shape, non-integer-valued
    provider = GR.ExpressionDERankingProvider()
    result = provider.run_pydeseq2_from_prepared_design(
        prepared_design=manifest, counts=tampered_counts
    )
    assert result.status == "failed"
    assert "non-integer" in str(result.reason)
    assert "pydeseq2" not in result.reason.lower()


def test_run_pydeseq2_from_prepared_design_rejects_a_flipped_contrast_sign(tmp_path):
    manifest, counts = _rb_prepared_fixture(tmp_path)
    tampered = dict(manifest)
    tampered["contrast_vector"] = [-v for v in manifest["contrast_vector"]]
    provider = GR.ExpressionDERankingProvider()
    result = provider.run_pydeseq2_from_prepared_design(prepared_design=tampered, counts=counts)
    assert result.status == "failed"
    assert "re-validation" in result.reason
    assert "pydeseq2" not in result.reason.lower()


def test_run_pydeseq2_from_prepared_design_rejects_an_all_zero_design_matrix(tmp_path):
    """The design_matrix AND its own design_matrix_sha256 are tampered
    consistently with each other, so the first-layer "hash matches the
    stored matrix" check alone could not catch it -- only the semantic
    recompute from sample_ids/donor_ids/group_values/comparison (which
    rebuilds the TRUE, non-zero design matrix and compares) can. This is
    exactly the "VACUOUS estimability check" scenario de_design_providers.
    py's own module docstring describes.
    """
    manifest, counts = _rb_prepared_fixture(tmp_path)
    tampered = dict(manifest)
    zeroed = [[0.0 for _ in row] for row in manifest["design_matrix"]]
    tampered["design_matrix"] = zeroed
    tampered["design_matrix_sha256"] = DD.array_content_sha256(np.array(zeroed, dtype=float))
    provider = GR.ExpressionDERankingProvider()
    result = provider.run_pydeseq2_from_prepared_design(prepared_design=tampered, counts=counts)
    assert result.status == "failed"
    assert "re-validation" in result.reason
    assert "pydeseq2" not in result.reason.lower()


def test_run_pydeseq2_from_prepared_design_passes_an_ndarray_contrast_not_a_list(
    tmp_path, monkeypatch
):
    """Per R2: pydeseq2 v0.5.4's own DeseqStats constructor branches on the
    TYPE of `contrast` (ndarray -> numeric contrast; list -> 3-item
    categorical contrast spec, IndexError for a 2-column design). This
    injects a minimal fake `pydeseq2` package as a spy at the import
    boundary -- it only records the type of the `contrast=` argument and
    then raises, so this proves what TYPE the calling code constructs and
    passes, not that real pydeseq2 accepts it (that is independently
    verified directly against the real v0.5.4 source in
    DE_BACKEND_EXECUTION_PACKET_v1.md section 5 -- a double establishes
    control flow here, never third-party library behaviour).
    """
    import types

    manifest, counts = _rb_prepared_fixture(tmp_path)
    captured = {}

    fake_pydeseq2 = types.ModuleType("pydeseq2")
    fake_pydeseq2.__version__ = "0.0.0-fake"
    fake_dds_mod = types.ModuleType("pydeseq2.dds")
    fake_ds_mod = types.ModuleType("pydeseq2.ds")

    class _FakeDeseqDataSet:
        def __init__(self, *, counts, metadata, design, n_cpus=None):
            pass

        def deseq2(self):
            pass

    class _StopAfterCapture(Exception):
        pass

    class _FakeDeseqStats:
        def __init__(self, dds, contrast, n_cpus=None):
            captured["contrast_type"] = type(contrast)
            captured["contrast_value"] = contrast
            raise _StopAfterCapture("spy stops here -- never a real fit")

    fake_dds_mod.DeseqDataSet = _FakeDeseqDataSet
    fake_ds_mod.DeseqStats = _FakeDeseqStats
    monkeypatch.setitem(sys.modules, "pydeseq2", fake_pydeseq2)
    monkeypatch.setitem(sys.modules, "pydeseq2.dds", fake_dds_mod)
    monkeypatch.setitem(sys.modules, "pydeseq2.ds", fake_ds_mod)

    provider = GR.ExpressionDERankingProvider()
    result = provider.run_pydeseq2_from_prepared_design(prepared_design=manifest, counts=counts)

    assert "contrast_type" in captured, "the spy was never reached"
    assert captured["contrast_type"] is np.ndarray, (
        f"contrast was passed as {captured['contrast_type']!r}, not np.ndarray -- "
        "pydeseq2 v0.5.4 would treat a plain list as a 3-item categorical "
        "contrast and raise IndexError for this 2-column design"
    )
    assert result.status == "failed"
    assert "pydeseq2 raised during fit/stats" in result.reason


def test_short_circuit_reports_self_run_s_own_reason_before_the_slower_recheck(tmp_path):
    """`run_pydeseq2_from_prepared_design`'s short-circuit (`if stub.status
    != "not_run": return stub`) means a structurally-invalid input is
    refused with `self.run`'s OWN, specific reason (here: the
    `comparison.group_column`/`sample_structure.group_column` mutual
    mismatch self.run checks directly) -- cheaply, before the slower R1
    re-validation (`validate_prepared_de_design` + the live-source
    `_source_binding_problems` re-fetch) is even attempted. Confirms the
    reason is the SPECIFIC one, not a generic fallback.
    """
    manifest, counts = _rb_prepared_fixture(tmp_path)
    tampered = dict(manifest)
    tampered["sample_structure"] = dict(manifest["sample_structure"])
    tampered["sample_structure"]["group_column"] = "disagrees_with_comparisons_own_group_column"
    provider = GR.ExpressionDERankingProvider()
    result = provider.run_pydeseq2_from_prepared_design(prepared_design=tampered, counts=counts)
    assert result.status == "failed"
    assert "must be the SAME column" in result.reason


def test_bounded_n_cpus_reads_slurm_allocation_not_all_available_cores(monkeypatch):
    """Per R4: pydeseq2's own DefaultInference defaults n_cpus=None to
    os.cpu_count()-style "every core on the node" (verified directly
    against the real v0.5.4 source: `utils.get_num_processes`'s own
    docstring says "Returns the maximum number of available cpus by
    default"). `_bounded_n_cpus` must instead read the actual SLURM
    allocation, and fall back to the SAFE default (1), never silently to
    "every core," when no SLURM allocation is visible.
    """
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
    monkeypatch.delenv("SLURM_CPUS_ON_NODE", raising=False)
    assert GR._bounded_n_cpus() == 1

    monkeypatch.setenv("SLURM_CPUS_ON_NODE", "8")
    assert GR._bounded_n_cpus() == 8

    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "4")
    # SLURM_CPUS_PER_TASK takes priority when both are set.
    assert GR._bounded_n_cpus() == 4


def test_mutation_kill_run_pydeseq2_short_circuit(tmp_path, monkeypatch):
    """Proves the short-circuit in
    `run_pydeseq2_from_prepared_design` (`if stub.status != "not_run": return
    stub`) is load-bearing: removing it lets a structurally-failed input fall
    through past `self.run`'s own check entirely. Per
    `review:20261015_142903` R1, this project's input is now ALSO
    caught by the slower, independent `validate_prepared_de_design`/
    `_source_binding_problems` re-check (both the declared
    `comparison.group_column` and `sample_structure.group_column` are
    independently compared against the live upstream manifest's own
    `condition_column`) -- so the mutant no longer lets the bad input
    through completely unrejected; it demonstrates the short-circuit's
    DISTINCT role is reporting `self.run`'s own specific mismatch reason
    cheaply, not reaching a correct-but-different, more expensive reason via
    the re-validation path instead.
    """
    manifest, counts = _rb_prepared_fixture(tmp_path)
    tampered = dict(manifest)
    tampered["comparison"] = dict(manifest["comparison"])
    tampered["comparison"]["group_column"] = "not_the_real_column"
    provider = GR.ExpressionDERankingProvider()

    real_run = GR.ExpressionDERankingProvider.run

    def _mutant_always_not_run(self, *args, **kwargs):
        result = real_run(self, *args, **kwargs)
        return _dc.replace(result, status="not_run", reason="mutated: pretend valid")

    monkeypatch.setattr(GR.ExpressionDERankingProvider, "run", _mutant_always_not_run)
    # With the mutant, the short-circuit's own condition (`!= "not_run"`) can
    # never fire, so `self.run`'s own specific "...must be the SAME column"
    # mismatch message is never produced -- proving THAT message's gate is
    # the short-circuit, specifically (both messages otherwise share the
    # substring "comparison.group_column", so this test keys off the ONE
    # phrase unique to self.run's own check, not the shared field name).
    result = provider.run_pydeseq2_from_prepared_design(prepared_design=tampered, counts=counts)
    assert "must be the SAME column" not in result.reason
    # The input is still correctly refused overall (R1's independent
    # live-source re-check also catches this same tampered field) -- the
    # mutant changes WHICH check reports it, not whether it is caught.
    assert result.status == "failed"
    assert "condition_column" in result.reason


def test_run_pydeseq2_from_prepared_design_known_effect_and_null_ranks_correctly(tmp_path):
    """The real numeric path: a known large positive effect
    (GENE_EFFECT_UP), a known large negative effect (GENE_EFFECT_DOWN), a
    genuinely null gene (GENE_NULL, group-mean-matched but NOT a constant --
    see the fixture's own docstring for why that distinction matters), and
    four stable background genes that should calibrate near-null too, run
    through the REAL `PreparedDEDesign`/backend path end to end,
    independently re-read via `read_expression_de_result` (never trusting
    the in-memory object alone). SKIPPED in this environment -- `pydeseq2`
    is not installed in `302006_env_clean` and no isolated DE-backend
    environment exists yet. This test is written now so it runs for real
    the moment one does; it is NOT executed this round.
    """
    pytest.importorskip("pydeseq2", reason="no isolated DE-backend environment exists yet")
    manifest, counts = _rb_prepared_fixture(tmp_path)
    provider = GR.ExpressionDERankingProvider()
    result = provider.run_pydeseq2_from_prepared_design(prepared_design=manifest, counts=counts)
    assert result.status == "complete", result.reason

    out_dir = str(tmp_path / "expr_de_result")
    os.makedirs(out_dir)
    written = GR.write_expression_de_result(out_dir, result)
    reloaded = GR.read_expression_de_result(os.path.join(out_dir, "EXPRESSION_DE_RESULT.json"))
    assert reloaded["status"] == "complete"

    genes = list(reloaded["ranking"]["genes"])
    scores = dict(zip(genes, reloaded["ranking"]["scores"]))
    for g in _RB_GENE_NAMES:
        assert g in scores, f"{g} missing from the ranking -- dropped by the backend?"

    # KNOWN injected effects: disease ~10x healthy on GENE_EFFECT_UP (large
    # POSITIVE shrunken log2FoldChange; disease is the comparison level) and
    # the mirror-image ~10x DOWN on GENE_EFFECT_DOWN (large NEGATIVE).
    assert scores["GENE_EFFECT_UP"] > 1.0
    assert scores["GENE_EFFECT_DOWN"] < -1.0
    # The genuinely null/background genes must each sit well inside the two
    # known effects, in EITHER direction -- not asserting exactly 0 for any
    # single one (a single draw is not required to), but none of them may
    # approach the magnitude of a real 10x effect.
    effect_magnitude = min(abs(scores["GENE_EFFECT_UP"]), abs(scores["GENE_EFFECT_DOWN"]))
    for g in ("GENE_NULL", "GENE_BG1", "GENE_BG2", "GENE_BG3", "GENE_BG4"):
        assert abs(scores[g]) < effect_magnitude / 2, (
            f"{g} scored {scores[g]!r}, too close to the known +-10x effect "
            "magnitude for a gene with no intended group difference"
        )

    assert written.params["shrinkage_method"] == "pydeseq2.DeseqStats.lfc_shrink"
    assert "raw_wald_stat" in written.params and "raw_lfc" in written.params
    # Direction/coefficient-name consistency (R2's "also needs checking in a
    # real integration test"): the raw Wald stat's sign must agree with the
    # shrunken LFC's sign for both known-effect genes, and the contrast name
    # used for shrinkage must be the SAME one recorded for the raw stats.
    assert written.params["contrast_name"] == manifest["contrast_name"]
    for g in ("GENE_EFFECT_UP", "GENE_EFFECT_DOWN"):
        same_sign = (written.params["raw_wald_stat"][g] > 0) == (scores[g] > 0)
        assert same_sign, f"{g}: raw Wald stat sign disagrees with shrunken LFC sign"
