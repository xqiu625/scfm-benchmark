"""Gate 2 scoring-only runtime config DRAFT. IT DOES NOT PRODUCE A RUNNABLE CONFIG.

Implements the written ruling in
`reports/manuscript/revision/methodology_review_2026-09-24/REVIEW.md`, and stops
exactly where that ruling stops.

The subsequent GATE2_PARAMETER_RULING.md approves all eleven methodological
choices. The production slot remains empty because capture eligibility, config
audit and a scoped pilot authorization are separate gates. Every unit here still
carries `runtime_eligible = false`; this module cannot produce a runnable config.

Filling a default to make the draft "work" is the one thing this module exists to
prevent. There is no code path that guesses a value.

Nothing here runs Scrublet, submits a job, selects a rate, or deletes a cell.
"""
import json
import types

from gate2_stress_design import (COHORT_USE_CLASS, DECISION_REFERENCE, DESIGN_VERSION,
                                 NEVER_POOLED, RESIDUAL_RATE_GRID, USE_NOT_RUN,
                                 USE_RESIDUAL_DIAGNOSTIC, USE_RESIDUAL_SCORING_ONLY,
                                 USE_STANDARD_SCORING, Gate2DesignError,
                                 assert_grid_is_frozen)
from gate2_stress_runner import (PRODUCTION_PARAMETER_RULING, REQUIRED_PARAMETER_NAMES,
                                 Gate2RunnerError)
from gate2_parameter_ruling import (METHODOLOGY_PARAMETER_RULING,
                                   PARAMETER_RULING_REFERENCE,
                                   PARAMETER_RULING_VERSION)

CONFIG_DRAFT_VERSION = 'gate2-runtime-config-draft/2'
RULING_REFERENCE = ('reports/manuscript/revision/methodology_review_2026-09-24/'
                    'REVIEW.md :: A. Gate 2 Q1-Q5 and the per-cohort table')

#: The ONE label the grid values may carry. The ruling forbids calling them a
#: prevalence, a bound, a vendor estimate or a deletion quota, so the vocabulary
#: is frozen to a single string and checked.
STRESS_SETTING_LABEL = 'assumed_algorithmic_stress_setting'
FORBIDDEN_LABELS = ('prevalence', 'bound', 'vendor_estimate', 'vendor_rate',
                    'deletion_quota', 'expected_prevalence', 'residual_rate_estimate')

#: The ruling's common conditions. A unit is eligible only if EVERY one is
#: satisfied, and none is satisfied by default.
COMMON_CONDITIONS = types.MappingProxyType({
    'C1_input': 'un-normalised RNA counts with source evidence; integrality is necessary '
                'but not sufficient; ADT, integrated/scaled, log and DecontX-corrected '
                'matrices are excluded as additive simulation input',
    'C2_unit': 'a fixed physical capture/array and raw-barcode namespace; sequencing lane, '
               'GSM, Sample_ID and donor become analysis keys only on mapping evidence; '
               'captures are never merged to raise the cell count',
    'C3_support': 'the complete basic cell-called/QC population of that stage, not '
                  'pre-selected by target lineage and not balanced-subsampled; parent-state '
                  'support must be reviewed and residual-simulation limits stated; a '
                  'limitation label alone does not establish eligibility',
    'C4_fixed_implementation': 'parameters, software version, randomness, threshold rule '
                               'and effective-dimension condition fixed in advance; illegal '
                               'input, infeasible dimension, threshold failure and identity '
                               'mismatch are unknown/NA, never False, 0 or singlet; no '
                               'per-sample PCA reduction to rescue a run',
    'C5_interpretation': 'report coverage and results by cohort, capture and input stage; '
                         'no top-r quota, no outcome-driven tuning, no best-grid-point '
                         'selection, no stability-implies-correctness inference',
})

#: Per cohort: the ruling's decision, evidence labels and stop conditions,
#: transcribed rather than re-derived. `can_design_scoring_experiment` is the
#: ruling's own word -- "conditional" means the DESIGN is permitted, not that any
#: input has been shown to qualify.
_COHORT_RULING = {
    'qiu': dict(
        can_design_scoring_experiment='conditional', use=USE_STANDARD_SCORING,
        purpose='standard scoring sensitivity on the complete capture',
        evidence_labels=('published_protocol_loading', 'chemistry_unresolved', STRESS_SETTING_LABEL),
        stop_conditions=(
            'GSM -> physical capture unproven',
            'raw-count source unproven',
            'any implementation that depends on chemistry/vendor rate while the chemistry conflict is unresolved',
        ),
        preconditions_outstanding=(
            ('C1', 'raw-count source for the GSM-level input not established'),
            ('C2', 'GSM -> physical capture not established'),
            ('C3', 'parent-state support at the deposited stage not established'),
        ),
    ),
    'chongqing': dict(
        can_design_scoring_experiment='conditional', use=USE_STANDARD_SCORING,
        purpose='standard scoring sensitivity; missing loading is not manufactured',
        evidence_labels=('loading_unknown', 'upstream_handling_unknown', STRESS_SETTING_LABEL),
        stop_conditions=(
            '"none reported" treated as confirmed absence of upstream removal',
            'capture or input stage unclear, so the result cannot be bounded',
        ),
        preconditions_outstanding=(
            ('C1', 'deposited filtered matrix; un-normalised RNA provenance not established'),
            ('C2', 'GSM -> physical capture not established'),
            ('C3', 'input stage and upstream handling unknown'),
        ),
    ),
    'cillo': dict(
        can_design_scoring_experiment='conditional', use=USE_RESIDUAL_DIAGNOSTIC,
        purpose='residual_input_diagnostic on the post-HTO input',
        evidence_labels=('upstream_filter_published', 'composition_assumed', 'residual_null_unvalidated', STRESS_SETTING_LABEL),
        stop_conditions=(
            'HTO sample conflated with donor',
            'unrestricted pairing read as a calibrated within-sample residual',
            'only a subset too small to represent parent states remains',
        ),
        preconditions_outstanding=(
            ('C1', 'post-HTO input; raw-count provenance not established'),
            ('C2', 'GSM -> capture rests on an empirically mapped h5ad suffix; 8 libraries unmapped'),
            ('C3', 'post-filter subset may not represent parent states'),
        ),
    ),
    'liu': dict(
        can_design_scoring_experiment='conditional', use=USE_RESIDUAL_DIAGNOSTIC,
        purpose='residual_input_diagnostic on the post-demultiplex input',
        evidence_labels=('sng_composition_measured_at_post_demux', 'loading_input_proxy', 'residual_null_unvalidated', STRESS_SETTING_LABEL),
        stop_conditions=(
            'the 55 lane file keys unproven as physical captures',
            'RNA and ADT not separated',
            'barcode selection stage unclear',
            'any statement of a same-donor rate estimate',
        ),
        preconditions_outstanding=(
            ('C1', 'RNA/ADT separation at the scored stage not established'),
            ('C2', 'the 55 lane keys not established as physical captures'),
            ('C3', 'post-demultiplex input; parent-state support not established'),
        ),
    ),
    'combat': dict(
        can_design_scoring_experiment='conditional', use=USE_RESIDUAL_DIAGNOSTIC,
        purpose='residual-input diagnostic after multi-stage selection',
        # `residual_null_unvalidated` is NOT in the ruling's COMBAT table row; it
        # comes from the paragraph setting the shared boundary for the three
        # upstream-filtered cohorts, which applies the unvalidated-residual-null
        # limit to Cillo, Liu and COMBAT together. Taking only the table row would
        # have let COMBAT carry a weaker label than the two it is grouped with --
        # a test caught exactly that.
        evidence_labels=('published_nominal_loading_50000', 'composition_assumed', 'multiple_selection_stages', 'residual_null_unvalidated', STRESS_SETTING_LABEL),
        stop_conditions=(
            'a prior assigned by vendor extrapolation or truncation',
            'Channel_ID -> capture unclear',
            'a score read as validated residual prevalence',
        ),
        preconditions_outstanding=(
            ('C1', 'raw-count provenance for the deposited container not established'),
            ('C2', 'Channel_ID -> physical capture not established'),
            ('C3', 'multi-stage selection; parent-state support not established'),
        ),
    ),
    'ren': dict(
        can_design_scoring_experiment='conditional', use=USE_RESIDUAL_SCORING_ONLY,
        purpose='residual_scoring_only, reported separately, never pooled',
        evidence_labels=('author_scrublet_then_manual_filter', 'counts_level_agreement_only', STRESS_SETTING_LABEL),
        stop_conditions=(
            'a second Scrublet pass presented as independent validation',
            'counts or analysis unit unconfirmed, leaving only the historical count description',
            'any barcode-level transfer without separate identity evidence',
        ),
        preconditions_outstanding=(
            ('C1', 'counts-level agreement only; barcode-level provenance absent'),
            ('C2', 'sampleID -> capture asserted by inventory, not verified here'),
            ('C3', "input is downstream of the authors' own Scrublet plus manual removal"),
        ),
    ),
    'stephenson': dict(
        can_design_scoring_experiment='no', use=USE_NOT_RUN,
        purpose='no approved parameters or run design',
        evidence_labels=('unmapped_to_gem_well', 'not_eligible_without_additional_evidence'),
        stop_conditions=(
            'the 143 records do not decompose into legal captures',
            'the grid must not be used to paper over a missing unit',
        ),
        preconditions_outstanding=(
            ('C2', 'no legal capture unit: the 143 records do not decompose'),
        ),
    ),
    'wilk': dict(
        can_design_scoring_experiment='no', use=USE_NOT_RUN,
        purpose='no approved parameters or run design',
        evidence_labels=('seq_well', 'deposit_filter_stage_unknown', 'not_eligible_without_additional_evidence'),
        stop_conditions=(
            'which side of the heuristic filter the deposit sits on is unknown',
            'per-array identity unknown',
            'a 10x rate model applied to Seq-Well',
        ),
        preconditions_outstanding=(
            ('C2', 'per-array identity unknown'),
            ('C3', 'which side of the heuristic filter the deposit sits on is unknown'),
        ),
    ),
}
COHORT_RULING = types.MappingProxyType(_COHORT_RULING)

#: The ruling keeps `r`, `s` and `e` in the interpretation layer. Naming them here
#: is what makes their absence from the run axes checkable.
INTERPRETATION_ONLY_AXES = ('r_upstream_doublet_fraction', 's_same_source_probability',
                            'e_removal_efficiency')

#: Rank and AUC invariance across the grid follow from the score identity when the
#: neighbour graph is fixed. They are an implementation check, never evidence of
#: biological robustness.
INVARIANCE_CHECKS = types.MappingProxyType({
    'rank_invariance_across_grid': 'implementation_check_only',
    'observed_vs_simulated_auc_invariance': 'implementation_check_only',
    'auc_interpretation': 'distribution separability, NOT diagnostic accuracy against a '
                          'truth label: the observed set contains unknown doublets',
})

#: Fixed once per (capture, seed) and reused across all eight grid points.
FIXED_PER_CAPTURE_SEED = ('input_cell_order', 'hvg_selection', 'simulated_parent_pairs',
                          'pca_representation', 'neighbour_graph')


def unresolved_parameters(ruling=None):
    """The parameters with no methodological ruling. Sorted, so a checklist is stable."""
    ruled = METHODOLOGY_PARAMETER_RULING if ruling is None else ruling
    return sorted(REQUIRED_PARAMETER_NAMES - set(ruled))


def assert_label_permitted(label):
    if label != STRESS_SETTING_LABEL:
        raise Gate2DesignError(
            f'grid values may only be labelled {STRESS_SETTING_LABEL!r}; {label!r} is '
            'refused. The ruling forbids presenting them as a prevalence, bound, vendor '
            'estimate or deletion quota.')
    return label


def assert_no_runtime_axis_beyond_the_prior(axes):
    """`e` (and `r`, `s`) are interpretation-layer only."""
    bad = [a for a in axes if a in INTERPRETATION_ONLY_AXES]
    if bad:
        raise Gate2DesignError(
            f'{bad} are interpretation-layer parameters; the ruling permits ONE runtime '
            'axis, the prior grid. Adding them creates unidentifiable combinations that '
            'produce the same residual.')
    return tuple(axes)


def emit_runnable_config(ruling=None):
    """Refuses while any required parameter is unruled. There is no override."""
    missing = unresolved_parameters(ruling)
    if missing:
        raise Gate2RunnerError(
            f'cannot emit a runnable config: {len(missing)} parameter(s) have no '
            f'methodological ruling: {missing}. Missing values cannot be filled '
            'from library defaults.')
    raise Gate2RunnerError(
        'a complete parameter ruling is necessary but not sufficient: no reviewed '
        'per-cohort eligibility artifact exists, so the C1-C3 input preconditions are '
        'unestablished for every cohort and no unit is runtime_eligible.')


def config_draft():
    """The draft. Every unit is runtime_eligible=false, and the reasons are listed."""
    assert_grid_is_frozen(RESIDUAL_RATE_GRID)
    missing = unresolved_parameters()
    cohorts = {}
    for name, r in COHORT_RULING.items():
        eligible_in_principle = r['can_design_scoring_experiment'] == 'conditional'
        blockers = []
        if not eligible_in_principle:
            blockers.append('cohort is not_run under the ruling')
        else:
            # Per-cohort only. An earlier draft cited the Stage 31 preflight's
            # 302/302 unknown captures as a blocker on all eight cohorts. That is
            # a different work line: Stage 31's 302 Sample_IDs are a donor-level
            # platelet analysis, not the input evidence for these eight cohorts.
            # Coupling them would let one cohort's input audit silently gate
            # seven others.
            if missing:
                blockers.append('global: %d runtime parameter(s) unruled' % len(missing))
            blockers.append('global: config audit and scoped pilot authorization not recorded')
            blockers.extend('%s: %s' % (c, w)
                            for c, w in r['preconditions_outstanding'])
        cohorts[name] = {
            'use_class': COHORT_USE_CLASS[name],
            'can_design_scoring_experiment': r['can_design_scoring_experiment'],
            'purpose': r['purpose'],
            'evidence_labels': list(r['evidence_labels']),
            'stop_conditions': list(r['stop_conditions']),
            'never_pooled': name in NEVER_POOLED,
            'grid': list(RESIDUAL_RATE_GRID) if eligible_in_principle else [],
            'grid_label': STRESS_SETTING_LABEL if eligible_in_principle else None,
            'runtime_eligible': False,
            'runtime_eligible_blockers': blockers,
        }
    return {
        'config_draft_version': CONFIG_DRAFT_VERSION,
        'design_version': DESIGN_VERSION,
        'decision_reference': DECISION_REFERENCE,
        'ruling_reference': RULING_REFERENCE,
        'parameter_ruling_reference': PARAMETER_RULING_REFERENCE,
        'parameter_ruling_version': PARAMETER_RULING_VERSION,
        'is_runnable_config': False,
        'any_unit_runtime_eligible': False,
        'grid': list(RESIDUAL_RATE_GRID),
        'grid_label': STRESS_SETTING_LABEL,
        'runtime_axes': ['expected_doublet_rate'],
        'interpretation_only_axes': list(INTERPRETATION_ONLY_AXES),
        'fixed_per_capture_seed': list(FIXED_PER_CAPTURE_SEED),
        'invariance_checks': dict(INVARIANCE_CHECKS),
        'common_conditions': dict(COMMON_CONDITIONS),
        'unresolved_parameters': missing,
        'n_unresolved_parameters': len(missing),
        'ruled_parameters': sorted(METHODOLOGY_PARAMETER_RULING),
        'methodologically_approved_parameters': dict(METHODOLOGY_PARAMETER_RULING),
        'production_configuration_activated': bool(PRODUCTION_PARAMETER_RULING),
        'cohorts': cohorts,
        'what_must_hold_before_a_runnable_config': [
            'every parameter in REQUIRED_PARAMETER_NAMES has a recorded methodological '
            'ruling with a value, scope and invalidation condition',
            'C2: a physical capture key and raw-barcode namespace are established for '
            'each eligible Gate 2 unit from that cohort\'s own evidence',
            'C1: the raw-count source is shown to be un-normalised RNA counts by '
            'provenance, not by integrality alone',
            'C3: the complete cell-called/QC population for the declared stage is '
            'available and parent-state support is reviewed; a limitation label alone '
            'does not establish eligibility',
            'C4: threshold rule and effective-dimension condition fixed in advance, with '
            'failures recorded as unknown/NA',
            'effective-dimension feasibility is assessed on each eligible Gate 2 capture '
            'using the complete input that would actually be scored',
            'at least one capture has reviewed C1-C3 evidence, the exact config passes '
            'audit, and the reviewer records a scoped small-pilot authorization',
        ],
    }


if __name__ == '__main__':          # pragma: no cover - inspection helper
    print(json.dumps(config_draft(), indent=2, sort_keys=True))
