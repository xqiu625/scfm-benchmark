"""Frozen design for the Gate 2 scoring stability experiment.

This module is the *design*, not a run. It encodes the methodological ruling
recorded in `docs/decisions.md` (commit acf9692) so that the grid, the seed, the
per-cohort permitted use and the fail-closed rules cannot be varied per channel,
per run, or after results are seen.

WHAT THIS EXPERIMENT IS
-----------------------
An **algorithmic stress test** of Scrublet's calibration against the expected
rate. `expected_doublet_rate` feeds Scrublet's score calibration, so the scores
themselves move with it; the grid exists to measure *how much*.

WHAT IT IS NOT
--------------
Not a prevalence estimate. No grid value is any cohort's doublet rate, and a
result at any grid point may not be read back as one. Scanpy documents
`expected_doublet_rate` as the experiment's *estimated* doublet rate; using it
as a stress axis demotes it from an estimate to a calibration-sensitivity
parameter, and that demotion has to be stated wherever a number from this
experiment is reported.

This authorises **scoring only**. Deleting a cell needs a separate ruling; no
such ruling exists, and this module deliberately provides nothing that could be
read as one.

WHY A DIRECT RESIDUAL GRID AND NOT `base rate x removal efficiency`
-------------------------------------------------------------------
The mechanistic form

    r_res = r[s + (1-e)(1-s)] / ((1-r) + r[s + (1-e)(1-s)])

is retained for the reporting layer only. `r`, `s` and `e` are not jointly
identifiable from the available evidence, so sweeping them as a 2-D runtime grid
would manufacture many parameter combinations that look precise and are in fact
equivalent — and it would add researcher degrees of freedom to a question that
already has too many. The runtime axis is the residual rate directly.

Vendor recovered-stage tables are background and a sanity check on the grid's
span. They are **not** convertible into a post-Gate-1 residual rate for any
capture: the mapping needs capture-specific measured loading, a cohort recovered
count, Gate 1's differential retention of singlets vs doublets, and the full
selection mechanism of upstream filtering. None of those is available, and a
single-valued conversion would dress an unidentifiable quantity as a
measurement.
"""
import hashlib
import json
import math
import types

#: Bump when the frozen design changes. A run whose recorded version differs
#: from this one did not use this design.
DESIGN_VERSION = 'gate2-stress-design/1'

#: The decision this module implements. Cited in every diagnostics row so a
#: number can be traced to the ruling that permitted it.
DECISION_REFERENCE = 'docs/decisions.md :: Gate 2 scoring stability experiment (2026-09-23)'

# ---------------------------------------------------------------------------
# The grid
# ---------------------------------------------------------------------------

#: The residual expected-rate stress grid. FROZEN: identical for every cohort
#: permitted to score, so a difference between cohorts is a difference in the
#: data and never in the parameterisation.
#:
#: A tuple, not a list, and never rebuilt from a config file — a per-channel
#: grid would let a cohort be scored on its own axis, which is exactly the
#: researcher degree of freedom this design removes.
RESIDUAL_RATE_GRID = (0.005, 0.01, 0.02, 0.04, 0.08, 0.12, 0.16, 0.20)

#: One seed for the whole experiment. Per-capture seeds are derived from it
#: deterministically by the runner; nothing draws from an unseeded generator.
DESIGN_SEED = 20260923

# ---------------------------------------------------------------------------
# Per-cohort permitted use
# ---------------------------------------------------------------------------

#: What a cohort's scores may be used for. These are **use classes**, not
#: quality grades: `residual_diagnostic` is not a weaker version of
#: `standard_scoring`, it is a different claim about what the number means.
USE_STANDARD_SCORING = 'standard_scoring'
USE_RESIDUAL_DIAGNOSTIC = 'residual_diagnostic'
USE_RESIDUAL_SCORING_ONLY = 'residual_scoring_only'
USE_NOT_RUN = 'not_run'

PERMITTED_USES = frozenset({
    USE_STANDARD_SCORING,
    USE_RESIDUAL_DIAGNOSTIC,
    USE_RESIDUAL_SCORING_ONLY,
    USE_NOT_RUN,
})

#: The per-cohort ruling. A cohort absent from this table is NOT implicitly
#: permitted — `use_class_for()` raises, because a new cohort arriving without a
#: ruling is exactly the case that must stop the run.
_COHORT_USE_CLASS = {
    # Standard scoring sensitivity.
    'qiu':        USE_STANDARD_SCORING,
    'chongqing':  USE_STANDARD_SCORING,

    # Upstream filtering already removed an unknown, non-recoverable fraction.
    # Scores describe what survived that filtering, not the original population.
    'cillo':      USE_RESIDUAL_DIAGNOSTIC,
    'liu':        USE_RESIDUAL_DIAGNOSTIC,
    'combat':     USE_RESIDUAL_DIAGNOSTIC,

    # Ren's input is downstream of the authors' own Scrublet pass plus an
    # undefined manual removal, and its provenance is count-level, not
    # barcode-level. Scored for diagnosis only: never pooled with the standard
    # cohorts and never used to remove a cell.
    'ren':        USE_RESIDUAL_SCORING_ONLY,

    # Stephenson: no legal GEM-well/capture analysis unit — public sample/pool
    # identity does not map to a capture. Wilk: Seq-Well, where the 10x rate
    # model does not apply in principle.
    'stephenson': USE_NOT_RUN,
    'wilk':       USE_NOT_RUN,
}

#: Read-only at runtime. A plain dict could be reassigned after import --
#: `COHORT_USE_CLASS['stephenson'] = USE_STANDARD_SCORING` would then make
#: assert_may_score() let it through, and no source-level mutation test would
#: see it. design_manifest() copies back to a plain dict for JSON.
COHORT_USE_CLASS = types.MappingProxyType(_COHORT_USE_CLASS)

#: What each use class does and does not license, carried into the diagnostics
#: so the constraint travels with the number rather than living only here.
_USE_CLASS_MEANING = {
    USE_STANDARD_SCORING: (
        'Scores describe this capture under the stress grid. Still not a prevalence '
        'estimate: the grid is a calibration axis, not a measured rate.'),
    USE_RESIDUAL_DIAGNOSTIC: (
        'Upstream filtering removed an unknown fraction before these cells. Scores describe '
        'the residual population only and may not be reported as the original rate, nor '
        'compared across cohorts as prevalence.'),
    USE_RESIDUAL_SCORING_ONLY: (
        'Diagnostic only. Never pooled with standard-scoring cohorts, never used to remove a '
        'cell, and never reported as a rate for this cohort.'),
    USE_NOT_RUN: (
        'No scoring experiment may be designed for this cohort under this design. The '
        'obstruction is structural, not a matter of tuning.'),
}

USE_CLASS_MEANING = types.MappingProxyType(_USE_CLASS_MEANING)

#: Cohorts whose scores may never be pooled with another cohort's, whatever the
#: grid point. Kept separate from the use class so a future cohort can be
#: unpoolable without being diagnostic-only.
NEVER_POOLED = frozenset({'ren'})

# ---------------------------------------------------------------------------
# Fail-closed rules
# ---------------------------------------------------------------------------

#: Conditions that make a capture UNKNOWN at a grid point. Every one of these
#: has produced a silent false result in this repository at least once: scanpy
#: *warns* rather than raising when automatic thresholding fails, then returns
#: finite scores and an all-False prediction, so a run that only checks for an
#: exception records a failure as a clean negative.
#:
#: An all-False prediction is NOT itself in this list — it can be legitimate.
#: The threshold is what must be read.
UNKNOWN_CONDITIONS = (
    'threshold_absent',           # no threshold in the capture's own metadata
    'threshold_non_finite',       # NaN/inf threshold
    'scrublet_raised',            # the call raised
    'calls_not_strict_boolean',   # not a bool or 0/1 int array of the right length
    'calls_contain_na',           # a single NA disqualifies the whole capture
    'scores_not_finite',          # NaN/inf scores
    'scores_wrong_shape',         # wrong length or not 1-D
    # Two checks a draft runner merely RECORDED while still letting the grid
    # point be `examined`. Recording a check and ignoring its result is not a
    # check, so each now has a status a run cannot bypass.
    'identity_mismatch',          # obs_names differ from the order frozen at
                                  # preparation, so the Jaccard's positional
                                  # alignment no longer holds and the point may
                                  # not contribute to it
    'parameter_echo_mismatch',    # scanpy reports having used a value other
                                  # than the one requested, so the row's
                                  # provenance does not describe the run
)

#: Layer 2. RECORDED AT EVERY GRID POINT, NO PASS/FAIL CUTOFF.
#:
#: These describe method *reliability*, which this repository keeps separate
#: from mechanical validity above. A first version of this module listed
#: `simulation_unimodal` among the mechanical UNKNOWN conditions -- but no
#: modality detector, statistic, threshold or implementation is defined
#: anywhere, so a runner would have been free to decide case by case when to
#: reject a result, bypassing the `method_unreliable` layer whose cutoff has
#: not been adjudicated.
#:
#: They are diagnostics until a pre-specified, testable, methodologically
#: reviewed cutoff exists. None of them may fail a capture today.
RELIABILITY_DIAGNOSTICS_REQUIRED = (
    'simulated_score_histogram',   # incl. modality; no cutoff selected
    'observed_score_histogram',
    'threshold_placement',         # where the threshold sits in both distributions
    'seed_stability',              # call agreement across reruns
    'tail_enrichment',
)

#: CAPTURE-level statuses, distinct from the per-grid-point mechanical ones
#: above. A capture that never reached the grid did not fail eight times; it
#: failed once, before the grid. Frozen here rather than invented in the runner:
#: `capture_preparation_failed` was previously written straight into a row's
#: `mechanical_status`, which added an unfrozen value to a frozen enum.
CAPTURE_PREPARATION_FAILED = 'capture_preparation_failed'
CAPTURE_STATUSES = (CAPTURE_PREPARATION_FAILED,)

#: Every value `mechanical_status` may take. Checked by the runner so a new
#: status cannot be introduced without appearing here first.
ALL_STATUSES = ('examined',) + UNKNOWN_CONDITIONS + CAPTURE_STATUSES


def assert_status_is_frozen(status):
    """Raise unless `status` is one of the frozen values."""
    if status not in ALL_STATUSES:
        raise Gate2DesignError(
            f'status {status!r} is not one of the frozen statuses {ALL_STATUSES!r}. '
            'A run may not introduce a status the design has not declared.')
    return status


#: Never auto-filled. `False` asserts *examined and not a doublet*; an unknown
#: capture asserts nothing. This constant exists so the prohibition is
#: greppable, not only documented.
UNKNOWN_FILL_VALUE = None


def freeze_before_results():
    """The claims that must hold before any result is looked at.

    Returned as data so a run can record them, and so a reader can see that the
    list is short and specific rather than a paragraph of intent.
    """
    return (
        'grid, seed, preprocessing, failure states and comparison metrics are fixed '
        'before any result is inspected',
        'every grid point retains its full threshold, score distribution, call count '
        'and failure reason',
        'a non-finite or absent threshold, or an interface error, is recorded UNKNOWN '
        'and never filled with False',
        'reliability diagnostics -- including simulated-distribution modality -- are '
        'recorded at every grid point and carry NO pass/fail cutoff; none of them fails '
        'a capture until a cutoff is adjudicated',
        'no grid point is chosen because its biological result, mortality association '
        'or deletion fraction looks better',
        'if the candidate set moves materially across the grid, the finding is '
        '"the method is unstable here", not a choice of rate',
        'this design authorises scoring only; removing a cell requires a separate ruling',
    )


# ---------------------------------------------------------------------------
# Accessors -- every one fails closed
# ---------------------------------------------------------------------------

class Gate2DesignError(RuntimeError):
    """A violation of the frozen design. Always fatal, never downgraded."""


def use_class_for(cohort):
    """The permitted use for a cohort. An unruled cohort raises.

    Absence is not permission: a cohort that appears without a ruling is the
    case that must stop the run, not the case that gets a default.
    """
    key = (cohort or '').strip().lower()
    if key not in COHORT_USE_CLASS:
        raise Gate2DesignError(
            f'cohort {cohort!r} has no recorded use class. A cohort may not be scored on a '
            f'default; add an explicit ruling to COHORT_USE_CLASS and cite the decision. '
            f'Ruled cohorts: {sorted(COHORT_USE_CLASS)}')
    return COHORT_USE_CLASS[key]


def assert_may_score(cohort):
    """Raise unless this cohort may be scored at all."""
    use = use_class_for(cohort)
    if use == USE_NOT_RUN:
        raise Gate2DesignError(
            f'cohort {cohort!r} is {USE_NOT_RUN}: {USE_CLASS_MEANING[USE_NOT_RUN]}')
    return use


def assert_grid_is_frozen(grid):
    """Raise unless `grid` is exactly the frozen grid.

    Compared elementwise against the tuple rather than by set or by length, so a
    reordering, a substitution or an extra point is all refused.
    """
    if tuple(grid) != RESIDUAL_RATE_GRID:
        raise Gate2DesignError(
            f'grid {tuple(grid)!r} is not the frozen grid {RESIDUAL_RATE_GRID!r}. '
            'The grid is identical for every cohort by design; a per-cohort grid would '
            'make a cohort difference indistinguishable from a parameterisation difference.')
    return RESIDUAL_RATE_GRID


def capture_seed(cohort, channel_id):
    """A deterministic per-capture seed. The RATE IS NOT AN INPUT.

    Common random numbers: every grid point for one capture draws the same
    simulated doublets, so a difference in the call set across the grid is the
    calibration effect of `expected_doublet_rate` and not Monte Carlo variation.

    A first version folded the rate into this payload, and its test asserted
    eight distinct seeds for the eight grid points of one capture. That would
    have made the sensitivity experiment uninterpretable: every call-set
    difference would have carried both effects at once, with no way to separate
    them.

    Different captures still get different streams, so two captures never share
    one by accident, and a rerun of any capture reproduces exactly.
    """
    if cohort is None or not str(cohort).strip():
        raise Gate2DesignError('capture_seed: cohort is required')
    if channel_id is None or not str(channel_id).strip():
        raise Gate2DesignError('capture_seed: channel_id is required')
    payload = json.dumps(
        [DESIGN_VERSION, DESIGN_SEED, str(cohort).strip().lower(), str(channel_id).strip()],
        separators=(',', ':'), sort_keys=False, ensure_ascii=True)
    digest = hashlib.sha256(payload.encode('utf-8')).digest()
    return int.from_bytes(digest[:8], 'big') % (2 ** 32 - 1)


def assert_rate_on_grid(rate):
    """Raise unless `rate` is one of the frozen grid values.

    Exact membership, not a tolerance: a rate that is nearly a grid point came
    from somewhere other than the grid.
    """
    if rate is None:
        raise Gate2DesignError('expected_doublet_rate is required, got None')
    if isinstance(rate, bool):
        raise Gate2DesignError(
            f'expected_doublet_rate={rate!r} is a bool; True would otherwise compare equal '
            'to 1.0 and a rate must never arrive as a flag')
    if not isinstance(rate, (int, float)):
        raise Gate2DesignError(
            f'expected_doublet_rate={rate!r} is {type(rate).__name__}, not a number')
    val = float(rate)
    if not math.isfinite(val):
        raise Gate2DesignError(f'expected_doublet_rate={rate!r} is not finite')
    if val not in RESIDUAL_RATE_GRID:
        raise Gate2DesignError(
            f'expected_doublet_rate={rate!r} is not on the frozen grid {RESIDUAL_RATE_GRID!r}. '
            'Runtime rates come from the grid and nowhere else.')
    return val


def design_manifest():
    """The frozen design as a recordable dict, for the run manifest.

    Emitted so that a result file states which design produced it, and so that a
    later reader can compare rather than trust.
    """
    return {
        'design_version': DESIGN_VERSION,
        'decision_reference': DECISION_REFERENCE,
        'residual_rate_grid': list(RESIDUAL_RATE_GRID),
        'grid_semantics': 'algorithmic stress-test axis; NOT a prevalence estimate',
        'design_seed': DESIGN_SEED,
        'cohort_use_class': dict(sorted(COHORT_USE_CLASS.items())),  # plain dict for JSON
        'never_pooled': sorted(NEVER_POOLED),
        'unknown_conditions': list(UNKNOWN_CONDITIONS),
        'capture_statuses': list(CAPTURE_STATUSES),
        'all_statuses': list(ALL_STATUSES),
        'reliability_diagnostics_required': list(RELIABILITY_DIAGNOSTICS_REQUIRED),
        'reliability_cutoffs': 'none adjudicated; diagnostics are recorded, never fail a capture',
        'seeding': 'common random numbers -- one seed per capture, shared by every grid point',
        'authorises': 'scoring only; cell deletion requires a separate ruling',
        'frozen_before_results': list(freeze_before_results()),
    }
