"""Gate 2 call-status model.

Two DIFFERENT things were previously carried in one `status` column, which is
how a capture whose automatic thresholding failed ended up recorded as a
capture full of singlets:

  CALL status     per CELL   -- what we concluded about this cell
  CAPTURE status  per CAPTURE -- how the scoring attempt executed

They are kept in separate columns here and must never be merged again.

The load-bearing rule: **a keep mask is not evidence of singlet status.** Only
`scrublet_noncandidate` asserts "this cell was examined by a working detector
and not called a doublet" -- nothing more. Every other non-doublet state means
"we do not know", and must be carried as a nullable call, never as `False`.

TWO LAYERS, and this module implements only the first:

    mechanical candidate   call_status == CALL_PREDICTED_DOUBLET
    final deletion authorization
                           mechanical candidate
                           + the cohort is eligible for this scoring at all
                           + an approved runtime rate
                           + the capture passed a reliability criterion
                           + a sensitivity analysis supports the call

**No final-authorization state exists in this codebase, and none should be
created until every one of those is settled.** `predicted_doublet` is a
NECESSARY but NOT SUFFICIENT condition for removing a cell. Nothing here
authorizes deletion, and no downstream consumer may reconstruct an
authorization from a boolean mask.
"""

# ---- per-cell call status -------------------------------------------------

CALL_PREDICTED_DOUBLET = 'predicted_doublet'
#: The capture scored successfully, the automatic threshold was valid, and this
#: cell did not exceed it. That is ALL it means.
#:
#: It does NOT assert the cell is a true singlet, does not confer any retention
#: authorization, and does not exclude a homotypic doublet -- Scrublet cannot
#: detect those by construction.
#:
#: The former name `auto_singlet` was withdrawn because "singlet" already means
#: something else in this repository: demuxlet/genotype SNG, i.e. one DONOR, not
#: one cell. A Liu barcode can be a demuxlet singlet AND a same-donor doublet at
#: the same time. After this rename, "singlet" in this codebase means the
#: demuxlet sense only.
CALL_SCRUBLET_NONCANDIDATE = 'scrublet_noncandidate'
CALL_TOO_SMALL = 'too_small'
CALL_THRESHOLD_FAILED = 'threshold_failed'
CALL_METHOD_UNRELIABLE = 'method_unreliable'
CALL_RUNTIME_FAILED = 'runtime_failed'

CALL_STATUSES = (
    CALL_PREDICTED_DOUBLET,
    CALL_SCRUBLET_NONCANDIDATE,
    CALL_TOO_SMALL,
    CALL_THRESHOLD_FAILED,
    CALL_METHOD_UNRELIABLE,
    CALL_RUNTIME_FAILED,
)

#: Statuses that make a cell a MECHANICAL CANDIDATE -- the algorithm's raw
#: output. This is not an authorization to remove anything. See the two-layer
#: note in the module docstring.
MECHANICAL_DOUBLET_CANDIDATE_STATUSES = frozenset({CALL_PREDICTED_DOUBLET})

#: Statuses that assert nothing about the cell. Their candidate-doublet value
#: must be null (pd.NA), never False -- False reads as "examined, not a
#: doublet", which is exactly the claim these states cannot support.
UNKNOWN_CALL_STATUSES = frozenset({
    CALL_TOO_SMALL,
    CALL_THRESHOLD_FAILED,
    CALL_METHOD_UNRELIABLE,
    CALL_RUNTIME_FAILED,
})

#: Assigned only by an explicit, separately recorded policy decision. No
#: reliability threshold has been chosen, so nothing in this codebase may
#: assign it; it exists so that decision has a place to land.
POLICY_ASSIGNED_ONLY = frozenset({CALL_METHOD_UNRELIABLE})


# ---- per-capture execution status -----------------------------------------

CAPTURE_SCORED = 'scored'
CAPTURE_TOO_SMALL = 'too_small'
CAPTURE_THRESHOLD_FAILED = 'threshold_failed'
CAPTURE_RUNTIME_FAILED = 'runtime_failed'

CAPTURE_STATUSES = (
    CAPTURE_SCORED,
    CAPTURE_TOO_SMALL,
    CAPTURE_THRESHOLD_FAILED,
    CAPTURE_RUNTIME_FAILED,
)

#: Capture execution status -> the call status every cell in it receives when
#: the capture did not produce usable calls.
CAPTURE_TO_UNIFORM_CALL = {
    CAPTURE_TOO_SMALL: CALL_TOO_SMALL,
    CAPTURE_THRESHOLD_FAILED: CALL_THRESHOLD_FAILED,
    CAPTURE_RUNTIME_FAILED: CALL_RUNTIME_FAILED,
}


def validate_call_statuses(values):
    """Raise on any value outside the enum. Used on checkpoint restore, where a
    stale or hand-edited file could otherwise reintroduce free-text states."""
    bad = sorted({str(v) for v in values} - set(CALL_STATUSES))
    if bad:
        raise ValueError(f'call_status values outside the enum: {bad}')


def validate_capture_status(value):
    if value not in CAPTURE_STATUSES:
        raise ValueError(f'capture_status {value!r} outside {CAPTURE_STATUSES}')
    return value


def mechanical_candidate_mask(call_status_values):
    """Mark the cells Scrublet mechanically called doublets.

    This is NOT a deletion decision and NOT an authorization. It answers only
    "did the detector call this cell a doublet", which is a necessary but not
    sufficient condition for removal. Cohort eligibility, an approved runtime
    rate, a capture reliability criterion and a sensitivity analysis are all
    still required, and none of them exists yet.

    The previous name was `deletion_mask`, returning `DELETION_ELIGIBLE`. Both
    asserted an authorization this layer cannot give. The rename is a
    deliberate BREAKING semantic change: the old names are gone rather than
    aliased, so any caller still using them fails loudly instead of quietly
    receiving a candidate mask it may read as permission.
    """
    return [str(v) in MECHANICAL_DOUBLET_CANDIDATE_STATUSES for v in call_status_values]


class CaptureInvariantError(RuntimeError):
    """A capture's status fields contradict each other. Always fatal."""


#: What each capture execution status permits. Validating fields one at a time
#: lets contradictory combinations through -- a checkpoint could claim
#: capture_status='too_small' while every cell said 'scrublet_noncandidate', and every
#: individual field would still be "valid". This table is the contract; it is
#: enforced on the fresh path, on save, and on restore, by the same function.
CAPTURE_INVARIANTS = {
    CAPTURE_SCORED: {
        'calls': frozenset({CALL_PREDICTED_DOUBLET, CALL_SCRUBLET_NONCANDIDATE}),
        'candidate': frozenset({0, 1}),
        'scores_finite': True,
        'threshold_finite': True,
    },
    CAPTURE_THRESHOLD_FAILED: {
        'calls': frozenset({CALL_THRESHOLD_FAILED}),
        'candidate': frozenset({-1}),
        # scores are retained as evidence when Scrublet produced them, so this
        # is the one status where finite and NaN are both legitimate -- but not
        # a mixture, which would mean a partially applied result.
        'scores_finite': None,
        'threshold_finite': False,
    },
    CAPTURE_TOO_SMALL: {
        'calls': frozenset({CALL_TOO_SMALL}),
        'candidate': frozenset({-1}),
        'scores_finite': False,
        'threshold_finite': False,
    },
    CAPTURE_RUNTIME_FAILED: {
        'calls': frozenset({CALL_RUNTIME_FAILED}),
        'candidate': frozenset({-1}),
        'scores_finite': False,
        'threshold_finite': False,
    },
}

#: call_status -> the only candidate value it may carry. Enforced per cell, so
#: `predicted_doublet` with candidate 0 (or the reverse) is impossible.
CALL_TO_CANDIDATE = {
    CALL_PREDICTED_DOUBLET: 1,
    CALL_SCRUBLET_NONCANDIDATE: 0,
    CALL_TOO_SMALL: -1,
    CALL_THRESHOLD_FAILED: -1,
    CALL_METHOD_UNRELIABLE: -1,
    CALL_RUNTIME_FAILED: -1,
}


def validate_capture_invariants(capture_status, call_status, candidate, scores,
                                threshold=None, diagnostics=None, *, where='capture'):
    """Enforce the whole contract at once. Called from the fresh path, from
    checkpoint save and from checkpoint restore, so all three agree by
    construction rather than by three separate partial checks."""
    import math

    import numpy as np

    validate_capture_status(capture_status)
    validate_call_statuses(call_status)
    rules = CAPTURE_INVARIANTS[capture_status]

    n = len(call_status)
    if n == 0:
        raise CaptureInvariantError(
            f'{where}: capture has zero cells; every check below would pass vacuously, '
            'so an empty capture is refused rather than recorded as scored')
    if not (len(candidate) == len(scores) == n):
        raise CaptureInvariantError(
            f'{where}: lengths disagree -- call_status={n} candidate={len(candidate)} '
            f'scores={len(scores)}')

    bad_calls = sorted(set(map(str, call_status)) - rules['calls'])
    if bad_calls:
        raise CaptureInvariantError(
            f'{where}: capture_status={capture_status!r} permits {sorted(rules["calls"])} '
            f'but found {bad_calls}')

    # Checked on the RAW object: a prior version did {int(c) for c in candidate},
    # which silently accepted the strings "0"/"1" and would coerce 2.5. dtype is
    # part of the contract, not an implementation detail.
    cand_arr = candidate if isinstance(candidate, np.ndarray) else None
    if cand_arr is None:
        raise CaptureInvariantError(
            f'{where}: candidate must be a numpy array, got {type(candidate).__name__}')
    if cand_arr.ndim != 1:
        raise CaptureInvariantError(f'{where}: candidate must be 1-D, got ndim={cand_arr.ndim}')
    if cand_arr.dtype != np.int8:
        raise CaptureInvariantError(
            f'{where}: candidate dtype must be int8, got {cand_arr.dtype}')

    bad_cand = sorted(set(cand_arr.tolist()) - rules['candidate'])
    if bad_cand:
        raise CaptureInvariantError(
            f'{where}: capture_status={capture_status!r} permits candidate '
            f'{sorted(rules["candidate"])} but found {bad_cand}')

    for i, (cs, cand) in enumerate(zip(call_status, candidate)):
        want = CALL_TO_CANDIDATE[str(cs)]
        if int(cand) != want:
            raise CaptureInvariantError(
                f'{where}: cell {i} has call_status={cs!r} but candidate={int(cand)} '
                f'(must be {want})')

    finite = [bool(math.isfinite(float(s))) for s in scores]
    if rules['scores_finite'] is True and not all(finite):
        raise CaptureInvariantError(f'{where}: capture_status={capture_status!r} requires all '
                                    'scores finite')
    if rules['scores_finite'] is False and any(finite):
        raise CaptureInvariantError(f'{where}: capture_status={capture_status!r} requires all '
                                    'scores NA/non-finite')
    if rules['scores_finite'] is None and any(finite) and not all(finite):
        raise CaptureInvariantError(
            f'{where}: capture_status={capture_status!r} scores are a mixture of finite and '
            'non-finite, which would mean a partially applied result')

    thr_finite = threshold is not None and math.isfinite(float(threshold))
    if rules['threshold_finite'] and not thr_finite:
        raise CaptureInvariantError(
            f'{where}: capture_status={capture_status!r} requires a finite threshold, got '
            f'{threshold!r}')
    if not rules['threshold_finite'] and thr_finite:
        raise CaptureInvariantError(
            f'{where}: capture_status={capture_status!r} must not carry a finite threshold, got '
            f'{threshold!r}')

    if diagnostics is not None and 'n_predicted_doublets' in diagnostics:
        claimed = int(diagnostics['n_predicted_doublets'])
        actual = int(sum(1 for c in candidate if int(c) == 1))
        if claimed != actual:
            raise CaptureInvariantError(
                f'{where}: diagnostics claim {claimed} predicted doublets but the array has '
                f'{actual}')
