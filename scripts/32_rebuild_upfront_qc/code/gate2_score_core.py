"""Gate 2 scoring core: score one capture, validate what Scrublet returned.

Replaces `scrublet_one_channel`'s `(keep_mask, score, status, diag)` tuple. The
mask was the defect: a caller could read `keep=True` as "singlet" when the real
meaning was "we could not tell". This returns a structured result whose only
non-doublet assertion is `scrublet_noncandidate`, and candidate calls are tri-state so
"unknown" has a representation.

What is deliberately NOT here: any manual threshold, any top-r fallback, any
default rate. If automatic thresholding does not succeed, the capture is
`threshold_failed` and every cell in it is unknown. Scanpy warns rather than
raising in that case, which is exactly how the previous code recorded a whole
Wilk capture as singlets, so the return value is validated rather than trusted.
"""
import warnings
from dataclasses import dataclass, field

import numpy as np

from gate2_diagnostics import DISTRIBUTION_FIELDS, score_distribution_summary
from gate2_status import (CALL_PREDICTED_DOUBLET, CALL_SCRUBLET_NONCANDIDATE,
                          CAPTURE_RUNTIME_FAILED,
                          CAPTURE_SCORED, CAPTURE_THRESHOLD_FAILED, CAPTURE_TOO_SMALL,
                          CAPTURE_TO_UNIFORM_CALL, validate_capture_invariants)

#: PROJECT RULE (qc_checklist 2), not a reliability boundary: below this a
#: bimodal threshold cannot be fitted meaningfully, so the capture is skipped
#: and marked -- never passed through as singlets.
MIN_CELLS_FOR_SCRUBLET = 70

#: Computational limits, not scientific ones. 2**53 is where float64 stops
#: being able to distinguish consecutive integers, so "this value is an
#: integer" ceases to be evidence. Both are fixed properties of the arithmetic
#: and must never be tuned to cohort results or outcomes.
FLOAT_EXACT_INT_MAX = 2 ** 53
INT_COUNT_MAX = 2 ** 53

#: Every parameter here is PASSED to the Scrublet call below, so "recorded"
#: and "used" are the same set by construction. An earlier version listed
#: parameters that were never passed, which made the fingerprint claim a
#: binding it did not have.
#:
#: Values are scanpy 1.11.3's own defaults, pinned rather than changed. Pinning
#: them is itself protective: if a future scanpy changes a default, the
#: fingerprint no longer silently matches.
SCRUBLET_PARAMS = {
    'sim_doublet_ratio': 2.0,
    'stdev_doublet_rate': 0.02,
    'synthetic_doublet_umi_subsampling': 1.0,
    'knn_dist_metric': 'euclidean',
    'normalize_variance': True,
    'log_transform': False,
    'mean_center': True,
    'n_prin_comps': 30,
    'n_neighbors': None,
    'threshold': None,          # None = automatic detection, which is the point
}


class RawCountError(RuntimeError):
    """The input is not the unnormalized raw-count matrix Scrublet requires."""


def validate_raw_counts(mat, *, where='capture'):
    """Scrublet's simulated-doublet null is built by ADDING observed cells, so
    the input must be raw counts. A normalized, log1p-ed, scaled or
    DecontX-corrected matrix runs to completion and produces a
    well-formed-but-meaningless result, complete with a valid checkpoint
    fingerprint -- nothing downstream could tell.

    dtype is not the test: counts are commonly stored as float. The test is
    that the values ARE integers, finite, non-negative, and small enough that
    "is an integer" is still a provable statement.

    The bounds below are COMPUTATIONAL safety limits. They are fixed properties
    of float64 and int64 arithmetic and must never be tuned to a cohort's
    results or to any outcome.
    """
    data = getattr(mat, 'data', None)
    if data is None:
        data = np.asarray(mat).ravel()
    data = np.asarray(data)

    if data.dtype == bool:
        raise RawCountError(
            f'{where}: matrix is boolean -- presence/absence, not counts')
    if data.dtype.kind == 'c':
        raise RawCountError(f'{where}: matrix is complex -- not counts')
    if data.dtype.kind not in 'fiu':
        raise RawCountError(f'{where}: matrix dtype {data.dtype} is not a real numeric type')

    if data.size == 0:
        # Not a claim that the matrix is malformed -- an all-zero capture is a
        # legal count matrix. It simply cannot be scored: Scrublet has no
        # expression to build a simulated-doublet null from.
        raise RawCountError(
            f'{where}: matrix has no non-zero counts, so Scrublet cannot be run on it '
            '(this is a runnability limit, not a claim that the matrix is invalid)')

    if data.dtype.kind == 'f':
        if not bool(np.all(np.isfinite(data))):
            raise RawCountError(f'{where}: matrix contains NaN or inf -- not raw counts')
        if bool(np.any(data < 0)):
            raise RawCountError(
                f'{where}: matrix contains negative values (min {float(data.min()):.4g}) -- '
                'this is scaled or centred data, not raw counts')
        # Above 2**53 every float64 is an integer, so "integer-valued" stops
        # being evidence of anything: a normalized value could land there and
        # pass. Refuse rather than pretend the check still holds.
        if bool(np.any(data > FLOAT_EXACT_INT_MAX)):
            raise RawCountError(
                f'{where}: float matrix has values above 2**53 (max {float(data.max()):.4g}); '
                'integer-valuedness is not provable there, so the raw-count check cannot pass')
        if not bool(np.all(data == np.rint(data))):
            frac = float(np.mean(data != np.rint(data)))
            raise RawCountError(
                f'{where}: {100 * frac:.1f}% of stored values are non-integer -- this looks '
                'normalized, log1p-transformed or DecontX-corrected, not raw counts')
    else:
        if data.dtype.kind == 'i' and bool(np.any(data < 0)):
            raise RawCountError(
                f'{where}: integer matrix contains negative values (min {int(data.min())})')
        if bool(np.any(data.astype(object) > INT_COUNT_MAX)):
            raise RawCountError(
                f'{where}: integer counts exceed {INT_COUNT_MAX} (max {int(data.max())}); '
                'beyond this, downstream float64 conversion inside Scrublet and sparse '
                'reductions are no longer exact')


@dataclass
class CaptureResult:
    capture_status: str
    call_status: list
    #: tri-state per cell: 1 doublet, 0 not-doublet, -1 unknown. Never bool --
    #: a bool cannot say "unknown", and False reads as "examined, not a doublet".
    candidate: np.ndarray
    scores: np.ndarray
    diagnostics: dict = field(default_factory=dict)

    def n_predicted_doublets(self):
        return int(np.sum(self.candidate == 1))


def _uniform(n, capture_status, scores=None, diagnostics=None):
    """Every cell gets the same unknown call. Candidate is -1 (unknown), never 0."""
    call = CAPTURE_TO_UNIFORM_CALL[capture_status]
    diag = diagnostics or {}
    res = CaptureResult(
        capture_status=capture_status,
        call_status=[call] * n,
        candidate=np.full(n, -1, dtype=np.int8),
        scores=np.full(n, np.nan) if scores is None else np.asarray(scores, dtype=float),
        diagnostics=diag,
    )
    validate_capture_invariants(res.capture_status, res.call_status, res.candidate,
                                res.scores, threshold=diag.get('threshold'),
                                diagnostics=diag, where=f'fresh:{capture_status}')
    return res


def _safe_repr_array(arr):
    """Describe an array for an error message WITHOUT touching its values.

    np.unique / sorting / truthiness on an unvalidated array is itself a way to
    crash: a pandas nullable BooleanArray containing pd.NA raises
    "boolean value of NA is ambiguous" from inside np.unique, and that
    exception escaped the validator entirely. Only shape and dtype are read
    here, and even that is guarded.
    """
    try:
        return f'type={type(arr).__name__} dtype={getattr(arr, "dtype", "?")} ' \
               f'shape={getattr(arr, "shape", "?")}'
    except Exception:
        return '<undescribable object>'


def _extract_predictions(obs_value, n_cells):
    """Return (preds_bool_or_None, problem_or_None). Never raises.

    Accepts exactly two shapes of input:
      * a real boolean array (python bool / np.bool_);
      * an integer-dtype array whose values are strictly {0, 1}.

    Everything else -- float 0.0/1.0, strings, categoricals, object arrays,
    pandas nullable booleans, anything containing NA -- is refused. A bare
    `.astype(bool)` would map 0.5, 2 and -1 to True, inventing calls out of
    whatever the column happened to hold.
    """
    try:
        # Reject pandas extension arrays before converting: an ExtensionArray
        # holding pd.NA cannot be inspected safely by numpy at all.
        if hasattr(obs_value, 'dtype') and not isinstance(getattr(obs_value, 'dtype', None),
                                                          np.dtype):
            return None, (f'predicted_doublet is a pandas extension type '
                          f'({_safe_repr_array(obs_value)}); a nullable or categorical dtype '
                          'cannot express a strict boolean call')
        raw = np.asarray(obs_value)
    except Exception as exc:
        return None, f'predicted_doublet could not be read as an array ({type(exc).__name__})'

    try:
        if raw.ndim != 1:
            return None, f'predicted_doublet is not 1-D ({_safe_repr_array(raw)})'
        if raw.shape[0] != n_cells:
            return None, f'predicted_doublet length {raw.shape[0]} != {n_cells} cells'
        if raw.dtype == bool:
            return raw.astype(bool), None
        if raw.dtype.kind in 'iu':
            # Safe: an integer array cannot contain NA, so comparison is total.
            if bool(np.all((raw == 0) | (raw == 1))):
                return raw.astype(bool), None
            return None, (f'predicted_doublet is integer but not strictly 0/1 '
                          f'({_safe_repr_array(raw)})')
        return None, (f'predicted_doublet is not a strict boolean or 0/1 integer '
                      f'({_safe_repr_array(raw)})')
    except Exception as exc:
        return None, (f'predicted_doublet could not be validated '
                      f'({type(exc).__name__}: {exc})')


def _safe_float(value):
    """float() that never raises; returns nan for anything unusable."""
    try:
        out = float(value)
    except Exception:
        return float('nan')
    return out


def _validate_scrublet_output(adata, n_cells):
    """Return (scores, predictions, threshold, sim_scores, problems).

    NEVER raises. Everything the caller would otherwise have to trust is
    checked here, and every check is itself guarded -- the validator running on
    a malformed return value must classify it, not crash on it.
    """
    problems = []
    try:
        obs = adata.obs
    except Exception as exc:
        return None, None, float('nan'), None, [f'adata.obs unreadable ({type(exc).__name__})']

    scores = None
    try:
        if 'doublet_score' not in obs:
            problems.append('doublet_score missing from obs')
        else:
            cand = np.asarray(obs['doublet_score'], dtype=float)
            if cand.ndim != 1:
                problems.append(f'doublet_score is not 1-D ({_safe_repr_array(cand)})')
            elif cand.shape[0] != n_cells:
                problems.append(f'doublet_score length {cand.shape[0]} != {n_cells} cells')
            elif not bool(np.all(np.isfinite(cand))):
                problems.append('doublet_score contains non-finite values')
            else:
                scores = cand
    except Exception as exc:
        problems.append(f'doublet_score could not be validated ({type(exc).__name__}: {exc})')

    preds = None
    try:
        if 'predicted_doublet' not in obs:
            problems.append('predicted_doublet missing from obs')
        else:
            preds, prob = _extract_predictions(obs['predicted_doublet'], n_cells)
            if prob:
                problems.append(prob)
    except Exception as exc:
        problems.append(f'predicted_doublet could not be read ({type(exc).__name__}: {exc})')

    threshold = float('nan')
    sim_scores = None
    try:
        uns = adata.uns.get('scrublet', {}) or {}
        threshold = _safe_float(uns.get('threshold', None))
        if not np.isfinite(threshold):
            # Scanpy WARNS ("Failed to automatically identify doublet score
            # threshold") and returns; it does not raise. Treating that as
            # success is the defect that recorded a Wilk capture as all-singlet.
            problems.append('automatic doublet score threshold missing or non-finite')
        sim = uns.get('doublet_scores_sim', None)
        if sim is not None:
            arr = np.asarray(sim, dtype=float)
            if arr.size:
                sim_scores = arr
    except Exception as exc:
        problems.append(f'scrublet uns could not be read ({type(exc).__name__}: {exc})')

    return scores, preds, threshold, sim_scores, problems


def score_one_capture(mat_t, cell_ids, *, expected_doublet_rate, seed,
                      donor_proportions=None, platform='10x', scrublet_fn=None,
                      dump_simulated_scores_to=None):
    """Score ONE capture's FULL cell population (never a pre-extracted lineage).

    `mat_t` is cells x genes, unnormalized raw counts. `expected_doublet_rate`
    comes from an authorized runtime config -- there is no default and no
    derivation from cell counts here.

    `scrublet_fn` is injected so tests can force success, a NaN threshold, or a
    runtime exception without weakening any validation below.

    `dump_simulated_scores_to` is a DIAGNOSTIC-ONLY escape hatch, default OFF.
    When given a path it writes that capture's full simulated score vector as
    a separate evidence artifact, so a future pilot can verify the summaries
    against the raw distribution. Production runs must leave it None: the full
    vectors are roughly twice the cell count per capture per rate per seed, and
    are not needed once the summaries are validated.
    """
    n_cells = mat_t.shape[0]
    if n_cells == 0:
        raise RawCountError('capture has zero cells; nothing to score')
    validate_raw_counts(mat_t, where='scrublet input')
    base_diag = {
        'n_cells_scored': int(n_cells),
        'platform': platform,
        'expected_doublet_rate': float(expected_doublet_rate),
        'seed': int(seed),
        'n_donors_in_capture': (int(np.sum(np.asarray(donor_proportions) > 0))
                                if donor_proportions is not None else 1),
        # `threshold` is the OPERATIVE threshold -- the one the accepted calls
        # rest on. It is set only for a capture that actually scored. A value
        # Scrublet reported for a capture that then failed for some other
        # reason is evidence, not an operative threshold, and goes in
        # `threshold_observed` so a failed capture can never look scored.
        'threshold': None,
        'threshold_observed': None,
        'obs_score_median': None,
        'obs_score_p95': None,
        'sim_score_median': None,
        'sim_score_p95': None,
        'sim_scores_available': False,
        'n_predicted_doublets': 0,
        # populated from gate2_diagnostics while both score arrays are still in
        # memory; they cannot be rebuilt afterwards
        **{f: None for f in DISTRIBUTION_FIELDS},
        'warnings': [],
        'problems': [],
    }

    if n_cells < MIN_CELLS_FOR_SCRUBLET:
        d = dict(base_diag, problems=[f'capture has {n_cells} cells, '
                                      f'below MIN_CELLS_FOR_SCRUBLET={MIN_CELLS_FOR_SCRUBLET}'])
        return _uniform(n_cells, CAPTURE_TOO_SMALL, diagnostics=d)

    if scrublet_fn is None:
        import scanpy as sc
        scrublet_fn = sc.pp.scrublet

    from anndata import AnnData
    adata = AnnData(X=mat_t)

    caught = []
    try:
        with warnings.catch_warnings(record=True) as wlist:
            warnings.simplefilter('always')
            # **SCRUBLET_PARAMS, not a subset: the fingerprint records exactly
            # what is passed, so a changed parameter invalidates checkpoints.
            scrublet_fn(adata, expected_doublet_rate=float(expected_doublet_rate),
                        random_state=int(seed), **SCRUBLET_PARAMS)
            caught = [str(w.message) for w in wlist]
    except Exception as exc:
        d = dict(base_diag, problems=[f'{type(exc).__name__}: {exc}'], warnings=caught)
        return _uniform(n_cells, CAPTURE_RUNTIME_FAILED, diagnostics=d)

    try:
        scores, preds, threshold, sim_scores, problems = _validate_scrublet_output(adata, n_cells)
    except Exception as exc:
        # _validate_scrublet_output is written not to raise; this is the belt to
        # its braces. An unexpected failure IN THE VALIDATOR becomes a classified
        # runtime failure with every cell unknown -- never a singlet or doublet.
        d = dict(base_diag, warnings=caught,
                 problems=[f'validator failed: {type(exc).__name__}: {exc}'])
        return _uniform(n_cells, CAPTURE_RUNTIME_FAILED, diagnostics=d)
    diag = dict(base_diag, warnings=caught, problems=problems)
    if threshold is not None and np.isfinite(threshold):
        diag['threshold_observed'] = float(threshold)
        if not problems:
            diag['threshold'] = float(threshold)
    if scores is not None and scores.shape[0] == n_cells and np.all(np.isfinite(scores)):
        diag['obs_score_median'] = float(np.median(scores))
        diag['obs_score_p95'] = float(np.percentile(scores, 95))
    if sim_scores is not None and np.all(np.isfinite(sim_scores)):
        diag['sim_scores_available'] = True
        diag['sim_score_median'] = float(np.median(sim_scores))
        diag['sim_score_p95'] = float(np.percentile(sim_scores, 95))
    else:
        # Recorded explicitly rather than left ambiguous: diagnostics that are
        # absent must be visibly absent.
        diag['sim_scores_available'] = False

    # Distribution evidence is computed even when the capture failed: for a
    # threshold_failed capture the shape of the two distributions IS the
    # evidence. Threshold-based fields stay null there, handled inside.
    diag.update(score_distribution_summary(
        scores, sim_scores, threshold=diag.get('threshold')))

    if dump_simulated_scores_to and sim_scores is not None:
        np.save(dump_simulated_scores_to, np.asarray(sim_scores, dtype=float))

    if problems:
        # Scores are kept when they exist -- they remain useful evidence -- but
        # no cell is called, and candidate stays unknown.
        keep_scores = scores if (scores is not None and scores.shape[0] == n_cells) else None
        return _uniform(n_cells, CAPTURE_THRESHOLD_FAILED, scores=keep_scores, diagnostics=diag)

    call_status = [CALL_PREDICTED_DOUBLET if p else CALL_SCRUBLET_NONCANDIDATE
                   for p in preds]
    candidate = preds.astype(np.int8)
    diag['n_predicted_doublets'] = int(candidate.sum())
    validate_capture_invariants(CAPTURE_SCORED, call_status, candidate, scores,
                                threshold=diag.get('threshold'), diagnostics=diag,
                                where='fresh:scored')
    return CaptureResult(capture_status=CAPTURE_SCORED, call_status=call_status,
                         candidate=candidate, scores=scores, diagnostics=diag)
