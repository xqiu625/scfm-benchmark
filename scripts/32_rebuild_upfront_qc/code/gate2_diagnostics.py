"""Gate 2 score-distribution diagnostics.

These are DIAGNOSTIC EVIDENCE, not a reliability criterion. Nothing here
decides whether a capture is trustworthy, nothing here assigns
`method_unreliable`, and no value produced here authorizes removing a cell.

They exist because the metrics needed to judge separability cannot be
reconstructed after a run: Scrublet's simulated scores live only inside the
call, and the previous implementation kept a median and a p95 before
discarding them. Everything below is therefore computed while both arrays are
still in memory.

WHAT THE NUMBERS DO NOT MEAN
  * `observed_simulated_auc` measures how separable the two score
    distributions are. It is not a doublet detection accuracy, and a high
    value does not authorize deletion.
  * `simulated_recall_at_threshold` is recall against SIMULATED doublets --
    synthetic profiles built by summing observed cells. It is an upper-bound-ish
    proxy and is NOT the true doublet recall, which is unmeasurable here.
  * percentiles and histograms summarize the distributions; they do not
    reconstruct them. A claim needing the full distribution needs the
    diagnostic-only capture dump, not these.

AUC DIRECTION, fixed and documented:
    simulated is the POSITIVE class
    AUC = P(sim > obs) + 0.5 * P(sim == obs)
So AUC = 1 means simulated doublets score strictly above every observed cell
(perfect separation), 0.5 means the distributions are indistinguishable, and
0 means the direction is inverted -- which would itself be a red flag.
"""
import math

import numpy as np

#: Persisted alongside the fields so a future change to their definitions
#: cannot be silently absorbed by an old checkpoint.
DIAGNOSTICS_SCHEMA_VERSION = 'gate2-diagnostics/1'

#: 101 points, p0..p100 inclusive.
PERCENTILE_POINTS = tuple(range(101))

#: Shared bins for the two distributions, so they can be compared directly.
HISTOGRAM_N_BINS = 50

DISTRIBUTION_FIELDS = (
    'diagnostics_schema_version',
    'n_observed', 'n_simulated',
    'observed_score_min', 'observed_score_max',
    'simulated_score_min', 'simulated_score_max',
    'simulated_recall_at_threshold',
    'observed_tail_fraction_at_threshold',
    'simulated_tail_fraction_at_threshold',
    'observed_simulated_auc',
    'observed_score_percentiles_p0_to_p100',
    'simulated_score_percentiles_p0_to_p100',
    'shared_histogram_bin_edges',
    'observed_histogram_counts',
    'simulated_histogram_counts',
    'histogram_range_degenerate',
    'distribution_problems',
)


def _plain(x):
    """numpy -> plain Python, so the JSON in the checkpoint holds no numpy
    objects and round-trips identically."""
    if isinstance(x, np.generic):
        return x.item()
    if isinstance(x, np.ndarray):
        return [_plain(v) for v in x]
    if isinstance(x, (list, tuple)):
        return [_plain(v) for v in x]
    return x


def _finite_1d(arr):
    """Return a finite 1-D float array, or None with a reason."""
    if arr is None:
        return None, 'absent'
    a = np.asarray(arr, dtype=float).ravel()
    if a.size == 0:
        return None, 'empty'
    if not bool(np.all(np.isfinite(a))):
        return None, 'contains non-finite values'
    return a, None


def observed_simulated_auc(observed, simulated):
    """P(sim > obs) + 0.5 * P(sim == obs), simulated as the positive class.

    Rank-based (Mann-Whitney U), so ties get exactly half weight and the result
    is exact rather than sampled.
    """
    n_obs, n_sim = observed.size, simulated.size
    if n_obs == 0 or n_sim == 0:
        return None
    combined = np.concatenate([simulated, observed])
    order = np.argsort(combined, kind='mergesort')
    ranks = np.empty(combined.size, dtype=float)
    sorted_vals = combined[order]
    i = 0
    while i < sorted_vals.size:
        j = i
        while j + 1 < sorted_vals.size and sorted_vals[j + 1] == sorted_vals[i]:
            j += 1
        # average rank over the tie group (ranks are 1-based)
        ranks[order[i:j + 1]] = 0.5 * ((i + 1) + (j + 1))
        i = j + 1
    rank_sum_sim = float(ranks[:n_sim].sum())
    u_sim = rank_sum_sim - n_sim * (n_sim + 1) / 2.0
    return float(u_sim / (n_sim * n_obs))


def _bin_edges(observed, simulated):
    """Shared edges over the JOINT finite range of both distributions.

    Scrublet scores are conventionally in [0, 1] but that is not enforced
    anywhere in this pipeline, so the range is taken from the data rather than
    assumed -- clipping to [0, 1] would silently hide an out-of-range score,
    which is exactly the kind of thing these diagnostics exist to reveal.
    """
    lo = float(min(observed.min(), simulated.min()))
    hi = float(max(observed.max(), simulated.max()))
    degenerate = not (hi > lo)
    if degenerate:
        # A constant distribution has no width. Give it a defined, symmetric
        # unit window rather than producing zero-width bins or dividing by zero.
        lo, hi = lo - 0.5, lo + 0.5
    return np.linspace(lo, hi, HISTOGRAM_N_BINS + 1), degenerate


def score_distribution_summary(observed, simulated, threshold=None):
    """Summarize both score distributions. Never raises.

    `threshold`-based fields are computed ONLY when the threshold is finite.
    Otherwise they are None -- never 0, which would read as a measured absence
    of tail mass rather than an unavailable measurement.

    Distribution-shape fields (AUC, percentiles, histograms) are computed
    whenever both arrays are usable, INCLUDING for a threshold_failed capture,
    where they are the main evidence available.
    """
    out = {f: None for f in DISTRIBUTION_FIELDS}
    out['diagnostics_schema_version'] = DIAGNOSTICS_SCHEMA_VERSION
    problems = []

    obs, obs_why = _finite_1d(observed)
    sim, sim_why = _finite_1d(simulated)
    if obs_why:
        problems.append(f'observed scores {obs_why}')
    if sim_why:
        problems.append(f'simulated scores {sim_why}')

    if obs is not None:
        out['n_observed'] = int(obs.size)
        out['observed_score_min'] = float(obs.min())
        out['observed_score_max'] = float(obs.max())
        out['observed_score_percentiles_p0_to_p100'] = _plain(
            np.percentile(obs, PERCENTILE_POINTS))
    if sim is not None:
        out['n_simulated'] = int(sim.size)
        out['simulated_score_min'] = float(sim.min())
        out['simulated_score_max'] = float(sim.max())
        out['simulated_score_percentiles_p0_to_p100'] = _plain(
            np.percentile(sim, PERCENTILE_POINTS))

    if obs is not None and sim is not None:
        out['observed_simulated_auc'] = observed_simulated_auc(obs, sim)
        edges, degenerate = _bin_edges(obs, sim)
        out['shared_histogram_bin_edges'] = _plain(edges)
        out['histogram_range_degenerate'] = bool(degenerate)
        out['observed_histogram_counts'] = _plain(np.histogram(obs, bins=edges)[0])
        out['simulated_histogram_counts'] = _plain(np.histogram(sim, bins=edges)[0])

    thr_ok = threshold is not None and isinstance(threshold, (int, float)) \
        and math.isfinite(float(threshold))
    if thr_ok:
        t = float(threshold)
        # `>=` matches the convention that a cell AT the threshold is called.
        if obs is not None:
            out['observed_tail_fraction_at_threshold'] = float(np.mean(obs >= t))
        if sim is not None:
            recall = float(np.mean(sim >= t))
            out['simulated_recall_at_threshold'] = recall
            out['simulated_tail_fraction_at_threshold'] = recall
    elif threshold is not None:
        problems.append('threshold is not finite; threshold-based metrics are null')

    out['distribution_problems'] = problems
    return out
