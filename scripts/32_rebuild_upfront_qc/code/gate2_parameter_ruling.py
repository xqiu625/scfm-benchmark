"""Methodological approval, not a runnable configuration or a pilot permit.

All eleven values were adjudicated together on 2026-09-24. See the linked
decision for scope, singleton allowed ranges, rationale and invalidation.
Stage 31 is outside this approval. Production activation remains gated on
capture C1-C3 evidence, config audit and a recorded scoped pilot authorization.
"""
from types import MappingProxyType

PARAMETER_RULING_VERSION = 'gate2-parameter-ruling/1'
PARAMETER_RULING_REFERENCE = (
    'reports/manuscript/revision/methodology_review_2026-09-24/'
    'GATE2_PARAMETER_RULING.md'
)

METHODOLOGY_PARAMETER_RULING = MappingProxyType({
    'sim_doublet_ratio': 2.0,
    'synthetic_doublet_umi_subsampling': 1.0,
    'log_transform': False,
    'n_prin_comps': 30,
    # A frozen rule: round(0.5 * sqrt(post-filter observed n)), resolved before
    # each capture's scoring call and held constant across all prior settings.
    'n_neighbors': None,
    # Grid priors are fixed scenarios, not estimated prevalence distributions.
    # Zero removes the prior-SD term only; it does not assert zero score error.
    'stdev_doublet_rate': 0.0,
    'knn_dist_metric': 'euclidean',
    'normalize_variance': True,
    'mean_center': True,
    'use_approx_neighbors': False,
    'get_doublet_neighbor_parents': False,
})


def assert_ruling_version(expected):
    """Refuse if a consumer was written against a different ruling version."""
    if expected != PARAMETER_RULING_VERSION:
        raise SystemExit(
            f'ruling version mismatch: caller expects {expected!r}, this module is '
            f'{PARAMETER_RULING_VERSION!r}. Re-read {PARAMETER_RULING_REFERENCE} '
            'before proceeding.')
    return PARAMETER_RULING_VERSION


def approved_for(required_names):
    """Every required parameter, or a refusal naming what is absent or extra.

    Consumers call this instead of copying the eleven values. A second copy is
    the drift defect this module exists to prevent -- and one was created and
    deleted on 2026-09-24, because a grep of two files was mistaken for a search
    of the code directory.
    """
    missing = sorted(set(required_names) - set(METHODOLOGY_PARAMETER_RULING))
    if missing:
        raise SystemExit(f'the ruling does not cover {missing}')
    extra = sorted(set(METHODOLOGY_PARAMETER_RULING) - set(required_names))
    if extra:
        raise SystemExit(f'the ruling covers names this path does not use: {extra}')
    return dict(METHODOLOGY_PARAMETER_RULING)
