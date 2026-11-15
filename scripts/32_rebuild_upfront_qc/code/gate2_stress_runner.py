"""Runner draft for the Gate 2 scoring-stability experiment. NOTHING RUNS HERE.

This module wires `gate2_stress_design` to scanpy. It is submitted for static
audit: there is no `__main__`, no path constant, no output directory, no config
loader and no job script. It cannot be invoked against cohort data as it stands.

Every DESIGN constant comes from the design module by import. Nothing is copied:
a copied grid can drift from the ruling, and a drifted copy is
indistinguishable from a deliberate change.

PARAMETER RULING -- PRODUCTION CONFIGURATION NOT ACTIVATED
------------------------------------------------
All eleven parameters have a methodological ruling in gate2_parameter_ruling.
`PRODUCTION_PARAMETER_RULING` stays empty until capture eligibility, config
audit and scoped pilot authorization are recorded. The default entry point
therefore still raises before copying a matrix. Explicit injected rulings are
for technical fixtures; supplying one does not authorize a cohort run.

COMMON RANDOM NUMBERS -- STRUCTURAL, NOT ASSERTED
-------------------------------------------------
scanpy exposes `scrublet_simulate_doublets()` separately from `pp.scrublet()`,
and the simulator does **not** take `expected_doublet_rate`:

    scrublet_simulate_doublets(adata, layer, sim_doublet_ratio,
                               synthetic_doublet_umi_subsampling, random_seed)

so the simulated doublets can be built ONCE per capture and reused across the
whole grid by passing `adata_sim=` to every call. `simulated_checksum` verifies
that simulation reuse. PCA and neighbours are still recomputed in each wrapper
call: the intended calibration-only comparison additionally requires direct
representation/graph identity evidence. Equal simulation hashes and seeds do
not establish that remaining config-audit condition.

The cost is that supplying `adata_sim` makes `pp.scrublet` SKIP its internal
preprocessing -- it only preprocesses when `adata_sim is None`. So this module
reproduces scanpy's own preprocessing, in scanpy's order, using scanpy's
functions. `prepare_capture()` is that reproduction and nothing else.

Measured on the pinned 1.11.3 (`inspect.getsource(sc.pp.scrublet)`): the
`log_transform`, `sim_doublet_ratio` and `synthetic_doublet_umi_subsampling`
kwargs are read only inside that `adata_sim is None` branch. Passing them to the
scoring call alongside `adata_sim` would be silently ignored while the row
claimed they had been applied, which is why they are grouped as SIMULATION and
PREPROCESSING parameters and applied here instead.

WHAT THIS RUNNER MAY NOT DO
---------------------------
No deletion mask. No retained matrix. No overwriting. No "best rate".

The entry point takes a matrix and an identifier, never a patient table, so it
never itself loads outcome. It cannot make a claim about what the caller loaded
before handing over the AnnData -- that is outside this module. What it does
guarantee is that no unapproved metadata column is propagated onward:
`minimise_input()` rebuilds the object from the matrix and cell identity, and
column NAMES are recorded without any value being read.
"""
import hashlib

import numpy as np

from gate2_diagnostics import DIAGNOSTICS_SCHEMA_VERSION, score_distribution_summary
from gate2_stress_design import (CAPTURE_PREPARATION_FAILED, DECISION_REFERENCE,
                                 DESIGN_VERSION, NEVER_POOLED,
                                 RELIABILITY_DIAGNOSTICS_REQUIRED, RESIDUAL_RATE_GRID,
                                 UNKNOWN_FILL_VALUE, Gate2DesignError, assert_grid_is_frozen,
                                 assert_may_score, assert_rate_on_grid,
                                 assert_status_is_frozen, capture_seed,
                                 design_manifest, use_class_for)

#: Scanpy's own scrublet preprocessing, reproduced in `prepare_capture`.
#:
#: This tuple records THIS RUNNER'S EXPECTED ORDER and nothing more. It cannot
#: detect that scanpy changed its own implementation -- nothing here reads
#: scanpy's source. What makes an upstream change visible is the version guard
#: below: the recorded contract version plus the resolved scanpy version pin the
#: pair that was actually executed, so a later reader can tell which
#: reproduction ran against which library.
SCANPY_PREPROCESSING_ORDER = (
    'filter_genes(min_cells=3)',
    'filter_cells(min_genes=3)',
    'save raw layer',
    'normalize_total',
    'log1p -> temporary layer',
    'highly_variable_genes(layer=log1p)',
    'subset to highly_variable',
    'scrublet_simulate_doublets(layer=raw)',
    'log1p both (when log_transform)',
    'normalize_total(both, target_sum=1e6)',
)

#: Bumped whenever SCANPY_PREPROCESSING_ORDER or the steps in `prepare_capture`
#: change. Written on EVERY row, including a capture-preparation failure row, so
#: no result exists without saying which reproduction produced it.
PREPROCESSING_CONTRACT_VERSION = 'gate2-preprocessing/1'

#: The scanpy versions this reproduction has been checked against. A real run on
#: any other version fails before the matrix is touched: the preprocessing here
#: mirrors scanpy's internals, and an unreviewed upgrade can silently change what
#: `adata_sim` was built from.
APPROVED_SCANPY_VERSIONS = ('1.11.3',)

#: THE complete set of parameters this path must have a ruling for. Frozen, so
#: "which parameters are open" is computed as `REQUIRED - ruled`, never supplied
#: by a caller.
#:
#: A previous draft let the caller pass `unadjudicated={}` and declare the
#: problem away; the tests did exactly that. That is the same defect as the
#: guard never being called, wearing a different hat.
REQUIRED_PARAMETER_NAMES = frozenset({
    'sim_doublet_ratio', 'synthetic_doublet_umi_subsampling',
    'log_transform',
    'n_prin_comps', 'n_neighbors', 'stdev_doublet_rate', 'knn_dist_metric',
    'normalize_variance', 'mean_center', 'use_approx_neighbors',
    'get_doublet_neighbor_parents',
})

#: Which interface each parameter belongs to. A simulation parameter passed to
#: the scoring call is a TypeError at best and a silently ignored kwarg at worst.
PARAMETER_GROUPS = {
    'simulation': ('sim_doublet_ratio', 'synthetic_doublet_umi_subsampling'),
    'preprocessing': ('log_transform',),
    'scoring': ('n_prin_comps', 'n_neighbors', 'stdev_doublet_rate', 'knn_dist_metric',
                'normalize_variance', 'mean_center', 'use_approx_neighbors',
                'get_doublet_neighbor_parents'),
}

# The groups must partition the required set exactly: no parameter unassigned,
# none in two groups. Checked at import so a later edit cannot drift.
_grouped = [n for names in PARAMETER_GROUPS.values() for n in names]
assert set(_grouped) == set(REQUIRED_PARAMETER_NAMES), (
    'PARAMETER_GROUPS does not cover REQUIRED_PARAMETER_NAMES exactly: '
    f'{set(_grouped) ^ set(REQUIRED_PARAMETER_NAMES)}')
assert len(_grouped) == len(set(_grouped)), 'a parameter appears in two groups'
del _grouped

#: Deployment slot, deliberately inactive. The eleven approved methodological
#: values live in gate2_parameter_ruling; an empty deployment slot no longer
#: means that the parameters are undecided. Promotion requires a reviewed
#: capture/config and a recorded pilot permit, not another parameter choice.
PRODUCTION_PARAMETER_RULING = {}


class Gate2RunnerError(RuntimeError):
    """A runner-level defect. Distinct from Gate2DesignError so an audit can
    tell a violated ruling from a broken implementation."""


class CapturePreparationError(Gate2RunnerError):
    """A capture is not preparable: a DATA feasibility problem, raised
    deliberately by an explicit check.

    Only this exception becomes a capture-level UNKNOWN. A previous draft
    converted every non-fatal Exception, which meant a TypeError, AttributeError
    or NameError -- ordinary code defects -- were recorded as data-quality
    findings.
    """


def validate_parameter_ruling(ruling):
    """Return the ruling, or raise naming exactly what is wrong.

    Open parameters are REQUIRED minus ruled. A caller cannot assert that a
    parameter is settled; it can only supply a ruling that covers it.
    """
    if ruling is None:
        if not PRODUCTION_PARAMETER_RULING:
            raise Gate2RunnerError(
                'production configuration is not activated: all eleven parameters '
                'are methodologically approved, but capture C1-C3, config audit '
                'and scoped pilot authorization must be recorded before activation')
        ruling = PRODUCTION_PARAMETER_RULING
    if not isinstance(ruling, dict):
        raise Gate2RunnerError(f'parameter ruling must be a mapping, got {type(ruling).__name__}')
    unknown = sorted(set(ruling) - REQUIRED_PARAMETER_NAMES)
    if unknown:
        raise Gate2RunnerError(
            f'parameter ruling contains names this path does not use: {unknown}')
    missing = sorted(REQUIRED_PARAMETER_NAMES - set(ruling))
    if missing:
        raise Gate2RunnerError(
            'these Scrublet parameters have no methodological ruling and this path '
            f'will not run until they do: {missing}. Recording scanpy\'s default is '
            'provenance, not approval.')
    return dict(ruling)


def parameters_for(group, ruling):
    """The kwargs for one interface, from the validated ruling only."""
    if group not in PARAMETER_GROUPS:
        raise Gate2RunnerError(f'unknown parameter group {group!r}')
    return {k: ruling[k] for k in PARAMETER_GROUPS[group]}


def scrublet_call_parameters(rate, seed, ruling, *, n_obs=None):
    """THE scoring kwargs. The row's provenance is this same object."""
    params = parameters_for('scoring', ruling)
    if params['n_neighbors'] is None:
        if n_obs is None or n_obs < 1:
            raise Gate2RunnerError('n_neighbors rule requires the post-filter observed cell count')
        resolved = round(0.5 * np.sqrt(n_obs))
        if resolved < 1:
            raise CapturePreparationError('resolved n_neighbors is below one')
        # Scanpy echoes a resolved integer. Passing/recording None made the old
        # equality check reject a correct library result as a parameter mismatch.
        params['n_neighbors'] = int(resolved)
    params.update(expected_doublet_rate=float(rate), random_state=int(seed))
    return params


#: Secondary-seed design does not exist yet, so seed stability is NOT computed.
#: Naming the gap is the record; inventing a second seed would be choosing one.
SEED_STABILITY_STATUS = 'not_computed_pending_secondary_seed_design'


#: The exceptions a scanpy/numpy/scipy call may raise ABOUT THE DATA.
#: Deliberately a short allow-list caught at each minimal call boundary, not
#: "every Exception except a few".
#:
#: `TypeError`, `AttributeError`, `NameError`, `KeyError`, `AssertionError` are
#: NOT here: an unexpected keyword, a missing attribute, an absent dict key or a
#: typo is a defect in this file, and filing it as a doublet-detection finding
#: would bury it. A previous draft caught every `Exception` bar four, and a test
#: asserted that an arbitrary `TypeError` counted as data degeneracy.
LIBRARY_DATA_EXCEPTIONS = (ValueError, IndexError, ZeroDivisionError,
                           FloatingPointError, np.linalg.LinAlgError)


def _step(stage, fn, *args, **kwargs):
    """Run one library call, converting only a data failure, and recording WHICH
    stage failed. Chained so the original traceback survives."""
    try:
        return fn(*args, **kwargs)
    except LIBRARY_DATA_EXCEPTIONS as exc:
        raise CapturePreparationError(f'{stage}: {type(exc).__name__}: {exc}') from exc


def dependency_identity(sc):
    """Who actually provided the preprocessing, resolved from the RUNNING module.

    A real scanpy is identified by its own `__version__`. A test double has no
    version and must declare `gate2_test_identity`; it is then recorded as a test
    double with `scanpy_version=None`, so no row produced under a fake can ever
    carry a production version string. A module with neither is refused rather
    than recorded as unknown.
    """
    version = getattr(sc, '__version__', None)
    if version is not None:
        return {'scanpy_version': str(version),
                'scanpy_identity_source': 'module.__version__',
                'preprocessing_contract_version': PREPROCESSING_CONTRACT_VERSION}
    declared = getattr(sc, 'gate2_test_identity', None)
    if declared is None:
        raise Gate2RunnerError(
            'the injected preprocessing module has no __version__ and does not '
            'declare gate2_test_identity; provenance cannot describe what ran')
    return {'scanpy_version': None,
            'scanpy_identity_source': 'test_double',
            'test_double_identity': str(declared),
            'preprocessing_contract_version': PREPROCESSING_CONTRACT_VERSION}


def assert_dependency_approved(identity):
    """Refuse an unreviewed scanpy before any matrix work.

    A test double is not version-checked -- it never claims a version -- but it
    is recorded as a double, which is the part that matters downstream.
    """
    if identity['scanpy_identity_source'] == 'test_double':
        return identity
    if identity['scanpy_version'] not in APPROVED_SCANPY_VERSIONS:
        raise Gate2RunnerError(
            f"scanpy {identity['scanpy_version']} is not in the approved set "
            f'{APPROVED_SCANPY_VERSIONS}. This module reproduces scanpy\'s own '
            'scrublet preprocessing step for step; an unreviewed version may '
            'reproduce something else.')
    return identity


def _checksum_matrix(x):
    """Stable checksum of a count matrix, dense or sparse."""
    h = hashlib.sha256()
    if hasattr(x, 'tocsr'):
        m = x.tocsr()
        for part in (m.indptr, m.indices, m.data):
            h.update(np.ascontiguousarray(part).tobytes())
        h.update(str(m.shape).encode())
    else:
        a = np.ascontiguousarray(x)
        h.update(a.tobytes())
        h.update(str(a.shape).encode())
    return h.hexdigest()


#: The only obs columns allowed past the runner boundary. Everything else the
#: caller loaded is dropped before scanpy sees it.
#:
#: This does NOT mean outcome was never loaded -- the caller hands this function
#: an already-constructed AnnData and the runner cannot control what happened
#: before that. What it means is that nothing outside this list is copied into
#: the object scanpy scores, or into anything scanpy returns.
APPROVED_OBS_COLUMNS = ()


def minimise_input(adata, cohort, channel_id):
    """Strip the capture down to the matrix and its cell identity.

    Returns `(minimal, retained, dropped)`. Column NAMES are recorded; values
    are never inspected. `adata.copy()` alone would carry the caller's entire
    obs into preprocessing and into scanpy's output.
    """
    import anndata as _ad
    obs_cols = list(getattr(getattr(adata, 'obs', None), 'columns', []))
    retained = [c for c in obs_cols if c in APPROVED_OBS_COLUMNS]
    dropped = [c for c in obs_cols if c not in APPROVED_OBS_COLUMNS]
    minimal = _ad.AnnData(X=adata.X.copy())
    minimal.obs_names = list(adata.obs_names)
    minimal.var_names = list(adata.var_names)
    for c in retained:
        minimal.obs[c] = adata.obs[c].values
    return minimal, retained, dropped


def prepare_capture(adata, cohort, channel_id, ruling, *, sc=None):
    """Preprocess one capture and simulate its doublets ONCE.

    Returns `(ad_obs, ad_sim, meta)`. The same `ad_sim` is
    then passed to every grid point, which is what makes the grid a calibration
    sweep rather than eight independent experiments.

    `ruling` is the SAME validated object the scoring call is built from, so the
    simulation and preprocessing choices recorded on each row cannot differ from
    the ones the simulator actually received. A previous draft derived the row's
    `simulation_parameters` from an injected object while the simulator kept
    using module constants -- the two were free to disagree.

    Raises `CapturePreparationError` for a data feasibility problem. Anything
    else propagates: a `TypeError` here is a bug in this file, not a finding
    about the capture.

    `sc` is injected so the draft is testable without scanpy present; the runner
    passes the real module.
    """
    if sc is None:                       # pragma: no cover - real runs inject it
        import scanpy as sc

    ruling = validate_parameter_ruling(ruling)
    sim_kwargs = parameters_for('simulation', ruling)
    prep = parameters_for('preprocessing', ruling)
    identity = assert_dependency_approved(dependency_identity(sc))
    assert_may_score(cohort)             # before any matrix work
    seed = capture_seed(cohort, channel_id)

    ad, retained, dropped = minimise_input(adata, cohort, channel_id)
    obs_names_checksum = hashlib.sha256(
        '\u0000'.join(map(str, ad.obs_names)).encode()).hexdigest()
    _step('filter_genes', sc.pp.filter_genes, ad, min_cells=3)
    _step('filter_cells', sc.pp.filter_cells, ad, min_genes=3)
    if ad.n_obs == 0 or ad.n_vars == 0:
        raise CapturePreparationError(
            f'{cohort}/{channel_id}: filtering emptied the capture '
            f'({ad.n_obs} cells x {ad.n_vars} genes)')

    ad.layers['raw'] = ad.X.copy()
    _step('normalize_total', sc.pp.normalize_total, ad)
    ad.layers['log1p'] = ad.X.copy()
    _step('log1p', sc.pp.log1p, ad, layer='log1p')
    _step('highly_variable_genes', sc.pp.highly_variable_genes, ad, layer='log1p')
    del ad.layers['log1p']
    ad = ad[:, ad.var['highly_variable']].copy()
    if ad.n_vars == 0:
        raise CapturePreparationError(
            f'{cohort}/{channel_id}: no highly variable genes')

    ad_sim = _step('scrublet_simulate_doublets', sc.pp.scrublet_simulate_doublets,
                   ad, layer='raw', random_seed=seed, **sim_kwargs)
    del ad.layers['raw']
    if prep['log_transform']:
        _step('log1p_observed', sc.pp.log1p, ad)
        _step('log1p_simulated', sc.pp.log1p, ad_sim)
    _step('normalize_observed', sc.pp.normalize_total, ad, target_sum=1e6)
    _step('normalize_simulated', sc.pp.normalize_total, ad_sim, target_sum=1e6)

    # Recomputed AFTER filtering: this is the order every grid point scores, and
    # the positional alignment the Jaccard depends on.
    post_filter_checksum = hashlib.sha256(
        '\u0000'.join(map(str, ad.obs_names)).encode()).hexdigest()
    return ad, ad_sim, {
        'simulated_checksum': _checksum_matrix(ad_sim.X),
        'seed': seed,
        # Recorded from the SAME dict that was passed to the simulator above.
        'simulation_parameters': dict(sim_kwargs),
        'preprocessing_parameters': dict(prep),
        # Resolved from the module that actually did the work, not from a comment.
        'dependency_identity': identity,
        'obs_names_checksum_pre_filter': obs_names_checksum,
        'obs_names_checksum': post_filter_checksum,
        'obs_columns_retained': retained,
        'obs_columns_dropped': dropped,
    }


#: The fields the pinned scanpy actually emits, measured on 1.11.3 rather than
#: read from documentation:
#:
#:     out.obs['predicted_doublet']                bool
#:     out.obs['doublet_score']                    float64
#:     out.uns['scrublet']['threshold']            float
#:     out.uns['scrublet']['doublet_scores_sim']   ndarray, n_obs * sim_ratio
#:     out.uns['scrublet']['parameters']           dict of what it actually used
#:
#: A first draft of this module read `uns['scrublet']['_calls']` and
#: `['_scores']`. Those keys do not exist in scanpy; they existed only in this
#: repository's test double, which had invented them. Every unit test passed
#: while a real successful run would have been classified
#: `calls_not_strict_boolean`.
SCANPY_OBS_CALLS = 'predicted_doublet'
SCANPY_OBS_SCORES = 'doublet_score'
SCANPY_UNS_KEY = 'scrublet'
SCANPY_UNS_THRESHOLD = 'threshold'
SCANPY_UNS_SIM_SCORES = 'doublet_scores_sim'
SCANPY_UNS_PARAMETERS = 'parameters'


def _read_score_vectors(out, n_cells):
    """Observed and simulated scores, read BEFORE and INDEPENDENTLY of the
    threshold.

    A first draft checked the threshold first and returned `observed=None` on
    failure, so a real threshold failure silently discarded the observed
    histogram, percentiles and AUC -- the very evidence a failed capture still
    has. The test that was supposed to cover this fed `threshold=None` straight
    to the diagnostics helper, bypassing this function, and so proved only that
    the helper tolerates a missing threshold.
    """
    uns = (getattr(out, 'uns', {}) or {}).get(SCANPY_UNS_KEY) or {}
    obs = getattr(out, 'obs', None)
    cols = list(getattr(obs, 'columns', [])) if obs is not None else []

    sim = uns.get(SCANPY_UNS_SIM_SCORES)
    sim_arr = np.asarray(sim, dtype=float) if sim is not None else None

    obs_arr, obs_problem = None, None
    if SCANPY_OBS_SCORES not in cols:
        obs_problem = 'scores_wrong_shape'
    else:
        cand = np.asarray(obs[SCANPY_OBS_SCORES], dtype=float)
        if cand.ndim != 1 or cand.shape[0] != n_cells:
            obs_problem = 'scores_wrong_shape'
        elif not np.all(np.isfinite(cand)):
            obs_problem = 'scores_not_finite'
            obs_arr = cand          # kept: the diagnostics helper drops non-finite itself
        else:
            obs_arr = cand
    return obs_arr, sim_arr, obs_problem


def _mechanical_status(out, n_cells):
    """Threshold and calls. Returns (status, calls, observed, simulated, threshold).

    The score vectors are read first and returned whatever the threshold does,
    so a threshold failure keeps its distribution evidence.
    """
    uns = (getattr(out, 'uns', {}) or {}).get(SCANPY_UNS_KEY) or {}
    obs = getattr(out, 'obs', None)
    cols = list(getattr(obs, 'columns', [])) if obs is not None else []
    obs_arr, sim_arr, obs_problem = _read_score_vectors(out, n_cells)

    def done(status, calls=None, thr=None):
        return assert_status_is_frozen(status), calls, obs_arr, sim_arr, thr

    if SCANPY_UNS_THRESHOLD not in uns:
        return done('threshold_absent')
    try:
        thr = float(uns[SCANPY_UNS_THRESHOLD])
    except (TypeError, ValueError):
        return done('threshold_non_finite')
    if not np.isfinite(thr):
        return done('threshold_non_finite')

    if SCANPY_OBS_CALLS not in cols:
        return done('calls_not_strict_boolean', thr=thr)
    raw = np.asarray(obs[SCANPY_OBS_CALLS])
    if raw.ndim != 1 or raw.shape[0] != n_cells:
        return done('calls_not_strict_boolean', thr=thr)
    if raw.dtype == object and any(x is None or x != x for x in raw):
        return done('calls_contain_na', thr=thr)
    if raw.dtype != bool and not (raw.dtype.kind in 'iu'
                                  and bool(np.all((raw == 0) | (raw == 1)))):
        return done('calls_not_strict_boolean', thr=thr)

    if obs_problem is not None:
        return done(obs_problem, thr=thr)
    return done('examined', calls=raw.astype(bool), thr=thr)




def score_capture_across_grid(cohort, channel_id, adata, *, scrublet_fn, sc=None,
                              grid=RESIDUAL_RATE_GRID, ruling=None):
    """One capture, all eight grid points, one simulated dataset.

    The production deployment slot remains empty pending capture C1-C3, config
    audit and a pilot permit, although all eleven methodological choices are
    recorded. A technical fixture may inject a complete explicit ruling; doing
    so never supplies cohort execution authorization.
    """
    if sc is None:                       # pragma: no cover - real runs inject it
        import scanpy as sc

    assert_grid_is_frozen(grid)
    # All four guards run BEFORE the matrix is copied or preprocessed. The
    # use-class check concerns work THIS function would do; it says nothing
    # about the caller's load, which has already happened. The dependency is
    # resolved here rather than inside the try, so a preparation failure row
    # still carries the identity of the module that failed.
    use_class = assert_may_score(cohort)
    ruling = validate_parameter_ruling(ruling)
    identity = assert_dependency_approved(dependency_identity(sc))

    try:
        ad_obs, ad_sim, meta = prepare_capture(adata, cohort, channel_id, ruling, sc=sc)
    except CapturePreparationError as exc:
        # ONLY the dedicated exception becomes a data finding. A previous draft
        # caught every non-fatal Exception, so a TypeError, AttributeError or
        # NameError -- ordinary code defects -- were filed as capture quality
        # problems. Unknown exceptions now propagate and stay visible as bugs.
        # Same shape as the success path: a caller must be able to unpack one
        # way. Returning a bare list here made the failure path a different API.
        return ([_preparation_failure_row(cohort, channel_id, use_class, rate,
                                          exc, identity)
                 for rate in grid], {}, {})

    # Two sidecars, deliberately separate. `calls` holds boolean call vectors and
    # is filled ONLY for examined points, because the Jaccard is defined only
    # there. `scores` holds the per-cell OBSERVED doublet scores and is filled
    # whenever scanpy produced them -- including when the threshold failed, which
    # is exactly when the distribution is the remaining evidence. Inferring
    # scores from calls is impossible and was never attempted.
    sidecars, score_sidecars, rows = {}, {}, []
    for rate in grid:
        assert_rate_on_grid(rate)
        call_kwargs = scrublet_call_parameters(rate, meta['seed'], ruling, n_obs=ad_obs.n_obs)
        row = {
            'design_version': DESIGN_VERSION,
            'decision_reference': DECISION_REFERENCE,
            'grid_semantics': design_manifest()['grid_semantics'],
            'diagnostics_schema_version': DIAGNOSTICS_SCHEMA_VERSION,
            'preprocessing_contract_version':
                meta['dependency_identity']['preprocessing_contract_version'],
            'scanpy_version': meta['dependency_identity']['scanpy_version'],
            'scanpy_identity_source':
                meta['dependency_identity']['scanpy_identity_source'],
            'cohort': cohort,
            'channel_id': channel_id,
            'use_class': use_class,
            'never_pooled': cohort.strip().lower() in NEVER_POOLED,
            'expected_doublet_rate': float(rate),
            'capture_seed': meta['seed'],
            'simulated_checksum': meta['simulated_checksum'],
            'obs_names_checksum': meta['obs_names_checksum'],
            'obs_columns_retained': list(meta['obs_columns_retained']),
            'obs_columns_dropped': list(meta['obs_columns_dropped']),
            'n_cells': int(ad_obs.n_obs),
            'seed_stability': SEED_STABILITY_STATUS,
            # The provenance IS the call. Same object, not a parallel record.
            'scrublet_parameters': dict(call_kwargs),
            'simulation_parameters': dict(meta['simulation_parameters']),
            'preprocessing_parameters': dict(meta['preprocessing_parameters']),
        }

        try:
            out = scrublet_fn(ad_obs, adata_sim=ad_sim, copy=True, **call_kwargs)
        except LIBRARY_DATA_EXCEPTIONS as exc:
            row.update(mechanical_status=assert_status_is_frozen('scrublet_raised'),
                       failure_detail=f'{type(exc).__name__}: {exc}',
                       n_called=UNKNOWN_FILL_VALUE, threshold=None,
                       scanpy_reported_parameters=None, identity_verified=None,
                       parameters_echo_agree=None)
            row.update(_reliability_block(None, None, None))
            rows.append(_finalize_row(row))
            continue

        # The order every grid point scores must be the order the Jaccard
        # assumes. Verified per call rather than asserted in a comment.
        row['identity_verified'] = (
            hashlib.sha256('\u0000'.join(map(str, out.obs_names)).encode()).hexdigest()
            == meta['obs_names_checksum'])

        uns = (getattr(out, 'uns', {}) or {}).get(SCANPY_UNS_KEY) or {}
        row['scanpy_reported_parameters'] = dict(uns.get(SCANPY_UNS_PARAMETERS) or {})
        row['parameters_echo_agree'] = _echo_agrees(call_kwargs,
                                                    row['scanpy_reported_parameters'])

        status, calls, obs_scores, sim_scores, thr = _mechanical_status(out, ad_obs.n_obs)
        # A threshold failure keeps identity-verified scores. An out-of-order
        # vector cannot be labelled with the input cell order in a sidecar.
        if obs_scores is not None and row['identity_verified']:
            score_sidecars[float(rate)] = obs_scores
        # Computed whatever the threshold did: a threshold failure keeps its
        # distribution evidence.
        row.update(_reliability_block(obs_scores, sim_scores, thr))

        # Two checks that a previous draft RECORDED and then ignored, letting the
        # grid point be `examined` and enter the Jaccard anyway. A check whose
        # result changes nothing is not a check.
        detail = None
        # Applied only to a point that would OTHERWISE be examined: overwriting
        # `threshold_absent` with `identity_mismatch` would hide the more
        # specific scanpy-level diagnosis. Both exclude the point either way,
        # and `identity_verified` is on the row regardless.
        if status != 'examined':
            pass
        elif row['identity_verified'] is False:
            # Positional alignment is the only thing making the Jaccard meaningful.
            status, detail = 'identity_mismatch', 'obs_names differ from preparation'
        elif not row['parameters_echo_agree']['agree']:
            # scanpy says it used a value other than the one requested, so the
            # row's provenance does not describe the run that produced it.
            disagreed = sorted(k for k in row['parameters_echo_agree']['compared']
                               if call_kwargs[k] != row['scanpy_reported_parameters'][k])
            status, detail = 'parameter_echo_mismatch', f'scanpy echoed {disagreed}'

        if status != 'examined':
            row.update(mechanical_status=assert_status_is_frozen(status),
                       failure_detail=detail,
                       n_called=UNKNOWN_FILL_VALUE, threshold=None)
        else:
            row.update(mechanical_status='examined', failure_detail=None,
                       n_called=int(calls.sum()), threshold=float(thr))
            sidecars[float(rate)] = calls    # ONLY examined points reach the Jaccard
        rows.append(_finalize_row(row))
    return rows, sidecars, score_sidecars


def _to_native(obj):
    """Coerce a row to natively JSON-serialisable types, RAISING on anything
    unrecognised.

    A previous test used `json.dumps(rows, default=str)`, which stringifies any
    object it cannot serialise -- so it would have passed on a row full of
    ndarrays and `np.float64`, and proved nothing about the checkpoint write that
    will actually consume these rows. This converts what it knows and refuses the
    rest, so an unexpected type is a failure here rather than a surprise later.
    """
    # NumPy scalars are checked FIRST, and this order is load-bearing:
    # `np.float64` is a SUBCLASS of Python `float`, so an `isinstance(obj, float)`
    # branch placed above this one returns it unconverted -- which is exactly
    # what happened here, and was only caught once a NumPy scalar was put into
    # the test ruling.
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        v = float(obj)
        return v if np.isfinite(v) else None      # JSON has no NaN/Infinity
    if obj is None or isinstance(obj, (str, bool, int)):
        return obj
    if isinstance(obj, float):
        # A plain float carrying NaN/Infinity would otherwise blow up at
        # `allow_nan=False` in whatever writes the checkpoint, long after the run.
        return obj if np.isfinite(obj) else None
    if isinstance(obj, np.ndarray):
        return [_to_native(x) for x in obj.tolist()]
    if isinstance(obj, dict):
        return {str(k): _to_native(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_native(x) for x in obj]
    raise Gate2RunnerError(
        f'row field of type {type(obj).__name__} is not JSON-serialisable and has '
        'no declared conversion')


#: Every row carries every one of these, on every exit path. Missing is an
#: error, not an implicit None: a failure row that quietly lacks a field makes
#: the failure and success schemas different, and only the reader downstream
#: finds out.
ROW_SCHEMA = (
    'design_version', 'decision_reference', 'grid_semantics',
    'diagnostics_schema_version', 'preprocessing_contract_version',
    'scanpy_version', 'scanpy_identity_source',
    'cohort', 'channel_id', 'use_class', 'never_pooled',
    'expected_doublet_rate', 'capture_seed', 'simulated_checksum',
    'obs_names_checksum', 'obs_columns_retained', 'obs_columns_dropped',
    'n_cells', 'seed_stability', 'scrublet_parameters', 'simulation_parameters',
    'preprocessing_parameters', 'scanpy_reported_parameters',
    'parameters_echo_agree', 'identity_verified', 'mechanical_status',
    'failure_detail', 'n_called', 'threshold', 'distribution_summary',
) + RELIABILITY_DIAGNOSTICS_REQUIRED


def _finalize_row(row):
    """The ONE exit for every row: schema, native types, strict JSON.

    A previous draft applied `_to_native` on the success path and appended the
    raw dict on the `scrublet_raised` path, so the failure rows -- the ones most
    likely to carry an odd value -- were the ones that skipped the conversion.
    Every path now returns through here, and the strict dump happens in
    PRODUCTION, not only in a test, so a contract violation surfaces at the row
    that caused it rather than at the checkpoint write.
    """
    import json as _json
    missing = [k for k in ROW_SCHEMA if k not in row]
    if missing:
        raise Gate2RunnerError(f'row is missing required fields: {missing}')
    extra = [k for k in row if k not in ROW_SCHEMA]
    if extra:
        raise Gate2RunnerError(f'row carries undeclared fields: {sorted(extra)}')
    native = _to_native(row)
    try:
        _json.dumps(native, allow_nan=False)
    except (TypeError, ValueError) as exc:
        raise Gate2RunnerError(f'row does not satisfy the JSON contract: {exc}') from exc
    return native


def _echo_agrees(requested, echoed):
    """Compare what we asked for against what scanpy says it used.

    scanpy echoes only some parameters, so the result names the boundary rather
    than pretending to a full check.
    """
    common = sorted(set(requested) & set(echoed))
    return {
        'compared': common,
        'not_echoed_by_scanpy': sorted(set(requested) - set(echoed)),
        'agree': all(requested[k] == echoed[k] for k in common),
    }


def _reliability_block(observed, simulated, threshold):
    """The reliability diagnostics, computed by the repository's existing
    module rather than by a second implementation.

    Each conceptual slot in RELIABILITY_DIAGNOSTICS_REQUIRED is mapped to the
    concrete fields the data contract already defines. A first draft set all
    five to None permanently while the report called them 'recorded'.
    """
    summary = score_distribution_summary(observed, simulated, threshold)
    return {
        'simulated_score_histogram': {
            'counts': summary.get('simulated_histogram_counts'),
            'bin_edges': summary.get('shared_histogram_bin_edges'),
            'percentiles_p0_to_p100': summary.get('simulated_score_percentiles_p0_to_p100'),
        },
        'observed_score_histogram': {
            'counts': summary.get('observed_histogram_counts'),
            'bin_edges': summary.get('shared_histogram_bin_edges'),
            'percentiles_p0_to_p100': summary.get('observed_score_percentiles_p0_to_p100'),
        },
        'threshold_placement': {k: summary.get(k) for k in summary
                                if 'threshold' in k or 'above' in k or 'tail' in k},
        'tail_enrichment': {'observed_simulated_auc': summary.get('observed_simulated_auc')},
        'seed_stability': SEED_STABILITY_STATUS,
        'distribution_summary': summary,
    }


def _preparation_failure_row(cohort, channel_id, use_class, rate, exc, identity):
    """Same schema as every other row. The nulls are written out explicitly
    rather than left to `_finalize_row` to fill: a field silently defaulted is a
    field nobody decided about."""
    row = {
        'design_version': DESIGN_VERSION,
        'decision_reference': DECISION_REFERENCE,
        'grid_semantics': design_manifest()['grid_semantics'],
        'diagnostics_schema_version': DIAGNOSTICS_SCHEMA_VERSION,
        'preprocessing_contract_version': PREPROCESSING_CONTRACT_VERSION,
        'scanpy_version': identity['scanpy_version'],
        'scanpy_identity_source': identity['scanpy_identity_source'],
        'cohort': cohort,
        'channel_id': channel_id,
        'use_class': use_class,
        'never_pooled': cohort.strip().lower() in NEVER_POOLED,
        'expected_doublet_rate': float(rate),
        'capture_seed': None,
        'simulated_checksum': None,
        'obs_names_checksum': None,
        'obs_columns_retained': None,
        'obs_columns_dropped': None,
        'n_cells': None,
        'seed_stability': SEED_STABILITY_STATUS,
        'scrublet_parameters': None,
        'simulation_parameters': None,
        'preprocessing_parameters': None,
        'scanpy_reported_parameters': None,
        'parameters_echo_agree': None,
        'identity_verified': None,
        'mechanical_status': assert_status_is_frozen(CAPTURE_PREPARATION_FAILED),
        'failure_detail': f'{type(exc).__name__}: {exc}',
        'n_called': UNKNOWN_FILL_VALUE,
        'threshold': None,
        'distribution_summary': None,
    }
    for diag in RELIABILITY_DIAGNOSTICS_REQUIRED:
        row.setdefault(diag, None)
    return _finalize_row(row)


def grid_stability(rows, sidecars=None):
    """How the CALL SET moves across the grid. NO BEST RATE.

    Counts are not sets: two rates can call the same NUMBER of cells and share
    almost none of them. A first draft reported only n_called min/max/span and
    called it stability, which it is not.

    Agreement is Jaccard over the boolean call vectors, aligned by position
    within one capture -- every grid point scores the same `ad_obs` in the same
    order, so position IS cell identity here; a future cross-capture comparison
    would have to align on cell ids explicitly.

    The full pairwise matrix is returned along with its min and max. There is
    deliberately no argmax, no ranking and no field a caller could read as a
    recommendation.
    """
    examined = [r for r in rows if r.get('mechanical_status') == 'examined']
    called = [r['n_called'] for r in examined]
    sidecars = sidecars or {}
    vecs = [(r['expected_doublet_rate'], sidecars.get(float(r['expected_doublet_rate'])))
            for r in examined
            if sidecars.get(float(r['expected_doublet_rate'])) is not None]

    pairwise, jac = {}, []
    for i in range(len(vecs)):
        for j in range(i + 1, len(vecs)):
            a, b = vecs[i][1], vecs[j][1]
            union = int(np.logical_or(a, b).sum())
            inter = int(np.logical_and(a, b).sum())
            # Both call sets empty: identical by definition, not undefined.
            v = 1.0 if union == 0 else inter / union
            pairwise[f'{vecs[i][0]}|{vecs[j][0]}'] = v
            jac.append(v)

    out = {
        'n_grid_points': len(rows),
        'n_examined': len(examined),
        'n_unknown': len(rows) - len(examined),
        'unknown_reasons': sorted({r['mechanical_status'] for r in rows
                                   if r.get('mechanical_status') != 'examined'}),
        'n_called_min': min(called) if called else None,
        'n_called_max': max(called) if called else None,
        'n_called_span': (max(called) - min(called)) if called else None,
        'call_set_jaccard_pairwise': pairwise,
        'call_set_jaccard_min': min(jac) if jac else None,
        'call_set_jaccard_max': max(jac) if jac else None,
        'call_set_alignment': 'positional within one capture; same ad_obs, same order',
        'identity_verified_all_points': all(r.get('identity_verified') is True
                                            for r in examined) if examined else None,
        'seed_stability': SEED_STABILITY_STATUS,
        'interpretation_note': (
            'A low Jaccard minimum means the candidate set moves with the rate, i.e. the '
            'method is unstable on this capture. It does NOT indicate which rate to use; '
            'no rate is selected from these numbers. A stable n_called with a low Jaccard '
            'is instability, not agreement.'),
    }
    sims = {r['simulated_checksum'] for r in rows}
    out['common_random_numbers_held'] = (len(sims) == 1 and None not in sims)
    if not out['common_random_numbers_held']:
        out['common_random_numbers_detail'] = (
            f'{len(sims)} distinct simulated datasets across the grid; the calibration '
            'effect cannot be separated from Monte Carlo variation')
    return out


def assert_not_pooled_with_others(rows):
    """Refuse to aggregate a NEVER_POOLED cohort with any other.

    Raises rather than dropping the offending rows: silently excluding them
    would produce a plausible aggregate from a forbidden input.
    """
    cohorts = {r['cohort'].strip().lower() for r in rows}
    offending = cohorts & NEVER_POOLED
    if offending and len(cohorts) > 1:
        raise Gate2RunnerError(
            f'cohort(s) {sorted(offending)} are never pooled, but this aggregation also '
            f'contains {sorted(cohorts - NEVER_POOLED)}')
    return sorted(cohorts)
