"""NON-RUNNABLE scoring config for GSM5102900_HC1, bound to measured evidence.

This module builds a config and refuses to make it runnable. Scoring requires a
separate scoped pilot-permit artifact that matches the config and input hashes
recorded here; no such permit exists, so `assert_permitted()` raises.

Everything below is bound to something that was measured, not assumed:
capture, the three deposit file hashes, the deploy commit, the preflight
manifest hash, the dependency versions the preflight actually ran under, the
eleven approved parameters imported from the ruling, and the `n_neighbors` /
`k_adj` the preflight resolved.

`graph_preflight_verified` is **false**, and the reason is recorded: the kNN
graph is built inside `calculate_doublet_scores`, so the score-free preflight
could not check it. A scoring run would be the first time that graph exists.
"""
import hashlib
import json
import os
import types

from gate2_parameter_ruling import (METHODOLOGY_PARAMETER_RULING,
                                    PARAMETER_RULING_REFERENCE,
                                    PARAMETER_RULING_VERSION, approved_for,
                                    assert_ruling_version)
from gate2_stress_design import RESIDUAL_RATE_GRID, capture_seed
from gate2_stress_runner import REQUIRED_PARAMETER_NAMES, Gate2RunnerError

CONFIG_VERSION = 'gate2-scoring-config-draft/GSM5102900_HC1/1'
EXPECT_RULING_VERSION = 'gate2-parameter-ruling/1'

CAPTURE = 'GSM5102900_HC1'
COHORT, CHANNEL = 'qiu', 'GSM5102900'
#: The commit the SCORE-FREE PREFLIGHT ran from. It is NOT a scoring deployment
#: commit -- the scoring entrypoint did not exist at 9aab7be. The permit binds the
#: scoring implementation separately, by file hash, because those files are newer.
PREFLIGHT_COMMIT = '9aab7befba0bd2a59ea1224aff95d261237611c6'
PREFLIGHT_JOB_ID = '29073092'
PREFLIGHT_MANIFEST_SHA256 = (
    'a78e5a01b4ff05f8a3901d77e348a3006a4dee2ed9e3fff8733fceb0ecc7a532')

DEPOSIT_DIR = '/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/GSE167363_qiu/extracted'
INPUT_SHA256 = types.MappingProxyType({
    'GSM5102900_HC1_barcodes.tsv.gz':
        'e6dbf9f587d9a6cf3ecff31cbf4d53b09c6faf9837c82d5c1c0a26227f377a02',
    'GSM5102900_HC1_features.tsv.gz':
        'c7408e82af270ca6308a8aa7ea659795a15ab76f2e3a85bd5db4cd3d53244199',
    'GSM5102900_HC1_matrix.mtx.gz':
        '4ea259bfc4de02cc9799df563779fce7c9bd4c889532d3ff1178e2caf545b9f8',
})

#: Measured by the preflight, not assumed. `n_neighbors` is the ruling's rule
#: resolved on THIS capture's post-filter N; `k_adj` is scanpy's own adjustment.
RESOLVED = types.MappingProxyType({
    'n_cells_effective': 7666,
    'n_usable_features_hvg': 1576,
    'n_simulated': 15332,
    'n_neighbors': 44,
    'k_adj': 132,
})

# ---- graph acceptance, MEASURED on the pinned scanpy, not assumed -----------
#
# scanpy 1.11.3 does NOT reach sklearn's `NearestNeighbors` on this path: it goes
# through `scanpy.Neighbors.compute_neighbors`, whose backend is a
# `KNeighborsTransformer` (and pynndescent is installed in this environment, so
# the backend is not even guaranteed to be sklearn). The graph the doublet score
# actually CONSUMES is the return of
# `_get_indices_distances_from_sparse_matrix(knn.distances, k_adj)` at
# scanpy/preprocessing/_scrublet/core.py:353, which is backend-independent. That
# is the capture point.
#
# Measured with a pinned-scanpy probe (300 obs / 600 sim, n_prin_comps=10), and
# re-derived by the `scanpy_integration` test so a scanpy change breaks the test
# rather than silently changing production:
#
#   constructions per prior : 1
#   n_neighbors argument    : k_adj
#   returned index columns  : k_adj - 1        (NOT k_adj, despite the helper's
#                                               docstring)
#   first column is self    : False            (likewise contrary to the docstring)
#
# Asserting self-inclusion or k_adj columns would have failed a CORRECT run.
GRAPH_CAPTURE_POINT = ('scanpy.preprocessing._scrublet.core.'
                       '_get_indices_distances_from_sparse_matrix')
GRAPH_CONSTRUCTIONS_PER_PRIOR = 1
GRAPH_COLUMNS_RULE = 'k_adj - 1'


def expected_graph_shape():
    """(rows, columns) the consumed graph must have for this capture."""
    rows = RESOLVED['n_cells_effective'] + RESOLVED['n_simulated']   # 22,998
    return (rows, RESOLVED['k_adj'] - 1)                             # (22998, 131)
DEPENDENCY_VERSIONS = types.MappingProxyType({
    'scanpy': '1.11.3', 'anndata': '0.11.4', 'numpy': '2.2.6',
    'scipy': '1.15.3', 'sklearn': '1.7.1', 'python': '3.10.18',
})
THRESHOLD_RULE = ('scanpy automatic histogram threshold; a failure is recorded as '
                  'unknown/NA and never as 0, False or singlet')
PRIOR_LABEL = 'assumed_algorithmic_stress_setting'

#: The permit that does not exist. A scoped pilot permit must name this config
#: hash and these input hashes; anything else is refused.
OUTPUT_ROOT = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/gate2_capture_scoring'

#: The files that ARE the scoring implementation. The permit binds their hashes,
#: so a permit cannot survive an edit to the entrypoint, the wrapper, the config
#: or any imported module.
IMPLEMENTATION_FILES = (
    'scripts/32_rebuild_upfront_qc/code/gate2_score_capture_GSM5102900_HC1.py',
    'scripts/32_rebuild_upfront_qc/code/gate2_scoring_config_GSM5102900_HC1.py',
    'scripts/32_rebuild_upfront_qc/code/gate2_parameter_ruling.py',
    'scripts/32_rebuild_upfront_qc/code/gate2_stress_design.py',
    'scripts/32_rebuild_upfront_qc/code/gate2_stress_runner.py',
    # TRANSITIVE: gate2_stress_runner imports and uses this. Omitting it let a
    # changed diagnostic module alter the output while the permit stayed valid.
    'scripts/32_rebuild_upfront_qc/code/gate2_diagnostics.py',
    'reports/manuscript/revision/gate2_capture_audit_2026-09-24/run_capture_scoring.sh',
)

PERMIT_PATH_ENV = 'GATE2_PILOT_PERMIT'

#: A legacy environment bypass existed here and was REMOVED. Production code that
#: honours an env var disabling its own authenticity check is not a guard: anyone
#: who can set the job environment can turn it off. Tests now inject a verifier
#: through the function signature instead, which the production caller never does.
#: The wrapper additionally refuses to run if this variable is set at all.
REMOVED_BYPASS_ENV = 'GATE2_PERMIT_TEST_BYPASS'
PERMIT_REQUIRED_FIELDS = ('permit_id', 'authorised_by', 'authorised_utc', 'capture',
                          'config_sha256', 'input_sha256', 'scope',
                          'implementation_sha256', 'implementation_commit',
                          'resource_ceiling', 'output_root')

#: Exact values, not merely "present". An arbitrary non-empty string previously
#: satisfied both of these.
REQUIRED_SCOPE = ('one job, one capture, eight approved priors, scoring only; '
                  'no outcome access, no cell deletion, no retry, no second job')
REQUIRED_RESOURCE_CEILING = types.MappingProxyType({
    'partition': 'short', 'time': '2:00:00', 'mem': '48G', 'cpus_per_task': 4,
    'max_jobs': 1,
})

#: Columns that must never appear in a scoring output for this capture.
FORBIDDEN_OUTPUT_FIELDS = ('Death', 'y', 'covid_severity', 'mild_vs_severe',
                           'severity_score', 'three_group', 'progression_group',
                           'outcome', 'mortality', 'survival')


def _config_body():
    """Everything the config hash covers. Excludes the output directory, which is
    DERIVED from that hash -- including it would be circular."""
    assert_ruling_version(EXPECT_RULING_VERSION)
    approved = approved_for(REQUIRED_PARAMETER_NAMES)
    return {
        'config_version': CONFIG_VERSION,
        'runnable': False,
        'capture': CAPTURE, 'cohort': COHORT, 'channel_id': CHANNEL,
        'one_complete_capture': True,
        'input': {'deposit_dir': DEPOSIT_DIR, 'sha256': dict(INPUT_SHA256),
                  'population': 'the complete 7,666-cell CellRanger-called deposit '
                                'matrix, NOT the 65-cell canonical-container subset'},
        'provenance': {'preflight_commit': PREFLIGHT_COMMIT,
                       'preflight_job_id': PREFLIGHT_JOB_ID,
                       'preflight_manifest_sha256': PREFLIGHT_MANIFEST_SHA256},
        'ruling': {'version': PARAMETER_RULING_VERSION,
                   'reference': PARAMETER_RULING_REFERENCE,
                   'imported_not_copied': True},
        'parameters': dict(approved, n_neighbors=RESOLVED['n_neighbors']),
        'resolved': dict(RESOLVED),
        'prior_grid': list(RESIDUAL_RATE_GRID), 'prior_label': PRIOR_LABEL,
        'seed': {'value': capture_seed(COHORT, CHANNEL),
                 'derivation': 'capture_seed(cohort, channel_id)'},
        'threshold_rule': THRESHOLD_RULE,
        'dependency_versions': dict(DEPENDENCY_VERSIONS),
        'output': {'root': OUTPUT_ROOT,
                   'derivation': 'OUTPUT_ROOT/<capture>__cfg<config sha[:12]>__'
                                 'impl<implementation digest[:12]>',
                   'new_directory_required': True, 'reuse_forbidden': True},
        # Bound by the config hash: these are reviewed ACCEPTANCE criteria, so a
        # permit must not survive someone relaxing them.
        'graph_acceptance': {
            'capture_point': GRAPH_CAPTURE_POINT,
            'constructions_per_prior': GRAPH_CONSTRUCTIONS_PER_PRIOR,
            'columns_rule': GRAPH_COLUMNS_RULE,
            'expected_shape': list(expected_graph_shape()),
            'enforced': ['construction count', 'n_neighbors == k_adj',
                         'shape', 'non-empty checksums',
                         'topology identical across priors'],
            'recorded_not_enforced': ['distance checksums -- invariance across '
                                      'priors is not claimed',
                                      'self-first column -- measured False on '
                                      'this path'],
        },
        'graph_preflight_verified': False,
        'graph_preflight_reason':
            'the kNN graph is built inside calculate_doublet_scores, so the '
            'score-free preflight could not construct or check it. A scoring run '
            'would be the first time this graph exists; it is unverified.',
        'permit': {'required': True, 'present': False,
                   'env_var': PERMIT_PATH_ENV,
                   'must_be_tracked_and_identical_to_head': True,
                   'must_match': ['capture', 'config_sha256', 'input_sha256',
                                  'implementation_sha256', 'output_root']},
    }


def scoring_config():
    """The config. Complete, bound to measured evidence, and NOT runnable."""
    cfg = _config_body()
    cfg['output'] = dict(cfg['output'], directory=output_directory())
    return cfg


def _git(*args):
    import subprocess
    return subprocess.run(['git', '-C', _repo_root(), *args],
                          capture_output=True, text=True)


def _head_commit():
    r = _git('rev-parse', 'HEAD')
    if r.returncode != 0:
        raise Gate2RunnerError('cannot resolve deployed HEAD')
    return r.stdout.strip()


def _assert_implementation_commit(want_commit):
    """The permit names the reviewed implementation commit, which must be an
    ANCESTOR of deployed HEAD -- not equal to it.

    Requiring equality was a deadlock, not a safeguard. A production permit must
    be tracked and identical to HEAD, so it can only arrive in a commit made AFTER
    the implementation commit; deployed HEAD is then the permit commit. Naming the
    permit commit's own hash inside that permit is impossible, because changing the
    permit changes the hash of the commit containing it. No legitimate permit could
    satisfy the old gate.

    Ancestry alone would be too weak: it permits any later commit to have edited
    the implementation. So every IMPLEMENTATION_FILES path must additionally be
    BYTE-IDENTICAL between the reviewed commit and deployed HEAD, compared by git
    blob id. That is what "the permit is for one reviewed implementation" means
    once the permit itself has to live in a later commit.
    """
    if not _git('cat-file', '-e', want_commit + '^{commit}').returncode == 0:
        raise Gate2RunnerError(
            f'permit implementation_commit {want_commit} is not a commit in this '
            'repository')
    head = _head_commit()
    if _git('merge-base', '--is-ancestor', want_commit, head).returncode != 0:
        raise Gate2RunnerError(
            f'permit implementation_commit {want_commit} is not an ancestor of '
            f'deployed HEAD {head}. A permit is for a reviewed implementation that '
            'this deployment actually descends from.')
    drifted = []
    for rel in IMPLEMENTATION_FILES:
        a = _git('rev-parse', f'{want_commit}:{rel}')
        b = _git('rev-parse', f'{head}:{rel}')
        if a.returncode != 0:
            raise Gate2RunnerError(
                f'{rel} does not exist at the reviewed implementation commit '
                f'{want_commit}; the permit names a commit that predates the '
                'implementation it claims to authorise')
        if b.returncode != 0:
            raise Gate2RunnerError(f'{rel} does not exist at deployed HEAD {head}')
        if a.stdout.strip() != b.stdout.strip():
            drifted.append(rel)
    if drifted:
        raise Gate2RunnerError(
            f'implementation changed between the reviewed commit {want_commit} and '
            f'deployed HEAD {head}: {drifted}. A permit does not carry over to '
            'edited code, even when the edit arrived in a later commit.')
    # implementation_sha256() digests the files ON DISK. The comparison above is
    # between two commits, so without this the digest could be satisfied by a
    # working tree that differs from deployed HEAD.
    dirty = [rel for rel in IMPLEMENTATION_FILES
             if _git('diff', '--quiet', 'HEAD', '--', rel).returncode != 0
             or _git('diff', '--cached', '--quiet', 'HEAD', '--', rel).returncode != 0]
    if dirty:
        raise Gate2RunnerError(
            f'implementation files differ from deployed HEAD on disk: {dirty}. '
            'The digest must describe the committed code, not a local edit.')
    return head


def _repo_root():
    return os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        '..', '..', '..'))


def _sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def implementation_sha256():
    """One digest over every file that IS the scoring implementation.

    A permit binds this, so editing the entrypoint, the wrapper, the config or an
    imported module invalidates it. The preflight commit cannot serve this role:
    the scoring entrypoint did not exist at that commit.
    """
    root = _repo_root()
    parts = []
    for rel in IMPLEMENTATION_FILES:
        p = os.path.join(root, rel)
        if not os.path.exists(p):
            raise Gate2RunnerError(f'implementation file missing: {rel}')
        parts.append(f'{rel}:{_sha256_file(p)}')
    return hashlib.sha256('\n'.join(parts).encode()).hexdigest()


def output_directory():
    """Deterministic, and different for any config or implementation change."""
    try:
        impl = implementation_sha256()[:12]
    except Gate2RunnerError:
        impl = 'noimpl'
    # The body, not scoring_config(), or this recurses through output_directory().
    # config_sha256 excludes output.directory, so this equals config_sha256().
    cfg = config_sha256(_config_body())[:12]
    return os.path.join(OUTPUT_ROOT, f'{CAPTURE}__cfg{cfg}__impl{impl}')


def config_sha256(config=None):
    """THE config hash -- one value, whether called on the body or the full config.

    `output.directory` is DERIVED from this hash and is therefore excluded from it.
    Without that exclusion the two call sites diverge: `config_sha256()` hashes the
    full config including the directory, while `output_directory()` must hash the
    body alone to avoid recursion. The directory would then be named after a hash
    that no published field equals, and the documented derivation
    `cfg<config sha[:12]>` would be false -- an auditor could not reproduce the
    output path from the permit's `config_sha256`.
    """
    cfg = dict(scoring_config() if config is None else config)
    cfg.pop('_config_sha256', None)
    out = cfg.get('output')
    if isinstance(out, dict) and 'directory' in out:
        cfg['output'] = {k: v for k, v in out.items() if k != 'directory'}
    blob = json.dumps(cfg, sort_keys=True, separators=(',', ':')).encode()
    return hashlib.sha256(blob).hexdigest()


def assert_inputs_unchanged():
    """Refuse if any deposit file drifted from the hash the evidence rests on."""
    for name, want in INPUT_SHA256.items():
        p = os.path.join(DEPOSIT_DIR, name)
        if not os.path.exists(p):
            raise Gate2RunnerError(f'input missing: {p}')
        h = hashlib.sha256()
        with open(p, 'rb') as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b''):
                h.update(chunk)
        got = h.hexdigest()
        if got != want:
            raise Gate2RunnerError(
                f'{name}: sha256 {got} != the audited {want}. The C1-C3 evidence '
                'does not describe this file.')
    return True


def _assert_permit_is_repository_artifact(path):
    """Production permits live in the deployed repo, tracked and identical to HEAD.

    Without this, `$GATE2_PILOT_PERMIT` pointing at any writable path is a permit
    factory: anyone who can read the config hashes can mint a matching JSON in
    /tmp. A test substitutes this function through `assert_permitted(
    repo_verifier=...)`; there is no environment variable that switches it off.
    """
    import subprocess
    root = _repo_root()
    real = os.path.realpath(path)
    if not real.startswith(os.path.realpath(root) + os.sep):
        raise Gate2RunnerError(
            f'permit {real} is outside the deployed repository {root}; a permit '
            'must be a tracked repository artifact, not an arbitrary file')
    rel = os.path.relpath(real, root)

    def git(*a):
        return subprocess.run(['git', '-C', root, *a], capture_output=True, text=True)

    if git('ls-files', '--error-unmatch', rel).returncode != 0:
        raise Gate2RunnerError(f'permit {rel} is not tracked')
    if git('diff', '--quiet', '--', rel).returncode != 0:
        raise Gate2RunnerError(f'permit {rel} has unstaged changes')
    if git('diff', '--cached', '--quiet', '--', rel).returncode != 0:
        raise Gate2RunnerError(f'permit {rel} has staged changes')
    return rel


def assert_permitted(permit_path=None, repo_verifier=None):
    """Refuse to score without a permit that matches THIS config and THESE inputs.

    `repo_verifier` exists so a test can substitute the repository-artifact check.
    Production calls `assert_permitted()` with no arguments and therefore always
    gets the real verifier; there is no runtime flag that weakens it.
    """
    if os.environ.get(REMOVED_BYPASS_ENV) is not None:
        raise Gate2RunnerError(
            f'{REMOVED_BYPASS_ENV} is set. That bypass was removed; its presence in '
            'the environment is treated as an attempt to disable the permit check.')
    verify = _assert_permit_is_repository_artifact if repo_verifier is None else repo_verifier
    cfg = scoring_config()
    if not cfg['graph_preflight_verified'] and not cfg['graph_preflight_reason']:
        raise Gate2RunnerError('graph limitation must be stated, not merely false')
    path = permit_path or os.environ.get(PERMIT_PATH_ENV)
    if not path:
        raise Gate2RunnerError(
            'no scoped pilot permit: scoring is not authorised. Set '
            f'{PERMIT_PATH_ENV} to a permit naming this config and these inputs.')
    if not os.path.exists(path):
        raise Gate2RunnerError(f'permit not found: {path}')
    verify(path)
    with open(path) as fh:
        permit = json.load(fh)
    missing = [f for f in PERMIT_REQUIRED_FIELDS if f not in permit]
    if missing:
        raise Gate2RunnerError(f'permit is missing required field(s): {missing}')
    if permit['capture'] != CAPTURE:
        raise Gate2RunnerError(
            f'permit is for capture {permit["capture"]!r}, this config is {CAPTURE!r}')
    want_cfg = config_sha256(cfg)
    if permit['config_sha256'] != want_cfg:
        raise Gate2RunnerError(
            f'permit config_sha256 {permit["config_sha256"]} != this config {want_cfg}. '
            'A permit does not carry over to a changed config.')
    if dict(permit['input_sha256']) != dict(INPUT_SHA256):
        raise Gate2RunnerError('permit input_sha256 does not match the audited inputs')
    want_impl = implementation_sha256()
    if permit['implementation_sha256'] != want_impl:
        raise Gate2RunnerError(
            f'permit implementation_sha256 {permit["implementation_sha256"]} != the '
            f'deployed implementation {want_impl}. A permit does not carry over to '
            'edited code.')
    if permit['output_root'] != OUTPUT_ROOT:
        raise Gate2RunnerError(
            f'permit output_root {permit["output_root"]!r} != {OUTPUT_ROOT!r}')
    if permit.get('scope') != REQUIRED_SCOPE:
        raise Gate2RunnerError(
            f'permit scope {permit.get("scope")!r} != the only authorised scope '
            f'{REQUIRED_SCOPE!r}')
    ceiling = permit.get('resource_ceiling')
    if ceiling != dict(REQUIRED_RESOURCE_CEILING):
        raise Gate2RunnerError(
            f'permit resource_ceiling {ceiling!r} != {dict(REQUIRED_RESOURCE_CEILING)!r}')
    want_commit = permit.get('implementation_commit')
    if not want_commit:
        raise Gate2RunnerError('permit states no implementation_commit')
    _assert_implementation_commit(want_commit)
    return permit


def scoped_ruling(permit):
    """The reviewed ruling, handed over ONLY after the permit gate succeeded.

    The global PRODUCTION_PARAMETER_RULING stays {} -- activating it would make
    every capture runnable, which no permit authorises.
    """
    if not permit:
        raise Gate2RunnerError('scoped ruling requested without a validated permit')
    return dict(METHODOLOGY_PARAMETER_RULING)


if __name__ == '__main__':          # pragma: no cover - inspection helper
    cfg = scoring_config()
    cfg['_config_sha256'] = config_sha256(cfg)
    print(json.dumps(cfg, indent=2, sort_keys=True))
