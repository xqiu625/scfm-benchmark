import ast
import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import socket
import subprocess
import sys
import time

root = Path(sys.argv[1]).resolve()
out = Path(sys.argv[2]).resolve()
out.mkdir(parents=True, exist_ok=False)
os.chdir(root)
frozen = json.loads((root / 'FROZEN_BYTES.json').read_text())
env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='4',
           OPENBLAS_NUM_THREADS='4', MKL_NUM_THREADS='4', NUMBA_NUM_THREADS='4')
env['PYTHONPYCACHEPREFIX'] = str(out / 'bytecode')
summary = {'started_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
           'host': socket.gethostname(), 'python': sys.version, 'worktree': str(root),
           'head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
           'phases': [], 'mutations': [], 'complete': False}


def hashes():
    return {p: hashlib.sha256((root / p).read_bytes()).hexdigest() for p in frozen['files']}


def verify(label):
    current = hashes()
    (out / (label + '_SHA256.json')).write_text(json.dumps(current, indent=2) + '\n')
    drift = [p for p, h in frozen['files'].items() if current[p] != h]
    if drift:
        raise RuntimeError(f'{label}: frozen bytes changed: {drift}')
    print(f'{label}: {len(current)} frozen file hashes MATCH', flush=True)


def save():
    (out / 'SUMMARY.json').write_text(json.dumps(summary, indent=2) + '\n')


def run(label, selection):
    cmd = [sys.executable, '-m', 'pytest', *selection, '-q', '-p', 'no:cacheprovider']
    start = time.monotonic()
    print(f'{label}: {shlex.join(cmd)}', flush=True)
    with (out / (label + '.log')).open('w') as log:
        log.write('$ ' + shlex.join(cmd) + '\n'); log.flush()
        result = subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT)
        log.write(f'\nEXIT_CODE={result.returncode}\n')
    record = {'label': label, 'command': cmd, 'exit_code': result.returncode,
              'elapsed_seconds': round(time.monotonic() - start, 2)}
    print(f'{label}: exit {result.returncode}, {record["elapsed_seconds"]} seconds', flush=True)
    return record


try:
    verify('BEFORE_TESTS')
    import scanpy, anndata, numpy, scipy, sklearn, pytest
    summary['versions'] = {m.__name__: m.__version__ for m in
                           (scanpy, anndata, numpy, scipy, sklearn, pytest)}
    for label, selection in [
        ('focused', ['tests/test_gate2_scoring_config.py', 'tests/test_gate2_scoring_wrapper.py']),
        ('scanpy_integration', ['tests/', '-m', 'scanpy_integration']),
        ('full_suite', ['tests/'])]:
        rec = run(label, selection)
        summary['phases'].append(rec); save()
        verify('AFTER_' + label.upper())
        if rec['exit_code'] != 0:
            raise RuntimeError(f'{label} failed; no later approval is possible')

    code = 'scripts/32_rebuild_upfront_qc/code/'
    config = code + 'gate2_scoring_config_GSM5102900_HC1.py'
    entry = code + 'gate2_score_capture_GSM5102900_HC1.py'
    mutations = [
        ('M4', entry, 'if len(per_call) != expected_n:', 'if False:',
         'tests/test_gate2_scoring_wrapper.py::test_zero_captures_is_a_failure_not_a_pass'),
        ('M5', entry, "if g['n_neighbors_arg'] != k_adj_expected:", 'if False:',
         'tests/test_gate2_scoring_wrapper.py::test_wrong_k_is_enforced_not_merely_recorded'),
        ('M6', entry, "if tuple(g['shape']) != tuple(expected_shape):", 'if False:',
         'tests/test_gate2_scoring_wrapper.py::test_wrong_shape_is_enforced'),
        ('M7', entry, 'if per_call and len(idx_shas) != 1:', 'if False:',
         'tests/test_gate2_scoring_wrapper.py::test_topology_differing_across_priors_is_enforced'),
        ('M8', entry, 'if per_call and len(idx_shas) != 1:',
         'if per_call and (len(idx_shas) != 1 or len(dist_shas) != 1):',
         'tests/test_gate2_scoring_wrapper.py::test_distances_are_recorded_but_not_claimed_invariant'),
        ('M9', config, "if _git('merge-base', '--is-ancestor', want_commit, head).returncode != 0:",
         'if head != want_commit:',
         'tests/test_gate2_scoring_config.py::test_the_permit_commit_gate_is_satisfiable_at_all'),
        ('M10', config, 'if drifted:', 'if False:',
         'tests/test_gate2_scoring_config.py::test_an_implementation_edited_after_the_reviewed_commit_is_refused'),
        ('M11', config, 'if dirty:', 'if False:',
         'tests/test_gate2_scoring_config.py::test_an_implementation_dirty_on_disk_is_refused'),
        ('M12', config, 'if a.returncode != 0:', 'if False:',
         'tests/test_gate2_scoring_config.py::test_a_commit_predating_the_implementation_is_refused'),
    ]
    for label, rel, old, new, test in mutations:
        path = root / rel
        baseline = path.read_bytes()
        assert hashlib.sha256(baseline).hexdigest() == frozen['files'][rel]
        source = baseline.decode()
        assert source.count(old) == 1, (label, source.count(old))
        mutated = source.replace(old, new, 1).encode()
        ast.parse(mutated)
        assert mutated != baseline
        try:
            path.write_bytes(mutated)
            assert path.read_bytes() == mutated
            rec = run(label, [test])
            rec.update(target=rel, test=test, baseline_sha256=hashlib.sha256(baseline).hexdigest(),
                       mutation_sha256=hashlib.sha256(mutated).hexdigest(), injected=True)
            rec['killed'] = rec['exit_code'] == 1 and ('FAILED ' + test) in (out / (label + '.log')).read_text()
        finally:
            path.write_bytes(baseline)
            assert path.read_bytes() == baseline
        rec['restored_byte_identical'] = True
        summary['mutations'].append(rec); save()
        verify('AFTER_' + label + '_RESTORE')
        if not rec['killed']:
            raise RuntimeError(f'{label} was not killed by its intended assertion')
    verify('FINAL')
    summary['complete'] = True
except BaseException as exc:
    summary['error'] = f'{type(exc).__name__}: {exc}'
    raise
finally:
    summary['ended_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    save()
