"""Gate 2 checkpoints: one schema for both runners, bound to a full fingerprint.

The previous checkpoints were keyed by a sanitized channel name and nothing
else. Resuming after any change -- a corrected rate, a different seed, new code
-- silently reused results computed under the old settings, and there was no
way to tell from the file that it had happened. Worse, the sanitizer mapped
`A-B` and `A_B` onto the same filename, so two captures could overwrite each
other.

Everything that can change a score is therefore hashed into a fingerprint, and
a restore with ANY mismatched field is refused. There is no warn-and-continue
path: a mismatch means the stored result does not answer the question now being
asked.

CONCURRENCY. Temp files are unique per call (pid + uuid4, created with O_EXCL),
so two writers can never share or delete each other's temp. The FINAL file is
last-writer-wins: `os.replace` is atomic, so a reader sees one complete
checkpoint or the other, never a mixture, but no lock arbitrates between two
writers holding DIFFERENT fingerprints for the same capture -- the later one
simply stands, and a subsequent run whose fingerprint matches neither will
refuse it. **Running two writers against one --out-dir is therefore a
precondition violation, not a supported mode.** A lock file would be the fix if
that ever becomes necessary; it is deliberately not added for a single-writer
runner.

Layout per capture: `<dir>/<collision-resistant key>.npz`, written to a
same-directory temp file, flushed and fsynced, then `os.replace`d, so a partial
file is never visible under the real name. Metadata travels as a JSON string,
not a pickled object array -- `allow_pickle` on a file that may be stale or
hand-edited is an unnecessary execution risk.
"""
import hashlib
import json
import os
import uuid

import numpy as np

from gate2_diagnostics import DIAGNOSTICS_SCHEMA_VERSION
from gate2_status import (CALL_STATUSES, validate_call_statuses,
                          validate_capture_invariants, validate_capture_status)

#: /1 -> /2: cell IDs were stored as a fixed 'U64' array, which silently
#: truncated anything longer and made the capture permanently unrestorable
#: (the fingerprint held the full checksum, so the stored IDs could never match
#: again) while the error blamed 'identity/order'.
#:
#: /2 -> /3: the call_status value domain changed. A /2 file contains the
#: literal 'auto_singlet', which is no longer a member of the enum, so
#: validate_call_statuses would reject it -- but with a message about an
#: invalid status rather than an outdated format. The fingerprint already
#: carries call_status_schema and would refuse these files anyway; the version
#: bump exists so the REASON reported is the true one.
#:
#: /3 -> /4: the diagnostics dict persisted inside the checkpoint gained the
#: score-distribution fields (AUC, 101 percentiles, shared-bin histograms,
#: threshold-based tail fractions). Diagnostics are restored verbatim on
#: resume, so a /3 file would come back missing them and a resumed run would
#: silently disagree with a fresh one on exactly the fields the sensitivity
#: analysis depends on. The diagnostics field set also carries its own
#: DIAGNOSTICS_SCHEMA_VERSION, so a later change to the field DEFINITIONS is
#: caught even when the checkpoint schema itself is unchanged.
#:
#: /1, /2 and /3 are all rejected outright. No migration path exists by design.
CHECKPOINT_SCHEMA_VERSION = 'gate2-checkpoint/4'
REJECTED_SCHEMA_VERSIONS = ('gate2-checkpoint/1', 'gate2-checkpoint/2',
                            'gate2-checkpoint/3')

#: Fields whose change must invalidate a stored capture. Order-independent; the
#: fingerprint is a hash of the sorted JSON.
FINGERPRINT_FIELDS = (
    'checkpoint_schema_version',
    'call_status_schema',
    'cohort',
    'channel_id_raw',
    'n_cells',
    'cell_id_checksum',
    'matrix_fingerprint',
    'config_row_checksum',
    'expected_doublet_rate',
    'donor_proportions_checksum',
    'donor_proportion_stage',
    'platform',
    'base_seed',
    'capture_seed',
    'scrublet_params',
    'software_versions',
    'git_commit',
    'code_fingerprint',
)


class CheckpointMismatch(RuntimeError):
    """Stored checkpoint does not match the current run. Always refuse."""


def _exact_unicode_dtype(values):
    """Unicode dtype wide enough for the longest value, so nothing truncates."""
    width = max((len(str(v)) for v in values), default=1)
    return f'U{max(width, 1)}'


def _cell_id_array(cell_ids):
    ids = [str(c) for c in cell_ids]
    return np.asarray(ids, dtype=_exact_unicode_dtype(ids))


def capture_key(cohort, channel_id_raw):
    """Collision-resistant per-capture filename stem.

    A sanitizing transform is not usable here: `A-B` and `A_B` both normalize to
    `A_B` and would share a file. The raw identity is hashed instead, and a
    readable but non-authoritative prefix is kept only for human browsing.
    """
    digest = hashlib.sha256(f'{cohort}\x00{channel_id_raw}'.encode()).hexdigest()[:24]
    readable = ''.join(c if c.isalnum() else '-' for c in str(channel_id_raw))[:40]
    return f'{readable}.{digest}'


class CellIdentityError(RuntimeError):
    """The cell-ID vector cannot carry identity. Always fatal."""


def validate_cell_ids(cell_ids, n_expected, *, cohort='?', channel='?'):
    """Every identity claim in this layer rests on these IDs, so they are
    checked before any checksum, scoring or checkpoint I/O.

    Identity is taken VERBATIM: no trimming, no case folding, no unicode
    normalization. Rewriting an ID to make it match would defeat the purpose,
    so near-collisions that differ only by surrounding whitespace or case are
    reported as errors rather than silently merged.
    """
    where = f'cohort={cohort!r} capture={channel!r}'
    if len(cell_ids) != n_expected:
        raise CellIdentityError(
            f'{where}: {len(cell_ids)} cell IDs but the matrix has {n_expected} cells')
    if n_expected == 0:
        raise CellIdentityError(f'{where}: no cells')

    seen, folded = set(), {}
    for i, raw in enumerate(cell_ids):
        if raw is None or (isinstance(raw, float) and raw != raw):   # None / NaN
            raise CellIdentityError(f'{where}: cell ID at position {i} is missing')
        cid = str(raw)
        if cid.strip() == '' or cid.lower() in ('nan', 'none', '<na>'):
            raise CellIdentityError(f'{where}: cell ID at position {i} is empty or null-like ({cid!r})')
        if cid in seen:
            raise CellIdentityError(f'{where}: duplicate cell ID {cid!r}')
        seen.add(cid)
        key = cid.strip().casefold()
        if key in folded:
            raise CellIdentityError(
                f'{where}: cell IDs {folded[key]!r} and {cid!r} collide after whitespace/case '
                'normalization. They are kept distinct rather than merged; resolve upstream.')
        folded[key] = cid
    return [str(c) for c in cell_ids]


def ordered_cell_checksum(cell_ids):
    """Checksum of cell identity AND order. Row count alone cannot detect a
    reordered matrix, which would silently misalign every score. Callers must
    have passed the IDs through validate_cell_ids first."""
    h = hashlib.sha256()
    for cid in cell_ids:
        h.update(str(cid).encode())
        h.update(b'\x00')
    return h.hexdigest()


def matrix_fingerprint(mat_csc_or_csr):
    """Shape, nnz and the content of the sparse arrays. Cheap relative to
    scoring, and the only way a changed matrix invalidates a checkpoint."""
    m = mat_csc_or_csr
    h = hashlib.sha256()
    h.update(str(m.shape).encode())
    h.update(str(m.nnz).encode())
    for attr in ('data', 'indices', 'indptr'):
        arr = getattr(m, attr, None)
        if arr is not None:
            h.update(np.ascontiguousarray(arr).tobytes())
    return h.hexdigest()


def donor_proportions_checksum(proportions):
    """Bind the donor mixture, rounded so float noise does not spuriously
    invalidate, but finely enough that a real composition change does."""
    if proportions is None:
        return 'none'
    arr = np.asarray(proportions, dtype=float)
    arr = np.sort(arr[arr > 0])[::-1]
    return hashlib.sha256(','.join(f'{v:.9f}' for v in arr).encode()).hexdigest()


def capture_seed(base_seed, cohort, channel_id_raw):
    """Deterministic per-capture seed. Explicit base seed, no implicit default;
    derived from identity so it is reproducible and recorded alongside the base."""
    digest = hashlib.sha256(f'{int(base_seed)}\x00{cohort}\x00{channel_id_raw}'.encode()).digest()
    return int.from_bytes(digest[:4], 'big') % (2 ** 31 - 1)


def code_fingerprint(paths):
    """Hash of the sources that decide the result, so an edited scorer
    invalidates checkpoints even without a commit."""
    h = hashlib.sha256()
    for p in sorted(paths):
        h.update(os.path.basename(p).encode())
        try:
            with open(p, 'rb') as fh:
                h.update(hashlib.sha256(fh.read()).digest())
        except OSError:
            h.update(b'<unreadable>')
    return h.hexdigest()


def software_versions():
    out = {}
    for name in ('scanpy', 'numpy', 'scipy', 'anndata', 'sklearn', 'skimage'):
        try:
            mod = __import__(name)
            out[name] = str(getattr(mod, '__version__', 'unknown'))
        except Exception:
            out[name] = 'absent'
    try:
        import scrublet
        out['scrublet'] = str(getattr(scrublet, '__version__', 'unknown'))
    except Exception:
        # scanpy>=1.9 vendors Scrublet; absence of the standalone package is
        # normal and is recorded rather than passed over.
        out['scrublet'] = 'vendored-by-scanpy-or-absent'
    return out


def build_fingerprint(**kwargs):
    missing = [f for f in FINGERPRINT_FIELDS if f not in kwargs]
    if missing:
        raise ValueError(f'fingerprint missing fields: {missing}')
    extra = sorted(set(kwargs) - set(FINGERPRINT_FIELDS))
    if extra:
        raise ValueError(f'fingerprint has unexpected fields: {extra}')
    return {k: kwargs[k] for k in FINGERPRINT_FIELDS}


def _fingerprint_digest(fp):
    return hashlib.sha256(
        json.dumps(fp, sort_keys=True, separators=(',', ':'), default=str).encode()).hexdigest()


def save_checkpoint(ckpt_dir, key, fingerprint, *, scores, call_status, candidate,
                    cell_ids, capture_status, diagnostics):
    """Atomic write. Temp file in the SAME directory (so os.replace stays atomic
    on one filesystem), flushed and fsynced before the rename."""
    os.makedirs(ckpt_dir, exist_ok=True)
    # The full contract, not field-by-field: the same function the fresh path and
    # the restore path use, so all three agree by construction.
    validate_capture_invariants(capture_status, call_status, candidate, scores,
                                threshold=(diagnostics or {}).get('threshold'),
                                diagnostics=diagnostics, where=f'save:{key}')

    meta = json.dumps({
        'fingerprint': fingerprint,
        'fingerprint_digest': _fingerprint_digest(fingerprint),
        'capture_status': capture_status,
        'diagnostics': diagnostics,
    }, sort_keys=True, default=str)

    final = os.path.join(ckpt_dir, f'{key}.npz')
    # PID alone collides between threads, retries and nested calls in ONE
    # process. uuid4 makes the path unique per call, so a failed save can
    # only ever delete its own temp.
    tmp = os.path.join(ckpt_dir, f'.{key}.{os.getpid()}.{uuid.uuid4().hex}.npz.tmp')
    replaced = False
    try:
        # O_EXCL: create-or-fail, so the path cannot be an existing file another
        # writer is using. uuid4 already makes that essentially impossible; the
        # flag removes the check-then-open window entirely.
        fd = os.open(tmp, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
        with os.fdopen(fd, 'wb') as fh:
            np.savez(
                fh,
                meta=np.array(meta),
                scores=np.asarray(scores, dtype=float),
                call_status=np.asarray([str(s) for s in call_status],
                                       dtype=_exact_unicode_dtype(CALL_STATUSES)),
                # candidate is tri-state: 1 doublet, 0 not-doublet, -1 unknown.
                # A plain bool array cannot express "unknown", which is the whole point.
                # NOT np.asarray(..., dtype=np.int8): that would launder strings,
            # floats and 2-D input into a valid-looking array. The invariant
            # check above already required a 1-D int8 array, so this is the
            # object itself.
            candidate=candidate,
                # dtype derived from the longest ID present -- never a fixed width.
                cell_ids=_cell_id_array(cell_ids),
            )
            fh.flush()
            os.fsync(fh.fileno())
        # Read the temp back and confirm the IDs survived byte-for-byte before
        # the file is allowed to become the checkpoint. A truncating dtype
        # would otherwise only surface much later, as a misleading
        # "identity/order differs" on every future restore.
        with np.load(tmp, allow_pickle=False) as _chk:
            if [str(c) for c in _chk['cell_ids']] != [str(c) for c in cell_ids]:
                raise CheckpointMismatch(
                    f'{final}: cell IDs did not survive serialization intact; refusing to '
                    'install a checkpoint that could never be restored')
        os.replace(tmp, final)
        replaced = True
        # The rename itself needs the directory entry flushed, or a crash can
        # lose a checkpoint that np.savez already fsynced. If THIS fails the
        # checkpoint is already visible and its contents are already durable,
        # so the write succeeded -- reporting failure here would tell the caller
        # nothing was written when something was. It is recorded as a warning
        # instead, and the file stands.
        try:
            dfd = os.open(ckpt_dir, os.O_RDONLY)
            try:
                os.fsync(dfd)
            finally:
                os.close(dfd)
        except OSError as exc:
            print(f'WARNING: {final} was written and replaced, but the directory fsync '
                  f'failed ({exc}); the checkpoint is valid but its rename may not '
                  'survive a power loss', flush=True)
    finally:
        # Only this call's temp file, named with our pid, is ever removed --
        # a concurrent writer's temp is not ours to delete.
        if not replaced and os.path.exists(tmp):
            os.unlink(tmp)
    return final


def load_checkpoint(ckpt_dir, key, fingerprint, *, expected_cell_ids):
    """Return the stored capture, or None if absent. Raise CheckpointMismatch on
    ANY disagreement -- a stale checkpoint is refused, never warned about."""
    path = os.path.join(ckpt_dir, f'{key}.npz')
    if not os.path.exists(path):
        return None
    try:
        with np.load(path, allow_pickle=False) as z:
            meta = json.loads(str(z['meta']))
            scores = z['scores']
            call_status = [str(s) for s in z['call_status']]
            candidate = z['candidate']
            cell_ids = [str(c) for c in z['cell_ids']]
    except Exception as exc:
        # Covers a truncated or half-written file: unreadable is refused, and the
        # caller rescores rather than proceeding on a partial result.
        raise CheckpointMismatch(f'{path}: unreadable or corrupt checkpoint ({exc})')

    stored_fp = meta.get('fingerprint', {})
    if meta.get('fingerprint_digest') != _fingerprint_digest(stored_fp):
        raise CheckpointMismatch(f'{path}: fingerprint digest does not match its own fingerprint')

    stored_schema = stored_fp.get('checkpoint_schema_version')
    if stored_schema != CHECKPOINT_SCHEMA_VERSION:
        why = (' (a superseded schema; no migration path exists)'
               if stored_schema in REJECTED_SCHEMA_VERSIONS else '')
        raise CheckpointMismatch(
            f'{path}: checkpoint_schema_version={stored_schema!r}{why}, this code writes '
            f'{CHECKPOINT_SCHEMA_VERSION!r}. Superseded checkpoints are rejected outright and '
            'are never migrated or read.')

    differing = [f for f in FINGERPRINT_FIELDS if stored_fp.get(f) != fingerprint.get(f)]
    if differing:
        raise CheckpointMismatch(
            f'{path}: fingerprint mismatch on {differing}. The stored result was computed under '
            'different inputs, configuration, seed, software or code and cannot be reused.')

    n = len(expected_cell_ids)
    if scores.ndim != 1:
        raise CheckpointMismatch(f'{path}: scores has ndim={scores.ndim}, expected 1-D')
    if scores.dtype.kind not in 'fi' or scores.dtype == bool:
        raise CheckpointMismatch(
            f'{path}: scores dtype {scores.dtype} is not numeric (object/string/bool refused)')
    if not np.all(np.isfinite(scores)) and meta.get('capture_status') == 'scored':
        raise CheckpointMismatch(
            f'{path}: capture_status=scored but scores contain non-finite values')
    if not (len(scores) == len(call_status) == len(candidate) == len(cell_ids) == n):
        raise CheckpointMismatch(
            f'{path}: restored array lengths '
            f'{(len(scores), len(call_status), len(candidate), len(cell_ids))} do not all equal '
            f'the expected {n}')
    if cell_ids != [str(c) for c in expected_cell_ids]:
        raise CheckpointMismatch(f'{path}: restored cell identity/order differs from the input')

    # The diagnostics carry their own version because the checkpoint schema can
    # stay put while a metric DEFINITION changes -- AUC direction, percentile
    # count, histogram binning. Writing that version without checking it on
    # restore left the contract open at exactly the point it was meant to close.
    stored_diag = meta.get('diagnostics') or {}
    stored_diag_version = stored_diag.get('diagnostics_schema_version')
    if stored_diag_version is None:
        raise CheckpointMismatch(
            f'{path}: diagnostics carry no diagnostics_schema_version. Written before the '
            'diagnostics contract existed; refused rather than assumed compatible.')
    if stored_diag_version != DIAGNOSTICS_SCHEMA_VERSION:
        raise CheckpointMismatch(
            f'{path}: diagnostics_schema_version={stored_diag_version!r}, this code writes '
            f'{DIAGNOSTICS_SCHEMA_VERSION!r}. The metric definitions differ; a stored summary '
            'computed under other definitions is not comparable and is never migrated.')

    if candidate.dtype != np.int8:
        raise CheckpointMismatch(
            f'{path}: candidate array dtype {candidate.dtype} invalid, expected int8')
    # The whole contract, so a stored file cannot claim one status while its
    # per-cell arrays say another -- the gap that let a corrupt or hand-edited
    # checkpoint restore "unknown" cells as singlets.
    try:
        validate_capture_invariants(
            meta.get('capture_status'), call_status, candidate, scores,
            threshold=(meta.get('diagnostics') or {}).get('threshold'),
            diagnostics=meta.get('diagnostics'), where=path)
    except Exception as exc:
        raise CheckpointMismatch(f'{path}: stored capture violates the status contract: {exc}')

    return {
        'scores': scores,
        'call_status': call_status,
        'candidate': candidate,
        'cell_ids': cell_ids,
        'capture_status': meta['capture_status'],
        'diagnostics': meta['diagnostics'],
    }
