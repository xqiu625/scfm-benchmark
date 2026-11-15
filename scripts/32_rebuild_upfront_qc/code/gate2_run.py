"""Gate 2 runner — one implementation for both cohort layouts.

Replaces gate2_scrublet_score.py (combined) and gate2_scrublet_split.py (split),
which had drifted into different checkpoint schemas, different prior handling
and different resume semantics. Layout is now an input-reading detail; status
model, checkpointing and diagnostics are identical either way.

Scoring only. Nothing is deleted here. `call_status == "predicted_doublet"`
marks a MECHANICAL CANDIDATE, which is necessary but not sufficient for
removal: cohort eligibility, an approved runtime rate, a capture reliability
criterion and a sensitivity analysis are all still required, and no
final-authorization artifact exists in this repository.

**A production run cannot currently start.** `--runtime-config` is mandatory and
no authorized configuration exists: every T7 row is runtime_eligible=false and
rate selection is unresolved. This is the intended behaviour, not a gap.

    python gate2_run.py --cohort <name> --layout combined|split \
        --input <dir> --out-dir <dir> --runtime-config <csv> --base-seed <int>
"""
import argparse
import csv
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from gate2_checkpoint import (CHECKPOINT_SCHEMA_VERSION, CellIdentityError, CheckpointMismatch,
                              build_fingerprint,
                              capture_key, capture_seed, code_fingerprint,
                              donor_proportions_checksum, load_checkpoint, matrix_fingerprint,
                              ordered_cell_checksum, save_checkpoint, software_versions,
                              validate_cell_ids)
from gate2_runtime_config import Gate2ConfigError, load_runtime_config
from gate2_score_core import SCRUBLET_PARAMS, score_one_capture
from gate2_status import CALL_STATUSES, mechanical_candidate_mask

#: Every module whose content can change a stored result. gate2_diagnostics.py
#: belongs here: editing a metric definition changes the persisted diagnostics
#: without touching any other file, and a checkpoint written under the old
#: definition would otherwise still restore.
CODE_FILES = [os.path.join(HERE, f) for f in
              ('gate2_run.py', 'gate2_score_core.py', 'gate2_checkpoint.py',
               'gate2_runtime_config.py', 'gate2_status.py', 'gate2_diagnostics.py')]


def _git_state():
    """Distinguish the three cases the previous version folded into 'unknown'.

    A dirty worktree matters: the commit alone would claim provenance the
    working files do not have. code_fingerprint covers the five gate2 modules,
    but not anything else that changed, so the dirty flag is recorded too.
    """
    import subprocess

    def _run(args):
        return subprocess.run(args, cwd=HERE, capture_output=True, text=True, timeout=10)

    try:
        head = _run(['git', 'rev-parse', 'HEAD'])
    except FileNotFoundError:
        return 'git-unavailable'
    except Exception as exc:
        return f'git-error:{type(exc).__name__}'
    if head.returncode != 0:
        return 'not-a-git-repository'
    commit = head.stdout.strip()
    if not commit:
        return 'not-a-git-repository'
    try:
        st = _run(['git', 'status', '--porcelain'])
        dirty = bool(st.stdout.strip()) if st.returncode == 0 else None
    except Exception:
        dirty = None
    if dirty is None:
        return f'{commit}+dirty-unknown'
    return f'{commit}+dirty' if dirty else commit


def derived_cell_key(cohort, channel_id, cell_id):
    """A DERIVED convenience column, not the authoritative identity.

    Canonical JSON array: injective for any content, including '|', commas,
    quotes, whitespace and non-ASCII. Consumers must still join on the three
    columns; nothing -- least of all a removal decision -- may depend on this
    string.
    """
    return json.dumps([str(cohort), str(channel_id), str(cell_id)],
                      ensure_ascii=False, separators=(',', ':'))


def _require_cell_id(cells, path, cohort):
    """cell_id must be present and named. The previous fallback to
    `cells.columns[0]` was a silent guess at an unnamed index or barcode
    column, and every identity guarantee in this layer is built on these
    values -- the same defect class as the prior fallbacks that were removed."""
    if 'cell_id' not in cells.columns:
        raise CellIdentityError(
            f'cohort={cohort!r}: {path} has no cell_id column (found: '
            f'{list(cells.columns)[:8]}). Gate 2 will not guess an index, barcode or '
            'first column; add an explicit cell_id upstream.')


def iter_captures_combined(input_dir, cohort):
    """One counts.mtx + cells.csv for the whole cohort."""
    import pandas as pd
    from scipy.io import mmread
    path = os.path.join(input_dir, 'cells.csv')
    mat = mmread(os.path.join(input_dir, 'counts.mtx')).tocsc()
    cells = pd.read_csv(path, low_memory=False)
    if mat.shape[1] != len(cells):
        raise RuntimeError(f'matrix has {mat.shape[1]} columns but cells.csv has {len(cells)} rows')
    _require_cell_id(cells, path, cohort)
    id_col = 'cell_id'
    for ch in sorted(cells.channel_id.astype(str).unique()):
        idx = np.where((cells.channel_id.astype(str) == ch).to_numpy())[0]
        sub = cells.iloc[idx]
        props = (sub['donor_id'].value_counts(normalize=True).to_numpy()
                 if 'donor_id' in sub.columns else None)
        sub_mat = mat[:, idx].T.tocsr()
        ids = validate_cell_ids(sub[id_col].tolist(), sub_mat.shape[0],
                                cohort=cohort, channel=ch)
        yield ch, sub_mat, ids, props


def iter_captures_split(input_dir, cohort):
    """split_per_channel/<channel>/{counts.mtx,cells.csv}."""
    import pandas as pd
    from scipy.io import mmread
    root = os.path.join(input_dir, 'split_per_channel')
    for ch in sorted(d for d in os.listdir(root) if os.path.isdir(os.path.join(root, d))):
        path = os.path.join(root, ch, 'cells.csv')
        cells = pd.read_csv(path, low_memory=False)
        if not len(cells):
            continue
        mat = mmread(os.path.join(root, ch, 'counts.mtx')).tocsc()
        if mat.shape[1] != len(cells):
            raise RuntimeError(f'{ch}: matrix/cells length mismatch')
        _require_cell_id(cells, path, cohort)
        props = (cells['donor_id'].value_counts(normalize=True).to_numpy()
                 if 'donor_id' in cells.columns else None)
        sub_mat = mat.T.tocsr()
        ids = validate_cell_ids(cells['cell_id'].tolist(), sub_mat.shape[0],
                                cohort=cohort, channel=str(ch))
        yield str(ch), sub_mat, ids, props


def run(cohort, layout, input_dir, out_dir, runtime_config_path, base_seed, scrublet_fn=None):
    os.makedirs(out_dir, exist_ok=True)
    ckpt_dir = os.path.join(out_dir, 'checkpoints')

    iterator = iter_captures_combined if layout == 'combined' else iter_captures_split
    captures = list(iterator(input_dir, cohort))
    channels = [c[0] for c in captures]

    # Fail-closed: raises if absent, incomplete, ineligible, or a T7 path.
    config = load_runtime_config(runtime_config_path, cohort, channels)

    versions = software_versions()
    commit = _git_state()
    code_fp = code_fingerprint(CODE_FILES)

    rows, diag_rows, restored = [], [], 0
    for ch, mat_t, cell_ids, props in captures:
        entry = config[ch]
        seed = capture_seed(base_seed, cohort, ch)
        fp = build_fingerprint(
            checkpoint_schema_version=CHECKPOINT_SCHEMA_VERSION,
            call_status_schema=list(CALL_STATUSES),
            cohort=cohort,
            channel_id_raw=ch,
            n_cells=int(len(cell_ids)),
            cell_id_checksum=ordered_cell_checksum(cell_ids),
            matrix_fingerprint=matrix_fingerprint(mat_t),
            config_row_checksum=entry['config_row_checksum'],
            expected_doublet_rate=float(entry['expected_doublet_rate']),
            donor_proportions_checksum=donor_proportions_checksum(props),
            donor_proportion_stage=entry['donor_proportion_stage'],
            platform=entry['platform'],
            base_seed=int(base_seed),
            capture_seed=int(seed),
            scrublet_params=SCRUBLET_PARAMS,
            software_versions=versions,
            git_commit=commit,
            code_fingerprint=code_fp,
        )
        key = capture_key(cohort, ch)

        stored = load_checkpoint(ckpt_dir, key, fp, expected_cell_ids=cell_ids)
        if stored is not None:
            restored += 1
            scores = stored['scores']
            call_status = stored['call_status']
            candidate = stored['candidate']
            capture_status = stored['capture_status']
            diagnostics = stored['diagnostics']
        else:
            res = score_one_capture(
                mat_t, cell_ids,
                expected_doublet_rate=entry['expected_doublet_rate'], seed=seed,
                donor_proportions=props, platform=entry['platform'], scrublet_fn=scrublet_fn)
            scores, call_status = res.scores, res.call_status
            candidate, capture_status = res.candidate, res.capture_status
            diagnostics = res.diagnostics
            save_checkpoint(ckpt_dir, key, fp, scores=scores, call_status=call_status,
                            candidate=candidate, cell_ids=cell_ids,
                            capture_status=capture_status, diagnostics=diagnostics)

        for cid, sc, st, cand in zip(cell_ids, scores, call_status, candidate):
            # AUTHORITATIVE IDENTITY is the three columns (cohort, channel_id,
            # cell_id). gate2_cell_key is a DERIVED convenience column only.
            # It is a canonical JSON array rather than a '|' join, because an
            # unescaped join is not injective: a channel containing '|' would
            # collide with a different (channel, cell) pair.
            rows.append({'gate2_cell_key': derived_cell_key(cohort, ch, cid),
                         'cohort': cohort, 'channel_id': ch, 'cell_id': cid,
                         'scrublet_score': sc, 'call_status': st,
                         # written as the tri-state it is; '' is unknown, so no
                         # consumer can read a blank as False
                         # nullable integer, NOT True/False: a blank is UNKNOWN and
                         # must never be coerced to False. See docs/data-contracts.md.
                         'candidate_doublet': {1: '1', 0: '0'}.get(int(cand), '')})
        # Diagnostics are rebuilt for EVERY capture on every run, restored ones
        # included -- the old split runner emitted only newly processed channels,
        # so a resumed run's summary silently covered part of the cohort.
        diag_rows.append({'cohort': cohort, 'channel_id': ch, 'capture_status': capture_status,
                          'capture_seed': seed, 'restored_from_checkpoint': stored is not None,
                          **{k: (json.dumps(v) if isinstance(v, (list, dict)) else v)
                             for k, v in diagnostics.items()}})

    if len(diag_rows) != len(captures):
        raise RuntimeError(f'diagnostics rows {len(diag_rows)} != captures {len(captures)}')

    cells_out = os.path.join(out_dir, f'{cohort}_gate2_calls.csv')
    with open(cells_out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator='\n')
        w.writeheader()
        w.writerows(rows)

    diag_out = os.path.join(out_dir, f'{cohort}_gate2_diagnostics.csv')
    fields = sorted({k for r in diag_rows for k in r})
    with open(diag_out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, lineterminator='\n', restval='')
        w.writeheader()
        w.writerows(diag_rows)

    n_candidates = sum(mechanical_candidate_mask(r['call_status'] for r in rows))
    print(f'{cohort}: {len(captures)} captures ({restored} restored), {len(rows)} cells')
    print(f'  mechanical doublet candidates: {n_candidates}')
    print('  NOTHING DELETED -- Gate 2 is a scoring pass.')
    print('  No cell has final deletion authorization.')
    return {'captures': len(captures), 'restored': restored, 'cells': len(rows),
            'mechanical_doublet_candidates': n_candidates,
            'calls_csv': cells_out, 'diagnostics_csv': diag_out}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cohort', required=True)
    ap.add_argument('--layout', required=True, choices=['combined', 'split'])
    ap.add_argument('--input', required=True)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--runtime-config', required=True,
                    help='authorized Gate 2 runtime rate config; no default, no fallback')
    ap.add_argument('--base-seed', required=True, type=int,
                    help='explicit base seed; per-capture seeds are derived from it')
    a = ap.parse_args()
    try:
        run(a.cohort, a.layout, a.input, a.out_dir, a.runtime_config, a.base_seed)
    except Gate2ConfigError as exc:
        print(f'GATE 2 REFUSED TO START: {exc}', file=sys.stderr)
        raise SystemExit(2)
    except CellIdentityError as exc:
        print(f'GATE 2 REFUSED TO START: {exc}', file=sys.stderr)
        raise SystemExit(4)
    except CheckpointMismatch as exc:
        print(f'GATE 2 REFUSED TO RESUME: {exc}', file=sys.stderr)
        raise SystemExit(3)


if __name__ == '__main__':
    main()
