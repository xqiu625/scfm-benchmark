"""Negative tests for the capture-provenance bundle.

Each test exists because the v1 bundle got that exact thing wrong. They are
written to FAIL if the defect returns, not to confirm the current output.

Importing the builder re-runs it, which regenerates the derived CSVs. That is
deliberate: a mutation to the builder is then visible in the artifacts these
tests read.
"""
import csv
import os
import re
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import build_capture_provenance as B  # noqa: E402


def load(name):
    with open(os.path.join(HERE, name)) as fh:
        return list(csv.DictReader(fh))


LONG = load('capture_provenance_long.csv')
SUMMARY = load('cohort_summary.csv')
RECON = load('donor_coverage_reconciliation.csv')
BOUNDARY = load('outcome_column_boundary.csv')


# --- 1. coverage-only donors must never be membership rows ------------------

def test_a_donor_without_membership_is_not_in_the_membership_table():
    """v1 appended 15 coverage placeholders to the long table, inflating
    n_membership_rows, n_biological_samples and the status distribution with rows
    that record no membership at all."""
    total_missing = 0
    for rec in RECON:
        cohort = rec['cohort']
        missing = [d for d in rec['canonical_without_membership'].split(';') if d]
        total_missing += len(missing)
        have = {B.norm(cohort, r['donor_key']) for r in LONG if r['cohort'] == cohort}
        for d in missing:
            assert B.norm(cohort, d) not in have, f'{cohort}/{d} leaked into the table'
    assert total_missing > 0, 'no coverage gap left to test against'


def test_no_membership_row_has_the_unresolved_membership_kind():
    """`membership_kind='unresolved'` was the marker the placeholders carried."""
    assert [r for r in LONG if r['membership_kind'] == 'unresolved'] == []


def test_every_membership_row_names_an_analysis_channel():
    for r in LONG:
        assert r['current_analysis_channel_id'].strip(), (
            f'{r["cohort"]}/{r["biological_sample_key"]} has no analysis channel, so '
            'it is not a membership')


# --- 2. an empty string is not a category ----------------------------------

def test_summary_counts_exclude_empty_keys():
    """v1 counted the empty channel of a placeholder row as one channel:
    Chongqing 133->134, Cillo 44->45, Wilk 23->24."""
    for s in SUMMARY:
        sub = [r for r in LONG if r['cohort'] == s['cohort']]
        for field, col in (('n_analysis_channels', 'current_analysis_channel_id'),
                           ('n_physical_captures', 'physical_capture_key'),
                           ('n_donors', 'donor_key'),
                           ('n_biological_samples', 'biological_sample_key')):
            with_empty = len({r[col] for r in sub})
            without = len({r[col] for r in sub if r[col].strip()})
            assert int(s[field]) == without, f'{s["cohort"]}.{field} counts an empty key'
            if with_empty != without:      # this cohort actually has empties
                assert int(s[field]) != with_empty


def test_stephenson_has_no_physical_captures_rather_than_one_empty_one():
    st = [s for s in SUMMARY if s['cohort'] == 'Stephenson 2021'][0]
    assert int(st['n_physical_captures']) == 0


# --- 3. a repository gloss is not `resolved` -------------------------------

def test_a_key_read_by_repository_annotation_is_not_resolved():
    """SCHEMA.md defines `resolved` as the DEPOSIT identifying the capture unit.
    Liu's lane token is read as a 10x chip well by this repository's annotation,
    which is not that."""
    gloss = [r for r in LONG
             if "annotation, not by the deposit" in r['ambiguity_reason']]
    assert gloss, 'expected at least one annotation-derived capture key'
    for r in gloss:
        assert r['physical_capture_status'] != 'resolved', (
            f'{r["cohort"]} claims deposit-resolved on a repository gloss')


def test_a_builder_derived_key_is_not_resolved():
    derived = [r for r in LONG if 'derived by' in r['ambiguity_reason']]
    assert derived, 'expected at least one builder-derived capture key'
    for r in derived:
        assert r['physical_capture_status'] != 'resolved'


def test_resolved_requires_a_capture_key():
    for r in LONG:
        if r['physical_capture_status'] == 'resolved':
            assert r['physical_capture_key'].strip()


# --- 4. an input holding outcome must not be called unread -----------------

KNOWN_OUTCOME_INPUTS = ('stephenson_patient_outcomes.csv', 'canonical_donor_labels.csv',
                        'ren_gse158055_channel_inventory.csv',
                        'qiu_gse167363_channel_inventory.csv')


def test_every_input_that_holds_outcome_records_it_as_excluded():
    """v1 claimed "no outcome column was read" while csv.DictReader materialised
    Outcome, event, outcome and disease_state."""
    by_name = {os.path.basename(b['input']): b for b in BOUNDARY}
    for name in KNOWN_OUTCOME_INPUTS:
        assert name in by_name, f'{name} is no longer declared as an input'
        assert by_name[name]['outcome_columns_excluded'].strip(), (
            f'{name} holds an outcome column but records none as excluded')


def test_no_input_reads_an_outcome_column():
    for b in BOUNDARY:
        for col in b['columns_read'].split(';'):
            assert not B.OUTCOME_PATTERN.search(col), (
                f'{b["input"]} reads outcome-like column {col!r}')


def test_the_parser_itself_receives_usecols():
    """Capture the kwargs the parser ACTUALLY received.

    Asserting on the returned columns proves nothing here: the data read is
    followed by `df[list(usecols)]`, so removing `usecols=` still yields only the
    requested keys while pandas has parsed every column. An earlier version of
    this test checked the returned keys and was reported as killing that mutation
    -- it did not; the mutation had also removed the projection, so the kill never
    isolated the parser argument.
    """
    path = os.path.join(B.T3, 'stephenson_patient_outcomes.csv')
    wanted = ['patient_id', 'sample_ids']
    calls = []

    def spy(p, **kw):
        calls.append(kw)
        return B.pd.read_csv(p, **kw)

    recs = B.read_columns(path, wanted, reader=spy)
    assert recs, 'no rows returned'
    header = [c for c in calls if c.get('nrows') == 0]
    data = [c for c in calls if c.get('nrows') != 0]
    assert len(header) == 1, f'expected exactly one header read, got {len(header)}'
    assert len(data) == 1, f'expected exactly one data read, got {len(data)}'
    assert data[0].get('usecols') == wanted, (
        f'the data read received usecols={data[0].get("usecols")!r}; without it the '
        'parser reads every column, outcome included, and the later projection '
        'hides that')


def test_read_columns_returns_only_the_requested_columns():
    """Weaker than the kwargs test above and kept only as a schema check on the
    returned records -- it cannot detect a dropped `usecols=`."""
    path = os.path.join(B.T3, 'stephenson_patient_outcomes.csv')
    wanted = ['patient_id', 'sample_ids']
    recs = B.read_columns(path, wanted)
    assert recs, 'no rows returned'
    for rec in recs[:5]:
        assert list(rec.keys()) == wanted
    present = list(B.pd.read_csv(path, nrows=0).columns)
    assert 'Outcome' in present, 'input no longer holds the column under test'


# --- COMBAT membership is observed, not intended ---------------------------

def test_combat_membership_is_not_claimed_as_intended_loading():
    """The pairs are donors observed among a channel's retained cells after the
    depositors' QC. No pool-design artifact lists the intended donors per channel,
    so absence of a pair is not evidence the donor was not loaded."""
    combat = [r for r in LONG if r['cohort'] == 'COMBAT 2022']
    assert combat
    for r in combat:
        assert r['membership_kind'] == 'observed_assigned', r['membership_kind']
        assert r['membership_kind'] != 'intended'


def test_every_observed_membership_states_its_limitation():
    obs = [r for r in LONG if r['membership_kind'].startswith('observed')]
    assert obs, 'expected observed memberships'
    for r in obs:
        assert r['ambiguity_reason'].strip(), (
            f'{r["cohort"]} records an observed membership with no stated limitation')


# --- cited evidence is frozen too -----------------------------------------

def test_a_cited_but_unparsed_evidence_file_is_still_checksummed():
    """The Stephenson rows rest on the T7 memo's SDRF reading. A cited file with no
    checksum can change without the bundle noticing."""
    with open(os.path.join(HERE, 'INPUT_CHECKSUMS.txt')) as fh:
        text = fh.read()
    assert 'T7_doublet_prior_evidence/memo.md' in text
    cited = [r for r in LONG if 'memo.md' in r['evidence_locator']]
    assert cited, 'expected rows citing the memo'


def test_requesting_an_outcome_column_is_refused():
    path = os.path.join(B.T3, 'stephenson_patient_outcomes.csv')
    with pytest.raises(SystemExit, match='refusing to read outcome-like'):
        B.read_columns(path, ['patient_id', 'Outcome'])


def test_the_outcome_pattern_covers_the_columns_actually_present():
    """Not a restatement of the constant: these are the real column names in the
    four inputs, so a narrowed pattern fails here."""
    for name in ('Outcome', 'outcome_unique_within_patient', 'event',
                 'event_definition', 'disease', 'outcome', 'disease_state'):
        assert B.OUTCOME_PATTERN.search(name), f'{name} would slip past the guard'


# --- 5. summaries must be rebuildable from the long table ------------------

def test_every_summary_number_is_rebuildable_from_the_long_table():
    for s in SUMMARY:
        sub = [r for r in LONG if r['cohort'] == s['cohort']]
        assert int(s['n_membership_rows']) == len(sub)
        assert int(s['n_donors']) == len({r['donor_key'] for r in sub if r['donor_key']})


def test_no_duplicate_membership_keys():
    seen = set()
    for r in LONG:
        k = (r['cohort'], r['donor_key'], r['biological_sample_key'],
             r['physical_capture_key'], r['current_analysis_channel_id'])
        assert k not in seen, f'duplicate membership {k}'
        seen.add(k)


# --- transactional publication and raw identity ----------------------------

def _load_extractor():
    """Import the extraction module without running it. Its HPCC paths are only
    touched by the extract functions, which these tests do not call."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        'extract_hpcc_membership', os.path.join(HERE, 'extract_hpcc_membership.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


X = _load_extractor()


def test_cheap_checksum_matches_the_repository_implementation():
    """Byte-for-byte port, verified against the original rather than assumed. An
    identity recorded by one algorithm and checked by another is not a check."""
    import importlib.util
    src = os.path.join(B.REPO, 'scripts/32_rebuild_upfront_qc/code/gate2_source_crosswalk.py')
    text = open(src).read()
    start = text.index('def cheap_checksum')
    end = text.index('\ndef ', start)
    ns = {'os': os, 'hashlib': __import__('hashlib'), 'pd': B.pd}
    exec(text[start:end], ns)                      # noqa: S102 - the original, run as-is
    probe = os.path.join(HERE, 'capture_provenance_long.csv')
    assert X.cheap_checksum(probe) == ns['cheap_checksum'](probe)
    assert importlib.util  # keep the import used


def _identity_probe(tmp_path, head, tail, size):
    p = tmp_path / 'probe.bin'
    body = head + b'\x00' * (size - len(head) - len(tail)) + tail
    p.write_bytes(body)
    return str(p)


def test_a_same_size_file_with_changed_head_fails_the_identity_check():
    """Comparing only the size would pass this. v1 compared only the size."""
    import tempfile, pathlib  # noqa: E401
    with tempfile.TemporaryDirectory() as d:
        tmp = pathlib.Path(d)
        size = (1 << 20) + 4096
        good = _identity_probe(tmp, b'HEAD-A' + b'x' * 100, b'TAIL' + b'y' * 100, size)
        trusted = X.cheap_checksum(good)
        os.replace(good, str(tmp / 'old'))
        bad = _identity_probe(tmp, b'HEAD-B' + b'x' * 100, b'TAIL' + b'y' * 100, size)
        assert os.path.getsize(bad) == size
        assert X.cheap_checksum(bad) != trusted, 'head change escaped the identity'
        with pytest.raises(SystemExit, match='does not match the trusted'):
            X.assert_raw_identity(bad, trusted)


def test_a_same_size_file_with_changed_tail_fails_the_identity_check():
    import tempfile, pathlib  # noqa: E401
    with tempfile.TemporaryDirectory() as d:
        tmp = pathlib.Path(d)
        size = (1 << 20) + 4096
        good = _identity_probe(tmp, b'HEAD' + b'x' * 100, b'TAIL-A' + b'y' * 100, size)
        trusted = X.cheap_checksum(good)
        os.replace(good, str(tmp / 'old'))
        bad = _identity_probe(tmp, b'HEAD' + b'x' * 100, b'TAIL-B' + b'y' * 100, size)
        assert os.path.getsize(bad) == size
        assert X.cheap_checksum(bad) != trusted, 'tail change escaped the identity'
        with pytest.raises(SystemExit, match='does not match the trusted'):
            X.assert_raw_identity(bad, trusted)


PUBLISHED = ('combat_donor_by_channel.csv', 'liu_donor_by_lane.csv',
             'RAW_INPUT_MANIFEST.txt')


def _seed(dest):
    """A miniature stand-in bundle, so a fault-injection test never touches the
    real extracts. Driving publish() into a failure with a guard mutated away
    would otherwise leave a genuinely empty extract behind, and the next reader
    fails for a reason that looks unrelated."""
    before = {}
    for i, n in enumerate(PUBLISHED):
        before[n] = f'original-{i}\n'.encode()
        with open(os.path.join(dest, n), 'wb') as fh:
            fh.write(before[n])
    return before


def _read_all(dest, names):
    return {n: open(os.path.join(dest, n), 'rb').read() for n in names}


def test_a_failure_part_way_through_leaves_every_previous_version_intact():
    """v2 overwrote the COMBAT csv, then the Liu csv, then the manifest, so a
    failure between any two left a new extract beside an old one."""
    import tempfile
    with tempfile.TemporaryDirectory() as dest:
        before = _seed(dest)
        real_replace = os.replace
        calls = {'n': 0}

        def flaky(src, dst):
            calls['n'] += 1
            if calls['n'] == 2:
                raise OSError('simulated failure publishing the second output')
            return real_replace(src, dst)

        os.replace = flaky
        try:
            with pytest.raises(OSError, match='simulated failure'):
                X.publish({n: before[n] + b'# new\n' for n in PUBLISHED}, dest=dest)
        finally:
            os.replace = real_replace
        assert not [f for f in os.listdir(dest) if f.endswith(X.INPROGRESS)], \
            'a .inprogress temp survived the failure'


def test_an_empty_payload_is_refused_before_anything_is_replaced():
    import tempfile
    with tempfile.TemporaryDirectory() as dest:
        before = _seed(dest)
        with pytest.raises(SystemExit, match='empty payload'):
            X.publish({PUBLISHED[0]: b'', PUBLISHED[1]: b'x', PUBLISHED[2]: b'y'},
                      dest=dest)
        assert _read_all(dest, PUBLISHED) == before, \
            'an output was replaced despite the refusal'
        assert not [f for f in os.listdir(dest) if f.endswith(X.INPROGRESS)]


def test_publish_replaces_every_output_and_leaves_no_temp_on_success():
    import tempfile
    with tempfile.TemporaryDirectory() as dest:
        _seed(dest)
        new = {n: f'new-{n}\n'.encode() for n in PUBLISHED}
        X.publish(new, dest=dest)
        assert _read_all(dest, PUBLISHED) == new
        assert not [f for f in os.listdir(dest) if f.endswith(X.INPROGRESS)]


def test_the_real_extracts_are_never_written_by_these_tests():
    """A guard on the guards: the fault-injection tests above must stay inside a
    temporary directory."""
    src = open(os.path.join(HERE, 'test_capture_provenance.py')).read()
    import ast as _ast
    tree = _ast.parse(src)
    for node in _ast.walk(tree):
        if (isinstance(node, _ast.Call) and getattr(node.func, 'attr', '') == 'publish'):
            kw = {k.arg for k in node.keywords}
            assert 'dest' in kw, 'a publish() call in these tests omits dest='
