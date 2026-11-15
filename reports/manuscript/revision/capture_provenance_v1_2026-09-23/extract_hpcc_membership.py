"""READ-ONLY extraction: HPCC raw sources -> the two membership extracts.

Run on HPCC from the repository root:

    python reports/.../extract_hpcc_membership.py

Publishes, beside itself, ALL-OR-NOTHING:
    combat_donor_by_channel.csv
    liu_donor_by_lane.csv
    RAW_INPUT_MANIFEST.txt

TRANSACTIONAL PUBLICATION
-------------------------
A previous version overwrote the COMBAT csv, then the Liu csv, then the manifest.
A failure between any two left a new extract beside an old one, or two new
extracts beside a stale manifest, or a truncated csv -- states that defeat the
point of a provenance bundle. Nothing is extracted straight to its final path
now: every output is built in full, written to a `.inprogress` sibling, verified,
and only then `os.replace`d. A failure leaves all three previous versions intact
and removes every temporary file.

RAW IDENTITY
------------
The COMBAT h5ad is 6.4 GB. Rather than re-hash it whole, the identity is the
`<size>:<sha256-of(size + first 1 MiB + last 1 MiB)>[:16]` form an earlier audit
already recorded for this exact path. `cheap_checksum()` below is a byte-for-byte
port of `gate2_source_crosswalk.py:62`, and the WHOLE identity is recomputed and
required to match -- not just the size. Comparing only the size would pass a
same-size file with replaced content.

OUTCOME BOUNDARY
----------------
COMBAT: `h5py` opens NAMED DATASETS, so only `obs/Channel_ID`, `obs/COMBAT_ID`
and the legacy `obs/__categories` entries for those two are materialised. The
deposited obs also holds clinical fields; none is opened.

Liu: the `.best` files carry genotype calls and no outcome field.

Note the boundary honestly: `cheap_checksum` and the Liu hashing read raw BYTES,
which include whatever bytes an outcome field occupies. What is guaranteed is
that no outcome column is parsed into a field, inspected, propagated or exported.

Nothing is written outside this directory. No cohort matrix is loaded, no cell is
deleted, no job is submitted, no parameter or rate is selected.
"""
import collections
import csv
import glob
import hashlib
import io
import os

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = '/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw'
COMBAT_H5AD = os.path.join(RAW, 'COMBAT/COMBAT-CITESeq-DATA.h5ad')
LIU_BEST_GLOB = os.path.join(RAW, 'GSE161918_liu/demuxlet/*_demuxed.best')

#: Recorded by an earlier audit for this exact path:
#:   T7_doublet_prior_evidence/gate2_source_crosswalk.csv :: column
#:   source_checksum, rows where cohort == 'combat'
COMBAT_TRUSTED_IDENTITY = '6409089483:d4435c68dd07c170'

COMBAT_OBS_CHANNEL = 'Channel_ID'
COMBAT_OBS_DONOR = 'COMBAT_ID'
INPROGRESS = '.inprogress'


def cheap_checksum(path, nbytes=1 << 20):
    """Byte-for-byte port of gate2_source_crosswalk.py:62.

    `<size>:<first 16 hex of sha256(str(size) + first nbytes + last nbytes)>`.
    Kept identical on purpose: an identity recorded by one script and verified by
    a different algorithm in another is not a verification.
    """
    sz = os.path.getsize(path)
    h = hashlib.sha256(str(sz).encode())
    with open(path, 'rb') as fh:
        h.update(fh.read(nbytes))
        if sz > nbytes:
            fh.seek(max(0, sz - nbytes))
            h.update(fh.read(nbytes))
    return f'{sz}:{h.hexdigest()[:16]}'


def assert_raw_identity(path, trusted):
    """Recompute the WHOLE identity and require equality."""
    actual = cheap_checksum(path)
    if actual != trusted:
        raise SystemExit(
            f'{path}: identity {actual} does not match the trusted {trusted}. '
            'The reused checksum does not describe this file.')
    return actual


def read_obs_categorical(obs, name):
    """One obs column, handling both anndata encodings. Only the named datasets
    are opened: modern stores `codes`/`categories` in a Group, legacy (which this
    deposit uses) stores codes as a plain int dataset with labels in
    `obs/__categories/<name>`."""
    node = obs[name]
    if isinstance(node, h5py.Group):
        codes, cats = node['codes'][:], node['categories'][:]
    else:
        codes, cats = node[:], obs['__categories'][name][:]
    cats = np.array([c.decode() if isinstance(c, bytes) else str(c) for c in cats])
    return cats[codes]


def csv_bytes(fieldnames, records):
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=fieldnames, lineterminator='\r\n')
    w.writeheader()
    w.writerows(records)
    return buf.getvalue().encode()


def extract_combat():
    """Returns (records, manifest_line). Writes nothing."""
    identity = assert_raw_identity(COMBAT_H5AD, COMBAT_TRUSTED_IDENTITY)
    with h5py.File(COMBAT_H5AD, 'r') as f:
        obs = f['obs']
        for req in (COMBAT_OBS_CHANNEL, COMBAT_OBS_DONOR):
            if req not in obs:
                raise SystemExit(f'COMBAT obs has no {req!r}')
        channel = read_obs_categorical(obs, COMBAT_OBS_CHANNEL)
        donor = read_obs_categorical(obs, COMBAT_OBS_DONOR)
    pairs = collections.Counter(zip(channel, donor))
    records = [{'channel_id': c, 'donor_id': d, 'n_cells': n}
               for (c, d), n in sorted(pairs.items())]
    print(f'COMBAT: {len(channel)} cells, {len(set(channel))} channels, '
          f'{len(set(donor))} donors, {len(records)} membership pairs')
    return records, (
        f'{identity}  {COMBAT_H5AD}  (whole identity recomputed via cheap_checksum '
        f'and matched; datasets read: obs/{COMBAT_OBS_CHANNEL}, obs/{COMBAT_OBS_DONOR})')


def extract_liu():
    """Returns (records, manifest_lines). Writes nothing."""
    paths = sorted(glob.glob(LIU_BEST_GLOB))
    if not paths:
        raise SystemExit(f'no .best files at {LIU_BEST_GLOB}')
    records, lines = [], []
    for p in paths:
        lane = os.path.basename(p).replace('_demuxed.best', '')
        h = hashlib.sha256()
        with open(p, 'rb') as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b''):
                h.update(chunk)
        lines.append(f'{h.hexdigest()}  {p}')
        with open(p) as fh:
            hdr = fh.readline().rstrip('\n').split('\t')
            i_type, i_sng = hdr.index('DROPLET.TYPE'), hdr.index('SNG.BEST.GUESS')
            c = collections.Counter()
            for line in fh:
                r = line.rstrip('\n').split('\t')
                if r[i_type] == 'SNG':
                    c[r[i_sng]] += 1
        for donor, n in sorted(c.items()):
            records.append({'lane': lane, 'donor': donor, 'n_sng_cells': n})
    print(f'Liu: {len(paths)} lanes, {len({r["donor"] for r in records})} SNG donors, '
          f'{len(records)} membership pairs')
    return records, lines


def publish(outputs, dest=None):
    """Write every output to a temp sibling, verify, then replace all of them.

    `outputs` maps final filename -> bytes. Nothing is replaced until every temp
    exists and re-reads byte-identical, so a failure part-way leaves all previous
    versions in place.

    `dest` defaults to this directory and is injectable so the fault-injection
    tests can exercise the failure paths in a temporary directory. They must not
    run against the real bundle: a test that drives publish() into a failure while
    a guard is mutated away would otherwise leave a genuinely empty extract on
    disk, and the next thing to read it fails for an unrelated-looking reason.
    """
    dest = HERE if dest is None else dest
    temps = {name: os.path.join(dest, name + INPROGRESS) for name in outputs}
    try:
        for name, payload in outputs.items():
            if not payload:
                raise SystemExit(f'{name}: refusing to publish an empty payload')
            with open(temps[name], 'wb') as fh:
                fh.write(payload)
                fh.flush()
                os.fsync(fh.fileno())
            with open(temps[name], 'rb') as fh:
                if fh.read() != payload:
                    raise SystemExit(f'{name}: temp file does not match its payload')
        for name in outputs:                       # only now is anything visible
            os.replace(temps[name], os.path.join(dest, name))
    finally:
        for t in temps.values():
            if os.path.exists(t):
                os.remove(t)


if __name__ == '__main__':
    combat_records, combat_line = extract_combat()
    liu_records, liu_lines = extract_liu()

    manifest = ['# Raw HPCC sources read by extract_hpcc_membership.py.',
                '# sha256 unless the line marks a recomputed <size>:<partial> identity.',
                combat_line] + liu_lines
    publish({
        'combat_donor_by_channel.csv': csv_bytes(
            ['channel_id', 'donor_id', 'n_cells'], combat_records),
        'liu_donor_by_lane.csv': csv_bytes(
            ['lane', 'donor', 'n_sng_cells'], liu_records),
        'RAW_INPUT_MANIFEST.txt': ('\n'.join(manifest) + '\n').encode(),
    })
    print(f'published 3 outputs atomically; raw inputs recorded: {len(manifest) - 2}')
