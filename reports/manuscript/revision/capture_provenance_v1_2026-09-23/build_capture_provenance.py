"""READ-ONLY: build the biological-unit x physical-capture membership long table.

Writes only into this directory. Modifies no builder, no DecontX output, no Gate 2
use class, no existing crosswalk. Submits nothing. Selects no parameter or rate.

Row = ONE MEMBERSHIP of a biological unit in a physical capture. A donor with no
membership in any capture is NOT a row here -- it belongs to
`donor_coverage_reconciliation.csv`. A previous version appended 15 such
coverage-only placeholders to this table, which inflated `n_membership_rows`,
`n_biological_samples` and the capture-status distribution with rows recording no
membership at all.

OUTCOME BOUNDARY -- STRICT, AND WHY IT HAD TO CHANGE
----------------------------------------------------
A previous version read every input with `csv.DictReader`, which materialises the
whole row. Four inputs carry outcome columns:

    stephenson_patient_outcomes.csv   Outcome, outcome_unique_within_patient
    canonical_donor_labels.csv        event, event_definition, disease
    ren_gse158055_channel_inventory   outcome
    qiu_gse167363_channel_inventory   disease_state

so "no outcome column was read" was FALSE: the values were parsed and held in
memory, and merely never indexed. Not referencing `r['Outcome']` is not the same
as not reading it -- the same defect this repository already recorded for
`pd.read_csv` followed by a column drop.

Every input is now read with `pandas.read_csv(usecols=...)`, which excludes the
other columns AT PARSE TIME. `read_columns()` records per input which columns were
requested and which present columns were excluded, and refuses any request for a
column matching the outcome vocabulary.

The claim that follows, and nothing stronger: **no outcome column is parsed into a
field, inspected, propagated, exported or analysed.** Checksumming opens each
input with `open(path, 'rb').read()`, so the raw BYTES of the whole file -- an
outcome field's bytes included -- do pass through a hash. Hashing bytes and
reading a value are different acts and are described separately here.
"""
import csv
import hashlib
import os
import re
import sys
from collections import Counter, defaultdict

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))
LIT = os.path.join(REPO, 'reports/manuscript/revision/literature_audit')
INV = os.path.join(REPO, 'scripts/32_rebuild_upfront_qc/channel_inventories')
T3 = os.path.join(LIT, 'T3_patient_sample_mapping')

FIELDS = ['cohort', 'accession', 'donor_key', 'biological_sample_key', 'draw_key',
          'physical_capture_key', 'current_analysis_channel_id', 'membership_kind',
          'physical_capture_status', 'analysis_matches_physical', 'evidence_source',
          'evidence_locator', 'evidence_grade', 'ambiguity_reason']

#: Any column whose name matches is an outcome/label field and must never be
#: requested. Deliberately broader than the four known names, so a new input with
#: a differently spelled label column is refused rather than silently read.
OUTCOME_PATTERN = re.compile(
    r'outcome|event|death|died|surviv|mortal|deceased|discharg|fatal|severity|disease',
    re.I)

INPUTS = {}
COLUMN_BOUNDARY = []


def read_columns(path, usecols, reader=None):
    """Read ONLY `usecols`, and record what was excluded at parse time.

    `reader` defaults to `pandas.read_csv` and exists so a test can capture the
    kwargs the parser actually received. Asserting on the RETURNED columns cannot
    do that: the data read is followed by a projection, so dropping `usecols=`
    still yields only the requested keys while the parser has read every column.
    """
    reader = pd.read_csv if reader is None else reader
    present = list(reader(path, nrows=0).columns)
    missing = [c for c in usecols if c not in present]
    if missing:
        raise SystemExit(f'{path}: requested columns absent: {missing}')
    forbidden = [c for c in usecols if OUTCOME_PATTERN.search(c)]
    if forbidden:
        raise SystemExit(
            f'{path}: refusing to read outcome-like column(s) {forbidden}. This '
            'table records capture provenance; a label column has no use here.')
    excluded = [c for c in present if c not in usecols]
    with open(path, 'rb') as fh:
        INPUTS[os.path.relpath(path, REPO)] = hashlib.sha256(fh.read()).hexdigest()
    COLUMN_BOUNDARY.append({
        'input': os.path.relpath(path, REPO),
        'n_columns_present': len(present),
        'columns_read': ';'.join(usecols),
        'columns_excluded_at_parse_time': ';'.join(excluded),
        'outcome_columns_excluded': ';'.join(c for c in excluded
                                            if OUTCOME_PATTERN.search(c)),
    })
    df = reader(path, usecols=usecols, dtype=str, keep_default_na=False)
    return df[list(usecols)].to_dict('records')


def row(**kw):
    unknown = set(kw) - set(FIELDS)
    if unknown:
        raise SystemExit(f'undeclared column(s): {sorted(unknown)}')
    r = {k: '' for k in FIELDS}
    r.update(kw)
    return r


rows = []

# --- Qiu -------------------------------------------------------------------
for r in read_columns(os.path.join(INV, 'qiu_gse167363_channel_inventory.csv'),
                      ['gsm', 'sample_name', 'donor', 'draw_time']):
    rows.append(row(
        cohort='Qiu 2021', accession='GSE167363',
        donor_key=r['donor'], biological_sample_key=r['sample_name'],
        draw_key=r['draw_time'], physical_capture_key=r['gsm'],
        current_analysis_channel_id=r['gsm'],
        membership_kind='intended', physical_capture_status='resolved',
        analysis_matches_physical='true',
        evidence_source='qiu_gse167363_channel_inventory.csv',
        evidence_locator='channel_inventories/qiu_gse167363_channel_inventory.csv :: '
                         'one row per gsm; columns donor, sample_name, draw_time',
        evidence_grade='measured_metadata'))

# --- Chongqing -------------------------------------------------------------
for r in read_columns(os.path.join(INV, 'chongqing_build_summary.csv'),
                      ['gsm', 'donor_id']):
    rows.append(row(
        cohort='Chongqing sepsis 2025', accession='GSE279452',
        donor_key=r['donor_id'], biological_sample_key=r['gsm'],
        draw_key='single_draw_no_variation', physical_capture_key=r['gsm'],
        current_analysis_channel_id=r['gsm'],
        membership_kind='intended', physical_capture_status='resolved',
        analysis_matches_physical='true',
        evidence_source='chongqing_build_summary.csv',
        evidence_locator='channel_inventories/chongqing_build_summary.csv :: '
                         'one row per gsm; column donor_id',
        evidence_grade='measured_metadata'))

# --- Ren -------------------------------------------------------------------
# 284 libraries across 196 donors: a donor contributes several libraries
# (different sample types and draws), so donor != capture here.
for r in read_columns(os.path.join(INV, 'ren_gse158055_channel_inventory.csv'),
                      ['gsm', 'sampleID', 'donor', 'sample_time']):
    rows.append(row(
        cohort='Ren 2021', accession='GSE158055',
        donor_key=r['donor'], biological_sample_key=r['sampleID'],
        draw_key=r['sample_time'], physical_capture_key=r['sampleID'],
        current_analysis_channel_id=r['sampleID'],
        membership_kind='intended', physical_capture_status='resolved',
        analysis_matches_physical='true',
        evidence_source='ren_gse158055_channel_inventory.csv',
        evidence_locator='channel_inventories/ren_gse158055_channel_inventory.csv :: '
                         'one row per gsm; columns sampleID, donor, sample_time',
        evidence_grade='measured_metadata'))

# --- Wilk ------------------------------------------------------------------
# The benchmark cohort is GSE174072 (Seq-Well). The T3 mapping audit covers
# GSE150728, a DIFFERENT deposit with different donor ids -- see FINDINGS.md.
# channel_id here is builder-derived, which is why it is never `resolved`.
for r in read_columns(os.path.join(INV, 'wilk_build_summary.csv'),
                      ['channel_id', 'donor_id']):
    rows.append(row(
        cohort='Wilk 2021 JEM', accession='GSE174072',
        donor_key=r['donor_id'], biological_sample_key=r['channel_id'],
        draw_key='', physical_capture_key=r['channel_id'],
        current_analysis_channel_id=r['channel_id'],
        membership_kind='intended', physical_capture_status='provisional',
        analysis_matches_physical='unknown',
        evidence_source='wilk_build_summary.csv',
        evidence_locator='channel_inventories/wilk_build_summary.csv :: '
                         'one row per channel_id; column donor_id',
        evidence_grade='derived',
        ambiguity_reason='Seq-Well array, not a 10x GEM well; the key is derived by '
                         'the builder rather than stated by the deposit, and the '
                         'deposit does not state that one channel_id is one array'))

# --- Cillo -----------------------------------------------------------------
# Hashtag-multiplexed: `members` is "donor@day" pairs separated by ";".
# GSM5465163 is `exact member-set match` with an EMPTY member set
# (n_donors_in_lib=0), which is why 40 exact-match libraries yield 39 captures.
for r in read_columns(os.path.join(INV, 'cillo_gse180578_channel_inventory.csv'),
                      ['gsm', 'members', 'h5ad_suffix', 'join_status']):
    members = [m for m in r['members'].split(';') if m.strip()]
    suffix_ok = bool(r['h5ad_suffix'].strip())
    status = 'resolved' if suffix_ok else 'unresolved'
    reason = '' if suffix_ok else (
        f"no h5ad_suffix ({r['join_status']}); the library's cells cannot be "
        'attributed to a single GEM well through the h5ad join key')
    for m in members:
        donor, _, day = m.partition('@')
        rows.append(row(
            cohort='Cillo 2021', accession='GSE180578',
            donor_key=donor, biological_sample_key=m, draw_key=day,
            physical_capture_key=r['gsm'] if suffix_ok else '',
            current_analysis_channel_id=r['gsm'],
            membership_kind='intended',
            physical_capture_status=status,
            analysis_matches_physical='true' if suffix_ok else 'unknown',
            evidence_source='cillo_gse180578_channel_inventory.csv',
            evidence_locator='channel_inventories/cillo_gse180578_channel_inventory.csv '
                             ':: one row per gsm; column members split on ";" then "@"; '
                             'h5ad_suffix / join_status',
            evidence_grade='measured_metadata',
            ambiguity_reason=reason))

# --- Stephenson ------------------------------------------------------------
# The SDRF resolves samples to SEQUENCE POOLS, not GEM wells. A sequencing lane is
# a flow-cell partition and a GEM well is a droplet-generation channel, so no well
# count follows from a lane count. Every row is unresolved, including the 83
# one-sample pools: lacking a pool prefix is not evidence of a dedicated well.
for r in read_columns(os.path.join(T3, 'stephenson_patient_outcomes.csv'),
                      ['patient_id', 'n_samples', 'sample_ids', 'pools',
                       'Collection_Days']):
    for s in [x for x in r['sample_ids'].split(';') if x.strip()]:
        pooled = '_' in s and s.split('_')[0].startswith('BGCV')
        pool = s.split('_')[0] if pooled else s
        rows.append(row(
            cohort='Stephenson 2021', accession='E-MTAB-10026',
            donor_key=r['patient_id'], biological_sample_key=s,
            draw_key=r['Collection_Days'],
            physical_capture_key='',
            current_analysis_channel_id=s,
            membership_kind='intended',
            physical_capture_status='unresolved',
            analysis_matches_physical='unknown',
            evidence_source='stephenson_patient_outcomes.csv + T7 memo.md SDRF reading',
            evidence_locator='T3_patient_sample_mapping/stephenson_patient_outcomes.csv '
                             ':: columns patient_id, sample_ids (";"-split), pools; '
                             'T7_doublet_prior_evidence/memo.md :: Stephenson section, '
                             'SDRF "sequence pool" = 98 values (83 of one, 15 of four)',
            evidence_grade='published',
            ambiguity_reason=(
                f'sequence pool {pool} was sequenced across two lanes with the lane '
                'assignment unpublished, and the pool-to-GEM-well relation is absent'
                if pooled else
                'one-sample sequence pool; the deposit does not state that it is a '
                'single GEM well, and absence of a pool prefix is not evidence of one')))

# --- COMBAT ----------------------------------------------------------------
# Channel_ID is the GEM well (70 = 10 Pool_ID x 7 sorted regions); COMBAT_ID is
# the donor. Produced by extract_hpcc_membership.py from the deposited h5ad,
# reading only those two obs datasets.
#
# `observed_assigned`, NOT `intended`. These pairs are donors observed among the
# RETAINED cells of a channel after the depositors' own QC. No pool-design
# artifact listing the intended donors per channel exists in this repository --
# `combat_channel_inventory.csv` carries only an `n_donors` count -- so a donor
# absent from a channel here cannot be called unloaded. Distinct from
# `observed_singlet`, which is genotype demultiplexing: this is the deposit's own
# per-cell donor assignment.
for r in read_columns(os.path.join(HERE, 'combat_donor_by_channel.csv'),
                      ['channel_id', 'donor_id']):
    rows.append(row(
        cohort='COMBAT 2022', accession='COMBAT-CITESeq',
        donor_key=r['donor_id'], biological_sample_key=r['donor_id'],
        draw_key='', physical_capture_key=r['channel_id'],
        current_analysis_channel_id=r['channel_id'],
        membership_kind='observed_assigned', physical_capture_status='resolved',
        analysis_matches_physical='true',
        evidence_source='COMBAT-CITESeq-DATA.h5ad obs',
        evidence_locator='combat_donor_by_channel.csv :: one row per '
                         '(obs/Channel_ID, obs/COMBAT_ID) pair; produced by '
                         'extract_hpcc_membership.py, raw identity in '
                         'RAW_INPUT_MANIFEST.txt',
        evidence_grade='measured_metadata',
        ambiguity_reason='donor observed among the channel\'s retained cells after '
                         'the depositors\' QC; no pool-design artifact states the '
                         'intended donors per channel, so absence of a pair is not '
                         'evidence that the donor was not loaded'))

# --- Liu -------------------------------------------------------------------
# OBSERVED singlet membership, and PROVISIONAL capture: the lane token comes from
# the deposit's GSM titles, but reading that token as one 10x chip well is this
# repository's gloss, not a statement in the deposit. `resolved` is defined as
# "the deposit identifies the capture unit directly", which this is not.
for r in read_columns(os.path.join(HERE, 'liu_donor_by_lane.csv'),
                      ['lane', 'donor']):
    rows.append(row(
        cohort='Liu 2021', accession='GSE161918',
        donor_key=r['donor'], biological_sample_key=r['donor'],
        draw_key='unresolved_mixed_pool',
        physical_capture_key=r['lane'], current_analysis_channel_id=r['lane'],
        membership_kind='observed_singlet',
        physical_capture_status='provisional',
        analysis_matches_physical='unknown',
        evidence_source='GSE161918 demuxlet .best',
        evidence_locator='liu_donor_by_lane.csv :: one row per (lane, '
                         'SNG.BEST.GUESS) where DROPLET.TYPE=="SNG"; produced by '
                         'extract_hpcc_membership.py, raw identity in '
                         'RAW_INPUT_MANIFEST.txt',
        evidence_grade='derived',
        ambiguity_reason='post-demultiplex singlet recovery, not intended pool '
                         'composition; and the lane token is read as a 10x chip '
                         "well by this repository's annotation, not by the deposit"))

# ---------------------------------------------------------------------------
# Canonical-donor coverage. Its OWN table, never membership rows.
# ---------------------------------------------------------------------------
CANON = read_columns(os.path.join(LIT, 'canonical_donor_labels.csv'),
                     ['cohort', 'accession', 'donor_id'])
COHORT_ACC = {'Qiu 2021': 'GSE167363', 'Chongqing sepsis 2025': 'GSE279452',
              'Ren 2021': 'GSE158055', 'Wilk 2021 JEM': 'GSE174072',
              'Cillo 2021': 'GSE180578', 'Stephenson 2021': 'E-MTAB-10026',
              'COMBAT 2022': 'COMBAT-CITESeq', 'Liu 2021': 'GSE161918'}


def norm(cohort, donor):
    """Cohort-specific normalisation, stated rather than applied silently."""
    if cohort == 'Liu 2021':
        return re.sub(r'_(B\d|T\d)$', '', donor)      # demuxlet lane/draw suffix
    if cohort == 'Chongqing sepsis 2025':
        return donor.replace(' ', '').upper()         # "pHC 01" vs "PHC01"
    return donor


NORM_RULE = {'Liu 2021': 'strip trailing _B<n>/_T<n> from the demuxlet genotype id',
             'Chongqing sepsis 2025': 'remove spaces and upper-case'}

recon = []
for cohort, acc in COHORT_ACC.items():
    have = {norm(cohort, r['donor_key']) for r in rows if r['cohort'] == cohort}
    canon_n = {norm(cohort, r['donor_id']): r['donor_id']
               for r in CANON if r['cohort'] == cohort}
    missing = sorted(d for n, d in canon_n.items() if n not in have)
    extra = sorted(have - set(canon_n))
    recon.append({'cohort': cohort, 'accession': acc,
                  'n_canonical_donors': len(canon_n),
                  'n_donors_with_membership': len({r['donor_key'] for r in rows
                                                   if r['cohort'] == cohort}),
                  'normalisation_rule': NORM_RULE.get(cohort, 'none'),
                  'n_canonical_without_membership': len(missing),
                  'canonical_without_membership': ';'.join(missing),
                  'n_membership_donors_not_canonical': len(extra),
                  'membership_donors_not_canonical': ';'.join(extra)})

# ---------------------------------------------------------------------------
# Invariants. Enumerations only -- nothing interpretive is asserted.
# ---------------------------------------------------------------------------
problems = []
for i, r in enumerate(rows):
    where = f'row {i} ({r["cohort"]}/{r["biological_sample_key"]})'
    for k in ('cohort', 'accession', 'donor_key', 'biological_sample_key',
              'current_analysis_channel_id'):
        if not str(r[k]).strip():
            problems.append(f'{where}: empty {k}')
    if r['physical_capture_status'] == 'resolved' and not r['physical_capture_key'].strip():
        problems.append(f'{where}: resolved with no physical_capture_key')
    if r['physical_capture_status'] != 'resolved' and not r['ambiguity_reason'].strip():
        problems.append(f'{where}: non-resolved with no ambiguity_reason')
    if '::' not in r['evidence_locator']:
        problems.append(f'{where}: evidence_locator names no column/row rule')
    if r['membership_kind'] == 'observed_singlet' and r['evidence_grade'] == 'measured_metadata':
        problems.append(f'{where}: observed membership graded as measured intended metadata')
    if r['physical_capture_status'] == 'unresolved' and r['analysis_matches_physical'] != 'unknown':
        problems.append(f'{where}: unresolved capture compared against the analysis channel')
    if r['membership_kind'] == 'unresolved':
        problems.append(f'{where}: coverage-only row must not be in the membership table')
    if r['membership_kind'].startswith('observed') and not r['ambiguity_reason'].strip():
        problems.append(f'{where}: observed membership with no stated limitation; '
                        'absence of a row must never read as absence from the capture')

for coh in {r['cohort'] for r in rows}:
    sub = [r for r in rows if r['cohort'] == coh]
    key = Counter((r['donor_key'], r['biological_sample_key'],
                   r['physical_capture_key'], r['current_analysis_channel_id'])
                  for r in sub)
    dups = {k: v for k, v in key.items() if v > 1}
    if dups:
        problems.append(f'{coh}: {len(dups)} duplicated membership keys')

EXPECT = {              # (donors, biological samples, membership rows)
    'Qiu 2021': (7, 12, 12),
    'Chongqing sepsis 2025': (133, 133, 133),
    'Ren 2021': (196, 284, 284),
    'Wilk 2021 JEM': (23, 23, 23),
    'Cillo 2021': (29, 70, 71),
    'Stephenson 2021': (130, 143, 143),
    'COMBAT 2022': (124, 124, 979),
    'Liu 2021': (50, 50, 877),
}
for coh, (nd, ns, nr) in EXPECT.items():
    sub = [r for r in rows if r['cohort'] == coh]
    got = (len({r['donor_key'] for r in sub}),
           len({r['biological_sample_key'] for r in sub}), len(sub))
    if got != (nd, ns, nr):
        problems.append(f'{coh}: expected donors/samples/rows {(nd, ns, nr)}, measured {got}')

if problems:
    print('INVARIANT FAILURES:', file=sys.stderr)
    for p in problems:
        print('  ' + p, file=sys.stderr)
    raise SystemExit(1)


def write_csv(name, fieldnames, records):
    with open(os.path.join(HERE, name), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(records)


# Files CITED as evidence but never parsed still belong in the frozen checksum
# list: the Stephenson rows rest on the T7 memo's SDRF reading, and an input whose
# checksum is absent can change without the bundle's provenance noticing.
for cited in (os.path.join(LIT, 'T7_doublet_prior_evidence/memo.md'),):
    with open(cited, 'rb') as fh:
        INPUTS[os.path.relpath(cited, REPO) + '  (cited, not parsed)'] = \
            hashlib.sha256(fh.read()).hexdigest()

write_csv('capture_provenance_long.csv', FIELDS, rows)
write_csv('donor_coverage_reconciliation.csv', list(recon[0]), recon)
write_csv('outcome_column_boundary.csv', list(COLUMN_BOUNDARY[0]), COLUMN_BOUNDARY)

# --- cohort summary --------------------------------------------------------
# An empty string is NOT a category. A previous version counted the empty
# `current_analysis_channel_id` of a coverage-only row as one channel, inflating
# Chongqing 133->134, Cillo 44->45 and Wilk 23->24. Those rows are gone, and the
# counts filter empties anyway so the defect cannot come back silently.
SUM = ['cohort', 'accession', 'n_donors', 'n_biological_samples', 'n_draws',
       'n_physical_captures', 'n_analysis_channels', 'n_membership_rows',
       'physical_capture_status', 'analysis_matches_physical',
       'max_donors_per_capture', 'max_captures_per_donor']


def nonempty(sub, field):
    return {r[field] for r in sub if str(r[field]).strip()}


summary = []
for coh in sorted({r['cohort'] for r in rows}):
    sub = [r for r in rows if r['cohort'] == coh]
    per_cap = Counter(r['physical_capture_key'] for r in sub
                      if str(r['physical_capture_key']).strip())
    per_don = defaultdict(set)
    for r in sub:
        if str(r['physical_capture_key']).strip():
            per_don[r['donor_key']].add(r['physical_capture_key'])
    summary.append({
        'cohort': coh,
        'accession': sorted({r['accession'] for r in sub})[0],
        'n_donors': len(nonempty(sub, 'donor_key')),
        'n_biological_samples': len(nonempty(sub, 'biological_sample_key')),
        'n_draws': len({(r['donor_key'], r['draw_key']) for r in sub
                        if str(r['draw_key']).strip()}),
        'n_physical_captures': len(nonempty(sub, 'physical_capture_key')),
        'n_analysis_channels': len(nonempty(sub, 'current_analysis_channel_id')),
        'n_membership_rows': len(sub),
        'physical_capture_status': '|'.join(f'{k}={v}' for k, v in sorted(
            Counter(r['physical_capture_status'] for r in sub).items())),
        'analysis_matches_physical': '|'.join(f'{k}={v}' for k, v in sorted(
            Counter(r['analysis_matches_physical'] for r in sub).items())),
        'max_donors_per_capture': max(per_cap.values()) if per_cap else 0,
        'max_captures_per_donor': max((len(v) for v in per_don.values()), default=0),
    })
write_csv('cohort_summary.csv', SUM, summary)

with open(os.path.join(HERE, 'INPUT_CHECKSUMS.txt'), 'w') as fh:
    for p in sorted(INPUTS):
        fh.write(f'{INPUTS[p]}  {p}\n')

print(f'membership rows: {len(rows)}')
for s in summary:
    print(f"  {s['cohort']:<24} donors={s['n_donors']:>4} samples={s['n_biological_samples']:>4} "
          f"captures={s['n_physical_captures']:>4} analysis_ch={s['n_analysis_channels']:>4} "
          f"rows={s['n_membership_rows']:>4}  {s['physical_capture_status']}")
print('coverage-only donors (reconciliation table only): '
      f'{sum(r["n_canonical_without_membership"] for r in recon)}')
print(f'inputs read with usecols: {len(COLUMN_BOUNDARY)}; inputs where an outcome '
      f'column was excluded at parse time: '
      f'{sum(1 for c in COLUMN_BOUNDARY if c["outcome_columns_excluded"])}')
print('all invariants held')
