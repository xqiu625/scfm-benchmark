"""READ-ONLY: render the sensitivity-design evidence-gap matrix AND revalidate it
against the T7 artifacts it claims to summarize.

Two jobs, kept separate on purpose:

  RENDER    the per-cohort narrative fields below are a MANUAL ADJUDICATION --
            a human read the sources and wrote down what they say. A script
            cannot derive "the deposit does not state which side of the filter
            it sits on" from a CSV.

  REVALIDATE every COUNTABLE claim the packet makes is re-derived here from the
            T7 artifacts and asserted. If T7 changes, this fails instead of
            letting the packet drift silently.

An earlier version defined the T7 path and never read it, so `regenerate with
this script` only reprinted the transcription. That gap is the failure mode this
split exists to prevent: regenerate is not revalidate.

Selects no rate, writes no runtime config, runs nothing, touches no checkpoint,
reads no outcome field.
"""
import argparse
import collections
import csv
import hashlib
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, '..', '..', '..'))
T7 = os.path.join(REPO, 'reports', 'manuscript', 'revision', 'literature_audit',
                  'T7_doublet_prior_evidence')
OUT = os.path.join(REPO, 'reports', 'manuscript', 'revision', 'sensitivity_design_review')
LIVING = os.path.join(REPO, 'docs', 'current-status.md')

#: Loading granularity, as an EXPLICIT per-cohort field. It was previously
#: inferred with `'per-capture' in f['loading']`, which matched Ren's
#: "... NOT per-capture" and filed Ren as capture-addressable -- contradicting
#: Q2 in the same document. Granularity is a classification, not something to
#: read out of prose.
GRAN_MEASURED = 'capture-specific measured'
GRAN_ADDRESSABLE = 'capture-addressable protocol value (NOT measured per capture)'
GRAN_UNASSIGNABLE = 'cohort-or-protocol-level, unassignable to a capture'
GRAN_NONE = 'none'
VALID_GRANULARITIES = (GRAN_MEASURED, GRAN_ADDRESSABLE, GRAN_UNASSIGNABLE, GRAN_NONE)

COHORTS = ['qiu', 'chongqing', 'cillo', 'liu', 'combat', 'ren', 'stephenson', 'wilk']

# ---------------------------------------------------------------------------
# Claims the packet makes that ARE countable. Each is re-derived from T7 below.
# ---------------------------------------------------------------------------
# (cohort, mapping_status) -> n. Keyed on the STATUS too: an earlier version
# counted only by cohort, so Stephenson's 143 could have flipped from
# unmapped_to_gem_well to mapped without the check failing while the README
# went on claiming 0/143.
EXPECTED_MAPPING = {
    ('qiu', 'mapped'): 12, ('chongqing', 'mapped'): 124, ('cillo', 'mapped'): 35,
    ('liu', 'mapped'): 55, ('combat', 'mapped'): 70, ('ren', 'mapped'): 284,
    ('stephenson', 'unmapped_to_gem_well'): 143, ('wilk', 'inherited'): 23,
}
EXPECTED_SOURCE_COUNTS = {'doublet_prior_capture_evidence.csv': 605,
                          'doublet_prior_source_audit.csv': 23,
                          'gate2_source_crosswalk.csv': 790}
EXPECTED_VENDOR_DOCS = 4
COMBAT_CELLS_LOADED = 50_000
EXPECTED_RATE_ROWS = 21
EXPECTED_PAIRING_ROWS = 58
EXPECTED_MEASURED_PAIRING = {'Liu': 55}          # the ONLY measured donor composition
EXPECTED_ASSUMED_PAIRING = {'COMBAT': 1, 'Cillo (2-donor libraries)': 1,
                            'Stephenson (Cambridge pools)': 1}
EXPECTED_REN_SAMPLES = 284

# ---------------------------------------------------------------------------
# Manual adjudication. `source` names where a human read it; `*_grade` says how
# strong it is. Nothing here is re-derivable by script.
# ---------------------------------------------------------------------------
FACTS = {
    'qiu': dict(
        platform='10x Chromium', capture_unit='GSM = GEM well',
        input_stage='Cell Ranger filtered triplet, before any doublet step',
        chemistry="CONFLICT: GEO says 3' v3.1, paper text says 5' UMI counting",
        chemistry_grade='published',
        loading='15,000 target cells per sample (protocol figure, uniform across samples)',
        loading_grade='published',
        loading_granularity=GRAN_ADDRESSABLE,
        recovered='unknown', recovered_grade='unknown',
        upstream_doublet='none reported', upstream_grade='published',
        donor_proportions='n/a (single donor per capture)', donor_grade='n/a',
        classification='standard_gate2_candidate',
        source='doublet_prior_capture_evidence.csv; J Leukoc Biol 2021 Methods',
        unidentifiable=('recovered count; which chemistry table applies; residual '
                        'prevalence of the post-Gate-1 matrix')),
    'chongqing': dict(
        platform='10x Chromium', capture_unit='GSM',
        input_stage='Cell Ranger filtered triplet, before any doublet step',
        chemistry="3' v3.1 and/or 5' VDJ, not resolvable per sample",
        chemistry_grade='published',
        loading='NOT REPORTED (no linked publication)', loading_grade='unknown',
        loading_granularity=GRAN_NONE,
        recovered='unknown', recovered_grade='unknown',
        upstream_doublet='none reported', upstream_grade='unknown',
        donor_proportions='n/a (single donor per capture)', donor_grade='n/a',
        classification='standard_gate2_candidate',
        source='doublet_prior_capture_evidence.csv; GEO extract protocols',
        unidentifiable=('loading; recovered count; chemistry per sample; whether "none '
                        'reported" reflects absence of a step or absence of a paper')),
    'cillo': dict(
        platform='10x Chromium Chip A', capture_unit='GSM',
        input_stage='Cell Ranger filtered, after HTO cross-sample exclusion',
        chemistry="5' v1", chemistry_grade='published',
        loading='NOT REPORTED', loading_grade='unknown',
        loading_granularity=GRAN_NONE,
        recovered='unknown', recovered_grade='unknown',
        upstream_doublet=('HTO double-positive exclusion, CROSS-SAMPLE ONLY; no '
                          'within-sample caller'),
        upstream_grade='published',
        donor_proportions=('equal mixture ASSUMED for the 2-donor libraries; per-library '
                           'composition never measured'),
        donor_grade='assumed',
        classification='upstream_doublet_filter_already_applied',
        source='doublet_prior_capture_evidence.csv; pooled-same-donor-pairing-scenarios-v3.csv',
        unidentifiable=('loading; HTO removal efficiency; per-library sample composition; '
                        'within-sample residual')),
    'liu': dict(
        platform='10x Chromium', capture_unit='lane',
        input_stage=("after demuxlet genotype demultiplexing, then this project's "
                     'inflection cell call'),
        chemistry="5' Next GEM v1.1", chemistry_grade='published',
        loading='NOT REPORTED', loading_grade='unknown',
        loading_granularity=GRAN_NONE,
        recovered='NOT OBTAINABLE (no raw droplet matrix deposited)',
        recovered_grade='unknown',
        upstream_doublet=('demuxlet v2 genotype demultiplexing (cross-donor only); '
                          'same-donor pairs untouched'),
        upstream_grade='published',
        donor_proportions=('MEASURED per lane (55 lanes) from demuxlet SNG composition -- '
                           'but measured AFTER DBL/AMB exclusion, so it is a STAGE-MISMATCHED '
                           'proxy for loading-stage p_d'),
        donor_grade='proxy',
        classification='upstream_doublet_filter_already_applied',
        source='liu_demuxlet_per_lane.csv; pooled-same-donor-pairing-scenarios-v3.csv',
        unidentifiable=('loading; recovered count; demuxlet removal efficiency; '
                        'loading-stage donor composition; same-donor residual')),
    'combat': dict(
        platform='10x Chromium Chip G', capture_unit='Channel_ID = GEM well',
        input_stage='post-study-QC AnnData, after demuxlet + Vireo',
        chemistry="5' scCITE-seq", chemistry_grade='published',
        loading='50,000 cells per channel (protocol figure, uniform across channels)',
        loading_grade='published',
        loading_granularity=GRAN_ADDRESSABLE,
        recovered='unknown', recovered_grade='unknown',
        upstream_doublet=('demuxlet v2 + Vireo genotype (cross-donor only); Scrublet score '
                          'retained as a QC covariate, NOT used to remove'),
        upstream_grade='published',
        donor_proportions='equal mixture ASSUMED, 14 donors; per-channel composition never measured',
        donor_grade='assumed',
        classification='upstream_doublet_filter_already_applied',
        source='doublet_prior_capture_evidence.csv; pooled-same-donor-pairing-scenarios-v3.csv',
        unidentifiable=('recovered count; per-channel donor composition; genotype removal '
                        'efficiency; same-donor residual; whether the vendor table extends '
                        'to 50,000 loaded')),
    'ren': dict(
        platform='10x Chromium', capture_unit='sampleID',
        input_stage=("after the authors' own Scrublet AND a further undefined manual removal"),
        chemistry="10x 3'", chemistry_grade='published',
        loading='cohort-level range only: 300-14,000 target cells per reaction',
        loading_grade='published',
        loading_granularity=GRAN_UNASSIGNABLE,
        recovered='unknown', recovered_grade='unknown',
        upstream_doublet=("authors' Scrublet removed 131,962 cells (7.6%); manual removal a "
                          'further 136,006 (8.5%), content undefined'),
        upstream_grade='measured',
        donor_proportions='n/a (single donor per capture)', donor_grade='n/a',
        classification='upstream_doublet_filter_already_applied + residual_scoring_only',
        source='ren_gse158055_capture_qc.csv (284 rows, Cell 2021 Table S1)',
        unidentifiable=('per-capture loading; what manual removal removed; residual after two '
                        'upstream passes; barcode-level correspondence (counts match 284/284, '
                        'identity never compared)')),
    'stephenson': dict(
        platform='10x Chromium', capture_unit='NOT RESOLVABLE to GEM well',
        input_stage='processed objects; per-capture matrices not resolvable',
        chemistry="CONFLICT: paper 5' NextGEM v1.1, SDRF annotates 10xV3",
        chemistry_grade='published',
        loading=('20,000-30,000 per sample (Newcastle/London); 2 lanes x 25,000 per pool of '
                 'four (Cambridge) -- protocol level, NOT attachable to a capture'),
        loading_grade='published',
        loading_granularity=GRAN_UNASSIGNABLE,
        recovered='unknown', recovered_grade='unknown',
        upstream_doublet=('SITE-DEPENDENT: Scrublet (Newcastle), emptyDrops + genotype tagging '
                          '(Cambridge), souporcell (London)'),
        upstream_grade='published',
        donor_proportions='equal mixture ASSUMED for Cambridge pools of four; never measured',
        donor_grade='assumed',
        classification='not_eligible_without_additional_evidence',
        source='doublet_prior_capture_evidence.csv (cohort-level row); capture_evidence_join_key.csv',
        unidentifiable=('capture identity itself -- there is no legitimate unit of analysis; '
                        'which site processed which capture; chemistry; everything downstream')),
    'wilk': dict(
        platform='Seq-Well S3 (NOT 10x)', capture_unit='converted directory basename',
        input_stage='author cell-count matrices; stage relative to the multiplet filter NOT stated',
        chemistry='Seq-Well array + ChemGenes barcoded beads', chemistry_grade='published',
        loading='15,000 cells per array, uniform protocol figure', loading_grade='published',
        loading_granularity=GRAN_UNASSIGNABLE,
        recovered='unknown', recovered_grade='unknown',
        upstream_doublet=('heuristic: cells with >75 genes per 100 UMIs removed as putative '
                          'multiplets; deposit does not say whether this was applied before '
                          'deposition'),
        upstream_grade='published',
        donor_proportions='n/a (single donor per array)', donor_grade='n/a',
        classification='not_eligible_without_additional_evidence',
        source='doublet_prior_capture_evidence.csv (cohort-level row)',
        unidentifiable=("which side of the authors' multiplet filter the deposit sits on; "
                        'whether ANY 10x-derived multiplet model applies to Seq-Well at all')),
}

COLUMNS = ['cohort', 'platform', 'capture_unit', 'input_stage',
           'chemistry', 'chemistry_evidence_grade',
           'loading_evidence', 'loading_evidence_grade', 'loading_granularity',
           'recovered_count', 'recovered_count_grade',
           'upstream_doublet_handling', 'upstream_evidence_grade',
           'donor_proportion_source', 'donor_proportion_grade',
           'donor_proportion_measured_at_stage',
           'capture_mapping_coverage', 'runnability_classification',
           'evidence_source', 'still_unidentifiable']


def _read(name):
    with open(os.path.join(T7, name), newline='', encoding='utf-8-sig') as fh:
        return list(csv.DictReader(fh))


def revalidate():
    """Re-derive every countable claim from T7. Raises on any disagreement."""
    checks = []

    jk = _read('capture_evidence_join_key.csv')
    cov = collections.Counter((r['cohort'], r['mapping_status']) for r in jk)
    assert dict(cov) == EXPECTED_MAPPING, f'mapping coverage drifted: {dict(cov)}'
    assert sum(cov.values()) == 746, sum(cov.values())
    checks.append('capture coverage AND mapping status per cohort match (746 total; '
                  'Stephenson 143 unmapped_to_gem_well, Wilk 23 inherited)')

    for name, n in EXPECTED_SOURCE_COUNTS.items():
        actual = len(_read(name))
        assert actual == n, f'{name}: {actual} rows, packet says {n}'
    checks.append('source-list row counts match (605 / 23 / 790)')

    rate = _read('rate-sensitivity-scenarios-v3.csv')
    assert len(rate) == EXPECTED_RATE_ROWS, len(rate)
    elig = {r['runtime_eligible'].strip().lower() for r in rate}
    assert elig == {'false'}, elig
    checks.append(f'{len(rate)}/{len(rate)} rate scenarios runtime_eligible=false')

    pair = _read('pooled-same-donor-pairing-scenarios-v3.csv')
    assert len(pair) == EXPECTED_PAIRING_ROWS, len(pair)
    elig = {r['runtime_eligible'].strip().lower() for r in pair}
    assert elig == {'false'}, elig
    checks.append(f'{len(pair)}/{len(pair)} pairing scenarios runtime_eligible=false')

    measured = collections.Counter(
        r['cohort'] for r in pair if r['source_of_proportions'].startswith('MEASURED'))
    assert dict(measured) == EXPECTED_MEASURED_PAIRING, dict(measured)

    # Cross-check against the SOURCE lane artifact, not just the count in the
    # scenario table. Agreeing on the number 55 is not the same as agreeing on
    # WHICH 55 lanes: an earlier version only counted rows, so a scenario row
    # naming a lane absent from the demuxlet output would have passed.
    lanes_src = _read('liu_demuxlet_per_lane.csv')
    ids_src = [r['lane'] for r in lanes_src]
    ids_pair = [r['lane'] for r in pair if r['cohort'] == 'Liu']
    assert len(ids_src) == len(set(ids_src)), 'duplicate lane in liu_demuxlet_per_lane.csv'
    assert len(ids_pair) == len(set(ids_pair)), 'duplicate Liu lane in the pairing table'
    assert set(ids_src) == set(ids_pair), (
        f'Liu lane IDs disagree: only in source {sorted(set(ids_src) - set(ids_pair))[:3]}, '
        f'only in scenarios {sorted(set(ids_pair) - set(ids_src))[:3]}')
    checks.append(f'Liu is the ONLY cohort with measured donor composition, and its '
                  f'{len(ids_pair)} scenario lanes match the demuxlet source lane IDs '
                  'one-for-one (no duplicates on either side)')

    assumed = collections.Counter(
        r['cohort'] for r in pair if r['source_of_proportions'].startswith('ASSUMED'))
    assert dict(assumed) == EXPECTED_ASSUMED_PAIRING, dict(assumed)
    checks.append('COMBAT / Cillo / Stephenson each carry exactly one assumed row')

    # the measured rows must still be flagged proxy in the matrix, because they
    # were measured at the post-demultiplex stage, not at loading
    stages = {r['proportions_measured_at_stage'] for r in pair
              if r['source_of_proportions'].startswith('MEASURED')}
    assert all('post-demultiplex' in s for s in stages), stages
    assert FACTS['liu']['donor_grade'] == 'proxy', \
        "Liu's measured composition is a stage-mismatched proxy for loading p_d"
    checks.append("Liu's measured composition is stage-mismatched -> graded proxy")

    ren = _read('ren_gse158055_capture_qc.csv')
    assert len(ren) == EXPECTED_REN_SAMPLES, len(ren)
    tot = sum(int(r['Ncells_manual_removal']) for r in ren if r['Ncells_manual_removal'].strip())
    assert tot == 1_462_702, tot
    checks.append(f'Ren: {len(ren)} samples, Ncells_manual_removal sums to {tot:,}')

    # Classification agreement, per cohort. `needle in whole_document` was too
    # weak: Stephenson and Wilk share a needle, and any other mention of it
    # anywhere in the file satisfied the check. Match within the line that names
    # the cohort instead.
    living_lines = open(LIVING).read().splitlines()
    for cohort, needle in [('Ren', 'residual_scoring_only'),
                           ('Stephenson', 'not_eligible_without_additional_evidence'),
                           ('Wilk', 'not_eligible_without_additional_evidence')]:
        assert needle in FACTS[cohort.lower()]['classification'], cohort
        hits = [ln for ln in living_lines
                if cohort in ln and needle in ln and ln.lstrip().startswith('|')]
        assert hits, (f'{LIVING}: no table row names {cohort} together with {needle!r}; '
                      'the classification table may have changed')
    checks.append('Ren / Stephenson / Wilk classifications matched on the cohort\'s own '
                  'table row in current-status.md')

    # Granularity is an explicit field, so a typo or a new cohort cannot fall
    # through. It was previously inferred by a substring test that matched Ren's
    # "... NOT per-capture" and filed Ren as capture-addressable -- contradicting
    # Q2 inside the same document, in the two places a reviewer reads first.
    bad_gran = {c: FACTS[c]['loading_granularity'] for c in COHORTS
                if FACTS[c]['loading_granularity'] not in VALID_GRANULARITIES}
    assert not bad_gran, f'unrecognized loading_granularity: {bad_gran}'
    assert all('loading_granularity' in FACTS[c] for c in COHORTS)
    assert not any(FACTS[c]['loading_granularity'] == GRAN_MEASURED for c in COHORTS), (
        'no cohort has capture-specific MEASURED loading; Q2 says so and the matrix must agree')
    assert FACTS['ren']['loading_granularity'] == GRAN_UNASSIGNABLE, (
        "Ren's loading is a cohort-level range and cannot be attached to a capture")
    checks.append('loading granularity is an explicit per-cohort field; none is '
                  '"capture-specific measured"; Ren is unassignable')

    # Recovered count, checked against the ARTIFACTS rather than against FACTS
    # (the previous version asserted its own hardcoded dict, which would still
    # pass if T7 later gained a count -- self-certifying, not validating).
    #
    # SCOPE. These are the two COHORT-LEVEL recovered-count fields. They are not
    # the only fields in T7 whose name suggests recovery:
    # tenx_multiplet_rate_tables.csv carries `cells_recovered`, populated and
    # numeric in all 54 rows -- but that is the VENDOR's own published
    # loaded -> recovered reference mapping, not a measurement of any cohort's
    # capture. Excluded deliberately, and named here so the exclusion is visible
    # rather than implicit.
    cw = _read('gate2_source_crosswalk.csv')
    assert 'n_input_barcodes' in cw[0], 'crosswalk schema changed'
    populated = [r for r in cw if r['n_input_barcodes'].strip()]
    assert not populated, f'{len(populated)} crosswalk rows now carry n_input_barcodes'

    # Sentinel-based, NOT float()-based. `~5000`, `5000 cells`, `estimated 5,000`
    # and `5k` are all recovered-count evidence but none parses as a float, so
    # catching ValueError and treating the value as prose would wave them
    # through. Anything outside the allowed sentinels trips the check instead.
    RECOVERED_SENTINELS = {
        'not obtainable (no raw droplet matrix deposited)',
        'not reported', 'not obtainable', 'na', 'n/a', 'unknown', '',
    }
    ev = _read('doublet_prior_capture_evidence.csv')
    suspect = [(r['cohort'], v) for r in ev
               if (v := r.get('estimated_or_recovered_cells', '').strip())
               and v.lower() not in RECOVERED_SENTINELS]
    assert not suspect, (
        f'estimated_or_recovered_cells holds a value outside the agreed sentinels, which may '
        f'be a recovered count and needs human review: {suspect[:3]}')
    checks.append('no cohort-level recovered count in either canonical field '
                  '(n_input_barcodes empty across 790 crosswalk rows; '
                  'estimated_or_recovered_cells holds only agreed sentinels across 605 '
                  'evidence rows; the vendor table\'s own cells_recovered column is '
                  'reference data and is excluded by design)')
    assert all(FACTS[c]['recovered_grade'] == 'unknown' for c in COHORTS), \
        'the matrix must agree with the artifacts it summarizes'

    # COMBAT's loading is outside every vendor table on file
    vendor = _read('tenx_multiplet_rate_tables.csv')
    docs = {r['document'] for r in vendor}
    assert len(docs) == EXPECTED_VENDOR_DOCS, docs
    max_loaded = max(float(r['cells_loaded']) for r in vendor if r['cells_loaded'].strip())
    assert COMBAT_CELLS_LOADED > max_loaded, (
        f'COMBAT {COMBAT_CELLS_LOADED} no longer exceeds the vendor maximum {max_loaded}')
    checks.append(f'COMBAT {COMBAT_CELLS_LOADED:,} loaded exceeds every vendor table '
                  f'(max {max_loaded:,.0f} across {len(docs)} documents)')
    return checks


HEADER_TEX = r'''\usepackage{fancyhdr}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small 302006 --- Gate 2 sensitivity-design review packet}
\fancyhead[R]{\small\thepage}
\fancyfoot[C]{\scriptsize Prepared 2026-09-22 \textbullet{} no rate selected \textbullet{} no runtime configuration exists \textbullet{} nothing here authorises deleting a cell}
\renewcommand{\headrulewidth}{0.4pt}
\setlength{\emergencystretch}{5em}

% --- TYPOGRAPHY ONLY (2026-09-24). Changes no value, no rate, no classification.
% The 2026-09-22 build ran text off the PHYSICAL page edge on three pages,
% measured with `pdftotext -bbox` against the 612pt page width:
%   page 2  xMax 613.6pt  scripts/32_rebuild_upfront_qc/code/build_...
%   page 3  xMax 615.1pt  not_eligible_without_a...
%   page 4  xMax 615.1pt  not_eligible_without_a...
% Both causes are unbreakable words, which no amount of \emergencystretch can fix:
% pandoc escapes `_` as `\_`, at which TeX will not break, and typewriter text is
% unhyphenated by default.
% Break after each underscore, adding NO character. hyphenat[htt] was tried and
% rejected: it fixed the width but inserted a real hyphen inside file paths and
% identifiers ("scripts/32_rebuild_upfront_-qc"), which a reviewer could copy as
% part of the path. A discretionary break with an empty pre-break has no such
% hazard, and snake_case gives TeX plenty of break points on its own.
\let\GtwoUnderscore\_
\renewcommand{\_}{\GtwoUnderscore\discretionary{}{}{}}
'''

PDF_NAME = 'sensitivity_design_review_packet.pdf'
#: Rendered here first. pandoc writing straight to the tracked path means a
#: failure after truncation leaves a corrupt frozen artifact; the real file
#: is only replaced once the render succeeded and the output looks like a PDF.
PDF_TMP = '_sensitivity_design_review_packet.inprogress.pdf'

#: The ONE definition of how the PDF is rendered. Previously the command lived
#: in two places -- printed here and transcribed into BUILD_PROVENANCE.md -- and
#: they had already drifted: the printed form was missing the title, subtitle
#: and date, so anyone rebuilding from it got a different title page.
PANDOC_ARGV = [
    'pandoc', '_packet_full.md', '-s', '--toc', '--toc-depth=2',
    '--pdf-engine=xelatex',
    '-V', 'geometry:margin=2.2cm', '-V', 'fontsize=10pt',
    '-V', 'mainfont=Helvetica Neue', '-V', 'monofont=Menlo',
    '-V', 'colorlinks=true', '-V', 'linkcolor=black', '-V', 'urlcolor=blue',
    '-V', 'title=Gate 2 doublet rate: sensitivity-design review packet',
    '-V', 'subtitle=Project 302006 — five-FM platelet mortality benchmark',
    '-V', 'date=Prepared 2026-09-22 · for independent methodological review',
    '-H', '_header.tex',
    '-o', PDF_TMP,
]

MATRIX_FIELDS = [
    ('platform', 'Platform'), ('capture_unit', 'Capture unit'),
    ('input_stage', 'Input stage'),
    ('chemistry', 'Chemistry'), ('chemistry_evidence_grade', '-- grade'),
    ('loading_evidence', 'Loading evidence'), ('loading_evidence_grade', '-- grade'),
    ('loading_granularity', '-- granularity'),
    ('recovered_count', 'Recovered count'), ('recovered_count_grade', '-- grade'),
    ('upstream_doublet_handling', 'Upstream doublet handling'),
    ('upstream_evidence_grade', '-- grade'),
    ('donor_proportion_source', 'Donor proportions'),
    ('donor_proportion_grade', '-- grade'),
    ('donor_proportion_measured_at_stage', '-- measured at stage'),
    ('capture_mapping_coverage', 'Capture mapping'),
    ('runnability_classification', 'Classification'),
    ('evidence_source', 'Evidence source'),
    ('still_unidentifiable', 'Still unidentifiable'),
]


def render_packet_markdown(rows):
    """README plus the matrix rendered inline, so the PDF is self-contained --
    a reviewer must not need the CSV open to read the packet."""
    md = ['## Cohort evidence-gap matrix, in full', '',
          'Rendered from `cohort_evidence_gap_matrix.csv` so this document is self-contained.',
          'The CSV is attached for machine use; nothing below requires it.', '',
          # Raw LaTeX: the matrix is the densest part of the packet and the only
          # place where a single value can exceed a column. \small buys margin
          # without touching any value.
          r'\begingroup\small', '']
    for r in rows:
        # The separator lengths are load-bearing, not decoration. With `|:---|:---|`
        # pandoc emits non-wrapping `l` columns, so a long value overflows instead
        # of wrapping -- which is what put `not_eligible_without_a...` 65.5pt past
        # the right margin on pages 3-4. Long separators make pandoc compute
        # relative widths and emit wrapping `p{}` columns.
        md += [f"### {r['cohort'].capitalize()}", '', '| Field | Value |',
               '|:' + '-' * 26 + '|:' + '-' * 78 + '|']
        for key, label in MATRIX_FIELDS:
            val = (r.get(key) or '').strip() or '--'
            md.append(f"| {label} | {val.replace('|', chr(92) + '|')} |")
        md.append('')
    md += [r'\endgroup', '']
    src = open(os.path.join(OUT, 'README.md')).read()
    anchor = '## The five questions'
    if src.count(anchor) != 1:
        raise RuntimeError('README anchor for the matrix insert not found exactly once')
    # the arrow glyph is absent from the PDF font; it is the only substitution
    return src.replace(anchor, '\n'.join(md) + '\n' + anchor).replace('\u2192', '->')


def build_pdf():
    """Render the PDF, validate it, then atomically replace the tracked file.

    Runs the one PANDOC_ARGV above, so there is no second copy of the command to
    drift. Returns (sha256, n_bytes, n_pages).

    Ordering matters and is easy to get wrong: the render, the validation and
    the rename all sit INSIDE the try, so the finally that cleans up runs after
    all of them. On success os.replace has already moved PDF_TMP away, so the
    cleanup finds nothing; on any failure PDF_TMP still exists and is removed,
    because a failing pandoc can truncate or partially write its output before
    returning non-zero and a corrupt .inprogress.pdf must not be left beside the
    frozen artifact. The tracked PDF is never opened until the replace.
    """
    tmp_pdf = os.path.join(OUT, PDF_TMP)
    final = os.path.join(OUT, PDF_NAME)
    try:
        try:
            subprocess.run(PANDOC_ARGV, cwd=OUT, check=True, capture_output=True, text=True)
        except FileNotFoundError:
            raise SystemExit('pandoc not found; the PDF step needs pandoc and xelatex')
        except subprocess.CalledProcessError as exc:
            raise SystemExit('pandoc failed; the existing PDF is untouched:\n'
                             + exc.stderr[:2000])

        try:
            data = open(tmp_pdf, 'rb').read()
        except OSError as exc:
            raise SystemExit(f'pandoc reported success but wrote no PDF ({exc}); '
                             'the existing PDF is untouched')
        if not data.startswith(b'%PDF-') or len(data) < 10_000:
            raise SystemExit(f'rendered output is not a plausible PDF ({len(data)} bytes); '
                             'the existing PDF is untouched')
        os.replace(tmp_pdf, final)
    finally:
        for tmp in ('_packet_full.md', '_header.tex', PDF_TMP):
            path = os.path.join(OUT, tmp)
            if os.path.exists(path):
                os.remove(path)

    pages = None
    try:
        info = subprocess.run(['pdfinfo', final], capture_output=True, text=True).stdout
        pages = next((l.split()[-1] for l in info.splitlines() if l.startswith('Pages')), None)
    except FileNotFoundError:
        pass
    return hashlib.sha256(data).hexdigest(), len(data), pages


def tool_versions():
    out = {}
    for name, argv in [('pandoc', ['pandoc', '--version']),
                       ('xetex', ['xelatex', '--version'])]:
        try:
            out[name] = subprocess.run(argv, capture_output=True, text=True
                                       ).stdout.splitlines()[0].strip()
        except Exception:
            out[name] = 'not found'
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pdf', action='store_true',
                    help='also render the reviewer PDF (needs pandoc + xelatex)')
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    print('REVALIDATING packet claims against T7 artifacts')
    for c in revalidate():
        print(f'  [ok] {c}')

    jk = _read('capture_evidence_join_key.csv')
    status = collections.Counter((r['cohort'], r['mapping_status']) for r in jk)
    pair = _read('pooled-same-donor-pairing-scenarios-v3.csv')
    pair_stage = {}
    for r in pair:
        key = r['cohort'].split(' (')[0].lower()
        pair_stage.setdefault(key, r['proportions_measured_at_stage'])

    rows = []
    for c in COHORTS:
        f = FACTS[c]
        cov = {s: n for (coh, s), n in status.items() if coh == c}
        rows.append({
            'cohort': c, 'platform': f['platform'], 'capture_unit': f['capture_unit'],
            'input_stage': f['input_stage'],
            'chemistry': f['chemistry'], 'chemistry_evidence_grade': f['chemistry_grade'],
            'loading_evidence': f['loading'], 'loading_evidence_grade': f['loading_grade'],
            # Read from the explicit field. Inferring it from f['loading'] with a
            # substring test matched Ren's "... NOT per-capture" and filed Ren as
            # capture-addressable, contradicting Q2 in the same document.
            'loading_granularity': f['loading_granularity'],
            'recovered_count': f['recovered'], 'recovered_count_grade': f['recovered_grade'],
            'upstream_doublet_handling': f['upstream_doublet'],
            'upstream_evidence_grade': f['upstream_grade'],
            'donor_proportion_source': f['donor_proportions'],
            'donor_proportion_grade': f['donor_grade'],
            'donor_proportion_measured_at_stage': pair_stage.get(c, 'n/a'),
            'capture_mapping_coverage': '; '.join(f'{k}={v}' for k, v in sorted(cov.items())),
            'runnability_classification': f['classification'],
            'evidence_source': f['source'],
            'still_unidentifiable': f['unidentifiable'],
        })

    path = os.path.join(OUT, 'cohort_evidence_gap_matrix.csv')
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    print(f'\nwrote {len(rows)} cohort rows -> {path}')

    md_path = os.path.join(OUT, '_packet_full.md')
    with open(md_path, 'w') as fh:
        fh.write(render_packet_markdown(rows))
    print(f'wrote the self-contained markdown -> {md_path}')

    # Per-page footer. The packet circulates outside the project, so the three
    # things it is NOT must be on every page, not only in the introduction.
    hdr_path = os.path.join(OUT, '_header.tex')
    with open(hdr_path, 'w') as fh:
        fh.write(HEADER_TEX)
    print(f'wrote the page furniture      -> {hdr_path}')
    if not args.pdf:
        print('\nRe-run with --pdf to render the reviewer PDF (needs pandoc + xelatex).')
        print('Intermediates left in place for inspection; --pdf removes them after use.')
        return

    sha, nbytes, pages = build_pdf()
    versions = tool_versions()
    print(f'\nbuilt {PDF_NAME}')
    print(f'  sha256 : {sha}')
    print(f'  size   : {nbytes:,} bytes')
    print(f'  pages  : {pages or "unknown (pdfinfo absent)"}')
    for k, v in versions.items():
        print(f'  {k:<7}: {v}')
    print('\nRecord these in BUILD_PROVENANCE.md, with the commit that ships them.')
    print('NOTE: builds are NOT byte-identical between runs -- measured, not assumed.')
    print('The PDF is a frozen, checksummed artifact, not a reproducible-to-the-byte one.')


if __name__ == '__main__':
    main()
