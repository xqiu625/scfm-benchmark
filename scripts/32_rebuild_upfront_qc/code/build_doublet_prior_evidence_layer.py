"""
Build the two derived tables of the doublet-prior EVIDENCE LAYER.

This is an evidence layer, not a configuration layer. Nothing here selects a
rate, and no Gate 2 runtime code may read these files. The boundary:

    source audit      -> where each claim came from
    capture evidence  -> mixed-level claims, each tagged with its level
    JOIN KEY          <- this script; identity mapping onto the 746 captures
    VERDICTS          <- this script; three orthogonal status fields
    -------- contract boundary; nothing above crosses it --------
    prior config      -> only values that pass docs/data-contracts.md
    sensitivity grid  -> rate grid, source, chemistry, model, limits

Why a join-key table is needed at all: the evidence table's `capture_id` is
NOT a single key. Ren/COMBAT/Liu match Gate 2 `channel_id` directly; Qiu,
Cillo and Chongqing match only on `repo_sample_id` (their `capture_id` is the
author's internal label, e.g. Qiu 'HC1' vs GEO 'GSM5102900'). Joining on
`capture_id` alone silently resolves 409 of 746 captures and drops the rest,
so the resolved key is materialized here rather than re-derived per consumer.

Two kinds of non-measurement are kept explicitly distinct, because collapsing
them is how a cohort-level protocol statement turns into a fake per-capture
measurement:
  cohort_protocol_inheritance  the protocol figure is uniform and genuinely
                               applies to every capture (Wilk: 15,000 cells
                               per array), inherited with capture_specific=false
  unmapped_to_gem_well         the deposit cannot be resolved to capture units
                               at all (Stephenson: 98 sequence pools, GEM wells
                               not identifiable)

Read-only with respect to HPCC: reads the committed crosswalk copy, runs no
job, touches no Gate 2 artifact.
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
T7 = os.path.join(HERE, '..', '..', '..', 'reports', 'manuscript', 'revision',
                  'literature_audit', 'T7_doublet_prior_evidence')
T7 = os.path.normpath(T7)

EVIDENCE = os.path.join(T7, 'doublet_prior_capture_evidence.csv')
CROSSWALK = os.path.join(T7, 'gate2_source_crosswalk.csv')

# Evidence-table cohort label -> repository cohort key
COHORT = {'Ren': 'ren', 'Cillo': 'cillo', 'COMBAT': 'combat', 'Liu': 'liu',
          'Qiu': 'qiu', 'Chongqing': 'chongqing', 'Stephenson': 'stephenson',
          'Wilk': 'wilk'}

# Which evidence column actually resolves to Gate 2 channel_id, per cohort.
# Verified 2026-09-21: each listed key joins 100% of that cohort's captures.
JOIN_ON = {'ren': 'capture_id', 'combat': 'capture_id', 'liu': 'capture_id',
           'qiu': 'repo_sample_id', 'cillo': 'repo_sample_id',
           'chongqing': 'repo_sample_id'}

# Cohorts with no capture-level evidence row, and why.
NO_CAPTURE_EVIDENCE = {
    'wilk': ('cohort_protocol_inheritance', 'cohort_protocol', 'false', 'inherited',
             '15,000 cells per Seq-Well array, uniform across the cohort, stated in the '
             'paper and in GEO -- a protocol figure that genuinely applies to every '
             'capture, NOT a per-capture measurement'),
    'stephenson': ('none', 'none', 'false', 'unmapped_to_gem_well',
                   'the deposit resolves to 98 sequence pools; GEM wells are not '
                   'identifiable, so no loading figure can be attached to a capture'),
}

# Three orthogonal status fields. A single status cannot express a cohort whose
# loading is known but whose upstream filtering leaves the residual unknown.
VERDICTS = [
    # cohort, loading_evidence_status, upstream_doublet_handling, rate_model_status, note
    ('ren', 'PARTIAL_OR_UNMAPPABLE_LOADING', 'AUTHOR_SCRUBLET_PLUS_MANUAL_FILTER',
     'RESIDUAL_RATE_UNKNOWN',
     "Loading is reported cohort-level only (300-14,000 target cells per reaction) and "
     "cannot be attached per capture. Our Gate 2 input is demonstrably POST author "
     "Scrublet and post manual removal: summing Ncells_manual_removal over the 284 "
     "samples gives 1,462,702, matching our per-capture called-cell counts sample by "
     "sample. The authors' expected_doublet_rate = 0.08 was their INPUT prior, not a "
     "measured removed fraction, so it must not be reused as our residual prior -- and "
     "the residual prevalence after their filtering is unknown."),
    ('wilk', 'LOADING_KNOWN', 'AUTHOR_HEURISTIC_MULTIPLET_FILTER',
     'PLATFORM_MODEL_APPLICABLE_WITH_UNCERTAINTY',
     "15,000 cells per array, uniform and stated. A Seq-Well loading-density model exists "
     "(Hughes 2020), so this is not 'no model available'; the uncertainty is cross-protocol "
     "transfer plus an unknown residual after the authors' heuristic filter (cells with "
     ">75 genes per 100 UMIs removed as putative multiplets)."),
    ('combat', 'LOADING_KNOWN', 'GENOTYPE_DEMUX',
     'PLATFORM_MODEL_EXTRAPOLATION_REQUIRED',
     "50,000 cells loaded per channel, 14 donors pooled. That load is beyond the range of "
     "the official 10x multiplet tables, so applying the vendor model requires explicit "
     "extrapolation -- it may not be extended linearly without stating so. Cross-donor "
     "doublets were removed by demuxlet v2 + Vireo; Scrublet score was retained only as a "
     "QC covariate."),
    ('qiu', 'LOADING_KNOWN', 'NONE_REPORTED',
     'PLATFORM_MODEL_APPLICABLE_WITH_UNCERTAINTY',
     "15,000 cells loaded per sample, protocol-level nominal loading, single donor per "
     "capture. No author doublet removal reported, so no upstream depletion to account "
     "for. Uncertainty is chemistry: GEO records 3' v3.1 while the paper text describes "
     "5' UMI counting -- unresolved, and the two have different multiplet tables."),
    ('stephenson', 'PARTIAL_OR_UNMAPPABLE_LOADING', 'MIXED_OR_SITE_SPECIFIC',
     'RESIDUAL_RATE_UNKNOWN',
     "Nominal loading is reported (20,000-30,000 per sample at Newcastle/London; 2 lanes x "
     "25,000 per pool of four at Cambridge) but cannot be mapped to GEM wells. Upstream "
     "handling differs by site: Scrublet per sample at Newcastle, emptyDrops plus genotype "
     "doublet tagging at Cambridge. Chemistry is also contested (paper 5' NextGEM v1.1 vs "
     "SDRF '10xV3')."),
    ('liu', 'LOADING_UNAVAILABLE', 'GENOTYPE_DEMUX', 'RESIDUAL_RATE_UNKNOWN',
     "No loading figure reported and no raw droplet matrix deposited; the per-lane count is "
     "a fixed 30,000-barcode whitelist, identical across all 55 lanes, carrying no "
     "capture-specific loading information. demuxlet v2 removed cross-donor doublets, so "
     "within-donor residual is unknown. Per-lane SNG/DBL/AMB are preserved in "
     "liu_demuxlet_per_lane.csv -- they classify the whitelist, which includes empty "
     "droplets, and are therefore NOT a doublet rate among real cells."),
    ('cillo', 'LOADING_UNAVAILABLE', 'HTO_FILTER', 'RESIDUAL_RATE_UNKNOWN',
     "Target cells loaded not reported. HTO double-positive exclusion removed cross-sample "
     "doublets only; no within-sample doublet call was made, so within-sample residual is "
     "unknown."),
    ('chongqing', 'LOADING_UNAVAILABLE', 'NONE_REPORTED', 'MODEL_UNAVAILABLE',
     "No linked publication and no loading figure in the GEO extract protocols. No upstream "
     "doublet handling reported. With neither a loading input nor upstream filtering to "
     "reason from, no rate model can be applied at all. NOTE: this assignment was made "
     "here rather than specified in review -- confirm before it is relied on."),
]

VALID = {
    'loading_evidence_status': {'LOADING_KNOWN', 'PARTIAL_OR_UNMAPPABLE_LOADING',
                                 'LOADING_UNAVAILABLE'},
    'upstream_doublet_handling': {'NONE_REPORTED', 'GENOTYPE_DEMUX', 'HTO_FILTER',
                                   'AUTHOR_SCRUBLET', 'AUTHOR_SCRUBLET_PLUS_MANUAL_FILTER',
                                   'AUTHOR_HEURISTIC_MULTIPLET_FILTER',
                                   'MIXED_OR_SITE_SPECIFIC'},
    'rate_model_status': {'DIRECT_RESIDUAL_ESTIMATE',
                           'PLATFORM_MODEL_APPLICABLE_WITH_UNCERTAINTY',
                           'PLATFORM_MODEL_EXTRAPOLATION_REQUIRED',
                           'RESIDUAL_RATE_UNKNOWN', 'MODEL_UNAVAILABLE'},
}


def load(path):
    with open(path, newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))


def main():
    ev = load(EVIDENCE)
    cw = load(CROSSWALK)

    gate2 = {}
    for r in cw:
        if r['in_gate2'] == 'True':
            gate2.setdefault(r['cohort'], []).append(r['channel_id'])

    # index evidence rows by (cohort, resolved key)
    idx = {}
    for r in ev:
        c = COHORT.get(r['cohort'])
        if c in JOIN_ON:
            idx[(c, r[JOIN_ON[c]])] = r

    rows = []
    for cohort in sorted(gate2):
        for ch in sorted(gate2[cohort]):
            if cohort in JOIN_ON:
                e = idx.get((cohort, ch))
                rows.append({
                    'cohort': cohort, 'gate2_channel_id': ch,
                    'evidence_capture_id': e['capture_id'] if e else '',
                    'evidence_repo_sample_id': e['repo_sample_id'] if e else '',
                    'join_column': JOIN_ON[cohort],
                    'mapping_method': 'capture_level_join' if e else 'none',
                    'evidence_granularity': 'capture' if e else 'none',
                    'capture_specific': 'true' if e else 'false',
                    'mapping_status': 'mapped' if e else 'unmatched',
                    'mapping_note': '' if e else 'no evidence row resolved for this capture'})
            else:
                method, gran, spec, status, note = NO_CAPTURE_EVIDENCE[cohort]
                rows.append({
                    'cohort': cohort, 'gate2_channel_id': ch,
                    'evidence_capture_id': '', 'evidence_repo_sample_id': '',
                    'join_column': '', 'mapping_method': method,
                    'evidence_granularity': gran, 'capture_specific': spec,
                    'mapping_status': status, 'mapping_note': note})

    out = os.path.join(T7, 'capture_evidence_join_key.csv')
    with open(out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        w.writeheader()
        w.writerows(rows)

    vout = os.path.join(T7, 'prior_evidence_verdicts.csv')
    with open(vout, 'w', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow(['cohort', 'gate2_captures', 'loading_evidence_status',
                    'upstream_doublet_handling', 'rate_model_status', 'note'])
        for c, le, ud, rm, note in sorted(VERDICTS):
            assert le in VALID['loading_evidence_status'], le
            assert ud in VALID['upstream_doublet_handling'], ud
            assert rm in VALID['rate_model_status'], rm
            w.writerow([c, len(gate2[c]), le, ud, rm, ' '.join(note.split())])

    # ---- validation
    n = len(rows)
    uniq = len({(r['cohort'], r['gate2_channel_id']) for r in rows})
    by_status = {}
    for r in rows:
        by_status[r['mapping_status']] = by_status.get(r['mapping_status'], 0) + 1
    total_gate2 = sum(len(v) for v in gate2.values())

    print(f'join key rows         : {n}')
    print(f'unique (cohort,chan)  : {uniq}')
    print(f'Gate 2 captures       : {total_gate2}')
    print(f'mapping_status        : {by_status}')
    assert n == uniq == total_gate2 == 746, 'capture coverage broken'
    assert by_status.get('mapped') == 580, by_status
    assert by_status.get('inherited') == 23, by_status
    assert by_status.get('unmapped_to_gem_well') == 143, by_status
    assert not any(v[3] == 'DIRECT_RESIDUAL_ESTIMATE' for v in VERDICTS), \
        'no cohort currently supports a direct residual estimate'
    print('\nall assertions passed')
    print(f'wrote {out}')
    print(f'wrote {vout}')


if __name__ == '__main__':
    main()
