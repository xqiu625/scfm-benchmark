"""Assemble the two outbound review packages from VERSIONED sources.

    python scripts/32_rebuild_upfront_qc/code/assemble_review_delivery.py

Delivery text used to exist only in the Desktop staging directory, so a commit
could not reconstruct or prove what a reviewer had actually been sent. Every
document that defines the review task, the AI-opinion boundary or the attachment
identity now lives in the repository, and this script is the only thing that
writes the staging directories.

It COPIES and CHECKSUMS. It renders no PDF, recomputes no matrix, and changes no
value. The Gate 2 PDF is rebuilt by `build_sensitivity_review_packet.py --pdf`;
the Stage 31 tables are frozen artifacts from 2026-09-22 that are not regenerated
here and never modified.

Nothing here sends anything. The handover logs stay unsent until a real send is
recorded.
"""
import hashlib
import os
import shutil

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    '..', '..', '..'))
STAGING = '/Users/xinruqiu/Desktop/302006/review-materials'
SD = os.path.join(REPO, 'reports/manuscript/revision/sensitivity_design_review')
S31 = os.path.join(REPO, 'reports/manuscript/revision/stage31_dimension_review')
T7 = os.path.join(REPO, 'reports/manuscript/revision/literature_audit/T7_doublet_prior_evidence')

GATE2_OUT = os.path.join(STAGING, 'gate2_packet_2026-09-24')
S31_OUT = os.path.join(STAGING, 'stage31_dimension_2026-09-22')

#: Attachments are NOT duplicated into the repository -- they already live in the
#: T7 evidence directory. The manifest records where each came from so the sent
#: copy is traceable to a repository path, which is what the checksum is for.
GATE2_ATTACHMENTS = [
    (SD, 'cohort_evidence_gap_matrix.csv'),
    (T7, 'doublet_prior_capture_evidence.csv'),
    (T7, 'capture_evidence_join_key.csv'),
    (T7, 'doublet_prior_source_audit.csv'),
    (T7, 'ren_gse158055_capture_qc.csv'),
    (T7, 'pooled-same-donor-pairing-scenarios-v3.csv'),
    (T7, 'liu_demuxlet_per_lane.csv'),
    (T7, 'tenx_multiplet_rate_tables.csv'),
]

#: The eight Stage 31 content files. The seven frozen ones are NOT in the
#: repository and are not regenerated: they are the 2026-09-22 artifacts already
#: staged, and this script must not touch them. Only the addendum and the cover
#: come from version control.
S31_FROZEN = ['README.md', 'RUN_MANIFEST.json', 'scrublet_dimension_parameter.md',
              'dimension_audit.csv', 'dimension_coverage_table.csv',
              'dimension_excluded_samples_long.csv', 'build_coverage_table.py']


def sha256(path):
    with open(path, 'rb') as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def copy_from_repo(src_dir, name, dest_dir, dest_name=None):
    src = os.path.join(src_dir, name)
    if not os.path.exists(src):
        raise SystemExit(f'versioned source missing: {src}')
    dst = os.path.join(dest_dir, dest_name or name)
    shutil.copyfile(src, dst)
    return os.path.relpath(src, REPO), sha256(dst)


def main():
    lines = []

    # ---- Gate 2 --------------------------------------------------------
    os.makedirs(os.path.join(GATE2_OUT, 'attachments'), exist_ok=True)
    provenance = []
    for name, dest in (('sensitivity_design_review_packet.pdf', None),
                       ('BUILD_PROVENANCE.md', None),
                       ('DELIVERY_README_2026-09-24.md', 'README.md')):
        rel, digest = copy_from_repo(SD, name, GATE2_OUT, dest)
        provenance.append((dest or name, rel, digest))
    for src_dir, name in GATE2_ATTACHMENTS:
        rel, digest = copy_from_repo(src_dir, name,
                                     os.path.join(GATE2_OUT, 'attachments'))
        provenance.append((f'attachments/{name}', rel, digest))

    with open(os.path.join(GATE2_OUT, 'SHA256SUMS.txt'), 'w') as fh:
        for name, _rel, digest in provenance:
            fh.write(f'{digest}  {name}\n')

    with open(os.path.join(GATE2_OUT, 'ATTACHMENT_MANIFEST.md'), 'w') as fh:
        fh.write('# Gate 2 package — what was sent, and where each file came from\n\n')
        fh.write('Assembled by `scripts/32_rebuild_upfront_qc/code/assemble_review_delivery.py`\n')
        fh.write('from versioned sources. Verify with `shasum -a 256 -c SHA256SUMS.txt`.\n\n')
        fh.write('| file in the package | repository source | SHA-256 |\n|---|---|---|\n')
        for name, rel, digest in provenance:
            fh.write(f'| `{name}` | `{rel}` | `{digest}` |\n')
        fh.write('\n**Handover status: not yet sent.**\n')
    print(f'Gate 2  -> {GATE2_OUT}  ({len(provenance)} files)')

    # ---- Stage 31 ------------------------------------------------------
    if not os.path.isdir(S31_OUT):
        raise SystemExit(f'the frozen Stage 31 bundle is not staged at {S31_OUT}')
    missing = [f for f in S31_FROZEN if not os.path.exists(os.path.join(S31_OUT, f))]
    if missing:
        raise SystemExit(f'frozen Stage 31 files absent: {missing}')

    s31 = []
    rel, digest = copy_from_repo(S31, 'ADDENDUM_2026-09-24.md', S31_OUT)
    s31.append(('ADDENDUM_2026-09-24.md', rel, digest))
    copy_from_repo(S31, 'READ_FIRST_2026-09-24.md', S31_OUT)
    for f in S31_FROZEN:                    # frozen: checksummed, never rewritten
        s31.append((f, '(frozen 2026-09-22 artifact, not in version control)',
                    sha256(os.path.join(S31_OUT, f))))

    with open(os.path.join(S31_OUT, 'SHA256SUMS.txt'), 'w') as fh:
        for name, _rel, digest in s31:
            fh.write(f'{digest}  {name}\n')
    print(f'Stage 31 -> {S31_OUT}  ({len(s31)} content files + cover)')
    print('\nBoth packages assembled. Neither has been sent; handover logs remain open.')
    for name, _rel, digest in s31:
        lines.append(name)


if __name__ == '__main__':
    main()
