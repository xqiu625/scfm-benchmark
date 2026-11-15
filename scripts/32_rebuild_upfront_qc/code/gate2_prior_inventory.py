"""
READ-ONLY provenance inventory for the Scrublet expected-rate (loading) prior.

Authorized scope (audit, 2026-09-21): inspect raw downloads, metadata, any
Cell Ranger / author summaries, and the existing builders. Do NOT select a
prior, modify calls, start Gate 2, or rebuild matrices.

Purpose. `loading = 0.008 x n_recovered / 1000` is defined over DROPLETS
RECOVERED FROM THE INSTRUMENT. Every prior currently in use counts something
else. This script asks, per cohort and per capture, which count levels are
actually obtainable and how each one is defined -- by measuring the artifacts,
not by reading builder docstrings (one of which was already found to be
factually wrong about Liu).

Columns are deliberately allowed to be NA. A level that cannot be obtained
must stay NA rather than be filled with the nearest available number; that
substitution is exactly the defect this inventory exists to expose.

Verdicts, per the audit:
  VALIDATED    definition matches the rate model, with source and checksum
  PROXY_ONLY   a related count exists but its definition does not match
  UNAVAILABLE  no defensible loading count

Output: gate2_prior_inventory.csv (per capture) + a per-cohort verdict summary.
"""
import glob
import gzip
import hashlib
import os

import h5py
import numpy as np
import pandas as pd

RAW = '/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw'
RES = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc'
OUT = f'{RES}/gate2_prior_inventory'

COLS = ['cohort', 'channel_id', 'platform',
        'n_input_barcodes', 'n_candidate_barcodes', 'n_called_pre_qc',
        'n_post_demux_singlets', 'n_unresolved_or_amb',
        'n_cross_donor_doublets_removed', 'n_gate1_retained',
        'source_artifact', 'source_stage', 'source_checksum',
        'count_definition', 'loading_prior_eligibility', 'reason_not_eligible']


def cheap_checksum(path, nbytes=1 << 20):
    """Size + head/tail digest. A full hash of hundreds of GB is not worth the
    IO here; this is enough to detect the artifact being swapped."""
    try:
        sz = os.path.getsize(path)
        h = hashlib.sha256(str(sz).encode())
        with open(path, 'rb') as fh:
            h.update(fh.read(nbytes))
            if sz > nbytes:
                fh.seek(max(0, sz - nbytes))
                h.update(fh.read(nbytes))
        return f'{sz}:{h.hexdigest()[:16]}'
    except OSError:
        return pd.NA


def gz_lines(path):
    with gzip.open(path, 'rt') as fh:
        return sum(1 for _ in fh)


def gate1_counts(cohort):
    """n_gate1_retained per capture, from whichever Gate 1 artifact exists."""
    f = f'{RES}/gate1_qc/{cohort}/cells.csv'
    if os.path.exists(f):
        c = pd.read_csv(f, usecols=['channel_id'], low_memory=False)
        return c.channel_id.astype(str).value_counts().to_dict()
    split = f'{RES}/gate1_qc/{cohort}/split_per_channel'
    if os.path.isdir(split):
        out = {}
        for ch in sorted(os.listdir(split)):
            p = f'{split}/{ch}/cells.csv'
            if os.path.exists(p):
                out[str(ch)] = sum(1 for _ in open(p)) - 1
        return out
    return {}


def rows_from_geo_triplets(cohort, pattern, platform, stage, defn):
    """Qiu / Chongqing: GEO barcodes/features/matrix per GSM."""
    g1 = gate1_counts(cohort)
    rows = []
    for bc in sorted(glob.glob(pattern)):
        base = os.path.basename(bc)
        # Gate 1 keys channels by the bare GSM accession; the filename carries
        # an extra sample suffix, so strip it or the join silently yields nothing
        ch = base.split('_barcodes')[0].split('.barcodes')[0].split('_')[0]
        n = gz_lines(bc)
        rows.append({
            'cohort': cohort, 'channel_id': ch, 'platform': platform,
            'n_input_barcodes': pd.NA, 'n_candidate_barcodes': pd.NA,
            'n_called_pre_qc': n, 'n_post_demux_singlets': pd.NA,
            'n_unresolved_or_amb': pd.NA, 'n_cross_donor_doublets_removed': pd.NA,
            'n_gate1_retained': g1.get(ch, pd.NA),
            'source_artifact': bc, 'source_stage': stage,
            'source_checksum': cheap_checksum(bc), 'count_definition': defn,
            'loading_prior_eligibility': 'UNAVAILABLE',
            'reason_not_eligible': 'depositor supplied a filtered/called cell matrix; '
                                    'no unfiltered droplet matrix and no Cell Ranger '
                                    'metrics_summary in the deposit'})
    return rows


def rows_liu():
    """Liu: per-lane filtered_feature_bc_matrix.h5 -- verified to be a fixed
    30,000-barcode demuxlet whitelist, not a Cell Ranger cell call."""
    g1 = gate1_counts('liu')
    demux = f'{RAW}/GSE161918_liu/singlets/liu_demuxlet_filter_summary.csv'
    dm = pd.read_csv(demux) if os.path.exists(demux) else None
    rows = []
    for h5 in sorted(glob.glob(f'{RAW}/GSE161918_liu/rna_h5/*.h5')):
        base = os.path.basename(h5)
        ch = base.split('_RNA_')[0].split('_', 1)[1] if '_' in base else base
        try:
            with h5py.File(h5, 'r') as f:
                n = f[list(f.keys())[0]]['barcodes'].shape[0]
        except Exception:
            n = pd.NA
        sng = amb = pd.NA
        if dm is not None:
            m = dm[dm.astype(str).apply(lambda r: ch in ' '.join(r.values), axis=1)]
            if len(m) == 1:
                for c in m.columns:
                    lc = c.lower()
                    if 'sng' in lc or 'singlet' in lc:
                        sng = m.iloc[0][c]
                    if 'amb' in lc or 'dbl' in lc or 'doublet' in lc:
                        amb = m.iloc[0][c]
        rows.append({
            'cohort': 'liu', 'channel_id': ch, 'platform': '10x',
            'n_input_barcodes': pd.NA, 'n_candidate_barcodes': n,
            'n_called_pre_qc': pd.NA, 'n_post_demux_singlets': sng,
            'n_unresolved_or_amb': amb, 'n_cross_donor_doublets_removed': pd.NA,
            'n_gate1_retained': g1.get(ch, pd.NA),
            'source_artifact': h5, 'source_stage': 'demuxlet whitelist h5',
            'source_checksum': cheap_checksum(h5),
            'count_definition': 'fixed 30,000-barcode per-lane demuxlet whitelist '
                                 '(verified identical across all 55 lanes) -- neither a '
                                 'droplet count nor a cell call',
            'loading_prior_eligibility': 'UNAVAILABLE',
            'reason_not_eligible': 'the only per-lane count is a fixed whitelist size, '
                                    'identical for every lane, so it carries no '
                                    'capture-specific loading information'})
    return rows


def rows_from_h5ad(cohort, path, obs_channel_key, platform, note):
    """COMBAT / Stephenson: one author-processed h5ad. layers['raw'] means raw
    COUNT VALUES, not raw droplets -- the cell set is already the author's."""
    g1 = gate1_counts(cohort)
    rows = []
    per_ch = {}
    try:
        with h5py.File(path, 'r') as f:
            o = f['obs'][obs_channel_key]
            if isinstance(o, h5py.Group):
                codes = o['codes'][:]
                cats = [c.decode() if isinstance(c, bytes) else str(c)
                        for c in o['categories'][:]]
                vals = np.asarray(cats, dtype=object)[codes]
            else:
                vals = np.asarray([v.decode() if isinstance(v, bytes) else str(v)
                                   for v in o[:]], dtype=object)
        u, cnt = np.unique(vals, return_counts=True)
        per_ch = dict(zip(u.tolist(), cnt.tolist()))
    except Exception as e:
        per_ch = {}
        print(f'  {cohort}: could not read obs/{obs_channel_key} ({e})')
    for ch, n in sorted(per_ch.items()):
        rows.append({
            'cohort': cohort, 'channel_id': str(ch), 'platform': platform,
            'n_input_barcodes': pd.NA, 'n_candidate_barcodes': pd.NA,
            'n_called_pre_qc': n, 'n_post_demux_singlets': pd.NA,
            'n_unresolved_or_amb': pd.NA, 'n_cross_donor_doublets_removed': pd.NA,
            'n_gate1_retained': g1.get(str(ch), pd.NA),
            'source_artifact': path, 'source_stage': 'author-processed h5ad',
            'source_checksum': cheap_checksum(path), 'count_definition': note,
            'loading_prior_eligibility': 'UNAVAILABLE',
            'reason_not_eligible': "layers['raw'] holds raw count VALUES, not raw "
                                    'droplets; the cell set is the depositor\'s own '
                                    'post-QC selection and its filtering is not recorded'})
    return rows


def rows_wilk():
    """Wilk: per-GSM Seq-Well cell.counts.matrices.rds."""
    g1 = gate1_counts('wilk')
    rows = []
    for rds in sorted(glob.glob(f'{RAW}/GSE174072_wilk/*cell.counts.matrices.rds.gz')):
        ch = os.path.basename(rds).split('_cell.counts')[0]
        rows.append({
            'cohort': 'wilk', 'channel_id': ch, 'platform': 'seqwell',
            'n_input_barcodes': pd.NA, 'n_candidate_barcodes': pd.NA,
            'n_called_pre_qc': pd.NA, 'n_post_demux_singlets': pd.NA,
            'n_unresolved_or_amb': pd.NA, 'n_cross_donor_doublets_removed': pd.NA,
            'n_gate1_retained': g1.get(ch, pd.NA),
            'source_artifact': rds, 'source_stage': 'author-processed Seq-Well RDS',
            'source_checksum': cheap_checksum(rds),
            'count_definition': 'Seq-Well cell x gene counts matrix as deposited',
            'loading_prior_eligibility': 'UNAVAILABLE',
            'reason_not_eligible': 'Seq-Well is not droplet loading -- the 10x '
                                    '0.8%/1,000 multiplet model does not apply, and no '
                                    'Seq-Well loading record exists for this cohort'})
    return rows


def rows_cillo():
    """Cillo: per-GSM GEX tar.gz holding a filtered barcodes/features/matrix
    triplet. Counted from the tar index so nothing is extracted."""
    import tarfile
    g1 = gate1_counts('cillo')
    rows = []
    for tp in sorted(glob.glob(f'{RAW}/GSE180578_cillo/*_GEX.tar.gz')):
        ch = os.path.basename(tp).split('_')[0]
        n = pd.NA
        try:
            with tarfile.open(tp, 'r:gz') as tf:
                for m in tf:
                    if m.name.endswith('barcodes.tsv'):
                        n = sum(1 for _ in tf.extractfile(m))
                        break
        except Exception:
            pass
        rows.append({
            'cohort': 'cillo', 'channel_id': ch, 'platform': '10x',
            'n_input_barcodes': pd.NA, 'n_candidate_barcodes': pd.NA,
            'n_called_pre_qc': n, 'n_post_demux_singlets': pd.NA,
            'n_unresolved_or_amb': pd.NA, 'n_cross_donor_doublets_removed': pd.NA,
            'n_gate1_retained': g1.get(ch, pd.NA),
            'source_artifact': tp, 'source_stage': 'GEO filtered matrix in GEX tar',
            'source_checksum': cheap_checksum(tp),
            'count_definition': 'barcodes.tsv of the deposited filtered matrix',
            'loading_prior_eligibility': 'UNAVAILABLE',
            'reason_not_eligible': 'depositor supplied a filtered/called cell matrix; '
                                    'no unfiltered droplet matrix and no Cell Ranger '
                                    'metrics_summary in the deposit'})
    return rows


def rows_ren():
    """Ren: source lives outside this project, under 302005/datasets. One
    pooled annotation file; every barcode already carries a celltype, i.e.
    these are the authors' called cells."""
    ann = ('/bigdata/godziklab/shared/Xinru/302005/datasets/302004data15/'
           'GSE158055_cell_annotation.csv.gz')
    g1 = gate1_counts('ren')
    rows = []
    if not os.path.exists(ann):
        return rows
    per = pd.read_csv(ann, usecols=['sampleID']).sampleID.astype(str).value_counts()
    ck = cheap_checksum(ann)
    for ch, n in per.sort_index().items():
        rows.append({
            'cohort': 'ren', 'channel_id': ch, 'platform': '10x',
            'n_input_barcodes': pd.NA, 'n_candidate_barcodes': pd.NA,
            'n_called_pre_qc': int(n), 'n_post_demux_singlets': pd.NA,
            'n_unresolved_or_amb': pd.NA, 'n_cross_donor_doublets_removed': pd.NA,
            'n_gate1_retained': g1.get(ch, pd.NA),
            'source_artifact': ann, 'source_stage': 'author cell annotation (pooled)',
            'source_checksum': ck,
            'count_definition': "barcodes in the authors' cell_annotation.csv; every "
                                 'barcode carries a celltype, so this is a called-cell set',
            'loading_prior_eligibility': 'UNAVAILABLE',
            'reason_not_eligible': 'only a called-cell annotation was deposited; no '
                                    'unfiltered droplet matrix and no per-library '
                                    'recovered-cell record'})
    return rows


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = []
    rows += rows_from_geo_triplets(
        'qiu', f'{RAW}/GSE167363_qiu/extracted/*barcodes.tsv.gz', '10x',
        'GEO filtered matrix', 'barcodes.tsv of the deposited filtered matrix')
    rows += rows_from_geo_triplets(
        'chongqing', f'{RAW}/chongqing_sepsis_GSE279452/extracted_full/*barcodes.tsv.gz',
        '10x', 'GEO filtered matrix', 'barcodes.tsv of the deposited filtered matrix')
    rows += rows_liu()
    rows += rows_wilk()
    rows += rows_from_h5ad('combat', f'{RAW}/COMBAT/COMBAT-CITESeq-DATA.h5ad',
                            'Channel_ID', '10x',
                            'cells in the depositor-assembled CITE-seq object')
    rows += rows_from_h5ad('stephenson',
                            f'{RAW}/Stephenson_EMTAB10026/covid_portal_210320_with_raw.h5ad',
                            'patient_id', '10x',
                            'cells in the depositor-assembled portal object')
    rows += rows_cillo()
    rows += rows_ren()

    df = pd.DataFrame(rows, columns=COLS)
    df.to_csv(f'{OUT}/gate2_prior_inventory.csv', index=False)

    print(f'\nwrote {len(df)} capture rows -> {OUT}/gate2_prior_inventory.csv\n')
    print('Per-cohort verdicts:')
    s = (df.groupby(['cohort', 'loading_prior_eligibility'])
           .size().rename('captures').reset_index())
    print(s.to_string(index=False))
    print('\nWhich count levels are obtainable (non-NA captures):')
    lv = ['n_input_barcodes', 'n_candidate_barcodes', 'n_called_pre_qc',
          'n_post_demux_singlets', 'n_unresolved_or_amb', 'n_gate1_retained']
    print(df.groupby('cohort')[lv].count().to_string())


if __name__ == '__main__':
    main()
