"""
READ-ONLY crosswalk: source artifacts -> Gate 2 captures.

Supersedes gate2_prior_inventory.py, which reported per-cohort verdicts but
left three join defects that made its NA/zero columns unreadable:

  1. Stephenson's obs/sample_id is an int16 dataset whose labels live in the
     legacy obs/__categories/ group (AnnData encoding 0.1.0). Reading it as a
     modern codes/categories Group yielded raw integers, so 130 bogus
     "captures" appeared instead of the real 143.
  2. Qiu/Chongqing channel_ids carried the deposit's sample suffix while
     Gate 1 keys on the bare GSM accession, so every row failed to join.
  3. Wilk's n_called_pre_qc was never read, because the deposit ships RDS.
     The converted/<channel>/barcodes.tsv produced by convert_wilk_rds_to_mtx.R
     holds the same count without needing R.

Each cohort's channel_id definition is taken from its builder, which is the
authoritative producer:
  qiu/chongqing  bare GSM accession        build_*_decontx_input.py
  cillo          GSM accession             build_cillo_decontx_input.py:220
  stephenson     obs['sample_id']          build_stephenson_decontx_input.py:84
  combat         obs['Channel_ID']         build_combat_decontx_input.py:92
  wilk           converted dir basename    build_wilk_decontx_input.py:131
  liu            lane id from the h5 name  build_liu_decontx_input.py
  ren            annotation sampleID       build_ren_decontx_input.py

The join is reported in both directions. A source artifact with no Gate 2
capture and a Gate 2 capture with no source artifact are both defects worth
seeing, so neither is silently dropped.

Writes gate2_source_crosswalk.csv and a reconciliation summary. Selects no
prior, modifies no calls, starts no Gate 2 run.
"""
import glob
import gzip
import hashlib
import os
import tarfile

import h5py
import numpy as np
import pandas as pd

RAW = '/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw'
REN_SRC = ('/bigdata/godziklab/shared/Xinru/302005/datasets/302004data15/'
           'GSE158055_cell_annotation.csv.gz')
RES = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc'
OUT = f'{RES}/gate2_prior_inventory'

COLS = ['cohort', 'channel_id', 'platform', 'in_gate2', 'match',
        'n_input_barcodes', 'n_candidate_barcodes', 'n_called_pre_qc',
        'n_post_demux_singlets', 'n_unresolved_or_amb',
        'n_cross_donor_doublets_removed', 'n_gate1_retained',
        'source_artifact', 'source_stage', 'source_checksum',
        'obs_field', 'producer_script', 'count_definition',
        'loading_prior_eligibility', 'reason_not_eligible']

NOT_ELIGIBLE_FILTERED = ('depositor supplied a filtered/called cell matrix; no unfiltered '
                          'droplet matrix and no Cell Ranger metrics_summary in the deposit')


def cheap_checksum(path, nbytes=1 << 20):
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


def read_obs_categorical(f, key):
    """Handle both AnnData encodings. The legacy 0.1.0 form stores codes as a
    plain integer dataset with labels in obs/__categories/<key>; reading that
    as a modern Group silently returns the integer codes instead."""
    o = f['obs'][key]
    if isinstance(o, h5py.Group):
        cats = [c.decode() if isinstance(c, bytes) else str(c) for c in o['categories'][:]]
        return np.asarray(cats, dtype=object)[o['codes'][:]]
    if '__categories' in f['obs'] and key in f['obs']['__categories']:
        cats = [c.decode() if isinstance(c, bytes) else str(c)
                for c in f['obs']['__categories'][key][:]]
        return np.asarray(cats, dtype=object)[o[:]]
    return np.asarray([v.decode() if isinstance(v, bytes) else str(v) for v in o[:]],
                      dtype=object)


def gate2_captures(cohort):
    """Authoritative Gate 2 capture set. Liu and Wilk run frozen cohort-specific
    rules whose output lives in gate2_input/, NOT gate1_qc/."""
    for p in (f'{RES}/gate2_input/{cohort}/cells.csv', f'{RES}/gate1_qc/{cohort}/cells.csv'):
        if os.path.exists(p):
            c = pd.read_csv(p, usecols=['channel_id'], low_memory=False)
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


def row(cohort, ch, platform, producer, stage, artifact, defn, obs_field=pd.NA,
        reason=NOT_ELIGIBLE_FILTERED, **counts):
    r = {c: pd.NA for c in COLS}
    r.update(cohort=cohort, channel_id=str(ch), platform=platform,
             source_artifact=artifact, source_stage=stage,
             source_checksum=cheap_checksum(artifact) if artifact else pd.NA,
             obs_field=obs_field, producer_script=producer, count_definition=defn,
             loading_prior_eligibility='UNAVAILABLE', reason_not_eligible=reason)
    r.update(counts)
    return r


def collect():
    rows = []

    # --- Qiu, Chongqing: GEO barcodes/features/matrix triplets, bare GSM key
    for cohort, pat, prod in [
            ('qiu', f'{RAW}/GSE167363_qiu/extracted/*barcodes.tsv.gz',
             'build_qiu_decontx_input.py'),
            ('chongqing', f'{RAW}/chongqing_sepsis_GSE279452/extracted_full/*barcodes.tsv.gz',
             'build_chongqing_decontx_input.py')]:
        for bc in sorted(glob.glob(pat)):
            base = os.path.basename(bc)
            ch = base.split('_barcodes')[0].split('.barcodes')[0].split('_')[0]
            with gzip.open(bc, 'rt') as fh:
                n = sum(1 for _ in fh)
            rows.append(row(cohort, ch, '10x', prod, 'GEO filtered matrix', bc,
                            'barcodes.tsv of the deposited filtered matrix',
                            n_called_pre_qc=n))

    # --- Cillo: per-GSM GEX tar, counted from the tar index (nothing extracted)
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
        rows.append(row('cillo', ch, '10x', 'build_cillo_decontx_input.py',
                        'GEO filtered matrix in GEX tar', tp,
                        'barcodes.tsv of the deposited filtered matrix',
                        n_called_pre_qc=n))

    # --- Liu: fixed 30,000-barcode demuxlet whitelist per lane.
    # Liu is the ONE cohort that publishes its demultiplexing breakdown, so it is
    # the only place n_cross_donor_doublets_removed is directly observed rather
    # than assumed -- which matters, because the residual-rate formula's premise
    # is precisely that cross-donor doublets were removed upstream.
    demux = f'{RAW}/GSE161918_liu/singlets/liu_demuxlet_filter_summary.csv'
    dm = (pd.read_csv(demux).set_index('lane') if os.path.exists(demux) else None)
    for h5 in sorted(glob.glob(f'{RAW}/GSE161918_liu/rna_h5/*.h5')):
        base = os.path.basename(h5)
        ch = base.split('_RNA_')[0].split('_', 1)[1]
        try:
            with h5py.File(h5, 'r') as f:
                n = f[list(f.keys())[0]]['barcodes'].shape[0]
        except Exception:
            n = pd.NA
        sng = dbl = amb = pd.NA
        if dm is not None and ch in dm.index:
            sng = int(dm.loc[ch, 'n_sng'])
            dbl = int(dm.loc[ch, 'n_dbl'])
            amb = int(dm.loc[ch, 'n_amb'])
        rows.append(row('liu', ch, '10x', 'build_liu_decontx_input.py',
                        'demuxlet whitelist h5 + liu_demuxlet_filter_summary.csv', h5,
                        'fixed 30,000-barcode per-lane demuxlet whitelist, identical '
                        'across all 55 lanes -- neither a droplet count nor a cell call; '
                        'the SNG/DBL/AMB split is demuxlet output over that whitelist, so '
                        'it is a genotype classification of candidate barcodes, not of '
                        'confirmed cells',
                        reason='the only per-lane count is a fixed whitelist size, identical '
                               'for every lane, so it carries no capture-specific loading '
                               'information',
                        n_candidate_barcodes=n, n_post_demux_singlets=sng,
                        n_cross_donor_doublets_removed=dbl, n_unresolved_or_amb=amb))

    # --- Wilk: Seq-Well; count from converted/, which convert_wilk_rds_to_mtx.R wrote
    for d in sorted(glob.glob(f'{RAW}/GSE174072_wilk/converted/*')):
        if not os.path.isdir(d):
            continue
        bc = f'{d}/barcodes.tsv'
        n = sum(1 for _ in open(bc)) if os.path.exists(bc) else pd.NA
        rows.append(row('wilk', os.path.basename(d), 'seqwell',
                        'convert_wilk_rds_to_mtx.R -> build_wilk_decontx_input.py',
                        'author-processed Seq-Well RDS, converted to mtx', bc,
                        'barcodes of the deposited Seq-Well cell x gene matrix',
                        reason='Seq-Well is not droplet loading -- the 10x 0.8%/1,000 '
                               'multiplet model does not apply, and no Seq-Well loading '
                               'record exists for this cohort',
                        n_called_pre_qc=n))

    # --- COMBAT, Stephenson: author-processed h5ad; layers['raw'] is raw VALUES
    for cohort, path, key, prod in [
            ('combat', f'{RAW}/COMBAT/COMBAT-CITESeq-DATA.h5ad', 'Channel_ID',
             'build_combat_decontx_input.py'),
            ('stephenson', f'{RAW}/Stephenson_EMTAB10026/covid_portal_210320_with_raw.h5ad',
             'sample_id', 'build_stephenson_decontx_input.py')]:
        with h5py.File(path, 'r') as f:
            vals = read_obs_categorical(f, key)
        u, cnt = np.unique(vals, return_counts=True)
        ck = cheap_checksum(path)
        for ch, n in zip(u.tolist(), cnt.tolist()):
            r = row(cohort, ch, '10x', prod, 'author-processed h5ad', None,
                    "cells in the depositor-assembled object; layers['raw'] holds raw "
                    'count VALUES, not raw droplets',
                    obs_field=f'obs/{key}',
                    reason="layers['raw'] holds raw count VALUES, not raw droplets; the "
                           "cell set is the depositor's own post-QC selection and its "
                           'filtering is not recorded',
                    n_called_pre_qc=int(n))
            r['source_artifact'] = path
            r['source_checksum'] = ck
            rows.append(r)

    # --- Ren: called-cell annotation, stored outside this project
    if os.path.exists(REN_SRC):
        per = pd.read_csv(REN_SRC, usecols=['sampleID']).sampleID.astype(str).value_counts()
        ck = cheap_checksum(REN_SRC)
        for ch, n in per.items():
            r = row('ren', ch, '10x', 'build_ren_decontx_input.py',
                    'author cell annotation (pooled)', None,
                    "barcodes in the authors' cell_annotation.csv; every barcode carries a "
                    'celltype, so this is a called-cell set',
                    obs_field='cell_annotation.sampleID',
                    reason='only a called-cell annotation was deposited; no unfiltered '
                           'droplet matrix and no per-library recovered-cell record',
                    n_called_pre_qc=int(n))
            r['source_artifact'] = REN_SRC
            r['source_checksum'] = ck
            rows.append(r)

    return pd.DataFrame(rows, columns=COLS)


def main():
    os.makedirs(OUT, exist_ok=True)
    df = collect()

    # reconcile against the authoritative Gate 2 capture universe, both directions
    extra = []
    for cohort in sorted(df.cohort.unique()):
        g2 = gate2_captures(cohort)
        m = df.cohort == cohort
        df.loc[m, 'in_gate2'] = df.loc[m, 'channel_id'].isin(g2)
        df.loc[m, 'n_gate1_retained'] = df.loc[m, 'channel_id'].map(g2)
        df.loc[m, 'match'] = np.where(df.loc[m, 'in_gate2'], 'matched', 'source_only')
        for ch in sorted(set(g2) - set(df.loc[m, 'channel_id'])):
            r = {c: pd.NA for c in COLS}
            r.update(cohort=cohort, channel_id=ch, in_gate2=True, match='gate2_only',
                     n_gate1_retained=g2[ch],
                     reason_not_eligible='Gate 2 capture with no located source artifact')
            extra.append(r)
    if extra:
        df = pd.concat([df, pd.DataFrame(extra, columns=COLS)], ignore_index=True)

    df = df.sort_values(['cohort', 'channel_id']).reset_index(drop=True)
    df.to_csv(f'{OUT}/gate2_source_crosswalk.csv', index=False)
    print(f'wrote {len(df)} rows -> {OUT}/gate2_source_crosswalk.csv\n')

    piv = (df.pivot_table(index='cohort', columns='match', values='channel_id',
                          aggfunc='count', fill_value=0))
    piv['gate2_total'] = [len(gate2_captures(c)) for c in piv.index]
    print('Source <-> Gate 2 reconciliation:')
    print(piv.to_string())
    print(f'\ntotal source artifacts : {int((df.match != "gate2_only").sum())}')
    print(f'total Gate 2 captures  : {int((df.in_gate2 == True).sum())}')

    print('\nCount levels obtainable, matched captures only (non-NA):')
    lv = ['n_input_barcodes', 'n_candidate_barcodes', 'n_called_pre_qc',
          'n_post_demux_singlets', 'n_cross_donor_doublets_removed',
          'n_unresolved_or_amb', 'n_gate1_retained']
    print(df[df.match == 'matched'].groupby('cohort')[lv].count().to_string())

    print('\nloading_prior_eligibility:',
          df.loading_prior_eligibility.value_counts(dropna=False).to_dict())


if __name__ == '__main__':
    main()
