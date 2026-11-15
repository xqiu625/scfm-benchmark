"""READ-ONLY Stage 31 preflight: can the redesigned workflow's inputs be identified?

    python scripts/31_donor_level_revision/code/stage31_input_preflight.py --out-dir <NEW dir>

Answers, for the 302 pre-subsample COVID `Sample_ID`s across nine `Data_NO`
groups, what the redesign needs before anything can be scored:

  * what the raw-count source is, and whether it is un-normalised RNA counts;
  * what the physical capture is, and whether a `Sample_ID` is one;
  * how a canonical cell ID decomposes into an id shape plus a nucleotide token;
  * whether every cell of the 60-cell analytic subset maps back to exactly one
    (capture, raw barcode) -- and it reports that it does NOT, because only a
    container-internal token candidate is available.

It produces NO doublet score, call or rate, selects NO `n_prin_comps`, and
submits nothing. A field that cannot be proven is written `unknown` -- never
inferred from string similarity or from a matching total cell count.

OUTCOME BOUNDARY
----------------
`h5py` opens NAMED DATASETS, so only the technical obs columns below are
materialised from a container whose obs also holds `Death`, `covid_severity`,
`mild_vs_severe`, `severity_score`, `three_group` and `progression_group`. None
of those is opened. The CSV inputs are read with `usecols`.

`disease` IS read from the patient list, because the COVID population cannot be
defined without it. It is a POPULATION SELECTOR, not an outcome: it selects
which patients are in scope and is recorded as such in the manifest. No value of
it reaches any output -- only the resulting id set does.

No output carries mortality, Death, survival, severity or any other clinical
value. The manifest separates `loaded` from `inspected` for every input.
"""
import argparse
import collections
import hashlib
import json
import os
import re
import sys
import time

import h5py
import numpy as np
import pandas as pd

RAW_H5AD = ("/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/"
            "human_platelet_covid_severity_raw_counts.h5ad")
PATIENT_LIST_CSV = ("/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/"
                    "31_donor_level_revision/results_patient_v2/qc_baseline/"
                    "patient_qc_distributions.csv")
SAMPLE_MAP_CSV = ("/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/"
                  "31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv")

#: Mirrors the production entry point, so the subset this preflight reconstructs
#: is the subset that entry point would build. Both constants are read from
#: doublet_detection_and_reaudit.py at run time and asserted equal, so a change
#: there fails this preflight instead of silently invalidating it.
CELLS_PER_PATIENT = 60
RNG_SEED = 20260910
ENTRY_POINT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'doublet_detection_and_reaudit.py')

#: The only obs datasets opened. Technical identity and technical QC only.
OBS_TECHNICAL = ('Sample_ID', 'Data_NO', 'Cell_ID', 'Cell_ID2', 'orig.ident',
                 'nCount_RNA', 'nFeature_RNA')
OBS_NEVER_OPENED = ('Death', 'covid_severity', 'mild_vs_severe', 'severity_score',
                    'three_group', 'progression_group', 'label', 'y')

UNKNOWN = 'unknown'
BOUNDARY = []


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b''):
            h.update(chunk)
    return h.hexdigest()


def cheap_identity(path, nbytes=1 << 20):
    """`<size>:<partial>` for files too large to hash whole, same form the
    repository already uses for the COMBAT container."""
    sz = os.path.getsize(path)
    h = hashlib.sha256(str(sz).encode())
    with open(path, 'rb') as fh:
        h.update(fh.read(nbytes))
        if sz > nbytes:
            fh.seek(max(0, sz - nbytes))
            h.update(fh.read(nbytes))
    return f'{sz}:{h.hexdigest()[:16]}'


def read_categorical(obs, name):
    node = obs[name]
    if isinstance(node, h5py.Group):
        codes, cats = node['codes'][:], node['categories'][:]
    else:
        if '__categories' not in obs or name not in obs['__categories']:
            return node[:]                       # plain numeric column
        codes, cats = node[:], obs['__categories'][name][:]
    cats = np.array([c.decode() if isinstance(c, bytes) else str(c) for c in cats])
    return cats[codes]


def read_csv_cols(path, usecols, why):
    present = list(pd.read_csv(path, nrows=0).columns)
    missing = [c for c in usecols if c not in present]
    if missing:
        raise SystemExit(f'{path}: missing {missing}')
    BOUNDARY.append({'input': path, 'kind': 'csv', 'columns_read': ';'.join(usecols),
                     'columns_excluded_at_parse_time':
                         ';'.join(c for c in present if c not in usecols),
                     'values_inspected': why})
    return pd.read_csv(path, usecols=usecols, dtype=str)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out-dir', required=True, help='NEW directory; must be empty')
    args = ap.parse_args()
    out = args.out_dir
    os.makedirs(out, exist_ok=True)
    if os.listdir(out):
        raise SystemExit(f'{out} is not empty; this preflight writes only to a new directory')

    src = open(ENTRY_POINT).read()
    for const, val in (('CELLS_PER_PATIENT', CELLS_PER_PATIENT), ('RNG_SEED', RNG_SEED)):
        if f'{const} = {val}' not in src:
            raise SystemExit(f'{const} no longer equals {val} in the entry point; '
                             'the subset this preflight reconstructs would not be the '
                             'one the entry point builds')

    # ---- population: COVID patients -> Sample_IDs --------------------------
    patients = read_csv_cols(PATIENT_LIST_CSV, ['patient_id', 'disease'],
                             'disease inspected as a POPULATION SELECTOR only')
    covid = set(patients.loc[patients['disease'] == 'COVID-19', 'patient_id'])
    smap = read_csv_cols(SAMPLE_MAP_CSV, ['sample_id', 'patient_id'],
                         'not inspected; used as an identity map')
    smap = smap[smap['patient_id'].isin(covid)]
    sample_to_patient = dict(zip(smap['sample_id'], smap['patient_id']))

    # ---- container: technical obs only -------------------------------------
    t0 = time.time()
    with h5py.File(RAW_H5AD, 'r') as f:
        obs = f['obs']
        opened = [c for c in OBS_TECHNICAL if c in obs]
        leaked = [c for c in OBS_NEVER_OPENED if c in opened]
        if leaked:
            raise SystemExit(f'refusing to open outcome columns {leaked}')
        cols = {c: read_categorical(obs, c) for c in opened}
        cell_id = np.array([x.decode() if isinstance(x, bytes) else str(x)
                            for x in obs['_index'][:]])
        n_obs, n_var = f['X'].shape
        x_dtype = str(f['X'].dtype)
        probe = f['X'][:200, :500]
        x_integral = bool(np.all(np.equal(np.mod(probe, 1), 0)))
        x_min = float(probe.min())
        has_layers = bool(list(f['layers'].keys())) if 'layers' in f else False
    BOUNDARY.append({'input': RAW_H5AD, 'kind': 'h5ad',
                     'columns_read': ';'.join(opened),
                     'columns_excluded_at_parse_time':
                         ';'.join(c for c in obs_keys_snapshot(RAW_H5AD) if c not in opened),
                     'values_inspected': 'technical identity/QC only; no outcome dataset opened'})

    df = pd.DataFrame({'cell_id': cell_id, **{c: cols[c] for c in opened}})
    df['patient_id'] = df['Sample_ID'].astype(str).map(sample_to_patient)
    full = df[df['patient_id'].notna()].copy()

    # ---- canonical cell-id decomposition, MEASURED -------------------------
    # There is NO single rule. A first pass assumed `<barcode>-<gem>_<suffix>`
    # and left 76% of ids unparsed, so the namespace is characterised instead of
    # assumed: the only element common to every id is the nucleotide run, and
    # everything around it is recorded as an observed SHAPE per Data_NO group.
    BC = re.compile(r'[ACGTN]{8,}')

    def barcode_of(cid):
        m = BC.search(cid or '')
        return m.group(0) if m else None

    def shape_of(cid):
        if not isinstance(cid, str):
            return UNKNOWN
        return re.sub(r'\d+', '<n>', BC.sub('<BC>', cid))

    # NOT a verified raw barcode. This token is extracted from the canonical
    # cell-id STRING in this same container. No raw matrix or barcode file was
    # opened, and the token was never joined against a 10x barcode universe, so
    # calling it a raw barcode would claim an identity mapping never performed.
    full['nucleotide_token_candidate'] = full['cell_id'].map(barcode_of)
    full['id_shape'] = full['cell_id'].map(shape_of)
    full['gemgroup'] = UNKNOWN
    full['suffix'] = UNKNOWN
    n_unparsed = int(full['nucleotide_token_candidate'].isna().sum())

    shapes_per_dno = full.groupby('Data_NO')['id_shape'].nunique()
    dno_per_shape = full.groupby('id_shape')['Data_NO'].nunique()
    suffix_rule = (
        f'NO single rule: {full["id_shape"].nunique()} distinct cell-id shapes across '
        f'{full["Data_NO"].nunique()} Data_NO groups; max {int(shapes_per_dno.max())} '
        'shapes within one group. Each group needs its own documented '
        'prefix/suffix rule before a barcode namespace can be fixed.')
    suffix_per_sample = shapes_per_dno            # reported under measured names below
    sample_per_suffix = dno_per_shape
    suffix_per_dno = shapes_per_dno

    shape_tbl = (full.groupby(['Data_NO', 'id_shape'])
                 .agg(n_cells=('cell_id', 'size'),
                      n_sample_ids=('Sample_ID', 'nunique'),
                      example=('cell_id', 'first'))
                 .reset_index().sort_values(['Data_NO', 'n_cells'],
                                            ascending=[True, False]))
    shape_tbl.to_csv(os.path.join(out, 'cell_id_namespace_by_data_no.csv'), index=False)

    # ---- the 60-cell analytic subset, reconstructed ------------------------
    rng = np.random.default_rng(RNG_SEED)
    keep = []
    for pid, grp in full.groupby('patient_id'):
        take = min(len(grp), CELLS_PER_PATIENT)
        keep.extend(rng.choice(grp.index.to_numpy(), size=take, replace=False).tolist())
    subset = full.loc[keep]

    # ---- per-Sample_ID inventory ------------------------------------------
    rows = []
    for sid, g in full.groupby('Sample_ID'):
        sub = subset[subset['Sample_ID'] == sid]
        bc_dup = int(g['nucleotide_token_candidate'].duplicated().sum()) if n_unparsed == 0 else -1
        rows.append({
            'sample_id': sid,
            'data_no': ';'.join(sorted(set(g['Data_NO']))),
            'patient_id': ';'.join(sorted(set(g['patient_id']))),
            'n_cells_pre_subsample': len(g),
            'n_cells_in_60_cell_subset': len(sub),
            'raw_count_source_path': RAW_H5AD,
            'raw_count_source_format': 'h5ad, dense X, no layers',
            'raw_count_dtype': x_dtype,
            'raw_counts_integral_probe': x_integral,
            'raw_counts_un_normalised': UNKNOWN,
            'matrix_origin_layer': 'X (the only matrix present)',
            'physical_capture_key': UNKNOWN,
            'physical_capture_evidence': 'not established by this preflight',
            'barcode_namespace': 'container-internal canonical cell_id; the nucleotide '
                                 'token is a candidate, not a verified raw barcode',
            'raw_barcode_identity_verified': False,
            'n_token_collisions_within_sample': bc_dup,
            'input_stage': 'canonical integrated container, pre-60-cell-subsample',
            'parent_states_retained': UNKNOWN,
            'n_cell_id_shapes': int(g['id_shape'].nunique()),
            'cell_id_shapes_seen': ';'.join(sorted(set(g['id_shape']))[:4]),
        })
    inv = pd.DataFrame(rows).sort_values('sample_id')
    inv.to_csv(os.path.join(out, 'source_inventory.csv'), index=False)

    # ---- cell-level crosswalk ---------------------------------------------
    cw = full[['cell_id', 'nucleotide_token_candidate', 'id_shape', 'Sample_ID', 'Data_NO',
               'patient_id']].copy()
    cw['in_60_cell_subset'] = cw.index.isin(set(keep))
    cw['physical_capture_key'] = UNKNOWN
    cw['raw_barcode_identity_verified'] = False
    cw.rename(columns={'Sample_ID': 'sample_id', 'Data_NO': 'data_no'}, inplace=True)
    cw.to_csv(os.path.join(out, 'capture_barcode_crosswalk.csv'), index=False)

    # ---- subset mapping ----------------------------------------------------
    per_cell_key = full.set_index('cell_id')[['Sample_ID', 'nucleotide_token_candidate']]
    # Only meaningful where a token was extractable at all.
    with_bc = full[full['nucleotide_token_candidate'].notna()]
    dup_key = int(with_bc.duplicated(subset=['Sample_ID', 'nucleotide_token_candidate']).sum())
    sm = []
    for sid, g in subset.groupby('Sample_ID'):
        gb = g[g['nucleotide_token_candidate'].notna()]
        uniq = int((~gb.duplicated(subset=['Sample_ID', 'nucleotide_token_candidate'])).sum())
        sm.append({'sample_id': sid,
                   'n_subset_cells': len(g),
                   'n_uniquely_mapped_to_sample_plus_token_within_container': uniq,
                   'n_ambiguous': len(g) - uniq,
                   'n_unmatched_in_full_input': int((~g['cell_id']
                                                     .isin(per_cell_key.index)).sum()),
                   'subset_link_basis': 'canonical cell_id -> the same pre-subsample '
                                        'container (NOT a raw source)',
                   'raw_source_link_established': 'no',
                   'physical_capture_resolved': 'no'})
    pd.DataFrame(sm).sort_values('sample_id').to_csv(
        os.path.join(out, 'subset_mapping_summary.csv'), index=False)

    manifest = {
        'preflight': 'stage31_input_preflight/1',
        'generated_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        'read_only': True,
        'produced_doublet_scores_calls_or_rates': False,
        'selected_n_prin_comps': False,
        'submitted_jobs': False,
        'entry_point': {'path': ENTRY_POINT, 'sha256': sha256(ENTRY_POINT),
                        'cells_per_patient': CELLS_PER_PATIENT, 'rng_seed': RNG_SEED},
        'inputs': {
            'raw_h5ad': {'path': RAW_H5AD, 'identity': cheap_identity(RAW_H5AD),
                         'shape': [int(n_obs), int(n_var)], 'x_dtype': x_dtype,
                         'x_integral_on_probe': x_integral, 'x_min_on_probe': x_min,
                         'has_layers': has_layers},
            'patient_list': {'path': PATIENT_LIST_CSV, 'sha256': sha256(PATIENT_LIST_CSV)},
            'sample_map': {'path': SAMPLE_MAP_CSV, 'sha256': sha256(SAMPLE_MAP_CSV)},
        },
        'outcome_boundary': {
            'known_outcome_values_loaded': [],
            'outcome_values_inspected': [],
            'obs_datasets_opened': opened,
            'obs_datasets_never_opened': list(OBS_NEVER_OPENED),
            'population_selector_inspected': 'patient_list.disease',
            'per_input': BOUNDARY,
        },
        'population': {
            'covid_patients': len(covid),
            'sample_ids_pre_subsample': int(full['Sample_ID'].nunique()),
            'data_no_groups': sorted(set(full['Data_NO'])),
            'n_data_no_groups': int(full['Data_NO'].nunique()),
            'cells_pre_subsample': int(len(full)),
            'cells_in_60_cell_subset': int(len(subset)),
            'container_sample_ids_all_diseases': int(df['Sample_ID'].nunique()),
        },
        'identity': {
            'cell_ids_unique': bool(full['cell_id'].is_unique),
            'cell_ids_unparsed': n_unparsed,
            'sample_plus_token_duplicates_within_container': dup_key,
            'suffix_rule_measured': suffix_rule,
            'n_cell_id_shapes_total': int(full['id_shape'].nunique()),
            'max_shapes_within_one_data_no': int(shapes_per_dno.max()),
            'max_data_no_groups_sharing_one_shape': int(dno_per_shape.max()),
            'cells_with_no_extractable_nucleotide_token': n_unparsed,
            'raw_barcode_identity_verified': False,
            'identity_scope': 'container-internal identity candidate only. The nucleotide '
                              'token is extracted from the canonical cell-id string in '
                              'this same h5ad; no raw matrix or barcode file was opened '
                              'and no join against a 10x barcode universe was performed. '
                              'The subset link is canonical cell_id -> the same '
                              'pre-subsample container, NOT a link to a raw source.',
            'sample_plus_token_duplicates_where_token_extractable': dup_key,
        },
        'unknown_fields': {
            'raw_barcode_identity': 'NOT verified. Only a container-internal nucleotide '
                                    'token candidate exists; the raw-input identity '
                                    'mapping has not been performed.',
            'physical_capture_key': 'not established for any Sample_ID',
            'raw_counts_un_normalised':
                'integral dtype is necessary, not sufficient; provenance not established',
            'parent_states_retained':
                'not established. The container is named for platelets and most cell ids '
                'carry no lineage token, but 6,030 carry a Megakaryocytes prefix and a '
                'handful carry T_cell/B_cell/Other_cell. Those tokens are inherited from '
                'upstream object names and are NOT evidence of the current lineage '
                'composition. Whether this container retains the parent states C3 '
                'requires is unresolved either way.',
            'subsample_drops_whole_sample_ids':
                'the 60-cell cap is applied PER PATIENT, so a patient with several '
                'Sample_IDs can lose one entirely: 004_TA has 16 cells pre-subsample and '
                '0 in the subset. That is the documented 302-vs-301 gap, measured.',
        },
        'elapsed_seconds': round(time.time() - t0, 1),
    }
    with open(os.path.join(out, 'PRECHECK_MANIFEST.json'), 'w') as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)

    print(f'sample_ids={manifest["population"]["sample_ids_pre_subsample"]} '
          f'data_no={manifest["population"]["n_data_no_groups"]} '
          f'cells={manifest["population"]["cells_pre_subsample"]} '
          f'subset={manifest["population"]["cells_in_60_cell_subset"]}')
    print(f'suffix rule: {suffix_rule}; unparsed cell ids: {n_unparsed}; '
          f'sample+barcode duplicates: {dup_key}')
    print(f'cell-id shapes: {full["id_shape"].nunique()} total, '
          f'max {int(shapes_per_dno.max())} within one Data_NO group')
    print(f'wrote 5 artifacts to {out}')


def obs_keys_snapshot(path):
    with h5py.File(path, 'r') as f:
        return sorted(k for k in f['obs'].keys() if not k.startswith('_')
                      and k != '__categories')


if __name__ == '__main__':
    main()
