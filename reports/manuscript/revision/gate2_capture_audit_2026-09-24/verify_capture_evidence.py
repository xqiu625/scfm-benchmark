"""READ-ONLY, fail-closed verification of the GSM5102900_HC1 capture evidence.

    python verify_capture_evidence.py [--out-dir DIR]

Recomputes every quantity CAPTURE_C1_C3_AUDIT.md asserts and refuses to finish if
any of them disagrees. The audit document is therefore replayable rather than a
set of numbers a reader must take on trust.

OUTCOME BOUNDARY
----------------
Reads three deposit files (barcodes / features / matrix) and exactly two named
datasets from the canonical container: `obs/_index` and `obs/Sample_ID`. `h5py`
opens named datasets, so the container's clinical columns are never materialised.
No outcome column or value is read, written or inspected. The capture is a
healthy control, so no mortality label exists for it in any case.

Submits nothing, scores nothing, deletes nothing, and leaves
PRODUCTION_PARAMETER_RULING untouched.
"""
import argparse
import gzip
import hashlib
import json
import os
import re
import sys
import time

import h5py
import numpy as np

CAPTURE = 'GSM5102900_HC1'
GSM = 'GSM5102900'
DEPOSIT_DIR = '/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/GSE167363_qiu/extracted'
CONTAINER = ('/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/'
             'human_platelet_covid_severity_raw_counts.h5ad')
GEO_SERIES_URL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167363'
GEO_SAMPLE_URL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5102900'
GEO_SUPPL_BASE = 'ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5102nnn/GSM5102900/suppl/'

#: Every value the audit document states. A mismatch stops the run.
#: Full 64-character SHA-256 of every input, ASSERTED. A previous version declared
#: a 16-character prefix and never compared any hash at all -- the manifest's
#: values were observations, not checks.
EXPECTED_SHA256 = {
    'GSM5102900_HC1_barcodes.tsv.gz':
        'e6dbf9f587d9a6cf3ecff31cbf4d53b09c6faf9837c82d5c1c0a26227f377a02',
    'GSM5102900_HC1_features.tsv.gz':
        'c7408e82af270ca6308a8aa7ea659795a15ab76f2e3a85bd5db4cd3d53244199',
    'GSM5102900_HC1_matrix.mtx.gz':
        '4ea259bfc4de02cc9799df563779fce7c9bd4c889532d3ff1178e2caf545b9f8',
}
EXPECTED_CONTAINER_SHA256 = (
    '1a0f37902cd33f410a71428ea45b3d1573fddec393ce27f07a93074d0ed536ba')
EXPECTED_CONTAINER_BYTES = 5027039977

EXPECTED = {
    'n_genes_declared': 33538,
    'n_cells_declared': 7666,
    'nnz_declared': 9903096,
    'n_features': 33538,
    'feature_types': {'Gene Expression'},
    'n_barcodes': 7666,
    'barcode_token_lengths': {16},
    'barcode_gem_suffixes': {'1'},
    'container_hc1_cells': 65,
    'container_tokens_in_deposit': 65,
    'genes_passing_min_cells_3': 17226,
    'cells_passing_min_genes_3': 7666,
}
SAMPLE_ID_IN_CONTAINER = 'HC1'


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out-dir', default=os.path.dirname(os.path.abspath(__file__)))
    args = ap.parse_args()
    t0 = time.time()
    fail = []

    def check(name, got, want):
        ok = got == want
        if not ok:
            fail.append(f'{name}: recomputed {got!r} but the audit states {want!r}')
        return {'value': got, 'expected': want, 'agrees': ok}

    paths = {k: os.path.join(DEPOSIT_DIR, f'{CAPTURE}_{k}')
             for k in ('barcodes.tsv.gz', 'features.tsv.gz', 'matrix.mtx.gz')}
    files = {}
    for k, p in paths.items():
        if not os.path.exists(p):
            raise SystemExit(f'input missing: {p}')
        digest = sha256(p)
        want = EXPECTED_SHA256[f'{CAPTURE}_{k}']
        if digest != want:
            raise SystemExit(f'{CAPTURE}_{k}: sha256 {digest} != expected {want}')
        files[k] = {'path': p, 'bytes': os.path.getsize(p), 'sha256': digest,
                    'sha256_asserted': True,
                    'geo_supplementary_url': GEO_SUPPL_BASE + f'{CAPTURE}_{k}'}
    csize = os.path.getsize(CONTAINER)
    if csize != EXPECTED_CONTAINER_BYTES:
        raise SystemExit(f'container: {csize} bytes != {EXPECTED_CONTAINER_BYTES}')
    cdigest = sha256(CONTAINER)
    if cdigest != EXPECTED_CONTAINER_SHA256:
        raise SystemExit(f'container: sha256 {cdigest} != {EXPECTED_CONTAINER_SHA256}')

    # ---- C1: the matrix, streamed (never materialised) ---------------------
    with gzip.open(paths['matrix.mtx.gz'], 'rt') as fh:
        banner = fh.readline().strip()
        fh.readline()
        ngene, ncell, nnz = (int(x) for x in fh.readline().split())
        gene_cells = np.zeros(ngene + 1, dtype=np.int32)
        cell_genes = np.zeros(ncell + 1, dtype=np.int32)
        n_lines, non_integer, vmin, vmax = 0, 0, None, None
        for line in fh:
            a = line.find(' ')
            b = line.find(' ', a + 1)
            g = int(line[:a]); c = int(line[a + 1:b]); v = float(line[b + 1:])
            gene_cells[g] += 1
            cell_genes[c] += 1
            if v != int(v):
                non_integer += 1
            vmin = v if vmin is None or v < vmin else vmin
            vmax = v if vmax is None or v > vmax else vmax
            n_lines += 1

    c1 = {
        'matrixmarket_banner': banner,
        'banner_declares_integer': 'integer' in banner,
        'n_genes': check('n_genes', ngene, EXPECTED['n_genes_declared']),
        'n_cells': check('n_cells', ncell, EXPECTED['n_cells_declared']),
        'nnz_declared': check('nnz_declared', nnz, EXPECTED['nnz_declared']),
        'nnz_streamed': check('nnz_streamed', n_lines, EXPECTED['nnz_declared']),
        'non_integer_values': check('non_integer_values', non_integer, 0),
        'value_min': vmin, 'value_max': vmax,
    }
    if not c1['banner_declares_integer']:
        fail.append('MatrixMarket banner does not declare integer')

    with gzip.open(paths['features.tsv.gz'], 'rt') as fh:
        feats = [l.rstrip('\n').split('\t') for l in fh if l.strip()]
    ftypes = {f[2] for f in feats if len(f) > 2}
    c1['n_features'] = check('n_features', len(feats), EXPECTED['n_features'])
    c1['feature_types'] = check('feature_types', ftypes, EXPECTED['feature_types'])

    # ---- C2: the barcode namespace ----------------------------------------
    with gzip.open(paths['barcodes.tsv.gz'], 'rt') as fh:
        bcs = [l.strip() for l in fh if l.strip()]
    tokens = [b.split('-')[0] for b in bcs]
    suffixes = {b.split('-')[1] for b in bcs if '-' in b}
    c2 = {
        'n_barcodes': check('n_barcodes', len(bcs), EXPECTED['n_barcodes']),
        'all_unique': check('barcodes_unique', len(set(bcs)), len(bcs)),
        'token_lengths': check('token_lengths', {len(t) for t in tokens},
                               EXPECTED['barcode_token_lengths']),
        'gem_suffixes': check('gem_suffixes', suffixes, EXPECTED['barcode_gem_suffixes']),
        'example_barcode': bcs[0],
    }

    # ---- C3: container subset, row-level join ------------------------------
    with h5py.File(CONTAINER, 'r') as f:
        obs = f['obs']
        opened = ['_index', 'Sample_ID']
        node = obs['Sample_ID']
        if isinstance(node, h5py.Group):
            codes, cats = node['codes'][:], node['categories'][:]
        else:
            codes, cats = node[:], obs['__categories']['Sample_ID'][:]
        cats = np.array([c.decode() if isinstance(c, bytes) else str(c) for c in cats])
        sample_id = cats[codes]
        cell_id = np.array([x.decode() if isinstance(x, bytes) else str(x)
                            for x in obs['_index'][:]])
    hc1 = cell_id[sample_id == SAMPLE_ID_IN_CONTAINER]
    BC = re.compile(r'[ACGTN]{8,}')
    hc1_tokens = [BC.search(c).group(0) for c in hc1]
    deposit_tokens = set(tokens)
    found = sum(1 for t in hc1_tokens if t in deposit_tokens)
    c3 = {
        'container_cells': check('container_hc1_cells', len(hc1),
                                 EXPECTED['container_hc1_cells']),
        'tokens_found_in_deposit': check('container_tokens_in_deposit', found,
                                         EXPECTED['container_tokens_in_deposit']),
        'container_is_subset': check('container_is_subset', found == len(hc1), True),
        'deposit_barcodes': len(deposit_tokens),
        'retained_fraction_percent': round(100.0 * len(hc1) / len(deposit_tokens), 4),
        'complete_declared_stage_population': len(deposit_tokens),
        'container_admissible_for_c3': False,
    }

    # ---- feasibility inputs (streamed; HVG deliberately NOT computed) ------
    feas = {
        'genes_passing_min_cells_3': check('genes_passing_min_cells_3',
                                           int((gene_cells[1:] >= 3).sum()),
                                           EXPECTED['genes_passing_min_cells_3']),
        'cells_passing_min_genes_3': check('cells_passing_min_genes_3',
                                           int((cell_genes[1:] >= 3).sum()),
                                           EXPECTED['cells_passing_min_genes_3']),
        'median_genes_per_cell': int(np.median(cell_genes[1:])),
        'min_genes_per_cell': int(cell_genes[1:].min()),
        'n_usable_features_hvg': None,
        'n_usable_features_status': 'NOT MEASURED - requires materialising the matrix; '
                                    'the login-node attempt was OOM-killed and no SLURM job '
                                    'is authorised',
        'strict_rule': 'n_prin_comps < min(n_cells_effective, n_usable_features)',
        'approved_n_prin_comps': 30,
        'feasibility_confirmed': False,
    }

    manifest = {
        'verification': 'gate2_capture_evidence/1',
        'capture': CAPTURE,
        'geo_sample': GSM,
        'generated_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        'read_only': True,
        'submitted_jobs': False,
        'produced_scores_calls_or_rates': False,
        'outcome_columns_read': [],
        'container_obs_datasets_opened': opened,
        'sources': {
            'deposit_files': files,
            'geo_series_url': GEO_SERIES_URL,
            'geo_sample_url': GEO_SAMPLE_URL,
            'canonical_container': {'path': CONTAINER, 'bytes': csize,
                                    'sha256': cdigest, 'sha256_asserted': True},
        },
        'C1_raw_count_provenance': c1,
        'C2_capture_and_barcode_namespace': c2,
        'C3_population_completeness': c3,
        'pc_feasibility': feas,
        'disagreements': fail,
        'elapsed_seconds': round(time.time() - t0, 1),
    }
    dest = os.path.join(args.out_dir, 'EVIDENCE_MANIFEST.json')
    with open(dest, 'w') as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True, default=str)

    for k, v in files.items():
        print(f'  {k:<18} {v["bytes"]:>12,} bytes  sha256 {v["sha256"]}')
    print(f'  C1 banner            : {banner}')
    print(f'  C1 dims/nnz          : {ngene} x {ncell}, nnz {nnz} (streamed {n_lines})')
    print(f'  C1 non-integer values: {non_integer}')
    print(f'  C1 feature types     : {sorted(ftypes)}')
    print(f'  C2 barcodes          : {len(bcs)} unique={len(set(bcs))} '
          f'lengths={sorted({len(t) for t in tokens})} suffixes={sorted(suffixes)}')
    print(f'  C3 container/deposit : {len(hc1)}/{len(deposit_tokens)} '
          f'= {c3["retained_fraction_percent"]}%  join {found}/{len(hc1)}')
    print(f'  feasibility          : N_eff={feas["cells_passing_min_genes_3"]["value"]} '
          f'G_HVG=NOT MEASURED -> PC=30 unconfirmed')
    print(f'  wrote {dest}')
    if fail:
        print('\nFAIL-CLOSED: recomputation disagrees with the audit document:', file=sys.stderr)
        for f_ in fail:
            print('  ' + f_, file=sys.stderr)
        raise SystemExit(1)
    print('\nall recomputed values agree with CAPTURE_C1_C3_AUDIT.md')


if __name__ == '__main__':
    main()
