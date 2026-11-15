"""READ-ONLY, ad-hoc C1-C3 style evidence check for GSM5102901_HC2 (Qiu, GSE167363),
proposed as the next Gate 2 capture candidate. Modeled on the already-reviewed
reports/manuscript/revision/gate2_capture_audit_2026-09-24/verify_capture_evidence.py
(which is hardcoded to GSM5102900_HC1 and is NOT modified here), but this is a
fresh, unreviewed script for candidate-evidence gathering only -- not part of
IMPLEMENTATION_FILES, not a permit input, not run against the scoring entry point.

OUTCOME BOUNDARY: reads three deposit files (barcodes/features/matrix) and exactly
two named datasets from the canonical container (obs/_index, obs/Sample_ID). No
outcome column read. HC2 is a healthy control; no mortality label exists for it.

Submits nothing, scores nothing, deletes nothing.

Writes gate2_hc2_candidate_evidence.json to --out-dir (refuses to overwrite an
existing non-empty directory -- a prior version wrote to a single fixed remote
path and would silently overwrite its own prior result on a rerun).
"""
import argparse
import gzip
import hashlib
import json
import os
import re

import h5py
import numpy as np

CAPTURE = 'GSM5102901_HC2'
GSM = 'GSM5102901'
DEPOSIT_DIR = '/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/GSE167363_qiu/extracted'
CONTAINER = ('/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/'
             'human_platelet_covid_severity_raw_counts.h5ad')
SAMPLE_ID_IN_CONTAINER = 'HC2'
BARCODE_RE = re.compile(r'^[ACGTN]{16}$')
CONTAINER_TOKEN_RE = re.compile(r'[ACGTN]{8,}')


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def decode_categorical(group_or_dataset, obs_group, key):
    """Decode an AnnData obs categorical, refusing to silently misrepresent a
    missing/out-of-range code as a real category.

    `categories[codes]` with a negative `code` (anndata's missing-value marker
    is -1) wraps around via Python/NumPy negative indexing and silently returns
    the LAST category -- which reads as a real value downstream. This counts and
    reports any out-of-range code instead of indexing through it blindly.
    """
    node = group_or_dataset
    if isinstance(node, h5py.Group):
        codes = node['codes'][:]
        cats = [c.decode() if isinstance(c, bytes) else str(c) for c in node['categories'][:]]
    elif '__categories' in obs_group and key in obs_group['__categories']:
        codes = node[:]
        cats = [c.decode() if isinstance(c, bytes) else str(c)
                for c in obs_group['__categories'][key][:]]
    else:
        # Not categorical at all -- plain string/object dataset.
        return np.asarray([v.decode() if isinstance(v, bytes) else str(v) for v in node[:]]), \
            {'categorical': False, 'n_out_of_range_codes': 0, 'n_missing_marker_codes': 0}

    codes = np.asarray(codes)
    n_cats = len(cats)
    missing_mask = codes == -1
    out_of_range_mask = (codes < -1) | (codes >= n_cats)
    bad = missing_mask | out_of_range_mask
    safe_codes = np.where(bad, 0, codes)  # placeholder index, never trusted where bad=True
    decoded = np.asarray(cats, dtype=object)[safe_codes]
    decoded = decoded.astype(object)
    decoded[bad] = None
    return decoded, {
        'categorical': True,
        'n_categories': n_cats,
        'n_missing_marker_codes': int(missing_mask.sum()),
        'n_out_of_range_codes': int(out_of_range_mask.sum()),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()
    if os.path.exists(args.out_dir) and os.listdir(args.out_dir):
        raise SystemExit(f'REFUSING: {args.out_dir} exists and is non-empty -- '
                          'use a new output directory so no prior result is overwritten')
    os.makedirs(args.out_dir, exist_ok=True)

    out = {'capture': CAPTURE}
    paths = {k: os.path.join(DEPOSIT_DIR, f'{CAPTURE}_{k}')
             for k in ('barcodes.tsv.gz', 'features.tsv.gz', 'matrix.mtx.gz')}
    files = {}
    for k, p in paths.items():
        if not os.path.exists(p):
            raise SystemExit(f'input missing: {p}')
        files[k] = {'path': p, 'bytes': os.path.getsize(p), 'sha256': sha256(p)}
    out['files'] = files

    # ---- C1: matrix header + integrality, streamed -------------------------
    with gzip.open(paths['matrix.mtx.gz'], 'rt') as fh:
        banner = fh.readline().strip()
        fh.readline()
        ngene, ncell, nnz = (int(x) for x in fh.readline().split())
        non_integer = 0
        min_val = None
        n_lines = 0
        for line in fh:
            n_lines += 1
            parts = line.split()
            v = float(parts[2])
            if v != int(v):
                non_integer += 1
            if min_val is None or v < min_val:
                min_val = v
    if n_lines != nnz:
        raise SystemExit(f'matrix nnz mismatch: header declares {nnz}, streamed {n_lines} rows')
    out['c1_matrix'] = {
        'banner': banner, 'n_genes_declared': ngene, 'n_cells_declared': ncell,
        'nnz_declared': nnz, 'nnz_observed': n_lines, 'non_integer_values': non_integer,
        'min_value': min_val,
    }

    with gzip.open(paths['features.tsv.gz'], 'rt') as fh:
        ftypes = set()
        n_features = 0
        for line in fh:
            n_features += 1
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3:
                ftypes.add(parts[2])
    out['c1_features'] = {'n_features': n_features, 'feature_types': sorted(ftypes)}

    # ---- C2: barcode namespace ----------------------------------------------
    with gzip.open(paths['barcodes.tsv.gz'], 'rt') as fh:
        bcs = [line.strip() for line in fh]
    lengths = set(len(b.split('-')[0]) for b in bcs)
    suffixes = set(b.split('-')[1] if '-' in b else '' for b in bcs)
    malformed = [b for b in bcs if not BARCODE_RE.match(b.split('-')[0])]
    out['c2_barcodes'] = {
        'n_barcodes': len(bcs), 'n_unique': len(set(bcs)),
        'token_lengths_observed': sorted(lengths), 'gem_suffixes_observed': sorted(suffixes),
        'n_malformed_tokens': len(malformed),
    }

    # ---- Matrix declared dims vs actual deposit record counts ---------------
    # A declared ncell/ngene that silently disagrees with the real barcode/feature
    # file would mean the matrix and the identity files describe different things.
    out['dims_consistency'] = {
        'matrix_ncell_declared': ncell, 'barcodes_file_count': len(bcs),
        'matrix_ncell_matches_barcodes_file': ncell == len(bcs),
        'matrix_ngene_declared': ngene, 'features_file_count': n_features,
        'matrix_ngene_matches_features_file': ngene == n_features,
    }
    if ncell != len(bcs) or ngene != n_features:
        raise SystemExit(
            f'REFUSING: matrix header dims ({ngene} genes x {ncell} cells) disagree with '
            f'the barcode/feature files ({n_features} features, {len(bcs)} barcodes)')

    # ---- C3: container join --------------------------------------------------
    with h5py.File(CONTAINER, 'r') as f:
        obs = f['obs']
        idx = obs['_index'][:]
        idx = np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in idx])
        sid, sid_decode_report = decode_categorical(obs['Sample_ID'], obs, 'Sample_ID')
        all_sample_ids = sorted({v for v in sid.tolist() if v is not None})
        mask = sid == SAMPLE_ID_IN_CONTAINER
        container_cells = idx[mask]
    out['c3_sample_id_decode'] = sid_decode_report
    deposit_tokens = set(b.split('-')[0] for b in bcs)
    # Container cell IDs are "<Sample_ID>_<16nt-barcode>" (e.g. "HC2_AACCACAGTTGTGTAC"),
    # NOT a bare/suffixed barcode -- a split-on-'-' extraction leaves the "HC2_" prefix
    # attached and never matches. Use the same nucleotide-token regex
    # verify_capture_evidence.py uses for GSM5102900_HC1's own C3 join.
    container_tokens = [CONTAINER_TOKEN_RE.search(c).group(0) for c in container_cells]
    n_container_tokens_duplicated = len(container_tokens) - len(set(container_tokens))
    malformed_container_tokens = [t for t in container_tokens if not BARCODE_RE.match(t)]
    found = sum(1 for t in container_tokens if t in deposit_tokens)
    out['c3_container'] = {
        'all_sample_id_values_in_container': all_sample_ids,
        'container_cells_sample_id_HC2': int(mask.sum()),
        'container_tokens_extracted': len(container_tokens),
        'container_tokens_duplicated': n_container_tokens_duplicated,
        'container_tokens_malformed_not_16nt': len(malformed_container_tokens),
        'container_tokens_found_in_deposit': found,
        'container_is_subset_of_deposit': found == len(container_tokens),
        'fraction_of_deposit_retained': found / len(bcs) if bcs else None,
    }
    if sid_decode_report['n_out_of_range_codes'] or sid_decode_report['n_missing_marker_codes']:
        out['c3_container']['WARNING'] = (
            f"{sid_decode_report['n_missing_marker_codes']} missing-marker and "
            f"{sid_decode_report['n_out_of_range_codes']} out-of-range Sample_ID codes were "
            "excluded from the HC2 mask, not silently mapped to a real category")

    out_path = os.path.join(args.out_dir, 'gate2_hc2_candidate_evidence.json')
    with open(out_path, 'w') as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
    print(json.dumps(out, indent=2, sort_keys=True))
    print('wrote', out_path)


if __name__ == '__main__':
    main()
