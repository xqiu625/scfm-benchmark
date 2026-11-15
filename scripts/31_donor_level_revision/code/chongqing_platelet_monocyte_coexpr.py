#!/usr/bin/env python3
"""
Fix-4 (code review 2026-09-04): platelet-monocyte aggregate/co-expression
proportion, per the sepsis_arm_expansion.md design doc's hard requirement --
not computed in the original integrate_chongqing_sepsis.py run. Added
post-hoc on the already-saved platelet subset (no need to redo the ~30 min
full integration+clustering run).

Definition: within the platelet cluster cells, the fraction whose raw counts
show monocyte-marker expression (CD14/LST1/FCN1, >=1 count in at least one of
the three) -- a standard proxy for platelet-monocyte aggregates/doublets in
droplet scRNA-seq (matches the design doc's "marker co-expression" framing,
same category of check as the whole-cohort PPBP/PF4/etc. platelet score,
just the reverse direction).

Merges a `pct_monocyte_coexpr` column into
chongqing_platelet_counts_per_donor.csv (does not overwrite other columns).
"""
import numpy as np
import pandas as pd
import scanpy as sc

PLATELET_H5AD = ('/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/'
                  'chongqing_sepsis_platelet_raw_counts.h5ad')
COUNTS_CSV = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
              '31_donor_level_revision/results_patient/'
              'chongqing_platelet_counts_per_donor.csv')

MONOCYTE_MARKERS = ['CD14', 'LST1', 'FCN1']


def main():
    print(f'Loading {PLATELET_H5AD}...', flush=True)
    a = sc.read_h5ad(PLATELET_H5AD)
    print(f'  {a.n_obs} cells x {a.n_vars} genes', flush=True)

    symbol_col = 'gene_symbol' if 'gene_symbol' in a.var.columns else None
    if symbol_col is None:
        raise SystemExit('ERROR: no gene_symbol column in var -- cannot look up markers')

    present = [m for m in MONOCYTE_MARKERS if m in a.var[symbol_col].to_numpy()]
    print(f'  monocyte markers present: {present} ({len(present)}/{len(MONOCYTE_MARKERS)})',
          flush=True)
    if not present:
        raise SystemExit('ERROR: none of CD14/LST1/FCN1 present in this gene panel')

    gene_idx = [np.where(a.var[symbol_col].to_numpy() == m)[0][0] for m in present]
    X_markers = a.X[:, gene_idx]
    if hasattr(X_markers, 'toarray'):
        X_markers = X_markers.toarray()
    has_coexpr = (X_markers >= 1).any(axis=1)

    donor = a.obs['Sample_ID'].to_numpy()
    df = pd.DataFrame({'sample_id': donor, 'has_coexpr': has_coexpr})
    per_donor = df.groupby('sample_id').agg(
        n_platelet_cells=('has_coexpr', 'size'),
        n_monocyte_coexpr=('has_coexpr', 'sum'),
    ).reset_index()
    per_donor['pct_monocyte_coexpr'] = (100 * per_donor['n_monocyte_coexpr'] /
                                         per_donor['n_platelet_cells'])

    overall_pct = 100 * has_coexpr.mean()
    print(f'\nOverall platelet-monocyte co-expression: {has_coexpr.sum()}/{len(has_coexpr)} '
          f'({overall_pct:.1f}%)', flush=True)

    existing = pd.read_csv(COUNTS_CSV)
    assert 'sample_id' in existing.columns or 'Sample_ID' in existing.columns, \
        f'unexpected columns in {COUNTS_CSV}: {existing.columns.tolist()}'
    id_col = 'sample_id' if 'sample_id' in existing.columns else 'Sample_ID'
    merged = existing.merge(
        per_donor[['sample_id', 'n_monocyte_coexpr', 'pct_monocyte_coexpr']],
        left_on=id_col, right_on='sample_id', how='left', suffixes=('', '_dup'))
    if 'sample_id_dup' in merged.columns:
        merged = merged.drop(columns=['sample_id_dup'])
    merged.to_csv(COUNTS_CSV, index=False)
    print(f'SAVED (merged pct_monocyte_coexpr) -> {COUNTS_CSV}', flush=True)


if __name__ == '__main__':
    main()
