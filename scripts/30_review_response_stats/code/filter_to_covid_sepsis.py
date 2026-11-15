#!/usr/bin/env python3
"""
Build a COVID+Sepsis-matched embedding directory (roadmap D).

The auxiliary tables (multiclassifier, embedding quality, activation, ablation, fusion)
currently run on the 5-disease pool (47,977 cells), while the headline tables use the
disease-matched COVID+Sepsis subset (42,095). Reviewers flagged cross-cohort comparison
as a correctness issue. This filters each FM's canonical h5ad to COVID+Sepsis and writes
a parallel embedding dir so the SAME auxiliary scripts produce subset-matched tables.

Preserves X (raw genes, needed for activation/ablation PCA) + obsm + obs.

In:  02_EMBEDDINGS/option_c/<fm>_covid_canonical.h5ad
Out: 02_EMBEDDINGS/option_c_covid_sepsis/<fm>_covid_canonical.h5ad
"""
import os
import glob
import anndata as ad

SRC = '/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c'
DST = '/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c_covid_sepsis'
KEEP_DISEASE = ['COVID-19', 'Sepsis']
KEEP_SEVERITY = ['fatal', 'mild', 'moderate', 'recovered', 'severe']


def main():
    os.makedirs(DST, exist_ok=True)
    files = sorted(glob.glob(os.path.join(SRC, '*_covid_canonical.h5ad')))
    print(f'{len(files)} FM files', flush=True)
    for f in files:
        name = os.path.basename(f)
        out = os.path.join(DST, name)
        a = ad.read_h5ad(f)
        n0 = a.n_obs
        mask = a.obs['Disease'].isin(KEEP_DISEASE) & a.obs['covid_severity'].isin(KEEP_SEVERITY)
        a = a[mask.values].copy()
        a.write_h5ad(out)
        print(f'{name}: {n0} -> {a.n_obs} cells  (donors={a.obs["Sample_ID"].nunique()})', flush=True)
    print(f'SAVED -> {DST}', flush=True)


if __name__ == '__main__':
    main()
