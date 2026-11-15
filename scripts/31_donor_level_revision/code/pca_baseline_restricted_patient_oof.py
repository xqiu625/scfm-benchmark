#!/usr/bin/env python3
"""
contamination_response_plan.md item 4: restricted-feature refit of the PCA-200
baseline -- "the one that answers the reviewer." Removes a standard,
pre-specified non-platelet-lineage marker panel (T/B/NK/monocyte-myeloid/DC/
erythroid/pan-leukocyte) from the gene space BEFORE PCA, refits under the
exact same protocol as pca_baseline_patient_oof_v2.py (same disease/severity
filter, same patient map, same CV harness, same seeds, same N_COMPONENTS),
and reports the ensemble-OOF AUC alongside the unrestricted PCA-200 (0.7423
per Table 3) via baseline_ladder_delong.py's comparison machinery.

If the headline AUC survives feature restriction, the outcome claim is
contamination-independent -- the whole objection is contained in one row.

Marker panel (excluded, ~46 genes): standard curated PBMC lineage markers for
T, B, NK, monocyte/myeloid, DC, erythroid, plus the pan-leukocyte marker
PTPRC. Deliberately does NOT exclude generic housekeeping/ribosomal genes
(RPLP2/RPL24/RPL18A etc, flagged separately in the R6 audit as ambient-like)
since those are not LINEAGE markers -- removing them would conflate two
different questions (does the AUC depend on cross-lineage content, vs. does
it depend on generically-abundant transcripts). Also does NOT exclude
platelet-native genes (PPBP, PF4, ITGA2B, GP1BA, CCL5, etc) -- those stay in.
"""
import argparse
import os

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from repeated_cv import DEFAULT_SEEDS, run_repeated_kfold

RAW = ('/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/'
       'human_platelet_covid_severity_raw_counts.h5ad')
PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2'

KEEP_DISEASE = ['COVID-19', 'Sepsis']
KEEP_SEVERITY = ['fatal', 'mild', 'moderate', 'recovered', 'severe']
N_COMPONENTS = 200

SEVERITY_FIX = {
    'GSM4995436': 'moderate',
    'GSM4995442': 'severe',
    'GSM4995432': 'moderate',
}

NON_PLATELET_LINEAGE_MARKERS = {
    "T_cell": ["CD3D", "CD3E", "CD3G", "CD2", "TRBC1", "TRBC2", "IL7R", "GZMA", "GZMB", "GZMK",
               "GZMH", "CD8A", "CD8B", "CD4", "CCL5", "IL32", "CD7"],
    "B_cell": ["CD79A", "CD79B", "MS4A1", "CD19", "IGHM", "IGHD", "CD74"],
    "NK_cell": ["GNLY", "NKG7", "KLRD1", "KLRB1", "PRF1", "NCAM1", "KLRF1"],
    "Monocyte_Myeloid": ["CD14", "LYZ", "FCN1", "S100A8", "S100A9", "S100A12", "FCER1G", "TYROBP",
                          "AIF1", "MNDA", "C1QA", "C1QB", "C1QC", "AXL", "VCAN"],
    "DC": ["FCER1A", "CLEC10A", "ITGAX", "CD1C"],
    "Erythroid": ["HBB", "HBA1", "HBA2", "SLC25A37", "ALAS2"],
    "Pan_leukocyte": ["PTPRC"],
}
EXCLUDE_GENES = sorted({g for genes in NON_PLATELET_LINEAGE_MARKERS.values() for g in genes})


def fit_predict(X_tr, y_tr, X_te):
    scaler = StandardScaler(with_mean=False).fit(X_tr)
    Xs_tr = scaler.transform(X_tr)
    Xs_te = scaler.transform(X_te)
    pca = PCA(n_components=N_COMPONENTS, random_state=42)
    Z_tr = pca.fit_transform(Xs_tr.toarray() if hasattr(Xs_tr, 'toarray') else Xs_tr)
    Z_te = pca.transform(Xs_te.toarray() if hasattr(Xs_te, 'toarray') else Xs_te)
    clf = LogisticRegression(C=1.0, solver='lbfgs', max_iter=2000,
                              class_weight='balanced', n_jobs=4)
    clf.fit(Z_tr, y_tr)
    p = clf.predict_proba(Z_te)[:, 1]
    return {'p_lr': p, 'p_rf': p}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n-seeds', type=int, default=1)
    args = ap.parse_args()
    seeds = DEFAULT_SEEDS[:args.n_seeds] if args.n_seeds > 1 else [42]

    os.makedirs(OUT_DIR, exist_ok=True)
    print('Loading raw counts...', flush=True)
    a = sc.read_h5ad(RAW)
    print(f'  {a.n_obs} cells x {a.n_vars} genes', flush=True)

    var_names = a.var_names.astype(str)
    exclude_mask = var_names.isin(EXCLUDE_GENES)
    found_excluded = sorted(set(var_names[exclude_mask]))
    print(f'  Excluding {exclude_mask.sum()} / {len(EXCLUDE_GENES)} requested non-platelet-lineage '
          f'markers found in the gene panel: {found_excluded}', flush=True)
    a = a[:, ~exclude_mask].copy()
    print(f'  After exclusion: {a.n_vars} genes remain', flush=True)

    disease = a.obs['Disease'].astype(str).values
    severity = a.obs['covid_severity'].astype(str).values
    sample_id = a.obs['Sample_ID'].astype(str).values
    cell_id = a.obs_names.values.astype(str)

    for sid, new_sev in SEVERITY_FIX.items():
        severity[sample_id == sid] = new_sev

    keep = np.isin(disease, KEEP_DISEASE) & np.isin(severity, KEEP_SEVERITY)
    print(f'  Disease+severity filter: {keep.sum()} / {a.n_obs} cells', flush=True)

    pmap = pd.read_csv(PATIENT_MAP)[['sample_id', 'patient_id']]
    df = pd.DataFrame({'cell_id': cell_id[keep], 'sample_id': sample_id[keep],
                        'severity': severity[keep]})
    df['_row'] = np.arange(len(df))
    n_before = len(df)
    df = df.merge(pmap, on='sample_id', how='inner')
    print(f'  v2 map filter: {len(df)} / {n_before} cells', flush=True)

    X = a.X[keep][df['_row'].to_numpy()]
    y = (df['severity'].to_numpy() == 'fatal').astype(int)
    groups = df['patient_id'].to_numpy()
    n_patients = len(np.unique(groups))
    print(f'  {X.shape[0]} cells, {n_patients} patients, {int(y.sum())} fatal', flush=True)

    oof_by_seed = run_repeated_kfold(X, y, groups, fit_predict_fn=fit_predict, seeds=seeds)

    base = dict(cell_id=df['cell_id'].to_numpy(), sample_id=df['sample_id'].to_numpy(),
                severity=df['severity'].to_numpy(), y=y)

    if args.n_seeds <= 1:
        out = pd.DataFrame({**base, **oof_by_seed[seeds[0]]})
        out_path = os.path.join(OUT_DIR, 'oof_PCA-200-restricted.csv')
        out.to_csv(out_path, index=False)
        print(f'SAVED -> {out_path}', flush=True)
    else:
        for seed in seeds:
            out = pd.DataFrame({**base, **oof_by_seed[seed]})
            out_path = os.path.join(OUT_DIR, f'oof_PCA-200-restricted_seed{seed}.csv')
            out.to_csv(out_path, index=False)
        print(f'SAVED {len(seeds)} seed files -> '
              f'{OUT_DIR}/oof_PCA-200-restricted_seed<42..{seeds[-1]}>.csv', flush=True)


if __name__ == '__main__':
    main()
