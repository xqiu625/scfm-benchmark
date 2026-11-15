#!/usr/bin/env python3
"""
subpopulation_biomarker_design.md section 7 ("the one thing to do immediately"):
a platelet-to-T-cell composition ratio, pre-specified by Qiu et al. 2024
(qiu2024deciphering, published AUC 0.754 for this exact ratio in a different,
multi-disease sample set) as a single comparator -- not a scan, not cluster
discovery. Tested here under this project's own patient-level protocol on
the SAME 3-cohort scope as the existing monocyte positive control (R4:
Ren + Stephenson + Qiu-sepsis, curated/markered cell-type labels), since
these are the only cohorts in this project with non-platelet cell-type
annotations available at all.

T-cell counts per patient (this script computes only COUNTS, not full
expression matrices -- much cheaper than build_monocyte_canonical.py, which
this reuses the annotation-reading logic from):
  Ren (GSE158055): majorType in {CD4, CD8} (curated atlas label; full raw
    annotation CSV read directly, no need to touch the 8GB count matrix).
  Stephenson (E-MTAB-10026): full_clustering in T_CELL_LABELS_STEPHENSON
    (curated atlas label; .obs read from the backed h5ad, no expression
    data loaded).
  Qiu sepsis (GSE167363): marker-threshold call (OTHER_LINEAGE["T"] score
    wins over MONO_MARKERS and the other three lineages) -- same method
    build_monocyte_canonical.py uses for its monocyte call, mirrored here.

Platelet counts per patient reuse the ALREADY-ESTABLISHED per-patient cell
counts from this project's main platelet pipeline (no new extraction
needed) -- see PATIENT_QC_CSV below (results_patient_v2/qc_baseline/
patient_qc_distributions.csv, n_cells column), the same file used
throughout this thread for the canonical-cohort filter.

Ratio definition: log2(n_platelet / n_T) per patient (raw ratio is heavily
right-skewed with small-N patients; log form matches how the analogous
count-based null model treats cell counts elsewhere in this project, see
Methods' "null model using only... cell count"). Fit through the SAME
repeated_cv.py harness (single-feature logistic regression per fold,
patient-level rows) as every other baseline-ladder rung this thread has
built.

CAVEAT flagged explicitly in the output: this is restricted to the
148-patient monocyte-control cohort (Ren+Stephenson+Qiu-sepsis), NOT the
full 240-patient within-COVID cohort every other baseline-ladder number in
this thread uses -- do not directly compare this AUC to PCA-200's 0.742
without noting the different N and patient set.
"""
import gzip
import os
import tarfile

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
from sklearn.linear_model import LogisticRegression

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from repeated_cv import DEFAULT_SEEDS, run_repeated_kfold

RAW = "/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/monocyte_source"
MAPPING_CSV = ("/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/"
               "31_donor_level_revision/results/donor_meta_with_patient_id.csv")
PATIENT_QC_CSV = ("/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/"
                   "31_donor_level_revision/results_patient_v2/qc_baseline/"
                   "patient_qc_distributions.csv")
OUT_DIR = "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/bib_revision_tier2/pla_t_ratio"
os.makedirs(OUT_DIR, exist_ok=True)

MONO_MARKERS = ["CD14", "LYZ", "FCN1", "S100A8", "S100A9"]
OTHER_LINEAGE = {
    "T": ["CD3D", "CD3E", "CD3G", "TRBC2"],
    "B": ["MS4A1", "CD79A", "CD79B"],
    "NK": ["GNLY", "NKG7", "KLRD1"],
    "Platelet": ["PPBP", "PF4", "ITGA2B"],
}

# VERIFIED against the actual full_clustering value_counts (2026-09-13,
# ad.read_h5ad(..., backed="r").obs["full_clustering"].value_counts()) --
# every label below is confirmed present. Deliberately EXCLUDES "NKT"
# (hybrid NK/T lineage, kept as its own category in this atlas rather than
# folded into T, so excluded here for a clean/conservative T-cell
# definition matching the CD3D/CD3E/CD3G/TRBC2 marker panel used elsewhere
# in this project) and "Platelets" (present in this atlas as its own label,
# but NOT used here -- platelet counts are taken from this project's
# already-established main pipeline instead, for consistency with every
# other platelet count in this thread; see main() below).
T_LABELS_STEPHENSON = [
    "CD4.Naive", "CD4.CM", "CD4.EM", "CD4.IL22", "CD4.Tfh", "CD4.Th1",
    "CD4.Th17", "CD4.Th2", "Treg", "CD4.Prolif",
    "CD8.Naive", "CD8.EM", "CD8.TE", "CD8.Prolif",
    "MAIT", "gdT",
]


def load_meta():
    df = pd.read_csv(MAPPING_CSV)
    return df.set_index("sample_id")


def score_lineage(X, gene_idx_map, genes):
    idx = [gene_idx_map[g] for g in genes if g in gene_idx_map]
    if not idx:
        return np.zeros(X.shape[0])
    sub = X[:, idx]
    if sp.issparse(sub):
        sub = sub.toarray()
    return sub.mean(axis=1).ravel()


def ren_tcell_counts(meta):
    print("=== Ren GSE158055: T-cell counts (annotation-only, no matrix read) ===", flush=True)
    target = meta[meta["Data_NO"] == "302004data15"]
    anno = pd.read_csv(f"{RAW}/ren/GSE158055_cell_annotation.csv.gz")
    is_t = anno["majorType"].isin(["CD4", "CD8"])
    anno_t = anno[is_t]
    joined = anno_t.join(target[["patient_id"]], on="sampleID")
    joined = joined.dropna(subset=["patient_id"])
    counts = joined.groupby("patient_id").size()
    print(f"  {is_t.sum()} T cells total in atlas, {len(joined)} in target cohort, "
          f"{counts.shape[0]} patients", flush=True)
    return counts


def stephenson_tcell_counts(meta):
    print("=== Stephenson E-MTAB-10026: T-cell counts (.obs only, no expression read) ===", flush=True)
    import anndata as ad
    target = meta[meta["Data_NO"] == "302004data17"]
    path = f"{RAW}/stephenson/covid_portal_210320_with_raw.h5ad"
    backed = ad.read_h5ad(path, backed="r")
    obs = backed.obs
    is_t = obs["full_clustering"].astype(str).isin(T_LABELS_STEPHENSON)
    target_sample_ids = set(target.index)
    in_cohort = obs["sample_id"].astype(str).isin(target_sample_ids)
    mask = is_t & in_cohort
    print(f"  {int(is_t.sum())} T cells total in atlas, {int(mask.sum())} in target cohort", flush=True)
    sub_obs = obs[mask]
    counts = sub_obs.groupby(sub_obs["patient_id"].astype(str)).size()
    print(f"  {counts.shape[0]} patients", flush=True)
    return counts


def qiu_sepsis_tcell_counts(meta):
    print("=== Qiu sepsis GSE167363: T-cell counts (marker scoring) ===", flush=True)
    QIU_GSM_TO_SAMPLE = {
        "P25_T0": "NS1_T0", "P25_T6": "NS1_T6",
        "P50_T0": "S1_T0", "P50_T6": "S1_T6",
        "NSES_T0": "NS2_T0", "NSES_T6": "NS2_T6",
        "S2_T0": "S2_T0", "S2_T6": "S2_T6",
        "S3_T0": "S3_T0", "S3_T6": "S3_T6",
    }
    d = f"{RAW}/qiu_sepsis"
    extracted = f"{d}/extracted"
    os.makedirs(extracted, exist_ok=True)
    if not os.listdir(extracted):
        with tarfile.open(f"{d}/GSE167363_RAW.tar") as t:
            t.extractall(extracted)
    files = os.listdir(extracted)
    counts_by_sample = {}
    for gsm_name, sample_id in QIU_GSM_TO_SAMPLE.items():
        mtx_f = [x for x in files if gsm_name in x and x.endswith("matrix.mtx.gz")]
        if not mtx_f:
            continue
        prefix = mtx_f[0].replace("matrix.mtx.gz", "")
        feat = pd.read_csv(f"{extracted}/{prefix}features.tsv.gz", header=None, sep="\t")
        gene_symbols = feat[1].values if feat.shape[1] > 1 else feat[0].values
        with gzip.open(f"{extracted}/{prefix}matrix.mtx.gz", "rb") as fh:
            mat = sio.mmread(fh).tocsr()
        bc = pd.read_csv(f"{extracted}/{prefix}barcodes.tsv.gz", header=None)[0].values
        if mat.shape[0] == len(gene_symbols):
            mat = mat.T.tocsr()
        gidx = {g: i for i, g in enumerate(gene_symbols)}
        scores = {k: score_lineage(mat, gidx, v) for k, v in
                  {**{"Mono": MONO_MARKERS}, **OTHER_LINEAGE}.items()}
        is_t = (scores["T"] > 0)
        for k, s in scores.items():
            if k != "T":
                is_t &= (scores["T"] >= s)
        counts_by_sample[sample_id] = int(is_t.sum())
        print(f"  {gsm_name} -> {sample_id}: {mat.shape[0]} cells, {int(is_t.sum())} marker-T", flush=True)

    s = pd.Series(counts_by_sample, name="n_t")
    joined = s.to_frame().join(meta[["patient_id"]], how="inner")
    counts = joined.groupby("patient_id")["n_t"].sum()
    print(f"  {counts.shape[0]} patients", flush=True)
    return counts


def main():
    meta = load_meta()
    ren_t = ren_tcell_counts(meta)
    steph_t = stephenson_tcell_counts(meta)
    qiu_t = qiu_sepsis_tcell_counts(meta)
    t_counts = pd.concat([ren_t, steph_t, qiu_t])
    t_counts = t_counts.groupby(t_counts.index).sum()
    print(f"\nTotal T-cell counts: {len(t_counts)} patients", flush=True)

    qc = pd.read_csv(PATIENT_QC_CSV, dtype={"patient_id": str})
    qc = qc.set_index("patient_id")
    t_counts.index = t_counts.index.astype(str)

    common = t_counts.index.intersection(qc.index)
    print(f"Patients with BOTH platelet and T-cell counts: {len(common)}", flush=True)
    missing = t_counts.index.difference(qc.index)
    if len(missing):
        print(f"  {len(missing)} patients have T-cell counts but no platelet QC row "
              f"(not in the canonical dataset) -- dropped: {list(missing)[:10]}", flush=True)

    df = pd.DataFrame({
        "n_platelet": qc.loc[common, "n_cells"],
        "n_t": t_counts.loc[common],
        "y": qc.loc[common, "y"],
        "disease": qc.loc[common, "disease"],
    })
    df = df[(df["n_platelet"] > 0) & (df["n_t"] > 0)]
    df["log_ratio"] = np.log2(df["n_platelet"] / df["n_t"])
    print(f"\nFinal analysis set: {len(df)} patients, {int(df['y'].sum())} fatal", flush=True)
    print(df.groupby("disease").size(), flush=True)

    out_path = os.path.join(OUT_DIR, "pla_t_ratio_per_patient.csv")
    df.to_csv(out_path)
    print(f"SAVED {out_path}", flush=True)

    if df["y"].nunique() < 2 or len(df) < 20:
        print("Too few patients / only one class -- stopping before CV fit.", flush=True)
        return

    X = df[["log_ratio"]].to_numpy()
    y = df["y"].to_numpy()
    groups = df.index.to_numpy()

    def fit_predict(X_tr, y_tr, X_te):
        clf = LogisticRegression(C=1.0, solver="lbfgs", max_iter=2000, class_weight="balanced")
        clf.fit(X_tr, y_tr)
        p = clf.predict_proba(X_te)[:, 1]
        return {"p_lr": p, "p_rf": p}

    seeds = DEFAULT_SEEDS
    oof_by_seed = run_repeated_kfold(X, y, groups, fit_predict_fn=fit_predict, seeds=seeds)

    from repeated_cv import repeated_patient_auc, summarize_repeated
    aucs, _ = repeated_patient_auc(oof_by_seed, groups, y, score_key="p_lr", agg="mean")
    summ = summarize_repeated(aucs)
    print(f"\nPla-T log-ratio, patient-level repeated CV (n={len(df)}, {int(y.sum())} fatal): "
          f"repeat-mean AUC = {summ['mean']:.4f} +/- {summ['sd']:.4f} (n_seeds={summ['n']})", flush=True)

    for seed in seeds:
        out = pd.DataFrame({"patient_id": groups, "y": y, **oof_by_seed[seed]})
        out.to_csv(os.path.join(OUT_DIR, f"oof_patient_PlaTRatio_seed{seed}.csv"), index=False)
    print(f"SAVED {len(seeds)} seed files -> {OUT_DIR}/oof_patient_PlaTRatio_seed<42..51>.csv", flush=True)


if __name__ == "__main__":
    main()
