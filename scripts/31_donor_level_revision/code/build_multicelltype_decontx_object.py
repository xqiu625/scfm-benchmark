#!/usr/bin/env python3
"""
Multi-cell-type DecontX rerun, step 2: combine the EXISTING canonical
platelet cells (148/240-patient subset overlapping Ren/Stephenson/
Qiu-sepsis) with genuine non-platelet-lineage cells from those same three
raw sources, so DecontX has real heterogeneous background to estimate an
ambient/contamination profile from -- directly answering
decontx_audit_review.md point B (platelet-only object under-identifies
contamination), as endorsed by reprocessing_protocol.md Stage 4 ("DecontX
on the full multi-cell-type object -- never on a platelet-only matrix").

Sources of non-platelet-lineage cells (reusing this project's already-
validated loading logic from build_monocyte_canonical.py / platelet_tcell_
ratio.py, generalized from monocyte-only/T-only to a full coarse-lineage
panel):
  Ren (GSE158055): CD4/CD8/B/Mono/NK/Macro/DC/Mega, col-index mapping
    verified empirically against the pre-existing ren_target_monocyte_
    cols.tsv (20,000/20,000 exact match, 2026-09-14).
  Stephenson (E-MTAB-10026): full_clustering categories mapped to coarse
    lineages, excluding Platelets/RBC/HSC_*/technical categories (Stephenson
    DOES have its own curated "Platelets" label -- deliberately NOT used
    here, to keep the platelet arm exactly the same cells already analyzed
    throughout this project, avoiding a second, inconsistent platelet
    source).
  Qiu sepsis (GSE167363): marker-argmax lineage call per cell (same 4-gene-
    set scoring as build_monocyte_canonical.py), keeping T/B/NK/Mono
    winners only.

Platelet cells: reused directly from the already-exported canonical-cohort
counts.mtx/cells.csv (decontx_io/), column-subset to the 148 patients that
overlap these 3 sources -- NOT re-extracted from the raw h5ad, so this is
provably the exact same cells already analyzed in the platelet-only
DecontX run.
"""
import gzip
import os
import tarfile

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

RAW = "/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/monocyte_source"
MAPPING_CSV = ("/bigdata/godziklab/shared/Xinru/302006/302006_updategit/"
               "scripts/31_donor_level_revision/results/donor_meta_with_patient_id.csv")
SAMPLE_MAP_CSV = ("/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/"
                   "31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv")
PATIENT_LIST_CSV = ("/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/"
                     "31_donor_level_revision/results_patient_v2/qc_baseline/patient_qc_distributions.csv")
PLATELET_IO_DIR = "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/bib_revision_tier2/contamination_audit/decontx_io"
OUT_DIR = "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/bib_revision_tier2/contamination_audit/decontx_multicelltype"
os.makedirs(OUT_DIR, exist_ok=True)

CELLS_PER_LINEAGE_PER_PATIENT_STEPH = 6
CELLS_PER_LINEAGE_PER_PATIENT_QIU = 6
RNG_SEED = 20260914

MONO_MARKERS = ["CD14", "LYZ", "FCN1", "S100A8", "S100A9"]
OTHER_LINEAGE = {
    "T": ["CD3D", "CD3E", "CD3G", "TRBC2"],
    "B": ["MS4A1", "CD79A", "CD79B"],
    "NK": ["GNLY", "NKG7", "KLRD1"],
    "Platelet": ["PPBP", "PF4", "ITGA2B"],
}

STEPHENSON_LINEAGE_MAP = {}
for lbl in ["CD4.Naive", "CD4.CM", "CD4.IL22", "CD4.Tfh", "CD4.EM", "CD4.Prolif", "CD4.Th1", "CD4.Th2", "CD4.Th17",
            "CD8.TE", "CD8.Naive", "CD8.EM", "CD8.Prolif", "gdT", "MAIT", "NKT", "Treg"]:
    STEPHENSON_LINEAGE_MAP[lbl] = "T"
for lbl in ["NK_16hi", "NK_56hi", "NK_prolif", "ILC1_3", "ILC2"]:
    STEPHENSON_LINEAGE_MAP[lbl] = "NK"
for lbl in ["B_naive", "B_switched_memory", "B_immature", "B_non-switched_memory", "B_exhausted",
            "Plasmablast", "Plasma_cell_IgG", "Plasma_cell_IgA", "Plasma_cell_IgM", "B_malignant"]:
    STEPHENSON_LINEAGE_MAP[lbl] = "B"
for lbl in ["CD14_mono", "CD83_CD14_mono", "CD16_mono", "C1_CD16_mono", "Mono_prolif"]:
    STEPHENSON_LINEAGE_MAP[lbl] = "Mono"
for lbl in ["pDC", "DC3", "DC2", "DC1", "DC_prolif", "ASDC"]:
    STEPHENSON_LINEAGE_MAP[lbl] = "DC"
# Deliberately excluded: Platelets, RBC, HSC_* (technical/irrelevant for this contrast)


def load_ren_multilineage():
    print("=== Ren multilineage ===", flush=True)
    features = pd.read_csv(f"{RAW}/ren/GSE158055_covid19_features.tsv.gz", header=None)[0].values
    n_genes = len(features)
    cols_df = pd.read_csv(f"{OUT_DIR}/ren_multilineage_cols.tsv", sep="\t").drop_duplicates("col")
    meta_df = pd.read_csv(f"{OUT_DIR}/ren_multilineage_meta.csv")
    n_cells = len(cols_df)
    col_to_local = pd.Series(np.arange(n_cells), index=cols_df["col"].values.astype(np.int32))

    trip = pd.read_csv(f"{OUT_DIR}/ren_multilineage_triplets.txt", sep=" ", header=None,
                        names=["row", "col", "val"],
                        dtype={"row": np.int32, "col": np.int32, "val": np.int32}, engine="c")
    local_c = col_to_local.reindex(trip["col"].values)
    keep = local_c.notna().values
    trip = trip[keep]
    local_c = local_c.values[keep].astype(np.int32)
    mat = sp.coo_matrix((trip["val"].values, (local_c, trip["row"].values - 1)),
                         shape=(n_cells, n_genes)).tocsr()

    obs = meta_df.drop_duplicates("col").set_index("col").reindex(cols_df["col"]).reset_index()
    obs = obs.rename(columns={"sampleID": "Sample_ID", "majorType": "lineage"})
    obs.index = [f"ren_ml_{cn}" for cn in obs["cellName"]]

    meta = pd.read_csv(MAPPING_CSV).set_index("sample_id")
    obs = obs.join(meta[["disease", "y", "patient_id"]], on="Sample_ID")
    obs["source_study"] = "Ren_GSE158055"

    adata = ad.AnnData(X=mat, obs=obs, var=pd.DataFrame(index=features))
    adata.var_names_make_unique()
    print(f"Ren multilineage adata: {adata.shape}, patients={adata.obs.patient_id.nunique()}", flush=True)
    return adata


def load_stephenson_multilineage():
    print("=== Stephenson multilineage ===", flush=True)
    meta = pd.read_csv(MAPPING_CSV).set_index("sample_id")
    target = meta[meta["Data_NO"] == "302004data17"]
    target_sample_ids = set(target.index.astype(str))

    path = f"{RAW}/stephenson/covid_portal_210320_with_raw.h5ad"
    backed = ad.read_h5ad(path, backed="r")
    obs = backed.obs
    lineage = obs["full_clustering"].astype(str).map(STEPHENSON_LINEAGE_MAP)
    in_cohort = obs["sample_id"].astype(str).isin(target_sample_ids)
    mask_candidates = lineage.notna() & in_cohort

    rng = np.random.default_rng(RNG_SEED)
    idx_pool = np.where(mask_candidates.values)[0]
    sub_obs = obs.iloc[idx_pool][["patient_id"]].copy()
    sub_obs["lineage"] = lineage.values[idx_pool]
    keep_local = []
    for (pid, lin), grp in sub_obs.groupby(["patient_id", "lineage"]):
        take = min(len(grp), CELLS_PER_LINEAGE_PER_PATIENT_STEPH)
        chosen = rng.choice(grp.index.to_numpy(), size=take, replace=False)
        keep_local.extend(chosen.tolist())
    keep_local = sorted(set(keep_local))
    print(f"Stephenson multilineage cells selected: {len(keep_local)}", flush=True)

    sub = backed[keep_local].to_memory()
    sub.X = sub.layers["raw"]
    del sub.layers["raw"]

    new_obs = pd.DataFrame(index=sub.obs.index)
    new_obs["Sample_ID"] = sub.obs["sample_id"].astype(str).values
    new_obs["patient_id"] = sub.obs["patient_id"].astype(str).values
    new_obs["lineage"] = obs["full_clustering"].astype(str).map(STEPHENSON_LINEAGE_MAP).loc[sub.obs.index].values
    new_obs = new_obs.join(target[["disease", "y"]], on="Sample_ID")
    new_obs["source_study"] = "Stephenson_EMTAB10026"
    new_obs.index = [f"steph_ml_{i}" for i in sub.obs.index]
    sub.obs = new_obs
    sub.var_names_make_unique()
    print(f"Stephenson multilineage adata: {sub.shape}, patients={sub.obs.patient_id.nunique()}", flush=True)
    return sub


QIU_GSM_TO_SAMPLE = {
    "P25_T0": "NS1_T0", "P25_T6": "NS1_T6",
    "P50_T0": "S1_T0", "P50_T6": "S1_T6",
    "NSES_T0": "NS2_T0", "NSES_T6": "NS2_T6",
    "S2_T0": "S2_T0", "S2_T6": "S2_T6",
    "S3_T0": "S3_T0", "S3_T6": "S3_T6",
}


def score_lineage(X, gene_idx_map, genes):
    idx = [gene_idx_map[g] for g in genes if g in gene_idx_map]
    if not idx:
        return np.zeros(X.shape[0])
    sub = X[:, idx]
    if sp.issparse(sub):
        sub = sub.toarray()
    return sub.mean(axis=1).ravel()


def load_qiu_multilineage(meta):
    print("=== Qiu sepsis multilineage ===", flush=True)
    d = f"{RAW}/qiu_sepsis"
    extracted = f"{d}/extracted"
    files = os.listdir(extracted)
    rng = np.random.default_rng(RNG_SEED)
    adatas = []
    for gsm_name, sample_id in QIU_GSM_TO_SAMPLE.items():
        mtx_f = [x for x in files if gsm_name in x and x.endswith("matrix.mtx.gz")]
        if not mtx_f:
            continue
        prefix = mtx_f[0].replace("matrix.mtx.gz", "")
        bc = pd.read_csv(f"{extracted}/{prefix}barcodes.tsv.gz", header=None)[0].values
        feat = pd.read_csv(f"{extracted}/{prefix}features.tsv.gz", header=None, sep="\t")
        gene_symbols = feat[1].values if feat.shape[1] > 1 else feat[0].values
        with gzip.open(f"{extracted}/{prefix}matrix.mtx.gz", "rb") as fh:
            mat = sio.mmread(fh).tocsr()
        if mat.shape[0] == len(gene_symbols):
            mat = mat.T.tocsr()
        gidx = {g: i for i, g in enumerate(gene_symbols)}
        scores = {k: score_lineage(mat, gidx, v) for k, v in
                   {**{"Mono": MONO_MARKERS}, **OTHER_LINEAGE}.items()}
        score_df = pd.DataFrame(scores)
        winner = score_df.idxmax(axis=1)
        has_signal = score_df.max(axis=1) > 0
        keep_mask = has_signal & winner.isin(["Mono", "T", "B", "NK"])
        keep_idx = np.where(keep_mask.values)[0]
        if len(keep_idx) == 0:
            continue
        take = min(len(keep_idx), CELLS_PER_LINEAGE_PER_PATIENT_QIU * 4)
        chosen = rng.choice(keep_idx, size=take, replace=False)
        sub = mat[chosen, :]
        obs = pd.DataFrame({"Sample_ID": sample_id, "lineage": winner.values[chosen]},
                            index=[f"qiu_ml_{gsm_name}_{i}" for i in chosen])
        one = ad.AnnData(X=sub, obs=obs, var=pd.DataFrame(index=gene_symbols))
        one.var_names_make_unique()
        adatas.append(one)
        print(f"  {gsm_name}: {len(chosen)} multilineage cells kept", flush=True)

    if not adatas:
        return None
    adata = ad.concat(adatas, join="outer", fill_value=0)
    adata.obs = adata.obs.join(meta[["disease", "y", "patient_id"]], on="Sample_ID")
    adata.obs["source_study"] = "Qiu_sepsis_GSE167363"
    adata.var_names_make_unique()
    print(f"Qiu multilineage adata: {adata.shape}, patients={adata.obs.patient_id.nunique()}", flush=True)
    return adata


def load_platelet_subset(target_patient_ids):
    print("=== Platelet cells (reused from existing platelet-only export) ===", flush=True)
    genes = pd.read_csv(f"{PLATELET_IO_DIR}/genes.csv")["gene"].tolist()
    cells = pd.read_csv(f"{PLATELET_IO_DIR}/cells.csv")
    raw = sio.mmread(f"{PLATELET_IO_DIR}/counts.mtx")
    raw_gc = (raw.tocsc() if raw.shape[0] == len(genes) else raw.T.tocsc())

    keep_mask = cells["patient_id"].astype(str).isin(target_patient_ids).to_numpy()
    print(f"Platelet cells in the 148-patient subset: {keep_mask.sum()} / {len(cells)}", flush=True)
    sub_mat = raw_gc[:, keep_mask].T.tocsr()  # cells x genes
    sub_obs = cells.loc[keep_mask].copy()
    sub_obs["lineage"] = "Platelet"
    sub_obs["source_study"] = "platelet_canonical_object"
    sub_obs.index = [f"plt_{i}" for i in sub_obs["cell_id"]]

    adata = ad.AnnData(X=sub_mat, obs=sub_obs.drop(columns=["cell_id"]), var=pd.DataFrame(index=genes))
    print(f"Platelet subset adata: {adata.shape}, patients={adata.obs.patient_id.nunique()}", flush=True)
    return adata


if __name__ == "__main__":
    patient_list = pd.read_csv(PATIENT_LIST_CSV, dtype={"patient_id": str})
    patient_list = patient_list[patient_list["disease"] == "COVID-19"]
    canonical_240 = set(patient_list["patient_id"].astype(str))
    meta = pd.read_csv(MAPPING_CSV).set_index("sample_id")

    ren = load_ren_multilineage()
    qiu = load_qiu_multilineage(meta)
    steph = load_stephenson_multilineage()

    multilineage_patients = set(ren.obs.patient_id.astype(str))
    if qiu is not None:
        multilineage_patients |= set(qiu.obs.patient_id.astype(str))
    multilineage_patients |= set(steph.obs.patient_id.astype(str))
    target_148 = multilineage_patients & canonical_240
    print(f"\nPatients with real multi-lineage data AND in the canonical 240-cohort: "
          f"{len(target_148)}", flush=True)

    platelet_adata = load_platelet_subset(target_148)

    ren = ren[ren.obs.patient_id.astype(str).isin(target_148)].copy()
    steph = steph[steph.obs.patient_id.astype(str).isin(target_148)].copy()
    if qiu is not None:
        qiu = qiu[qiu.obs.patient_id.astype(str).isin(target_148)].copy()
    print(f"After restricting to target_148: Ren={ren.n_obs}, Stephenson={steph.n_obs}, "
          f"Qiu={qiu.n_obs if qiu is not None else 0}, Platelet={platelet_adata.n_obs}", flush=True)

    parts = [a for a in (ren, qiu, steph, platelet_adata) if a is not None]
    combined = ad.concat(parts, join="inner", fill_value=0)
    combined.obs["is_platelet"] = (combined.obs["lineage"] == "Platelet")
    combined.write(f"{OUT_DIR}/multicelltype_raw.h5ad")

    print("\n=== Multi-cell-type DecontX object ===", flush=True)
    print(f"Total: {combined.shape[0]} cells x {combined.shape[1]} genes (gene intersection)", flush=True)
    print(combined.obs.groupby(["source_study", "lineage"]).size(), flush=True)
    print(f"TOTAL patients: {combined.obs.patient_id.nunique()}", flush=True)

    X = combined.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    X = X.tocsc().astype(np.float64)
    sio.mmwrite(f"{OUT_DIR}/multicelltype_counts.mtx", X.T.tocsr())
    pd.Series(combined.var_names.astype(str)).to_csv(f"{OUT_DIR}/multicelltype_genes.csv", index=False, header=["gene"])
    combined.obs[["patient_id", "Sample_ID", "lineage", "is_platelet", "source_study", "y"]].to_csv(
        f"{OUT_DIR}/multicelltype_cells.csv", index_label="cell_id")
    print("SAVED multicelltype_counts.mtx, multicelltype_genes.csv, multicelltype_cells.csv", flush=True)
