#!/usr/bin/env python3
"""Job F Phase 2: build the monocyte canonical object from Ren, Stephenson, and
Qiu-sepsis source data (the platelet cohort's 3 largest/most useful contributing
studies with either curated cell-type annotations or a tractable path to one).

Provenance, flagged explicitly per-cell via obs['monocyte_source_provenance']:
  'curated_atlas_label'   Ren (majorType=='Mono') and Stephenson (full_clustering
                          monocyte categories) -- labels come from the original
                          published atlas annotation.
  'marker_gene_threshold' Qiu sepsis has NO cell-type annotation on GEO; monocytes
                          are identified by scoring CD14/LYZ/FCN1/S100A8/S100A9
                          against competing lineage marker sets (T/B/NK/Platelet)
                          per cell and calling monocyte where the myeloid score
                          wins and is positive. This is NOT the same standard of
                          evidence as the curated labels -- do not blend silently
                          in any report referencing monocyte counts by source.

Patient identity: reuses Job H's donor_meta_with_patient_id.csv patient_id for
Ren and Qiu-sepsis (validated correct -- draw-suffix stripping). For Stephenson,
OVERRIDES with the SDRF Characteristics[individual] (CVxxxx) for BGCV-coded
Sample_IDs, because Job H's patient_id used the sequencing POOL id (BGCVxx) for
these -- a confirmed bug (flagged separately to the coordinator, out of scope to
fix here). Non-BGCV Stephenson rows (AP*, MH*, newcastle*) reuse patient_id as-is.
"""
import gzip
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.io as sio
import anndata as ad

RAW = "/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/monocyte_source"
OUT_DIR = "/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c_monocyte"
MAPPING_CSV = ("/bigdata/godziklab/shared/Xinru/302006/302006_updategit/"
               "scripts/31_donor_level_revision/results/donor_meta_with_patient_id.csv")
SDRF = f"{RAW}/stephenson/E-MTAB-10026.sdrf.txt"

MONO_MARKERS = ["CD14", "LYZ", "FCN1", "S100A8", "S100A9"]
OTHER_LINEAGE = {
    "T": ["CD3D", "CD3E", "CD3G", "TRBC2"],
    "B": ["MS4A1", "CD79A", "CD79B"],
    "NK": ["GNLY", "NKG7", "KLRD1"],
    "Platelet": ["PPBP", "PF4", "ITGA2B"],
}


def load_meta():
    df = pd.read_csv(MAPPING_CSV)
    return df.set_index("sample_id")


def load_stephenson_individual_map():
    sdrf = pd.read_csv(SDRF, sep="\t")
    sdrf = sdrf[["Source Name", "Characteristics[individual]", "Characteristics[sequence pool]"]]
    sdrf.columns = ["source_name", "individual", "pool"]
    return sdrf.set_index("source_name")


# --------------------------------------------------------------------- Ren
def build_ren(meta):
    print("=== Ren GSE158055 ===", flush=True)
    target = meta[meta["Data_NO"] == "302004data15"]

    features = pd.read_csv(f"{RAW}/ren/GSE158055_covid19_features.tsv.gz", header=None)[0].values
    n_genes = len(features)
    print(f"n_genes={n_genes}")

    cols_df = pd.read_csv(f"{RAW}/ren/ren_target_monocyte_cols.tsv", sep="\t")
    # cols_df: col (1-based global cell index), cellName
    cols_df = cols_df.drop_duplicates("col").reset_index(drop=True)
    n_cells = len(cols_df)
    print(f"n_target_monocyte_cells={n_cells}")
    col_to_local = pd.Series(np.arange(n_cells), index=cols_df["col"].values.astype(np.int32))

    trip = pd.read_csv(f"{RAW}/ren/ren_monocyte_triplets.txt", sep=" ", header=None,
                        names=["row", "col", "val"],
                        dtype={"row": np.int32, "col": np.int32, "val": np.int32},
                        engine="c")
    local_c = col_to_local.reindex(trip["col"].values)
    keep = local_c.notna().values
    trip = trip[keep]
    local_c = local_c.values[keep].astype(np.int32)
    print(f"n_nonzero_entries={len(trip)}")
    mat = sp.coo_matrix((trip["val"].values, (local_c, trip["row"].values - 1)),
                         shape=(n_cells, n_genes)).tocsr()

    anno = pd.read_csv(f"{RAW}/ren/GSE158055_cell_annotation.csv.gz").set_index("cellName")
    cell_anno = anno.loc[cols_df["cellName"].values]

    obs = pd.DataFrame({
        "Sample_ID": cell_anno["sampleID"].values,
        "orig_celltype": cell_anno["celltype"].values,
    }, index=[f"ren_{cn}" for cn in cols_df["cellName"].values])
    obs = obs.join(target[["disease", "y", "Data_NO", "patient_id"]], on="Sample_ID")
    obs["monocyte_source_provenance"] = "curated_atlas_label"
    obs["source_study"] = "Ren_GSE158055"

    adata = ad.AnnData(X=mat, obs=obs, var=pd.DataFrame(index=features))
    adata.var_names_make_unique()
    print(f"Ren adata: {adata.shape}, fatal cells={int((adata.obs.y==1).sum())}, "
          f"donors={adata.obs.patient_id.nunique()}", flush=True)
    return adata


# --------------------------------------------------------------- Qiu sepsis
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


def build_qiu_sepsis(meta):
    print("=== Qiu sepsis GSE167363 ===", flush=True)
    import tarfile
    d = f"{RAW}/qiu_sepsis"
    extracted = f"{d}/extracted"
    os.makedirs(extracted, exist_ok=True)
    if not os.listdir(extracted):
        with tarfile.open(f"{d}/GSE167363_RAW.tar") as t:
            t.extractall(extracted)

    files = os.listdir(extracted)
    adatas = []
    for gsm_name, sample_id in QIU_GSM_TO_SAMPLE.items():
        mtx_f = [x for x in files if gsm_name in x and x.endswith("matrix.mtx.gz")]
        if not mtx_f:
            print(f"  MISSING files for {gsm_name}, skipping")
            continue
        prefix = mtx_f[0].replace("matrix.mtx.gz", "")
        bc = pd.read_csv(f"{extracted}/{prefix}barcodes.tsv.gz", header=None)[0].values
        feat = pd.read_csv(f"{extracted}/{prefix}features.tsv.gz", header=None, sep="\t")
        gene_symbols = feat[1].values if feat.shape[1] > 1 else feat[0].values
        with gzip.open(f"{extracted}/{prefix}matrix.mtx.gz", "rb") as fh:
            mat = sio.mmread(fh).tocsr()
        # standard 10x convention: genes x cells
        if mat.shape[0] == len(gene_symbols):
            mat = mat.T.tocsr()
        assert mat.shape[0] == len(bc), f"{gsm_name}: shape mismatch {mat.shape} vs {len(bc)} barcodes"

        gidx = {g: i for i, g in enumerate(gene_symbols)}
        mono_score = score_lineage(mat, gidx, MONO_MARKERS)
        other_scores = {k: score_lineage(mat, gidx, v) for k, v in OTHER_LINEAGE.items()}
        is_mono = (mono_score > 0)
        for k, s in other_scores.items():
            is_mono &= (mono_score >= s)
        print(f"  {gsm_name} -> {sample_id}: {mat.shape[0]} cells, "
              f"{int(is_mono.sum())} marker-monocyte", flush=True)
        if is_mono.sum() == 0:
            continue
        sub = mat[is_mono, :]
        obs = pd.DataFrame({
            "Sample_ID": sample_id,
        }, index=[f"qiu_{gsm_name}_{i}" for i in np.where(is_mono)[0]])
        one = ad.AnnData(X=sub, obs=obs, var=pd.DataFrame(index=gene_symbols))
        one.var_names_make_unique()
        adatas.append(one)

    if not adatas:
        return None
    adata = ad.concat(adatas, join="outer", fill_value=0)
    adata.obs = adata.obs.join(meta[["disease", "y", "Data_NO", "patient_id"]], on="Sample_ID")
    adata.obs["monocyte_source_provenance"] = "marker_gene_threshold"
    adata.obs["source_study"] = "Qiu_sepsis_GSE167363"
    adata.var_names_make_unique()
    print(f"Qiu sepsis adata: {adata.shape}, fatal cells={int((adata.obs.y==1).sum())}, "
          f"donors={adata.obs.patient_id.nunique()}", flush=True)
    return adata


MONO_LABELS_STEPHENSON = ["C1_CD16_mono", "CD14_mono", "CD16_mono", "CD83_CD14_mono", "Mono_prolif"]


def build_stephenson(meta):
    print("=== Stephenson E-MTAB-10026 ===", flush=True)
    target = meta[meta["Data_NO"] == "302004data17"]
    path = f"{RAW}/stephenson/covid_portal_210320_with_raw.h5ad"
    backed = ad.read_h5ad(path, backed="r")
    obs = backed.obs
    is_mono = obs["full_clustering"].astype(str).isin(MONO_LABELS_STEPHENSON)
    target_sample_ids = set(target.index)
    in_cohort = obs["sample_id"].astype(str).isin(target_sample_ids)
    mask = (is_mono & in_cohort).values
    print(f"Stephenson: {int(mask.sum())} monocyte cells in {int(in_cohort.sum())} "
          f"target-sample cells (of {backed.n_obs} total in atlas)", flush=True)

    sub = backed[mask].to_memory()
    sub.X = sub.layers["raw"]
    del sub.layers["raw"]

    new_obs = pd.DataFrame(index=sub.obs.index)
    new_obs["Sample_ID"] = sub.obs["sample_id"].astype(str).values
    new_obs["orig_celltype"] = sub.obs["full_clustering"].astype(str).values
    # patient_id here is the h5ad's OWN field (Characteristics[individual]-equivalent,
    # verified e.g. sample_id "BGCV02_CV0052" -> patient_id "CV0052"), NOT Job H's
    # buggy patient_id column (which used the sequencing pool "BGCV02" for these rows).
    new_obs["patient_id"] = sub.obs["patient_id"].astype(str).values
    new_obs = new_obs.join(target[["disease", "y", "Data_NO"]], on="Sample_ID")
    new_obs["monocyte_source_provenance"] = "curated_atlas_label"
    new_obs["source_study"] = "Stephenson_EMTAB10026"
    new_obs.index = [f"steph_{i}" for i in sub.obs.index]
    sub.obs = new_obs
    sub.var_names_make_unique()
    print(f"Stephenson adata: {sub.shape}, fatal cells={int((sub.obs.y==1).sum())}, "
          f"donors={sub.obs.patient_id.nunique()}", flush=True)
    return sub


def combine_and_report(ren, qiu, steph):
    parts = [a for a in (ren, qiu, steph) if a is not None]
    combined = ad.concat(parts, join="inner", fill_value=0)
    combined.obs["Disease"] = combined.obs["disease"]
    combined.obs["covid_severity"] = np.where(combined.obs["y"] == 1, "fatal", "survived")
    combined.write(f"{OUT_DIR}/monocyte_canonical_raw.h5ad")

    print("\n=== Job F Phase 2: monocyte canonical object (checkpoint 2) ===", flush=True)
    print(f"Total: {combined.shape[0]} cells x {combined.shape[1]} genes "
          f"(gene intersection across Ren/Stephenson/Qiu-sepsis)", flush=True)
    for study, sub in combined.obs.groupby("source_study"):
        n_fatal_pt = sub[sub.y == 1].patient_id.nunique()
        n_surv_pt = sub[sub.y == 0].patient_id.nunique()
        print(f"  {study}: {len(sub)} cells, {sub.patient_id.nunique()} patients "
              f"({n_fatal_pt} fatal / {n_surv_pt} survived), "
              f"provenance={sub.monocyte_source_provenance.iloc[0]}", flush=True)
    print(f"TOTAL patients: {combined.obs.patient_id.nunique()}, "
          f"fatal={combined.obs[combined.obs.y==1].patient_id.nunique()}, "
          f"survived={combined.obs[combined.obs.y==0].patient_id.nunique()}", flush=True)
    go = (combined.shape[0] >= 20000 and
          (combined.obs.y == 1).any() and (combined.obs.y == 0).any())
    print(f"GO criteria (>=20k cells, both classes represented): "
          f"{'GO' if go else 'NO-GO'}", flush=True)


if __name__ == "__main__":
    os.makedirs(OUT_DIR, exist_ok=True)
    meta = load_meta()
    ren = build_ren(meta)
    ren.write(f"{OUT_DIR}/ren_monocyte_raw.h5ad")
    qiu = build_qiu_sepsis(meta)
    if qiu is not None:
        qiu.write(f"{OUT_DIR}/qiu_sepsis_monocyte_raw.h5ad")
    steph = build_stephenson(meta)
    steph.write(f"{OUT_DIR}/stephenson_monocyte_raw.h5ad")
    combine_and_report(ren, qiu, steph)
    print("DONE Phase 2", flush=True)
