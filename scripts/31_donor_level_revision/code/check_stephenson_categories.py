import anndata as ad
backed = ad.read_h5ad(
    "/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/monocyte_source/stephenson/covid_portal_210320_with_raw.h5ad",
    backed="r")
print(backed.obs["full_clustering"].value_counts())
print(backed.obs.columns.tolist())
