
cd /bigdata/godziklab/shared/Xinru/302006/scfm/UCE

python eval_single_anndata.py \
    --adata_path "/bigdata/godziklab/shared/Xinru/302006/data/raw/HCT116_platelets_raw_counts.h5ad" \
    --dir "/bigdata/godziklab/shared/Xinru/302006/embeddings/uce/" \
    --species "human" \
    --batch_size 200

python eval_single_anndata.py \
    --adata_path "/bigdata/godziklab/shared/Xinru/302006/data/raw/HCT116_retn_raw_counts.h5ad" \
    --dir "/bigdata/godziklab/shared/Xinru/302006/embeddings/uce/" \
    --species "human" \
    --batch_size 200

python eval_single_anndata.py \
    --adata_path "/bigdata/godziklab/shared/Xinru/302006/data/raw/HEK293T_platelets_raw_counts.h5ad" \
    --dir "/bigdata/godziklab/shared/Xinru/302006/embeddings/uce/" \
    --species "human" \
    --batch_size 200

python eval_single_anndata.py \
    --adata_path "/bigdata/godziklab/shared/Xinru/302006/data/raw/HEK293T_retn_raw_counts.h5ad" \
    --dir "/bigdata/godziklab/shared/Xinru/302006/embeddings/uce/" \
    --species "human" \
    --batch_size 200
