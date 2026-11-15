#!/bin/bash
#SBATCH --job-name=uce_mouse_epithelial
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --time=24:00:00
#SBATCH --mem=250GB
#SBATCH --gres=gpu:1
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/logs/uce_mouse_epithelial_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/logs/uce_mouse_epithelial_%j.err


echo "=== UCE mouse macrophage ==="
echo "Start: $(date)"

source activate /bigdata/godziklab/shared/Xinru/302006_env

cd /bigdata/godziklab/shared/Xinru/302006/scfm/UCE

python eval_single_anndata.py \
    --adata_path "/bigdata/godziklab/shared/Xinru/302006/data/raw/mouse_epithelial_retnla_ko_wt.h5ad" \
    --dir "/bigdata/godziklab/shared/Xinru/302006/embeddings/uce_embeddings/" \
    --species "mouse" \
    --batch_size 200

echo "✅ UCE mouse epithelial completed"