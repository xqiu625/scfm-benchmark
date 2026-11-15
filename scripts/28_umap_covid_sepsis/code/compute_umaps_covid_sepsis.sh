#!/bin/bash
#SBATCH --job-name=umap_covid_sepsis
#SBATCH --partition=epyc
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/28_umap_covid_sepsis/logs/umap_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/28_umap_covid_sepsis/logs/umap_%j.err
#SBATCH --requeue

# PATH-prepend activation (avoids 'source activate' silent-failure in batch shells)
export PATH=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean/bin:$PATH

REPO=/bigdata/godziklab/shared/Xinru/302006/302006_updategit

cd "$REPO"

mkdir -p "$REPO/scripts/28_umap_covid_sepsis/logs"
mkdir -p /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/option_c_primary_benchmark/covid_sepsis_only/umap

python3 -u "$REPO/scripts/28_umap_covid_sepsis/code/compute_umaps_covid_sepsis.py"
