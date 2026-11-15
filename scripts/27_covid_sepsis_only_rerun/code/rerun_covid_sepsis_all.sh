#!/bin/bash
#SBATCH --job-name=covid_sepsis_all
#SBATCH --partition=epyc
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/27_covid_sepsis_only_rerun/logs/rerun_all_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/27_covid_sepsis_only_rerun/logs/rerun_all_%j.err
#SBATCH --requeue

# PATH-prepend activation (avoids 'source activate' silent-failure in batch shells)
export PATH=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean/bin:$PATH

REPO=/bigdata/godziklab/shared/Xinru/302006/302006_updategit

cd "$REPO"

mkdir -p "$REPO/scripts/27_covid_sepsis_only_rerun/logs"
mkdir -p /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/option_c_primary_benchmark/covid_sepsis_only

python3 -u "$REPO/scripts/27_covid_sepsis_only_rerun/code/rerun_covid_sepsis_all.py"
