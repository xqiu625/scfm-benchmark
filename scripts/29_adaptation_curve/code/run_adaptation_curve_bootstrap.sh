#!/bin/bash
#SBATCH --job-name=0006_adapt_curve_boot
#SBATCH --partition=epyc
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/29_adaptation_curve/logs/adapt_boot_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/29_adaptation_curve/logs/adapt_boot_%j.err
#SBATCH --requeue
#SBATCH --mail-user=xqiu014@ucr.edu
#SBATCH --mail-type=END,FAIL

# Adaptation curve, bootstrap-corrected (2026-05-22).
# Adds 200-iter patient-level bootstrap on test predictions to match the
# cross-disease analysis protocol (cross_disease_fatal_vs_survived.py).
# Output -> 05_RESULTS/.../adaptation_curve_bootstrap/ (separate dir; legacy
# adaptation_curve/ output is preserved).

export PATH=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean/bin:$PATH
export PYTHONNOUSERSITE=1

REPO=/bigdata/godziklab/shared/Xinru/302006/302006_updategit
cd "$REPO"

mkdir -p "$REPO/scripts/29_adaptation_curve/logs"
mkdir -p /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/option_c_primary_benchmark/covid_sepsis_only/adaptation_curve_bootstrap

python3 -u "$REPO/scripts/29_adaptation_curve/code/run_adaptation_curve_bootstrap.py"
