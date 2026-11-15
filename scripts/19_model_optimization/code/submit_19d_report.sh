#!/bin/bash
#SBATCH --job-name=19d_report
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/19d_report_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/19d_report_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xqiu014@ucr.edu

echo "=================================================="
echo "ACC 19D: Optimization Report"
echo "Job ID: $SLURM_JOB_ID"
echo "Start Time: $(date)"
echo "=================================================="

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

export ENV_PATH="/bigdata/godziklab/shared/Xinru/302006/state_env"
export PATH="${ENV_PATH}/bin:$PATH"
export CONDA_PREFIX="${ENV_PATH}"
export PYTHONNOUSERSITE=1

mkdir -p /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/19_model_optimization/figures

python scripts/19_model_optimization/code/optimization_report.py

echo "Exit Code: $?"
echo "End Time: $(date)"
