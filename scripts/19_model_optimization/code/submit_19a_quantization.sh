#!/bin/bash
#SBATCH --job-name=19a_quant
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:a100:1
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/19a_quantization_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/19a_quantization_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xqiu014@ucr.edu

echo "=================================================="
echo "ACC 19A: Quantization Benchmark"
echo "Job ID: $SLURM_JOB_ID"
echo "Start Time: $(date)"
echo "=================================================="

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

export ENV_PATH="/bigdata/godziklab/shared/Xinru/302006/state_env"
export PATH="${ENV_PATH}/bin:$PATH"
export CONDA_PREFIX="${ENV_PATH}"
export PYTHONNOUSERSITE=1

mkdir -p /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/19_model_optimization/quantization

python scripts/19_model_optimization/code/quantization_benchmark.py

echo "Exit Code: $?"
echo "End Time: $(date)"
