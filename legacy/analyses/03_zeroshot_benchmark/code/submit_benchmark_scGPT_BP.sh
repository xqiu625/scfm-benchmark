#!/bin/bash
#SBATCH --job-name=bench_scGPT_BP
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/bench_scGPT_BP_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/bench_scGPT_BP_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

# Zero-shot Benchmark: scGPT Best Practices
# Tasks: binary and 6-class classification

echo "=========================================="
echo "ZERO-SHOT BENCHMARK: scGPT_BP"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start time: $(date)"
echo ""

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate /bigdata/godziklab/shared/Xinru/302006/state_env

# Navigate to code directory
cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit/analyses/03_zeroshot_benchmark/code

# Run benchmark for scGPT_BP - 3-class only (binary & 6-class done)
# Using quick mode: LogReg + RF only
python benchmark_single_model.py --model scGPT_BP --tasks 3-class --quick

echo ""
echo "=========================================="
echo "COMPLETE"
echo "End time: $(date)"
echo "=========================================="
