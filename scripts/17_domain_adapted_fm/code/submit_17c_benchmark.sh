#!/bin/bash
#SBATCH --job-name=17c_benchmark
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/17c_benchmark_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/17c_benchmark_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xqiu014@ucr.edu

# ============================================================================
# ACC 17C: Benchmark Domain-Adapted Models
# ============================================================================
# Step 1: Extract embeddings from adapted models
# Step 2: Run unified 7-model benchmark (5 original + 2 adapted)
# ============================================================================

echo "=================================================="
echo "ACC 17C: Benchmark Domain-Adapted Models"
echo "Job ID: $SLURM_JOB_ID"
echo "Start Time: $(date)"
echo "=================================================="

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

export ENV_PATH="/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean"
export PATH="${ENV_PATH}/bin:$PATH"
export CONDA_PREFIX="${ENV_PATH}"
export PYTHONNOUSERSITE=1

mkdir -p /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/17_domain_adapted_fm/benchmark

echo "Python: $(which python)"

# Step 1: Extract embeddings (requires GPU env -- run separately if needed)
# python scripts/17_domain_adapted_fm/code/extract_adapted_embeddings.py --model both

# Step 2: Benchmark all 7 models
echo ""
echo "Running 7-model benchmark..."
python scripts/17_domain_adapted_fm/code/benchmark_adapted_models.py

EXIT_CODE=$?

echo ""
echo "=================================================="
echo "Exit Code: $EXIT_CODE"
echo "End Time: $(date)"
echo "=================================================="

exit $EXIT_CODE
