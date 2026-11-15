#!/bin/bash
#SBATCH --job-name=scgpt_rerun
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/scgpt_rerun_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/scgpt_rerun_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xqiu014@ucr.edu

# ============================================================================
# Re-run scGPT benchmarks with cleaned embeddings (Orion data removed)
# ============================================================================

echo "=================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start Time: $(date)"
echo "=================================================="

# Setup environment
cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

export ENV_PATH="/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean"
export PATH="${ENV_PATH}/bin:$PATH"
export CONDA_PREFIX="${ENV_PATH}"
export PYTHONNOUSERSITE=1

echo "Python: $(which python)"

# Run finetuned benchmark (severity prediction)
echo ""
echo "=========================================="
echo "1. Finetuned Benchmark (Severity)"
echo "=========================================="
python scripts/rerun_scgpt_only/run_scgpt_finetuned.py

# Run differential diagnosis
echo ""
echo "=========================================="
echo "2. Differential Diagnosis (5-class)"
echo "=========================================="
python scripts/rerun_scgpt_only/run_scgpt_differential_diagnosis.py

echo ""
echo "=================================================="
echo "End Time: $(date)"
echo "=================================================="

# Copy results to analyses folders
cp -v /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/scgpt_rerun/*.csv \
      /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/02_embedding_classification/results/ 2>/dev/null || true

cp -v /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/scgpt_rerun/*.png \
      /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/15_differential_diagnosis/figures/ 2>/dev/null || true

echo "Job complete!"
