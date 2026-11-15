#!/bin/bash
#SBATCH --job-name=20_drug_val
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/20_drug_validation_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/20_drug_validation_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xqiu014@ucr.edu

# ============================================================================
# ACC 20: Drug Target Validation
# ============================================================================
# NOTE: API calls need internet. Run initial queries interactively first
# to cache results, then batch analysis can run offline.
# ============================================================================

echo "=================================================="
echo "ACC 20: Drug Target Validation"
echo "Job ID: $SLURM_JOB_ID"
echo "Start Time: $(date)"
echo "=================================================="

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

export ENV_PATH="/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean"
export PATH="${ENV_PATH}/bin:$PATH"
export CONDA_PREFIX="${ENV_PATH}"
export PYTHONNOUSERSITE=1

mkdir -p /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/20_drug_target_validation/{figures,cache}

echo "Python: $(which python)"

# Run validation (use --skip-api if running offline)
python analyses/20_drug_target_validation/code/validate_therapeutic_targets.py

EXIT_CODE=$?

echo ""
echo "=================================================="
echo "Exit Code: $EXIT_CODE"
echo "End Time: $(date)"
echo "=================================================="

exit $EXIT_CODE
