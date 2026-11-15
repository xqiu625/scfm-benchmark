#!/bin/bash
#SBATCH --job-name=0006_cs_gather
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/0006_cs_gather_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/0006_cs_gather_%j.err
#SBATCH --mail-user=xqiu014@ucr.edu
#SBATCH --mail-type=END,FAIL

# Gather step for the array-parallel stage-27 re-run. Submit with
#   --dependency=afterok:<array_jobid>
# so it runs only after all per-FM + PCA tasks succeed. Concatenates the per-FM
# partials (EMBEDDINGS order) into mortality_lr_rf.csv / calibration.csv / lobo.csv,
# then runs the cross-FM stages: DeLong (stage2) + per-cohort OOF (stage5).
set -euo pipefail

ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
if [ ! -x "$ENV_DIR/bin/python" ]; then echo "ERROR: env not found: $ENV_DIR"; exit 1; fi
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="$ENV_DIR"
export PYTHONNOUSERSITE=1
echo "Gather | Python: $(which python) | Start: $(date)"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
python scripts/27_covid_sepsis_only_rerun/code/rerun_covid_sepsis_all.py --mode gather

echo "Gather done: $(date)"
