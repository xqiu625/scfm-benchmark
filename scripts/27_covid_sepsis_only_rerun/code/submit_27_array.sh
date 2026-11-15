#!/bin/bash
#SBATCH --job-name=0006_cs_array
#SBATCH --partition=batch
#SBATCH --array=0-7
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/0006_cs_array_%A_%a.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/0006_cs_array_%A_%a.err
#SBATCH --mail-user=xqiu014@ucr.edu
#SBATCH --mail-type=END,FAIL

# Array-parallel COVID+Sepsis-only re-run (stage 27).
# Mapping: SLURM_ARRAY_TASK_ID 0-6 index directly into the EMBEDDINGS list in
#   rerun_covid_sepsis_all.py (0=STATE, 1=UCE, 2=TranscriptFormer, 3=scPRINT-2,
#   4=STACK, 5=Geneformer, 6=scGPT); task 7 = stage3 PCA baseline.
#   tasks 0-6 : per-FM stage1 (OOF + calibration) + stage4 (LOBO) -> _partials/
#   task 7    : stage3 PCA baseline -> pca_baseline.csv
# Cross-FM steps (DeLong, per-cohort) + partial concat run in the gather job.
set -euo pipefail

ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
if [ ! -x "$ENV_DIR/bin/python" ]; then echo "ERROR: env not found: $ENV_DIR"; exit 1; fi
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="$ENV_DIR"
export PYTHONNOUSERSITE=1
echo "Task $SLURM_ARRAY_TASK_ID | Python: $(which python) | Start: $(date)"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
SCRIPT=scripts/27_covid_sepsis_only_rerun/code/rerun_covid_sepsis_all.py

if [ "$SLURM_ARRAY_TASK_ID" -lt 7 ]; then
    python "$SCRIPT" --mode per_fm --index "$SLURM_ARRAY_TASK_ID"
else
    python "$SCRIPT" --mode pca
fi

echo "Task $SLURM_ARRAY_TASK_ID done: $(date)"
