#!/bin/bash
# =============================================================================
# Re-run 3 BUGFIXED jobs (Round 3)
#
# Bugs fixed:
#   ACC 02c — ensemble_models.py: scGPT barcode '-batch_N' suffix stripped
#             so cell alignment finds 47k+ common cells instead of 0
#   ACC 11  — interpretability_analysis.py: multiclass AUC scoring changed
#             from 'roc_auc' to 'roc_auc_ovr' for 3-class tasks
#   ACC 16  — run_scgpt_bp_benchmark.py: classes with <n_splits samples
#             filtered before CV (Flu=0, Non_covid=0 caused crash)
#
# Usage:
#   cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
#   bash submit_bugfix_rerun.sh
#
# Author: Xinru Qiu
# Date: April 2026
# =============================================================================

LOG_DIR="/bigdata/godziklab/shared/Xinru/302006/06_LOGS/bugfix_rerun"
mkdir -p "$LOG_DIR"

echo "============================================================"
echo "SUBMITTING 3 BUGFIX JOBS — ROUND 3"
echo "Date: $(date)"
echo "Log dir: $LOG_DIR"
echo "============================================================"

# -----------------------------------------------------------------
# ACC 02c: Ensemble Models — barcode alignment fix
# -----------------------------------------------------------------
JOB_02c=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=02c_ensemble_v3
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/bugfix_rerun/02c_ensemble_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/bugfix_rerun/02c_ensemble_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 02c: Ensemble Models (bugfix v3 — barcode alignment) ==="
echo "Start: $(date)"
python analyses/02_embedding_classification/code/ensemble_models.py
echo "Completed: $(date)"
SLURM
)
echo "ACC 02c submitted: $JOB_02c"

# -----------------------------------------------------------------
# ACC 11: Interpretability — multiclass AUC fix
# -----------------------------------------------------------------
JOB_11=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=11_interp_v3
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/bugfix_rerun/11_interp_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/bugfix_rerun/11_interp_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 11: Interpretability Analysis (bugfix v3 — multiclass AUC) ==="
echo "Start: $(date)"
python analyses/11_interpretability/code/interpretability_analysis.py
echo "Completed: $(date)"
SLURM
)
echo "ACC 11 submitted: $JOB_11"

# -----------------------------------------------------------------
# ACC 16: scGPT Best Practices — empty class filter fix
# -----------------------------------------------------------------
JOB_16=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=16_scgpt_bp_v3
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=8:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/bugfix_rerun/16_scgpt_bp_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/bugfix_rerun/16_scgpt_bp_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 16: scGPT Best Practices (bugfix v3 — empty class filter) ==="
echo "Start: $(date)"
python analyses/16_embedding_generation_5models/code/run_scgpt_bp_benchmark.py
echo "Completed: $(date)"
SLURM
)
echo "ACC 16 submitted: $JOB_16"

# -----------------------------------------------------------------
# Summary
# -----------------------------------------------------------------
echo ""
echo "============================================================"
echo "ROUND 3 — BUGFIX JOBS SUBMITTED"
echo "============================================================"
echo "  ACC 02c (ensemble, 12h):   $JOB_02c"
echo "  ACC 11  (interp, 12h):     $JOB_11"
echo "  ACC 16  (scGPT BP, 8h):    $JOB_16"
echo "============================================================"
echo ""
echo "Monitor: squeue -u xinruqiu"
echo "Logs:    ls -lt $LOG_DIR/"
echo "============================================================"
