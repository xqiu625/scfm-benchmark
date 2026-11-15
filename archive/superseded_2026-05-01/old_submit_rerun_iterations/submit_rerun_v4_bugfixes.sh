#!/bin/bash
# =============================================================================
# Re-run analyses affected by critical bug fixes (v4)
#
# Fixes applied in this round:
#   #1.2  Bootstrap OOB evaluation      → ACC 08
#   #1.3  DeLong formula corrected       → ACC 08
#   #1.5  Silent raw-expression fallback → ACC 06
#   #1.6  PCA inside CV loop            → ACC 08, 09, 15
#   #1.7  Cell alignment assertion       → ACC 18
#   #1.8  Attention fusion groups bug    → ACC 18
#
# Usage:
#   1. Push updated code to HPCC:
#      rsync -avz --exclude='.git' --exclude='__pycache__' \
#        302006_updategit/ \
#        cluster.hpcc.ucr.edu:/bigdata/godziklab/shared/Xinru/302006/302006_updategit/
#
#   2. SSH into HPCC and submit:
#      cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
#      bash submit_rerun_v4_bugfixes.sh
#
# Author: Xinru Qiu
# Date: April 2026
# =============================================================================

BASE_DIR="/bigdata/godziklab/shared/Xinru/302006"
CODE_DIR="${BASE_DIR}/302006_updategit"
LOG_DIR="${BASE_DIR}/06_LOGS/patientCV_rerun_v4"
ENV_CLEAN="/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean"
ENV_STATE="${BASE_DIR}/state_env"

mkdir -p "$LOG_DIR"

echo "============================================================"
echo "ROUND 4 — BUG FIX RE-RUNS"
echo "Date: $(date)"
echo "Log dir: $LOG_DIR"
echo "============================================================"
echo ""
echo "Fixes in this round:"
echo "  #1.2  Bootstrap → OOB evaluation (ACC 08)"
echo "  #1.3  DeLong → correct Hanley-McNeil (ACC 08)"
echo "  #1.5  Silent fallback → raise error (ACC 06)"
echo "  #1.6  PCA → inside CV loop (ACC 08, 09, 15)"
echo "  #1.7  Cell alignment → severity assertion (ACC 18)"
echo "  #1.8  Attention fusion → groups passed to cv.split (ACC 18)"
echo ""

# -----------------------------------------------------------------
# ACC 06: Cross-Disease Generalization (fix #1.5)
# Now raises error if no embedding key found (was silently using adata.X)
# -----------------------------------------------------------------
JOB_06=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=06_crossdisease_v4
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=4:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/06_crossdisease_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/06_crossdisease_%j.err

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 06: Cross-Disease Generalization (fix #1.5: no silent fallback) ==="
echo "Start: $(date)"
python analyses/06_batch_generalization/code/cross_disease_generalization.py
echo "Exit code: $?"
echo "Completed: $(date)"
SLURM
)
echo "ACC 06 submitted: $JOB_06"

# -----------------------------------------------------------------
# ACC 08: Statistical Significance (fixes #1.2, #1.3, #1.6)
# Bootstrap now uses OOB, DeLong formula corrected, PCA inside CV
# -----------------------------------------------------------------
JOB_08=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=08_stats_v4
#SBATCH --partition=batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/08_stats_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/08_stats_%j.err

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 08: Statistical Significance (fixes #1.2 OOB, #1.3 DeLong, #1.6 PCA-in-CV) ==="
echo "Start: $(date)"
python analyses/08_statistical_significance/code/statistical_significance_fast.py \
    --models STATE UCE TranscriptFormer scGPT Geneformer \
    --tasks binary 3-class \
    --n-bootstrap 100 --pca-dims 256
echo "Exit code: $?"
echo "Completed: $(date)"
SLURM
)
echo "ACC 08 submitted: $JOB_08"

# -----------------------------------------------------------------
# ACC 09: Embedding Ablation (fix #1.6)
# PCA now fit inside each CV fold, not on full data
# -----------------------------------------------------------------
JOB_09=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=09_ablation_v4
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=8:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/09_ablation_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/09_ablation_%j.err

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 09: Embedding Ablation (fix #1.6: PCA inside CV) ==="
echo "Start: $(date)"
python analyses/09_embedding_ablation/code/embedding_ablations_fast.py \
    --tasks binary 3-class
echo "Exit code: $?"
echo "Completed: $(date)"
SLURM
)
echo "ACC 09 submitted: $JOB_09"

# -----------------------------------------------------------------
# ACC 15: Differential Diagnosis (fix #1.6)
# PCA now fit inside each CV fold for scGPT denoising
# -----------------------------------------------------------------
JOB_15=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=15_diffdiag_v4
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=128gb
#SBATCH --time=8:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/15_diffdiag_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/15_diffdiag_%j.err

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 15: Differential Diagnosis (fix #1.6: PCA inside CV) ==="
echo "Start: $(date)"
python analyses/15_differential_diagnosis/code/run_15_differential_diagnosis.py
echo "Exit code: $?"
echo "Completed: $(date)"
SLURM
)
echo "ACC 15 submitted: $JOB_15"

# -----------------------------------------------------------------
# ACC 18: Multimodal Fusion (fixes #1.7, #1.8) — GPU
# Groups now passed to attention_fusion cv.split + alignment assertion
# -----------------------------------------------------------------
JOB_18=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=18_fusion_v4
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/18_fusion_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/18_fusion_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xqiu014@ucr.edu

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

export ENV_PATH="/bigdata/godziklab/shared/Xinru/302006/state_env"
export PATH="${ENV_PATH}/bin:$PATH"
export CONDA_PREFIX="${ENV_PATH}"
export PYTHONNOUSERSITE=1

echo "=== ACC 18: Multimodal Fusion (fixes #1.7 alignment, #1.8 groups) ==="
echo "Start: $(date)"
python -c "import torch; print(f'PyTorch {torch.__version__}, CUDA: {torch.cuda.is_available()}')"

python analyses/18_multimodal_fusion/code/multimodal_fusion.py \
    --strategy all --attention-epochs 100 --common-dim 256
echo "Exit code: $?"
echo "Completed: $(date)"
SLURM
)
echo "ACC 18 submitted: $JOB_18"

# -----------------------------------------------------------------
# Summary
# -----------------------------------------------------------------
echo ""
echo "============================================================"
echo "ROUND 4 — ALL BUG FIX JOBS SUBMITTED"
echo "============================================================"
echo "  ACC 06  (cross-disease, fix #1.5):        $JOB_06"
echo "  ACC 08  (stats, fixes #1.2/#1.3/#1.6):    $JOB_08"
echo "  ACC 09  (ablation, fix #1.6):             $JOB_09"
echo "  ACC 15  (diff diag, fix #1.6):            $JOB_15"
echo "  ACC 18  (fusion GPU, fixes #1.7/#1.8):    $JOB_18"
echo "============================================================"
echo ""
echo "Monitor:  squeue -u xinruqiu"
echo "Logs:     ls -lt $LOG_DIR/"
echo "============================================================"
