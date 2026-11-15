#!/bin/bash
# =============================================================================
# Round 4 — Fix all remaining timeout and bug issues
#
# Changes from v3:
#   ACC 02 & 03: Added --quick (LR+RF only, no SVM/MLP). Time: 24h.
#   ACC 10 STATE/UCE/TF: Increased to 24h. GF and scGPT already completed.
#   ACC 06: Bugfix rerun (severity column + disease name matching).
#
# Already completed (not resubmitted):
#   ACC 04 PCA:      running (job 23605292, 48h)
#   ACC 09 Ablation: completed
#   ACC 10 GF:       completed (AUC-ROC 0.778)
#   ACC 10 scGPT:    completed (AUC-ROC 0.630)
#   ACC 08 Stats:    completed
#   ACC 18 Fusion:   completed
#
# IMPORTANT: git pull before running!
#
# Usage:
#   cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
#   git pull
#   bash submit_rerun_v4_fixes.sh
#
# Author: Xinru Qiu
# Date: April 2026
# =============================================================================

LOG_DIR="/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4"
mkdir -p "$LOG_DIR"

# Verify bugfix is present before submitting
if ! grep -q 'covid_severity' analyses/06_batch_generalization/code/cross_disease_generalization.py; then
    echo "ERROR: ACC 06 bugfix not found. Did you git pull?"
    exit 1
fi

echo "============================================================"
echo "ROUND 4 — ALL FIXES (14 jobs)"
echo "Date: $(date)"
echo "Log dir: $LOG_DIR"
echo "============================================================"

# =================================================================
# ACC 02: Core Severity — --quick (LR+RF only), 24h per model
# =================================================================

JOB_02_STATE=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=02_STATE_v4
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/02_STATE_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/02_STATE_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 02: STATE — severity benchmark (quick, patient-level CV) v4 ==="
echo "Start: $(date)"
python analyses/02_embedding_classification/code/benchmark_single_model.py \
    --model STATE --tasks binary 3-class 6-class --quick
echo "Completed: $(date)"
SLURM
)
echo "ACC 02 STATE:  $JOB_02_STATE"

JOB_02_UCE=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=02_UCE_v4
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/02_UCE_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/02_UCE_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 02: UCE — severity benchmark (quick, patient-level CV) v4 ==="
echo "Start: $(date)"
python analyses/02_embedding_classification/code/benchmark_single_model.py \
    --model UCE --tasks binary 3-class 6-class --quick
echo "Completed: $(date)"
SLURM
)
echo "ACC 02 UCE:    $JOB_02_UCE"

JOB_02_scGPT=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=02_scGPT_v4
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/02_scGPT_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/02_scGPT_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 02: scGPT — severity benchmark (quick, patient-level CV) v4 ==="
echo "Start: $(date)"
python analyses/02_embedding_classification/code/benchmark_single_model.py \
    --model scGPT --tasks binary 3-class 6-class --quick
echo "Completed: $(date)"
SLURM
)
echo "ACC 02 scGPT:  $JOB_02_scGPT"

JOB_02_GF=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=02_GF_v4
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/02_GF_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/02_GF_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 02: Geneformer — severity benchmark (quick, patient-level CV) v4 ==="
echo "Start: $(date)"
python analyses/02_embedding_classification/code/benchmark_single_model.py \
    --model Geneformer --tasks binary 3-class 6-class --quick
echo "Completed: $(date)"
SLURM
)
echo "ACC 02 GF:     $JOB_02_GF"

JOB_02_TF=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=02_TF_v4
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/02_TF_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/02_TF_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 02: TranscriptFormer — severity benchmark (quick, patient-level CV) v4 ==="
echo "Start: $(date)"
python analyses/02_embedding_classification/code/benchmark_single_model.py \
    --model TranscriptFormer --tasks binary 3-class 6-class --quick
echo "Completed: $(date)"
SLURM
)
echo "ACC 02 TF:     $JOB_02_TF"

# =================================================================
# ACC 03: Zero-shot Benchmark — --quick (LR+RF only), 24h per model
# =================================================================

JOB_03_STATE=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=03_STATE_v4
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/03_STATE_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/03_STATE_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 03: STATE — zero-shot benchmark (quick, patient-level CV) v4 ==="
echo "Start: $(date)"
python analyses/03_zeroshot_benchmark/code/benchmark_single_model.py \
    --model STATE --tasks binary 3-class 6-class --quick
echo "Completed: $(date)"
SLURM
)
echo "ACC 03 STATE:  $JOB_03_STATE"

JOB_03_UCE=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=03_UCE_v4
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/03_UCE_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/03_UCE_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 03: UCE — zero-shot benchmark (quick, patient-level CV) v4 ==="
echo "Start: $(date)"
python analyses/03_zeroshot_benchmark/code/benchmark_single_model.py \
    --model UCE --tasks binary 3-class 6-class --quick
echo "Completed: $(date)"
SLURM
)
echo "ACC 03 UCE:    $JOB_03_UCE"

JOB_03_scGPT=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=03_scGPT_v4
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/03_scGPT_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/03_scGPT_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 03: scGPT — zero-shot benchmark (quick, patient-level CV) v4 ==="
echo "Start: $(date)"
python analyses/03_zeroshot_benchmark/code/benchmark_single_model.py \
    --model scGPT --tasks binary 3-class 6-class --quick
echo "Completed: $(date)"
SLURM
)
echo "ACC 03 scGPT:  $JOB_03_scGPT"

JOB_03_GF=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=03_GF_v4
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/03_GF_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/03_GF_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 03: Geneformer — zero-shot benchmark (quick, patient-level CV) v4 ==="
echo "Start: $(date)"
python analyses/03_zeroshot_benchmark/code/benchmark_single_model.py \
    --model Geneformer --tasks binary 3-class 6-class --quick
echo "Completed: $(date)"
SLURM
)
echo "ACC 03 GF:     $JOB_03_GF"

JOB_03_TF=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=03_TF_v4
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/03_TF_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/03_TF_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 03: TranscriptFormer — zero-shot benchmark (quick, patient-level CV) v4 ==="
echo "Start: $(date)"
python analyses/03_zeroshot_benchmark/code/benchmark_single_model.py \
    --model TranscriptFormer --tasks binary 3-class 6-class --quick
echo "Completed: $(date)"
SLURM
)
echo "ACC 03 TF:     $JOB_03_TF"

# =================================================================
# ACC 10: Clinical Metrics — STATE, UCE, TF only (GF/scGPT done), 24h
# =================================================================

JOB_10_STATE=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=10_STATE_v4
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/10_STATE_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/10_STATE_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 10: STATE — clinical metrics (patient-level CV) v4 ==="
echo "Start: $(date)"
python analyses/10_clinical_metrics/code/clinical_metrics_extended.py --models STATE
echo "Completed: $(date)"
SLURM
)
echo "ACC 10 STATE:  $JOB_10_STATE"

JOB_10_UCE=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=10_UCE_v4
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/10_UCE_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/10_UCE_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 10: UCE — clinical metrics (patient-level CV) v4 ==="
echo "Start: $(date)"
python analyses/10_clinical_metrics/code/clinical_metrics_extended.py --models UCE
echo "Completed: $(date)"
SLURM
)
echo "ACC 10 UCE:    $JOB_10_UCE"

JOB_10_TF=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=10_TF_v4
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/10_TF_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/10_TF_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 10: TranscriptFormer — clinical metrics (patient-level CV) v4 ==="
echo "Start: $(date)"
python analyses/10_clinical_metrics/code/clinical_metrics_extended.py --models TranscriptFormer
echo "Completed: $(date)"
SLURM
)
echo "ACC 10 TF:     $JOB_10_TF"

# =================================================================
# ACC 06: Cross-Disease Generalization — bugfix rerun
# =================================================================

JOB_06=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=06_crossdis_v4
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=4:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/06_crossdis_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v4/06_crossdis_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 06: Cross-Disease Generalization (bugfix v4) ==="
echo "Start: $(date)"
python analyses/06_batch_generalization/code/cross_disease_generalization.py
echo "Completed: $(date)"
SLURM
)
echo "ACC 06 cross-disease: $JOB_06"

# -----------------------------------------------------------------
# Summary
# -----------------------------------------------------------------
echo ""
echo "============================================================"
echo "ROUND 4 — 14 JOBS SUBMITTED"
echo "============================================================"
echo "ACC 02 (severity, --quick, 24h each):"
echo "  STATE:            $JOB_02_STATE"
echo "  UCE:              $JOB_02_UCE"
echo "  scGPT:            $JOB_02_scGPT"
echo "  Geneformer:       $JOB_02_GF"
echo "  TranscriptFormer: $JOB_02_TF"
echo ""
echo "ACC 03 (zero-shot, --quick, 24h each):"
echo "  STATE:            $JOB_03_STATE"
echo "  UCE:              $JOB_03_UCE"
echo "  scGPT:            $JOB_03_scGPT"
echo "  Geneformer:       $JOB_03_GF"
echo "  TranscriptFormer: $JOB_03_TF"
echo ""
echo "ACC 10 (clinical, 24h each):"
echo "  STATE:            $JOB_10_STATE"
echo "  UCE:              $JOB_10_UCE"
echo "  TranscriptFormer: $JOB_10_TF"
echo ""
echo "ACC 06 (cross-disease, bugfix, 4h):"
echo "  cross-disease:    $JOB_06"
echo "============================================================"
echo ""
echo "Already completed (not resubmitted):"
echo "  ACC 04 PCA:      running (job 23605292, 48h)"
echo "  ACC 09 Ablation: completed"
echo "  ACC 10 GF:       completed (AUC-ROC 0.778)"
echo "  ACC 10 scGPT:    completed (AUC-ROC 0.630)"
echo "  ACC 08 Stats:    completed"
echo "  ACC 18 Fusion:   completed"
echo ""
echo "Monitor: squeue -u xinruqiu"
echo "Logs:    ls -lt $LOG_DIR/"
echo "============================================================"
