#!/bin/bash
# =============================================================================
# Re-run FAILED jobs from patient-level CV batch (v2)
#
# Round 1 results (April 21-22):
#   COMPLETED: ACC 08, 09, 11b, 15, 17, 19 (+ 02c crashed with bug, now fixed)
#   TIMED OUT: ACC 02, 02b, 03, 04, 05, 07, 10, 16
#   MISSING:   ACC 18, 22 (GPU — may still be running from round 1)
#
# Fixes in this round:
#   - Time limits increased (12-24h for heavy jobs)
#   - Bug fix: ensemble_models.py — embedding key lookup now covers all FMs
#   - Bug fix: interpretability_analysis.py — groups index mismatch resolved
#
# Usage:
#   cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
#   bash submit_rerun_v2_failed.sh
#
# Author: Xinru Qiu
# Date: April 2026
# =============================================================================

BASE_DIR="/bigdata/godziklab/shared/Xinru/302006"
CODE_DIR="${BASE_DIR}/302006_updategit"
LOG_DIR="${BASE_DIR}/06_LOGS/patientCV_rerun_v2"
ENV_CLEAN="/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean"
ENV_STATE="${BASE_DIR}/state_env"

mkdir -p "$LOG_DIR"

echo "============================================================"
echo "SUBMITTING FAILED JOBS — ROUND 2 (longer time limits + bug fixes)"
echo "Date: $(date)"
echo "Log dir: $LOG_DIR"
echo "============================================================"

# -----------------------------------------------------------------
# ACC 02: Core Severity Benchmark — was 2h, now 12h
# -----------------------------------------------------------------
JOB_02=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=02_severity_v2
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v2/02_severity_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v2/02_severity_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 02: Core Severity Benchmark (patient-level CV) v2 ==="
echo "Start: $(date)"

for MODEL in STATE UCE scGPT Geneformer TranscriptFormer; do
    echo "--- $MODEL ---"
    python analyses/02_embedding_classification/code/benchmark_single_model.py \
        --model "$MODEL" --tasks binary 3-class 6-class
done

echo "Completed: $(date)"
SLURM
)
echo "ACC 02 submitted: $JOB_02"

# -----------------------------------------------------------------
# ACC 02b: Full benchmark — was 6h, now 24h
# -----------------------------------------------------------------
JOB_02b=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=02b_full_v2
#SBATCH --partition=batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --time=24:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v2/02b_full_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v2/02b_full_%j.err

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 02b: Full 5-classifier benchmark (patient-level CV) v2 ==="
echo "Start: $(date)"
python analyses/02_embedding_classification/code/covid_severity_benchmark.py --quick
echo "Completed: $(date)"
SLURM
)
echo "ACC 02b submitted: $JOB_02b"

# -----------------------------------------------------------------
# ACC 02c: Ensemble Models — BUG FIXED (embedding key lookup)
# -----------------------------------------------------------------
JOB_02c=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=02c_ensemble_v2
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v2/02c_ensemble_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v2/02c_ensemble_%j.err

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 02c: Ensemble Models (patient-level CV) v2 — bug fixed ==="
python analyses/02_embedding_classification/code/ensemble_models.py
echo "Completed: $(date)"
SLURM
)
echo "ACC 02c submitted: $JOB_02c"

# -----------------------------------------------------------------
# ACC 03: Zero-shot Benchmark — was 2h, now 12h
# -----------------------------------------------------------------
JOB_03=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=03_zeroshot_v2
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v2/03_zeroshot_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v2/03_zeroshot_%j.err

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 03: Zero-shot Benchmark (patient-level CV) v2 ==="
for MODEL in STATE UCE scGPT Geneformer TranscriptFormer; do
    echo "--- $MODEL ---"
    python analyses/03_zeroshot_benchmark/code/benchmark_single_model.py \
        --model "$MODEL" --tasks binary 3-class 6-class
done
echo "Completed: $(date)"
SLURM
)
echo "ACC 03 submitted: $JOB_03"

# -----------------------------------------------------------------
# ACC 04: PCA Baseline — was 2h, now 12h
# -----------------------------------------------------------------
JOB_04=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=04_pca_v2
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v2/04_pca_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v2/04_pca_%j.err

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 04: PCA Baseline (patient-level CV) v2 ==="
python analyses/04_pca_baseline/code/benchmark_pca_baseline.py
echo "Completed: $(date)"
SLURM
)
echo "ACC 04 submitted: $JOB_04"

# -----------------------------------------------------------------
# ACC 05: Sample Efficiency — was 4h, now 24h
# -----------------------------------------------------------------
JOB_05=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=05_sample_eff_v2
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v2/05_sample_eff_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v2/05_sample_eff_%j.err

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 05: Sample Efficiency (patient-level CV) v2 ==="
python analyses/05_sample_efficiency/code/subsampling_robustness.py \
    --tasks binary --n-repeats 5
echo "Completed: $(date)"
SLURM
)
echo "ACC 05 submitted: $JOB_05"

# -----------------------------------------------------------------
# ACC 07: Gene Dropout — was 6h, now 24h
# Got through STATE+UCE+TF in ~6h, need ~10h total for all 5
# -----------------------------------------------------------------
JOB_07=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=07_dropout_v2
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v2/07_dropout_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v2/07_dropout_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 07: Gene Dropout Robustness (patient-level CV) v2 ==="
python analyses/07_gene_dropout/code/gene_dropout_robustness.py \
    --tasks binary 3-class --n-repeats 5
echo "Completed: $(date)"
SLURM
)
echo "ACC 07 submitted: $JOB_07"

# -----------------------------------------------------------------
# ACC 10: Clinical Metrics — was 4h, now 12h
# -----------------------------------------------------------------
JOB_10=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=10_clinical_v2
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v2/10_clinical_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v2/10_clinical_%j.err

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 10: Clinical Metrics (patient-level CV) v2 ==="
python analyses/10_clinical_metrics/code/clinical_metrics_extended.py
echo "Completed: $(date)"
SLURM
)
echo "ACC 10 submitted: $JOB_10"

# -----------------------------------------------------------------
# ACC 11a: Interpretability — BUG FIXED (groups index mismatch)
# -----------------------------------------------------------------
JOB_11=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=11_interp_v2
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v2/11_interp_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v2/11_interp_%j.err

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 11a: Interpretability Analysis (patient-level CV) v2 — bug fixed ==="
python analyses/11_interpretability/code/interpretability_analysis.py
echo "Completed: $(date)"
SLURM
)
echo "ACC 11 submitted: $JOB_11"

# -----------------------------------------------------------------
# ACC 16: scGPT Best Practices — was 2h, now 8h
# -----------------------------------------------------------------
JOB_16=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=16_scgpt_bp_v2
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=8:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v2/16_scgpt_bp_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v2/16_scgpt_bp_%j.err

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 16: scGPT Best Practices Benchmark (patient-level CV) v2 ==="
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
echo "ROUND 2 — FAILED JOBS RESUBMITTED"
echo "============================================================"
echo "  ACC 02  (severity, 12h):     $JOB_02"
echo "  ACC 02b (full bench, 24h):   $JOB_02b"
echo "  ACC 02c (ensemble, fixed):   $JOB_02c"
echo "  ACC 03  (zeroshot, 12h):     $JOB_03"
echo "  ACC 04  (PCA baseline, 12h): $JOB_04"
echo "  ACC 05  (sample eff, 24h):   $JOB_05"
echo "  ACC 07  (gene dropout, 24h): $JOB_07"
echo "  ACC 10  (clinical, 12h):     $JOB_10"
echo "  ACC 11  (interp, fixed):     $JOB_11"
echo "  ACC 16  (scGPT BP, 8h):      $JOB_16"
echo "============================================================"
echo ""
echo "NOT resubmitted (already completed in round 1):"
echo "  ACC 08  (stats)       — completed"
echo "  ACC 09  (ablation)    — completed"
echo "  ACC 15  (diff diag)   — completed"
echo "  ACC 17  (adapted)     — completed"
echo "  ACC 19  (quantization)— completed"
echo "  ACC 18  (fusion GPU)  — check: squeue -u xinruqiu"
echo "  ACC 22  (ESM GPU)     — check: squeue -u xinruqiu"
echo ""
echo "Monitor: squeue -u xinruqiu"
echo "Logs:    ls -lt $LOG_DIR/"
echo "============================================================"
