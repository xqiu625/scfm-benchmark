#!/bin/bash
# =============================================================================
# Re-run ACC 02, 03, 04 — Round 3
#
# ACC 02 & 03 & 10: v2 ran all 5 models in one job (12h). STATE binary alone
#   took ~4h, so jobs timed out. Fix: one job per model, 12h each.
# ACC 04: PCA baselines finished but Raw_LogReg on 12k features ran 9h+
#   without finishing. Fix: resubmit with 48h.
# ACC 06: Cross-disease produced 0 results — two bugs:
#   (a) severity column lookup missed 'covid_severity' (used 'Unknown' fallback)
#   (b) disease matching used 'COVID' but data has 'COVID-19'
# ACC 09: Convergence warnings — max_iter increased from 500 to 2000.
#
# Total: 18 jobs (5×ACC02 + 5×ACC03 + 1×ACC04 + 5×ACC10 + 1×ACC06 + 1×ACC09).
#
# Usage:
#   cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
#   bash submit_rerun_v3_split.sh
#
# Author: Xinru Qiu
# Date: April 2026
# =============================================================================

LOG_DIR="/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3"
mkdir -p "$LOG_DIR"

echo "============================================================"
echo "SUBMITTING ACC 02 & ACC 03 — SPLIT BY MODEL (ROUND 3)"
echo "Date: $(date)"
echo "Log dir: $LOG_DIR"
echo "============================================================"

# =================================================================
# ACC 02: Core Severity Benchmark — 1 job per model
# =================================================================

# --- ACC 02: STATE ---
JOB_02_STATE=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=02_STATE_v3
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/02_STATE_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/02_STATE_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 02: STATE — severity benchmark (patient-level CV) v3 ==="
echo "Start: $(date)"
python analyses/02_embedding_classification/code/benchmark_single_model.py \
    --model STATE --tasks binary 3-class 6-class
echo "Completed: $(date)"
SLURM
)
echo "ACC 02 STATE:  $JOB_02_STATE"

# --- ACC 02: UCE ---
JOB_02_UCE=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=02_UCE_v3
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/02_UCE_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/02_UCE_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 02: UCE — severity benchmark (patient-level CV) v3 ==="
echo "Start: $(date)"
python analyses/02_embedding_classification/code/benchmark_single_model.py \
    --model UCE --tasks binary 3-class 6-class
echo "Completed: $(date)"
SLURM
)
echo "ACC 02 UCE:    $JOB_02_UCE"

# --- ACC 02: scGPT ---
JOB_02_scGPT=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=02_scGPT_v3
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/02_scGPT_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/02_scGPT_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 02: scGPT — severity benchmark (patient-level CV) v3 ==="
echo "Start: $(date)"
python analyses/02_embedding_classification/code/benchmark_single_model.py \
    --model scGPT --tasks binary 3-class 6-class
echo "Completed: $(date)"
SLURM
)
echo "ACC 02 scGPT:  $JOB_02_scGPT"

# --- ACC 02: Geneformer ---
JOB_02_GF=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=02_GF_v3
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/02_GF_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/02_GF_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 02: Geneformer — severity benchmark (patient-level CV) v3 ==="
echo "Start: $(date)"
python analyses/02_embedding_classification/code/benchmark_single_model.py \
    --model Geneformer --tasks binary 3-class 6-class
echo "Completed: $(date)"
SLURM
)
echo "ACC 02 GF:     $JOB_02_GF"

# --- ACC 02: TranscriptFormer ---
JOB_02_TF=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=02_TF_v3
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/02_TF_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/02_TF_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 02: TranscriptFormer — severity benchmark (patient-level CV) v3 ==="
echo "Start: $(date)"
python analyses/02_embedding_classification/code/benchmark_single_model.py \
    --model TranscriptFormer --tasks binary 3-class 6-class
echo "Completed: $(date)"
SLURM
)
echo "ACC 02 TF:     $JOB_02_TF"

# =================================================================
# ACC 03: Zero-shot Benchmark — 1 job per model
# =================================================================

# --- ACC 03: STATE ---
JOB_03_STATE=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=03_STATE_v3
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/03_STATE_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/03_STATE_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 03: STATE — zero-shot benchmark (patient-level CV) v3 ==="
echo "Start: $(date)"
python analyses/03_zeroshot_benchmark/code/benchmark_single_model.py \
    --model STATE --tasks binary 3-class 6-class
echo "Completed: $(date)"
SLURM
)
echo "ACC 03 STATE:  $JOB_03_STATE"

# --- ACC 03: UCE ---
JOB_03_UCE=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=03_UCE_v3
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/03_UCE_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/03_UCE_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 03: UCE — zero-shot benchmark (patient-level CV) v3 ==="
echo "Start: $(date)"
python analyses/03_zeroshot_benchmark/code/benchmark_single_model.py \
    --model UCE --tasks binary 3-class 6-class
echo "Completed: $(date)"
SLURM
)
echo "ACC 03 UCE:    $JOB_03_UCE"

# --- ACC 03: scGPT ---
JOB_03_scGPT=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=03_scGPT_v3
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/03_scGPT_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/03_scGPT_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 03: scGPT — zero-shot benchmark (patient-level CV) v3 ==="
echo "Start: $(date)"
python analyses/03_zeroshot_benchmark/code/benchmark_single_model.py \
    --model scGPT --tasks binary 3-class 6-class
echo "Completed: $(date)"
SLURM
)
echo "ACC 03 scGPT:  $JOB_03_scGPT"

# --- ACC 03: Geneformer ---
JOB_03_GF=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=03_GF_v3
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/03_GF_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/03_GF_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 03: Geneformer — zero-shot benchmark (patient-level CV) v3 ==="
echo "Start: $(date)"
python analyses/03_zeroshot_benchmark/code/benchmark_single_model.py \
    --model Geneformer --tasks binary 3-class 6-class
echo "Completed: $(date)"
SLURM
)
echo "ACC 03 GF:     $JOB_03_GF"

# --- ACC 03: TranscriptFormer ---
JOB_03_TF=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=03_TF_v3
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/03_TF_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/03_TF_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 03: TranscriptFormer — zero-shot benchmark (patient-level CV) v3 ==="
echo "Start: $(date)"
python analyses/03_zeroshot_benchmark/code/benchmark_single_model.py \
    --model TranscriptFormer --tasks binary 3-class 6-class
echo "Completed: $(date)"
SLURM
)
echo "ACC 03 TF:     $JOB_03_TF"

# =================================================================
# ACC 04: PCA Baseline — resubmit with 48h
# Raw_LogReg on 12,395 features ran 9h+ without finishing at 12h.
# =================================================================
JOB_04=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=04_pca_v3
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=48:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/04_pca_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/04_pca_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 04: PCA Baseline (patient-level CV) v3 — 48h ==="
echo "Start: $(date)"
python analyses/04_pca_baseline/code/benchmark_pca_baseline.py
echo "Completed: $(date)"
SLURM
)
echo "ACC 04 PCA:    $JOB_04"

# =================================================================
# ACC 10: Clinical Metrics — 1 job per model
# Each model takes ~4h for binary alone; 5 models can't fit in 12h.
# =================================================================

# --- ACC 10: STATE ---
JOB_10_STATE=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=10_STATE_v3
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/10_STATE_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/10_STATE_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 10: STATE — clinical metrics (patient-level CV) v3 ==="
echo "Start: $(date)"
python analyses/10_clinical_metrics/code/clinical_metrics_extended.py --models STATE
echo "Completed: $(date)"
SLURM
)
echo "ACC 10 STATE:  $JOB_10_STATE"

# --- ACC 10: UCE ---
JOB_10_UCE=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=10_UCE_v3
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/10_UCE_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/10_UCE_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 10: UCE — clinical metrics (patient-level CV) v3 ==="
echo "Start: $(date)"
python analyses/10_clinical_metrics/code/clinical_metrics_extended.py --models UCE
echo "Completed: $(date)"
SLURM
)
echo "ACC 10 UCE:    $JOB_10_UCE"

# --- ACC 10: scGPT ---
JOB_10_scGPT=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=10_scGPT_v3
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/10_scGPT_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/10_scGPT_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 10: scGPT — clinical metrics (patient-level CV) v3 ==="
echo "Start: $(date)"
python analyses/10_clinical_metrics/code/clinical_metrics_extended.py --models scGPT
echo "Completed: $(date)"
SLURM
)
echo "ACC 10 scGPT:  $JOB_10_scGPT"

# --- ACC 10: Geneformer ---
JOB_10_GF=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=10_GF_v3
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/10_GF_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/10_GF_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 10: Geneformer — clinical metrics (patient-level CV) v3 ==="
echo "Start: $(date)"
python analyses/10_clinical_metrics/code/clinical_metrics_extended.py --models Geneformer
echo "Completed: $(date)"
SLURM
)
echo "ACC 10 GF:     $JOB_10_GF"

# --- ACC 10: TranscriptFormer ---
JOB_10_TF=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=10_TF_v3
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/10_TF_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/10_TF_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 10: TranscriptFormer — clinical metrics (patient-level CV) v3 ==="
echo "Start: $(date)"
python analyses/10_clinical_metrics/code/clinical_metrics_extended.py --models TranscriptFormer
echo "Completed: $(date)"
SLURM
)
echo "ACC 10 TF:     $JOB_10_TF"

# =================================================================
# ACC 06: Cross-Disease Generalization — bugfix
# Two bugs fixed: severity column lookup + disease name matching
# =================================================================
JOB_06=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=06_crossdis_v3
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=4:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/06_crossdis_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/06_crossdis_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 06: Cross-Disease Generalization (bugfix v3) ==="
echo "Start: $(date)"
python analyses/06_batch_generalization/code/cross_disease_generalization.py
echo "Completed: $(date)"
SLURM
)
echo "ACC 06 cross-disease: $JOB_06"

# =================================================================
# ACC 09: Embedding Ablation — max_iter 500→2000
# =================================================================
JOB_09=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=09_ablation_v3
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/09_ablation_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun_v3/09_ablation_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 09: Embedding Ablation (max_iter fix) ==="
echo "Start: $(date)"
python analyses/09_embedding_ablation/code/embedding_ablations_fast.py
echo "Completed: $(date)"
SLURM
)
echo "ACC 09 ablation: $JOB_09"

# -----------------------------------------------------------------
# Summary
# -----------------------------------------------------------------
echo ""
echo "============================================================"
echo "ROUND 3 — 18 JOBS"
echo "============================================================"
echo "ACC 02 (severity, 12h each):"
echo "  STATE:            $JOB_02_STATE"
echo "  UCE:              $JOB_02_UCE"
echo "  scGPT:            $JOB_02_scGPT"
echo "  Geneformer:       $JOB_02_GF"
echo "  TranscriptFormer: $JOB_02_TF"
echo ""
echo "ACC 03 (zero-shot, 12h each):"
echo "  STATE:            $JOB_03_STATE"
echo "  UCE:              $JOB_03_UCE"
echo "  scGPT:            $JOB_03_scGPT"
echo "  Geneformer:       $JOB_03_GF"
echo "  TranscriptFormer: $JOB_03_TF"
echo ""
echo "ACC 04 (PCA baseline, 48h):"
echo "  PCA+Raw:          $JOB_04"
echo ""
echo "ACC 10 (clinical metrics, 12h each):"
echo "  STATE:            $JOB_10_STATE"
echo "  UCE:              $JOB_10_UCE"
echo "  scGPT:            $JOB_10_scGPT"
echo "  Geneformer:       $JOB_10_GF"
echo "  TranscriptFormer: $JOB_10_TF"
echo ""
echo "ACC 06 (cross-disease, bugfix, 4h):"
echo "  cross-disease:    $JOB_06"
echo ""
echo "ACC 09 (ablation, max_iter fix, 12h):"
echo "  ablation:         $JOB_09"
echo "============================================================"
echo ""
echo "Monitor: squeue -u xinruqiu"
echo "Logs:    ls -lt $LOG_DIR/"
echo "============================================================"
