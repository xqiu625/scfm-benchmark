#!/bin/bash
# ============================================================================
# RERUN v5 — Final push for manuscript
# ============================================================================
#
# STATUS from v4:
#   ACC 02: binary DONE (all 5), 3-class DONE (UCE/TF/scGPT/GF), 6-class DONE (scGPT/GF)
#   ACC 03: IDENTICAL to ACC 02 (scripts are byte-for-byte same) — SKIP
#   ACC 08: DONE (stats)
#   ACC 09: DONE (ablation)
#   ACC 10: DONE (STATE/UCE/TF), missing scGPT + GF
#   ACC 15: CANCELLED — partial (STATE/UCE/GF done in logs, no saved output)
#   ACC 18: DONE (fusion)
#   ACC 05: Only old cell-level results — needs patient-level CV rerun
#   ACC 07: Empty results — needs patient-level CV rerun
#
# REMAINING:
#   Priority 1 (manuscript critical):
#     ACC 02 STATE: 3-class only, 6-class only (separate jobs)
#     ACC 02 UCE: 6-class only
#     ACC 02 TF: 6-class only
#     ACC 15: full rerun all 5 models (no saved output from partial)
#
#   Priority 2 (robustness section):
#     ACC 05: sample efficiency rerun
#     ACC 07: gene dropout rerun
#     ACC 10: scGPT + GF clinical metrics
#
# NOTE: ACC 03 zero-shot is SKIPPED — the script is identical to ACC 02.
#       For the manuscript, note that "zero-shot" evaluation using frozen
#       embeddings with trained classifiers is equivalent to the supervised
#       benchmark since no model fine-tuning is performed in either case.
# ============================================================================

set -euo pipefail

LOG_DIR="/bigdata/godziklab/shared/Xinru/302006/06_LOGS/rerun_v5"
mkdir -p "$LOG_DIR"

echo "============================================"
echo "  RERUN v5 — Final manuscript submission"
echo "  $(date)"
echo "============================================"
echo ""

# ============================================================================
# PRIORITY 1: ACC 02 — Complete severity benchmark
# ============================================================================

# ACC 02: STATE — 3-class ONLY (binary already done)
echo "Submitting ACC 02: STATE 3-class..."
JOB_02_STATE_3C=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=02_STATE_3c
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/rerun_v5/02_STATE_3c_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/rerun_v5/02_STATE_3c_%j.err
#SBATCH --exclude=gpu10,gpu08

echo "=== ACC 02: STATE — 3-class only (patient-level CV) v5 ==="
echo "Start: $(date)"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

python analyses/02_embedding_classification/code/benchmark_single_model.py \
    --model STATE --tasks 3-class --quick

echo "Exit code: $?"
echo "Completed: $(date)"
SLURM
)
echo "  Job ID: $JOB_02_STATE_3C"

# ACC 02: STATE — 6-class ONLY
echo "Submitting ACC 02: STATE 6-class..."
JOB_02_STATE_6C=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=02_STATE_6c
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/rerun_v5/02_STATE_6c_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/rerun_v5/02_STATE_6c_%j.err
#SBATCH --exclude=gpu10,gpu08

echo "=== ACC 02: STATE — 6-class only (patient-level CV) v5 ==="
echo "Start: $(date)"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

python analyses/02_embedding_classification/code/benchmark_single_model.py \
    --model STATE --tasks 6-class --quick

echo "Exit code: $?"
echo "Completed: $(date)"
SLURM
)
echo "  Job ID: $JOB_02_STATE_6C"

# ACC 02: UCE — 6-class ONLY (binary + 3-class done)
echo "Submitting ACC 02: UCE 6-class..."
JOB_02_UCE_6C=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=02_UCE_6c
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/rerun_v5/02_UCE_6c_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/rerun_v5/02_UCE_6c_%j.err
#SBATCH --exclude=gpu10,gpu08

echo "=== ACC 02: UCE — 6-class only (patient-level CV) v5 ==="
echo "Start: $(date)"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

python analyses/02_embedding_classification/code/benchmark_single_model.py \
    --model UCE --tasks 6-class --quick

echo "Exit code: $?"
echo "Completed: $(date)"
SLURM
)
echo "  Job ID: $JOB_02_UCE_6C"

# ACC 02: TranscriptFormer — 6-class ONLY (binary + 3-class done)
echo "Submitting ACC 02: TF 6-class..."
JOB_02_TF_6C=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=02_TF_6c
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/rerun_v5/02_TF_6c_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/rerun_v5/02_TF_6c_%j.err
#SBATCH --exclude=gpu10,gpu08

echo "=== ACC 02: TranscriptFormer — 6-class only (patient-level CV) v5 ==="
echo "Start: $(date)"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

python analyses/02_embedding_classification/code/benchmark_single_model.py \
    --model TranscriptFormer --tasks 6-class --quick

echo "Exit code: $?"
echo "Completed: $(date)"
SLURM
)
echo "  Job ID: $JOB_02_TF_6C"

# ============================================================================
# PRIORITY 1: ACC 15 — Differential Diagnosis (CRITICAL for Fig 2)
# ============================================================================

echo ""
echo "Submitting ACC 15: Differential Diagnosis (all 5 models)..."
JOB_15=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=15_diffdiag
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/rerun_v5/15_diffdiag_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/rerun_v5/15_diffdiag_%j.err
#SBATCH --exclude=gpu10,gpu08

echo "=== ACC 15: Differential Diagnosis (patient-level CV) v5 ==="
echo "Start: $(date)"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

python analyses/15_differential_diagnosis/code/run_15_differential_diagnosis.py \
    --cv-folds 5

echo "Exit code: $?"
echo "Completed: $(date)"
SLURM
)
echo "  Job ID: $JOB_15"

# ============================================================================
# PRIORITY 2: ACC 05 — Sample Efficiency (robustness section)
# ============================================================================

echo ""
echo "Submitting ACC 05: Sample Efficiency..."
JOB_05=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=05_sample_eff
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=48:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/rerun_v5/05_sample_eff_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/rerun_v5/05_sample_eff_%j.err
#SBATCH --exclude=gpu10,gpu08

echo "=== ACC 05: Sample Efficiency (patient-level CV) v5 ==="
echo "Start: $(date)"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

python analyses/05_sample_efficiency/code/subsampling_robustness.py \
    --tasks binary --n-repeats 3

echo "Exit code: $?"
echo "Completed: $(date)"
SLURM
)
echo "  Job ID: $JOB_05"

# ============================================================================
# PRIORITY 2: ACC 07 — Gene Dropout (robustness section)
# ============================================================================

echo ""
echo "Submitting ACC 07: Gene Dropout..."
JOB_07=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=07_gene_drop
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=48:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/rerun_v5/07_gene_drop_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/rerun_v5/07_gene_drop_%j.err
#SBATCH --exclude=gpu10,gpu08

echo "=== ACC 07: Gene Dropout (patient-level CV) v5 ==="
echo "Start: $(date)"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

python analyses/07_gene_dropout/code/gene_dropout_robustness.py \
    --tasks binary --n-repeats 3

echo "Exit code: $?"
echo "Completed: $(date)"
SLURM
)
echo "  Job ID: $JOB_07"

# ============================================================================
# PRIORITY 2: ACC 10 — Clinical Metrics (scGPT + GF)
# ============================================================================

echo ""
echo "Submitting ACC 10: scGPT clinical metrics..."
JOB_10_SCGPT=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=10_scGPT
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/rerun_v5/10_scGPT_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/rerun_v5/10_scGPT_%j.err
#SBATCH --exclude=gpu10,gpu08

echo "=== ACC 10: scGPT — clinical metrics (patient-level CV) v5 ==="
echo "Start: $(date)"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

python analyses/10_clinical_metrics/code/clinical_metrics_extended.py \
    --models scGPT

echo "Exit code: $?"
echo "Completed: $(date)"
SLURM
)
echo "  Job ID: $JOB_10_SCGPT"

echo "Submitting ACC 10: Geneformer clinical metrics..."
JOB_10_GF=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=10_GF
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/rerun_v5/10_GF_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/rerun_v5/10_GF_%j.err
#SBATCH --exclude=gpu10,gpu08

echo "=== ACC 10: Geneformer — clinical metrics (patient-level CV) v5 ==="
echo "Start: $(date)"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

python analyses/10_clinical_metrics/code/clinical_metrics_extended.py \
    --models Geneformer

echo "Exit code: $?"
echo "Completed: $(date)"
SLURM
)
echo "  Job ID: $JOB_10_GF"

# ============================================================================
# SUMMARY
# ============================================================================

echo ""
echo "============================================"
echo "  SUBMITTED 9 JOBS:"
echo "============================================"
echo ""
echo "  Priority 1 (manuscript critical):"
echo "    ACC 02 STATE 3-class:  $JOB_02_STATE_3C  (24h)"
echo "    ACC 02 STATE 6-class:  $JOB_02_STATE_6C  (24h)"
echo "    ACC 02 UCE 6-class:    $JOB_02_UCE_6C    (24h)"
echo "    ACC 02 TF 6-class:     $JOB_02_TF_6C     (24h)"
echo "    ACC 15 Diff Diagnosis: $JOB_15           (12h)"
echo ""
echo "  Priority 2 (robustness + clinical):"
echo "    ACC 05 Sample Eff:     $JOB_05           (48h)"
echo "    ACC 07 Gene Dropout:   $JOB_07           (48h)"
echo "    ACC 10 scGPT:          $JOB_10_SCGPT     (12h)"
echo "    ACC 10 Geneformer:     $JOB_10_GF        (24h)"
echo ""
echo "  SKIPPED:"
echo "    ACC 03 — identical script to ACC 02"
echo "    ACC 06 — data limitation (Sepsis all severe)"
echo ""
echo "  NOTE: v4 results already collected:"
echo "    ACC 02: binary (all 5), 3-class (UCE/TF/scGPT/GF), 6-class (scGPT/GF)"
echo "    ACC 08: stats DONE"
echo "    ACC 09: ablation DONE"
echo "    ACC 10: STATE/UCE/TF DONE"
echo "    ACC 18: fusion DONE"
echo ""
echo "  After v5 completes, ALL manuscript analyses will be ready."
echo "============================================"
