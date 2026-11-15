#!/bin/bash
#SBATCH --job-name=cross_disease_fatal
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=168:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/cross_disease_fatal_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/cross_disease_fatal_%j.err
#SBATCH --mail-user=xqiu014@ucr.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

echo "=== Cross-Disease Generalization: fatal vs survived ==="
echo "Start: $(date)"
echo "Host:  $(hostname)"

# Bypass conda activate (which fails silently in SLURM batch shells without
# proper conda init). PATH-prepend the env's bin directly.
ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
if [ ! -x "$ENV_DIR/bin/python" ]; then
    echo "ERROR: $ENV_DIR/bin/python not found or not executable"
    exit 1
fi
export PATH="$ENV_DIR/bin:$PATH"
export CONDA_PREFIX="$ENV_DIR"
export CONDA_DEFAULT_ENV=302006_env_clean
export PYTHONNOUSERSITE=1

# Verify
echo "Python: $(which python)"
python --version
python -c "import sklearn, scanpy, anndata; print(f'sklearn={sklearn.__version__} scanpy={scanpy.__version__} anndata={anndata.__version__}')"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

# Defaults verified against actual obs columns (inspect_sepsis_outcome_columns.py
# run on 2026-04-30):
#   COVID:  covid_severity has {severe, mild, fatal, moderate, recovered, healthy}
#   Sepsis: covid_severity has {fatal: 2613, severe: 809}  ← 'severe' = ICU survived
#
# UNIFIED IN-HOSPITAL MORTALITY (2026-05-01 decision):
#   - fatal = 1 (both diseases)
#   - survived = 0 includes COVID-severe (ICU survivors) and sepsis-severe (also ICU survivors)
#   - healthy and unknown EXCLUDED (not part of mortality prediction)
#
# PROVENANCE PIN (do not confuse the two invocations):
#   * Manuscript tab:crossdisease (claim C002) is backed by the 2026-05-04 run
#     invoked with:  --covid-survived-tokens mild,moderate,recovered,healthy
#     (COVID-severe EXCLUDED from the survived set; evidence CSV
#     cross_disease_fatal_vs_survived_20260504_152836.csv, n_train=21286;
#     scGPT's row comes from the 2026-04-30 run).
#   * The invocation BELOW (current default: --covid-survived-tokens
#     mild,moderate,recovered,severe — severe INCLUDED as ICU-survivor) backs
#     the post-2026-05-01 unified-label rerun outputs under
#     05_RESULTS/cross_disease_fatal_vs_survived/ (e.g. the 2026-05-22
#     severe-included CSV, n_train=38673) — NOT the manuscript table.
#
# scGPT excluded from default model list per project decision (anti-correlated
# cross-disease, manuscript-grade negative result documented separately).
SEPSIS_LABEL_COL="${SEPSIS_LABEL_COL:-covid_severity}"
SEPSIS_FATAL_TOKENS="${SEPSIS_FATAL_TOKENS:-fatal}"
SEPSIS_SURVIVED_TOKENS="${SEPSIS_SURVIVED_TOKENS:-severe}"
MODELS="${MODELS:-STATE UCE TranscriptFormer scPRINT2 STACK Geneformer}"  # Option C 6-FM (Geneformer added 2026-05-04 via native scpeft_env)

echo "Sepsis label col: $SEPSIS_LABEL_COL"
echo "Sepsis fatal tokens: $SEPSIS_FATAL_TOKENS"
echo "Models: $MODELS"

python scripts/06_batch_generalization/code/cross_disease_fatal_vs_survived.py \
    --models $MODELS \
    --covid-label-col covid_severity \
    --covid-fatal-tokens fatal \
    --covid-survived-tokens mild,moderate,recovered,severe \
    --sepsis-label-col "$SEPSIS_LABEL_COL" \
    --sepsis-fatal-tokens "$SEPSIS_FATAL_TOKENS" \
    --sepsis-survived-tokens "$SEPSIS_SURVIVED_TOKENS" \
    --group-col Sample_ID \
    --n-bootstrap 200

echo "Exit code: $?"
echo "Completed: $(date)"
