#!/bin/bash
#SBATCH --partition=batch
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --job-name=99_audit
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/logs/99_audit_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/logs/99_audit_%j.err
#SBATCH --mail-type=END,FAIL

set -euo pipefail

PROJECT_ROOT=/bigdata/godziklab/shared/Xinru/302006
LOGS_DIR="${PROJECT_ROOT}/logs"

cd "${PROJECT_ROOT}"
mkdir -p "${LOGS_DIR}" scripts/99_audit/results

echo "Python: $(which python) ($(python --version 2>&1))"
echo "Start : $(date)"

python scripts/99_audit/code/run_audit.py --repo-root . \
    2>&1 | tee scripts/99_audit/results/run_log.txt

echo "End   : $(date)"
