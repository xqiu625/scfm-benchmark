#!/bin/bash
#SBATCH --job-name=jobR_winner
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/jobR_winner_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/jobR_winner_%j.err

set -euo pipefail
echo "=== Job R: winner-corrected max-statistic test + fragility + consensus ==="
date

ENV_DIR=/bigdata/godziklab/shared/Xinru/302006/state_env
export PATH="$ENV_DIR/bin:$PATH"
export PYTHONNOUSERSITE=1

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
python scripts/31_donor_level_revision/code/jobR_winner_test_ensemble.py

echo "DONE_JOBR"
date
