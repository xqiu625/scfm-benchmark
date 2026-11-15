#!/bin/bash
#SBATCH --job-name=jobS2_r4primary
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/jobS2_r4primary_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/jobS2_r4primary_%j.err

set -euo pipefail
echo "=== Round-4 primary-test correction: permutation null + recentered bootstrap + LOFO + SE check ==="
date

ENV_DIR=/bigdata/godziklab/shared/Xinru/302006/state_env
export PATH="$ENV_DIR/bin:$PATH"
export PYTHONNOUSERSITE=1

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
python scripts/31_donor_level_revision/code/jobS2_round4_primary.py

echo "DONE_JOBS2_R4PRIMARY"
date
