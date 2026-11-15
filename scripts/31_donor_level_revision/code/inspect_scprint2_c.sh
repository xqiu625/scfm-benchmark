#!/bin/bash
#SBATCH --job-name=inspect_scprint2_c
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/inspect_scprint2_c_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/inspect_scprint2_c_%j.err

set -euo pipefail
ENV_DIR=/bigdata/godziklab/shared/Xinru/302006/scprint2_env
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export PYTHONNOUSERSITE=1
export HF_HOME=/bigdata/godziklab/shared/Xinru/.huggingface
cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

python3 -u scripts/31_donor_level_revision/code/inspect_scprint2_c.py
