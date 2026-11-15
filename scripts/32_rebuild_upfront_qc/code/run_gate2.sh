#!/bin/bash
# Submits the unified Gate 2 runner. It refuses to start without an authorized
# --runtime-config; no such config exists yet, so this will exit 2 by design.
#SBATCH --job-name=gate2
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/gate2_%x_%j.log
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --time=02:00:00
ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PATH="$ENV_DIR/bin:$PATH"
export PYTHONNOUSERSITE=1
python3 /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/32_rebuild_upfront_qc/code/gate2_run.py "$@"
