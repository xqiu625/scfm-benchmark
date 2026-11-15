#!/bin/bash
#SBATCH --job-name=gate1_apply
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/gate1_%x_%j.log
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --time=02:00:00
ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PATH="$ENV_DIR/bin:$PATH"
export PYTHONNOUSERSITE=1
python3 /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/32_rebuild_upfront_qc/code/qc_gate1_apply.py \
  --cohort-dir "$1" --out-dir "$2" --layout "$3"
