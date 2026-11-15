#!/bin/bash
#SBATCH --job-name=0006_gf_v1_zeroshot
#SBATCH --partition=short_gpu
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu04,gpu13,gpu14
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=2:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/gf_v1_zeroshot_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/gf_v1_zeroshot_%j.err
#SBATCH --mail-user=xqiu014@ucr.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail
echo "=== Zero-shot Geneformer-V1-10M control (roadmap J) ==="
echo "Start: $(date)"; echo "Host: $(hostname)"; nvidia-smi || true

ENV_DIR=/bigdata/godziklab/shared/Xinru/302006/scpeft_env
if [ ! -x "$ENV_DIR/bin/python" ]; then echo "ERROR: scpeft_env not found"; exit 1; fi
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="$ENV_DIR"
export PYTHONNOUSERSITE=1
export PYTHONPATH=/bigdata/godziklab/shared/Xinru/302006/03_MODELS/scPEFT:${PYTHONPATH:-}
echo "Python: $(which python)  ($(python --version 2>&1))"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
python scripts/30_review_response_stats/code/geneformer_v1_zeroshot.py
echo "Completed: $(date)"
