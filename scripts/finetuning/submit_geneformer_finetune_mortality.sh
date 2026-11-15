#!/bin/bash
#SBATCH --job-name=0006_gf_ft_mortality
#SBATCH --partition=short_gpu
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu04,gpu13,gpu14
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=2:00:00
#SBATCH --requeue
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/gf_ft_surv_fatal_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/gf_ft_surv_fatal_%j.err
#SBATCH --mail-user=xqiu014@ucr.edu
#SBATCH --mail-type=END,FAIL,REQUEUE

set -euo pipefail

echo "=== Geneformer LoRA fine-tune: in-hospital mortality (Option C) ==="
echo "Start: $(date)"
echo "Host: $(hostname)"
nvidia-smi || true

ENV_DIR=/bigdata/godziklab/shared/Xinru/302006/scpeft_env
if [ ! -x "$ENV_DIR/bin/python" ]; then
    echo "ERROR: $ENV_DIR/bin/python not found. Set up scpeft_env first."
    exit 1
fi
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="$ENV_DIR"
export PYTHONNOUSERSITE=1
export PYTHONPATH=/bigdata/godziklab/shared/Xinru/302006/03_MODELS/scPEFT:${PYTHONPATH:-}
echo "Python: $(which python)  ($(python --version 2>&1))"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

# Idempotent skip-if-done
RESULT_DIR=/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/geneformer_peft
if ls "$RESULT_DIR"/*mortality*final*.json >/dev/null 2>&1; then
    echo "Already complete: $(ls $RESULT_DIR/*mortality*final*.json | head -1)"
    exit 0
fi

# Note: the existing geneformer_peft_finetune.py looks for a Geneformer-ready
# h5ad at $DATA_DIR/geneformer_ready/covid_geneformer_ready.h5ad. If you've
# tokenized a new Option-C version under that path, this will pick it up
# automatically. Otherwise it falls back to the raw file with ensembl_id.
#
# --resume is always passed so that a fresh start (no prior ckpts) is fine
# (the resume logic checks for files; absence triggers a normal start). On
# subsequent runs (after short_gpu 2h preemption), resume picks up from
# ckpt_latest.pt of whichever fold was in flight.
python scripts/finetuning/geneformer_peft_finetune.py \
    --task mortality \
    --epochs 25 \
    --lr 0.0008 \
    --resume

echo "Completed: $(date)"
ls -lh "$RESULT_DIR/" | tail -10
