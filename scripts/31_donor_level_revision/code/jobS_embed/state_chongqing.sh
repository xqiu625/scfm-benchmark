#!/bin/bash
#SBATCH --job-name=jobS_state
#SBATCH --partition=short_gpu
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu04,gpu13,gpu14
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=2:00:00
#SBATCH --requeue
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/jobS_state_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/jobS_state_%j.err

set -euo pipefail
ENV_DIR=/bigdata/godziklab/shared/Xinru/302006/state_env
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="$ENV_DIR"
export PYTHONNOUSERSITE=1
export TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1

INPUT=/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/chongqing_sepsis_platelet_raw_counts.h5ad
OUTDIR=/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c_chongqing
OUTPUT="${OUTDIR}/state_chongqing_canonical.h5ad"
mkdir -p "$OUTDIR"
echo "Input: $INPUT"; echo "Output: $OUTPUT"; date

state emb transform --model-folder /bigdata/godziklab/shared/Xinru/302006/03_MODELS/SE-600M \
  --input "$INPUT" --output "$OUTPUT" --embed-key X_state --batch-size 64
ls -lh "$OUTPUT"; date
