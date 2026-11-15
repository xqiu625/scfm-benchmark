#!/bin/bash
#SBATCH --job-name=jobS_scprint2
#SBATCH --partition=short_gpu
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu04,gpu13,gpu14
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=2:00:00
#SBATCH --requeue
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/jobS_scprint2_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/jobS_scprint2_%j.err

set -euo pipefail
ENV_DIR=/bigdata/godziklab/shared/Xinru/302006/scprint2_env
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="$ENV_DIR"
export PYTHONNOUSERSITE=1
export HF_HOME=/bigdata/godziklab/shared/Xinru/.huggingface

INPUT=/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/chongqing_sepsis_platelet_raw_counts.h5ad
OUTDIR=/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c_chongqing
OUTPUT="${OUTDIR}/scprint2_chongqing_canonical.h5ad"
CKPT="${HF_HOME}/scprint2_ckpts/small-v2.ckpt"
[ -f "$CKPT" ] || CKPT=$(find "$HF_HOME/hub" -name "small-v2.ckpt" 2>/dev/null | head -1)
mkdir -p "$OUTDIR"
echo "Input: $INPUT"; echo "Output: $OUTPUT"; date

scprint2 embed --adata "$INPUT" --ckpt_path "$CKPT" --species "NCBITaxon:9606" --output_filename "$OUTPUT"
ls -lh "$OUTPUT"; date
