#!/bin/bash
#SBATCH --job-name=jobS_scgpt
#SBATCH --partition=short_gpu
#SBATCH --gres=gpu:1
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu04,gpu13,gpu14
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=2:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/jobS_scgpt_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/jobS_scgpt_%j.err

set -euo pipefail
BASE=/bigdata/godziklab/shared/Xinru/302006
REPO=$BASE/302006_updategit
SIF=$BASE/07_CONTAINERS/scgpt_0.1.7.sif
MODEL=$BASE/models/scGPT_human/scGPT_human
INPUT=$BASE/01_DATA/raw/chongqing_sepsis_platelet_raw_counts_symbol.h5ad
OUTPUT=$BASE/02_EMBEDDINGS/option_c_chongqing/scgpt_chongqing_canonical.h5ad
SCRIPT=$REPO/scripts/foundation_models/scGPT/embed_scgpt_option_c.py
mkdir -p "$BASE/02_EMBEDDINGS/option_c_chongqing"

export SINGULARITYENV_PYTHONNOUSERSITE=1
module load singularity 2>/dev/null || echo "using system singularity"
echo "Input: $INPUT"; echo "Output: $OUTPUT"; date

singularity exec --nv --bind "$BASE" --env PYTHONNOUSERSITE=1 "$SIF" \
    python "$SCRIPT" --input "$INPUT" --model "$MODEL" --output "$OUTPUT" --batch-size 32
date
