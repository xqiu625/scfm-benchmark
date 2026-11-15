#!/bin/bash
#SBATCH --job-name=jobS_uce
#SBATCH --partition=short_gpu
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu04,gpu13,gpu14
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128gb
#SBATCH --time=2:00:00
#SBATCH --requeue
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/jobS_uce_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/jobS_uce_%j.err

set -euo pipefail
ENV_DIR=/bigdata/godziklab/shared/Xinru/302006/state_env
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="$ENV_DIR"
export PYTHONNOUSERSITE=1

INPUT=/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/chongqing_sepsis_platelet_raw_counts_symbol.h5ad
OUTDIR=/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c_chongqing
OUTPUT="${OUTDIR}/uce_chongqing_canonical.h5ad"
mkdir -p "$OUTDIR"

UCE_DIR=/bigdata/godziklab/shared/Xinru/302006/03_MODELS/uce/UCE
UCE_OUTDIR="${OUTDIR}/uce_run"
mkdir -p "$UCE_OUTDIR"
cd "$UCE_DIR"
echo "Input: $INPUT"; echo "Output: $OUTPUT"; date
python eval_single_anndata.py --adata_path "$INPUT" --dir "$UCE_OUTDIR" --model_loc "$UCE_DIR/model_files/4layer_model.torch" --species human --batch_size 32

PRODUCED=$(find "$UCE_OUTDIR" -name "*_uce_adata.h5ad" -newer "$INPUT" | head -1)
if [ -n "$PRODUCED" ]; then
    mv "$PRODUCED" "$OUTPUT"
    rm -rf "$UCE_OUTDIR"
    echo "Moved $PRODUCED -> $OUTPUT"
else
    echo "ERROR: UCE did not produce *_uce_adata.h5ad"; ls "$UCE_OUTDIR"; exit 1
fi
date
