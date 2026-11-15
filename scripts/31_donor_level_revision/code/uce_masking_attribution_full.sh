#!/bin/bash
#SBATCH --job-name=uce_mask_full
#SBATCH --partition=short_gpu
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu04,gpu13,gpu14
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=2:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/uce_mask_full_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/uce_mask_full_%j.err

set -euo pipefail
ENV_DIR=/bigdata/godziklab/shared/Xinru/302006/state_env
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export PYTHONNOUSERSITE=1
cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

SCRIPT_SELF=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/code/uce_masking_attribution_full.sh
FINAL_OUTPUT=/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/gene_attribution_gsea/uce_masking/uce_masking_gene_attribution_ranked.csv

python3 -u scripts/31_donor_level_revision/code/uce_masking_attribution.py

if [ -f "$FINAL_OUTPUT" ]; then
    echo "COMPLETE -- final output exists at $FINAL_OUTPUT, no further chaining."
else
    echo "Not yet complete -- chaining next short_gpu allocation (dependency: afterany:$SLURM_JOB_ID)"
    sbatch --dependency=afterany:$SLURM_JOB_ID "$SCRIPT_SELF"
fi
