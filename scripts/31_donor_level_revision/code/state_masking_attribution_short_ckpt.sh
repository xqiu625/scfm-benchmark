#!/bin/bash
# Checkpointed short_gpu chain: runs state_masking_attribution.py within a
# 2h short_gpu allocation. The python script saves a checkpoint and exits
# cleanly (code 0) if it approaches TIME_BUDGET_SECONDS before finishing.
# This wrapper then checks whether the FINAL output exists; if not, it
# resubmits itself with --dependency=afterany so the next link in the chain
# starts as soon as this allocation ends, without any manual babysitting.
# Chosen over waiting for one uninterrupted multi-hour `gpu`-partition
# allocation because short_gpu has repeatedly shown shorter queue waits in
# this session's own observations (more "mix"/partially-free nodes).
#SBATCH --job-name=state_mask_ckpt
#SBATCH --partition=short_gpu
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu04,gpu13,gpu14
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/state_mask_ckpt_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/state_mask_ckpt_%j.err

set -euo pipefail
ENV_DIR=/bigdata/godziklab/shared/Xinru/302006/state_env
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export PYTHONNOUSERSITE=1
export TIME_BUDGET_SECONDS=5400
cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

SCRIPT_SELF=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/code/state_masking_attribution_short_ckpt.sh
FINAL_OUTPUT=/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/gene_attribution_gsea/state_masking/state_masking_gene_attribution_ranked.csv

python3 -u scripts/31_donor_level_revision/code/state_masking_attribution.py

if [ -f "$FINAL_OUTPUT" ]; then
    echo "COMPLETE -- final output exists at $FINAL_OUTPUT, no further chaining."
else
    echo "Not yet complete -- chaining next short_gpu allocation (dependency: afterany:$SLURM_JOB_ID)"
    sbatch --dependency=afterany:$SLURM_JOB_ID "$SCRIPT_SELF"
fi
