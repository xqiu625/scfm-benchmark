#!/bin/bash
#SBATCH --job-name=0006_scgpt_optc
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu04,gpu13,gpu14
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/0006_scgpt_optc_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/0006_scgpt_optc_%j.err
#SBATCH --mail-user=xqiu014@ucr.edu
#SBATCH --mail-type=END,FAIL

# scGPT Option-C embedding on the canonical COVID input (via scgpt_0.1.7.sif container).
# Usage:  sbatch submit_scgpt_option_c.sh [N_CELLS]
#   N_CELLS > 0  -> smoke test on that many cells, writes *_smoketest.h5ad
#   (omitted)    -> full run -> 02_EMBEDDINGS/option_c/scgpt_covid_canonical.h5ad
set -euo pipefail

BASE=/bigdata/godziklab/shared/Xinru/302006
REPO=$BASE/302006_updategit
SIF=$BASE/07_CONTAINERS/scgpt_0.1.7.sif
MODEL=$BASE/models/scGPT_human/scGPT_human
INPUT=$BASE/01_DATA/raw/human_platelet_covid_severity_raw_counts.h5ad
SCRIPT=$REPO/scripts/foundation_models/scGPT/embed_scgpt_option_c.py

N_CELLS="${1:-0}"
if [ "$N_CELLS" -gt 0 ]; then
    OUTPUT=$BASE/02_EMBEDDINGS/option_c/scgpt_covid_canonical_smoketest.h5ad
    EXTRA="--n-cells $N_CELLS"
    echo "=== SMOKE TEST: $N_CELLS cells ==="
else
    OUTPUT=$BASE/02_EMBEDDINGS/option_c/scgpt_covid_canonical.h5ad
    EXTRA=""
    echo "=== FULL RUN ==="
fi

echo "Host: $(hostname) | Start: $(date)"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true
module load singularity 2>/dev/null || echo "using system singularity"
mkdir -p "$BASE/02_EMBEDDINGS/option_c" "$BASE/06_LOGS"

# Make the container python ignore the host's ~/.local site-packages — otherwise it
# imports the user-site pyarrow (libarrow.so.2100, needs GLIBC_2.28) and fails against
# the container's older libc. Classic HPCC trap (see docs techcontext / hpcc_snippets).
export SINGULARITYENV_PYTHONNOUSERSITE=1

singularity exec --nv --bind "$BASE" --env PYTHONNOUSERSITE=1 "$SIF" \
    python "$SCRIPT" \
        --input "$INPUT" \
        --model "$MODEL" \
        --output "$OUTPUT" \
        --batch-size 32 $EXTRA

echo "Done: $(date)  -> $OUTPUT"
