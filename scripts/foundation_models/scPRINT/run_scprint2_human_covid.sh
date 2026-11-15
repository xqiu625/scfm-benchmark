#!/bin/bash
#SBATCH --job-name=scprint2_human_covid
#SBATCH --partition=short_gpu
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu04,gpu13,gpu14
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=2:00:00
#SBATCH --requeue
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/scprint2_human_covid_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/scprint2_human_covid_%j.err
#SBATCH --mail-user=xqiu014@ucr.edu
#SBATCH --mail-type=END,FAIL,REQUEUE

set -euo pipefail

echo "=== scPRINT-2 human COVID embeddings (Option C, canonical raw input) ==="
echo "Start: $(date)"
echo "Host: $(hostname)"
echo "Job:  $SLURM_JOB_ID  ${SLURM_RESTART_COUNT:+restart=$SLURM_RESTART_COUNT}"
nvidia-smi || true

ENV_DIR=/bigdata/godziklab/shared/Xinru/302006/scprint2_env
if [ ! -x "$ENV_DIR/bin/python" ]; then
    echo "ERROR: $ENV_DIR/bin/python not found. Run install_scprint2.sh first."
    exit 1
fi
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="$ENV_DIR"
export CONDA_DEFAULT_ENV=scprint2_env
export PYTHONNOUSERSITE=1
export HF_HOME=/bigdata/godziklab/shared/Xinru/.huggingface
echo "Python: $(which python)  ($(python --version 2>&1))"

INPUT="${INPUT:-/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/human_platelet_covid_severity_raw_counts.h5ad}"
OUTDIR=/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c
OUTPUT="${OUTDIR}/scprint2_covid_canonical.h5ad"
CKPT="${HF_HOME}/scprint2_ckpts/small-v2.ckpt"

# Fall back to scanning HF cache if local-dir download wasn't used
if [ ! -f "$CKPT" ]; then
    CKPT=$(find "$HF_HOME/hub" -name "small-v2.ckpt" 2>/dev/null | head -1)
fi
if [ ! -f "$CKPT" ]; then
    echo "ERROR: small-v2.ckpt not found. Run install_scprint2.sh first."
    exit 1
fi

if [ ! -f "$INPUT" ]; then
    echo "ERROR: canonical raw input not found: $INPUT"
    exit 2
fi

mkdir -p "$OUTDIR"

# Idempotent skip-if-done
if python -c "import anndata,sys; a=anndata.read_h5ad('$OUTPUT'); sys.exit(0 if any(k.startswith('X_scprint') for k in a.obsm) else 1)" 2>/dev/null; then
    echo "Output already complete: $OUTPUT. Exit."
    exit 0
fi

echo "Input:  $INPUT"
echo "Ckpt:   $CKPT"
echo "Output: $OUTPUT"

# Use scprint2 CLI — much cleaner than custom Python wrapper
scprint2 embed \
    --adata "$INPUT" \
    --ckpt_path "$CKPT" \
    --species "NCBITaxon:9606" \
    --output_filename "$OUTPUT"

echo "Output: $OUTPUT"
ls -lh "$OUTPUT"
echo "Completed: $(date)"
