#!/bin/bash
#SBATCH --job-name=scprint_human_covid
#SBATCH --partition=short_gpu
#SBATCH --exclude=gpu10,gpu08
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=2:00:00
#SBATCH --requeue
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/scprint_human_covid_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/scprint_human_covid_%j.err
#SBATCH --mail-user=xqiu014@ucr.edu
#SBATCH --mail-type=END,FAIL,REQUEUE

set -euo pipefail

echo "=== scPRINT human COVID embeddings (short_gpu, chunked checkpoint) ==="
echo "Start: $(date)"
echo "Host: $(hostname)"
echo "Job:  $SLURM_JOB_ID  ${SLURM_RESTART_COUNT:+restart=$SLURM_RESTART_COUNT}"
nvidia-smi || true

ENV_DIR=/bigdata/godziklab/shared/Xinru/302006/scprint_env
if [ ! -x "$ENV_DIR/bin/python" ]; then
    echo "ERROR: $ENV_DIR/bin/python not found"; exit 1
fi
export PATH="$ENV_DIR/bin:$PATH"
# Prepend conda's lib/ so we use its libstdc++ (CXXABI_1.3.15+), not Rocky 8's old one.
# Required for scprint → lamindb → django → _sqlite3 → libicui18n.so.78
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="$ENV_DIR"
export CONDA_DEFAULT_ENV=scprint_env
export PYTHONNOUSERSITE=1
export HF_HOME=/bigdata/godziklab/shared/Xinru/.huggingface
echo "Python: $(which python)  ($(python --version 2>&1))"
echo "libstdc++ in use: $(ldd $ENV_DIR/bin/python | grep libstdc++ || echo 'not linked directly')"

REPO=/bigdata/godziklab/shared/Xinru/302006/302006_updategit
INPUT=/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/state/covid_state.h5ad
OUTDIR=/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/scprint
OUTPUT="${OUTDIR}/covid_scprint.h5ad"
CHECKPOINT_DIR="${OUTDIR}/checkpoints/covid_scprint"
mkdir -p "$OUTDIR" "$CHECKPOINT_DIR"

# Skip-if-done: if final h5ad exists and has X_scprint, exit cleanly (idempotent)
if [ -f "$OUTPUT" ]; then
    HAS_EMB=$(python -c "
import anndata as ad
try:
    a = ad.read_h5ad('$OUTPUT')
    print('yes' if 'X_scprint' in a.obsm else 'no')
except Exception:
    print('no')
" 2>/dev/null)
    if [ "$HAS_EMB" = "yes" ]; then
        echo "Output already complete: $OUTPUT (X_scprint present). Nothing to do."
        exit 0
    fi
fi

python "$REPO/scripts/foundation_models/scPRINT/generate_scprint_embeddings.py" \
    --input  "$INPUT" \
    --output "$OUTPUT" \
    --ckpt   auto \
    --batch-size 64 \
    --max-len 2000 \
    --checkpoint-dir "$CHECKPOINT_DIR" \
    --chunk-size 5000

echo "Output: $OUTPUT"
ls -lh "$OUTPUT"
echo "Completed: $(date)"
