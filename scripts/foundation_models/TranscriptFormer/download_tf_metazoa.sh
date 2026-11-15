#!/bin/bash
#SBATCH --job-name=tf_metazoa_dl
#SBATCH --partition=batch
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/tf_metazoa_dl_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/tf_metazoa_dl_%j.err
#SBATCH --mail-user=xqiu014@ucr.edu
#SBATCH --mail-type=END,FAIL

# Download TF-Metazoa weights for cross-species TranscriptFormer (mouse + invertebrate inference).
# TF-Metazoa was trained on 1.5B years of evolution across 12 species including mouse.
# Source: TranscriptFormer GitHub (czi-ai/transcriptformer).

set -euo pipefail
echo "=== TF-Metazoa weights download ==="
echo "Start: $(date)"

# Activate transcriptformer env to use its CLI for download
ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/transcriptformer
if [ -x "$ENV_DIR/bin/python" ]; then
    export PATH="$ENV_DIR/bin:$PATH"
    export PYTHONNOUSERSITE=1
    echo "Python: $(which python)"
fi

REPO_DIR=/bigdata/godziklab/shared/Xinru/302006/302006_updategit
CHECKPOINT_DIR="$REPO_DIR/checkpoints"
mkdir -p "$CHECKPOINT_DIR"

# Clean up empty placeholders from prior failed run
rm -rf "$CHECKPOINT_DIR/tf_metazoa" 2>/dev/null

# The transcriptformer CLI uses positional model name (matches the tf_sapiens
# download pattern in install_transcriptformer_conda.sh):
#   transcriptformer download tf-sapiens   (existing, works)
#   transcriptformer download tf-metazoa   (new, expected)
#
# It writes to ./checkpoints/<model_name>/ relative to CWD.

cd "$REPO_DIR"
echo "Working dir: $(pwd)"
echo "Running: transcriptformer download tf-metazoa"
transcriptformer download tf-metazoa

echo ""
echo "=== Download status ==="
ls -lh "$CHECKPOINT_DIR/tf_metazoa" 2>/dev/null || ls -lh "$CHECKPOINT_DIR"
echo ""
echo "Done: $(date)"
