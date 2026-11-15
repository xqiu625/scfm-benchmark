#!/bin/bash
#SBATCH --job-name=tf_regen
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/transcriptformer_regenerate_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/transcriptformer_regenerate_%j.log
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00

echo "============================================================"
echo "TRANSCRIPTFORMER EMBEDDING RE-GENERATION"
echo "============================================================"
echo "REASON: Previous embeddings contain 18% synthetic fallback data"
echo "============================================================"
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'N/A')"
echo "============================================================"

# Activate environment
source ~/.bashrc
conda activate /bigdata/godziklab/shared/Xinru/302006/state_env

# Verify environment
echo "Python: $(which python)"
python -c "import torch; print(f'PyTorch: {torch.__version__}, CUDA: {torch.cuda.is_available()}')"

# Navigate to TranscriptFormer directory
cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/foundation_models/TranscriptFormer

# Download model if not present
CHECKPOINT_DIR="./checkpoints/tf_sapiens"
if [ ! -d "$CHECKPOINT_DIR" ]; then
    echo ""
    echo "Downloading TranscriptFormer tf-sapiens model..."
    transcriptformer download tf-sapiens
    echo "Download complete."
else
    echo "Model checkpoint already exists at: $CHECKPOINT_DIR"
fi

# Verify model exists
if [ ! -d "$CHECKPOINT_DIR" ]; then
    echo "ERROR: Model download failed. Cannot proceed."
    exit 1
fi

# Backup old embeddings first
BACKUP_DIR="/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/transcriptformer/backup_$(date +%Y%m%d_%H%M%S)"
echo "Backing up old embeddings to: $BACKUP_DIR"
mkdir -p "$BACKUP_DIR"
cp /bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/transcriptformer/*.h5ad "$BACKUP_DIR/" 2>/dev/null || echo "No existing embeddings to backup"

# Run embedding generation (COVID only - the contaminated one)
echo ""
echo "Generating new embeddings..."
python generate_transcriptformer_embeddings.py --dataset covid --batch-size 16

echo ""
echo "============================================================"
echo "TRANSCRIPTFORMER RE-GENERATION COMPLETE"
echo "Completed at: $(date)"
echo "============================================================"
echo ""
echo "NEXT STEPS:"
echo "1. Run embedding integrity check to verify new embeddings are clean"
echo "2. Re-run analyses that used TranscriptFormer"
echo "============================================================"
