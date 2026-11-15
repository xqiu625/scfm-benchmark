#!/bin/bash
#SBATCH --job-name=gene_dropout
#SBATCH --partition=batch
#SBATCH --mem=64gb
#SBATCH --cpus-per-task=8
#SBATCH --time=6:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/gene_dropout_%j.log
#SBATCH --exclude=gpu10,gpu08

echo "============================================================"
echo "GENE DROPOUT ROBUSTNESS ANALYSIS"
echo "Start time: $(date)"
echo "Node: $(hostname)"
echo "============================================================"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

# PATH-prepend pattern (avoid silent-fail in batch shells without conda init)
ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export PYTHONNOUSERSITE=1

# 6-FM Option C re-run (2026-05-05): STATE STACK Geneformer scPRINT2 TranscriptFormer UCE
python scripts/07_gene_dropout/code/gene_dropout_robustness.py --tasks binary --models STATE STACK Geneformer scPRINT2 TranscriptFormer UCE

echo "============================================================"
echo "Completed at: $(date)"
echo "============================================================"
