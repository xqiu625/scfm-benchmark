#!/bin/bash
#SBATCH --job-name=bench_tf
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/benchmark_transcriptformer_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/benchmark_transcriptformer_%j.log
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256gb
#SBATCH --time=72:00:00

echo "============================================================"
echo "TRANSCRIPTFORMER BENCHMARK (Clean Embeddings)"
echo "============================================================"
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "============================================================"
echo ""
echo "Running full benchmark suite for TranscriptFormer"
echo "Embeddings: covid_transcriptformer.h5ad (regenerated Jan 2, 2026)"
echo "============================================================"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/analysis

PYTHON=/bigdata/godziklab/shared/Xinru/302006/state_env/bin/python
echo "Using Python: $PYTHON"

# 1. Fine-tuned classification (binary + 3-class)
echo ""
echo "============================================================"
echo "1. FINE-TUNED CLASSIFICATION"
echo "============================================================"
$PYTHON benchmark_single_model.py --model TranscriptFormer --task binary
$PYTHON benchmark_single_model.py --model TranscriptFormer --task 3-class

# 2. Clinical metrics
echo ""
echo "============================================================"
echo "2. CLINICAL METRICS"
echo "============================================================"
$PYTHON clinical_metrics.py --models TranscriptFormer --tasks binary 3-class

# 3. Embedding ablations
echo ""
echo "============================================================"
echo "3. EMBEDDING ABLATIONS"
echo "============================================================"
$PYTHON embedding_ablations.py --models TranscriptFormer --tasks binary 3-class

# 4. Statistical significance
echo ""
echo "============================================================"
echo "4. STATISTICAL SIGNIFICANCE"
echo "============================================================"
$PYTHON statistical_significance.py --models TranscriptFormer --n-bootstrap 500 --tasks binary 3-class

echo ""
echo "============================================================"
echo "TRANSCRIPTFORMER BENCHMARK COMPLETE"
echo "Completed at: $(date)"
echo "============================================================"
echo ""
echo "Results saved to: /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/"
echo "============================================================"
