#!/bin/bash
# =============================================================================
# Script: download_results_from_hpcc.sh
# Purpose: Download results and figures from HPCC to local scripts/ stage folders
# Usage: Run this LOCALLY (on your Mac)
# =============================================================================

HPCC_USER="xinruqiu"
HPCC_HOST="cluster.hpcc.ucr.edu"
HPCC_BASE="/bigdata/godziklab/shared/Xinru/302006"
LOCAL_BASE="/Users/xinruqiu/Desktop/302006/302006_updategit/scripts"

echo "=============================================="
echo "Downloading results from HPCC to local"
echo "=============================================="

# -----------------------------------------------------------------------------
# 01_embedding_integrity
# -----------------------------------------------------------------------------
echo "[01] Embedding Integrity..."
scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/05_RESULTS/embedding_integrity_report.csv \
    ${LOCAL_BASE}/01_embedding_integrity/results/ 2>/dev/null && echo "  ✓ results"

# -----------------------------------------------------------------------------
# 02_embedding_classification
# -----------------------------------------------------------------------------
echo "[02] Fine-tuned Benchmark..."

# Results
scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/05_RESULTS/covid_severity_benchmark/benchmark_STATE_20251212_164146.csv \
    ${LOCAL_BASE}/02_embedding_classification/results/benchmark_STATE.csv 2>/dev/null
scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/05_RESULTS/covid_severity_benchmark/benchmark_UCE_20251212_032014.csv \
    ${LOCAL_BASE}/02_embedding_classification/results/benchmark_UCE.csv 2>/dev/null
scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/05_RESULTS/covid_severity_benchmark/benchmark_scGPT_20251212_013208.csv \
    ${LOCAL_BASE}/02_embedding_classification/results/benchmark_scGPT.csv 2>/dev/null
scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/05_RESULTS/covid_severity_benchmark/benchmark_Geneformer_20251212_171449.csv \
    ${LOCAL_BASE}/02_embedding_classification/results/benchmark_Geneformer.csv 2>/dev/null
scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/05_RESULTS/covid_severity_benchmark/benchmark_TranscriptFormer_20251212_114102.csv \
    ${LOCAL_BASE}/02_embedding_classification/results/benchmark_TranscriptFormer.csv 2>/dev/null
echo "  ✓ results"

# Figures
scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/04_FIGURES/fig_model_rankings.{png,pdf} \
    ${LOCAL_BASE}/02_embedding_classification/figures/ 2>/dev/null
scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/04_FIGURES/fig_roc_curves.{png,pdf} \
    ${LOCAL_BASE}/02_embedding_classification/figures/ 2>/dev/null
scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/04_FIGURES/fig_confusion_matrix.{png,pdf} \
    ${LOCAL_BASE}/02_embedding_classification/figures/ 2>/dev/null
echo "  ✓ figures"

# -----------------------------------------------------------------------------
# 03_zeroshot_benchmark
# -----------------------------------------------------------------------------
echo "[03] Zero-shot Benchmark..."
scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/04_FIGURES/fig_zeroshot_vs_finetuned.{png,pdf} \
    ${LOCAL_BASE}/03_zeroshot_benchmark/figures/ 2>/dev/null && echo "  ✓ figures"

# -----------------------------------------------------------------------------
# 05_sample_efficiency
# -----------------------------------------------------------------------------
echo "[05] Sample Efficiency..."
scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/05_RESULTS/covid_severity_benchmark/subsampling_robustness_20251227_215149.csv \
    ${LOCAL_BASE}/05_sample_efficiency/results/subsampling_robustness.csv 2>/dev/null
scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/05_RESULTS/covid_severity_benchmark/subsampling_robustness_20251227_215149.json \
    ${LOCAL_BASE}/05_sample_efficiency/results/subsampling_robustness.json 2>/dev/null
echo "  ✓ results"
scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/04_FIGURES/fig_sample_efficiency.{png,pdf} \
    ${LOCAL_BASE}/05_sample_efficiency/figures/ 2>/dev/null && echo "  ✓ figures"

# -----------------------------------------------------------------------------
# 06_batch_generalization
# -----------------------------------------------------------------------------
echo "[06] Batch Generalization..."
scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/05_RESULTS/covid_severity_benchmark/batch_shift_generalization_20251226_182340.csv \
    ${LOCAL_BASE}/06_batch_generalization/results/batch_shift_generalization.csv 2>/dev/null
scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/05_RESULTS/covid_severity_benchmark/batch_shift_generalization_20251226_182340.json \
    ${LOCAL_BASE}/06_batch_generalization/results/batch_shift_generalization.json 2>/dev/null
echo "  ✓ results"
scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/04_FIGURES/fig_batch_generalization.{png,pdf} \
    ${LOCAL_BASE}/06_batch_generalization/figures/ 2>/dev/null && echo "  ✓ figures"

# -----------------------------------------------------------------------------
# 07_gene_dropout
# -----------------------------------------------------------------------------
echo "[07] Gene Dropout..."
scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/04_FIGURES/fig_gene_dropout.{png,pdf} \
    ${LOCAL_BASE}/07_gene_dropout/figures/ 2>/dev/null && echo "  ✓ figures"

# -----------------------------------------------------------------------------
# 08_statistical_significance
# -----------------------------------------------------------------------------
echo "[08] Statistical Significance..."

# Bootstrap CI results
for suffix in 125319 125332 125348 125356; do
    scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/05_RESULTS/covid_severity_benchmark/bootstrap_ci_fast_20260109_${suffix}.csv \
        ${LOCAL_BASE}/08_statistical_significance/results/ 2>/dev/null
    scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/05_RESULTS/covid_severity_benchmark/delong_tests_fast_20260109_${suffix}.csv \
        ${LOCAL_BASE}/08_statistical_significance/results/ 2>/dev/null
    scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/05_RESULTS/covid_severity_benchmark/fold_aucs_fast_20260109_${suffix}.csv \
        ${LOCAL_BASE}/08_statistical_significance/results/ 2>/dev/null
done
echo "  ✓ results"

scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/04_FIGURES/fig_statistical_significance.{png,pdf} \
    ${LOCAL_BASE}/08_statistical_significance/figures/ 2>/dev/null && echo "  ✓ figures"

# -----------------------------------------------------------------------------
# 09_embedding_ablation
# -----------------------------------------------------------------------------
echo "[09] Embedding Ablation..."
scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/05_RESULTS/covid_severity_benchmark/embedding_ablations_fast_20260110_015800.csv \
    ${LOCAL_BASE}/09_embedding_ablation/results/embedding_ablations.csv 2>/dev/null
scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/05_RESULTS/covid_severity_benchmark/embedding_ablations_fast_20260110_015800.json \
    ${LOCAL_BASE}/09_embedding_ablation/results/embedding_ablations.json 2>/dev/null
echo "  ✓ results"

scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/04_FIGURES/fig_robustness_summary.{png,pdf} \
    ${LOCAL_BASE}/09_embedding_ablation/figures/ 2>/dev/null && echo "  ✓ figures"

# -----------------------------------------------------------------------------
# 10_clinical_metrics
# -----------------------------------------------------------------------------
echo "[10] Clinical Metrics..."
scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/05_RESULTS/covid_severity_benchmark/clinical_metrics_20260102_224502.csv \
    ${LOCAL_BASE}/10_clinical_metrics/results/clinical_metrics.csv 2>/dev/null
scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/05_RESULTS/covid_severity_benchmark/clinical_metrics_20260102_224502.json \
    ${LOCAL_BASE}/10_clinical_metrics/results/clinical_metrics.json 2>/dev/null
echo "  ✓ results"

scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/04_FIGURES/fig_clinical_metrics.{png,pdf} \
    ${LOCAL_BASE}/10_clinical_metrics/figures/ 2>/dev/null && echo "  ✓ figures"

# -----------------------------------------------------------------------------
# 11_interpretability
# -----------------------------------------------------------------------------
echo "[11] Interpretability..."
scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/05_RESULTS/gene_programs/disease_programs_STATE_20251231_014342.csv \
    ${LOCAL_BASE}/11_interpretability/results/disease_programs_STATE.csv 2>/dev/null
scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/05_RESULTS/gene_programs/disease_programs_UCE_20251231_014342.csv \
    ${LOCAL_BASE}/11_interpretability/results/disease_programs_UCE.csv 2>/dev/null
scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/05_RESULTS/gene_programs/program_comparison_20251231_014342.csv \
    ${LOCAL_BASE}/11_interpretability/results/program_comparison.csv 2>/dev/null
echo "  ✓ results"

scp ${HPCC_USER}@${HPCC_HOST}:${HPCC_BASE}/04_FIGURES/fig_biomarkers.{png,pdf} \
    ${LOCAL_BASE}/11_interpretability/figures/ 2>/dev/null && echo "  ✓ figures"

# -----------------------------------------------------------------------------
# 12_umap_visualization
# -----------------------------------------------------------------------------
echo "[12] UMAP Visualization..."
scp ${HPCC_USER}@${HPCC_HOST}:"${HPCC_BASE}/05_RESULTS/covid_severity_benchmark/figures/umap_5models_*.{png,pdf}" \
    ${LOCAL_BASE}/12_umap_visualization/figures/ 2>/dev/null
scp ${HPCC_USER}@${HPCC_HOST}:"${HPCC_BASE}/05_RESULTS/covid_severity_benchmark/figures/state_severity_gradient.{png,pdf}" \
    ${LOCAL_BASE}/12_umap_visualization/figures/ 2>/dev/null
scp ${HPCC_USER}@${HPCC_HOST}:"${HPCC_BASE}/05_RESULTS/covid_severity_benchmark/figures/uce_severity_gradient.{png,pdf}" \
    ${LOCAL_BASE}/12_umap_visualization/figures/ 2>/dev/null
echo "  ✓ figures"

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
echo ""
echo "=============================================="
echo "Download complete!"
echo "=============================================="
echo ""
echo "Verify with:"
echo "  find ${LOCAL_BASE} -name '*.csv' | wc -l"
echo "  find ${LOCAL_BASE} -name '*.png' | wc -l"
