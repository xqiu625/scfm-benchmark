#!/bin/bash
#SBATCH --job-name=0006_subset_aux
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/subset_aux_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/subset_aux_%j.err
#SBATCH --mail-user=xqiu014@ucr.edu
#SBATCH --mail-type=END,FAIL
set -euo pipefail
echo "=== auxiliary tables on COVID+Sepsis matched cohort (roadmap D) ==="; echo "Start: $(date)"
ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PATH="$ENV_DIR/bin:$PATH"; export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="$ENV_DIR"; export PYTHONNOUSERSITE=1
cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

C=scripts/24_option_c_canonical_input/code
EMB=/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c_covid_sepsis
R=/bigdata/godziklab/shared/Xinru/302006/05_RESULTS

echo "--- 1. filter to COVID+Sepsis ---"
python scripts/30_review_response_stats/code/filter_to_covid_sepsis.py

echo "--- 2. multiclassifier (subset) ---"
python $C/multiclassifier_option_c.py --embedding-dir $EMB --output-dir $R/option_c_multiclassifier_subset

echo "--- 3. embedding quality / scIB (subset) ---"
python $C/embedding_quality_option_c.py --embedding-dir $EMB --output-dir $R/option_c_embedding_quality_subset

echo "--- 4. activation axis (subset) ---"
python $C/score_activation_markers.py --embedding-dir $EMB \
    --raw-input /bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/human_platelet_covid_severity_raw_counts.h5ad \
    --output-dir $R/option_c_activation_markers_subset

echo "--- 5. ablation / matched-dim (subset, incl scGPT) ---"
python $C/ablation_option_c.py --embedding-dir $EMB --output-dir $R/option_c_ablation_subset

echo "--- 6. fusion (subset) ---"
python $C/fusion_option_c.py --embedding-dir $EMB --output-dir $R/option_c_fusion_subset

echo "Completed: $(date)"
