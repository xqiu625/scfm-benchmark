#!/bin/bash
#SBATCH --job-name=0006_scgpt_aux
#SBATCH --partition=batch
#SBATCH --array=0-4
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=30-00:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/0006_scgpt_aux_%A_%a.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/0006_scgpt_aux_%A_%a.err
#SBATCH --mail-user=xqiu014@ucr.edu
#SBATCH --mail-type=END,FAIL

# scGPT auxiliary-evaluation top-up (stage 30).
#
# The scGPT Option C canonical embedding
# (02_EMBEDDINGS/option_c/scgpt_covid_canonical.h5ad, obsm X_scgpt, 47,977 cells)
# landed 2026-06-08, after the auxiliary analyses below were first run, so scGPT
# shows "---" in five manuscript tables. Each array task recomputes exactly one
# missing scGPT evaluation under the identical protocol as the other FMs:
#
#   task 0  adaptation curve (COVID+Sepsis, bootstrap)  -> tab:adaptation_curve
#   task 1  sample efficiency (5-disease pool, 5-fold)  -> tab:sample_eff
#   task 2  gene dropout robustness (5-disease pool)    -> tab:dropout
#   task 3  embedding ablation (5-disease pool)         -> tab:ablation
#   task 4  multimodal fusion refit, primary 5-FM set   -> tab:fusion
#
# All tasks are CPU-only sklearn jobs on frozen embeddings (no GPU).
# See scripts/30_review_response_stats/SCGPT_AUXILIARY_RERUNS.md for details.
#
# Submit:   sbatch scripts/30_review_response_stats/code/submit_scgpt_auxiliary_reruns.sh
# Single:   sbatch --array=2 scripts/30_review_response_stats/code/submit_scgpt_auxiliary_reruns.sh

set -euo pipefail

ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
if [ ! -x "$ENV_DIR/bin/python" ]; then echo "ERROR: env not found: $ENV_DIR"; exit 1; fi
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="$ENV_DIR"
export PYTHONNOUSERSITE=1
echo "Task $SLURM_ARRAY_TASK_ID | Python: $(which python) | Start: $(date)"

REPO=/bigdata/godziklab/shared/Xinru/302006/302006_updategit
cd "$REPO"

HPCC=/bigdata/godziklab/shared/Xinru/302006
EMB=$HPCC/02_EMBEDDINGS/option_c
RES=$HPCC/05_RESULTS

# Pre-flight: the whole point of this batch is the 2026-06-08 scGPT embedding.
SCGPT_H5AD=$EMB/scgpt_covid_canonical.h5ad
if [ ! -f "$SCGPT_H5AD" ]; then
    echo "ERROR: scGPT Option C embedding not found: $SCGPT_H5AD"; exit 1
fi

C24=scripts/24_option_c_canonical_input/code

case "$SLURM_ARRAY_TASK_ID" in
  0)
    echo "=== [0] Adaptation curve (bootstrap), scGPT only -> tab:adaptation_curve ==="
    # Merges scGPT rows into the existing all_rows.csv / summary.csv in-place
    # (writes scGPT.csv; other methods' rows are preserved).
    python3 -u scripts/29_adaptation_curve/code/run_adaptation_curve_bootstrap.py \
        --models scGPT
    ;;
  1)
    echo "=== [1] Sample efficiency, scGPT only -> tab:sample_eff ==="
    # Separate output dir so the canonical primary-benchmark CSVs are untouched;
    # append the scgpt rows to $RES/option_c_primary_benchmark/sample_efficiency.csv
    # when filling the table.
    mkdir -p "$RES/option_c_primary_benchmark/scgpt_aux"
    python3 -u $C24/primary_benchmark_option_c.py \
        --embedding-dir "$EMB" \
        --output-dir    "$RES/option_c_primary_benchmark/scgpt_aux" \
        --models scgpt --analyses sample_efficiency
    ;;
  2)
    echo "=== [2] Gene dropout robustness, scGPT only -> tab:dropout ==="
    # Same command as the original 6-FM run (Job 23851383), restricted to scGPT.
    # Timestamped outputs; nothing is overwritten.
    python3 -u scripts/07_gene_dropout/code/gene_dropout_robustness.py \
        --tasks binary --models scGPT
    ;;
  3)
    echo "=== [3] Embedding ablation, scGPT only -> tab:ablation ==="
    # Writes per_fm_scgpt.csv and merges scgpt rows into the canonical
    # ablation_summary.csv in-place (other FMs' rows are preserved).
    python3 -u $C24/ablation_option_c.py \
        --embedding-dir "$EMB" \
        --output-dir    "$RES/option_c_ablation" \
        --models scgpt
    ;;
  4)
    echo "=== [4] Multimodal fusion refit, primary 5-FM set incl. scGPT -> tab:fusion ==="
    # Fusion cannot add a single model: the combination strategies use every FM,
    # so the primary 5-FM set (STATE, UCE, Geneformer, scPRINT-2, scGPT) is refit
    # together on a fresh Cell_ID intersection. New output dir; the legacy 6-FM
    # results in $RES/option_c_fusion are preserved.
    mkdir -p "$RES/option_c_fusion_primary5"
    python3 -u $C24/fusion_option_c.py \
        --embedding-dir "$EMB" \
        --output-dir    "$RES/option_c_fusion_primary5" \
        --models geneformer scgpt scprint2 state uce
    ;;
  *)
    echo "ERROR: unexpected SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"; exit 1
    ;;
esac

echo "Task $SLURM_ARRAY_TASK_ID done: $(date)"
