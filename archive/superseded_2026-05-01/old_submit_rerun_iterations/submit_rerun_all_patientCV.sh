#!/bin/bash
# =============================================================================
# Re-run ALL analyses with patient-level CV fix (StratifiedGroupKFold)
#
# Usage:
#   1. Push updated code to HPCC:
#      rsync -avz --exclude='.git' 302006_updategit/ \
#        cluster.hpcc.ucr.edu:/bigdata/godziklab/shared/Xinru/302006/302006_updategit/
#
#   2. SSH into HPCC:
#      ssh cluster.hpcc.ucr.edu
#
#   3. Submit all jobs:
#      cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
#      bash submit_rerun_all_patientCV.sh
#
#   4. Or submit individual jobs (see below)
#
# Author: Xinru Qiu
# Date: April 2026
# =============================================================================

BASE_DIR="/bigdata/godziklab/shared/Xinru/302006"
CODE_DIR="${BASE_DIR}/302006_updategit"
LOG_DIR="${BASE_DIR}/06_LOGS/patientCV_rerun"
ENV_CLEAN="/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean"
ENV_STATE="${BASE_DIR}/state_env"

mkdir -p "$LOG_DIR"

echo "============================================================"
echo "SUBMITTING ALL PATIENT-LEVEL CV RE-RUNS"
echo "Date: $(date)"
echo "Log dir: $LOG_DIR"
echo "============================================================"

# -----------------------------------------------------------------
# ACC 02: Core Severity Benchmark (benchmark_single_model.py)
# ~30 min, CPU only
# -----------------------------------------------------------------
JOB_02=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=02_severity_patientCV
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=2:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/02_severity_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/02_severity_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 02: Core Severity Benchmark (patient-level CV) ==="
echo "Start: $(date)"

for MODEL in STATE UCE scGPT Geneformer TranscriptFormer; do
    echo "--- $MODEL ---"
    python analyses/02_embedding_classification/code/benchmark_single_model.py \
        --model "$MODEL" --tasks binary 3-class 6-class
done

echo "Completed: $(date)"
SLURM
)
echo "ACC 02 submitted: $JOB_02"

# -----------------------------------------------------------------
# ACC 02b: Full benchmark (covid_severity_benchmark.py)
# All 5 classifiers x 4 models — heavier
# -----------------------------------------------------------------
JOB_02b=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=02b_full_benchmark
#SBATCH --partition=batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --time=6:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/02b_full_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/02b_full_%j.err

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 02b: Full 5-classifier benchmark (patient-level CV) ==="
echo "Start: $(date)"
python analyses/02_embedding_classification/code/covid_severity_benchmark.py --quick
echo "Completed: $(date)"
SLURM
)
echo "ACC 02b submitted: $JOB_02b"

# -----------------------------------------------------------------
# ACC 02c: Ensemble Models
# -----------------------------------------------------------------
JOB_02c=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=02c_ensemble
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=4:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/02c_ensemble_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/02c_ensemble_%j.err

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 02c: Ensemble Models (patient-level CV) ==="
python analyses/02_embedding_classification/code/ensemble_models.py
echo "Completed: $(date)"
SLURM
)
echo "ACC 02c submitted: $JOB_02c"

# -----------------------------------------------------------------
# ACC 03: Zero-shot Benchmark (same code as 02)
# -----------------------------------------------------------------
JOB_03=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=03_zeroshot_patientCV
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=2:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/03_zeroshot_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/03_zeroshot_%j.err

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 03: Zero-shot Benchmark (patient-level CV) ==="
for MODEL in STATE UCE scGPT Geneformer TranscriptFormer; do
    echo "--- $MODEL ---"
    python analyses/03_zeroshot_benchmark/code/benchmark_single_model.py \
        --model "$MODEL" --tasks binary 3-class 6-class
done
echo "Completed: $(date)"
SLURM
)
echo "ACC 03 submitted: $JOB_03"

# -----------------------------------------------------------------
# ACC 04: PCA Baseline
# -----------------------------------------------------------------
JOB_04=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=04_pca_patientCV
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=2:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/04_pca_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/04_pca_%j.err

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 04: PCA Baseline (patient-level CV) ==="
python analyses/04_pca_baseline/code/benchmark_pca_baseline.py
echo "Completed: $(date)"
SLURM
)
echo "ACC 04 submitted: $JOB_04"

# -----------------------------------------------------------------
# ACC 05: Sample Efficiency
# -----------------------------------------------------------------
JOB_05=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=05_sample_eff_patientCV
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=4:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/05_sample_eff_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/05_sample_eff_%j.err

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 05: Sample Efficiency (patient-level CV) ==="
python analyses/05_sample_efficiency/code/subsampling_robustness.py \
    --tasks binary --n-repeats 5
echo "Completed: $(date)"
SLURM
)
echo "ACC 05 submitted: $JOB_05"

# -----------------------------------------------------------------
# ACC 07: Gene Dropout Robustness
# -----------------------------------------------------------------
JOB_07=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=07_dropout_patientCV
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=6:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/07_dropout_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/07_dropout_%j.err
#SBATCH --exclude=gpu10,gpu08

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 07: Gene Dropout Robustness (patient-level CV) ==="
python analyses/07_gene_dropout/code/gene_dropout_robustness.py \
    --tasks binary 3-class --n-repeats 5
echo "Completed: $(date)"
SLURM
)
echo "ACC 07 submitted: $JOB_07"

# -----------------------------------------------------------------
# ACC 08: Statistical Significance (parallel bootstrap)
# -----------------------------------------------------------------
JOB_08=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=08_stats_patientCV
#SBATCH --partition=batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --time=8:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/08_stats_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/08_stats_%j.err

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 08a: Statistical Significance Fast (patient-level CV) ==="
python analyses/08_statistical_significance/code/statistical_significance_fast.py \
    --models STATE UCE TranscriptFormer scGPT Geneformer \
    --tasks binary 3-class \
    --n-bootstrap 100 --pca-dims 256

echo ""
echo "=== ACC 08b: Per-Fold Analysis (patient-level CV) ==="
python analyses/08_statistical_significance/code/per_fold_analysis.py

echo "Completed: $(date)"
SLURM
)
echo "ACC 08 submitted: $JOB_08"

# -----------------------------------------------------------------
# ACC 09: Embedding Ablation
# -----------------------------------------------------------------
JOB_09=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=09_ablation_patientCV
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=6:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/09_ablation_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/09_ablation_%j.err

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 09: Embedding Ablation (patient-level CV) ==="
python analyses/09_embedding_ablation/code/embedding_ablations_fast.py \
    --tasks binary 3-class
echo "Completed: $(date)"
SLURM
)
echo "ACC 09 submitted: $JOB_09"

# -----------------------------------------------------------------
# ACC 10: Clinical Metrics
# -----------------------------------------------------------------
JOB_10=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=10_clinical_patientCV
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=4:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/10_clinical_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/10_clinical_%j.err

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 10: Clinical Metrics (patient-level CV) ==="
python analyses/10_clinical_metrics/code/clinical_metrics_extended.py
echo "Completed: $(date)"
SLURM
)
echo "ACC 10 submitted: $JOB_10"

# -----------------------------------------------------------------
# ACC 11: Interpretability
# -----------------------------------------------------------------
JOB_11=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=11_interp_patientCV
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=4:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/11_interp_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/11_interp_%j.err

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 11a: Interpretability Analysis (patient-level CV) ==="
python analyses/11_interpretability/code/interpretability_analysis.py

echo ""
echo "=== ACC 11b: Gene Programs NMF (patient-level CV) ==="
python analyses/11_interpretability/code/gene_programs_nmf.py

echo "Completed: $(date)"
SLURM
)
echo "ACC 11 submitted: $JOB_11"

# -----------------------------------------------------------------
# ACC 15: Differential Diagnosis
# -----------------------------------------------------------------
JOB_15=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=15_diffdiag_patientCV
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=128gb
#SBATCH --time=8:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/15_diffdiag_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/15_diffdiag_%j.err

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 15a: 5-Model Differential Diagnosis (patient-level CV) ==="
python analyses/15_differential_diagnosis/code/run_15_differential_diagnosis.py

echo ""
echo "=== ACC 15b: scGPT Clean Differential Diagnosis ==="
python analyses/15_differential_diagnosis/code/run_scgpt_differential_diagnosis.py

echo ""
echo "=== ACC 15c: scGPT Clean Finetuned ==="
python analyses/15_differential_diagnosis/code/run_scgpt_finetuned.py

echo "Completed: $(date)"
SLURM
)
echo "ACC 15 submitted: $JOB_15"

# -----------------------------------------------------------------
# ACC 16: scGPT Best Practices Benchmark
# -----------------------------------------------------------------
JOB_16=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=16_scgpt_bp_patientCV
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=2:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/16_scgpt_bp_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/16_scgpt_bp_%j.err

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 16: scGPT Best Practices Benchmark (patient-level CV) ==="
python analyses/16_embedding_generation_5models/code/run_scgpt_bp_benchmark.py
echo "Completed: $(date)"
SLURM
)
echo "ACC 16 submitted: $JOB_16"

# -----------------------------------------------------------------
# ACC 17: Domain-Adapted FM Benchmark
# -----------------------------------------------------------------
JOB_17=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=17_adapted_patientCV
#SBATCH --partition=batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=6:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/17_adapted_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/17_adapted_%j.err

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
source activate /rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PYTHONNOUSERSITE=1

echo "=== ACC 17: Domain-Adapted FM Benchmark (patient-level CV) ==="
python analyses/17_domain_adapted_fm/code/benchmark_adapted_models.py
echo "Completed: $(date)"
SLURM
)
echo "ACC 17 submitted: $JOB_17"

# -----------------------------------------------------------------
# ACC 18: Multimodal Fusion (requires GPU for attention fusion)
# -----------------------------------------------------------------
JOB_18=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=18_fusion_patientCV
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/18_fusion_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/18_fusion_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xqiu014@ucr.edu

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

export ENV_PATH="/bigdata/godziklab/shared/Xinru/302006/state_env"
export PATH="${ENV_PATH}/bin:$PATH"
export CONDA_PREFIX="${ENV_PATH}"
export PYTHONNOUSERSITE=1

echo "=== ACC 18: Multimodal Fusion (patient-level CV) ==="
python -c "import torch; print(f'PyTorch {torch.__version__}, CUDA: {torch.cuda.is_available()}')"

python analyses/18_multimodal_fusion/code/multimodal_fusion.py \
    --strategy all --attention-epochs 100 --common-dim 256
echo "Completed: $(date)"
SLURM
)
echo "ACC 18 submitted: $JOB_18"

# -----------------------------------------------------------------
# ACC 19: Quantization Benchmark (requires GPU + STATE model)
# -----------------------------------------------------------------
JOB_19=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=19_quant_patientCV
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/19_quant_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/19_quant_%j.err

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

export ENV_PATH="/bigdata/godziklab/shared/Xinru/302006/state_env"
export PATH="${ENV_PATH}/bin:$PATH"
export CONDA_PREFIX="${ENV_PATH}"
export PYTHONNOUSERSITE=1

echo "=== ACC 19: Quantization Benchmark (patient-level CV) ==="
python analyses/19_model_optimization/code/quantization_benchmark.py
echo "Completed: $(date)"
SLURM
)
echo "ACC 19 submitted: $JOB_19"

# -----------------------------------------------------------------
# ACC 22: ESM Protein Multimodal (requires GPU for attention)
# -----------------------------------------------------------------
JOB_22=$(sbatch --parsable <<'SLURM'
#!/bin/bash
#SBATCH --job-name=22_esm_patientCV
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=128G
#SBATCH --time=8:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/22_esm_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/patientCV_rerun/22_esm_%j.err

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

export ENV_PATH="/bigdata/godziklab/shared/Xinru/302006/state_env"
export PATH="${ENV_PATH}/bin:$PATH"
export CONDA_PREFIX="${ENV_PATH}"
export PYTHONNOUSERSITE=1

echo "=== ACC 22: 5-Modality Fusion with ESM (patient-level CV) ==="
python analyses/22_esm_protein_multimodal/code/step5_multimodal_fusion_5model.py \
    --strategy all --attention-epochs 100 --common-dim 256
echo "Completed: $(date)"
SLURM
)
echo "ACC 22 submitted: $JOB_22"

# -----------------------------------------------------------------
# Summary
# -----------------------------------------------------------------
echo ""
echo "============================================================"
echo "ALL JOBS SUBMITTED"
echo "============================================================"
echo "  ACC 02  (severity):        $JOB_02"
echo "  ACC 02b (full benchmark):  $JOB_02b"
echo "  ACC 02c (ensemble):        $JOB_02c"
echo "  ACC 03  (zeroshot):        $JOB_03"
echo "  ACC 04  (PCA baseline):    $JOB_04"
echo "  ACC 05  (sample eff):      $JOB_05"
echo "  ACC 07  (gene dropout):    $JOB_07"
echo "  ACC 08  (stat signif):     $JOB_08"
echo "  ACC 09  (ablation):        $JOB_09"
echo "  ACC 10  (clinical):        $JOB_10"
echo "  ACC 11  (interpretability):$JOB_11"
echo "  ACC 15  (diff diagnosis):  $JOB_15"
echo "  ACC 16  (scGPT BP):        $JOB_16"
echo "  ACC 17  (adapted FM):      $JOB_17"
echo "  ACC 18  (fusion, GPU):     $JOB_18"
echo "  ACC 19  (quantization):    $JOB_19"
echo "  ACC 22  (ESM fusion):      $JOB_22"
echo "============================================================"
echo ""
echo "Monitor with:"
echo "  squeue -u xinruqiu"
echo "  watch -n 30 'squeue -u xinruqiu'"
echo ""
echo "Check logs:"
echo "  ls -lt $LOG_DIR/"
echo "  tail -f $LOG_DIR/02_severity_*.log"
echo "============================================================"
