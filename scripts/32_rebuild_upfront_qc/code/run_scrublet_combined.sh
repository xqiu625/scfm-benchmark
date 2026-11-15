# RETIRED: scrublet_filter_combined.py is superseded by gate2_run.py.
# Kept so an old submission fails loudly instead of running old semantics.
#!/bin/bash
#SBATCH --job-name=scrublet_combined
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/scrublet_combined_%j.log
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00

ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export PYTHONNOUSERSITE=1

python3 /bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/32_rebuild_upfront_qc/code/scrublet_filter_combined.py \
  --cohort-dir "$1" --out-dir "$2"
