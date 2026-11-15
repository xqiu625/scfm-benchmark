#!/bin/bash
#SBATCH --job-name=sample_eff_6fm
#SBATCH --partition=batch
#SBATCH --mem=64gb
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/sample_eff_6fm_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/sample_eff_6fm_%j.err
#SBATCH --mail-user=xqiu014@ucr.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail
echo "=== Sample efficiency 6-FM Option C re-run ==="
echo "Start: $(date)"

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export PYTHONNOUSERSITE=1

python scripts/05_sample_efficiency/code/subsampling_robustness.py --tasks binary

echo "Completed: $(date)"
