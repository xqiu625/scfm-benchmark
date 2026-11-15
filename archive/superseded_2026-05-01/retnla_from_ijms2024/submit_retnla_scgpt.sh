#!/bin/bash
#SBATCH --job-name=retnla_scgpt
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/logs/retnla_scgpt_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/logs/retnla_scgpt_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xinru.qiu@ucr.edu

echo "RETNLA scGPT Processing Job"
echo "=========================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "GPU: $CUDA_VISIBLE_DEVICES"
echo "Start: $(date)"

cd /bigdata/godziklab/shared/Xinru/302006
source activate /bigdata/godziklab/shared/Xinru/302006_env


# Create logs directory if it doesn't exist
mkdir -p logs

# Load environment
module load singularity 2>/dev/null || echo "Using system singularity"

# Set GPU environment
export CUDA_VISIBLE_DEVICES=${SLURM_LOCALID:-0}
echo "CUDA devices: $CUDA_VISIBLE_DEVICES"

# Run RETNLA scGPT processing
echo "Starting RETNLA scGPT processing..."
python process_retnla_scgpt.py

exit_code=$?

echo "Job completed with exit code: $exit_code"
echo "End: $(date)"

# Check results
if [ $exit_code -eq 0 ]; then
    echo "SUCCESS: RETNLA scGPT processing completed"
    echo "Checking results..."
    ls -la /bigdata/godziklab/shared/Xinru/302006/retnla/scgpt_*
else
    echo "FAILED: RETNLA scGPT processing failed"
fi

exit $exit_code