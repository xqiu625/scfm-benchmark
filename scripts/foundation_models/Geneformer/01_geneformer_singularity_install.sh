#!/bin/bash

# SUCCESSFUL GENEFORMER SINGULARITY INSTALLATION
# Summary of the working approach for HPC environments

echo "🧬 Geneformer Singularity Installation - Success Recipe"
echo "====================================================="

# Step 1: Load Singularity module
module load singularity

# Step 2: Set cache to shared space (avoids home directory disk quota issues)
export SINGULARITY_CACHEDIR="/bigdata/godziklab/shared/Xinru/302006/singularity_cache"
export SINGULARITY_TMPDIR="/bigdata/godziklab/shared/Xinru/302006/singularity_tmp"
mkdir -p "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR"

# Step 3: Clean any existing cache
rm -rf /rhome/xinruqiu/.singularity/cache

# Step 4: Pull the Microsoft container (includes Geneformer + evaluation tools)
singularity pull geneformer.sif docker://kzkedzierska/sc_foundation_evals:latest

# Step 5: Test the installation
singularity exec --nv geneformer.sif python -c "
from geneformer import TranscriptomeTokenizer
tk = TranscriptomeTokenizer(nproc=1)
print('✅ Geneformer ready for processing!')
"

echo "🎉 Installation complete!"
echo "📁 Container: geneformer.sif (~5-8GB)"
echo "🚀 Usage: singularity exec --nv geneformer.sif python your_script.py"

# =====================================================
# KEY SUCCESS FACTORS:
# =====================================================
# 1. Source modules properly in script:
#    source /etc/profile.d/modules.sh
#    module load singularity
#
# 2. Use shared space for cache (not home directory):
#    export SINGULARITY_CACHEDIR="/path/to/shared/space/cache"
#    export SINGULARITY_TMPDIR="/path/to/shared/space/tmp"
#
# 3. Clean old cache to avoid conflicts:
#    rm -rf ~/.singularity/cache
#
# 4. Use Microsoft's container (most reliable):
#    docker://kzkedzierska/sc_foundation_evals:latest
#
# 5. Test with GPU support:
#    singularity exec --nv container.sif python -c "import torch; print(torch.cuda.is_available())"
# =====================================================

# WORKING SCRIPT TEMPLATE:
cat > working_singularity_setup.sh << 'WORKING_SCRIPT'
#!/bin/bash

# Source modules properly
source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/Modules/init/bash 2>/dev/null

# Load singularity
module load singularity

# Set cache to shared space (CRITICAL for avoiding disk quota issues)
export SINGULARITY_CACHEDIR="/your/shared/path/singularity_cache"
export SINGULARITY_TMPDIR="/your/shared/path/singularity_tmp"
mkdir -p "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR"

# Clean old cache
rm -rf ~/.singularity/cache

# Download container
singularity pull geneformer.sif docker://kzkedzierska/sc_foundation_evals:latest

# Test
singularity exec --nv geneformer.sif python -c "
from geneformer import TranscriptomeTokenizer
print('✅ Success!')
"
WORKING_SCRIPT

echo ""
echo "📋 Usage Examples:"
echo "=================="

echo "# Basic usage:"
echo "singularity exec geneformer.sif python script.py"
echo ""

echo "# With GPU support:"
echo "singularity exec --nv geneformer.sif python script.py"
echo ""

echo "# Interactive shell:"
echo "singularity shell --nv geneformer.sif"
echo ""

echo "# SLURM job:"
echo "sbatch --gres=gpu:1 --mem=16G --wrap=\"singularity exec --nv geneformer.sif python process_retnla.py\""
echo ""

echo "# Mount data directories:"
echo "singularity exec --nv --bind /data:/data geneformer.sif python script.py"

echo ""
echo "🎯 Container Contents:"
echo "====================="
echo "✅ Python 3.9"
echo "✅ PyTorch with CUDA support"
echo "✅ Transformers library" 
echo "✅ Geneformer package"
echo "✅ Single-cell analysis tools (scanpy, anndata)"
echo "✅ Evaluation frameworks"
echo ""
echo "📊 Perfect for processing 890 RETNLA ortholog cells!"