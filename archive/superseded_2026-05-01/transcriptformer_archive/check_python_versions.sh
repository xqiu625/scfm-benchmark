#!/bin/bash
# Check available Python versions on cluster

echo "=========================================="
echo "CHECKING PYTHON VERSIONS"
echo "=========================================="

echo ""
echo "Current default Python:"
python3 --version

echo ""
echo "Available Python modules:"
module avail python 2>&1 | grep -i python

echo ""
echo "Checking conda environments:"
which conda && conda env list

echo ""
echo "=========================================="
echo "RECOMMENDATION:"
echo "=========================================="
echo ""
echo "TranscriptFormer requires Python 3.11+"
echo ""
echo "If you have a Python 3.11+ module:"
echo "  module load python/3.11"
echo "  bash scripts/foundation_models/TranscriptFormer/install_transcriptformer.sh"
echo ""
echo "If using conda:"
echo "  conda create -n transcriptformer python=3.11"
echo "  conda activate transcriptformer"
echo "  # Then install packages manually"
echo ""
echo "If neither available:"
echo "  Use Singularity container (I can create the definition)"
