# Complete scGPT Embedding Errors and Solutions Guide

This document catalogs all errors encountered during scGPT embedding generation and their tested solutions, based on real troubleshooting for foundation model comparison studies.

## Error 1: Conda Environment Creation Failures

### Error Message:
```
LibMambaUnsatisfiableError: Encountered problems while solving:
- package cuda-compiler-12.4.1-hbad6d8a_0 is excluded by strict repo priority
- package cudatoolkit-11.8.0-h4ba93d1_13 is excluded by strict repo priority
```

### Root Cause:
HPC systems with strict repository priorities blocking CUDA and compiler packages from the official scGPT environment YAML file.

### Solutions Tried:
1. ❌ **Force installation**: `conda env create -f protocol_env.yml --force` (command not supported)
2. ❌ **Channel priority modification**: Temporary conda config changes still blocked by system policies
3. ❌ **Manual environment creation**: Dependency conflicts with system-installed packages
4. ✅ **Singularity container**: Downloaded pre-built container bypassing all dependency issues

### Working Solution:
```bash
# Set singularity cache to shared directory
export SINGULARITY_CACHEDIR=/bigdata/godziklab/shared/Xinru/302006/singularity_cache
mkdir -p $SINGULARITY_CACHEDIR

# Download working container
cd /bigdata/godziklab/shared/Xinru/302006/
singularity pull scgpt_0.1.7.sif docker://xueerchen/scgpt:0.1.7
```

## Error 2: Library Compatibility Issues (GLIBCXX)

### Error Message:
```
ImportError: /lib64/libstdc++.so.6: version `GLIBCXX_3.4.29' not found
```

### Root Cause:
Newer pandas/scipy packages compiled with newer GCC requiring GLIBCXX versions not available on older HPC systems.

### Solutions Tried:
1. ❌ **Downgrade pandas**: `pip install pandas==1.5.3` still had GLIBCXX requirements
2. ❌ **Install libstdcxx-ng**: `conda install libstdcxx-ng=12` helped but didn't resolve all issues
3. ❌ **Multiple package downgrades**: Created dependency conflicts
4. ✅ **Pre-built container**: Contained compatible library versions

### Working Solution:
Use the tested Singularity container which has pre-compiled compatible libraries.

## Error 3: Missing IPython Dependency

### Error Message:
```
ModuleNotFoundError: No module named 'IPython'
```

### Root Cause:
scGPT utils expecting IPython for certain operations but not listed in core dependencies.

### Solution:
```bash
pip install ipython
```

## Error 4: Gene Column AssertionError

### Error Message:
```
AssertionError: 
File "/opt/conda/lib/python3.9/site-packages/scgpt/tasks/cell_emb.py", line 189, in embed_data
    assert gene_col in adata.var
```

### Root Cause:
scGPT's `embed_data` function expects gene names in `adata.var[gene_col]` column, not in the index. The function has a hard-coded assertion checking for this column.

### Investigation:
```python
# Our data structure (WRONG for scGPT):
adata.var.index = ['DDX11L2', 'MIR1302-2HG', ...]  # Genes in index
adata.var.columns = ['n_cells', 'mt', 'highly_variable', ...]  # No gene_name column

# What scGPT expects (CORRECT):
adata.var['gene_name'] = ['DDX11L2', 'MIR1302-2HG', ...]  # Genes in column
```

### Solutions Tried:
1. ❌ **No gene_col parameter**: `embed_data(adata, model_dir)` - still failed, default gene_col expected
2. ❌ **gene_col=None**: Still triggered assertion error
3. ❌ **Using existing columns**: Tried `feature_name`, none existed
4. ✅ **Add gene_name column**: Create the required column from index

### Working Solution:
```python
def prepare_data_for_scgpt(adata):
    """Prepare AnnData for scGPT"""
    adata_copy = adata.copy()
    # Add required gene_name column that scGPT expects
    adata_copy.var['gene_name'] = adata_copy.var.index.astype(str)
    return adata_copy

# Usage:
adata = prepare_data_for_scgpt(adata)
embed_data(adata, model_dir, gene_col='gene_name', ...)
```

### Evidence from Benchmark Code:
Found in perturbation benchmark code that scGPT always expects:
```python
adata.var["gene_name"].values.tolist()  # Always accessing this column
```

## Error 5: Flash Attention Data Type Error

### Error Message:
```
AssertionError: qkv.dtype in [torch.float16, torch.bfloat16]
File "/opt/conda/lib/python3.9/site-packages/flash_attn/flash_attention.py", line 36
    assert qkv.dtype in [torch.float16, torch.bfloat16]
```

### Root Cause:
Flash attention (CUDA optimization) requires specific half-precision data types but scGPT was passing float32 tensors. H100 GPU compatibility issue with PyTorch version in container.

### Investigation:
```python
# Error occurred during model inference:
# 1. Gene matching successful: "match 23312/36500 genes in vocabulary"
# 2. Model loading successful: Got to embedding stage
# 3. Flash attention failed: Data type assertion
```

### GPU Compatibility Issue:
```
NVIDIA H100 NVL with CUDA capability sm_90 is not compatible with the current PyTorch installation.
The current PyTorch install supports CUDA capabilities sm_37 sm_50 sm_60 sm_61 sm_70 sm_75 sm_80 sm_86
```

### Solutions Tried:
1. ❌ **Convert to float16**: Data conversion broke other parts of pipeline
2. ❌ **Disable mixed precision**: No clear API to disable
3. ✅ **Use CPU device**: Bypassed GPU-specific flash attention
4. ✅ **Environment variable**: `os.environ['SCGPT_DISABLE_FLASH_ATTN'] = '1'`

### Working Solution:
```python
import os
os.environ['SCGPT_DISABLE_FLASH_ATTN'] = '1'  # Disable flash attention

# OR use CPU device
embed_data(adata, model_dir, device='cpu', ...)
```

## Error 6: CUDA Capability Warnings

### Error Message:
```
UserWarning: NVIDIA H100 NVL with CUDA capability sm_90 is not compatible with the current PyTorch installation.
```

### Root Cause:
PyTorch version in container (compiled for older CUDA architectures) doesn't support H100's compute capability 9.0.

### Solution:
Not a blocking error - just a warning. Use CPU device to avoid GPU compatibility issues:
```python
device = 'cpu'  # Avoid GPU compatibility issues
```

## Error 7: Disk Quota Exceeded (Singularity Cache)

### Error Message:
```
FATAL: disk quota exceeded
write /rhome/xinruqiu/.singularity/cache/blob/oci-put-blob: disk quota exceeded
```

### Root Cause:
Singularity cache in home directory exceeded quota when downloading large container.

### Solution:
```bash
# Redirect cache to shared directory with more space
export SINGULARITY_CACHEDIR=/bigdata/godziklab/shared/Xinru/302006/singularity_cache
mkdir -p $SINGULARITY_CACHEDIR
```

## Complete Working Solution Summary

### Final Working Script Structure:
```python
#!/usr/bin/env python3
import os
os.environ['SCGPT_DISABLE_FLASH_ATTN'] = '1'  # Fix flash attention

import scgpt as scg
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

def prepare_data_for_scgpt(adata):
    """Fix gene column issue"""
    adata_copy = adata.copy()
    adata_copy.var['gene_name'] = adata_copy.var.index.astype(str)  # Add required column
    adata_copy.X = adata_copy.X.astype('float32')  # Ensure proper dtype
    return adata_copy

def generate_scgpt_embeddings(input_path, output_path, model_dir):
    adata = sc.read_h5ad(input_path)
    adata = prepare_data_for_scgpt(adata)  # Fix data structure
    
    from scgpt.tasks import embed_data
    ref_embed_adata = embed_data(
        adata,
        str(model_dir),
        gene_col='gene_name',  # Use the column we created
        obs_to_save=['cell_type', 'batch'],
        batch_size=16,  # Conservative batch size
        device='cpu'    # Avoid GPU compatibility issues
    )
    
    ref_embed_adata.write_h5ad(output_path)
    return True
```

### Environment Setup:
```bash
# Use Singularity container (most reliable)
export SINGULARITY_CACHEDIR=/bigdata/godziklab/shared/Xinru/302006/singularity_cache
mkdir -p $SINGULARITY_CACHEDIR
cd /bigdata/godziklab/shared/Xinru/302006/
singularity pull scgpt_0.1.7.sif docker://xueerchen/scgpt:0.1.7

# Run with container
singularity exec --nv scgpt_0.1.7.sif python generate_scgpt_embeddings_final.py
```

## Key Lessons Learned

### 1. Data Structure Requirements
- scGPT requires `adata.var['gene_name']` column, not just genes in index
- This is poorly documented but critical for `embed_data` function
- Check benchmark code for real-world usage patterns

### 2. Hardware Compatibility
- Modern GPUs (H100) may not be compatible with older PyTorch in containers
- CPU fallback is reliable for embedding generation
- Flash attention can be disabled via environment variables

### 3. Dependency Management
- HPC systems with strict policies make custom environments difficult
- Pre-built containers bypass most dependency issues
- Singularity cache location matters for disk quota

### 4. Debugging Strategy
- Start with smallest possible test case (50 cells)
- Add `gene_name` column to fix assertion errors
- Use CPU device to avoid GPU compatibility issues
- Check benchmark code for real implementation patterns

### 5. Foundation Model Comparison Context
- These solutions enable reliable scGPT embedding generation
- Critical for UCE vs scGPT vs Geneformer comparisons
- CPU performance acceptable for embedding generation (inference only)

## Alternative Approaches (If Above Fails)

### 1. Direct Model Loading
```python
from scgpt.model import TransformerModel
from scgpt.tokenizer.gene_tokenizer import GeneVocab

# Load model components directly
vocab = GeneVocab.from_file(vocab_file)
model = TransformerModel.load_from_checkpoint(model_file)
# Custom preprocessing and inference pipeline
```

### 2. HuggingFace Integration
```python
from transformers import AutoModel
# Use HuggingFace Hub version if available
```

### 3. Docker Alternative
```bash
docker pull xueerchen/scgpt:0.1.7
docker run --gpus all -v /data:/data xueerchen/scgpt:0.1.7
```

This comprehensive guide should help anyone facing similar scGPT embedding generation challenges in HPC environments or for foundation model comparison studies.