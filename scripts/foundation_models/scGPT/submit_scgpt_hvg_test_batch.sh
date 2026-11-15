#!/bin/bash
# ===================================================
# Process First 10 HVG-Filtered Platelet Chunks with scGPT - Test Batch
# Using HVG-optimized chunks for 5-10x processing speedup
# ===================================================

echo "Processing First 10 HVG-Filtered Platelet Chunks with scGPT - Test Batch"
echo "========================================================================"

BASE_DIR="/bigdata/godziklab/shared/Xinru/302006"
CHUNKS_DIR="$BASE_DIR/data/gpu_chunks_hvg"  # HVG-filtered chunks location
MODEL_DIR="$BASE_DIR/models/scGPT_human/scGPT_human"
JOBS_DIR="$BASE_DIR/jobs/scgpt_hvg_test"
LOGS_DIR="$BASE_DIR/logs/scgpt_hvg_test"

# Create directories
mkdir -p "$JOBS_DIR"
mkdir -p "$LOGS_DIR"

echo "HVG-filtered chunks directory: $CHUNKS_DIR"
echo "Model directory: $MODEL_DIR"

# Find first 10 HVG-filtered chunks
echo "Finding first 10 HVG-filtered chunks..."
hvg_chunks=($(ls -1 $CHUNKS_DIR/*_hvg.h5ad 2>/dev/null | head -10))

if [ ${#hvg_chunks[@]} -eq 0 ]; then
    echo "ERROR: No HVG-filtered chunks found in $CHUNKS_DIR"
    echo "Looking for files with pattern: *_hvg.h5ad"
    ls -la "$CHUNKS_DIR" | head -20
    exit 1
fi

echo "Found ${#hvg_chunks[@]} HVG-filtered chunks to process:"
for chunk in "${hvg_chunks[@]}"; do
    echo "  $(basename $chunk)"
done

# Cancel any existing test jobs
echo "Canceling existing scGPT HVG test jobs..."
scancel --name=scgpt_hvg_test_ 2>/dev/null || true

job_ids=()

# Process each HVG chunk
for i in "${!hvg_chunks[@]}"; do
    chunk_path="${hvg_chunks[$i]}"
    chunk_name=$(basename "$chunk_path" .h5ad)
    output_name="${chunk_name}_scgpt.h5ad"
    output_path="$CHUNKS_DIR/$output_name"
    
    # Skip if already processed
    if [ -f "$output_path" ]; then
        echo "Already processed: $output_name"
        continue
    fi
    
    job_num=$((i + 1))
    job_name="scgpt_hvg_test_${job_num}"
    job_script="$JOBS_DIR/${job_name}.sh"
    
    echo "Creating job $job_num/10: $(basename $chunk_path)"
    
    cat > "$job_script" << JOBEOF
#!/bin/bash
#SBATCH --job-name=${job_name}
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=512GB
#SBATCH --time=24:00:00
#SBATCH --output=${LOGS_DIR}/${job_name}_%j.out
#SBATCH --error=${LOGS_DIR}/${job_name}_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xinru.qiu@ucr.edu

echo "HVG-Filtered scGPT Processing - Test Batch Job ${job_num}/10"
echo "========================================================"
echo "Job ID: \$SLURM_JOB_ID"
echo "Node: \$SLURM_NODELIST"
echo "Memory: 512GB"
echo "Time: 24 hours"
echo "Start: \$(date)"
echo "Input: $(basename $chunk_path)"
echo "Output: $output_name"

cd $BASE_DIR

# Load environment
echo "Loading environment..."
module load singularity 2>/dev/null || echo "Using system singularity"
source ~/.bashrc 2>/dev/null || true
mamba activate scgpt || conda activate scgpt || echo "Using system python"

echo "Environment verification:"
echo "Python: \$(which python)"
echo "Conda env: \$CONDA_DEFAULT_ENV"
echo "Singularity: \$(which singularity)"

# Check input file
if [ ! -f "$chunk_path" ]; then
    echo "ERROR: Input file not found: $chunk_path"
    exit 1
fi

# Check if output exists
if [ -f "$output_path" ]; then
    echo "Output already exists: $output_path"
    exit 0
fi

# Process HVG chunk with scGPT
echo "Processing HVG-filtered chunk with scGPT..."

singularity exec scgpt_0.1.7.sif python << 'PYTHON_EOF'
import sys
sys.path.append('$BASE_DIR/302006_updategit/scripts/scGPT')

import scgpt as scg
import scanpy as sc
import numpy as np
import torch
import json
import time
import gc
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

from scgpt.model import TransformerModel
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.utils import set_seed

# File paths
input_file = "$chunk_path"
output_file = "$output_path"
model_dir = "$MODEL_DIR"

print(f"Processing HVG-filtered chunk ${job_num}/10")
print(f"Input: {input_file}")
print(f"Output: {output_file}")
print(f"Model: {model_dir}")

def process_hvg_chunk():
    """Process HVG-filtered platelet chunk with scGPT"""
    
    start_time = time.time()
    
    # Load HVG-filtered chunk
    print("Loading HVG-filtered chunk...")
    adata = sc.read_h5ad(input_file)
    print(f"Loaded HVG chunk: {adata.shape}")
    print(f"HVG genes: {adata.n_vars:,} (down from ~38,000)")
    
    # Load scGPT model components
    model_dir_path = Path(model_dir)
    vocab_file = model_dir_path / "vocab.json"
    model_file = model_dir_path / "best_model.pt"
    args_file = model_dir_path / "args.json"
    
    print(f"Loading vocabulary from {vocab_file}")
    vocab = GeneVocab.from_file(vocab_file)
    
    # Add special tokens
    special_tokens = ["<pad>", "<cls>", "<eoc>"]
    for s in special_tokens:
        if s not in vocab:
            vocab.append_token(s)
    vocab.set_default_index(vocab["<pad>"])
    
    print(f"Vocabulary size: {len(vocab):,}")
    
    # Match genes to vocabulary
    gene_names = adata.var.index.tolist()
    gene_ids = []
    valid_indices = []
    
    print("Matching HVG genes to scGPT vocabulary...")
    for idx, gene in enumerate(gene_names):
        if gene in vocab:
            gene_ids.append(vocab[gene])
            valid_indices.append(idx)
    
    vocab_match_rate = len(gene_ids) / len(gene_names) * 100
    print(f"HVG genes matched to vocab: {len(gene_ids):,}/{len(gene_names):,} ({vocab_match_rate:.1f}%)")
    
    # Filter to vocabulary genes
    adata_vocab = adata[:, valid_indices].copy()
    print(f"Using {adata_vocab.n_vars:,} vocabulary-matched HVG genes")
    
    # Get expression matrix
    X = adata_vocab.X
    if hasattr(X, 'toarray'):
        X = X.toarray()
    
    print(f"Expression matrix: {X.shape}")
    print(f"Expression range: [{np.min(X):.3f}, {np.max(X):.3f}]")
    
    # Binning (optimized for HVG genes)
    print(f"Binning {X.shape[1]:,} HVG genes...")
    X_binned = np.zeros_like(X, dtype=int)
    
    for i in range(X.shape[1]):
        if i % 200 == 0 and i > 0:
            print(f"   Gene {i:,}/{X.shape[1]:,}")
        
        gene_expr = X[:, i]
        if np.max(gene_expr) > 0:
            non_zero = gene_expr[gene_expr > 0]
            if len(non_zero) > 50:
                # Use quantile-based binning
                bins = np.quantile(non_zero, np.linspace(0, 1, 51))
                bins = np.unique(bins)
                if len(bins) > 1:
                    X_binned[:, i] = np.digitize(gene_expr, bins)
                else:
                    X_binned[:, i] = (gene_expr > 0).astype(int)
            else:
                # Simple thresholding for low-expression genes
                max_val = np.max(gene_expr)
                if max_val > 0:
                    X_binned[:, i] = np.clip((gene_expr / max_val * 50).astype(int), 0, 50)
    
    print(f"Binning complete. Range: [{np.min(X_binned)}, {np.max(X_binned)}]")
    
    # Load scGPT model
    print("Loading scGPT model...")
    with open(args_file, "r") as f:
        args = json.load(f)
    
    model = TransformerModel(
        len(vocab),
        args.get("embsize", 512),
        args.get("nheads", 8),
        args.get("d_hid", 512),
        args.get("nlayers", 12),
        vocab=vocab,
        dropout=0.0,
        pad_token="<pad>",
        pad_value=0
    )
    
    # Load model weights
    state_dict = torch.load(model_file, map_location='cpu')
    model.load_state_dict(state_dict, strict=False)
    model.eval()
    
    print("Model loaded successfully")
    
    # Process cells efficiently
    n_cells = X_binned.shape[0]
    embeddings = []
    
    print(f"Processing {n_cells:,} cells with HVG-optimized pipeline...")
    
    with torch.no_grad():
        for i in range(n_cells):
            if i % 100 == 0:
                elapsed = time.time() - start_time
                if i > 0:
                    rate = i / elapsed
                    remaining = n_cells - i
                    est_remaining = remaining / rate
                    print(f"   Cell {i+1:,}/{n_cells:,} | {elapsed:.1f}s | {rate:.1f}/s | Est: {est_remaining:.1f}s")
                else:
                    print(f"   Starting cell processing...")
            
            try:
                # Prepare tensors
                gene_tensor = torch.tensor(gene_ids, dtype=torch.long).unsqueeze(0)
                value_tensor = torch.tensor(X_binned[i:i+1], dtype=torch.float)
                padding_mask = torch.zeros_like(value_tensor, dtype=torch.bool)
                
                # Generate embedding
                output = model(gene_tensor, value_tensor, src_key_padding_mask=padding_mask)
                
                # Extract cell embedding
                if hasattr(output, 'cell_emb') and output.cell_emb is not None:
                    emb = output.cell_emb.cpu().numpy()
                elif isinstance(output, dict) and 'cell_emb' in output:
                    emb = output['cell_emb'].cpu().numpy()
                else:
                    # Fallback: use mean of sequence embeddings
                    emb = output.mean(dim=1).cpu().numpy()
                
                embeddings.append(emb)
                
                # Clean up tensors
                del gene_tensor, value_tensor, padding_mask, output
                
            except Exception as e:
                print(f"   CRITICAL ERROR: Cell {i} failed: {e}")
                raise RuntimeError(
                    f"Cell {i} embedding generation failed. "
                    f"Cannot use synthetic data. Original error: {e}"
                )
            
            # Periodic cleanup
            if i % 50 == 0:
                gc.collect()
    
    # Combine embeddings
    final_embeddings = np.vstack(embeddings)
    processing_time = time.time() - start_time
    
    print(f"Generated embeddings: {final_embeddings.shape} in {processing_time:.2f}s")
    print(f"Processing speed: {n_cells/processing_time:.1f} cells/second")
    print(f"Embedding stats: mean={np.mean(final_embeddings):.3f}, std={np.std(final_embeddings):.3f}")
    
    # Create result AnnData
    result_adata = sc.AnnData(X=final_embeddings, obs=adata.obs.copy())
    result_adata.uns['scgpt_processing'] = {
        'model': 'scGPT_human',
        'embedding_dim': final_embeddings.shape[1],
        'hvg_genes_input': adata.n_vars,
        'vocab_genes_used': len(gene_ids),
        'vocab_match_rate': vocab_match_rate,
        'processing_time_seconds': processing_time,
        'cells_per_second': n_cells / processing_time,
        'hvg_optimized': True,
        'test_batch': True,
        'chunk_number': '${job_num}/10'
    }
    
    # Save result
    result_adata.write_h5ad(output_file)
    print(f"Saved HVG-processed result: {output_file}")
    
    # Cleanup
    del model, adata, adata_vocab, embeddings, result_adata
    gc.collect()
    
    return True

# Execute processing
try:
    success = process_hvg_chunk()
    
    if success:
        print(f"SUCCESS: HVG chunk ${job_num}/10 processed")
        
        # Verify result
        result = sc.read_h5ad(output_file)
        print(f"Final result: {result.shape}")
        
        if 'scgpt_processing' in result.uns:
            info = result.uns['scgpt_processing']
            print(f"Processing info:")
            print(f"  - HVG genes: {info.get('hvg_genes_input', 'N/A'):,}")
            print(f"  - Vocab genes: {info.get('vocab_genes_used', 'N/A'):,}")
            print(f"  - Match rate: {info.get('vocab_match_rate', 'N/A'):.1f}%")
            print(f"  - Speed: {info.get('cells_per_second', 'N/A'):.1f} cells/sec")
        
        print("HVG-optimized scGPT processing completed successfully!")
        
    else:
        print(f"FAILED: HVG chunk ${job_num}/10 processing failed")
        sys.exit(1)
        
except Exception as e:
    print(f"ERROR: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

PYTHON_EOF

echo "Completed: \$(date)"
echo "HVG-optimized scGPT chunk ${job_num}/10 finished"
JOBEOF

    chmod +x "$job_script"
    
    # Submit job
    job_id=$(sbatch --parsable "$job_script")
    job_ids+=($job_id)
    
    echo "  Submitted job $job_num/10: Job ID $job_id"
    sleep 0.3
done

echo
echo "HVG-Filtered scGPT Test Batch Summary:"
echo "======================================"
echo "Jobs submitted: ${#job_ids[@]}/10"
echo "Job IDs: ${job_ids[*]}"
echo
echo "Processing specifications:"
echo "• HVG genes: ~3,000 (down from ~38,000)"
echo "• Expected speedup: 5-10x"
echo "• Memory: 200GB per job"
echo "• Time: 4 hours per job"
echo "• Method: scGPT embeddings on HVG-filtered data"
echo
echo "Monitor progress:"
echo "squeue -u \$USER"
echo "watch -n 10 'echo \"Test jobs: \$(squeue -u \$USER --name=scgpt_hvg_test_ | wc -l), Completed: \$(ls -1 $CHUNKS_DIR/*_hvg_scgpt.h5ad 2>/dev/null | wc -l)/10\"'"
echo
echo "Expected benefits from HVG optimization:"
echo "• 5-10x faster processing time"
echo "• Reduced memory requirements"
echo "• Maintained biological signal"
echo "• Scalable to full 188 chunks"
echo
echo "Test batch submitted successfully!"