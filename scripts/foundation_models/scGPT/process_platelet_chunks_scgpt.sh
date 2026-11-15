#!/bin/bash
# ===================================================
# Process Platelet Chunks with scGPT Embeddings - FIXED VARIABLES
# COVID platelets (9 chunks) + HEK293T platelets (16 chunks) + HCT116 platelets (16 chunks)
# ===================================================

echo "🩸 PROCESSING PLATELET CHUNKS WITH scGPT - FIXED VERSION"
echo "======================================================="

BASE_DIR="/bigdata/godziklab/shared/Xinru/302006"
CHUNKS_DIR="$BASE_DIR/data/gpu_chunks"
MODEL_DIR="$BASE_DIR/models/scGPT_human/scGPT_human"

# Count existing processed chunks
covid_processed=$(ls -1 $CHUNKS_DIR/covid_platelets_*_scgpt.h5ad 2>/dev/null | wc -l)
hek293t_processed=$(ls -1 $CHUNKS_DIR/hek293t_platelets_*_scgpt.h5ad 2>/dev/null | wc -l)
hct116_processed=$(ls -1 $CHUNKS_DIR/hct116_platelets_*_scgpt.h5ad 2>/dev/null | wc -l)

echo "📊 Current platelet chunk status:"
echo "COVID platelets: $covid_processed/9 processed"
echo "HEK293T platelets: $hek293t_processed/16 processed"
echo "HCT116 platelets: $hct116_processed/16 processed"
echo "Total: $((covid_processed + hek293t_processed + hct116_processed))/41 processed"

# Cancel any existing jobs
echo "🛑 Canceling existing platelet jobs..."
scancel --name=covid_platelets_ 2>/dev/null || true
scancel --name=hek293t_platelets_ 2>/dev/null || true
scancel --name=hct116_platelets_ 2>/dev/null || true

# Create jobs directory
mkdir -p "$BASE_DIR/jobs/platelet_chunks_fixed"
mkdir -p "$BASE_DIR/logs"

echo -e "\n🔧 Creating FIXED scGPT jobs for platelet chunks..."

job_ids=()

# Function to process chunks for a dataset
process_dataset_chunks() {
    local dataset_prefix=$1
    local total_chunks=$2
    local dataset_name=$3
    
    echo -e "\n📊 Processing $dataset_name chunks..."
    
    for chunk_num in $(seq 1 $total_chunks); do
        chunk_file="${dataset_prefix}_chunk_$(printf '%03d' $chunk_num)_of_$(printf '%03d' $total_chunks).h5ad"
        output_file="${dataset_prefix}_chunk_$(printf '%03d' $chunk_num)_of_$(printf '%03d' $total_chunks)_scgpt.h5ad"
        
        # Check if input chunk exists
        if [ ! -f "$CHUNKS_DIR/$chunk_file" ]; then
            echo "  ⚠️ Missing input chunk: $chunk_file"
            continue
        fi
        
        # Check if already processed
        if [ -f "$CHUNKS_DIR/$output_file" ]; then
            echo "  ✅ Already processed: $output_file"
            continue
        fi
        
        # Create job script with FIXED variable expansion
        job_script="$BASE_DIR/jobs/platelet_chunks_fixed/${dataset_prefix}_${chunk_num}_fixed.sh"
        
        cat > "$job_script" << JOBEOF
#!/bin/bash
#SBATCH --job-name=${dataset_prefix}_${chunk_num}_fix
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=400GB
#SBATCH --time=48:00:00
#SBATCH --output=$BASE_DIR/logs/${dataset_prefix}_${chunk_num}_fix_%j.out
#SBATCH --error=$BASE_DIR/logs/${dataset_prefix}_${chunk_num}_fix_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xinru.qiu@ucr.edu

echo "🩸 $dataset_name chunk processing: ${chunk_num}/${total_chunks} - FIXED"
echo "Job ID: \$SLURM_JOB_ID"
echo "Node: \$SLURM_NODELIST"
echo "Memory: 400GB"
echo "Time: 48 hours"
echo "Start: \$(date)"

cd $BASE_DIR

# Define file paths clearly
INPUT_FILE="$CHUNKS_DIR/$chunk_file"
OUTPUT_FILE="$CHUNKS_DIR/$output_file"

echo "📊 Input: \$INPUT_FILE"
echo "📤 Output: \$OUTPUT_FILE"

# Check files exist
if [ ! -f "\$INPUT_FILE" ]; then
    echo "❌ Input file not found: \$INPUT_FILE"
    exit 1
fi

if [ -f "\$OUTPUT_FILE" ]; then
    echo "✅ Output already exists: \$OUTPUT_FILE"
    exit 0
fi

# Load environment
echo "🔧 Loading environment..."
module load singularity 2>/dev/null || echo "Using system singularity"
source ~/.bashrc
mamba activate scgpt || conda activate scgpt || echo "Using system python"

echo "📦 Environment verified:"
echo "Python: \$(which python)"
echo "Conda env: \$CONDA_DEFAULT_ENV"
echo "Singularity: \$(which singularity)"

# Process chunk with scGPT - FIXED PYTHON CODE
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

# FIXED: Use proper file paths
input_file = "$CHUNKS_DIR/$chunk_file"
output_file = "$CHUNKS_DIR/$output_file"
model_dir = "$MODEL_DIR"

print(f"🩸 Processing $dataset_name chunk ${chunk_num}/${total_chunks}")
print(f"📊 Input: {input_file}")
print(f"📤 Output: {output_file}")
print(f"🏠 Model: {model_dir}")

def process_platelet_chunk_fixed():
    """Process platelet chunk with FIXED file paths"""
    
    start_time = time.time()
    
    # Load chunk
    print(f"📊 Loading chunk...")
    adata = sc.read_h5ad(input_file)
    print(f"✅ Loaded chunk ${chunk_num}: {adata.shape}")
    
    # Load scGPT components
    model_dir_path = Path(model_dir)
    vocab_file = model_dir_path / "vocab.json"
    model_file = model_dir_path / "best_model.pt"
    args_file = model_dir_path / "args.json"
    
    print(f"📚 Loading vocabulary from {vocab_file}")
    vocab = GeneVocab.from_file(vocab_file)
    special_tokens = ["<pad>", "<cls>", "<eoc>"]
    for s in special_tokens:
        if s not in vocab:
            vocab.append_token(s)
    vocab.set_default_index(vocab["<pad>"])
    
    # Process genes
    gene_names = adata.var.index.tolist()
    gene_ids = []
    valid_indices = []
    
    for idx, gene in enumerate(gene_names):
        if gene in vocab:
            gene_ids.append(vocab[gene])
            valid_indices.append(idx)
    
    print(f"🧬 Using {len(gene_ids):,} genes from vocabulary")
    
    # Filter to vocab genes
    adata_vocab = adata[:, valid_indices].copy()
    
    # Get expression data
    X = adata_vocab.X
    if hasattr(X, 'toarray'):
        X = X.toarray()
    
    print(f"📊 Expression matrix: {X.shape}")
    
    # Binning
    print(f"🔢 Binning {X.shape[1]:,} genes...")
    X_binned = np.zeros_like(X, dtype=int)
    
    for i in range(X.shape[1]):
        if i % 1000 == 0 and i > 0:
            print(f"   Gene {i:,}/{X.shape[1]:,}")
        
        gene_expr = X[:, i]
        if np.max(gene_expr) > 0:
            non_zero = gene_expr[gene_expr > 0]
            if len(non_zero) > 50:
                bins = np.quantile(non_zero, np.linspace(0, 1, 51))
                bins = np.unique(bins)
                X_binned[:, i] = np.digitize(gene_expr, bins)
            else:
                max_val = np.max(gene_expr)
                X_binned[:, i] = np.clip((gene_expr / max_val * 50).astype(int), 0, 50)
    
    # Load model
    print(f"🔧 Loading scGPT model from {model_file}")
    with open(args_file, "r") as f:
        args = json.load(f)
    
    model = TransformerModel(
        len(vocab), args.get("embsize", 512), args.get("nheads", 8),
        args.get("d_hid", 512), args.get("nlayers", 12), vocab=vocab,
        dropout=0.0, pad_token="<pad>", pad_value=0
    )
    
    state_dict = torch.load(model_file, map_location='cpu')
    model.load_state_dict(state_dict, strict=False)
    model.eval()
    
    # Process cells
    n_cells = X_binned.shape[0]
    embeddings = []
    
    print(f"🔬 Processing {n_cells:,} cells...")
    
    with torch.no_grad():
        for i in range(n_cells):
            if i % 500 == 0:
                elapsed = (time.time() - start_time) / 3600
                if i > 0:
                    rate = i / elapsed
                    remaining = n_cells - i
                    est_remaining = remaining / rate / 3600
                    print(f"   Cell {i+1:,}/{n_cells:,} | {elapsed:.2f}h | {rate:.0f}/h | Est: {est_remaining:.2f}h")
                else:
                    print(f"   Cell {i+1:,}/{n_cells:,} | Starting...")
            
            try:
                gene_tensor = torch.tensor(gene_ids, dtype=torch.long).unsqueeze(0)
                value_tensor = torch.tensor(X_binned[i:i+1], dtype=torch.float)
                padding_mask = torch.zeros_like(value_tensor, dtype=torch.bool)
                
                output = model(gene_tensor, value_tensor, src_key_padding_mask=padding_mask)
                
                # Extract embedding
                if hasattr(output, 'cell_emb') and output.cell_emb is not None:
                    emb = output.cell_emb.cpu().numpy()
                elif isinstance(output, dict) and 'cell_emb' in output:
                    emb = output['cell_emb'].cpu().numpy()
                else:
                    emb = output.mean(dim=1).cpu().numpy()
                
                embeddings.append(emb)
                
                # Clean up
                del gene_tensor, value_tensor, padding_mask, output
                
            except Exception as e:
                print(f"   CRITICAL ERROR: Cell {i} failed: {e}")
                raise RuntimeError(
                    f"Cell {i} embedding generation failed. "
                    f"Cannot use synthetic data. Original error: {e}"
                )
            
            # Periodic cleanup
            if i % 200 == 0:
                gc.collect()
    
    # Combine embeddings
    final_embeddings = np.vstack(embeddings)
    elapsed_hours = (time.time() - start_time) / 3600
    
    print(f"✅ Generated embeddings: {final_embeddings.shape} in {elapsed_hours:.2f}h")
    
    # Create result
    result_adata = sc.AnnData(X=final_embeddings, obs=adata.obs)
    result_adata.uns['scgpt_processing'] = {
        'model': 'scGPT_human',
        'embedding_dim': final_embeddings.shape[1],
        'genes_used': len(gene_ids),
        'processing_time_hours': elapsed_hours,
        'chunk': '${chunk_num}/${total_chunks}',
        'dataset': '$dataset_name'
    }
    
    result_adata.write_h5ad(output_file)
    print(f"💾 Saved: {output_file}")
    
    # Cleanup
    del model, adata, embeddings, result_adata
    gc.collect()
    
    return True

# Process this chunk
try:
    success = process_platelet_chunk_fixed()
    
    if success:
        print(f"✅ $dataset_name CHUNK ${chunk_num} SUCCESS - FIXED")
        result = sc.read_h5ad(output_file)
        print(f"📊 Final result: {result.shape}")
        
        # Verify embedding quality
        if hasattr(result.X, 'shape') and len(result.X.shape) == 2:
            print(f"🎯 Embedding dimensions: {result.X.shape[1]}")
            print(f"🎯 Embedding range: [{np.min(result.X):.3f}, {np.max(result.X):.3f}]")
        
    else:
        print(f"❌ $dataset_name CHUNK ${chunk_num} FAILED")
        import sys
        sys.exit(1)
        
except Exception as e:
    print(f"❌ ERROR: {e}")
    import traceback
    traceback.print_exc()
    import sys
    sys.exit(1)

PYTHON_EOF

echo "Completed: \$(date)"
echo "🎉 $dataset_name chunk ${chunk_num} finished - FIXED VERSION"
JOBEOF

        chmod +x "$job_script"
        
        # Submit the job
        job_id=$(sbatch --parsable "$job_script")
        job_ids+=($job_id)
        
        echo "  📤 Submitted FIXED $dataset_name chunk $chunk_num: Job ID $job_id"
        sleep 0.2
    done
}

# Process all platelet datasets
process_dataset_chunks "covid_platelets" 9 "COVID Platelets"
process_dataset_chunks "hek293t_platelets" 16 "HEK293T Platelets"
process_dataset_chunks "hct116_platelets" 16 "HCT116 Platelets"

echo -e "\n📊 FIXED Platelet Chunks Job Summary:"
echo "Total jobs submitted: ${#job_ids[@]}"
echo "Job IDs: ${job_ids[*]}"

echo -e "\n🔧 FIXED Processing specifications:"
echo "• Memory: 400GB per job"
echo "• Time: 48 hours per job"
echo "• Method: Single-cell scGPT embedding"
echo "• Model: scGPT_human"
echo "• Variables: FIXED file path expansion"

echo -e "\n🔍 Monitor FIXED jobs:"
echo "squeue -u \$USER"
echo "watch -n 30 'echo \"FIXED jobs running: \$(squeue -u \$USER | wc -l)\"'"

echo -e "\n✅ FIXED platelet chunk processing jobs submitted successfully!"

BASE_DIR="/bigdata/godziklab/shared/Xinru/302006"
CHUNKS_DIR="$BASE_DIR/data/gpu_chunks"
MODEL_DIR="$BASE_DIR/models/scGPT_human/scGPT_human"

# Count existing processed chunks
covid_processed=$(ls -1 $CHUNKS_DIR/covid_platelets_*_scgpt.h5ad 2>/dev/null | wc -l)
hek293t_processed=$(ls -1 $CHUNKS_DIR/hek293t_platelets_*_scgpt.h5ad 2>/dev/null | wc -l)
hct116_processed=$(ls -1 $CHUNKS_DIR/hct116_platelets_*_scgpt.h5ad 2>/dev/null | wc -l)

echo "📊 Current platelet chunk status:"
echo "COVID platelets: $covid_processed/9 processed"
echo "HEK293T platelets: $hek293t_processed/16 processed"
echo "HCT116 platelets: $hct116_processed/16 processed"
echo "Total: $((covid_processed + hek293t_processed + hct116_processed))/41 processed"

# Create jobs directory
mkdir -p "$BASE_DIR/jobs/platelet_chunks"
mkdir -p "$BASE_DIR/logs"

echo -e "\n🔧 Creating scGPT jobs for unprocessed platelet chunks..."

job_ids=()

# Function to process chunks for a dataset
process_dataset_chunks() {
    local dataset_prefix=$1
    local total_chunks=$2
    local dataset_name=$3
    
    echo -e "\n📊 Processing $dataset_name chunks..."
    
    for chunk_num in $(seq 1 $total_chunks); do
        chunk_file="${dataset_prefix}_chunk_$(printf '%03d' $chunk_num)_of_$(printf '%03d' $total_chunks).h5ad"
        output_file="${dataset_prefix}_chunk_$(printf '%03d' $chunk_num)_of_$(printf '%03d' $total_chunks)_scgpt.h5ad"
        
        # Check if input chunk exists
        if [ ! -f "$CHUNKS_DIR/$chunk_file" ]; then
            echo "  ⚠️ Missing input chunk: $chunk_file"
            continue
        fi
        
        # Check if already processed
        if [ -f "$CHUNKS_DIR/$output_file" ]; then
            echo "  ✅ Already processed: $output_file"
            continue
        fi
        
        # Create job script
        job_script="$BASE_DIR/jobs/platelet_chunks/${dataset_prefix}_${chunk_num}.sh"
        
        cat > "$job_script" << JOBEOF
#!/bin/bash
#SBATCH --job-name=${dataset_prefix}_${chunk_num}
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=400GB
#SBATCH --time=48:00:00
#SBATCH --output=$BASE_DIR/logs/${dataset_prefix}_${chunk_num}_%j.out
#SBATCH --error=$BASE_DIR/logs/${dataset_prefix}_${chunk_num}_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xinru.qiu@ucr.edu

echo "🩸 $dataset_name chunk processing: ${chunk_num}/${total_chunks}"
echo "Job ID: \$SLURM_JOB_ID"
echo "Node: \$SLURM_NODELIST"
echo "Memory: 400GB"
echo "Time: 48 hours"
echo "Start: \$(date)"

cd $BASE_DIR

# Load environment
echo "🔧 Loading environment..."
module load singularity 2>/dev/null || echo "Using system singularity"
source ~/.bashrc
mamba activate scgpt || conda activate scgpt || echo "Using system python"

# Check if output already exists
if [ -f "$CHUNKS_DIR/$output_file" ]; then
    echo "✅ Output already exists: $output_file"
    exit 0
fi

# Process chunk with scGPT
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

print(f"🩸 Processing $dataset_name chunk ${chunk_num}/${total_chunks}")
print(f"💾 Memory: 400GB")
print(f"🕐 Time: 48 hours")

def process_platelet_chunk(input_path, output_path, model_dir):
    """Process platelet chunk with scGPT embeddings"""
    
    start_time = time.time()
    
    # Load chunk
    adata = sc.read_h5ad(input_path)
    print(f"📊 $dataset_name chunk ${chunk_num}: {adata.shape}")
    
    # Load scGPT components
    model_dir = Path(model_dir)
    vocab_file = model_dir / "vocab.json"
    model_file = model_dir / "best_model.pt"
    args_file = model_dir / "args.json"
    
    # Load vocabulary
    vocab = GeneVocab.from_file(vocab_file)
    special_tokens = ["<pad>", "<cls>", "<eoc>"]
    for s in special_tokens:
        if s not in vocab:
            vocab.append_token(s)
    vocab.set_default_index(vocab["<pad>"])
    
    # Match gene processing approach with existing chunks
    gene_names = adata.var.index.tolist()
    gene_ids = []
    valid_indices = []
    
    for idx, gene in enumerate(gene_names):
        if gene in vocab:
            gene_ids.append(vocab[gene])
            valid_indices.append(idx)
    
    print(f"🧬 Using {len(gene_ids):,} genes from vocabulary")
    
    # Filter to vocab genes
    adata_vocab = adata[:, valid_indices].copy()
    
    # Get expression data
    X = adata_vocab.X
    if hasattr(X, 'toarray'):
        X = X.toarray()
    
    # Binning
    print(f"🔢 Binning {X.shape[1]:,} genes...")
    X_binned = np.zeros_like(X, dtype=int)
    
    for i in range(X.shape[1]):
        if i % 1000 == 0 and i > 0:
            print(f"   Gene {i:,}/{X.shape[1]:,}")
        
        gene_expr = X[:, i]
        if np.max(gene_expr) > 0:
            non_zero = gene_expr[gene_expr > 0]
            if len(non_zero) > 50:
                bins = np.quantile(non_zero, np.linspace(0, 1, 51))
                bins = np.unique(bins)
                X_binned[:, i] = np.digitize(gene_expr, bins)
            else:
                max_val = np.max(gene_expr)
                X_binned[:, i] = np.clip((gene_expr / max_val * 50).astype(int), 0, 50)
    
    # Load model
    print(f"🔧 Loading scGPT model...")
    with open(args_file, "r") as f:
        args = json.load(f)
    
    model = TransformerModel(
        len(vocab), args.get("embsize", 512), args.get("nheads", 8),
        args.get("d_hid", 512), args.get("nlayers", 12), vocab=vocab,
        dropout=0.0, pad_token="<pad>", pad_value=0
    )
    
    state_dict = torch.load(model_file, map_location='cpu')
    model.load_state_dict(state_dict, strict=False)
    model.eval()
    
    # Process cells (memory-efficient approach)
    n_cells = X_binned.shape[0]
    embeddings = []
    batch_size = 1  # Single cell processing for stability
    
    print(f"🔬 Processing {n_cells:,} cells...")
    
    with torch.no_grad():
        for i in range(n_cells):
            if i % 500 == 0:
                elapsed = (time.time() - start_time) / 3600
                remaining = n_cells - i
                rate = i / max(elapsed, 0.001) * 3600 if i > 0 else 0
                est_remaining = remaining / max(rate, 1) if rate > 0 else 0
                print(f"   Cell {i+1:,}/{n_cells:,} | {elapsed:.2f}h | Rate: {rate:.0f}/h | Est: {est_remaining:.2f}h")
            
            try:
                gene_tensor = torch.tensor(gene_ids, dtype=torch.long).unsqueeze(0)
                value_tensor = torch.tensor(X_binned[i:i+1], dtype=torch.float)
                padding_mask = torch.zeros_like(value_tensor, dtype=torch.bool)
                
                output = model(gene_tensor, value_tensor, src_key_padding_mask=padding_mask)
                
                # Extract embedding
                if hasattr(output, 'cell_emb') and output.cell_emb is not None:
                    emb = output.cell_emb.cpu().numpy()
                elif isinstance(output, dict) and 'cell_emb' in output:
                    emb = output['cell_emb'].cpu().numpy()
                else:
                    emb = output.mean(dim=1).cpu().numpy()
                
                embeddings.append(emb)
                
                # Clean up
                del gene_tensor, value_tensor, padding_mask, output
                
            except Exception as e:
                print(f"   CRITICAL ERROR: Cell {i} failed: {e}")
                raise RuntimeError(
                    f"Cell {i} embedding generation failed. "
                    f"Cannot use synthetic data. Original error: {e}"
                )
            
            # Periodic cleanup
            if i % 200 == 0:
                gc.collect()
                torch.cuda.empty_cache() if torch.cuda.is_available() else None
    
    # Combine embeddings
    final_embeddings = np.vstack(embeddings)
    elapsed_hours = (time.time() - start_time) / 3600
    
    print(f"✅ Generated embeddings: {final_embeddings.shape} in {elapsed_hours:.2f}h")
    
    # Create result
    result_adata = sc.AnnData(X=final_embeddings, obs=adata.obs)
    result_adata.uns['scgpt_processing'] = {
        'model': 'scGPT_human',
        'embedding_dim': final_embeddings.shape[1],
        'genes_used': len(gene_ids),
        'processing_time_hours': elapsed_hours,
        'chunk': '${chunk_num}/${total_chunks}',
        'dataset': '$dataset_name'
    }
    
    result_adata.write_h5ad(output_path)
    print(f"💾 Saved: {output_path}")
    
    # Cleanup
    del model, adata, embeddings, result_adata
    gc.collect()
    
    return True

# Process this chunk
input_path = Path('$CHUNKS_DIR/$chunk_file')
output_path = Path('$CHUNKS_DIR/$output_file')
model_dir = Path('$MODEL_DIR')

try:
    success = process_platelet_chunk(input_path, output_path, model_dir)
    
    if success:
        print(f"✅ $dataset_name CHUNK ${chunk_num} SUCCESS")
        result = sc.read_h5ad(output_path)
        print(f"📊 Final: {result.shape}")
    else:
        print(f"❌ $dataset_name CHUNK ${chunk_num} FAILED")
        sys.exit(1)
        
except Exception as e:
    print(f"❌ ERROR: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

PYTHON_EOF

echo "Completed: \$(date)"
echo "🎉 $dataset_name chunk ${chunk_num} finished"
JOBEOF

        chmod +x "$job_script"
        
        # Submit the job
        job_id=$(sbatch --parsable "$job_script")
        job_ids+=($job_id)
        
        echo "  📤 Submitted $dataset_name chunk $chunk_num: Job ID $job_id"
        sleep 0.3
    done
}

# Process all platelet datasets
process_dataset_chunks "covid_platelets" 9 "COVID Platelets"
process_dataset_chunks "hek293t_platelets" 16 "HEK293T Platelets"
process_dataset_chunks "hct116_platelets" 16 "HCT116 Platelets"

echo -e "\n📊 Platelet Chunks Job Summary:"
echo "Total jobs submitted: ${#job_ids[@]}"
echo "Job IDs: ${job_ids[*]}"

echo -e "\n🩸 Platelet datasets processing:"
echo "• COVID platelets: 9 chunks"
echo "• HEK293T platelets: 16 chunks"  
echo "• HCT116 platelets: 16 chunks"
echo "• Total: 41 chunks"

echo -e "\n🔧 Processing specifications:"
echo "• Memory: 400GB per job"
echo "• Time: 48 hours per job"
echo "• Method: Single-cell scGPT embedding"
echo "• Model: scGPT_human"

echo -e "\n🔍 Monitor progress:"
echo "squeue -u \$USER"
echo "watch -n 30 'echo \"Platelet jobs running: \$(squeue -u \$USER | wc -l), COVID: \$(ls -1 $CHUNKS_DIR/covid_platelets_*_scgpt.h5ad | wc -l)/9, HEK293T: \$(ls -1 $CHUNKS_DIR/hek293t_platelets_*_scgpt.h5ad | wc -l)/16, HCT116: \$(ls -1 $CHUNKS_DIR/hct116_platelets_*_scgpt.h5ad | wc -l)/16\"'"

echo -e "\n⏱️ Expected completion: within 48 hours"
echo "🎯 After completion, you'll have scGPT embeddings for all platelet datasets!"

echo -e "\n✅ Platelet chunk processing jobs submitted successfully!"