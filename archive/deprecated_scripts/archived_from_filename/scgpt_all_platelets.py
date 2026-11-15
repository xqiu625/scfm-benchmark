#!/usr/bin/env python3
"""
scGPT All Platelets Processing
Process all platelet datasets through scGPT to generate embeddings
Uses proven memory-efficient chunked processing approach
"""

import pandas as pd
import numpy as np
import scanpy as sc
from pathlib import Path
import subprocess
import os
import math

class scGPTPlateletsProcessor:
    def __init__(self, base_path="/bigdata/godziklab/shared/Xinru/302006"):
        self.base_path = Path(base_path)
        self.raw_path = self.base_path / "data" / "raw"  # Raw data location
        self.scgpt_path = self.base_path / "embeddings" / "scgpt_full"
        self.scgpt_path.mkdir(parents=True, exist_ok=True)
        self.chunks_path = self.base_path / "data" / "gpu_chunks"
        self.chunks_path.mkdir(parents=True, exist_ok=True)
        
        # Platelet datasets to process (CORRECTED PATHS)
        self.platelet_datasets = {
            'hek293t_platelets': {
                'file': 'hek293t_platelets.h5ad',  # Raw data file
                'expected_cells': 88434,
                'description': 'HEK293T platelets'
            },
            'hct116_platelets': {
                'file': 'hct116_platelets.h5ad',   # Raw data file
                'expected_cells': 89738,
                'description': 'HCT116 platelets'
            },
            'covid_platelets': {
                'file': 'human_platelet_covid_severity.h5ad',  # Raw data file
                'expected_cells': 47977,
                'description': 'COVID platelets with severity progression'
            }
        }
        
        # scGPT processing parameters (proven from HEK293T epithelial success)
        self.processing_params = {
            'chunk_size': 5832,  # Proven chunk size
            'max_genes': 3000,   # Most variable genes
            'memory_gb': 400,    # RAM allocation
            'time_limit': '12:00:00',  # 12 hours
            'partition': 'highmem'
        }
        
    def analyze_platelet_datasets(self):
        """Analyze all platelet datasets for processing requirements"""
        print("🔍 ANALYZING PLATELET DATASETS FOR scGPT PROCESSING")
        print("📁 Data source: /bigdata/godziklab/shared/Xinru/302006/data/raw/")
        print("=" * 60)
        
        analysis_results = {}
        
        for dataset_name, dataset_info in self.platelet_datasets.items():
            print(f"\n📊 Analyzing {dataset_name}:")
            
            filepath = self.raw_path / dataset_info['file']  # Use raw_path
            
            if not filepath.exists():
                print(f"   ❌ File not found: {filepath}")
                analysis_results[dataset_name] = {'status': 'file_not_found'}
                continue
            
            try:
                # Load and analyze dataset
                adata = sc.read_h5ad(filepath)
                
                actual_cells = adata.n_obs
                actual_genes = adata.n_vars
                
                # Calculate chunking requirements
                chunk_size = self.processing_params['chunk_size']
                n_chunks = math.ceil(actual_cells / chunk_size)
                
                print(f"   ✅ {dataset_info['description']}")
                print(f"      Cells: {actual_cells:,} (expected: {dataset_info['expected_cells']:,})")
                print(f"      Genes: {actual_genes:,}")
                print(f"      Required chunks: {n_chunks}")
                print(f"      Chunk size: ~{chunk_size:,} cells each")
                
                # Check if already processed
                output_file = self.scgpt_path / f"{dataset_name}_scgpt_full.h5ad"
                already_processed = output_file.exists()
                
                if already_processed:
                    try:
                        processed_adata = sc.read_h5ad(output_file)
                        print(f"      🎯 Already processed: {processed_adata.shape}")
                        if 'X_scgpt' in processed_adata.obsm:
                            print(f"      ✅ scGPT embeddings: {processed_adata.obsm['X_scgpt'].shape}")
                    except:
                        print(f"      ⚠️ Output file exists but corrupted")
                        already_processed = False
                
                analysis_results[dataset_name] = {
                    'status': 'ready',
                    'cells': actual_cells,
                    'genes': actual_genes,
                    'chunks_needed': n_chunks,
                    'already_processed': already_processed,
                    'file_path': filepath,
                    'output_path': output_file
                }
                
            except Exception as e:
                print(f"   ❌ Error analyzing {dataset_name}: {e}")
                analysis_results[dataset_name] = {'status': 'error', 'error': str(e)}
        
        return analysis_results
    
    def create_chunks(self, dataset_name, dataset_info):
        """Create chunks for a specific platelet dataset"""
        print(f"\n📦 Creating chunks for {dataset_name}...")
        
        filepath = dataset_info['file_path']
        adata = sc.read_h5ad(filepath)
        
        chunk_size = self.processing_params['chunk_size']
        n_chunks = dataset_info['chunks_needed']
        
        print(f"   📊 Dataset: {adata.shape}")
        print(f"   🔢 Creating {n_chunks} chunks of ~{chunk_size} cells each")
        
        chunk_files = []
        
        for chunk_idx in range(n_chunks):
            start_idx = chunk_idx * chunk_size
            end_idx = min((chunk_idx + 1) * chunk_size, adata.n_obs)
            
            # Extract chunk
            chunk_adata = adata[start_idx:end_idx].copy()
            
            # Save chunk
            chunk_filename = f"{dataset_name}_chunk_{chunk_idx+1:03d}_of_{n_chunks:03d}.h5ad"
            chunk_filepath = self.chunks_path / chunk_filename
            
            chunk_adata.write_h5ad(chunk_filepath)
            chunk_files.append(chunk_filepath)
            
            print(f"   ✅ Chunk {chunk_idx+1}/{n_chunks}: {chunk_adata.n_obs} cells -> {chunk_filename}")
        
        print(f"   🎯 Created {len(chunk_files)} chunks for {dataset_name}")
        return chunk_files
    
    def create_processing_script(self, dataset_name, chunk_files):
        """Create SLURM batch script for processing chunks with proper environment setup"""
        print(f"\n📝 Creating processing script for {dataset_name}...")
        
        script_content = f"""#!/bin/bash
#SBATCH --job-name=scgpt_{dataset_name}
#SBATCH --partition={self.processing_params['partition']}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem={self.processing_params['memory_gb']}G
#SBATCH --time={self.processing_params['time_limit']}
#SBATCH --array=1-{len(chunk_files)}
#SBATCH --output={self.base_path}/logs/scgpt_{dataset_name}_%A_%a.out
#SBATCH --error={self.base_path}/logs/scgpt_{dataset_name}_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xinru.qiu@ucr.edu

echo "🚀 Starting scGPT processing for {dataset_name}"
echo "📊 Job: $SLURM_JOB_ID, Array task: $SLURM_ARRAY_TASK_ID"
echo "🖥️ Node: $SLURMD_NODENAME"
echo "💾 Memory: {self.processing_params['memory_gb']}G"
echo "⏰ Time limit: {self.processing_params['time_limit']}"
echo "Start: $(date)"

# Navigate to working directory
cd {self.base_path}

# Create logs directory if it doesn't exist
mkdir -p {self.base_path}/logs

# Load required modules and environment
echo "🔧 Loading environment..."
module load singularity 2>/dev/null || echo "Using system singularity"
source ~/.bashrc

# Activate conda environment
echo "🐍 Activating conda environment..."
mamba activate scgpt || conda activate scgpt || echo "Using system python"

echo "📦 Environment verification:"
echo "Python: $(which python)"
echo "Conda env: $CONDA_DEFAULT_ENV"
echo "Singularity: $(which singularity)"

# Set up CUDA (force CPU processing for stability)
export CUDA_VISIBLE_DEVICES=""
echo "🖥️ CUDA devices: $CUDA_VISIBLE_DEVICES (CPU mode)"

# Get chunk file based on array task ID
chunk_files=({' '.join([str(f) for f in chunk_files])})
chunk_file="${{chunk_files[$SLURM_ARRAY_TASK_ID-1]}}"
chunk_basename=$(basename "$chunk_file" .h5ad)
output_file="{self.chunks_path}/${{chunk_basename}}_scgpt.h5ad"

echo "🔄 Processing chunk: $chunk_file"
echo "📤 Output will be: $output_file"

# Check if output already exists
if [ -f "$output_file" ]; then
    echo "✅ Output already exists: $output_file"
    echo "Skipping processing."
    exit 0
fi

# Verify input file exists
if [ ! -f "$chunk_file" ]; then
    echo "❌ Input chunk file not found: $chunk_file"
    exit 1
fi

# Run scGPT processing using singularity
echo "🚀 Starting scGPT processing with singularity..."

singularity exec {self.base_path}/scgpt_0.1.7.sif python << 'PYTHON_EOF'
import sys
sys.path.append('{self.base_path}/302006_updategit/scripts/scGPT')

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

def process_platelet_chunk(input_path, output_path, model_dir):
    # Process platelet chunk through scGPT with memory-efficient approach
    
    start_time = time.time()
    
    try:
        # Load chunk
        print(f"📊 Loading chunk: {{input_path}}")
        adata = sc.read_h5ad(input_path)
        print(f"✅ Loaded chunk: {{adata.shape}}")
        
        # Model paths
        model_dir = Path(model_dir)
        vocab_file = model_dir / "vocab.json"
        model_file = model_dir / "best_model.pt"
        args_file = model_dir / "args.json"
        
        print(f"🔧 Loading scGPT components...")
        print(f"   Vocab: {{vocab_file}}")
        print(f"   Model: {{model_file}}")
        print(f"   Args: {{args_file}}")
        
        # Load vocabulary
        vocab = GeneVocab.from_file(vocab_file)
        special_tokens = ["<pad>", "<cls>", "<eoc>"]
        for s in special_tokens:
            if s not in vocab:
                vocab.append_token(s)
        vocab.set_default_index(vocab["<pad>"])
        
        # Gene selection - use most variable genes (proven approach)
        max_genes = {self.processing_params['max_genes']}
        if adata.n_vars > max_genes:
            print(f"🧬 Selecting {{max_genes}} most variable genes from {{adata.n_vars}}")
            sc.pp.highly_variable_genes(adata, n_top_genes=max_genes, flavor='seurat_v3')
            adata = adata[:, adata.var.highly_variable].copy()
            print(f"✅ Reduced to {{adata.shape}}")
        
        # Map genes to vocabulary
        gene_names = adata.var.index.tolist()
        gene_ids = []
        valid_indices = []
        
        for idx, gene in enumerate(gene_names):
            if gene in vocab:
                gene_ids.append(vocab[gene])
                valid_indices.append(idx)
        
        print(f"🧬 Using {{len(gene_ids):,}} genes in vocabulary")
        
        if len(gene_ids) == 0:
            print("❌ No genes found in vocabulary!")
            return False
        
        # Filter to vocabulary genes
        adata_vocab = adata[:, valid_indices].copy()
        
        # Get expression data
        X = adata_vocab.X
        if hasattr(X, 'toarray'):
            X = X.toarray()
        
        print(f"📊 Expression matrix: {{X.shape}}")
        
        # Binning (memory efficient)
        print(f"🔢 Binning expression values...")
        X_binned = np.zeros_like(X, dtype=int)
        
        for i in range(X.shape[1]):
            if i % 500 == 0:
                print(f"   Gene {{i+1:,}}/{{X.shape[1]:,}}")
            
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
        
        print(f"✅ Model loaded successfully")
        
        # Process cells with single-cell batches (memory efficient)
        n_cells = X_binned.shape[0]
        embeddings = []
        
        print(f"🔬 Processing {{n_cells:,}} cells (single-cell batches)...")
        
        with torch.no_grad():
            for i in range(n_cells):
                if i % 100 == 0:
                    elapsed = (time.time() - start_time) / 3600
                    remaining = n_cells - i
                    est_remaining = (elapsed / max(i, 1)) * remaining
                    print(f"   Cell {{i+1:,}}/{{n_cells:,}} | {{elapsed:.2f}}h | Est: {{est_remaining:.2f}}h")
                
                cell_data = X_binned[i:i+1]
                
                try:
                    gene_tensor = torch.tensor(gene_ids, dtype=torch.long).unsqueeze(0)
                    value_tensor = torch.tensor(cell_data, dtype=torch.float)
                    padding_mask = torch.zeros_like(value_tensor, dtype=torch.bool)
                    
                    output = model(gene_tensor, value_tensor, src_key_padding_mask=padding_mask)
                    
                    # Extract embeddings
                    if hasattr(output, 'cell_emb') and output.cell_emb is not None:
                        emb = output.cell_emb.numpy()
                    elif isinstance(output, dict) and 'cell_emb' in output:
                        emb = output['cell_emb'].numpy()
                    else:
                        emb = output.mean(dim=1).numpy()
                    
                    embeddings.append(emb)
                    
                    # Clean up tensors
                    del gene_tensor, value_tensor, padding_mask, output
                    
                except Exception as e:
                    print(f"CRITICAL ERROR: Cell {i} failed: {e}")
                    raise RuntimeError(
                        f"Cell {i} embedding generation failed. "
                        f"Cannot use synthetic data. Original error: {e}"
                    )
                
                # Memory cleanup every 100 cells
                if i % 100 == 0:
                    gc.collect()
        
        # Combine embeddings
        final_embeddings = np.vstack(embeddings)
        elapsed_hours = (time.time() - start_time) / 3600
        
        print(f"✅ Generated embeddings: {{final_embeddings.shape}} in {{elapsed_hours:.2f}}h")
        
        # Save result
        result_adata = sc.AnnData(X=final_embeddings, obs=adata.obs.copy())
        result_adata.write_h5ad(output_path)
        
        print(f"💾 Saved to: {{output_path}}")
        
        # Cleanup
        del model, adata, embeddings, result_adata
        gc.collect()
        
        return True
        
    except Exception as e:
        print(f"❌ Error processing chunk: {{e}}")
        import traceback
        traceback.print_exc()
        return False

# Main processing
input_path = Path('$chunk_file')
output_path = Path('$output_file')
model_dir = Path('{self.base_path}/models/scGPT_human/scGPT_human')

print(f"🎯 Processing {dataset_name} chunk $SLURM_ARRAY_TASK_ID")
print(f"📊 Input: {{input_path}}")
print(f"📤 Output: {{output_path}}")
print(f"🔧 Model: {{model_dir}}")

try:
    success = process_platelet_chunk(input_path, output_path, model_dir)
    
    if success:
        print(f"✅ CHUNK SUCCESS")
        # Verify output
        result = sc.read_h5ad(output_path)
        print(f"📊 Final result: {{result.shape}}")
        
        if hasattr(result, 'X') and result.X is not None:
            print(f"📊 Embeddings shape: {{result.X.shape}}")
        else:
            print("⚠️ No embeddings matrix found")
            
    else:
        print(f"❌ CHUNK FAILED")
        sys.exit(1)
        
except Exception as e:
    print(f"❌ MAIN ERROR: {{e}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

PYTHON_EOF

echo "Completed: $(date)"
echo "🎉 {dataset_name} chunk $SLURM_ARRAY_TASK_ID finished"
"""
        
        script_path = self.base_path / f"process_{dataset_name}_scgpt.sh"
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        # Make executable
        script_path.chmod(0o755)
        
        print(f"   ✅ Created script: {script_path}")
        print(f"   📋 Features added:")
        print(f"      • Proper SLURM resource allocation")
        print(f"      • Conda environment activation")  
        print(f"      • Singularity container execution")
        print(f"      • Email notifications")
        print(f"      • Comprehensive logging")
        print(f"      • Memory-efficient processing")
        
        return script_path
    
    def submit_processing_jobs(self, dataset_name, script_path):
        """Submit SLURM processing jobs"""
        print(f"\n🚀 Submitting processing jobs for {dataset_name}...")
        
        try:
            # Submit job array
            result = subprocess.run(
                ['sbatch', str(script_path)],
                capture_output=True,
                text=True,
                cwd=str(self.base_path)
            )
            
            if result.returncode == 0:
                job_id = result.stdout.strip().split()[-1]
                print(f"   ✅ Submitted job array: {job_id}")
                print(f"   📊 Monitor with: squeue -u $USER | grep {dataset_name}")
                return job_id
            else:
                print(f"   ❌ Failed to submit: {result.stderr}")
                return None
                
        except Exception as e:
            print(f"   ❌ Error submitting job: {e}")
            return None
    
    def create_combination_script(self, dataset_name, dataset_info):
        """Create script to combine processed chunks"""
        print(f"\n🔗 Creating combination script for {dataset_name}...")
        
        n_chunks = dataset_info['chunks_needed']
        output_file = dataset_info['output_path']
        
        combination_script = f"""#!/usr/bin/env python3
'''
Combine {dataset_name} scGPT processed chunks
'''

import scanpy as sc
import pandas as pd
from pathlib import Path
import numpy as np

print('🔗 COMBINING {dataset_name.upper()} scGPT CHUNKS')
print('=' * 50)

chunks_path = Path('{self.chunks_path}')
output_path = Path('{output_file}')

# Find all processed chunks
chunk_files = []
for i in range(1, {n_chunks + 1}):
    chunk_file = chunks_path / f'{dataset_name}_chunk_{{i:03d}}_of_{n_chunks:03d}_scgpt.h5ad'
    if chunk_file.exists():
        chunk_files.append(chunk_file)
    else:
        print(f'❌ Missing chunk: {{chunk_file}}')

print(f'📊 Found {{len(chunk_files)}}/{n_chunks} processed chunks')

if len(chunk_files) != {n_chunks}:
    print('⚠️ Not all chunks processed yet')
    exit(1)

# Load and combine chunks
print('🔄 Loading chunks...')
chunk_list = []

for i, chunk_file in enumerate(chunk_files):
    print(f'   Loading chunk {{i+1}}/{{len(chunk_files)}}...')
    chunk_adata = sc.read_h5ad(chunk_file)
    
    if 'X_scgpt' not in chunk_adata.obsm:
        print(f'❌ Missing scGPT embeddings in chunk {{i+1}}')
        exit(1)
    
    chunk_list.append(chunk_adata)
    print(f'   ✅ Chunk {{i+1}}: {{chunk_adata.shape}}, embeddings: {{chunk_adata.obsm["X_scgpt"].shape}}')

# Concatenate chunks
print('🔗 Combining chunks...')
combined_adata = sc.concat(chunk_list, axis=0, join='outer')

print(f'✅ Combined dataset: {{combined_adata.shape}}')
print(f'✅ scGPT embeddings: {{combined_adata.obsm["X_scgpt"].shape}}')

# Save combined dataset
combined_adata.write_h5ad(output_path)
print(f'💾 Saved to: {{output_path}}')

print('🎉 {dataset_name} scGPT processing complete!')
"""
        
        script_path = self.base_path / f"combine_{dataset_name}_scgpt.py"
        with open(script_path, 'w') as f:
            f.write(combination_script)
        
        script_path.chmod(0o755)
        
        print(f"   ✅ Created combination script: {script_path}")
        return script_path
    
    def process_all_platelets(self):
        """Process all platelet datasets through scGPT"""
        print("🚀 scGPT ALL PLATELETS PROCESSING")
        print("🎯 Generating embeddings for HEK293T, HCT116, and COVID platelets")
        print("=" * 70)
        
        # Analyze datasets
        analysis_results = self.analyze_platelet_datasets()
        
        # Process each dataset
        job_ids = {}
        combination_scripts = {}
        
        for dataset_name, dataset_info in analysis_results.items():
            if dataset_info['status'] != 'ready':
                print(f"\n⚠️ Skipping {dataset_name}: {dataset_info['status']}")
                continue
            
            if dataset_info['already_processed']:
                print(f"\n✅ {dataset_name} already processed, skipping")
                continue
            
            print(f"\n🔄 PROCESSING {dataset_name.upper()}")
            print("-" * 50)
            
            # Create chunks
            chunk_files = self.create_chunks(dataset_name, dataset_info)
            
            # Create processing script
            script_path = self.create_processing_script(dataset_name, chunk_files)
            
            # Submit jobs
            job_id = self.submit_processing_jobs(dataset_name, script_path)
            if job_id:
                job_ids[dataset_name] = job_id
            
            # Create combination script
            combination_script = self.create_combination_script(dataset_name, dataset_info)
            combination_scripts[dataset_name] = combination_script
        
        # Summary
        print(f"\n" + "=" * 70)
        print("✅ scGPT PLATELET PROCESSING SETUP COMPLETE")
        print("=" * 70)
        
        print(f"\n📊 PROCESSING JOBS SUBMITTED:")
        for dataset_name, job_id in job_ids.items():
            print(f"   • {dataset_name}: Job {job_id}")
        
        print(f"\n🔗 COMBINATION SCRIPTS CREATED:")
        for dataset_name, script_path in combination_scripts.items():
            print(f"   • {dataset_name}: {script_path}")
        
        print(f"\n📋 MONITORING COMMANDS:")
        print(f"   • Check job status: squeue -u $USER")
        print(f"   • Check progress: ls {self.chunks_path}/*_scgpt.h5ad | wc -l")
        print(f"   • Combine chunks: python combine_<dataset>_scgpt.py")
        
        return job_ids, combination_scripts

if __name__ == "__main__":
    # Initialize processor
    processor = scGPTPlateletsProcessor()
    
    # Process all platelet datasets
    job_ids, scripts = processor.process_all_platelets()
    
    print(f"\n🎉 SETUP COMPLETE!")
    print("📋 Next steps:")
    print("   1. Monitor job progress with squeue")
    print("   2. Wait for all chunks to complete")
    print("   3. Run combination scripts to create final embeddings")
    print("   4. Proceed with cross-species transfer learning analysis")