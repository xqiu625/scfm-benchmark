#!/usr/bin/env python3
"""
Run RETNLA scGPT processing via Singularity container
Uses the same container approach that works for platelet processing
"""

import os
import subprocess
import sys
import logging
import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Configuration
BASE_DIR = "/bigdata/godziklab/shared/Xinru/302006"
DATA_DIR = f"{BASE_DIR}/data/raw"
OUTPUT_DIR = f"{BASE_DIR}/results/scgpt/retnla"
TEMP_DIR = f"{BASE_DIR}/temp/scgpt_retnla"
CONTAINER = f"{BASE_DIR}/scgpt_0.1.7.sif"
MODEL_DIR = f"{BASE_DIR}/models/scGPT_human/scGPT_human"

# Ortholog mapping
ORTHOLOG_FILE = f"{BASE_DIR}/data/external/genome_wide_orthologs/genome_wide_human_mouse_orthologs_complete_20250722_174716.tsv"

def prepare_retnla_data():
    """Prepare and combine RETNLA datasets with ortholog conversion"""
    
    logger.info("Preparing RETNLA datasets...")
    
    # Load ortholog mapping
    if not os.path.exists(ORTHOLOG_FILE):
        raise FileNotFoundError(f"Ortholog mapping file not found: {ORTHOLOG_FILE}")
    
    ortholog_df = pd.read_csv(ORTHOLOG_FILE, sep='\t')
    mouse_to_human = dict(zip(
        ortholog_df['mouse_gene_symbol'], 
        ortholog_df['human_gene_symbol']
    ))
    logger.info(f"Loaded {len(mouse_to_human)} mouse->human ortholog mappings")
    
    datasets = []
    
    # Load HCT116
    hct116_path = f"{DATA_DIR}/HCT116_retn_raw_counts.h5ad"
    if os.path.exists(hct116_path):
        hct116 = sc.read_h5ad(hct116_path)
        hct116.obs['dataset'] = 'hct116'
        hct116.obs['species'] = 'human'
        hct116.obs['gene_target'] = 'RETN'
        datasets.append(hct116)
        logger.info(f"HCT116: {hct116.shape}")
    
    # Load HEK293T
    hek293t_path = f"{DATA_DIR}/HEK293T_retn_raw_counts.h5ad"
    if os.path.exists(hek293t_path):
        hek293t = sc.read_h5ad(hek293t_path)
        hek293t.obs['dataset'] = 'hek293t'
        hek293t.obs['species'] = 'human'
        hek293t.obs['gene_target'] = 'RETN'
        datasets.append(hek293t)
        logger.info(f"HEK293T: {hek293t.shape}")
    
    # Load mouse with ortholog conversion
    mouse_path = f"{DATA_DIR}/mouse_epithelial_retnla_ko_wt.h5ad"
    if os.path.exists(mouse_path):
        mouse = sc.read_h5ad(mouse_path)
        
        # Convert mouse genes to human orthologs
        mouse_genes = mouse.var.index.tolist()
        logger.info(f"Mouse dataset has {len(mouse_genes)} genes")
        
        ortholog_mapping_dict = {}
        for mouse_gene in mouse_genes:
            if mouse_gene in mouse_to_human:
                human_ortholog = mouse_to_human[mouse_gene]
                if human_ortholog in ortholog_mapping_dict:
                    # Keep gene with higher expression
                    existing_mouse_gene = ortholog_mapping_dict[human_ortholog]
                    existing_expr = mouse[:, existing_mouse_gene].X.mean()
                    current_expr = mouse[:, mouse_gene].X.mean()
                    if current_expr > existing_expr:
                        ortholog_mapping_dict[human_ortholog] = mouse_gene
                else:
                    ortholog_mapping_dict[human_ortholog] = mouse_gene
        
        # Filter to ortholog genes
        mouse_genes_with_orthologs = list(ortholog_mapping_dict.values())
        ortholog_genes = list(ortholog_mapping_dict.keys())
        
        mouse_ortholog_mask = mouse.var.index.isin(mouse_genes_with_orthologs)
        mouse_filtered = mouse[:, mouse_ortholog_mask].copy()
        mouse_filtered = mouse_filtered[:, mouse_genes_with_orthologs]
        
        # Update gene names
        mouse_filtered.var.index = ortholog_genes
        mouse_filtered.obs['dataset'] = 'mouse'
        mouse_filtered.obs['species'] = 'mouse'
        mouse_filtered.obs['gene_target'] = 'RETNLA'
        
        datasets.append(mouse_filtered)
        logger.info(f"Mouse (converted): {mouse_filtered.shape}")
        logger.info(f"Found orthologs for {len(ortholog_genes)}/{len(mouse_genes)} mouse genes")
    
    # Combine datasets
    combined = sc.concat(datasets, join='outer', index_unique='-')
    logger.info(f"Combined RETNLA dataset: {combined.shape}")
    
    # Basic preprocessing
    sc.pp.filter_cells(combined, min_genes=200)
    sc.pp.filter_genes(combined, min_cells=3)
    logger.info(f"After filtering: {combined.shape}")
    
    return combined

def save_for_container(adata, temp_dir):
    """Save data for container processing"""
    
    os.makedirs(temp_dir, exist_ok=True)
    
    # Save combined dataset
    adata_path = f"{temp_dir}/retnla_combined.h5ad"
    adata.write_h5ad(adata_path)
    
    logger.info(f"Saved combined dataset: {adata_path}")
    logger.info(f"Dataset info: {adata.shape[0]} cells, {adata.shape[1]} genes")
    
    return adata_path

def create_container_script():
    """Create Python script to run inside container"""
    
    script_content = f'''
import sys
sys.path.append('{BASE_DIR}/302006_updategit/scripts/scGPT')

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

print("Processing RETNLA data with scGPT...")

def process_retnla_container():
    """Process RETNLA data inside container"""
    
    start_time = time.time()
    
    # Load data
    adata = sc.read_h5ad("retnla_combined.h5ad")
    print(f"Loaded RETNLA data: {{adata.shape}}")
    
    # Load scGPT components
    model_dir = Path("{MODEL_DIR}")
    vocab_file = model_dir / "vocab.json"
    model_file = model_dir / "best_model.pt"
    args_file = model_dir / "args.json"
    
    # Load vocabulary
    print("Loading vocabulary...")
    vocab = GeneVocab.from_file(vocab_file)
    special_tokens = ["<pad>", "<cls>", "<eoc>"]
    for s in special_tokens:
        if s not in vocab:
            vocab.append_token(s)
    vocab.set_default_index(vocab["<pad>"])
    
    # Match genes
    gene_names = adata.var.index.tolist()
    gene_ids = []
    valid_indices = []
    
    for idx, gene in enumerate(gene_names):
        if gene in vocab:
            gene_ids.append(vocab[gene])
            valid_indices.append(idx)
    
    print(f"Using {{len(gene_ids):,}} genes from vocabulary")
    
    # Filter to vocab genes
    adata_vocab = adata[:, valid_indices].copy()
    
    # Get expression data
    X = adata_vocab.X
    if hasattr(X, 'toarray'):
        X = X.toarray()
    
    # Binning
    print(f"Binning {{X.shape[1]:,}} genes...")
    X_binned = np.zeros_like(X, dtype=int)
    
    for i in range(X.shape[1]):
        if i % 1000 == 0 and i > 0:
            print(f"   Gene {{i:,}}/{{X.shape[1]:,}}")
        
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
    print("Loading scGPT model...")
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
    
    print(f"Processing {{n_cells:,}} cells...")
    
    with torch.no_grad():
        for i in range(n_cells):
            if i % 100 == 0:
                elapsed = (time.time() - start_time) / 3600
                if i > 0:
                    rate = i / elapsed * 3600
                    remaining = n_cells - i
                    est_remaining = remaining / rate if rate > 0 else 0
                    print(f"   Cell {{i+1:,}}/{{n_cells:,}} | {{elapsed:.2f}}h | {{rate:.0f}}/h | {{est_remaining:.2f}}h")
            
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
            if i % 50 == 0:
                gc.collect()
    
    # Combine embeddings
    final_embeddings = np.vstack(embeddings)
    elapsed_hours = (time.time() - start_time) / 3600
    
    print(f"Generated embeddings: {{final_embeddings.shape}} in {{elapsed_hours:.2f}}h")
    
    # Create result
    result_adata = sc.AnnData(X=final_embeddings, obs=adata.obs)
    result_adata.uns['scgpt_processing'] = {{
        'model': 'scGPT_human',
        'embedding_dim': final_embeddings.shape[1],
        'genes_used': len(gene_ids),
        'processing_time_hours': elapsed_hours,
        'datasets': ['hct116', 'hek293t', 'mouse'],
        'ortholog_conversion': True
    }}
    
    # Save result
    result_adata.write_h5ad("retnla_scgpt_result.h5ad")
    np.save("retnla_scgpt_embeddings.npy", final_embeddings)
    
    print("RETNLA scGPT processing completed successfully!")
    print(f"Result saved: {{final_embeddings.shape}}")
    
    return True

# Run processing
try:
    success = process_retnla_container()
    if success:
        print("SUCCESS: RETNLA scGPT processing completed")
    else:
        print("FAILED: RETNLA scGPT processing")
        sys.exit(1)
except Exception as e:
    print(f"ERROR: {{e}}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
'''
    
    script_path = f"{TEMP_DIR}/process_retnla_scgpt.py"
    with open(script_path, 'w') as f:
        f.write(script_content)
    
    return script_path

def run_container_processing(data_path, script_path):
    """Run processing via Singularity container"""
    
    logger.info("Running scGPT processing via Singularity container...")
    
    # Change to temp directory
    os.chdir(TEMP_DIR)
    
    # Run container
    cmd = [
        'singularity', 'exec',
        '--bind', f"{MODEL_DIR}:{MODEL_DIR}",
        CONTAINER,
        'python', script_path
    ]
    
    logger.info(f"Running command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)  # 2 hour timeout
        
        logger.info(f"Container exit code: {result.returncode}")
        logger.info(f"STDOUT:\n{result.stdout}")
        
        if result.stderr:
            logger.warning(f"STDERR:\n{result.stderr}")
        
        if result.returncode != 0:
            raise subprocess.CalledProcessError(result.returncode, cmd, result.stdout, result.stderr)
        
        return True
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Container processing failed: {e}")
        return False
    except subprocess.TimeoutExpired:
        logger.error("Container processing timed out")
        return False

def move_results():
    """Move results to final output directory"""
    
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Move result files
    result_files = [
        ("retnla_scgpt_result.h5ad", "retnla_scgpt_embeddings.h5ad"),
        ("retnla_scgpt_embeddings.npy", "retnla_scgpt_embeddings.npy")
    ]
    
    for temp_file, final_file in result_files:
        temp_path = f"{TEMP_DIR}/{temp_file}"
        final_path = f"{OUTPUT_DIR}/{final_file}"
        
        if os.path.exists(temp_path):
            subprocess.run(['mv', temp_path, final_path], check=True)
            logger.info(f"Moved: {final_path}")
        else:
            logger.warning(f"Result file not found: {temp_path}")

def main():
    """Main execution"""
    
    logger.info("Starting RETNLA scGPT processing via container...")
    
    try:
        # Prepare data
        combined_adata = prepare_retnla_data()
        data_path = save_for_container(combined_adata, TEMP_DIR)
        
        # Create container script
        script_path = create_container_script()
        
        # Run processing
        success = run_container_processing(data_path, script_path)
        
        if success:
            # Move results
            move_results()
            logger.info("RETNLA scGPT processing completed successfully!")
            logger.info(f"Results saved to: {OUTPUT_DIR}")
        else:
            logger.error("Container processing failed")
            return False
        
    except Exception as e:
        logger.error(f"Processing failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    return True

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)