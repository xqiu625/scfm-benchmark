#!/usr/bin/env python3
"""
Fixed Geneformer embedding generation script
Addresses CUDA CUBLAS initialization issues
"""

import os
import pickle
import torch
import numpy as np
from datasets import load_from_disk
from transformers import BertModel, BertConfig
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def initialize_cuda():
    """Initialize CUDA environment properly"""
    try:
        # Force CUDA initialization
        if torch.cuda.is_available():
            torch.cuda.init()
            torch.cuda.empty_cache()
            
            # Initialize cuBLAS handle explicitly
            device = torch.device("cuda")
            dummy_tensor = torch.randn(2, 2, device=device)
            _ = torch.mm(dummy_tensor, dummy_tensor)  # Force cuBLAS initialization
            
            logger.info(f"CUDA initialized successfully on {torch.cuda.get_device_name()}")
            return device
        else:
            logger.warning("CUDA not available, using CPU")
            return torch.device("cpu")
            
    except Exception as e:
        logger.warning(f"CUDA initialization failed: {e}, falling back to CPU")
        return torch.device("cpu")

def process_with_fallback_strategies(model, input_tensor, device):
    """Try multiple processing strategies with fallbacks"""
    
    strategies = [
        ("Standard GPU processing", lambda: standard_processing(model, input_tensor, device)),
        ("Mixed precision", lambda: mixed_precision_processing(model, input_tensor, device)),
        ("CPU processing", lambda: cpu_processing(model, input_tensor)),
        ("Chunked processing", lambda: chunked_processing(model, input_tensor, device, chunk_size=1))
    ]
    
    for strategy_name, strategy_func in strategies:
        try:
            logger.info(f"Trying {strategy_name}...")
            embeddings = strategy_func()
            logger.info(f"{strategy_name} successful! Shape: {embeddings.shape}")
            return embeddings
            
        except Exception as e:
            logger.warning(f"{strategy_name} failed: {e}")
            continue
    
    raise RuntimeError("All processing strategies failed")

def standard_processing(model, input_tensor, device):
    """Standard GPU processing"""
    input_tensor = input_tensor.to(device)
    with torch.no_grad():
        outputs = model(input_tensor)
        return outputs.last_hidden_state.mean(dim=1).cpu().numpy()

def mixed_precision_processing(model, input_tensor, device):
    """Mixed precision processing to reduce memory pressure"""
    input_tensor = input_tensor.to(device)
    with torch.no_grad():
        with torch.amp.autocast('cuda'):
            outputs = model(input_tensor)
            return outputs.last_hidden_state.mean(dim=1).float().cpu().numpy()

def cpu_processing(model, input_tensor):
    """CPU-only processing as fallback"""
    model = model.cpu()
    input_tensor = input_tensor.cpu()
    with torch.no_grad():
        outputs = model(input_tensor)
        return outputs.last_hidden_state.mean(dim=1).numpy()

def chunked_processing(model, input_tensor, device, chunk_size=1):
    """Process samples in small chunks to avoid CUDA issues"""
    embeddings_list = []
    
    for i in range(0, input_tensor.shape[0], chunk_size):
        chunk = input_tensor[i:i+chunk_size].to(device)
        
        with torch.no_grad():
            outputs = model(chunk)
            chunk_embeddings = outputs.last_hidden_state.mean(dim=1).cpu().numpy()
            embeddings_list.append(chunk_embeddings)
    
    return np.vstack(embeddings_list)

def generate_embeddings():
    """Main embedding generation with robust error handling"""
    
    dataset_path = "/data/input/retnla_ortholog_tokenized.dataset"
    output_dir = "/data/output"
    
    logger.info("Loading tokenized dataset...")
    dataset = load_from_disk(dataset_path)
    logger.info(f"Loaded {len(dataset)} cells")
    
    # Initialize CUDA properly
    device = initialize_cuda()
    
    # Load model
    logger.info("Loading Geneformer model...")
    model_name = "ctheodoris/Geneformer"
    model = BertModel.from_pretrained(model_name, trust_remote_code=False)
    model = model.to(device)
    model.eval()
    logger.info(f"Model loaded on {device}")
    
    # Process in batches with robust error handling
    batch_size = 8  # Smaller batch size to avoid memory issues
    all_embeddings = []
    cell_metadata = []
    
    total_batches = (len(dataset) + batch_size - 1) // batch_size
    
    for batch_idx in range(total_batches):
        start_idx = batch_idx * batch_size
        end_idx = min((batch_idx + 1) * batch_size, len(dataset))
        
        logger.info(f"Processing batch {batch_idx + 1}/{total_batches} (cells {start_idx}-{end_idx-1})")
        
        try:
            # Get batch data
            batch_dict = dataset[start_idx:end_idx]
            input_ids = batch_dict['input_ids']
            
            # Convert to tensor
            input_tensor = torch.tensor(input_ids, dtype=torch.long)
            
            # Process with fallback strategies
            batch_embeddings = process_with_fallback_strategies(model, input_tensor, device)
            
            all_embeddings.append(batch_embeddings)
            
            # Store metadata
            for i in range(len(batch_dict['species'])):
                cell_metadata.append({
                    'species': batch_dict['species'][i],
                    'gene_target': batch_dict['gene_target'][i],
                    'treatment_status': batch_dict['treatment_status'][i],
                    'dataset': batch_dict['dataset'][i],
                    'cell_index': start_idx + i
                })
            
            logger.info(f"Batch {batch_idx + 1} completed successfully")
            
        except Exception as e:
            logger.error(f"CRITICAL: Batch {batch_idx + 1} failed: {e}")
            logger.error(f"Cells {start_idx} to {end_idx-1} could not be processed.")
            logger.error("Aborting to prevent synthetic data contamination.")
            raise RuntimeError(
                f"Batch {batch_idx + 1} failed after all fallback strategies. "
                f"Cannot generate embeddings for cells {start_idx}-{end_idx-1}. "
                f"Original error: {e}"
            )
    
    # Combine all embeddings
    final_embeddings = np.vstack(all_embeddings)
    logger.info(f"Generated embeddings shape: {final_embeddings.shape}")
    
    # Save results
    os.makedirs(output_dir, exist_ok=True)
    
    # Save embeddings
    embeddings_file = os.path.join(output_dir, "retnla_ortholog_geneformer_embeddings.npy")
    np.save(embeddings_file, final_embeddings)
    
    embeddings_pkl = os.path.join(output_dir, "retnla_ortholog_geneformer_embeddings.pkl")
    with open(embeddings_pkl, 'wb') as f:
        pickle.dump(final_embeddings, f)
    
    # Save metadata
    metadata_file = os.path.join(output_dir, "cell_metadata.csv")
    import pandas as pd
    pd.DataFrame(cell_metadata).to_csv(metadata_file, index=False)
    
    # Save processing info
    processing_info = {
        'total_cells': len(dataset),
        'embedding_dim': final_embeddings.shape[1],
        'batch_size': batch_size,
        'device_used': str(device),
        'model_name': model_name
    }
    
    with open(os.path.join(output_dir, "processing_metadata.pkl"), 'wb') as f:
        pickle.dump(processing_info, f)
    
    # Create summary file
    with open(os.path.join(output_dir, "embedding_info.txt"), 'w') as f:
        f.write(f"Geneformer Embeddings Generated Successfully\n")
        f.write(f"Total cells: {len(dataset)}\n")
        f.write(f"Embedding dimensions: {final_embeddings.shape[1]}\n")
        f.write(f"Final embedding shape: {final_embeddings.shape}\n")
        f.write(f"Processing device: {device}\n")
        f.write(f"Model: {model_name}\n")
    
    logger.info("All output files saved successfully!")
    
    files_created = [
        "retnla_ortholog_geneformer_embeddings.npy",
        "retnla_ortholog_geneformer_embeddings.pkl", 
        "cell_metadata.csv",
        "processing_metadata.pkl",
        "embedding_info.txt"
    ]
    
    logger.info("Files created:")
    for file in files_created:
        logger.info(f"  - {file}")
    
    return final_embeddings

if __name__ == "__main__":
    try:
        embeddings = generate_embeddings()
        logger.info("Manual embedding extraction completed successfully!")
        
    except Exception as e:
        logger.error(f"Embedding generation failed: {e}")
        exit(1)





# # Replace the container script with the fixed version
# cp generate_geneformer_embeddings_fixed.py /bigdata/godziklab/shared/Xinru/302006/temp/geneformer_embeddings/generate_embeddings.py
# 
# # Run the fixed version
# singularity exec --nv \
#   --bind /bigdata/godziklab/shared/Xinru/302006/temp/geneformer_embeddings/container_input:/data/input \
#   --bind /bigdata/godziklab/shared/Xinru/302006/temp/geneformer_embeddings/container_output:/data/output \
#   --bind /bigdata/godziklab/shared/Xinru/302006/temp/geneformer_embeddings/generate_embeddings.py:/opt/generate_embeddings.py \
#   /bigdata/godziklab/shared/Xinru/302006/scfm/geneformer.sif python /opt/generate_embeddings.py