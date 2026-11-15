#!/usr/bin/env python3
"""
Generate Geneformer Embeddings - Chunked Version

This script processes data in chunks to work within short_gpu time limits.
It saves checkpoints after each chunk so progress isn't lost.

Usage:
    python generate_geneformer_embeddings_chunked.py --dataset covid --chunk 0
    python generate_geneformer_embeddings_chunked.py --dataset covid --chunk 1
    ...
    python generate_geneformer_embeddings_chunked.py --dataset covid --merge

Author: Project 302006
Date: December 2025
"""

import os
import sys
import subprocess
from pathlib import Path
import logging
from datetime import datetime
import shutil
import argparse
import numpy as np

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)


class GeneformerChunkedEmbedding:
    """Generate Geneformer embeddings in chunks to work within time limits"""

    # Chunk size - process this many cells at a time
    CHUNK_SIZE = 5000  # ~5K cells per chunk, ~10-15 min on GPU

    def __init__(self, base_dir="/bigdata/godziklab/shared/Xinru/302006"):
        self.base_dir = Path(base_dir)
        self.input_dir = self.base_dir / "01_DATA/geneformer_ready"
        self.output_dir = self.base_dir / "02_EMBEDDINGS/geneformer"
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Chunk output directory
        self.chunk_dir = self.output_dir / "chunks"
        self.chunk_dir.mkdir(parents=True, exist_ok=True)

        # Singularity container path
        self.container_path = self.base_dir / "07_CONTAINERS/geneformer.sif"

        # Temp directory for container operations
        self.temp_dir = self.base_dir / "temp/geneformer_embeddings"
        self.temp_dir.mkdir(parents=True, exist_ok=True)

        # Check if container exists
        if not self.container_path.exists():
            logger.error(f"Geneformer container not found: {self.container_path}")
            raise FileNotFoundError(f"Container not found: {self.container_path}")

        logger.info(f"Container: {self.container_path}")
        logger.info(f"Input directory: {self.input_dir}")
        logger.info(f"Output directory: {self.output_dir}")
        logger.info(f"Chunk size: {self.CHUNK_SIZE} cells")

    def get_dataset_info(self, dataset_name):
        """Get info about dataset size"""
        # Check tokenized dataset
        tokenized_path = self.temp_dir / "container_output" / f"{dataset_name}_tokenized.dataset"
        if tokenized_path.exists():
            from datasets import load_from_disk
            dataset = load_from_disk(str(tokenized_path))
            return len(dataset)
        return None

    def create_tokenization_script(self, dataset_name):
        """Create Python script for tokenization inside container"""

        script = f'''#!/usr/bin/env python3
"""Tokenization script for {dataset_name} - runs inside Geneformer container"""

import sys
sys.path.insert(0, '/opt/Geneformer')

from geneformer import TranscriptomeTokenizer
import scanpy as sc
import anndata as ad
import numpy as np
import scipy.sparse as sp
import os

print("Loading data...")
adata = sc.read_h5ad("/data/input/{dataset_name}_geneformer_ready.h5ad")
print(f"Loaded {{adata.n_obs}} cells, {{adata.n_vars}} genes")

# Check for ensembl_id
if 'ensembl_id' not in adata.var.columns:
    print("ERROR: Missing ensembl_id column!")
    sys.exit(1)

print(f"var_names (first 5): {{adata.var_names[:5].tolist()}}")
print(f"ensembl_id (first 5): {{adata.var['ensembl_id'][:5].tolist()}}")

# Convert to dense array to avoid h5py fancy indexing issues
print("Converting matrix format for compatibility...")
if sp.issparse(adata.X):
    print(f"  Original format: {{type(adata.X).__name__}}")
    adata.X = np.array(adata.X.todense())
    print(f"  Converted to dense array: {{adata.X.shape}}")
else:
    adata.X = np.array(adata.X)
    print(f"  Matrix format: dense array")

# Ensure n_counts is present
if 'n_counts' not in adata.obs.columns:
    if 'total_counts' in adata.obs.columns:
        adata.obs['n_counts'] = adata.obs['total_counts']
    else:
        adata.obs['n_counts'] = np.array(adata.X.sum(axis=1)).flatten()
    print(f"  Added n_counts column")

# Save to temp location as loom file (avoids h5py fancy indexing issues)
temp_dir = "/data/temp/{dataset_name}_only"
os.makedirs(temp_dir, exist_ok=True)
temp_loom = f"{{temp_dir}}/{dataset_name}_for_tokenizer.loom"

# Write as loom file
print(f"Saving to {{temp_loom}}...")
adata.var['ensembl_id'] = adata.var_names.tolist()
adata.write_loom(temp_loom)

# Initialize tokenizer - V1 API
print("Initializing TranscriptomeTokenizer...")
try:
    tk = TranscriptomeTokenizer(
        custom_attr_name_dict={{"cell_type": "cell_type"}},
        nproc=4,
        model_input_size=4096,
        special_token=True
    )
    print("Using Geneformer V2 API")
except TypeError:
    tk = TranscriptomeTokenizer(
        custom_attr_name_dict={{"cell_type": "cell_type"}},
        nproc=4
    )
    print("Using Geneformer V1 API")

# Tokenize
print("Tokenizing data...")
tk.tokenize_data(
    data_directory=temp_dir,
    output_directory="/data/output",
    output_prefix="{dataset_name}_tokenized",
    file_format="loom"
)

print("Tokenization complete!")
print(f"Total cells tokenized: {{adata.n_obs}}")
'''
        return script

    def create_chunked_embedding_script(self, dataset_name, chunk_idx, chunk_size):
        """Create Python script for embedding extraction of a specific chunk"""

        script = f'''#!/usr/bin/env python3
"""Chunked embedding extraction script for {dataset_name} - chunk {chunk_idx}"""

import sys
sys.path.insert(0, '/opt/Geneformer')

import torch
import numpy as np
from datasets import load_from_disk

print("=" * 60)
print("GENEFORMER EMBEDDING EXTRACTION: {dataset_name} - CHUNK {chunk_idx}")
print("=" * 60)

# Try GPU first, fall back to CPU
if torch.cuda.is_available():
    device = torch.device("cuda")
    print(f"Using GPU: {{torch.cuda.get_device_name(0)}}")
else:
    device = torch.device("cpu")
    print("Using CPU for embedding extraction")

# Load tokenized dataset
tokenized_path = "/data/output/{dataset_name}_tokenized.dataset"
print(f"Loading tokenized data from: {{tokenized_path}}")

try:
    dataset = load_from_disk(tokenized_path)
    total_cells = len(dataset)
    print(f"Total cells in dataset: {{total_cells}}")
except Exception as e:
    print(f"ERROR loading tokenized data: {{e}}")
    sys.exit(1)

# Calculate chunk boundaries
chunk_size = {chunk_size}
chunk_idx = {chunk_idx}
start_idx = chunk_idx * chunk_size
end_idx = min(start_idx + chunk_size, total_cells)

if start_idx >= total_cells:
    print(f"Chunk {{chunk_idx}} is beyond dataset size. Nothing to process.")
    sys.exit(0)

print(f"Processing chunk {{chunk_idx}}: cells {{start_idx}} to {{end_idx}} ({{end_idx - start_idx}} cells)")

# Load model
from transformers import BertModel

print("Loading model from HuggingFace (ctheodoris/Geneformer)...")
model = BertModel.from_pretrained("ctheodoris/Geneformer")
model = model.to(device)
model.eval()

vocab_size = model.config.vocab_size
print(f"Model vocabulary size: {{vocab_size}}")

# Get max length for this chunk
chunk_data = dataset.select(range(start_idx, end_idx))

# Check token ID range and fix if needed
print("Checking token IDs...")
all_token_ids = [tid for seq in chunk_data['input_ids'] for tid in seq]
max_token_id = max(all_token_ids) if all_token_ids else 0
min_token_id = min(all_token_ids) if all_token_ids else 0
print(f"Token ID range: {{min_token_id}} - {{max_token_id}}")
print(f"Vocab size: {{vocab_size}}")

if max_token_id >= vocab_size:
    print(f"WARNING: Token IDs exceed vocab! Will clip to valid range.")

max_len = min(max(len(seq) for seq in chunk_data['input_ids']), 2048)
print(f"Max sequence length in chunk: {{max_len}}")

# Process chunk - use larger batch on GPU
batch_size = 16 if device.type == "cuda" else 4
all_embeddings = []

print(f"Extracting embeddings (batch_size={{batch_size}})...")
import time
start_time = time.time()

for i in range(0, len(chunk_data), batch_size):
    batch_data = chunk_data[i:i+batch_size]
    batch_seqs = batch_data['input_ids']

    # Pad sequences and clip token IDs
    padded_seqs = []
    attention_masks = []
    for seq in batch_seqs:
        seq = list(seq)[:max_len]
        # IMPORTANT: Clip token IDs to valid vocabulary range
        seq = [min(max(tid, 0), vocab_size - 1) for tid in seq]
        pad_len = max_len - len(seq)
        padded_seq = seq + [0] * pad_len
        mask = [1] * len(seq) + [0] * pad_len
        padded_seqs.append(padded_seq)
        attention_masks.append(mask)

    input_ids = torch.tensor(padded_seqs, dtype=torch.long).to(device)
    attention_mask = torch.tensor(attention_masks, dtype=torch.long).to(device)

    try:
        with torch.no_grad():
            outputs = model(input_ids, attention_mask=attention_mask)
            hidden_states = outputs.last_hidden_state
            mask_expanded = attention_mask.unsqueeze(-1).expand(hidden_states.size()).float()
            sum_embeddings = torch.sum(hidden_states * mask_expanded, dim=1)
            sum_mask = mask_expanded.sum(dim=1).clamp(min=1e-9)
            batch_emb = (sum_embeddings / sum_mask).cpu().numpy()
            all_embeddings.append(batch_emb)
    except RuntimeError as e:
        print(f"GPU error at batch {{i}}: {{e}}")
        print("Falling back to CPU for this batch...")
        model_cpu = model.cpu()
        input_ids_cpu = input_ids.cpu()
        attention_mask_cpu = attention_mask.cpu()
        with torch.no_grad():
            outputs = model_cpu(input_ids_cpu, attention_mask=attention_mask_cpu)
            hidden_states = outputs.last_hidden_state
            mask_expanded = attention_mask_cpu.unsqueeze(-1).expand(hidden_states.size()).float()
            sum_embeddings = torch.sum(hidden_states * mask_expanded, dim=1)
            sum_mask = mask_expanded.sum(dim=1).clamp(min=1e-9)
            batch_emb = (sum_embeddings / sum_mask).numpy()
            all_embeddings.append(batch_emb)
        model = model_cpu.to(device)

    if (i // batch_size) % 20 == 0 or i == 0:
        elapsed = time.time() - start_time
        cells_done = i + len(batch_seqs)
        rate = cells_done / elapsed if elapsed > 0 else 0
        eta = (len(chunk_data) - cells_done) / rate if rate > 0 else 0
        print(f"Processed {{cells_done}}/{{len(chunk_data)}} cells | Rate: {{rate:.1f}} cells/s | ETA: {{eta:.0f}}s")

embeddings = np.vstack(all_embeddings)
print(f"Chunk embeddings shape: {{embeddings.shape}}")
print(f"Total time: {{time.time() - start_time:.1f}}s")

# Save chunk embeddings
output_file = f"/data/output/{dataset_name}_chunk_{{chunk_idx}}.npy"
np.save(output_file, embeddings)
print(f"Saved to {{output_file}}")

# Also save metadata
metadata = {{
    'chunk_idx': chunk_idx,
    'start_idx': start_idx,
    'end_idx': end_idx,
    'n_cells': embeddings.shape[0],
    'embedding_dim': embeddings.shape[1]
}}
import json
with open(f"/data/output/{dataset_name}_chunk_{{chunk_idx}}_meta.json", 'w') as f:
    json.dump(metadata, f)

print("Done!")
'''
        return script

    def run_in_container(self, script_content, script_name, dataset_name):
        """Run a Python script inside the Singularity container"""

        container_input = self.temp_dir / "container_input"
        container_output = self.temp_dir / "container_output"
        container_temp = self.temp_dir / "container_temp"

        container_input.mkdir(parents=True, exist_ok=True)
        container_output.mkdir(parents=True, exist_ok=True)
        container_temp.mkdir(parents=True, exist_ok=True)

        # Copy input file if needed
        input_file = self.input_dir / f"{dataset_name}_geneformer_ready.h5ad"
        if input_file.exists():
            dest_file = container_input / f"{dataset_name}_geneformer_ready.h5ad"
            if not dest_file.exists():
                shutil.copy(input_file, dest_file)

        # Save script
        script_path = self.temp_dir / script_name
        with open(script_path, 'w') as f:
            f.write(script_content)

        # Run Singularity
        cmd = [
            'singularity', 'exec',
            '--nv',
            '--bind', f"{container_input}:/data/input",
            '--bind', f"{container_output}:/data/output",
            '--bind', f"{container_temp}:/data/temp",
            '--bind', f"{script_path}:/opt/run_script.py",
            str(self.container_path),
            'python', '/opt/run_script.py'
        ]

        logger.info(f"Running: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=7200  # 2 hour timeout
            )

            logger.info(f"STDOUT:\n{result.stdout}")
            if result.stderr:
                logger.warning(f"STDERR:\n{result.stderr}")

            if result.returncode != 0:
                raise RuntimeError(f"Container command failed with code {result.returncode}")

            return container_output

        except subprocess.TimeoutExpired:
            logger.error("Container operation timed out")
            raise

    def tokenize_dataset(self, dataset_name):
        """Tokenize dataset (run once)"""
        logger.info(f"\n{'='*80}")
        logger.info(f"TOKENIZING: {dataset_name}")
        logger.info(f"{'='*80}")

        # Check if already tokenized
        tokenized_path = self.temp_dir / "container_output" / f"{dataset_name}_tokenized.dataset"
        if tokenized_path.exists():
            logger.info(f"Dataset already tokenized: {tokenized_path}")
            return True

        # Clear temp
        container_temp = self.temp_dir / "container_temp"
        if container_temp.exists():
            shutil.rmtree(container_temp)

        tokenize_script = self.create_tokenization_script(dataset_name)
        self.run_in_container(tokenize_script, "tokenize.py", dataset_name)
        logger.info(f"✓ {dataset_name} tokenized!")
        return True

    def process_chunk(self, dataset_name, chunk_idx):
        """Process a single chunk"""
        logger.info(f"\n{'='*80}")
        logger.info(f"PROCESSING CHUNK {chunk_idx}: {dataset_name}")
        logger.info(f"{'='*80}")

        # Check if chunk already exists
        chunk_file = self.temp_dir / "container_output" / f"{dataset_name}_chunk_{chunk_idx}.npy"
        if chunk_file.exists():
            logger.info(f"Chunk {chunk_idx} already processed: {chunk_file}")
            return True

        embed_script = self.create_chunked_embedding_script(dataset_name, chunk_idx, self.CHUNK_SIZE)
        self.run_in_container(embed_script, f"embed_chunk_{chunk_idx}.py", dataset_name)

        # Copy chunk to output dir
        if chunk_file.exists():
            shutil.copy(chunk_file, self.chunk_dir / f"{dataset_name}_chunk_{chunk_idx}.npy")
            meta_file = self.temp_dir / "container_output" / f"{dataset_name}_chunk_{chunk_idx}_meta.json"
            if meta_file.exists():
                shutil.copy(meta_file, self.chunk_dir / f"{dataset_name}_chunk_{chunk_idx}_meta.json")

        logger.info(f"✓ Chunk {chunk_idx} processed!")
        return True

    def merge_chunks(self, dataset_name):
        """Merge all chunks into final embeddings"""
        import scanpy as sc
        import json

        logger.info(f"\n{'='*80}")
        logger.info(f"MERGING CHUNKS: {dataset_name}")
        logger.info(f"{'='*80}")

        # Find all chunk files
        chunk_files = sorted(self.chunk_dir.glob(f"{dataset_name}_chunk_*.npy"))
        if not chunk_files:
            # Check container output
            chunk_files = sorted((self.temp_dir / "container_output").glob(f"{dataset_name}_chunk_*.npy"))

        if not chunk_files:
            logger.error(f"No chunk files found for {dataset_name}")
            return False

        logger.info(f"Found {len(chunk_files)} chunk files")

        # Load and concatenate
        all_embeddings = []
        for chunk_file in chunk_files:
            emb = np.load(chunk_file)
            all_embeddings.append(emb)
            logger.info(f"  {chunk_file.name}: {emb.shape}")

        embeddings = np.vstack(all_embeddings)
        logger.info(f"Total embeddings shape: {embeddings.shape}")

        # Load original data and add embeddings
        input_file = self.input_dir / f"{dataset_name}_geneformer_ready.h5ad"
        adata = sc.read_h5ad(input_file)

        if embeddings.shape[0] != adata.n_obs:
            logger.error(f"Embedding count ({embeddings.shape[0]}) != cell count ({adata.n_obs})")
            return False

        adata.obsm['X_geneformer'] = embeddings
        output_h5ad = self.output_dir / f"{dataset_name}_geneformer.h5ad"
        adata.write_h5ad(output_h5ad)

        logger.info(f"✓ Saved: {output_h5ad}")
        logger.info(f"  Cells: {adata.n_obs}")
        logger.info(f"  Embedding dimensions: {embeddings.shape[1]}")

        return True

    def get_num_chunks(self, dataset_name):
        """Calculate number of chunks needed for dataset"""
        # Try to get from tokenized dataset
        tokenized_path = self.temp_dir / "container_output" / f"{dataset_name}_tokenized.dataset"
        if tokenized_path.exists():
            from datasets import load_from_disk
            dataset = load_from_disk(str(tokenized_path))
            n_cells = len(dataset)
        else:
            # Estimate from input file
            import scanpy as sc
            input_file = self.input_dir / f"{dataset_name}_geneformer_ready.h5ad"
            adata = sc.read_h5ad(input_file)
            n_cells = adata.n_obs

        n_chunks = (n_cells + self.CHUNK_SIZE - 1) // self.CHUNK_SIZE
        return n_chunks, n_cells


def main():
    parser = argparse.ArgumentParser(description='Generate Geneformer embeddings in chunks')
    parser.add_argument('--dataset', type=str, required=True,
                        choices=['covid', 'hek293t', 'hct116'],
                        help='Dataset to process')
    parser.add_argument('--chunk', type=int, default=None,
                        help='Chunk index to process (0-indexed)')
    parser.add_argument('--tokenize', action='store_true',
                        help='Only run tokenization')
    parser.add_argument('--merge', action='store_true',
                        help='Merge all chunks into final file')
    parser.add_argument('--info', action='store_true',
                        help='Show info about chunks needed')
    args = parser.parse_args()

    embedder = GeneformerChunkedEmbedding()

    if args.info:
        n_chunks, n_cells = embedder.get_num_chunks(args.dataset)
        print(f"\nDataset: {args.dataset}")
        print(f"Total cells: {n_cells}")
        print(f"Chunk size: {embedder.CHUNK_SIZE}")
        print(f"Number of chunks: {n_chunks}")
        print(f"\nTo process all chunks, run:")
        print(f"  # Step 1: Tokenize (once)")
        print(f"  python {sys.argv[0]} --dataset {args.dataset} --tokenize")
        print(f"  # Step 2: Process each chunk")
        for i in range(n_chunks):
            print(f"  python {sys.argv[0]} --dataset {args.dataset} --chunk {i}")
        print(f"  # Step 3: Merge chunks")
        print(f"  python {sys.argv[0]} --dataset {args.dataset} --merge")
        return 0

    if args.tokenize:
        embedder.tokenize_dataset(args.dataset)
        return 0

    if args.merge:
        embedder.merge_chunks(args.dataset)
        return 0

    if args.chunk is not None:
        # First ensure tokenization is done
        embedder.tokenize_dataset(args.dataset)
        # Then process the chunk
        embedder.process_chunk(args.dataset, args.chunk)
        return 0

    # Default: show usage
    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
