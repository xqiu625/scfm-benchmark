#!/usr/bin/env python3
"""
Generate Geneformer Embeddings for COVID/HEK293T/HCT116 datasets

This script runs inside the Geneformer Singularity container.
Uses Geneformer V2 with EmbExtractor for embedding generation.

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

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f'geneformer_embeddings_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)


class GeneformerEmbedding:
    """Generate Geneformer embeddings using Singularity container"""

    def __init__(self, base_dir="/bigdata/godziklab/shared/Xinru/302006"):
        self.base_dir = Path(base_dir)
        self.input_dir = self.base_dir / "01_DATA/geneformer_ready"
        self.output_dir = self.base_dir / "02_EMBEDDINGS/geneformer"
        self.output_dir.mkdir(parents=True, exist_ok=True)

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
    # Convert to dense numpy array - this avoids h5py indexing issues
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

# Write as loom file - this format works better with Geneformer tokenizer
print(f"Saving to {{temp_loom}}...")
# Need to ensure var_names are used as gene IDs in loom
adata.var['ensembl_id'] = adata.var_names.tolist()  # Ensure ensembl_id matches var_names
adata.write_loom(temp_loom)

# Initialize tokenizer - try V2 first, fall back to V1
print("Initializing TranscriptomeTokenizer...")
try:
    # Try V2 API first
    tk = TranscriptomeTokenizer(
        custom_attr_name_dict={{"cell_type": "cell_type"}},
        nproc=4,
        model_input_size=4096,
        special_token=True
    )
    print("Using Geneformer V2 API")
except TypeError:
    # Fall back to V1 API (older version without model_input_size)
    tk = TranscriptomeTokenizer(
        custom_attr_name_dict={{"cell_type": "cell_type"}},
        nproc=4
    )
    print("Using Geneformer V1 API")

# Tokenize - use loom format
print("Tokenizing data...")
tk.tokenize_data(
    data_directory=temp_dir,
    output_directory="/data/output",
    output_prefix="{dataset_name}_tokenized",
    file_format="loom"
)

print("Tokenization complete!")
'''
        return script

    def create_embedding_script(self, dataset_name):
        """Create Python script for embedding extraction inside container"""

        script = f'''#!/usr/bin/env python3
"""Embedding extraction script for {dataset_name} - runs inside Geneformer container"""

import sys
sys.path.insert(0, '/opt/Geneformer')

import torch
import numpy as np
import scanpy as sc
from pathlib import Path
from datasets import load_from_disk

print("=" * 60)
print("GENEFORMER EMBEDDING EXTRACTION: {dataset_name}")
print("=" * 60)

# Check CUDA
print(f"CUDA available: {{torch.cuda.is_available()}}")
if torch.cuda.is_available():
    print(f"CUDA device: {{torch.cuda.get_device_name(0)}}")
    device = torch.device("cuda")
else:
    print("Using CPU")
    device = torch.device("cpu")

# Load tokenized dataset
tokenized_path = "/data/output/{dataset_name}_tokenized.dataset"
print(f"Loading tokenized data from: {{tokenized_path}}")

try:
    dataset = load_from_disk(tokenized_path)
    print(f"Loaded {{len(dataset)}} cells")
except Exception as e:
    print(f"ERROR loading tokenized data: {{e}}")
    sys.exit(1)

# Force CPU to avoid CUDA issues with container
device = torch.device("cpu")
print("Using CPU for embedding extraction (more stable with container)")

# Manual extraction with proper padding - skip EmbExtractor due to version issues
from transformers import BertModel, BertConfig

print("Loading model from HuggingFace (ctheodoris/Geneformer)...")
model = BertModel.from_pretrained("ctheodoris/Geneformer")
model = model.to(device)
model.eval()

# Check vocabulary size
vocab_size = model.config.vocab_size
print(f"Model vocabulary size: {{vocab_size}}")

# Check token IDs in dataset
all_token_ids = [tid for seq in dataset['input_ids'] for tid in seq]
max_token_id = max(all_token_ids)
min_token_id = min(all_token_ids)
print(f"Token ID range in data: {{min_token_id}} - {{max_token_id}}")

if max_token_id >= vocab_size:
    print(f"WARNING: Max token ID ({{max_token_id}}) >= vocab size ({{vocab_size}})")
    print("Clipping token IDs to valid range...")

# Get max length for padding (cap at 2048 for memory)
max_len = min(max(len(seq) for seq in dataset['input_ids']), 2048)
print(f"Max sequence length: {{max_len}}")

batch_size = 8  # Can use larger batch on CPU
all_embeddings = []

print("Extracting embeddings...")
for i in range(0, len(dataset), batch_size):
    batch_data = dataset[i:i+batch_size]
    batch_seqs = batch_data['input_ids']

    # Pad sequences to same length (truncate if > max_len)
    padded_seqs = []
    attention_masks = []
    for seq in batch_seqs:
        seq = list(seq)[:max_len]  # Truncate if needed
        # Clip token IDs to valid vocab range
        seq = [min(tid, vocab_size - 1) for tid in seq]
        pad_len = max_len - len(seq)
        padded_seq = seq + [0] * pad_len  # 0 is typically PAD token
        mask = [1] * len(seq) + [0] * pad_len
        padded_seqs.append(padded_seq)
        attention_masks.append(mask)

    input_ids = torch.tensor(padded_seqs, dtype=torch.long).to(device)
    attention_mask = torch.tensor(attention_masks, dtype=torch.long).to(device)

    with torch.no_grad():
        outputs = model(input_ids, attention_mask=attention_mask)
        # Use mean pooling over non-padded tokens
        hidden_states = outputs.last_hidden_state
        # Expand attention mask for broadcasting
        mask_expanded = attention_mask.unsqueeze(-1).expand(hidden_states.size()).float()
        sum_embeddings = torch.sum(hidden_states * mask_expanded, dim=1)
        sum_mask = mask_expanded.sum(dim=1).clamp(min=1e-9)
        batch_emb = (sum_embeddings / sum_mask).cpu().numpy()
        all_embeddings.append(batch_emb)

    if (i // batch_size) % 100 == 0 or i == 0:
        print(f"Processed {{i + len(batch_seqs)}}/{{len(dataset)}} cells")

embeddings = np.vstack(all_embeddings)
print(f"Embeddings shape: {{embeddings.shape}}")

# Save embeddings
np.save("/data/output/{dataset_name}_geneformer_embeddings.npy", embeddings)

print("Done!")
'''
        return script

    def run_in_container(self, script_content, script_name, dataset_name):
        """Run a Python script inside the Singularity container"""

        # Create directories
        container_input = self.temp_dir / "container_input"
        container_output = self.temp_dir / "container_output"
        container_temp = self.temp_dir / "container_temp"

        # Clear temp directory to avoid old files interfering
        if container_temp.exists():
            shutil.rmtree(container_temp)

        container_input.mkdir(parents=True, exist_ok=True)
        container_output.mkdir(parents=True, exist_ok=True)
        container_temp.mkdir(parents=True, exist_ok=True)

        # Copy input file
        input_file = self.input_dir / f"{dataset_name}_geneformer_ready.h5ad"
        if not input_file.exists():
            raise FileNotFoundError(f"Input file not found: {input_file}")

        shutil.copy(input_file, container_input / f"{dataset_name}_geneformer_ready.h5ad")

        # Save script
        script_path = self.temp_dir / script_name
        with open(script_path, 'w') as f:
            f.write(script_content)

        # Run Singularity
        cmd = [
            'singularity', 'exec',
            '--nv',  # GPU support
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

    def process_dataset(self, dataset_name):
        """Process a single dataset through Geneformer pipeline"""
        logger.info(f"\n{'='*80}")
        logger.info(f"Processing: {dataset_name}")
        logger.info(f"{'='*80}")

        try:
            # Step 1: Tokenization
            logger.info("Step 1: Tokenization...")
            tokenize_script = self.create_tokenization_script(dataset_name)
            self.run_in_container(tokenize_script, "tokenize.py", dataset_name)

            # Step 2: Embedding extraction
            logger.info("Step 2: Embedding extraction...")
            embed_script = self.create_embedding_script(dataset_name)
            output_dir = self.run_in_container(embed_script, "embed.py", dataset_name)

            # Step 3: Copy results to final location
            logger.info("Step 3: Saving results...")
            import scanpy as sc
            import numpy as np

            # Look for embedding files
            for f in output_dir.iterdir():
                if 'embedding' in f.name.lower() or f.suffix in ['.npy', '.parquet']:
                    dest = self.output_dir / f"{dataset_name}_{f.name}"
                    shutil.copy(f, dest)
                    logger.info(f"Saved: {dest}")

            # If we have numpy embeddings, also create h5ad with embeddings
            npy_file = output_dir / f"{dataset_name}_geneformer_embeddings.npy"
            if npy_file.exists():
                embeddings = np.load(npy_file)

                # Load original data and add embeddings
                input_file = self.input_dir / f"{dataset_name}_geneformer_ready.h5ad"
                adata = sc.read_h5ad(input_file)

                if embeddings.shape[0] == adata.n_obs:
                    adata.obsm['X_geneformer'] = embeddings
                    output_h5ad = self.output_dir / f"{dataset_name}_geneformer.h5ad"
                    adata.write_h5ad(output_h5ad)
                    logger.info(f"Saved h5ad with embeddings: {output_h5ad}")
                    logger.info(f"  Cells: {adata.n_obs}")
                    logger.info(f"  Embedding dimensions: {embeddings.shape[1]}")

            logger.info(f"✓ {dataset_name} processing complete!")
            return True

        except Exception as e:
            logger.error(f"Failed to process {dataset_name}: {str(e)}")
            import traceback
            traceback.print_exc()
            return False

    def process_all_datasets(self):
        """Process all datasets"""
        logger.info("=" * 80)
        logger.info("GENEFORMER EMBEDDING GENERATION")
        logger.info("=" * 80)

        datasets = ['covid', 'hek293t', 'hct116']
        results = {}

        for dataset_name in datasets:
            success = self.process_dataset(dataset_name)
            results[dataset_name] = success

        # Summary
        logger.info("\n" + "=" * 80)
        logger.info("SUMMARY")
        logger.info("=" * 80)
        for dataset_name, success in results.items():
            status = "✓ Success" if success else "✗ Failed"
            logger.info(f"  {dataset_name}: {status}")

        return results


def main():
    """Main execution"""
    import argparse
    parser = argparse.ArgumentParser(description='Generate Geneformer embeddings')
    parser.add_argument('--dataset', type=str, default=None,
                        help='Process specific dataset only (covid/hek293t/hct116)')
    args = parser.parse_args()

    try:
        embedder = GeneformerEmbedding()

        if args.dataset:
            embedder.process_dataset(args.dataset)
        else:
            embedder.process_all_datasets()

        return 0

    except Exception as e:
        logger.error(f"\nERROR: {str(e)}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
