#!/usr/bin/env python3
"""
Process RETNLA Ortholog Cells with Geneformer via Singularity Container

This script processes 890 RETNLA ortholog cells (mouse + human) through the 
Geneformer model using the successfully deployed Singularity container.

Container: /bigdata/godziklab/shared/Xinru/302006/scfm/geneformer.sif
Data files:
- /bigdata/godziklab/shared/Xinru/302006/data/raw/HCT116_retn_expressing_vs_targeted.h5ad
- /bigdata/godziklab/shared/Xinru/302006/data/raw/HEK293T_retn_expressing_vs_targeted.h5ad  
- /bigdata/godziklab/shared/Xinru/302006/data/raw/mouse_epithelial_retnla_ko_wt.h5ad
"""

import os
import sys
import subprocess
import pandas as pd
import scanpy as sc
import numpy as np
from pathlib import Path
import logging

# Ortholog mapping file
ORTHOLOG_FILE = "/bigdata/godziklab/shared/Xinru/302006/data/external/genome_wide_orthologs/genome_wide_human_mouse_orthologs_complete_20250722_174716.tsv"

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Project paths
PROJECT_BASE = "/bigdata/godziklab/shared/Xinru/302006"
DATA_DIR = f"{PROJECT_BASE}/data/raw"
CONTAINER_PATH = f"{PROJECT_BASE}/scfm/geneformer.sif"
OUTPUT_DIR = f"{PROJECT_BASE}/results/geneformer"
TEMP_DIR = f"{PROJECT_BASE}/temp/geneformer"

# Create output directories
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(TEMP_DIR, exist_ok=True)

class GeneformerProcessor:
    def __init__(self):
        self.container_path = CONTAINER_PATH
        self.output_dir = OUTPUT_DIR
        self.temp_dir = TEMP_DIR
        
        # Data file paths
        self.data_files = {
            'hct116': f"{DATA_DIR}/HCT116_retn_expressing_vs_targeted.h5ad",
            'hek293t': f"{DATA_DIR}/HEK293T_retn_expressing_vs_targeted.h5ad", 
            'mouse': f"{DATA_DIR}/mouse_epithelial_retnla_ko_wt.h5ad"
        }
        
        # Load ortholog mapping
        self.ortholog_mapping = self.load_ortholog_mapping()
        
        # Verify container exists
        if not os.path.exists(self.container_path):
            raise FileNotFoundError(f"Geneformer container not found: {self.container_path}")
            
        logger.info(f"Geneformer container found: {self.container_path}")
    
    def load_ortholog_mapping(self):
        """Load genome-wide human-mouse ortholog mapping"""
        if not os.path.exists(ORTHOLOG_FILE):
            raise FileNotFoundError(f"Ortholog mapping file not found: {ORTHOLOG_FILE}")
        
        logger.info(f"Loading ortholog mapping from: {ORTHOLOG_FILE}")
        ortholog_df = pd.read_csv(ORTHOLOG_FILE, sep='\t')
        
        # Create mouse->human gene symbol mapping
        mouse_to_human = dict(zip(
            ortholog_df['mouse_gene_symbol'], 
            ortholog_df['human_gene_symbol']
        ))
        
        logger.info(f"Loaded {len(mouse_to_human)} mouse->human ortholog mappings")
        logger.info(f"Sample mappings: {dict(list(mouse_to_human.items())[:5])}")
        
        return mouse_to_human
    
    def convert_mouse_genes_to_human_orthologs(self, mouse_adata):
        """Convert mouse gene symbols to human orthologs for cross-species analysis"""
        
        logger.info("Converting mouse genes to human orthologs...")
        
        # Get current mouse gene names
        mouse_genes = mouse_adata.var.index.tolist()
        logger.info(f"Mouse dataset has {len(mouse_genes)} genes")
        
        # Find orthologs and handle duplicates
        ortholog_mapping_dict = {}
        genes_without_orthologs = []
        
        for mouse_gene in mouse_genes:
            if mouse_gene in self.ortholog_mapping:
                human_ortholog = self.ortholog_mapping[mouse_gene]
                
                # Handle duplicates by keeping the gene with highest mean expression
                if human_ortholog in ortholog_mapping_dict:
                    existing_mouse_gene = ortholog_mapping_dict[human_ortholog]
                    
                    # Get mean expression for both genes
                    existing_expr = mouse_adata[:, existing_mouse_gene].X.mean()
                    current_expr = mouse_adata[:, mouse_gene].X.mean()
                    
                    if current_expr > existing_expr:
                        logger.info(f"Duplicate ortholog '{human_ortholog}': replacing '{existing_mouse_gene}' (mean={existing_expr:.3f}) with '{mouse_gene}' (mean={current_expr:.3f})")
                        ortholog_mapping_dict[human_ortholog] = mouse_gene
                    else:
                        logger.info(f"Duplicate ortholog '{human_ortholog}': keeping '{existing_mouse_gene}' (mean={existing_expr:.3f}) over '{mouse_gene}' (mean={current_expr:.3f})")
                else:
                    ortholog_mapping_dict[human_ortholog] = mouse_gene
            else:
                genes_without_orthologs.append(mouse_gene)
        
        # Extract final mappings
        mouse_genes_with_orthologs = list(ortholog_mapping_dict.values())
        ortholog_genes = list(ortholog_mapping_dict.keys())
        
        logger.info(f"Found orthologs for {len(ortholog_genes)}/{len(mouse_genes)} mouse genes")
        logger.info(f"Genes without orthologs: {len(genes_without_orthologs)}")
        logger.info(f"Duplicate orthologs resolved: {len(mouse_genes) - len(genes_without_orthologs) - len(ortholog_genes)}")
        
        if len(genes_without_orthologs) > 0:
            logger.info(f"First 10 genes without orthologs: {genes_without_orthologs[:10]}")
        
        # Create filtered mouse dataset with only unique ortholog genes
        mouse_ortholog_mask = mouse_adata.var.index.isin(mouse_genes_with_orthologs)
        mouse_adata_filtered = mouse_adata[:, mouse_ortholog_mask].copy()
        
        # Reorder to match ortholog mapping
        mouse_adata_filtered = mouse_adata_filtered[:, mouse_genes_with_orthologs]
        
        # Update gene names to human orthologs
        mouse_adata_filtered.var.index = ortholog_genes
        
        # Add ortholog conversion info to var
        mouse_adata_filtered.var['original_mouse_symbol'] = mouse_genes_with_orthologs
        mouse_adata_filtered.var['human_ortholog_symbol'] = ortholog_genes
        
        logger.info(f"Mouse dataset after ortholog conversion: {mouse_adata_filtered.n_obs} cells, {mouse_adata_filtered.n_vars} genes")
        
        # Special validation for RETNLA -> RETNLB/RETN conversion
        if 'Retnla' in mouse_genes_with_orthologs:
            retnla_idx = mouse_genes_with_orthologs.index('Retnla')
            human_ortholog = ortholog_genes[retnla_idx]
            logger.info(f"RETNLA ortholog validation: mouse 'Retnla' -> human '{human_ortholog}'")
        
        return mouse_adata_filtered
    def verify_data_files(self):
        """Verify all input data files exist and get cell counts"""
        total_cells = 0
        file_info = {}
        
        for dataset, filepath in self.data_files.items():
            if not os.path.exists(filepath):
                raise FileNotFoundError(f"Data file not found: {filepath}")
            
            # Load and inspect the data
            adata = sc.read_h5ad(filepath)
            n_cells = adata.n_obs
            total_cells += n_cells
            
            file_info[dataset] = {
                'path': filepath,
                'n_cells': n_cells,
                'n_genes': adata.n_vars,
                'obs_columns': list(adata.obs.columns)
            }
            
            logger.info(f"{dataset}: {n_cells} cells, {adata.n_vars} genes")
            logger.info(f"  Obs columns: {list(adata.obs.columns)}")
            
        logger.info(f"Total RETNLA ortholog cells: {total_cells}")
        return file_info, total_cells
    
    def prepare_combined_dataset(self, file_info):
        """Combine all datasets with proper metadata for cross-species analysis"""
        combined_data = []
        
        for dataset, info in file_info.items():
            logger.info(f"Loading {dataset} data...")
            adata = sc.read_h5ad(info['path'])
            
            # Add dataset and species metadata
            if dataset == 'mouse':
                # Convert mouse genes to human orthologs FIRST
                adata = self.convert_mouse_genes_to_human_orthologs(adata)
                
                adata.obs['species'] = 'mouse'
                adata.obs['gene_target'] = 'RETNLA'
                # Map condition to standardized names
                if 'condition' in adata.obs.columns:
                    adata.obs['treatment_status'] = adata.obs['condition']
                elif 'treatment' in adata.obs.columns:
                    adata.obs['treatment_status'] = adata.obs['treatment']
                elif 'group' in adata.obs.columns:
                    adata.obs['treatment_status'] = adata.obs['group']
                else:
                    logger.warning(f"No condition/treatment/group column found in {dataset}")
                    logger.warning(f"Available columns: {list(adata.obs.columns)}")
            else:
                adata.obs['species'] = 'human'
                adata.obs['cell_line'] = dataset.upper()
                adata.obs['gene_target'] = 'RETN'  # Both HCT116 and HEK293T target RETN
                
                # Map condition/treatment to standardized format
                if 'condition' in adata.obs.columns:
                    adata.obs['treatment_status'] = adata.obs['condition']
                elif 'treatment' in adata.obs.columns:
                    adata.obs['treatment_status'] = adata.obs['treatment']
                elif 'group' in adata.obs.columns:
                    adata.obs['treatment_status'] = adata.obs['group']
                elif 'cell_group' in adata.obs.columns:
                    adata.obs['treatment_status'] = adata.obs['cell_group']
                else:
                    logger.warning(f"No condition/treatment/group/cell_group column found in {dataset}")
                    logger.warning(f"Available columns: {list(adata.obs.columns)}")
            
            # Add dataset identifier
            adata.obs['dataset'] = dataset
            
            combined_data.append(adata)
        
        # Concatenate all datasets
        logger.info("Combining datasets...")
        combined_adata = sc.concat(combined_data, join='outer', index_unique='_')
        
        # Log combined dataset info
        logger.info(f"Combined dataset: {combined_adata.n_obs} cells, {combined_adata.n_vars} genes")
        logger.info(f"Species distribution: {combined_adata.obs['species'].value_counts().to_dict()}")
        
        if 'treatment_status' in combined_adata.obs.columns:
            logger.info(f"Treatment status distribution: {combined_adata.obs['treatment_status'].value_counts().to_dict()}")
        
        return combined_adata
    
    def prepare_geneformer_format(self, combined_adata):
        """Prepare data in format expected by Geneformer TranscriptomeTokenizer"""
        
        # Add required Ensembl ID column for Geneformer
        # For human orthologs, we need to map gene symbols to Ensembl IDs
        logger.info("Adding Ensembl ID annotations for Geneformer...")
        
        # Load human Ensembl ID mapping from ortholog file
        ortholog_df = pd.read_csv(ORTHOLOG_FILE, sep='\t')
        symbol_to_ensembl = dict(zip(
            ortholog_df['human_gene_symbol'], 
            ortholog_df['human_ensembl_id']
        ))
        
        # Map gene symbols to Ensembl IDs
        ensembl_ids = []
        genes_without_ensembl = []
        
        for gene_symbol in combined_adata.var.index:
            if gene_symbol in symbol_to_ensembl:
                ensembl_ids.append(symbol_to_ensembl[gene_symbol])
            else:
                # For genes without Ensembl mapping, use the symbol itself
                ensembl_ids.append(gene_symbol)
                genes_without_ensembl.append(gene_symbol)
        
        # Add ensembl_id column to var
        combined_adata.var['ensembl_id'] = ensembl_ids
        
        logger.info(f"Mapped {len(combined_adata.var) - len(genes_without_ensembl)}/{len(combined_adata.var)} genes to Ensembl IDs")
        if len(genes_without_ensembl) > 0:
            logger.info(f"Genes without Ensembl mapping: {len(genes_without_ensembl)} (using gene symbols)")
            logger.info(f"First 10 genes without Ensembl IDs: {genes_without_ensembl[:10]}")
        
        # Add required cell metadata columns for Geneformer
        logger.info("Adding required cell metadata columns for Geneformer...")
        
        # Add n_counts (total UMI counts per cell)
        if 'n_counts' not in combined_adata.obs.columns:
            if 'total_counts' in combined_adata.obs.columns:
                combined_adata.obs['n_counts'] = combined_adata.obs['total_counts']
            else:
                # Calculate total counts per cell
                combined_adata.obs['n_counts'] = np.array(combined_adata.X.sum(axis=1)).flatten()
        
        # Add filter_pass (quality control filter)
        if 'filter_pass' not in combined_adata.obs.columns:
            # Assume all cells pass filter since they're already in our curated dataset
            combined_adata.obs['filter_pass'] = 1
        
        # Add cell_type if not present (using treatment_status as proxy)
        if 'cell_type' not in combined_adata.obs.columns:
            if 'treatment_status' in combined_adata.obs.columns:
                combined_adata.obs['cell_type'] = combined_adata.obs['treatment_status']
            else:
                combined_adata.obs['cell_type'] = 'unknown'
        
        logger.info("Required Geneformer columns added:")
        logger.info(f"  n_counts: {combined_adata.obs['n_counts'].describe()}")
        logger.info(f"  filter_pass: {combined_adata.obs['filter_pass'].value_counts().to_dict()}")
        logger.info(f"  cell_type distribution: {combined_adata.obs['cell_type'].value_counts().to_dict()}")
        
        # Save combined dataset for Geneformer processing
        combined_path = f"{self.temp_dir}/retnla_ortholog_combined.h5ad"
        combined_adata.write_h5ad(combined_path)
        logger.info(f"Combined dataset saved: {combined_path}")
        
        # Create metadata CSV for Geneformer
        metadata_df = combined_adata.obs.copy()
        metadata_df.index.name = 'cell_id'
        metadata_path = f"{self.temp_dir}/retnla_ortholog_metadata.csv"
        metadata_df.to_csv(metadata_path)
        logger.info(f"Metadata saved: {metadata_path}")
        
        return combined_path, metadata_path
    
    def run_geneformer_tokenizer(self, data_path, metadata_path):
        """Run TranscriptomeTokenizer via Singularity container"""
        
        # Create tokenizer script to run inside container
        tokenizer_script = f"""
import sys
sys.path.append('/opt/Geneformer')

from geneformer import TranscriptomeTokenizer
import scanpy as sc
import pandas as pd
import os

# Load data
print("Loading combined RETNLA ortholog dataset...")
combined_adata = sc.read_h5ad("/data/input/retnla_ortholog_combined.h5ad")
print(f"Loaded {{combined_adata.n_obs}} cells, {{combined_adata.n_vars}} genes")

# Initialize tokenizer
print("Initializing TranscriptomeTokenizer...")
tk = TranscriptomeTokenizer(
    custom_attr_name_dict={{
        "species": "species",
        "gene_target": "gene_target", 
        "treatment_status": "treatment_status",
        "dataset": "dataset"
    }},
    nproc=4  # 4 processes as tested successful
)

# Tokenize data
print("Tokenizing RETNLA ortholog cells...")
tk.tokenize_data(
    data_directory="/data/input",
    output_directory="/data/output", 
    output_prefix="retnla_ortholog_tokenized",
    file_format="h5ad"
)

print("Tokenization complete!")
"""
        
        # Save tokenizer script
        script_path = f"{self.temp_dir}/run_tokenizer.py"
        with open(script_path, 'w') as f:
            f.write(tokenizer_script)
        
        # Create input/output directories for container
        container_input = f"{self.temp_dir}/container_input"
        container_output = f"{self.temp_dir}/container_output"
        os.makedirs(container_input, exist_ok=True)
        os.makedirs(container_output, exist_ok=True)
        
        # Copy files to container input
        subprocess.run([
            'cp', data_path, f"{container_input}/retnla_ortholog_combined.h5ad"
        ], check=True)
        
        subprocess.run([
            'cp', metadata_path, f"{container_input}/retnla_ortholog_metadata.csv" 
        ], check=True)
        
        # Run tokenization via Singularity
        logger.info("Running Geneformer TranscriptomeTokenizer via Singularity...")
        
        singularity_cmd = [
            'singularity', 'exec',
            '--nv',  # GPU support
            '--bind', f"{container_input}:/data/input",
            '--bind', f"{container_output}:/data/output",
            '--bind', f"{script_path}:/opt/run_tokenizer.py",
            self.container_path,
            'python', '/opt/run_tokenizer.py'
        ]
        
        logger.info(f"Running command: {' '.join(singularity_cmd)}")
        
        try:
            result = subprocess.run(
                singularity_cmd,
                capture_output=True,
                text=True,
                check=True,
                timeout=3600  # 1 hour timeout
            )
            
            logger.info("Tokenization completed successfully!")
            logger.info(f"STDOUT:\n{result.stdout}")
            
            if result.stderr:
                logger.warning(f"STDERR:\n{result.stderr}")
                
        except subprocess.CalledProcessError as e:
            logger.error(f"Tokenization failed with return code {e.returncode}")
            logger.error(f"STDOUT:\n{e.stdout}")
            logger.error(f"STDERR:\n{e.stderr}")
            raise
        except subprocess.TimeoutExpired:
            logger.error("Tokenization timed out after 1 hour")
            raise
        
        return container_output
    
    def verify_tokenization_output(self, output_dir):
        """Verify tokenization completed successfully"""
        
        # Look for tokenized output files
        output_files = os.listdir(output_dir)
        tokenized_files = [f for f in output_files if 'tokenized' in f.lower()]
        
        if not tokenized_files:
            raise RuntimeError("No tokenized files found in output directory")
        
        logger.info(f"Tokenization output files: {tokenized_files}")
        
        # Move tokenized files to final output directory
        for filename in tokenized_files:
            src = os.path.join(output_dir, filename)
            dst = os.path.join(self.output_dir, filename)
            subprocess.run(['mv', src, dst], check=True)
            logger.info(f"Moved tokenized file: {dst}")
        
        return [os.path.join(self.output_dir, f) for f in tokenized_files]
    
    def process_all_data(self):
        """Main processing pipeline for RETNLA ortholog cells"""
        
        logger.info("Starting RETNLA ortholog Geneformer processing...")
        
        # Step 1: Verify all data files exist
        logger.info("Step 1: Verifying data files...")
        file_info, total_cells = self.verify_data_files()
        
        # Step 2: Combine datasets with proper metadata
        logger.info("Step 2: Combining datasets...")
        combined_adata = self.prepare_combined_dataset(file_info)
        
        # Step 3: Prepare Geneformer format
        logger.info("Step 3: Preparing Geneformer format...")
        data_path, metadata_path = self.prepare_geneformer_format(combined_adata)
        
        # Step 4: Run tokenization
        logger.info("Step 4: Running TranscriptomeTokenizer...")
        output_dir = self.run_geneformer_tokenizer(data_path, metadata_path)
        
        # Step 5: Verify and organize output
        logger.info("Step 5: Verifying tokenization output...")
        tokenized_files = self.verify_tokenization_output(output_dir)
        
        logger.info("RETNLA ortholog Geneformer processing completed successfully!")
        logger.info(f"Tokenized files: {tokenized_files}")
        
        return tokenized_files

def main():
    """Main execution function"""
    
    try:
        processor = GeneformerProcessor()
        tokenized_files = processor.process_all_data()
        
        print("\n" + "="*50)
        print("GENEFORMER PROCESSING COMPLETE")
        print("="*50)
        print(f"Container: {CONTAINER_PATH}")
        print(f"Output directory: {OUTPUT_DIR}")
        print(f"Tokenized files: {len(tokenized_files)}")
        for f in tokenized_files:
            print(f"  - {f}")
        print("\nNext step: Generate embeddings using generate_geneformer_embeddings.py")
        print("="*50)
        
    except Exception as e:
        logger.error(f"Processing failed: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()