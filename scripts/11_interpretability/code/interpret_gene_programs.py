#!/usr/bin/env python3
"""
Interpret Gene Programs using Gene-Program Correlation

Since NMF is run on embeddings (not genes directly), we need to:
1. Load program scores (W matrix) for each cell
2. Load original gene expression data
3. Correlate gene expression with program scores
4. Extract top correlated genes per program

This gives us biologically interpretable gene lists for each program.

Usage:
    python interpret_gene_programs.py --model STATE
    python interpret_gene_programs.py --model UCE --top-genes 50

Author: Project 302006
Date: December 2025
"""

import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
from datetime import datetime
from scipy import stats
from scipy.sparse import issparse
import argparse
import logging
import json
import warnings
warnings.filterwarnings('ignore')

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Paths
BASE_DIR = Path("/bigdata/godziklab/shared/Xinru/302006")
DATA_DIR = BASE_DIR / "01_DATA"
EMB_DIR = BASE_DIR / "02_EMBEDDINGS"
RESULTS_DIR = BASE_DIR / "05_RESULTS" / "gene_programs"
OUTPUT_DIR = BASE_DIR / "05_RESULTS" / "gene_programs_interpretation"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_program_scores(model_name='STATE'):
    """Load NMF program loadings (W matrix) for each cell."""
    # Find the most recent loadings file (W matrix from NMF)
    score_files = list(RESULTS_DIR.glob(f"program_loadings_{model_name}_*.npy"))
    if not score_files:
        raise FileNotFoundError(f"No program loadings found for {model_name}")

    latest_file = max(score_files, key=lambda x: x.stat().st_mtime)
    logger.info(f"Loading program loadings from: {latest_file}")

    W = np.load(latest_file)  # Shape: (n_cells, n_programs)
    logger.info(f"  Program loadings shape: {W.shape} (cells × programs)")

    return W, latest_file


def load_disease_programs(model_name='STATE'):
    """Load disease-associated program rankings."""
    disease_files = list(RESULTS_DIR.glob(f"disease_programs_{model_name}_*.csv"))
    if not disease_files:
        return None

    latest_file = max(disease_files, key=lambda x: x.stat().st_mtime)
    df = pd.read_csv(latest_file)
    logger.info(f"Loaded disease programs from: {latest_file}")

    return df


def load_expression_data(model_name='STATE'):
    """Load original gene expression data.

    We need the raw/normalized counts to correlate with program scores.
    """
    # Try COVID data first (main analysis dataset)
    if model_name == 'STATE':
        # STATE uses the original COVID h5ad
        data_path = DATA_DIR / "covid" / "human_platelet_covid_severity.h5ad"
    else:
        # UCE uses its own processed file
        data_path = EMB_DIR / "uce" / "human_platelet_covid_severity_uce_adata.h5ad"

    if not data_path.exists():
        # Try alternative paths
        alt_paths = [
            DATA_DIR / "human_platelet_covid_severity.h5ad",
            EMB_DIR / "state" / "covid_state.h5ad",
        ]
        for alt in alt_paths:
            if alt.exists():
                data_path = alt
                break
        else:
            raise FileNotFoundError(f"No expression data found for {model_name}")

    logger.info(f"Loading expression data from: {data_path}")
    adata = sc.read_h5ad(data_path)
    logger.info(f"  Expression shape: {adata.shape} (cells × genes)")
    logger.info(f"  Gene names: {adata.var_names[:5].tolist()}...")

    return adata


def correlate_genes_with_programs(W, adata, n_top=30, batch_size=1000):
    """Correlate gene expression with program scores.

    For each program, compute Spearman correlation with each gene.
    Return top correlated genes per program.

    Args:
        W: Program scores matrix (n_cells, n_programs)
        adata: AnnData with expression
        n_top: Number of top genes per program
        batch_size: Process genes in batches for memory efficiency
    """
    n_cells, n_programs = W.shape
    n_genes = adata.n_vars
    gene_names = adata.var_names.tolist()

    # Verify cell count matches
    if n_cells != adata.n_obs:
        logger.warning(f"Cell count mismatch: W has {n_cells}, adata has {adata.n_obs}")
        # Use minimum
        n_cells = min(n_cells, adata.n_obs)
        W = W[:n_cells]

    logger.info(f"Computing correlations: {n_programs} programs × {n_genes} genes")

    # Get expression matrix
    X = adata.X[:n_cells]
    if issparse(X):
        X = X.toarray()

    # Store correlations for each program
    program_genes = {}

    for prog_idx in range(n_programs):
        prog_scores = W[:, prog_idx]

        # Compute correlations in batches
        correlations = np.zeros(n_genes)
        pvalues = np.zeros(n_genes)

        for start in range(0, n_genes, batch_size):
            end = min(start + batch_size, n_genes)
            gene_batch = X[:, start:end]

            for i, gene_idx in enumerate(range(start, end)):
                gene_expr = gene_batch[:, i]
                # Skip genes with no variance
                if np.std(gene_expr) == 0:
                    correlations[gene_idx] = 0
                    pvalues[gene_idx] = 1.0
                else:
                    corr, pval = stats.spearmanr(prog_scores, gene_expr)
                    correlations[gene_idx] = corr if not np.isnan(corr) else 0
                    pvalues[gene_idx] = pval if not np.isnan(pval) else 1.0

        # Get top positively correlated genes
        top_pos_idx = np.argsort(correlations)[::-1][:n_top]
        # Get top negatively correlated genes
        top_neg_idx = np.argsort(correlations)[:n_top]

        program_genes[prog_idx] = {
            'top_positive': {
                'genes': [gene_names[i] for i in top_pos_idx],
                'correlations': correlations[top_pos_idx].tolist(),
                'pvalues': pvalues[top_pos_idx].tolist()
            },
            'top_negative': {
                'genes': [gene_names[i] for i in top_neg_idx],
                'correlations': correlations[top_neg_idx].tolist(),
                'pvalues': pvalues[top_neg_idx].tolist()
            }
        }

        if prog_idx % 10 == 0:
            logger.info(f"  Processed program {prog_idx}/{n_programs}")
            top_gene = gene_names[top_pos_idx[0]]
            top_corr = correlations[top_pos_idx[0]]
            logger.info(f"    Top gene: {top_gene} (r={top_corr:.3f})")

    return program_genes


def prepare_geneagent_input(program_genes, disease_programs, model_name, n_top=30):
    """Prepare gene sets for GeneAgent analysis."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Prioritize disease-associated programs
    if disease_programs is not None:
        priority_order = disease_programs.sort_values('abs_coef', ascending=False)['program'].tolist()
    else:
        priority_order = list(program_genes.keys())

    # Prepare output
    geneagent_input = {
        'model': model_name,
        'timestamp': timestamp,
        'n_programs': len(program_genes),
        'n_top_genes': n_top,
        'method': 'spearman_correlation',
        'programs': {}
    }

    # Text file for batch processing
    gene_lists_file = OUTPUT_DIR / f"gene_lists_{model_name}_{timestamp}.txt"
    gene_json_file = OUTPUT_DIR / f"geneagent_input_{model_name}_{timestamp}.json"

    with open(gene_lists_file, 'w') as f:
        for rank, prog_id in enumerate(priority_order):
            if prog_id not in program_genes:
                continue

            top_pos = program_genes[prog_id]['top_positive']
            top_neg = program_genes[prog_id]['top_negative']

            # Get disease direction if available
            direction = "UNKNOWN"
            coef = 0.0
            if disease_programs is not None:
                prog_row = disease_programs[disease_programs['program'] == prog_id]
                if len(prog_row) > 0:
                    direction = prog_row['direction'].values[0]
                    coef = prog_row['coefficient'].values[0]

            geneagent_input['programs'][str(prog_id)] = {
                'rank': rank + 1,
                'direction': direction,
                'coefficient': float(coef),
                'top_positive_genes': top_pos['genes'],
                'top_positive_correlations': top_pos['correlations'],
                'top_negative_genes': top_neg['genes'],
                'top_negative_correlations': top_neg['correlations']
            }

            # Write gene list (positively correlated genes for pathway analysis)
            f.write(f"# Program {prog_id} ({direction}, coef={coef:.4f})\n")
            f.write(f"# Top positively correlated genes:\n")
            f.write(",".join(top_pos['genes'][:n_top]) + "\n")
            f.write(f"# Top negatively correlated genes:\n")
            f.write(",".join(top_neg['genes'][:n_top]) + "\n\n")

    # Save JSON
    with open(gene_json_file, 'w') as f:
        json.dump(geneagent_input, f, indent=2)

    logger.info(f"Saved gene lists to: {gene_lists_file}")
    logger.info(f"Saved JSON input to: {gene_json_file}")

    return gene_lists_file, gene_json_file


def run_enrichment_analysis(program_genes, disease_programs, model_name, n_top=30):
    """Run traditional enrichment analysis using gseapy."""
    try:
        import gseapy as gp
    except ImportError:
        logger.warning("gseapy not installed. Skipping enrichment analysis.")
        logger.info("Install with: pip install gseapy")
        return None

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    enrichment_results = {}

    # Prioritize disease-associated programs
    if disease_programs is not None:
        priority_order = disease_programs.sort_values('abs_coef', ascending=False)['program'].tolist()[:20]
    else:
        priority_order = list(program_genes.keys())[:20]

    logger.info(f"Running GO/KEGG enrichment for top {len(priority_order)} programs...")

    for prog_id in priority_order:
        if prog_id not in program_genes:
            continue

        # Use positively correlated genes
        genes = program_genes[prog_id]['top_positive']['genes'][:n_top]

        try:
            # Run enrichment
            enr = gp.enrichr(
                gene_list=genes,
                gene_sets=['GO_Biological_Process_2023', 'KEGG_2021_Human'],
                organism='human',
                outdir=None,
                cutoff=0.05
            )

            if enr.results is not None and len(enr.results) > 0:
                top_terms = enr.results.head(10)[['Term', 'Adjusted P-value', 'Genes']].to_dict('records')
                enrichment_results[str(prog_id)] = {
                    'top_genes': genes[:10],
                    'enrichment': top_terms
                }
                logger.info(f"  Program {prog_id}: {enr.results.iloc[0]['Term'][:60]}...")

        except Exception as e:
            logger.warning(f"  Program {prog_id}: Enrichment failed - {e}")
            continue

    # Save results
    if enrichment_results:
        output_file = OUTPUT_DIR / f"enrichment_results_{model_name}_{timestamp}.json"
        with open(output_file, 'w') as f:
            json.dump(enrichment_results, f, indent=2)
        logger.info(f"Saved enrichment results to: {output_file}")

    return enrichment_results


def print_summary(program_genes, disease_programs, model_name):
    """Print summary of top disease-associated programs."""

    print("\n" + "=" * 80)
    print(f"TOP DISEASE-ASSOCIATED GENE PROGRAMS ({model_name})")
    print("=" * 80)

    if disease_programs is None:
        print("No disease programs available.")
        return

    # Top 5 SEVERE and NON-SEVERE programs
    severe = disease_programs[disease_programs['direction'] == 'SEVERE'].head(5)
    nonsevere = disease_programs[disease_programs['direction'] == 'NON-SEVERE'].head(5)

    print("\n🔴 TOP SEVERE-ASSOCIATED PROGRAMS:")
    print("-" * 60)
    for _, row in severe.iterrows():
        prog_id = row['program']
        coef = row['coefficient']
        if prog_id in program_genes:
            top_genes = program_genes[prog_id]['top_positive']['genes'][:5]
            print(f"  Program {prog_id} (coef={coef:+.3f}): {', '.join(top_genes)}")

    print("\n🟢 TOP NON-SEVERE-ASSOCIATED PROGRAMS:")
    print("-" * 60)
    for _, row in nonsevere.iterrows():
        prog_id = row['program']
        coef = row['coefficient']
        if prog_id in program_genes:
            top_genes = program_genes[prog_id]['top_positive']['genes'][:5]
            print(f"  Program {prog_id} (coef={coef:+.3f}): {', '.join(top_genes)}")

    print("\n" + "=" * 80)


def main():
    parser = argparse.ArgumentParser(description='Interpret gene programs via correlation')
    parser.add_argument('--model', type=str, default='STATE', choices=['STATE', 'UCE'],
                        help='Foundation model')
    parser.add_argument('--top-genes', type=int, default=30,
                        help='Number of top genes per program')
    parser.add_argument('--run-enrichment', action='store_true',
                        help='Run gseapy enrichment analysis')
    args = parser.parse_args()

    logger.info("=" * 80)
    logger.info(f"GENE PROGRAM INTERPRETATION: {args.model}")
    logger.info("Method: Gene-Program Spearman Correlation")
    logger.info("=" * 80)

    # Load program scores (W matrix)
    W, score_file = load_program_scores(args.model)

    # Load disease program rankings
    disease_programs = load_disease_programs(args.model)

    # Load expression data
    adata = load_expression_data(args.model)

    # Correlate genes with programs
    logger.info("\nComputing gene-program correlations...")
    program_genes = correlate_genes_with_programs(
        W, adata, n_top=args.top_genes
    )

    # Prepare GeneAgent input
    gene_lists_file, gene_json_file = prepare_geneagent_input(
        program_genes, disease_programs, args.model, n_top=args.top_genes
    )

    # Run enrichment analysis if requested
    if args.run_enrichment:
        run_enrichment_analysis(program_genes, disease_programs, args.model, n_top=args.top_genes)

    # Print summary
    print_summary(program_genes, disease_programs, args.model)

    logger.info("\n" + "=" * 80)
    logger.info("DONE")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
