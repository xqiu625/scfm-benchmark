#!/usr/bin/env python3
"""
ACC 20: Drug Target Validation Pipeline

Validate FM-nominated therapeutic targets from ACC 11 perturbation analysis
against known drug databases (ChEMBL, OpenTargets, DrugBank).

Pipeline:
1. Load perturbation recovery rankings from ACC 11 (STATE + UCE)
2. Build ground truth from drug databases
3. Compute precision@k and recall@k
4. Fisher's exact test for enrichment
5. Literature validation for top candidates

Usage:
    python validate_therapeutic_targets.py
    python validate_therapeutic_targets.py --skip-api  # Use cached results only

Author: Xinru Qiu
Date: February 2026
"""

import argparse
import logging
from pathlib import Path
import json
import time
import numpy as np
import pandas as pd
import requests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# Import precision_at_k metrics
import sys
sys.path.insert(0, str(Path(__file__).parent))
from precision_at_k import evaluate_ranking, precision_at_k, fisher_exact_enrichment

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

BASE_DIR = Path("/bigdata/godziklab/shared/Xinru/302006")
RESULTS_DIR = BASE_DIR / "05_RESULTS/20_drug_target_validation"
FIGURES_DIR = RESULTS_DIR / "figures"
CACHE_DIR = RESULTS_DIR / "cache"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
CACHE_DIR.mkdir(parents=True, exist_ok=True)


# =============================================================================
# Ground Truth: Known Drug Targets
# =============================================================================

# Curated anti-inflammatory / anti-COVID / anti-thrombosis drug targets
KNOWN_TARGETS = {
    # Anti-inflammatory targets
    'IL6', 'IL1B', 'TNF', 'JAK1', 'JAK2', 'JAK3', 'TYK2',
    'STAT3', 'STAT1', 'NFKB1', 'RELA', 'PTGS2', 'PTGS1',
    'PDE4B', 'PDE4D', 'ALOX5', 'LTA4H', 'MAPK14', 'MAPK1',
    'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'MTOR', 'AKT1',
    'BTK', 'SYK', 'IRAK4', 'TLR4', 'MYD88', 'NLRP3',

    # Anti-COVID targets (repurposed)
    'ACE2', 'TMPRSS2', 'CTSL', 'FURIN', 'BSG', 'NRP1',
    'PLpro', 'RDRP', 'CASP1', 'CASP3',

    # Anti-thrombotic targets (relevant to platelet biology)
    'ITGA2B', 'ITGB3', 'GP1BA', 'GP1BB', 'GP5', 'GP9',
    'P2RY12', 'TBXA2R', 'F2R', 'F2RL3',
    'PTGIR', 'PLA2G4A', 'TBXAS1',
    'SELP', 'SELPLG', 'ITGAL', 'ITGB2',

    # Kinase targets (broad anti-inflammatory)
    'ABL1', 'SRC', 'FYN', 'LCK', 'LYN',
    'EGFR', 'ERBB2', 'KIT', 'PDGFRA', 'FGFR1',
    'VEGFA', 'KDR', 'FLT1',
}


def search_chembl_target(gene_name, cache_dir=None):
    """Search ChEMBL for target info. Caches results."""
    if cache_dir:
        cache_file = cache_dir / f"chembl_{gene_name}.json"
        if cache_file.exists():
            with open(cache_file) as f:
                return json.load(f)

    result = {
        'gene': gene_name,
        'chembl_target_id': None,
        'target_type': None,
        'n_compounds': 0,
        'n_approved_drugs': 0,
        'has_activity_data': False,
    }

    try:
        url = f"https://www.ebi.ac.uk/chembl/api/data/target/search.json"
        params = {'q': gene_name, 'limit': 5}
        resp = requests.get(url, params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()

        if data.get('targets'):
            for target in data['targets']:
                if 'Homo sapiens' in str(target.get('organism', '')):
                    result['chembl_target_id'] = target.get('target_chembl_id')
                    result['target_type'] = target.get('target_type')
                    result['pref_name'] = target.get('pref_name')
                    break

        if result['chembl_target_id']:
            # Count activities
            act_url = f"https://www.ebi.ac.uk/chembl/api/data/activity.json"
            act_params = {
                'target_chembl_id': result['chembl_target_id'],
                'limit': 1
            }
            act_resp = requests.get(act_url, params=act_params, timeout=30)
            if act_resp.status_code == 200:
                act_data = act_resp.json()
                result['n_compounds'] = act_data.get('page_meta', {}).get('total_count', 0)
                result['has_activity_data'] = result['n_compounds'] > 0

        time.sleep(0.5)  # Rate limiting

    except Exception as e:
        result['error'] = str(e)

    if cache_dir:
        with open(cache_dir / f"chembl_{gene_name}.json", 'w') as f:
            json.dump(result, f, indent=2)

    return result


def search_opentargets(gene_name, cache_dir=None):
    """Search OpenTargets for disease associations."""
    if cache_dir:
        cache_file = cache_dir / f"opentargets_{gene_name}.json"
        if cache_file.exists():
            with open(cache_file) as f:
                return json.load(f)

    result = {
        'gene': gene_name,
        'covid_association': False,
        'inflammation_association': False,
        'thrombosis_association': False,
        'overall_score': 0.0,
    }

    try:
        # OpenTargets GraphQL API
        url = "https://api.platform.opentargets.org/api/v4/graphql"
        query = """
        query {
            search(queryString: "%s", entityNames: ["target"], page: {size: 1, index: 0}) {
                hits {
                    id
                    name
                    entity
                }
            }
        }
        """ % gene_name

        resp = requests.post(url, json={'query': query}, timeout=30)
        if resp.status_code == 200:
            data = resp.json()
            hits = data.get('data', {}).get('search', {}).get('hits', [])
            if hits:
                result['opentargets_id'] = hits[0].get('id')
                result['found'] = True

        time.sleep(0.3)

    except Exception as e:
        result['error'] = str(e)

    if cache_dir:
        with open(cache_dir / f"opentargets_{gene_name}.json", 'w') as f:
            json.dump(result, f, indent=2)

    return result


def extract_gene_name(perturbation):
    """Extract clean gene name from perturbation string.

    ACC 11 perturbation names have dual-guide suffixes:
      ICMT_P1P2-1|ICMT_P1P2-2  -> ICMT
      SLC28A1_P1-1|SLC28A1_P1-2 -> SLC28A1
      KRTAP21-3_ENST00000444335.1-1|... -> KRTAP21-3
      non-targeting_02896|... -> non-targeting
    """
    if pd.isna(perturbation):
        return None
    gene = str(perturbation).split('_')[0]
    if gene.startswith('non-targeting'):
        return 'non-targeting'
    return gene.upper()


def load_perturbation_rankings():
    """Load perturbation recovery rankings from ACC 11.

    Returns:
        rankings: dict mapping model name -> list of clean gene symbols (ranked)
        total_gene_universe: int, total unique genes across all perturbation files
    """
    import glob
    rankings = {}
    all_unique_genes = set()

    pert_dir = BASE_DIR / "05_RESULTS/perturbation_recovery"

    for model in ['STATE', 'UCE']:
        # Prefer all_perturbations (full ranking ~2500 genes) over therapeutic_candidates (top 100)
        glob_patterns = [
            str(pert_dir / f"{model.lower()}_all_perturbations_*.csv"),
            str(pert_dir / f"{model.lower()}_therapeutic_candidates_*.csv"),
            str(pert_dir / f"{model.lower()}_recovery_rankings*.csv"),
            str(pert_dir / f"recovery_rankings_{model}*.csv"),
            str(pert_dir / f"therapeutic_candidates_{model.lower()}*.csv"),
        ]

        for pattern in glob_patterns:
            matches = sorted(glob.glob(pattern))
            if matches:
                # Use the most recent file (last in sorted order)
                f = Path(matches[-1])
                df = pd.read_csv(f)
                gene_col = next((c for c in df.columns if 'gene' in c.lower()), df.columns[0])
                raw_names = df[gene_col].tolist()

                # Clean perturbation guide suffixes to get base gene symbols
                clean_names = []
                seen = set()
                for name in raw_names:
                    gene = extract_gene_name(name)
                    if gene and gene != 'non-targeting' and gene not in seen:
                        clean_names.append(gene)
                        seen.add(gene)

                rankings[model] = clean_names
                all_unique_genes.update(clean_names)
                logger.info(f"Loaded {model} rankings: {len(clean_names)} unique genes "
                           f"(from {len(raw_names)} perturbations) from {f.name}")
                logger.info(f"  Top-5: {clean_names[:5]}")
                break

    if not rankings:
        logger.warning("No perturbation rankings found. Using default candidate list.")
        default_candidates = [
            'ICMT', 'DUSP11', 'LILRB2', 'B4GALT1', 'PINK1',
            'PIK3R2', 'RAMP3', 'SLC28A1', 'ATP2B3', 'ESD',
            'CST7', 'GZMA', 'PRF1', 'NKG7', 'CCL5',
            'IFITM1', 'ISG15', 'IFI6', 'MX1', 'OAS1',
        ]
        rankings['default'] = default_candidates
        all_unique_genes.update(default_candidates)

    total_universe = len(all_unique_genes) if all_unique_genes else 2907  # fallback: X-Atlas has ~2907 targets
    logger.info(f"Total gene universe: {total_universe} unique genes")

    return rankings, total_universe


def validate_targets(rankings, total_gene_universe=2907, skip_api=False):
    """Run full validation pipeline."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    logger.info("=" * 60)
    logger.info("ACC 20: Drug Target Validation")
    logger.info(f"Known targets: {len(KNOWN_TARGETS)} | Gene universe: {total_gene_universe}")
    logger.info("=" * 60)

    all_results = {}

    # Diagnostic: check overlap between KNOWN_TARGETS and the ranked gene lists
    for model_name, ranked_genes in rankings.items():
        ranked_set = set(g.upper() for g in ranked_genes)
        known_upper = set(g.upper() for g in KNOWN_TARGETS)
        testable = known_upper & ranked_set
        logger.info(f"\n  Diagnostic ({model_name}): {len(testable)}/{len(known_upper)} "
                   f"known targets are in the {len(ranked_set)}-gene ranking")
        if testable:
            logger.info(f"    Testable known targets: {sorted(testable)}")
        else:
            logger.info(f"    NO known targets found in ranking — enrichment test is vacuous")
            # Check which known targets exist in the full perturbation universe
            # by looking at the all_perturbations files
            import glob
            pert_dir = BASE_DIR / "05_RESULTS/perturbation_recovery"
            for pattern in [str(pert_dir / f"{model_name.lower()}_all_perturbations_*.csv")]:
                matches = sorted(glob.glob(pattern))
                if matches:
                    df_full = pd.read_csv(Path(matches[-1]))
                    col = next((c for c in df_full.columns if 'gene' in c.lower()), df_full.columns[0])
                    full_genes = set()
                    for name in df_full[col]:
                        gene = extract_gene_name(name)
                        if gene and gene != 'non-targeting':
                            full_genes.add(gene)
                    testable_full = known_upper & full_genes
                    logger.info(f"    In FULL perturbation library ({len(full_genes)} genes): "
                               f"{len(testable_full)} known targets present")
                    if testable_full:
                        logger.info(f"    Present known targets: {sorted(testable_full)}")
                    break

    for model_name, ranked_genes in rankings.items():
        logger.info(f"\n--- {model_name} Rankings ({len(ranked_genes)} genes) ---")

        # Evaluate against known targets
        # Use K values up to the ranking size
        ks = [k for k in [5, 10, 20, 50, 100, 200, 500, 1000] if k <= len(ranked_genes)]
        eval_results = evaluate_ranking(
            ranked_genes, KNOWN_TARGETS, ks=ks,
            total_gene_universe=total_gene_universe
        )

        logger.info(f"\n  Precision@K results:")
        for k, metrics in eval_results['metrics_by_k'].items():
            sig_str = "***" if metrics['fisher_p_value'] < 0.001 else \
                      "**" if metrics['fisher_p_value'] < 0.01 else \
                      "*" if metrics['fisher_p_value'] < 0.05 else ""
            logger.info(
                f"    P@{k:3d}: {metrics['precision']:.3f} "
                f"(hits={metrics['hits']}, "
                f"random={metrics['random_baseline_mean']:.3f}, "
                f"fold={metrics['fold_enrichment']:.1f}x, "
                f"Fisher p={metrics['fisher_p_value']:.2e}) {sig_str}"
            )
        logger.info(f"  Average Precision: {eval_results['average_precision']:.4f}")

        # Validate top-10 against databases
        top_10 = ranked_genes[:10]
        logger.info(f"\n  Top-10 candidate validation:")

        candidate_info = []
        for gene in top_10:
            gene_upper = gene.upper()
            in_known = gene_upper in KNOWN_TARGETS

            chembl_info = {'gene': gene_upper}
            ot_info = {'gene': gene_upper}

            if not skip_api:
                chembl_info = search_chembl_target(gene_upper, CACHE_DIR)
                ot_info = search_opentargets(gene_upper, CACHE_DIR)

            # Classify validation status
            if in_known:
                status = "validated"
            elif chembl_info.get('has_activity_data', False):
                status = "plausible"
            else:
                status = "novel"

            candidate_info.append({
                'gene': gene_upper,
                'rank': top_10.index(gene) + 1,
                'validation_status': status,
                'in_known_targets': in_known,
                'chembl_id': chembl_info.get('chembl_target_id'),
                'chembl_compounds': chembl_info.get('n_compounds', 0),
            })

            status_icon = {"validated": "[V]", "plausible": "[P]", "novel": "[N]"}[status]
            logger.info(f"    {status_icon} {gene_upper}: {status} "
                       f"(ChEMBL compounds: {chembl_info.get('n_compounds', 'N/A')})")

        all_results[model_name] = {
            'evaluation': eval_results,
            'top_10_candidates': candidate_info,
            'n_validated': sum(1 for c in candidate_info if c['validation_status'] == 'validated'),
            'n_plausible': sum(1 for c in candidate_info if c['validation_status'] == 'plausible'),
            'n_novel': sum(1 for c in candidate_info if c['validation_status'] == 'novel'),
        }

    # Save results
    with open(RESULTS_DIR / f"validation_results_{timestamp}.json", 'w') as f:
        json.dump(all_results, f, indent=2, default=str)

    # Create summary CSV
    rows = []
    for model, res in all_results.items():
        for k, metrics in res['evaluation']['metrics_by_k'].items():
            rows.append({
                'model': model,
                'k': k,
                'precision': metrics['precision'],
                'recall': metrics['recall'],
                'fisher_p_value': metrics['fisher_p_value'],
                'fold_enrichment': metrics['fold_enrichment'],
                'hits': metrics['hits'],
            })
    pd.DataFrame(rows).to_csv(RESULTS_DIR / f"precision_at_k_{timestamp}.csv", index=False)

    # Plot precision@K
    plot_precision_at_k(all_results)

    # Summary
    logger.info(f"\n{'='*60}")
    logger.info("VALIDATION SUMMARY")
    logger.info(f"{'='*60}")
    for model, res in all_results.items():
        logger.info(f"  {model}: {res['n_validated']} validated, "
                    f"{res['n_plausible']} plausible, {res['n_novel']} novel (top-10)")

    return all_results


def plot_precision_at_k(all_results):
    """Plot precision@K curves for all models."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    colors = {'STATE': '#e74c3c', 'UCE': '#3498db', 'default': '#95a5a6'}

    # Get K values from actual results
    first_res = next(iter(all_results.values()))
    ks = sorted(first_res['evaluation']['metrics_by_k'].keys())

    # Precision@K
    ax = axes[0]
    for model, res in all_results.items():
        valid_ks = [k for k in ks if k in res['evaluation']['metrics_by_k']]
        valid_precs = [res['evaluation']['metrics_by_k'][k]['precision'] for k in valid_ks]
        ax.plot(valid_ks, valid_precs, 'o-', label=model, color=colors.get(model, '#333'), linewidth=2)

    # Random baseline
    if all_results:
        first = next(iter(all_results.values()))
        random_precs = [first['evaluation']['metrics_by_k'].get(k, {}).get('random_baseline_mean', 0) for k in ks]
        valid_ks = [k for k in ks if k in first['evaluation']['metrics_by_k']]
        valid_random = [first['evaluation']['metrics_by_k'][k]['random_baseline_mean'] for k in valid_ks]
        ax.plot(valid_ks, valid_random, '--', color='gray', label='Random', linewidth=1)

    ax.set_xlabel('K')
    ax.set_ylabel('Precision@K')
    ax.set_title('Precision@K: FM Targets vs Known Drugs')
    ax.legend()
    ax.set_xscale('log')

    # Validation status pie chart
    ax = axes[1]
    if all_results:
        first_model = next(iter(all_results.values()))
        sizes = [first_model['n_validated'], first_model['n_plausible'], first_model['n_novel']]
        labels = ['Validated', 'Plausible', 'Novel']
        colors_pie = ['#2ecc71', '#f39c12', '#e74c3c']
        ax.pie(sizes, labels=labels, colors=colors_pie, autopct='%1.0f%%',
              startangle=90, textprops={'fontsize': 12})
        ax.set_title('Top-10 Candidate Validation Status')

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'precision_at_k.png', dpi=150, bbox_inches='tight')
    plt.savefig(FIGURES_DIR / 'precision_at_k.pdf', bbox_inches='tight')
    plt.close()
    logger.info(f"Saved figures to {FIGURES_DIR}")


def main():
    parser = argparse.ArgumentParser(description='Drug Target Validation')
    parser.add_argument('--skip-api', action='store_true', help='Skip API calls, use cached results')
    args = parser.parse_args()

    rankings, total_gene_universe = load_perturbation_rankings()
    validate_targets(rankings, total_gene_universe=total_gene_universe, skip_api=args.skip_api)


if __name__ == '__main__':
    main()
