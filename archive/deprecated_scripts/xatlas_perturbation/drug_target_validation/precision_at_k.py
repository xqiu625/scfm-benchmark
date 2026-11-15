#!/usr/bin/env python3
"""
ACC 20: Precision@K and Recall@K Metrics

Statistical evaluation of FM-nominated drug targets against ground truth databases.

Functions:
    - precision_at_k: fraction of top-k predictions that are known targets
    - recall_at_k: fraction of known targets found in top-k predictions
    - fisher_exact_enrichment: Fisher's exact test for target enrichment
    - average_precision: area under precision-recall curve

Author: Xinru Qiu
Date: February 2026
"""

import numpy as np
from scipy import stats


def precision_at_k(ranked_genes, known_targets, k):
    """
    Compute precision@k: fraction of top-k that are known targets.

    Args:
        ranked_genes: list of gene names ranked by perturbation score
        known_targets: set of known drug target gene names
        k: number of top genes to evaluate

    Returns:
        precision: float in [0, 1]
    """
    top_k = set(g.upper() for g in ranked_genes[:k])
    known = set(g.upper() for g in known_targets)
    hits = top_k & known
    return len(hits) / k if k > 0 else 0.0


def recall_at_k(ranked_genes, known_targets, k):
    """
    Compute recall@k: fraction of known targets found in top-k.

    Args:
        ranked_genes: list of gene names ranked by perturbation score
        known_targets: set of known drug target gene names
        k: number of top genes to evaluate

    Returns:
        recall: float in [0, 1]
    """
    top_k = set(g.upper() for g in ranked_genes[:k])
    known = set(g.upper() for g in known_targets)
    if len(known) == 0:
        return 0.0
    hits = top_k & known
    return len(hits) / len(known)


def average_precision(ranked_genes, known_targets, max_k=None):
    """
    Compute average precision (AP) over the full ranking.

    Args:
        ranked_genes: list of gene names ranked by perturbation score
        known_targets: set of known drug target gene names
        max_k: optional maximum rank to evaluate (default: all)

    Returns:
        ap: float, average precision score
    """
    known = set(g.upper() for g in known_targets)
    if len(known) == 0:
        return 0.0

    if max_k is None:
        max_k = len(ranked_genes)

    n_relevant = 0
    sum_precision = 0.0

    for i, gene in enumerate(ranked_genes[:max_k]):
        if gene.upper() in known:
            n_relevant += 1
            sum_precision += n_relevant / (i + 1)

    return sum_precision / len(known)


def fisher_exact_enrichment(ranked_genes, known_targets, k, total_genes):
    """
    Fisher's exact test for enrichment of known targets in top-k.

    Tests whether top-k ranked genes contain more known targets than
    expected by chance.

    Args:
        ranked_genes: list of gene names ranked by perturbation score
        known_targets: set of known drug target gene names
        k: number of top genes to evaluate
        total_genes: total number of genes in the ranking

    Returns:
        dict with odds_ratio, p_value, and contingency table
    """
    top_k = set(g.upper() for g in ranked_genes[:k])
    known = set(g.upper() for g in known_targets)
    all_genes = set(g.upper() for g in ranked_genes)

    # Contingency table
    #                 Known Target    Not Target
    # In Top-k            a              b
    # Not in Top-k        c              d

    a = len(top_k & known)  # known targets in top-k
    b = k - a               # non-targets in top-k
    c = len(known) - a      # known targets not in top-k
    d = total_genes - k - c  # non-targets not in top-k

    # Ensure non-negative values
    a = max(0, a)
    b = max(0, b)
    c = max(0, c)
    d = max(0, d)

    table = np.array([[a, b], [c, d]])
    odds_ratio, p_value = stats.fisher_exact(table, alternative='greater')

    return {
        'odds_ratio': float(odds_ratio),
        'p_value': float(p_value),
        'hits_in_top_k': a,
        'total_known_targets': len(known),
        'k': k,
        'total_genes': total_genes,
        'contingency_table': table.tolist(),
    }


def random_baseline_precision(n_known_targets, total_genes, k, n_simulations=10000):
    """
    Compute expected precision@k under random ranking.

    Args:
        n_known_targets: number of known targets in gene set
        total_genes: total number of genes
        k: number of top genes
        n_simulations: number of random permutations

    Returns:
        dict with mean, std, and percentiles of random precision@k
    """
    precisions = []
    for _ in range(n_simulations):
        random_ranking = np.random.permutation(total_genes)
        # Count how many of the first k are among the n_known_targets
        hits = np.sum(random_ranking[:k] < n_known_targets)
        precisions.append(hits / k)

    return {
        'mean': float(np.mean(precisions)),
        'std': float(np.std(precisions)),
        'p5': float(np.percentile(precisions, 5)),
        'p95': float(np.percentile(precisions, 95)),
    }


def evaluate_ranking(ranked_genes, known_targets, ks=[5, 10, 20, 50, 100],
                     total_gene_universe=None):
    """
    Comprehensive evaluation of a gene ranking.

    Args:
        ranked_genes: list of gene names
        known_targets: set of known target gene names
        ks: list of k values to evaluate
        total_gene_universe: total number of genes in the full perturbation
            library (for Fisher's exact test and random baseline).
            If None, defaults to len(ranked_genes).

    Returns:
        dict with all metrics at each k
    """
    total_ranked = len(ranked_genes)
    total_genes = total_gene_universe if total_gene_universe else total_ranked
    results = {
        'total_ranked': total_ranked,
        'total_gene_universe': total_genes,
        'total_known_targets': len(known_targets),
        'metrics_by_k': {},
    }

    for k in ks:
        if k > total_ranked:
            continue

        prec = precision_at_k(ranked_genes, known_targets, k)
        rec = recall_at_k(ranked_genes, known_targets, k)
        fisher = fisher_exact_enrichment(ranked_genes, known_targets, k, total_genes)
        random = random_baseline_precision(len(known_targets), total_genes, k, n_simulations=1000)

        results['metrics_by_k'][k] = {
            'precision': prec,
            'recall': rec,
            'fisher_odds_ratio': fisher['odds_ratio'],
            'fisher_p_value': fisher['p_value'],
            'hits': fisher['hits_in_top_k'],
            'random_baseline_mean': random['mean'],
            'random_baseline_std': random['std'],
            'fold_enrichment': prec / random['mean'] if random['mean'] > 0 else float('inf'),
        }

    results['average_precision'] = average_precision(ranked_genes, known_targets)

    return results
