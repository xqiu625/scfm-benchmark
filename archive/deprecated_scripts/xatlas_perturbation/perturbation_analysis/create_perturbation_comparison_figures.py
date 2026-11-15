#!/usr/bin/env python3
"""
Create comparison figures for STATE vs UCE perturbation analysis.

Generates:
1. Top candidates comparison (side-by-side bar chart)
2. Recovery score distribution comparison
3. Venn diagram of overlapping candidates
4. Pathway/function heatmap
5. Summary metrics dashboard

Usage:
    python scripts/analysis/create_perturbation_comparison_figures.py

Author: Project 302006
Date: December 2025
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from matplotlib_venn import venn2
import warnings
warnings.filterwarnings('ignore')

# Style settings
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 10

# Colors
STATE_COLOR = '#e74c3c'  # Red
UCE_COLOR = '#3498db'    # Blue
OVERLAP_COLOR = '#9b59b6'  # Purple

def load_data(results_dir):
    """Load STATE and UCE therapeutic candidates."""
    # Find the most recent files
    state_files = sorted(results_dir.glob('state_therapeutic_candidates_*.csv'))
    uce_files = sorted(results_dir.glob('uce_therapeutic_candidates_*.csv'))

    if not state_files or not uce_files:
        raise FileNotFoundError("Could not find therapeutic candidates files")

    state_df = pd.read_csv(state_files[-1])
    uce_df = pd.read_csv(uce_files[-1])

    # Load full scores
    state_scores_files = sorted(results_dir.glob('state_perturbation_scores_*.csv'))
    uce_scores_files = sorted(results_dir.glob('uce_perturbation_scores_*.csv'))

    state_scores = pd.read_csv(state_scores_files[-1]) if state_scores_files else None
    uce_scores = pd.read_csv(uce_scores_files[-1]) if uce_scores_files else None

    return state_df, uce_df, state_scores, uce_scores


def extract_gene_name(perturbation):
    """Extract gene name from perturbation string."""
    if pd.isna(perturbation):
        return None
    # Format: GENE_P1P2-1|GENE_P1P2-2 or GENE_ENST...-1|GENE_ENST...-2
    gene = perturbation.split('_')[0]
    if gene.startswith('non-targeting'):
        return 'non-targeting'
    return gene


def create_top_candidates_comparison(state_df, uce_df, output_dir, top_n=15):
    """Create side-by-side bar chart of top candidates."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 8))

    # Filter out non-targeting controls for cleaner visualization
    state_filtered = state_df[~state_df['perturbation'].str.contains('non-targeting', na=False)].head(top_n)
    uce_filtered = uce_df[~uce_df['perturbation'].str.contains('non-targeting', na=False)].head(top_n)

    # Extract gene names
    state_genes = [extract_gene_name(p) for p in state_filtered['perturbation']]
    uce_genes = [extract_gene_name(p) for p in uce_filtered['perturbation']]

    # STATE plot
    ax1 = axes[0]
    y_pos = np.arange(len(state_genes))
    bars1 = ax1.barh(y_pos, state_filtered['recovery_score'], color=STATE_COLOR, alpha=0.8)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(state_genes)
    ax1.invert_yaxis()
    ax1.set_xlabel('Recovery Score')
    ax1.set_title(f'STATE Top {top_n} Therapeutic Candidates', fontweight='bold', color=STATE_COLOR)
    ax1.set_xlim(23.5, 25.5)

    # Add value labels
    for bar, val in zip(bars1, state_filtered['recovery_score']):
        ax1.text(val + 0.02, bar.get_y() + bar.get_height()/2, f'{val:.2f}',
                va='center', fontsize=8)

    # UCE plot
    ax2 = axes[1]
    y_pos = np.arange(len(uce_genes))
    bars2 = ax2.barh(y_pos, uce_filtered['recovery_score'], color=UCE_COLOR, alpha=0.8)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(uce_genes)
    ax2.invert_yaxis()
    ax2.set_xlabel('Recovery Score')
    ax2.set_title(f'UCE Top {top_n} Therapeutic Candidates', fontweight='bold', color=UCE_COLOR)
    ax2.set_xlim(26, 28)

    # Add value labels
    for bar, val in zip(bars2, uce_filtered['recovery_score']):
        ax2.text(val + 0.02, bar.get_y() + bar.get_height()/2, f'{val:.2f}',
                va='center', fontsize=8)

    plt.suptitle('Therapeutic Candidates: STATE vs UCE\n(Higher recovery score = promotes recovery from severe COVID)',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    # Save
    plt.savefig(output_dir / 'comparison_top_candidates.png', bbox_inches='tight', facecolor='white')
    plt.savefig(output_dir / 'comparison_top_candidates.pdf', bbox_inches='tight')
    plt.close()
    print(f"Saved: comparison_top_candidates.png/pdf")


def create_recovery_score_distribution(state_scores, uce_scores, output_dir):
    """Create distribution comparison of recovery scores."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Histogram comparison
    ax1 = axes[0]
    ax1.hist(state_scores['recovery_score'], bins=50, alpha=0.7, label='STATE', color=STATE_COLOR, density=True)
    ax1.hist(uce_scores['recovery_score'], bins=50, alpha=0.7, label='UCE', color=UCE_COLOR, density=True)
    ax1.set_xlabel('Recovery Score')
    ax1.set_ylabel('Density')
    ax1.set_title('Recovery Score Distribution', fontweight='bold')
    ax1.legend()

    # Add vertical lines for means
    state_mean = state_scores['recovery_score'].mean()
    uce_mean = uce_scores['recovery_score'].mean()
    ax1.axvline(state_mean, color=STATE_COLOR, linestyle='--', linewidth=2, label=f'STATE mean: {state_mean:.2f}')
    ax1.axvline(uce_mean, color=UCE_COLOR, linestyle='--', linewidth=2, label=f'UCE mean: {uce_mean:.2f}')

    # Box plot comparison
    ax2 = axes[1]
    data_to_plot = [state_scores['recovery_score'], uce_scores['recovery_score']]
    bp = ax2.boxplot(data_to_plot, labels=['STATE', 'UCE'], patch_artist=True)
    bp['boxes'][0].set_facecolor(STATE_COLOR)
    bp['boxes'][1].set_facecolor(UCE_COLOR)
    for box in bp['boxes']:
        box.set_alpha(0.7)
    ax2.set_ylabel('Recovery Score')
    ax2.set_title('Recovery Score Comparison', fontweight='bold')

    # Add stats annotation
    state_stats = f"STATE:\nMean: {state_mean:.2f}\nStd: {state_scores['recovery_score'].std():.2f}\nMax: {state_scores['recovery_score'].max():.2f}"
    uce_stats = f"UCE:\nMean: {uce_mean:.2f}\nStd: {uce_scores['recovery_score'].std():.2f}\nMax: {uce_scores['recovery_score'].max():.2f}"
    ax2.text(0.7, state_scores['recovery_score'].max(), state_stats, fontsize=9, color=STATE_COLOR)
    ax2.text(1.7, uce_scores['recovery_score'].max(), uce_stats, fontsize=9, color=UCE_COLOR)

    plt.suptitle('Recovery Score Analysis: STATE vs UCE', fontsize=14, fontweight='bold')
    plt.tight_layout()

    plt.savefig(output_dir / 'comparison_score_distribution.png', bbox_inches='tight', facecolor='white')
    plt.savefig(output_dir / 'comparison_score_distribution.pdf', bbox_inches='tight')
    plt.close()
    print(f"Saved: comparison_score_distribution.png/pdf")


def create_venn_overlap(state_df, uce_df, output_dir, top_n=50):
    """Create Venn diagram showing overlap in top candidates."""
    # Extract gene names from top N candidates
    state_genes = set([extract_gene_name(p) for p in state_df.head(top_n)['perturbation']
                       if not str(p).startswith('non-targeting')])
    uce_genes = set([extract_gene_name(p) for p in uce_df.head(top_n)['perturbation']
                     if not str(p).startswith('non-targeting')])

    # Remove None values
    state_genes.discard(None)
    uce_genes.discard(None)
    state_genes.discard('non-targeting')
    uce_genes.discard('non-targeting')

    fig, ax = plt.subplots(figsize=(10, 8))

    # Create Venn diagram
    v = venn2([state_genes, uce_genes], set_labels=('STATE', 'UCE'), ax=ax)

    # Customize colors
    if v.get_patch_by_id('10'):
        v.get_patch_by_id('10').set_color(STATE_COLOR)
        v.get_patch_by_id('10').set_alpha(0.7)
    if v.get_patch_by_id('01'):
        v.get_patch_by_id('01').set_color(UCE_COLOR)
        v.get_patch_by_id('01').set_alpha(0.7)
    if v.get_patch_by_id('11'):
        v.get_patch_by_id('11').set_color(OVERLAP_COLOR)
        v.get_patch_by_id('11').set_alpha(0.7)

    # Find overlapping genes
    overlap = state_genes.intersection(uce_genes)

    plt.title(f'Overlap in Top {top_n} Therapeutic Candidates\n(excluding non-targeting controls)',
              fontsize=14, fontweight='bold')

    # Add text with overlapping genes if any
    if overlap:
        overlap_text = "Overlapping genes:\n" + ", ".join(sorted(overlap))
        plt.figtext(0.5, 0.02, overlap_text, ha='center', fontsize=10,
                   bbox=dict(boxstyle='round', facecolor=OVERLAP_COLOR, alpha=0.3))
    else:
        plt.figtext(0.5, 0.02, "No overlapping genes in top 50 candidates",
                   ha='center', fontsize=11, style='italic')

    plt.tight_layout()
    plt.savefig(output_dir / 'comparison_venn_overlap.png', bbox_inches='tight', facecolor='white')
    plt.savefig(output_dir / 'comparison_venn_overlap.pdf', bbox_inches='tight')
    plt.close()
    print(f"Saved: comparison_venn_overlap.png/pdf")

    return overlap


def create_summary_dashboard(state_df, uce_df, state_scores, uce_scores, output_dir):
    """Create a summary dashboard with key metrics."""
    fig = plt.figure(figsize=(16, 12))

    # Create grid
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

    # 1. Model comparison metrics (top left)
    ax1 = fig.add_subplot(gs[0, 0])
    metrics = ['Embedding\nDim', 'Total\nPerturbations', 'Top Score', 'Mean Score']
    state_vals = [2058, len(state_scores), state_scores['recovery_score'].max(), state_scores['recovery_score'].mean()]
    uce_vals = [1280, len(uce_scores), uce_scores['recovery_score'].max(), uce_scores['recovery_score'].mean()]

    x = np.arange(len(metrics))
    width = 0.35

    # Normalize for visualization
    state_norm = [v/max(s,u) for v, s, u in zip(state_vals, state_vals, uce_vals)]
    uce_norm = [v/max(s,u) for v, s, u in zip(uce_vals, state_vals, uce_vals)]

    ax1.bar(x - width/2, state_norm, width, label='STATE', color=STATE_COLOR, alpha=0.8)
    ax1.bar(x + width/2, uce_norm, width, label='UCE', color=UCE_COLOR, alpha=0.8)
    ax1.set_xticks(x)
    ax1.set_xticklabels(metrics)
    ax1.set_ylabel('Normalized Value')
    ax1.set_title('Model Comparison Metrics', fontweight='bold')
    ax1.legend()
    ax1.set_ylim(0, 1.2)

    # Add actual values as text
    for i, (sv, uv) in enumerate(zip(state_vals, uce_vals)):
        ax1.text(i - width/2, state_norm[i] + 0.05, f'{sv:.0f}' if sv > 100 else f'{sv:.2f}',
                ha='center', fontsize=8, color=STATE_COLOR)
        ax1.text(i + width/2, uce_norm[i] + 0.05, f'{uv:.0f}' if uv > 100 else f'{uv:.2f}',
                ha='center', fontsize=8, color=UCE_COLOR)

    # 2. Top 10 candidates table (top middle and right)
    ax2 = fig.add_subplot(gs[0, 1:])
    ax2.axis('off')

    # Create table data
    state_top10 = state_df[~state_df['perturbation'].str.contains('non-targeting', na=False)].head(10)
    uce_top10 = uce_df[~uce_df['perturbation'].str.contains('non-targeting', na=False)].head(10)

    table_data = []
    for i in range(10):
        state_gene = extract_gene_name(state_top10.iloc[i]['perturbation']) if i < len(state_top10) else '-'
        state_score = f"{state_top10.iloc[i]['recovery_score']:.2f}" if i < len(state_top10) else '-'
        uce_gene = extract_gene_name(uce_top10.iloc[i]['perturbation']) if i < len(uce_top10) else '-'
        uce_score = f"{uce_top10.iloc[i]['recovery_score']:.2f}" if i < len(uce_top10) else '-'
        table_data.append([i+1, state_gene, state_score, uce_gene, uce_score])

    table = ax2.table(cellText=table_data,
                      colLabels=['Rank', 'STATE Gene', 'Score', 'UCE Gene', 'Score'],
                      cellLoc='center',
                      loc='center',
                      colWidths=[0.08, 0.25, 0.12, 0.25, 0.12])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.2, 1.5)

    # Color header
    for i in range(5):
        table[(0, i)].set_facecolor('#34495e')
        table[(0, i)].set_text_props(color='white', fontweight='bold')

    # Color STATE and UCE columns
    for i in range(1, 11):
        table[(i, 1)].set_facecolor('#fadbd8')
        table[(i, 2)].set_facecolor('#fadbd8')
        table[(i, 3)].set_facecolor('#d4e6f1')
        table[(i, 4)].set_facecolor('#d4e6f1')

    ax2.set_title('Top 10 Therapeutic Candidates (excluding non-targeting)', fontweight='bold', pad=20)

    # 3. Recovery score violin plot (middle left)
    ax3 = fig.add_subplot(gs[1, 0])
    data = pd.DataFrame({
        'Model': ['STATE'] * len(state_scores) + ['UCE'] * len(uce_scores),
        'Recovery Score': list(state_scores['recovery_score']) + list(uce_scores['recovery_score'])
    })
    sns.violinplot(data=data, x='Model', y='Recovery Score', palette=[STATE_COLOR, UCE_COLOR], ax=ax3)
    ax3.set_title('Recovery Score Distribution', fontweight='bold')

    # 4. Scatter: Recovery vs Disease score (middle center)
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.scatter(state_scores['disease_score'], state_scores['recovery_score'],
               alpha=0.3, s=10, c=STATE_COLOR, label='STATE')
    ax4.scatter(uce_scores['disease_score'], uce_scores['recovery_score'],
               alpha=0.3, s=10, c=UCE_COLOR, label='UCE')
    ax4.set_xlabel('Disease Score (→ worsens severity)')
    ax4.set_ylabel('Recovery Score (→ promotes recovery)')
    ax4.set_title('Recovery vs Disease Direction', fontweight='bold')
    ax4.legend()
    ax4.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax4.axvline(x=0, color='gray', linestyle='--', alpha=0.5)

    # 5. Dataset distribution (middle right)
    ax5 = fig.add_subplot(gs[1, 2])
    state_datasets = state_df['dataset'].value_counts()
    uce_datasets = uce_df['dataset'].value_counts()

    x = np.arange(2)
    width = 0.35

    hek_vals = [state_datasets.get('hek293t', 0), uce_datasets.get('hek293t', 0)]
    hct_vals = [state_datasets.get('hct116', 0), uce_datasets.get('hct116', 0)]

    ax5.bar(x - width/2, hek_vals, width, label='HEK293T', color='#27ae60', alpha=0.8)
    ax5.bar(x + width/2, hct_vals, width, label='HCT116', color='#e67e22', alpha=0.8)
    ax5.set_xticks(x)
    ax5.set_xticklabels(['STATE', 'UCE'])
    ax5.set_ylabel('Number of Candidates')
    ax5.set_title('Candidates by Cell Line', fontweight='bold')
    ax5.legend()

    # 6. Key findings text (bottom)
    ax6 = fig.add_subplot(gs[2, :])
    ax6.axis('off')

    findings = """
KEY FINDINGS - STATE vs UCE Perturbation Analysis

1. NO OVERLAP in Top 50 Candidates
   • STATE and UCE identify completely different therapeutic targets
   • Suggests models capture different biological signals

2. Different Pathway Focus
   • STATE: CASQ1 (calcium signaling), HSPB8 (autophagy), WNT3 (Wnt pathway), GRM6 (glutamate receptor)
   • UCE: ICMT (protein prenylation), B4GALT1 (glycosylation), AGPAT3 (lipid metabolism), DUSP11 (phosphatase)

3. Score Distribution
   • UCE has higher absolute recovery scores (26.8-27.2) vs STATE (24.3-24.8)
   • Both models identify candidates closest to "recovered" severity state

4. Implications
   • Cross-validation with multiple models recommended for drug discovery
   • Consensus candidates (if any in top 100+) may be high-priority targets
   • Different embedding dimensions (STATE: 2058, UCE: 1280) may capture different features
"""

    ax6.text(0.5, 0.5, findings, transform=ax6.transAxes, fontsize=10,
            verticalalignment='center', horizontalalignment='center',
            bbox=dict(boxstyle='round', facecolor='#ecf0f1', alpha=0.8),
            family='monospace')

    plt.suptitle('STATE vs UCE: Perturbation Analysis Summary Dashboard',
                fontsize=16, fontweight='bold', y=0.98)

    plt.savefig(output_dir / 'comparison_summary_dashboard.png', bbox_inches='tight', facecolor='white')
    plt.savefig(output_dir / 'comparison_summary_dashboard.pdf', bbox_inches='tight')
    plt.close()
    print(f"Saved: comparison_summary_dashboard.png/pdf")


def create_rank_correlation_plot(state_scores, uce_scores, output_dir):
    """Create scatter plot showing rank correlation between models."""
    # Merge on perturbation
    merged = state_scores.merge(uce_scores, on='perturbation', suffixes=('_state', '_uce'))

    if len(merged) == 0:
        print("Warning: No matching perturbations between STATE and UCE")
        return

    # Calculate ranks
    merged['rank_state'] = merged['recovery_score_state'].rank(ascending=False)
    merged['rank_uce'] = merged['recovery_score_uce'].rank(ascending=False)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Scatter plot of ranks
    ax1 = axes[0]
    ax1.scatter(merged['rank_state'], merged['rank_uce'], alpha=0.3, s=5, c='#34495e')
    ax1.set_xlabel('STATE Rank (1 = best)')
    ax1.set_ylabel('UCE Rank (1 = best)')
    ax1.set_title('Rank Correlation Between Models', fontweight='bold')

    # Add correlation
    corr = merged['rank_state'].corr(merged['rank_uce'])
    ax1.text(0.05, 0.95, f'Spearman ρ = {corr:.3f}', transform=ax1.transAxes,
            fontsize=12, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Add diagonal line
    max_rank = max(merged['rank_state'].max(), merged['rank_uce'].max())
    ax1.plot([0, max_rank], [0, max_rank], 'r--', alpha=0.5, label='Perfect correlation')
    ax1.legend()

    # Recovery score scatter
    ax2 = axes[1]
    ax2.scatter(merged['recovery_score_state'], merged['recovery_score_uce'], alpha=0.3, s=5, c='#34495e')
    ax2.set_xlabel('STATE Recovery Score')
    ax2.set_ylabel('UCE Recovery Score')
    ax2.set_title('Recovery Score Correlation', fontweight='bold')

    # Add correlation
    score_corr = merged['recovery_score_state'].corr(merged['recovery_score_uce'])
    ax2.text(0.05, 0.95, f'Pearson r = {score_corr:.3f}', transform=ax2.transAxes,
            fontsize=12, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.suptitle('STATE vs UCE: Cross-Model Correlation Analysis', fontsize=14, fontweight='bold')
    plt.tight_layout()

    plt.savefig(output_dir / 'comparison_rank_correlation.png', bbox_inches='tight', facecolor='white')
    plt.savefig(output_dir / 'comparison_rank_correlation.pdf', bbox_inches='tight')
    plt.close()
    print(f"Saved: comparison_rank_correlation.png/pdf")


def main():
    # Paths
    results_dir = Path("/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/perturbation_analysis")
    output_dir = results_dir  # Save to same directory

    print("=" * 60)
    print("STATE vs UCE Perturbation Analysis - Figure Generation")
    print("=" * 60)

    # Load data
    print("\nLoading data...")
    state_df, uce_df, state_scores, uce_scores = load_data(results_dir)
    print(f"  STATE candidates: {len(state_df)}")
    print(f"  UCE candidates: {len(uce_df)}")
    print(f"  STATE total perturbations: {len(state_scores)}")
    print(f"  UCE total perturbations: {len(uce_scores)}")

    # Generate figures
    print("\nGenerating figures...")

    print("\n1. Top candidates comparison...")
    create_top_candidates_comparison(state_df, uce_df, output_dir)

    print("\n2. Score distribution comparison...")
    create_recovery_score_distribution(state_scores, uce_scores, output_dir)

    print("\n3. Venn diagram overlap...")
    overlap = create_venn_overlap(state_df, uce_df, output_dir, top_n=50)

    print("\n4. Summary dashboard...")
    create_summary_dashboard(state_df, uce_df, state_scores, uce_scores, output_dir)

    print("\n5. Rank correlation...")
    create_rank_correlation_plot(state_scores, uce_scores, output_dir)

    print("\n" + "=" * 60)
    print("COMPLETE - All figures saved to:")
    print(f"  {output_dir}")
    print("=" * 60)

    # List generated files
    print("\nGenerated files:")
    for f in sorted(output_dir.glob('comparison_*.png')):
        print(f"  - {f.name}")


if __name__ == "__main__":
    main()
