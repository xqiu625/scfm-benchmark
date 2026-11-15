#!/usr/bin/env python3
"""
Plot the adaptation curve from scripts/29_adaptation_curve summary.csv.

Two-panel figure:
  Panel A: AUC vs n_sepsis added to COVID training (log+0 x-axis), 6 FMs + PCA-200
           from-scratch baseline. Shaded band = +/-1 SE across 10 seeds.
  Panel B: Three-quantity bar chart per method: Y-intercept (n=0), AUC@n=100,
           asymptote (n=all). PCA-200 reference line on each group.

Reads:
  <DATA_DIR>/summary.csv (or path passed via --data)
Writes:
  <OUT_DIR>/figure_adaptation_curve.pdf
  <OUT_DIR>/figure_adaptation_curve.png
"""
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Okabe-Ito + extras, color-blind safe; reuse across manuscript
FM_COLORS = {
    'STATE':            '#0072B2',
    'UCE':              '#E69F00',
    'TranscriptFormer': '#009E73',
    'scPRINT-2':        '#CC79A7',
    'STACK':            '#56B4E9',
    'Geneformer':       '#D55E00',
    'PCA-200+LR':       '#555555',
}
ORDER = ['scPRINT-2', 'STATE', 'TranscriptFormer', 'Geneformer', 'UCE', 'STACK', 'PCA-200+LR']


def load(data_dir):
    df = pd.read_csv(os.path.join(data_dir, 'summary.csv'))
    df['std']    = df['std'].fillna(0.0)
    df['se']     = df['std'] / np.sqrt(df['count'].clip(lower=1))
    df['is_all'] = (df['n_sepsis_req'] == -1)
    return df


def x_pos(n):
    """Plot n on a log scale, with n=0 shown to the left of n=1."""
    return np.where(n == 0, 0.5, n.astype(float))


def panel_a(ax, df):
    for fm in ORDER:
        d = df[df['method'] == fm].copy()
        if d.empty:
            continue
        d_curve = d[~d['is_all']].sort_values('n_sepsis')
        d_asym  = d[d['is_all']]
        c = FM_COLORS[fm]
        ls = '--' if fm == 'PCA-200+LR' else '-'
        xs = x_pos(d_curve['n_sepsis'].values)
        ys = d_curve['mean'].values
        se = d_curve['se'].values
        ax.plot(xs, ys, color=c, linestyle=ls, linewidth=1.7,
                marker='o', markersize=4, label=fm, zorder=3 if fm == 'PCA-200+LR' else 4)
        ax.fill_between(xs, ys - se, ys + se, color=c, alpha=0.15, linewidth=0, zorder=2)
        if not d_asym.empty:
            ax.axhline(d_asym['mean'].iloc[0], color=c, linestyle=':',
                       linewidth=0.8, alpha=0.5, zorder=1)
    ax.set_xscale('log')
    ax.set_xlim(0.4, 1500)
    ax.set_xticks([0.5, 10, 50, 100, 500, 1000])
    ax.set_xticklabels(['0', '10', '50', '100', '500', '1000'])
    ax.set_xlabel('Sepsis cells added to COVID training pool')
    ax.set_ylabel('AUC on held-out sepsis test set')
    ax.set_ylim(0.60, 1.00)
    ax.set_title('A. Adaptation curve: COVID + n sepsis cells $\\rightarrow$ sepsis test',
                 loc='left', fontsize=11)
    ax.legend(loc='lower right', fontsize=8, frameon=False, ncol=1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def panel_b(ax, df):
    fms = [m for m in ORDER if m != 'PCA-200+LR']
    pca = df[df['method'] == 'PCA-200+LR']
    pca_y0    = pca.loc[pca['n_sepsis_req'] == 0, 'mean'].iloc[0]
    pca_100   = pca.loc[pca['n_sepsis_req'] == 100, 'mean'].iloc[0]
    pca_asymp = pca.loc[pca['is_all'], 'mean'].iloc[0]

    width = 0.26
    x = np.arange(len(fms))
    for i, kind in enumerate(['Y-intercept (n=0)', 'AUC at n=100', 'Asymptote (all)']):
        vals, ses = [], []
        for fm in fms:
            d = df[df['method'] == fm]
            if   i == 0: row = d[d['n_sepsis_req'] == 0]
            elif i == 1: row = d[d['n_sepsis_req'] == 100]
            else:        row = d[d['is_all']]
            vals.append(row['mean'].iloc[0])
            ses.append(row['se'].iloc[0])
        offset = (i - 1) * width
        colors = [FM_COLORS[fm] for fm in fms]
        alpha = [0.45, 0.75, 1.0][i]
        ax.bar(x + offset, vals, width, yerr=ses, color=colors, alpha=alpha,
               edgecolor='black', linewidth=0.5, capsize=2, label=kind)

    # PCA reference lines (one per metric)
    for y, label, ls in [(pca_y0,    f'PCA Y-int = {pca_y0:.3f}',    '--'),
                         (pca_100,   f'PCA @n=100 = {pca_100:.3f}',  '-.'),
                         (pca_asymp, f'PCA asymp = {pca_asymp:.3f}', ':')]:
        ax.axhline(y, color='#555555', linestyle=ls, linewidth=0.9, alpha=0.7)
        ax.text(len(fms) - 0.45, y, label, fontsize=7,
                ha='right', va='bottom', color='#444444')

    ax.set_xticks(x)
    ax.set_xticklabels(fms, rotation=20, ha='right', fontsize=9)
    ax.set_ylim(0.60, 1.00)
    ax.set_ylabel('AUC on held-out sepsis test set')
    ax.set_title('B. Three quantities per FM (PCA-200+LR baseline as reference lines)',
                 loc='left', fontsize=11)
    ax.legend(loc='lower left', fontsize=8, frameon=False, ncol=3, columnspacing=1.0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data',    default='/Users/xinruqiu/Downloads/biohub_position_figs/_data/adaptation_curve')
    ap.add_argument('--out-dir', default='/Users/xinruqiu/Desktop/Files/302006/302006_updategit/scripts/29_adaptation_curve/figures')
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    df = load(args.data)
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 5),
                                   gridspec_kw={'width_ratios': [1.0, 1.1], 'wspace': 0.25})
    panel_a(axA, df)
    panel_b(axB, df)
    fig.suptitle('Adaptation curve: how much sepsis data before each model is useful?',
                 fontsize=12, y=1.00)
    pdf_path = os.path.join(args.out_dir, 'figure_adaptation_curve.pdf')
    png_path = os.path.join(args.out_dir, 'figure_adaptation_curve.png')
    fig.savefig(pdf_path, bbox_inches='tight')
    fig.savefig(png_path, bbox_inches='tight', dpi=200)
    print(f'Wrote {pdf_path}')
    print(f'Wrote {png_path}')


if __name__ == '__main__':
    main()
