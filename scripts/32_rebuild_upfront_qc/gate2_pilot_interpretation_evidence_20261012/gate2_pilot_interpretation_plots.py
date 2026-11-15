"""READ-ONLY plotting of the Gate 2 GSM5102900_HC1 pilot (job 29882644) per-prior
score distributions and automatic thresholds.

Reads the self-contained, tracked `gate2_pilot_histogram_snapshot.json` that
`gate2_pilot_interpretation_analysis.py` writes (a verbatim, sha256-bound copy of
SCORING_MANIFEST.json's per-prior histogram/threshold fields) -- NOT a live,
uncommitted HPCC manifest path. A prior version of this script hardcoded the live
manifest path, so "reproducible from tracked files" was false; this version takes
that path explicitly via --snapshot and fails loudly if it is missing, rather than
silently depending on an unversioned external location.

Plots ONLY the histogram counts/bin edges and threshold values already recorded in
the snapshot for each of the 8 frozen priors. Does not recompute or refit any
histogram or threshold, does not touch the four original output files, does not
read any outcome field.

IMPORTANT CAVEATS, stated on every plot rather than left implicit:
  - observed (blue) and simulated (orange) use INDEPENDENT log-scaled y-axes.
    Bar heights are not comparable across colors by eye; only each color's own
    trend across the x-axis is meaningful on that axis.
  - each prior panel has ITS OWN x-axis scale. expected_doublet_rate is an input
    to scanpy's per-cell score formula, so score values (and therefore threshold
    values) from different priors are not on one shared numeric scale -- a
    threshold at 0.45 in one panel and 0.54 in another are not two cut points on
    the same axis.

Writes one PNG per prior plus one combined overview PNG to --out-dir (refuses to
overwrite an existing non-empty directory).
"""
import argparse
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--snapshot', required=True,
                     help='path to a gate2_pilot_histogram_snapshot.json written by '
                          'gate2_pilot_interpretation_analysis.py')
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()
    if os.path.exists(args.out_dir) and os.listdir(args.out_dir):
        raise SystemExit(f'REFUSING: {args.out_dir} exists and is non-empty -- '
                          'use a new output directory so no prior result is overwritten')
    os.makedirs(args.out_dir, exist_ok=True)

    snapshot = json.load(open(args.snapshot))
    by_prior = snapshot['by_prior']
    rows = sorted(by_prior.values(), key=lambda r: r['expected_doublet_rate'])

    for row in rows:
        if not row['observed_simulated_bin_edges_identical']:
            raise SystemExit(
                f"prior {row['expected_doublet_rate']}: observed and simulated bin edges "
                "differ -- this script assumes they are shared and verified identical by "
                "the analysis script; refusing to plot them as if sharing one x-axis")

    CAVEAT = ('observed/simulated: independent log y-axes, not height-comparable by eye.\n'
              'x-axis scale is PRIOR-SPECIFIC; not comparable across panels/priors.')

    fig, axes = plt.subplots(4, 2, figsize=(11, 15))
    for ax, row in zip(axes.flat, rows):
        edges = row['bin_edges']
        centers = [(edges[i] + edges[i + 1]) / 2 for i in range(len(edges) - 1)]
        width = (edges[1] - edges[0]) * 0.9
        ax.bar(centers, row['observed_counts'], width=width, alpha=0.6, label='observed',
               color='tab:blue')
        ax2 = ax.twinx()
        ax2.bar(centers, row['simulated_counts'], width=width, alpha=0.4, label='simulated',
                color='tab:orange')
        ax.axvline(row['threshold'], color='red', linestyle='--', linewidth=1,
                   label=f"threshold {row['threshold']:.4f}")
        ax.set_yscale('log')
        ax2.set_yscale('log')
        rate = row['expected_doublet_rate']
        n_called = row['n_called']
        ax.set_title(f'prior={rate}  n_called={n_called}  ({n_called/row["n_cells"]*100:.2f}% of {row["n_cells"]})',
                     fontsize=9)
        ax.set_xlabel('Scrublet doublet score (this prior\'s own scale)', fontsize=7)
        ax.set_ylabel('observed count (log)', color='tab:blue', fontsize=8)
        ax2.set_ylabel('simulated count (log)', color='tab:orange', fontsize=8)
        ax.tick_params(axis='both', labelsize=7)
        ax2.tick_params(axis='both', labelsize=7)
        ax.legend(fontsize=6, loc='upper right')

    fig.suptitle('GSM5102900_HC1 pilot (job 29882644): observed vs simulated score\n'
                  'histograms and automatic thresholds, 8 frozen priors\n'
                  'plotted from a tracked, sha256-bound snapshot; no refit\n' + CAVEAT,
                  fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    overview_path = os.path.join(args.out_dir, 'gate2_pilot_per_prior_histograms_overview.png')
    fig.savefig(overview_path, dpi=150)
    plt.close(fig)

    for row in rows:
        edges = row['bin_edges']
        centers = [(edges[i] + edges[i + 1]) / 2 for i in range(len(edges) - 1)]
        width = (edges[1] - edges[0]) * 0.9
        fig, ax = plt.subplots(figsize=(6, 4.6))
        ax.bar(centers, row['observed_counts'], width=width, alpha=0.6, label='observed',
               color='tab:blue')
        ax2 = ax.twinx()
        ax2.bar(centers, row['simulated_counts'], width=width, alpha=0.4, label='simulated',
                color='tab:orange')
        ax.axvline(row['threshold'], color='red', linestyle='--', linewidth=1.2,
                   label=f"threshold {row['threshold']:.4f}")
        ax.set_yscale('log')
        ax2.set_yscale('log')
        rate = row['expected_doublet_rate']
        n_called = row['n_called']
        ax.set_title(f'GSM5102900_HC1 pilot, prior={rate}\n'
                     f'n_called={n_called} ({n_called/row["n_cells"]*100:.2f}% of {row["n_cells"]})\n'
                     + CAVEAT, fontsize=8)
        ax.set_xlabel("Scrublet doublet score (this prior's own scale)")
        ax.set_ylabel('observed count (log)', color='tab:blue')
        ax2.set_ylabel('simulated count (log)', color='tab:orange')
        ax.legend(fontsize=8, loc='upper right')
        fig.tight_layout()
        fname = f'gate2_pilot_prior_{str(rate).replace(".", "p")}_histogram.png'
        fig.savefig(os.path.join(args.out_dir, fname), dpi=150)
        plt.close(fig)
        print('wrote', fname)
    print('wrote', os.path.basename(overview_path))


if __name__ == '__main__':
    main()
