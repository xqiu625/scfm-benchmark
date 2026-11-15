#!/usr/bin/env python3
"""
BIB revision memo item I (2026-09-13 memo, Tier 2/3): consensus mortality-
attribution gene set. Reuses the already-computed cross-model masking-
attribution comparison (results_gene_attribution_gsea/cross_model/
consensus_ranking.csv, from compare_masking_attribution_5models.py -- see
project_302006_gsea_pathway_comparison memory) -- no new attribution compute,
just a k-of-5 consensus extraction on data that already exists.

Writes consensus_gene_set.csv: for each gene ranked by all 5 models, its
5-character sign pattern, how many models agree on a single sign (k), and
(for k>=4, the "consensus core") the mean cross-model rank -- ready to (a)
report directly, and (b) hand to score_consensus_signature.py (Tier 3) as
the gene list for a scoreable signature.
"""
import csv
import os

IN_PATH = ('scripts/31_donor_level_revision/results_gene_attribution_gsea/'
           'cross_model/consensus_ranking.csv')
OUT_DIR = 'scripts/31_donor_level_revision/results_gene_attribution_gsea/cross_model'
OUT_PATH = os.path.join(OUT_DIR, 'consensus_gene_set.csv')

rows = list(csv.DictReader(open(IN_PATH)))
n = len(rows)
print(f'{n} genes ranked by all 5 models', flush=True)

for r in rows:
    signs = r['signs']
    n_pos = signs.count('+')
    n_neg = signs.count('-')
    r['k_agreeing'] = max(n_pos, n_neg)
    r['consensus_sign'] = '+' if n_pos > n_neg else ('-' if n_neg > n_pos else '0')
    r['is_unanimous'] = signs in ('+++++', '-----')

by_k = {}
for r in rows:
    by_k.setdefault(int(r['k_agreeing']), []).append(r)

for k in sorted(by_k, reverse=True):
    print(f'  k={k}: {len(by_k[k])} genes', flush=True)

# "consensus core" = k>=4 (agree in all but at most one model)
core = [r for r in rows if int(r['k_agreeing']) >= 4]
core.sort(key=lambda r: float(r['mean_rank']))
print(f'\nConsensus core (k>=4 of 5 models agree on sign): {len(core)} genes', flush=True)
print(f'  positive: {sum(1 for r in core if r["consensus_sign"] == "+")}', flush=True)
print(f'  negative: {sum(1 for r in core if r["consensus_sign"] == "-")}', flush=True)

fieldnames = list(rows[0].keys())
with open(OUT_PATH, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(sorted(rows, key=lambda r: (-int(r['k_agreeing']), float(r['mean_rank']))))
print(f'\nSAVED {OUT_PATH}', flush=True)

print('\nTop 20 consensus-core genes by mean rank:', flush=True)
for r in core[:20]:
    print(f'  {r["gene"]:10s} k={r["k_agreeing"]} sign={r["consensus_sign"]} '
          f'mean_rank={float(r["mean_rank"]):.1f} signs={r["signs"]}', flush=True)
