import h5py
import numpy as np
from collections import defaultdict

f = h5py.File('/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/COMBAT/COMBAT-CITESeq-DATA.h5ad', 'r')
obs = f['obs']

def read_cat(name):
    grp = obs[name]
    if isinstance(grp, h5py.Group):
        codes = grp['codes'][:]
        cats = grp['categories'][:]
        cats = np.array([c.decode() if isinstance(c, bytes) else c for c in cats])
        return cats[codes]
    else:
        arr = grp[:]
        return np.array([c.decode() if isinstance(c, bytes) else c for c in arr])

combat_id = read_cat('COMBAT_ID')
pool_id = read_cat('Pool_ID')
channel_id = read_cat('Channel_ID')

print('n cells:', len(combat_id))
print('n unique donors:', len(set(combat_id)))
print('n unique pools:', len(set(pool_id)))
print('n unique channels:', len(set(channel_id)))

donor_pools = defaultdict(set)
for d, p in zip(combat_id, pool_id):
    donor_pools[d].add(p)

multi_pool_donors = {d: p for d, p in donor_pools.items() if len(p) > 1}
print('donors spanning >1 pool:', len(multi_pool_donors))
for d, p in list(multi_pool_donors.items())[:20]:
    print(' ', d, sorted(p))

pool_list = sorted(set(pool_id))
parent = {p: p for p in pool_list}
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x
def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb:
        parent[ra] = rb

for d, pools in donor_pools.items():
    pools = list(pools)
    for i in range(1, len(pools)):
        union(pools[0], pools[i])

blocks = defaultdict(set)
for p in pool_list:
    blocks[find(p)].add(p)
print('n connected pool-blocks:', len(blocks))
for root, members in blocks.items():
    print(' block root', root, 'pools:', sorted(members))
