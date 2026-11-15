from scipy.io import mmread
import numpy as np
mat = mmread("/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/fm_embeddings/subsampled_input/decontx_counts.mtx").tocsc()
sums = np.asarray(mat.sum(axis=0)).ravel()
print("n cells:", mat.shape[1])
print("n cells with sum==0:", (sums==0).sum())
print("n cells with sum<10:", (sums<10).sum())
print("min/max/median sum:", sums.min(), sums.max(), np.median(sums))
print("lowest 20 sums:", np.sort(sums)[:20])
