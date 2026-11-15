"""Does the high-mito platelet population retain coherent platelet identity and
usable expression content? Evidence for the Gate-1 mito decision (qc_checklist
0.4 standard: multiple consistent markers, non-platelet signal, complexity,
contamination -- never a single marker, never retention counts alone)."""
import numpy as np, pandas as pd
from scipy.io import mmread

D = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc'
qc = pd.read_csv(f'{D}/qc_flags_trial/cillo/qc_flags_per_cell.csv', low_memory=False)
mat = mmread(f'{D}/cillo_decontx_input/counts.mtx').tocsc()
genes = pd.read_csv(f'{D}/cillo_decontx_input/genes.csv')['gene'].astype(str).tolist()
gidx = {g: i for i, g in enumerate(genes)}

PLT = ['PF4','PPBP','ITGA2B','GP9','GP1BA','TUBB1','NRGN','CAVIN2','SDPR','CLU','TAGLN2','MYL9','SPARC','TMSB4X']
OTHER = {'Monocyte':['LYZ','S100A8','S100A9','CD14','FCN1'],
         'T':['CD3D','CD3E','IL32','TRAC'],
         'B':['MS4A1','CD79A','CD79B'],
         'NK':['NKG7','GNLY','KLRD1'],
         'RBC':['HBB','HBA1','HBA2']}
plt_i = [gidx[g] for g in PLT if g in gidx]
print(f'platelet markers present: {len([g for g in PLT if g in gidx])}/{len(PLT)}')

is_plt = (qc.harmonized_lineage == 'Platelet').to_numpy()
hi = is_plt & (qc.pct_counts_mt > 20).to_numpy()
lo = is_plt & (qc.pct_counts_mt <= 20).to_numpy()
mono = (qc.harmonized_lineage == 'Monocyte_Myeloid').to_numpy()

def cpm(cols):
    sub = mat[:, cols]
    tot = np.asarray(sub.sum(axis=0)).ravel()
    return sub, np.where(tot > 0, tot, 1)

def summarize(name, mask):
    cols = np.where(mask)[0]
    sub, tot = cpm(cols)
    pm = sub[plt_i, :]
    n_marks = np.asarray((pm > 0).sum(axis=0)).ravel()
    plt_frac = 100 * np.asarray(pm.sum(axis=0)).ravel() / tot
    dense = np.asarray(sub.todense()) if len(cols) < 5000 else None
    top1 = (100 * np.asarray(sub.max(axis=0).todense()).ravel() / tot)
    row = {'n_cells': len(cols),
           'median_platelet_markers_detected': float(np.median(n_marks)),
           'pct_cells_with_>=5_markers': round(100*np.mean(n_marks >= 5), 1),
           'median_%counts_in_platelet_markers': round(float(np.median(plt_frac)), 1),
           'median_%counts_in_top1_gene': round(float(np.median(top1)), 1)}
    for lin, gl in OTHER.items():
        idx = [gidx[g] for g in gl if g in gidx]
        f = 100 * np.asarray(mat[idx][:, cols].sum(axis=0)).ravel() / tot
        row[f'median_%counts_{lin}'] = round(float(np.median(f)), 2)
    return pd.Series(row, name=name)

out = pd.DataFrame([summarize('platelet_HIGH_mito(>20)', hi),
                    summarize('platelet_low_mito(<=20)', lo),
                    summarize('monocyte_reference', mono)])
pd.set_option('display.width', 200)
print('\n=== identity & content: high-mito platelets vs low-mito vs monocyte reference ===')
print(out.T.to_string())

# DecontX per-cell contamination (model output -- evidence, not ground truth)
try:
    dec = pd.read_csv(f'{D}/cillo_decontx_input/cillo_decontx_contamination_scores.csv')
    m = qc[['cell_id']].merge(dec[['cell_id','decontX_contamination']], on='cell_id', how='left')
    c = m.decontX_contamination.to_numpy()
    print('\n=== DecontX contamination (model-derived evidence, not ground truth) ===')
    for nm, msk in [('high-mito platelets', hi), ('low-mito platelets', lo), ('monocytes', mono)]:
        v = c[msk]; v = v[np.isfinite(v)]
        if len(v): print(f'  {nm:22s} median={np.median(v):.3f}  p75={np.percentile(v,75):.3f}  n={len(v)}')
except Exception as e:
    print(f'\n(contamination scores unavailable: {e})')
