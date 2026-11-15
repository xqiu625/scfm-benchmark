"""Materialize the Gate-2 inputs for the two cohorts with cohort-specific rules:
 Liu  -> barcode-rank CALLED cells only (not the full candidate set)
 Wilk -> Seq-Well adaptive QC retained cells (not the rejected 10x floors)
Both keep the FULL lineage population, as Scrublet requires."""
import numpy as np, pandas as pd, os
from scipy.io import mmread, mmwrite
R='/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc'

for cohort, flagfile, keepexpr, src in [
    ('liu',  f'{R}/liu_cellcall/liu_cells_with_callflag.csv', 'called_cell',        f'{R}/gate1_qc/liu'),
    ('wilk', f'{R}/wilk_seqwell_qc/wilk_cells_seqwell_flags.csv', 'not_sw_delete',  f'{R}/gate1_qc/wilk')]:
    out=f'{R}/gate2_input/{cohort}'; os.makedirs(out, exist_ok=True)
    if os.path.exists(f'{out}/DONE'):
        print(f'{cohort}: already materialized'); continue
    flags=pd.read_csv(flagfile, low_memory=False)
    # Gate-1 full flagged frame has ALL cells in the same order as the ORIGINAL
    # cohort matrix, so index alignment is against that, not the retained subset.
    orig_cells=pd.read_csv(f'{src}/qc_flags_per_cell.csv', low_memory=False)
    assert len(flags)==len(orig_cells), f'{cohort}: {len(flags)} vs {len(orig_cells)}'
    if keepexpr=='called_cell':
        keep=flags['called_cell'].to_numpy().astype(bool)
    else:
        keep=(~flags['sw_would_delete'].to_numpy().astype(bool))
    src_mtx=(f'{R}/liu_decontx_input/counts.mtx' if cohort=='liu'
             else f'{R}/wilk_decontx_input/counts.mtx')
    mat=mmread(src_mtx).tocsc()
    assert mat.shape[1]==len(flags), f'{cohort}: matrix {mat.shape[1]} vs flags {len(flags)}'
    mmwrite(f'{out}/counts.mtx', mat[:, keep])
    flags[keep].reset_index(drop=True).to_csv(f'{out}/cells.csv', index=False)
    # prior_counts.csv is NO LONGER WRITTEN.
    #
    # It recorded flags.groupby('channel_id').size() -- the number of cells
    # ENTERING Gate 1 -- under the undefined name `n_pre_qc`, and that value was
    # then fed to a loading-rate formula defined over droplets RECOVERED FROM
    # THE INSTRUMENT. The two are different levels of the count ladder, and the
    # ambiguous column name is what let the mismatch pass unnoticed through six
    # cohorts' worth of invalid scoring.
    #
    # The replacement is docs/data-contracts.md's runtime configuration, which
    # requires an explicit `loading_prior_column` and a decision reference. No
    # authorized configuration exists yet; see docs/execution-plan.md P0.
    #
    # Existing copies under gate2_input/{liu,wilk}/ are left in place rather
    # than deleted -- they are inputs to already-quarantined runs and part of
    # that audit trail -- but nothing may consume them.
    open(f'{out}/DONE','w').write('ok\n')
    print(f'{cohort}: {int(keep.sum())} / {len(keep)} cells kept -> {out}')
