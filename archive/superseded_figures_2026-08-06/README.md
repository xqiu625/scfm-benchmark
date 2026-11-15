# Superseded figure scripts (archived 2026-08-06)

These four files were moved out of `scripts/26_figures/code/` during a
pre-submission hygiene pass. `generate_fig8_lora.py` contains a fabricated
per-fold scatter (invented fold-level numbers rather than data-traced ones)
and must not be used for any manuscript figure. The three shell runners
(`run_all_figures.sh`, `run_generate_main_figures.sh`, `run_fig7_umap_grid.sh`)
are broken: they invoke scripts that do not exist in this clone
(`generate_main_figures.py`, `generate_fig7_umap_grid.py`), so they fail or
silently skip figures. The canonical figure generator is
`scripts/26_figures/code/regenerate_indexed_figures.py`
(`regenerate_main_figures_5fm.py` is retained but deprecated — see the banner
at its top).
