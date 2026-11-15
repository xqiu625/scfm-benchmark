# data/

Stable data layers (architecture §2). Large binaries (`*.h5ad`, `*.h5`, `*.pkl`,
`*.npy`, …) are gitignored — the canonical copies live on HPCC under
`configs/hpcc.yaml:data_root`. These folders hold small mappings/manifests and
serve as documented mount points.

| Layer | Purpose | Edit? |
|---|---|---|
| `raw/` | immutable original data (downloads, GEO pulls) | **never edit by hand** |
| `interim/` | temporary transformed data | regenerable |
| `processed/` | reusable cleaned data (e.g. `covid_sepsis_mortality.h5ad`) | regenerable |
| `external/` | references, annotations, ortholog maps | curated |

Each major file should have a contract in `docs/data-contracts.md`.
