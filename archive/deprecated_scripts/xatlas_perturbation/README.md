# X-Atlas / Perturbation Analysis Archive

Archived on 2026-04-16.

**Reason**: X-Atlas/Orion perturbation data uses HCT116 and HEK293T cell lines, which are not biologically comparable to the platelet transcriptomes used in the main analysis. Removed from manuscript scope.

## Contents

- `data_download/` - X-Atlas/Orion dataset download scripts (HCT116 + HEK293T, ~8M cells)
- `embedding_generation/` - STATE/UCE embedding generation for perturbation cells
- `perturbation_analysis/` - Perturbation recovery analysis, severity mapping, comparison figures
- `cell_line_extraction/` - UCE cell-type extraction from HCT116/HEK293T
- `drug_discovery/` - Full drug discovery pipeline (protein structure, docking, repurposing)
- `drug_target_validation/` - ACC 20: Precision@K, therapeutic target validation
- `figures/` - Perturbation landscape, therapeutic ranking, drug target figures

## Note
These analyses may be useful for a future paper focused on perturbation-based drug target discovery using a more appropriate cell type comparison.
