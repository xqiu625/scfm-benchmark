#!/bin/bash
#SBATCH --job-name=ren_ml_awk
#SBATCH --partition=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/ren_ml_awk_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/logs/ren_ml_awk_%j.err

RAW=/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/monocyte_source/ren
OUT_DIR=/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/bib_revision_tier2/contamination_audit/decontx_multicelltype

awk 'FNR==NR{keep[$1]=1;next} ($2 in keep)' \
    <(tail -n +2 "$OUT_DIR/ren_multilineage_cols.tsv") \
    <(zcat "$RAW/GSE158055_covid19_counts.mtx.gz" | tail -n +4) \
    > "$OUT_DIR/ren_multilineage_triplets.txt"

wc -l "$OUT_DIR/ren_multilineage_triplets.txt"
echo AWK_FILTER_DONE
