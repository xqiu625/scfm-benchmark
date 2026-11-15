"""Parses the hash-verified signatures source document's own markdown tables into gene
lists, so SCHULTE_SCHREPPING_CANDIDATE_MARKERS_v1.json is reproducible from its source rather
than hand-copied. Refuses to run if the source document's content has changed.
"""
import re
import json
import hashlib

SRC = "codex_output/evidence/mono_20261017_151151_signatures/mono_signatures_source.md"
SRC_SHA256 = "e55ed44f61172ea2378f7a5267c5ed0351cc29db0a630fcfb5ddad129cb86dff"


def load_source(path):
    data = open(path, "rb").read()
    actual = hashlib.sha256(data).hexdigest()
    if actual != SRC_SHA256:
        raise RuntimeError(f"source document changed: expected {SRC_SHA256}, got {actual}")
    return data.decode("utf-8")


def section(text, start_marker, end_marker):
    i = text.index(start_marker)
    j = text.index(end_marker, i + len(start_marker))
    return text[i:j]


def md_table_rows(block):
    rows = []
    for line in block.splitlines():
        line = line.strip()
        if not line.startswith("|"):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if all(set(c) <= {"-", ":"} for c in cells):
            continue
        rows.append(cells)
    return rows


def gene_col(rows, col_name_candidates):
    header = rows[0]
    idx = next((header.index(c) for c in col_name_candidates if c in header), None)
    if idx is None:
        raise ValueError(f"no column named any of {col_name_candidates} in {header}")
    return [r[idx].strip("`").strip() for r in rows[1:] if r[idx].strip("`").strip()]


def extract(text):
    hla_union = re.findall(
        r"`([^`]+)`", section(text, "### 1b. HLA-II 并集", "## 2.")
    )
    nonclassical = sorted(set(gene_col(
        md_table_rows(section(text, "### 4a. Schulte-Schrepping Table S4（non-classical", "### 4b.")),
        ["gene"],
    )))
    s100ahi_c2_c0 = gene_col(
        md_table_rows(section(text, "### 6a. HLA-DRlo S100Ahi", "### 6b.")), ["gene"]
    )
    cd163hi_c2_c2 = gene_col(
        md_table_rows(section(text, "### 6b. HLA-DRlo CD163hi", "### 6c.")), ["gene"]
    )
    s100ahi_c1_c3 = gene_col(
        md_table_rows(section(text, "### 6c. Cohort 1 severe", "## 7.")), ["gene"]
    )
    mild_c2_c1 = gene_col(
        md_table_rows(section(text, "## 5. Mild 相关", "## 6.")), ["gene"]
    )
    isg_c1_c2_single_donor = gene_col(
        md_table_rows(section(text, "#### 2b-i. Cohort 1 monocyte cluster 2", "#### 2b-ii.")),
        ["gene"],
    )
    isg_c2_c2 = sorted(gene_col(
        md_table_rows(section(text, "#### 2b-iii. Cohort 2 monocyte cluster 2", "## 3.")),
        ["gene"],
    ))
    s100_alarmins = sorted(set(gene_col(
        md_table_rows(section(text, "## 3. S100 alarmins", "## 4.")), ["gene"]
    )))
    return {
        "HLA_II_UNION_ORIGINAL_11_RESEARCHER_DERIVED_v1": sorted(hla_union),
        "NONCLASSICAL_MONOCYTE_CD16_C1Q_UNION_PAPER_CLUSTER_MARKER_v1": nonclassical,
        "SEVERE_HLA_DRLO_S100AHI_COHORT2_CLUSTER0_PAPER_CLUSTER_MARKER_v1": s100ahi_c2_c0,
        "SEVERE_HLA_DRLO_CD163HI_COHORT2_CLUSTER2_PAPER_CLUSTER_MARKER_v1": sorted(cd163hi_c2_c2),
        "SEVERE_HLA_DRLO_S100AHI_COHORT1_CLUSTER3_PAPER_CLUSTER_MARKER_v1": s100ahi_c1_c3,
        "MILD_HLA_DRHI_CD83HI_COHORT2_CLUSTER1_PAPER_CLUSTER_MARKER_v1": sorted(mild_c2_c1),
        "ISG_COHORT1_CLUSTER2_SINGLE_DONOR_CAVEAT_PAPER_CLUSTER_MARKER_v1": sorted(isg_c1_c2_single_donor),
        "ISG_COHORT2_CLUSTER2_PAPER_CLUSTER_MARKER_v1": isg_c2_c2,
        "S100_ALARMINS_UNION_PAPER_CLUSTER_MARKER_v1": s100_alarmins,
    }


if __name__ == "__main__":
    extracted = extract(load_source(SRC))
    for name, genes in extracted.items():
        print(f"{name}: n={len(genes)}")
    print(json.dumps(extracted, indent=2))
