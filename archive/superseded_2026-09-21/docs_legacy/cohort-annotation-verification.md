# 队列标注核查清单（Cohort Annotation Verification Checklist）

> **Archived 2026-09-21.** The verification was completed on 2026-08-25 and
> its durable conclusions were incorporated into `docs/cohort-sources.md`.
> Keep this record as provenance, not as an active checklist.

**Last updated:** 2026-08-25 — **RESOLVED, 情形 A confirmed.** See "核查结果" section
at the end of this file for the full Check 1–3 output and the conclusion.

**背景:** 2026-08-25 的 cohort rename pass 把稿件中的内部队列编号（302004dataXX /
302005dataXX）全部替换为文献队列名，并确认了两件事：

1. 302005data02 = Combes et al.（Nature 2021，COVID-19 研究，PubMed 记录无 sepsis），
   但旧稿件曾称其为 "sepsis cohort"（已按作者确认改标 COVID）。
2. 302004data18 = Jiang et al.（JCI Insight 2020，PMID 32554932），其队列为
   pneumonia + sepsis 患者（早期 ARDS vs 未发展 ARDS 的 sepsis 对照），但其在
   整合数据中的 `Disease` 标注未知。

本清单的目的：在 HPCC 上对 canonical h5ad 做三项只读检查，确定 (a) 现有疾病
标注是否与文献一致（即 cross-disease 等分析是否需要重跑），(b) data18 是否属于
疾病匹配子集。所有检查均为只读，不改动任何数据。

---

## 检查对象

任一 unified-input canonical h5ad 均可（obs 元数据在各文件间一致）：

```
02_EMBEDDINGS/option_c/state_covid_canonical.h5ad
# 或整合后的原始输入对象（47,977 cells × 30,127 genes）
```

---

## Check 1 — Data_NO × Disease 交叉表（核心检查）

```python
import anndata as ad
adata = ad.read_h5ad("02_EMBEDDINGS/option_c/state_covid_canonical.h5ad", backed="r")

# 细胞级
print(adata.obs.groupby("Data_NO")["Disease"].value_counts())

# 供体级（每 donor 一行）
donors = adata.obs.drop_duplicates("Sample_ID")
print(donors.groupby("Data_NO")["Disease"].value_counts())
```

**预期输出（与当前稿件一致）:**

| Data_NO | 文献 | 预期 Disease 标注 |
|---|---|---|
| 302004data04 | Schulte-Schrepping (Cell 2020) | COVID-19 |
| 302004data15 | Ren (Cell 2021) | COVID-19 |
| 302004data16 | Liu (Cell 2021) | COVID-19 |
| 302004data17 | Stephenson (Nat Med 2021) | COVID-19 |
| 302005data01 | Bernardes (Immunity 2020) | **COVID-19** |
| 302005data02 | Combes (Nature 2021) | **COVID-19** |
| 302005data03 | Qiu (J Leukoc Biol 2021) | Sepsis |
| 302004data18 | Jiang (JCI Insight 2020) | **待定 — 见 Check 2** |

**判定:**

- 若 302005data01/02 全部 = COVID-19 → 疾病标注本来就正确，**无需重跑任何分析**
  （2026-08-25 的修改只是修正了正文措辞，数值不受影响）。
- 若 302005data01/02 出现 Sepsis 标注 → cross-disease 双向迁移（tab:crossdisease、
  Figure 3a）、adaptation curve（tab:adaptation_curve、Figure 5）、severe-excluded
  bootstrap 全部被污染，relabel 后需重跑；disease_composition 表的计数也要改。

---

## Check 2 — data18（Jiang）的标注与子集归属

```python
d18 = adata.obs[adata.obs["Data_NO"] == "302004data18"]
print(d18["Disease"].value_counts())
print(d18["covid_severity"].value_counts())
print("donors:", d18["Sample_ID"].nunique(), "cells:", len(d18))
```

**情形 A：data18 标注为 `Sepsis`（在 322-donor 疾病匹配子集内）**

- 17 个 sepsis donors = Qiu（GSE167363 标称 12 samples）+ Jiang（5 donors）= 17，
  数字吻合。分析数值**依然不需要重跑**（sepsis 就是 sepsis，标注是对的）。
- 但需要修叙述（我可以来做）:
  - main.tex per-cohort 段 "the Qiu et al. sepsis cohort" 改为 Qiu + Jiang 两个来源；
  - supplementary 两个 per-cohort 表的 "Sepsis cohort" 分组加 Jiang;
  - 核查 LOBO / per-cohort 表为何只有 7 行（data18 是否因 donor 太少被剔出
    LOBO？如是，需在文中说明）。

**情形 B：data18 标注为 `Non_covid` / `Lung_infect`（在子集之外）**

- 一切维持现状，分析数值不变、文稿不用动。
- 仅修 `docs/cohort-sources.md` 中 Jiang 一行的疾病描述
  （ARDS → pneumonia + sepsis with early ARDS）。
- **如果你决定把 data18 的 5 个供体改标 Sepsis 并纳入主任务**，则需要重跑：
  疾病匹配子集构建（322 → 327 donors）、cross-disease 双向迁移、adaptation
  curve、LOBO 队列集，以及所有 headline 表（tab:performance、tab:delong、
  tab:calibration、tab:batch 等）。这是一个数据决策，请先告知再执行。

---

## Check 3 — sepsis 供体的队列分布（回答 17 donors 的构成）

```python
sepsis = adata.obs[adata.obs["Disease"] == "Sepsis"]
print(sepsis.drop_duplicates("Sample_ID").groupby("Data_NO")["covid_severity"]
      .value_counts())
print("total sepsis donors:", sepsis["Sample_ID"].nunique())
```

**预期（当前稿件口径）:** 合计 17 donors = 4 fatal + 13 severe。
若输出显示 sepsis 供体只来自 302005data03 → 支持情形 B 之外的另一种可能
（即 Qiu 一家即 17 donors，GEO 的 12 samples 只是上传口径差异）；
若同时来自 302005data03 和 302004data18 → 支持情形 A。

---

## Check 4（附带）— `docs/cohort-sources.md` 待修条目

无论上面结果如何，这份文档（被 Data Availability 引用）有既有错误，建议修正：

1. sepsis 队列误记为 "Reyes et al., J Leukoc Biol 2021"。GEO 核实 GSE167363
   的 PubMed 记录为 **Qiu et al. 2021（PMID 34558746）**；Reyes et al. 是
   Nat Med 2020（在稿中仅为背景引用 reyes2020immune）。
2. "Structural caveat" 段称主任务来自 4 个队列（Wilk, Lee, Stephenson,
   Reyes），与现行 7 队列结构（Schulte-Schrepping/Ren/Liu/Stephenson/
   Bernardes/Combes/Qiu）不符，且人名也错。
3. Jiang 一行的疾病栏（ARDS）按 Check 2 的结果更正。
4. 该文档的 Data_NO ↔ 文献映射应与新命名（作者姓）保持一致。

---

## 反馈方式

把 Check 1–3 的完整输出贴回来即可。我将据此：

- 确认"无需重跑"或列出需要重跑的分析清单；
- 视结果修订 main.tex / supplementary.tex 的队列叙述（情形 A）或仅修订
  `docs/cohort-sources.md`（情形 B）。

---

## 核查结果（2026-08-25，跑在 `state_covid_canonical.h5ad`，backed 模式）

**Check 1 — Data_NO × Disease crosstab：** 302005data01（Bernardes）100% `COVID-19`
（6,030 cells / 28 donors + 482 HC）；302005data02（Combes）100% `COVID-19`/`Lung_infect`
（4,582 + 1,749 cells / 16 + 8 donors），**两者均无 `Sepsis` 标注**。→ 疾病标注本来就
正确，author 的 2026-08-25 修改是对的，**headline 数值全部不需要重跑**。

**Check 2 — data18（Jiang）：** 476 cells / 7 donors，`Disease=Sepsis`，
`covid_severity=severe`（全部 7 donor，0 fatal）。data18 确实存在于 canonical h5ad 中
且已标注为 Sepsis。

**Check 3 — sepsis donor 构成：** 302005data03（Qiu）= 10 donors（4 fatal + 6 severe）+
2 HC；302004data18（Jiang）= 7 donors（0 fatal + 7 severe）。合计 **17 = 4 fatal + 13
severe**，与稿件既有的 "17 sepsis donors (4 fatal, 13 survived)" 描述精确吻合。
→ **情形 A 成立**：sepsis 供体确实同时来自 data03 和 data18。

**超出原清单范围的新发现：** 稿件 LOBO/per-cohort 段落的 "seven cohorts" 措辞低估了
实际的 322-donor headline 队列——该队列跨 **11** 个 Data_NO 分组，不是 7 个。除了已列出
的 7 个（Schulte-Schrepping/Ren/Liu/Stephenson/Bernardes/Combes/Qiu，共 290 donors），
还有 4 个队列（Wilk 302004data01、Lee 302004data02、Arunachalam 302004data05、Jiang
302004data18，共 32 donors / 6,090 cells）**每个都是 0 fatal**（全部 severe/moderate，
即全部落在 y=0 存活类），因此 AUC 在这些队列上无法定义，无法单独成为一折 LOBO 或
per-cohort-AUC 行——这正是它们不出现在 Table batch / Figure 4 / Table per_cohort_auc /
FigureS4 中的真正原因（而不是"donor 太少被剔除"）。供体/细胞数精确核对：
290+32=322，36,005+6,090=42,095。**这不是标注错误，是既有正确设计的一个未被文字披露的
事实**——已在 main.tex 的 LOBO 段落 + per-cohort-OOF 段落加入披露句，并新增
`jiang2020ards` 参考文献（Jiang et al. 2020, JCI Insight, PMID 32554932，
doi:10.1172/jci.insight.135678，WebFetch 核实）。**没有任何表格/图的数值改动。**

**Check 4（`docs/cohort-sources.md`）：** 两处既有错误已修正——sepsis 队列
的第一作者由 "Reyes et al." 更正为 "Qiu et al. 2021"（PMID 34558746；Reyes et al. 2020
Nat Med 是稿件中另一处作为背景引用的不同文献，未受影响）；"4-cohort structural caveat"
段落改写为准确的 11-cohort / 7-with-fatal 结构说明。

**结论：本清单已完全解决，无需进一步 HPCC 检查或任何数值重跑。**
