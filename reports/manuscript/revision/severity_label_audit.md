# Severity/结局标签质检报告（2026-09-03）

范围：322 Sample_IDs（COVID+Sepsis 子集）的 `covid_severity` 与 mortality 标签 y。
比对基准：audit2 从各仓库抠出的权威记录（`literature_audit/T3_patient_sample_mapping/patient_sample_map.csv`）
+ GEO SOFT 一手核查（Qiu）。当前磁盘映射 = v2 草稿（259 patients；Stephenson 已修正、
Combes 池未剔、Wilk/Lee 待并）。

## ✅ 全绿项（可直接信任）

1. **内部一致性**：322/322 样本 severity 与 y 单值；y == (severity=='fatal') 零例外；
   按当前映射分组后**无任何病人跨 draw 的 y 冲突**（v1 时的 3 个"mixed-label"实为
   Stephenson 池分组假象，v2 规则下消失）。
2. **Ren（GSE158055）**：59/59 样本与 GEO 记录匹配；我们的 9 个 fatal =
   仓库记录的 9 个 deceased，零错配。
3. **Stephenson（E-MTAB-10026）**：102/102 匹配；我们的 8 个 fatal patients 与
   仓库 8 例 Death 完全一致（CV0052, CV0180, CV0198, MH8919329, MH9143320,
   MH9143325, MH9179821, newcastle004v2）。
4. **Qiu sepsis（GSE167363）**：2 个 fatal 病人确认（NS ES/NS LS 均为 non-survivor；
   GEO 的 Sex/age 属性错绑不影响结局标签）。
5. **分布与稿件 Table 1 完全一致**（COVID 29 fatal/100 severe/35 moderate/122 mild/
   19 recovered；sepsis 4 fatal/13 severe，Sample_ID 口径）。

## 🟡 良性备注（无需动作）

- Bernardes 3 个病人（004/005/012）跨访视 severity 组合为 moderate/severe/recovered
  混合——draw 日状态，非病人结局；y 全 0 不受影响。印证 severity 在该队列为
  draw-level，论文写作时勿当作病人级结局描述。

## 🔴 无法从仓库验证，需查原文 supplement（建议 Task 4 给 claude science）

| 队列 | 我们的 fatal sid | 问题 |
|---|---|---|
| Combes (GSE163668) | 2（Patient 29、Patient 47） | GEO 只有 covid_status/表型（POSITIVE/SEVERE 等），**无 mortality 字段**；Patient 47 记为 SEVERE——需对 Combes 2021 Nature supplement 核实这两人是否 died in hospital；注意 Patient 29 ≠ Donor 29 |
| Liu (GSE161918) | 4 | 仓库"报了死亡但未对时"；且我们的 HGR 标签来源待确认（sc lanes 为池化）——需对 Liu Cell supplement 核实这 4 个病人的结局与标签出处 |
| Schulte-Schrepping (EGA) | 4 | controlled access，需论文 Table S1 核对 |
| Bernardes | 2（patient 合并后 1） | bulk GSE161777 有 patient number；结局需对 Immunity supplement |

- Wilk / Lee / Arunachalam / Jiang 在我们子集中 **0 fatal**，mortality 标签无需验证
  （**2026-09-07 更新：GSE151263 (Jiang) 因原文无 per-patient 结局、survivor 标签为
  设计推断不可验证，已从 mortality 任务全部剔除（用户决策 B）**——由此所有保留
  标签均可溯源）
  （Lee 的 nCoV9 severe / nCoV10 mild 为同一人采血日 NEWS 评分差异，预期行为）。

## 结论

- **已验证 21/33 fatal sid（64%）**：Ren 9 + Stephenson 8 + Qiu 4，与仓库记录
  100% 一致，零错配。
- **待验证 12/33**：Schulte-Schrepping 4 + Liu 4 + Bernardes 2 + Combes 2，
  需用各原文 supplement 关闭（建议作为 Task 4 发给 claude science）。
- y 的构造逻辑（severity=='fatal'）与跨 draw 一致性：全部通过。

---

## 扩展：非 fatal severity 等级质检（2026-09-03）

### 可逐样本比对的部分

- **Combes**：仓库分组（POSITIVE/SEVERE vs POSITIVE/M/M）与我们的标签比对，
  13 个单病人文库中 **11 一致，2 个存疑**：
  - `GSM4995432`：我们标 severe，仓库 M/M → **疑似错标**（或两口径定义差异），
    已加入 Task 4 待 supplement 裁定。
  - `GSM4995436`（我们标 fatal）：仓库分组 M/M（轻/中症组）——fatal 病人出现在
    M/M 组不寻常，需 supplement 核实（原本就在 fatal 待验证清单）。
  - 池化文库 3 个：仓库分组字段为 "mixed (pooled sample)"，不可比对（已剔除队列）。
- **Lee**：11/11 一致，包括最 tricky 的 nCoV_9（severe）/nCoV_10（mild）
  同一病人对——draw 日 NEWS 评分，仓库与我们的标签完全吻合。
- **Qiu survivors / Jiang**：severe 标注与队列定义自洽（ICU 出院存活 /
  Berlin ARDS + 机械通气）。

### 语义层：各队列 severity 不是同一把尺子（写作与 Limitations 用）

| 队列 | "severe" 的定义 |
|---|---|
| Qiu (sepsis) | qSOFA≥2，结局二分（survivor/non-survivor） |
| Jiang (sepsis-arm) | Berlin ARDS + 全部机械通气 |
| Ren | moderate/severe × progression/convalescence 分期 |
| Stephenson | WHO 式 ordinal（mild= ward 不需吸氧） |
| Liu | NIH 指南（mild/moderate/severe/critical） |
| Schulte-Schrepping | WHO 2–4=mild / 5–7=severe |
| Bernardes | 改良 WHO + CRP/IL-6/ferritin |
| Combes | 是否机械通气 / ICU≤1 天 |
| Wilk | Berlin ARDS，ventilated vs not |
| Lee | NEWS <5 mild / ≥5 severe（采血当日） |
| Arunachalam | 未报（GEO 仅 Moderate/Severe 二值） |

→ 11 个队列至少 6 种不同的 severity 标尺。我们的整合把它们映射到统一的
mild/moderate/severe/recovered/fatal——**这正是主任务只用 mortality 二分类
（severe 归 survived）的根本原因**，Methods 里应保留这句 justification；
细粒度 severity 跨队列不可比，不做、不解释。

### 待 Task 4 补充核实（非 fatal 部分）

- Combes `GSM4995432`（severe vs M/M 冲突）。
- Liu / Schulte-Schrepping / Bernardes / Stephenson：per-patient severity 等级
  对 supplement/SDRF（NIH / WHO / 改良 WHO 各表）。
- Ren：moderate/severe × stage 映射到我们的四级标签的规则核对。

---

## Task 4 / 4b 闭环（2026-09-03，`literature_audit/T4_outcome_severity_verification/`）

**33/33 fatal sid 全部裁定完毕：31 正确，2 错误（已定位更正方案）。**

### 🔴 发现的标签错误（Combes，2 个，必须改）

| 样本 | COMET | 论文记录（Supplementary Table 1, p.3） | 我们的标签 | 更正 |
|---|---|---|---|---|
| GSM4995436 (Patient 29) | 1029 | "POS Mild/Moderate 85 2 2 Viral No 0 Yes"——无 ICU、0 通气天数、**Discharged=Yes（存活）** | fatal（324 细胞） | severity→moderate，**y 1→0** |
| GSM4995442 (Patient 47) | 1047 | **Discharged=Yes（存活）** | fatal（159 细胞） | severity→severe，**y 1→0** |
| GSM4995432 (Patient 14) | 1014 | "POS Mild/Moderate … Yes 0 Yes"——M/M 组 | severe（y=0） | severity→moderate（y 不变） |

→ COVID fatal：27 sid → 25 sid；fatal patients：29 → **27**（COVID 25 + sepsis 2）。
483 个细胞从 fatal 改标 survived。错误来源：2024 整合时 Combes 的 severity 分组
被误映到 fatal。

### ✅ 确认无误的（含出处）

- **Liu** 4/4：HGR0000071 / HGR0000110 / HGR0000141 / HGR0000430 均 outcome=Deceased
  （Table S6B mmc6.xlsx）——同时**关闭了 HGR 标签出处问题**：确为真实 subject 级标签。
- **Schulte-Schrepping** 4 sid / 3 patients：C19-CB-0013（2 draws）、0016、0021
  均 Outcome=deceased（Table S1 mmc1.xlsx）。
- **Bernardes**：Patient 002（002_TA/TA2 均 y=1）确认 died（Recovered=no，
  Critical）。Patient 003 论文亦 died 但**不在我们子集**（无需处理）。
- Ren 9、Stephenson 8、Qiu 4：此前已验。

### 4b 非 fatal 等级：权威表已提取（比对留待 v2 重标后做）

- Liu：33 subjects = 25 Critical / 5 Severe / 3 Moderate（`liu_pbmc_subjects.csv`），
  severity 为 subject 级常量。
- Schulte-Schrepping：WHO 按 visit 记录（多次采样者逐 visit），
  `schulte_cohort1_patients.csv`。
- Bernardes：modified WHO 为 **patient 级**（修正了"visit 级"的假设），
  `bernardes_patients.csv`。
- Stephenson：`Status_on_day_collection`（draw 级）+ `Worst_Clinical_Status`
  （patient 级），不在 SDRF，来自作者 processed h5ad（`stephenson_sample_severity.csv`）。
- Ren：GEO 三字段（severity 两级 × phase × outcome）→ 四级标签的映射规则已记录。
