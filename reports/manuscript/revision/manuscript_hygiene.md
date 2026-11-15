# 稿件卫生规则与数据漏洞登记（P2 写作期执行）

> 2026-09-07 用户三条规则：①全文禁用"修订史"词汇，发表稿读起来就是第一版；
> ②GSE151263 保留在数据集、仅 mortality 分析剔除；③词汇与数据漏洞全面筛查。

## 一、禁用词汇表（manuscript-facing 文本，含 drafts/*.tex、main.tex、图注、PPT）

**修订史/版本史类（一律禁止）**：legacy、original（历史义）、pre-revision、
earlier/previous version、old(er)、new(er)（版本义）、updated、revised、
corrected/correction/fix(ed)、bug、mislabel、error（指我们的错误）、audit、
v1/v2（分析版本义）、"第几版"、"更正"、"修正"、"审计"、"旧/新（版本义）"。

**允许的中性替代**：
- "legacy sepsis test set" → "the 5-patient GSE167363 sepsis test set"
- "expanded pool" → "the GSE279452 sepsis pool" / "the 91-patient sepsis pool"
- "pre-revision cohort pool" → "the five-disease platelet pool"
- "corrected labels" → 不单独提及；直接用正确数字（Methods 只陈述标签定义）
- "new cohort" → "the GSE279452 sepsis cohort"

**历史义合法例外**：引用文献的年份与版本（"the 2023 release of Census"）、
数据集自身的版本（Genecorpus-30M vs 104M、V1-10M vs V2-316M 模型版本）、
"first" 用于文献新颖性（"to our knowledge, the largest…"）。

## 二、数据漏洞登记（Jiang 类问题的完整清单与状态）

| # | 问题 | 状态 | 论文处理 |
|---|---|---|---|
| 1 | GSE151263 无 per-patient 结局，survivor 标签为设计推断 | **已处理**：保留在数据集；一切 mortality 分析剔除（Job J 重跑 reverse transfer + adaptation curve） | 数据集段一句事实陈述；Table 1 该行标注 |
| 2 | Combes 3 个多病人池化文库（GSM4995425/26/30） | 已剔除（2,597 cells） | Methods 一句（池化文库不进分析） |
| 3 | Combes 3 个错标（2 fatal→survived、1 severe→moderate） | 已更正（读时 SEVERITY_FIX） | 不提及，直接用正确标签 |
| 4 | 6 队列多次采样当独立 donor（Stephenson/Liu/Combes/Bernardes/Ren/Qiu） | 已按 patient_id 归并（322→252） | Methods 一句归并规则 |
| 5 | Stephenson 池/病人规则反转（BGCV=10x 池） | 已按 CVxxxx 归并；3 个"混合标签"实为池内不同病人 | 同 #4 |
| 6 | Liu sc 10x lanes 为多病人池（timepoint=mixed） | 标签来自作者 demultiplex 对象；4 fatal 已对 supplement 确认 Deceased | Table 1 注明"标签来自作者 processed object" |
| 7 | Ren `-1/-2` 后缀（重复 vs 采样） | 保守归并为同病人；结局全库唯一（已验） | 同 #4 |
| 8 | Bernardes 访视级 severity 波动（moderate/severe/recovered 混合） | draw 级状态非结局；y 全 0 不受影响 | 不写；severity 不当病人级结局描述 |
| 9 | Lee severity = 采血日 NEWS（nCoV9/10 同人不同级） | 预期行为；y 全 0 不受影响 | 同 #8 |
| 10 | Wilk 第 8 个样本缺失（论文 8 samples/7 patients，GEO 7 libraries/6 patients） | 与 GEO 一致（7 sid→6 patients） | Table 1 脚注 "GEO holds 7 of 8 reported samples" |
| 11 | GSE167363 GEO 的 Sex/age 与论文错绑（NS ES/LS 交换） | 结局不受影响（两病人均 fatal）；更正邮件已拟 | 不提及（可选：你发 GEO 更正） |
| 12 | GSE279452 mort28 缺失 2 例 / `hosp_days` 语义（draw-day vs 总住院） | 缺失 2 例已从标签集剔除；`hosp_days` 按"入院→采血"理解（字段与采血日期同列） | Table 1 时点列照实写 "median 16 d (range 3–121)" |
| 13 | TF/STACK 的 4-patient 测试集数字（0.897/0.858） | 已标注 test_set 不可比，不进任何正文图 | 不进正文；如提，只作 supplementary 标注 |
| 14 | scPRINT-2 `small-v2.ckpt` 预训练语料无文档 | provenance 无法核实 | Methods 模型表加一句 provenance 注 |
| 15 | Chongqing 28 天死亡 vs in-hospital 口径差 | 已定为 mort28 主口径 | Methods 终点定义一句差异说明 |

## 三、词汇以外的一致性检查（P2 组装前终检）

- 全文搜索禁用词表（第一节）+ "legacy/original/earlier/previous/old/new version/
  corrected/fixed/audit/v1/v2/pre-revision"，命中=0（main.tex 与图注）。
- 每个数字与 `results_patient_v2/` CSV 重算一致（claims-vs-CSV 终检）。
- 术语唯一（Terminology Ledger）：patient-level OOF AUC、repeated-CV、
  severe-as-survived、GSE accession 全文一致；无 "donor"（除文献引用）。
