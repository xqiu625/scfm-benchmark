# P1-D — Cohort characteristics 表骨架（含采样时间）

> **已被调研结果取代（2026-09-02）**：正式表 =
> `revision/literature_audit/T1T2_fm_corpora_and_cohort_windows/cohort_sampling_windows.csv`
> （13 列 × 11 队列，含 sampling window 原文、design/timepoints、platform、
> severity 定义、fatal-case timing、status、source_urls；证据行 93 条在
> `cohort_sampling_evidence.csv`）。时间线图：
> `literature_audit/T1T2_fm_corpora_and_cohort_windows/cohort_sampling_timeline.png`（可作补充图）。
> 本文件只保留计数口径说明和遗留事项。

## 计数口径（2026-09-02 决策）

表的计数列报 **patients（括号附 Sample_IDs）**；正文队列描述按疾病组报病人数
（COVID-19 240 / sepsis 12 / healthy 70）；2024 IJMS 整合文引用保留但**不引
其任何数字**（仅支撑 platelet 重要性 + 数据来源一句）。

- patients/Sample_IDs 归并规则：`patient_overlap_finding.md`（内部文档）。
- ⚠ 238 仍为上界：Wilk（1 donor 采 9d/11d 两次）、Lee（3 donors 各采两次）
  待 GEO 核实后可能再减（见 finding 文档 ⚠⚠ 节）。
- healthy 71 为 donor 口径，patient 归并核查在 HPCC Job H 第 6 步。

## 遗留事项（调研未闭合，表中写 "not reported"）

- Stephenson：无数字采样窗口（仅 "before inclusion to interventional trials"）；
  per-donor days-from-onset 需 ArrayExpress SDRF，可取但尚未取。
- Lee：完全无采样时点（论文与 GEO 均无）；per-donor day-of-illness 或在
  supplement，未核。
- Bernardes：scRNA-seq accession 原文写 "to be confirmed"；processed 在
  FastGenomics，bulk 在 GSE161777。
- Schulte-Schrepping：EGA controlled access，只有论文聚合窗口（cohort1 3–20d /
  cohort2 2–25d；图注另称 4–13d，两处矛盾照实并录）。
- Arunachalam：Methods 在 Science supplement（未取到）；platform chemistry 与
  severity 定义 UNVERIFIED；GEO 有 per-sample days_since_symptom_onset（2–16d，
  median 9）。
-  fatal-case timing：Jiang/Combes/Schulte-Schrepping/Wilk/Lee/Arunachalam 未报；
  Liu 报了死亡但未对时；Ren 的 GEO metadata 给出 16 个 deceased 样本确切采血天
  （day 1–64）——表中可用。

## Limitations 措辞

直接用 memo §2.4 的成句（`literature_audit/T1T2_fm_corpora_and_cohort_windows/memo.md`），
核心是：11 队列 4 种时间锚（sepsis 诊断 / 插管 / 入院 / 症状发作），窗口
6h–126 天，两个队列未报；mortality 任务不控制该异质性，cross-cohort 效应量
按此解读。平叙、不辩护。
