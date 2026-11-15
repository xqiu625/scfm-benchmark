# literature_audit — 外部文献调研结果存档（claude science 执行）

两个 Task 的一手调研结果，全部附原文引用与证据行。调研要求见
`../research_brief.md`（含防幻觉规则：一手来源、逐字引用、查不到写 not reported）。

## T1T2_fm_corpora_and_cohort_windows/（Task 1 + 2，2026-09-02 回）

| 文件 | 内容 |
|---|---|
| `memo.md` | 总 memo：P1-E verdict（"essentially absent" 不可支持，改为 ~0.1% Census 上限 + 注释错位 + STATE 癌细胞系）+ P1-D 采样窗口结论（4 种时间锚，6h–126d）+ gap audit |
| `fm_corpora.csv` | Task 1 主表：5 个 FM 的语料名/规模/类型/platelet 与 megakaryocyte 是否报告 |
| `fm_corpora_evidence.csv` | Task 1 证据行（逐字引用 + 链接） |
| `census_platelet_megakaryocyte.csv` | CELLxGENE Census 各 LTS release 的 platelet/MK-lineage 计数与占比（可复现） |
| `cohort_sampling_windows.csv` | **Task 2 主表（13 列 × 11 队列）**：accession、采样窗口原文、design/timepoints、平台、severity 定义、fatal-case timing、status、source_urls → 新 Table 1 的直接来源 |
| `cohort_sampling_evidence.csv` | Task 2 证据行（93 行） |
| `cohort_sampling_timeline.png` | 11 队列采样窗口时间线图（可作补充图） |
| `fetch_log.csv` | 全部 68 个来源文件的获取清单（provenance） |

## T3_patient_sample_mapping/（Task 3，2026-09-03 回）

| 文件 | 内容 |
|---|---|
| `memo.md` | 总 memo：9 队列 verdict 表；两条错误规则（Stephenson 池/病人反转、Combes 混池）；Liu sc 池化警示；Ren/Bernardes/SS 确认；Wilk/Lee 重复对 |
| `mapping_verdicts.csv` | 每队列一行 × 15 字段（规则、verdict、病人/样本数、draw 含义、结局可用性） |
| `patient_sample_map.csv` | **806 行**：9 个仓库每个 sample 的 patient key、draw、时间、记录结局 —— severity 标签质检与 v2 映射的比对基准 |
| `mapping_evidence.csv` | Task 3 证据行（15/15 逐字核验通过） |
| `stephenson_patient_outcomes.csv` | Stephenson 130 病人的逐人结局（8 例 Death 名单） |
| `stephenson_pool_composition.csv` | 15 个 10x 池 × 成员构成 |
| `stephenson_pool_vs_patient.png` | 池 vs 病人分组示意（mixed-label 由来） |

相关内部文档：`../patient_overlap_finding.md`（泄漏发现与 v2 归并规则）、
`../severity_label_audit.md`（结局标签质检）、`../P1D_cohort_table_skeleton.md`
（cohort 表口径）。

## T7_doublet_prior_evidence/（2026-09-21）

**证据层，非配置层。Gate 2 运行代码禁止读取本目录任何文件。**
完整说明见 `T7_doublet_prior_evidence/README.md`。

| 文件 | 内容 |
|---|---|
| `README.md` | 层级边界、连接键为何不可省、746 capture 的三类覆盖、Ren/Liu 两条易失警示 |
| `memo.md` | 原始调研 survey |
| `doublet_prior_source_audit.csv` | 23 个来源索引（标注 primary/secondary） |
| `doublet_prior_capture_evidence.csv` | 605 行，**混合层级**（603 capture + 2 cohort-level），各自标 `level` |
| `prior_evidence_verdicts.csv` | 8 队列 × 三个正交字段（loading / upstream / rate_model）+ 理由 |
| `capture_evidence_join_key.csv` | **746 行**身份映射：580 mapped、23 cohort_protocol_inheritance（Wilk，`capture_specific=false`）、143 unmapped_to_gem_well（Stephenson）。mapped = **capture-specific evidence linkage**，非直接测量 |
| `doublet-rate-adjudication-v3.md` | stage-alignment 裁定 + 逐队列推理 + **prospective** manuscript 模板（尚未执行，禁止直接使用） |
| `rate-sensitivity-scenarios-v3.csv` | 21 行候选 rate 场景，**21/21 `runtime_eligible=false`** |
| `pooled-same-donor-pairing-scenarios-v3.csv` | 58 行 Σp² 建模配对概率，**58/58 `runtime_eligible=false`** |
| `gate2_source_crosswalk.csv` | 790 源 artifact → 746 capture |
| `geo_capture_cell_counts.csv` | 仓库侧实测 barcode 计数 |
| `ren_gse158055_capture_qc.csv` | Ren Table S1 逐样本 QC 阶梯（证明我们的输入在作者 Scrublet 下游） |
| `liu_demuxlet_per_lane.csv` | Liu **逐 lane** SNG/DBL/AMB —— 是 30,000 白名单上的基因型分类，非真实细胞 doublet 率 |
| `tenx_multiplet_rate_tables.csv` | 10x 官方 multiplet 表（COMBAT 的 50,000 已超出其范围） |

派生表由 `scripts/32_rebuild_upfront_qc/code/build_doublet_prior_evidence_layer.py` 生成。

**本目录任何表都不是 runtime prior configuration。** Gate 2 代码不得直接或经 wrapper 读取；
表中出现数值、列可读、或某行标注 `middle`，均不构成使用授权。将任一 scenario 提升为 runtime
parameter，需要在 `docs/decisions.md` 中作出单独、显式、可审计的决策。

科学结论未变：无队列为 `DIRECT_RESIDUAL_ESTIMATE`，无队列具备可识别的 post-Gate-1 residual
doublet prior；部分队列有 nominal loading 证据，但那与送入 Scrublet 的矩阵中的 residual
prevalence 不是同一阶段。**Gate 2 的 rate 选择仍未解决。**
