# Gate 2 / Stage 31 方法学审阅与书面裁决

> 后续状态：同日 [十一参数裁决](GATE2_PARAMETER_RULING.md) 已完成 11/11 方法学批准，并明确新的启动条件。下文“十一参数未裁决”是本报告形成时的历史状态；其他适用边界继续有效。Stage 31 固定维度仍未选择。

Reviewer / role：Codex，受项目所有者委托作本次方法学审阅。

Date：2026-09-24，America/Los_Angeles。

Conflicts / prior involvement：AI 审阅者；此前参与过交付材料核验，已读过项目现有意见与部分历史结果摘要。本次判断基于算法、输入定义和证据，而非拟保留的细胞数；不声称盲法、外部人工审阅或预注册。此前缺失的 `20260923_111811.md` 未被重构；本文是新的、完整可查的裁决。

采用 academic-paper-reviewer 的方法学审阅框架，按 Codex adapter 在当前会话完成；没有虚构多个独立审阅者。先记录方法学发现，再给出综合裁决。用户已明确要求审阅及修改，故同步修订项目的现行设计文档；提交审阅的 PDF、CSV 和旧 Stage 31 七文件保留为证据快照。

## 裁决摘要

**Gate 2：批准下述 scoring/sensitivity-only 设计；具体 capture 是否能评分仍取决于原始计数、分析单位及父细胞状态证据。不批准 prevalence 估计或细胞删除。**

**Stage 31：有条件批准先在完整、子采样前的物理 capture 输入上评分，再映射回分析子集的流程重设。拒绝继续用旧 30–60-cell/Sample_ID 设计生成科学证据；现有证据不足以选择新的固定 PCA 维度。**

本次审阅解决了“由谁给出书面方法学意见”的任务。后续工作应落实本文列明的输入与实现条件；无需再以交付排版或缺失旧 AI 意见为由循环审阅。配置、HPCC pilot、detection、outcome analysis 和删除的执行权限仍按项目原有边界处理；本次没有执行这些操作。

## 已审阅的证据

- Gate 2 12 页 PDF：SHA-256 `1ea5507daf5775fbdf8da1c2482ec0a52e511309b3a663968680ab3602343c5a`，及随附 8 个 CSV。使用其中的证据等级；队列的实测标签来自已提供的审计材料，不表示本次重新测量了原始数据。
- Stage 31 `ADDENDUM_2026-09-24.md`：实际 SHA-256 `e6f859b7e12ddf4531fc882d048ae3edf07263e578f3b14b115215dc9a382b98`。`20260924_195645.md` 所写 `6a1a5f…` 已不是当前副本的校验值。
- Stage 31 原七文件、现行 decisions/current-status/execution-plan，以及评分入口的相关定义。原 301/161 维度分布仅用于判断旧设计的可行性。
- 原始 Scrublet 论文、作者软件说明、Scanpy 1.11.x 官方文档及源码，见文末。查阅 1.11.x 文档不替代执行时对实际 1.11.3 安装的核对。

本次没有使用未审阅的 `capture_provenance_v1_2026-09-23/` 草稿。新的 PDF 为 12 页，页面文字边界检查无越界；现有 SHA 清单全部通过。无需再重复旧覆盖表的计算。

## 方法学发现

1. **MAJOR：残余 prevalence 不可由当前材料识别。** 上游过滤可能同时删除 singlet、same-source doublet 和 cross-source doublet；仅给一个 cross-source removal efficiency 不能描述整个选择过程。
2. **MAJOR：过滤后的模拟分布有额外限制。** 用 post-demux/HTO 剩余细胞任意配对，会模拟出已被上游大量排除的跨来源组合。把 prior 调成 residual 值不能纠正这种类别构成差异，也不能把分数解释为同供体 doublet 概率。
3. **MAJOR：prior 网格的排序稳定性不是独立验证。** 固定邻居图后，改变 prior 对分数施加单调变换；排序和基于排序的 observed-versus-simulated AUC 理应不变。它们主要检查实现一致性。
4. **MAJOR：旧 Stage 31 的 parent-state 支持不足。** 先按患者取最多 60 个细胞，再按 Sample_ID 建图，可能丢失可模拟的父状态与局部结构。降低 PCA 维度只改变计算可行性，不能补回丢失的数据。这里没有宣称文献存在普适的“60 个细胞必定失败”阈值。
5. **MINOR：包中 Q4 的十四倍残余率表述需修正。** `s` 从 `1/14` 到 `1` 的比为 14，但归一化后的 residual rate 比一般不是 14，详见 Q4。

Scrublet 论文要求主要父细胞状态在 singlet 数据中有所代表，并指出少数 multiplet 假设及缺失父状态的限制。上述第 2、4 项是将这些条件应用于本项目的判断；完整 capture 输入仍需验证，不能因细胞数变大就宣称有效。[原始论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC6625319/)

## A. Gate 2：Q1–Q5

### Q1：缺少 loading/recovered 时的敏感性范围

**裁决：采用固定的 direct residual expected-rate stress grid：**

```text
0.005, 0.01, 0.02, 0.04, 0.08, 0.12, 0.16, 0.20
```

这里的数值是传给 `expected_doublet_rate` 的算法情景值，统一标签为 `assumed_algorithmic_stress_setting`。它们不是 cohort prevalence、可信区间、经验下界/上界或 vendor 外推结果。0.5%–20% 提供跨度 40 倍的、事先确定的压力范围；低段近倍增、高段保留若干压力点。采用现有候选网格以限制新增选择自由度；不声称八个点是唯一最优设计或经文献校准。

所有合格 capture 使用同一完整网格，按 cohort/input stage 分层解释。高端 12%–20% 仅表示算法压力，不保证真实高负载下两细胞模拟仍足以代表 higher-order multiplets。真实 residual 也没有被证明落在该区间内。

不根据 call counts、保留率、临床关联或运行失败改变端点、选“最佳率”或缩窄网格。缺少 loading 本身不禁止这种评分设计；缺少合法分析单位或计数输入则禁止运行。

### Q2：vendor recovered-stage rate 到 residual 的转换

**裁决：现有证据不支持用于本项目运行的数值转换。** Qiu 的 nominal loading 可描述协议目标，不能当作逐 capture 实测。其 chemistry 冲突对采用 vendor 表构成障碍；在不调用 vendor 模型的算法压力实验中，可保留冲突标签，不必人为挑一张 chemistry 表。COMBAT 的 50,000/channel 也不能通过截顶或线性外推赋予 residual rate。

审阅者推导：设上游 doublet fraction 为 `r`，两细胞独立配对时 same-source 概率为 `s`；三个类别经全部后续步骤的条件存活概率分别为 `a_S`（singlet）、`a_W`（within-source doublet）、`a_X`（cross-source doublet）。则

```text
D = r [s a_W + (1-s) a_X]
r_residual = D / [(1-r) a_S + D]
```

该式是两细胞类别模型中的质量守恒，并非本项目的估计器。上游率和存活概率未知时，无法从最终总细胞数唯一求解。多阶段处理应在各类别的连续条件存活概率中表示；不能将几个边际 removal fractions 简单相减。完全不约束这些概率时，残余率没有比 `[0,1]` 更有信息的普适界；分母为零时 residual 未定义。

### Q3：removal efficiency `e`

**裁决：`e∈[0,1]` 仅进入解释性模型或明确标注的情景推导，不作为本次 scoring 的额外运行轴。没有依据缩窄其区间。**

包内公式对应 `a_S=a_W=1, a_X=1-e` 的特殊情况。`e=1` 只有在 `r,s` 固定、singlet 和 within-source doublet 均未被额外删除时，才给出该模型内的最小 residual。它不是实际 residual 的无条件下界。

必须标明：配对是否近似独立、source 是 donor 还是 HTO sample、loading composition 是否已知、移除是否依赖 donor pair/cell type/UMI、AMB 与不确定 barcode 的处理、额外 QC 对三类的差异选择、higher-order multiplet 是否可忽略。若需要研究不同机制，应报告 `a_S,a_W,a_X` 及阶段定义；不能把全部过程压缩成一个已测量的 `e`。

### Q4：donor-proportion proxy

**裁决：Liu 的 post-demultiplex SNG composition 可作“存活者组成”的描述性测量，作为 loading composition 时只能标 `proxy`；偏差方向不可识别。**

审阅者反例：只考虑 singlet 选择，loading 为 `(0.5,0.5)`、两供体存活概率为 `(1,0.25)`，则观察份额变成 `(0.8,0.2)`，`Σp²` 从 0.50 上升到 0.68。反过来，loading 为 `(0.8,0.2)`、存活概率为 `(0.25,1)`，观察份额为 `(0.5,0.5)`，`Σp²` 从 0.68 降到 0.50。因此不能统一声称 proxy 保守或偏高。

Liu 每 lane 的 SNG/DBL/AMB 汇总计数能量化分类规模，不能识别每个未分配事件的 donor 构成；DBL barcode 也不等于一个可任意分配的 singlet。没有 donor-pair 信息或有依据的选择模型，不批准数值“校正”。

在已知 `K` 个实际来源且独立两细胞配对的模型内，`1/K≤s≤1`。COMBAT 的 `1/14`、Cillo 两个 HTO source 的 `1/2`、Stephenson 四来源的 `1/4` 是各自假设下的等份下界，不是实测 loading 值；Cillo 的 source 不应自动改写为 donor。

修正包内 Q4：当 `e=1`、`0<r<1` 时，`r_res(s=1)/r_res(s=1/14)=14-13r`，仅在 `r→0` 时趋近 14。该修正不改变任何冻结 CSV；此处作为原 PDF 的方法学勘误。

### Q5：二维机制网格还是 residual grid

**裁决：运行设计采用 Q1 的一个 prior 轴；`r,s,e` 留在解释层。** 多组不可识别机制参数可给出同一 residual，重复运行这些组合不会增加可识别性。

同一个 capture、同一次随机重复内，必须固定输入/细胞顺序、HVG、模拟 parent pairs、PCA 表示、邻居图及所有其他参数。已有 runner draft 复用模拟矩阵的做法方向正确；同一个 seed 本身不等于已经证明邻居图一致，需在实际实现中核对。

设 `q` 是平滑后的 simulated-neighbor fraction，`λ=n_sim/n_obs`，`ρ` 是网格 prior。依据 Scrublet 分数式可推得：

```text
Lρ(q) = (qρ/λ) / [(1-ρ)(1-q) + qρ/λ]
odds(Lρ) = odds(q) × odds(ρ) / λ
```

因此固定 `q,λ` 时，网格上排序应一致。自动 histogram threshold 可因非线性尺度变化而改变；应报告 threshold 成功/失败、观察与模拟分布、call-set Jaccard、状态和覆盖，而非把排序稳定当成生物学稳健性。自动 threshold 随网格变化本身是压力实验结果，不使整个实验失效。[作者实现](https://github.com/swolock/scrublet/blob/master/src/scrublet/scrublet.py)

若显示 simulated-versus-observed AUC，应称为分布区分度；observed 中含有未知 doublet，该 AUC 不是具有真值标签的诊断准确率。`detected_fraction / detectable_fraction` 在上述模型失配下也不能升级成实测 prevalence。稳定性并不消除系统性偏差。

## A. 逐队列返回

表中“conditional”指允许设计，且必须满足共同条件 C1–C5 与本行条件；并非已证明输入合格。G 表示 Q1 的完整八点网格。所有 `ρ` 的证据等级均为 `assumed`。

| Cohort | can_design_scoring_experiment | 参数及用途 | 证据标签 | 专属失效/停止条件 |
|---|---|---|---|---|
| Qiu | conditional | G；完整 capture 的标准评分敏感性 | `published_protocol_loading`、`chemistry_unresolved`、`assumed_algorithmic_stress_setting` | GSM→physical capture 或原始计数来源不能证实时停止该单元；若实现依赖 chemistry/vendor rate，须先解决 chemistry 冲突 |
| Chongqing | conditional | G；标准评分敏感性；loading 缺失无需补造 | `loading_unknown`、`upstream_handling_unknown`、`assumed_algorithmic_stress_setting` | 不得把“none reported”当作确认无上游筛除；capture/输入阶段不明而无法限定解释时停止 |
| Cillo | conditional | G；post-HTO 输入的 `residual_input_diagnostic` | `upstream_filter_published`、`composition_assumed`、`residual_null_unvalidated` | HTO sample 与 donor 混同，或将 unrestricted pairing 解释成已校准的 within-sample residual 时拒绝该解释；若仅剩不足以代表父状态的子集，停止科学检测 |
| Liu | conditional | G；post-demux 输入的 `residual_input_diagnostic` | SNG composition 为 `measured_at_post_demux`；loading input 为 `proxy`；`residual_null_unvalidated` | 55 lanes 的文件键未证实为物理 capture 时停止相应单元；RNA/ADT 未分离、barcode 选择阶段不明则不批准评分；不得称 same-donor rate estimate |
| COMBAT | conditional | G；多阶段筛除后的 residual-input 诊断 | `published_nominal_loading_50000`、`composition_assumed`、`multiple_selection_stages` | vendor 外推/截顶赋 prior、Channel_ID→capture 不明、把 score 当作已验证 residual prevalence 时停止相应实现或解释 |
| Ren | conditional | G；保持 `residual_scoring_only`，独立报告 | `author_scrublet_then_manual_filter`、`counts_level_agreement_only` | 不得用再次 Scrublet 作为独立验证；若不能确认计数与分析单元，只保留历史计数描述；任何 barcode 级转移必须另证 identity |
| Stephenson | no | 当前没有获批参数或运行设计 | `unmapped_to_gem_well`、`not_eligible_without_additional_evidence` | 143 个记录尚不能分解为合法 capture；mapping 与处理站点未解决前停止，不用 grid 掩盖单位缺失 |
| Wilk | no（当前包） | 当前没有获批参数或运行设计 | `Seq_Well`、`deposit_filter_stage_unknown`、`not_eligible_without_additional_evidence` | 沉积数据位于 heuristic filter 哪一侧及 per-array identity 不明时停止；禁止套用 10x rate |

**Wilk 的限制范围：** Seq-Well 本身不是 Scrublet 的原理禁区。若来源、array 单位、RNA counts 和父状态支持补齐，可另行设计明确标为 arbitrary prior stress 的 per-array 诊断；不要求先发明 occupancy model 才能讨论纯算法压力。若要估计真实 prevalence，才需平台专属 loading/capture/selection 模型及校准。本次现有输入证据不足，所以维持不合格；不是以“非 10x”代替科学论证。

**上游过滤三队列的共同边界：** unrestricted pooled simulation 可作为存活输入的表达相似度诊断。若目标改为真实 same-source residual detection，需要与来源/选择机制相符的模拟分布及验证；仅改变 `ρ` 不够。本次不批准把当前诊断输出直接纳入患者间 residual prevalence 比较。

### 共同条件与失效规则

- **C1 输入：** 有来源证据的未归一化 RNA counts；整数值检查只是必要条件，不能证明数据从未归一化/取整。排除 ADT、integrated/scaled、log、DecontX-corrected 数据作为加和模拟输入。
- **C2 单位：** 固定物理 capture/array 与原始 barcode 命名空间。技术测序 lane、GSM、Sample_ID、donor 只有经证据映射后才可用作分析键；不把不同 capture 合并来增加细胞数。
- **C3 支持：** 输入包含该已知处理阶段的完整、基本 cell-called/QC 后群体，且未先按目标 lineage 或分析平衡数抽样；若只能用已筛除的输入，明确标记父状态和残余模拟分布的限制。
- **C4 固定实现：** 参数、软件版本、随机性、阈值规则和有效维度条件先固定。非法输入、维度不可行、threshold 失败或 identity mismatch 均为 unknown/NA；不降低单个样本 PCA 来救运行，不把失败写为零或 singlet。
- **C5 解释：** 按 cohort、capture、input stage 报告覆盖和结果；无 top-r 配额、无 outcome 调参、无选“最合意”网格点、无稳定即正确的推论。多 seed 检查若未做，明确不支持跨 seed 稳健性结论。

需要停止的是违反 C1–C4 的运行单元或违反 C5 的推论。网格间差异大、诊断不稳定是应保留的结果，不是删除不利样本或重选参数的理由。未定义的 reliability cutoff 保持未定义；本次不任意设 AUC、Jaccard 或最低 N 作为科学可信门槛。

**Overall Gate 2 decision：Scoring/sensitivity only; no deletion authorization。** 合格单元的设计获上述条件性批准；真实 residual prevalence 和删除规则仍证据不足。

## B. Stage 31 workflow redesign

### 六项书面裁决

1. **Analysis unit：有条件批准按实际物理 capture/array 建模。** donor、患者、Sample_ID、测序 lane 都不自动等于 capture。多个技术文件若来自同一物理反应，可在证明 barcode/gene 空间和去重方式后重建同一输入；多个 capture 不按患者混合。一个 capture 内多个 donor 的模拟范围取决于上游处理阶段，遵守 Gate 2 的 residual-null 限制。
2. **Raw-count source：使用完整的、子采样前的同 capture RNA 群体。** “完整”指保留该已知阶段可用细胞和父状态，不是把空液滴都送入 Scrublet。基本 cell calling/QC 应明示；不先只选 platelet、不用 balanced 60-cell 子集或为 DecontX 抽样拼出的背景。原始整数计数的来源需与内容检查共同成立。[Scanpy 输入接口](https://scanpy.readthedocs.io/en/1.11.x/api/generated/scanpy.pp.scrublet.html)
3. **Capture mapping：须有可追溯 crosswalk。** 至少记录 cohort、source artifact、physical capture、原 barcode、canonical cell ID、Sample_ID、input stage；明确 suffix/前缀转换规则和冲突检查。不能用同样的细胞总数或未经解释的字符串交集证明 identity。source inventory 按项目记录的 302 个 pre-subsample Sample_ID、九个 Data_NO 组交代；301 是旧抽样后审计总体，不能替代 302。
4. **Score-to-subset mapping：按身份作确定性连接，禁止重新评分。** 每个分析子集细胞应唯一连接到一个原始 barcode/capture 及其 status；row order、nearest neighbor 或表达相似度不能替代 identity。未匹配、原输入被 QC 排除、未评分和失败分别保留状态与 NA；不填 False/0。
5. **Cross-sample comparability：固定方法有助于复现，但不校准跨 capture 分数。** 输入组成、细胞数、测序深度、上游选择和阈值均会改变结果。按 capture 先报告；患者描述性汇总以已评分 eligible 细胞为分母，并同时显示总数与未评分覆盖。它不能自动成为跨患者真实 doublet-rate biomarker。作者建议逐样本检测以保留原实验组成，这不授权把任意 Sample_ID 当作实验单位。[作者 best practices](https://github.com/swolock/scrublet)
6. **Failure conditions：** 原始来源无法证明、capture 键不唯一、不同 capture 混合、counts 已变换、缺失关键父状态、单元维度不足、子集身份连接有歧义或覆盖被隐性丢弃，都阻止相应科学检测/转移。无足够可靠的完整输入与映射时，将 Stage 31 doublet-rate outcome 分支排除出证据链；这不要求删除这些患者或细胞。

### PCA 五问的当前回答

| 问题 | 裁决 |
|---|---|
| Fixed `n_prin_comps` | **不选择。现有证据不足以为新流程选择一个固定整数。** |
| Decision population | 新流程通过来源与 capture 验证的完整输入总体；不是旧 all-301 或 eligible-161 |
| Coverage requirement | 新总体须完整记录每个预定单元的 measurement/status；旧 `<=29` 不能替代新维度审计 |
| Maximum acceptable exclusion | **0% 未记录或为清除错误而追认的静默排除。** 允许无法评分者保持 unknown 并披露；未批准任何固定比例的科学分析排除 |
| Invalidation | 总体/输入阶段/预处理改变；维度等于或超过有效最小维度；按样本动态降维；根据临床结果或运行后 coverage 挑值 |

这里批准的是流程方向，尚不是具体 raw-count 来源和 capture crosswalk 已验收。新 dimension audit 应在完整输入上、在 detection/outcome 之前完成。数值选择须结合实际有效维度和结构证据；默认 30 或最大可行 29 都不能直接作为校准依据。固定维度也不会单独赋予跨样本可比性。

**Overall Stage 31 decision：Conditionally approved（workflow redesign）；Insufficient evidence to select a value（fixed PCA）。**

## 综合修订要求与执行边界

编辑性判断：方案可推进为受条件约束的诊断研究。尚不能推进的结论是 cohort residual prevalence、经验证的患者间 doublet-rate 比较，以及依据这些输出删除细胞。本文的限制针对这些具体推论，不否定已正确完成的技术审计。

本次已要求同步修正现行文档中的三处实质问题：把供体 residual 目标写得过强；把固定 prior 下的排序一致性当作验证；runbook 仍教人从 recovered counts/观察 donor shares 算 runtime prior 并给 Wilk 借用 10x 模型。历史证据快照与 append-only 决策保留，新增本裁决说明优先关系。

下一份实现差异应落实本文的分析单位、输入、模拟与映射条件，显式列出所用参数及其适用性，接受一次有范围的静态核对。现有 stress runner 的 `PRODUCTION_PARAMETER_RULING` 仍为空；本文的八点 prior 设计不等于另外 11 个运行参数自动获批。Stage 31 当前完整输入证据尚未验收，所以本次没有凭空填生产参数或绕过 guard。

## 核验与来源

仅以虚构数字核验了质量守恒式的 `e=0/1` 极限、`14-13r` 比值、score odds 恒等式和 prior 网格上的排序性质；未在 cohort 数据上运行检测。修订对象为决策/设计/说明文档，未改算法或现有检测结果。

1. Wolock, Lopez & Klein (2019), *Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data*. DOI: [10.1016/j.cels.2018.11.005](https://doi.org/10.1016/j.cels.2018.11.005)，[全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC6625319/)。核对了模型假设与缺失父状态的限制。
2. [Scrublet 作者仓库](https://github.com/swolock/scrublet)及[分数/阈值实现](https://github.com/swolock/scrublet/blob/master/src/scrublet/scrublet.py)。核对了逐样本建议、Bayesian score 与自动 threshold 机制；以上 odds/敏感性推论由本审阅推导。
3. [Scanpy 1.11.x Scrublet 文档](https://scanpy.readthedocs.io/en/1.11.x/api/generated/scanpy.pp.scrublet.html)及[wrapper 源码](https://github.com/scverse/scanpy/blob/1.11.x/src/scanpy/preprocessing/_scrublet/__init__.py)。核对输入处理、外部 simulated 对象接口及参数含义。

访问日期均为 2026-09-24。部分 cohort 原始论文页面和旧 vendor PDF 在本次网络访问中不可读取，其队列具体事实沿用随包 source audit 的标注，没有将它们伪报为本次逐条外部重新核实。最终裁决按这些条件性输入作出。
