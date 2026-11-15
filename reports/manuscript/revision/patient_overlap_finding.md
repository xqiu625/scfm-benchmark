# 发现记录：多次采样导致的病人级泄漏（2026-09-02）

> **内部工作文档，不进论文。** 2026-09-02 决策：论文里不写审计/更正故事；
> Methods 用一句话陈述「同一病人多次采样归并为 patient_id 分析单位」，
> 队列直接按疾病组报病人数（COVID-19 240 / sepsis 12 / healthy 70）。
> 本文件仅作为归并规则与影响的技术记录，供重跑与答疑用。
>
> **2026-09-03 Task-3 核查（`literature_audit/T3_patient_sample_mapping/`）已推翻本文档 v1 的两条规则
> （Stephenson 规则反了；Combes 含多病人池）并确认 GSE167363 GEO 标题有误
> （NS ES/NS LS 与论文互换）——归并规则以文末「Task-3 修订（v2）」为准。
> 238 这个数随 v2 重建而更新，当前应视为待修订值。**

## 起因与结论

Adam 问 "When the samples were taken during the disease?"。顺着这个问题核查
`null_model_donor_meta.csv`（Job A 产物，含 Data_NO）发现：**6/11 个队列把同一病人
的多次采样当作不同 `Sample_ID`**，而全部分析（stage-27 CV、cross-disease、
adaptation curve、LOBO、donor-level）都以 `Sample_ID` 为分组/聚合单位。

**322 Sample_IDs → 238 独立病人（COVID 305→226、sepsis 17→12；fatal 33→29，
其中 sepsis fatal 4→2）。**

## 逐队列证据与归并规则（已应用于 `results/donor_meta_with_patient_id.csv`）

| Data_NO | 队列 | 模式 | 规则 | Sample_IDs→patients | fatal sid→pid | 置信 |
|---|---|---|---|---|---|---|
| 302005data03 | Qiu 2021 (sepsis) | `NS1_T0/NS1_T6` | 去 `_T\d+` | 10→5 | 4→2 | **确证**（论文 T0/T+6h within-subject 设计） |
| 302004data16 | Liu 2021 | `HGR0000042_T0/_T1/_T2` | 去 `_T\d+` | 52→31 | 4→4 | 高（time-resolved 设计） |
| 302004data04 | Schulte-Schrepping | `C19-CB-0001_d7/_d16` | 去 `_d\d+` | 23→18 | 4→3 | 高（dN=采样天） |
| 302004data17 | Stephenson | `BGCV01_CV0025/CV0144` | BGCVnn 为病人 | 102→69 | 8→8 | 高（BGCV=病人，CV=访视样本） |
| 302005data01 | Bernardes | `001_TA/TA2/TB/TC` | 去 `_T[A-Z]\d*` | 28→13 | 2→1 | 高（longitudinal；TA2 语义待 GEO 确认） |
| 302004data15 | Ren | `S-M026-1/-2` | 去 `-\d+`（保守归并） | 59→54 | 9→9 | **待验证**（-1/-2 可能是技术重复；无论哪种都同病人） |
| 其余 5 队列 | Wilk/Lee/Arunachalam/Jiang/Combes | 无碰撞 | 不变 | 各不变 | | 见 ⚠ |

⚠ 残留疑点：Combes（302005data02）16 个 ID 全是 GSM 号，**GSM↔病人的映射
需查 GEO metadata**（同一病人的多个 sample 是不同 GSM）——已加入调研清单。

⚠⚠ 2026-09-02 调研（`literature_audit/T1T2_fm_corpora_and_cohort_windows/cohort_sampling_windows.csv`）新增两条线索，
**当前 238 的 patient 数仍可能是上界**：
- Wilk（302004data01，7 IDs）：原文报道 1 个 donor 在 9d 和 11d 各采一次
  → 7 IDs 可能是 6 patients。
- Lee（302004data02，11 IDs `nCoV_1..11`）：原文报道 3 个 donors 各采两次
  → 11 IDs 可能对应 8 patients。
两组的重复采样本规则未命中（ID 形态不带时点后缀），必须查 GEO sample 表
确认 patient 归属后再定。Task 3 其余开放项：Stephenson BGCV 分组 + 病人级
结局（3 个 mixed-label 病人）、Combes GSM 映射、Ren `-1/-2` 语义、
Bernardes `TA2` 语义。

## 影响分析

1. **CV 泄漏**：`StratifiedGroupKFold(groups=Sample_ID)` 可把同一病人的 T0/T6
   分进 train 和 test。波及 84 个多余 Sample_IDs / ~50 个多 draw 病人。
   方向：AUC 偏乐观（draw 间高度相关），幅度待重跑量化。
2. **sepsis 臂**：12 病人 / 2 fatal。reverse transfer（sepsis→COVID）与
   adaptation curve 的 n 进一步减半；adaptation curve 当初 "5/17 donors 留出"
   的划分也可能切开了同一病人的 draw → 必须按 patient 重跑。
3. **donor-level 聚合**（P1-B）：T0/T6 应归并为一个 patient 分数，n=238。
4. **null model**（Job A）：分组同理需按 patient 重跑；预期结论方向不变。
5. LOBO 按队列留出，不受影响；per-cohort AUC 不受分组定义影响。

## 修订中的处理（不改数据，改分组与表述）

- 新增 `patient_id` 列作为一切 CV 分组与 patient 聚合的唯一单位；
  `Sample_ID` 保留为采样单元描述。**论文中的可见改动仅此一句 Methods 陈述。**
- 受影响分析全部重跑（HPCC **Job H**，见 runbook）：stage-27 风格 CV（重出
  OOF 预测）→ cross-disease transfer → adaptation curve → stage-31 patient-level
  与 null model。
- 正文队列口径：直接按疾病组报病人数（COVID-19 240 / sepsis 12 / healthy 70）；
  不出现 "322 donors"，也不写 322→238 的更正叙事。
- sepsis 写 "12 patients (2 fatal)"；reverse-transfer / adaptation curve 标
  descriptive only；不做 within-sepsis CV（2 fatal 数学上不支持）。
- Limitations 平叙两句：多次采样归并规则 + sepsis 小 n。不展开、不辩护。

---

## Task-3 修订（v2，2026-09-03，依据 audit2 + GEO SOFT 一手核实）

### Qiu 2021 (GSE167363)：GEO 元数据有误——两个 non-survivor 的性别/年龄绑定互换

一手核实（GEO family SOFT + BioSample efetch + PMC8629881 全文 + 作者确认 2026-09-03）：
- **论文/作者口径（正确）**：NS LS = P50 = 女 65–70（入组 24h 内死亡，
  APACHE 38 / SOFA 16）；NS ES = P34 = 男 90–95（30 天内死亡，
  APACHE 18 / SOFA 11）。S1 = 女 65–70（survivor，与 NS LS 同性别同年龄段，
  但为不同病人——Table 1 中两列并存）。
- **GEO 当前记录（错误）**：标题 "NS LS T0/T6"（GSM5102902/3，BioSample
  SAMN18039640/39）绑定 Sex=Male, age 90–95；标题 "NS ES_T0/T6"
  （GSM5511351/2，SAMN20695241/…）绑定 Sex=Female, age 65–70。
  GSM characteristics 与 BioSample attributes 两层一致地错 → 提交时两个
  non-survivor 的 demographics 被交叉绑定（2020 批与 2021 批同错，说明
  错误在提交源表）。
- 附带瑕疵：标题分隔符不一（"NS LS T0" 空格 vs "NS ES_T0" 下划线）；
  两批 characteristic 键名不同（"donor:" vs "donor id:"）。
- **对 benchmark 无实质影响**：两名患者均为 fatal（y=1 不变）；T0/T6 配对
  按 donor id 一致，patient 分组不受影响。影响面：cohort 表 fatal-case
  timing 的归属（"24h 内死亡"对应 NS LS = 女 65–70，不可挂错 GSM）；
  以及任何按 GEO demographics 做的 per-patient 描述。
- **更正路径（数据所有者）**：联系 GEO（geo@ncbi.nlm.nih.gov）更正
  GSM5102902/3 与 GSM5511351/2 的 Sex/age 属性（保持标题不动——标题与论文
  一致且已被下游引用）。

### Stephenson 2021：v1 规则反了——BGCVnn 是 10x 测序池，CVxxxx 才是病人

- ArrayExpress sample table：`Characteristics[individual]`=CV0052，
  `Characteristics[sequence pool]`=BGCV02。15 个池 × 每池 4 个不同病人；
  13/15 个池内 Outcome 不一 → v1 的「3 个 mixed-label 病人」实为池内不同
  病人，**按正确规则后 mixed-label 消失**（全部 130 病人 Outcome 唯一）。
- 规则 v2：patient = 尾部 id（CVxxxx / MH* / newcastle* / AP*）；同一 CV id
  出现在不同 BGCV 前缀下 = 同一病人的重采（D0/D28，全库 12 对）。
  302004data17 的 102 IDs 需按此重数（v1 估的 69 作废，预计 ~90）。
- 全库 8 例死亡：CV0052, CV0180, CV0198, MH8919329, MH9143320, MH9143325,
  MH9179821, newcastle004v2（`audit2/stephenson_patient_outcomes.csv`）。

### Combes 2021：3 个「donor」是多病人池，应剔除

- 我们的 16 个 GSM 中：**GSM4995425 = Patients 1+2+3、GSM4995426 = Patients
  5+6、GSM4995430 = Patients 50+51**（GEO `pooled: TRUE`，covid_status/
  phenotype = "mixed (pooled sample)"，需 demultiplex 才能归属）。
  这三个目前在 benchmark 里被当作 3 个 survivor（y=0）、合计 2,597 细胞。
- **处理：剔除这 3 个文库**（302005data02: 16→13 单病人文库），并在
  cohort 表注明 pooled libraries excluded。注意 trap：`Patient 29` ≠ `Donor 29`。

### Liu 2021：HGR 命名只适用于 bulk 记录；sc 文库是多病人池

- GSE161918 的 10x lanes 标题为 "Batch N PBMC 10x Lane N …"，timepoint 字面值
  = `mixed`，无 HGR id。**我们 h5ad 里的 HGRxxx_Tn 标签不可能来自 GEO 标题**，
  必须确认其来源（推测为作者 demultiplex 的 processed object）→ 待核。
  规则 v1（52→31）在标签来源确认前标记为「条件成立」。

### Wilk 2020：7→6（重复对属实，映射未闭合）

- 重复对：C1 在 9d（吸氧）与 11d（插管后）各采一次（`covid_555_1/_2`）。
- GEO 只有 7 个 COVID 文库 = 6 病人（论文称 8 samples / 7 patients，缺 1 个）；
  `covid_5xx` ↔ `C1–C7` 映射两方均未给出 → 我们的 C 标签来源待核。

### Lee 2020：11→8（人口学重复对， severity 为采血日状态）

- 论文：C3/C6/C7 三人各采两次。GEO 无 patient 字段，但 11 个 nCoV 样本中
  恰有三对相邻同年龄同性（nCoV 3&4 67F、7&8 73F、9&10 61M）→ 11→8。
- **nCoV 9 标 severe / nCoV 10 标 mild 几乎确定是同一人**：severity 按采血日
  NEWS 评分，会随病程变化——draw 级 severity ≠ 结局。该队列 y 全 0，不受影响。

### 确认无误的规则

- Ren（59→54）：suffix 是样本索引而非 draw 索引，但 stem 与 `patient`
  characteristic 1:1；权威 draw 键 = patient + sampling day；结局全库唯一。
- Bernardes（28→13）：TA/TA2/TB/TB2/TC/TE = 访视码（day 1/3/8/11/14/20，
  名义值 ±1 天；`TD`/`rec` 未定义）。
- Schulte-Schrepping（23→18）：精确计数论证成立（cohort-1: 27 samples /
  18 subjects）；`_dN` 不在 EGA 记录中，天锚=症状发作系推断。
- Jiang（7）：单次横断面（插管 24h 内），7 IDs 均为独立病人。

### v2 行动项

**最终计数（2026-09-04，agent v2 diff 报告为权威）：COVID-19 240 patients
（25 fatal，Task-4 更正后）/ sepsis 12（2 fatal）/ healthy 70；sepsis 扩展臂
（Chongqing）另计 91 donors / 17 fatal。**

1. **重建 `donor_meta_with_patient_id.csv`（v2）**：Stephenson 用尾部 id 规则；
   Combes 剔除 3 个池化文库（细胞同步剔除，下游一切 n 随之更新）；Wilk 7→6、
   Lee 11→8（在映射核实前先按人口学对归并并标注）；Liu 待标签来源核实。
2. **Job H 更新**：第 1 步 verify 改用 v2 映射；第 5 步 within-COVID 主表的
   n（patients/fatal）以 v2 为准；Job H 第 6 步 healthy 归并同步用 v2 规则。
3. **Liu/Wilk 标签来源**：查 canonical h5ad 构建脚本（`scripts/data/` 整合流程）
   确认 HGR 与 C 标签的出处；若为作者 processed object，引用其 demultiplex 说明。
4. （可选，数据所有者行动）向 GEO 申请更正 GSE167363 两个 non-survivor 的标题。
