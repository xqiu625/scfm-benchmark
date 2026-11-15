# HPCC 代码审查报告（2026-09-04）

范围：`scripts/31_donor_level_revision/code/` 全部在跑/待跑脚本（Job H2/S/F），
对照 CV Protocol v2、P0 estimand、Task-3/4 审计结论。

## ✅ 通过（协议合规，无需动作）

| 文件 | 审查点 |
|---|---|
| `repeated_cv.py` | OOF 逻辑正确（每行每 seed 恰好留出一次，NaN 断言）；per-seed AUC + 同 seed 配对 Δ + t-CI/Wilcoxon 符合协议 |
| `rerun_covid_sepsis_mortality_patient_v2.py/.sh` | SEVERITY_FIX 读时应用（不改共享 h5ad）✓；v2 inner-merge 剔池 ✓；`.sh` 传 `--n-seeds 10` ✓；RF 固定 random_state=42，隔离 fold 方差 ✓ |
| `pca_baseline_patient_oof_v2.py` | **无泄漏**（scaler/PCA 仅 fit train）✓；补了 Disease 过滤（原脚本只按 severity，会拉进 Flu/SLE 等）✓；p_rf=p_lr 重复列已注明 LR-only ✓ |
| `build_patient_map_v2.py` | Combes 剔池 / Wilk C1 A\|B / Lee 三对逻辑正确；sample_id 唯一性成立（merge 不会放大细胞数）；diff 报告齐全。**但见 Fix-1** |
| `integrate_chongqing_sepsis.py` | **关键 catch**：用 `group=='Adult Sepsis'`（93−2=91）而非 `cohort=adult`（108 会混入 8 ICU + 7 健康对照）；QC 阈值全部记录；LIMITATION 诚实声明（原 R 脚本不可恢复，用文献默认参数）；cluster+marker 注释而非裸计数 ✓ |
| `jobS_within_sepsis_cv.py` / `jobS_pca_baseline.py` | groups=Sample_ID（Chongqing 一人一一样本）✓；默认 10 seeds ✓；PCA 同样无泄漏 ✓ |
| `jobS_embed/*.sh` | 各 FM 用各自 env；Geneformer 路径名（`geneformer-12L-30M`）为 legacy 命名，但**与 COVID 侧同一 checkpoint**，embedding 空间一致 ✓ |
| `build_monocyte_canonical.py` | provenance 逐细胞标记（curated vs marker 阈值）✓；Stephenson patient_id 用 h5ad 自身字段（CV 规则正确）✓；**QIU_GSM_TO_SAMPLE 经 GEO 文件名核实正确**（P25=GSM5102902/3 non-survived→NS1，P50=GSM5102904/5 survived→S1，NSES=GSM5511351/2 non-survived→NS2）✓ |
| `jobS_cross_disease.py` | StandardScaler 仅 fit train ✓；test 侧 patient-clustered bootstrap（200 次）✓；COVID 侧 patient_id 零缺失断言 ✓ |

## 🔴 必须修（3 项）

**Fix-1：v2 映射文件的 y 未应用 Task-4 标签更正。**
`build_patient_map_v2.py` 没有 SEVERITY_FIX → `donor_meta_with_patient_id_v2.csv`
里 GSM4995436/442 仍 y=1（应为 0）。影响：计数报告多 2 个 fatal；
`null_model_cellcount.py` 从 map 读 y → null model 也带错。
修法：build 脚本加同一个 fix dict（GSM4995436→moderate、GSM4995442→severe、
GSM4995432→moderate）并重算 y，重跑生成 v2 map。

**Fix-2：transfer 脚本的 AUC 是 cell-weighted，与 P0 estimand 不一致。**
`jobS_cross_disease.py`、`cross_disease_fatal_vs_survived_patient_v2.py` 的
AUC 在细胞级计算（bootstrap 按 patient 聚类，但点估计是 cell-weighted）。
P0 规定一切评分为 patient-level 聚合。修法：transfer 也按 patient 聚合
（每 patient 细胞分均值 → patient-level AUC）为主口径；cell-weighted 保留为
secondary 并标注。

**Fix-3：transfer 的 severe 协议与 P0 决定不符。**
两个 transfer 脚本的 COVID 侧均为 severe-EXCLUDED（旧 transfer 协议）。
P0 已定 **severe=survived 为主协议**。修法：severe-included 作为主口径跑一遍，
severe-excluded 保留为敏感性；两个版本同报、标主次，不写解释。

## 🟡 应该修（2 项）

4. `integrate_chongqing_sepsis.py` 缺 **platelet–monocyte 聚合体比例**报告
   （设计文档明确要求）：在 platelet 簇中统计单核 marker（CD14/LST1/FCN1）
   共表达比例并输出到 per-donor 表。
5. 数字同步：Chongqing 实际为 **91 donors / 17 fatal**（93 个 Adult Sepsis 减
   2 个 mort28 缺失），把 runbook Job S、`sepsis_arm_expansion.md`、一页纸里的
   "108 patients / 20 fatal" 全部改过来。

## 📝 备注（2 项）

6. Geneformer 目录名 legacy，两侧同 checkpoint，勿在论文中写错版本号
   （统一写 Geneformer-V2-316M）。
7. Qiu sepsis 的 NS1/NS2 与论文 ES/LS 顺序相反（NS1=P25=男 90–95=论文 NS ES，
   <30 天；NS2=NSES=女 65–70=论文 NS LS，<24h）——y 不受影响，cohort 表
   fatal timing 归属按此写。
