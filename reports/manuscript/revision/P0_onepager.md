# P0 一页纸 — 发 Adam（终稿，2026-09-08，全-null 上界框架）

> **状态：Jobs A–S 全部完成，全部数字经源 CSV 复核；外部评审 round 3–17 逐条对账、
> 全部关闭（round-18 为 main.tex 一致性 review 触发的单条更正，见 §5）；
> BIB 投稿定位已由用户拍板（§7）；Key Points 板块已正式起草（§7）。**
> 本文档是 P2 全文重写的唯一权威依据——数字与本文档冲突的地方，以本文档为准。
>
> **头条：预声明主对比（最佳 FM vs from-scratch 基线）未通过**——
> permutation P=0.181，recentered bootstrap P=0.073，25/25
> leave-one-fatal-out refits 全 P>0.05；这是唯一的预声明主检验。
>
> **round-12 到 round-15 曾在此打开一条"互补性"支线**（控制 PCA-200+
> 队列后 FM 是否仍额外贡献信息，与"FM 排序是否赢 PCA-200"这一优越性
> 问题分开），历经四轮修正后一度报告为"三个模型稳健显著、一个视比较
> 对象而定、UCE 始终不显著"。**round-16 判定：这个发现主要是一个 harness
> bug，不是真实信号，已用真实计算撤销，不是靠论证——详见 §3a round-16
> 段。** 诊断：round-12–15 用来算这个增量的 OOF 检验框架把一个概率尺度
> （范围窄）的冻结分数喂进未标准化的 L2（C=0.3）逻辑回归重新拟合——
> 单特征重拟合本该几乎无损，实测却丢了 0.07–0.12 AUC；已用门控单元
> 测试核实（`round16_harness_unittest.py`：把 PCA-200 自己的分数套进
> 同一套 harness，精确复现了 round-15 报告的异常损失）并修正（每折内
> 单独 `StandardScaler`，无信息泄漏）。**用修正后的 harness 完整重跑
> 三个比较对象 × 五个 FM，配合全文一贯的 BH 标准：三个比较对象下全部
> 0/5 通过**，唯一残留信号只是 STATE 边缘、未经多重校正的 p=0.025–0.045，
> 达不到本文任何其他地方采用的显著性标准。**"in discrimination"限定词
> 随之收回**——不是因为怀疑检验力不对称这个问题本身，而是因为认真核实
> 后，那个不对称背后并没有一个真实存在的互补信号，只有一个此前未被
> 发现的 harness bug。deviance permutation 检验（round-12 原始的
> winner-corrected P<0.0001）从未经过这个 bug（in-sample 构造，不涉及
> OOF pooling），但已按 round-9 处理已撤回统计量的先例降级为补充交叉
> 检验，不再承担标题结论，因为它现在冗余于一个已经算清楚的、更值得
> 信任的检验族。本文因而是**全-null + 量化上界框架**，与 round-12 之前
> 的原始表述一致——五轮修正探索了一个真问题（不同检验力的统计量给出
> 不同答案），最终结论是探索本身没有找到一个真实的额外信号，而不是
> 找到了一个又被撤回了。
>
> 唯一能分辨的结果，且其内容是关于失败方而非优势方：**双雄（STATE/Geneformer）
> 与 UCE 之间 ~0.12 AUC 的差距**（confirmatory，range-statistic permutation
> P=0.0112，通过 15 对 BH-FDR）——UCE 是全文唯一在任何 null 检验中都未显著
> 超过 count+disease null 的方法（这是"UCE 不适配"的发现，不是"双雄有优势"的
> 发现；语料占比归因假设已被语料测量本身推翻，见 §5）。
>
> **标题已从 "out-of-distribution" 改为 "underrepresented"**（0.12% 的 Census
> 释出版本里是数以万计、有标注的真实存在，不是缺失）。**QC/队列增量价值的
> 预声明主检验是 PCA-200 单项 nested LR（round-10 修正后 P<0.0001）**，其余
> 5 个 FM 各自读数为 exploratory（见 §3a）。**round-9 撤回了一个数学上不可能
> 的"六模型 max-statistic" p 值**（corrected P 小于其自身 family 最小
> uncorrected p，与 round-4 被弃用的 resample-then-max 同一类 bug 的第二次
> 出现）——已整段删除，新增 Methods QC 标准检查防止此类错误再次通过评审，
> 详见 §3a。
>
> **round-10 修正了 nested-LR 检验族的参考分布，"0/6 通过 BH"改为"6/6 通过
> BH"**——round-9 上一轮留作"不影响结论方向"的渐近-vs-经验参考分布差异，
> round-10 指出该差异的方向被想反了（渐近 p 保护的是"显著"结论，不保护
> "不显著"结论，而 nested-LR 族里恰恰有三个"不显著"判定在渐近参考下做出）。
> round-10 猜测的具体机制（"表示分数在置换标签下被重新拟合"）经代码级核实
> 是错的——该分数确认全程冻结（object id + 内容 hash 双重核验）；真正原因
> 是 chi²(1) 渐近参考本身就不适配一个 L2 正则（C=0.3）惩罚过的统计量（未
> 惩罚版本的经验零分布确实紧贴 chi²(1)，惩罚版本不贴）。修正：模型设定不变
> （仍是全文统一使用的 C=0.3 as-fit 设定），只把参考分布从 chi²(1) 换成同一
> 惩罚统计量自身的经验置换零分布（10,000 次，与全文其余两个置换检验同一
> 惯例）——6 个模型中 5 个（STATE/Geneformer/scPRINT-2/PCA-200/scGPT）撞到
> 10,000 次置换的分辨率下限（0 次超出，P<0.0001），**UCE 仍显著但明显最弱**
> （P=0.0022，观测统计量 1.43 vs 其余五者 2.09–5.76）——BH 结论从"0/6 通过"
> 变为"**6/6 通过**"，详见 §3a。"UCE 在任何检验中都未显著超过 null"这一
> 贯穿全文的表述需要收窄为"UCE 在三项检验中两项未过（方法-vs-null 检验
> p=0.157、增量价值 AUC-赛跑口径 p=0.072），在修正后的增量价值 nested-LR
> 口径下以全场最弱的优势（P=0.0022 vs 其余五者 P<0.0001）勉强通过"，不能再
> 无条件地说"任何检验都未通过"。
>
> **round-11：三处文本一致性修复，均已核实、非新计算。** (1) round-10 的
> UCE 措辞收窄只改了 §3a，abstract 与 §5 里各留了一句"never cleared the
> null model in any test"/"在任何检验中都未显著超过...的方法"——已分别
> 改为"the weakest method throughout...the only one that failed to
> clear [the count+disease null model]"（摘要）和明确限定为
> count+disease null 检验（§5）。(2) nested-LR 的 P 值单独摆出会让读者
> 误读为大效应——已在 §3a 建六模型对照表，P 与 AUC 增量成对报告，加一段
> "存在性 vs 幅度"的分工说明，并写明 P<0.0001 是 10,000 次置换 0 次超出
> 的分辨率下限（P<1/10,001），不是算出的具体值。(3) 综述表只有 3 行会让
> 读者先注意到缺席的知名 cell-level scFM benchmark——已加一段点名
> Kedzierska et al. 2025 Genome Biology（WebSearch 核实为真实文献），
> 明确它支撑的是"姊妹领域已确立同型结论"这一论点，不是本表的"病人级
> 无 selection 校正"论点。§7 Key Points 已从散文提升为正式 5 条板块。
>
> **round-12 item 1（真实新计算）+ round-13（同一发现的正确框架 + 缺失的
> 数字）：头条对比换一个问题后翻转为显著，已核实非矛盾、已并列报告，
> 已算出幅度、已改正归因。** round-12 用同一套 winner-corrected 置换机制
> 重跑头条对比（`round12_headline_deviance_test.py`，job 28229781，
> 4×2500=10,000 次患者级标签置换，置换零分布=每次置换下 5 个 FM 里的
> 最大 deviance 统计量）：reduced 模型=QC+队列+PCA-200 冻结分数，
> full 模型=reduced+某 FM 冻结分数——**置换零分布最大值仅 3.819，STATE
> 观测值 4.328 超出，0/10,000 次超出，经验 P<0.0001**。
>
> **round-13 指出原有的"更敏感统计量"框架会招致质疑（如果只是更强版本，
> 为什么论文报告较弱的那个作为 primary？），并给出了正确框架：这是两个
> 不同的问题，不是同一问题的两种检验力。** ΔAUC 问的是**优越性**——
> FM 的分数排序是否比 PCA-200 好，答案否（P=0.181，仍是预声明、部署相关
> 主检验）；deviance 问的是**互补性**——控制 PCA-200+队列后 FM 是否还
> 额外贡献信息，答案是（P<0.0001）。两个高度相关（r=0.71–0.90）的表示
> 完全可以互补而不优越，不需要"谁更接近真相"的论证。**round-13 还指出
> 结果不该只归因于 STATE**：五个 FM 的 deviance 统计量逐一对照零分布
> 百分位——STATE 超过最大值、Geneformer/scPRINT-2/scGPT 超过 99th pct、
> UCE 超过 95th pct——**全部五个表示都携带可检测的互补信息**，这与
> §3a"STATE 与 Geneformer 不可区分（p=0.60）"放在一起才不矛盾。
> **round-13 还要求算出此前只是断言的幅度**：从已拟合的同一批模型直接
> 读出 in-sample nested ΔAUC——STATE +0.0279（与 round-13 自己基于
> deviance 统计量位置模拟出的预期 +0.025 几乎精确吻合）、Geneformer
> +0.0227、scPRINT-2 +0.0214、scGPT +0.0086、UCE +0.0033——约为主对比
> 失败优势（+0.064）的一半，明确低于本文分辨率论证在此样本量下能确立
> 的规模（~0.10 AUC @25 events），两个结果因而是同一发现的两种测法，
> 不是张力。
>
> "no foundation model significantly outperformed the baseline"统一
> 加注"in discrimination"，摘要/§2/§3a/§5/Key Point 1 均已同步；摘要
> 字数从 279（round-12 item 3 省字）回升到 328（round-12 item 1 新增内容
> + 修复一处自相矛盾），是诚实报告一个真实、已量化发现所需的代价，不是
> 疏漏。**round-13 item 4 另要求 Key Point 1 瘦身**——此前塞进了四个
> 统计量、两个从句，~65 词，标题反而被埋在句子中间，而 BIB Key Points
> 是被扫读而非精读的；互补性细节摘要已完整承担，KP1 收窄回单句头条即可。
>
> **round-14：round-13 自己提出的分辨率框架用错了标尺，且原句自相矛盾；
> in-sample 幅度换成 genuine OOF。** round-13 说"+0.0279 明确低于本文
> 分辨率论证能确立的规模（~0.10 AUC@25 events）"——round-14 指出这句
> 自相矛盾（deviance 检验已经以 P<0.0001 探测到它，"探测到的东西低于
> 可分辨下限"不能同时成立）且标尺用错了（~0.10 是两个**独立**模型对比
> 的分辨率，这里的对比是**嵌套**的——reduced/full 共享 13/14 个特征，
> 预测相关性远高于两个独立表示，方差小得多，可分辨规模也小得多）。
> **修正分两步**：(1) 把此前的 in-sample 幅度换成 genuine OOF——
> `round14_nested_auc_oof.py`，把 reduced/full 协变量组合模型重新纳入
> 真正的 10-seed×5-fold OOF 协议（representation 分数仍冻结），STATE
> 从 in-sample +0.0279 变为 **OOF +0.0391**（95% CI [+0.0334,+0.0449]，
> 独立的配对 Wilcoxon p=0.0020，五个 FM 全部 p=0.0020——10-seed 配对
> 检验能给出的最小可能 p 值，意味着全部 10 个 seed 同号，是一致效应
> 不是巧合）；**非预期方向**——OOF 值全部比 in-sample **更大**（不是
> 更小，与 round-14 自己的预期相反），已如实记录、给出一个待验证的
> 可能机制（reduced 模型因 8 个 cohort dummy 在仅 25 事件下的 OOF
> 过拟合更严重），不作为定论。(2) 用这一组 OOF 配对 Wilcoxon 结果本身
> 作为"此嵌套对比在这一样本量下可分辨"的**直接经验证据**，取代此前
> 依赖 round-13 模拟表格的论证——不再需要假设一个相关系数再推算分辨率，
> 数据本身已经给出答案。摘要、§2（新增 "OOF nested ΔAUC" 到禁裸 AUC
> 术语表）、§3a、§5、Key Point 1（新增复合从句）均已同步；摘要 345→364 词。
> **（round-15 更正：这里描述的"配对 Wilcoxon"证据本身被 round-15 判定为
> 无效检验，见下方 round-15 段——本段作为对账记录保留，不删改历史，但
> 引用其结论时请以 round-15 的修正版本为准。）**
>
> **round-15：round-14 用来支撑"OOF 可分辨"的配对 Wilcoxon 本身是错的
> 检验，且更深层的问题是比较对象本身受损——两项都已用真实计算解决，
> 不是靠论证。** (1) 10-seed 配对 Wilcoxon 反映的是同一批 240 人在不同
> CV 划分下的一致性，不是病人抽样层面的可分辨性——与 round-3 已经排除
> 的"seed 间 sd 当误差项"同一类错误，五个模型全部报告同一个"10 个 seed
> 全部同号"的最小可能 p 值（2/2¹⁰），不携带更多信息。**改用 exact
> DeLong**（`round15_delong_and_undamaged.py`——项目里生成
> `incremental_value_test.csv` 的原始 exact DeLong 实现已不在盘上，
> 本脚本复用另一处（`24_option_c_canonical_input`）保留的同一份标准
> 实现，运行前用两项独立核验确认正确：AUC 点估计与 sklearn 精确一致，
> DeLong SE 与患者级 bootstrap SE 高度吻合，ratio=0.997）后：STATE
> p=0.00215、Geneformer p=0.00180、scPRINT-2 p=0.01684 三者稳健显著
> （BH 内部亦通过）；scGPT p=0.06367、UCE p=0.18083 不显著。(2) round-15
> 进一步指出 reduced 比较对象（QC+队列+PCA-200）本身在 OOF 口径下严重
> 受损（AUC=0.635，比 PCA-200 单独低 ~0.09），提出"增量是否只是在修复
> 这一损伤"的判别实验：换用 PCA-200 单独、PCA-200+log1p(n_cells) 两个
> 未受损的比较对象重跑同一构造——**增量不但没有缩小反而变大**（STATE
> +0.037→+0.088/+0.075，Geneformer/scPRINT-2 同向变大），且 scGPT 在
> 两个未受损版本下都转为显著（p=0.014/0.019）——若增量主要是修复过
> 拟合，换成几乎没有过拟合可修复的干净比较对象后增量应收缩趋零，而
> 不是变大；结果因而支持而非推翻互补性解读。**三版本汇总**：STATE/
> Geneformer/scPRINT-2 在全部三个比较对象下都显著；scGPT 在两个未受损
> 版本下显著、退化版本下边缘不显著；UCE 在任何版本下都不显著——一个
> 经受住压力测试的、内部一致的结论。(3) deviance 检验与新的 DeLong
> 检验共享同一个（已知曾受损的）reduced 模型，round-15 指出"方法学
> 独立"这一措辞不准确，已软化为"两者方向一致，且共享比较对象，真正
> 检验结论的是判别实验本身"。(4) **一个 round-15 未要求、但由本次计算
> 顺带暴露的新问题，如实标注、本轮不展开**：round-13 的 deviance-统计量
> 百分位表（"全部五个超过 95th/99th 百分位或最大值"）完全基于 **in-sample**
> 似然比构造，从未像上方 ΔAUC 幅度那样被重新验证为 OOF 版本；它对 UCE/
> scGPT 的"显著"判定与本段 genuine OOF DeLong 的"不显著"判定不一致——
> 大概率是 in-sample 统计量比 OOF 更容易呈现虚高显著性的又一实例，但
> deviance 置换检验本身（含其零分布）从未被重新用 OOF 版本跑过，留作
> 未闭合项，不影响 winner-corrected P<0.0001（"至少一个模型显著"）这一
> 结论，因为 STATE 本身已被上方 DeLong 独立、稳健地确认显著。摘要、§2、
> §3a、§5、Key Point 1 均已同步；摘要 364→361 词。
> **（round-16 更正：round-15 这里"STATE 本身已被上方 DeLong 独立、
> 稳健地确认显著"这句话已不成立——round-16 发现 DeLong 检验用的
> harness 本身有 bug，修正后 STATE 只在三个比较对象下都邊缘、未经多重
> 校正显著，BH 后全部不过。见下方 round-16 段。）**
>
> **round-16：round-15 留下的"未受损比较对象"竟然互相不一致（比它们
> 各自包含的 PCA-200 分数本身丢了 0.07–0.12 AUC，且丢失量不随特征数
> 单调），round-15 的整个"互补性经受住压力测试"结论其实是同一个
> harness bug 的产物，已用真实计算撤销。** round-16 论证一个单特征
> 逻辑回归重拟合在数学上不该丢这么多 AUC（单特征模型是其输入的单调
> 变换），并给出诊断：PCA-200 的冻结分数是概率尺度（窄范围、std 小），
> L2 penalty 在 C=0.3 下作用于**未标准化**的系数——小尺度特征需要大
> 系数才有判别力，大系数恰恰被罚得最狠，斜率被压扁，折间/seed 间截距
> 噪声主导了 pooled OOF 排序。**已用门控单元测试直接核实**
> （`round16_harness_unittest.py`）：把 PCA-200 自己的冻结分数套进
> round-12–15 的同一套 harness 重新拟合，OOF AUC 从原始 0.7423 掉到
> 0.6569——与 round-15 报告的"PCA-200 单独"reduced AUC 精确一致，确认
> 复现的是同一个 bug。**修正**：每折内单独 fit `StandardScaler`（只用
> 训练集，无泄漏），单测里把丢失从 0.0854 降到 0.0182，通过门控阈值。
> **用修正后的 harness 完整重跑三个比较对象 × 五个 FM**
> （`round16_delong_scaled.py`）：三个 reduced OOF AUC 现在都落在
> PCA-200 自身 OOF AUC（0.725）附近的合理范围（0.710–0.724），不再有
> round-15 那种不合理的大幅跌落；**配合全文一贯的 BH 标准，逐一核对每
> 个排名（不是"排名1不过就自动全不过"，那不是 BH 的通用性质），三个
> 比较对象下全部 0/5 通过**——完整 rank-by-rank 表见 §3a。名义
> 显著性也只有 STATE 边缘（p=0.025–0.045，三个比较对象都在但都不过
> 0.05 的多重校正），Geneformer 有一处边缘（p=0.045），其余全部
> p>0.26。**round-12 到 round-15 建立的"全部五个/三个稳健"发现主要
> 是这个 bug 的产物，不是真实的生物学信号**——deviance permutation
> 检验（从未经过这个 bug）已按 round-9 处理已撤回统计量的先例降级为
> 补充交叉检验，不再承担标题结论；"in discrimination"限定词收回，
> 摘要/§2/§3a/§5/Key Point 1 头条恢复为不带限定词的原始表述。新增
> 第四个 Methods QC 标准检查：任何"reduced/full 嵌套对比"harness 投入
> 使用前必须先过一个门控单元测试（把比较对象自身包含的已知分数单独
> 重跑一遍，OOF AUC 必须在原分数 ~0.02–0.04 以内）。摘要 361→328 词
> ——这次的收缩来自一个真实结论的撤回，不是文字精简；详见 §3a
> round-16 段的完整数字、诊断与三个比较对象表格。
>
> **round-17：round-16 撤回的 harness bug 会不会也污染了预声明主检验——
> 已直接验证，答案是不会，机制归因已修正。** round-16 的门控只针对
> OOF pooling harness；预声明增量价值主检验（PCA-200 nested LR）是
> in-sample 单次拟合，没有跨折 pooling，理论上不受同一个机制影响，
> round-17 要求不要只靠论证、要真的跑一遍验证，并顺带核实一个相关问题：
> round-10 曾把这个检验族坍缩的置换零分布笼统归因为"L2 正则"，round-16
> 在另一个 harness 上更精确地定位到"L2 作用于未标准化的窄尺度特征"——
> 这是不是同一个机制？**已用真实计算核实两件事**
> （`round17_pca200_primary_scaling_check.py`）：(1) 标准化后观测 LR
> 统计量从 4.4772 跳到 18.2530，置换零分布中位数从 0.084 移到
> 0.462——几乎精确落在 chi²(1) 的 0.455 上，**确认机制归因应该更新为
> round-16 的精确版本**；(2) 但两个版本的经验置换 P 完全相同（2000 次
> 置换下都是 0.0005）——**因为观测统计量与它自己的置换零分布用同一套
> 代码拟合，扭曲对两边对称地起作用，被参考分布吸收了，不像 round-16
> 的 OOF 情形那样直接污染一个被拿去跟其他未受扭曲数值比较的绝对
> AUC**，这正是 round-17 自己的预判，现已证实。**预声明主检验的
> P<0.0001 不需要重跑替换**，只更新了机制归因文字。同时修正了摘要与
> §5 里违反本文档自身写作规则的过程叙述（"after correcting a scaling
> defect"这类话不该出现在正文/摘要里，只该留在本文档这本账簿里），
> 摘要 328→308 词；修正了 §3a 一处不成立的 BH 一般性推理（"排名1不过
> 就自动全不过"不是 BH 的通用性质，虽然这次三个配置确实是逐一核对后
> 全部不过，已换成完整 rank-by-rank 表）；标注清楚 round-10 家族"6/6
> 通过"与 round-16 家族"0/5 通过"问的是两个不同比较对象，不矛盾，并在
> §5 加入 round-17 提出的一句合并陈述作为全文论点。
>
> **用户澄清（2026-09-09，非评审轮次，用户直接提出）：STATE 只用了 SE
> 组件，未用 ST 组件——已在 §2 加一条评测范围声明。** 这不是"没测全"：
> ST 预测扰动下的状态转变，本任务（观测转录组→死亡结局）没有扰动输入，
> ST 无自然角色；另外四个模型也全部只用 frozen embedding，未动用各自的
> 微调/生成式能力——单独给 STATE 加 ST 会破坏 apples-to-apples 比较。
> 不改变任何已报告数字，是一条需要在 P2 正文 Methods/Discussion 落地的
> scope 声明，防守一个 STATE/Arc Institute 相关审稿人大概率会问的问题。
>
> 写作规则：数字如实摆，正文不叙述更正过程；本文档本身作为对账记录保留
> 必要的溯源标注（哪一轮评审、哪一次重算），但不做辩护性解释。

---

## 1. 队列（论文口径）

- **COVID-19：240 patients（25 fatal）**。
- **sepsis：91 patients（17 fatal，GSE279452，28 天死亡口径）**是唯一进入
  within-sepsis 分析的臂。另有 GSE167363（5 patients，2 fatal）和
  GSE151263/Jiang（7 patients，无结局记录）——三者合计 103 个"sepsis 相关"
  病人，但只有 91+5=96 个有可验证结局，用于逆向 transfer 训练池；GSE151263
  的 7 人不进入任何带标签分析（资格在数据组装期即定，独立于任何结局相关
  结果，不出现在任何结果表/图中）。
- healthy：70 patients（对照，不进 mortality 任务）。
- 多次采样归并为 `patient_id`，Methods 一句话陈述规则。GSE279452 采血中位
  16 天（3–121 天，fatal 更早：median 10 vs 18.5 天）；旧队列时点用 4 种锚
  （6h–126d）。
- 2024 IJMS platelet 整合文仅用于 platelet 生物学重要性论证，不引其任何数字
  （donor 数、细胞数、其分析结果）。

## 2. Estimand（主指标，唯一）

- **评测范围声明（用户 2026-09-09 提出，此前文档从未写明）：全文五个模型
  一律只用 frozen embedding（representation-only），不使用任何模型的
  生成式/微调/任务特化能力。** 对 STATE 这一点需要专门说明，因为它由两个
  组件构成——**SE（State Embedding，产生细胞级 embedding，与
  Geneformer/scGPT/UCE/scPRINT-2 做的是同一类事）和 ST（State
  Transition，预测细胞在某个扰动下会转变到什么状态，用于 perturbation
  response 预测）**——全文只用了 SE 产生的 embedding，未使用 ST。**这不是
  "没测全"，是本任务设计本身用不上 ST**：ST 需要一个扰动输入（"如果给这个
  细胞施加某种药物/基因敲低，它会变成什么状态"），而本文的任务是"给定
  观测到的、未受任何实验扰动的血小板转录组，预测病人是否死亡"——没有扰动
  可言，ST 在这一设计里没有自然可扮演的角色，把它用上需要把任务重新定义
  成完全不同的问题（预测某种干预如何改变病人转录组轨迹），不是本文的
  estimand。**公平性也要求这样做**：另外四个模型在全文里同样只用
  frozen-embedding 用法，没有一个被允许动用其微调/生成式能力（如
  Geneformer 的 in-context 微调、scGPT 的生成式功能）——若单独给 STATE
  加上 ST，是让它比其余四个模型多一件武器，破坏 apples-to-apples 比较，
  不是让评测更完整。**本文因而评的是"预训练 embedding 在这一分类任务上的
  可迁移性"，不是"每个模型的全部能力"**——这本来就是下方 estimand 的既定
  范围，只是此前从未针对 STATE 的双组件结构专门写明；若 STATE/Arc
  Institute 相关审稿人提出此问题，这条声明直接给出答案，属于 P2 正文
  Methods/Discussion 需要落地的一句 scope 声明，不改变任何已报告数字。
- **主指标**：within-COVID patient-level OOF ROC-AUC（每 patient 的细胞预测
  概率取均值聚合到 patient 后计算），**seed-ensemble OOF 口径**（10 个
  repeated-CV seed 的分数按 patient 取均值）。
- **CV 协议**：10 seeds × 5-fold `StratifiedGroupKFold`（patient 分组）；
  报告口径 = repeat-mean±sd。
- **预声明主对比**：最佳 FM vs from-scratch 基线，判别力（ΔAUC）上
  winner-corrected max-statistic **permutation** test（患者级标签置换
  10,000 次；recentered bootstrap 作交叉验证）——**不是**简单的
  resample-then-max（该方法已证实无效，见 §3a）。这是唯一的预声明主检验，
  部署相关（判别力是最终要用的量），不因下条而改变。其余两两对比一律
  exploratory 名义 p。
- **round-12 新增的 exploratory 二级读数（round-13 重新框定为互补性，非
  "更敏感的同一问题"）**：同一 240 人换用 deviance 统计量（C=0.3 惩罚
  nested-LR，winner-corrected 经验置换零分布，与 §3a QC 增量价值检验同一
  机制）——问的是**互补性**（"控制 QC+队列+PCA-200 后 FM 的表示是否还
  额外贡献信息"），与预声明主检验问的**优越性**（"FM 的分数排序是否
  比 PCA-200 更好"）是两个不同的问题，不是同一问题的两种检验力，见 §3a。
- **Pairwise SE**：**exact DeLong**（placement-value 协方差）为准；bootstrap
  仅保留给没有解析解的 max/range 统计量。
- **Cross-disease**：固定 split + 200 次 patient bootstrap。**Adaptation
  curve**：10-seed 训练池子采样。
- **Secondary 术语**（禁裸 "AUC"）：cell-weighted OOF AUC、per-fold mean、
  拼接 pooled LOBO AUC、seed-ensemble OOF AUC、**OOF nested ΔAUC**（reduced-
  vs-full 协变量组合模型在同一 10-seed×5-fold 协议内重新拟合后的配对
  OOF AUC 差，显著性用 exact DeLong 计算——历经 round-14 误用跨-seed
  Wilcoxon、round-15 改为 patient-level DeLong、round-16 发现并修正
  harness 的标准化缺陷三轮修正，**最终结论是配合全文 BH 标准后 0/5 通过，
  不是一个可报告的正向发现**，完整数字与诊断见 §3a round-16 段）。

## 3. 最终数字

### 3a. within-COVID（240 patients / 25 fatal，主结果）

| 模型 | repeat-mean±sd |
|---|---|
| Geneformer | **0.786±0.026** |
| STATE | **0.778±0.033**（与 Geneformer 不可区分，p=0.60） |
| scPRINT-2 | 0.740±0.027 |
| PCA-200+LR | 0.725±0.033 |
| scGPT | 0.720±0.047 |
| UCE | 0.661±0.046 |

**基线三件套**（三者都是独立基线，不合并成一句"非 XX 假象"）：

| 基线 | AUC | 说明 |
|---|---|---|
| null（count + disease） | 0.554 | 仅细胞数 + 疾病标签 |
| cohort-only（诚实 out-of-fold） | 0.486±0.038 | 折内 Laplace-smoothed cohort 死亡率作为唯一特征；范围 0.420–0.557 |
| QC + cell count（as-fit） | 0.618 | QC 3 特征 + log1p(n_cells)，无需任何符号选择；见下方 QC 段脚注：QC 3 特征单独（无 count）as-fit 为 0.376，条件于 depth 后才转为 0.618，是抑制模式（suppression），非异常 |

**预声明主对比：未通过。** 观测最大优势（Geneformer − PCA-200 = +0.064）不
超过挑选程序在零假设下自然产生的量级——**permutation P=0.181**（10,000 次
患者级标签置换，frozen seed-ensemble 分数上的 scoring-step null）；
**recentered bootstrap P=0.073**，方向一致；**25/25 leave-one-fatal-out
refits 全部 P>0.05**（范围 0.141–0.288）——null 结果稳定，非单一死亡病人
驱动。（一个更早的方法——resample-then-max 后不重新居中——给出过 P=0.0045，
但该方法被证实无效：它重新估计了 winner's-curse 放大后的效应而非将其移除，
模拟显示其在真零假设下 Type I error≈0.36；已弃用，不再引用。）

---

## round-16：判定性结论——round 12–15 的"互补性"发现主要是一个 harness bug，
## 修正后配合全文一致的 BH 标准，结论不成立

**这是本轮（round-16）的判定性结果，取代下方 round-12 到 round-15 全部
关于"互补性"的具体数字结论——下方内容作为对账记录完整保留（不删改
历史），但阅读时请以本节为准。**

**round-16 指出三个比较对象都比它们各自包含的 PCA-200 分数本身丢失了
0.07–0.12 AUC，且丢失量不随特征数单调**（PCA-200 单独 1 特征丢
0.068，PCA-200+log1p(n_cells) 2 特征反而丢 0.119，比 13 特征的退化版本
丢得更多）——round-16 论证这在数学上不可能只由"过拟合"解释：单特征
逻辑回归是其输入的单调变换，重新拟合不应该丢这么多。**已用一个门控
单元测试直接核实，不是靠论证**（`round16_harness_unittest.py`）：把
PCA-200 自己的冻结 OOF 集成分数（原始 AUC=0.7423）套进 round-12–15 用的
同一套 harness（10-seed×5-fold、C=0.3、**未标准化**）重新拟合，OOF AUC
掉到 0.6569——与 round-15 报告的"PCA-200 单独"reduced AUC 精确一致，
确认复现的是同一个 harness、同一个 bug。**诊断（round-16，已核实）**：
PCA-200 的冻结分数是概率尺度（约 0.04–0.76，std=0.146），L2 penalty
在 C=0.3 下作用于**未标准化**的系数——一个只有这么窄尺度的单特征想要
产生足够的判别力就需要一个较大的系数，而较大系数恰恰被 L2 惩罚得最
狠，斜率被压得很扁，折间/seed 间的截距噪声便主导了 pooled OOF 排序。
**修正**（round-16 建议的方案二，标准做法）：`StandardScaler` 只在每折
训练集上 fit、同时变换训练与测试（无信息泄漏），单测里把丢失从 0.0854
降到 0.0182，通过门控阈值（<0.04）。

**用修正后的 harness 完整重跑三个比较对象 × 五个 FM**
（`round16_delong_scaled.py`，其余构造——冻结 representation 分数、
exact DeLong、10-seed×5-fold OOF——与 round-15 完全一致，只改了这一处
标准化）：

| 比较对象 | reduced OOF AUC | STATE Δ(p) | Geneformer Δ(p) | scPRINT-2 Δ(p) | scGPT Δ(p) | UCE Δ(p) |
|---|---|---|---|---|---|---|
| QC+队列+PCA-200（13 特征） | 0.7098 | +0.0244 (0.0305) | +0.0331 (0.0826) | −0.0017 (0.891) | +0.0004 (0.967) | −0.0002 (0.906) |
| PCA-200 单独（1 特征） | 0.7241 | +0.0517 (0.0454) | +0.0486 (0.0559) | +0.0192 (0.402) | +0.0097 (0.496) | −0.0073 (0.306) |
| PCA-200+log1p(n_cells)（2 特征） | 0.7131 | +0.0601 (0.0252) | +0.0564 (0.0448) | +0.0277 (0.264) | +0.0134 (0.361) | −0.0078 (0.271) |

三个 reduced OOF AUC 现在都落在 PCA-200 自身 OOF AUC（0.725，见上方主
结果表）附近的合理范围内，不再有 round-15 那种"比原始分数低 0.09–0.12"
的不合理跌幅——harness 修正后的行为符合直觉，是又一个确认修正正确的
信号。

**BH（α=0.05，五模型内部，与全文其余部分同一标准）：三个比较对象下
全部 0/5 通过——逐一核对每个排名，不是"第一名不过就自动全不过"这种
不成立的一般性说法**（BH 是 step-up 程序：从最大排名往回找第一个
p(i)≤(i/5)×0.05 成立的排名 i，拒绝 1..i；排名 1 不过不代表更靠后的
排名不能过，round-17 item 3 举了一个 rank-1 不过但 BH 仍拒绝全部 5 个
的反例，提醒不要把"排名1不过"当通用性质来讲）：

| 比较对象 | rank1 p (阈值0.01) | rank2 p (阈值0.02) | rank3 p (阈值0.03) | rank4 p (阈值0.04) | rank5 p (阈值0.05) |
|---|---|---|---|---|---|
| QC+队列+PCA-200 | 0.0305 fail | 0.0826 fail | 0.8912 fail | 0.9061 fail | 0.9669 fail |
| PCA-200 单独 | 0.0454 fail | 0.0559 fail | 0.3058 fail | 0.4020 fail | 0.4964 fail |
| PCA-200+log1p(n_cells) | 0.0252 fail | 0.0448 fail | 0.2636 fail | 0.2710 fail | 0.3608 fail |

三个比较对象里每一个排名在每一处都不满足阈值，逐一核对后确认 0/5 通过
——这次确实是全部排名都失败，不是靠"排名1失败即全部失败"这条不成立
的捷径推出来的。未经多重校正的名义显著性也只有 STATE 在三个比较
对象里都低于 0.05（p=0.025–0.045，边缘）、Geneformer 在一个比较对象
里低于 0.05（p=0.045，另两个 0.056/0.083，边缘不显著）；scPRINT-2、
scGPT、UCE 在任何比较对象下都不显著（p=0.26–0.97）。

**结论：round-12 到 round-15 建立起来的"互补性"发现——"全部五个模型"
（round-13/14）、"至少三个稳健显著"（round-15）——主要是同一个 harness
bug 的产物，不是真实的生物学信号。** 用修正后的 harness、配合全文
一贯的 BH 标准，这个额外的互补性问题回到了 null：唯一残留的信号是
STATE 在三个比较对象下都邊缘、未经多重校正显著（p=0.025–0.045），这
达不到本文任何其他地方采用的显著性标准，不能作为一个发现来报告。

**deviance permutation 检验（round-12 原始的 winner-corrected
P<0.0001）如何处理**：该检验从未经过这次发现问题的 harness（它是单次
in-sample 拟合，不涉及 10-seed×5-fold OOF pooling，不受同一 bug 影响），
但 round-15 已经指出它与 genuine OOF DeLong 结果在 UCE/scGPT 上不一致
（in-sample"显著" vs. OOF"不显著"），而现在 OOF 侧的结论本身又经历了
一次实质性下修——round-16 建议要么把 deviance 检验重新做成 OOF 版本，
要么按 round-9 处理已撤回的六模型 max-statistic 的先例，将其降级为
补充交叉检验、不再承担标题结论。**采纳降级方案**（原因同 round-9：
它现在冗余于一个已经算清楚的、更值得信任的检验族——这里是"互补性回到
null"这个结论本身——继续用一个从未 OOF 验证过的构造去支撑一个标题级
的 P<0.0001，是本文档在其余每一处都不会接受的做法）：winner-corrected
deviance P<0.0001 保留在下方 round-12 的历史记录中，标注为**未经 OOF
验证的补充数字，不再是任何标题/摘要/Key Point 的依据**。

**"in discrimination"限定词随之收回**：这一限定词的存在理由是"判别力
未显著，但互补性显著"——既然互补性本身在修正后回到 null，就不再有
需要与判别力分开限定的第二个结论。摘要/§2/§3a/§5/Key Point 1 的头条
表述**恢复为不带限定词的"no foundation model significantly
outperformed the baseline"**，与 round-12 之前的原始表述一致——不是
因为撤回了 round-12 的动机（那个动机——检验力不对称——依然是个真问题），
而是因为认真核实后，检验力不对称背后并没有一个真实存在的互补信号，
只有一个此前未被发现的 harness bug。

**新增 Methods QC 标准检查（第四个，与 round-9 的 floor-check invariant
并列）**：任何"reduced/full 嵌套对比"类型的 **OOF pooling** harness，
投入使用前必须先过一个门控单元测试——把比较对象自身包含的某个已知
分数单独作为 frozen 输入重新跑一遍这套 harness，其 OOF AUC 必须在原
分数 AUC 的 ~0.02–0.04 以内，不过这一测试的 harness 不得产出任何
显著性结论。`round16_harness_unittest.py` 保留作为这一检查的参考实现。

**round-17 item 1 明确范围：这一门控只适用于跨折 pooling 的 OOF
harness，不适用于本文其余的 in-sample 嵌套检验（round-8/9/10 的 QC
增量价值 nested-LR 族，含预声明主检验）。** 原因已用真实计算核实
（见上方"存在性 vs 幅度"段）：in-sample 构造里，观测统计量与它自己的
置换零分布用同一套代码拟合，未标准化造成的扭曲对两边对称地起作用，
经验置换 P 因而不受影响（本轮诊断中标准化前后经验 P 完全相同，都是
0.0005）——扭曲只让检验偏保守，不产出错误的"显著"结论。OOF pooling
harness 不满足这个对称性：它比较的是一个绝对 AUC **数值**（比如
"reduced=0.635"），这个数值会被单独拿去和文档其他地方未受扭曲的
参照点（比如 PCA-200 自身的 0.725）比较，扭曲因而直接污染了这个
比较，没有对称性可以吸收它——这是两类 harness 之间真实存在、已验证
的差别，不是回避核查的借口。

---

**以下 round-12 至 round-15 的历史记录原样保留，供对账追溯使用，
其具体数字结论已被上方 round-16 取代，不再代表本文档的当前立场：**

**round-12 item 1（round-13 重新框定）：头条对比换一个统计量后结论翻转为
显著——已核实为真实发现，不是 bug，与 ΔAUC 头条并列报告，ΔAUC 仍是
预声明主检验。** **round-13 修正框架**：这两个检验回答的不是"同一个问题、
两种检验力"，而是**两个不同的问题**——ΔAUC winner-corrected permutation
问的是"最佳 FM 的分数本身是否比 PCA-200 的分数排序更好"（**优越性**，
superiority）；deviance 检验问的是"在 PCA-200 的分数之上，FM 的分数是否
还**额外贡献**信息"（**互补性**，complementarity）。两个高度相关
（r=0.71–0.90）的表示完全可以互补而不优越——FM 可以携带 PCA-200 没有的
信息，同时不在独立排序上赢过 PCA-200。在这个框架下两个结果都是对的，
不需要"谁更敏感、谁更接近真相"的论证；这也解释了一个此前框架下显得
奇怪之处——如果 deviance 只是同一件事的更强版本，审稿人会问为什么论文
报告的是检验力更弱的那个作为 primary；在优越性/互补性框架下答案是
现成的：判别力是部署相关的量，互补性不是——这正是 §2 早已给出的答案。

**用同一套 winner-corrected 置换机制**（`round12_headline_deviance_test.py`，
job 28229781，**4 个 array task 各 2,500 次、共 10,000 次患者级标签置换，
置换零分布 = 每次置换下 5 个 FM 里的最大统计量，跨 4 个 task 合并——
"10,000 次"指的是 10,000 个 per-permutation 最大值，不是 5 个 FM 的统计量
被合并成一个池；round-13 item 4 要求明确这一点，因为两次更早的
max-statistic 构造出过错**）**重新检验头条对比，统计量换成 deviance**：
reduced 模型 = QC(3)+log1p(n_cells)+cohort one-hot(8)+PCA-200 冻结 OOF
分数（13 特征，C=0.3），full 模型 = reduced + 某 FM 的冻结 OOF 分数
（14 特征）；per-FM 统计量 = 2×(ll_full−ll_reduced)；winner 统计量 = 5 个
FM 里的最大值，PCA-200 本身在 reduced 模型内、不作为候选。

**round-13 item 3：五个模型逐一对照置换零分布的百分位（中位数 0.188、
95th pct 0.790、99th pct 1.230、最大值 3.819），不只报告 argmax。**
winner-corrected 的 max-of-five 零分布随机占优于任一单模型自己的边际零
分布，因此下表的相对位置是保守估计（每个模型自己的边际 p 会更小）：

| 模型 | deviance 统计量 | 相对于 winner 零分布的位置 |
|---|---|---|
| STATE | 4.3284 | 超过零分布最大值（3.819） |
| Geneformer | 3.2227 | 超过 99th pct（1.230） |
| scPRINT-2 | 3.1843 | 超过 99th pct |
| scGPT | 1.5365 | 超过 99th pct |
| UCE | 0.8749 | 超过 95th pct（0.790），低于 99th pct |

**结论改写（round-13）：不是"STATE 的表示携带额外信息"，而是"全部五个
被测表示都在 PCA-200+队列之上携带可检测的额外信息"**——单独归因于 STATE
会与 §3a 上方"STATE 与 Geneformer 不可区分（p=0.60）"放在一起显得矛盾，
已改正为面向全部五者的陈述；winner-corrected P<0.0001（"至少一个表示
显著"，只需要 argmax=STATE 一个成立）保留作为标题推断（headline
inference），不受下方任何修正影响。**round-15 更新："全部五个"这一具体
措辞需要修正**：上表基于的是 in-sample（C=0.3 惩罚 log-likelihood）
deviance 统计量与其自身置换零分布的百分位比较，而下方 exact DeLong 是
genuine OOF 的直接患者级检验——两者对 UCE（95th pct vs p=0.18）与
scGPT（99th pct vs p=0.06）的结论不一致：in-sample 版本说"显著"，OOF
版本说"不显著"。**这是 in-sample 统计量比 OOF 更容易呈现虚高显著性的
又一实例**（与 round-14 自己最初的直觉方向一致，只是这次影响的是显著性
判定本身，不是幅度大小）——deviance 置换检验的观测统计量与零分布**都**
是 in-sample 构造，从未像下方 ΔAUC 幅度那样被重新用 OOF 版本验证过；
本轮不新增这项计算（不在 round-15 明确要求的范围内），如实标注为留待
后续处理的未闭合项，不影响 winner-corrected P<0.0001 这一"至少一个显著"
的结论（因为下方 OOF DeLong 独立确认了 STATE 本身稳健显著），但"全部
五个都携带可检测信息"这一更强表述应收窄为"至少三个（STATE、
Geneformer、scPRINT-2）稳健显著，一个（scGPT）在更干净的比较对象下
显著，UCE 在任何构造下都不显著"——与下方 DeLong 结果一致，见下段。

**round-13 item 2 起先只算出 in-sample 幅度，round-14 item 2 要求换成
genuine OOF——已重算，数字发生了非预期方向的变化，如实报告。**
存在性（deviance 显著）不等于幅度（判别力提升多少），此前只是断言
"提升量小到不足以移动判别力"而没有算出真实数字；round-13 用同一批
reduced/full 拟合算出了 in-sample nested ΔAUC（STATE +0.0279 等），但
in-sample 值在全文语境下是未加限定的例外——§2 明确禁止裸 "AUC"，全篇
其余每一个 AUC 都是 OOF，round-14 指出这一个不是、且摘要里没有标注这一
点，是审稿人会重点核对的地方。

**round-14 曾用配对 Wilcoxon（跨 10 个 CV seed）作 OOF 显著性证据——
round-15 item 1 指出这是错的检验，与全文早前 round-3 已经排除的那类
错误同形（seed 间±sd 反映的是同一批 240 人在不同 CV 划分下的分区噪声，
不是病人抽样层面的不确定性，round-3 已因此把"跨 10-seed 的 sd"排除
出误差项候选）。10-seed 配对 Wilcoxon 的可能最小 p 值是 2/2¹⁰≈0.00195，
五个模型全部报告这个值，说明这个统计量只记录了一件事——十次 CV 划分
方向一致——不携带任何关于病人抽样层面可分辨性的信息，既不能当第二套
独立显著性证据，也不能当"此对比在此样本量下可分辨"的经验证明。**已按
§2 早已声明的 pairwise SE 标准改正**：对 reduced/full 的 seed-ensemble
OOF 预测向量（患者级、10 个 seed 取均值聚合为一个分数，与 model_vecs 的
构造惯例一致）跑 exact DeLong（`round15_delong_and_undamaged.py`——
项目里此前生成 `incremental_value_test.csv` 的 exact DeLong 实现已不在
盘上，本脚本原样复用 `scripts/24_option_c_canonical_input` 里保留的那份
实现，运行前先用两项独立核验确认正确：AUC 点估计与 sklearn 精确一致；
DeLong SE 与患者级 paired bootstrap SE 高度吻合，ratio=0.997）：

| 模型 | full OOF AUC | reduced OOF AUC | ΔAUC | DeLong SE | z | DeLong p | 预测相关 r |
|---|---|---|---|---|---|---|---|
| STATE | 0.6716 | 0.6348 | +0.0368 | 0.0120 | 3.069 | **0.00215** | 0.979 |
| Geneformer | 0.6646 | 0.6348 | +0.0298 | 0.0095 | 3.121 | **0.00180** | 0.987 |
| scPRINT-2 | 0.6633 | 0.6348 | +0.0285 | 0.0119 | 2.390 | **0.01684** | 0.979 |
| scGPT | 0.6482 | 0.6348 | +0.0134 | 0.0072 | 1.854 | 0.06367 | 0.992 |
| UCE | 0.6437 | 0.6348 | +0.0089 | 0.0067 | 1.338 | 0.18083 | 0.993 |

**BH（α=0.05，五模型内部）：3/5 通过**（Geneformer/STATE/scPRINT-2，
排序后逐一对比阈值 [0.01,0.02,0.03,0.04,0.05] 全部通过；scGPT 0.06367
在原始与 BH 口径下均不过；UCE 0.18083 同样不过）——不是 round-14 曾报告
的"全部五个显著"，是三个稳健显著、一个边缘不显著（scGPT）、一个明确
不显著（UCE）。**测得的预测相关 r=0.979–0.993**，落在 round-14 模拟表格
假设的 rho=0.97–0.99 区间内，二者互相印证：嵌套对比的真实相关性确实
在这个高区间，不是凭空假设。

**round-15 item 2：更根本的问题——reduced 比较对象本身严重受损
（QC+队列+PCA-200 的 OOF AUC=0.6348，比 PCA-200 单独的 OOF AUC 低
~0.09–0.11），"增量"有多少是修复这一损伤、而非真实互补信息？已用一个
判别实验直接检验，不是靠论证。** 8 个 cohort dummy 在仅 25 个事件下的
过拟合倾向此前已确认（边际贡献仅 +0.009），round-15 指出这解释了
round-14 自己"OOF 值反而大于 in-sample"这一非预期发现的确切方向：
in-sample 时 reduced 模型未受损（没有可修复的），增量小；OOF 时 reduced
模型受损（有大量可修复的），增量变大——这正是"修复"机制会预测的模式，
把它从一个有趣的巧合变成了对增量解读的真实威胁。**判别实验**：换用
两个未受损的 reduced 比较对象重跑同一构造（representation 分数依旧
冻结，同一 10-seed×5-fold OOF 协议，exact DeLong）——

| reduced 版本 | reduced OOF AUC | STATE Δ (p) | Geneformer Δ (p) | scPRINT-2 Δ (p) | scGPT Δ (p) | UCE Δ (p) |
|---|---|---|---|---|---|---|
| PCA-200 单独（1 特征） | 0.6569 | +0.0878 (0.00007) | +0.0690 (0.00085) | +0.0659 (0.01661) | +0.0346 (0.01398) | +0.0123 (0.51866) |
| PCA-200+log1p(n_cells)（2 特征） | 0.6063 | +0.0753 (0.00066) | +0.0566 (0.00032) | +0.0638 (0.00542) | +0.0214 (0.01943) | +0.0171 (0.08083) |

**结果与"修复"假说的方向相反——增量在两个更干净的对比对象下不但没有
缩小，反而变大**（STATE 从退化对比的 +0.0368 涨到 +0.0878/+0.0753），
且新增 scGPT 在两个未受损版本下都转为显著（p=0.014/0.019）——若增量
主要来自修复 8 个 cohort dummy 的过拟合，换成一个只有 1–2 个干净特征、
几乎没有过拟合可修复的 reduced 模型后，增量应该收缩趋零，而不是变大。
**三个版本汇总：STATE/Geneformer/scPRINT-2 在全部三个 reduced 版本下都
显著；scGPT 在两个未受损版本下显著、在退化版本下边缘不显著（p=0.064）；
UCE 在全部三个版本下都不显著**（p=0.081–0.519）——这是一个经受住压力
测试的、内部一致的结论：至少三个、很可能四个 FM 携带真实的互补信息，
UCE 始终是唯一的例外，与全文其余部分对 UCE 的定位完全一致，不是修复
产生的假象。

**round-15 item 3：deviance 检验与 DeLong 检验共享同一个（此前已知受损
的）reduced 模型，"方法学独立"这一措辞不准确，已软化。** 二者在**统计量**
上独立，但共享**比较对象**——一个损坏的 reduced 模型会同时抬高两者的
似然比与 OOF AUC。准确的表述是：两个检验方向一致，且共享同一个比较
对象，因此真正检验结论的不是二者的一致性，而是上方的未受损比较对象
敏感性分析——分析结果支持结论成立，但支持它的证据是敏感性分析本身，
不是"两种独立方法都同意"这句话。

**round-14 item 1：分辨率对比此前用错了标尺，且原句在字面上自我矛盾，
已改正（结论不变，仅证据来源从已撤回的 Wilcoxon 换成上方 DeLong 结果）。**
全文引用的"~0.10 AUC @25 events"分辨率来自**两个不同模型**之间对比
（如主对比 Geneformer vs PCA-200，分数相关 r≈0.8），而这里的"nested"
对比——reduced 与 full 模型共享几乎全部特征——预测相关性远高于两个独立
模型（测得 r=0.979–0.993，见上表），方差因而小得多，可分辨的最小效应量
也小得多；用两模型对比的分辨率下限去衡量一个嵌套对比的效应量，是错的
标尺。**原句更严重的问题是自我矛盾**：这个增量已经被检测出来了（无论是
deviance 的 P<0.0001 还是上方 DeLong 的 p=0.002–0.017），"检测到的东西
低于可分辨下限"在逻辑上不能同时成立。**直接经验证据**：三个模型在
patient-level exact DeLong 下稳健显著、一个在未受损比较对象下显著，
这就是"这一嵌套对比在此样本量下可分辨"的直接证明。正确的表述：**嵌套
对比与其比较对象共享几乎全部模型，预测高度相关（r≈0.98–0.99），因而在
这一样本量下是可分辨的；标准对比（两个独立表示）预测相关性低得多
（r≈0.8），需要 ~0.10 AUC 才可分辨——这正是同一批数据能同时支持"探测到
显著互补信号"与"无法确立显著优越性"这两个结论的机制解释，不是张力**。

**"no foundation model significantly outperformed the baseline"这一
贯穿全文/摘要/Key Points 的表述需要"in discrimination"这一限定词**——
判别力（部署相关的量）上未显著超过，不代表表示内容上没有可检测的互补
信息；两者都如实报告，不取更有利的一个，与 5 个 FM 增量价值两种口径的
处理方式一致。

**round-13 item 4：floor-check invariant 在此处空转，如实标注。**
Methods QC 标准检查（`assert corrected_p >= min(uncorrected_ps)`）应用
于本检验时无信息——winner-corrected P 与 STATE 自身在这一 family 里的
位置都撞在 10,000 次置换的分辨率下限（0/10,000 次超出），检查在此分辨
率下恒为通过，不构成对这一具体构造的验证；已在下方 Methods QC 段单独
标注，不与另外三个真正起到区分作用的案例并列。

**两两排名本身仍有信息（exploratory）**：pairwise SE 以 exact DeLong 为准
（与 bootstrap 高度吻合，未系统性偏差）：Geneformer−PCA p=0.036、
STATE−PCA p=0.048，均不过 winner 校正；scPRINT-2/scGPT 与 PCA 不可分
（p=0.65/0.91）；UCE 低于 PCA 不显著（p=0.081）。

**唯一在两种方差估计下都成立的分离：双雄 vs UCE。** 差距 +0.118/+0.130
（exact DeLong z≈4.3–5.0，p<0.0001）。已升级为 **confirmatory**：
range-statistic permutation（与主检验同批标签置换，nperm=10,000，
rng_seed=20260908）——观测 range=0.1302，**P=0.0112**；同时通过 15 对
pairwise 的 BH-FDR（q=0.030）。分数相关 r=0.71–0.90；consensus（5-FM 均值）
AUC=0.768，描述性，不作推荐。

**各方法 vs null（count+disease）配对检验（exact DeLong）**：STATE
p=0.00135、Geneformer p=0.00063、scPRINT-2 p=0.0113、scGPT p=0.0232、
PCA-200 p=0.0126 均显著；**UCE p=0.157 不显著**——这是 UCE 未通过的三项
检验之一（另两项见下方增量价值段：AUC-赛跑口径 p=0.072 未通过，
round-10 修正后的 nested-LR 口径 p=0.0022 **勉强通过**，是六者中唯一
未撞到置换分辨率下限的）。"UCE 在任何检验中都未显著超过 null"这一表述
（round-9 及更早版本）经 round-10 修正后不再准确，收窄为"UCE 在本文
全部相关检验中始终是表现最弱、最常未过检验的方法"，Discussion 需按此
措辞陈述。

**QC 增量价值检验（结论：非假象，通过增量价值检验得出；round-8 机制定案：
0.376 是抑制模式，不是折不稳，且与 §3b 的 sepsis PCA-200=0.428 不是同一
机制）。** QC-only（仅 3 特征，不含 count）as-fit AUC=0.376，是随机的错误
一侧 0.124；QC 3 特征**条件于** library depth（加 log1p(n_cells)）后
as-fit AUC=0.618，方向正确，比"仅 depth"高 +0.064。这是教科书式的统计
**抑制（suppression）模式**——%mito 与基因数这类 QC 特征只在**条件于**
测序深度后才携带结局信息方向正确的信号，这在生物学上不意外（两者本身与
depth 机械相关）。**这是该基线特有、已定案的机制，不再与 §3b 的 sepsis
PCA-200=0.428 并作"两个低于随机的 out-of-fold AUC 是模式而非巧合"的联合
陈述**——那条联合陈述已撤回，0.428 的机制未知、独立处理。去掉 8 个
cohort dummy 单独看 QC+log1p(n_cells)：AUC=0.618，加回 cohort 后
0.618→0.627，**cohort 变量对方向性贡献仅 +0.009**——cohort 组成不是
0.376→0.618 这一步的机制，抑制模式（depth 条件化）才是。properly-specified
增量价值检验（基线 = QC 3 特征
+ log1p(n_cells) + cohort one-hot [8 dummy, Data_NO]，L2 正则 C=0.3，同
10-seed×5-fold 患者分组协议）：基线 **AUC=0.627，as-fit，从未做过符号翻转**
——与 Job K 最初那个更简单的 QC-only-无-cohort 模型（0.376→0.624 才是
真的翻了符号）是两个不同的基线，此处不存在"看到结果朝哪个方向倒再决定
方向"的问题。**已核实（round-7 诊断）**：4 个协变量的单变量 Spearman 关联
全部正向（log_counts ρ=+0.092、log_genes ρ=+0.095、pct_mito ρ=+0.064、
log_n_cells ρ=+0.076，n=240 下均未达名义显著）；50 次拟合中多变量系数
方向与单变量一致（log_counts/log_n_cells 100% 为正，log_genes 72%、
pct_mito 64% 为正）——**无拟合/正则化伪影**，符号本就一致，模型是按其
应有方向工作的（cohort dummy 的贡献已在上方基线表脚注给出）。

**预声明增量价值主检验（round-8 item 1 定案选择的检验，round-10 item 1
修正参考分布后的数字）：PCA-200 单独一项，nested likelihood-ratio 检验
（协变量-only vs 协变量+PCA-200，同折拟合），P<0.0001（10,000 次经验
置换零分布，0 次超出观测统计量 4.4772；旧的 chi²(1) 渐近参考给出的
p=0.0344 已核实是参考分布失配导致的过度保守估计，详见下方 round-10 段）。**
选它为主检验的理由不变，与是否达到显著无关：(a) 不涉及模型选择——PCA-200
只有一个，不像 5 个 FM 需要挑；(b) 它就是被问的问题本身："在 QC+计数+队列
之上，platelet 表示是否提供增量信息"；(c) 作为唯一预声明检验，不需要多重
校正，单独站得住。**"learnable from platelet transcriptomes...非 QC/队列
假象"的结论引用这一个预声明检验，不再引用六个未校正检验的并集。**

**5 个 FM 各自的增量价值数字仍为 exploratory（不预声明，两种口径如实并录，
不取更有利的一个）；round-11 item 2 要求两种口径必须成对报告、不得单独
引用某一列的 P 值——理由见下方"存在性 vs 幅度"段。** 完整六模型对照表
（AUC 增量取自 `incremental_value_test.csv`，与上方 §3a 基线表 0.627 一致；
nested-LR P 取自 round-10 修正后的经验置换零分布，`round10_nested_lr_
permutation_null.py`，job 28227832）：

| 模型 | AUC 增量（vs 0.627 基线） | exact DeLong p | nested-LR P（round-10 修正后） |
|---|---|---|---|
| STATE | +0.065 | 0.00015 | P<0.0001 |
| Geneformer | +0.053 | 0.00014 | P<0.0001 |
| PCA-200 | +0.050 | 0.0034 | P<0.0001 |
| scPRINT-2 | +0.046 | 0.0054 | P<0.0001 |
| scGPT | +0.026 | 0.0063 | P<0.0001 |
| UCE | +0.016 | 0.072（不显著） | P=0.0022 |

**置换 p 值惯例（round-11 item 2 要求明确陈述）**：10,000 次置换、0 次
超出观测统计量，可报告的界是估计器本身的分辨率下限 P<1/10,001——本文档
统一写作"P<0.0001"，指的是这个下限，不是一个算出的具体值；更高精度需要
更多置换次数，本节结论（6/6 通过 BH）不依赖这一位数的精确度。

**存在性 vs 幅度（round-11 item 2 核心提醒，防止读者误读；round-17 item 1
修正了机制归因，用真实计算而非论证）**：nested-LR 检验之所以把 5 个
模型都压到置换分辨率下限，round-10 当时把机制归因为笼统的"C=0.3 的
L2 正则把一个真正无信息的新增系数压向零"。**round-16 在完全不同的
（OOF pooling）harness 上定位到了更精确的机制——L2 作用于一个
**未标准化的窄尺度概率特征**，而不是"L2 正则"本身——round-17 指出这
应该是同一个现象，并要求直接验证：标准化是否会让这里的坍缩零分布
显著移向 chi²(1)。** 已直接验证（`round17_pca200_primary_scaling_check.py`，
PCA-200 单独一项，2000 次置换）：**不标准化**（round-10 原始构造）——
观测 LR=4.4772（精确复现），置换零分布中位数=0.084（vs chi²(1)
中位数 0.455），渐近 chi²(1) p=0.0344（精确复现 round-10 报告的
"旧参考给出 p=0.0344"）；**标准化后**——观测 LR 跳到 18.2530（同一份
数据、同一个信号，标准化后系数不再被人为压扁），置换零分布中位数=
0.462，**几乎精确落在 chi²(1) 的 0.455 上**，渐近 chi²(1) p=0.00002。
**机制归因现在有直接实验支持，不再是论证**：确实是"L2 作用于未标准化
窄尺度特征"这一具体机制，不是"L2 正则"这个笼统说法。**但这不影响任何
已报告的显著性结论**——两个版本的经验置换 P 完全相同（2000 次置换下
都是 0.0005，撞在这一诊断跑的分辨率下限）：因为观测统计量与它自己的
置换零分布用的是同一套（可能被扭曲的）拟合代码，扭曲对两边对称地起
作用，参考分布把它"吸收"了，round-17 自己的判断——"这与 round-16 的
OOF 情形不同，那里的扭曲直接污染了一个被拿去和其他未受扭曲的参照点
比较的绝对 AUC 数值；这里扭曲只降低检验力（偏保守），检验本身仍然
显著"——已被证实成立。**预声明主检验（PCA-200 nested LR，10,000 次
置换 P<0.0001）不需要重跑替换，结论稳健；仅机制归因文字已按上述更新。**
nested-LR 仍是一个强有力的**存在性**检验（"增量信息是否存在"），但对
**幅度**（"增量信息有多大"）几乎没有区分力——**这正是本文其余部分建立
"分辨率"框架所针对的量**。把"P<0.0001"单独放在"+0.026 AUC"旁边会让读者
读出一个大效应；正确的分工陈述是：nested-LR 族确立了 platelet 表示在
QC+计数+队列之上确实携带增量信息这一**存在性结论**；DeLong AUC 赛跑把这
一信息**量化**在 +0.016 到 +0.065 之间；本文的分辨率论证（25 events 分辨
~0.10 AUC）适用于这个量级，不适用于 nested-LR 的 P 值本身——不能从
"P<0.0001"推出"增量信息很大"。

**BH 校正结论反转**：对六模型（5 FM+PCA-200）分别在两种口径内部做 BH
（α=0.05）——AUC 赛跑口径仍 5/6 通过（仅 UCE 不过，未变）；**nested LR
口径由"0/6 通过"改为"6/6 通过"**（六个 p 排序后逐一对比 BH 阈值
[0.00833, 0.01667, 0.025, 0.03333, 0.04167, 0.05]：五个 P<0.0001 全部
远低于各自阈值，UCE P=0.0022 亦低于其阈值 0.05，全部通过）。两种口径
现在仅在 UCE 一个模型上分歧（AUC-赛跑不显著 p=0.072，nested-LR 勉强显著
P=0.0022），不再是"是否存在增量价值"这一问题本身的分歧——这一变化的
完整诊断过程见下方 round-10 段。**round-17 item 1 明确标注：这里的
"6/6 通过"与后面 round-16 段的"0/5 通过"问的是两个不同的比较对象，不
矛盾，不是同一问题的两次计数**——这里（round-8/9/10 家族）问的是"表示
是否比裸 QC+计数+队列（不含 PCA-200）多提供信息"；round-16 家族问的是
"FM 是否比 PCA-200+队列（含 PCA-200）多提供信息"。合并读法见 §5：
platelet 表示携带超出 QC+计数+队列的信息，但没有一个基础模型表示携带
超出一个 from-scratch PCA 本身的信息。

**round-10 item 1：nested-LR 参考分布修正（重大更新，反转了上一段的
BH 结论）。** round-9 曾把"nested-LR 检验族的经验置换零分布比 chi²(1)
渐近参考更集中于零"记录为"不影响结论方向的 round-10 候选项"，理由是
"渐近 p 是更保守的方向"。**round-10 external review 指出这个方向判断是
反的**：保守（偏大）的 p 值只保护"显著"结论，不保护"不显著"结论——而
nested-LR 族里恰恰有 UCE（p=0.231）、scGPT（p=0.148）两个"不显著"判定，
以及依赖 0/6-vs-5/6 不对称性所建立的"需要单独预声明主检验"这一论证结构，
都是用渐近参考做出的。round-10 同时提出了一个具体机制猜测：置换循环里
代表分数（representation score）可能被在置换标签下重新拟合，导致该
特征在置换零假设下接近常数、统计量整体坍缩而非只是变窄。**独立代码级
核实（`round10_diagnose_null.py`）确认这个具体猜测是错的**：代表分数
向量在 2,000 次置换的诊断跑中，object id 与内容 hash 全程不变——它从
未被重新拟合，从一开始就是正确冻结的（这也是 round8_item1 脚本原本的
设计，round-8/9 都没有改动过这一点）。真正原因用 C=1e6（近乎无惩罚）
重跑同一构造后确认：**无惩罚版本的经验零分布紧贴 chi²(1)**（STATE：
中位数 0.468 vs chi²(1) 中位数 0.455；95th pct 4.10 vs 3.841），而
**C=0.3（全文统一使用的 as-fit 惩罚强度）版本的经验零分布明显更集中于
零**（STATE：中位数 0.082，95th pct 0.687）——chi²(1) 本来就不是一个
L2 正则化统计量的正确渐近参考，这是一个参考分布失配问题，不是置换实现
的 bug。**修正方案（在看到全部结果前声明，round-10 明确要求）**：不
改变模型设定本身（继续用全文统一的 C=0.3 as-fit 惩罚，不像 round-10
建议的备选方案那样去掉惩罚——去惩罚会引入一个新的、事后才决定的建模
选择，改参考分布不会），只把参考分布从 chi²(1) 换成该惩罚统计量自身
的经验置换零分布（label permutation，频率与全文其余两个置换检验一致：
10,000 次/模型）。六模型的完整重算结果与由此导致的 BH 结论反转见上段。

**关于 UCE 的连带更新**：UCE 现在在 round-10 修正后的 nested-LR 检验中
以全场最弱的优势通过（P=0.0022 vs 其余五者的 P<0.0001，观测统计量
1.43 vs 2.09–5.76），但在同一增量价值问题的 AUC-赛跑口径下仍不显著
（p=0.072），在完全独立的 vs-null exact DeLong 检验中也仍不显著
（p=0.157）。**"UCE 在任何检验中都未显著超过 null"这一贯穿 §3a/§5/
abstract 的表述不再准确**，已在上方"各方法 vs null"段收窄措辞；双雄
vs UCE 的 confirmatory 分离（range-statistic P=0.0112）本身不受
round-10 影响（该检验与 nested-LR 族完全独立，未使用 chi²(1) 参考）。

**round-8 的六模型 max-statistic permutation 已撤回（round-9 item 1）**：
该段原报告置换 P=0.0001（`round8_item1_max_lr_permutation.py`，HPCC job
28214868），round-9 外部复核指出这违反 selection correction 的数学下界——
被选择的六个 uncorrected p 中最小者是 STATE 的渐近 p=0.0164，corrected P
不可能比这更小。独立诊断复现（state_env，重构 STATE 与 PCA-200 各自的
边际置换零分布，300 perms，两者观测值均可精确复现：STATE LR_stat=5.7596/
渐近 p=0.0164、PCA-200 LR_stat=4.4772/渐近 p=0.0344）确认：(a) 原脚本的
六模型 max-statistic 循环存在真实、未定位到具体代码行的 bug——即便用两个
argmax 候选各自的经验置换边际 p（≈0.003，300 perms 下 0 次超出，
1/301）作为更宽松的下界，0.0001 仍违反下界约 30 倍；(b) 附带发现、留待
后续 round 处理：L2 正则（C=0.3）下 nested-LR 统计量的经验置换零分布
显著比 chi2(1) 渐近参考更集中于零（STATE 置换零分布中位数 0.083、95th
pct 0.739，对比 chi2(1) 参考中位数 0.455、95th pct 3.841，约小 5 倍）——
意味着本节引用的渐近 p（STATE 0.0164、PCA-200 0.0344 等）本身可能不是
六个 nested-LR 检验的正确参考分布，经验置换很可能给出更小、而非更大的
p 值；**这一判断本身是对的，但"不影响不显著结论的方向"这一判断被
round-10 指出是反的（保守 p 只保护显著结论）——已在下方"round-10 item 1"
段完整修正并重算全部六个模型，不再是未解决项。**

依据 round-9 review 的明确建议（"其冗余于已预声明的主检验，且即便正确
计算出的值也只是 marginal（~0.04–0.06），不如直接删除"），**整段撤回、
不做"修正后保留"**——"learnable from platelet transcriptomes"的结论
继续且仅由上方 PCA-200 单项预声明检验（round-10 修正后 P<0.0001）支撑。
（round-9 对"即便修正也只是 marginal ~0.04–0.06"的估计后来被 round-10
的正确参考分布证明过于保守——真实值远小于此，但这不改变"整段撤回"这个
决定本身：即便一个正确构造的 max-statistic 版本会显著，它仍然冗余于
已经显著的预声明单项主检验，不需要为此恢复那段代码。）

**Methods QC 标准检查（round-9 item 2 新增，一行代码）**：任何
selection-corrected p 值必须 ≥ 其所在 family 中最小的 uncorrected p——
`assert corrected_p >= min(uncorrected_ps)`。retroactive 核对全文三个
selection-corrected 统计量：主检验（uncorrected 0.036 → corrected
0.181，通过）、range statistic（uncorrected <0.0001 → corrected 0.0112，
通过）、六模型 max-statistic nested-LR（uncorrected 0.0164 → corrected
0.0001，**违反，已撤回**，见上）——违规仅此一处，已按上述处理。
**round-12 头条-deviance winner-corrected 检验（round-13 item 4 标注）：
该检查在此处空转，不计入"通过"案例。** winner-corrected P 与 STATE 自身
在其 family 里的 uncorrected 位置都撞在 10,000 次置换的分辨率下限
（0/10,000 次超出，两者都只能报告 P<0.0001），检查因而恒真，不构成对
这一具体构造是否正确的验证——不同于上面三个案例，这里"通过"不提供
信息，如实标注为空转，而非与真正起到过区分作用的案例并列计数。

**QC 附带检查（与上述撤回无关，独立成立予以保留）**：按 cohort 的
QC-结局关联多数 p>0.05（1/6 可测 cohort 达名义 p=0.031，未过多重校正）；
折间 cohort 组成均衡（`fold_cohort_balance_seed42.csv`）。

**cohort-only 基线的诚实计算（机制已定案：add-one 平滑的规模反转，非
CV 噪声稀释）。** 一版手算曾用 LOBO 折规则估计 cohort-only 预测器下限
≥0.579（"34 名排除者全部得最低分"的简化 concordance 论证，后证实基于
未做交叉验证的 in-sample 假设，已弃用）。**实际算出来的、与全文同协议**
（10-seed×5-fold patient-grouped，折内 Laplace-smoothed cohort mortality
rate 作为唯一特征，out-of-fold 评分）**的诚实数字是 0.486±0.038**——接近
随机。**机制**：Laplace add-one 平滑 (d+1)/(n+2) 按 cohort 规模系统性
反转风险排序——一个 6 人零死亡 cohort 得分 0.125，高于一个 90 人、8 例
真实死亡的 cohort（0.098）（完整 9-cohort n/死亡/平滑分见
`item4_per_cohort_table.csv`）。这是**该特定估计器的构造性缺陷**，不是
"CV 噪声把 cohort 间方差稀释掉了"。**0.486 因而对"信号是否为 cohort
假象"这一问题保持沉默，既不支持也不反驳**——真正回答该问题的证据是
带 cohort one-hot 的增量价值检验（见下段，结论仍然成立：控制 cohort 后
platelet 表示仍普遍提供显著增量）。

### 3b. within-sepsis（91 patients / 17 fatal，GSE279452，28 天死亡）

UCE 0.548±0.038、scPRINT-2 0.523±0.049、Geneformer 0.515±0.043、
STATE 0.511±0.040、scGPT 0.449±0.052、PCA-200 0.428±0.082。**完全 null**：
患者级配对检验 + 多重校正后无任何两两对比显著（唯一名义对比 UCE>scGPT
p=0.047，不过校正；所有 FM vs PCA 的 95% CI 均含 0，最小名义 p=0.079 为
UCE）——不仅全部 ≈ 随机，任何模型间排序在此样本量下都不可辩护。**功效
边界**：17 fatal events 下 80% 功效仅能分辨 AUC≳0.71，此 null 是"小信号
上界"，不是"无信号证明"。文献（Task 6a）：sepsis 血小板结局信号此前仅
报道于早期窗（0–6h / <72h）——我们是晚时点（中位 16 天）的首次测量；早晚
窗对比是研究间比较（不同队列/平台），不是同一队列内信号消失。

**年龄结构**（数据定案，精确核对到 91-patient 分析集本身）：91-patient
within-sepsis 分析集中 90 人有年龄记录（1 人 `Skin-S152` 未匹配到年龄），
**全部已知年龄者均为成人**（19–93 岁，中位 67，零儿科）。93 是 GSE279452
GEO "Adult Sepsis" 候选组原始人数（91 = 93 减 2 个 mort28 缺失），92 是该
候选组中成功匹配到年龄的人数（1 人年龄未匹配，与上述 91 分析集缺年龄的
同一人）——93/92/91 三个数字指代不同层级，非笔误。GSE279452 全系列另有
26 个儿科样本（0.1–15 岁），**无一进入 93-candidate 或 91-analysis 任一
层级**。

### 3c. monocyte 阳性对照（in-distribution，148 patients / 19 fatal）

STATE 0.696、scGPT 0.677、Geneformer 0.674、scPRINT-2 0.672、UCE 0.648、
PCA-200 0.654。**无任何 FM 显著赢 PCA**（STATE +0.043 p=0.23、Geneformer
p=0.57、scPRINT-2 p=0.68、UCE p=0.80、scGPT p=0.94）——in-distribution 上
from-scratch 基线统计上就够用。mono < platelet 倒挂（全体，含 PCA −0.072）
经 COVID-only 子集确认成立 → 任务属性，非 FM 属性。唯一方向性例外：UCE 是
全部四种配置（seed × disease-pool）下唯一偏好 monocyte 的模型（对 PCA 仍
不显著）。

**与主结果一致**：pipeline 在任何设计格（含 within-COVID 主格、含
in-distribution 的 monocyte 对照）都未分辨出 FM 优势——platelet 的 null
不能归咎于 OOD，因为 in-distribution 同样 null；该对照也无法证明 pipeline
在任何地方能分辨 FM 优势。

### 3d. cross-disease（patient-level，severe=survived 主协议）

- **COVID→Sepsis（测试池 = 91 GSE279452 patients/17 fatal）**：全部
  0.64–0.68，CI 全重叠（UCE-LR 0.684 [0.543, 0.805]；STATE 0.653；
  scGPT 0.650；scPRINT-2 0.647；Geneformer 0.644）。**这不是"转移弱"，是
  "转移不可测"**：§3b 已证明该 91-patient 池内死亡率即便用同域训练也不可
  学习（全部方法≈随机、无任何配对对比显著、17 fatal events 功效边界
  AUC≳0.71）——把模型迁移到一个本身检测不出信号的任务上，测到的还是那个
  任务的 null，不是迁移能力的证据。这同时解释了正向/逆向不对称：逆向在
  信号确实存在的 COVID 队列上测（可达 0.79），正向在什么都学不到的队列上
  测——不对称是两个目标任务可学习性不同的属性，不是迁移方向本身的属性。
- **Sepsis→COVID（训练池 = 96 patients = 91 GSE279452 + 5 GSE167363；
  测试 = 240 patients/25 fatal）**：点估计整体更高（STATE-LR 0.7209
  [0.621, 0.809]、RF 0.7113；Geneformer-LR 0.7118 [0.612, 0.814]、
  RF 0.6358；UCE 0.693/0.607、scGPT 0.642/0.622、scPRINT-2 0.619/0.517），
  CI 全重叠。**排名不可辩护**：severe-excluded 敏感性下头两名反转
  （Geneformer 0.775 > STATE 0.765）；训练池从 91→96（+5 病人）使 UCE
  移动 +0.086——此尺度下 transfer 估计本身不可辨识，这本身就是 no-transfer
  结论的直接证据。正向"全弱"与逆向"点估计持续更高但排名不稳"的不对称
  如实并录，不作解释。
- **adaptation curve**（Jiang/GSE151263 剔除后重跑，32-patient 检验集
  [原 34]，96-patient sepsis 训练池）：**34→32 的机制（round-9 item 4
  澄清）**：Jiang 的 7 名病人在 `run_adaptation_curve_chongqing.py` 里于
  原始 embedding 载入阶段即被过滤（`data_no != JIANG_DATA_NO`，
  `load_platelet_embedding()`），发生在 `make_sepsis_holdout()` 对
  disease/severity/patient_id 做分层切分**之前**——即他们从未进入被按结局
  分层的病人池，因为他们本就没有可验证的结局记录，无法被分层或打分。
  他们的移除只是把上游总病人池缩小了 7 人；`make_sepsis_holdout()` 内
  测试集大小是该池的固定比例（`TEST_FRACTION=1/3`，
  `n_test = round(n_pool * 1/3)`），池缩小 7 人使测试集按比例缩小
  7×⅓≈2.3、四舍五入为 2——是训练/测试比例重算的结果，不是两名特定
  "无结局"病人被字面剔除。用满 sepsis 数据时全部方法 0.40–0.62（零射击
  更低）。**该检验集规模下不可分辨，不是 null 结果**：32-patient 检验集在
  该阳性率下约 3 死亡/29 存活，Hanley-McNeil SE≈0.178、95% CI 宽度≈0.70
  （[0.15, 0.85] 内均与随机不可区分）；观测的 0.40–0.62 落差 ≈1.2 SE，与
  "单一真值+噪声"一致。措辞统一为"该检验集规模下不可分辨"，不用
  "no rescue"——与 §3b 对 17-event 功效边界的处理方式一致。

### 3e. 其余次要轴

**校准**：STATE 最佳（LR：Brier 0.118/ECE 0.168；RF：STATE-RF Brier
0.095/ECE 0.079，RF 全体 ECE 0.054–0.080）。UCE 最抗嵌入维度置零。

**拼接 pooled LOBO**（206 patients/25 fatal，全部 fatal 保留）：前三并列
（STATE 0.753 / scPRINT-2 0.751 / Geneformer 0.748，差 ≤0.005），UCE
0.681、scGPT **0.638（五者中最低）**。scGPT 队列间最稳（σ=0.123）、UCE
最不稳（0.152）——但低方差围绕一个差均值不构成"稳定性"推荐理由（见 §6
部署指引）。单 seed 排名不可靠：STATE 10 次 repeat 范围 0.711–0.815，主
估计一律取 repeat-mean。

**LOBO 折规则**：held-out cohort 需 ≥1 fatal 且 ≥1 survivor，否则只进
训练池（34 名排除者全为全存活 cohort，25 fatal 全部保留在 206 人测试集
中，患病率从 10.4%→12.1%）。这一规则本身是对结局变量的选择——LOBO 度量
的是"对含死亡病例的 cohort 的泛化"，全存活 cohort 只贡献训练信号、从未
被真正留出验证；不应描述为"泛化到任意 cohort"。（这里的 34 是 LOBO 排除的
零死亡 cohort 病人数，与 §3d adaptation curve 检验集的 34→32 是两件不同的
"34"——后者是 Jiang 剔除后按比例重新分层切分的结果，7 人×⅓≈2.3，见 §3d。）

## 4. 标题与 Abstract 定稿

**标题**：*Stress-testing single-cell foundation models on underrepresented
platelets: a patient-level mortality benchmark*

（round-7 更正：原 "out-of-distribution" 在语料测量落地后不成立——三个可
测语料的模型 platelet 占比均 ~0.12%（scGPT 0.125%、UCE 0.120%、scPRINT-2
0.120–0.125%），是 Census 释出版本里数以万计、有标注的真实存在，不是
"absent"。"underrepresented"/"minority cell type" 是与测量结果一致的
表述；OOD 一词在正文中仅作为内部设计格标签保留——monocyte in-distribution
vs. platelet 对照——不再作为对 corpus 内容的科学断言使用。）

**Abstract（终稿，308 词，原 401 词——已删语料共变归因句、精简语料成分
细节至一句、UCE 重新定位、多处合并句子；round-8 item 3：310 词版意外删掉了
"recentered bootstrap P=0.073"，只留 P=0.181，重现 round 5 要求同时给出两个
数字才能避免的挑樱桃不对称——已恢复为 "P = 0.181; cross-check P = 0.073"
四词括注，314 词；round-12 item 3：Key Points 已正式成篇后，摘要内与其
重复的细节已删（corpus "~0.12%"具体数字、monocyte 结果句）——314→279 词；
round-12 item 1 新增必要限定与 deviance 结果，279→328 词；round-13 把
"更敏感统计量"改写为互补性框架、STATE-only 改为全部五个模型、把未算出的
"too small a signal"换成真实幅度，328→345 词；**round-14：修正了分辨率
标尺自相矛盾（原句同时说"检测到"又说"低于可分辨下限"）、把 in-sample
幅度换成 genuine OOF 配对结果（`round14_nested_auc_oof.py`，五模型
Wilcoxon p 全部 0.0020），同时删掉第三句冗余的控制组枚举与结尾句对第四句
的重复陈述，345→364 词**——每一次涨字数都对应一个真实的修正（模型数量、
标尺错误、in-sample→OOF），删的都是真正冗余的部分，不是注水也不是选择性
省略；**round-15：跨-seed Wilcoxon（与全文早前已排除的"seed 间 sd 当误差
项"同一类错误）换成 patient-level exact DeLong，"全部五个模型"改为
"三个稳健、一个视比较对象而定、UCE 始终不显著"，并加入判别实验发现
（增量随比较对象变干净而变大，不支持"只是修复过拟合"的替代解释），
364→361 词；**round-16：判定性结果——harness 的标准化缺陷修正后，配合
全文一贯的 BH 标准，互补性发现在三个比较对象下全部 0/5 通过，不再是
一个可报告的发现。"in discrimination"限定词随之收回，头条句恢复为无
限定的"no foundation model significantly outperformed the baseline"；
互补性段落改写为一句诚实的负结果陈述（investigated...no model's
increment survived），而不是删掉不提——这是一次真实的方法学纠错过程，
不是可以悄悄抹去的插曲**，361→328 词（净减 33 词，这次的收缩来自一个
真实结论的撤回，不是文字精简）；**round-17 item 2 指出"after correcting
a scaling defect"这半句本身违反了本文档自己的写作规则——正文不叙述
更正过程，更正过程只留在本文档（对账记录）里，摘要作为最终会进入
manuscript 正文的文本，不应该讲述内部调试情节**（一个读者读到"after
correcting a scaling defect"会追问"那其余结果是不是也在未修正版本上
跑的"——这正是 round-17 item 1 要处理的问题，但不该由摘要的措辞主动
引出）；已删除该半句，只保留干净的负结果陈述，328→308 词。

> We benchmarked five single-cell foundation models (STATE, UCE, Geneformer, scPRINT-2, scGPT) on human platelets, a minority cell type in the public atlases these models draw their pretraining corpora from. We used in-hospital and 28-day mortality as a patient-level task in 240 COVID-19 and 91 sepsis patients. All models received identical raw input and were evaluated under repeated patient-level cross-validation, every prediction aggregated to the patient before scoring, against a from-scratch PCA + logistic-regression baseline, a null model, and an in-distribution monocyte positive control. Mortality was learnable from platelet transcriptomes (patient-level AUC up to 0.79), but under a pre-specified, selection-corrected test no foundation model significantly outperformed the baseline (label-permutation P = 0.181; cross-check P = 0.073; stable across all 25 leave-one-fatal-out refits). A separate question — whether a model's score adds information beyond the baseline's even without ranking patients better on its own — was examined with a nested-model comparison; no model's increment survived the same multiple-comparison correction used throughout (Benjamini-Hochberg, 0/5 in every comparator configuration tested). Within sepsis, no model exceeded chance at late sampling windows (17 events, resolving only AUC ≳ 0.71), consistent with published platelet-outcome signals being reported only in early disease windows. Cross-disease transfer was rank-indeterminate in the direction where signal exists and uninformative in the direction where it does not; the adaptation curve has the same resolution limit. The one result that survived selection-corrected and multiple-comparison testing was that two models (STATE, Geneformer, statistically indistinguishable from each other) separated from a third (UCE), the weakest method throughout — the only one that failed to clear the null model built from cell counts and disease identity. At realistic clinical event counts, pretraining bought no detectable advantage over a strong from-scratch baseline on this underrepresented cell type; the benchmark quantifies the resolution achieved (~0.10 AUC at 25 events) and the scale needed to detect smaller effects.

**摘要数字核对**：语料共变归因句已从摘要**删除**（round-7 item 1：三个可测
语料模型 platelet 占比几乎相同（~0.12%），其 within-COVID AUC 却横跨
0.079——全部模型间总差距的 63%——语料占比在可测的三个模型上根本不排序；
而双雄（真正撑起这个差距的两个模型）恰好是语料**不可测**的两个，归因句在
两端都没有测量支撑）。该假设保留在 Discussion，措辞为"未经检验的候选属性
之一"，不作为摘要结论。UCE 重新定位（round-7 item 5）：双雄-UCE 分离的
内容是"两个模型区分于一个连 null 都过不了的模型"，即关于 UCE **不适配**
此任务的陈述，不是关于双雄优势的陈述——语料假设因而是在解释一个失败模式，
不是解释一个优势。QC 增量价值检验在全部六个代表（5 FM+PCA-200）上显著
（round-10 修正后 nested-LR 口径：五个 P<0.0001，UCE P=0.0022；DeLong
AUC 赛跑口径：五个 p<0.01，UCE p=0.072 不显著，见 §3a），故 "learnable
from platelet transcriptomes" 不降级；但方法-vs-null 配对检验显示
**UCE 是唯一未显著超过 count+disease null 的方法**（p=0.157，这是与
增量价值检验完全独立的另一个检验，不受 round-10 影响），Discussion 需
明确陈述。

## 5. 一句话故事线（P2 重写围绕它）

在所有设计格里——包括唯一曾看似成立的 within-COVID——预训练相对
from-scratch 基线**在判别力上**的增值都不可分辨：主检验 permutation
P=0.181（recentered bootstrap P=0.073），25/25 LOFO refits 全 P>0.05。
一个额外的互补性问题（控制 PCA-200+队列后 FM 是否仍额外贡献信息，与
"FM 排序是否赢 PCA-200"这一优越性问题分开）也检验过：配合全文一贯的
BH 标准，三个比较对象下全部 0/5 通过，不构成第二个发现（完整数字、
构造与诊断见 §3a）。没有需要额外说明的第二种"存在但小"的信号。能分辨的
模型间差距只有一个——双雄 vs UCE 超 +0.11（confirmatory：range-statistic
permutation + BH q=0.030）。**这个结果关于 UCE，不关于双雄**（round-7 item 5 更正）：
UCE 是全文六个方法里唯一未能在 count+disease null 检验中显著超过该 null
的方法（p=0.157，且 UCE 是 6 个 null 对比里唯一 BH 校正后仍不过的一个；
round-11 item 1 澄清：原"在任何检验中"措辞对 nested-LR 检验不再准确，
已收窄为明确指代这一个检验）；
双雄-UCE 分离的真实内容是"两个模型区分于一个连 null 都过不了的模型"，
是一个关于 UCE **不适配此任务**的发现，不是关于双雄有预训练优势的发现。

**语料占比归因假设已被语料测量本身推翻，降级为未经检验的候选属性
（round-7 item 1）**：三个可测语料的模型（scPRINT-2/scGPT/UCE）platelet
占比几乎相同（~0.12–0.125%），但它们的 within-COVID AUC 横跨 0.079——占
全部五模型间总差距（0.125）的 63%——语料占比在唯一能测的三个模型上根本
不排序模型表现。而双雄（真正撑起 UCE 差距的两个模型）恰好是语料**完全
不可测**的两个（Geneformer/STATE 均未公开语料细胞类型构成）。归因句因而
在两端都没有测量支撑：能测的地方不排序，会排序的地方测不了。血液/免疫
（而非 platelet 特异）语料占比这一更贴切的量仍未以可比单位测过全部五个
模型（Geneformer 只有数据集数、UCE 只有样本比例、scGPT 只有细胞数——
三种不同单位）。语料假设保留在 Discussion，措辞为 12 个候选属性之一、
未经检验、正在解释 UCE 的失败模式而非双雄的优势，不作为正文结论。

**round-18 更正（2026-09-11，main.tex 一致性 review 发现）**：上文"双雄
（真正撑起 UCE 差距的两个模型）恰好是语料完全不可测的两个（Geneformer/
STATE 均未公开语料细胞类型构成）"不准确。model_architecture_analysis.md
（2026-09-06 核实）显示 STATE 的语料 platelet 状态实为**已界定**：其
scBaseCount 部分（约占语料 2/3）在注释层与 QC 层双重确定 ABSENT，CxG 部分
与 scPRINT-2/scGPT/UCE 共享同一 Census 上限（≤0.1%），合计将 STATE 的
platelet 占比上界压至 <0.05%——这与 UCE 同属"可测"范畴，不是"不可测"。
真正语料构成完全无法界定的只有 Geneformer 一个模型（文档层已穷尽，全部
出版材料均无细胞类型级构成）。main.tex 已按此改写（Introduction、Methods
§ssec:methods_models、Discussion §ssec:discussion_property 三处同步）：
"三个可数值界定的模型...STATE 的语料界定更低（低于 0.05%）...仅
Geneformer 未公开任何细胞类型级语料构成"。上文原表述保留作为历史记录，
本条为其更正。

两条件解释框架 [疾病贴近 × 语料覆盖] 仍不作为正文解释：条件一只能由
基线自身表现度量，循环不可证伪。

**本文的贡献因而是标尺而非冠军榜**：在临床现实事件数下给"预训练能买多少"
量化出上界（25 events 分辨 ~0.10 AUC；0.03–0.05 需 ~10³ patients），并
提供任何新模型可即插即测的协议——对下表这样一个已系统核实、多数未做
selection/multiplicity 校正的领域，这是更可用的东西。

**综述内容表（round-9 item 3 新增，round-10 item 2 补充搜索策略与 EGSP
范围标注，实地核实——补足"benchmarking review"定位所需的 review 成分，
把 §5 原先"未复现阳性结果的领域"这一无引用断言换成可核对的三项）**：

**搜索策略（round-10 明确要求，如实陈述——这是"我们能找到的三项"，
不是穷尽性系统综述）**：PubMed（经 Google 检索 `site:pubmed.ncbi.nlm.nih.gov`）
+ bioRxiv/medRxiv，查询串 `"single-cell foundation model" patient-level
clinical outcome benchmark survival mortality 2025 2026`，检索时间
2026-09-09，纳入标准=对≥1个单细胞基础模型（预训练 transformer，直接跑
单细胞输入或将其单细胞预训练 embedding 迁移到 patient-level/bulk 结局
预测）在病人级临床结局任务上的系统评测。该查询召回 4 条不重复记录，
2 条纳入（EGSP、9-模型 benchmark）、2 条排除——scSurvival（Cancer
Discovery 2026，PMC）：attention-MIL Cox 回归，特征提取器是从头训练的
VAE，不基于任何预训练单细胞基础模型，排除；Chen 等《Leveraging
Hematologic Single-Cell Measurements...》（J Appl Lab Med 2025）：数据是
血液分析仪 cell population data（CPD），非 scRNA-seq，亦非基础模型，
排除。第三项 TissueFormer 并非来自此次通用查询——它在本项目更早的
文献核查（§3b 引用的 Task 6a）中已被识别，此次用其名称做定向复核确认
仍可检索到（bioRxiv 2025.08.17.670735），故其纳入路径与前两项不同，
如实标注。

| 研究 | scFM(s) | 单细胞输入方式 | 任务/结局 | 队列规模 | 事件数 | from-scratch 基线 | 病人级聚合 | 多重/selection 校正 |
|---|---|---|---|---|---|---|---|---|
| EGSP（*Leveraging single-cell foundation models for accurate survival outcome prediction*, **Bioinformatics Advances 2026（已同行评审）**，PMC13032892） | scFoundation | **embedding 迁移**（scFoundation 单细胞预训练 embedding 迁移到 bulk RNA-seq，非直接跑单细胞输入——按上方纳入标准的第二支路径纳入，非直接单细胞评测） | 泛癌 overall survival（TCGA，25 个癌种） | 8,776 样本（按癌种 57–1,079 不等） | 3,095 death / 5,681 censored（合并计） | 有（mRNA+Clin、Ridge-Cox、Lasso-Cox） | 有（bulk 衍生的病人级 embedding） | 无（Wilcoxon signed-rank，未见多重校正） |
| *Empirical Evaluation of Single-Cell Foundation Models for Predicting Cancer Outcomes*（**bioRxiv 预印本**，2025.10.31.685892, PMC12637420） | 9 个模型（Geneformer V1/V2、scGPT ×2、scFoundation、SCimilarity、CellPLM） | **直接**（单细胞输入） | 6 项癌症临床/生物学结局任务 | 未按任务披露具体人数（原文称"typically small"），合计 1,170 个 supervised 实验 | 未按任务披露 | 有（HVG、PCA、scVI 三种） | 有（majority-vote / pseudo-bulk / MIL 三种聚合） | 无（9 模型 × 6 任务未做多重校正） |
| TissueFormer（**bioRxiv 预印本**，2025.08.17.670735, PMC13466306） | Murine-adapted Geneformer | **直接**（单细胞输入） | COVID-19 严重度三分类（对照/轻/重） | 376 donors（COMBAT 85、Ren 185、Stephenson 106，三队列合并） | 各队列 severe 组 46/83/28（最小分层 28） | 有（random-init TissueFormer，"与预训练模型表现不可区分"） | 有（按 donor 分折） | 无 |

三项均未报告 selection/multiplicity 校正，且除 EGSP 的合并事件数（3,095）外，
其余两项均未按具体任务/分层披露事件数，意味着读者无法核对其比较是否有
足够事件数来分辨所报告的差异——这正是本文引入"分辨率"框架（在报告
是否显著之前先给出该样本量能分辨的效应量下限，如本文 17 events→AUC≳0.71、
25 events→~0.10 AUC）而非只报"是否显著"要解决的问题。该表同时提供本文
定位为"五个已发表模型的 benchmarking review"（§7 定案）所需的 review
成分，并把本段的资源比较从修辞断言变成可核对的数据。

**round-11 item 3：三行表看起来太窄，需要点名它自己排除了什么（已核实
引用真实存在）。** 上表按纳入标准（病人级临床结局）只收 3 项，但细胞级
scFM benchmarking 这一更成熟的姊妹领域确实存在且被广泛引用，读者会先
注意到它们缺席，再注意到本文的纳入标准。代表性的一项：**Kedzierska,
Crawford, Amini 等，*Zero-shot evaluation reveals limitations of
single-cell foundation models*，Genome Biology 26:101（2025），
doi:10.1186/s13059-025-03574-x（WebSearch 核实，早期版本为 bioRxiv
2023.10.16.561085）**——在聚类与 batch integration 等细胞级任务上，
零样本 Geneformer/scGPT 被 scVI、HVG 选择、均值表达等更简单的方法超越；
扰动预测与细胞动态重建方向也有同形状结论的相关 benchmark。**点名它服务
两个不同的论点，需要在表述上分开**：(a) **gap 论点**（本表 3 项支撑）——
病人级临床结局的 scFM 评测本质上不存在有效事件数披露、也几乎不做
selection/multiplicity 校正；(b) **姊妹领域已确立同型结论**（Kedzierska
等支撑，不需要更多行）——细胞级 scFM benchmarking 已经成熟，且已发现
"预训练不敌精心选择的 from-scratch 基线"这一同形态结论，本文是这一结论
在病人级的类比，不是孤立的 null 结果。Kedzierska 等的建议之一
（评测集应与全部未来预训练数据隔离，因为无法排除污染）也直接命中本文
corpus 段的局限——5 个模型中 2 个完全未公开语料构成，platelet 占比与
本文队列重叠与否均无法核查——引用它把这段从"本文的局限"变成"一个已被
文献确认的结构性问题的一个实例"。**两项为 bioRxiv
预印本（检索时 2026-09-09 尚未见同行评审版本），EGSP 已同行评审——如
实标注于表内，供读者自行判断证据权重。**

monocyte"不行"不是模型的失败：倒挂全体成立（含 PCA），是任务属性。它
照亮的是另一面——即使 in-distribution，也没有 FM 显著赢基线，所以
platelet 的 null 不是 OOD 的特例，而是这套 pipeline 在各处的一致读数。

QC 一役的最终形态：QC 3 特征单独 as-fit 反向（0.376）是条件于 library
depth 的经典抑制模式（suppression），条件后（QC+count=0.618）方向转正，
不是折结构不稳，也不是异常；增量价值检验（QC+计数+队列之上加表示）是
唯一能定案的实验——这也回答了 Adam 第一封邮件的终极版本（"几乎可以从
细胞数预测结局"）：**预声明主检验**（PCA-200 nested LR，round-10 修正
参考分布后 P<0.0001，无需模型选择）显示控制 QC+计数+队列后，platelet
表示仍提供显著增量信息，且这一结论现在比 round-8/9 定案时更强而非更弱
——round-10 修正后 5 个 FM 的 exploratory 读数在 nested-LR 口径下也全部
显著（UCE 最弱，P=0.0022），仅在 DeLong AUC 赛跑口径下 UCE 仍不显著
（见 §3a）。

**round-17 item 1 提出的合并陈述**（这个检验与 round-16 的互补性检验
问的是两个不同比较对象，round-10 问"表示是否比裸 QC+计数+队列多提供
信息"[比较对象不含 PCA-200]，round-16 问"FM 是否比 PCA-200+队列多提供
信息"[比较对象含 PCA-200]，二者可以同时成立，合在一起是本文更完整的
论点）——**platelet 表示携带超出 QC、细胞数与队列的结局信息；但没有
任何一个基础模型表示携带超出一个 from-scratch PCA（用同一批数据算出）
的结局信息**。这是本文一句话论点，两条证据链都支持它，不需要读者自己
拼合两处相隔数页的 BH 计数。

## 6. P2 编辑原则（BIB 密度控制）

正文一个故事 + 四张主表（within-COVID / within-sepsis / cross-disease /
in-dist-vs-OOD 对照）；审计、v2 映射、cohort 细节、模型三维档案、冲突协议
全部进 supplement；Results 压缩至 6–7 节、Discussion 至 4 节。

**部署指引**：全文结论 = 在任何设计格（含 in-distribution 的 monocyte
对照）都未分辨出任一 FM 显著优于 from-scratch 基线，部署指引必须与此一致。
**对这一细胞类型、这一事件数量级的病人级结局预测，可辩护的建议是使用
from-scratch PCA+LR 基线**——与全部测试过的 FM 统计不可区分（in-distribution
与 OOD 皆然），且运行成本低几个数量级。校准观察（STATE-RF Brier 0.095/
ECE 0.079，RF 全体 ECE 0.054–0.080）保留为"若仍选择使用某 FM，选哪个"的
附加说明，不作为使用 FM 本身的推荐理由。**不推荐 scGPT 的"稳定性"**：其
拼接 LOBO 是五者中最低（0.638），低方差围绕一个差均值不构成推荐理由。

## 7. BIB 投稿风险核对（round-8 item 4，实地核实，仅供人工决策——本 agent
未联系期刊）

来源：`academic.oup.com/bib/pages/msprep_submission`（Manuscript Preparation，
2026-09-08 网页核实）。

**(a) 摘要/关键点/关键词要求**：官方指引未给出明确的摘要字数上限，只要求
"a short abstract"；**同时要求最多 5 条 Key Points**和最多 6 个
keywords，外加约 30 词的作者简介。当前 314 词摘要相对"short"的措辞可能
偏长。**round-11 item 4 要求把这 5 条从散文提升为真正的板块（BIB 明确
要求的必备元素，目前仍缺）——正式稿如下**，投稿前直接誊入 cover
letter/摘要页：

**Key Points（≤5 条，BIB 必备板块，round-9 排序建议已采纳：③提到②之前，
因为它最反直觉、最容易被读者记住；round-10/11 后措辞已同步更新）**

1. Five published single-cell foundation models were benchmarked on platelet mortality prediction; none significantly outperformed a from-scratch PCA+logistic-regression baseline under a pre-specified, selection-corrected test (permutation P=0.181).
2. A pretraining-corpus-composition hypothesis for the one confirmatory finding was overturned by the corpus measurement itself — the three models whose corpora can be measured show no ordering by platelet content.
3. The one result that survives correction is not that two models have a pretraining advantage — it is that a third (UCE) is the weakest performer throughout, the only one that fails to clear the count+disease null model.
4. The same null held on an in-distribution positive control (monocytes), showing platelet underrepresentation is not the explanation — the benchmark's pipeline finds no foundation-model advantage anywhere it looked.
5. The paper's contribution is a resolution-quantified benchmarking protocol (~0.10 AUC resolvable at 25 events), not a leaderboard — a distinction that matters because comparable patient-level scFM evaluations in the literature report no correction for the comparisons they make, and none reports event counts at the granularity of those comparisons.

**(b) 范围风险（更严重）**：BIB 官方 author guidance **明确将自己定位为
review 期刊，声明不发表 "pure original research"**，并建议此类稿件另投
他处。实践中 BIB 确实经常发表 benchmarking 类研究，本文大概率能存活，但
这是一个在投入大量修改后才发现的 desk-rejection 风险，值得投稿前低成本
解决。

**决定（用户，2026-09-09）：走 (ii)**——本文正式定位为"对五个已发表单细胞
基础模型的 benchmarking review"，非原创研究。这也确实更贴近本文的实际
性质（五个已发表模型 + 标准化评测协议，不提出新模型/新方法）。落实：
cover letter 需明确使用"benchmarking review"措辞；Introduction 开篇一句
话应同样定位为对已发表模型的系统性评测而非新方法论文；标题已含
"benchmark"一词，不需改动。此决定不改变任何数字或统计结论，仅为投稿
框定，P2 重写时在 cover letter + Introduction 落实即可，不需要在
本文档中进一步展开。
