# 本次修订新增分析：口径与结果快照（2026-10-19）

本文件记录本次修订采用的分析口径和结果快照，包含输入、输出与日期。分析是在既有结果与
审阅意见基础上形成的，属于**事后诊断**；它们不是原研究的事前预设（pre-registered plan），
也不得倒写为原研究预设。某些实现细节（如去重判据与代表选择规则）在运行前已确定，但不应
据此把全部分析说成"未知结果时预先固定"。

## 0. 背景与边界

- 论文既有主分析（3,919 个共同可评估任务）**保持不变**，仍是主要基准。
- 六 OS 队列敏感性、METABRIC 独立流程复制均不变。
- 新增内容只有三块：① 非恒定 LASSO 分层诊断；② 验证预测去重敏感性；
  ③ 预测级审计软件入口。不新增队列、不做 nested CV、不重训模型（除非另有说明）。

## 1. 非恒定 LASSO 分层诊断

- **脚本**：`scripts/47_nonconstant_diagnostic.R`
- **输入**：`results/tables/20_rule_performance_long.tsv`、`20_model_metrics_long.tsv`；
  METABRIC 段另读 `results/metabric/31_metabric_predictions.rds` 与 31/32 表。
- **基准集合**：`argmax_val` 的 3,919 个共同可评估任务（与主分析一致）。
- **分层变量**：LASSO 的**验证预测是否恒定**（`20_model_metrics_long.tsv` 的 `constant`
  列，定义 = 验证预测分数只有一个唯一值）。不按测试表现筛选。
- **参照**：固定 LASSO（规则 `none`）；稠密对照 = 固定 Ridge。
- **汇总层级**：任务 → (队列×划分) 平均 → 队列内跨划分平均 → 八队列等权。
- **输出**：`47_nonconstant_overall.tsv`、`47_nonconstant_by_cohort.tsv`、
  `47_nonconstant_by_cell.tsv`、`47_nonconstant_os6.tsv`、`47_nonconstant_metabric.tsv`、
  `logs/47_nonconstant_diagnostic.log`。
- **结果快照（全精度复核）**：
  - 全体 3,919：Gv=0.06101166985575、Gt=0.02365691698918、Orel=0.03735475286657；
    argmax 测试 C=0.60789200715545、Ridge=0.61427437335451、argmax−Ridge=−0.00638236619906。
  - 非恒定 3,057：Gv=0.04943319959182、Gt=0.01092472551194、Orel=0.03850847407988；
    argmax 测试 C=0.61037058730147、Ridge=0.61723413038314、argmax−Ridge=−0.00686354308167；
    覆盖全部 8 队列、40 个 cohort×split 单元，八队列方向一致（均负）。
  - 常数任务 862 个：COADREAD 239、BRCA 222、STAD 197、LUAD 96、MESO 72、PAAD 26、
    KIRC 9、LGG 1；仅覆盖 33 个 cell，不作队列级宏平均强推断。
- **口径限制**：事后条件分层，不是因果证明，也不替换主分析集。

## 2. 验证预测去重敏感性

- **脚本**：`scripts/48_dedup_sensitivity.R`
- **输入**：`results/main/19_*_100.rds`（八队列×五划分）、`20_model_metrics_long.tsv`、
  `20_rule_performance_long.tsv`、`20_rule_members_long.tsv`。
- **判据**：以**验证预测向量精确相等**为同一组（不使用容差、不排序去重、不剔除常数）。
- **选代表**：同组内按候选名称**固定字典序**取最小者为代表。
- **参照**：固定 LASSO / Ridge 输出始终保留，即使与池中其他候选重复。
- **新池规则**：validation argmax、top-five 集成（唯一候选不足 5 个时用 top-min(5, K)）。
- **测试信息**：不得用于分组或选代表；本实现只用验证预测。
- **输出**：`48_dedup_summary.tsv`、`48_dedup_by_cohort.tsv`、`48_dedup_coverage.tsv`、
  `48_dedup_per_task.tsv`、`48_dedup_performance_long.tsv`、`logs/48_dedup_sensitivity.log`。
- **结果快照（三层级）**：
  - argmax 获胜者在 68/3,919（1.7%）任务中改变；Orel 原池/去重池均 +0.0374（四位小数）。
  - top-five 成员集合在 699/3,919（17.8%）任务中改变；Orel 由 +0.0225 变为 +0.0211；
    top-five−argmax 由 +0.0048 变为 +0.0052。
  - argmax−Ridge 保持 −0.0064。
  - 分析池平均候选 16.9643（唯一 14.1056）；原始池 17.4437（唯一 14.5606）；
    1,879 个任务原始池恰好多一个候选；含重复对的任务占 72.9%（原始池 73.1%）；
    唯一候选不足 5 个的任务 6/3,919（0.15%，均在 MESO）。
  - `48_dedup_coverage.tsv` 的 `n_tasks` 为队列任务数（BRCA 495、COAD 495、KIRC 472、
    LGG 477、LUAD 495、MESO 495、PAAD 495、STAD 495），不是平均候选数。

## 3. 预测级审计软件入口

- **实现**：`R/audit_predictions.R`（导出 `audit_predictions()`；不重训模型）。
- **一致性核对脚本**：`scripts/49_audit_prediction_api.R`
- **输入**：与 20 相同的 40 个冻结预测；逐患者风险分数、时间、事件与候选标识。
- **输出**：`49_audit_api_check.tsv`（一致性）、`49_audit_api_reference_switch.tsv`
  （更换参照）、`49_audit_api_dedup.tsv`（去重演示）、`logs/49_audit_api.log`。
- **判定**：与原流水线逐任务比较，Cv/Ct/Gv/Gt/Orel 的最大绝对差须 < 1e-10；
  argmax_val/top5_ensemble 的获胜者与成员须完全一致；记录数须等于预期、键须唯一、
  API 的 (task,rule) 须全部能在主分析表中找到，否则脚本以非零状态停止。
- **本轮修复（大修后评审）**：`top_k` 统一注册并实际计算；`linear_only/tree_only`
  可请求且未知规则名报错；默认规则不含 `oracle_test`；`by_cohort` 复用
  "队列×划分 → 队列"的中间层；校验事件长度与 0/1 编码、候选名唯一、ID 唯一且
  验证/测试不重叠；缺失按任务级共同有限掩码处理并记录丢弃数；汇总只在共同
  可评估集合上计算，保证汇总层 `Orel = Gv − Gt`。
- **复用来源映射**：见 `R/audit_predictions.R` 文件头（ci_direct、.fixed_prio、
  ens_pred、rules_one、三层级汇总口径）。

## 4. 汇总与推断口径纪律（不随结果调整）

- 主基准仍为八队列全任务；六 OS 为终点一致的敏感性；METABRIC 为独立流程复制。
- 不因哪组 p 值更小或排序更漂亮切换主集。
- 新增分析一律标注为本次修订的探索性/事后诊断。
- `Orel/Gv` 比例在分母很小或非正时不赋予稳定意义；显示值与全精度的 0.0001 差异
  由表注说明，不为对齐显示值而篡改 Orel。

## 5. 基线条件下的搜索规模曲线（第 8 包评审新增）

- **脚本**：`scripts/50_search_budget_audit.R`；绘图 `scripts/51_search_budget_figure.R`。
- **输入**：`20_model_metrics_long.tsv`、`20_rule_performance_long.tsv`、
  `31_metabric_model_metrics.tsv`、`31_metabric_rule_performance.tsv`。
- **设计**：参考锚定搜索。每个任务固定起点（`Lasso` 或 `Ridge`，起点始终在池中），
  其余 K−1 个候选从剩余 N−1 个中**均匀无放回**抽取；只按**验证** C-index 选优，
  并列用原完整池的 `make_prio(names, seed)`（不在子池内重新生成随机数）。
  K ∈ {1,2,4,8,12,full}；full = 该任务自己的完整分析池，**不统一冒称 K=18**。
- **期望值精确计算**：顺序按验证 C 降序、以固定优先级破并得到 i=1..N，起点位置 r。
  K≥2 时 p_r = C(N−r, K−1)/C(N−1, K−1)；i<r 时 p_i = C(N−i−1, K−2)/C(N−1, K−1)；
  i>r 时为 0。K=1 取起点。加权求和即期望 Cv/Ct，相减得 Gv/Gt/Orel。
  另给出"子池内均匀随机选一个模型"的期望 C 作为对照（非风险分数集成）。
- **汇总口径**：TCGA 任务 → (队列×划分) → 队列 → 八队列等权；六 OS 同法；
  METABRIC 任务 → 签名内跨 seed 平均 → 六签名等权。
- **输出**：`50_search_budget_task.tsv`、`_by_cohort.tsv`、`_summary.tsv`、`_os6.tsv`、
  `_metabric.tsv`、`_check.tsv`；图 `fig5_search_budget.*`、`figs2_*`、`figs3_*`；
  图数据 `figdata_fig5.tsv`、`figdata_figs2.tsv`、`figdata_figs3.tsv`。
- **核对（`50_search_budget_check.tsv`，全部通过）**：
  - 概率和为一（最大偏差 5.99e-15）；K=1 取起点、增益恰为 0；
  - K=full 与旧 `argmax_val` 逐任务一致，TCGA 与 METABRIC 的最大绝对差均为 0；
  - 公式对 N=2..8 的全部起点位置与 K 与穷举子池一致（最大偏差 3.33e-16）；
  - 期望验证表现随 K 非下降：7,838 个 任务×起点 序列中 0 例违反；
  - 47,028 行任务级记录、3,919 个任务、7,838 个 任务×起点。
- **结果快照（八队列等权 / 六签名等权）**：
  - TCGA、LASSO 起点：期望测试 C 0.5842→0.6075(K=12)→0.6079(full)；
    Gv=+0.0610、Gt=+0.0237、Orel=+0.0374（与主稿完全一致）。
  - TCGA、Ridge 起点：Gv=+0.0314、Gt=−0.0064、Orel=+0.0377（与主稿一致）。
  - 六 OS：LASSO 起点 Gt=+0.0292；Ridge 起点 Gt=−0.0057。
  - METABRIC：LASSO 起点 Gv=+0.0151/Gt=+0.0067/Orel=+0.0084；Ridge 起点
    Gv=+0.0079/Gt=+0.0002/Orel=+0.0077（与主稿一致）。
  - 均匀随机对照：TCGA full 池期望测试 C 0.5952，低于验证选优的 0.6079，
    更低于固定 Ridge 起点的 0.6143。
- **属性**：已见结果后的修订性探索，不是事前预设；不按测试表现选"最优 K"；
  每池至多约 18 个配置，不外推到大型算法库。

## 6. 真实数据软件示例（第 8 包评审新增）

- **导出脚本**：`inst/examples/export_audit_example.R`（从本地冻结预测导出，去标识化）。
- **示例数据**：`inst/extdata/audit_example_{val,test,survival}.tsv`。
  选择规则固定且与效果无关：签名按字典序最小（EndoPredict_11）、seed 最小（20261014）；
  18 个候选、验证 494、测试 495 位患者；ID 改写为 V####/T####；只含风险分数与生存结局。
- **运行示例**：`inst/examples/audit_predictions_metabric.R`，只调用**安装后的包**，
  不 source 研究脚本、不重训模型。
- **实测**：本机（R 4.6.1，Windows）一次审计 wall time 约 0.5 秒（含读取文件）。
  示例同时演示更换固定参照后 Cv/Ct 不变、增益与 Orel 改变。

## 7. 第 8 包评审的收尾修复

- `audit_predictions()`：无可评估记录时给出明确错误（不再对空表做 `aggregate`）。
- 汇总口径文字由"共同可评估集合"改为准确的 **within each rule**：每条规则只在
  其自身可评估任务上汇总，不同规则集合可不同，因此宏平均不是规则间配对比较；
  文档建议比较两条规则时取共同任务键。
- 家族规则：`linear_only`/`tree_only` 改为**非空即返回**（单候选时返回该候选），
  与原 `rules_one()` 一致。
- `49`：预期覆盖键改为**直接从冻结 `20_rule_performance_long.tsv` 构造**，
  并做**双向差集**（预期缺 API / API 多出预期），单个 RDS 整体缺失不再掩盖缺口。
- `46`：仓库白名单补 `48_dedup_coverage.tsv` 与全部 `50_*` 表。
- `43`：新增图文件同步——按手稿图号把 5 张主图（含新 Fig. 5）从 `results/figures`
  复制到 `submission_cmpb/figures`，缺文件或数量与图注不一致时**直接失败**（此前该目录
  靠手工维护，Fig. 5 一度缺失）。
- `44`：新增逐图文件存在性断言（Fig. 1–5 各有 PDF 与 PNG）与"恰好 5 组图"断言；
  通过项由 34 增至 40。
- `46`：导出统计不再把 `.git` 计入文件数与体积（此前因提交对象而虚高）。

