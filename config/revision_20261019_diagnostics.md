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
