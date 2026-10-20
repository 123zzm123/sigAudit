# Uno/IPCW 敏感性分析的实现配置（事后记录）

**状态说明。** 本文件记录脚本 `52_uno_sensitivity.R` 与 `55_uno_validation.R` 实际执行的设置，
写于分析完成之后，用于让读者核对“做了什么”，**不是预注册、不是分析前冻结的方案**。
`config/frozen_analysis_plan.md` 与 `config/metabric_frozen_config.md` 仍然只覆盖主分析。

## 固定设置

| 项 | 取值 |
|---|---|
| 截断时间 tau | 60 个月（固定，不随队列、任务、结果调整） |
| 删失分布 G | 训练折内、对删失指示 (1 − event) 做 Kaplan–Meier；右连续取值；t 小于首个时间点时取 1 |
| 权重 | 1 / G(T_i)^2，只对事件方 i 加权 |
| 可比对 | T_i < T_j（严格）、delta_i = 1、T_i < tau |
| 事件时间并列 | 两个事件同一时刻：不构成可比对（`survival::concordance` 记为 tied outcome，同样不进入分母） |
| 事件–删失同刻 | 不构成可比对（要求 T_i < T_j）；`survival::concordance` 把删失者视为更晚，会保留该对——这是两者分母差异的主要来源 |
| 分数并列 | 计 0.5（两种实现一致） |
| 事件时间并列按 0.5 | **不采用**；仅在 S43 中作为对照约定列出 |
| 事件方支持 | G(T_i) < 0.05 的事件不进入分母 |
| 队列×划分支持 | G_train(tau) >= 0.05 且训练随访上限 >= tau；不满足则整个划分不进入分析 |
| 最少事件数 | 每个任务 tau 前测试事件 >= 10 |
| 参评估计量 | 固定 Lasso、固定 Ridge、固定 Enet0.5、argmax_val、top5_ensemble |
| 选择规则 | 不变：验证集 Harrell C，并列按 `_lib_analysis.R::make_prio()` |
| 汇总层级 | 任务 → 队列×划分均值 → 队列内各划分等权 → 六队列等权；METABRIC 为 任务 → 签名×划分均值 → 签名内五划分等权 → 六签名等权 |

## 实现核对

1. **同定义核对（主证据）**：`scripts/55_uno_validation.py` 用独立的 Python 实现重写上述定义，
   在 7 组合成数据集（分数全等、完全正序、完全反序、分数并列、事件时间并列、tau 边界、
   随机半数删失）与 1 个真实任务上逐值比对 `_lib_uno.R`，共 15 个 C 值比较，
   最大绝对差 1.5e-14。结果见 `results/tables/55_uno_validation.tsv`。
2. **删失分布来源差异（参考）**：训练折 G 与 `survival::concordance(timewt = "n/G2", ymax = 60)`
   在 180 个任务上的差异，中位 0.003、最大 0.058（STAD 的某一划分）。两者不是同一定义：
   后者用分析数据估计 G，且 `ymax` 在 tau 处截尾而不是按 T_i < tau 排除事件。

## 未做的事

未用 Uno 重新选择模型；未在 Uno 与 Harrell 之间做差的差；未搜索 tau；未做时间依赖 AUC、
Brier 或校准；未把该分析扩展到 LGG、PAAD 两个进展终点队列（事件构成未核实）。
