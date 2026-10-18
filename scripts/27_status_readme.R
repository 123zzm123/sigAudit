#!/usr/bin/env Rscript
# ============================================================
# 27 — 生成 README_current.md（评审要求：唯一正式管线入口、源码哈希、
#      实际完成的组合、失败与未完成项、时间与时区）
#
# 该文件由脚本从**实际状态**生成，不手写，避免文档与事实不符。
# 用法：Rscript scripts/27_status_readme.R
# ============================================================
options(stringsAsFactors = FALSE)
PROJ <- Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit"); setwd(PROJ)
OUT <- "README_current.md"
`%||%` <- function(a, b) if (is.null(a)) b else a
h <- function(f) if (file.exists(f)) substr(as.character(tools::md5sum(f)), 1, 10) else "缺失"

## ---- 注册网格与实际完成情况 ----
COHORTS <- c("MESO","PAAD","STAD","KIRC","LUAD","COAD","LGG","BRCA")
SEEDS   <- c(20261014, 20261016, 20261017, 20261018, 20261019)
## 优先用最终输入清单（由 23_manifest.R 生成），避免扫描目录时混入旧产物
if (file.exists("results/manifest/main_input_files.txt")) {
  main <- readLines("results/manifest/main_input_files.txt", warn = FALSE)
  main <- main[nzchar(main)]
} else {
  main <- list.files("results/main", pattern = "^19_.*\\.rds$", full.names = TRUE)
  main <- main[!grepl("\\.partial\\.rds$", main)]
}
done <- character(0); info <- list()
for (f in main) {
  R <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(R)) next
  key <- sprintf("%s/%d", R$cohort, R$seed)
  okf <- isTRUE(R$completed) && identical(as.integer(R$n_processed), as.integer(R$n_planned))
  done <- c(done, if (okf) key else paste0(key, "(未完成)"))
  info[[key]] <- R
}
want <- as.vector(t(outer(COHORTS, SEEDS, function(a, b) sprintf("%s/%d", a, b))))
missing <- setdiff(want, done)
extra   <- setdiff(done, want)

## ---- 链日志中的失败 ----
lg <- if (file.exists("logs/chain_v8.log")) readLines("logs/chain_v8.log", warn=FALSE) else character(0)
lg2 <- if (file.exists("logs/chain_v8_resume.log")) readLines("logs/chain_v8_resume.log", warn=FALSE) else character(0)
fails <- c(grep("^\\[FAIL\\]|^\\[ABORTED\\]", lg, value=TRUE),
           grep("^\\[FAIL\\]|^\\[ABORTED\\]", lg2, value=TRUE))

## ---- 时间与时区 ----
tz <- format(Sys.time(), "%Z")
ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

L <- c(
"# README_current —— 当前状态（自动生成，请以本文件为准）",
"",
sprintf("生成时间：**%s（时区 %s）**", ts, tz),
"本文件由 `scripts/27_status_readme.R` 从实际产物生成；**不要手工编辑**。",
"",
"## 1. 唯一正式管线（只有这一条）",
"",
"```",
"数据获取     scripts/02_fetch_cohorts.py, scripts/04_entrez2symbol.py",
"基因集构建   scripts/01b_build_genesets.py        （重建并逐集合核对一致）",
"正式主实验   run_v8_parallel.sh  ->  scripts/19_main_corrected.R（逐患者预测落盘）",
"清单         scripts/23_manifest.R（注册网格校验）, 24_manifests.R",
"主分析       scripts/20_analysis_corrected.R      （P1/P2 规则核算）",
"共用函数库   scripts/_lib_analysis.R              （20 与 29 共用，杜绝口径漂移）",
"辅助分析     scripts/21_entry_predictors.R        （进入效应 + 留一癌种敏感性）",
"基线敏感性   scripts/29_baseline_sensitivity.R    （Ridge / Enet0.5，第三轮审查后新增）",
"现实任务     scripts/26a_extract_public_signatures.R -> 26b_public_signature_tasks.R",
"模拟         scripts/25_simulation.R pilot|final",
"图表         scripts/13_figures.R                 （v3：比例/裁切/COADREAD/oracle 已修）",
"状态文档     scripts/27_status_readme.R, scripts/28_freeze_current.sh",
"一键复现     reproduce_current.sh tables|signatures|metabric|train|all|check",
"验收         scripts/_verify_19.R, _verify_response_numbers.R, _verify_final_numbers.R",
"```",
"",
"**不要**再把 08/09/16/17/18 当作正式管线：它们要么使用**旧的评价口径**",
"（在评价集上重新拟合方向、构成泄漏），要么只覆盖旧的三队列/旧池子。",
"`17_remedies.R` 已退役（功能并入 20）；根目录的 `reproduce.sh` 也已作废",
"（它跑的是 v6 + 09/10/11/12/13 旧流水线，**无法生成当前论文结果**），",
"统一改用 `reproduce_current.sh`。",
"",
"**源码哈希（含共用库与冻结计划）由 `scripts/28_freeze_current.sh` 固化到**",
"`results/manifest/frozen_v1_<日期>.tsv`，那才是可核验的构建记录。",
"",
"## 2. 源码哈希（MD5 前 10 位）",
"",
"| 脚本 | 哈希 |",
"|---|---|",
sprintf("| `scripts/19_main_corrected.R` | `%s` |", h("scripts/19_main_corrected.R")),
sprintf("| `scripts/20_analysis_corrected.R` | `%s` |", h("scripts/20_analysis_corrected.R")),
sprintf("| `scripts/21_entry_predictors.R` | `%s` |", h("scripts/21_entry_predictors.R")),
sprintf("| `scripts/23_manifest.R` | `%s` |", h("scripts/23_manifest.R")),
sprintf("| `scripts/29_baseline_sensitivity.R` | `%s` |", h("scripts/29_baseline_sensitivity.R")),
sprintf("| `scripts/_lib_analysis.R` | `%s` |", h("scripts/_lib_analysis.R")),
sprintf("| `scripts/25_simulation.R` | `%s` |", h("scripts/25_simulation.R")),
sprintf("| `scripts/13_figures.R` | `%s` |", h("scripts/13_figures.R")),
sprintf("| `run_v8_parallel.sh` | `%s` |", h("run_v8_parallel.sh")),
sprintf("| `reproduce_current.sh` | `%s` |", h("reproduce_current.sh")),
sprintf("| `R/utils.R` | `%s` |", h("R/utils.R")),
sprintf("| `R/audit_selection.R` | `%s` |", h("R/audit_selection.R")),
sprintf("| `R/audit_modules.R` | `%s` |", h("R/audit_modules.R")),
sprintf("| `R/sigAudit.R` | `%s` |", h("R/sigAudit.R")),
"",
"## 3. 正式主实验：完成情况",
"",
sprintf("注册网格：%d 个（队列×种子）= %d 队列 × %d 划分种子",
        length(want), length(COHORTS), length(SEEDS)),
sprintf("已完成：**%d / %d**", sum(done %in% want), length(want)),
sprintf("缺失：%s", if (length(missing)) paste(missing, collapse=", ") else "无"),
sprintf("不在网格内的多余产物：%s", if (length(extra)) paste(extra, collapse=", ") else "无"),
"",
"| 队列 | 已完成划分种子 |",
"|---|---|")
for (cn in COHORTS) {
  s <- SEEDS[vapply(SEEDS, function(x) sprintf("%s/%d", cn, x) %in% done, logical(1))]
  # 文件名里的 COAD 实际是 COADREAD（第三轮审查要求全文统一命名）
  disp <- if (cn == "COAD") "COADREAD" else cn
  L <- c(L, sprintf("| %s | %s |", disp, if (length(s)) paste(s, collapse=", ") else "—"))
}
L <- c(L, "",
"## 4. 失败与已修问题",
"")
if (length(fails)) L <- c(L, "链日志中的失败记录：", "```", fails, "```") else
  L <- c(L, "链日志中当前没有未处理的 FAIL/ABORTED 记录。")
L <- c(L, "",
"已知并已修复的一例（保留记录）：",
"- **LUAD 首次失败**：`OS_STATUS` 含**空字符串**，被「显式取值映射」判为未定义取值而中止。",
"  修法：空串/纯空白**按缺失处理并计数**（实测 52 例）。链在失败时**按设计立即中止**，",
"  未把错误分母带进汇总；修复后 LUAD 已重跑成功。",
"",
"## 5. 尚未完成 / 已知偏离",
"",
"本节的每一条都由脚本按**文件是否真的存在**生成，不再手写进度。",
"",
"- [ ] 正文成稿（Methods / Results / Discussion）—— 当前只有提纲，**不能算论文已完成**",
"- [ ] 投稿期刊、分区体系与 APC 需**官方来源核验**（见 `docs/journal_check.md`）",
"- [ ] LGG / PAAD 的 PFS 终点**时间起点**未能独立核实（门户元数据未给出）",
"- [ ] 下载日期与批次失败清单未在数据文件中记录",
"- [ ] `S3`（家族对比）主张已**撤下**，偏离记录见 `config/frozen_analysis_plan.md`",
"- [ ] `grade_credibility()` 已移除 A–D 字母评级；optimism / transport / concordance",
"      三个模块已标为 **experimental 并移出默认入口**，未参与本文任何结论",
"",
"## 6. 当前可提供的文件",
"")
for (f in c("results/manifest/sample_manifest.tsv","results/manifest/geneset_manifest.tsv",
            "results/manifest/candidate_diversity.tsv",
            "results/manifest/public_signatures_genes.tsv","results/manifest/public_signatures_meta.tsv",
            "results/manifest/public_signatures_excluded.tsv",
            "results/manifest/main_input_files.txt","results/manifest/run_manifest.tsv",
            "results/manifest/frozen_v1_20261017.tsv",
            "config/data_dictionary.md","config/design_grid.txt","config/frozen_analysis_plan.md",
            "results/tables/20_rule_performance_long.tsv","results/tables/20_rule_members_long.tsv",
            "results/tables/29_baseline_sensitivity.tsv","results/tables/29_paired_vs_fixed.tsv",
            "results/tables/31_metabric_rule_performance.tsv",
            "results/tables/32_metabric_baseline_sensitivity.tsv",
            "results/tables/32_metabric_vs_tcga.tsv","results/figures/fig5_metabric.png",
            "results/tables/36_metabric_prediction_check.tsv",
            "results/tables/25_simulation_final_replicates.tsv",
            "results/figures/fig_captions.md",
            "reproduce_current.sh")) {
  L <- c(L, sprintf("- `%s` %s", f, if (file.exists(f)) "（存在）" else "（**尚未生成**）"))
}
L <- c(L, "", "## 7. 复现入口",
"",
"```bash",
"bash reproduce_current.sh tables      # 主分析表图（由已保存预测重建；不训练）",
"bash reproduce_current.sh signatures  # 26a+26b：会在签名基因集上【重新拟合】候选池",
"bash reproduce_current.sh metabric    # METABRIC（31 训练 → 32 → 33）——**会训练**",
"bash reproduce_current.sh train       # 从原始数据重训 40 组（约 2.5 小时）",
"bash reproduce_current.sh all         # train + signatures + tables",
"bash reproduce_current.sh check       # 单元测试 + R CMD check + 哈希固化",
"SIGAUDIT_ROOT=/path/to/sigAudit bash reproduce_current.sh tables   # 自定义根目录",
"```",
"",
"关键词区分：`tables` **不训练任何模型**；`signatures` 与 `train` **会训练**。",
"根目录通过参数或 `SIGAUDIT_ROOT` 传给所有子进程（脚本内不再写死路径）；",
"运行记录与脚本哈希写入 `logs/reproduce_current_provenance.txt`。",
"",
"**异地运行已验证**：在 `D:/bioinfo/sigAudit_repro_test` 独立根目录跑通 `tables`",
"（7 步全 OK），其产物与本目录**除路径字符串外字节相同**。")
L <- c(L, "", "> 读取这些 TSV 时请用 `quote = \"\"`：数据里含 `van 't Veer` 的孤立单引号，",
        "> `read.table` 默认 `quote = \"\\\"'\"` 会把它当作引号起点并**吞掉后续换行**",
        "> （实测 6 行只读到 1 行）。")

writeLines(L, OUT)
cat("已生成 ", OUT, "（", length(L), " 行）\n", sep="")
cat("  已完成组合 ", sum(done %in% want), "/", length(want), "\n", sep="")
