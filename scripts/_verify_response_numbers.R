# 读取辅助：必须 quote="" —— 数据里含 van 't Veer 的孤立单引号，
# read.table 默认 quote="\"' 会把它当引号起点并吞掉后续换行（实测只读到 1 行）。
read.table_plain <- function(f, ...) read.table(f, ..., quote="", comment.char="")
# 核对第二轮回应文档"第十六节"中引用的每一个数字。
ok <- 0; bad <- 0
ck <- function(cond, lab) {
  v <- isTRUE(cond)
  cat(sprintf("  [%s] %s\n", if (v) "OK" else "**错**", lab))
  if (v) ok <<- ok + 1 else bad <<- bad + 1
}
has <- function(lines, pat) any(grepl(pat, lines, fixed = TRUE))

## ---- 16.1 R 包 ----
cat("=== 16.1 R 包 ===\n")
L <- readLines("D:/bioinfo/sigAudit/results/check/00check.log", warn = FALSE)
ck(has(L, "Status: OK"), "check 日志: Status: OK")
ck(has(L, "checking examples ... OK"), "checking examples ... OK")
ck(has(L, "checking tests ... OK"), "checking tests ... OK")
ck(length(list.files("D:/bioinfo/sigAudit/man", pattern = "[.]Rd$")) == 10, "man/ 共 10 个 Rd")
ck(has(readLines("D:/bioinfo/sigAudit/man/audit_selection.Rd", warn = FALSE), "n_background"),
   "audit_selection.Rd 含 n_background")
ck(has(readLines("D:/bioinfo/sigAudit/man/audit_rules.Rd", warn = FALSE), "argmax_nonreference"),
   "audit_rules.Rd 含 argmax_nonreference")
sa <- readLines("D:/bioinfo/sigAudit/R/sigAudit.R", warn = FALSE)
ck(!any(grepl("grade <- if", sa, fixed = TRUE)), "sigAudit.R 已无 A-D 分级赋值")
ck(has(sa, "no letter grade"), "print.sigAudit_summary 明确声明不产生等级")

## ---- 16.2 模块 21 ----
cat("\n=== 16.2 模块 21 ===\n")
e21 <- "D:/bioinfo/sigAudit/results/tables/21_entry_effect.tsv"
if (!file.exists(e21)) {
  cat("  [INFO] 21 的正式输出尚未产生（需 40 组跑完后运行 21）；文中已如此标注\n")
} else {
  e <- read.table_plain(e21, sep = "\t", header = TRUE)
  mc <- e$mc_se[is.finite(e$mc_se)]
  ck(length(mc) > 0, "存在 mc_se 列")
  cat(sprintf("  [INFO] MCSE 实测范围 %.5f ~ %.5f\n", min(mc), max(mc)))
}

## ---- 16.3 模拟 ----
cat("\n=== 16.3 模拟 ===\n")
s <- read.table_plain("D:/bioinfo/sigAudit/results/tables/25_simulation_final.tsv", sep = "\t", header = TRUE)
ck(nrow(s) == 8, "共 8 个条件")
ck(setequal(unique(s$n_rep), c(100, 65, 69)), "有效重复数为 100 / 65 / 69")
r1 <- s$Orel[s$K == 6  & s$signal == "strong" & s$rho == 0]
r2 <- s$Orel[s$K == 20 & s$signal == "strong" & s$rho == 0]
ck(abs(round(r2 / r1, 1) - 2.3) < 0.06, sprintf("K=6->20 的 Orel 比值 %.2f（文中 2.3 倍）", r2 / r1))
w  <- s$Orel[s$K == 6 & s$signal == "weak" & s$rho == 0]
ck(abs(round(w / r1, 1) - 7.8) < 0.06, sprintf("弱/强信号 Orel 比值 %.2f（文中应为 7.8 倍）", w / r1))
neg <- s$dCt_ens[s$K == 6 & s$signal == "strong"]
ck(all(neg < 0), sprintf("K=6 强信号的集成 dCt 均为负: %s", paste(round(neg, 4), collapse = ", ")))
ck(all(abs(s$dCt_ens[s$K == 20 & s$signal == "strong"]) < 0.001),
   "K=20 强信号的集成 dCt 约等于零")
ck(has(readLines("D:/bioinfo/sigAudit/scripts/25_simulation.R", warn = FALSE),
       "分数层面"), "脚本内明确标注了相关性因子是分数层面实现")

## ---- 16.4 发表签名模块 ----
cat("\n=== 16.4 发表签名模块 ===\n")
m <- read.table_plain("D:/bioinfo/sigAudit/results/manifest/public_signatures_meta.tsv", sep = "\t", header = TRUE)
ck(nrow(m) == 6, "纳入 6 个签名")
ck(m$n_gene_unique[m$signature == "GENE70_MammaPrint"] == 51, "GENE70 唯一基因 = 51")
ck(m$n_probe_or_entry[m$signature == "GENE70_MammaPrint"] == 70, "GENE70 探针数 = 70")
ck(m$n_gene_unique[m$signature == "GENE76_Rotterdam"] == 60 &&
   m$n_probe_or_entry[m$signature == "GENE76_Rotterdam"] == 76, "GENE76 = 60 唯一 / 76 探针")
ck(m$n_gene_unique[m$signature == "GGI"] == 99 &&
   m$n_probe_or_entry[m$signature == "GGI"] == 128, "GGI = 99 唯一 / 128 探针")
ex <- read.table_plain("D:/bioinfo/sigAudit/results/manifest/public_signatures_excluded.tsv",
                 sep = "\t", header = TRUE)
ck(nrow(ex) == 7, "排除 7 个对象并记录原因")
a <- read.table_plain("D:/bioinfo/sigAudit/results/tables/26_signature_summary.tsv", sep = "\t", header = TRUE)
a$pct <- 100 * a$Orel / a$Gv
ck(abs(mean(a$pct) - 67.1) < 0.2, sprintf("平均 Orel/Gv = %.1f%%（文中 67.1%%）", mean(a$pct)))
ck(abs(min(a$pct) - 54.7) < 0.2, sprintf("最小 = %.1f%%（文中 54.7%%）", min(a$pct)))
ck(abs(max(a$pct) - 80.5) < 0.2, sprintf("最大 = %.1f%%（文中 80.5%%）", max(a$pct)))
ck(nrow(read.table_plain("D:/bioinfo/sigAudit/results/tables/26_signature_tasks.tsv", sep = "\t", header = TRUE)) == 240,
   "签名任务单元数 = 240（6 签名 × 8 队列 × 5 划分）")

## ---- 16.5 LUAD ----
cat("\n=== 16.5 LUAD 修复 ===\n")
f <- "D:/bioinfo/sigAudit/results/main/19_LUAD_20261014_100.rds"
if (!file.exists(f)) {
  cat("  [INFO] LUAD 产物尚未生成（重跑中）\n")
} else {
  R <- readRDS(f)
  ck(isTRUE(R$completed), "LUAD 产物 completed = TRUE")
  ck(length(R$sample_ids) == 497, sprintf("样本数 %d（文中 497）", length(R$sample_ids)))
  g <- R$genesets[[1]]
  ck(sum(g$e_tr) + sum(g$eva) + sum(g$ete) == 180, "总事件数 180")
}

cat(sprintf("\n核对结果: %d 项一致, %d 项不一致\n", ok, bad))
if (bad > 0) quit(status = 1)
