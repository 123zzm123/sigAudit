#!/usr/bin/env Rscript
# ============================================================
# 真实数据示例：对已保存的 METABRIC 预测运行预测级审计
#
# 只调用【安装后的 sigAudit 包】audit_predictions()，不 source 研究脚本、
# 不重新训练任何模型。输入是 inst/extdata 中一个真实的
# (签名 × 划分) 任务：18 个候选风险分数 + 验证/测试生存结局。
# 数据由 inst/examples/export_audit_example.R 从本地冻结预测导出（见该文件）；
# 选择规则固定为"签名与种子按字典序最小"，与效果无关。
#
# 用法：
#   Rscript inst/examples/audit_predictions_metabric.R
#   （需先安装包：R CMD INSTALL .  或  remotes::install_local(".")）
# ============================================================
suppressMessages(library(sigAudit))

## ---- 定位示例数据：优先用已安装包，其次回退到源码树 ----
dir <- system.file("extdata", package = "sigAudit")
if (!nzchar(dir)) dir <- "inst/extdata"          # 从源码树直接运行时
rd <- function(f) read.delim(file.path(dir, f), check.names = FALSE)

val <- rd("audit_example_val.tsv"); rownames(val) <- val$id; val$id <- NULL
tst <- rd("audit_example_test.tsv"); rownames(tst) <- tst$id; tst$id <- NULL
sur <- rd("audit_example_survival.tsv")
sv <- sur[sur$set == "validation", ]; st <- sur[sur$set == "test", ]
stopifnot(identical(rownames(val), sv$id), identical(rownames(tst), st$id),
          identical(colnames(val), colnames(tst)))

task <- list(
  name = "METABRIC|EndoPredict_11",
  val_pred = as.matrix(val), test_pred = as.matrix(tst),
  time_val = sv$time, event_val = sv$event,
  time_test = st$time, event_test = st$event)

## ---- 一次调用：候选诊断 + 固定参照比较 + 规则结果 ----
t0 <- Sys.time()
res <- audit_predictions(task, reference = "Lasso",
                         rules = c("none", "argmax_val", "top5_ensemble"))
dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

cat("== sigAudit on one real METABRIC task ==\n")
cat(sprintf("candidates: %d   n_val: %d   n_test: %d\n\n",
            res$diagnostics$n_candidates, nrow(val), nrow(tst)))
print(res)
cat(sprintf("\naudit wall time: %.2f s\n", dt))
cat("  (one audit of the requested rules on the already-loaded task;\n")
cat("   excludes file reading and model fitting)\n")

## ---- 更换固定参照：哪些 Cv/Ct 不变，哪些增益改变 ----
resR <- audit_predictions(task, reference = "Ridge",
                          rules = c("none", "argmax_val", "top5_ensemble"))
cat("\n-- the same prediction, a different fixed reference --\n")
cat(sprintf("%-14s %10s %10s %10s %10s\n", "rule", "Cv(Lasso)", "Ct(Lasso)", "Cv(Ridge)", "Ct(Ridge)"))
a <- res$rules; b <- resR$rules
for (i in seq_len(nrow(a))) {
  j <- which(b$rule == a$rule[i])
  cat(sprintf("%-14s %10.4f %10.4f %10.4f %10.4f\n",
              a$rule[i], a$Cv[i], a$Ct[i], b$Cv[j], b$Ct[j]))
}
cat("\nCv/Ct belong to the selected candidate, so argmax_val and top5_ensemble keep the same\n")
cat("Cv/Ct while the candidate pool is unchanged; `none` *is* the current reference, so its\n")
cat("Cv/Ct move with the reference. Only the gains and Orel are reference-adjusted:\n")
cat(sprintf("  argmax_val: Orel(Lasso) = %+.4f   Orel(Ridge) = %+.4f\n",
            a$Orel[a$rule == "argmax_val"], b$Orel[b$rule == "argmax_val"]))
cat(sprintf("  top5_ensemble: Orel(Lasso) = %+.4f   Orel(Ridge) = %+.4f\n",
            a$Orel[a$rule == "top5_ensemble"], b$Orel[b$rule == "top5_ensemble"]))
cat("\nNote: a single task is one unit; this example demonstrates the software path, not an\n")
cat("effect estimate (see the manuscript for the cohort-level analysis).\n")
