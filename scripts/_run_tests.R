#!/usr/bin/env Rscript
# 统一测试入口：source R/ 下全部源文件后，用 testthat 跑 tests/testthat。
# 任何 failed / error 都会以非零码退出（评审要求：断言式，不能只 cat 打印）。
# 注意：文件匹配用 "[.]R$" 而不是反斜杠转义，避免不同 shell 下的转义差异。
suppressMessages(library(testthat))
src <- sort(list.files("D:/bioinfo/sigAudit/R", pattern = "[.]R$", full.names = TRUE))
cat("加载源文件:\n"); for (f in src) cat("  ", basename(f), "\n", sep = "")
for (f in src) source(f)
res <- test_dir("D:/bioinfo/sigAudit/tests/testthat", reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("\n---- 汇总 ----\n  test_that 块 %d | passed %d | failed %d | error %d | warning %d | skipped %d\n",
            nrow(df), sum(df$passed), sum(df$failed), sum(df$error), sum(df$warning), sum(df$skipped)))
if (sum(df$failed) + sum(df$error) > 0) quit(status = 1)
cat("全部通过（退出码 0）\n")
