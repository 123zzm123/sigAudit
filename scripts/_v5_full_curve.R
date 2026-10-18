#!/usr/bin/env Rscript
# 用 v5 已落盘的表，把完整的 k=2..19 曲线补出来。
# 动机：之前只打印了 k=2,3,5,8,12,16,19，把"跳升在 k=16"当成了结论。
#       v5 池子第 14 位是 Ridge_1se —— 若跳升其实在 k=14，则整个
#       "树集成导致跳升"的解释是抽样打印造成的错觉。
PROJ <- "D:/bioinfo/sigAudit"
V5 <- c("Lasso","Ridge","Lasso_1se","Enet0.5","Enet0.5_1se","Enet0.9","Enet0.7",
        "Enet0.3","Enet0.1","Enet0.9_1se","Enet0.7_1se","Enet0.3_1se","Enet0.1_1se",
        "Ridge_1se","UniCox+Lasso","RSF100","GBM100","GBM250","RSF250")
files <- list.files(file.path(PROJ,"results/tables"), pattern="^08_[A-Z]+_100\\.tsv$", full.names=TRUE)
cat(sprintf("%-6s %5s | %s\n", "队列", "n", paste(sprintf("k=%-2d", 2:19), collapse=" ")))
res <- list()
for (f in files) {
  coh <- sub("^08_(.+)_100\\.tsv$","\\1", basename(f))
  D <- read.table(f, header=TRUE, sep="\t", check.names=FALSE)
  cols <- paste0("bias_k", 2:19)
  cols <- cols[cols %in% names(D)]
  m <- colMeans(D[, cols, drop=FALSE], na.rm=TRUE)
  n <- nrow(D)
  res[[coh]] <- m
  cat(sprintf("%-6s %5d | %s\n", coh, n, paste(sprintf("%+.4f", m), collapse=" ")))
}
cat("\n各队列的【相邻增量】(后一 k 减前一 k)，找跳升点：\n")
for (coh in names(res)) {
  m <- res[[coh]]
  d <- diff(as.numeric(m))
  j <- which.max(d)
  cat(sprintf("  %-6s  最大增量出现在 k=%d->%d : %+.4f   该处模型 = %s\n",
              coh, 1+j, 2+j, d[j], V5[2+j]))
}
cat("\nv5 池子逐位对照：\n")
for (i in seq_along(V5)) cat(sprintf("  k=%-2d  %s\n", i, V5[i]))
