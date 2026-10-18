#!/usr/bin/env Rscript
# 交叉验证：包里的 audit_selection() 是否复现 10_analysis.R 的独立实现？
# 若两者数值一致，说明方法实现无分歧（同一估计量在两个代码路径下等价）。
suppressMessages(source("D:/bioinfo/sigAudit/R/audit_modules.R"))
PROJ <- "D:/bioinfo/sigAudit"
POOL_NAMES <- c("Lasso","Ridge","Enet0.9","Enet0.7","Enet0.5","Enet0.3","Enet0.1",
                "Ridge_1se","UniCox+Lasso","LassoTop100","LassoTop50","LassoTop25",
                "LassoTop10","Enet0.5Top25","RSF100","GBM100","GBM250","RSF250")

D <- read.table(file.path(PROJ,"results/tables/09_LGG_100.tsv"),
                header=TRUE, sep="\t", check.names=FALSE)
D <- D[D$model %in% POOL_NAMES, ]
tab <- table(D$geneset); keep <- names(tab)[tab == length(POOL_NAMES)]
D <- D[D$geneset %in% keep, ]
mk <- function(v) { M <- tapply(D[[v]], list(D$geneset, D$model), function(x) x[1])
  M <- M[keep, POOL_NAMES, drop=FALSE]
  matrix(as.numeric(M), nrow=nrow(M), dimnames=dimnames(M)) }
V <- mk("val_ci"); Tt <- mk("test_ci")
ok <- complete.cases(V) & complete.cases(Tt)
V <- V[ok,,drop=FALSE]; Tt <- Tt[ok,,drop=FALSE]
cat("基因集数:", nrow(V), " 候选数:", ncol(V), "\n\n")

r <- audit_selection(V, Tt, reference = "Lasso")
print(r)

cat("\n--- 与 10_analysis.R 独立实现对照 ---\n")
cat("  10_analysis: Ridge_1se 进入效应 +0.0099, r = +0.680 (n=289)\n")
cat(sprintf("  audit_selection: Ridge_1se 进入效应 %+.4f, r = %+.3f (n=%d)\n",
    r$entry_summary$mean[r$entry_summary$target=="Ridge_1se"],
    r$cor_optimism$estimate, nrow(r$entry)))
a <- r$entry_summary
cat("\n--- 排名前 5 ---\n"); print(head(a, 5), row.names=FALSE)
