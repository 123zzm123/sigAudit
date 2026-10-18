#!/usr/bin/env Rscript
# 交叉验证：包里的 audit_rules() 是否复现 12_rules.R 的 LGG 结果？
# 12_rules.R (LGG): none 0.773/0.648 | argmax_val 0.797/0.650 +0.0023 |
#                   runner_up 0.792/0.651 +0.0035 | linear_only 0.794/0.644 -0.0036 |
#                   tree_only 0.778/0.660 +0.0121 | oracle 0.771/0.682 +0.0348
suppressMessages(source("D:/bioinfo/sigAudit/R/audit_modules.R"))
PROJ <- "D:/bioinfo/sigAudit"
POOL_NAMES <- c("Lasso","Ridge","Enet0.9","Enet0.7","Enet0.5","Enet0.3","Enet0.1",
                "Ridge_1se","UniCox+Lasso","LassoTop100","LassoTop50","LassoTop25",
                "LassoTop10","Enet0.5Top25","RSF100","GBM100","GBM250","RSF250")
FAM <- setNames(ifelse(POOL_NAMES %in% c("RSF100","GBM100","GBM250","RSF250"),"tree","linear"), POOL_NAMES)
D <- read.table(file.path(PROJ,"results/tables/09_LGG_100.tsv"), header=TRUE, sep="\t", check.names=FALSE)
D <- D[D$model %in% POOL_NAMES, ]
tab <- table(D$geneset); keep <- names(tab)[tab == length(POOL_NAMES)]
D <- D[D$geneset %in% keep, ]
mk <- function(v) { M <- tapply(D[[v]], list(D$geneset, D$model), function(x) x[1])
  M <- M[keep, POOL_NAMES, drop=FALSE]
  matrix(as.numeric(M), nrow=nrow(M), dimnames=dimnames(M)) }
V <- mk("val_ci"); Tt <- mk("test_ci")
ok <- complete.cases(V) & complete.cases(Tt)
r <- audit_rules(V[ok,,drop=FALSE], Tt[ok,,drop=FALSE], reference="Lasso", family=FAM)
print(r)
