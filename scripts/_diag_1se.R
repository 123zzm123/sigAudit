#!/usr/bin/env Rscript
# 诊断：lambda.1se 是否在小样本高维 Cox 里退化成空模型（系数全 0 -> C-index=0.500）
# 若退化比例很高，则 v5 模型池里 6/19 个变体是"死模型"，剂量-反应曲线被污染。
suppressMessages({library(survival); library(glmnet); library(jsonlite)})
options(stringsAsFactors = FALSE)
PROJ <- "D:/bioinfo/sigAudit"

expr <- as.matrix(read.table("D:/bioinfo/BRCA_TLS/data/cancers/LGG_broad_expr.tsv",
                             header=TRUE, sep="\t", row.names=1, check.names=FALSE))
if (nrow(expr) < ncol(expr)) expr <- t(expr)
expr <- log2(expr + 1); expr <- expr[apply(expr,1,sd)>0, , drop=FALSE]
sv <- read.table("D:/bioinfo/BRCA_TLS/data/cancers/LGG_broad_surv.tsv",
                 header=TRUE, sep="\t", row.names=1, check.names=FALSE)
tcol <- grep("_MONTHS$", names(sv), value=TRUE)[1]
ecol <- sub("_MONTHS$", "_STATUS", tcol)
sv$t <- suppressWarnings(as.numeric(sv[[tcol]]))
sv$e <- ifelse(grepl("^1", sv[[ecol]]), 1L, 0L)
sv <- sv[!is.na(sv$t) & sv$t > 0, ]
colnames(expr) <- substr(colnames(expr),1,12); expr <- expr[,!duplicated(colnames(expr))]
com <- intersect(colnames(expr), rownames(sv)); expr <- expr[,com]; sv <- sv[com,]

set.seed(20261014)
ev <- which(sv$e==1); ce <- which(sv$e==0)
tr <- sort(c(sample(ev, floor(.5*length(ev))), sample(ce, floor(.5*length(ce)))))
rest <- setdiff(seq_along(sv$e), tr)
re <- rest[sv$e[rest]==1]; rc <- rest[sv$e[rest]==0]
va <- sort(c(sample(re, floor(.5*length(re))), sample(rc, floor(.5*length(rc)))))
te <- setdiff(rest, va)
Z <- t(scale(t(expr))); Z[!is.finite(Z)] <- 0
G <- function(idx) t(Z[, colnames(Z)[idx], drop=FALSE])
Xtr <- G(tr); Xva <- G(va); Xte <- G(te)
ttr<-sv$t[tr]; etr<-sv$e[tr]; tva<-sv$t[va]; eva<-sv$e[va]; tte<-sv$t[te]; ete<-sv$e[te]
ci <- function(t,e,s){ s<-as.numeric(s); ok<-is.finite(s)&!is.na(t)
  if(sum(ok)<15) return(NA_real_); as.numeric(summary(coxph(Surv(t[ok],e[ok])~s[ok]))$concordance[1]) }

gs <- fromJSON(file.path(PROJ,"data/signatures/msigdb_genesets.json"))
gs <- gs[sapply(gs,length)>=10 & sapply(gs,length)<=500]
set.seed(20261015); pick <- sample(names(gs), 10)

cat(sprintf("%-10s %-4s %8s %8s %8s %8s\n","alpha","lam","n_nzero","val_ci","test_ci","val_ci_min"))
agg <- list()
for (nm in pick) {
  g <- intersect(gs[[nm]], rownames(Z)); if (length(g)<10) next
  a <- Xtr[,g,drop=FALSE]; b <- Xva[,g,drop=FALSE]; c2 <- Xte[,g,drop=FALSE]
  y <- Surv(ttr, etr)
  for (al in c(1, 0.5, 0.1, 0)) {
    set.seed(1); cv <- tryCatch(cv.glmnet(a, y, family="cox", alpha=al, nfolds=5), error=function(e) NULL)
    if (is.null(cv)) next
    for (lm in c("min","1se")) {
      l <- if (lm=="min") cv$lambda.min else cv$lambda.1se
      co <- as.numeric(coef(cv, s=l))
      nz <- sum(co != 0)
      pv <- as.numeric(predict(cv, b, s=l, type="link"))
      pt <- as.numeric(predict(cv, c2, s=l, type="link"))
      key <- sprintf("a=%.1f_%s", al, lm)
      agg[[key]] <- c(agg[[key]], c(nz, ci(tva,eva,pv), ci(tte,ete,pt)))
    }
  }
}
cat("\n===== 10 个基因集汇总 =====\n")
cat(sprintf("%-14s %10s %10s %10s %10s %10s\n","模型","平均n_nzero","零系数占比","平均val","平均test","val=0.5占比"))
for (k in names(agg)) {
  m <- matrix(agg[[k]], ncol=3, byrow=TRUE)
  cat(sprintf("%-14s %10.1f %10.0f%% %10.3f %10.3f %10.0f%%\n", k,
      mean(m[,1]), 100*mean(m[,1]==0), mean(m[,2],na.rm=TRUE), mean(m[,3],na.rm=TRUE),
      100*mean(abs(m[,2]-0.5)<1e-9, na.rm=TRUE)))
}
