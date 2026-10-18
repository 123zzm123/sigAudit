#!/usr/bin/env Rscript
# ============================================================
# 排查：LGG 里"所有模型"的 val−test C-index 落差都在 0.125~0.131，
#       这既可能是这一组 val/test 恰好难易不同（可接受，设计里已用
#       "增益之差"抵消），也可能是划分/索引错位之类会推翻全部结论的 bug。
#
# 判据：同一个 tr 训练出来的模型，对 rest 的 255 个样本做 200 次【随机二分】，
#       看"两半 C-index 之差"的分布。若观测到的固定划分落差落在这个分布内，
#       说明只是划分运气；若远远超出，则确有问题。
# ============================================================
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
cat("tr", length(tr), "va", length(va), "te", length(te), "\n")
cat("事件率: tr", round(mean(sv$e[tr]),3), " va", round(mean(sv$e[va]),3),
    " te", round(mean(sv$e[te]),3), "\n")
cat("随访中位数(月): tr", round(median(sv$t[tr]),1), " va", round(median(sv$t[va]),1),
    " te", round(median(sv$t[te]),1), "\n\n")

Z <- t(scale(t(expr))); Z[!is.finite(Z)] <- 0
G <- function(idx) t(Z[, colnames(Z)[idx], drop=FALSE])
ci <- function(t,e,s){ s<-as.numeric(s); ok<-is.finite(s)&!is.na(t)
  if(sum(ok)<15) return(NA_real_); as.numeric(summary(coxph(Surv(t[ok],e[ok])~s[ok]))$concordance[1]) }

gs <- fromJSON(file.path(PROJ,"data/signatures/msigdb_genesets.json"))
gs <- gs[sapply(gs,length)>=10 & sapply(gs,length)<=500]
set.seed(20261015); pick <- sample(names(gs), 25)

Xtr <- G(tr); Xre <- G(rest)
ttr<-sv$t[tr]; etr<-sv$e[tr]
tre<-sv$t[rest]; ere<-sv$e[rest]
i_va <- match(va, rest); i_te <- match(te, rest)   # 在 rest 内的下标（必须一致，用来防错位）

obs <- numeric(0); nullsd <- numeric(0); nullall <- numeric(0)
for (k in seq_along(pick)) {
  g <- intersect(gs[[pick[k]]], rownames(Z)); if (length(g) < 10) next
  set.seed(1)
  cv <- tryCatch(cv.glmnet(Xtr[,g,drop=FALSE], Surv(ttr,etr), family="cox",
                           alpha=1, nfolds=5), error=function(e) NULL)
  if (is.null(cv)) next
  sre <- as.numeric(predict(cv, Xre[,g,drop=FALSE], s=cv$lambda.min, type="link"))
  # 观测落差：固定划分的 va / te
  c_va <- ci(tre[i_va], ere[i_va], sre[i_va])
  c_te <- ci(tre[i_te], ere[i_te], sre[i_te])
  obs <- c(obs, c_va - c_te)
  # 零分布：把 rest 随机二分 200 次，看"两半之差"的尺度
  for (b in 1:200) {
    perm <- sample(length(rest)); h1 <- perm[1:length(va)]; h2 <- perm[(length(va)+1):length(rest)]
    d <- ci(tre[h1], ere[h1], sre[h1]) - ci(tre[h2], ere[h2], sre[h2])
    if (is.finite(d)) { nullall <- c(nullall, d); nullsd <- c(nullsd, d) }
  }
}
cat("=== 固定划分下的 val−test 落差 ===\n")
cat(sprintf("  基因集数 %d, 平均 %.4f, 中位 %.4f, 范围 [%.4f, %.4f]\n",
            length(obs), mean(obs), median(obs), min(obs), max(obs)))
cat("\n=== 随机二分下的『两半之差』零分布（同一样本量 127/128）===\n")
cat(sprintf("  样本数 %d, 均值 %.4f, SD %.4f\n", length(nullall), mean(nullall), sd(nullall)))
cat(sprintf("  2.5%%~97.5%% 区间 [%.4f, %.4f]\n",
            quantile(nullall,.025), quantile(nullall,.975)))
z <- mean(obs) / sd(nullall)
cat(sprintf("\n  观测落差折合标准差数 z = %.2f\n", z))
cat(if (abs(z) < 3) "  => 落在划分运气的正常范围内，设计正确。\n" else
    "  => 超出正常范围，必须进一步排查！\n")
