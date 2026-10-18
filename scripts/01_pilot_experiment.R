#!/usr/bin/env Rscript
# ============================================================
# 01 — 实验：量化"模型选择规则"造成的性能虚高
#
# 研究问题：标准 101 组合流程按「训练+测试均值」选模型，
#           是否系统性选出过拟合模型？虚高多少？
#
# 对照（同一数据/基因集/模型池，只改选法）：
#   A  mean   : 按 (Train+Test)/2 最高选 —— 流行做法
#   B  innercv: 只用训练集内部 CV 选      —— 诚实做法
#   C  oracle : 按 Test 最高选            —— 上界（实际不可用）
#
# 性能优化（v2）：
#   * RSF ntree 300->150, GBM n.trees 300->150
#   * 内部 CV 折数 5->3
#   * 每 25 个基因集增量存盘（v1 因只在最后存盘，被杀后丢了 25 个结果）
# ============================================================
suppressMessages({library(survival); library(glmnet); library(randomForestSRC); library(gbm); library(jsonlite)})
options(stringsAsFactors = FALSE)
PROJ <- "D:/bioinfo/sigAudit"
LOGDIR <- file.path(PROJ, "logs"); dir.create(LOGDIR, showWarnings = FALSE, recursive = TRUE)
ARGS <- commandArgs(trailingOnly = TRUE)
N_GS   <- if (length(ARGS) >= 1) as.integer(ARGS[1]) else 200
COHORT <- if (length(ARGS) >= 2) ARGS[2] else "LGG"
TAG    <- paste0("01_", COHORT, "_", N_GS)

# 日志路径必须带队列标签 —— 之前不同队列共用一个 01_pilot.log，
# 手动测试会把正在跑的任务日志覆盖掉（实际发生过）。
# 注意: 必须放在 ARGS 解析之后，否则 COHORT/N_GS 尚未定义。
LOG <- file.path(LOGDIR, paste0(TAG, ".log")); con <- file(LOG, "wt")
say <- function(...) { m <- paste0(...); cat(m,"\n",sep=""); cat(m,"\n",sep="",file=con); flush(con) }

# ---- 数据 ----
say("载入 ", COHORT, " 队列 ...")
# 注意：R 的 else 必须与前一行的 } 或表达式在同一行，换行会报
# "unexpected 'else'" —— 本项目已踩过一次同类语法坑
if (COHORT == "LGG") {
  cf <- "D:/bioinfo/BRCA_TLS/data/cancers/LGG_broad_expr.tsv"
  sf <- "D:/bioinfo/BRCA_TLS/data/cancers/LGG_broad_surv.tsv"
} else {
  cf <- sprintf("%s/data/cohorts/%s_expr.tsv", PROJ, COHORT)
  sf <- sprintf("%s/data/cohorts/%s_surv.tsv", PROJ, COHORT)
}
expr <- as.matrix(read.table(cf, header=TRUE, sep="\t", row.names=1, check.names=FALSE))
if (nrow(expr) < ncol(expr)) expr <- t(expr)          # -> 基因 × 样本
expr <- log2(expr + 1)
expr <- expr[apply(expr,1,sd)>0, , drop=FALSE]
sv <- read.table(sf, header=TRUE, sep="\t", row.names=1, check.names=FALSE)
tcol <- grep("_MONTHS$", names(sv), value=TRUE)[1]
ecol <- sub("_MONTHS$", "_STATUS", tcol)
sv$t <- suppressWarnings(as.numeric(sv[[tcol]]))
sv$e <- ifelse(grepl("^1", sv[[ecol]]), 1L, ifelse(grepl("^0", sv[[ecol]]), 0L, NA_integer_))
sv <- sv[!is.na(sv$t) & !is.na(sv$e) & sv$t > 0, ]
colnames(expr) <- substr(colnames(expr),1,12)          # 标准步骤：统一 patient 级
expr <- expr[, !duplicated(colnames(expr)), drop=FALSE]
com <- intersect(colnames(expr), rownames(sv))
expr <- expr[, com]; sv <- sv[com, ]
say("  ", COHORT, ": ", ncol(expr), " 样本 × ", nrow(expr), " 基因, 事件 ", sum(sv$e),
    "  终点 ", tcol)

d <- data.frame(sampleID=com, t=sv$t, e=sv$e)
set.seed(20261014)
idx_ev <- which(d$e==1); idx_cen <- which(d$e==0)
tr <- sort(c(sample(idx_ev, floor(.7*length(idx_ev))), sample(idx_cen, floor(.7*length(idx_cen)))))
te <- setdiff(seq_len(nrow(d)), tr)
say("  训练 ", length(tr), " / 测试 ", length(te))

Z <- t(scale(t(expr))); Z[!is.finite(Z)] <- 0
Etr <- t(Z[, d$sampleID[tr]]); Ete <- t(Z[, d$sampleID[te]])
ttr <- d$t[tr]; etr <- d$e[tr]; tte <- d$t[te]; ete <- d$e[te]

# 小队列（如 MESO 测试集仅 26 例）会被过严的门槛整批跳过；放宽到 15/5，
# 但保留最低要求以免 C-index 完全不稳定
ci <- function(t,e,s) { ok <- !is.na(s) & is.finite(s)
  if (sum(ok)<15 || sum(e[ok])<5) return(NA_real_)
  as.numeric(summary(coxph(Surv(t[ok],e[ok]) ~ s[ok]))$concordance[1]) }

NTREE <- 150; NFOLD_INNER <- 3          # << 优化参数

fit_models <- function(Xtr, ttr, etr, Xte) {
  y <- Surv(ttr, etr); res <- list()
  for (a in c(1, 0.5, 0)) {              # Lasso / Enet / Ridge
    nm <- c("Lasso","Enet","Ridge")[match(a, c(1,0.5,0))]
    set.seed(1); cv <- tryCatch(cv.glmnet(Xtr, y, family="cox", alpha=a, nfolds=5), error=function(e) NULL)
    if (!is.null(cv)) res[[nm]] <- list(
      tr=as.numeric(predict(cv, Xtr, s=cv$lambda.min, type="link")),
      te=as.numeric(predict(cv, Xte, s=cv$lambda.min, type="link")))
  }
  set.seed(1); rf <- tryCatch(rfsrc(Surv(ttr,etr) ~ ., data=data.frame(Surv(ttr,etr), Xtr),
      ntree=NTREE, nodesize=5, splitrule="logrank"), error=function(e) NULL)
  if (!is.null(rf)) res[["RSF"]] <- list(tr=as.numeric(rf$predicted),
                                         te=as.numeric(predict(rf, data.frame(Xte))$predicted))
  set.seed(1); gb <- tryCatch(gbm(Surv(ttr,etr) ~ ., data=data.frame(Surv(ttr,etr), Xtr),
      distribution="coxph", n.trees=NTREE, shrinkage=0.01, interaction.depth=3,
      n.minobsinnode=10, verbose=FALSE), error=function(e) NULL)
  if (!is.null(gb)) res[["GBM"]] <- list(tr=as.numeric(predict(gb, data.frame(Xtr), n.trees=NTREE, type="link")),
                                         te=as.numeric(predict(gb, data.frame(Xte), n.trees=NTREE, type="link")))
  up <- apply(Xtr, 2, function(x){ s<-tryCatch(summary(coxph(Surv(ttr,etr)~x))$coef, error=function(e) NULL)
    if (is.null(s)) NA_real_ else s[1,"Pr(>|z|)"] })
  keep <- which(!is.na(up) & up < 0.05)
  if (length(keep) >= 5) {
    set.seed(1); cv6 <- tryCatch(cv.glmnet(Xtr[,keep,drop=FALSE], y, family="cox", alpha=1, nfolds=5), error=function(e) NULL)
    if (!is.null(cv6)) res[["UniCox+Lasso"]] <- list(
      tr=as.numeric(predict(cv6, Xtr[,keep,drop=FALSE], s=cv6$lambda.min, type="link")),
      te=as.numeric(predict(cv6, Xte[,keep,drop=FALSE], s=cv6$lambda.min, type="link")))
  }
  res
}

inner_cv <- function(Xtr, ttr, etr) {
  set.seed(7); folds <- sample(rep(1:NFOLD_INNER, length.out=length(ttr)))
  sc <- rep(NA_real_, length(ttr))
  for (k in 1:NFOLD_INNER) {
    i <- folds==k
    if (sum(etr[i])<5 || sum(etr[!i])<10) next
    m <- fit_models(Xtr[!i,,drop=FALSE], ttr[!i], etr[!i], Xtr[i,,drop=FALSE])
    if (!length(m)) next
    cavg <- sapply(m, function(x) ci(ttr[!i], etr[!i], x$tr))
    sc[i] <- m[[names(which.max(cavg))]]$te
  }
  ci(ttr, etr, sc)
}

gs <- fromJSON(file.path(PROJ,"data/signatures/msigdb_genesets.json"))
gs <- gs[sapply(gs, length) >= 10 & sapply(gs, length) <= 500]
say("可用基因集: ", length(gs))
set.seed(20261015); pick <- sample(names(gs), min(N_GS, length(gs)))
say("本次评估 ", length(pick), " 个基因集 (队列 ", COHORT, ")")

HDR <- c("geneset","n_gene","sel_mean","c_mean","sel_oracle","c_oracle","c_innercv","max_train","max_test")
OUT <- file.path(PROJ, "results/tables", paste0(TAG, ".tsv"))
write.table(as.data.frame(matrix(NA,0,length(HDR),dimnames=list(NULL,HDR))), OUT,
            sep="\t", quote=FALSE, row.names=FALSE)

res <- list(); t0 <- Sys.time()
for (i in seq_along(pick)) {
  nm <- pick[i]; g <- intersect(gs[[nm]], rownames(Z))
  if (length(g) < 10) next
  Xtr <- Etr[, g, drop=FALSE]; Xte <- Ete[, g, drop=FALSE]
  ms <- tryCatch(fit_models(Xtr, ttr, etr, Xte), error=function(e) NULL)
  if (is.null(ms) || length(ms) < 3) next
  tab <- data.frame(model=names(ms),
                    Train=sapply(ms, function(x) ci(ttr,etr,x$tr)),
                    Test =sapply(ms, function(x) ci(tte,ete,x$te)))
  tab <- tab[!is.na(tab$Test), ]; if (nrow(tab) < 3) next
  tab$Mean <- (tab$Train + tab$Test)/2
  selA <- tab$model[which.max(tab$Mean)]; selC <- tab$model[which.max(tab$Test)]
  cB <- tryCatch(inner_cv(Xtr,ttr,etr), error=function(e) NA_real_)
  res[[nm]] <- data.frame(geneset=nm, n_gene=length(g), sel_mean=selA, c_mean=tab$Test[tab$model==selA][1],
                          sel_oracle=selC, c_oracle=tab$Test[tab$model==selC][1], c_innercv=cB,
                          max_train=max(tab$Train), max_test=max(tab$Test))
  # 增量存盘
  if (i %% 10 == 0) {
    R0 <- do.call(rbind, res)
    write.table(R0, OUT, sep="\t", quote=FALSE, row.names=FALSE)
    say("  进度 ", i, "/", length(pick), "  已完成 ", nrow(R0),
        "  用时 ", round(as.numeric(difftime(Sys.time(),t0,units="min")),1), " 分钟")
  }
}
R <- do.call(rbind, res)
write.table(R, OUT, sep="\t", quote=FALSE, row.names=FALSE)

say("")
say("========== 结果 (", COHORT, ") ==========")
say("完成基因集: ", nrow(R))
say(sprintf("选法A(训练+测试均值) 选出模型的测试集 C-index = %.3f", mean(R$c_mean, na.rm=TRUE)))
say(sprintf("选法B(训练集内部CV)  的测试集 C-index       = %.3f", mean(R$c_innercv, na.rm=TRUE)))
say(sprintf("选法C(oracle 上界)   的测试集 C-index       = %.3f", mean(R$c_oracle, na.rm=TRUE)))
say(sprintf("训练集最高 C-index 平均 = %.3f  ← 论文常报的数", mean(R$max_train, na.rm=TRUE)))
say(sprintf("测试集最高 C-index 平均 = %.3f  ← 真实上界",     mean(R$max_test, na.rm=TRUE)))
say(sprintf("A 相对 B(诚实) 的虚高 = %.3f", mean(R$c_mean - R$c_innercv, na.rm=TRUE)))
close(con)
