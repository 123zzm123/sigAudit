#!/usr/bin/env Rscript
# ============================================================
# 05 — 实验 v2：三分法，直接测量「模型选择偏倚」
#
# v1 的设计缺陷（已在 01 中暴露）：
#   * 对比的 A(均值选) 与 B(内部CV) 是两条不同规则，差值不纯
#   * "oracle" 也在测试集上算，本身同样乐观
#
# v2 设计（干净）：
#   数据 -> 训练 / 验证 / 测试  (50/25/25)
#   标准做法：用【训练+验证】挑选模型，报告【验证集】C-index  <- 论文里的数
#   真实性能：同一模型在【完全未参与任何选择的测试集】上的 C-index
#   >>> inflation = 报告值 - 真实值   <- 这就是纯粹的模型选择偏倚
#
# 另外记录一个无选择的诚实基线（固定用 Lasso），用于回答
# 「不挑模型会损失多少」。
# ============================================================
suppressMessages({library(survival); library(glmnet); library(randomForestSRC); library(gbm); library(jsonlite)})
options(stringsAsFactors = FALSE)
PROJ <- "D:/bioinfo/sigAudit"
LOGDIR <- file.path(PROJ, "logs"); dir.create(LOGDIR, showWarnings = FALSE, recursive = TRUE)

ARGS <- commandArgs(trailingOnly = TRUE)
N_GS   <- if (length(ARGS) >= 1) as.integer(ARGS[1]) else 100
COHORT <- if (length(ARGS) >= 2) ARGS[2] else "LGG"
TAG    <- sprintf("05_%s_%d", COHORT, N_GS)
LOG    <- file.path(LOGDIR, paste0(TAG, ".log")); con <- file(LOG, "wt")
say <- function(...) { m <- paste0(...); cat(m,"\n",sep=""); cat(m,"\n",sep="",file=con); flush(con) }

# ---- 数据 ----
if (COHORT == "LGG") {
  cf <- "D:/bioinfo/BRCA_TLS/data/cancers/LGG_broad_expr.tsv"
  sf <- "D:/bioinfo/BRCA_TLS/data/cancers/LGG_broad_surv.tsv"
} else {
  cf <- sprintf("%s/data/cohorts/%s_expr.tsv", PROJ, COHORT)
  sf <- sprintf("%s/data/cohorts/%s_surv.tsv", PROJ, COHORT)
}
say("载入 ", COHORT)
expr <- as.matrix(read.table(cf, header=TRUE, sep="\t", row.names=1, check.names=FALSE))
if (nrow(expr) < ncol(expr)) expr <- t(expr)
expr <- log2(expr + 1); expr <- expr[apply(expr,1,sd)>0, , drop=FALSE]
sv <- read.table(sf, header=TRUE, sep="\t", row.names=1, check.names=FALSE)
tcol <- grep("_MONTHS$", names(sv), value=TRUE)[1]
ecol <- sub("_MONTHS$", "_STATUS", tcol)
sv$t <- suppressWarnings(as.numeric(sv[[tcol]]))
sv$e <- ifelse(grepl("^1", sv[[ecol]]), 1L, ifelse(grepl("^0", sv[[ecol]]), 0L, NA_integer_))
sv <- sv[!is.na(sv$t) & !is.na(sv$e) & sv$t > 0, ]
colnames(expr) <- substr(colnames(expr),1,12)
expr <- expr[, !duplicated(colnames(expr)), drop=FALSE]
com <- intersect(colnames(expr), rownames(sv)); expr <- expr[,com]; sv <- sv[com,]
say("  ", ncol(expr), " 样本 × ", nrow(expr), " 基因, 事件 ", sum(sv$e), " (", tcol, ")")

# ---- 三分法 50/25/25（按事件分层）----
set.seed(20261014)
ev <- which(sv$e==1); ce <- which(sv$e==0)
tr <- sort(c(sample(ev, floor(.5*length(ev))), sample(ce, floor(.5*length(ce)))))
rest <- setdiff(seq_along(sv$e), tr)
rest_ev <- rest[sv$e[rest]==1]; rest_ce <- rest[sv$e[rest]==0]
va <- sort(c(sample(rest_ev, floor(.5*length(rest_ev))), sample(rest_ce, floor(.5*length(rest_ce)))))
te <- setdiff(rest, va)
say("  训练 ", length(tr), " / 验证 ", length(va), " / 测试 ", length(te))
say("  事件: ", sum(sv$e[tr]), " / ", sum(sv$e[va]), " / ", sum(sv$e[te]))

Z <- t(scale(t(expr))); Z[!is.finite(Z)] <- 0
G <- function(idx) t(Z[, colnames(Z)[idx]])          # -> 样本 × 基因

# 防御性 C-index：长度不符/类型异常一律返回 NA，而不是让 sum() 变 NA 后崩掉
ci <- function(t, e, s) {
  s <- tryCatch(as.numeric(s), error = function(x) NULL)
  if (is.null(s) || length(s) != length(t)) return(NA_real_)
  ok <- is.finite(s) & !is.na(t) & !is.na(e)
  if (sum(ok) < 15 || sum(e[ok]) < 5) return(NA_real_)
  r <- tryCatch(as.numeric(summary(coxph(Surv(t[ok], e[ok]) ~ s[ok]))$concordance[1]),
                error = function(x) NA_real_)
  r
}

NTREE <- 150
fit_models <- function(Xa, ta, ea, Xb) {
  y <- Surv(ta, ea); res <- list()
  for (a in c(1, 0.5, 0)) {
    nm <- c("Lasso","Enet","Ridge")[match(a, c(1,0.5,0))]
    set.seed(1); cv <- tryCatch(cv.glmnet(Xa, y, family="cox", alpha=a, nfolds=5), error=function(e) NULL)
    if (!is.null(cv)) res[[nm]] <- list(fit=cv,
      pa=as.numeric(predict(cv, Xa, s=cv$lambda.min, type="link")),
      pb=as.numeric(predict(cv, Xb, s=cv$lambda.min, type="link")))
  }
  set.seed(1); rf <- tryCatch(rfsrc(Surv(ta,ea) ~ ., data=data.frame(Surv(ta,ea), Xa),
      ntree=NTREE, nodesize=5, splitrule="logrank"), error=function(e) NULL)
  if (!is.null(rf)) res[["RSF"]] <- list(pa=as.numeric(rf$predicted),
                                         pb=as.numeric(predict(rf, data.frame(Xb))$predicted))
  set.seed(1); gb <- tryCatch(gbm(Surv(ta,ea) ~ ., data=data.frame(Surv(ta,ea), Xa),
      distribution="coxph", n.trees=NTREE, shrinkage=0.01, interaction.depth=3,
      n.minobsinnode=10, verbose=FALSE), error=function(e) NULL)
  if (!is.null(gb)) res[["GBM"]] <- list(pa=as.numeric(predict(gb, data.frame(Xa), n.trees=NTREE, type="link")),
                                         pb=as.numeric(predict(gb, data.frame(Xb), n.trees=NTREE, type="link")))
  up <- apply(Xa, 2, function(x){ s<-tryCatch(summary(coxph(Surv(ta,ea)~x))$coef, error=function(e) NULL)
    if (is.null(s)) NA_real_ else s[1,"Pr(>|z|)"] })
  k <- which(!is.na(up) & up<0.05)
  if (length(k)>=5) {
    set.seed(1); cv6 <- tryCatch(cv.glmnet(Xa[,k,drop=FALSE], y, family="cox", alpha=1, nfolds=5), error=function(e) NULL)
    if (!is.null(cv6)) res[["UniCox+Lasso"]] <- list(
      pa=as.numeric(predict(cv6, Xa[,k,drop=FALSE], s=cv6$lambda.min, type="link")),
      pb=as.numeric(predict(cv6, Xb[,k,drop=FALSE], s=cv6$lambda.min, type="link")))
  }
  res
}

gs <- fromJSON(file.path(PROJ,"data/signatures/msigdb_genesets.json"))
gs <- gs[sapply(gs,length)>=10 & sapply(gs,length)<=500]
set.seed(20261015); pick <- sample(names(gs), min(N_GS, length(gs)))
say("基因集 ", length(pick), " 个")

HDR <- c("geneset","n_gene","sel","val_reported","true_test","inflation",
         "lasso_val","lasso_test","oracle_test","gap_train_val")
OUT <- file.path(PROJ,"results/tables",paste0(TAG,".tsv"))
write.table(as.data.frame(matrix(NA,0,length(HDR),dimnames=list(NULL,HDR))),OUT,
            sep="\t",quote=FALSE,row.names=FALSE)

Xtr <- G(tr); Xva <- G(va); Xte <- G(te)
ttr<-sv$t[tr]; etr<-sv$e[tr]; tva<-sv$t[va]; eva<-sv$e[va]; tte<-sv$t[te]; ete<-sv$e[te]

res <- list(); t0 <- Sys.time()
for (i in seq_along(pick)) {
  nm <- pick[i]; g <- intersect(gs[[nm]], rownames(Z))
  if (length(g) < 10) next
  a <- Xtr[,g,drop=FALSE]; b <- Xva[,g,drop=FALSE]; c2 <- Xte[,g,drop=FALSE]
  ms <- tryCatch(fit_models(a, ttr, etr, b), error=function(e) NULL)
  if (is.null(ms) || length(ms) < 3) next
  vci <- sapply(ms, function(x) ci(tva, eva, x$pb))          # 验证集（用于选模型）
  if (sum(!is.na(vci)) < 3) next
  sel <- names(which.max(vci))                                # 标准做法：按验证集选
  reported <- vci[[sel]]                                      # 报告的就是这个数
  # 真实性能：同一模型在未碰过的测试集上
  set.seed(1)
  # 用训练集重新拟合被选中的模型，再预测测试集（避免用验证集信息）
  ms_plain <- fit_models(a, ttr, etr, c2)
  true_test <- if (sel %in% names(ms_plain)) ci(tte, ete, ms_plain[[sel]]$pb) else NA_real_
  lv <- if (!is.null(ms_plain[["Lasso"]])) ci(tva, eva, ms_plain[["Lasso"]]$pb) else NA_real_
  lt <- if (!is.null(ms_plain[["Lasso"]])) ci(tte, ete, ms_plain[["Lasso"]]$pb) else NA_real_
  orac <- max(sapply(ms_plain, function(x) ci(tte, ete, x$pb)), na.rm=TRUE)
  res[[nm]] <- data.frame(geneset=nm, n_gene=length(g), sel=sel,
                          val_reported=reported, true_test=true_test,
                          inflation=reported-true_test,
                          lasso_val=lv, lasso_test=lt, oracle_test=orac,
                          gap_train_val=max(sapply(ms, function(x) ci(ttr,etr,x$pa)),na.rm=TRUE)-reported)
  if (i %% 10 == 0) {
    write.table(do.call(rbind,res), OUT, sep="\t", quote=FALSE, row.names=FALSE)
    say("  进度 ", i, "/", length(pick), "  用时 ",
        round(as.numeric(difftime(Sys.time(),t0,units="min")),1), " 分钟")
  }
}
R <- do.call(rbind, res)
write.table(R, OUT, sep="\t", quote=FALSE, row.names=FALSE)

say(""); say("========== v2 结果 (", COHORT, ") ==========")
say("完成 ", nrow(R), " 个基因集")
say(sprintf("报告值(验证集) 平均     = %.3f", mean(R$val_reported, na.rm=TRUE)))
say(sprintf("真实性能(测试集) 平均   = %.3f", mean(R$true_test, na.rm=TRUE)))
say(sprintf(">>> 模型选择偏倚        = %.3f", mean(R$inflation, na.rm=TRUE)))
say(sprintf("无选择基线(Lasso 测试)  = %.3f", mean(R$lasso_test, na.rm=TRUE)))
ok <- !is.na(R$inflation)
if (sum(ok)>=10) { tt <- t.test(R$inflation[ok])
  say(sprintf("配对 t 检验: 95%% CI [%.3f, %.3f], p = %.2e", tt$conf.int[1], tt$conf.int[2], tt$p.value)) }
close(con)
