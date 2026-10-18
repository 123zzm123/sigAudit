#!/usr/bin/env Rscript
# ============================================================
# 07 — 实验 v4：用「增益之差」测量【纯】模型选择偏倚
#
# v2 的三处问题（已在 05 的结果中暴露）：
#   P1  true_test 只用【训练集】重拟合 -> 真实性能被低估 -> 虚高被高估。
#       文献标准做法是用【训练+验证】重拟合最终模型。
#   P2  没有量化 P1 这个口径差异本身有多大。
#   P3  缺少"为什么某些队列虚高、某些不虚高"的解释变量。
#
# v3 设计：
#   三分法 训练50 / 验证25 / 测试25（按事件分层）
#   标准做法：用【训练+验证】挑模型，报告【验证集】C-index   <- 论文里的数
#   真实性能A：同一模型用【训练+验证】重拟合 -> 测【测试集】   <- 文献口径
#   真实性能B：同一模型用【训练】重拟合    -> 测【测试集】   <- v2 口径, 供对比
#   inflation = 报告值 − 真实性能A
#
#   解释变量（用于回答"何时虚高"）：
#     disp_val  = 各模型验证集 C-index 的标准差   （模型间可分辨程度）
#     gain_sel  = 选中模型真实A − 固定Lasso真实A  （挑模型带来的真实收益）
#     n_events_test / n_events_val
# ============================================================
suppressMessages({library(survival); library(glmnet); library(randomForestSRC); library(gbm); library(jsonlite)})
options(stringsAsFactors = FALSE)
PROJ <- "D:/bioinfo/sigAudit"
LOGDIR <- file.path(PROJ, "logs"); dir.create(LOGDIR, showWarnings = FALSE, recursive = TRUE)

ARGS <- commandArgs(trailingOnly = TRUE)
N_GS   <- if (length(ARGS) >= 1) as.integer(ARGS[1]) else 100
COHORT <- if (length(ARGS) >= 2) ARGS[2] else "LGG"
TAG    <- sprintf("07_%s_%d", COHORT, N_GS)
LOG    <- file.path(LOGDIR, paste0(TAG, ".log")); con <- file(LOG, "wt")
say <- function(...) { m <- paste0(...); cat(m,"\n",sep=""); cat(m,"\n",sep="",file=con); flush(con) }

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

set.seed(20261014)
ev <- which(sv$e==1); ce <- which(sv$e==0)
tr <- sort(c(sample(ev, floor(.5*length(ev))), sample(ce, floor(.5*length(ce)))))
rest <- setdiff(seq_along(sv$e), tr)
re <- rest[sv$e[rest]==1]; rc <- rest[sv$e[rest]==0]
va <- sort(c(sample(re, floor(.5*length(re))), sample(rc, floor(.5*length(rc)))))
te <- setdiff(rest, va)
say("  训练 ", length(tr), " / 验证 ", length(va), " / 测试 ", length(te),
    "  事件 ", sum(sv$e[tr]), "/", sum(sv$e[va]), "/", sum(sv$e[te]))

Z <- t(scale(t(expr))); Z[!is.finite(Z)] <- 0
G <- function(idx) t(Z[, colnames(Z)[idx]])

ci <- function(t, e, s) {
  s <- tryCatch(as.numeric(s), error = function(x) NULL)
  if (is.null(s) || length(s) != length(t)) return(NA_real_)
  ok <- is.finite(s) & !is.na(t) & !is.na(e)
  if (sum(ok) < 15 || sum(e[ok]) < 5) return(NA_real_)
  tryCatch(as.numeric(summary(coxph(Surv(t[ok], e[ok]) ~ s[ok]))$concordance[1]),
           error = function(x) NA_real_)
}

NTREE <- 150
# 在 (Xa,ta,ea) 上拟合模型池，并对 Xb 预测
fit_pool <- function(Xa, ta, ea, Xb) {
  y <- Surv(ta, ea); res <- list()
  for (a in c(1, 0.5, 0)) {
    nm <- c("Lasso","Enet","Ridge")[match(a, c(1,0.5,0))]
    set.seed(1); cv <- tryCatch(cv.glmnet(Xa, y, family="cox", alpha=a, nfolds=5), error=function(e) NULL)
    if (!is.null(cv)) res[[nm]] <- as.numeric(predict(cv, Xb, s=cv$lambda.min, type="link"))
  }
  set.seed(1); rf <- tryCatch(rfsrc(Surv(ta,ea) ~ ., data=data.frame(Surv(ta,ea), Xa),
      ntree=NTREE, nodesize=5, splitrule="logrank"), error=function(e) NULL)
  if (!is.null(rf)) res[["RSF"]] <- as.numeric(predict(rf, data.frame(Xb))$predicted)
  set.seed(1); gb <- tryCatch(gbm(Surv(ta,ea) ~ ., data=data.frame(Surv(ta,ea), Xa),
      distribution="coxph", n.trees=NTREE, shrinkage=0.01, interaction.depth=3,
      n.minobsinnode=10, verbose=FALSE), error=function(e) NULL)
  if (!is.null(gb)) res[["GBM"]] <- as.numeric(predict(gb, data.frame(Xb), n.trees=NTREE, type="link"))
  up <- apply(Xa, 2, function(x){ s<-tryCatch(summary(coxph(Surv(ta,ea)~x))$coef, error=function(e) NULL)
    if (is.null(s)) NA_real_ else s[1,"Pr(>|z|)"] })
  k <- which(!is.na(up) & up<0.05)
  if (length(k)>=5) {
    set.seed(1); cv6 <- tryCatch(cv.glmnet(Xa[,k,drop=FALSE], y, family="cox", alpha=1, nfolds=5), error=function(e) NULL)
    if (!is.null(cv6)) res[["UniCox+Lasso"]] <- as.numeric(predict(cv6, Xb[,k,drop=FALSE], s=cv6$lambda.min, type="link"))
  }
  res
}

gs <- fromJSON(file.path(PROJ,"data/signatures/msigdb_genesets.json"))
gs <- gs[sapply(gs,length)>=10 & sapply(gs,length)<=500]
set.seed(20261015); pick <- sample(names(gs), min(N_GS, length(gs)))
say("基因集 ", length(pick), " 个")

HDR <- c("geneset","n_gene","sel","val_reported","val_lasso","val_premium",
         "true_tv","true_t","test_lasso","test_premium","sel_bias",
         "inflation_raw","disp_val","oracle_test")
OUT <- file.path(PROJ,"results/tables",paste0(TAG,".tsv"))
write.table(as.data.frame(matrix(NA,0,length(HDR),dimnames=list(NULL,HDR))),OUT,
            sep="\t",quote=FALSE,row.names=FALSE)

Xtr <- G(tr); Xva <- G(va); Xte <- G(te)
ttr<-sv$t[tr]; etr<-sv$e[tr]; tva<-sv$t[va]; eva<-sv$e[va]; tte<-sv$t[te]; ete<-sv$e[te]
tv <- c(ttr, tva); ev <- c(etr, eva); Xtv <- rbind(Xtr, Xva)      # 训练+验证

res <- list(); t0 <- Sys.time()
for (i in seq_along(pick)) {
  nm <- pick[i]; g <- intersect(gs[[nm]], rownames(Z))
  if (length(g) < 10) next
  a <- Xtr[,g,drop=FALSE]; b <- Xva[,g,drop=FALSE]; c2 <- Xte[,g,drop=FALSE]; ab <- Xtv[,g,drop=FALSE]
  pv <- tryCatch(fit_pool(a, ttr, etr, b), error=function(e) NULL)     # 训练拟合 -> 预测验证
  if (is.null(pv) || length(pv) < 3) next
  vci <- sapply(pv, function(s) ci(tva, eva, s))
  if (sum(!is.na(vci)) < 3) next
  sel <- names(which.max(vci)); reported <- vci[[sel]]
  disp <- sd(vci, na.rm = TRUE)                                       # 模型间可分辨程度
  # 真实性能：用 训练+验证 重拟合（文献口径） vs 只用训练（v2 口径）
  ptv <- tryCatch(fit_pool(ab, tv, ev, c2), error=function(e) NULL)
  pt  <- tryCatch(fit_pool(a,  ttr, etr, c2), error=function(e) NULL)
  true_tv <- if (!is.null(ptv[[sel]])) ci(tte, ete, ptv[[sel]]) else NA_real_
  true_t  <- if (!is.null(pt[[sel]]))  ci(tte, ete, pt[[sel]])  else NA_real_
  # 固定基线 Lasso 在两个集合上的表现 -> 用"增益之差"抵消集合难度差异
  val_lasso  <- if (!is.null(pv[["Lasso"]]))  ci(tva, eva, pv[["Lasso"]])  else NA_real_
  test_lasso <- if (!is.null(ptv[["Lasso"]])) ci(tte, ete, ptv[["Lasso"]]) else NA_real_
  val_prem   <- reported - val_lasso          # 验证集上"选模型"相对固定基线的增益
  test_prem  <- true_tv  - test_lasso         # 测试集上的同一增益（真实值）
  sel_bias   <- val_prem - test_prem          # <<< 纯粹的选择偏倚，集合难度已抵消
  orac <- if (!is.null(ptv)) max(sapply(ptv, function(s) ci(tte, ete, s)), na.rm=TRUE) else NA_real_
  res[[nm]] <- data.frame(geneset=nm, n_gene=length(g), sel=sel, val_reported=reported,
                          val_lasso=val_lasso, val_premium=val_prem,
                          true_tv=true_tv, true_t=true_t, test_lasso=test_lasso,
                          test_premium=test_prem, sel_bias=sel_bias,
                          inflation_raw=reported-true_tv,
                          disp_val=disp, oracle_test=orac)
  if (i %% 10 == 0) {
    write.table(do.call(rbind,res), OUT, sep="\t", quote=FALSE, row.names=FALSE)
    say("  进度 ", i, "/", length(pick), "  用时 ",
        round(as.numeric(difftime(Sys.time(),t0,units="min")),1), " 分钟")
  }
}
R <- do.call(rbind, res)
write.table(R, OUT, sep="\t", quote=FALSE, row.names=FALSE)

say(""); say("========== v3 结果 (", COHORT, ") ==========")
say("完成 ", nrow(R))
say(sprintf("报告值(验证集)            = %.3f", mean(R$val_reported, na.rm=TRUE)))
say(sprintf("真实性能(训练+验证重拟合) = %.3f", mean(R$true_tv, na.rm=TRUE)))
say(sprintf("真实性能(仅训练重拟合)    = %.3f", mean(R$true_t, na.rm=TRUE)))
say(sprintf("原始口径偏倚(验证-测试)    = %.3f  (含集合难度混淆)", mean(R$inflation_raw, na.rm=TRUE)))
say(sprintf("验证集上选模型的增益       = %.3f", mean(R$val_premium, na.rm=TRUE)))
say(sprintf("测试集上同一增益(真实)     = %.3f", mean(R$test_premium, na.rm=TRUE)))
say(sprintf(">>> 纯选择偏倚(增益之差)   = %.3f", mean(R$sel_bias, na.rm=TRUE)))
say(sprintf("无选择基线(固定Lasso 测试) = %.3f", mean(R$test_lasso, na.rm=TRUE)))
ok <- !is.na(R$sel_bias)
if (sum(ok)>=10) { tt <- t.test(R$sel_bias[ok])
  say(sprintf("配对 t: 95%% CI [%.3f, %.3f], p = %.2e", tt$conf.int[1], tt$conf.int[2], tt$p.value)) }
close(con)
