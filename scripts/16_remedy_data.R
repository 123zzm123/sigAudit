#!/usr/bin/env Rscript
# ============================================================
# 16 — 对策实验（数据生成）：存下【逐样本预测】，让"减偏倚对策"可以被离线检验
#
# 动机：09 只存了每个 (签名 × 候选) 的一个 C-index 数值。但检验对策时，
#       很多规则需要在【同一批模型】上对验证集重新划分 / bootstrap ——
#       例如"重复二分验证集取最差表现"、"bootstrap 多数投票选模型"、
#       "把 top-k 候选的预测平均成集成"。这些都无法从单个 C-index 复原，
#       必须保留逐样本预测。
#
# 因此本脚本与 09 唯一的区别是：把 predict 出来的风险分原样存盘（RDS），
# 不做任何 C-index 计算。所有指标留给 17_remedies.R 离线算。
#
# 输出：results/remedies/16_<COHORT>_<N>.rds
#   一个 list，每个元素对应一个基因集：
#     name, models(候选名), val_pred(样本×候选), test_pred(样本×候选),
#     tva, eva, tte, ete, n_gene
# ============================================================
suppressMessages({library(survival); library(glmnet); library(randomForestSRC); library(gbm); library(jsonlite)})
options(stringsAsFactors = FALSE)
`%||%` <- function(a, b) if (is.null(a)) b else a

PROJ <- "D:/bioinfo/sigAudit"
LOGDIR <- file.path(PROJ, "logs"); dir.create(LOGDIR, showWarnings = FALSE, recursive = TRUE)
OUTDIR <- file.path(PROJ, "results/remedies"); dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

ARGS <- commandArgs(trailingOnly = TRUE)
N_GS   <- if (length(ARGS) >= 1) as.integer(ARGS[1]) else 100
COHORT <- if (length(ARGS) >= 2) ARGS[2] else "LGG"
SEED_SPLIT <- if (length(ARGS) >= 3) as.integer(ARGS[3]) else 20261014L
TAG <- sprintf("16_%s_%d", COHORT, N_GS)
LOG <- file.path(LOGDIR, paste0(TAG, ".log")); con <- file(LOG, "wt")
say <- function(...) { m <- paste0(...); cat(m,"\n",sep=""); cat(m,"\n",sep="",file=con); flush(con) }

# ---------------- 数据 ----------------
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
if (!ecol %in% names(sv)) stop("缺少事件列: ", ecol)
sv$t <- suppressWarnings(as.numeric(sv[[tcol]]))
sv$e <- ifelse(grepl("^1", sv[[ecol]]), 1L, ifelse(grepl("^0", sv[[ecol]]), 0L, NA_integer_))
sv <- sv[!is.na(sv$t) & !is.na(sv$e) & sv$t > 0, ]
colnames(expr) <- substr(colnames(expr),1,12)
expr <- expr[, !duplicated(colnames(expr)), drop=FALSE]
com <- intersect(colnames(expr), rownames(sv)); expr <- expr[,com,drop=FALSE]; sv <- sv[com,]
say("  ", ncol(expr), " 样本 × ", nrow(expr), " 基因, 事件 ", sum(sv$e), " (", tcol, ")")

# ---------------- 三分法（与 09 完全同一划分）----------------
set.seed(SEED_SPLIT)
say("  划分种子 ", SEED_SPLIT)
ev <- which(sv$e==1); ce <- which(sv$e==0)
tr <- sort(c(sample(ev, floor(.5*length(ev))), sample(ce, floor(.5*length(ce)))))
rest <- setdiff(seq_along(sv$e), tr)
re <- rest[sv$e[rest]==1]; rc <- rest[sv$e[rest]==0]
va <- sort(c(sample(re, floor(.5*length(re))), sample(rc, floor(.5*length(rc)))))
te <- setdiff(rest, va)
say("  训练 ", length(tr), " / 验证 ", length(va), " / 测试 ", length(te),
    "  事件 ", sum(sv$e[tr]), "/", sum(sv$e[va]), "/", sum(sv$e[te]))

Z <- t(scale(t(expr))); Z[!is.finite(Z)] <- 0
G <- function(idx) t(Z[, colnames(Z)[idx], drop=FALSE])

# ---------------- 模型池（与 09 逐字一致）----------------
POOL <- list(
  list(n="Lasso",        k="lin", a=1.0, lam="min"),
  list(n="Ridge",        k="lin", a=0.0, lam="min"),
  list(n="Enet0.9",      k="lin", a=0.9, lam="min"),
  list(n="Enet0.7",      k="lin", a=0.7, lam="min"),
  list(n="Enet0.5",      k="lin", a=0.5, lam="min"),
  list(n="Enet0.3",      k="lin", a=0.3, lam="min"),
  list(n="Enet0.1",      k="lin", a=0.1, lam="min"),
  list(n="Ridge_1se",    k="lin", a=0.0, lam="1se"),
  list(n="UniCox+Lasso", k="uni", lam="min"),
  list(n="LassoTop100",  k="lin", a=1.0, lam="min", pre=100),
  list(n="LassoTop50",   k="lin", a=1.0, lam="min", pre=50),
  list(n="LassoTop25",   k="lin", a=1.0, lam="min", pre=25),
  list(n="LassoTop10",   k="lin", a=1.0, lam="min", pre=10),
  list(n="Enet0.5Top25", k="lin", a=0.5, lam="min", pre=25),
  list(n="RSF100",       k="rsf", ntree=100),
  list(n="GBM100",       k="gbm", nt=100),
  list(n="GBM250",       k="gbm", nt=250),
  list(n="RSF250",       k="rsf", ntree=250)
)
POOL_NAMES <- vapply(POOL, function(m) m$n, character(1))
say("模型池 ", length(POOL), " 个变体")

# 与 09 的 fit_pool 相同，但返回逐样本预测（本来就有），不做 C-index
fit_pool <- function(Xa, ta, ea, Xb, Xc) {
  y <- Surv(ta, ea); out <- list()
  for (m in POOL) {
    nm <- m$n; s <- NULL
    Aa <- Xa; Ab <- Xb; Ac <- Xc
    if (!is.null(m$pre)) {
      zz <- apply(Aa, 2, function(x){
        q <- tryCatch(summary(coxph(Surv(ta,ea) ~ x))$coef, error=function(e) NULL)
        if (is.null(q)) NA_real_ else abs(q[1, "z"]) })
      zz[!is.finite(zz)] <- -Inf
      kk <- min(m$pre, ncol(Aa))
      if (kk < 5 || sum(zz > -Inf) < 5) next
      ix <- order(zz, decreasing = TRUE)[seq_len(kk)]
      Aa <- Aa[, ix, drop=FALSE]; Ab <- Ab[, ix, drop=FALSE]; Ac <- Ac[, ix, drop=FALSE]
    }
    if (m$k == "lin") {
      set.seed(1)
      cv <- tryCatch(cv.glmnet(Aa, y, family="cox", alpha=m$a, nfolds=5,
                               cox.ties = "breslow"), error=function(e) NULL)
      if (!is.null(cv)) {
        l <- if (m$lam == "min") cv$lambda.min else cv$lambda.1se
        s <- list(pb = as.numeric(predict(cv, Ab, s=l, type="link")),
                  pc = as.numeric(predict(cv, Ac, s=l, type="link")))
      }
    } else if (m$k == "uni") {
      up <- apply(Aa, 2, function(x){
        z <- tryCatch(summary(coxph(Surv(ta,ea) ~ x))$coef, error=function(e) NULL)
        if (is.null(z)) NA_real_ else z[1, "Pr(>|z|)"] })
      k <- which(!is.na(up) & up < 0.05)
      if (length(k) >= 5) {
        set.seed(1)
        cv6 <- tryCatch(cv.glmnet(Aa[,k,drop=FALSE], y, family="cox", alpha=1, nfolds=5,
                                  cox.ties = "breslow"), error=function(e) NULL)
        if (!is.null(cv6)) s <- list(
          pb = as.numeric(predict(cv6, Ab[,k,drop=FALSE], s=cv6$lambda.min, type="link")),
          pc = as.numeric(predict(cv6, Ac[,k,drop=FALSE], s=cv6$lambda.min, type="link")))
      }
    } else if (m$k == "rsf") {
      d <- data.frame(.time = ta, .status = ea, Aa, check.names = FALSE)
      set.seed(1)
      rf <- tryCatch(rfsrc(Surv(.time, .status) ~ ., data = d, ntree = m$ntree,
                           nodesize = m$ns %||% 5, splitrule = "logrank"),
                     error = function(e) NULL)
      if (!is.null(rf)) s <- list(
        pb = as.numeric(predict(rf, data.frame(Ab, check.names = FALSE))$predicted),
        pc = as.numeric(predict(rf, data.frame(Ac, check.names = FALSE))$predicted))
    } else if (m$k == "gbm") {
      set.seed(1)
      gb <- tryCatch(gbm.fit(x = Aa, y = Surv(ta, ea), distribution = "coxph",
                             n.trees = m$nt, shrinkage = m$sh %||% 0.01,
                             interaction.depth = m$d %||% 3, n.minobsinnode = 10,
                             verbose = FALSE), error = function(e) NULL)
      if (!is.null(gb)) s <- list(
        pb = as.numeric(predict(gb, Ab, n.trees = m$nt, type = "link")),
        pc = as.numeric(predict(gb, Ac, n.trees = m$nt, type = "link")))
    }
    if (!is.null(s)) out[[nm]] <- s
  }
  out
}

# ---------------- 基因集 ----------------
gs <- fromJSON(file.path(PROJ,"data/signatures/msigdb_genesets.json"))
gs <- gs[sapply(gs,length)>=10 & sapply(gs,length)<=500]
set.seed(20261015); pick <- sample(names(gs), min(N_GS, length(gs)))
say("基因集 ", length(pick), " 个")

OUT <- file.path(OUTDIR, paste0(TAG, ".rds"))
res <- list(); t0 <- Sys.time(); n_skip <- 0
Xtr <- G(tr); Xva <- G(va); Xte <- G(te)
ttr<-sv$t[tr]; etr<-sv$e[tr]; tva<-sv$t[va]; eva<-sv$e[va]; tte<-sv$t[te]; ete<-sv$e[te]

for (i in seq_along(pick)) {
  nm <- pick[i]; g <- intersect(gs[[nm]], rownames(Z))
  if (length(g) < 10) { n_skip <- n_skip + 1; next }
  pm <- tryCatch(fit_pool(Xtr[,g,drop=FALSE], ttr, etr,
                          Xva[,g,drop=FALSE], Xte[,g,drop=FALSE]), error=function(e) NULL)
  if (is.null(pm) || length(pm) < 3) { n_skip <- n_skip + 1; next }
  mods <- POOL_NAMES[POOL_NAMES %in% names(pm)]
  if (!("Lasso" %in% mods)) { n_skip <- n_skip + 1; next }
  P <- sapply(mods, function(x) pm[[x]]$pb)
  Q <- sapply(mods, function(x) pm[[x]]$pc)
  if (is.null(dim(P))) { P <- matrix(P, ncol=1, dimnames=list(NULL, mods))
                         Q <- matrix(Q, ncol=1, dimnames=list(NULL, mods)) }
  colnames(P) <- mods; colnames(Q) <- mods
  res[[nm]] <- list(name=nm, n_gene=length(g), models=mods,
                    val_pred=P, test_pred=Q, tva=tva, eva=eva, tte=tte, ete=ete)
  if (i %% 10 == 0) {
    saveRDS(res, OUT)
    say("  进度 ", i, "/", length(pick), "  已完成 ", length(res),
        "  跳过 ", n_skip, "  用时 ",
        round(as.numeric(difftime(Sys.time(), t0, units="min")),1), " 分钟")
    gc(verbose = FALSE)
  }
}
if (length(res)) {
  saveRDS(res, OUT)
  say(""); say("========== 16 逐样本预测 (", COHORT, ") ==========")
  say("完成 ", length(res), " 个基因集, 跳过 ", n_skip)
  say("验证集样本 ", length(tva), " 测试集样本 ", length(tte))
  say("含 Lasso 的基因集: ", sum(sapply(res, function(z) "Lasso" %in% z$models)))
  say("可用候选数分布:"); print(table(sapply(res, function(z) length(z$models))))
} else say("没有完成的基因集！")
close(con)
