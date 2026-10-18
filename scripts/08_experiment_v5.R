#!/usr/bin/env Rscript
# ============================================================
# 08 — 实验 v5：模型池大小 → 选择偏倚的剂量-反应
#
# 研究问题：真实论文用的是「101 组合」式的大模型池。"在众多候选里挑
#           验证集上最好的那个"会抬高报告值，而这个抬高的幅度是否随
#           候选数量增长？
#
# 设计（v5 终版，已通过 4 项代码审查修正）：
#   1) 三分法：训练 50% / 验证 25% / 测试 25%（按事件分层）
#   2) 模型池 19 个变体（线性靠 alpha×lambda 扩展；RSF/GBM 各 1 个）
#   3) 【同一批模型】同时预测验证集与测试集 —— 省一半算力，且让
#      "验证增益 vs 测试增益"的比较严格来自同一模型（v3 已验证
#      "训练拟合 vs 训练+验证拟合"的差异仅 0.004，可忽略）
#   4) 嵌套子池 k=2..K：bias(k) = [val增益] − [test增益]
#      val增益  = val(sel_k)  − val(ref)
#      test增益 = test(sel_k) − test(ref)      ref = Lasso(alpha=1)
#      → 集合难度差异在此相减中抵消，剩下的就是纯选择偏倚
#   5) 每 10 个基因集增量存盘；子池全 NA 时安全跳过
# ============================================================
suppressMessages({library(survival); library(glmnet); library(randomForestSRC); library(gbm); library(jsonlite)})
options(stringsAsFactors = FALSE)
`%||%` <- function(a, b) if (is.null(a)) b else a

PROJ <- "D:/bioinfo/sigAudit"
LOGDIR <- file.path(PROJ, "logs"); dir.create(LOGDIR, showWarnings = FALSE, recursive = TRUE)

ARGS <- commandArgs(trailingOnly = TRUE)
N_GS   <- if (length(ARGS) >= 1) as.integer(ARGS[1]) else 100
COHORT <- if (length(ARGS) >= 2) ARGS[2] else "LGG"
TAG    <- sprintf("08_%s_%d", COHORT, N_GS)
LOG    <- file.path(LOGDIR, paste0(TAG, ".log")); con <- file(LOG, "wt")
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
if (nrow(expr) < ncol(expr)) expr <- t(expr)          # -> 基因 × 样本
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

# ---------------- 三分法 ----------------
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
G <- function(idx) t(Z[, colnames(Z)[idx], drop=FALSE])     # -> 样本 × 基因

# 防御性 C-index：长度不符/值为 NA 一律返回 NA
ci <- function(t, e, s) {
  s <- tryCatch(as.numeric(s), error=function(x) NULL)
  if (is.null(s) || length(s) != length(t)) return(NA_real_)
  ok <- is.finite(s) & !is.na(t) & !is.na(e)
  if (sum(ok) < 15 || sum(e[ok]) < 5) return(NA_real_)
  tryCatch(as.numeric(summary(coxph(Surv(t[ok],e[ok]) ~ s[ok]))$concordance[1]),
           error = function(x) NA_real_)
}

# ---------------- 模型池（19 个变体，顺序=常用优先）----------------
POOL <- list(
  list(n="Lasso",        k="lin", a=1.0, lam="min"),   # 同时作为参照 ref
  list(n="Ridge",        k="lin", a=0.0, lam="min"),
  list(n="Lasso_1se",    k="lin", a=1.0, lam="1se"),
  list(n="Enet0.5",      k="lin", a=0.5, lam="min"),
  list(n="Enet0.5_1se",  k="lin", a=0.5, lam="1se"),
  list(n="Enet0.9",      k="lin", a=0.9, lam="min"),
  list(n="Enet0.7",      k="lin", a=0.7, lam="min"),
  list(n="Enet0.3",      k="lin", a=0.3, lam="min"),
  list(n="Enet0.1",      k="lin", a=0.1, lam="min"),
  list(n="Enet0.9_1se",  k="lin", a=0.9, lam="1se"),
  list(n="Enet0.7_1se",  k="lin", a=0.7, lam="1se"),
  list(n="Enet0.3_1se",  k="lin", a=0.3, lam="1se"),
  list(n="Enet0.1_1se",  k="lin", a=0.1, lam="1se"),
  list(n="Ridge_1se",    k="lin", a=0.0, lam="1se"),
  list(n="UniCox+Lasso", k="uni", lam="min"),
  list(n="RSF100",       k="rsf", ntree=100),
  list(n="GBM100",       k="gbm", nt=100),
  list(n="GBM250",       k="gbm", nt=250),
  list(n="RSF250",       k="rsf", ntree=250)
)
KMAX <- length(POOL)
POOL_NAMES <- vapply(POOL, function(m) m$n, character(1))
say("模型池 ", KMAX, " 个变体")

# 在 (Xa,ta,ea) 上一次性拟合全池，并对 Xb / Xc 两个集合分别预测
fit_pool <- function(Xa, ta, ea, Xb, Xc) {
  y <- Surv(ta, ea); out <- list()
  for (m in POOL) {
    nm <- m$n; s <- NULL
    if (m$k == "lin") {
      set.seed(1)
      cv <- tryCatch(cv.glmnet(Xa, y, family="cox", alpha=m$a, nfolds=5),
                     error=function(e) NULL)
      if (!is.null(cv)) {
        l <- if (m$lam == "min") cv$lambda.min else cv$lambda.1se
        s <- list(pb = as.numeric(predict(cv, Xb, s=l, type="link")),
                  pc = as.numeric(predict(cv, Xc, s=l, type="link")))
      }
    } else if (m$k == "uni") {
      up <- apply(Xa, 2, function(x){
        z <- tryCatch(summary(coxph(Surv(ta,ea) ~ x))$coef, error=function(e) NULL)
        if (is.null(z)) NA_real_ else z[1, "Pr(>|z|)"] })
      k <- which(!is.na(up) & up < 0.05)
      if (length(k) >= 5) {
        set.seed(1)
        cv6 <- tryCatch(cv.glmnet(Xa[,k,drop=FALSE], y, family="cox", alpha=1, nfolds=5),
                        error=function(e) NULL)
        if (!is.null(cv6)) s <- list(
          pb = as.numeric(predict(cv6, Xb[,k,drop=FALSE], s=cv6$lambda.min, type="link")),
          pc = as.numeric(predict(cv6, Xc[,k,drop=FALSE], s=cv6$lambda.min, type="link")))
      }
    } else if (m$k == "rsf") {
      # randomForestSRC 3.9.0：公式写 Surv(ta,ea) 会报 "Survival formula
      # incorrectly specified"，x=/y= 接口也已移除（y 会部分匹配 ytry）。
      # 正确做法：时间/状态作为两列放进 data.frame。
      d <- data.frame(.time = ta, .status = ea, Xa, check.names = FALSE)
      set.seed(1)
      rf <- tryCatch(rfsrc(Surv(.time, .status) ~ ., data = d, ntree = m$ntree,
                           nodesize = m$ns %||% 5, splitrule = "logrank"),
                     error = function(e) NULL)
      if (!is.null(rf)) s <- list(
        pb = as.numeric(predict(rf, data.frame(Xb, check.names = FALSE))$predicted),
        pc = as.numeric(predict(rf, data.frame(Xc, check.names = FALSE))$predicted))
    } else if (m$k == "gbm") {
      # gbm 的公式接口不可用；gbm.fit + predict 可用
      set.seed(1)
      gb <- tryCatch(gbm.fit(x = Xa, y = Surv(ta, ea), distribution = "coxph",
                             n.trees = m$nt, shrinkage = m$sh %||% 0.01,
                             interaction.depth = m$d %||% 3, n.minobsinnode = 10,
                             verbose = FALSE), error = function(e) NULL)
      if (!is.null(gb)) s <- list(
        pb = as.numeric(predict(gb, Xb, n.trees = m$nt, type = "link")),
        pc = as.numeric(predict(gb, Xc, n.trees = m$nt, type = "link")))
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

OUT <- file.path(PROJ,"results/tables", paste0(TAG,".tsv"))
res <- list(); t0 <- Sys.time(); n_skip <- 0

Xtr <- G(tr); Xva <- G(va); Xte <- G(te)
ttr<-sv$t[tr]; etr<-sv$e[tr]; tva<-sv$t[va]; eva<-sv$e[va]; tte<-sv$t[te]; ete<-sv$e[te]

for (i in seq_along(pick)) {
  nm <- pick[i]; g <- intersect(gs[[nm]], rownames(Z))
  if (length(g) < 10) { n_skip <- n_skip + 1; next }
  a <- Xtr[,g,drop=FALSE]; b <- Xva[,g,drop=FALSE]; c2 <- Xte[,g,drop=FALSE]
  # 同一批（训练集拟合的）模型，分别预测验证集与测试集
  pm <- tryCatch(fit_pool(a, ttr, etr, b, c2), error=function(e) NULL)
  if (is.null(pm) || length(pm) < 3) { n_skip <- n_skip + 1; next }
  avail <- POOL_NAMES[POOL_NAMES %in% names(pm)]
  if (length(avail) < 3 || !("Lasso" %in% avail)) { n_skip <- n_skip + 1; next }
  vv <- sapply(avail, function(x) ci(tva, eva, pm[[x]]$pb))
  tt <- sapply(avail, function(x) ci(tte, ete, pm[[x]]$pc))
  vref <- vv[["Lasso"]]; tref <- tt[["Lasso"]]
  if (is.na(vref) || is.na(tref)) { n_skip <- n_skip + 1; next }
  K <- length(avail)
  bias_k <- setNames(rep(NA_real_, KMAX), paste0("bias_k", 1:KMAX))
  for (k in 2:K) {
    sub <- avail[1:k]
    vs <- vv[sub]
    if (all(is.na(vs))) next                       # 子池全 NA 时安全跳过
    vsel <- sub[which.max(vs)]
    if (is.na(vv[[vsel]]) || is.na(tt[[vsel]])) next
    bias_k[[paste0("bias_k", k)]] <- (vv[[vsel]] - vref) - (tt[[vsel]] - tref)
  }
  vs_all <- vv[avail]
  if (all(is.na(vs_all))) { n_skip <- n_skip + 1; next }
  sel_max <- avail[which.max(vs_all)]
  row <- data.frame(geneset = nm, n_gene = length(g), K_avail = K, stringsAsFactors = FALSE)
  for (k in 2:KMAX) row[[paste0("bias_k", k)]] <- bias_k[[paste0("bias_k", k)]]
  row$val_prem_maxK  <- vv[[sel_max]] - vref
  row$test_prem_maxK <- tt[[sel_max]] - tref
  row$bias_maxK      <- row$val_prem_maxK - row$test_prem_maxK
  row$n_models_val   <- sum(!is.na(vv[avail]))
  res[[nm]] <- row
  if (i %% 10 == 0) {
    write.table(do.call(rbind, res), OUT, sep="\t", quote=FALSE, row.names=FALSE)
    say("  进度 ", i, "/", length(pick), "  已完成 ", length(res),
        "  跳过 ", n_skip, "  用时 ",
        round(as.numeric(difftime(Sys.time(), t0, units="min")),1), " 分钟")
  }
}
if (length(res)) {
  R <- do.call(rbind, res)
  write.table(R, OUT, sep="\t", quote=FALSE, row.names=FALSE)
  say(""); say("========== v5 剂量-反应结果 (", COHORT, ") ==========")
  say("完成 ", nrow(R), " 个基因集, 跳过 ", n_skip)
  say(sprintf("最大池 (K=%d) 的纯选择偏倚 = %+.4f", KMAX, mean(R$bias_maxK, na.rm=TRUE)))
  say("各子池 k 的平均选择偏倚:")
  for (k in c(2, 3, 5, 8, 12, 16, KMAX)) {
    v <- R[[paste0("bias_k", k)]]
    if (is.null(v)) next
    n <- sum(!is.na(v)); if (n < 10) next
    say(sprintf("  k=%2d  偏倚 = %+.4f  (n=%d)", k, mean(v, na.rm=TRUE), n))
  }
  ok <- !is.na(R$bias_maxK)
  if (sum(ok) >= 10) {
    t0t <- t.test(R$bias_maxK[ok])
    say(sprintf("最大池偏倚 95%% CI [%+.4f, %+.4f], p = %.2e",
                t0t$conf.int[1], t0t$conf.int[2], t0t$p.value))
  }
} else say("没有完成的基因集！")
close(con)
