#!/usr/bin/env Rscript
# ============================================================
# 09 — 实验 v6：模型级明细输出（供 10 做任意子集 / 任意顺序的离线分析）
#
# 为什么需要它：
#   08（v5）只保存了「前缀子池」的偏倚曲线（bias_k2..bias_k19），它能
#   说明 "池子越大偏倚越大，且在 k=16 处跳升"，但无法区分两个互相竞争
#   的解释：
#     H1 剂量假说：偏倚来自池子【变大】本身
#     H2 成分假说：偏倚来自池子里【混进了树集成类算法】
#   若 H2 成立，把 RSF/GBM 挪到池子第 2~3 位，跳升应当【提前】出现；
#   若 H1 成立，跳升位置应当【跟随 k，而不跟随算法类型】。
#   判别这两者需要每个模型在验证/测试集上的 C-index 明细。
#   拿到明细后，任何子集、任何排列都是离线零成本计算，不必重跑。
#
#   v5 还有一个致命缺陷已在本版修正（详见下方 POOL 注释）：lambda.1se 在
#   稀疏惩罚路径上退化成空模型，使 6/19 个候选"死"掉。v6 的池子保证每个
#   候选都有活性且互相可区分。
#
# 与 08 的一致性：同一份数据、同一种标准化、同一个划分种子
# （set.seed(20261014)）、同一套分位点。模型池已升级为 v6。
#
# 输出：results/tables/09_<COHORT>_<N>.tsv  （长表）
#   geneset, n_gene, model, val_ci, test_ci
# ============================================================
suppressMessages({library(survival); library(glmnet); library(randomForestSRC); library(gbm); library(jsonlite)})
options(stringsAsFactors = FALSE)
`%||%` <- function(a, b) if (is.null(a)) b else a

PROJ <- "D:/bioinfo/sigAudit"
LOGDIR <- file.path(PROJ, "logs"); dir.create(LOGDIR, showWarnings = FALSE, recursive = TRUE)

ARGS <- commandArgs(trailingOnly = TRUE)
N_GS   <- if (length(ARGS) >= 1) as.integer(ARGS[1]) else 100
COHORT <- if (length(ARGS) >= 2) ARGS[2] else "LGG"
# 第 3 个参数：划分随机种子。所有主结果用 20261014；换一个种子 = 一次
# 独立的划分复现，用于检验"每个候选的进入效应排名"是否稳定（审稿人必问）。
# 第 4 个参数：输出文件名后缀，避免覆盖主结果（主结果空串）。
SEED_SPLIT <- if (length(ARGS) >= 3) as.integer(ARGS[3]) else 20261014L
SFX        <- if (length(ARGS) >= 4) ARGS[4] else ""
TAG    <- sprintf("09%s_%s_%d", SFX, COHORT, N_GS)
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

# ---------------- 三分法（与 08 完全同一划分）----------------
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
G <- function(idx) t(Z[, colnames(Z)[idx], drop=FALSE])     # -> 样本 × 基因

# 防御性 C-index（与 08 同一实现）
ci <- function(t, e, s) {
  s <- tryCatch(as.numeric(s), error=function(x) NULL)
  if (is.null(s) || length(s) != length(t)) return(NA_real_)
  ok <- is.finite(s) & !is.na(t) & !is.na(e)
  if (sum(ok) < 15 || sum(e[ok]) < 5) return(NA_real_)
  tryCatch(as.numeric(summary(coxph(Surv(t[ok],e[ok]) ~ s[ok]))$concordance[1]),
           error = function(x) NA_real_)
}

# ---------------- 模型池 v6：18 个【全部有活性且互相可区分】的候选 ----------------
# 为什么改（v5 的致命缺陷）：
#   诊断（scripts/_diag_1se.R）显示 lambda.1se 在稀疏惩罚路径上 100% 退化成
#   系数全零的空模型（每折仅约 19 个事件，1 个 SE 裕度足以跨过整条非零路径），
#   于是 Lasso_1se / Enet*_1se 六个变体恒等于 C-index 0.500，永远不会被选中。
#   v5 的 "k=2..15 剂量增长" 因此有 6/15 是死格子 —— 曲线形状被污染。
#   （这也解释了 v5 里 bias_k2 与 bias_k3 分毫不差：第 3 位是死模型。）
#
# v6 的原则：k 轴上每一格都必须是一个【真实、有活性、且与其他格可区分】的候选。
#   - 删掉全部 1se 退化解（Ridge_1se 除外：Ridge 不稀疏，实测 0% 退化）
#   - 线性家族用 alpha 阶梯 + 「单变量预筛 top-K 后再惩罚」两条正交轴扩展
#   - 树集成家族全部置于末尾，使 "线性剂量段" 与 "成分进入点" 清晰可分
POOL <- list(
  # --- 线性/正则化家族（第 1..14 位）---
  list(n="Lasso",        k="lin", a=1.0, lam="min"),   # 参照 ref
  list(n="Ridge",        k="lin", a=0.0, lam="min"),
  list(n="Enet0.9",      k="lin", a=0.9, lam="min"),
  list(n="Enet0.7",      k="lin", a=0.7, lam="min"),
  list(n="Enet0.5",      k="lin", a=0.5, lam="min"),
  list(n="Enet0.3",      k="lin", a=0.3, lam="min"),
  list(n="Enet0.1",      k="lin", a=0.1, lam="min"),
  list(n="Ridge_1se",    k="lin", a=0.0, lam="1se"),   # 非稀疏，无退化
  list(n="UniCox+Lasso", k="uni", lam="min"),          # 单变量 p<0.05 预筛
  list(n="LassoTop100",  k="lin", a=1.0, lam="min", pre=100),
  list(n="LassoTop50",   k="lin", a=1.0, lam="min", pre=50),
  list(n="LassoTop25",   k="lin", a=1.0, lam="min", pre=25),
  list(n="LassoTop10",   k="lin", a=1.0, lam="min", pre=10),
  list(n="Enet0.5Top25", k="lin", a=0.5, lam="min", pre=25),
  # --- 树集成家族（第 15..18 位）---
  list(n="RSF100",       k="rsf", ntree=100),
  list(n="GBM100",       k="gbm", nt=100),
  list(n="GBM250",       k="gbm", nt=250),
  list(n="RSF250",       k="rsf", ntree=250)
)
POOL_FAMILY <- vapply(POOL, function(m) if (m$k %in% c("rsf","gbm")) "tree" else "linear", character(1))
K_LIN <- sum(POOL_FAMILY == "linear")   # 线性家族占据的前缀长度
KMAX <- length(POOL)
POOL_NAMES <- vapply(POOL, function(m) m$n, character(1))
say("模型池 ", KMAX, " 个变体")

fit_pool <- function(Xa, ta, ea, Xb, Xc) {
  y <- Surv(ta, ea); out <- list()
  for (m in POOL) {
    nm <- m$n; s <- NULL
    # 每个模型各自持有一份输入副本：预筛会按基因取子集，
    # 若在这里直接改写 Xa/Xb/Xc 会泄漏到后面的模型（变异 bug）。
    Aa <- Xa; Ab <- Xb; Ac <- Xc
    if (!is.null(m$pre)) {                       # 单变量预筛 top-K（只用训练集排名）
      zz <- apply(Aa, 2, function(x){
        q <- tryCatch(summary(coxph(Surv(ta,ea) ~ x))$coef, error=function(e) NULL)
        if (is.null(q)) NA_real_ else abs(q[1, "z"]) })
      zz[!is.finite(zz)] <- -Inf
      kk <- min(m$pre, ncol(Aa))
      if (kk < 5 || sum(zz > -Inf) < 5) { out[[nm]] <- NULL; next }
      ix <- order(zz, decreasing = TRUE)[seq_len(kk)]
      Aa <- Aa[, ix, drop=FALSE]; Ab <- Ab[, ix, drop=FALSE]; Ac <- Ac[, ix, drop=FALSE]
    }
    if (m$k == "lin") {
      set.seed(1)
      cv <- tryCatch(cv.glmnet(Aa, y, family="cox", alpha=m$a, nfolds=5),
                     error=function(e) NULL)
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
        cv6 <- tryCatch(cv.glmnet(Aa[,k,drop=FALSE], y, family="cox", alpha=1, nfolds=5),
                        error=function(e) NULL)
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

OUT <- file.path(PROJ,"results/tables", paste0(TAG,".tsv"))
res <- list(); t0 <- Sys.time(); n_skip <- 0

Xtr <- G(tr); Xva <- G(va); Xte <- G(te)
ttr<-sv$t[tr]; etr<-sv$e[tr]; tva<-sv$t[va]; eva<-sv$e[va]; tte<-sv$t[te]; ete<-sv$e[te]

for (i in seq_along(pick)) {
  nm <- pick[i]; g <- intersect(gs[[nm]], rownames(Z))
  if (length(g) < 10) { n_skip <- n_skip + 1; next }
  a <- Xtr[,g,drop=FALSE]; b <- Xva[,g,drop=FALSE]; c2 <- Xte[,g,drop=FALSE]
  pm <- tryCatch(fit_pool(a, ttr, etr, b, c2), error=function(e) NULL)
  if (is.null(pm) || length(pm) < 3) { n_skip <- n_skip + 1; next }
  avail <- POOL_NAMES[POOL_NAMES %in% names(pm)]
  if (!("Lasso" %in% avail)) { n_skip <- n_skip + 1; next }
  vv <- sapply(avail, function(x) ci(tva, eva, pm[[x]]$pb))
  tt <- sapply(avail, function(x) ci(tte, ete, pm[[x]]$pc))
  if (all(is.na(vv)) || all(is.na(tt))) { n_skip <- n_skip + 1; next }
  res[[nm]] <- data.frame(geneset = nm, n_gene = length(g), model = avail,
                          val_ci = as.numeric(vv), test_ci = as.numeric(tt),
                          stringsAsFactors = FALSE)
  if (i %% 10 == 0) {
    write.table(do.call(rbind, res), OUT, sep="\t", quote=FALSE, row.names=FALSE)
    say("  进度 ", i, "/", length(pick), "  已完成 ", length(res),
        "  跳过 ", n_skip, "  用时 ",
        round(as.numeric(difftime(Sys.time(), t0, units="min")),1), " 分钟")
    gc(verbose = FALSE)
  }
}
if (length(res)) {
  R <- do.call(rbind, res)
  write.table(R, OUT, sep="\t", quote=FALSE, row.names=FALSE)
  say(""); say("========== 09 模型级明细 (", COHORT, ") ==========")
  say("完成 ", length(res), " 个基因集, 跳过 ", n_skip, "  明细行 ", nrow(R))
  # 自检：参照模型必须每行都在
  say("含 Lasso 参照的基因集: ", length(unique(R$geneset[R$model=="Lasso"])))
  say("各模型平均 验证/测试 C-index、相对 Lasso 的增益、以及退化率:")
  say("  (退化 = val C-index 恰为 0.500，即模型输出常数风险分)")
  refmap_v <- setNames(R$val_ci [R$model=="Lasso"], R$geneset[R$model=="Lasso"])
  refmap_t <- setNames(R$test_ci[R$model=="Lasso"], R$geneset[R$model=="Lasso"])
  for (m in POOL_NAMES) {
    s <- R[R$model == m, ]
    if (!nrow(s)) { say(sprintf("  %-14s  ★缺失★", m)); next }
    rv <- refmap_v[s$geneset]; rt <- refmap_t[s$geneset]
    deg <- 100 * mean(abs(s$val_ci - 0.5) < 1e-9, na.rm = TRUE)
    say(sprintf("  %-14s  val %.3f  test %.3f  Δval %+.4f  Δtest %+.4f  退化 %3.0f%%  (n=%d)",
                m, mean(s$val_ci, na.rm=TRUE), mean(s$test_ci, na.rm=TRUE),
                mean(s$val_ci - rv, na.rm=TRUE), mean(s$test_ci - rt, na.rm=TRUE), deg, nrow(s)))
  }
  # 家族内"可区分度"：每个基因集上家族成员 val C-index 的标准差，再取均值。
  # 若该值接近 0，说明该家族成员彼此雷同 —— 增加成员数量不会改变选择结果，
  # 曲线就会是平的（这正是 v5 里 k=2..15 平坦的机制）。
  Wide <- do.call(rbind, lapply(split(R, R$geneset), function(d) {
    v <- setNames(d$val_ci, d$model)
    c(lin = sd(v[names(v) %in% POOL_NAMES[POOL_FAMILY=="linear"]], na.rm=TRUE),
      tree= sd(v[names(v) %in% POOL_NAMES[POOL_FAMILY=="tree"]],   na.rm=TRUE),
      all = sd(v, na.rm=TRUE)) }))
  say(sprintf("家族内可区分度（成员间 val C-index 的 SD，均值）: 线性 %.4f  树集成 %.4f  全体 %.4f",
      mean(Wide[,"lin"], na.rm=TRUE), mean(Wide[,"tree"], na.rm=TRUE), mean(Wide[,"all"], na.rm=TRUE)))
  say(sprintf("线性家族成员数 %d，树集成成员数 %d，K_MAX = %d",
      K_LIN, KMAX - K_LIN, KMAX))
} else say("没有完成的基因集！")
close(con)
