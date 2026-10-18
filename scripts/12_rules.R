#!/usr/bin/env Rscript
# ============================================================
# 12 — 决策规则的收益核算：多算法挑选到底买到了什么？
#
# 之前所有分析都在量【偏倚】（报告值 − 真实值）。但偏倚小不等于好：
# 一个"挑最差模型"的规则偏倚也为 0。评价规则的唯一诚实指标是
# 【绝对真实性能】——在完全未参与选择的测试集上的 C-index。
#
# 于是把各种"选模型规则"放在同一把尺子上（均为零额外算力，全部离线）：
#   none   不做任何选择（参照，偏倚恒为 0）
#   argmax_val   标准做法：在全部候选里挑验证集最优        <- 文献主流
#   runner_up    挑验证集次优（规避 winner's curse 的朴素对策）
#   linear_only  限制在线性家族里挑（领域俗语"别用复杂模型"）
#   tree_only    限制在树集成家族里挑（反向对照）
#   oracle_test  用测试集挑（不可实现，仅作性能上界参照）
#
# 对每条规则报告：平均真实（测试集）C-index、平均报告值（验证集）、偏倚，
# 以及与 none 的【配对】真实性能差（这是"选模型是否真有收益"的答案）。
# ============================================================
suppressMessages(source("D:/bioinfo/sigAudit/R/audit_modules.R"))
options(stringsAsFactors = FALSE)
PROJ <- "D:/bioinfo/sigAudit"

POOL_NAMES <- c("Lasso","Ridge","Enet0.9","Enet0.7","Enet0.5","Enet0.3","Enet0.1",
                "Ridge_1se","UniCox+Lasso","LassoTop100","LassoTop50","LassoTop25",
                "LassoTop10","Enet0.5Top25","RSF100","GBM100","GBM250","RSF250")
TREES <- c("RSF100","GBM100","GBM250","RSF250")

LOGF <- file.path(PROJ,"logs","12_rules.log"); con <- file(LOGF,"wt")
say <- function(...) { m <- paste0(...); cat(m,"\n",sep=""); cat(m,"\n",sep="",file=con); flush(con) }

files <- list.files(file.path(PROJ,"results/tables"), pattern="^09_.*_[0-9]+\\.tsv$", full.names=TRUE)
ALL <- list(); SUM <- list()

for (f in files) {
  coh <- sub("^09_(.+)_[0-9]+\\.tsv$","\\1", basename(f))
  D <- tryCatch(read.table(f, header=TRUE, sep="\t", check.names=FALSE), error=function(e) NULL)
  if (is.null(D) || !nrow(D)) next
  D <- D[D$model %in% POOL_NAMES, ]; if (!nrow(D)) next
  all_gs <- unique(D$geneset); n_all <- length(all_gs)
  cnt <- tapply(D$geneset, D$model, function(g) length(unique(g)))
  pnames <- POOL_NAMES[POOL_NAMES %in% names(cnt)[cnt/n_all >= 0.90]]
  if (!("Lasso" %in% pnames) || length(pnames) < 5) { say("[skip] ", coh, " 候选不足"); next }
  sub <- D[D$model %in% pnames, ]
  tab <- table(sub$geneset); keep <- names(tab)[tab == length(pnames)]
  if (length(keep) < 30) { say("[skip] ", coh, " 完整基因集 <30（可能增量写入中）"); next }
  mk <- function(v) { M <- tapply(sub[[v]], list(sub$geneset, sub$model), function(x) x[1])
    M <- M[keep, pnames, drop=FALSE]
    matrix(as.numeric(M), nrow=nrow(M), dimnames=dimnames(M)) }
  V <- mk("val_ci"); Tt <- mk("test_ci")
  ok <- complete.cases(V) & complete.cases(Tt)
  V <- V[ok,,drop=FALSE]; Tt <- Tt[ok,,drop=FALSE]
  if (nrow(V) < 30) next
  n <- nrow(V)
  LIN <- setdiff(pnames, TREES)

  # 每条规则 -> 每个基因集选中的候选名
  pick_argmax <- function(cols) apply(V[, cols, drop=FALSE], 1,
                                      function(z) cols[which.max(z)])
  pick_rank <- function(cols, r) t(apply(V[, cols, drop=FALSE], 1, function(z) {
      o <- cols[order(-z, seq_along(z))]; o[min(r, length(o))] }))
  RULES <- list(
    none  = rep("Lasso", n),
    argmax_val  = pick_argmax(pnames),
    runner_up   = pick_rank(pnames, 2),
    linear_only = pick_argmax(LIN),
    tree_only   = if (length(intersect(pnames,TREES)) >= 2) pick_argmax(intersect(pnames,TREES)) else NULL,
    oracle_test = apply(Tt, 1, function(z) pnames[which.max(z)])
  )

  rows <- do.call(rbind, lapply(names(RULES), function(rn) {
    sel <- RULES[[rn]]; if (is.null(sel)) return(NULL)
    i <- seq_len(n)
    # cbind(整数, 字符) 会生成字符矩阵，按名索引会越界 —— 必须转成列号
    j <- match(sel, colnames(V))
    if (anyNA(j)) return(NULL)
    rep_v <- V[cbind(i, j)]; tru <- Tt[cbind(i, j)]
    base  <- Tt[, "Lasso"]
    d <- tru - base; ok2 <- is.finite(d)
    q <- t.test(d[ok2]); w <- suppressWarnings(wilcox.test(d[ok2]))
    data.frame(cohort=coh, rule=rn, n=n,
               reported=mean(rep_v, na.rm=TRUE), true_test=mean(tru, na.rm=TRUE),
               bias=mean(rep_v - tru, na.rm=TRUE),
               gain_vs_ref=mean(d[ok2]), ci_lo=q$conf.int[1], ci_hi=q$conf.int[2],
               t_p=q$p.value, wilcox_p=w$p.value)
  }))
  SUM[[coh]] <- rows
  ALL[[coh]] <- data.frame(cohort=coh, geneset=rownames(V), rule=NA)
  say(""); say("===== ", coh, " : ", n, " 基因集 =====")
  for (i in seq_len(nrow(rows))) say(sprintf(
    "  %-12s 报告值 %.3f  真实值 %.3f  偏倚 %+.4f   相对 Lasso 的真实增益 %+.4f [%+.4f,%+.4f]  p=%.2e",
    rows$rule[i], rows$reported[i], rows$true_test[i], rows$bias[i],
    rows$gain_vs_ref[i], rows$ci_lo[i], rows$ci_hi[i], rows$t_p[i]))
}

if (length(SUM)) {
  S <- do.call(rbind, SUM)
  write.table(S, file.path(PROJ,"results/tables","12_rules.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)
  say(""); say("========== 跨队列汇总（各规则的真实收益均值）==========")
  ag <- aggregate(cbind(reported,true_test,bias,gain_vs_ref) ~ rule, data=S, FUN=mean)
  ag <- ag[order(-ag$true_test), ]
  for (i in seq_len(nrow(ag))) say(sprintf(
    "  %-12s 报告值 %.3f  真实值 %.3f  偏倚 %+.4f  相对 Lasso 增益 %+.4f",
    ag$rule[i], ag$reported[i], ag$true_test[i], ag$bias[i], ag$gain_vs_ref[i]))
  say(""); say("  按队列看 argmax_val 的真实收益（选模型到底买到了什么）：")
  s2 <- S[S$rule=="argmax_val", ]
  for (i in seq_len(nrow(s2))) say(sprintf(
    "    %-6s 真实增益 %+.4f [%+.4f,%+.4f]  p=%.2e   %s",
    s2$cohort[i], s2$gain_vs_ref[i], s2$ci_lo[i], s2$ci_hi[i], s2$t_p[i],
    if (s2$ci_lo[i] > 0) "显著为正" else if (s2$ci_hi[i] < 0) "显著为负" else "与不做选择无差异"))
}
close(con)
