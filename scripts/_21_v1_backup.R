#!/usr/bin/env Rscript
# ============================================================
# 21 — 进入效应与候选层面预测因子（修正版）
#
# 评审意见落地：
#  1. 旧 10_analysis.R 的"进入效应"用的是【固定 rest 顺序的前缀背景池】，
#     并非对该规模下所有子池平均，因此不得声称与位置/组成天然正交。
#     这里改为：对每个 (目标 t, 规模 k) 随机抽 B 个背景子池，
#     求边际变化的平均与分散程度 —— 称为【特定规模下的平均边际进入效应】。
#     不称 Shapley 值（未做任何公理化分解）。
#  2. 预测因子相关性：明确是【同一分析样本内的关联】，不是机制证明；
#     并补一个【留一癌种】评估，检查关联是否由单一癌种驱动。
#  3. 所有 C-index 用直接评价口径（与 20 相同），不使用旧写法。
#
# 输入：results/main/19_*.rds
# 输出：results/tables/21_entry_effect.tsv, 21_predictors.tsv, 21_loo.tsv
# ============================================================
suppressMessages(library(survival))
options(stringsAsFactors = FALSE)
PROJ <- "D:/bioinfo/sigAudit"; TAB <- file.path(PROJ, "results/tables")
TREES <- c("RSF100","GBM100","GBM250","RSF250")
B_BG <- 20                       # 每个 (目标, 规模) 的随机背景子池数
KS   <- c(2, 3, 5, 8, 12, 18)    # 评估的池规模（超出实际候选数则跳过）

LOGF <- file.path(PROJ,"logs","21_entry_predictors.log"); con <- file(LOGF,"wt")
say <- function(...) { m <- paste0(...); cat(m,"\n",sep=""); cat(m,"\n",sep="",file=con); flush(con) }

ci_new <- function(t, e, s) {
  s <- tryCatch(as.numeric(s), error=function(x) NULL)
  if (is.null(s) || length(s) != length(t)) return(NA_real_)
  ok <- is.finite(s) & !is.na(t) & !is.na(e)
  if (sum(ok) < 10) return(NA_real_)
  tryCatch(as.numeric(concordance(Surv(t[ok],e[ok]) ~ s[ok], reverse=TRUE, timewt="n")$concordance),
           error=function(x) NA_real_)
}

fs <- list.files(file.path(PROJ,"results/main"), pattern="^19_.*\\.rds$", full.names=TRUE)
if (!length(fs)) stop("没有 19_*.rds")
say("找到 ", length(fs), " 个文件；随机背景池数 B = ", B_BG)

ENT <- list(); FEA <- list()
for (f in fs) {
  R <- tryCatch(readRDS(f), error=function(e) NULL)
  if (is.null(R) || !length(R$genesets)) next
  coh <- R$cohort; seed <- R$seed; G <- R$genesets
  cnt <- table(unlist(lapply(G, function(z) colnames(z$val_pred))))
  keep_mods <- names(cnt)[cnt/length(G) >= 0.90]
  use <- Filter(function(z) all(keep_mods %in% colnames(z$val_pred)), G)
  if (length(use) < 30 || !("Lasso" %in% keep_mods)) { say("[skip] ", coh, " ", seed); next }
  K <- length(keep_mods)
  # ---- 预计算每个基因集的 Cv / Ct 矩阵（修正口径）----
  Vm <- Tm <- NULL; GSn <- character(0)
  Vlist <- Tlist <- list()
  for (z in use) {
    Pv <- z$val_pred[, keep_mods, drop=FALSE]; Pt <- z$test_pred[, keep_mods, drop=FALSE]
    Vlist[[z$name]] <- vapply(keep_mods, function(m) ci_new(z$tva, z$eva, Pv[,m]), numeric(1))
    Tlist[[z$name]] <- vapply(keep_mods, function(m) ci_new(z$tte, z$ete, Pt[,m]), numeric(1))
  }
  okgs <- names(Vlist)[vapply(Vlist, function(v) all(is.finite(v)), logical(1)) &
                       vapply(Tlist, function(v) all(is.finite(v)), logical(1))]
  Vm <- do.call(rbind, Vlist[okgs]); Tm <- do.call(rbind, Tlist[okgs])
  ref <- "Lasso"
  # 子池 S 的（跨基因集平均）乐观度
  bias_S <- function(S) {
    j <- match(S, keep_mods)
    cv <- Vm[, j, drop=FALSE]; ct <- Tm[, j, drop=FALSE]
    win <- max.col(cv, ties.method="first")
    vapply(seq_len(nrow(Vm)), function(i) (cv[i, win[i]] - Vm[i, ref]) - (ct[i, win[i]] - Tm[i, ref]),
           numeric(1)) |> mean(na.rm=TRUE)
  }
  cand <- setdiff(keep_mods, ref)
  set.seed(seed)
  rows <- list()
  for (k in KS) {
    if (k > K) next
    for (tg in cand) {
      pool <- setdiff(cand, tg)
      if (length(pool) < k - 1) next
      bb <- numeric(0)
      for (b in seq_len(B_BG)) {
        S <- sample(pool, k - 1)
        bb <- c(bb, bias_S(c(S, tg)) - bias_S(S))
      }
      rows[[length(rows)+1]] <- data.frame(cohort=coh, seed=seed, target=tg, k=k,
        entry=mean(bb, na.rm=TRUE), sd_bg=sd(bb, na.rm=TRUE), n_bg=length(bb))
    }
  }
  E <- do.call(rbind, rows); ENT[[paste(coh,seed)]] <- E
  # ---- 候选层面特征（用于预测因子分析）----
  dval <- sweep(Vm, 1, Vm[, ref], "-"); dtest <- sweep(Tm, 1, Tm[, ref], "-")
  win <- apply(Vm, 1, function(z) keep_mods[which.max(z)])
  agg <- aggregate(entry ~ target, data=E, FUN=function(x) mean(x, na.rm=TRUE))
  fe <- do.call(rbind, lapply(cand, function(m) {
    a <- dval[, m]; b <- dtest[, m]; q <- a >= quantile(a, .75, na.rm=TRUE)
    data.frame(cohort=coh, seed=seed, target=m, family=if (m %in% TREES) "tree" else "linear",
               entry_mean = agg$entry[match(m, agg$target)],
               mean_dval=mean(a,na.rm=TRUE), sd_dval=sd(a,na.rm=TRUE),
               winrate=mean(win==m,na.rm=TRUE), optimism=mean(a-b,na.rm=TRUE),
               tail_opt=mean((a-b)[q],na.rm=TRUE))
  }))
  FEA[[paste(coh,seed)]] <- fe
  say(sprintf("  %-5s seed=%d  %2d 基因集 × %2d 候选   %s",
              coh, seed, length(okgs), K, R$endpoint))
}

if (length(ENT)) {
  E <- do.call(rbind, ENT); write.table(E, file.path(TAB,"21_entry_effect.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  say(""); say("========== 平均边际进入效应（对随机背景池平均）==========")
  # 先按(队列×划分×目标×规模)已是值；对划分平均，再对队列等权
  a1 <- aggregate(entry ~ cohort + target + k, data=E, FUN=mean, na.rm=TRUE)
  a2 <- aggregate(entry ~ target + k, data=a1, FUN=mean, na.rm=TRUE)
  for (k in sort(unique(a2$k))) {
    s <- a2[a2$k==k, ]; s <- s[order(-s$entry), ]
    say(sprintf("  k=%2d : 最高 %s (%+.4f)   最低 %s (%+.4f)   中位 %+.4f",
        k, s$target[1], s$entry[1], s$target[nrow(s)], s$entry[nrow(s)], median(s$entry)))
  }
  say("")
  say("  注：这是【特定规模下的平均边际进入效应】，背景池随机抽取；")
  say("      不同规模之间不可直接比较（池组成同时变化）。不构成 Shapley 分解。")
}
if (length(FEA)) {
  F <- do.call(rbind, FEA); write.table(F, file.path(TAB,"21_predictors.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  say(""); say("========== 候选层面特征与进入效应的关联 ==========")
  say("  同一分析样本内的关联，不是机制证明。")
  fb <- aggregate(cbind(entry_mean, mean_dval, sd_dval, winrate, optimism, tail_opt) ~ cohort + target + family,
                  data=F, FUN=mean, na.rm=TRUE)
  for (v in c("mean_dval","sd_dval","winrate","optimism","tail_opt")) {
    ct <- cor.test(fb[[v]], fb$entry_mean)
    say(sprintf("    %-12s r = %+.3f   p = %.2e   (n=%d)", v, ct$estimate, ct$p.value, nrow(fb)))
  }
  say(""); say("  留一癌种检验（去掉一个队列后重算相关）：")
  loo <- do.call(rbind, lapply(unique(fb$cohort), function(cn) {
    s <- fb[fb$cohort != cn, ]
    data.frame(held_out=cn, n=nrow(s),
      r_mean_dval=cor(s$mean_dval, s$entry_mean), r_sd_dval=cor(s$sd_dval, s$entry_mean),
      r_winrate=cor(s$winrate, s$entry_mean), r_tail_opt=cor(s$tail_opt, s$entry_mean))
  }))
  write.table(loo, file.path(TAB,"21_loo.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  for (i in seq_len(nrow(loo))) say(sprintf(
    "    留下 %-6s : mean_dval %+.2f  sd_dval %+.2f  winrate %+.2f  tail_opt %+.2f",
    loo$held_out[i], loo$r_mean_dval[i], loo$r_sd_dval[i], loo$r_winrate[i], loo$r_tail_opt[i]))
  say("    （若去掉某一癌种后相关大幅下降甚至翻号，说明关联由该癌种驱动，不能当普遍规律）")
}
close(con)
