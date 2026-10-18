#!/usr/bin/env Rscript
# ============================================================
# 21 — 进入效应与候选层面关联（第二轮评审意见落地）
#
# 评审指出的 5 点，逐条修正：
#  1. 【规模定义】原实现 KS 里含 k=18，但背景池排除了 Lasso（实际最多 17），
#     故 k=18 永远算不出。现在：**按最终可搜索池定义 K**（= 全部保留候选，
#     含参照，与 20 的主规则一致），KS 覆盖到 K 本身。
#  2. 【候选宇宙一致】原实现 winrate 从含 Lasso 的矩阵算，而进入效应排除 Lasso。
#     现在：进入效应、胜率、特征全部使用**同一个候选宇宙**。
#  3. 【并列处理】原实现 `bias_S(c(S,tg))` 总把目标放最后，配合
#     `max.col(ties.method="first")` 使并列永远偏向背景模型。
#     现在：使用与 20 相同的**固定随机优先级**破并，与输入顺序无关。
#  4. 【Monte Carlo】报告背景平均的标准误，并明确 `sd_bg` 是
#     **背景组成分散**，不是患者抽样置信区间。
#  5. 【多 k 平均】不再把多个 k 平均成一个 entry_mean；
#     如保留汇总值，明确标注为"人为等权摘要"。
#
# 另按评审要求：
#  · 留一分析正名为 leave-one-cohort-out **sensitivity of association**，
#    不再称"留一验证"（没有训练预测器、没有留出误差）。
#  · `tail_opt` 使用测试结局，**只能作事后解释变量，不能作事前预警特征**。
#  · 不同规模之间不做"纯数量因果效应"的强解释。
#
# 输入：results/manifest/main_input_files.txt
# 输出：results/tables/21_entry_effect.tsv, 21_entry_summary.tsv,
#       21_predictors.tsv, 21_loo_sensitivity.tsv
# ============================================================
suppressMessages(library(survival))
options(stringsAsFactors = FALSE)
PROJ <- Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit"); TAB <- file.path(PROJ, "results/tables")
TREES <- c("RSF100","GBM100","GBM250","RSF250")
B_BG  <- 50          # 每个 (目标, 规模) 的背景子池数（比 v1 的 20 更充分）

LOGF <- file.path(PROJ,"logs","21_entry_predictors.log"); con <- file(LOGF,"wt")
say <- function(...) { m <- paste0(...); cat(m,"\n",sep=""); cat(m,"\n",sep="",file=con); flush(con) }

ci_direct <- function(t, e, s) {
  s <- tryCatch(as.numeric(s), error=function(x) NULL)
  if (is.null(s) || length(s) != length(t)) return(NA_real_)
  ok <- is.finite(s) & !is.na(t) & !is.na(e)
  if (sum(ok) < 10) return(NA_real_)
  tryCatch(as.numeric(concordance(Surv(t[ok],e[ok]) ~ s[ok], reverse=TRUE, timewt="n")$concordance),
           error=function(x) NA_real_)
}
pick_max <- function(v, prio) {          # 与 20 相同的、与输入顺序无关的破并
  m <- suppressWarnings(max(v, na.rm=TRUE))
  if (!is.finite(m)) return(NA_integer_)
  cand <- which(is.finite(v) & v == m); cand[which.max(prio[cand])]
}

mani <- file.path(PROJ, "results/manifest/main_input_files.txt")
if (!file.exists(mani)) stop("缺少 manifest，请先跑 scripts/23_manifest.R")
fs <- readLines(mani, warn=FALSE); fs <- fs[nzchar(fs)]
say("按 manifest 读入 ", length(fs), " 个（队列×种子）产物；每个(目标,规模)的背景池数 B = ", B_BG)

ENT <- list(); FEA <- list()
for (f in fs) {
  R <- tryCatch(readRDS(f), error=function(e) NULL)
  if (is.null(R) || !length(R$genesets)) next
  coh <- R$cohort; seed <- R$seed
  cnt <- table(unlist(lapply(R$genesets, function(z) colnames(z$val_pred))))
  pool_all <- names(cnt)[cnt/length(R$genesets) >= 0.90]        # 可搜索池（含参照）
  if (!("Lasso" %in% pool_all) || length(pool_all) < 4) { say("[skip] ", coh); next }
  set.seed(seed); prio <- runif(length(pool_all)); names(prio) <- pool_all

  ## 与 20 一致的唯一任务清单
  Vlist <- list(); Tlist <- list()
  for (z in R$genesets) {
    if (!all(pool_all %in% colnames(z$val_pred))) next
    P <- z$val_pred[, pool_all, drop=FALSE]; Q <- z$test_pred[, pool_all, drop=FALSE]
    if (!all(is.finite(P)) || !all(is.finite(Q))) next
    V <- vapply(pool_all, function(m) ci_direct(z$tva, z$eva, P[,m]), numeric(1))
    Tt<- vapply(pool_all, function(m) ci_direct(z$tte, z$ete, Q[,m]), numeric(1))
    if (!all(is.finite(V)) || !all(is.finite(Tt))) next
    Vlist[[z$name]] <- V; Tlist[[z$name]] <- Tt
  }
  if (length(Vlist) < 30) { say("[skip] ", coh, " ", seed, " 任务 <30"); next }
  Vm <- do.call(rbind, Vlist); Tm <- do.call(rbind, Tlist)
  ref <- "Lasso"; K <- length(pool_all)
  # 子池 S 的（跨任务平均）乐观度：Orel = (Cv(胜者)−Cv(参照)) − (Ct(胜者)−Ct(参照))
  bias_S <- function(S) {
    j <- match(S, pool_all)
    cv <- Vm[, j, drop=FALSE]; ct <- Tm[, j, drop=FALSE]
    win <- vapply(seq_len(nrow(cv)), function(i) pick_max(cv[i,], prio[S]), integer(1))
    val <- vapply(seq_len(nrow(cv)), function(i)
      (cv[i, win[i]] - Vm[i, ref]) - (ct[i, win[i]] - Tm[i, ref]), numeric(1))
    mean(val, na.rm=TRUE)
  }
  cand <- setdiff(pool_all, ref)                    # 可作为"目标"的候选
  # 规模集合：覆盖到 K 本身（含参照的全池）；全池时目标加入后即整池
  KS <- sort(unique(pmin(c(2,3,5,8,12,K), K))); KS <- KS[KS >= 2]
  set.seed(seed + 7)
  rows <- list()
  for (k in KS) {
    for (tg in cand) {
      bg <- setdiff(pool_all, tg)
      if (length(bg) < k - 1) next
      if (k - 1 == length(bg)) {                    # 全池情形：只有一个背景集合，不假装抽样
        d <- bias_S(c(bg, tg)) - bias_S(bg)
        rows[[length(rows)+1]] <- data.frame(cohort=coh, seed=seed, target=tg, k=k,
          entry=d, sd_bg=NA_real_, mc_se=NA_real_, n_bg=1, full_pool=TRUE)
      } else {
        bb <- numeric(B_BG)
        for (b in seq_len(B_BG)) { S <- sample(bg, k - 1); bb[b] <- bias_S(c(S, tg)) - bias_S(S) }
        rows[[length(rows)+1]] <- data.frame(cohort=coh, seed=seed, target=tg, k=k,
          entry=mean(bb), sd_bg=sd(bb), mc_se=sd(bb)/sqrt(B_BG), n_bg=B_BG, full_pool=FALSE)
      }
    }
  }
  if (!length(rows)) { say("[skip] ", coh, " 无进入效应"); next }
  E <- do.call(rbind, rows); ENT[[paste(coh,seed)]] <- E

  ## 候选层面特征（与进入效应同一候选宇宙）
  dv <- sweep(Vm, 1, Vm[, ref], "-"); dt <- sweep(Tm, 1, Tm[, ref], "-")
  win <- apply(Vm, 1, function(z) pool_all[pick_max(z, prio)])
  full <- E[E$full_pool, ]                     # 用全池情形的进入效应作特征关联
  FEA[[paste(coh,seed)]] <- do.call(rbind, lapply(cand, function(m) {
    a <- dv[, m]; b <- dt[, m]; q <- a >= quantile(a, .75, na.rm=TRUE)
    z <- full[full$target == m, ]
    data.frame(cohort=coh, seed=seed, target=m,
      family = if (m %in% TREES) "tree" else "linear",
      entry_full = if (nrow(z)) z$entry[1] else NA_real_,
      mean_dval = mean(a, na.rm=TRUE), sd_dval = sd(a, na.rm=TRUE),
      winrate   = mean(win == m, na.rm=TRUE),
      optimism  = mean(a - b, na.rm=TRUE),
      tail_opt  = mean((a - b)[q], na.rm=TRUE))
  }))
  say(sprintf("  %-5s seed=%d  %2d 任务 × %2d 候选   K=%d  规模集合: %s",
              coh, seed, length(Vlist), K, K, paste(KS, collapse=",")))
}

if (length(ENT)) {
  E <- do.call(rbind, ENT)
  F <- do.call(rbind, FEA)
  Ff <- F[, c("cohort","seed","target","family")]; Ff$key <- paste(Ff$cohort, Ff$seed, Ff$target)
  E$key <- paste(E$cohort, E$seed, E$target)
  E$family <- Ff$family[match(E$key, Ff$key)]
  E$key <- NULL
  write.table(E, file.path(TAB,"21_entry_effect.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  a1 <- aggregate(entry ~ cohort + target + k + family, data=E, FUN=function(x) mean(x, na.rm=TRUE))
  a2 <- aggregate(entry ~ target + k + family, data=a1, FUN=function(x) mean(x, na.rm=TRUE))
  write.table(a2, file.path(TAB,"21_entry_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

  say(""); say("========== 平均边际进入效应（特定规模下、对随机背景池平均）==========")
  say("  k 为【加入目标后的池规模】；K 为可搜索池（含参照）。")
  for (k in sort(unique(a2$k))) {
    s <- a2[a2$k==k, ]; s <- s[order(-s$entry), ]
    say(sprintf("  k=%2d : 最高 %-14s %+.4f   最低 %-14s %+.4f   中位 %+.4f",
                k, s$target[1], s$entry[1], s$target[nrow(s)], s$entry[nrow(s)], median(s$entry)))
    mc <- E[E$k==k, ]
    mc <- mc[is.finite(mc$mc_se), ]
    if (nrow(mc)) say(sprintf("         背景抽样的 Monte Carlo 标准误（中位）: %.5f", median(mc$mc_se)))
  }
  say("")
  say("  注：① 这是【特定规模下、对随机背景子池平均】的边际进入效应，")
  say("      背景池组成随规模同时变化，因此【不能解释成与组成完全分离的纯数量因果效应】；")
  say("      ② 不构成 Shapley 分解（未做任何公理化分配）；")
  say("      ③ k=K 时背景池只有一个（全池），此时 n_bg=1、不假装抽样；")
  say("      ④ sd_bg 是【背景组成】的分散，不是患者抽样的置信区间。")

  ## ---------- 候选层面关联（明确标注为同一分析样本内的关联）----------
  Fagg <- aggregate(cbind(entry_full, mean_dval, sd_dval, winrate, optimism, tail_opt) ~
                      cohort + target + family, data=F, FUN=mean, na.rm=TRUE)
  write.table(Fagg, file.path(TAB,"21_predictors.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  say(""); say("========== 候选层面特征与【全池情形】进入效应的关联 ==========")
  say("  同一分析样本内的关联，不是机制证明；候选及其特征共享验证信息。")
  for (v in c("mean_dval","sd_dval","winrate","optimism","tail_opt")) {
    ct <- suppressWarnings(cor.test(Fagg[[v]], Fagg$entry_full))
    note <- if (v == "tail_opt") "  ← 使用测试结局，只能作事后解释变量" else ""
    say(sprintf("    %-12s r = %+.3f  p = %.2e  (n=%d)%s", v, ct$estimate, ct$p.value, nrow(Fagg), note))
  }
  say("")
  say("  留一癌种敏感性（leave-one-cohort-out sensitivity of association）——")
  say("  仅回答「关联是否被某个癌种驱动」，**不是**留出验证（没有训练预测器、没有留出误差）：")
  loo <- do.call(rbind, lapply(unique(Fagg$cohort), function(cn) {
    s <- Fagg[Fagg$cohort != cn, ]
    data.frame(held_out_cohort=cn, n_cohorts=length(unique(s$cohort)), n=nrow(s),
      r_mean_dval=cor(s$mean_dval, s$entry_full), r_sd_dval=cor(s$sd_dval, s$entry_full),
      r_winrate=cor(s$winrate, s$entry_full), r_tail_opt=cor(s$tail_opt, s$entry_full))
  }))
  write.table(loo, file.path(TAB,"21_loo_sensitivity.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  for (i in seq_len(nrow(loo))) say(sprintf(
    "    去掉 %-6s 后: mean_dval %+.2f  sd_dval %+.2f  winrate %+.2f  tail_opt %+.2f",
    loo$held_out_cohort[i], loo$r_mean_dval[i], loo$r_sd_dval[i], loo$r_winrate[i], loo$r_tail_opt[i]))
} else say("没有可用的进入效应结果")
close(con)
