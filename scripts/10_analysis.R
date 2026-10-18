#!/usr/bin/env Rscript
# ============================================================
# 10 — 离线分析：从 09 的模型级明细做【任意子集 / 任意排列】的偏倚分析
#
# 三个互相竞争的假说：
#   H1 剂量假说：偏倚来自池子【变大】本身
#   H2 成分假说：偏倚来自池子【含树集成类算法】
#   H3 单点假说：偏倚来自池子里【存在一个乐观度高的候选】
#
# 核心恒等式（整个设计的基础）：
#   bias = (val_sel − test_sel) − (val_ref − test_ref)
#   展开后 = 「被选模型的乐观度」−「参照模型的乐观度」。
#   两侧各自减同一个参照，因此【参照 Lasso 不必在子集内】，
#   从而可以在同样 k 下自由比较任意模型组合。
#
# 主估计量 ——「进入效应」entry_effect(target, p)：
#   把 target 放进一个已有 (p−1) 个模型、且 target 原本不在其中的子池：
#     Δ = bias(子池 ∪ {target}) − bias(子池)
#   对 p = 2..K 逐个测。
#   - H1 预测：Δ 随 p 单调增大（越靠后加越有效）——位置依赖
#   - H2 预测：只有树模型的 Δ 是大正值
#   - H3 预测：同一 target 在任何 p 处 Δ 都大（位置无关），
#              且 Δ 与该候选自身的乐观度 (Δval−Δtest) 成正比
#
# 输入：results/tables/09_<COHORT>_<N>.tsv
# 输出：10_model_summary.tsv, 10_entry_effect.tsv, 10_dose_curves.tsv,
#       10_tests.tsv, logs/10_analysis.log
# ============================================================
options(stringsAsFactors = FALSE)
PROJ <- "D:/bioinfo/sigAudit"

POOL_NAMES <- c("Lasso","Ridge","Enet0.9","Enet0.7","Enet0.5","Enet0.3","Enet0.1",
                "Ridge_1se","UniCox+Lasso","LassoTop100","LassoTop50","LassoTop25",
                "LassoTop10","Enet0.5Top25","RSF100","GBM100","GBM250","RSF250")
FAM <- setNames(ifelse(POOL_NAMES %in% c("RSF100","GBM100","GBM250","RSF250"),
                       "tree", "linear"), POOL_NAMES)

say2 <- local({ con <- file(file.path(PROJ,"logs","10_analysis.log"), "wt")
  function(...) { m <- paste0(...); cat(m,"\n",sep=""); cat(m,"\n",sep="",file=con); flush(con) } })
say2("10_analysis  ", format(Sys.time()))

files <- list.files(file.path(PROJ,"results/tables"), pattern="^09_.*_[0-9]+\\.tsv$", full.names=TRUE)
if (!length(files)) stop("没有 09_*.tsv，先跑 scripts/09_model_level.R")
say2("找到 ", length(files), " 个 09 明细表")

summ <- list(); dose <- list(); tests <- list(); entry <- list(); keepinfo <- list()

for (f in files) {
  coh <- sub("^09_(.+)_[0-9]+\\.tsv$", "\\1", basename(f))
  # 增量写入中的文件可能被截断 -> 读取失败就跳过，而不是整个分析崩掉
  D <- tryCatch(read.table(f, header=TRUE, sep="\t", check.names=FALSE), error=function(e) NULL)
  if (is.null(D) || !nrow(D) || !all(c("geneset","model","val_ci","test_ci") %in% names(D))) {
    say2("[skip] ", coh, " 读取失败或列不符（可能正在写入）"); next }
  D <- D[D$model %in% POOL_NAMES, ]
  if (!nrow(D)) { say2("[skip] ", coh, " 无可用行"); next }

  # ---- 可用性过滤：小队列里个别模型（如 UniCox+Lasso 需要 ≥5 个显著基因）
  #      可能在某些基因集上拟合失败。要求"18 个全齐"会让小队列全被丢掉。
  #      规则：先按【该队列内可用率 ≥90%】筛模型，再取"这些模型全齐"的基因集。
  all_gs <- unique(D$geneset); n_all <- length(all_gs)
  cnt <- tapply(D$geneset, D$model, function(g) length(unique(g)))   # 每个模型覆盖多少基因集
  rate <- cnt / n_all
  pnames <- POOL_NAMES[POOL_NAMES %in% names(rate)[rate >= 0.90]]
  dropped <- setdiff(POOL_NAMES, pnames)
  if (!("Lasso" %in% pnames)) { say2("[skip] ", coh, " 参照 Lasso 不可用"); next }
  sub <- D[D$model %in% pnames, ]
  tab <- table(sub$geneset); keep <- names(tab)[tab == length(pnames)]
  if (!length(keep)) { say2("[skip] ", coh, " 无全齐基因集（可用模型 ", length(pnames), "）"); next }
  mk <- function(v) { M <- tapply(sub[[v]], list(sub$geneset, sub$model), function(x) x[1])
    M <- M[keep, pnames, drop=FALSE]
    matrix(as.numeric(M), nrow=nrow(M), dimnames=dimnames(M)) }
  V <- mk("val_ci"); Tt <- mk("test_ci")
  ok <- complete.cases(V) & complete.cases(Tt)
  V <- V[ok,,drop=FALSE]; Tt <- Tt[ok,,drop=FALSE]
  n_gs <- nrow(V); K <- length(pnames)
  keepinfo[[coh]] <- data.frame(cohort=coh, n_genesets=n_gs, n_models=K,
    dropped=if (length(dropped)) paste(dropped, collapse=",") else "-")
  say2(""); say2("===== ", coh, " : ", n_gs, " 个基因集, ", K, " 个可用候选",
       if (length(dropped)) paste0("  (因可用率<90% 剔除: ", paste(dropped, collapse=","), ")") else "", " =====")
  if (n_gs < 10 || K < 5) { say2("  样本或候选太少，跳过"); next }
  ref_v <- V[,"Lasso"]; ref_t <- Tt[,"Lasso"]

  # ---- 逐模型：与参照相比的乐观度 ----
  sm <- data.frame(model=pnames, family=FAM[pnames],
                   val=colMeans(V), test=colMeans(Tt), row.names=NULL)
  sm$dval <- sm$val - mean(ref_v); sm$dtest <- sm$test - mean(ref_t)
  sm$optimism <- sm$dval - sm$dtest
  summ[[coh]] <- sm
  for (i in seq_len(nrow(sm))) say2(sprintf(
    "    %-14s %-7s val %.3f test %.3f  Δval %+.4f Δtest %+.4f  乐观度 %+.4f",
    sm$model[i], sm$family[i], sm$val[i], sm$test[i], sm$dval[i], sm$dtest[i], sm$optimism[i]))

  # ---- 偏倚：子集内按 val 选模型，量它"乐观度超出参照"多少 ----
  bias_at <- function(order, k) {
    s2 <- order[seq_len(k)]
    M <- V[, s2, drop=FALSE]
    # max.col 逐行返回 argmax 的【列号】；which.max(M) 返回的是整矩阵线性下标，是错的
    idx <- max.col(M, ties.method="first")
    vapply(seq_len(n_gs), function(i) {
      m <- s2[idx[i]]; (V[i,m] - ref_v[i]) - (Tt[i,m] - ref_t[i]) }, numeric(1))
  }

  # ---- 1) 典型排列下的剂量曲线 ----
  interleave <- function(a,b) { o <- character(0); i <- 1; j <- 1
    while (i<=length(a) || j<=length(b)) {
      if (i<=length(a)) { o <- c(o,a[i]); i <- i+1 }
      if (j<=length(b)) { o <- c(o,b[j]); j <- j+1 } }; o }
  TRE <- pnames[FAM[pnames]=="tree"]; LIN <- pnames[FAM[pnames]=="linear"]
  ORDERS <- list(A_canonical=pnames)
  if (length(TRE) && length(LIN)) ORDERS <- c(ORDERS, list(
    B_tree_first = c(TRE, LIN), C_interleaved = interleave(LIN, TRE),
    D_tree_at3   = c(LIN[1:2], TRE, LIN[3:length(LIN)])))
  for (on in names(ORDERS)) {
    o <- ORDERS[[on]]
    dose[[paste0(coh,"|",on)]] <- do.call(rbind, lapply(2:length(o), function(k) data.frame(
      cohort=coh, order=on, k=k, model_at_k=o[k], n_tree=sum(o[seq_len(k)] %in% TRE),
      bias=mean(bias_at(o,k)), n=n_gs)))
  }
  bc <- dose[[paste0(coh,"|A_canonical")]]
  d <- diff(bc$bias)                 # d[j] = bias(k=j+1 -> j+2)，进入者是位置 j+2
  jm <- which.max(d)
  say2(sprintf("  典型排列：偏倚 %+.4f (k=2) -> %+.4f (k=%d)；最大单步增量在 k=%d->%d，「%s」进入 %+.4f",
               bc$bias[1], bc$bias[nrow(bc)], K, jm+1, jm+2, bc$model_at_k[jm+1], d[jm]))

  # ---- 2) 主估计量：进入效应（是否位置无关）----
  for (tg in setdiff(pnames, "Lasso")) {
    rest <- setdiff(pnames, tg)
    for (p in 2:K) {
      o <- append(rest, tg, after = p-1)
      dd <- bias_at(o, p) - bias_at(rest[seq_len(p-1)], p-1)
      dd <- dd[is.finite(dd)]
      if (length(dd) < 10) next
      q <- t.test(dd)
      entry[[paste0(coh,"|",tg,"|",p)]] <- data.frame(
        cohort=coh, target=tg, family=FAM[tg], pos=p, n=length(dd),
        mean=mean(dd), ci_lo=q$conf.int[1], ci_hi=q$conf.int[2], p=q$p.value)
    }
  }

  # ---- 3) 关键配对检验 ----
  tt <- function(a,b,lab) { dd <- a-b; ok <- is.finite(dd); if (sum(ok)<10) return(NULL)
    q <- t.test(dd[ok]); w <- suppressWarnings(wilcox.test(dd[ok]))
    data.frame(cohort=coh, contrast=lab, n=sum(ok), mean_diff=mean(dd[ok]),
               ci_lo=q$conf.int[1], ci_hi=q$conf.int[2], t_p=q$p.value, wilcox_p=w$p.value) }
  by_opt <- order(-sm$optimism)
  worst  <- sm$model[by_opt][1]      # 最乐观的候选
  pw <- match(worst, pnames)
  # 家族对比必须【两侧同 k】，且两侧都【排除参照 Lasso】：
  #   若线性子池含 Lasso，它可以直接"选中参照自己"使偏倚恰为 0，
  #   而树子池没有这个选项 -> 不对称，会系统性夸大树的危害。
  LINc <- setdiff(LIN, "Lasso")
  kk <- min(4, length(TRE), length(LINc))
  if (kk >= 2) tests[[paste0(coh,"_1")]] <- tt(bias_at(TRE, kk), bias_at(LINc, kk),
      sprintf("同样 k=%d 且均不含参照：树集子池 − 线性子池", kk))
  tests[[paste0(coh,"_2")]] <- tt(bias_at(pnames,K), bias_at(LIN,length(LIN)),
      sprintf("全池(%d) − 纯线性(%d)", K, length(LIN)))
  if (pw >= 2 && pw <= K) {
    tests[[paste0(coh,"_3")]] <- tt(bias_at(pnames,pw), bias_at(pnames,pw-1),
        sprintf("最乐观候选「%s」进入 (k=%d->%d)", worst, pw-1, pw))
  }
  if (length(TRE) >= 1) {
    pwt <- match(TRE[1], pnames)
    if (!is.na(pwt) && pwt >= 2) tests[[paste0(coh,"_4")]] <- tt(
      bias_at(pnames,pwt), bias_at(pnames,pwt-1),
      sprintf("首个树集成「%s」进入 (k=%d->%d)", TRE[1], pwt-1, pwt))
  }
}

# ---------------- 输出 ----------------
w <- function(x, nm) if (length(x)) write.table(do.call(rbind,x), file.path(PROJ,"results/tables",nm),
                                                sep="\t", quote=FALSE, row.names=FALSE)
w(dose,"10_dose_curves.tsv"); w(tests,"10_tests.tsv"); w(entry,"10_entry_effect.tsv")
w(keepinfo,"10_keepinfo.tsv")
if (length(summ)) write.table(do.call(rbind, lapply(names(summ), function(cn) cbind(cohort=cn, summ[[cn]]))),
                              file.path(PROJ,"results/tables","10_model_summary.tsv"),
                              sep="\t", quote=FALSE, row.names=FALSE)

if (length(entry)) {
  EE <- do.call(rbind, entry)
  say2(""); say2("========== 主估计量：各候选的『进入效应』==========")
  say2("  进入效应 = 把该候选加入已有 p-1 个候选的子池后，偏倚的变化量（对 p 取平均）")
  say2("  H1(剂量)->随 p 单调增大；H2(成分)->仅树为大正值；H3(单点)->位置无关且∝自身乐观度")
  ag <- aggregate(mean ~ cohort + target + family, data=EE, FUN=mean)
  sdv <- aggregate(mean ~ cohort + target, data=EE, FUN=sd)
  ag$sd_over_pos <- sdv$mean[match(paste(ag$cohort,ag$target), paste(sdv$cohort,sdv$target))]
  for (cn in unique(ag$cohort)) {
    a <- ag[ag$cohort==cn, ]; a <- a[order(-a$mean), ]
    say2("  --- ", cn, " ---")
    for (i in seq_len(nrow(a))) say2(sprintf(
      "    %-14s %-7s 平均进入效应 %+.4f   跨位置SD %.4f", a$target[i], a$family[i], a$mean[i], a$sd_over_pos[i]))
  }
  mrg <- merge(EE, do.call(rbind, lapply(names(summ), function(cn) cbind(cohort=cn, summ[[cn]]))),
               by.x=c("cohort","target"), by.y=c("cohort","model"))
  say2(""); say2("  进入效应 vs 该候选自身乐观度：")
  for (cn in unique(mrg$cohort)) { s <- mrg[mrg$cohort==cn, ]
    ct <- cor.test(s$optimism, s$mean, use="complete.obs")
    say2(sprintf("    %-6s  Pearson r = %+.3f  p = %.2e  (n=%d)", cn, ct$estimate, ct$p.value, nrow(s))) }
  say2(""); say2("  跨队列合并：")
  ct <- cor.test(mrg$optimism, mrg$mean, use="complete.obs")
  say2(sprintf("    ALL    Pearson r = %+.3f  p = %.2e  (n=%d)", ct$estimate, ct$p.value, nrow(mrg)))
}
if (length(tests)) {
  TS <- do.call(rbind, tests)
  say2(""); say2("========== 关键配对检验 ==========")
  for (i in seq_len(nrow(TS))) say2(sprintf(
    "  %-6s %-38s n=%3d 差 %+.4f [%+.4f,%+.4f] t-p=%.2e wilcox-p=%.2e",
    TS$cohort[i], TS$contrast[i], TS$n[i], TS$mean_diff[i],
    TS$ci_lo[i], TS$ci_hi[i], TS$t_p[i], TS$wilcox_p[i]))
}
say2(""); say2("每个队列最终使用的基因集数 / 候选数：")
if (length(keepinfo)) { KI <- do.call(rbind, keepinfo)
  for (i in seq_len(nrow(KI))) say2(sprintf("  %-6s n=%3d  K=%2d  剔除: %s",
    KI$cohort[i], KI$n_genesets[i], KI$n_models[i], KI$dropped[i])) }
