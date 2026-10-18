#!/usr/bin/env Rscript
# ============================================================
# 20 — 修正版分析：所有指标从 19 落盘的【逐患者预测】离线计算
#
# 与旧分析（10/11/12/18）的关键差别：
#  1. C-index 改为【直接评价】：concordance(Surv ~ s, reverse=TRUE, timewt="n")，
#     不在评价集上重新拟合，不做任何方向学习。
#  2. 同时计算【旧写法】的 C-index，用于量化"评价修正"本身造成多大差异
#     （评审要求：不能预判修正后的偏倚必然增大或减小，必须实测）。
#  3. 术语按评审意见统一：
#       Gv   = Cv(r) − Cv(b)          验证增益
#       Gt   = Ct(r) − Ct(b)          独立测试增益
#       Orel = Gv − Gt                【参照校正的选择乐观度】（不称"纯偏倚"）
#       Oraw = Cv(r) − Ct(r)          原始落差（另一指标，不可与 Orel 混称）
#  4. 汇总顺序：先按（队列 × 划分）内部对基因集求平均，再对队列等权汇总。
#  5. 报告每个检验的实际基因集数、可用候选数、以及可比对生存对数。
#
# 输入：results/main/19_*.rds
# 输出：results/tables/20_rules.tsv, 20_summary.tsv, 20_eval_fix.tsv, 20_coverage.tsv
# ============================================================
suppressMessages(library(survival))
options(stringsAsFactors = FALSE)
PROJ <- "D:/bioinfo/sigAudit"; TAB <- file.path(PROJ, "results/tables")
TREES <- c("RSF100", "GBM100", "GBM250", "RSF250")
# 重采样次数（对策规则）。这些规则每个基因集都要跑几十次 C-index，故取值适中。
B_VOTE <- 20
B_HALF <- 10

LOGF <- file.path(PROJ, "logs", "20_analysis_corrected.log"); con <- file(LOGF, "wt")
say <- function(...) { m <- paste0(...); cat(m, "\n", sep = ""); cat(m, "\n", sep = "", file = con); flush(con) }

## ---------- 两种 C-index 口径 ----------
ci_new <- function(t, e, s) {           # 直接评价，方向固定：高值=高风险
  s <- tryCatch(as.numeric(s), error = function(x) NULL)
  if (is.null(s) || length(s) != length(t)) return(NA_real_)
  ok <- is.finite(s) & !is.na(t) & !is.na(e)
  if (sum(ok) < 10) return(NA_real_)
  tryCatch(as.numeric(concordance(Surv(t[ok], e[ok]) ~ s[ok], reverse = TRUE, timewt = "n")$concordance),
           error = function(x) NA_real_)
}
ci_old <- function(t, e, s) {           # 旧写法：在评价集上重新拟合（泄漏方向）
  s <- tryCatch(as.numeric(s), error = function(x) NULL)
  if (is.null(s) || length(s) != length(t)) return(NA_real_)
  ok <- is.finite(s) & !is.na(t) & !is.na(e)
  if (sum(ok) < 10) return(NA_real_)
  tryCatch(as.numeric(summary(coxph(Surv(t[ok], e[ok]) ~ s[ok]))$concordance[1]),
           error = function(x) NA_real_)
}
# 可比对生存对数（评审要求报告）。
# 不能用 concordance()$count —— 它返回的是【向量】，塞进 data.frame 会被循环成多行幻影行。
# 直接按定义数：有序对 (i,j) 满足 t_i < t_j 且 i 发生了事件。
n_pairs <- function(t, e) {
  ok <- !is.na(t) & !is.na(e) & is.finite(t) & is.finite(e)
  t <- t[ok]; e <- e[ok]
  if (length(t) < 2) return(0)
  sum(outer(t, t, "<") & (e == 1))
}

## ---------- 规则实现（只用验证集信息选模型；测试集仅用于评价）----------
# 注意：V 是【命名数值向量】（每个候选一个 C-index），不是矩阵，不能写 V[, cols]
AMAX <- function(v, cols) { if (!length(cols) || all(is.na(v[cols]))) return(NA_character_)
  cols[which.max(v[cols])] }

ens_pred <- function(Pv, Pt, cols) {     # 用【验证集】统计量标准化后平均，固定应用到测试
  cols <- cols[!is.na(cols)]; if (!length(cols)) return(list(v = NULL, t = NULL))
  Mv <- Pv[, cols, drop = FALSE]; Mt <- Pt[, cols, drop = FALSE]
  mu <- colMeans(Mv); sdv <- apply(Mv, 2, sd); sdv[!is.finite(sdv) | sdv == 0] <- 1
  list(v = rowMeans(sweep(sweep(Mv, 2, mu, "-"), 2, sdv, "/")),
       t = rowMeans(sweep(sweep(Mt, 2, mu, "-"), 2, sdv, "/")))
}

rules_one <- function(z, ci = ci_new) {
  Pv <- z$val_pred; Pt <- z$test_pred; mods <- colnames(Pv)
  tva <- z$tva; eva <- z$eva; tte <- z$tte; ete <- z$ete
  V <- vapply(mods, function(m) ci(tva, eva, Pv[, m]), numeric(1)); names(V) <- mods
  T <- vapply(mods, function(m) ci(tte, ete, Pt[, m]), numeric(1)); names(T) <- mods
  cand <- setdiff(mods, "Lasso"); out <- list()
  put <- function(rn, cv, ct) out[[rn]] <<- c(cv, ct)
  put("none", V["Lasso"], T["Lasso"])
  j <- AMAX(V, cand); if (!is.na(j)) put("argmax_val", V[j], T[j])
  # 次优
  if (length(cand) >= 2) { o <- cand[order(-V[cand])]; j2 <- o[2]; put("runner_up", V[j2], T[j2]) }
  # --- 教科书式对策：bootstrap 验证集、多数投票选模型 ---
  # 只对【验证集】重采样；方向固定（高值=高风险），不做任何方向学习。
  set.seed(1)
  wins <- setNames(numeric(length(cand)), cand)
  for (b in seq_len(B_VOTE)) {
    ib <- sample(length(tva), replace = TRUE)
    if (sum(eva[ib]) < 5) next
    vb <- vapply(cand, function(m) ci(tva[ib], eva[ib], Pv[ib, m]), numeric(1))
    if (all(is.na(vb))) next
    w <- which.max(vb); wins[w] <- wins[w] + 1
  }
  if (sum(wins) > 0) { jm <- names(wins)[which.max(wins)]; put("majority_vote", V[jm], T[jm]) }
  # --- 随机二分验证集：取较差一半 / 均值减标准差 ---
  set.seed(1)
  acc <- setNames(vector("list", length(cand)), cand)
  for (b in seq_len(B_HALF)) {
    pm <- sample(length(tva)); h1 <- pm[1:floor(length(tva)/2)]; h2 <- setdiff(pm, h1)
    for (m in cand) {
      x1 <- ci(tva[h1], eva[h1], Pv[h1, m]); x2 <- ci(tva[h2], eva[h2], Pv[h2, m])
      if (is.finite(x1) && is.finite(x2)) acc[[m]] <- c(acc[[m]], x1, x2)
    }
  }
  hmin <- vapply(acc, function(z) if (length(z) >= 4) { k2 <- length(z) %/% 2 * 2
      mean(pmin(z[seq(1, k2, 2)], z[seq(2, k2, 2)])) } else NA_real_, numeric(1))
  hpen <- vapply(acc, function(z) if (length(z) >= 4) mean(z) - sd(z) else NA_real_, numeric(1))
  if (any(is.finite(hmin))) { jw <- cand[which.max(replace(hmin, !is.finite(hmin), -Inf))]
    put("worst_half", V[jw], T[jw]) }
  if (any(is.finite(hpen))) { jp <- cand[which.max(replace(hpen, !is.finite(hpen), -Inf))]
    put("penalized", V[jp], T[jp]) }
  lin <- setdiff(cand, TREES); tre <- intersect(cand, TREES)
  jl <- AMAX(V, lin); if (!is.na(jl)) put("linear_only", V[jl], T[jl])
  jt <- AMAX(V, tre); if (!is.na(jt)) put("tree_only", V[jt], T[jt])
  for (k in c(2, 3, 5, 10)) {
    tk <- names(sort(V[cand], decreasing = TRUE))[seq_len(min(k, length(cand)))]
    pr <- ens_pred(Pv, Pt, tk)
    if (!is.null(pr$v)) put(sprintf("top%d_ensemble", k), ci(tva, eva, pr$v), ci(tte, ete, pr$t))
  }
  pr <- ens_pred(Pv, Pt, cand)
  if (!is.null(pr$v)) put("all_ensemble", ci(tva, eva, pr$v), ci(tte, ete, pr$t))
  jo <- names(which.max(T[cand])); put("oracle_test", V[jo], T[jo])   # 仅作描述性参照
  out
}

## ---------- 读入 ----------
fs <- list.files(file.path(PROJ, "results/main"), pattern = "^19_.*\\.rds$", full.names = TRUE)
if (!length(fs)) stop("没有 19_*.rds，先跑 scripts/19_main_corrected.R")
say("找到 ", length(fs), " 个修正版主实验文件")
ROW <- list(); COV <- list(); FIX <- list(); P2 <- list()

for (f in fs) {
  R <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(R) || !length(R$genesets)) { say("[skip] ", basename(f)); next }
  coh <- R$cohort; seed <- R$seed; G <- R$genesets
  # 可用候选：在该（队列,划分）里出现率 >=90%
  cnt <- table(unlist(lapply(G, function(z) colnames(z$val_pred))))
  keep_mods <- names(cnt)[cnt / length(G) >= 0.90]
  use <- Filter(function(z) all(keep_mods %in% colnames(z$val_pred)), G)
  if (length(use) < 30) { say("[skip] ", coh, " seed", seed, " 完整基因集 <30"); next }
  # 流程链：计划抽样 -> 实际拟合成功 -> 可用候选筛选后共同完整任务
  # （评审要求正文区分这四者，且所有主图表用同一流程表）
  COV[[paste(coh, seed)]] <- data.frame(cohort = coh, seed = seed, endpoint = R$endpoint,
    n_planned = 100L, n_fitted = length(G),
    n_genesets = length(use), n_candidates = length(keep_mods),
    n_val = R$n_val, n_test = R$n_test,
    pairs_val = n_pairs(use[[1]]$tva, use[[1]]$eva),
    pairs_test = n_pairs(use[[1]]$tte, use[[1]]$ete))
  say(""); say("===== ", coh, " seed ", seed, " (", R$endpoint, ") : ",
      length(use), " 基因集 × ", length(keep_mods), " 候选 =====")

  # 逐基因集 × 逐规则
  A <- list(); Bo <- list()
  for (z in use) { z$val_pred <- z$val_pred[, keep_mods, drop=FALSE]
                   z$test_pred <- z$test_pred[, keep_mods, drop=FALSE]
    rn <- rules_one(z, ci_new); ro <- rules_one(z, ci_old)
    for (k in names(rn)) A[[length(A)+1]] <- data.frame(cohort=coh, seed=seed, geneset=z$name,
                                                        rule=k, cv=rn[[k]][1], ct=rn[[k]][2])
    for (k in names(ro)) Bo[[length(Bo)+1]] <- data.frame(cohort=coh, seed=seed, geneset=z$name,
                                                          rule=k, cv=ro[[k]][1], ct=ro[[k]][2])
  }
  A <- do.call(rbind, A); Bo <- do.call(rbind, Bo)
  # 逐规则：配对（相对 none）的 Gv / Gt / Orel，先在本（队列×划分）内对基因集平均
  a <- do.call(rbind, lapply(split(A, A$rule), function(d) {
    b <- A[A$rule == "none", ]; ii <- match(d$geneset, b$geneset)
    gv <- mean(d$cv - b$cv[ii], na.rm=TRUE); gt <- mean(d$ct - b$ct[ii], na.rm=TRUE)
    data.frame(cohort=coh, seed=seed, rule=d$rule[1], n=nrow(d),
               Cv=mean(d$cv, na.rm=TRUE), Ct=mean(d$ct, na.rm=TRUE),
               Gv=gv, Gt=gt, Orel=gv-gt, Oraw=mean(d$cv-d$ct, na.rm=TRUE)) }))
  ao <- do.call(rbind, lapply(split(Bo, Bo$rule), function(d) {
    b <- Bo[Bo$rule == "none", ]; ii <- match(d$geneset, b$geneset)
    data.frame(cohort=coh, seed=seed, rule=d$rule[1],
               Orel=mean((d$cv - b$cv[ii]) - (d$ct - b$ct[ii]), na.rm=TRUE),
               Oraw=mean(d$cv - d$ct, na.rm=TRUE)) }))
  ROW[[paste(coh, seed)]] <- a
  # ---- 主要比较 P2：相对 argmax_val 的配对差 ----
  # 参照改为【同队列、同划分、同基因集】的标准做法，而不是只跟固定 Lasso 比。
  ref_r <- A[A$rule == "argmax_val", ]
  P2[[paste(coh, seed)]] <- do.call(rbind, lapply(split(A, A$rule), function(d) {
    ii <- match(d$geneset, ref_r$geneset)
    dt <- d$ct - ref_r$ct[ii]
    do <- (d$cv - d$ct) - (ref_r$cv[ii] - ref_r$ct[ii])
    data.frame(cohort=coh, seed=seed, rule=d$rule[1],
               dCt=mean(dt, na.rm=TRUE), dOrel=mean(do, na.rm=TRUE),
               n_gs=sum(is.finite(dt) & is.finite(do)))
  }))
  FIX[[paste(coh, seed)]] <- merge(a[, c("cohort","seed","rule","Orel","Oraw","Cv","Ct")],
                                   ao, by=c("cohort","seed","rule"), suffixes=c("_new","_old"))
  for (i in seq_len(nrow(a))) if (a$rule[i] %in% c("none","argmax_val","linear_only","runner_up","tree_only")) say(sprintf(
    "  %-14s Cv %.4f  Ct %.4f   Gv %+.4f  Gt %+.4f  Orel %+.4f  Oraw %+.4f",
    a$rule[i], a$Cv[i], a$Ct[i], a$Gv[i], a$Gt[i], a$Orel[i], a$Oraw[i]))
}

if (!length(ROW)) stop("没有任何可用的（队列×划分）结果")
R <- do.call(rbind, ROW)
write.table(R, file.path(TAB, "20_rules.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(do.call(rbind, COV), file.path(TAB, "20_coverage.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
if (length(FIX)) write.table(do.call(rbind, FIX), file.path(TAB, "20_eval_fix.tsv"),
                             sep="\t", quote=FALSE, row.names=FALSE)

## ---------- 汇总：先按(队列×划分)得值，再对队列等权 ----------
say(""); say("==================== 汇总（先按队列×划分，再对队列等权）====================")
say("  划分种子: ", paste(sort(unique(R$seed)), collapse=", "), "   队列: ", length(unique(R$cohort)))
say("")
# 每个(队列,划分)一行已在 R 中；对「队列」先平均（跨划分），再对队列等权
byc <- aggregate(cbind(Cv, Ct, Gv, Gt, Orel, Oraw) ~ cohort + rule, data=R, FUN=mean, na.rm=TRUE)
agg <- aggregate(cbind(Cv, Ct, Gv, Gt, Orel, Oraw) ~ rule, data=byc, FUN=mean, na.rm=TRUE)
sd_ <- aggregate(Orel ~ rule, data=byc, FUN=function(x) diff(range(x, na.rm=TRUE)))
agg$Orel_range <- sd_$Orel[match(agg$rule, sd_$rule)]
agg <- agg[order(-agg$Ct), ]
say("  ", sprintf("%-14s %8s %8s %9s %9s %9s", "rule","Cv","Ct","Gv","Gt","Orel"))
for (i in seq_len(nrow(agg))) say("  ", sprintf("%-14s %8.4f %8.4f %+9.4f %+9.4f %+9.4f",
  agg$rule[i], agg$Cv[i], agg$Ct[i], agg$Gv[i], agg$Gt[i], agg$Orel[i]))

r1 <- agg[agg$rule=="argmax_val", ]; r0 <- agg[agg$rule=="none", ]
say("")
say(sprintf("  ★ 标准做法（argmax_val）：验证增益 Gv = %+.4f，独立测试增益 Gt = %+.4f",
            r1$Gv, r1$Gt))
say(sprintf("    → 参照校正的选择乐观度 Orel = %+.4f，占验证增益的 %.1f%%",
            r1$Orel, 100 * r1$Orel / r1$Gv))
say(sprintf("    （原始落差 Oraw 为 %+.4f —— 与 Orel 是不同指标，不可混称）", r1$Oraw))
say("")
say("  各规则相对参照（none）的配对比较：")
for (i in seq_len(nrow(agg))) say(sprintf(
  "    %-14s Gt %+.4f   Orel %+.4f（跨队列范围 %+.4f）  %s",
  agg$rule[i], agg$Gt[i], agg$Orel[i], agg$Orel_range[i],
  if (agg$rule[i]=="oracle_test") "（事后最大值，非总体上界，仅描述性参照）" else ""))

## ---------- 评价修正的影响（评审要求实测）----------
if (length(FIX)) {
  F <- do.call(rbind, FIX)
  byc2 <- aggregate(cbind(Orel_new, Orel_old, Oraw_new, Oraw_old) ~ cohort + rule, data=F, FUN=mean, na.rm=TRUE)
  a2 <- aggregate(cbind(Orel_new, Orel_old, Oraw_new, Oraw_old) ~ rule, data=byc2, FUN=mean, na.rm=TRUE)
  say(""); say("========== 评价函数修正的影响（同一批预测，两种口径对比）==========")
  say("  ", sprintf("%-14s %12s %12s %10s", "rule", "Orel(修正后)", "Orel(旧写法)", "差"))
  for (i in seq_len(nrow(a2))) say("  ", sprintf("%-14s %+12.4f %+12.4f %+10.4f",
    a2$rule[i], a2$Orel_new[i], a2$Orel_old[i], a2$Orel_new[i] - a2$Orel_old[i]))
  d <- a2[a2$rule=="argmax_val", ]
  say(sprintf("  → 标准做法处，修正前后的 Orel 差 %+.4f", d$Orel_new - d$Orel_old))
}
write.table(agg, file.path(TAB, "20_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## ---------- P2：各规则 vs 标准选优 argmax_val（评审指定的主要比较）----------
if (length(P2)) {
  P <- do.call(rbind, P2)
  write.table(P, file.path(TAB, "20_paired_vs_argmax.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  # 推断单位 = （队列 × 划分），不是基因集：约百个基因集共享患者，不视为独立样本
  say(""); say("========== 各规则 vs 标准做法 argmax_val（配对；推断单位 = 队列×划分）==========")
  say("  先在本（队列×划分）内对基因集配对平均，再对（队列×划分）单元做检验。")
  say("  差异不显著 != 等效或非劣；点估计接近零 != 已证实消除偏倚。")
  say("  ", sprintf("%-14s %10s %10s %8s %18s %s", "rule", "dCt", "dOrel", "n单元", "dCt 95%CI", "p"))
  u <- aggregate(cbind(dCt, dOrel) ~ rule, data=P, FUN=mean, na.rm=TRUE)
  u <- u[order(-u$dCt), ]
  for (i in seq_len(nrow(u))) {
    x <- P$dCt[P$rule == u$rule[i]]; x <- x[is.finite(x)]
    tt <- if (length(x) >= 4) t.test(x) else NULL
    say("  ", sprintf("%-14s %+10.4f %+10.4f %8d %18s %s", u$rule[i], u$dCt[i], u$dOrel[i], length(x),
        if (is.null(tt)) "n/a" else sprintf("[%+.4f,%+.4f]", tt$conf.int[1], tt$conf.int[2]),
        if (is.null(tt)) "" else sprintf("p=%.2e", tt$p.value)))
  }
  say("")
  say("  注：只报告效应量与区间，不自动生成“偏倚已消除”之类的结论标签。")
}
say(""); say("每个（队列×划分）的实际规模与可比对生存对数：")
CV <- do.call(rbind, COV)
say("  ", sprintf("%-6s %6s %6s %6s %6s %6s %6s %10s %10s",
       "cohort","seed","计划","拟合","完整","候选","测试n","可比对(验)","可比对(测)"))
for (i in seq_len(nrow(CV))) say("  ", sprintf("%-6s %6d %6d %6d %6d %6d %6d %10.0f %10.0f",
  CV$cohort[i], CV$seed[i], CV$n_planned[i], CV$n_fitted[i], CV$n_genesets[i],
  CV$n_candidates[i], CV$n_test[i], CV$pairs_val[i], CV$pairs_test[i]))
say("")
say(sprintf("  合计：计划 %d 个（队列×划分）任务；实际拟合成功 %d；共同完整分析任务 %d（%.1f%%）",
            sum(CV$n_planned), sum(CV$n_fitted), sum(CV$n_genesets),
            100 * sum(CV$n_genesets) / sum(CV$n_planned)))
say("  注意：'完整任务'指该基因集拥有【全部保留候选】且所有候选 C-index 均非缺失；")
say("        不同（队列×划分）的候选数可能不同，故不能用单一的 8×99×18 表述。")
close(con)
