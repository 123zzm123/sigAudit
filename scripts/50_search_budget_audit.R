#!/usr/bin/env Rscript
# ============================================================
# 50 — 基线条件下的"搜索规模曲线"（本次修订新增，2026-10-19，事后探索）
#
# 问题：当搜索从一个固定模型（起点）逐步扩展候选池时，新增候选带来的是
#      留出收益还是主要抬高验证表现？答案是否取决于起点是 LASSO 还是 Ridge？
#
# 设计（参考锚定搜索，reference-anchored search）：
#   · 使用原有 3,919 个共同可评估任务及其【实际分析候选池】，不重新筛选；
#   · 起点固定为 Lasso / Ridge（两个起点分别计算，每个任务每个起点都算）；
#   · 含 N 个候选的任务，搜索规模 K 表示子池大小：起点始终在池中，其余 K−1 个
#     从剩余 N−1 个候选中【均匀无放回】抽取；
#   · 选优只按【验证】C-index；并列优先级沿用原完整池预生成的 make_prio(mods, seed)，
#     不在每个子池重新生成随机数；测试 C-index 只用于评价；
#   · K ∈ {1, 2, 4, 8, 12, full}；full = 该任务自己的完整池。
#
# 随机子池的期望表现【精确计算】，不做蒙特卡洛：
#   按验证 C 从高到低、以固定优先级破并得严格顺序 i=1..N，起点位置为 r。
#   对 K>=2：p_r = C(N-r, K-1) / C(N-1, K-1)；
#            p_i = C(N-i-1, K-2) / C(N-1, K-1)，i < r；i > r 时为 0。
#   K=1 时直接取起点（概率 1）。期望 Cv/Ct = 各候选 C 的上述概率加权和。
#   同时给出"子池内均匀随机选一个模型"的期望 C（对照），
#   E_random = [C_start + (K-1)/(N-1) * sum_{j≠start} C_j] / K。
#
# 口径：TCGA 任务 -> (队列×划分) 平均 -> 队列内平均 -> 队列等权；
#      METABRIC 任务 -> 签名内跨 seed 平均 -> 六签名等权。
# 属性：这是【已见结果后】的修订性探索分析，不是事前预设；不按测试表现选择
#      "最优 K"，不给每个 K 配探索性 p 值。
#
# 输入：results/tables/20_model_metrics_long.tsv, 20_rule_performance_long.tsv,
#       31_metabric_model_metrics.tsv, 31_metabric_rule_performance.tsv
# 输出：50_search_budget_task.tsv / _by_cohort.tsv / _summary.tsv /
#       _os6.tsv / _metabric.tsv / _metabric_task.tsv / _check.tsv + logs/50_search_budget.log
# ============================================================
options(stringsAsFactors = FALSE)
PROJ <- Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
TAB  <- file.path(PROJ, "results/tables"); LOGD <- file.path(PROJ, "logs")
dir.create(LOGD, showWarnings = FALSE, recursive = TRUE)
con <- file(file.path(LOGD, "50_search_budget.log"), "wt")
say <- function(...) { m <- paste0(...); cat(m,"\n",sep=""); cat(m,"\n",sep="",file=con); flush(con) }
rd <- function(f) read.table(f, sep="\t", header=TRUE, quote="", comment.char="", stringsAsFactors=FALSE)
.lib <- file.path(PROJ, "scripts/_lib_analysis.R"); if (!file.exists(.lib)) stop("缺少 _lib_analysis.R")
source(.lib)   # make_prio / order_rule（并列口径与原流水线完全一致）
KS <- c(1L, 2L, 4L, 8L, 12L)
STARTS <- c("Lasso", "Ridge")

## ---- 精确子池选择概率：起点位置 r，池大小 N，搜索规模 K ----
prob_pick <- function(N, r, K) {
  p <- numeric(N)
  if (K == 1L) { p[r] <- 1; return(p) }
  den <- lchoose(N - 1, K - 1)
  for (i in seq_len(N)) {
    if (i < r) {
      a <- N - i - 1; b <- K - 2L
      p[i] <- if (a < b || b < 0) 0 else exp(lchoose(a, b) - den)
    } else if (i == r) {
      p[i] <- if (N - r < K - 1L) 0 else exp(lchoose(N - r, K - 1L) - den)
    }
  }
  p
}

## ---- 单个任务的搜索规模曲线（给定候选的验证/测试 C 与固定顺序）----
budget_one <- function(v, t, start, Ks = KS, label = "") {
  N <- length(v)
  prio <- attr(v, "prio")
  ord <- order(-v, -prio)                      # 严格顺序（与 order_rule 一致）
  r <- match(start, names(v)[ord])
  if (is.na(r)) return(NULL)
  vv <- v[ord]; tt <- t[ord]
  Kset <- sort(unique(c(Ks[Ks <= N], N)))
  out <- lapply(Kset, function(K) {
    p <- prob_pick(N, r, K)
    data.frame(K = K, is_full = as.integer(K == N),
               Cv_start = vv[r], Ct_start = tt[r],
               eCv = sum(p * vv), eCt = sum(p * tt),
               eRandCv = (vv[r] + (K - 1)/(N - 1) * (sum(vv) - vv[r])) / K,
               eRandCt = (tt[r] + (K - 1)/(N - 1) * (sum(tt) - tt[r])) / K,
               p_sum = sum(p), stringsAsFactors = FALSE)
  })
  D <- do.call(rbind, out)
  D$Gv <- D$eCv - D$Cv_start
  D$Gt <- D$eCt - D$Ct_start
  D$Orel <- D$Gv - D$Gt
  D
}

## ============================================================
## A. TCGA：3,919 个共同可评估任务
## ============================================================
MM <- rd(file.path(TAB,"20_model_metrics_long.tsv"))
MM$task <- paste(MM$cohort, MM$seed, MM$geneset, sep="|")
RP <- rd(file.path(TAB,"20_rule_performance_long.tsv"))
RP$task <- paste(RP$cohort, RP$seed, RP$geneset, sep="|")
am <- RP[RP$rule=="argmax_val", c("task","cv","ct")]
names(am)[2:3] <- c("am_cv","am_ct")

say("TCGA 候选指标表：", nrow(MM), " 行，", length(unique(MM$task)), " 个任务")

TASKS <- split(MM, MM$task)
TASK_ROWS <- list()
for (tk in names(TASKS)) {
  g <- TASKS[[tk]]
  mods <- sort(unique(g$model))
  if (length(mods) < 2 || !all(STARTS %in% mods)) next
  v <- setNames(g$val_ci[match(mods, g$model)], mods)
  t <- setNames(g$test_ci[match(mods, g$model)], mods)
  if (any(!is.finite(v)) || any(!is.finite(t))) next
  prio <- make_prio(mods, g$seed[1])
  v <- v[mods]; attr(v, "prio") <- prio[mods]
  for (st in STARTS) {
    B <- budget_one(v, t, st)
    if (is.null(B)) next
    TASK_ROWS[[length(TASK_ROWS)+1]] <- data.frame(
      source="TCGA", cohort=g$cohort[1], geneset=g$geneset[1], seed=g$seed[1],
      start=st, N=length(mods), B, stringsAsFactors=FALSE)
  }
}
TD <- do.call(rbind, TASK_ROWS); rownames(TD) <- NULL
## full 在所有任务上归为同一个 K 层级（K 列等于各自池大小 N，若按 K 分组会把
## full 拆成 N=15/16/17/18 多组，破坏队列等权）。Kgrp=0 表示 full。
TD$Kgrp <- ifelse(TD$is_full == 1L, 0L, TD$K)
## 直接写全精度（不预先 round）：显示时再舍入，避免建立第二套数据版本。
write.table(TD[, setdiff(colnames(TD), "Kgrp")],
            file.path(TAB,"50_search_budget_task.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
say("任务×起点×K 记录：", nrow(TD), "；任务 ", length(unique(paste(TD$cohort,TD$seed,TD$geneset))))

## ---- 三层级汇总：任务 -> (队列×划分) -> 队列 -> 队列等权 ----
MC <- c("eCv","eCt","Gv","Gt","Orel","eRandCv","eRandCt")
addlab <- function(D) {
  D$K_label <- ifelse(D$is_full == 1L, "full", as.character(D$Kgrp))
  D$K <- ifelse(D$is_full == 1L, NA_integer_, D$Kgrp)
  D[order(D$start, ifelse(D$Kgrp == 0L, 999L, D$Kgrp)), ]
}
agg3 <- function(df) {
  s1 <- aggregate(df[MC], by=df[c("cohort","seed","start","Kgrp","is_full")], FUN=mean, na.rm=TRUE)
  s2 <- aggregate(s1[MC], by=s1[c("cohort","start","Kgrp","is_full")], FUN=mean, na.rm=TRUE)
  aggregate(s2[MC], by=s2[c("start","Kgrp","is_full")], FUN=mean, na.rm=TRUE)
}
SUM <- addlab(agg3(TD)); rownames(SUM) <- NULL
nt <- aggregate(list(n_tasks=TD$eCv), by=TD[c("start","Kgrp","is_full")], FUN=function(x) sum(is.finite(x)))
SUM$n_tasks <- nt$n_tasks[match(paste(SUM$start,SUM$Kgrp,SUM$is_full), paste(nt$start,nt$Kgrp,nt$is_full))]
SUM$n_cohort <- 8L
SUM <- SUM[, c("start","K","K_label","is_full",MC,"n_tasks","n_cohort")]
write.table(SUM, file.path(TAB,"50_search_budget_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## ---- 队列级 ----
s1c <- aggregate(TD[MC], by=TD[c("cohort","seed","start","Kgrp","is_full")], FUN=mean, na.rm=TRUE)
COH <- aggregate(s1c[MC], by=s1c[c("cohort","start","Kgrp","is_full")], FUN=mean, na.rm=TRUE)
COH <- addlab(COH); COH <- COH[order(COH$start, ifelse(COH$Kgrp==0L,999L,COH$Kgrp), COH$cohort), ]; rownames(COH) <- NULL
COH <- COH[, c("cohort","start","K","K_label","is_full",MC)]
write.table(COH, file.path(TAB,"50_search_budget_by_cohort.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## ---- 六 OS 队列 ----
OS6 <- c("BRCA","COAD","KIRC","LUAD","MESO","STAD")
D6 <- TD[TD$cohort %in% OS6, ]
SUM6 <- addlab(agg3(D6)); rownames(SUM6) <- NULL
nt6 <- aggregate(list(n_tasks=D6$eCv), by=D6[c("start","Kgrp","is_full")], FUN=function(x) sum(is.finite(x)))
SUM6$n_tasks <- nt6$n_tasks[match(paste(SUM6$start,SUM6$Kgrp,SUM6$is_full), paste(nt6$start,nt6$Kgrp,nt6$is_full))]
SUM6$n_cohort <- 6L
SUM6 <- SUM6[, c("start","K","K_label","is_full",MC,"n_tasks","n_cohort")]
write.table(SUM6, file.path(TAB,"50_search_budget_os6.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## ============================================================
## B. METABRIC：签名内跨 seed 平均 -> 六签名等权
## ============================================================
MMb <- rd(file.path(TAB,"31_metabric_model_metrics.tsv"))
RPb <- rd(file.path(TAB,"31_metabric_rule_performance.tsv"))
mb_rows <- list()
for (key in split(MMb, paste(MMb$signature, MMb$seed, sep="|"))) {
  mods <- sort(unique(key$model))
  if (length(mods) < 2 || !all(STARTS %in% mods)) next
  v <- setNames(key$cv[match(mods, key$model)], mods)
  t <- setNames(key$ct[match(mods, key$model)], mods)
  if (any(!is.finite(v)) || any(!is.finite(t))) next
  prio <- make_prio(mods, key$seed[1])
  attr(v, "prio") <- prio[mods]
  for (st in STARTS) {
    B <- budget_one(v, t, st)
    if (is.null(B)) next
    mb_rows[[length(mb_rows)+1]] <- data.frame(
      signature=key$signature[1], seed=key$seed[1], start=st, N=length(mods), B,
      stringsAsFactors=FALSE)
  }
}
MBt <- do.call(rbind, mb_rows); rownames(MBt) <- NULL
MBt$Kgrp <- ifelse(MBt$is_full == 1L, 0L, MBt$K)
## 签名内跨 seed 平均 -> 六签名等权（"ALL" 行）
s1b <- aggregate(MBt[MC], by=MBt[c("signature","start","Kgrp","is_full")], FUN=mean, na.rm=TRUE)
ALLb <- aggregate(s1b[MC], by=s1b[c("start","Kgrp","is_full")], FUN=mean, na.rm=TRUE)
ALLb <- data.frame(signature="ALL", ALLb, stringsAsFactors=FALSE)
MB <- rbind(ALLb, s1b)
MB <- addlab(MB)
MB <- MB[order(MB$signature != "ALL", MB$start, ifelse(MB$Kgrp==0L,999L,MB$Kgrp), MB$signature), ]
MB$n_tasks <- NA_integer_
for (i in seq_len(nrow(MB))) {
  if (MB$signature[i] == "ALL") {
    MB$n_tasks[i] <- sum(MBt$start==MB$start[i] & MBt$Kgrp==MB$Kgrp[i] & MBt$is_full==MB$is_full[i])
  } else {
    MB$n_tasks[i] <- sum(MBt$signature==MB$signature[i] & MBt$start==MB$start[i] &
                           MBt$Kgrp==MB$Kgrp[i] & MBt$is_full==MB$is_full[i])
  }
}
MB <- MB[, c("signature","start","K","K_label","is_full",MC,"n_tasks")]
write.table(MB, file.path(TAB,"50_search_budget_metabric.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## 与 TCGA 逐任务表对称：写出 METABRIC 的任务×起点×K 明细（360 行）
MBT_OUT <- MBt[, c("signature","seed","start","N","K","is_full",
                   "Cv_start","Ct_start","eCv","eCt","eRandCv","eRandCt",
                   "Gv","Gt","Orel","p_sum")]
MBT_OUT <- MBT_OUT[order(MBT_OUT$signature, MBT_OUT$seed, MBT_OUT$start,
                         ifelse(MBT_OUT$is_full == 1L, 999L, MBT_OUT$K)), ]
write.table(MBT_OUT, file.path(TAB,"50_search_budget_metabric_task.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

## METABRIC 的 K=full 与旧 argmax 一致性核对
MBt$key <- paste(MBt$signature, MBt$seed, sep="|")
amb <- RPb[RPb$rule=="argmax_val", c("signature","seed","cv","ct")]
amb$key <- paste(amb$signature, amb$seed, sep="|")
FUb <- merge(MBt[MBt$is_full==1, c("key","start","eCv","eCt")], amb[,c("key","cv","ct")], by="key", all.x=TRUE)
mb_full_dev_Cv <- max(abs(FUb$eCv - FUb$cv), na.rm=TRUE)
mb_full_dev_Ct <- max(abs(FUb$eCt - FUb$ct), na.rm=TRUE)
say(sprintf("METABRIC full 与旧 argmax：Cv max|dev|=%.2e，Ct max|dev|=%.2e", mb_full_dev_Cv, mb_full_dev_Ct))

## ============================================================
## C. 核对
## ============================================================
## C1. 每个任务的概率和为一、K=1 取起点、K=full 与旧 argmax 一致
pdev <- max(abs(TD$p_sum - 1), na.rm=TRUE)
k1 <- TD[TD$K == 1, ]
k1_dev_Cv <- max(abs(k1$eCv - k1$Cv_start)); k1_dev_Ct <- max(abs(k1$eCt - k1$Ct_start))
## 旧 argmax（20 表）：与 K=full 的期望 Cv/Ct 逐任务比较
TD$task <- paste(TD$cohort, TD$seed, TD$geneset, sep="|")
FU <- merge(TD[TD$is_full == 1, c("task","start","eCv","eCt")], am, by="task", all.x=TRUE)
full_dev_Cv <- max(abs(FU$eCv - FU$am_cv), na.rm=TRUE)
full_dev_Ct <- max(abs(FU$eCt - FU$am_ct), na.rm=TRUE)

## C2. 小规模穷举一致性：N=2..8，所有起点位置 r 与 K，用枚举子池验证概率公式
exh_max <- 0
for (N in 2:8) for (r in 1:N) for (K in 2:N) {
  p <- prob_pick(N, r, K)
  cmb <- combn(N, K)
  cnt <- numeric(N); nsel <- 0
  for (j in seq_len(ncol(cmb))) { s <- cmb[, j]; if (!(r %in% s)) next
    nsel <- nsel + 1; cnt[min(s)] <- cnt[min(s)] + 1 }
  pe <- cnt / nsel
  exh_max <- max(exh_max, abs(p - pe))
}

## C3. 期望验证表现随 K 非下降（起点始终在池中）
TD$key <- paste(TD$task, TD$start, sep="@@")
mono_viol <- 0; mono_tot <- 0
for (k in unique(TD$key)) {
  z <- TD[TD$key==k, ]; z <- z[order(z$K), ]
  if (nrow(z) < 2) next
  mono_tot <- mono_tot + 1
  if (any(diff(z$eCv) < -1e-12)) mono_viol <- mono_viol + 1
}

## C4. K=1 时 Gv/Gt/Orel 为零
k1_dev_G <- max(abs(c(k1$Gv, k1$Gt, k1$Orel)))

## C5. 任务数与任务键集合：直接与预期规模比较（不是无条件 TRUE）
n_tcga_rows <- nrow(TD)
n_tcga_tasks <- length(unique(TD$task))
n_tcga_keys <- length(unique(TD$key))
n_mb_rows <- nrow(MBt)
n_mb_tasks <- length(unique(paste(MBt$signature, MBt$seed)))
n_mb_keys <- length(unique(paste(MBt$signature, MBt$seed, MBt$start)))
keys_tcga_ok <- setequal(unique(TD$task), unique(am$task))
keys_mb_ok <- setequal(unique(MBt$key), unique(amb$key))

CHK <- data.frame(
  item=c("p_sum_max_abs_dev","k1_eCv_equals_start","k1_eCt_equals_start","k1_gains_zero",
         "full_Cv_equals_argmax","full_Ct_equals_argmax",
         "metabric_full_Cv_equals_argmax","metabric_full_Ct_equals_argmax",
         "exhaustive_prob_max_abs_dev",
         "monotone_eCv_violations","monotone_tasks_checked",
         "n_task_rows","n_tasks","n_tasks_with_both_starts",
         "n_metabric_task_rows","n_metabric_tasks","n_metabric_with_both_starts",
         "full_keys_match_argmax","metabric_full_keys_match_argmax"),
  value=c(pdev, k1_dev_Cv, k1_dev_Ct, k1_dev_G, full_dev_Cv, full_dev_Ct,
          mb_full_dev_Cv, mb_full_dev_Ct, exh_max,
          mono_viol, mono_tot, n_tcga_rows, n_tcga_tasks, n_tcga_keys,
          n_mb_rows, n_mb_tasks, n_mb_keys,
          keys_tcga_ok, keys_mb_ok),
  pass=c(pdev < 1e-12, k1_dev_Cv < 1e-15, k1_dev_Ct < 1e-15, k1_dev_G < 1e-15,
         full_dev_Cv < 1e-10, full_dev_Ct < 1e-10, mb_full_dev_Cv < 1e-10,
         mb_full_dev_Ct < 1e-10, exh_max < 1e-12,
         mono_viol == 0, mono_tot == 7838L,
         n_tcga_rows == 47028L, n_tcga_tasks == 3919L, n_tcga_keys == 7838L,
         n_mb_rows == 360L, n_mb_tasks == 30L, n_mb_keys == 60L,
         isTRUE(keys_tcga_ok), isTRUE(keys_mb_ok)))
write.table(CHK, file.path(TAB,"50_search_budget_check.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

say(""); say("================ 核对 ================")
for (i in seq_len(nrow(CHK))) say(sprintf("  %-28s %s   %s", CHK$item[i],
  format(CHK$value[i], digits=6), if (isTRUE(CHK$pass[i])) "OK" else "FAIL"))
say(sprintf("  full 与旧 argmax：Cv max|dev|=%.2e，Ct max|dev|=%.2e",
            full_dev_Cv, full_dev_Ct))

## ============================================================
## D. 文本摘要
## ============================================================
say(""); say("================ TCGA（八队列等权）：起点 × K ================")
say("  ", sprintf("%-6s %5s %8s %8s %8s %8s %8s %8s", "start","K","eCv","eCt","Gv","Gt","Orel","randCt"))
for (i in seq_len(nrow(SUM))) say("  ", sprintf("%-6s %5s %8.4f %8.4f %+8.4f %+8.4f %+8.4f %8.4f",
  SUM$start[i], SUM$K_label[i], SUM$eCv[i], SUM$eCt[i],
  SUM$Gv[i], SUM$Gt[i], SUM$Orel[i], SUM$eRandCt[i]))
say("")
say("  六 OS 子集：")
for (i in seq_len(nrow(SUM6))) say("  ", sprintf("%-6s %5s %8.4f %8.4f %+8.4f %+8.4f %+8.4f",
  SUM6$start[i], SUM6$K_label[i], SUM6$eCv[i], SUM6$eCt[i],
  SUM6$Gv[i], SUM6$Gt[i], SUM6$Orel[i]))
say("")
say("================ METABRIC（六签名等权，ALL 行）================")
mball <- MB[MB$signature=="ALL", ]
for (i in seq_len(nrow(mball))) say("  ", sprintf("%-6s %5s %8.4f %8.4f %+8.4f %+8.4f %+8.4f",
  mball$start[i], mball$K_label[i], mball$eCv[i], mball$eCt[i],
  mball$Gv[i], mball$Gt[i], mball$Orel[i]))
say("")
say("说明：这是参考锚定搜索（起点始终在池中）；K=full 时两个起点回到同一完整池，")
say("      选中模型相同。曲线是【已见结果后】的修订性探索，不是事前预设；")
say("      测试表现没有单调保证，本文不据此选'最优 K'。")
if (any(!CHK$pass)) stop("50 核对未全部通过：", paste(CHK$item[!CHK$pass], collapse=", "))
close(con)
cat("50 完成。\n")
