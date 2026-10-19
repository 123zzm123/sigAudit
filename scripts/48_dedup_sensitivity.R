#!/usr/bin/env Rscript
# ============================================================
# 48 — 候选去重敏感性（本次修订新增，2026-10-19）
#
# 问题：候选池里存在大量重复（验证预测完全相同的候选对），这部分冗余是否
#      影响"搜索/集成相对固定参照"的实际结论？
#
# 设计（先冻结，再运行；不使用测试预测或测试结局做任何分组/选代表）：
#   · 判据：以【验证预测向量完全相同】（精确相等）为同一组；
#   · 选代表：同组内按候选名称【固定字典序】取最小的名字做代表；
#   · 参照：固定 Lasso / Ridge 的输出始终保留，即使它与其他候选重复；
#   · 新池重算：validation argmax 与 top-five（唯一候选不足 5 个时用 top-min(5,K)）；
#   · 原池与新池在【同一任务】上配对比较；不重训、不使用测试信息构造池。
#
# 口径与主分析一致：任务 -> (队列×划分) -> 队列 -> 八队列等权。
# ============================================================
suppressMessages({library(survival); library(glmnet); library(randomForestSRC); library(gbm)})
options(stringsAsFactors = FALSE)
PROJ <- Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
TAB  <- file.path(PROJ, "results/tables"); LOGD <- file.path(PROJ, "logs")
dir.create(LOGD, showWarnings = FALSE, recursive = TRUE)
con <- file(file.path(LOGD, "48_dedup_sensitivity.log"), "wt")
say <- function(...) { m <- paste0(...); cat(m,"\n",sep=""); cat(m,"\n",sep="",file=con); flush(con) }
rd <- function(f) read.table(f, sep="\t", header=TRUE, quote="", comment.char="", stringsAsFactors=FALSE)
.lib <- file.path(PROJ, "scripts/_lib_analysis.R"); if (!file.exists(.lib)) stop("缺少 _lib_analysis.R")
source(.lib)   # ci_direct / pick_max / order_rule / make_prio / ens_pred

# 候选指标（含 constant）与规则性能：用于对齐任务集合与提供 C 值
MM <- rd(file.path(TAB,"20_model_metrics_long.tsv"))
RP <- rd(file.path(TAB,"20_rule_performance_long.tsv"))
RPm <- rd(file.path(TAB,"20_rule_members_long.tsv"))
MM$task <- paste(MM$cohort,MM$seed,MM$geneset,sep="|")
RP$task <- paste(RP$cohort,RP$seed,RP$geneset,sep="|")
RPm$task <- paste(RPm$cohort,RPm$seed,RPm$geneset,sep="|")
# 查找表：task|model -> val_ci/test_ci
idx <- paste(MM$task, MM$model, sep="@@")
VAL <- setNames(MM$val_ci, idx); TST <- setNames(MM$test_ci, idx)
# 原始 argmax 获胜者（用于一致性核对）
am_win <- setNames(RPm$winner[RPm$rule=="argmax_val"], RPm$task[RPm$rule=="argmax_val"])

mani <- readLines(file.path(PROJ,"results/manifest/main_input_files.txt"), warn=FALSE); mani <- mani[nzchar(mani)]
say("读入 ", length(mani), " 个（队列×种子）产物")

PERF <- list(); FEAT <- list()
for (f in mani) {
  R <- tryCatch(readRDS(f), error=function(e) NULL)
  if (is.null(R) || !length(R$genesets)) next
  coh <- R$cohort; seed <- R$seed; G <- R$genesets
  cnt <- table(unlist(lapply(G, function(z) colnames(z$val_pred))))
  keep_mods <- names(cnt)[cnt/length(G) >= 0.90]
  if (!("Lasso" %in% keep_mods)) { say("[skip] ", coh, seed, " 参照缺失"); next }
  prio <- make_prio(keep_mods, seed)
  for (z in G) {
    if (!all(keep_mods %in% colnames(z$val_pred))) next
    Pv <- z$val_pred[, keep_mods, drop=FALSE]; Pt <- z$test_pred[, keep_mods, drop=FALSE]
    if (!all(is.finite(Pv)) || !all(is.finite(Pt))) next
    tk <- paste(coh, seed, z$name, sep="|")
    v <- VAL[paste(tk, keep_mods, sep="@@")]; t <- TST[paste(tk, keep_mods, sep="@@")]
    if (any(!is.finite(v)) || any(!is.finite(t))) next
    names(v) <- names(t) <- keep_mods
    # ---- 原始候选池（未做候选覆盖筛选）与筛选后分析池的对照 ----
    Pv_raw <- as.matrix(z$val_pred)
    n_cand_raw <- ncol(Pv_raw)
    n_unique_raw <- length(unique(apply(Pv_raw, 2, function(col) paste(format(col, digits=17), collapse=","))))
    # ---- 去重：按验证预测精确相等分组，组内按名字字典序取最小者为代表 ----
    key <- apply(Pv, 2, function(col) paste(format(col, digits=17), collapse=","))
    grp <- split(keep_mods, key)
    rep_names <- vapply(grp, function(g) sort(g)[1], character(1))
    dup_pool <- sort(unique(unname(rep_names)))
    K <- length(keep_mods); Ku <- length(dup_pool)
    n_dup_group <- sum(vapply(grp, length, integer(1)) > 1)
    rep_map <- paste(paste(rep_names, vapply(grp, function(g) paste(sort(g), collapse=","), character(1)), sep="<-"), collapse=";")
    # ---- 原池 argmax / top5（与主分析同一 prio 破并）----
    am_o <- keep_mods[pick_max(v[keep_mods], prio[keep_mods])]
    ord_o <- order_rule(v[keep_mods], prio[keep_mods])
    tk_o <- keep_mods[ord_o[seq_len(min(5, length(ord_o)))]]
    # ---- 去重池 argmax / top5（同一 prio 的对应子集，破并规则不变）----
    am_d <- dup_pool[pick_max(v[dup_pool], prio[dup_pool])]
    ord_d <- order_rule(v[dup_pool], prio[dup_pool])
    K5 <- min(5, Ku)
    tk_d <- dup_pool[ord_d[seq_len(K5)]]
    # ---- 集成（验证侧标准化 -> 取平均）----
    eo <- ens_pred(Pv, Pt, tk_o); ed <- ens_pred(Pv, Pt, tk_d)
    # ---- 固定参照 ----
    las_v <- v[["Lasso"]]; las_t <- t[["Lasso"]]
    has_ridge <- "Ridge" %in% keep_mods
    rg_v <- if (has_ridge) v[["Ridge"]] else NA_real_
    rg_t <- if (has_ridge) t[["Ridge"]] else NA_real_

    mk <- function(rule, cvv, ctt) data.frame(cohort=coh, seed=seed, geneset=z$name, rule=rule,
        Cv=cvv, Ct=ctt, Gv=cvv-las_v, Gt=ctt-las_t, Orel=(cvv-las_v)-(ctt-las_t),
        dCt_ridge=ctt-rg_t, stringsAsFactors=FALSE)
    PERF[[length(PERF)+1]] <- rbind(
      mk("argmax_orig",  v[[am_o]], t[[am_o]]),
      mk("argmax_dedup", v[[am_d]], t[[am_d]]),
      mk("top5_orig",  ci_direct(z$tva, z$eva, eo$v), ci_direct(z$tte, z$ete, eo$t)),
      mk("top5_dedup", ci_direct(z$tva, z$eva, ed$v), ci_direct(z$tte, z$ete, ed$t))
    )
    FEAT[[length(FEAT)+1]] <- data.frame(cohort=coh, seed=seed, geneset=z$name,
      n_cand=K, n_unique=Ku, n_dup_group=n_dup_group,
      n_cand_raw=n_cand_raw, n_unique_raw=n_unique_raw,
      argmax_changed=as.integer(!identical(am_o, am_d)),
      am_o=am_o, am_d=am_d,
      top5_changed=as.integer(!identical(sort(tk_o), sort(tk_d))),
      n_top5_members_changed=length(setdiff(tk_o, tk_d)),
      K5=K5, mem_o=paste(tk_o, collapse="|"), mem_d=paste(tk_d, collapse="|"),
      rep_map=rep_map, stringsAsFactors=FALSE)
  }
}
T <- do.call(rbind, FEAT)
Perf <- do.call(rbind, PERF)
rownames(Perf) <- NULL; rownames(T) <- NULL
T$task <- paste(T$cohort, T$seed, T$geneset, sep="|")
Perf$task <- paste(Perf$cohort, Perf$seed, Perf$geneset, sep="|")

## ---- 一致性核对：本脚本重算的 argmax_orig 是否等于主分析 20 的获胜者 ----
mine <- setNames(T$am_o, T$task)
chk <- data.frame(task=T$task, mine=mine[T$task], main=am_win[T$task], stringsAsFactors=FALSE)
agree <- sum(chk$mine == chk$main, na.rm=TRUE); tot <- sum(!is.na(chk$main))
say("一致性核对：本脚本原池 argmax 与主分析 20 获胜者一致 ", agree, "/", tot)

## ---- 三层级汇总 ----
MC <- c("Cv","Ct","Gv","Gt","Orel","dCt_ridge")
strat <- function(df) {
  s1 <- aggregate(df[MC], by=df[c("cohort","seed","rule")], FUN=mean, na.rm=TRUE)
  s2 <- aggregate(s1[MC], by=s1[c("cohort","rule")], FUN=mean, na.rm=TRUE)
  aggregate(s2[MC], by=s2["rule"], FUN=mean, na.rm=TRUE)
}
SUM <- strat(Perf)
write.table(SUM, file.path(TAB,"48_dedup_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## 队列级
s1c <- aggregate(Perf[MC], by=Perf[c("cohort","seed","rule")], FUN=mean, na.rm=TRUE)
COH <- aggregate(s1c[MC], by=s1c[c("cohort","rule")], FUN=mean, na.rm=TRUE)
write.table(COH, file.path(TAB,"48_dedup_by_cohort.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## 覆盖/成员变化（队列级：计数用求和，均值用任务级平均，不能套同一个 mean）
cnt_by_cohort <- aggregate(list(n_tasks=rep(1L, nrow(T))),
                           by=list(cohort=T$cohort), FUN=sum)
means_by_cohort <- aggregate(list(
                       mean_n_cand=T$n_cand, mean_n_unique=T$n_unique,
                       mean_n_cand_raw=T$n_cand_raw, mean_n_unique_raw=T$n_unique_raw,
                       mean_dup_group=T$n_dup_group,
                       frac_any_dup=as.integer(T$n_unique < T$n_cand),
                       frac_any_dup_raw=as.integer(T$n_unique_raw < T$n_cand_raw),
                       frac_argmax_changed=T$argmax_changed,
                       frac_top5_changed=T$top5_changed,
                       mean_top5_members_changed=T$n_top5_members_changed,
                       frac_K_lt5=as.integer(T$n_unique<5)),
                  by=list(cohort=T$cohort), FUN=mean)
cov2 <- merge(cnt_by_cohort, means_by_cohort, by="cohort", sort=FALSE)
cov2 <- cov2[match(unique(T$cohort), cov2$cohort), ]
write.table(cov2, file.path(TAB,"48_dedup_coverage.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(T, file.path(TAB,"48_dedup_per_task.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
## 逐任务性能（原池 vs 去重池），供独立复算，不只是诊断标记
write.table(Perf, file.path(TAB,"48_dedup_performance_long.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## ---- 文本摘要 ----
say(""); say("================ 原池 vs 去重池（三层级，八队列等权）================")
say("  ", sprintf("%-14s %8s %8s %9s %9s %9s %11s", "rule","Cv","Ct","Gv","Gt","Orel","Ct−Ridge"))
for (i in seq_len(nrow(SUM))) say("  ", sprintf("%-14s %8.4f %8.4f %+9.4f %+9.4f %+9.4f %+11.4f",
  SUM$rule[i], SUM$Cv[i], SUM$Ct[i], SUM$Gv[i], SUM$Gt[i], SUM$Orel[i], SUM$dCt_ridge[i]))
say("")
n_tot <- nrow(T)
say(sprintf("  任务总数 %d；分析池候选 %.4f（唯一 %.4f）；原始池候选 %.4f（唯一 %.4f）",
    n_tot, mean(T$n_cand), mean(T$n_unique), mean(T$n_cand_raw), mean(T$n_unique_raw)))
say(sprintf("  原始池比分析池恰好多一个候选的任务 %d（%.1f%%）",
    sum(T$n_cand_raw == T$n_cand + 1L), 100*mean(T$n_cand_raw == T$n_cand + 1L)))
say(sprintf("  含至少一对重复候选的任务：分析池 %d/%d = %.3f%%；原始池 %.3f%%",
    sum(T$n_unique < T$n_cand), n_tot, 100*mean(T$n_unique < T$n_cand),
    100*mean(T$n_unique_raw < T$n_cand_raw)))
say(sprintf("  平均重复组数 %.2f；唯一候选不足 5 个的任务 %d/%d = %.3f%%",
    mean(T$n_dup_group), sum(T$n_unique<5), n_tot, 100*sum(T$n_unique<5)/n_tot))
say(sprintf("  argmax 获胜者改变：%d/%d = %.3f%%（按任务计数）",
    sum(T$argmax_changed), n_tot, 100*sum(T$argmax_changed)/n_tot))
say(sprintf("  top5 成员集合改变：%d/%d = %.3f%%（按任务计数；队列等权宏平均 %.2f%%）",
    sum(T$top5_changed), n_tot, 100*sum(T$top5_changed)/n_tot, 100*mean(cov2$frac_top5_changed)))
say("")
say("  关键对照：")
g <- function(r) SUM[SUM$rule==r, ]
say(sprintf("    argmax 相对 Lasso：原池 Orel=%+.4f，去重池 Orel=%+.4f",
            g("argmax_orig")$Orel, g("argmax_dedup")$Orel))
say(sprintf("    top5  相对 Lasso：原池 Orel=%+.4f，去重池 Orel=%+.4f",
            g("top5_orig")$Orel, g("top5_dedup")$Orel))
say(sprintf("    top5−argmax：原池 %+.4f，去重池 %+.4f",
            g("top5_orig")$Ct - g("argmax_orig")$Ct, g("top5_dedup")$Ct - g("argmax_dedup")$Ct))
say(sprintf("    argmax−Ridge：原池 %+.4f，去重池 %+.4f",
            g("argmax_orig")$dCt_ridge, g("argmax_dedup")$dCt_ridge))
say("")
say("说明：去重只删除验证预测完全相同的候选，不改变参照输出；")
say("      若结论改变，原样报告，不选择性呈现。")
close(con)
cat("48 完成。\n")
