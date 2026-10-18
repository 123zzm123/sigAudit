#!/usr/bin/env Rscript
# ============================================================
# 20 — 修正版分析（第二轮评审意见落地）
#
# 主要比较（预先指定）：
#   P1  argmax_val 相对 none（固定 Lasso 参照）的 Gv / Gt / Orel
#   P2  top5_ensemble 相对 argmax_val 的配对差（推断单位 = 队列）
#
# 相对上一版的实质修正：
#  【B3】argmax_val / 集成 / runner_up **允许 Lasso 参与**（原版 setdiff 掉参照，
#        实际回答的是"强制选一个非参照"的收益）。原版口径保留为
#        argmax_nonreference，仅作敏感性分析。
#  【B4】先生成【唯一的任务可评估清单】：同一（队列,划分,基因集）上，
#        规则与参照的全部所需指标在【同一配对集合】上计算；
#        不再分别对 Gv、Gt 取均值（那会在缺失分布不同时失去配对含义）。
#        候选预测出现非有限值 -> 该任务整体记为失败，不静默缩小患者集。
#  【B2】P2 的推断单位改为【队列】：先把每队列各划分的配对差取平均，
#        得到队列级效应再等权汇总；不把共享患者的重复划分当独立样本。
#  【I6】输出每队列跨划分的 Gv/Gt/Orel 全部取值与范围。
#  【I7】评价修正的影响分解：报告两种口径下"选模是否改变"的比例。
#  【I8】可比对生存对数取 survival::concordance 的 count 字段之和，不用自写判据。
#  【I3】落盘长表 model_metrics_long / rule_metrics_long（含获胜者、成员）。
#  【并列】用与候选输入顺序无关的固定随机优先级破并（pick_max / order_rule）。
#
# 输入：results/manifest/main_input_files.txt（由 23_manifest.R 生成，不靠目录扫描）
# ============================================================
suppressMessages(library(survival))
options(stringsAsFactors = FALSE)
PROJ <- Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit"); TAB <- file.path(PROJ, "results/tables")
dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
TREES <- c("RSF100","GBM100","GBM250","RSF250")
# 重采样次数（只影响三条【次要】规则：majority_vote / worst_half / penalized）。
# 20/10 -> 12/6 是【运行成本控制】，主比较 P1/P2 完全不受影响；
# 该改动记录于此，并在 README_current 中如实说明。
B_VOTE <- 12; B_HALF <- 6

LOGF <- file.path(PROJ, "logs", "20_analysis_corrected.log"); con <- file(LOGF, "wt")
say <- function(...) { m <- paste0(...); cat(m,"\n",sep=""); cat(m,"\n",sep="",file=con); flush(con) }

## ---------- 规则与评价函数：全部来自共用库 ----------
## 之所以改成 source()：第三轮审查指出 20 脚本与 R 包对同名规则口径不一致
## （linear_only 是否含参照）、n_pairs 算错。把函数集中在 _lib_analysis.R，
## 让 20（主分析）与 29（固定基线敏感性）共用，杜绝脚本之间的口径漂移。
.lib <- file.path(PROJ, "scripts/_lib_analysis.R")
if (!file.exists(.lib)) stop("缺少共用库 scripts/_lib_analysis.R")
source(.lib)   # 提供 ci_direct / ci_old / n_pairs（已修）/ pick_max / order_rule /
               # make_prio / ens_pred / rules_one / strat_summary / cohort_paired

## ---------- 读入（只读 manifest 列出的文件）----------
mani_f <- file.path(PROJ, "results/manifest/main_input_files.txt")
if (!file.exists(mani_f)) stop("缺少 manifest，请先跑 scripts/23_manifest.R")
fs <- readLines(mani_f, warn = FALSE); fs <- fs[nzchar(fs)]
say("按 manifest 读入 ", length(fs), " 个（队列×种子）产物")

ROW <- list(); COV <- list(); FIX <- list(); P2 <- list()
MML <- list(); RML <- list(); RPL <- list()   # RML=成员清单；RPL=逐任务性能

for (f in fs) {
  R <- tryCatch(readRDS(f), error=function(e) NULL)
  if (is.null(R) || !length(R$genesets)) { say("[skip] ", basename(f)); next }
  coh <- R$cohort; seed <- R$seed; G <- R$genesets
  cnt <- table(unlist(lapply(G, function(z) colnames(z$val_pred))))
  keep_mods <- names(cnt)[cnt/length(G) >= 0.90]
  if (!("Lasso" %in% keep_mods)) { say("[skip] ", coh, " ", seed, " 参照缺失"); next }
  # 规范优先级：只依赖 (候选名字集合, seed)，与输入列顺序无关
  prio0 <- make_prio(keep_mods, seed)

  ## 【B4】唯一任务可评估清单
  tasks <- list(); pm_long <- list()
  for (z in G) {
    if (!all(keep_mods %in% colnames(z$val_pred))) next
    Pv <- z$val_pred[, keep_mods, drop=FALSE]; Pt <- z$test_pred[, keep_mods, drop=FALSE]
    if (!all(is.finite(Pv)) || !all(is.finite(Pt))) next      # 非有限值 -> 任务失败
    V <- vapply(keep_mods, function(m) ci_direct(z$tva, z$eva, Pv[,m]), numeric(1))
    Tt<- vapply(keep_mods, function(m) ci_direct(z$tte, z$ete, Pt[,m]), numeric(1))
    if (!all(is.finite(V)) || !all(is.finite(Tt))) next
    tasks[[length(tasks)+1]] <- list(name=z$name, val_pred=Pv, test_pred=Pt,
      tva=z$tva, eva=z$eva, tte=z$tte, ete=z$ete, V=V, Tt=Tt,
      t_tr=z$t_tr, e_tr=z$e_tr, train_pred=z$train_pred[, keep_mods, drop=FALSE])
    pm_long[[length(pm_long)+1]] <- data.frame(cohort=coh, seed=seed, geneset=z$name,
      model=keep_mods,
      train_ci = vapply(keep_mods, function(m) ci_direct(z$t_tr, z$e_tr, z$train_pred[,m]), numeric(1)),
      val_ci = V, test_ci = Tt,
      constant = vapply(keep_mods, function(m) length(unique(Pv[,m])) == 1, logical(1)))
  }
  if (length(tasks) < 30) { say("[skip] ", coh, " seed", seed, " 可评估任务 <30"); next }
  MML[[paste(coh,seed)]] <- do.call(rbind, pm_long)
  COV[[paste(coh,seed)]] <- data.frame(cohort=coh, seed=seed, endpoint=R$endpoint,
    study_id=R$study_id, n_planned=R$n_planned, n_processed=R$n_processed,
    n_fitted=length(G), n_tasks=length(tasks), n_candidates=length(keep_mods),
    n_train=length(R$train_ids), n_val=length(R$val_ids), n_test=length(R$test_ids),
    ev_train=sum(tasks[[1]]$e_tr), ev_val=sum(tasks[[1]]$eva), ev_test=sum(tasks[[1]]$ete),
    pairs_val=n_pairs(tasks[[1]]$tva, tasks[[1]]$eva),
    pairs_test=n_pairs(tasks[[1]]$tte, tasks[[1]]$ete))
  say(""); say("===== ", coh, " seed ", seed, " (", R$endpoint, ") : ", length(tasks),
      " 个可评估任务 × ", length(keep_mods), " 候选 =====")

  A <- list(); RL <- list(); same <- 0; tot <- 0
  eval_fix <- numeric(0); fixed_old_eval <- numeric(0)
  for (z2 in tasks) {
    rn <- rules_one(z2, prio0, ci_direct, B_VOTE = B_VOTE, B_HALF = B_HALF)
    for (k in names(rn$out))
      A[[length(A)+1]] <- data.frame(geneset=z2$name, rule=k, cv=rn$out[[k]][1], ct=rn$out[[k]][2])
    for (k in names(rn$rec)) RL[[length(RL)+1]] <- cbind(cohort=coh, seed=seed, geneset=z2$name, rn$rec[[k]])
    # 逐任务规则性能（第三轮审查要求：规则长表必须含每任务 Cv/Ct 与参照，而不是只有成员列）
    RPL[[length(RPL)+1]] <- data.frame(cohort=coh, seed=seed, geneset=z2$name,
      rule=names(rn$out),
      cv=vapply(rn$out, function(x) x[1], numeric(1)),
      ct=vapply(rn$out, function(x) x[2], numeric(1)),
      stringsAsFactors=FALSE)
    # 【I7】只对 argmax_val 做口径对照，避免整池跑两遍（原实现每个任务要跑约 1300 次 C-index）
    mods <- colnames(z2$val_pred)
    Vo <- vapply(mods, function(m) ci_old(z2$tva, z2$eva, z2$val_pred[,m]), numeric(1)); names(Vo) <- mods
    wn <- rn$rec$argmax_val$winner
    wo <- if (all(!is.finite(Vo))) NA_character_ else mods[pick_max(Vo, prio0[mods])]
    tot <- tot + 1
    if (identical(wn, wo)) same <- same + 1
    # 固定【新口径选出的获胜者】，仅替换测试评价器 -> 分离出纯评价器效应
    t_new <- ci_direct(z2$tte, z2$ete, z2$test_pred[, wn])
    t_old <- ci_old   (z2$tte, z2$ete, z2$test_pred[, wn])
    eval_fix <- c(eval_fix, t_old - t_new)
    # 参照侧同样处理，便于构成 Orel 的口径差
    r_n <- ci_direct(z2$tte, z2$ete, z2$test_pred[, "Lasso"])
    r_o <- ci_old   (z2$tte, z2$ete, z2$test_pred[, "Lasso"])
    fixed_old_eval <- c(fixed_old_eval, (t_old - r_o) - (t_new - r_n))
  }
  A <- do.call(rbind, A)
  RML[[paste(coh,seed)]] <- do.call(rbind, RL)

  ## 【B4】同一配对集合上先算 Gv/Gt/Orel
  b <- A[A$rule=="none", ]; ii <- match(A$geneset, b$geneset)
  A$Gv <- A$cv - b$cv[ii]; A$Gt <- A$ct - b$ct[ii]; A$Orel <- A$Gv - A$Gt
  ok <- is.finite(A$Gv) & is.finite(A$Gt) & is.finite(A$Orel)
  a <- do.call(rbind, lapply(split(A[ok,], A$rule[ok]), function(d)
    data.frame(cohort=coh, seed=seed, rule=d$rule[1], n_tasks=nrow(d),
               Cv=mean(d$cv), Ct=mean(d$ct), Gv=mean(d$Gv), Gt=mean(d$Gt),
               Orel=mean(d$Orel), Oraw=mean(d$cv - d$ct))))
  ROW[[paste(coh,seed)]] <- a

  ## P2：相对 argmax_val 的配对差（同一任务集）
  rf <- A[A$rule=="argmax_val", ]; ii2 <- match(A$geneset, rf$geneset)
  dm <- is.finite(A$ct) & is.finite(rf$ct[ii2]) & is.finite(A$cv) & is.finite(rf$cv[ii2])
  A$dCt <- A$ct - rf$ct[ii2]; A$dOrel <- (A$cv - A$ct) - (rf$cv[ii2] - rf$ct[ii2])
  P2[[paste(coh,seed)]] <- do.call(rbind, lapply(split(A[dm,], A$rule[dm]), function(d)
    data.frame(cohort=coh, seed=seed, rule=d$rule[1], n_tasks=nrow(d),
               dCt=mean(d$dCt), dOrel=mean(d$dOrel))))
  FIX[[paste(coh,seed)]] <- data.frame(cohort=coh, seed=seed,
    pct_selection_changed = 100*(1 - same/tot), n=tot,
    dCt_evaluator_only = mean(eval_fix, na.rm=TRUE),          # 固定获胜者，仅换测试评价器
    dOrel_evaluator_only = mean(fixed_old_eval, na.rm=TRUE))  # 同上，且参照侧同口径处理
  say(sprintf("  【I7】两种评价口径下 argmax_val 选中同一候选的任务比例: %.1f%%", 100*same/tot))
}

if (!length(ROW)) stop("没有可用的（队列×划分）结果")
R <- do.call(rbind, ROW)
write.table(R, file.path(TAB,"20_rules.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(do.call(rbind, COV), file.path(TAB,"20_coverage.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
if (length(MML)) write.table(do.call(rbind, MML), file.path(TAB,"20_model_metrics_long.tsv"),
                             sep="\t", quote=FALSE, row.names=FALSE)
if (length(RML)) write.table(do.call(rbind, RML), file.path(TAB,"20_rule_members_long.tsv"),
                             sep="\t", quote=FALSE, row.names=FALSE)
if (length(RPL)) {
  RP <- do.call(rbind, RPL)
  # 在同一配对集合上给出相对参照（none = 固定 Lasso）的 Gv/Gt/Orel
  b <- RP[RP$rule == "none", ]; ii <- match(paste(RP$cohort,RP$seed,RP$geneset),
                                            paste(b$cohort,b$seed,b$geneset))
  RP$ref_cv <- b$cv[ii]; RP$ref_ct <- b$ct[ii]
  RP$Gv <- RP$cv - RP$ref_cv; RP$Gt <- RP$ct - RP$ref_ct; RP$Orel <- RP$Gv - RP$Gt
  RP$Oraw <- RP$cv - RP$ct
  RP$status <- ifelse(is.finite(RP$cv) & is.finite(RP$ct), "ok", "invalid")
  write.table(RP, file.path(TAB,"20_rule_performance_long.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)
}
if (length(FIX)) write.table(do.call(rbind, FIX), file.path(TAB,"20_eval_fix.tsv"),
                             sep="\t", quote=FALSE, row.names=FALSE)
if (length(P2))  write.table(do.call(rbind, P2), file.path(TAB,"20_paired_vs_argmax.tsv"),
                             sep="\t", quote=FALSE, row.names=FALSE)

## ---------- 汇总：先队列内跨划分平均 -> 再队列等权 ----------
say(""); say("============ 汇总：先队列内跨划分平均，再对队列等权 ============")
byc <- aggregate(cbind(Cv,Ct,Gv,Gt,Orel,Oraw) ~ cohort + rule, data=R, FUN=mean, na.rm=TRUE)
agg <- aggregate(cbind(Cv,Ct,Gv,Gt,Orel,Oraw) ~ rule, data=byc, FUN=mean, na.rm=TRUE)
rng <- aggregate(Orel ~ rule, data=byc, FUN=function(x) diff(range(x, na.rm=TRUE)))
agg$Orel_range <- rng$Orel[match(agg$rule, rng$rule)]
agg <- agg[order(-agg$Ct), ]
say("  ", sprintf("%-22s %8s %8s %9s %9s %9s %10s", "rule","Cv","Ct","Gv","Gt","Orel","跨队列范围"))
for (i in seq_len(nrow(agg))) say("  ", sprintf("%-22s %8.4f %8.4f %+9.4f %+9.4f %+9.4f %10.4f",
  agg$rule[i], agg$Cv[i], agg$Ct[i], agg$Gv[i], agg$Gt[i], agg$Orel[i], agg$Orel_range[i]))
r1 <- agg[agg$rule=="argmax_val",]
say("")
say(sprintf("  ★ P1：argmax_val（全池，含参照） Gv=%+.4f  Gt=%+.4f  Orel=%+.4f",
            r1$Gv, r1$Gt, r1$Orel))
if (is.finite(r1$Gv) && r1$Gv > 0.02) say(sprintf("     比例：未在独立测试集保留的部分占验证增益 %.1f%%（仅在 Gv 明显为正时解释）",
                                                  100*r1$Orel/r1$Gv))
say("     绝对增益差为主信息；比例随基线改变而变化，不作为唯一卖点。")
write.table(agg, file.path(TAB,"20_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## ---------- 【I6】跨划分稳定性 ----------
say(""); say("============ 【I6】每队列跨划分稳定性（argmax_val）============")
ST <- do.call(rbind, lapply(split(R, R$cohort), function(d) {
  s <- d[d$rule=="argmax_val", ]
  data.frame(cohort=d$cohort[1], n_seed=nrow(s),
             Orel_mean=mean(s$Orel), Orel_min=min(s$Orel), Orel_max=max(s$Orel),
             Gt_mean=mean(s$Gt), Gt_min=min(s$Gt), Gt_max=max(s$Gt),
             values=paste(sprintf("%+.4f", s$Orel), collapse=" ")) }))
for (i in seq_len(nrow(ST))) say(sprintf("  %-8s %d 个划分  Orel 各值: %s   均值 %+.4f  极差 %.4f",
  ST$cohort[i], ST$n_seed[i], ST$values[i], ST$Orel_mean[i], ST$Orel_max[i]-ST$Orel_min[i]))
write.table(ST, file.path(TAB,"20_cross_split_stability.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## ---------- 【B2】P2：推断单位 = 队列 ----------
if (length(P2)) {
  P <- do.call(rbind, P2)
  coh_lv <- aggregate(cbind(dCt,dOrel) ~ cohort + rule, data=P, FUN=mean, na.rm=TRUE)
  u <- aggregate(cbind(dCt,dOrel) ~ rule, data=coh_lv, FUN=mean, na.rm=TRUE)
  u <- u[order(-u$dCt), ]
  say(""); say("============ 【B2】各规则 vs argmax_val（配对；推断单位 = 队列）============")
  say("  先队列内跨划分平均，再对队列等权。差异不显著 != 等效或非劣。")
  say("  ", sprintf("%-22s %10s %10s %7s %20s %s", "rule","dCt","dOrel","n队列","dCt 95%CI","p"))
  for (i in seq_len(nrow(u))) {
    x <- coh_lv$dCt[coh_lv$rule==u$rule[i]]; x <- x[is.finite(x)]
    tt <- if (length(x) >= 4) t.test(x) else NULL
    say("  ", sprintf("%-22s %+10.4f %+10.4f %7d %20s %s", u$rule[i], u$dCt[i], u$dOrel[i], length(x),
        if (is.null(tt)) "n/a" else sprintf("[%+.4f,%+.4f]", tt$conf.int[1], tt$conf.int[2]),
        if (is.null(tt)) "" else sprintf("p=%.2e", tt$p.value)))
  }
  say("  注：以癌种队列为重复单位的探索性基准汇总，不是对全部癌症总体的代表性推断。")
}
if (length(FIX)) {
  F <- do.call(rbind, FIX)
  say(""); say("============ 【I7】评价修正的影响 ============")
  say(sprintf("  ① 选模是否会变：两种口径下 argmax_val 选中同一候选的任务比例 %.1f%%（即改变 %.1f%%）",
              100 - mean(F$pct_selection_changed), mean(F$pct_selection_changed)))
  say(sprintf("  ② 固定获胜者后【仅替换测试评价器】：ΔCt = %+.4f，ΔOrel = %+.4f",
              mean(F$dCt_evaluator_only), mean(F$dOrel_evaluator_only)))
  say("  因此新旧口径的总差异 = 「选模变化」+「评价器变化」两部分；")
  say("  上一版把总差异当作纯评价器效应，现按评审要求分解为 ①②。")
}
say(""); say("规模与流程链：")
CV <- do.call(rbind, COV)
say("  ", sprintf("%-6s %6s %6s %6s %6s %6s %6s %6s", "队列","种子","计划","处理","拟合","任务","候选","测试n"))
for (i in seq_len(nrow(CV))) say("  ", sprintf("%-6s %6d %6d %6d %6d %6d %6d %6d",
  CV$cohort[i], CV$seed[i], CV$n_planned[i], CV$n_processed[i], CV$n_fitted[i],
  CV$n_tasks[i], CV$n_candidates[i], CV$n_test[i]))
say(sprintf("  合计：计划 %d，处理 %d，拟合成功 %d，可评估任务 %d",
            sum(CV$n_planned), sum(CV$n_processed), sum(CV$n_fitted), sum(CV$n_tasks)))
close(con)
