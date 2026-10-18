#!/usr/bin/env Rscript
# ============================================================
# 29 — 固定基线敏感性分析（第三轮审查后新增；第四轮审查后修正汇总层级）
#
# 为什么必须有：审查在**完全相同的 3919 个任务**上，把规则成员表与模型长表
# 连接、保留 argmax 原有获胜者、仅把核算参照从固定 Lasso 换成固定 Ridge，
# 发现 **固定 Ridge 的测试 C 高于 argmax，也高于 top-5 集成**。
#
# ── 第四轮审查指出的口径错误（本版已修）────────────────────
#  旧版在汇总规则与固定候选时**丢掉了 seed 层级**：直接从任务级按 cohort 平均，
#  再对队列做宏平均。而主分析要求 任务 -> (队列×划分) -> 队列 -> 队列等权。
#  由于各（队列×划分）的完整任务数不同（KIRC 1365-1485、LGG 1410-1455，
#  其余恒为 1485），两种层级给出**不同的权重**，因而结果**不是精度差异**而是口径差异。
#  实测：固定 Ridge 测试 C 在错口径下 0.6143268，正确层级下 **0.6142744**（差 5.2e-05）。
#  另外旧版把 `n_tasks` 写成 nrow(J)（= 3919×15 = 58785 的合计），
#  在每行重复，被当成单规则样本量 —— 本版改为**该规则在该参照下的实际任务数**。
#
# 本脚本**不重训任何模型**，只对已保存的逐任务指标做重新核算。
# 参照 = 固定候选（single-model）。所有固定候选都会被汇总，不挑有利基线。
#
# 输入：20_rule_performance_long.tsv / 20_model_metrics_long.tsv / 20_rule_members_long.tsv
# 输出：29_fixed_candidate_summary.tsv / 29_baseline_sensitivity.tsv /
#       29_paired_vs_fixed.tsv / 29_per_cohort_vs_ridge.tsv / 29_verification.tsv
# ============================================================
suppressMessages(library(survival))
options(stringsAsFactors = FALSE)
PROJ <- Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
TAB  <- file.path(PROJ, "results/tables")
LOGF <- file.path(PROJ, "logs", "29_baseline_sensitivity.log")
dir.create(dirname(LOGF), showWarnings = FALSE, recursive = TRUE)
con <- file(LOGF, "wt")
say <- function(...) { m <- paste0(...); cat(m,"\n",sep=""); cat(m,"\n",sep="",file=con); flush(con) }
rd <- function(f) read.table(f, sep="\t", header=TRUE, quote="", comment.char="", stringsAsFactors=FALSE)

F_RP <- file.path(TAB, "20_rule_performance_long.tsv")
F_MM <- file.path(TAB, "20_model_metrics_long.tsv")
F_RM <- file.path(TAB, "20_rule_members_long.tsv")
for (f in c(F_RP, F_MM, F_RM)) if (!file.exists(f)) stop("缺少输入: ", f)

RP <- rd(F_RP); MM <- rd(F_MM)
RP$task <- paste(RP$cohort, RP$seed, RP$geneset, sep="|")
MM$task <- paste(MM$cohort, MM$seed, MM$geneset, sep="|")
say("读入规则性能记录 ", nrow(RP), " 条（任务 ", length(unique(RP$task)), " 个 × 规则 ",
    length(unique(RP$rule)), " 条）")
say("读入模型指标记录 ", nrow(MM), " 条（候选 ", length(unique(MM$model)), " 个）")

## ============================================================
## 汇总层级（与主分析 20 完全一致）
##   任务 -> (队列×划分) 平均 -> 队列平均 -> 队列等权
## ============================================================
KEY1 <- c("cohort","seed")     # 第一层
KEY2 <- "cohort"               # 第二层
strat_rule <- function(df, cols, rule_col = "rule", extra = character(0)) {
  by1 <- c(KEY1, rule_col, extra); by2 <- c(KEY2, rule_col, extra)
  s1 <- aggregate(df[cols], by = df[by1], FUN = mean, na.rm = TRUE)
  s2 <- aggregate(s1[cols], by = s1[by2], FUN = mean, na.rm = TRUE)
  aggregate(s2[cols], by = s2[c(rule_col, extra)], FUN = mean, na.rm = TRUE)
}
## 中间层（队列×划分），配对差与 t 区间要用
strat_seed <- function(df, cols, rule_col = "rule", extra = character(0)) {
  aggregate(df[cols], by = df[c(KEY1, rule_col, extra)], FUN = mean, na.rm = TRUE)
}
strat_cohort <- function(s1, cols, rule_col = "rule", extra = character(0)) {
  aggregate(s1[cols], by = s1[c(KEY2, rule_col, extra)], FUN = mean, na.rm = TRUE)
}

## ============================================================
## A. 所有固定候选的汇总性能（同样用三层级；不挑有利基线）
## ============================================================
cand_rows <- do.call(rbind, lapply(split(MM, MM$model), function(d) {
  d <- d[is.finite(d$val_ci) & is.finite(d$test_ci), ]
  if (!nrow(d)) return(NULL)
  s1 <- aggregate(cbind(val_ci, test_ci, constant) ~ cohort + seed, data = d, FUN = mean, na.rm = TRUE)
  s2 <- aggregate(cbind(val_ci, test_ci, constant) ~ cohort, data = s1, FUN = mean, na.rm = TRUE)
  data.frame(model = d$model[1], n_tasks = nrow(d), n_cohort = nrow(s2),
             Cv = mean(s2$val_ci), Ct = mean(s2$test_ci),
             Cv_min = min(s2$val_ci), Cv_max = max(s2$val_ci),
             Ct_min = min(s2$test_ci), Ct_max = max(s2$test_ci),
             constant_rate_val = mean(s1$constant, na.rm = TRUE))
}))
cand_rows <- cand_rows[order(-cand_rows$Ct), ]
write.table(cand_rows, file.path(TAB, "29_fixed_candidate_summary.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
say("")
say("============ A. 所有固定候选的汇总性能（三层级，队列等权）============")
say("  ", sprintf("%-14s %7s %8s %8s %9s %8s", "候选","n任务","验证C","测试C","验证常数率","跨队列Ct范围"))
for (i in seq_len(nrow(cand_rows))) say("  ", sprintf("%-14s %7d %8.4f %8.4f %8.2f%% %5.4f~%.4f",
  cand_rows$model[i], cand_rows$n_tasks[i], cand_rows$Cv[i], cand_rows$Ct[i],
  100*cand_rows$constant_rate_val[i], cand_rows$Ct_min[i], cand_rows$Ct_max[i]))
say("  ⚠️ 只有 n任务 = 3919 的候选之间才可直接比较；任务数更少的候选（如 UniCox+Lasso）")
say("     是在**不同的任务子集**上评价的，其数值偏高不能解释为'更好'。")

say("")
say("  各队列的验证预测常数率：")
cs <- aggregate(constant ~ cohort + model, data = MM, FUN = mean, na.rm = TRUE)
for (m in c("Lasso","Ridge","Enet0.5")) {
  z <- cs[cs$model == m, ]; z <- z[order(-z$constant), ]
  say(sprintf("    %-8s %s", m, paste(sprintf("%s %.1f%%", z$cohort, 100*z$constant), collapse="  ")))
}

## ============================================================
## B/C. 参照敏感性
## ============================================================
mmv <- MM[, c("task","model","val_ci","test_ci")]
rule_tab <- RP[, c("task","cohort","seed","rule","cv","ct","status")]

refs <- c("Lasso", "Ridge", "Enet0.5")
BS <- list(); PA <- list(); Jall <- list()
for (rf in refs) {
  R <- mmv[mmv$model == rf & is.finite(mmv$val_ci) & is.finite(mmv$test_ci), ]
  if (!nrow(R)) { say("[跳过] 参照 ", rf, " 无可用任务"); next }
  names(R)[names(R) == "val_ci"] <- "ref_cv"; names(R)[names(R) == "test_ci"] <- "ref_ct"
  ## 配对完整任务集：该规则与参照在该任务上都可用
  J <- merge(rule_tab, R[, c("task","ref_cv","ref_ct")], by = "task")
  J <- J[is.finite(J$cv) & is.finite(J$ct), ]
  J$Gv <- J$cv - J$ref_cv; J$Gt <- J$ct - J$ref_ct; J$Orel <- J$Gv - J$Gt
  J$dCt <- J$ct - J$ref_ct; J$dOrel <- (J$cv - J$ct) - (J$ref_cv - J$ref_ct)
  J$reference <- rf
  Jall[[rf]] <- J

  ag <- strat_rule(J, c("cv","ct","Gv","Gt","Orel"))
  ## 【修正】每行的 n_tasks = 该规则在该参照下的实际配对数（不再重复写成合计）
  ntask <- aggregate(list(n_tasks = J$cv), by = list(rule = J$rule), FUN = function(x)
    sum(is.finite(x)))
  ag$n_tasks <- ntask$n_tasks[match(ag$rule, ntask$rule)]
  ag$reference <- rf
  BS[[rf]] <- ag

  ## 配对差（推断单位 = 队列；先队列内跨划分平均，再队列等权）
  s1 <- strat_seed(J, c("dCt","dOrel","Gv","Gt","Orel"))
  s2 <- strat_cohort(s1, c("dCt","dOrel","Gv","Gt","Orel"))
  pa <- do.call(rbind, lapply(split(s2, s2$rule), function(d) {
    x <- d$dCt[is.finite(d$dCt)]
    tt <- if (length(x) >= 4) tryCatch(t.test(x), error=function(e) NULL) else NULL
    data.frame(reference = rf, rule = d$rule[1], n_cohort = length(x),
               dCt = mean(x),
               lo = if (is.null(tt)) NA_real_ else tt$conf.int[1],
               hi = if (is.null(tt)) NA_real_ else tt$conf.int[2],
               p  = if (is.null(tt)) NA_real_ else tt$p.value,
               dOrel = mean(d$dOrel, na.rm = TRUE))
  }))
  PA[[rf]] <- pa
}
BS <- do.call(rbind, BS); PA <- do.call(rbind, PA)
BS <- BS[order(BS$reference, -BS$ct), ]
write.table(BS, file.path(TAB, "29_baseline_sensitivity.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
PA <- PA[order(PA$reference, -PA$dCt), ]
write.table(PA, file.path(TAB, "29_paired_vs_fixed.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

say("")
say("============ B. 每个（参照 × 规则）的核算结果（三层级）============")
for (rf in unique(BS$reference)) {
  z <- BS[BS$reference == rf, ]
  say(sprintf("  --- 参照 = 固定 %s ---", rf))
  say("  ", sprintf("%-22s %8s %8s %8s %9s %9s %9s", "规则","n任务","Cv","Ct","Gv","Gt","Orel"))
  for (i in seq_len(nrow(z))) say("  ", sprintf("%-22s %8d %8.4f %8.4f %+9.4f %+9.4f %+9.4f",
    z$rule[i], z$n_tasks[i], z$cv[i], z$ct[i], z$Gv[i], z$Gt[i], z$Orel[i]))
}

say("")
say("============ C. 配对差（推断单位 = 队列）============")
for (rf in unique(PA$reference)) {
  z <- PA[PA$reference == rf, ]
  say(sprintf("  --- 相对固定 %s 的 ΔCt ---", rf))
  for (i in seq_len(nrow(z))) say("  ", sprintf("%-22s %+9.4f [%+.4f,%+.4f] p=%-9s  n队列=%d",
    z$rule[i], z$dCt[i], z$lo[i], z$hi[i],
    if (is.na(z$p[i])) "n/a" else sprintf("%.2e", z$p[i]), z$n_cohort[i]))
}

## ============================================================
## C2. 逐队列：argmax / top5 相对固定 Ridge（三层级）
## ============================================================
Jr <- Jall[["Ridge"]]
Jr <- Jr[Jr$rule %in% c("argmax_val","top5_ensemble"), c("cohort","seed","rule","dCt")]
s1 <- aggregate(dCt ~ cohort + seed + rule, data = Jr, FUN = mean, na.rm = TRUE)
CL <- aggregate(dCt ~ cohort + rule, data = s1, FUN = mean, na.rm = TRUE)
say("")
say("============ C2. 逐队列：相对固定 Ridge 的平均测试差 ============")
say("  ", sprintf("%-8s %12s %12s", "队列", "argmax−Ridge", "top5−Ridge"))
cohs <- sort(unique(CL$cohort))
for (cc in cohs) {
  a <- CL$dCt[CL$cohort == cc & CL$rule == "argmax_val"]
  b <- CL$dCt[CL$cohort == cc & CL$rule == "top5_ensemble"]
  say("  ", sprintf("%-8s %+12.4f %+12.4f", cc, a, b))
}
say(sprintf("  低于 Ridge 的队列数：argmax %d/%d，top5 %d/%d",
            sum(CL$dCt[CL$rule=="argmax_val"] < 0), length(cohs),
            sum(CL$dCt[CL$rule=="top5_ensemble"] < 0), length(cohs)))
write.table(CL, file.path(TAB, "29_per_cohort_vs_ridge.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## ============================================================
## D. 与审查报告数字逐项对照
## ============================================================
getd <- function(rf, rl) PA$dCt[PA$reference == rf & PA$rule == rl]
getb <- function(rf, rl, col) BS[[col]][BS$reference == rf & BS$rule == rl]
exp_tab <- data.frame(
  item = c("固定Ridge测试C",
           "argmax − Ridge", "top5 − Ridge", "Lasso − Ridge",
           "argmax|Ridge: Gv","argmax|Ridge: Gt","argmax|Ridge: Orel",
           "argmax|Enet0.5: Gv","argmax|Enet0.5: Gt","argmax|Enet0.5: Orel"),
  reported = c(0.6142743734, -0.006382, -0.001575, -0.030039,
               0.031366, -0.006382, 0.037748,
               0.056677, 0.018840, 0.037836),
  computed = c(cand_rows$Ct[cand_rows$model == "Ridge"],
               getd("Ridge","argmax_val"), getd("Ridge","top5_ensemble"), getd("Ridge","none"),
               getb("Ridge","argmax_val","Gv"), getb("Ridge","argmax_val","Gt"), getb("Ridge","argmax_val","Orel"),
               getb("Enet0.5","argmax_val","Gv"), getb("Enet0.5","argmax_val","Gt"), getb("Enet0.5","argmax_val","Orel")),
  stringsAsFactors = FALSE)
exp_tab$diff <- exp_tab$computed - exp_tab$reported
exp_tab$pass <- abs(exp_tab$diff) < 1e-6        # 【修正】从 1e-4 收紧到 1e-6
write.table(exp_tab, file.path(TAB, "29_verification.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
say("")
say("============ D. 与审查报告数字逐项对照（容差 1e-6）============")
say("  ", sprintf("%-26s %14s %14s %12s %s", "项目","审查报告","本次复算","差值","判定"))
for (i in seq_len(nrow(exp_tab))) say("  ", sprintf("%-26s %14.8f %14.8f %+12.2e %s",
  exp_tab$item[i], exp_tab$reported[i], exp_tab$computed[i], exp_tab$diff[i],
  if (exp_tab$pass[i]) "一致" else "**不一致**"))
say("")
say("核对结果: ", sum(exp_tab$pass), "/", nrow(exp_tab), " 项在 1e-6 容差内一致")
say("")
say("★ 结论（供正文使用）：")
say("  · 绝对 Orel 在 Lasso / Ridge / Enet0.5 三种参照下都约为 0.037–0.038 —— 稳健。")
say("  · 但 '未兑现比例' = Orel/Gv 随参照显著变化；Ridge 参照下 Gt 为负，比例超过 100%。")
say("  · 因此正文以【绝对增益 + 参照依赖性】为主，61.2% 降为解释性次要数字。")
say("  · 表述纪律：Lasso 的零系数解是 **L1 稀疏性**的表现，不能称为'弱正则化'的证据；")
say("    与 Ridge 的对比应写成'稀疏 Lasso 基线'对'稠密 Ridge 比较模型'。")
close(con)
