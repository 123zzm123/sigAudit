#!/usr/bin/env Rscript
# ============================================================
# 36 — 从【逐患者预测】独立重算 METABRIC 的 C-index，并与提交的表核对
#
# 为什么需要（第七轮评审）：正文与复现说明声称结果可由逐患者预测重建，
# 但 32 读的是既有的 31 指标表。本脚本补上这条链：只读预测 RDS，用与主分析
# 【同一个】ci_direct / ens_pred 重新计算，再与 31 的两张指标表逐条对比。
#
# 【第八轮修正】旧版在匹配不到指标行时用 `if (is.na(i)) next` 静默跳过，
# 且把"实际检查了多少条"当成"应该检查多少条"，因此缺行也可能输出"全部通过"。
# 现在先按冻结设计构造【期望键集合】并核对齐全性，缺任何一条键就报出并停止。
#   期望：6 签名 × 5 种子 = 30 个任务
#         每任务 18 个候选 -> 540 条候选指标
#         每任务 7 条被核验规则 -> 210 条规则指标
#
# 【第九轮修正】上一版的期望集合仍由 PRED 自身派生（lapply(PRED, ...)），
# 若 RDS 本身少一个任务或某任务候选名单被缩减，期望会跟着缩水而仍然"通过"。
# 现在改为先核对 PRED 与冻结设计（6 签名 × 5 种子 × 18 模型）完全一致，
# 再由设计派生期望规模。实算条数改用循环计数器 k，不再用预分配向量的长度。
#   本项只服务当前冻结版本，不做旧版兼容或补值。
#
# 检查项：
#   1. 每个任务：样本 ID 唯一、训练/验证/测试互不重叠、并集覆盖全部纳入样本
#   2. 预测中的时间与事件与 METABRIC_surv.tsv 逐人一致
#   3. 训练、验证、测试预测均为有限值
#   4. 540 条候选指标与 31_metabric_model_metrics.tsv 一致
#   5. 210 条主要规则指标与 31_metabric_rule_performance.tsv 一致
#   6. 每个任务的 argmax 记录获胜者等于验证侧 C-index 最大者（30 条）
#
# 输出：results/tables/36_metabric_prediction_check.tsv
# 用法：Rscript scripts/36_verify_metabric_predictions.R
# ============================================================
suppressMessages(library(survival))
options(stringsAsFactors = FALSE)
PROJ <- Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
TAB  <- file.path(PROJ, "results/tables")
source(file.path(PROJ, "scripts/_lib_analysis.R"))   # ci_direct / ens_pred

F_PRED <- file.path(PROJ, "results/metabric/31_metabric_predictions.rds")
F_MM   <- file.path(TAB, "31_metabric_model_metrics.tsv")
F_RP   <- file.path(TAB, "31_metabric_rule_performance.tsv")
F_SURV <- file.path(PROJ, "data/cohorts/METABRIC_surv.tsv")
for (f in c(F_PRED, F_MM, F_RP, F_SURV)) if (!file.exists(f)) stop("缺少 ", f)

rd <- function(f) read.table(f, sep="\t", header=TRUE, quote="", comment.char="",
                             stringsAsFactors=FALSE)
PRED <- readRDS(F_PRED)
MM   <- rd(F_MM); RP <- rd(F_RP); SV <- rd(F_SURV)
sv_t <- setNames(as.numeric(SV$OS_MONTHS), SV$sampleId)
sv_e <- setNames(ifelse(SV$OS_STATUS == "1:DECEASED", 1L, 0L), SV$sampleId)

cat("预测文件含", length(PRED), "个任务\n")

## ============================================================
## 0a. 上游确认：PRED 必须与冻结设计逐项一致
##     （第九轮）期望规模由此派生，而不是由可能不完整的输入自行决定
## ============================================================
EXPECTED_SIGS <- c(
  "EndoPredict_11", "GENE70_MammaPrint", "GENE76_Rotterdam",
  "GGI", "OncotypeDX_21", "PAM50"
)
EXPECTED_SEEDS <- c(20261014, 20261016, 20261017, 20261018, 20261019)
EXPECTED_MODELS <- c(
  "Lasso", "Ridge", "Enet0.9", "Enet0.7", "Enet0.5",
  "Enet0.3", "Enet0.1", "Ridge_1se", "UniCox+Lasso",
  "LassoTop100", "LassoTop50", "LassoTop25", "LassoTop10",
  "Enet0.5Top25", "RSF100", "GBM100", "GBM250", "RSF250"
)
RULES <- c("none", "argmax_val", "top2_ensemble", "top3_ensemble",
           "top5_ensemble", "top10_ensemble", "all_ensemble")

expected_tasks <- with(
  expand.grid(signature = EXPECTED_SIGS, seed = EXPECTED_SEEDS,
              stringsAsFactors = FALSE),
  paste(signature, seed, sep = "|")
)
actual_tasks <- vapply(PRED, function(p)
  paste(p$signature, p$seed, sep = "|"), character(1))

design_ok <- identical(length(actual_tasks), length(expected_tasks)) &&
             !anyDuplicated(actual_tasks) &&
             setequal(actual_tasks, expected_tasks) &&
             all(vapply(PRED, function(p)
               length(p$model) == length(EXPECTED_MODELS) &&
                 !anyDuplicated(p$model) &&
                 setequal(p$model, EXPECTED_MODELS), logical(1)))

cat("\n=== 0a. 上游设计一致性 ===\n")
cat(sprintf("  任务：实际 %d / 期望 %d；目标集合一致 %s\n",
            length(actual_tasks), length(expected_tasks),
            if (setequal(actual_tasks, expected_tasks)) "是" else "否"))
if (!design_ok) {
  cat("  [中止] 预测 RDS 与冻结设计不一致：\n")
  cat("    缺少任务：\n"); print(head(setdiff(expected_tasks, actual_tasks), 20))
  cat("    多出任务：\n"); print(head(setdiff(actual_tasks, expected_tasks), 20))
  cat("    重复任务：\n"); print(head(actual_tasks[duplicated(actual_tasks)], 20))
  bad <- vapply(PRED, function(p)
    length(p$model) != length(EXPECTED_MODELS) ||
      anyDuplicated(p$model) || !setequal(p$model, EXPECTED_MODELS), logical(1))
  cat("    候选名单不符的任务：", sum(bad), "个\n")
  if (any(bad)) print(head(actual_tasks[bad], 20))
  cat(sprintf("  依据：%d 签名 × %d 种子 = %d 任务，每任务 %d 个候选，%d 条规则\n",
              length(EXPECTED_SIGS), length(EXPECTED_SEEDS),
              length(expected_tasks), length(EXPECTED_MODELS), length(RULES)))
  quit(status = 2)
}

## ============================================================
## 0b. 期望键集合的齐全性核对（缺项直接停止）
## ============================================================
MM$key <- paste(MM$signature, MM$seed, MM$model, sep="|")
RP$key <- paste(RP$signature, RP$seed, RP$rule, sep="|")

EXP_CAND <- unlist(lapply(strsplit(expected_tasks, "|", fixed = TRUE), function(t)
  paste(t[1], t[2], EXPECTED_MODELS, sep="|")))
EXP_RULE <- unlist(lapply(strsplit(expected_tasks, "|", fixed = TRUE), function(t)
  paste(t[1], t[2], RULES, sep="|")))

chk_keys <- function(exp, have, what, file) {
  dup <- unique(have[duplicated(have)])
  miss <- setdiff(exp, have)
  if (length(dup)) { cat("  [中止]", what, "存在重复键：\n"); print(head(dup, 10)); quit(status = 2) }
  if (length(miss)) {
    cat(sprintf("  [中止] %s 缺少 %d / %d 条期望指标（%s）\n", what, length(miss), length(exp), file))
    print(head(miss, 20)); quit(status = 2)
  }
  cat(sprintf("  %-22s 期望 %d 条，全部存在且唯一\n", what, length(exp)))
}
cat("\n=== 0. 期望键齐全性 ===\n")
chk_keys(EXP_CAND, MM$key, "候选指标", basename(F_MM))
chk_keys(EXP_RULE, RP$key, "规则指标", basename(F_RP))

## ============================================================
## 1-3. 结构检查
## ============================================================
struct <- do.call(rbind, lapply(PRED, function(p) {
  ids <- c(p$train_ids, p$val_ids, p$test_ids)
  data.frame(signature = p$signature, seed = p$seed,
             n_train = length(p$train_ids), n_val = length(p$val_ids),
             n_test = length(p$test_ids),
             unique_ids = length(unique(ids)) == length(ids),
             disjoint = length(intersect(p$train_ids, p$val_ids)) == 0 &&
                        length(intersect(p$train_ids, p$test_ids)) == 0 &&
                        length(intersect(p$val_ids,  p$test_ids)) == 0,
             covers_all = setequal(ids, names(sv_t)),
             time_match = all(p$t_train == sv_t[p$train_ids]) &&
                          all(p$t_val   == sv_t[p$val_ids]) &&
                          all(p$t_test  == sv_t[p$test_ids]),
             event_match = all(p$e_train == sv_e[p$train_ids]) &&
                           all(p$e_val   == sv_e[p$val_ids]) &&
                           all(p$e_test  == sv_e[p$test_ids]),
             ## 第八轮：训练/验证/测试三类预测都检查
             finite = all(is.finite(p$train_pred)) &&
                      all(is.finite(p$val_pred)) &&
                      all(is.finite(p$test_pred)),
             stringsAsFactors = FALSE)
}))
cat("\n=== 1-3. 结构检查 ===\n")
cat(sprintf("  ID 唯一 %d/%d；分区不重叠 %d/%d；覆盖全部 %d/%d；时间一致 %d/%d；事件一致 %d/%d；预测有限(含训练) %d/%d\n",
            sum(struct$unique_ids), nrow(struct), sum(struct$disjoint), nrow(struct),
            sum(struct$covers_all), nrow(struct), sum(struct$time_match), nrow(struct),
            sum(struct$event_match), nrow(struct), sum(struct$finite), nrow(struct)))

## ============================================================
## 4. 候选 C-index 重算（键已保证齐全）
## ============================================================
d_cand <- numeric(length(EXP_CAND))
k <- 0L
for (p in PRED) for (m in p$model) {
  k <- k + 1L
  i <- match(paste(p$signature, p$seed, m, sep="|"), MM$key)
  cv <- ci_direct(p$t_val,  p$e_val,  p$val_pred[, m])
  ct <- ci_direct(p$t_test, p$e_test, p$test_pred[, m])
  d_cand[k] <- max(abs(c(cv - MM$cv[i], ct - MM$ct[i])))
}
if (k != length(EXP_CAND))
  stop(sprintf("候选实算 %d 条，与期望 %d 条不符", k, length(EXP_CAND)))
cat(sprintf("\n=== 4. 候选指标 ===\n  重算 %d 条（期望 %d），最大绝对差 %.3e\n",
            k, length(EXP_CAND), max(d_cand)))

## ============================================================
## 5. 主要规则重算
## ============================================================
d_rule <- numeric(length(EXP_RULE)); n_rule <- 0L
for (p in PRED) {
  mods <- p$model
  V <- vapply(mods, function(m) ci_direct(p$t_val, p$e_val, p$val_pred[, m]), numeric(1))
  names(V) <- mods
  for (r in RULES) {
    n_rule <- n_rule + 1L
    i <- match(paste(p$signature, p$seed, r, sep="|"), RP$key)
    mem <- strsplit(RP$members[i], "|", fixed = TRUE)[[1]]
    if (r == "none") {
      cvt <- V[["Lasso"]]; ctt <- ci_direct(p$t_test, p$e_test, p$test_pred[, "Lasso"])
    } else if (grepl("ensemble", r)) {
      pr <- ens_pred(p$val_pred, p$test_pred, mem)
      cvt <- ci_direct(p$t_val,  p$e_val,  pr$v)
      ctt <- ci_direct(p$t_test, p$e_test, pr$t)
    } else {
      w <- RP$winner[i]
      cvt <- V[[w]]; ctt <- ci_direct(p$t_test, p$e_test, p$test_pred[, w])
    }
    d_rule[n_rule] <- max(abs(c(cvt - RP$cv[i], ctt - RP$ct[i])))
  }
}
cat(sprintf("\n=== 5. 规则指标 ===\n  重算 %d 条（期望 %d），最大绝对差 %.3e\n",
            n_rule, length(EXP_RULE), max(d_rule)))
cat("  说明：集成按记录成员重建，验证的是这些成员组成的预测与指标；\n")
cat("        不重新执行 majority_vote / worst_half / penalized 的随机破并与重采样选择。\n")

## ============================================================
## 6. argmax 获胜者核对（期望 = 任务数）
## ============================================================
win_ok <- vapply(PRED, function(p) {
  i <- match(paste(p$signature, p$seed, "argmax_val", sep="|"), RP$key)
  V <- vapply(p$model, function(m) ci_direct(p$t_val, p$e_val, p$val_pred[, m]), numeric(1))
  names(V) <- p$model
  abs(V[[RP$winner[i]]] - max(V)) < 1e-12
}, logical(1))
cat(sprintf("\n=== 6. argmax 获胜者 ===\n  %d / %d\n", sum(win_ok), length(win_ok)))

## ============================================================
## 写出（expected 用【期望集合】的规模，不是实际完成数）
## ============================================================
n_task <- length(expected_tasks)
out <- data.frame(
  item = c("任务数", "样本ID唯一", "分区不重叠", "覆盖全部纳入样本",
           "时间与生存表一致", "事件与生存表一致", "预测有限（含训练）",
           "候选指标期望条数", "候选指标实算条数", "候选C-index最大绝对差",
           "规则指标期望条数", "规则指标实算条数", "规则C-index最大绝对差",
           "argmax获胜者正确数"),
  value = c(n_task, sum(struct$unique_ids), sum(struct$disjoint), sum(struct$covers_all),
            sum(struct$time_match), sum(struct$event_match), sum(struct$finite),
            length(EXP_CAND), k, signif(max(d_cand), 6),
            length(EXP_RULE), n_rule, signif(max(d_rule), 6), sum(win_ok)),
  expected = c(n_task, n_task, n_task, n_task, n_task, n_task, n_task,
               length(EXP_CAND), length(EXP_CAND), 1e-10,
               length(EXP_RULE), length(EXP_RULE), 1e-10, n_task),
  stringsAsFactors = FALSE)
write.table(out, file.path(TAB, "36_metabric_prediction_check.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

ok <- all(struct$unique_ids) && all(struct$disjoint) && all(struct$covers_all) &&
      all(struct$time_match) && all(struct$event_match) && all(struct$finite) &&
      k == length(EXP_CAND) && max(d_cand) < 1e-10 &&
      n_rule == length(EXP_RULE) && max(d_rule) < 1e-10 &&
      sum(win_ok) == n_task
cat("\n", if (ok) "全部通过：METABRIC 的 C-index 可由随包预测独立重算。"
         else "**存在不一致，请检查**", "\n", sep = "")
cat("结果写入 results/tables/36_metabric_prediction_check.tsv\n")
if (!ok) quit(status = 1)
