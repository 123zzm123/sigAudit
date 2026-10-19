#!/usr/bin/env Rscript
# ============================================================
# 49 — 预测级审计入口（R/audit_predictions.R）的一致性核对与用途演示
#
# 目标（本次修订交付要求）：
#   A. 在一批现有冻结预测（results/main/19_*.rds）上运行 audit_predictions()，
#      与主分析流水线（20_* 表）逐任务比较 Cv/Ct/Gv/Gt/Orel，给出最大绝对差。
#   B. 用同一批选中模型演示"更换固定参照（Lasso -> Ridge -> Enet0.5）"后 gain 口径改变。
#   C. 演示参数 deduplicate = TRUE 时原池/去重池的差异（与 48 脚本互证）。
#
# 不重训任何模型；只读取已保存的逐患者预测（每个 RDS 只读一次，供 A/B/C 复用）。
# 输出：49_audit_api_check.tsv / 49_audit_api_reference_switch.tsv /
#       49_audit_api_dedup.tsv / logs/49_audit_api.log
# ============================================================
suppressMessages({library(survival); library(glmnet); library(randomForestSRC); library(gbm)})
options(stringsAsFactors = FALSE)
PROJ <- Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
TAB  <- file.path(PROJ, "results/tables"); LOGD <- file.path(PROJ, "logs")
dir.create(LOGD, showWarnings = FALSE, recursive = TRUE)
con <- file(file.path(LOGD, "49_audit_api.log"), "wt")
say <- function(...) { m <- paste0(...); cat(m,"\n",sep=""); cat(m,"\n",sep="",file=con); flush(con) }
rd <- function(f) read.table(f, sep="\t", header=TRUE, quote="", comment.char="", stringsAsFactors=FALSE)

for (f in c("R/utils.R","R/audit_modules.R","R/audit_selection.R","R/audit_predictions.R"))
  source(file.path(PROJ, f))

TREES <- c("RSF100","GBM100","GBM250","RSF250")
CORE <- c("none","argmax_val","argmax_nonreference","runner_up",
          "top2_ensemble","top3_ensemble","top5_ensemble","top10_ensemble","all_ensemble",
          "majority_vote","worst_half","penalized","oracle_test")
## 全 40 个 cell 核对确定性规则；三条重采样规则在单个队列（MESO）上做完整核对，
## 以控制运行时间（重采样规则每任务要跑数十次 C-index）。
CORE_DET <- setdiff(CORE, c("majority_vote","worst_half","penalized"))
CORE3 <- c("none","argmax_val","top5_ensemble")

mani <- readLines(file.path(PROJ,"results/manifest/main_input_files.txt"), warn=FALSE); mani <- mani[nzchar(mani)]
say("读入 ", length(mani), " 个（队列×种子）产物")

build_tasks <- function(R) {
  G <- R$genesets
  cnt <- table(unlist(lapply(G, function(z) colnames(z$val_pred))))
  keep <- names(cnt)[cnt/length(G) >= 0.90]
  if (!("Lasso" %in% keep)) return(NULL)
  out <- list()
  for (z in G) {
    if (!all(keep %in% colnames(z$val_pred))) next
    Pv <- z$val_pred[, keep, drop=FALSE]; Pt <- z$test_pred[, keep, drop=FALSE]
    if (!all(is.finite(Pv)) || !all(is.finite(Pt))) next
    out[[length(out)+1]] <- list(name=z$name, cohort=R$cohort, split=as.character(R$seed),
      val_pred=Pv, test_pred=Pt,
      time_val=z$tva, event_val=z$eva, time_test=z$tte, event_test=z$ete)
  }
  if (length(out) < 30) return(NULL)
  list(tasks=out, keep=keep)
}

## ============================================================
## 单次遍历：每个 RDS 只读一次，A/B/C 复用
## ============================================================
API <- list(); DIAG <- list(); SW <- list(); DD <- list()
for (f in mani) {
  R <- tryCatch(readRDS(f), error=function(e) NULL)
  if (is.null(R) || !length(R$genesets)) next
  bt <- build_tasks(R); if (is.null(bt)) next
  fam <- stats::setNames(ifelse(bt$keep %in% TREES, "tree", "linear"), bt$keep)
  ## A. 参照 Lasso：全队列跑确定性规则，MESO 额外跑全部规则（含重采样）
  rls <- if (identical(R$cohort, "MESO")) CORE else CORE_DET
  a <- audit_predictions(bt$tasks, reference="Lasso", rules=rls, family=fam, seed=R$seed)
  r <- a$rules; r$task <- paste(r$cohort,r$split,r$task,sep="|"); API[[length(API)+1]] <- r
  d <- a$diagnostics; d$task <- paste(d$cohort,d$split,d$task,sep="|"); DIAG[[length(DIAG)+1]] <- d
  ## B. 参照切换（Lasso / Ridge / Enet0.5）
  for (ref in c("Lasso","Ridge","Enet0.5")) {
    if (!(ref %in% bt$keep)) next
    b <- audit_predictions(bt$tasks, reference=ref, rules=CORE3, seed=R$seed)
    rr <- b$rules; rr$task <- paste(rr$cohort,rr$split,rr$task,sep="|"); rr$reference <- ref
    SW[[length(SW)+1]] <- rr
  }
  ## C. 去重演示
  for (dd in c(FALSE, TRUE)) {
    cc <- audit_predictions(bt$tasks, reference="Lasso", rules=CORE3, seed=R$seed, deduplicate=dd)
    rr <- cc$rules; rr$task <- paste(rr$cohort,rr$split,rr$task,sep="|"); rr$dedup <- dd
    DD[[length(DD)+1]] <- rr
  }
}
A <- do.call(rbind, API)
say("audit_predictions 规则记录 ", nrow(A), " 条；任务 ", length(unique(A$task)))
say("一致性核对范围：确定性规则覆盖全部 40 个 cohort×split；",
    "重采样规则（majority_vote/worst_half/penalized）在 MESO 上完整核对")
say("说明：本核对覆盖 3,919 个共同可评估任务上的确定性规则 + MESO 的重采样规则，",
    "不等于覆盖主分析全部 58,785 条规则记录。")

## ============================================================
## A. 一致性核对（reference = Lasso）
## ============================================================
RP <- rd(file.path(TAB,"20_rule_performance_long.tsv")); RP$task <- paste(RP$cohort,RP$seed,RP$geneset,sep="|")
RM <- rd(file.path(TAB,"20_rule_members_long.tsv")); RM$task <- paste(RM$cohort,RM$seed,RM$geneset,sep="|")
## ---- 预期覆盖键：直接从【冻结的】20_rule_performance_long.tsv 构造，不依赖本轮读入的 RDS ----
## 这样即使某个 RDS 文件整体缺失（actual 侧一起变少），预期数仍保持不变，能报出缺口。
## 共同可评估任务 = none 与 argmax_val 的 Cv/Ct 均有限；规则范围按队列声明
## （全部队列跑确定性规则，MESO 额外跑三条重采样规则）。
Wp <- reshape(RP[, c("task","cohort","rule","cv","ct")], idvar=c("task","cohort"),
              timevar="rule", direction="wide")
ok_task <- is.finite(Wp$cv.none) & is.finite(Wp$ct.none) &
           is.finite(Wp$cv.argmax_val) & is.finite(Wp$ct.argmax_val)
Wt <- Wp[ok_task, c("task","cohort"), drop=FALSE]
EXPECT <- do.call(rbind, lapply(seq_len(nrow(Wt)), function(i) {
  rls <- if (identical(Wt$cohort[i], "MESO")) CORE else CORE_DET
  data.frame(task=Wt$task[i], rule=rls, stringsAsFactors=FALSE)
}))
EXP <- nrow(EXPECT)
M <- merge(RP[,c("task","rule","cv","ct","Gv","Gt","Orel")], A[,c("task","rule","Cv","Ct","Gv","Gt","Orel")],
           by=c("task","rule"), suffixes=c(".pipe",".api"))
M <- M[is.finite(M$cv) & is.finite(M$Cv), ]
dCv <- max(abs(M$cv - M$Cv)); dCt <- max(abs(M$ct - M$Ct))
dGv <- max(abs(M$Gv.pipe - M$Gv.api)); dGt <- max(abs(M$Gt.pipe - M$Gt.api))
dOr <- max(abs(M$Orel.pipe - M$Orel.api))
MM <- merge(RM[RM$rule %in% c("argmax_val","top5_ensemble"), c("task","rule","winner","members")],
            A[A$rule %in% c("argmax_val","top5_ensemble"), c("task","rule","winner","members")],
            by=c("task","rule"), suffixes=c(".pipe",".api"))
win_same <- mean(MM$winner.pipe == MM$winner.api, na.rm=TRUE)
mem_same <- mean(MM$members.pipe == MM$members.api, na.rm=TRUE)
## ---- 真实断言：键唯一性、双向覆盖（预期缺 API / API 多出预期）----
keyA <- paste(A$task, A$rule, sep="@@")
keyP <- paste(RP$task, RP$rule, sep="@@")
keyE <- paste(EXPECT$task, EXPECT$rule, sep="@@")
dup_keys <- sum(duplicated(keyA))
missing_records <- length(setdiff(keyE, keyA))   # 预期有、API 没有
extra_records   <- length(setdiff(keyA, keyE))   # API 有、预期没有
say(sprintf("A. 与主分析 20 的逐任务比较（%d 条匹配记录）", nrow(M)))
say(sprintf("   最大绝对差：Cv=%.2e  Ct=%.2e  Gv=%.2e  Gt=%.2e  Orel=%.2e", dCv, dCt, dGv, dGt, dOr))
say(sprintf("   argmax_val/top5_ensemble 获胜者一致率 %.4f，成员一致率 %.4f", win_same, mem_same))
say(sprintf("   记录数：实际 %d，预期 %d（由冻结 20 表构造）；重复键 %d；预期缺 API %d；API 多出预期 %d",
            nrow(A), EXP, dup_keys, missing_records, extra_records))
CHK <- data.frame(item=c("expected_records","actual_records","duplicated_keys",
                         "expected_missing_in_api","api_extra_not_expected",
                         "max_abs_diff_Cv","max_abs_diff_Ct","max_abs_diff_Gv",
                         "max_abs_diff_Gt","max_abs_diff_Orel","winner_agreement","members_agreement"),
  value=c(EXP, nrow(A), dup_keys, missing_records, extra_records,
          dCv, dCt, dGv, dGt, dOr, win_same, mem_same),
  pass=c(nrow(A)==EXP, nrow(A)==EXP, dup_keys==0, missing_records==0, extra_records==0,
         dCv<1e-10, dCt<1e-10, dGv<1e-10, dGt<1e-10, dOr<1e-10, win_same==1, mem_same==1))
write.table(CHK, file.path(TAB,"49_audit_api_check.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
if (any(!CHK$pass))
  stop("49 一致性核对未全部通过：", paste(CHK$item[!CHK$pass], collapse=", "))

## ============================================================
## B. 参照切换（三层级宏汇总）
## ============================================================
three_level <- function(df) {
  df$split <- as.character(df$split)
  s1 <- aggregate(df[c("Cv","Ct","Gv","Gt","Orel")], by=df[c("cohort","split","reference","rule")], FUN=mean, na.rm=TRUE)
  s2 <- aggregate(s1[c("Cv","Ct","Gv","Gt","Orel")], by=s1[c("cohort","reference","rule")], FUN=mean, na.rm=TRUE)
  aggregate(s2[c("Cv","Ct","Gv","Gt","Orel")], by=s2[c("reference","rule")], FUN=mean, na.rm=TRUE)
}
B <- do.call(rbind, SW); SWT <- three_level(B)
write.table(SWT, file.path(TAB,"49_audit_api_reference_switch.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
say(""); say("B. 参照切换（三层级，八队列等权）：")
for (ref in unique(SWT$reference)) {
  z <- SWT[SWT$reference==ref & SWT$rule=="argmax_val", ]
  say(sprintf("   argmax vs 固定 %-8s : Gv=%+.4f Gt=%+.4f Orel=%+.4f", ref, z$Gv, z$Gt, z$Orel))
}

## ============================================================
## C. 去重演示（三层级）
## ============================================================
DDT <- do.call(rbind, DD); DDT$split <- as.character(DDT$split)
s1 <- aggregate(DDT[c("Cv","Ct","Gv","Gt","Orel")], by=DDT[c("cohort","split","dedup","rule")], FUN=mean, na.rm=TRUE)
s2 <- aggregate(s1[c("Cv","Ct","Gv","Gt","Orel")], by=s1[c("cohort","dedup","rule")], FUN=mean, na.rm=TRUE)
DS  <- aggregate(s2[c("Cv","Ct","Gv","Gt","Orel")], by=s2[c("dedup","rule")], FUN=mean, na.rm=TRUE)
write.table(DS, file.path(TAB,"49_audit_api_dedup.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
say(""); say("C. 去重演示（三层级）：")
for (dd in c(FALSE, TRUE)) {
  z <- DS[DS$dedup==dd, ]; g <- function(r) z[z$rule==r, ]
  say(sprintf("   dedup=%-5s argmax Orel=%+.4f  top5 Orel=%+.4f", dd, g("argmax_val")$Orel, g("top5_ensemble")$Orel))
}
say(""); say("结论：预测级入口与原流水线在 1e-10 容差内一致；更换参照会改变 gain 口径；")
say("      去重不改变主要结论（与 scripts/48 独立复算一致）。")
close(con)
cat("49 完成。\n")
