#!/usr/bin/env Rscript
# ============================================================
# 47 — 非恒定 LASSO 分层诊断（本次修订新增，2026-10-19）
#
# 目的：回答"论文的主要发现是否只是由 LASSO 在验证集上给出恒定（全零系数）
#      解所驱动"。以 argmax_val 的 3,919 个共同可评估任务为基准，按
#      【LASSO 的验证预测是否恒定】分层，不按测试表现筛选。
#
# 口径（严格与主分析 20/29 一致）：
#   任务 -> (队列×划分) 平均 -> 队列内跨划分平均 -> 八队列等权平均。
#   恒定判据 = 验证预测分数是否只有一个唯一值（沿用 20_model_metrics_long 的 constant 列）。
#   LASSO 参照 = 固定 Lasso（规则 none）；dense 参照 = 固定 Ridge。
#   非恒定子集是【事后条件分层】，不是新的主分析集，也不是因果证明。
#
# 输入：results/tables/20_rule_performance_long.tsv, 20_model_metrics_long.tsv,
#       及（可选）results/metabric/31_metabric_predictions.rds + 31/32 表
# 输出：47_nonconstant_overall.tsv / _by_cohort.tsv / _by_cell.tsv /
#        _os6.tsv / _metabric.tsv / _summary.txt
# ============================================================
suppressMessages({library(survival); library(glmnet); library(randomForestSRC); library(gbm)})
options(stringsAsFactors = FALSE)
PROJ <- Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
TAB  <- file.path(PROJ, "results/tables")
LOGD <- file.path(PROJ, "logs"); dir.create(LOGD, showWarnings = FALSE, recursive = TRUE)
con <- file(file.path(LOGD, "47_nonconstant_diagnostic.log"), "wt")
say <- function(...) { m <- paste0(...); cat(m,"\n",sep=""); cat(m,"\n",sep="",file=con); flush(con) }
rd <- function(f) read.table(f, sep="\t", header=TRUE, quote="", comment.char="", stringsAsFactors=FALSE)
F1 <- file.path(TAB,"20_rule_performance_long.tsv"); F2 <- file.path(TAB,"20_model_metrics_long.tsv")
for (f in c(F1,F2)) if (!file.exists(f)) stop("缺少输入: ", f)
RP <- rd(F1); MM <- rd(F2)
RP$task <- paste(RP$cohort, RP$seed, RP$geneset, sep="|")
MM$task <- paste(MM$cohort, MM$seed, MM$geneset, sep="|")

## ---------- 组装任务级基准表 ----------
am <- RP[RP$rule=="argmax_val", c("task","cohort","seed","cv","ct")]
names(am)[4:5] <- c("am_cv","am_ct")
t5 <- RP[RP$rule=="top5_ensemble", c("task","cv","ct")]
names(t5)[2:3] <- c("t5_cv","t5_ct")
ls <- RP[RP$rule=="none", c("task","cv","ct")]
names(ls)[2:3] <- c("lasso_cv","lasso_ct")
mmL <- MM[MM$model=="Lasso", c("task","constant")]
names(mmL)[2] <- "lasso_constant"
mmR <- MM[MM$model=="Ridge", c("task","val_ci","test_ci")]
names(mmR)[2:3] <- c("ridge_cv","ridge_ct")

D <- Reduce(function(a,b) merge(a,b,by="task"), list(am,t5,ls,mmL,mmR))
D$gv   <- D$am_cv - D$lasso_cv
D$gt   <- D$am_ct - D$lasso_ct
D$orel <- D$gv - D$gt
D$am_minus_ridge   <- D$am_ct - D$ridge_ct
D$t5_minus_am      <- D$t5_ct - D$am_ct
D$t5_minus_ridge   <- D$t5_ct - D$ridge_ct
D$gv5  <- D$t5_cv - D$lasso_cv
D$gt5  <- D$t5_ct - D$lasso_ct
say("共同可评估任务（argmax 基准）: ", nrow(D))
say("  其中 LASSO 验证预测恒定: ", sum(D$lasso_constant), "；非恒定: ", sum(!D$lasso_constant))

## ---------- 三层级汇总（任务 -> 队列×划分 -> 队列 -> 队列等权）----------
MC <- c("gv","gt","orel","lasso_ct","am_ct","ridge_ct","t5_ct",
        "am_minus_ridge","t5_minus_am","t5_minus_ridge","gv5","gt5")
agg3 <- function(df) {
  s1 <- aggregate(df[MC], by=df[c("cohort","seed")], FUN=mean, na.rm=TRUE)
  s2 <- aggregate(s1[MC], by=s1["cohort"], FUN=mean, na.rm=TRUE)
  colMeans(s2[MC])
}
agg_cohort <- function(df) {
  s1 <- aggregate(df[MC], by=df[c("cohort","seed")], FUN=mean, na.rm=TRUE)
  aggregate(s1[MC], by=s1["cohort"], FUN=mean, na.rm=TRUE)
}
mk_row <- function(tag, df) {
  a <- agg3(df)
  data.frame(stratum=tag, n_tasks=nrow(df), n_cohort=length(unique(df$cohort)),
             Gv=a[["gv"]], Gt=a[["gt"]], Orel=a[["orel"]], Orel_over_Gv=a[["orel"]]/a[["gv"]],
             Lasso_Ct=a[["lasso_ct"]], argmax_Ct=a[["am_ct"]], Ridge_Ct=a[["ridge_ct"]],
             top5_Ct=a[["t5_ct"]], argmax_minus_Ridge=a[["am_minus_ridge"]],
             top5_minus_argmax=a[["t5_minus_am"]], top5_minus_Ridge=a[["t5_minus_ridge"]],
             Gv_top5=a[["gv5"]], Gt_top5=a[["gt5"]], stringsAsFactors=FALSE)
}
OV <- rbind(mk_row("all_tasks", D), mk_row("nonconstant_Lasso", D[!D$lasso_constant, ]))
write.table(OV, file.path(TAB,"47_nonconstant_overall.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## ---------- 队列级（all / nonconstant 并列）----------
CL <- rbind(cbind(stratum="all_tasks", agg_cohort(D)),
            cbind(stratum="nonconstant_Lasso", agg_cohort(D[!D$lasso_constant, ])))
cnt <- aggregate(list(n_tasks=D$gv), by=list(cohort=D$cohort), FUN=function(x) sum(is.finite(x)))
cc  <- aggregate(list(n_constant=D$lasso_constant), by=list(cohort=D$cohort), FUN=sum)
cnt <- merge(cnt, cc, by="cohort")
cnt$n_nonconstant <- cnt$n_tasks - cnt$n_constant
CL <- merge(CL, cnt[,c("cohort","n_tasks","n_constant","n_nonconstant")], by="cohort", all.x=TRUE)
write.table(CL, file.path(TAB,"47_nonconstant_by_cohort.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## ---------- 队列×划分覆盖数与均值 ----------
s1_all <- aggregate(D[MC], by=D[c("cohort","seed")], FUN=mean, na.rm=TRUE)
n_all  <- aggregate(list(n_tasks=D$gv), by=D[c("cohort","seed")], FUN=function(x) sum(is.finite(x)))
n_con  <- aggregate(list(n_constant=D$lasso_constant), by=D[c("cohort","seed")], FUN=sum)
CI <- Reduce(function(a,b) merge(a,b,by=c("cohort","seed")), list(s1_all,n_all,n_con))
CI$n_nonconstant <- CI$n_tasks - CI$n_constant
write.table(CI, file.path(TAB,"47_nonconstant_by_cell.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## ---------- 六 OS 队列 ----------
OS6 <- c("BRCA","COAD","KIRC","LUAD","MESO","STAD")
D6 <- D[D$cohort %in% OS6, ]
OV6 <- rbind(mk_row("all_tasks_OS6", D6), mk_row("nonconstant_Lasso_OS6", D6[!D6$lasso_constant, ]))
write.table(OV6, file.path(TAB,"47_nonconstant_os6.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## ---------- 逐队列 argmax−Ridge 方向 ----------
dir_tab <- do.call(rbind, lapply(split(CL, CL$stratum), function(g)
  data.frame(stratum=g$stratum[1],
             n_cohort_negative=sum(g$argmax_minus_ridge < 0, na.rm=TRUE),
             n_cohort=length(unique(g$cohort)))))

## ============================================================
## METABRIC 描述（若有预测：由验证预测直接判定 Lasso 是否恒定）
## ============================================================
MB <- NULL
prf <- file.path(PROJ,"results/metabric/31_metabric_predictions.rds")
mmf <- file.path(TAB,"31_metabric_model_metrics.tsv")
rpf <- file.path(TAB,"31_metabric_rule_performance.tsv")
if (file.exists(prf) && file.exists(mmf) && file.exists(rpf)) {
  P <- readRDS(prf)
  cm <- do.call(rbind, lapply(P, function(z)
    data.frame(signature=z$signature, seed=z$seed,
               lasso_constant=("Lasso" %in% colnames(z$val_pred)) &&
                 length(unique(z$val_pred[,"Lasso"]))==1)))
  mm <- rd(mmf); colnames(mm)[colnames(mm)=="cv"] <- "ridge_cv"; colnames(mm)[colnames(mm)=="ct"] <- "ridge_ct"
  mmR <- mm[mm$model=="Ridge", c("signature","seed","ridge_cv","ridge_ct")]
  rp <- rd(rpf)
  am <- rp[rp$rule=="argmax_val", c("signature","seed","cv","ct")]; names(am)[3:4] <- c("am_cv","am_ct")
  t5 <- rp[rp$rule=="top5_ensemble", c("signature","seed","cv","ct")]; names(t5)[3:4] <- c("t5_cv","t5_ct")
  ls <- rp[rp$rule=="none", c("signature","seed","cv","ct")]; names(ls)[3:4] <- c("lasso_cv","lasso_ct")
  M <- Reduce(function(a,b) merge(a,b,by=c("signature","seed")), list(am,t5,ls,cm,mmR))
  M$gv <- M$am_cv-M$lasso_cv; M$gt <- M$am_ct-M$lasso_ct; M$orel <- M$gv-M$gt
  M$gv5 <- M$t5_cv-M$lasso_cv; M$gt5 <- M$t5_ct-M$lasso_ct
  M$am_minus_ridge <- M$am_ct-M$ridge_ct
  M$t5_minus_am <- M$t5_ct-M$am_ct; M$t5_minus_ridge <- M$t5_ct-M$ridge_ct
  agg2 <- function(df) {  # 任务 -> 签名平均 -> 六签名等权
    s1 <- aggregate(df[MC], by=df["signature"], FUN=mean, na.rm=TRUE)
    colMeans(s1[MC])
  }
  mk2 <- function(tag, df) { a <- agg2(df)
    data.frame(stratum=tag, n_tasks=nrow(df), n_signature=length(unique(df$signature)),
               Gv=a[["gv"]], Gt=a[["gt"]], Orel=a[["orel"]], Orel_over_Gv=a[["orel"]]/a[["gv"]],
               Lasso_Ct=a[["lasso_ct"]], argmax_Ct=a[["am_ct"]], Ridge_Ct=a[["ridge_ct"]],
               top5_Ct=a[["t5_ct"]], argmax_minus_Ridge=a[["am_minus_ridge"]],
               top5_minus_argmax=a[["t5_minus_am"]], top5_minus_Ridge=a[["t5_minus_ridge"]],
               stringsAsFactors=FALSE) }
  MB <- rbind(mk2("all_tasks_METABRIC", M), mk2("nonconstant_Lasso_METABRIC", M[!M$lasso_constant, ]))
  write.table(MB, file.path(TAB,"47_nonconstant_metabric.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
} else {
  say("[注] METABRIC 预测或指标表缺失，跳过 METABRIC 描述")
}

## ============================================================
## 文本摘要
## ============================================================
fmt <- function(x) sprintf("%+.4f", x)
say("")
say("================ 全体 vs 非恒定 LASSO 子集 ================")
say("  ", sprintf("%-22s %7s %9s %9s %9s %9s %9s %9s %11s", "分层","n任务","Gv","Gt","Orel",
                "Orel/Gv","argmaxC","RidgeC","argmax−Ridge"))
for (i in seq_len(nrow(OV))) say("  ", sprintf("%-22s %7d %s %s %s %9.1f%% %9.4f %9.4f %s",
  OV$stratum[i], OV$n_tasks[i], fmt(OV$Gv[i]), fmt(OV$Gt[i]), fmt(OV$Orel[i]),
  100*OV$Orel_over_Gv[i], OV$argmax_Ct[i], OV$Ridge_Ct[i], fmt(OV$argmax_minus_Ridge[i])))
say("")
say("  top5 − argmax / top5 − Ridge（对照）：全体 ", fmt(OV$top5_minus_argmax[1])," / ",fmt(OV$top5_minus_Ridge[1]),
    "；非恒定 ", fmt(OV$top5_minus_argmax[2])," / ",fmt(OV$top5_minus_Ridge[2]))
say("")
say("  逐队列方向（argmax−Ridge）：", paste(sprintf("%s:%s", dir_tab$stratum,
    sprintf("%d/%d 为负", dir_tab$n_cohort_negative, dir_tab$n_cohort)), collapse="  "))
say("")
say("  常数任务 862 个的队列分布（仅计数，不做宏平均强推断）：")
cs <- cnt[order(-cnt$n_constant), ]
for (i in seq_len(nrow(cs))) say(sprintf("    %-8s 常数 %4d / %4d  (%.1f%%)", cs$cohort[i],
    cs$n_constant[i], cs$n_tasks[i], 100*cs$n_constant[i]/cs$n_tasks[i]))
cell_cov <- length(unique(paste(CI$cohort[CI$n_nonconstant>0], CI$seed[CI$n_nonconstant>0])))
say(sprintf("  非恒定子集覆盖 %d 个 cohort×split 单元；常数组覆盖 %d 个单元",
            cell_cov, sum(CI$n_constant>0)))
say("")
say("================ 六 OS 队列 ================")
for (i in seq_len(nrow(OV6))) say("  ", sprintf("%-24s %7d Gv=%s Gt=%s Orel=%s", OV6$stratum[i],
  OV6$n_tasks[i], fmt(OV6$Gv[i]), fmt(OV6$Gt[i]), fmt(OV6$Orel[i])))
if (!is.null(MB)) {
  say("")
  say("================ METABRIC 描述 ================")
  for (i in seq_len(nrow(MB))) say("  ", sprintf("%-28s %4d Gv=%s Gt=%s Orel=%s argmax−Ridge=%s",
    MB$stratum[i], MB$n_tasks[i], fmt(MB$Gv[i]), fmt(MB$Gt[i]), fmt(MB$Orel[i]), fmt(MB$argmax_minus_Ridge[i])))
}
say("")
say("说明：非恒定子集是事后条件分层，不是新的主分析集，也不证明因果；")
say("      LASSO 恒定解是 L1 稀疏性的合法结果，不能称为'弱正则化'。")
close(con)
cat("47 完成。输出写入 ", TAB, "\n")
