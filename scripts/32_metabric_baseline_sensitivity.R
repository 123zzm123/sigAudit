#!/usr/bin/env Rscript
# ============================================================
# 32 — METABRIC：参照敏感性分析（与 29 同构，但按签名汇总）
#
# 回答评审关心的那个问题：**"收益依赖基线"的判断能否在独立来源队列上再现？**
#
# 汇总层级：任务（签名×划分） -> 签名内跨 5 个划分平均 -> 6 个签名等权
#           （5 次划分共享患者、6 个签名也共享患者，因此不做显著性检验）
#
# 输入：31_metabric_model_metrics.tsv（逐候选 Cv/Ct）
#       31_metabric_rule_performance.tsv（逐规则 Cv/Ct）
# 输出：32_metabric_baseline_sensitivity.tsv
#       32_metabric_fixed_candidates.tsv
#       32_metabric_per_signature_vs_ridge.tsv
#       32_metabric_vs_tcga.tsv
# ============================================================
suppressMessages(library(survival))
options(stringsAsFactors = FALSE)
PROJ <- Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
TAB  <- file.path(PROJ, "results/tables")
con <- file(file.path(PROJ, "logs", "32_metabric_baseline_sensitivity.log"), "wt")
say <- function(...) { m <- paste0(...); cat(m,"\n",sep=""); cat(m,"\n",sep="",file=con); flush(con) }
rd <- function(f) read.table(f, sep="\t", header=TRUE, quote="", comment.char="", stringsAsFactors=FALSE)

MM <- rd(file.path(TAB, "31_metabric_model_metrics.tsv"))
RP <- rd(file.path(TAB, "31_metabric_rule_performance.tsv"))
say("模型指标 ", nrow(MM), " 条；规则性能 ", nrow(RP), " 条")
say("签名 ", length(unique(MM$signature)), "；划分 ", length(unique(MM$seed)),
    "；候选 ", length(unique(MM$model)))

## 任务键
MM$task <- paste(MM$signature, MM$seed, sep="|")
RP$task <- paste(RP$signature, RP$seed, sep="|")

## 层级：签名内跨划分 -> 签名等权
strat <- function(df, cols, extra="rule") {
  s1 <- aggregate(df[cols], by = df[c("signature","seed",extra)], FUN = mean, na.rm = TRUE)
  s2 <- aggregate(s1[cols], by = s1[c("signature",extra)], FUN = mean, na.rm = TRUE)
  aggregate(s2[cols], by = s2[extra], FUN = mean, na.rm = TRUE)
}

## ---------- A. 所有固定候选的汇总 ----------
fixed <- do.call(rbind, lapply(split(MM, MM$model), function(d) {
  s1 <- aggregate(cbind(cv, ct) ~ signature + seed, data = d, FUN = mean, na.rm = TRUE)
  s2 <- aggregate(cbind(cv, ct) ~ signature, data = s1, FUN = mean, na.rm = TRUE)
  data.frame(model = d$model[1], n_tasks = nrow(d), n_sig = nrow(s2),
             Cv = mean(s2$cv), Ct = mean(s2$ct))
}))
fixed <- fixed[order(-fixed$Ct), ]
write.table(fixed, file.path(TAB, "32_metabric_fixed_candidates.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
say("")
say("===== A. 固定候选的汇总性能（签名等权）=====")
say("  ", sprintf("%-14s %7s %8s %8s", "候选","n任务","验证C","测试C"))
for (i in seq_len(nrow(fixed))) say("  ", sprintf("%-14s %7d %8.4f %8.4f",
  fixed$model[i], fixed$n_tasks[i], fixed$Cv[i], fixed$Ct[i]))

## ---------- B. 参照敏感性 ----------
refs <- c("Lasso", "Ridge", "Enet0.5")
BS <- list(); PA <- list(); Jall <- list()
for (rf in refs) {
  R <- MM[MM$model == rf, c("task","cv","ct")]
  names(R)[2:3] <- c("ref_cv","ref_ct")
  J <- merge(RP, R, by = "task")   # RP 自带 signature/seed，避免列名冲突
  J <- J[is.finite(J$cv) & is.finite(J$ct), ]
  J$Gv <- J$cv - J$ref_cv; J$Gt <- J$ct - J$ref_ct; J$Orel <- J$Gv - J$Gt
  J$dCt <- J$ct - J$ref_ct
  Jall[[rf]] <- J
  ag <- strat(J, c("cv","ct","Gv","Gt","Orel"))
  ag$reference <- rf
  ntask <- aggregate(list(n = J$cv), by = list(rule = J$rule), FUN = function(x) sum(is.finite(x)))
  ag$n_tasks <- ntask$n[match(ag$rule, ntask$rule)]
  BS[[rf]] <- ag
  ## 逐签名的配对差（不跨签名合并做检验）
  pl <- aggregate(dCt ~ signature + rule, data = J, FUN = mean, na.rm = TRUE)
  PA[[rf]] <- pl
}
BS <- do.call(rbind, BS)
BS <- BS[order(BS$reference, -BS$ct), ]
write.table(BS, file.path(TAB, "32_metabric_baseline_sensitivity.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

say("")
say("===== B. 每个（参照 × 规则）=====")
for (rf in refs) {
  z <- BS[BS$reference == rf, ]
  say(sprintf("  --- 参照 = 固定 %s ---", rf))
  say("  ", sprintf("%-22s %8s %8s %9s %9s %9s", "规则","n任务","Ct","Gv","Gt","Orel"))
  for (i in seq_len(nrow(z))) say("  ", sprintf("%-22s %8d %8.4f %+9.4f %+9.4f %+9.4f",
    z$rule[i], z$n_tasks[i], z$ct[i], z$Gv[i], z$Gt[i], z$Orel[i]))
}

## ---------- C. 逐签名：argmax/top5 相对固定 Ridge ----------
Jr <- Jall[["Ridge"]]
Jr <- Jr[Jr$rule %in% c("argmax_val","top5_ensemble"), c("signature","rule","dCt")]
PS <- aggregate(dCt ~ signature + rule, data = Jr, FUN = mean, na.rm = TRUE)
wide <- reshape(PS, idvar = "signature", timevar = "rule", direction = "wide")
## ⚠️ reshape 按 timevar 的【字母序】产生列：argmax_val 排在 top5_ensemble 之前。
##    旧实现按【位置】赋名字，结果两列被写反（见第四轮补充意见）。
##    现在按【名字匹配】取列，不再依赖位置。
nmw <- names(wide)
i_arg <- grep("argmax_val",   nmw); i_top <- grep("top5_ensemble", nmw)
if (length(i_arg) != 1 || length(i_top) != 1)
  stop("per-signature reshape 的列名无法唯一匹配：", paste(nmw, collapse=", "))
wide <- data.frame(signature = wide$signature,
                   argmax_minus_Ridge = wide[[i_arg]],
                   top5_minus_Ridge   = wide[[i_top]],
                   stringsAsFactors = FALSE)
wide <- wide[order(wide$argmax_minus_Ridge), ]
write.table(wide, file.path(TAB, "32_metabric_per_signature_vs_ridge.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
say("")
say("===== C. 逐签名：相对固定 Ridge 的测试差 =====")
say("  ", sprintf("%-20s %16s %16s", "签名", "argmax−Ridge", "top5−Ridge"))
for (i in seq_len(nrow(wide))) say("  ", sprintf("%-20s %+16.4f %+16.4f",
  wide$signature[i], wide$argmax_minus_Ridge[i], wide$top5_minus_Ridge[i]))
say(sprintf("  argmax 低于 Ridge 的签名数: %d / %d",
            sum(wide$argmax_minus_Ridge < 0), nrow(wide)))
say(sprintf("  top5   低于 Ridge 的签名数: %d / %d",
            sum(wide$top5_minus_Ridge < 0), nrow(wide)))

## ---------- D. 与 TCGA 主结果的对照 ----------
get <- function(rf, rl, col) BS[[col]][BS$reference == rf & BS$rule == rl]
tcga <- data.frame(
  quantity = c("argmax|Lasso: Gv","argmax|Lasso: Gt","argmax|Lasso: Orel",
               "argmax|Ridge: Gv","argmax|Ridge: Gt","argmax|Ridge: Orel",
               "argmax|Enet0.5: Orel",
               "固定Ridge测试C","固定Lasso测试C",
               "top5−argmax (Ct)"),
  TCGA   = c(0.061012, 0.023657, 0.037355,
             0.031366, -0.006382, 0.037748, 0.037836,
             0.614274, 0.584235, 0.004808),
  METABRIC = c(get("Lasso","argmax_val","Gv"), get("Lasso","argmax_val","Gt"),
               get("Lasso","argmax_val","Orel"),
               get("Ridge","argmax_val","Gv"), get("Ridge","argmax_val","Gt"),
               get("Ridge","argmax_val","Orel"),
               get("Enet0.5","argmax_val","Orel"),
               fixed$Ct[fixed$model=="Ridge"], fixed$Ct[fixed$model=="Lasso"],
               get("Lasso","top5_ensemble","ct") - get("Lasso","argmax_val","ct")),
  stringsAsFactors = FALSE)
## 补上 top5 − 固定Ridge（TCGA 为 −0.0016，METABRIC 需据实填入）
top5_r <- mean(wide$top5_minus_Ridge)
tcga <- rbind(tcga, data.frame(quantity="top5 − Ridge (Ct)", TCGA=-0.0016,
                               METABRIC=top5_r, stringsAsFactors=FALSE))
tcga$direction_agrees <- sign(tcga$TCGA) == sign(tcga$METABRIC)
write.table(tcga, file.path(TAB, "32_metabric_vs_tcga.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
say("")
say("===== D. 与 TCGA 主结果的方向对照 =====")
say("  ", sprintf("%-24s %12s %12s %s", "量","TCGA","METABRIC","方向一致"))
for (i in seq_len(nrow(tcga))) say("  ", sprintf("%-24s %+12.4f %+12.4f %s",
  tcga$quantity[i], tcga$TCGA[i], tcga$METABRIC[i],
  if (tcga$direction_agrees[i]) "是" else "**否**"))
say("")
say("  方向一致的项数: ", sum(tcga$direction_agrees), " / ", nrow(tcga))
say("")
say("⚠️ 6 个签名共享患者、5 次划分也共享患者 -> 以上均为描述性汇总，不做显著性检验。")
close(con)
