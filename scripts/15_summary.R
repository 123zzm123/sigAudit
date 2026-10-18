#!/usr/bin/env Rscript
# ============================================================
# 15 — 终版数字汇总：把论文要引用的所有关键数字一次性打印出来。
# 所有数字都来自已落盘的表，不做任何重新拟合。
# 论文写作时以本脚本的输出为准（避免手抄旧数字）。
# ============================================================
options(stringsAsFactors = FALSE)
PROJ <- "D:/bioinfo/sigAudit"; TAB <- file.path(PROJ, "results/tables")
say <- function(...) cat(paste0(...), "\n", sep = "")
hdr <- function(s) { say(""); say(strrep("=", 68)); say(s); say(strrep("=", 68)) }

## ---------- 1. 队列与语料 ----------
hdr("1. 队列覆盖")
EXPECT <- c("BRCA", "COAD", "KIRC", "LGG", "LUAD", "MESO", "PAAD", "STAD")
files <- list.files(TAB, pattern = "^09_[A-Z]+_100\\.tsv$", full.names = TRUE)
got <- character(0); low <- character(0)
for (f in files) {
  coh <- sub("^09_(.+)_100\\.tsv$", "\\1", basename(f))
  D <- read.table(f, header = TRUE, sep = "\t", check.names = FALSE)
  n_gs <- length(unique(D$geneset))
  nm <- length(unique(D$model))
  got <- c(got, coh)
  if (n_gs < 50) low <- c(low, sprintf("%s(%d)", coh, n_gs))
  say(sprintf("  %-6s 基因集 %3d   可用候选 %2d   行 %d%s", coh, n_gs, nm, nrow(D),
              if (n_gs < 50) "   ★ 基因集数异常少" else ""))
}
# 静默失败防线：这类脚本对"基因名对不上"是安静跳过的，只打印"完成了 N 个"，
# 看上去像正常结束。所以必须显式检查【期望的队列是否都在】+【每个队列的基因集数】。
miss <- setdiff(EXPECT, got)
if (length(miss)) say("  ★★ 缺失队列: ", paste(miss, collapse = ", "),
                      "  —— 极可能是静默跳过（例如基因名对不上），务必排查")
if (length(low)) say("  ★★ 以下队列基因集数异常少: ", paste(low, collapse = ", "))
if (!length(miss) && !length(low)) say("  ✅ 8 个队列齐备且基因集数正常")

## ---------- 2. 决策规则核算 ----------
hdr("2. 决策规则核算（核心数字）")
S <- read.table(file.path(TAB, "12_rules.tsv"), header = TRUE, sep = "\t")
say(sprintf("  队列数 %d，规则 %d 条", length(unique(S$cohort)), length(unique(S$rule))))
a <- aggregate(cbind(reported, true_test, bias, gain_vs_ref) ~ rule, data = S, FUN = mean)
a <- a[order(-a$true_test), ]
for (i in seq_len(nrow(a))) say(sprintf(
  "  %-12s 报告值 %.4f  真实值 %.4f  偏倚 %+.4f  相对参照的真实增益 %+.4f",
  a$rule[i], a$reported[i], a$true_test[i], a$bias[i], a$gain_vs_ref[i]))
r0 <- a$reported[a$rule == "none"]; r1 <- a$reported[a$rule == "argmax_val"]
t1 <- a$gain_vs_ref[a$rule == "argmax_val"]
say(sprintf("  >>> 报告增益 %+.4f  vs  真实增益 %+.4f   ->  虚高占报告增益的 %.0f%%",
            r1 - r0, t1, 100 * ((r1 - r0) - t1) / (r1 - r0)))
say("")
say("  argmax_val 逐队列的真实增益：")
s2 <- S[S$rule == "argmax_val", ]; s2 <- s2[order(-s2$gain_vs_ref), ]
for (i in seq_len(nrow(s2))) say(sprintf(
  "    %-6s %+.4f [%+.4f,%+.4f]  p=%.2e   %s",
  s2$cohort[i], s2$gain_vs_ref[i], s2$ci_lo[i], s2$ci_hi[i], s2$t_p[i],
  if (is.na(s2$t_p[i])) "n/a" else if (s2$ci_lo[i] > 0) "显著为正" else "与不做选择无差异"))
say("")
say("  三个民间对策 vs 标准做法（真实增益差）：")
for (r in c("runner_up", "linear_only", "tree_only")) {
  x <- S[S$rule == r, ]; y <- S[S$rule == "argmax_val", ]
  x <- x[match(y$cohort, x$cohort), ]
  d <- x$true_test - y$true_test
  say(sprintf("    %-12s 真实值差 %+.4f （负 = 比标准做法更差）", r, mean(d, na.rm = TRUE)))
}

## ---------- 3. 进入效应 ----------
hdr("3. 进入效应（主估计量）")
E <- read.table(file.path(TAB, "11_candidate_features.tsv"), header = TRUE, sep = "\t")
say(sprintf("  队列数 %d，候选层面样本 %d", length(unique(E$cohort)), nrow(E)))
say("  进入效应排名第一的候选（各队列）：")
for (cn in unique(E$cohort)) {
  s <- E[E$cohort == cn, ]
  w <- s$target[which.max(s$entry_mean)]; o <- s$target[which.max(s$optimism)]
  say(sprintf("    %-6s 危害最大 %-13s (%+.4f) | 最乐观 %-13s (%+.4f) | %s",
              cn, w, max(s$entry_mean), o, max(s$optimism),
              if (w == o) "一致" else "★不一致★"))
}
say("")
say("  候选层面指标 vs 进入效应的相关（合并全部队列）：")
lbl <- c(mean_dval = "mean dval   (验证集可观测)", sd_dval = "sd of dval  (验证集可观测)",
         winrate = "win rate    (验证集可观测)", mean_dtest = "mean dtest  (需测试集)",
         cor_vt = "cor(dval,dtest) (需测试集)", optimism = "mean optimism  (需测试集)",
         tail_opt = "upper-tail optimism (需测试集)")
for (v in names(lbl)) { ct <- cor.test(E[[v]], E$entry_mean)
  say(sprintf("    %-28s r = %+.3f  p = %.2e", lbl[v], ct$estimate, ct$p.value)) }
cen <- function(x) x - ave(x, E$cohort, FUN = function(z) mean(z, na.rm = TRUE))
ec <- cen(E$entry_mean)
say("  队列内中心化后（剔除队列间差异）：")
for (v in names(lbl)) { ct <- suppressWarnings(cor.test(cen(E[[v]]), ec))
  say(sprintf("    %-28s r = %+.3f  p = %.2e", lbl[v], ct$estimate, ct$p.value)) }
say("")
say("  OLS 用【仅验证集可观测】指标预测进入效应：")
fit <- lm(entry_mean ~ mean_dval + sd_dval + winrate, data = E)
say(sprintf("    R^2 = %.3f   校正 R^2 = %.3f   n = %d",
            summary(fit)$r.squared, summary(fit)$adj.r.squared, nrow(E)))
say(sprintf("    加上需测试集的指标后 R^2 = %.3f",
            summary(lm(entry_mean ~ mean_dval + sd_dval + winrate +
                         mean_dtest + cor_vt + optimism + tail_opt, data = E))$r.squared))

## ---------- 4. 家族对比 ----------
hdr("4. 家族对比（同样候选数，均不含参照）")
T0 <- read.table(file.path(TAB, "10_tests.tsv"), header = TRUE, sep = "\t")
T0 <- T0[grepl("^同样 k=", T0$contrast), ]
for (i in seq_len(nrow(T0))) say(sprintf(
  "  %-6s 树集子池 − 线性子池 = %+.4f [%+.4f,%+.4f]  p=%.2e   n=%d   %s",
  T0$cohort[i], T0$mean_diff[i], T0$ci_lo[i], T0$ci_hi[i], T0$t_p[i], T0$n[i],
  if (is.na(T0$mean_diff[i])) "n/a" else if (T0$mean_diff[i] > 0) "树更糟" else "树更好"))
say("")
say("  同上，改用全池 − 纯线性（候选数不同）：")
T1 <- read.table(file.path(TAB, "10_tests.tsv"), header = TRUE, sep = "\t")
T1 <- T1[grepl("^全池", T1$contrast), ]
for (i in seq_len(nrow(T1))) say(sprintf(
  "  %-6s %+.4f [%+.4f,%+.4f]  p=%.2e", T1$cohort[i], T1$mean_diff[i], T1$ci_lo[i], T1$ci_hi[i], T1$t_p[i]))

## ---------- 5. 剂量曲线 ----------
hdr("5. 剂量曲线（典型排列，各队列）")
D <- read.table(file.path(TAB, "10_dose_curves.tsv"), header = TRUE, sep = "\t")
for (cn in unique(D$cohort)) {
  d <- D[D$cohort == cn & D$order == "A_canonical", ]; d <- d[order(d$k), ]
  dd <- diff(d$bias); j <- which.max(dd)
  say(sprintf("  %-6s k=2 %.4f -> k=%d %.4f ; 最大单步增量在 k=%d->%d 「%s」 %+.4f",
              cn, d$bias[1], max(d$k), d$bias[nrow(d)],
              d$k[j], d$k[j + 1], d$model_at_k[j + 1], dd[j]))
}
say("")
say("完成后请核对：本脚本所有数字应与 logs/10_analysis.log、11_predictors.log、12_rules.log 一致。")
