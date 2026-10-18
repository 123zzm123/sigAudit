#!/usr/bin/env Rscript
# ============================================================
# 13 — 论文主图（v3，按"收益核算 + 参照依赖性"主线；第三轮审查意见落地）
#
# 相对 v2 的修改（全部来自第三轮审查）：
#   Fig 1  色块按【真实的 50 / 25 / 25 比例】绘制（旧版视觉上训练段没占一半）；
#          删除未经证实的 "pre-registered" 用语，改为"冻结的分析划分"。
#   Fig 2  左侧 "ALL (equal weight)" 不再被裁切；COAD 统一改为 COADREAD；
#          图注说明两点之间的连线只是 Gv 与 Gt 的【差】，不是置信区间。
#   Fig 3  图注明确灰线是【跨划分范围】，不是 95% 置信区间，
#          也不要把跨划分散布与患者抽样不确定性混为一谈。
#   Fig 4  oracle 单列为灰色【描述性参照】（不可实施），与可实施规则之间加分隔线；
#          主面板突出预先指定的 P2（top5_ensemble）；新增固定 Ridge 的事后敏感性对照；
#          长规则名不再裁切（缩短标签 + 放大左边距）。
#
# 同时输出【绘图数据】（results/figures/figdata_*.tsv）与完整图注，使"新图可复现"。
#
# 颜色仍用已跑过验证脚本的两个分类色（CVD ΔE：protan 24.7、tritan 32.7；正常 33.6，全部 PASS）。
# 输入：results/tables/20_rules.tsv, 20_paired_vs_argmax.tsv,
#       29_paired_vs_fixed.tsv（若有）
# 用法：Rscript scripts/13_figures.R
# ============================================================
options(stringsAsFactors = FALSE)
PROJ <- Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit"); TAB <- file.path(PROJ, "results/tables")
FIG  <- file.path(PROJ, "results/figures"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
CAP  <- file.path(FIG, "fig_captions.md")
rd <- function(f) read.table(f, sep = "\t", header = TRUE, quote = "", comment.char = "")

C_GV   <- "#2a78d6"   # 验证增益（分类槽 1）
C_GT   <- "#eb6834"   # 独立测试增益（分类槽 2）
INK    <- "#1a1a19"; INK2 <- "#6b6b66"; GRID <- "#e6e5e0"

## COAD 实际为 COADREAD（第三轮审查要求全文统一）
ren <- function(x) ifelse(x == "COAD", "COADREAD", x)

redraw <- function(name, f, w = 7, h = 5) {
  pdf(file.path(FIG, paste0(name, ".pdf")), width = w, height = h); f(); dev.off()
  png(file.path(FIG, paste0(name, ".png")), width = w, height = h, units = "in", res = 300)
  f(); dev.off()
  cat("  ->", name, "\n")
}
blank <- function(xlim, ylim, xlab, main, ylabels, yat) {
  plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = "", main = main,
       yaxt = "n", bty = "n", cex.main = 1.05, col.axis = INK2, col.lab = INK2)
  abline(v = pretty(xlim), col = GRID, lwd = 1)
  axis(2, at = yat, labels = ylabels, las = 1, tick = FALSE, col.axis = INK, cex.axis = 0.8)
}

CAPS <- c("# 图注（与 PDF/PNG 一一对应）", "")

## ---------------- Fig 1：设计与口径示意（按真实比例）----------------
redraw("fig1_design", function() {
  par(mar = c(1, 1, 3.2, 1), family = "sans")
  plot(NA, xlim = c(0, 10), ylim = c(0, 6), axes = FALSE, xlab = "", ylab = "",
       main = "Fig 1  Design: split, selection rule, and the two gains")
  ## 【修正】色块宽度严格等于样本比例 50 / 25 / 25
  ys <- 4.55; hs <- 0.85; x0 <- 0.4; W <- 9.2
  fr <- c(0.50, 0.25, 0.25)
  ws <- W * fr
  xs <- x0 + c(0, cumsum(ws)[1:2])
  lab <- c("train 50%", "validation 25%", "held-out test 25%")
  cl  <- c("#d9d8d2", C_GV, C_GT)
  for (i in 1:3) {
    rect(xs[i], ys, xs[i] + ws[i], ys + hs, col = cl[i], border = "white", lwd = 2)
    text(xs[i] + ws[i]/2, ys + hs/2, lab[i], col = if (i == 1) INK else "#ffffff", cex = 0.85)
  }
  text(x0, ys + hs + 0.33, "one of five frozen analysis splits (seeds fixed before the revised run)",
       adj = 0, col = INK2, cex = 0.78)
  text(x0, 3.6, "Selection: argmax of the validation C-index over the candidate pool (reference allowed)",
       adj = 0, col = INK, cex = 0.82)
  text(x0, 2.9, "Gv  =  Cv(selected) - Cv(reference)            validation gain",
       adj = 0, col = C_GV, cex = 0.9)
  text(x0, 2.3, "Gt  =  Ct(selected) - Ct(reference)            held-out test gain",
       adj = 0, col = C_GT, cex = 0.9)
  text(x0, 1.5, "Orel = Gv - Gt        reference-adjusted selection optimism",
       adj = 0, col = INK, cex = 0.9)
  text(x0, 0.75, "Notes: (i) block widths are proportional to the sample fractions;",
       adj = 0, col = INK2, cex = 0.72)
  text(x0, 0.35, "(ii) Ct is a held-out test estimate, not a known truth; (iii) the reference is a fixed Lasso unless stated otherwise.",
       adj = 0, col = INK2, cex = 0.72)
})
CAPS <- c(CAPS,
  "**Fig 1.** Design of the audit. Each of eight TCGA cohorts is split into train (50%), validation (25%) and a held-out test set (25%); block widths are drawn in proportion to these fractions. The selection rule picks the candidate with the highest validation C-index from the whole pool, with the reference (fixed Lasso) allowed to participate. Gv and Gt are the selection's gains in C-index over the reference on the validation and held-out test sets respectively, and Orel = Gv - Gt is the reference-adjusted selection optimism. Ct is an estimate on a held-out set, not a known truth; the seeds were fixed before the revised run rather than publicly pre-registered.", "")

## ---------------- 读入 ----------------
fr <- file.path(TAB, "20_rules.tsv")
if (!file.exists(fr)) stop("缺少 20_rules.tsv —— 请先跑 scripts/20_analysis_corrected.R")
R <- rd(fr); R$cohort <- ren(R$cohort)
R <- R[R$rule == "argmax_val", ]
if (!nrow(R)) stop("20_rules.tsv 里没有 argmax_val 行")

## ---------------- Fig 2：收益分解 ----------------
byc <- aggregate(cbind(Gv, Gt, Orel) ~ cohort, data = R, FUN = function(x) mean(x, na.rm = TRUE))
byc <- byc[order(byc$Gv), ]
allg <- c(mean(byc$Gv), mean(byc$Gt))
FD2 <- rbind(data.frame(cohort = byc$cohort, group = "cohort", Gv = byc$Gv, Gt = byc$Gt, Orel = byc$Orel),
             data.frame(cohort = "ALL (equal weight)", group = "overall",
                        Gv = allg[1], Gt = allg[2], Orel = allg[1] - allg[2]))
write.table(FD2, file.path(FIG, "figdata_fig2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
redraw("fig2_gain_decomposition", function() {
  par(mar = c(4.2, 8.2, 3.2, 1.6), family = "sans")   # 【修正】加宽左边距，ALL 标签不再被裁切
  y <- seq_len(nrow(byc))
  yl <- c(byc$cohort, "ALL (equal weight)")
  yv_gv <- c(byc$Gv, allg[1]); yv_gt <- c(byc$Gt, allg[2])
  yy <- c(y, nrow(byc) + 1.4)
  xl <- range(c(yv_gv, yv_gt, 0)) * 1.15
  blank(xl, c(0.4, max(yy) + 0.7), "gain in C-index relative to the reference",
        "Fig 2  Reported (validation) vs held-out test gain", yl, yy)
  abline(v = 0, col = INK2, lwd = 1.5)
  segments(yv_gv, yy, yv_gt, yy, col = "#c9c8c2", lwd = 2)
  points(yv_gv, yy, pch = 16, cex = 1.5, col = C_GV)
  points(yv_gt, yy, pch = 16, cex = 1.5, col = C_GT)
  segments(0, max(yy) + 0.5, 0, 0.6, col = INK2, lwd = 1.5, lty = 3)
  abline(h = nrow(byc) + 0.7, col = GRID, lwd = 1, lty = 2)
  legend("bottomright", c("Gv (validation)", "Gt (held-out test)"), pch = 16, pt.cex = 1.4,
         col = c(C_GV, C_GT), bty = "n", text.col = INK, cex = 0.85)
  mtext("grey connector = the difference Gv - Gt (Orel); it is NOT a confidence interval",
        side = 3, line = 0.2, adj = 0.02, col = INK2, cex = 0.72)
}, w = 8.4, h = 5)
CAPS <- c(CAPS,
  "**Fig 2.** Decomposition of the reported gain for the validation-set argmax rule, per cohort and overall (equal weight over the eight cohorts; COAD is COADREAD). Blue = Gv (gain on the validation set); orange = Gt (gain on the held-out test set). Both are relative to the fixed Lasso reference. The grey connector joining the two points is the arithmetic difference Gv - Gt (i.e. Orel) and is **not** a confidence interval; no uncertainty is displayed in this panel (see Fig 3 for cross-split dispersion).", "")

## ---------------- Fig 3：跨划分稳定性（单系列，不加图例）----------------
bycs <- aggregate(Orel ~ cohort + seed, data = R, FUN = function(x) mean(x, na.rm = TRUE))
FD3 <- do.call(rbind, lapply(split(bycs, bycs$cohort), function(d)
  data.frame(cohort = d$cohort[1], n_split = nrow(d), Orel_mean = mean(d$Orel),
             Orel_min = min(d$Orel), Orel_max = max(d$Orel),
             values = paste(sprintf("%+.4f", d$Orel), collapse = " "))))
write.table(FD3, file.path(FIG, "figdata_fig3.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
redraw("fig3_cross_split", function() {
  par(mar = c(4.4, 3.4, 3.4, 1.4), family = "sans")
  coh <- sort(unique(bycs$cohort))
  xl <- range(c(bycs$Orel, 0)) * 1.2
  plot(NA, xlim = xl, ylim = c(0.4, length(coh) + 0.6), axes = FALSE, xlab = "Orel per split",
       ylab = "", main = "Fig 3  Cross-split stability of Orel (argmax on validation)")
  abline(v = pretty(xl), col = GRID, lwd = 1)
  abline(v = 0, col = INK2, lwd = 1.5)
  axis(1, col.axis = INK2, col.lab = INK2); axis(2, at = seq_along(coh), labels = coh,
       las = 1, tick = FALSE, col.axis = INK, cex.axis = 0.85)
  for (i in seq_along(coh)) {
    v <- bycs$Orel[bycs$cohort == coh[i]]
    segments(min(v), i, max(v), i, col = "#c9c8c2", lwd = 1.5)
    points(v, rep(i, length(v)), pch = 16, cex = 1.1, col = C_GV)
    points(mean(v), i, pch = 18, cex = 1.8, col = INK)
  }
  legend("bottomright", c("each split", "cohort mean"), pch = c(16, 18),
         col = c(C_GV, INK), bty = "n", text.col = INK, cex = 0.85)
  mtext("grey line = range across the five splits, NOT a 95% CI; it is not patient-sampling uncertainty",
        side = 3, line = 0.2, adj = 0.02, col = INK2, cex = 0.7)
}, w = 7, h = 4.6)
CAPS <- c(CAPS,
  "**Fig 3.** Cross-split stability of Orel for the validation-set argmax rule. Each blue point is the Orel of one of the five analysis splits (same patients, different partitions); the black diamond is the cohort mean. The grey line spans the minimum to the maximum across the five splits: it is the **range across splits**, **not** a 95% confidence interval and **not** patient-sampling uncertainty.", "")

## ---------------- Fig 4：规则配对比较 ----------------
fp <- file.path(TAB, "20_paired_vs_argmax.tsv")
if (!file.exists(fp)) stop("缺少 20_paired_vs_argmax.tsv")
P <- rd(fp)
cl <- aggregate(dCt ~ cohort + rule, data = P, FUN = function(x) mean(x, na.rm = TRUE))
mk <- function(cl, lab) {
  out <- do.call(rbind, lapply(split(cl, cl$rule), function(d) {
    x <- d$dCt[is.finite(d$dCt)]
    if (length(x) < 4) return(NULL)
    tt <- t.test(x)
    data.frame(rule = d$rule[1], n = length(x), mean = mean(x),
               lo = tt$conf.int[1], hi = tt$conf.int[2])
  }))
  if (is.null(out) || !nrow(out)) return(NULL)
  out$kind <- lab; out
}
ag_impl <- mk(cl[cl$rule != "oracle_test", ], "implementable")
ag_orac <- mk(cl[cl$rule == "oracle_test", ], "descriptive")

## 固定 Ridge 的事后敏感性对照（来自 29；符号取反：Ridge − argmax）
fr29 <- file.path(TAB, "29_paired_vs_fixed.tsv")
ag_ridge <- NULL
if (file.exists(fr29)) {
  B <- rd(fr29)
  z <- B[B$reference == "Ridge" & B$rule == "argmax_val", ]
  if (nrow(z)) ag_ridge <- data.frame(rule = "fixed Ridge", n = z$n_cohort[1],
    mean = -z$dCt[1], lo = -z$hi[1], hi = -z$lo[1], kind = "descriptive")
}

if (is.null(ag_impl) || !nrow(ag_impl)) {
  cat("  (Fig 4 跳过：可用队列 < 4)\n")
} else {
  # 主面板：可实施规则（按效果排序）；下方分隔后放描述性参照
  ag_impl <- ag_impl[order(ag_impl$mean), ]
  ref_blk <- rbind(ag_orac, ag_ridge)
  if (!is.null(ref_blk)) ref_blk <- ref_blk[order(ref_blk$mean), ]
  AG <- rbind(ag_impl, ref_blk)
  # 缩短标签，避免裁切
  short <- c(top5_ensemble = "top5 ensemble (P2)", top3_ensemble = "top3 ensemble",
             top2_ensemble = "top2 ensemble", top10_ensemble = "top10 ensemble",
             all_ensemble = "all ensemble", linear_only = "linear candidates",
             tree_only = "tree candidates", runner_up = "runner-up",
             majority_vote = "majority vote", worst_half = "worst-half",
             penalized = "penalised", none = "no selection", argmax_nonreference = "argmax (no ref)",
             oracle_test = "oracle (test-set max)", `fixed Ridge` = "fixed Ridge")
  AG$lab <- ifelse(AG$rule %in% names(short), short[AG$rule], AG$rule)
  AG$lab[is.na(AG$lab)] <- AG$rule[is.na(AG$lab)]
  write.table(AG, file.path(FIG, "figdata_fig4.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  redraw("fig4_rule_pairing", function() {
    par(mar = c(4.4, 9.4, 3.6, 1.4), family = "sans")
    nI <- nrow(ag_impl); nR <- nrow(ref_blk)
    yI <- rev(seq_len(nI)); yR <- if (nR) rev(nI + 1.2 + seq_len(nR)) else numeric(0)
    yy <- c(yI, yR)
    xl <- range(c(AG$lo, AG$hi, 0), na.rm = TRUE) * 1.15
    blank(xl, c(0.3, max(yy) + 0.7), "paired difference in held-out test C-index vs argmax_val",
          "Fig 4  Does any implemented rule beat selecting the validation-best candidate?",
          AG$lab, yy)
    abline(v = 0, col = INK2, lwd = 1.5)
    isI <- seq_len(nI)
    segments(AG$lo[isI], yI, AG$hi[isI], yI, col = "#c9c8c2", lwd = 2)
    segments(AG$lo[-isI], yR, AG$hi[-isI], yR, col = "#dedcd6", lwd = 2)
    points(AG$mean[isI], yI, pch = 16, cex = 1.4, col = C_GT)
    points(AG$mean[-isI], yR, pch = 17, cex = 1.4, col = INK2)
    # 预先指定的 P2 单独标出
    j <- which(AG$rule[isI] == "top5_ensemble")
    if (length(j)) points(AG$mean[isI][j], yI[j], pch = 1, cex = 2.3, lwd = 1.6, col = INK)
    if (nR) { abline(h = nI + 0.6, col = INK2, lwd = 1, lty = 2)
      text(xl[1] + 0.02*diff(xl), nI + 0.85, "descriptive reference (not implementable / post-hoc)",
           adj = 0, col = INK2, cex = 0.72) }
    legend("bottomright", c("implemented rule", "descriptive reference", "pre-specified P2"),
           pch = c(16, 17, 1), pt.cex = c(1.3, 1.3, 1.9),
           col = c(C_GT, INK2, INK), bty = "n", text.col = INK, cex = 0.8)
    mtext("negative = worse on the held-out test set; bars are 95% CI (unit = cohort)",
          side = 3, line = 0.2, adj = 0.02, col = INK2, cex = 0.72)
  }, w = 8.6, h = 5.6)
  CAPS <- c(CAPS,
    "**Fig 4.** Paired difference in held-out test C-index between each rule and the validation-set argmax rule (unit of inference = cohort; bars are 95% CIs from an equal-weight t-interval over the eight cohorts). To the right of 0 the rule is better than argmax on the held-out test set. Rules in the upper block are implementable; the lower block is descriptive only — `oracle (test-set max)` is the post-hoc maximum on the test set and the `fixed Ridge` row is a post-hoc sensitivity comparison added after review (Ridge minus argmax). The circled point is the pre-specified comparison P2 (`top5 ensemble`).", "")
}

writeLines(CAPS, CAP)
cat("图表输出:", FIG, "\n")
cat("图注:", CAP, "\n")
cat("绘图数据: figdata_fig2.tsv / figdata_fig3.tsv / figdata_fig4.tsv\n")
cat("提示：颜色选择已通过可视化验证脚本（见脚本头注释），请勿在下游随意改色。\n")
