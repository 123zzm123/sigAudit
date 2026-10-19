#!/usr/bin/env Rscript
# ============================================================
# 51 — 搜索规模曲线图（主图 Fig. 5 + 补充图 S2–S5）
#
# 输入：results/tables/50_search_budget_summary.tsv（八队列等权）
#       50_search_budget_by_cohort.tsv（队列级）
#       50_search_budget_metabric.tsv（METABRIC，含 ALL 与逐签名）
#       50_search_budget_metabric_task.tsv（逐任务明细）
# 输出：results/figures/fig5_search_budget.{pdf,png}
#       results/figures/figs2_search_budget_cohorts_lasso.{pdf,png}
#       results/figures/figs3_search_budget_cohorts_ridge.{pdf,png}
#       results/figures/figs4_search_budget_signatures_lasso.{pdf,png}
#       results/figures/figs5_search_budget_signatures_ridge.{pdf,png}
#       results/figures/figdata_search_budget_{fig5,figs2,figs3,figs4,figs5}.tsv
#       results/figures/fig_captions_search_budget.md
#
# 注意：文件命名与手稿图号存在历史偏移（手稿 Fig. 4 = fig5_metabric，
#   手稿 Fig. 3 = fig4_rule_pairing；fig3_cross_split 是补充图 S1）。
#   本脚本的搜索规模曲线是手稿 Fig. 5，文件名统一带 search_budget，
#   避免与脚本 33 的 fig5_metabric / figdata_fig5.tsv / fig_captions_fig5.md 冲突。
#
# 图义：横轴 = 搜索规模 K（1,2,4,8,12,full），纵轴 = 相对固定起点的增益。
#   蓝 = Gv（验证增益），橙 = Gt（留出增益），两者之差即 Orel。零线表示固定起点。
#   起点在 K 全程都在池中；K = full 时两个起点回到同一完整池。
# ============================================================
options(stringsAsFactors = FALSE)
PROJ <- Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
TAB <- file.path(PROJ, "results/tables"); FIG <- file.path(PROJ, "results/figures")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) read.table(f, sep="\t", header=TRUE, quote="", comment.char="", stringsAsFactors=FALSE)

C_GV <- "#2a78d6"; C_GT <- "#eb6834"; INK <- "#1a1a19"; INK2 <- "#6b6b66"; GRID <- "#e6e5e0"
KL <- c("1", "2", "4", "8", "12", "full")   # 固定横轴刻度

redraw <- function(name, f, w, h) {
  pdf(file.path(FIG, paste0(name, ".pdf")), width = w, height = h); f(); dev.off()
  png(file.path(FIG, paste0(name, ".png")), width = w, height = h, units = "in", res = 300)
  f(); dev.off(); cat("  ->", name, "\n")
}

## 取某个起点、按固定 K 顺序排好
sel <- function(D, st) {
  D <- D[D$start == st, ]
  D$K_label <- factor(D$K_label, levels = KL)
  D[order(D$K_label), ]
}

## 单个面板：一条起点曲线上的 Gv / Gt
curve_panel <- function(z, main, show_ylab = TRUE, letter = NULL) {
  yv <- range(c(0, z$Gv, z$Gt), na.rm = TRUE)
  yv <- yv + c(-1, 1) * 0.12 * diff(yv)
  par(mar = c(4.0, if (show_ylab) 4.2 else 2.2, 3.4, 1.2), family = "sans")
  plot(NA, xlim = c(0.8, 6.2), ylim = yv, axes = FALSE, xlab = "search size K",
       ylab = if (show_ylab) "gain in C-index vs. the fixed start" else "",
       main = main, cex.main = 0.92, col.axis = INK2, col.lab = INK2)
  abline(h = pretty(yv), col = GRID, lwd = 1)
  abline(h = 0, col = INK2, lwd = 1.4, lty = 3)
  axis(1, at = 1:6, labels = KL, col.axis = INK, col.lab = INK2, cex.axis = 0.8)
  axis(2, col.axis = INK2, cex.axis = 0.75)
  points(1:6, z$Gv, type = "b", pch = 16, cex = 1.0, lwd = 1.8, col = C_GV)
  points(1:6, z$Gt, type = "b", pch = 17, cex = 1.0, lwd = 1.8, col = C_GT)
  if (!is.null(letter))
    mtext(letter, side = 3, adj = 0, line = 1.3, cex = 1.05, font = 2, col = INK)
}

## 四面板共用图例（独立一行，避免看起来只属于某个面板）
common_legend <- function() {
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", c("Gv (validation gain)", "Gt (held-out test gain)"),
         pch = c(16, 17), col = c(C_GV, C_GT), bty = "n", horiz = TRUE,
         text.col = INK, cex = 0.9, xpd = NA)
}

## ---------------- 读入 ----------------
SUM  <- rd(file.path(TAB,"50_search_budget_summary.tsv"))
COH  <- rd(file.path(TAB,"50_search_budget_by_cohort.tsv"))
MB   <- rd(file.path(TAB,"50_search_budget_metabric.tsv"))

## ---------------- Fig 5：四个面板 A–D ----------------
TC_las <- sel(SUM, "Lasso"); TC_rid <- sel(SUM, "Ridge")
MB_las <- sel(MB[MB$signature == "ALL", ], "Lasso")
MB_rid <- sel(MB[MB$signature == "ALL", ], "Ridge")
FD5 <- rbind(data.frame(source="TCGA",     SUM[, c("start","K_label","eCv","eCt","Gv","Gt","Orel","eRandCt")]),
             data.frame(source="METABRIC", MB[MB$signature=="ALL", c("start","K_label","eCv","eCt","Gv","Gt","Orel","eRandCt")]))
write.table(FD5, file.path(FIG,"figdata_search_budget_fig5.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

redraw("fig5_search_budget", function() {
  layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE), heights = c(1, 1, 0.17))
  par(oma = c(0.4, 0.6, 1.9, 0.6))
  curve_panel(TC_las, "TCGA, fixed LASSO start", TRUE, "A")
  curve_panel(TC_rid, "TCGA, fixed ridge start", TRUE, "B")
  curve_panel(MB_las, "METABRIC, fixed LASSO start", TRUE, "C")
  curve_panel(MB_rid, "METABRIC, fixed ridge start", TRUE, "D")
  common_legend()
  mtext("Reference-anchored search budgets", side = 3, line = 0.2, outer = TRUE,
        cex = 1.05, font = 2, col = INK)
}, w = 8.6, h = 7.6)

## ---------------- 补充图 S2 / S3：八队列，两个起点 ----------------
coh_names <- sort(unique(COH$cohort))
S2 <- COH[COH$start == "Lasso", c("cohort","start","K_label","eCv","eCt","Gv","Gt","Orel","eRandCt")]
S3 <- COH[COH$start == "Ridge", c("cohort","start","K_label","eCv","eCt","Gv","Gt","Orel","eRandCt")]
write.table(S2, file.path(FIG,"figdata_search_budget_figs2.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(S3, file.path(FIG,"figdata_search_budget_figs3.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

cohort_fig <- function(name, st, letter_title) {
  redraw(name, function() {
    par(mfrow = c(3, 3), oma = c(0.6, 0.6, 2.8, 0.6))
    for (cn in coh_names) curve_panel(sel(COH[COH$cohort == cn, ], st), cn, TRUE)
    curve_panel(sel(SUM, st), "ALL (TCGA, equal weight)", TRUE)
    mtext(letter_title, side = 3, line = 0.4, outer = TRUE, cex = 0.95, col = INK)
  }, w = 10.2, h = 9.9)
}
cohort_fig("figs2_search_budget_cohorts_lasso", "Lasso",
           "Fig. S2  Search-budget curves per TCGA cohort, fixed-LASSO start")
cohort_fig("figs3_search_budget_cohorts_ridge", "Ridge",
           "Fig. S3  Search-budget curves per TCGA cohort, fixed-ridge start")

## ---------------- 补充图 S4 / S5：六签名，两个起点 ----------------
sigs <- sort(unique(MB$signature[MB$signature != "ALL"]))
## 显式写图数据，避免依赖文件名推导
S4 <- MB[MB$signature != "ALL" & MB$start == "Lasso",
         c("signature","start","K_label","eCv","eCt","Gv","Gt","Orel","eRandCt")]
S5 <- MB[MB$signature != "ALL" & MB$start == "Ridge",
         c("signature","start","K_label","eCv","eCt","Gv","Gt","Orel","eRandCt")]
write.table(S4, file.path(FIG,"figdata_search_budget_figs4.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(S5, file.path(FIG,"figdata_search_budget_figs5.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

sig_fig <- function(name, st, letter_title) {
  redraw(name, function() {
    par(mfrow = c(2, 4), oma = c(0.6, 0.6, 2.8, 0.6))
    for (sg in sigs) curve_panel(sel(MB[MB$signature == sg, ], st), sg, TRUE)
    curve_panel(sel(MB[MB$signature == "ALL", ], st), "ALL (6 signatures)", TRUE)
    plot.new()
    mtext(letter_title, side = 3, line = 0.4, outer = TRUE, cex = 0.95, col = INK)
  }, w = 11.0, h = 7.2)
}
sig_fig("figs4_search_budget_signatures_lasso", "Lasso",
        "Fig. S4  Search-budget curves per METABRIC signature, fixed-LASSO start")
sig_fig("figs5_search_budget_signatures_ridge", "Ridge",
        "Fig. S5  Search-budget curves per METABRIC signature, fixed-ridge start")

## ---------------- 图注 ----------------
CAPS <- c(
"# 图注：搜索规模曲线（Fig. 5 / Figs. S2–S5）", "",
"**Fig. 5. Reference-anchored search-budget curves.** For each of the 3,919 jointly evaluable TCGA tasks (A, B) and each of the 30 METABRIC tasks (C, D), a fixed starting model (fixed LASSO, A and C; fixed ridge, B and D) was kept in the pool and the remaining K−1 candidates were drawn uniformly without replacement; the plotted values are the exact expectation, over all such sub-pools, of the validation gain (blue, `Gv`) and the held-out test gain (orange, `Gt`) over that same fixed start, for the candidate with the highest validation C-index in the sub-pool. The dashed horizontal line marks the fixed start, so its gain is zero; the vertical distance between the two curves is `Orel`. The points are the discrete budgets K = 1, 2, 4, 8, 12 and the task's own complete analysed pool (`full`); they are equally spaced on the axis only for display, and the connecting lines are drawn to help track a series, not to imply interpolation or a linear K scale. Each panel has its own vertical scale, so magnitudes are not directly comparable across panels. The expectation is combinatorial, not sampled, and the curves carry no Monte Carlo error. This is a post-hoc analysis added in revision and is not extrapolated to large algorithm libraries.", "",
"**Fig. S2.** The same curves for each of the eight TCGA cohorts separately and their equal-weight average (bottom-right), for the fixed-LASSO start. Colours and vertical axis as in Fig. 5; each small panel has its own vertical scale.", "",
"**Fig. S3.** The same curves for each of the eight TCGA cohorts separately and their equal-weight average (bottom-right), for the fixed-ridge start. This start is the one that produces negative held-out gains in every cohort; the panel makes the cross-cohort consistency visible. Colours and vertical axis as in Fig. 5.", "",
"**Fig. S4.** The same curves for each of the six METABRIC signatures separately and their equal-weight average, for the fixed-LASSO start. Colours and vertical axis as in Fig. 5.", "",
"**Fig. S5.** The same curves for each of the six METABRIC signatures separately and their equal-weight average, for the fixed-ridge start. Unlike the corresponding TCGA start (Fig. S3), the METABRIC ridge-anchored held-out gains stay small and positive and are not uniform across signatures. Colours and vertical axis as in Fig. 5.", "")
writeLines(CAPS, file.path(FIG,"fig_captions_search_budget.md"))
cat("图与绘图数据输出:", FIG, "\n")
