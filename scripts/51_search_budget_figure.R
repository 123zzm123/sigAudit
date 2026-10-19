#!/usr/bin/env Rscript
# ============================================================
# 51 — 搜索规模曲线图（主图 Fig. 5 + 补充图 S2/S3）
#
# 输入：results/tables/50_search_budget_summary.tsv（八队列等权）
#       50_search_budget_by_cohort.tsv（队列级）
#       50_search_budget_metabric.tsv（METABRIC，含 ALL 与逐签名）
# 输出：results/figures/fig5_search_budget.{pdf,png}
#       results/figures/figs2_search_budget_cohorts.{pdf,png}
#       results/figures/figs3_search_budget_signatures.{pdf,png}
#       results/figures/figdata_search_budget_fig5.tsv / _figs2.tsv / _figs3.tsv
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
curve_panel <- function(z, main, show_ylab = TRUE, show_legend = FALSE) {
  yv <- range(c(0, z$Gv, z$Gt), na.rm = TRUE)
  yv <- yv + c(-1, 1) * 0.12 * diff(yv)
  par(mar = c(4.0, if (show_ylab) 4.2 else 2.2, 3.2, 1.2), family = "sans")
  plot(NA, xlim = c(0.8, 6.2), ylim = yv, axes = FALSE, xlab = "search size K",
       ylab = if (show_ylab) "gain in C-index vs. the fixed start" else "",
       main = main, cex.main = 0.95, col.axis = INK2, col.lab = INK2)
  abline(h = pretty(yv), col = GRID, lwd = 1)
  abline(h = 0, col = INK2, lwd = 1.4, lty = 3)
  axis(1, at = 1:6, labels = KL, col.axis = INK, col.lab = INK2, cex.axis = 0.8)
  axis(2, col.axis = INK2, cex.axis = 0.75)
  points(1:6, z$Gv, type = "b", pch = 16, cex = 1.0, lwd = 1.8, col = C_GV)
  points(1:6, z$Gt, type = "b", pch = 17, cex = 1.0, lwd = 1.8, col = C_GT)
  if (show_legend)
    legend("topleft", c("Gv (validation)", "Gt (held-out test)"), pch = c(16, 17),
           col = c(C_GV, C_GT), bty = "n", text.col = INK, cex = 0.78)
}

## ---------------- 读入 ----------------
SUM  <- rd(file.path(TAB,"50_search_budget_summary.tsv"))
COH  <- rd(file.path(TAB,"50_search_budget_by_cohort.tsv"))
MB   <- rd(file.path(TAB,"50_search_budget_metabric.tsv"))

## ---------------- Fig 5：四个面板 ----------------
TC_las <- sel(SUM, "Lasso"); TC_rid <- sel(SUM, "Ridge")
MB_las <- sel(MB[MB$signature == "ALL", ], "Lasso")
MB_rid <- sel(MB[MB$signature == "ALL", ], "Ridge")
FD5 <- rbind(data.frame(source="TCGA",     SUM[, c("start","K_label","eCv","eCt","Gv","Gt","Orel","eRandCt")]),
             data.frame(source="METABRIC", MB[MB$signature=="ALL", c("start","K_label","eCv","eCt","Gv","Gt","Orel","eRandCt")]))
write.table(FD5, file.path(FIG,"figdata_search_budget_fig5.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

redraw("fig5_search_budget", function() {
  par(mfrow = c(2, 2), oma = c(0.6, 0.6, 2.6, 0.6))
  curve_panel(TC_las, "TCGA — start = fixed LASSO", TRUE)
  curve_panel(TC_rid, "TCGA — start = fixed ridge", TRUE)
  curve_panel(MB_las, "METABRIC — start = fixed LASSO", TRUE, show_legend = TRUE)
  curve_panel(MB_rid, "METABRIC — start = fixed ridge", TRUE)
  mtext("Fig 5  Reference-anchored search-budget curves: validation gain vs held-out gain as the pool grows",
        side = 3, line = 0.4, outer = TRUE, cex = 1.0, col = INK)
}, w = 8.4, h = 7.2)

## ---------------- Supplementary Fig. S2：八条队列曲线（LASSO 起点）----------------
coh_names <- sort(unique(COH$cohort))
FDS2 <- COH[COH$start == "Lasso", c("cohort","start","K_label","eCv","eCt","Gv","Gt","Orel","eRandCt")]
write.table(FDS2, file.path(FIG,"figdata_search_budget_figs2.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
redraw("figs2_search_budget_cohorts", function() {
  par(mfrow = c(3, 3), oma = c(0.6, 0.6, 2.8, 0.6))
  for (cn in coh_names) curve_panel(sel(COH[COH$cohort == cn, ], "Lasso"), cn, TRUE)
  curve_panel(TC_las, "ALL (TCGA, equal weight)", TRUE, show_legend = TRUE)
  mtext("Fig. S2  Search-budget curves per TCGA cohort, fixed-LASSO start",
        side = 3, line = 0.4, outer = TRUE, cex = 0.95, col = INK)
}, w = 10.2, h = 9.9)

## ---------------- Supplementary Fig. S3：六条签名曲线（LASSO 起点）----------------
sigs <- sort(unique(MB$signature[MB$signature != "ALL"]))
FDS3 <- MB[MB$signature != "ALL" & MB$start == "Lasso", c("signature","start","K_label","eCv","eCt","Gv","Gt","Orel","eRandCt")]
write.table(FDS3, file.path(FIG,"figdata_search_budget_figs3.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
redraw("figs3_search_budget_signatures", function() {
  par(mfrow = c(2, 4), oma = c(0.6, 0.6, 2.8, 0.6))
  for (sg in sigs) curve_panel(sel(MB[MB$signature == sg, ], "Lasso"), sg, TRUE)
  curve_panel(MB_las, "ALL (6 signatures)", TRUE)
  plot.new()
  mtext("Fig. S3  Search-budget curves per METABRIC signature, fixed-LASSO start",
        side = 3, line = 0.4, outer = TRUE, cex = 0.95, col = INK)
}, w = 11.0, h = 7.2)

## ---------------- 图注 ----------------
CAPS <- c(
"# 图注：搜索规模曲线（Fig. 5 / Fig. S2 / Fig. S3）", "",
"**Fig. 5. Reference-anchored search-budget curves.** For each of the 3,919 jointly evaluable TCGA tasks (top row) and each of the 30 METABRIC tasks (bottom row), a fixed starting model (fixed LASSO, left; fixed ridge, right) was kept in the pool and the remaining K-1 candidates were drawn uniformly without replacement; the plotted values are the exact expectation, over all such sub-pools, of the validation gain (blue, Gv) and the held-out test gain (orange, Gt) over that same fixed start, for the candidate with the highest validation C-index in the sub-pool. The dashed horizontal line marks the fixed start, so its gain is zero; the vertical distance between the two curves is Orel. `full` is each task's own complete analysed pool. The expectation is combinatorial, not sampled. This is a post-hoc analysis added in revision and is not extrapolated to large algorithm libraries.", "",
"**Fig. S2.** The same curves for each of the eight TCGA cohorts separately and their equal-weight average (bottom-right), for the fixed-LASSO start. Colours and vertical axis as in Fig. 5.", "",
"**Fig. S3.** The same curves for each of the six METABRIC signatures separately and their equal-weight average, for the fixed-LASSO start. Colours and vertical axis as in Fig. 5.", "")
writeLines(CAPS, file.path(FIG,"fig_captions_search_budget.md"))
cat("图与绘图数据输出:", FIG, "\n")
