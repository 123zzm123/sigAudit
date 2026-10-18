#!/usr/bin/env Rscript
# ============================================================
# 33 — Fig. 5：独立队列（METABRIC）的参照依赖性与逐签名乐观度
#
# 两个面板：
#   A  三种固定参照下的 Orel，TCGA vs METABRIC（分组柱）
#      —— 绝对乐观度在两种数据来源下都稳定，且 METABRIC 更小
#   B  METABRIC 逐签名 Orel（argmax vs 固定 LASSO）
#      —— 六个签名全部为正，范围 +0.0028 ~ +0.0170
#
# 规范同 13_figures.R：单一坐标轴、细笔画、文字用墨色、>=2 系列必有图例、
# PDF 用 ASCII 标签。颜色通过 validate_palette.js 验证过。
# 输入：32_metabric_baseline_sensitivity.tsv / 31_metabric_model_metrics.tsv
#       31_metabric_rule_performance.tsv / 29_baseline_sensitivity.tsv（TCGA 对照）
# 输出：results/figures/fig5_metabric.{pdf,png} + figdata_fig5.tsv + fig_captions_fig5.md
# ============================================================
options(stringsAsFactors = FALSE)
PROJ <- Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
TAB  <- file.path(PROJ, "results/tables"); FIG <- file.path(PROJ, "results/figures")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

C_A <- "#2a78d6"; C_B <- "#eb6834"
INK <- "#1a1a19"; INK2 <- "#6b6b66"; GRID <- "#e6e5e0"
rd <- function(f) read.table(f, sep="\t", header=TRUE, quote="", comment.char="",
                             stringsAsFactors=FALSE)

MB  <- rd(file.path(TAB, "32_metabric_baseline_sensitivity.tsv"))
TC  <- rd(file.path(TAB, "29_baseline_sensitivity.tsv"))
MM  <- rd(file.path(TAB, "31_metabric_model_metrics.tsv"))
RP  <- rd(file.path(TAB, "31_metabric_rule_performance.tsv"))

## ---- 面板 A 数据 ----
refs <- c("Lasso", "Ridge", "Enet0.5")
tc <- sapply(refs, function(r) TC$Orel[TC$reference == r & TC$rule == "argmax_val"])
mb <- sapply(refs, function(r) MB$Orel[MB$reference == r & MB$rule == "argmax_val"])

## ---- 面板 B 数据：逐签名 Orel（argmax vs 固定 LASSO，同一任务配对）----
MM$task <- paste(MM$signature, MM$seed, sep="|")
RP$task <- paste(RP$signature, RP$seed, sep="|")
L0 <- MM[MM$model == "Lasso", c("task", "cv", "ct")]; names(L0)[2:3] <- c("l_cv", "l_ct")
A  <- merge(RP[RP$rule == "argmax_val", c("task", "signature", "cv", "ct")], L0, by = "task")
A$Orel <- (A$cv - A$l_cv) - (A$ct - A$l_ct)
per <- aggregate(Orel ~ signature, data = A, FUN = mean)
per <- per[order(per$Orel), ]

draw <- function() {
  par(mfrow = c(1, 2), mar = c(4.6, 5.2, 3.8, 1.2), family = "sans")

  ## ---------- A ----------
  yl <- range(c(tc, mb, 0)) * 1.28
  bp <- barplot(rbind(tc, mb), beside = TRUE, col = c(C_A, C_B), border = NA,
                names.arg = refs, ylim = yl, las = 1,
                ylab = "Orel (reference-adjusted selection optimism)",
                main = "A  Optimism under three fixed references",
                cex.names = 0.95, col.axis = INK2, col.lab = INK2, cex.main = 1.0)
  grid(nx = NA, ny = NULL, col = GRID, lwd = 1); abline(h = 0, col = INK2, lwd = 1.2)
  text(bp[1, ], tc, sprintf("%.3f", tc), pos = 3, cex = 0.72, col = INK)
  text(bp[2, ], mb, sprintf("%.3f", mb), pos = 3, cex = 0.72, col = INK)
  legend("topright", c("TCGA (8 cohorts)", "METABRIC (6 signatures)"),
         fill = c(C_A, C_B), border = NA, bty = "n", text.col = INK, cex = 0.82)
  mtext("absolute optimism is stable across references in both sources",
        side = 3, line = 0.15, adj = 0.02, col = INK2, cex = 0.66)

  ## ---------- B ----------
  par(mar = c(4.6, 9.2, 3.8, 1.2))
  y <- seq_len(nrow(per))
  xl <- c(0, max(per$Orel) * 1.25)
  plot(NA, xlim = xl, ylim = c(0.4, nrow(per) + 0.6), yaxt = "n", bty = "n",
       xlab = "Orel (argmax vs fixed LASSO)", ylab = "",
       main = "B  METABRIC: optimism by signature",
       col.axis = INK2, col.lab = INK2, cex.main = 1.0)
  abline(v = pretty(xl), col = GRID, lwd = 1)
  axis(2, at = y, labels = per$signature, las = 1, tick = FALSE, col.axis = INK, cex.axis = 0.8)
  rect(0, y - 0.32, per$Orel, y + 0.32, col = C_B, border = NA)
  text(per$Orel, y, sprintf("%.4f", per$Orel), pos = 4, cex = 0.7, col = INK, xpd = NA)
  mtext(sprintf("positive for %d/%d signatures (range %.4f-%.4f)",
                sum(per$Orel > 0), nrow(per), min(per$Orel), max(per$Orel)),
        side = 3, line = 0.15, adj = 0.02, col = INK2, cex = 0.66)
}
pdf(file.path(FIG, "fig5_metabric.pdf"), width = 12.5, height = 5.0); draw(); dev.off()
png(file.path(FIG, "fig5_metabric.png"), width = 12.5, height = 5.0, units = "in", res = 300)
draw(); dev.off()

out <- rbind(
  data.frame(panel = "A", x = refs, TCGA = tc, METABRIC = mb, stringsAsFactors = FALSE),
  data.frame(panel = "B", x = per$signature, TCGA = NA_real_, METABRIC = per$Orel,
             stringsAsFactors = FALSE))
out$note <- c(rep("Orel under each reference (argmax)", 3), rep("per-signature Orel", nrow(per)))
write.table(out, file.path(FIG, "figdata_fig5.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## 图注以正文为唯一版本（第九轮）：此处逐字输出正文 §3.10 之后的 Fig. 5 图注
writeLines(c(
  "**Fig. 5. Independent-cohort check of the reference dependence.** (A) `Orel` for validation-set",
  "selection under three fixed references, in the TCGA benchmark (blue) and in METABRIC (orange).",
  "In both sources the three values are close to one another; the magnitude is smaller in",
  "METABRIC (about 0.008 against 0.037–0.038 in TCGA). (B) `Orel` by signature in METABRIC",
  "(validation-set selection against the fixed LASSO), after legacy symbols were mapped to current",
  "HGNC symbols and every signature was used at all unique gene symbols extracted from the",
  "`genefu` source objects. Optimism is positive for all six signatures. All values are",
  "descriptive: the six signatures and the five splits share the same patients, so the 30 tasks",
  "are not independent replicates and no p-values are reported."
), file.path(FIG, "fig_captions_fig5.md"))

cat("-> fig5_metabric (pdf+png) + figdata_fig5.tsv + fig_captions_fig5.md\n")
cat(sprintf("   A: TCGA %s | METABRIC %s\n",
            paste(sprintf("%.3f", tc), collapse="/"), paste(sprintf("%.3f", mb), collapse="/")))
cat(sprintf("   B: Orel 全为正 %d/%d，范围 %.4f–%.4f\n",
            sum(per$Orel > 0), nrow(per), min(per$Orel), max(per$Orel)))
