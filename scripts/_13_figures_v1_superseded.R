#!/usr/bin/env Rscript
# ============================================================
# 13 - 论文图表（base graphics，英文标签以避免 Windows 中文乱码）
# 输出 results/figures/*.pdf 与 *.png
# 数据来源全部是已落盘的表，可随时重跑：
#   12_rules.tsv           -> Fig 2 代价核算（主图）
#   08_<C>_100.tsv         -> Fig 3 "打印点错觉"（v5 完整曲线）
#   11_candidate_features  -> Fig 4 进入效应排名 / Fig 5 机制散点
#   10_tests.tsv           -> Fig 6 家族对比
# ============================================================
options(stringsAsFactors = FALSE)
PROJ <- "D:/bioinfo/sigAudit"
TAB  <- file.path(PROJ, "results/tables")
FIG  <- file.path(PROJ, "results/figures"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

OPEN <- function(name, w = 7, h = 5) {
  pdf(file.path(FIG, paste0(name, ".pdf")), width = w, height = h) }
redraw <- function(name, f, w = 7, h = 5) {   # PDF 与 PNG 用同一段绘图代码
  OPEN(name, w, h); f(); dev.off()
  png(file.path(FIG, paste0(name, ".png")), width = w, height = h, units = "in", res = 300)
  f(); dev.off()
  cat("  ->", name, "\n")
}
COL <- c(reported = "#c0392b", true = "#2471a3")
LINCOL <- c(linear = "#2e86c1", tree = "#c0392b")

## ---------- Fig 2：代价核算（主图）----------
f <- file.path(TAB, "12_rules.tsv")
if (file.exists(f)) {
  S <- read.table(f, header = TRUE, sep = "\t")
  ord <- c("none","argmax_val","runner_up","linear_only","tree_only","oracle_test")
  lab <- c(none="no selection (ref)", argmax_val="argmax (val)", runner_up="runner-up",
           linear_only="linear only", tree_only="tree only", oracle_test="oracle (test)")
  cohs <- unique(S$cohort); cohs <- c(cohs, "ALL (mean)")
  redraw("fig2_rule_accounting", function() {
    par(mfrow = c(2, 3), mar = c(4.2, 4.2, 3, 1), mgp = c(2.4, .8, 0))
    for (cn in cohs) {
      d <- if (cn == "ALL (mean)")
        aggregate(cbind(reported, true_test, bias, gain_vs_ref) ~ rule, data = S, FUN = mean)
        else S[S$cohort == cn, ]
      d <- d[match(ord, d$rule), ]; d <- d[!is.na(d$rule), ]
      y <- seq_len(nrow(d))
      plot(NA, xlim = range(c(d$reported, d$true_test)) + c(-.01, .01),
           ylim = c(.4, nrow(d) + .6), yaxt = "n", xlab = "C-index", ylab = "",
           main = cn, cex.main = 1.1)
      axis(2, at = y, labels = lab[d$rule], las = 1, cex.axis = .75)
      segments(d$reported, y, d$true_test, y, col = "grey55", lwd = 2)
      points(d$reported, y, pch = 16, col = COL["reported"], cex = 1.3)
      points(d$true_test, y, pch = 17, col = COL["true"], cex = 1.3)
    }
    mtext("Fig 2  Reported vs true C-index by model-selection rule", outer = TRUE, line = -1.2, font = 2)
    mtext("red circle = reported (validation set);  blue triangle = true (held-out test set)",
          outer = TRUE, line = -2.4, cex = .85, col = "grey25")
  }, w = 12, h = 6.5)
}

## ---------- Fig 3：v5 "打印点错觉" ----------
V5 <- c("Lasso","Ridge","Lasso_1se","Enet0.5","Enet0.5_1se","Enet0.9","Enet0.7",
        "Enet0.3","Enet0.1","Enet0.9_1se","Enet0.7_1se","Enet0.3_1se","Enet0.1_1se",
        "Ridge_1se","UniCox+Lasso","RSF100","GBM100","GBM250","RSF250")
fs <- list.files(TAB, pattern = "^08_[A-Z]+_100\\.tsv$", full.names = TRUE)
if (length(fs)) {
  redraw("fig3_printing_point_illusion", function() {
    par(mfrow = c(1, 3), mar = c(4.2, 4.2, 3, 1), mgp = c(2.4, .8, 0))
    for (f in fs) {
      coh <- sub("^08_(.+)_100\\.tsv$", "\\1", basename(f))
      D <- read.table(f, header = TRUE, sep = "\t")
      cols <- paste0("bias_k", 2:19); cols <- cols[cols %in% names(D)]
      m <- colMeans(D[, cols, drop = FALSE], na.rm = TRUE)
      k <- 2:(1 + length(m))
      plot(k, m, type = "o", pch = 16, col = "#2471a3", xlab = "pool size k", ylab = "selection bias (C-index)",
           main = sprintf("%s (n=%d)", coh, nrow(D)), cex.main = 1.1)
      # 只印 7 个点会看到的样子
      isp <- k %in% c(2, 3, 5, 8, 12, 16, 19)
      points(k[isp], m[isp], pch = 1, col = "#c0392b", cex = 1.7, lwd = 2)
      abline(v = 14, lty = 3, col = "grey40")
      mtext("Ridge_1se enters", side = 3, at = 14, line = -.2, cex = .7, col = "grey30")
      if (k[1] == 2) legend("topleft", c("full curve", "sparse printout (k=2,3,5,8,12,16,19)"),
                            pch = c(16, 1), col = c("#2471a3", "#c0392b"), bty = "n", cex = .75)
    }
    mtext("Fig 3  Sampling only 7 of the 18 k values creates a phantom jump at k=16", outer = TRUE, line = -1.4, font = 2)
  }, w = 13, h = 4.6)
}

## ---------- Fig 4：进入效应排名 ----------
f <- file.path(TAB, "11_candidate_features.tsv")
if (file.exists(f)) {
  FT <- read.table(f, header = TRUE, sep = "\t")
  cohs <- unique(FT$cohort)
  redraw("fig4_entry_effects", function() {
    par(mfrow = c(1, length(cohs)), mar = c(9.5, 4.2, 3, .6), mgp = c(2.4, .8, 0))
    for (cn in cohs) {
      s <- FT[FT$cohort == cn, ]; s <- s[order(s$entry_mean), ]
      bp <- barplot(s$entry_mean, names.arg = NA, las = 2, col = LINCOL[s$family],
                    border = NA, ylab = "entry effect (dC-index)",
                    main = sprintf("%s (n=%d)", cn, s$n_gs[1]),
                    ylim = range(c(0, s$entry_mean)) * c(1.35, 1.05), cex.main = 1.1)
      text(bp, par("usr")[3] - 1e-3, labels = s$target, srt = 60, adj = 1,
           xpd = TRUE, cex = .62)
      abline(h = 0, col = "grey40")
      if (cn == cohs[1]) legend("topleft", c("linear", "tree ensemble"), fill = LINCOL,
                                border = NA, bty = "n", cex = .8)
    }
    mtext("Fig 4  Per-candidate entry effect (change in bias when the candidate joins the pool)",
          outer = TRUE, line = -1.6, font = 2)
  }, w = 4 * length(cohs), h = 5.6)

  ## ---------- Fig 5：机制散点（2x3）----------
  P <- c("mean_dval","sd_dval","winrate","cor_vt","optimism","tail_opt")
  LBL <- c(mean_dval="mean dval (observable)", sd_dval="sd of dval (observable)",
           winrate="win rate (observable)", cor_vt="cor(dval,dtest) (hidden)",
           optimism="mean optimism (hidden)", tail_opt="upper-tail optimism (hidden)")
  cohcol <- setNames(rainbow(length(cohs)), cohs)
  redraw("fig5_mechanism", function() {
    par(mfrow = c(2, 3), mar = c(4.4, 4.4, 2.6, 1), mgp = c(2.4, .8, 0))
    for (v in P) {
      plot(FT[[v]], FT$entry_mean, pch = 16, col = cohcol[FT$cohort], cex = 1.05,
           xlab = LBL[v], ylab = "entry effect", main = "")
      abline(lm(entry_mean ~ FT[[v]], data = FT), col = "grey30", lwd = 1.6)
      ct <- cor.test(FT[[v]], FT$entry_mean)
      legend("topleft", sprintf("r = %+.2f\np = %.1e", ct$estimate, ct$p.value),
             bty = "n", cex = .85)
    }
    mtext("Fig 5  What predicts a candidate's contribution to selection bias?",
          outer = TRUE, line = -1.3, font = 2)
  }, w = 11, h = 7)
}

## ---------- Fig 6：家族对比 ----------
f <- file.path(TAB, "10_tests.tsv")
if (file.exists(f)) {
  T0 <- read.table(f, header = TRUE, sep = "\t")
  T0 <- T0[grepl("^同样 k=", T0$contrast), ]
  if (nrow(T0)) {
    redraw("fig6_family_contrast", function() {
      s <- T0[order(T0$mean_diff), ]
      y <- seq_len(nrow(s))
      par(mar = c(4.4, 4.6, 3, 1), mgp = c(2.6, .8, 0))
      plot(s$mean_diff, y, xlim = range(c(s$ci_lo, s$ci_hi, 0)) * 1.1,
           pch = 16, cex = 1.4, yaxt = "n", xlab = "tree-only - linear-only (same k, dC-index)",
           ylab = "", main = "Fig 6  Does the algorithm family drive the bias?")
      axis(2, at = y, labels = s$cohort, las = 1)
      segments(s$ci_lo, y, s$ci_hi, y, lwd = 2)
      abline(v = 0, lty = 2, col = "grey40")
      text(s$mean_diff, y + .28, sprintf("p = %.1e", s$t_p), cex = .72, col = "grey25")
    }, w = 7, h = 4.5)
  }
}
cat("图表输出目录:", FIG, "\n")
