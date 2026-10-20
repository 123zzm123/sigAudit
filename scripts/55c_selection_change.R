#!/usr/bin/env Rscript
# ============================================================
# 55c —— 53 组 B 并列规则的差异量化（旧规则只在这里出现）
#
# 53 本身只使用 make_prio()，与 A/C 组一致，不复现任何旧规则。
# 本脚本从已保存的预测重新计算两种并列规则下的赢家与 top-five，
# 写出 53_selection_change.tsv（S40），并核对现行规则的结果与
# 53_brca_signature_tasks.tsv 完全一致。
#
# 旧规则指 26b 的写法：set.seed(seed); runif(length(mods)) 按 mods 当前顺序
# 分配随机数。它只用于回答“改规则会不会改变结论”，不参与任何报告数值。
# ============================================================
suppressMessages(library(survival))
options(stringsAsFactors = FALSE)

PROJ <- Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
TAB  <- file.path(PROJ, "results/tables")
source(file.path(PROJ, "scripts/_lib_analysis.R"))

PRED <- readRDS(file.path(PROJ, "results/main/53_brca_signature_predictions.rds"))
B <- read.table(file.path(TAB, "53_brca_signature_tasks.tsv"), sep = "\t", header = TRUE,
                quote = "", stringsAsFactors = FALSE)

SEL <- list()
for (p in PRED) {
  mods <- colnames(p$val_pred)
  V <- vapply(mods, function(m) ci_direct(p$t_val, p$e_val, p$val_pred[, m]), numeric(1)); names(V) <- mods
  Tt <- vapply(mods, function(m) ci_direct(p$t_test, p$e_test, p$test_pred[, m]), numeric(1)); names(Tt) <- mods
  if (!all(is.finite(V)) || !all(is.finite(Tt))) next

  prio_new <- make_prio(mods, p$seed)
  set.seed(p$seed); prio_old <- runif(length(mods)); names(prio_old) <- mods

  wn <- mods[pick_max(V, prio_new[mods])]; wo <- mods[pick_max(V, prio_old[mods])]
  t5n <- mods[order_rule(V, prio_new[mods])[seq_len(min(5, length(mods)))]]
  t5o <- mods[order_rule(V, prio_old[mods])[seq_len(min(5, length(mods)))]]
  ct5n <- ci_direct(p$t_test, p$e_test, ens_pred(p$val_pred, p$test_pred, t5n)$t)
  ct5o <- ci_direct(p$t_test, p$e_test, ens_pred(p$val_pred, p$test_pred, t5o)$t)

  SEL[[length(SEL) + 1]] <- data.frame(signature = p$signature, seed = p$seed,
    winner_old = wo, winner_new = wn, winner_changed = wo != wn,
    top5_old = paste(t5o, collapse = "|"), top5_new = paste(t5n, collapse = "|"),
    top5_changed = !setequal(t5o, t5n), dCt_argmax = Tt[[wn]] - Tt[[wo]],
    dCt_top5 = ct5n - ct5o, dGv = (V[[wn]] - V[[wo]]),
    n_tied_max = sum(V == max(V)))
}

SEL <- do.call(rbind, SEL)
write.table(SEL, file.path(TAB, "53_selection_change.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## 现行规则与 53 的表逐值对齐（防止两处实现漂移）
mx <- 0
for (p in PRED) {
  j <- which(B$signature == p$signature & B$seed == p$seed)
  if (!length(j)) next
  mods <- colnames(p$val_pred)
  V <- vapply(mods, function(m) ci_direct(p$t_val, p$e_val, p$val_pred[, m]), numeric(1)); names(V) <- mods
  prio <- make_prio(mods, p$seed)
  w <- mods[pick_max(V, prio[mods])]
  t5 <- mods[order_rule(V, prio[mods])[seq_len(min(5, length(mods)))]]
  ct5 <- ci_direct(p$t_test, p$e_test, ens_pred(p$val_pred, p$test_pred, t5)$t)
  mx <- max(mx, abs(ci_direct(p$t_test, p$e_test, p$test_pred[, w]) - B$Ct_argmax[j]),
            abs(ct5 - B$Ct_top5[j]))
}

cat("任务：", nrow(SEL), "\n")
cat("旧规则改变赢家：", sum(SEL$winner_changed), "；改变 top-five 成员：", sum(SEL$top5_changed), "\n")
cat("现行规则与 53 表最大绝对差：", format(mx, scientific = TRUE), "\n")
cat("赢家标签变化但测试 C 差为 0 的任务：",
    sum(SEL$winner_changed & abs(SEL$dCt_argmax) < 1e-12), "\n")
cat("最大 |dCt_top5|：", format(max(abs(SEL$dCt_top5)), digits = 6), "\n")
