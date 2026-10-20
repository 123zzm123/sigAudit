#!/usr/bin/env Rscript
# ============================================================
# 52 —— Uno/IPCW 一致性敏感性（修订期探索性补充）
#
# 只替换“测试集评价器”，不改变任何选择规则：
#   * 验证集选择仍用 Harrell C（与原分析完全一致）
#   * 测试集同时报告 Harrell C 与 Uno 的 IPCW C（tau = 60 月）
# 五个估计量：固定 Lasso（参照）、固定 Ridge、固定 Enet0.5、
#   argmax_val、top5_ensemble
#
# 汇总层级（与原主分析一致，2026-10-20 第二轮修正）：
#   任务 -> 队列×seed 均值 -> 队列内对各 seed 等权 -> 各队列等权
#   METABRIC：任务 -> 签名×seed 均值 -> 签名内五划分等权 -> 六签名等权
#
# 支持检查在**队列×seed**层面执行（G_train(tau) >= 0.05 且训练随访覆盖 tau）；
# 未通过的 seed 整体不进入分析，而不是在任务内静默丢弃事件。
# ============================================================
suppressMessages({ library(survival) })
options(stringsAsFactors = FALSE)

PROJ <- Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
TAB  <- file.path(PROJ, "results/tables")
con  <- file(file.path(PROJ, "logs", "52_uno_sensitivity.log"), "wt")
say  <- function(...) { m <- paste0(...); cat(m, "\n", sep = ""); cat(m, "\n", sep = "", file = con); flush(con) }

source(file.path(PROJ, "scripts/_lib_analysis.R"))   # ci_direct / pick_max / order_rule / make_prio / ens_pred
source(file.path(PROJ, "scripts/_lib_uno.R"))        # gfun_train / pair_prep / uno_from_pairs / har_trunc

TAU    <- 60
G_MIN  <- 0.05
MIN_EV <- 10
OS_COH <- c("BRCA", "COAD", "KIRC", "LUAD", "MESO", "STAD")

fs <- sort(list.files(file.path(PROJ, "results/main"),
                      pattern = "^19_[A-Z]+_[0-9]{8}_100\\.rds$", full.names = TRUE))

TASKS <- list(); SUPP <- list(); CHECK <- list(); ATTR <- list()

for (f in fs) {
  R <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(R) || !length(R$genesets)) next
  coh <- R$cohort
  if (!(coh %in% OS_COH)) { say("[跳过] ", coh, "（", R$endpoint, "，非 OS）"); next }
  G <- R$genesets
  cnt <- table(unlist(lapply(G, function(z) colnames(z$val_pred))))
  keep_mods <- names(cnt)[cnt / length(G) >= 0.90]
  if (!("Lasso" %in% keep_mods)) { say("[跳过] ", coh, " 无参照"); next }
  prio0 <- make_prio(keep_mods, R$seed)
  z0 <- G[[1]]
  Gf <- gfun_train(z0$t_tr, z0$e_tr)
  g_tau <- unname(Gf(TAU))
  max_fu <- max(z0$t_tr)
  supported <- is.finite(g_tau) && g_tau >= G_MIN && max_fu >= TAU

  n_eval <- 0; n_ev_excl <- 0
  if (!supported) {
    say("[不支持] ", coh, " seed ", R$seed, "：G_train(60)=", round(g_tau, 4), "，训练随访上限=", round(max_fu, 1))
  } else {
    say(""); say("===== ", coh, " seed ", R$seed, " (", R$endpoint, ") G_train(60)=",
        sprintf("%.4f", g_tau), " =====")
    for (z in G) {
      if (!all(keep_mods %in% colnames(z$val_pred))) next
      Pv <- z$val_pred[, keep_mods, drop = FALSE]; Pt <- z$test_pred[, keep_mods, drop = FALSE]
      if (!all(is.finite(Pv)) || !all(is.finite(Pt))) next
      V <- vapply(keep_mods, function(m) ci_direct(z$tva, z$eva, Pv[, m]), numeric(1)); names(V) <- keep_mods
      if (!all(is.finite(V))) next
      n_eval <- n_eval + 1

      pp <- pair_prep(z$tte, z$ete, TAU, Gf, G_MIN)
      if (pp$n_ev < MIN_EV) { n_ev_excl <- n_ev_excl + 1; next }

      winner <- keep_mods[pick_max(V, prio0[keep_mods])]
      top5   <- keep_mods[order_rule(V, prio0[keep_mods])[seq_len(min(5, length(keep_mods)))]]
      ens    <- ens_pred(Pv, Pt, top5)

      est <- list(
        lasso  = as.numeric(Pt[, "Lasso"]),
        ridge  = if ("Ridge" %in% keep_mods) as.numeric(Pt[, "Ridge"]) else NULL,
        enet05 = if ("Enet0.5" %in% keep_mods) as.numeric(Pt[, "Enet0.5"]) else NULL,
        argmax = as.numeric(Pt[, winner]),
        top5   = if (!is.null(ens)) ens$t else NULL
      )
      u <- vapply(est, function(s) if (is.null(s)) NA_real_ else uno_from_pairs(s, pp), numeric(1))
      h <- vapply(est, function(s) if (is.null(s)) NA_real_ else ci_direct(z$tte, z$ete, s), numeric(1))
      h60 <- vapply(est, function(s) if (is.null(s)) NA_real_ else har_trunc(z$tte, z$ete, s, TAU), numeric(1))

      TASKS[[length(TASKS) + 1]] <- data.frame(
        cohort = coh, seed = R$seed, geneset = z$name, n_test = length(z$tte),
        n_event_tau = pp$n_ev, cens_frac = mean(z$ete == 0),
        G_tau = g_tau, winner = winner,
        uno_lasso = u[["lasso"]], uno_ridge = u[["ridge"]], uno_enet05 = u[["enet05"]],
        uno_argmax = u[["argmax"]], uno_top5 = u[["top5"]],
        har_lasso = h[["lasso"]], har_ridge = h[["ridge"]], har_enet05 = h[["enet05"]],
        har_argmax = h[["argmax"]], har_top5 = h[["top5"]],
        har60_lasso = h60[["lasso"]], har60_ridge = h60[["ridge"]], har60_enet05 = h60[["enet05"]],
        har60_argmax = h60[["argmax"]], har60_top5 = h60[["top5"]])

      ## 删失分布估计来源的差异（train-G vs 分析数据 G），只作参考
      if (n_eval <= 3) {
        for (nm in c("lasso", "argmax")) {
          cc <- tryCatch(as.numeric(concordance(Surv(z$tte, z$ete) ~ est[[nm]], reverse = TRUE,
                                               timewt = "n/G2", ymax = TAU)$concordance), error = function(e) NA_real_)
          CHECK[[length(CHECK) + 1]] <- data.frame(cohort = coh, seed = R$seed, geneset = z$name,
            estimator = nm, uno_trainG = u[[nm]], uno_analysisG = cc, diff = u[[nm]] - cc)
        }
      }
    }
  }

  ATTR[[length(ATTR) + 1]] <- data.frame(cohort = coh, seed = R$seed,
    n_planned = R$n_planned, n_gene_sets = length(G), n_jointly_evaluable = n_eval,
    n_excluded_event_rule = n_ev_excl, n_analysed = n_eval - n_ev_excl)
  SUPP[[length(SUPP) + 1]] <- data.frame(cohort = coh, seed = R$seed,
    n_test = length(z0$tte), n_train = length(z0$t_tr),
    ev_test = sum(z0$ete), ev_test_before_tau = sum(z0$ete == 1 & z0$tte < TAU),
    cens_test = mean(z0$ete == 0), G_train_tau = g_tau,
    max_followup_train = max_fu, supported = supported)
}

## ---------- METABRIC ----------
mb <- tryCatch(readRDS(file.path(PROJ, "results/metabric/31_metabric_predictions.rds")), error = function(e) NULL)
if (!is.null(mb)) {
  say(""); say("===== METABRIC（六签名 × 五划分）=====")
  for (z in mb) {
    mods <- z$model
    Gf <- gfun_train(z$t_train, z$e_train)
    g_tau <- unname(Gf(TAU)); max_fu <- max(z$t_train)
    supported <- is.finite(g_tau) && g_tau >= G_MIN && max_fu >= TAU
    SUPP[[length(SUPP) + 1]] <- data.frame(cohort = "METABRIC", seed = z$seed,
      n_test = length(z$t_test), n_train = length(z$t_train),
      ev_test = sum(z$e_test), ev_test_before_tau = sum(z$e_test == 1 & z$t_test < TAU),
      cens_test = mean(z$e_test == 0), G_train_tau = g_tau,
      max_followup_train = max_fu, supported = supported)
    if (!supported) next
    V <- vapply(mods, function(m) ci_direct(z$t_val, z$e_val, z$val_pred[, m]), numeric(1)); names(V) <- mods
    if (!all(is.finite(V))) next
    prio <- make_prio(mods, z$seed)
    winner <- mods[pick_max(V, prio[mods])]
    top5 <- mods[order_rule(V, prio[mods])[seq_len(min(5, length(mods)))]]
    ens <- ens_pred(z$val_pred, z$test_pred, top5)
    pp <- pair_prep(z$t_test, z$e_test, TAU, Gf, G_MIN)
    if (pp$n_ev < MIN_EV) next
    est <- list(lasso = as.numeric(z$test_pred[, "Lasso"]),
                ridge = as.numeric(z$test_pred[, "Ridge"]),
                enet05 = as.numeric(z$test_pred[, "Enet0.5"]),
                argmax = as.numeric(z$test_pred[, winner]),
                top5 = ens$t)
    u <- vapply(est, function(s) uno_from_pairs(s, pp), numeric(1))
    h <- vapply(est, function(s) ci_direct(z$t_test, z$e_test, s), numeric(1))
    h60 <- vapply(est, function(s) har_trunc(z$t_test, z$e_test, s, TAU), numeric(1))
    TASKS[[length(TASKS) + 1]] <- data.frame(
      cohort = "METABRIC", seed = z$seed, geneset = z$signature, n_test = length(z$t_test),
      n_event_tau = pp$n_ev, cens_frac = mean(z$e_test == 0), G_tau = g_tau, winner = winner,
      uno_lasso = u[["lasso"]], uno_ridge = u[["ridge"]], uno_enet05 = u[["enet05"]],
      uno_argmax = u[["argmax"]], uno_top5 = u[["top5"]],
      har_lasso = h[["lasso"]], har_ridge = h[["ridge"]], har_enet05 = h[["enet05"]],
      har_argmax = h[["argmax"]], har_top5 = h[["top5"]],
      har60_lasso = h60[["lasso"]], har60_ridge = h60[["ridge"]], har60_enet05 = h60[["enet05"]],
      har60_argmax = h60[["argmax"]], har60_top5 = h60[["top5"]])
  }
  ATTR[[length(ATTR) + 1]] <- data.frame(cohort = "METABRIC", seed = NA_integer_,
    n_planned = length(mb), n_gene_sets = length(mb), n_jointly_evaluable = length(mb),
    n_excluded_event_rule = 0, n_analysed = sum(sapply(TASKS, function(d) d$cohort[1] == "METABRIC")))
}

T <- do.call(rbind, TASKS)
S <- do.call(rbind, SUPP)
A <- do.call(rbind, ATTR)

## ---------- 逐任务差异（字段名即比较对象） ----------
T$d_uno_argmax_lasso <- T$uno_argmax - T$uno_lasso
T$d_uno_top5_lasso   <- T$uno_top5   - T$uno_lasso
T$d_uno_ridge_lasso  <- T$uno_ridge  - T$uno_lasso
T$d_uno_enet_lasso   <- T$uno_enet05 - T$uno_lasso
T$d_uno_argmax_ridge <- T$uno_argmax - T$uno_ridge
T$d_uno_top5_ridge   <- T$uno_top5   - T$uno_ridge
T$d_uno_top5_argmax  <- T$uno_top5   - T$uno_argmax
T$d_uno_argmax_enet  <- T$uno_argmax - T$uno_enet05
T$d_uno_top5_enet    <- T$uno_top5   - T$uno_enet05
T$d_har_argmax_lasso <- T$har_argmax - T$har_lasso
T$d_har_top5_lasso   <- T$har_top5   - T$har_lasso
T$d_har_argmax_ridge <- T$har_argmax - T$har_ridge
T$d_har60_argmax_lasso <- T$har60_argmax - T$har60_lasso
T$d_har60_top5_lasso   <- T$har60_top5   - T$har60_lasso

write.table(T, file.path(TAB, "52_uno_task.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(S, file.path(TAB, "52_uno_support.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(A, file.path(TAB, "52_uno_attrition.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---------- 分层汇总 ----------
dcols <- c("d_uno_argmax_lasso", "d_uno_top5_lasso", "d_uno_ridge_lasso", "d_uno_enet_lasso",
           "d_uno_argmax_ridge", "d_uno_top5_ridge", "d_uno_top5_argmax",
           "d_uno_argmax_enet", "d_uno_top5_enet",
           "d_har_argmax_lasso", "d_har_top5_lasso", "d_har_argmax_ridge",
           "d_har60_argmax_lasso", "d_har60_top5_lasso")

tc <- T[T$cohort != "METABRIC", ]
cell <- aggregate(tc[dcols], by = tc[c("cohort", "seed")], FUN = function(x) mean(x, na.rm = TRUE))
cell$n_task <- as.numeric(table(tc[c("cohort", "seed")])[cbind(cell$cohort, as.character(cell$seed))])
coh  <- aggregate(cell[dcols], by = cell["cohort"], FUN = function(x) mean(x, na.rm = TRUE))
coh$n_cell <- as.numeric(table(cell$cohort)[coh$cohort])
coh$n_task <- as.numeric(table(tc$cohort)[coh$cohort])
write.table(cell, file.path(TAB, "52_uno_by_cell.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(coh,  file.path(TAB, "52_uno_by_cohort.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

mbt <- T[T$cohort == "METABRIC", ]
mcell <- aggregate(mbt[dcols], by = mbt[c("geneset", "seed")], FUN = function(x) mean(x, na.rm = TRUE))
msig  <- aggregate(mcell[dcols], by = mcell["geneset"], FUN = function(x) mean(x, na.rm = TRUE))
msig$n_task <- as.numeric(table(mbt$geneset)[msig$geneset])
write.table(msig, file.path(TAB, "52_uno_by_signature.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

mkrow <- function(scope, n_unit, ntask, vals) {
  out <- data.frame(scope = scope, n_unit = n_unit, n_task = ntask, check.names = FALSE)
  for (k in names(vals)) out[[k]] <- as.numeric(vals[[k]])
  out
}
os <- coh[coh$cohort %in% OS_COH, ]
agg <- mkrow("6 OS cohorts (macro)", nrow(os), sum(os$n_task), as.list(colMeans(os[dcols], na.rm = TRUE)))
agg <- rbind(agg, mkrow("METABRIC (macro)", nrow(msig), sum(msig$n_task), as.list(colMeans(msig[dcols], na.rm = TRUE))))
write.table(agg, file.path(TAB, "52_uno_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

if (length(CHECK)) write.table(do.call(rbind, CHECK), file.path(TAB, "52_uno_check.tsv"),
                               sep = "\t", quote = FALSE, row.names = FALSE)

## ---------- 报告 ----------
fmt <- function(v) sprintf("%+.6f", v)
say("")
say("================ Uno/IPCW 敏感性（tau = ", TAU, " 月；分层汇总）================")
say(sprintf("  %-26s %5s %6s %12s %12s %12s %12s", "范围", "单元", "任务",
            "argmax-Lasso", "top5-Lasso", "argmax-Ridge", "top5-Ridge"))
for (i in seq_len(nrow(agg))) {
  say(sprintf("  %-26s %5d %6d %12s %12s %12s %12s", agg$scope[i], agg$n_unit[i], agg$n_task[i],
      fmt(agg$d_uno_argmax_lasso[i]), fmt(agg$d_uno_top5_lasso[i]),
      fmt(agg$d_uno_argmax_ridge[i]), fmt(agg$d_uno_top5_ridge[i])))
}
say("")
say(sprintf("  支持范围：G_train(60) %.6f–%.6f；tau 前测试事件 %d–%d（TCGA）",
    min(S$G_train_tau[S$cohort != "METABRIC"]), max(S$G_train_tau),
    min(S$ev_test_before_tau[S$cohort != "METABRIC"]), max(S$ev_test_before_tau[S$cohort != "METABRIC"])))
say(sprintf("  METABRIC：G_train(60) %.6f–%.6f；tau 前事件 %d–%d",
    min(S$G_train_tau[S$cohort == "METABRIC"]), max(S$G_train_tau[S$cohort == "METABRIC"]),
    min(S$ev_test_before_tau[S$cohort == "METABRIC"]), max(S$ev_test_before_tau[S$cohort == "METABRIC"])))
say("  说明：验证集仍用 Harrell C 选择；Uno 只换测试评价器；未用 Uno 重新选模，未跨指标相减。")
close(con)
