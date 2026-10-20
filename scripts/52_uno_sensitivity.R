#!/usr/bin/env Rscript
# ============================================================
# 52 —— Uno/IPCW 一致性敏感性（修订期探索性补充）
#
# 目的：只替换“测试集评价器”，不改变任何选择规则。
#   * 验证集选择仍用 Harrell C（与原分析完全一致）
#   * 测试集同时报告 Harrell C 与 Uno 的 IPCW C（tau = 60 月）
# 对比的五个估计量：固定 Lasso（参照）、固定 Ridge、固定 Enet0.5、
#   argmax_val（验证集最高分的候选）、top5_ensemble（验证集前五的标准化均值）
#
# 冻结配置（运行前写入 config/revision_20261020_uno.md）：
#   * tau = 60 月（5 年临床窗口），不按测试表现调整
#   * 删失分布 G 由**训练折**估计：Surv(t_tr, 1-e_tr) 的 KM，在测试事件时刻取值
#   * 权重 1/G(T_i)^2；可比对要求 T_i < T_j 且 T_i < tau 且 δ_i = 1
#   * 可行性：G_tr(tau) >= 0.05 且该任务 tau 前事件 >= 10；否则记为不可评估
#   * 同一 tau 下另给截断 Harrell（ymax = tau），用于区分“换指标”与“截断”两个来源
#   * 只用 6 个 OS 队列（LGG/PAAD 为 PFS，事件构成未核实，不纳入）
#
# 不做：不用 Uno 重新选模；不把 Harrell 的 Gv 与 Uno 的 Gt 相减；不扫多个 tau。
# ============================================================
suppressMessages({ library(survival) })
options(stringsAsFactors = FALSE)

PROJ <- Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
TAB  <- file.path(PROJ, "results/tables")
LOG  <- file.path(PROJ, "logs/52_uno_sensitivity.log")
con  <- file(file.path(PROJ, "logs", "52_uno_sensitivity.log"), "wt")
say  <- function(...) { m <- paste0(...); cat(m, "\n", sep = ""); cat(m, "\n", sep = "", file = con); flush(con) }

source(file.path(PROJ, "scripts/_lib_analysis.R"))   # ci_direct / pick_max / order_rule / make_prio / ens_pred

TAU      <- 60          # 月
G_MIN    <- 0.05        # 训练侧删失分布在该任务 tau 处的支持下限
MIN_EV   <- 10          # tau 前最少事件数
OS_COH   <- c("BRCA", "COAD", "KIRC", "LUAD", "MESO", "STAD")

## ---------- IPCW 工具 ----------
gfun_train <- function(t_tr, e_tr) {
  sf <- survfit(Surv(t_tr, 1 - e_tr) ~ 1)
  tm <- sf$time; sv <- sf$surv
  function(x) { i <- findInterval(x, tm); ifelse(i == 0, 1, sv[pmax(i, 1)]) }
}

## 预先算好“事件 i（T_i<tau）与更晚的 j”的可比对集，供同一任务的所有模型复用
pair_prep <- function(t, e, tau, G) {
  Gi <- G(t)
  ev <- which(e == 1 & t < tau & is.finite(Gi) & Gi >= G_MIN)
  if (!length(ev)) return(list(i = integer(0), j = integer(0), w = numeric(0), n_ev = 0L))
  ii <- integer(0); jj <- integer(0); ww <- numeric(0)
  for (i in ev) {
    j <- which(t > t[i])
    if (!length(j)) next
    ii <- c(ii, rep(i, length(j))); jj <- c(jj, j); ww <- c(ww, rep(1 / Gi[i]^2, length(j)))
  }
  list(i = ii, j = jj, w = ww, n_ev = length(ev))
}

uno_from_pairs <- function(s, pp) {
  if (!length(pp$w)) return(NA_real_)
  d <- s[pp$i] - s[pp$j]
  sum(pp$w * (d > 0) + 0.5 * pp$w * (d == 0)) / sum(pp$w)
}

har_trunc <- function(t, e, s, tau) {
  ok <- is.finite(s) & is.finite(t) & !is.na(e)
  if (sum(ok) < 10) return(NA_real_)
  tryCatch(as.numeric(concordance(Surv(t[ok], e[ok]) ~ s[ok], reverse = TRUE,
                                  timewt = "n", ymax = tau)$concordance), error = function(x) NA_real_)
}

## ---------- 读入 ----------
fs <- list.files(file.path(PROJ, "results/main"), pattern = "^19_[A-Z]+_[0-9]{8}_100\\.rds$", full.names = TRUE)
fs <- sort(fs)
say("找到 ", length(fs), " 个队列产物；Uno 只纳入 6 个 OS 队列：", paste(OS_COH, collapse = ", "))

TASKS <- list(); SUPP <- list(); CHECK <- list(); FROZEN <- list()

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
  Gf <- gfun_train(G[[1]]$t_tr, G[[1]]$e_tr)

  say(""); say("===== ", coh, " seed ", R$seed, " (", R$endpoint, ") =====")
  n_task <- 0
  for (z in G) {
    if (!all(keep_mods %in% colnames(z$val_pred))) next
    Pv <- z$val_pred[, keep_mods, drop = FALSE]; Pt <- z$test_pred[, keep_mods, drop = FALSE]
    if (!all(is.finite(Pv)) || !all(is.finite(Pt))) next
    V <- vapply(keep_mods, function(m) ci_direct(z$tva, z$eva, Pv[, m]), numeric(1)); names(V) <- keep_mods
    if (!all(is.finite(V))) next
    winner <- keep_mods[pick_max(V, prio0[keep_mods])]
    top5   <- keep_mods[order_rule(V, prio0[keep_mods])[seq_len(min(5, length(keep_mods)))]]
    ens    <- ens_pred(Pv, Pt, top5)

    pp <- pair_prep(z$tte, z$ete, TAU, Gf)
    if (pp$n_ev < MIN_EV) next

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
      G_tau = unname(Gf(TAU)), winner = winner,
      uno_lasso = u[["lasso"]], uno_ridge = u[["ridge"]], uno_enet05 = u[["enet05"]],
      uno_argmax = u[["argmax"]], uno_top5 = u[["top5"]],
      har_lasso = h[["lasso"]], har_ridge = h[["ridge"]], har_enet05 = h[["enet05"]],
      har_argmax = h[["argmax"]], har_top5 = h[["top5"]],
      har60_lasso = h60[["lasso"]], har60_ridge = h60[["ridge"]], har60_enet05 = h60[["enet05"]],
      har60_argmax = h60[["argmax"]], har60_top5 = h60[["top5"]])
    n_task <- n_task + 1

    ## 交叉核对：手写 IPCW 与 survival::concordance(timewt="n/G2", ymax) 的前 3 个任务
    if (n_task <= 3) {
      for (nm in c("lasso", "argmax")) {
        cc <- tryCatch(as.numeric(concordance(Surv(z$tte, z$ete) ~ est[[nm]], reverse = TRUE,
                                             timewt = "n/G2", ymax = TAU)$concordance), error = function(e) NA_real_)
        CHECK[[length(CHECK) + 1]] <- data.frame(cohort = coh, seed = R$seed, geneset = z$name,
          estimator = nm, uno_manual_trainG = u[[nm]], uno_survival_dataG = cc,
          diff = u[[nm]] - cc)
      }
    }
  }
  say("  可评估任务（tau 前事件 >= ", MIN_EV, "）：", n_task, " / ", length(G))
  SUPP[[length(SUPP) + 1]] <- data.frame(cohort = coh, endpoint = R$endpoint, seed = R$seed,
    n_test = length(G[[1]]$tte), ev_test = sum(G[[1]]$ete),
    ev_test_before_tau = sum(G[[1]]$ete == 1 & G[[1]]$tte < TAU),
    cens_test = mean(G[[1]]$ete == 0), G_train_tau = unname(Gf(TAU)),
    supported = unname(Gf(TAU)) >= G_MIN)
}

## ---------- METABRIC ----------
mb <- tryCatch(readRDS(file.path(PROJ, "results/metabric/31_metabric_predictions.rds")), error = function(e) NULL)
if (!is.null(mb)) {
  say(""); say("===== METABRIC（六签名 × 五划分）=====")
  n_task <- 0
  for (z in mb) {
    mods <- z$model
    V <- vapply(mods, function(m) ci_direct(z$t_val, z$e_val, z$val_pred[, m]), numeric(1)); names(V) <- mods
    if (!all(is.finite(V))) next
    prio <- make_prio(mods, z$seed)
    winner <- mods[pick_max(V, prio[mods])]
    top5 <- mods[order_rule(V, prio[mods])[seq_len(min(5, length(mods)))]]
    ens <- ens_pred(z$val_pred, z$test_pred, top5)
    Gf <- gfun_train(z$t_train, z$e_train)
    pp <- pair_prep(z$t_test, z$e_test, TAU, Gf)
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
      n_event_tau = pp$n_ev, cens_frac = mean(z$e_test == 0), G_tau = unname(Gf(TAU)),
      winner = winner,
      uno_lasso = u[["lasso"]], uno_ridge = u[["ridge"]], uno_enet05 = u[["enet05"]],
      uno_argmax = u[["argmax"]], uno_top5 = u[["top5"]],
      har_lasso = h[["lasso"]], har_ridge = h[["ridge"]], har_enet05 = h[["enet05"]],
      har_argmax = h[["argmax"]], har_top5 = h[["top5"]],
      har60_lasso = h60[["lasso"]], har60_ridge = h60[["ridge"]], har60_enet05 = h60[["enet05"]],
      har60_argmax = h60[["argmax"]], har60_top5 = h60[["top5"]])
    n_task <- n_task + 1
  }
  say("  可评估任务：", n_task, " / ", length(mb))
  SUPP[[length(SUPP) + 1]] <- data.frame(cohort = "METABRIC", endpoint = "OS", seed = NA_integer_,
    n_test = length(mb[[1]]$t_test), ev_test = sum(mb[[1]]$e_test),
    ev_test_before_tau = sum(mb[[1]]$e_test == 1 & mb[[1]]$t_test < TAU),
    cens_test = mean(mb[[1]]$e_test == 0),
    G_train_tau = unname(gfun_train(mb[[1]]$t_train, mb[[1]]$e_train)(TAU)), supported = TRUE)
}

nms <- lapply(TASKS, names)
if (length(unique(lapply(nms, paste, collapse = "|"))) > 1) {
  print(table(vapply(nms, function(x) paste(x, collapse = "|"), character(1))))
  stop("TASKS 列名不一致")
}
T <- do.call(rbind, TASKS)
S <- do.call(rbind, SUPP)

## ---------- 逐任务差异 ----------
T$d_uno_argmax_lasso  <- T$uno_argmax - T$uno_lasso
T$d_uno_top5_lasso    <- T$uno_top5   - T$uno_lasso
T$d_uno_ridge_lasso   <- T$uno_ridge  - T$uno_lasso
T$d_uno_enet_lasso    <- T$uno_enet05 - T$uno_lasso
T$d_uno_argmax_ridge  <- T$uno_argmax - T$uno_ridge
T$d_uno_top5_ridge    <- T$uno_top5   - T$uno_ridge
T$d_uno_top5_argmax   <- T$uno_top5   - T$uno_argmax
T$d_har_argmax_lasso  <- T$har_argmax - T$har_lasso
T$d_har_top5_lasso    <- T$har_top5   - T$har_lasso
T$d_har_argmax_ridge  <- T$har_argmax - T$har_ridge
T$d_har60_argmax_lasso<- T$har60_argmax - T$har60_lasso
T$d_har60_top5_lasso  <- T$har60_top5   - T$har60_lasso

write.table(T, file.path(TAB, "52_uno_task.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(S, file.path(TAB, "52_uno_support.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---------- 队列层均值 → 6 OS 宏平均 ----------
dcols <- c("d_uno_argmax_lasso", "d_uno_top5_lasso", "d_uno_ridge_lasso", "d_uno_enet_lasso",
           "d_uno_argmax_ridge", "d_uno_top5_ridge", "d_uno_top5_argmax",
           "d_har_argmax_lasso", "d_har_top5_lasso", "d_har_argmax_ridge",
           "d_har60_argmax_lasso", "d_har60_top5_lasso")
cco <- aggregate(T[dcols], by = T["cohort"], FUN = function(x) mean(x, na.rm = TRUE))
cco$n_task <- as.numeric(table(T$cohort)[cco$cohort])
write.table(cco, file.path(TAB, "52_uno_by_cohort.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

os <- cco[cco$cohort %in% OS_COH, ]
mkrow <- function(scope, ncoh, ntask, vals) {
  out <- data.frame(scope = scope, n_cohort = ncoh, n_task = ntask, check.names = FALSE)
  for (k in names(vals)) out[[k]] <- as.numeric(vals[[k]])
  out
}
agg <- mkrow("6 OS cohorts (macro)", nrow(os), sum(os$n_task), as.list(colMeans(os[dcols], na.rm = TRUE)))
mbrow <- cco[cco$cohort == "METABRIC", ]
if (nrow(mbrow)) {
  agg <- rbind(agg, mkrow("METABRIC (macro)", 1, mbrow$n_task, as.list(unlist(mbrow[dcols]))))
}
write.table(agg, file.path(TAB, "52_uno_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

if (length(CHECK)) write.table(do.call(rbind, CHECK), file.path(TAB, "52_uno_check.tsv"),
                               sep = "\t", quote = FALSE, row.names = FALSE)

## ---------- 报告 ----------
fmt <- function(v) sprintf("%+.4f", v)
say("")
say("================ Uno/IPCW 敏感性（tau = ", TAU, " 月）================")
say(sprintf("  %-34s %6s %7s %9s %9s %9s %9s", "对比（Uno）", "队列", "任务",
            "argmax-Lasso", "top5-Lasso", "argmax-Ridge", "top5-Ridge"))
for (i in seq_len(nrow(agg))) {
  say(sprintf("  %-34s %6d %7d %9s %9s %9s %9s", agg$scope[i], agg$n_cohort[i], agg$n_task[i],
      fmt(agg$d_uno_argmax_lasso[i]), fmt(agg$d_uno_top5_lasso[i]),
      fmt(agg$d_uno_argmax_ridge[i]), fmt(agg$d_uno_top5_ridge[i])))
}
say("")
say("  同 tau 的截断 Harrell（用于区分“换指标”与“截断”）：")
for (i in seq_len(nrow(agg))) {
  say(sprintf("  %-34s argmax-Lasso %s   top5-Lasso %s", agg$scope[i],
      fmt(agg$d_har60_argmax_lasso[i]), fmt(agg$d_har60_top5_lasso[i])))
}
if (length(CHECK)) {
  ch <- do.call(rbind, CHECK)
  say("")
  say(sprintf("  实现核对：手写(train-G) vs survival(n/G2, data-G)  n=%d  mean diff=%+.5f  max|diff|=%.5f",
      nrow(ch), mean(ch$diff, na.rm = TRUE), max(abs(ch$diff), na.rm = TRUE)))
}
say("")
say("  说明：验证集仍用 Harrell C 选择；Uno 只替换测试集评价器。")
say("  tau=60 月固定；删失分布由训练折估计；未用 Uno 重新选模，未跨指标相减。")
close(con)
