# ============================================================
# _lib_analysis.R —— 主分析与基线敏感性共用的函数库
#
# 存在的理由（第三轮审查意见）：
#  * 审查指出 20 脚本与 R 包对同名规则给出不同定义（linear_only 是否含参照），
#    并指出"所有接口共用同一函数"是必须的。把规则/评价函数集中到一份库文件，
#    让 20（主分析）与 29（固定基线敏感性）共用，消除脚本之间的口径漂移。
#  * 审查指出 n_pairs 算错：Harrell C 的有效可比对数为
#        concordant + discordant + tied.x
#    不能把 tied.y / tied.xy 一并计入。本库已修正。
#
# 本文件不产生任何输出，只定义函数。被 20 与 29 以 source() 载入。
# ============================================================
suppressMessages(library(survival))

TREES <- c("RSF100","GBM100","GBM250","RSF250")

## ---------- C-index：直接评价，方向固定 ----------
## 用 survival::concordance(reverse=TRUE, timewt="n")，不做 Cox 重拟合，
## 因此不会从评价集自身的结局反推风险方向。
ci_direct <- function(t, e, s) {
  s <- tryCatch(as.numeric(s), error=function(x) NULL)
  if (is.null(s) || length(s) != length(t)) return(NA_real_)
  ok <- is.finite(s) & !is.na(t) & !is.na(e)
  if (sum(ok) < 10) return(NA_real_)
  tryCatch(as.numeric(concordance(Surv(t[ok],e[ok]) ~ s[ok], reverse=TRUE, timewt="n")$concordance),
           error=function(x) NA_real_)
}

## 仅用于量化"评价口径修正"的影响（旧口径 = Cox 重拟合拿 concordance[1]）
ci_old <- function(t, e, s) {
  s <- tryCatch(as.numeric(s), error=function(x) NULL)
  if (is.null(s) || length(s) != length(t)) return(NA_real_)
  ok <- is.finite(s) & !is.na(t) & !is.na(e)
  if (sum(ok) < 10) return(NA_real_)
  tryCatch(as.numeric(summary(coxph(Surv(t[ok],e[ok]) ~ s[ok]))$concordance[1]),
           error=function(x) NA_real_)
}

## ---------- 可比对生存对数（已按审查意见修正）----------
## Harrell C 的分母 = concordant + discordant + tied.x。
## 说明：tied.x 是"风险评分并列但事件时间可排序"的对，计入分母；
##       tied.y（事件时间并列）与 tied.xy（两者都并列）不属于可比对。
n_pairs <- function(t, e) {
  ok <- !is.na(t) & !is.na(e) & is.finite(t) & is.finite(e)
  if (sum(ok) < 2) return(0)
  cd <- tryCatch(concordance(Surv(t[ok],e[ok]) ~ seq_along(t[ok]), timewt="n"),
                 error=function(x) NULL)
  if (is.null(cd) || is.null(cd$count)) return(NA_real_)
  cts <- cd$count
  if (is.null(names(cts))) return(as.numeric(sum(cts, na.rm=TRUE)))
  need <- c("concordant","discordant","tied.x")
  if (all(need %in% names(cts)))
    return(as.numeric(unname(cts["concordant"]) + unname(cts["discordant"]) + unname(cts["tied.x"])))
  # 兜底：只有 tied.x 可用时按名字取，绝不把其余 tied.* 一并相加
  tx <- if ("tied.x" %in% names(cts)) unname(cts["tied.x"]) else 0
  as.numeric(unname(cts["concordant"]) + unname(cts["discordant"]) + tx)
}

## ---------- 与输入列顺序【无关】的并列处理 ----------
## 语义：prio 是一个"名字 -> 优先级"的映射，调用方保证同一候选集合下 prio 值唯一。
## 关键：候选身份由【名字】决定，不由列位置决定。交换输入列顺序不改变任何结果。
pick_max <- function(v, prio) {
  m <- suppressWarnings(max(v, na.rm=TRUE))
  if (!is.finite(m)) return(NA_integer_)
  cand <- which(is.finite(v) & v == m)
  cand[which.max(prio[cand])]
}
order_rule <- function(v, prio) order(-v, -prio, na.last = NA)

## 生成规范优先级：先对候选名排序，再按排序后的规范顺序分配随机数。
## 这样 prio 只依赖 (候选名字集合, seed)，与输入列顺序无关。
make_prio <- function(mods, seed) {
  mods <- sort(unique(mods))
  set.seed(seed)
  p <- runif(length(mods))
  stats::setNames(p, mods)
}

## ---------- 集成与规则 ----------
## 均值/标准差只用【验证集预测】估计并固定应用到测试集（禁止用测试侧信息标准化）
ens_pred <- function(Pv, Pt, cols) {
  cols <- cols[!is.na(cols)]; if (!length(cols)) return(NULL)
  Mv <- Pv[, cols, drop=FALSE]; Mt <- Pt[, cols, drop=FALSE]
  mu <- colMeans(Mv); sdv <- apply(Mv, 2, sd)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  list(v = rowMeans(sweep(sweep(Mv, 2, mu, "-"), 2, sdv, "/")),
       t = rowMeans(sweep(sweep(Mt, 2, mu, "-"), 2, sdv, "/")))
}

## 规则集合。ci_eval 提供"测试侧评价器"以便做旧口径对照（默认与验证侧同一函数）。
## 返回 out: rule -> c(Cv, Ct)；rec: rule -> 获胜者与成员。
##
## 【命名与口径，按第三轮审查意见明确】
##   linear_only = 全部【线性】候选内选优，【包含参照 Lasso】（setdiff(mods, TREES)）
##   tree_only   = 全部【树模型】候选内选优（intersect(mods, TREES)）
##   -> 二者【不是】同规模、同任务集的对称家族比较。原 S3（家族对照）主张已撤下，
##      偏离记录见 README_current。不要把 linear_only − tree_only 解释为家族效应。
rules_one <- function(z, prio, ci_eval = ci_direct, B_VOTE = 12, B_HALF = 6) {
  Pv <- z$val_pred; Pt <- z$test_pred; mods <- colnames(Pv)
  tva <- z$tva; eva <- z$eva; tte <- z$tte; ete <- z$ete
  V <- vapply(mods, function(m) ci_direct(tva, eva, Pv[,m]), numeric(1)); names(V) <- mods
  T <- vapply(mods, function(m) ci_eval(tte, ete, Pt[,m]), numeric(1)); names(T) <- mods
  lin <- setdiff(mods, TREES); tre <- intersect(mods, TREES)
  out <- list(); rec <- list()
  put <- function(rn, cv, ct, winner, members) {
    out[[rn]] <<- c(cv, ct)
    rec[[rn]] <<- data.frame(rule=rn, winner=winner, members=paste(members, collapse="|"))
  }
  if ("Lasso" %in% mods) put("none", V[["Lasso"]], T[["Lasso"]], "Lasso", "Lasso")
  # 主要规则：全池，允许 Lasso 参与
  if (length(mods)) {
    w <- mods[pick_max(V[mods], prio[mods])]; put("argmax_val", V[w], T[w], w, w)
    ord <- order_rule(V[mods], prio[mods])
    if (length(ord) >= 2) { w2 <- mods[ord[2]]; put("runner_up", V[w2], T[w2], w2, w2) }
  }
  # 敏感性：原版口径（排除参照）
  cand0 <- setdiff(mods, "Lasso")
  if (length(cand0)) { w0 <- cand0[pick_max(V[cand0], prio[cand0])]
    put("argmax_nonreference", V[w0], T[w0], w0, w0) }
  # 家族（非对称，仅作描述性；S3 主张已撤下）
  if (length(lin)) { wl <- lin[pick_max(V[lin], prio[lin])]; put("linear_only", V[wl], T[wl], wl, wl) }
  if (length(tre)) { wt <- tre[pick_max(V[tre], prio[tre])]; put("tree_only", V[wt], T[wt], wt, wt) }
  # top-k / 全池集成（全池，含 Lasso）
  for (k in c(2,3,5,10)) {
    tk <- mods[order_rule(V[mods], prio[mods])[seq_len(min(k, length(mods)))]]
    pr <- ens_pred(Pv, Pt, tk)
    if (!is.null(pr)) put(sprintf("top%d_ensemble", k), ci_direct(tva, eva, pr$v),
                          ci_eval(tte, ete, pr$t), "ensemble", tk)
  }
  pr <- ens_pred(Pv, Pt, mods)
  if (!is.null(pr)) put("all_ensemble", ci_direct(tva, eva, pr$v), ci_eval(tte, ete, pr$t), "ensemble", mods)
  # 验证集重采样类（方向固定，不做方向学习）
  set.seed(1); wins <- stats::setNames(numeric(length(mods)), mods)
  for (b in seq_len(B_VOTE)) {
    ib <- sample(length(tva), replace=TRUE); if (sum(eva[ib]) < 5) next
    vb <- vapply(mods, function(m) ci_direct(tva[ib], eva[ib], Pv[ib,m]), numeric(1))
    if (all(!is.finite(vb))) next
    jj <- pick_max(vb, prio[mods]); wins[jj] <- wins[jj] + 1
  }
  # 【修正】最终票数并列也用 prio 破并，不再用 which.max(位置优先)
  if (sum(wins) > 0) { wm <- mods[pick_max(wins[mods], prio[mods])]; put("majority_vote", V[wm], T[wm], wm, wm) }
  set.seed(1); acc <- stats::setNames(vector("list", length(mods)), mods)
  for (b in seq_len(B_HALF)) {
    pm <- sample(length(tva)); h1 <- pm[1:floor(length(tva)/2)]; h2 <- setdiff(pm, h1)
    for (m in mods) { x1 <- ci_direct(tva[h1], eva[h1], Pv[h1,m]); x2 <- ci_direct(tva[h2], eva[h2], Pv[h2,m])
      if (is.finite(x1) && is.finite(x2)) acc[[m]] <- c(acc[[m]], x1, x2) }
  }
  hmin <- vapply(acc, function(z2) if (length(z2)>=4) { k2 <- length(z2)%/%2*2
      mean(pmin(z2[seq(1,k2,2)], z2[seq(2,k2,2)])) } else NA_real_, numeric(1))
  hpen <- vapply(acc, function(z2) if (length(z2)>=4) mean(z2)-sd(z2) else NA_real_, numeric(1))
  if (any(is.finite(hmin))) { ww <- mods[pick_max(hmin, prio[mods])]; put("worst_half", V[ww], T[ww], ww, ww) }
  if (any(is.finite(hpen))) { wp <- mods[pick_max(hpen, prio[mods])]; put("penalized", V[wp], T[wp], wp, wp) }
  # 描述性参照（用测试数据选优，不是总体上界，不可实施）
  if (length(mods)) { wo <- mods[pick_max(T[mods], prio[mods])]; put("oracle_test", V[wo], T[wo], wo, wo) }
  list(out=out, rec=rec, V=V, T=T)
}

## ---------- 分层汇总（先队列内跨划分平均 -> 再对队列等权）----------
## 返回各规则在队列层级上的均值，并给出跨队列范围。
strat_summary <- function(df, value_cols, rule_col = "rule", cohort_col = "cohort") {
  by <- aggregate(df[value_cols], by = df[c(cohort_col, rule_col)], FUN = mean, na.rm = TRUE)
  agg <- aggregate(by[value_cols], by = by[rule_col], FUN = mean, na.rm = TRUE)
  agg
}

## 队列层级配对差 + 探索性 t 区间（推断单位 = 队列）
cohort_paired <- function(df, rule_col = "rule", cohort_col = "cohort",
                          val_col = "dCt", n_min = 4) {
  lv <- aggregate(df[val_col], by = df[c(cohort_col, rule_col)], FUN = mean, na.rm = TRUE)
  out <- do.call(rbind, lapply(split(lv, lv[[rule_col]]), function(d) {
    x <- d[[val_col]]; x <- x[is.finite(x)]
    tt <- if (length(x) >= n_min) tryCatch(t.test(x), error=function(e) NULL) else NULL
    data.frame(rule = d[[rule_col]][1], n_cohort = length(x), mean = mean(x),
               lo = if (is.null(tt)) NA_real_ else tt$conf.int[1],
               hi = if (is.null(tt)) NA_real_ else tt$conf.int[2],
               p  = if (is.null(tt)) NA_real_ else tt$p.value)
  }))
  out[order(-out$mean), ]
}
