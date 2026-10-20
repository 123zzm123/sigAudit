# ============================================================
# _lib_uno.R —— Uno/IPCW 一致性分析的共用函数
#
# 定义（冻结，脚本 52 与 55 共用同一份实现）：
#   * 删失分布 G 由训练折估计，Kaplan–Meier 作用于删失指示 (1 - e)，
#     取右连续阶梯：G(t) = S(t_k), k = max{k: t_k <= t}；t 小于首个时间点取 1
#   * 可比对：T_i < T_j（严格），delta_i = 1，T_i < tau
#   * 权重 1/G(T_i)^2
#   * 风险分数并列按 0.5 计；时间并列不进入可比对
#   * tau 为固定临床窗（60 月），不按队列/任务调整
# ============================================================
suppressMessages(library(survival))

gfun_train <- function(t_tr, e_tr) {
  sf <- survfit(Surv(t_tr, 1 - e_tr) ~ 1)
  tm <- sf$time; sv <- sf$surv
  function(x) { i <- findInterval(x, tm); ifelse(i == 0, 1, sv[pmax(i, 1)]) }
}

## 预先算好“事件 i（T_i<tau）与更晚的 j”的可比对集，供同一任务的所有模型复用
pair_prep <- function(t, e, tau, G, g_min = 0.05) {
  Gi <- G(t)
  ev <- which(e == 1 & t < tau & is.finite(Gi) & Gi >= g_min)
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
