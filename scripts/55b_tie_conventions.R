#!/usr/bin/env Rscript
# ============================================================
# 55b —— 时间并列与分数并列的处理约定（小样本固定例子）
#
# 目的：把三种“并列”分开写清楚，避免把 Harrell/Uno 的差异概括成
# “所有时间并列都计 0.5”：
#   (a) 两个事件同一时刻
#   (b) 一个事件与一个删失同一时刻
#   (c) 两个风险分数相同
#
# 三种约定：
#   strict     本分析：T_i < T_j 才可比；事件并列不比对；分数并列 0.5
#   survival   survival::concordance(timewt = "n")：事件–删失同刻视为事件在前；
#              事件时间并列记为 tied outcome，不进入分母；分数并列 0.5
#   textbook   把事件时间并列按 0.5 计入（本分析不采用，仅作对照）
#
# 输出：results/tables/55_tie_conventions.tsv
# ============================================================
suppressMessages(library(survival))
options(stringsAsFactors = FALSE)

PROJ <- Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
TAB  <- file.path(PROJ, "results/tables")

## 严格约定（与 _lib_uno.R 相同）
strict_c <- function(t, e, s) {
  num <- den <- 0
  for (i in seq_along(t)) {
    if (e[i] != 1) next
    for (j in seq_along(t)) {
      if (t[j] > t[i]) {
        den <- den + 1
        num <- num + (s[i] > s[j]) + 0.5 * (s[i] == s[j])
      }
    }
  }
  if (den == 0) return(c(NA, 0))
  c(num / den, den)
}

## 教科书变体：事件时间并列按 0.5 计入
textbook_c <- function(t, e, s) {
  num <- den <- 0
  for (i in seq_along(t)) {
    if (e[i] != 1) next
    for (j in seq_along(t)) {
      if (t[j] > t[i]) {
        den <- den + 1
        num <- num + (s[i] > s[j]) + 0.5 * (s[i] == s[j])
      } else if (t[j] == t[i] && e[j] == 1 && j > i) {
        den <- den + 1
        num <- num + 0.5
      }
    }
  }
  if (den == 0) return(c(NA, 0))
  c(num / den, den)
}

surv_c <- function(t, e, s) {
  cc <- tryCatch(concordance(Surv(t, e) ~ s, reverse = TRUE, timewt = "n"),
                 error = function(x) NULL)
  if (is.null(cc)) return(c(NA, NA))
  c(as.numeric(cc$concordance), sum(cc$count[c("concordant", "discordant", "tied.x")]))
}

CASES <- list(
  list(name = "event_vs_censor_same_time", t = c(10, 10), e = c(1, 0), s = c(1, 0)),
  list(name = "two_events_same_time",      t = c(10, 10), e = c(1, 1), s = c(1, 0)),
  list(name = "mixed_seven_observations",  t = c(10, 10, 10, 20, 20, 30, 40),
       e = c(1, 1, 0, 1, 0, 1, 0), s = c(5, 4, 3, 2, 1, 0, -1))
)

rows <- list()
for (cs in CASES) {
  a <- strict_c(cs$t, cs$e, cs$s)
  b <- surv_c(cs$t, cs$e, cs$s)
  d <- textbook_c(cs$t, cs$e, cs$s)
  fmt <- function(x) ifelse(is.na(x), "NA", sprintf("%.6f", x))
  rows[[length(rows) + 1]] <- data.frame(
    dataset = cs$name, n = length(cs$t),
    convention = "strict (this analysis)",
    C = fmt(a[1]), n_pairs = a[2],
    note = "T_i < T_j only; tied event times and event-censor ties excluded; tied scores 0.5")
  rows[[length(rows) + 1]] <- data.frame(
    dataset = cs$name, n = length(cs$t),
    convention = "survival::concordance (timewt = n)",
    C = fmt(b[1]), n_pairs = b[2],
    note = "event-censor ties kept (censored considered later); tied outcomes dropped; tied scores 0.5")
  rows[[length(rows) + 1]] <- data.frame(
    dataset = cs$name, n = length(cs$t),
    convention = "textbook variant (not used)",
    C = fmt(d[1]), n_pairs = d[2],
    note = "tied event times credited 0.5; shown only to separate the conventions")
}
OUT <- do.call(rbind, rows)
write.table(OUT, file.path(TAB, "55_tie_conventions.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
print(OUT, row.names = FALSE)
