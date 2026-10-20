#!/usr/bin/env Rscript
# ============================================================
# 55 —— Uno/IPCW 的同定义实现核对（新增自写估计器必须做的检查）
#
# 产出：
#   results/tables/55_uno_cases.tsv       合成用例 + 一个真实任务的输入（长表）
#   results/tables/55_uno_r_values.tsv    R 版（_lib_uno.R）算出的 C
#   随后由 scripts/55_uno_validation.py 用独立的 Python 实现重算并逐值比对，
#   写出 results/tables/55_uno_validation.tsv。
#
# 用例直接覆盖公式中的约定：分数全等、完全正序/反序、分数并列、
# 事件时间并列、事件与删失同刻、tau 边界（恰好等于 tau 的事件不计入）。
# ============================================================
suppressMessages({ library(survival) })
options(stringsAsFactors = FALSE)
PROJ <- Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
TAB  <- file.path(PROJ, "results/tables")
source(file.path(PROJ, "scripts/_lib_uno.R"))
source(file.path(PROJ, "scripts/_lib_analysis.R"))

TAU <- 60
CASES <- list()

add_case <- function(name, tr_t, tr_e, te_t, te_e, scores, note = "") {
  CASES[[name]] <<- list(tr_t = tr_t, tr_e = tr_e, te_t = te_t, te_e = te_e,
                         scores = scores, note = note)
}

## 1. 分数全等 -> 0.5
set.seed(1)
tr_t <- rexp(200, 0.02); tr_e <- rbinom(200, 1, 0.6)
te_t <- rexp(100, 0.02); te_e <- rbinom(100, 1, 0.6)
add_case("all_equal_scores", tr_t, tr_e, te_t, te_e,
         list(s1 = rep(0, 100)), "all risk scores identical; C must be 0.5")

## 2. 完全正序：风险越高事件越早 -> 1
te_t2 <- sort(rexp(100, 0.03)); te_e2 <- rep(1L, 100)
add_case("perfect_concordance", tr_t, tr_e, te_t2, te_e2,
         list(s1 = -te_t2), "earlier event has higher score; C must be 1")

## 3. 完全反序 -> 0
add_case("perfect_discordance", tr_t, tr_e, te_t2, te_e2,
         list(s1 = te_t2), "earlier event has lower score; C must be 0")

## 4. 分数并列：一半相同
set.seed(2)
te_t4 <- rexp(120, 0.02); te_e4 <- rbinom(120, 1, 0.7)
s4 <- round(rnorm(120), 0)          # 大量并列
add_case("tied_scores", tr_t, tr_e, te_t4, te_e4, list(s1 = s4),
         "tied scores contribute 0.5")

## 5. 事件时间并列 + 事件/删失同刻
te_t5 <- c(10, 10, 10, 20, 20, 30, 40, 50, 60, 60, 70, 80)
te_e5 <- c(1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1)
tr_t5 <- c(5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105, 115, 125, 135, 145)
tr_e5 <- c(1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1)
add_case("tied_times", tr_t5, tr_e5, te_t5, te_e5,
         list(s1 = -te_t5, s2 = te_t5), "tied event times are not comparable; strict T_i < T_j")

## 6. tau 边界：一个事件恰好等于 tau（不计入），一个刚好小于
te_t6 <- c(30, 59.999, 60, 61, 90, 120, 150, 180, 210, 240, 270, 300)
te_e6 <- rep(1L, 12)
add_case("tau_boundary", tr_t, tr_e, te_t6, te_e6,
         list(s1 = -te_t6), "event at exactly tau is excluded; strictly below is included")

## 7. 一半删失、随机分数
set.seed(3)
te_t7 <- rexp(150, 0.02); te_e7 <- rbinom(150, 1, 0.5)
add_case("random_half_censored", tr_t, tr_e, te_t7, te_e7,
         list(s1 = rnorm(150), s2 = rnorm(150), s3 = rnorm(150)),
         "three score vectors on the same test set")

## 8. 真实任务：BRCA，seed 20261014，第一个基因集
R <- readRDS(file.path(PROJ, "results/main/19_BRCA_20261014_100.rds"))
G <- R$genesets
cnt <- table(unlist(lapply(G, function(z) colnames(z$val_pred))))
keep <- names(cnt)[cnt / length(G) >= 0.90]
z <- G[[1]]
Pv <- z$val_pred[, keep, drop = FALSE]; Pt <- z$test_pred[, keep, drop = FALSE]
V <- vapply(keep, function(m) ci_direct(z$tva, z$eva, Pv[, m]), numeric(1)); names(V) <- keep
prio <- make_prio(keep, R$seed)
win <- keep[pick_max(V, prio[keep])]
top5 <- keep[order_rule(V, prio[keep])[seq_len(min(5, length(keep)))]]
ens <- ens_pred(Pv, Pt, top5)
add_case("real_BRCA_task",
         z$t_tr, z$e_tr, z$tte, z$ete,
         list(lasso = as.numeric(Pt[, "Lasso"]), ridge = as.numeric(Pt[, "Ridge"]),
              enet05 = as.numeric(Pt[, "Enet0.5"]), argmax = as.numeric(Pt[, win]),
              top5 = as.numeric(ens$t)),
         paste0("BRCA seed 20261014, geneset ", z$name, "; winner ", win))

## ---------- 写出长表 ----------
rows <- list(); params <- list()
for (nm in names(CASES)) {
  c1 <- CASES[[nm]]
  for (sp in c("train", "test")) {
    tt <- if (sp == "train") c1$tr_t else c1$te_t
    ee <- if (sp == "train") c1$tr_e else c1$te_e
    df <- data.frame(dataset = nm, split = sp, id = seq_along(tt), time = tt, event = ee)
    for (mn in names(c1$scores)) {
      df[[mn]] <- if (sp == "test") c1$scores[[mn]] else NA_real_
    }
    rows[[length(rows) + 1]] <- df
  }
  params[[length(params) + 1]] <- data.frame(dataset = nm, tau = TAU, models = paste(names(c1$scores), collapse = ","),
                                             note = c1$note)
}
allcol <- unique(unlist(lapply(rows, names)))
rows <- lapply(rows, function(d) { for (k in setdiff(allcol, names(d))) d[[k]] <- NA_real_; d[allcol] })
write.table(do.call(rbind, rows), file.path(TAB, "55_uno_cases.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(do.call(rbind, params), file.path(TAB, "55_uno_params.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ---------- R 版数值 ----------
out <- list()
for (nm in names(CASES)) {
  c1 <- CASES[[nm]]
  Gf <- gfun_train(c1$tr_t, c1$tr_e)
  pp <- pair_prep(c1$te_t, c1$te_e, TAU, Gf, 0.05)
  for (mn in names(c1$scores)) {
    out[[length(out) + 1]] <- data.frame(dataset = nm, model = mn,
      n_pairs = length(pp$w), n_events = pp$n_ev,
      uno_r = uno_from_pairs(c1$scores[[mn]], pp))
  }
}
write.table(do.call(rbind, out), file.path(TAB, "55_uno_r_values.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat("cases:", length(CASES), " R values:", length(out), "\n")
