#!/usr/bin/env Rscript
# ============================================================
# 53 —— 同癌种、同签名对照（修订期探索性补充）
#
# 三组任务（任务 = 特征集 × 划分）：
#   A  TCGA-BRCA × 99 个 MSigDB 基因集 × 5 划分（读 19_BRCA 缓存）
#   B  TCGA-BRCA × 6 个已发表签名 × 5 划分（本脚本重训；若已保存预测则直接复用）
#   C  METABRIC  × 6 个已发表签名 × 5 划分（读 31 缓存）
#
# 2026-10-20 第二轮修正：
#   * B 组并列打破改用与 A/C 相同的 make_prio()（对候选名排序后分配随机数）；
#     并记录旧规则（按 mods 顺序 runif）下的赢家/top-five 变化与数值差
#   * 复用 results/main/53_brca_signature_predictions.rds，不因规则修正而重训
#   * lasso_constant 统一按**验证**预测定义（与主文一致）
#   * 基因保留量按 seed 检查一致性后再汇总
# ============================================================
suppressMessages({ library(survival); library(glmnet); library(randomForestSRC); library(gbm) })
options(stringsAsFactors = FALSE)
`%||%` <- function(a, b) if (is.null(a)) b else a

PROJ <- Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
TAB  <- file.path(PROJ, "results/tables")
con  <- file(file.path(PROJ, "logs/53_same_signature.log"), "wt")
say  <- function(...) { m <- paste0(...); cat(m, "\n", sep = ""); cat(m, "\n", sep = "", file = con); flush(con) }

source(file.path(PROJ, "scripts/_lib_analysis.R"))
source(file.path(PROJ, "scripts/_lib_signature.R"))

SEEDS <- c(20261014, 20261016, 20261017, 20261018, 20261019)
MAX_MISS <- 0.20
PRED_F <- file.path(PROJ, "results/main/53_brca_signature_predictions.rds")

## ---------- B：预测（复用或重训） ----------
if (file.exists(PRED_F) && Sys.getenv("SIG53_REFIT", "0") != "1") {
  PRED <- readRDS(PRED_F)
  say("复用已保存的 BRCA 签名任务预测：", length(PRED), " 个任务")
} else {
  say("===== 组 B：TCGA-BRCA 六签名（重训）=====")
  expr <- as.matrix(read.table(file.path(PROJ, "data/cohorts/BRCA_expr.tsv"),
                               header = TRUE, sep = "\t", row.names = 1, check.names = FALSE))
  if (nrow(expr) < ncol(expr)) expr <- t(expr)
  expr <- log2(expr + 1)
  sv <- read.table(file.path(PROJ, "data/cohorts/BRCA_surv.tsv"), header = TRUE, sep = "\t",
                   row.names = 1, check.names = FALSE)
  raw <- as.character(sv[["OS_STATUS"]])
  sv$t <- suppressWarnings(as.numeric(sv[["OS_MONTHS"]]))
  sv$e <- ifelse(raw %in% c("1:DECEASED"), 1L, ifelse(raw %in% c("0:LIVING"), 0L, NA_integer_))
  sv <- sv[!is.na(sv$t) & !is.na(sv$e) & sv$t > 0, ]
  colnames(expr) <- substr(colnames(expr), 1, 12)
  expr <- expr[, !duplicated(colnames(expr)), drop = FALSE]
  com <- intersect(colnames(expr), rownames(sv)); expr <- expr[, com, drop = FALSE]; sv <- sv[com, ]
  LF <- file.path(PROJ, "results/manifest/public_signatures_genes_mapped.tsv")
  L <- read.table(LF, sep = "\t", header = TRUE, quote = "", check.names = FALSE)
  SIGS <- lapply(split(L, L$signature), function(d) unique(d$gene_used))

  PRED <- list()
  for (sd in SEEDS) {
    set.seed(sd)
    ev <- which(sv$e == 1); ce <- which(sv$e == 0)
    tr <- sort(c(sample(ev, floor(.5 * length(ev))), sample(ce, floor(.5 * length(ce)))))
    rest <- setdiff(seq_along(sv$e), tr)
    re <- rest[sv$e[rest] == 1]; rc <- rest[sv$e[rest] == 0]
    va <- sort(c(sample(re, floor(.5 * length(re))), sample(rc, floor(.5 * length(rc)))))
    te <- setdiff(rest, va)
    for (sg in names(SIGS)) {
      g <- intersect(SIGS[[sg]], rownames(expr))
      if (length(g) < 8) next
      sub <- expr[g, , drop = FALSE]
      trm <- sub[, tr, drop = FALSE]; miss <- rowSums(is.na(trm)) / ncol(trm)
      ex <- sub[miss <= MAX_MISS, , drop = FALSE]
      if (nrow(ex) < 8) next
      med <- apply(ex[, tr, drop = FALSE], 1, function(v) median(v, na.rm = TRUE))
      for (j in seq_len(nrow(ex))) { v <- ex[j, ]; k <- is.na(v); if (any(k)) { v[k] <- med[j]; ex[j, ] <- v } }
      mu <- rowMeans(ex[, tr, drop = FALSE]); sdd <- apply(ex[, tr, drop = FALSE], 1, sd)
      kp <- is.finite(sdd) & sdd > 0; ex <- ex[kp, , drop = FALSE]; mu <- mu[kp]; sdd <- sdd[kp]
      Z <- (ex - mu) / sdd; Z[!is.finite(Z)] <- 0
      G <- function(idx) t(Z[, colnames(Z)[idx], drop = FALSE])
      pm <- tryCatch(fit_pool(G(tr), sv$t[tr], sv$e[tr], G(va), G(te)), error = function(e) NULL)
      if (is.null(pm) || length(pm) < 4) next
      mods <- POOL_NAMES[POOL_NAMES %in% names(pm)]
      if (!("Lasso" %in% mods)) next
      mk <- function(key) { M <- sapply(mods, function(x) pm[[x]][[key]]); if (is.null(dim(M))) M <- matrix(M, ncol = 1)
        colnames(M) <- mods; M }
      PRED[[paste(sg, sd)]] <- list(signature = sg, seed = sd,
        train_ids = colnames(Z)[tr], val_ids = colnames(Z)[va], test_ids = colnames(Z)[te],
        t_train = sv$t[tr], e_train = sv$e[tr], t_val = sv$t[va], e_val = sv$e[va],
        t_test = sv$t[te], e_test = sv$e[te],
        val_pred = mk("pb"), test_pred = mk("pc"), genes = rownames(Z))
    }
  }
  saveRDS(PRED, PRED_F)
  say("已保存预测：", PRED_F)
}

## ---------- B：指标（两种并列规则对照） ----------
say(""); say("===== 组 B 指标 =====")
BT <- list()
for (p in PRED) {
  mods <- colnames(p$val_pred)
  V <- vapply(mods, function(m) ci_direct(p$t_val, p$e_val, p$val_pred[, m]), numeric(1)); names(V) <- mods
  Tt <- vapply(mods, function(m) ci_direct(p$t_test, p$e_test, p$test_pred[, m]), numeric(1)); names(Tt) <- mods
  if (!all(is.finite(V)) || !all(is.finite(Tt))) next

  prio <- make_prio(mods, p$seed)                    # 与 A/C 组同一规则
  wn <- mods[pick_max(V, prio[mods])]
  t5n <- mods[order_rule(V, prio[mods])[seq_len(min(5, length(mods)))]]
  ens_n <- ens_pred(p$val_pred, p$test_pred, t5n)
  ct5n <- ci_direct(p$t_test, p$e_test, ens_n$t)

  r5 <- if ("Ridge" %in% mods) Tt[["Ridge"]] else NA_real_
  e5 <- if ("Enet0.5" %in% mods) Tt[["Enet0.5"]] else NA_real_
  BT[[length(BT) + 1]] <- data.frame(
    group = "TCGA-BRCA / six signatures", signature = p$signature, seed = p$seed,
    n_gene = length(p$genes), n_gene_final = length(p$genes), n_candidate = length(mods),
    n_test = length(p$t_test), ev_test = sum(p$e_test),
    lasso_constant_val = length(unique(p$val_pred[, "Lasso"])) == 1,
    Ct_lasso = Tt[["Lasso"]], Ct_ridge = r5, Ct_enet05 = e5,
    Ct_argmax = Tt[[wn]], Ct_top5 = ct5n,
    Gv = V[[wn]] - V[["Lasso"]], Gt = Tt[[wn]] - Tt[["Lasso"]],
    Orel = (V[[wn]] - V[["Lasso"]]) - (Tt[[wn]] - Tt[["Lasso"]]),
    d_argmax_ridge = Tt[[wn]] - r5, d_top5_ridge = ct5n - r5,
    d_top5_argmax = ct5n - Tt[[wn]], winner = wn,
    genes = paste(p$genes, collapse = "|"))
}
B <- do.call(rbind, BT)
write.table(B, file.path(TAB, "53_brca_signature_tasks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
say("  任务：", nrow(B), "（并列差异见 55c，本脚本不含旧规则）")

## 划分一致性：B 组与 A 组同 seed 的患者集合
f19 <- file.path(PROJ, "results/main", sprintf("19_BRCA_%d_100.rds", SEEDS[1]))
chk_ids <- list()
for (sd in SEEDS) {
  f <- file.path(PROJ, "results/main", sprintf("19_BRCA_%d_100.rds", sd))
  if (!file.exists(f)) next
  R <- readRDS(f)
  b <- PRED[[grep(sprintf(" %d$", sd), names(PRED))[1]]]
  if (is.null(b)) next
  chk_ids[[length(chk_ids) + 1]] <- data.frame(seed = sd,
    train_same = identical(sort(b$train_ids), sort(R$train_ids)),
    val_same = identical(sort(b$val_ids), sort(R$val_ids)),
    test_same = identical(sort(b$test_ids), sort(R$test_ids)))
}
IDS <- do.call(rbind, chk_ids)
write.table(IDS, file.path(TAB, "53_split_identity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
say("  与主分析划分完全一致：", all(IDS$train_same & IDS$val_same & IDS$test_same))

## ---------- A：TCGA-BRCA × MSigDB ----------
say(""); say("===== 组 A：TCGA-BRCA MSigDB（读缓存）=====")
AT <- list()
for (f in list.files(file.path(PROJ, "results/main"), pattern = "^19_BRCA_[0-9]{8}_100\\.rds$", full.names = TRUE)) {
  R <- readRDS(f); G <- R$genesets
  cnt <- table(unlist(lapply(G, function(z) colnames(z$val_pred))))
  keep <- names(cnt)[cnt / length(G) >= 0.90]
  prio <- make_prio(keep, R$seed)
  for (z in G) {
    if (!all(keep %in% colnames(z$val_pred))) next
    Pv <- z$val_pred[, keep, drop = FALSE]; Pt <- z$test_pred[, keep, drop = FALSE]
    if (!all(is.finite(Pv)) || !all(is.finite(Pt))) next
    V <- vapply(keep, function(m) ci_direct(z$tva, z$eva, Pv[, m]), numeric(1)); names(V) <- keep
    Tt <- vapply(keep, function(m) ci_direct(z$tte, z$ete, Pt[, m]), numeric(1)); names(Tt) <- keep
    if (!all(is.finite(V)) || !all(is.finite(Tt))) next
    w <- keep[pick_max(V, prio[keep])]
    top5 <- keep[order_rule(V, prio[keep])[seq_len(min(5, length(keep)))]]
    ens <- ens_pred(Pv, Pt, top5)
    ct5 <- ci_direct(z$tte, z$ete, ens$t)
    r5 <- Tt[["Ridge"]]; e5 <- if ("Enet0.5" %in% keep) Tt[["Enet0.5"]] else NA_real_
    AT[[length(AT) + 1]] <- data.frame(
      group = "TCGA-BRCA / MSigDB gene sets", signature = z$name, seed = R$seed,
      n_gene = NA_integer_, n_gene_final = NA_integer_, n_candidate = length(keep),
      n_test = length(z$tte), ev_test = sum(z$ete),
      lasso_constant_val = length(unique(Pv[, "Lasso"])) == 1,
      Ct_lasso = Tt[["Lasso"]], Ct_ridge = r5, Ct_enet05 = e5,
      Ct_argmax = Tt[[w]], Ct_top5 = ct5,
      Gv = V[[w]] - V[["Lasso"]], Gt = Tt[[w]] - Tt[["Lasso"]],
      Orel = (V[[w]] - V[["Lasso"]]) - (Tt[[w]] - Tt[["Lasso"]]),
      d_argmax_ridge = Tt[[w]] - r5, d_top5_ridge = ct5 - r5,
      d_top5_argmax = ct5 - Tt[[w]], winner = w, genes = NA_character_)
  }
}
A <- do.call(rbind, AT)
say("  组 A 任务：", nrow(A))

## ---------- C：METABRIC ----------
say(""); say("===== 组 C：METABRIC 六签名（读 31 缓存）=====")
mb <- readRDS(file.path(PROJ, "results/metabric/31_metabric_predictions.rds"))
CT <- list(); GENEC <- list()
for (z in mb) {
  mods <- z$model
  V <- vapply(mods, function(m) ci_direct(z$t_val, z$e_val, z$val_pred[, m]), numeric(1)); names(V) <- mods
  Tt <- vapply(mods, function(m) ci_direct(z$t_test, z$e_test, z$test_pred[, m]), numeric(1)); names(Tt) <- mods
  if (!all(is.finite(V)) || !all(is.finite(Tt))) next
  prio <- make_prio(mods, z$seed)
  w <- mods[pick_max(V, prio[mods])]
  top5 <- mods[order_rule(V, prio[mods])[seq_len(min(5, length(mods)))]]
  ens <- ens_pred(z$val_pred, z$test_pred, top5)
  ct5 <- ci_direct(z$t_test, z$e_test, ens$t)
  CT[[length(CT) + 1]] <- data.frame(
    group = "METABRIC / six signatures", signature = z$signature, seed = z$seed,
    n_gene = length(z$genes), n_gene_final = length(z$genes), n_candidate = length(mods),
    n_test = length(z$t_test), ev_test = sum(z$e_test),
    lasso_constant_val = length(unique(z$val_pred[, "Lasso"])) == 1,
    Ct_lasso = Tt[["Lasso"]], Ct_ridge = Tt[["Ridge"]],
    Ct_enet05 = if ("Enet0.5" %in% mods) Tt[["Enet0.5"]] else NA_real_,
    Ct_argmax = Tt[[w]], Ct_top5 = ct5,
    Gv = V[[w]] - V[["Lasso"]], Gt = Tt[[w]] - Tt[["Lasso"]],
    Orel = (V[[w]] - V[["Lasso"]]) - (Tt[[w]] - Tt[["Lasso"]]),
    d_argmax_ridge = Tt[[w]] - Tt[["Ridge"]], d_top5_ridge = ct5 - Tt[["Ridge"]],
    d_top5_argmax = ct5 - Tt[[w]], winner = w,
    genes = paste(z$genes, collapse = "|"))
  GENEC[[z$signature]] <- z$genes
}
C <- do.call(rbind, CT)
say("  组 C 任务：", nrow(C))

ALL <- rbind(A, B, C)
write.table(ALL, file.path(TAB, "53_same_signature_tasks_all.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---------- 汇总 ----------
mc <- c("Ct_lasso", "Ct_ridge", "Ct_enet05", "Ct_argmax", "Ct_top5",
        "Gv", "Gt", "Orel", "d_argmax_ridge", "d_top5_ridge", "d_top5_argmax")
grp <- do.call(rbind, lapply(split(ALL, ALL$group), function(d) {
  out <- data.frame(group = d$group[1], n_task = nrow(d), n_test = mean(d$n_test),
                    ev_test = mean(d$ev_test),
                    lasso_constant_pct = 100 * mean(d$lasso_constant_val, na.rm = TRUE))
  for (k in mc) out[[k]] <- mean(d[[k]], na.rm = TRUE)
  out
}))
write.table(grp, file.path(TAB, "53_same_signature_groups.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

sub <- ALL[ALL$group != "TCGA-BRCA / MSigDB gene sets", ]
sigsum <- do.call(rbind, lapply(split(sub, list(sub$group, sub$signature), drop = TRUE), function(d) {
  data.frame(group = d$group[1], signature = d$signature[1], n_task = nrow(d),
             n_gene_final = mean(d$n_gene_final, na.rm = TRUE),
             Ct_lasso = mean(d$Ct_lasso, na.rm = TRUE), Ct_ridge = mean(d$Ct_ridge, na.rm = TRUE),
             Ct_argmax = mean(d$Ct_argmax, na.rm = TRUE), Ct_top5 = mean(d$Ct_top5, na.rm = TRUE),
             d_argmax_ridge = mean(d$d_argmax_ridge, na.rm = TRUE),
             d_top5_ridge = mean(d$d_top5_ridge, na.rm = TRUE),
             d_top5_argmax = mean(d$d_top5_argmax, na.rm = TRUE),
             Gv = mean(d$Gv, na.rm = TRUE), Gt = mean(d$Gt, na.rm = TRUE), Orel = mean(d$Orel, na.rm = TRUE))
}))
sigsum <- sigsum[order(sigsum$group, sigsum$signature), ]
write.table(sigsum, file.path(TAB, "53_same_signature_by_signature.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## 基因保留：先检查各 seed 是否一致
ov <- do.call(rbind, lapply(unique(B$signature), function(sg) {
  gl <- Filter(Negate(is.null), lapply(SEEDS, function(sd) {
    r <- B$genes[B$signature == sg & B$seed == sd]
    if (!length(r)) return(NULL)
    sort(unique(strsplit(r[1], "|", fixed = TRUE)[[1]]))
  }))
  gb <- gl[[1]]
  same_b <- all(vapply(gl, identical, logical(1), gb))
  cl <- lapply(unique(C$seed[C$signature == sg]), function(sd) {
    r <- C$genes[C$signature == sg & C$seed == sd]
    sort(unique(strsplit(r[1], "|", fixed = TRUE)[[1]]))
  })
  gc <- cl[[1]]
  same_c <- all(vapply(cl, identical, logical(1), gc))
  data.frame(signature = sg, n_tcga_brca = length(gb), n_metabric = length(gc),
             n_shared = length(intersect(gb, gc)),
             genes_same_across_seeds = same_b && same_c)
}))
write.table(ov, file.path(TAB, "53_signature_gene_overlap.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

say("")
say("================ 三组对照（Harrell C，任务层均值）================")
for (i in seq_len(nrow(grp))) {
  say(sprintf("  %-34s %5d  Ct_L=%.4f  argmax-R=%+.4f  top5-R=%+.4f  Gt=%+.4f",
      grp$group[i], grp$n_task[i], grp$Ct_lasso[i], grp$d_argmax_ridge[i], grp$d_top5_ridge[i], grp$Gt[i]))
}
say("")
say("  并列规则差异见 55c（本脚本只使用 make_prio）")
say("  基因集跨 seed 是否一致：", all(ov$genes_same_across_seeds))
close(con)
