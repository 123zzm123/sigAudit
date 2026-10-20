#!/usr/bin/env Rscript
# ============================================================
# 53 —— 同癌种、同签名对照（修订期探索性补充）
#
# 问题：TCGA 里“固定 ridge 略优于验证集选择”的排名，在换成更接近的任务条件后
#       是否仍然存在？本模块只回答“观察到的差异是否在更可比的任务下仍出现”，
#       不声称证明生物学或平台原因。
#
# 三组任务（单位：任务 = 特征集 × 划分）：
#   A  TCGA-BRCA × 99 个 MSigDB 基因集 × 5 划分（主分析原样本，读 19_BRCA 缓存）
#   B  TCGA-BRCA × 6 个已发表签名 × 5 划分（本脚本重训，候选池与 26b 完全一致）
#   C  METABRIC  × 6 个已发表签名 × 5 划分（读 31 缓存）
#
# 每个任务报告：固定参照的测试 Ct、argmax−ridge、top5−ridge、top5−argmax，
# 以及相对固定 Lasso 的 Gv / Gt / Orel。签名不是独立外部队列，只作描述。
#
# 冻结：候选池、划分种子、标准化与缺失处理与主分析一致；不因结果方向更换签名、
#       不做亚型筛选、不把六签名当六个外部验证集。
# ============================================================
suppressMessages({ library(survival); library(glmnet); library(randomForestSRC); library(gbm) })
options(stringsAsFactors = FALSE)
`%||%` <- function(a, b) if (is.null(a)) b else a

PROJ <- Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
TAB  <- file.path(PROJ, "results/tables")
con  <- file(file.path(PROJ, "logs/53_same_signature.log"), "wt")
say  <- function(...) { m <- paste0(...); cat(m, "\n", sep = ""); cat(m, "\n", sep = "", file = con); flush(con) }

source(file.path(PROJ, "scripts/_lib_analysis.R"))    # ci_direct / pick_max / order_rule / ens_pred
source(file.path(PROJ, "scripts/_lib_signature.R"))   # POOL / POOL_NAMES / fit_pool（与 26b 同一实现）

SEEDS <- c(20261014, 20261016, 20261017, 20261018, 20261019)
MAX_MISS <- 0.20

## ---------- B：TCGA-BRCA × 六签名（重训 30 个任务） ----------
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
say("  BRCA 表达矩阵：", nrow(expr), " 基因 × ", ncol(expr), " 患者；事件 ", sum(sv$e))

LF <- file.path(PROJ, "results/manifest/public_signatures_genes_mapped.tsv")
if (!file.exists(LF)) stop("缺少签名清单")
L <- read.table(LF, sep = "\t", header = TRUE, quote = "", check.names = FALSE)
SIGS <- lapply(split(L, L$signature), function(d) unique(d$gene_used))

BT <- list(); PRED <- list()
for (sd in SEEDS) {
  set.seed(sd)
  ev <- which(sv$e == 1); ce <- which(sv$e == 0)
  tr <- sort(c(sample(ev, floor(.5 * length(ev))), sample(ce, floor(.5 * length(ce)))))
  rest <- setdiff(seq_along(sv$e), tr)
  re <- rest[sv$e[rest] == 1]; rc <- rest[sv$e[rest] == 0]
  va <- sort(c(sample(re, floor(.5 * length(re))), sample(rc, floor(.5 * length(rc)))))
  te <- setdiff(rest, va)
  say("  seed ", sd, ": train ", length(tr), " / val ", length(va), " / test ", length(te),
      "；测试事件 ", sum(sv$e[te]))
  for (sg in names(SIGS)) {
    g <- intersect(SIGS[[sg]], rownames(expr))
    if (length(g) < 8) { say("    跳过 ", sg, "：交集 ", length(g)); next }
    sub <- expr[g, , drop = FALSE]
    trm <- sub[, tr, drop = FALSE]; miss <- rowSums(is.na(trm)) / ncol(trm)
    ex <- sub[miss <= MAX_MISS, , drop = FALSE]
    if (nrow(ex) < 8) { say("    跳过 ", sg, "：缺失过滤后 <8"); next }
    med <- apply(ex[, tr, drop = FALSE], 1, function(v) median(v, na.rm = TRUE))
    for (j in seq_len(nrow(ex))) { v <- ex[j, ]; k <- is.na(v); if (any(k)) { v[k] <- med[j]; ex[j, ] <- v } }
    mu <- rowMeans(ex[, tr, drop = FALSE]); sdd <- apply(ex[, tr, drop = FALSE], 1, sd)
    kp <- is.finite(sdd) & sdd > 0; ex <- ex[kp, , drop = FALSE]; mu <- mu[kp]; sdd <- sdd[kp]
    Z <- (ex - mu) / sdd; Z[!is.finite(Z)] <- 0
    G <- function(idx) t(Z[, colnames(Z)[idx], drop = FALSE])
    pm <- tryCatch(fit_pool(G(tr), sv$t[tr], sv$e[tr], G(va), G(te)), error = function(e) NULL)
    if (is.null(pm) || length(pm) < 4) { say("    跳过 ", sg, "：候选 <4"); next }
    mods <- POOL_NAMES[POOL_NAMES %in% names(pm)]
    if (!("Lasso" %in% mods)) next
    mk <- function(key) { M <- sapply(mods, function(x) pm[[x]][[key]]); if (is.null(dim(M))) M <- matrix(M, ncol = 1)
      colnames(M) <- mods; M }
    Pv <- mk("pb"); Pt <- mk("pc")
    V <- vapply(mods, function(m) ci_direct(sv$t[va], sv$e[va], Pv[, m]), numeric(1)); names(V) <- mods
    Tt <- vapply(mods, function(m) ci_direct(sv$t[te], sv$e[te], Pt[, m]), numeric(1)); names(Tt) <- mods
    if (!all(is.finite(V)) || !all(is.finite(Tt))) { say("    跳过 ", sg, "：C-index 非法"); next }
    set.seed(sd); prio <- runif(length(mods)); names(prio) <- mods     # 与 26b 完全一致的并列打破规则
    w <- mods[pick_max(V, prio[mods])]
    top5 <- mods[order_rule(V, prio[mods])[seq_len(min(5, length(mods)))]]
    ens <- ens_pred(Pv, Pt, top5)
    ct5 <- ci_direct(sv$t[te], sv$e[te], ens$t)
    r5 <- if ("Ridge" %in% mods) Tt[["Ridge"]] else NA_real_
    e5 <- if ("Enet0.5" %in% mods) Tt[["Enet0.5"]] else NA_real_
    BT[[length(BT) + 1]] <- data.frame(
      group = "TCGA-BRCA / six signatures", signature = sg, seed = sd,
      n_gene = length(g), n_gene_final = nrow(Z), n_candidate = length(mods),
      n_test = length(te), ev_test = sum(sv$e[te]),
      lasso_constant = length(unique(Pt[, "Lasso"])) == 1,
      Ct_lasso = Tt[["Lasso"]], Ct_ridge = r5, Ct_enet05 = e5,
      Ct_argmax = Tt[[w]], Ct_top5 = ct5,
      Gv = V[[w]] - V[["Lasso"]], Gt = Tt[[w]] - Tt[["Lasso"]],
      Orel = (V[[w]] - V[["Lasso"]]) - (Tt[[w]] - Tt[["Lasso"]]),
      d_argmax_ridge = Tt[[w]] - r5, d_top5_ridge = ct5 - r5,
      d_top5_argmax = ct5 - Tt[[w]], winner = w,
      genes = paste(rownames(Z), collapse = "|"))
    PRED[[paste(sg, sd)]] <- list(signature = sg, seed = sd, train_ids = colnames(Z)[tr],
      val_ids = colnames(Z)[va], test_ids = colnames(Z)[te],
      t_train = sv$t[tr], e_train = sv$e[tr], t_val = sv$t[va], e_val = sv$e[va],
      t_test = sv$t[te], e_test = sv$e[te],
      val_pred = Pv, test_pred = Pt, genes = rownames(Z))
  }
}
B <- do.call(rbind, BT)
saveRDS(PRED, file.path(PROJ, "results/main/53_brca_signature_predictions.rds"))
write.table(B, file.path(TAB, "53_brca_signature_tasks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
say("  组 B 可评估任务：", nrow(B))

## ---------- A：TCGA-BRCA × MSigDB（读缓存） ----------
say(""); say("===== 组 A：TCGA-BRCA MSigDB（读 19_BRCA 缓存）=====")
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
      lasso_constant = length(unique(Pt[, "Lasso"])) == 1,
      Ct_lasso = Tt[["Lasso"]], Ct_ridge = r5, Ct_enet05 = e5,
      Ct_argmax = Tt[[w]], Ct_top5 = ct5,
      Gv = V[[w]] - V[["Lasso"]], Gt = Tt[[w]] - Tt[["Lasso"]],
      Orel = (V[[w]] - V[["Lasso"]]) - (Tt[[w]] - Tt[["Lasso"]]),
      d_argmax_ridge = Tt[[w]] - r5, d_top5_ridge = ct5 - r5,
      d_top5_argmax = ct5 - Tt[[w]], winner = w, genes = NA_character_)
  }
}
A <- do.call(rbind, AT)
say("  组 A 可评估任务：", nrow(A))

## ---------- C：METABRIC × 六签名（读缓存） ----------
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
    lasso_constant = length(unique(z$test_pred[, "Lasso"])) == 1,
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
say("  组 C 可评估任务：", nrow(C))

ALL <- rbind(A, B, C)
write.table(ALL, file.path(TAB, "53_same_signature_tasks_all.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---------- 汇总 ----------
mc <- c("Ct_lasso", "Ct_ridge", "Ct_enet05", "Ct_argmax", "Ct_top5",
        "Gv", "Gt", "Orel", "d_argmax_ridge", "d_top5_ridge", "d_top5_argmax")
grp <- do.call(rbind, lapply(split(ALL, ALL$group), function(d) {
  out <- data.frame(group = d$group[1], n_task = nrow(d), n_test = mean(d$n_test),
                    ev_test = mean(d$ev_test),
                    lasso_constant_pct = 100 * mean(d$lasso_constant, na.rm = TRUE))
  for (k in mc) out[[k]] <- mean(d[[k]], na.rm = TRUE)
  out
}))
write.table(grp, file.path(TAB, "53_same_signature_groups.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

sigsum <- do.call(rbind, lapply(split(ALL[ALL$group != "TCGA-BRCA / MSigDB gene sets", ],
                                     list(ALL$group[ALL$group != "TCGA-BRCA / MSigDB gene sets"],
                                          ALL$signature[ALL$group != "TCGA-BRCA / MSigDB gene sets"]),
                                     drop = TRUE), function(d) {
  out <- data.frame(group = d$group[1], signature = d$signature[1], n_task = nrow(d),
                    n_gene_final = mean(d$n_gene_final, na.rm = TRUE),
                    Ct_lasso = mean(d$Ct_lasso, na.rm = TRUE), Ct_ridge = mean(d$Ct_ridge, na.rm = TRUE),
                    Ct_argmax = mean(d$Ct_argmax, na.rm = TRUE), Ct_top5 = mean(d$Ct_top5, na.rm = TRUE),
                    d_argmax_ridge = mean(d$d_argmax_ridge, na.rm = TRUE),
                    d_top5_ridge = mean(d$d_top5_ridge, na.rm = TRUE),
                    d_top5_argmax = mean(d$d_top5_argmax, na.rm = TRUE),
                    Gv = mean(d$Gv, na.rm = TRUE), Gt = mean(d$Gt, na.rm = TRUE), Orel = mean(d$Orel, na.rm = TRUE))
  out
}))
sigsum <- sigsum[order(sigsum$group, sigsum$signature), ]
write.table(sigsum, file.path(TAB, "53_same_signature_by_signature.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## 基因保留与交集
ov <- do.call(rbind, lapply(names(SIGS), function(sg) {
  gb <- unique(unlist(strsplit(B$genes[B$signature == sg][1], "|", fixed = TRUE)))
  gc <- GENEC[[sg]]
  data.frame(signature = sg, n_source = length(unique(SIGS[[sg]])),
             n_tcga_brca = length(gb), n_metabric = length(gc),
             n_shared = length(intersect(gb, gc)))
}))
write.table(ov, file.path(TAB, "53_signature_gene_overlap.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

say("")
say("================ 三组对照（Harrell C，任务层均值）================")
say(sprintf("  %-34s %6s %8s %8s %8s %8s %11s %11s", "组", "任务", "Ct_Lasso", "Ct_Ridge",
            "argmax−R", "top5−R", "Gv", "Gt"))
for (i in seq_len(nrow(grp))) {
  say(sprintf("  %-34s %6d %8.4f %8.4f %+8.4f %+8.4f %+11.4f %+11.4f",
      grp$group[i], grp$n_task[i], grp$Ct_lasso[i], grp$Ct_ridge[i],
      grp$d_argmax_ridge[i], grp$d_top5_ridge[i], grp$Gv[i], grp$Gt[i]))
}
say("")
say("  签名层（B/C）描述：")
for (i in seq_len(nrow(sigsum))) {
  say(sprintf("  %-12s %-20s argmax−R %+.4f  top5−R %+.4f  Gt %+.4f", sigsum$group[i],
      sigsum$signature[i], sigsum$d_argmax_ridge[i], sigsum$d_top5_ridge[i], sigsum$Gt[i]))
}
say("")
say("  说明：签名不是独立外部队列；六签名只作描述，不做六个外部验证的推断。")
close(con)
