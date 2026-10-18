## ============================================================
## 关于 cox.ties：本文件所有 cv.glmnet 都显式传 cox.ties = "breslow"。
## 原因：glmnet 5.0 的默认是 breslow，而 5.1 起默认将改为 efron（包会主动警告）。
## 若不显式指定，升级 glmnet 会让全部结果静默变化、论文无法复现。
## 已验证：在 glmnet 5.0 下显式 breslow 与"不指定时的默认"逐位完全相同
## （lambda、系数、预测差均为 0），因此钉住它不改变任何已有结果。
## ============================================================
## audit_modules.R — 各审计模块
## ============================================================

#' 模块 1：乐观度校正（bootstrap）——**实验性**
#'
#' \strong{注意：本模块未参与论文的任何结论，属实验性功能。}
#' 已知局限：缺失过滤、中位数填补与零方差剔除在**整个队列**上估计，
#' 且没有在每次 bootstrap 重采样中重新估计；因此它的乐观度估计
#' 不是嵌套重采样的严格实现。
#'
#' 表观 C-index（在全部数据上拟合并评估）总是偏乐观。本模块用 bootstrap
#' 估计乐观度并给出校正值，即模型在"新样本"上的期望表现。
#'
#' @param expr 表达矩阵（基因 × 样本 或 样本 × 基因，自动识别）
#' @param surv 生存数据框，行名为 patient id
#' @param time_col,event_col 生存时间/事件列名
#' @param genes 签名基因
#' @param method 学习器："lasso" | "ridge" | "enet"
#' @param B bootstrap 次数，默认 200
#' @param seed 随机种子，默认 20261014（**表观模型与 bootstrap 的 CV 折都用它**，
#'   保证函数输出完全受该 seed 控制）
#' @return 含 apparent / optimism / corrected 的列表
#' @export
audit_optimism <- function(expr, surv, time_col, event_col, genes,
                           method = "lasso", B = 200, seed = 20261014) {
  d <- .align_cohort(expr, surv, time_col, event_col, genes)
  Z <- .zscore_genes(d$expr, genes)
  g <- intersect(genes, rownames(Z))
  if (length(g) < 3) stop("fewer than 3 usable signature genes")
  Z <- Z[g, , drop = FALSE]; X <- t(Z)
  t <- d$surv$t; e <- d$surv$e

  alpha <- switch(method, lasso = 1, ridge = 0, enet = 0.5, 1)
  fit_once <- function(Xa, ta, ea) {
    cv <- tryCatch(glmnet::cv.glmnet(Xa, survival::Surv(ta, ea), family = "cox",
                                     alpha = alpha, nfolds = 5,
                                     cox.ties = "breslow"), error = function(e) NULL)
    if (is.null(cv)) return(NULL)
    as.numeric(predict(cv, Xa, s = cv$lambda.min, type = "link"))
  }

  # 【修正】表观模型此前在 set.seed(seed) 之前拟合，其 CV 折用的是当时的全局随机状态，
  # 导致同一 seed 不能完全复现输出。现在把表观拟合也置于同一个 seed 之下。
  set.seed(seed)
  sc <- fit_once(X, t, e)
  apparent <- .c_index(t, e, if (is.null(sc)) rep(NA_real_, length(t)) else sc)

  set.seed(seed)
  opt <- numeric(0)
  for (i in seq_len(B)) {
    b <- sample(length(t), replace = TRUE)
    Xb <- X[b, , drop = FALSE]; tb <- t[b]; eb <- e[b]
    if (sum(eb) < 10) next
    cv <- tryCatch(glmnet::cv.glmnet(Xb, survival::Surv(tb, eb), family = "cox",
                                     alpha = alpha, nfolds = 5,
                                     cox.ties = "breslow"), error = function(e) NULL)
    if (is.null(cv)) next
    s_boot <- as.numeric(predict(cv, Xb, s = cv$lambda.min, type = "link"))
    s_orig <- as.numeric(predict(cv, X,  s = cv$lambda.min, type = "link"))
    c_boot <- .c_index(tb, eb, s_boot)
    c_orig <- .c_index(t,  e,  s_orig)
    if (!is.na(c_boot) && !is.na(c_orig)) opt <- c(opt, c_boot - c_orig)
  }
  optimism <- if (length(opt)) mean(opt) else NA_real_
  structure(list(module = "optimism", method = method, n = length(t),
                 events = sum(e), n_gene = length(g), B = length(opt),
                 apparent = apparent, optimism = optimism,
                 corrected = apparent - optimism),
            class = "sigAudit_module")
}

#' 模块 2：跨队列可迁移性——**实验性**
#'
#' \strong{注意：本模块未参与论文的任何结论，属实验性功能。}
#'
#' 训练队列拟合 -> 外部队列代入（**禁止在任何一步使用外部数据拟合**）。
#' 基因层 z-score 的均值/标准差、缺失过滤与零方差判定**全部只在训练队列上估计**，
#' 再原样应用到外部队列。
#'
#' @section 已修正（第三轮审查）:
#' 旧实现对外部队列**再次调用** `.zscore_genes()`，即用外部队列自身的均值、标准差
#' 与可用基因决定变换和特征集合 —— 这与本函数"禁止用外部数据拟合"的注释直接矛盾。
#' 现在改为 `.fit_zscore()`（仅训练侧） + `.apply_zscore()`（训练参数应用到两侧）。
#' 同时 `seed` 现在真正用于 CV 折（此前该参数未被使用，输出不受其控制）。
#'
#' @param train_expr,train_surv 训练队列
#' @param train_time_col,train_event_col 训练队列的生存时间/事件列名
#' @param test_expr,test_surv 外部队列
#' @param test_time_col,test_event_col 外部队列的生存时间/事件列名
#' @param genes 签名基因
#' @param method 学习器："lasso" | "ridge" | "enet"
#' @param seed 随机种子（CV 折）
#' @return 含 external C-index、落差、缺失基因的列表
#' @export
audit_transport <- function(train_expr, train_surv, train_time_col, train_event_col,
                            test_expr, test_surv, test_time_col, test_event_col,
                            genes, method = "lasso", seed = 20261014) {
  dtr <- .align_cohort(train_expr, train_surv, train_time_col, train_event_col, genes)
  dte <- .align_cohort(test_expr,  test_surv,  test_time_col,  test_event_col, genes)
  # 【修正】标准化参数只在训练队列上拟合，再原样应用到两个队列
  par <- .fit_zscore(dtr$expr, genes)
  Ztr <- .apply_zscore(dtr$expr, par, genes)
  Zte <- .apply_zscore(dte$expr, par, genes)
  g <- intersect(genes, intersect(rownames(Ztr), rownames(Zte)))
  missing_in_test <- setdiff(intersect(genes, rownames(Ztr)), rownames(Zte))
  if (length(g) < 3) stop("fewer than 3 signature genes shared by both cohorts")

  alpha <- switch(method, lasso = 1, ridge = 0, enet = 0.5, 1)
  Xtr <- t(Ztr[g, , drop = FALSE])
  set.seed(seed)                       # 【修正】seed 此前未被使用
  cv <- glmnet::cv.glmnet(Xtr, survival::Surv(dtr$surv$t, dtr$surv$e),
                          family = "cox", alpha = alpha, nfolds = 5,
                          cox.ties = "breslow")
  apparent <- .c_index(dtr$surv$t, dtr$surv$e,
                       as.numeric(predict(cv, Xtr, s = cv$lambda.min, type = "link")))

  Xte <- t(Zte[g, , drop = FALSE])
  ext <- .c_index(dte$surv$t, dte$surv$e,
                  as.numeric(predict(cv, Xte, s = cv$lambda.min, type = "link")))
  structure(list(module = "transportability", method = method,
                 n_gene = length(g), missing_in_external = missing_in_test,
                 apparent = apparent, external = ext, gap = apparent - ext,
                 n_train = nrow(Xtr), n_external = nrow(Xte),
                 events_external = sum(dte$surv$e),
                 standardisation = "train-only (.fit_zscore on train cohort)"),
            class = "sigAudit_module")
}

#' 模块 3：mRNA-蛋白一致性——**实验性**
#'
#' \strong{注意：本模块未参与论文的任何结论，属实验性功能。}
#'
#' 对签名基因中同时有 mRNA 与蛋白定量的部分，计算 Spearman 相关。
#' mRNA 有信号但蛋白不一致的基因，其生物学意义值得怀疑。
#'
#' @param expr mRNA 矩阵
#' @param protein 蛋白矩阵（样本 × 蛋白 或 蛋白 × 样本）
#' @param genes 签名基因
#' @param min_n 计算相关所需的最少共同样本数，默认 30
#' @return 每个基因的 rho/p，及汇总
#' @export
audit_concordance <- function(expr, protein, genes, min_n = 30) {
  em <- as.matrix(expr)
  if (!.is_genes_by_samples(em, genes)) em <- t(em)
  pm <- as.matrix(protein)
  # 统一成【样本 × 蛋白】：用标签匹配判方向，避免 nrow>ncol 这类脆弱判据
  if (.is_genes_by_samples(pm, genes)) pm <- t(pm)
  colnames(em) <- .harmonise_ids(colnames(em))
  rownames(pm) <- .harmonise_ids(rownames(pm))

  g <- intersect(genes, intersect(rownames(em), colnames(pm)))
  com <- intersect(colnames(em), rownames(pm))
  if (!length(g) || length(com) < min_n)
    return(structure(list(module = "concordance", n_gene_covered = length(g),
                          n_sample = length(com), table = NULL,
                          message = "insufficient shared genes or samples"),
                     class = "sigAudit_module"))

  tab <- do.call(rbind, lapply(g, function(x) {
    m <- as.numeric(em[x, com]); p <- as.numeric(pm[com, x])
    ok <- !is.na(m) & !is.na(p)
    if (sum(ok) < min_n) return(NULL)
    ct <- suppressWarnings(stats::cor.test(m[ok], p[ok], method = "spearman"))
    data.frame(gene = x, n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
  }))
  if (is.null(tab) || !nrow(tab)) tab <- NULL
  structure(list(module = "concordance", n_gene_covered = length(g),
                 n_sample = length(com), table = tab,
                 n_concordant = if (is.null(tab)) NA_integer_ else sum(tab$p < 0.05 & tab$rho > 0),
                 n_discordant = if (is.null(tab)) NA_integer_ else sum(!(tab$p < 0.05 & tab$rho > 0))),
            class = "sigAudit_module")
}


#' 打印审计模块结果
#' @param x 一个 sigAudit_module 对象
#' @param ... 其它参数（未使用）
#' @return 不可见地返回 x
#' @export
print.sigAudit_module <- function(x, ...) {
  cat(sprintf("<sigAudit module: %s>\n", x$module))
  # 出错时模块对象只有 module + error 两个字段。若不先拦下来，
  # 下面各个分支会把 NULL 喂给 sprintf -> 返回 character(0) -> 整行静默消失，
  # 使用者看到的是一片空白而不是错误原因。
  if (!is.null(x$error)) {
    cat("  failed: ", x$error, "\n", sep = "")
    return(invisible(x))
  }
  if (x$module == "optimism")
    cat(sprintf("  n=%d  events=%d  genes=%d  bootstrap=%d\n  apparent=%.3f  optimism=%.3f  corrected=%.3f\n",
                x$n, x$events, x$n_gene, x$B, x$apparent, x$optimism, x$corrected))
  else if (x$module == "transportability")
    cat(sprintf("  genes=%d (missing in external: %s)\n  apparent=%.3f  external=%.3f  gap=%.3f\n  n_train=%d n_external=%d\n",
                x$n_gene, if (length(x$missing_in_external)) paste(x$missing_in_external, collapse=",") else "none",
                x$apparent, x$external, x$gap, x$n_train, x$n_external))
  else if (x$module == "concordance") {
    cat(sprintf("  genes covered=%d  samples=%d\n", x$n_gene_covered, x$n_sample))
    if (!is.null(x$table)) {
      cat(sprintf("  concordant=%d  discordant=%d\n", x$n_concordant, x$n_discordant))
      print(x$table, row.names = FALSE)
    } else cat("  ", x$message %||% "", "\n", sep = "")
  }
  else if (x$module == "selection") {
    cat(sprintf("  signatures=%d  candidates=%d  reference=%s\n", x$n, x$n_candidate, x$reference))
    cat(sprintf("  pool-level Orel (reference-adjusted selection optimism) = %+.4f\n", x$pool_bias))
    d <- x$dose
    if (!is.null(d) && nrow(d)) {
      dd <- diff(d$orel)
      j <- if (length(dd)) which.max(dd) else NA_integer_
      cat(sprintf("  prefix pool curve: k=2 %+.4f -> k=%d %+.4f\n",
                  d$orel[1], d$k[nrow(d)], d$orel[nrow(d)]))
      if (!is.na(j)) cat(sprintf("    largest single step at k=%d->%d (%s) %+.4f\n",
                                 d$k[j], d$k[j + 1], d$candidate_added[j + 1], dd[j]))
      cat("    note: the prefix curve uses a FIXED background pool, not an average over subsets\n")
    }
    a <- x$entry_summary
    if (!is.null(a) && nrow(a)) {
      cat("  marginal entry effect (random background pools; averaged over k):\n")
      ag <- stats::aggregate(entry ~ target, data = a, FUN = function(z) mean(z, na.rm = TRUE))
      ag <- ag[order(-ag$entry), ]
      for (i in seq_len(nrow(ag)))
        cat(sprintf("    %-16s %+.4f\n", ag$target[i], ag$entry[i]))
      cat("    (see $entry for per-k values, sd_bg and Monte Carlo SE;\n")
      cat("     sd_bg is background-composition spread, NOT a patient-sampling CI)\n")
    }
    for (n in x$notes %||% character(0)) cat("  note: ", n, "\n", sep = "")
  }
  else if (x$module == "rules") {
    t <- x$table
    cat(sprintf("  signatures=%d  reference=%s\n", x$n, x$reference))
    cat(sprintf("  %-20s %9s %9s %9s %9s %9s\n", "rule", "Cv", "Ct", "Gv", "Gt", "Orel"))
    for (i in seq_len(nrow(t))) cat(sprintf(
      "  %-20s %9.4f %9.4f %+9.4f %+9.4f %+9.4f\n",
      t$rule[i], t$reported[i], t$test[i], t$Gv[i], t$Gt[i], t$Orel[i]))
    cat("  (Cv = validation-set value; Ct = held-out test estimate;\n")
    cat("   Orel = reference-adjusted selection optimism; Oraw is a different quantity)\n")
    for (n in x$notes %||% character(0)) cat("  note: ", n, "\n", sep = "")
  }
  invisible(x)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
