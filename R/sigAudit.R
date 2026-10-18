## ============================================================
## sigAudit.R — 主入口与可信度报告
## ============================================================

#' 审计一个 omics 预后签名的可信度
#'
#' \strong{本包经过验证的核心能力是 \code{\link{audit_selection}} 与
#' \code{\link{audit_rules}}（收益核算与规则比较）。}
#' 本函数串联的 optimism / transportability / concordance 三个模块
#' **均标记为实验性、已移出默认入口**，不出现在论文的任何结论中。
#' 若要继续使用它们，必须显式传 \code{experimental = TRUE}。
#'
#' @param genes 签名基因（字符向量）
#' @param train_expr,train_surv 训练队列表达与生存
#' @param train_time_col,train_event_col 生存列名
#' @param external 可选的列表：每个元素含 expr/surv/time_col/event_col，用作外部队列
#' @param protein 可选的蛋白矩阵，用于一致性检查
#' @param method 学习器："lasso" | "ridge" | "enet"
#' @param B bootstrap 次数
#' @param experimental 是否运行三个实验性模块（默认 FALSE：跳过并给出说明）
#' @param verbose 是否打印进度
#' @return `sigAudit_report` 对象
#' @export
sigAudit <- function(genes,
                     train_expr, train_surv, train_time_col, train_event_col,
                     external = NULL, protein = NULL,
                     method = "lasso", B = 200, experimental = FALSE, verbose = TRUE) {
  stopifnot(is.character(genes), length(genes) >= 3)
  res <- list(genes = genes, method = method, n_gene = length(genes))

  if (!experimental) {
    # 说明：本包要求 R 代码为纯 ASCII（R CMD check 的 portable-package 规则），
    # 所以消息文本用英文；中文说明保留在上面的 roxygen 注释里（注释不受此限制）。
    warning("The optimism / transportability / concordance modules are experimental and ",
            "have been removed from the default entry point; nothing was run. ",
            "For the validated core, call audit_selection() and audit_rules() directly. ",
            "Pass experimental = TRUE to run the legacy modules anyway.",
            call. = FALSE)
    res$experimental_skipped <- TRUE
    res$optimism <- NULL
    res$transport <- list()
    res$concordance <- NULL
    return(structure(res, class = "sigAudit_report"))
  }

  if (verbose) message("[1/3] optimism correction (experimental) ...")
  res$optimism <- tryCatch(
    audit_optimism(train_expr, train_surv, train_time_col, train_event_col,
                   genes, method = method, B = B),
    error = function(e) structure(list(module="optimism", error=conditionMessage(e)),
                                  class = "sigAudit_module"))

  res$transport <- list()
  if (!is.null(external)) {
    for (i in seq_along(external)) {
      e <- external[[i]]
      if (verbose) message(sprintf("[2/3] external transportability -- cohort %d/%d ...", i, length(external)))
      res$transport[[e$name %||% paste0("cohort", i)]] <- tryCatch(
        audit_transport(train_expr, train_surv, train_time_col, train_event_col,
                        e$expr, e$surv, e$time_col, e$event_col, genes, method = method),
        error = function(err) structure(list(module="transportability",
                                             error=conditionMessage(err)),
                                        class = "sigAudit_module"))
    }
  }

  if (!is.null(protein)) {
    if (verbose) message("[3/3] mRNA-protein concordance (experimental) ...")
    res$concordance <- tryCatch(
      audit_concordance(train_expr, protein, genes),
      error = function(e) structure(list(module="concordance", error=conditionMessage(e)),
                                    class = "sigAudit_module"))
  }
  structure(res, class = "sigAudit_report")
}

#' 汇总各项审计估计值与数据覆盖（**不再给出 A–D 总评级**）
#'
#' 为什么移除字母分级：原实现从 0 分起、模块失败或缺失只写 notes 而不加分，
#' 因此**一份没有任何有效证据的报告也能得到最高等级 A**；而且 0.05 / 0.60 /
#' 0.10 / 0.55 等阈值声称"经真实数据标定"却没有标定证据；
#' 外部队列越多、可扣分的机会越多，等级在不同研究之间也不可比。
#'
#' 现在改为直接报告：① 每个模块的**估计值**；② 每个模块的**数据覆盖**
#' （有多少样本/基因/任务真正参与了计算）；③ 失败原因。
#' **不做跨模块加总，也不产生等级标签。**
#'
#' @param report 由 \code{sigAudit()} 返回的 \code{sigAudit_report} 对象；
#'   也可手动附加 \code{report$selection}（\code{audit_selection()} 的结果）。
#' @return 含 estimates / coverage / notes 的列表
#' @export
grade_credibility <- function(report) {
  est <- list(); cov <- list(); notes <- character(0)

  o <- report$optimism
  if (!is.null(o) && is.null(o$error)) {
    est$apparent <- o$apparent; est$optimism <- o$optimism; est$corrected <- o$corrected
    cov$optimism <- sprintf("n=%s, events=%s, genes=%s, bootstrap=%s",
                            o$n, o$events, o$n_gene, o$B)
  } else if (!is.null(o) && !is.null(o$error)) {
    notes <- c(notes, sprintf("optimism module failed: %s", o$error))
  }

  for (nm in names(report$transport)) {
    t <- report$transport[[nm]]
    if (is.null(t)) next
    if (!is.null(t$error)) { notes <- c(notes, sprintf("%s module failed: %s", nm, t$error)); next }
    est[[paste0("external_", nm)]] <- t$external
    cov[[paste0("transport_", nm)]] <- sprintf(
      "n_train=%s, n_external=%s, events_external=%s, genes=%s, missing_in_external=%s",
      t$n_train, t$n_external, t$events_external, t$n_gene, length(t$missing_in_external))
  }

  co <- report$concordance
  if (!is.null(co) && is.null(co$error)) {
    est$concordance_n_concordant  <- co$n_concordant
    est$concordance_n_discordant  <- co$n_discordant
    cov$concordance <- sprintf("genes_covered=%s, samples=%s", co$n_gene_covered, co$n_sample)
  } else if (!is.null(co) && !is.null(co$error)) {
    notes <- c(notes, sprintf("concordance module failed: %s", co$error))
  }

  # 模型选择偏倚：报告池级偏倚（可复现的量）与有效任务数，不做等级加总。
  se <- report$selection
  if (!is.null(se) && is.null(se$error) && !is.null(se$pool_bias)) {
    est$pool_bias <- se$pool_bias
    cov$selection <- sprintf("signatures=%s, candidates=%s", se$n, se$n_candidate)
    if (!is.null(se$entry_summary)) {
      notes <- c(notes,
        "entry effects are marginal changes under a fixed background-pool design;",
        "their per-candidate ranking is only moderately reproducible across data splits",
        "and should not be used to name a single offending candidate.")
    }
  } else if (!is.null(report$selection) && !is.null(report$selection$error)) {
    notes <- c(notes, sprintf("selection-bias module failed: %s", report$selection$error))
  }

  structure(list(estimates = est, coverage = cov, notes = notes),
            class = "sigAudit_summary")
}

#' 打印审计汇总
#' @param x 一个 sigAudit_summary 对象
#' @param ... 其它参数（未使用）
#' @return 不可见地返回 x
#' @export
print.sigAudit_summary <- function(x, ...) {
  cat("--- estimates ---\n")
  if (!length(x$estimates)) cat("  (none)\n") else
    for (nm in names(x$estimates))
      cat(sprintf("  %-28s %s\n", nm,
                  if (is.numeric(x$estimates[[nm]]))
                    sprintf("%.4f", x$estimates[[nm]]) else as.character(x$estimates[[nm]])))
  cat("--- coverage (what actually entered the computation) ---\n")
  if (!length(x$coverage)) cat("  (none)\n") else
    for (nm in names(x$coverage)) cat(sprintf("  %-28s %s\n", nm, x$coverage[[nm]]))
  if (length(x$notes)) { cat("--- notes ---\n"); for (n in x$notes) cat("  - ", n, "\n", sep = "") }
  cat("\n(no letter grade is produced: module failures or absence of evidence\n",
      " must not yield a better label, and the thresholds were never calibrated)\n", sep = "")
  invisible(x)
}

#' 打印审计报告
#' @param x 一个 sigAudit_report 对象
#' @param ... 其它参数（未使用）
#' @return 不可见地返回 x
#' @export
print.sigAudit_report <- function(x, ...) {
  cat("===============================================\n")
  cat(" sigAudit credibility report\n")
  cat("===============================================\n")
  cat(sprintf("signature genes: %d   learner: %s\n\n", x$n_gene, x$method))
  for (nm in c("optimism", "concordance")) {
    if (!is.null(x[[nm]])) { print(x[[nm]]); cat("\n") }
  }
  if (length(x$transport)) {
    cat("--- external transportability ---\n")
    for (nm in names(x$transport)) {
      cat(sprintf("[%s]\n", nm)); print(x$transport[[nm]])
    }
  }
  g <- grade_credibility(x)
  cat("\n-----------------------------------------------\n")
  print(g)
  cat("-----------------------------------------------\n")
  invisible(x)
}
