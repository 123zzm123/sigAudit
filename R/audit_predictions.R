## ============================================================
## audit_predictions.R — 预测级审计入口
##
## 2026-10-19 新增。动机（本次修订）：原 audit_rules() 只接受 C-index 矩阵，
## 无法从各候选**已经生成的验证/测试风险分数**直接重建集成或审计选择程序。
## 本函数把研究脚本（scripts/_lib_analysis.R 的 rules_one / ens_pred /
## make_prio、scripts/20/29 的三层级汇总口径）中已被验证的逻辑提取为一个
## 可复用的预测级入口，使外部用户用自己的风险分数即可运行同一套审计。
##
## 复用来源映射：
##   · C-index            <- scripts/_lib_analysis.R::ci_direct（固定方向，无重拟合）
##   · 并列处理           <- R/utils.R::.fixed_prio（与候选名字绑定，与列顺序无关）
##   · 集成标准化         <- scripts/_lib_analysis.R::ens_pred（验证侧均值/标准差）
##   · 规则集合           <- scripts/_lib_analysis.R::rules_one
##   · 三层级汇总         <- scripts/29_baseline_sensitivity.R 的 strat_rule 口径
##
## 2026-10-19 修订（大修后评审）：统一规则注册与实际实现（top_k 真正生效、
## 家族规则可用、未知规则报错）；修复 by_cohort 的加权口径；补齐输入契约
## （事件长度/取值、候选名唯一、ID 唯一且验证/测试不重叠）；按任务级共同
## 掩码处理缺失并记录丢弃数；汇总只在共同可评估集合上进行，保证
## Orel = Gv − Gt 在汇总层仍成立；默认规则不含 oracle_test。
## ============================================================

## ---- 内部：固定方向的 Harrell C-index（与 scripts/_lib_analysis.R 的 ci_direct 一致）----
#' @keywords internal
#' @noRd
.ci_direct <- function(t, e, s, min_n = 10L) {
  s <- tryCatch(as.numeric(s), error = function(x) NULL)
  if (is.null(s) || length(s) != length(t)) return(NA_real_)
  ok <- is.finite(s) & !is.na(t) & !is.na(e)
  if (sum(ok) < min_n) return(NA_real_)
  tryCatch(as.numeric(survival::concordance(
      survival::Surv(t[ok], e[ok]) ~ s[ok], reverse = TRUE, timewt = "n")$concordance),
    error = function(x) NA_real_)
}

## ---- 内部：验证侧标准化 -> 取平均（与 ens_pred 一致；sd=0 时置 1）----
#' @keywords internal
#' @noRd
.ens_pred_local <- function(Pv, Pt, cols) {
  cols <- cols[!is.na(cols)]
  if (!length(cols)) return(NULL)
  Mv <- Pv[, cols, drop = FALSE]; Mt <- Pt[, cols, drop = FALSE]
  mu <- colMeans(Mv); sdv <- apply(Mv, 2, stats::sd)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  list(v = rowMeans(sweep(sweep(Mv, 2, mu, "-"), 2, sdv, "/")),
       t = rowMeans(sweep(sweep(Mt, 2, mu, "-"), 2, sdv, "/")))
}

## ---- 内部：判断是否为"单个任务"（含 val_pred 的 list）----
#' @keywords internal
#' @noRd
.is_single_task <- function(x) {
  is.list(x) && !is.null(x$val_pred) && !is.null(x$test_pred)
}

## ---- 内部：按任务级共同掩码清理一侧（同一患者集合上评估所有候选）----
## 我们不在每个候选上各自删除缺失（那会破坏配对比较）。这里定义共同患者掩码：
## 只有"随访时间有限 且 该患者在**所有**候选上的分数都有限"才保留；丢弃数会被记录。
#' @keywords internal
#' @noRd
.mask_side <- function(P, t, e, id) {
  P <- as.matrix(P)
  ok <- is.finite(t) & rowSums(!is.finite(P)) == 0L
  list(P = P[ok, , drop = FALSE], t = t[ok], e = e[ok],
       id = if (is.null(id)) NULL else id[ok],
       n_dropped = as.integer(sum(!ok)))
}

## ---- 内部：校验并规范化一个任务 ----
#' @keywords internal
#' @noRd
.norm_task <- function(z, i = NULL) {
  need <- c("val_pred", "test_pred", "time_val", "event_val", "time_test", "event_test")
  miss <- setdiff(need, names(z))
  if (length(miss)) stop("task ", i %||% 1L, " is missing: ", paste(miss, collapse = ", "))
  Pv <- as.matrix(z$val_pred); Pt <- as.matrix(z$test_pred)
  if (!is.numeric(Pv) || !is.numeric(Pt)) stop("val_pred/test_pred must be numeric matrices")
  if (is.null(colnames(Pv)) || is.null(colnames(Pt)))
    stop("val_pred/test_pred must have column names (candidate identifiers)")
  if (anyDuplicated(colnames(Pv))) stop("val_pred has duplicated candidate names")
  if (anyDuplicated(colnames(Pt))) stop("test_pred has duplicated candidate names")
  if (!identical(colnames(Pv), colnames(Pt)))
    stop("val_pred and test_pred must have identical candidate column names in the same order")
  tv <- as.numeric(z$time_val); ev <- as.numeric(z$event_val)
  tt <- as.numeric(z$time_test); et <- as.numeric(z$event_test)
  if (length(tv) != nrow(Pv)) stop("length(time_val) != nrow(val_pred)")
  if (length(ev) != nrow(Pv)) stop("length(event_val) != nrow(val_pred)")
  if (length(tt) != nrow(Pt)) stop("length(time_test) != nrow(test_pred)")
  if (length(et) != nrow(Pt)) stop("length(event_test) != nrow(test_pred)")
  if (anyNA(ev) || !all(ev %in% c(0, 1))) stop("event_val must be 0/1 with no missing values")
  if (anyNA(et) || !all(et %in% c(0, 1))) stop("event_test must be 0/1 with no missing values")
  vid <- z$val_ids; tid <- z$test_ids
  if (!is.null(vid)) {
    if (length(vid) != nrow(Pv)) stop("length(val_ids) != nrow(val_pred)")
    if (anyNA(vid) || any(!nzchar(as.character(vid)))) stop("val_ids must be non-missing and non-empty")
    if (anyDuplicated(vid)) stop("val_ids must be unique")
  }
  if (!is.null(tid)) {
    if (length(tid) != nrow(Pt)) stop("length(test_ids) != nrow(test_pred)")
    if (anyNA(tid) || any(!nzchar(as.character(tid)))) stop("test_ids must be non-missing and non-empty")
    if (anyDuplicated(tid)) stop("test_ids must be unique")
  }
  if (!is.null(vid) && !is.null(tid) &&
      length(intersect(as.character(vid), as.character(tid))))
    stop("val_ids and test_ids must not overlap")
  mv <- .mask_side(Pv, tv, ev, vid)
  mt <- .mask_side(Pt, tt, et, tid)
  ## 更新时间/事件后重新写回清理结果；事件已在上面校验
  list(val_pred = mv$P, test_pred = mt$P, time_val = mv$t, event_val = mv$e,
       time_test = mt$t, event_test = mt$e,
       n_drop_val = mv$n_dropped, n_drop_test = mt$n_dropped,
       name = z$name %||% paste0("task", i %||% 1L),
       cohort = z$cohort %||% NA_character_, split = z$split %||% NA_character_,
       val_ids = mv$id, test_ids = mt$id)
}

## 帮助：%||%（包内 utils 已有，这里兜底以免依赖顺序）
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

#' 预测级审计：由各候选的验证/测试风险分数核算搜索增益
#'
#' 给定每个任务上各候选**已经生成的**验证集与测试集风险分数、结局以及患者标识，
#' 本函数重算候选与规则的 C-index、参照校正的选择乐观度 \code{Orel = Gv - Gt}、
#' 固定参照对照、常数/重复诊断，并按研究口径给出分层汇总。
#'
#' \strong{这不是一个训练器}：它不拟合任何模型，也不下载数据；它只审计已有预测，
#' 因此任何产生风险分数的流程都能直接接入。
#'
#' @param tasks 单个任务对象，或任务对象列表。每个任务为一个 list，至少含：
#'   \describe{
#'     \item{val_pred,test_pred}{风险分数矩阵，行 = 患者（验证集与测试集人数可以不同），
#'       列 = 候选；两侧的候选列名必须一致且唯一（顺序也须一致）}
#'     \item{time_val,event_val}{验证集随访时间与事件指示（0/1，无缺失）}
#'     \item{time_test,event_test}{测试集随访时间与事件指示（0/1，无缺失）}
#'   }
#'   可选：\code{cohort} 与 \code{split}（分组标签，用于三层级汇总）、
#'   \code{name}、\code{val_ids} 与 \code{test_ids}（患者标识；若提供则须唯一、
#'   非缺失，且验证/测试两侧不得重叠）。行顺序由调用者保证对齐；本函数不会
#'   重新排序患者。默认要求预测与结局在任务级**共同对齐且有限**：若存在缺失，
#'   只在"随访时间有限且该患者在所有候选上分数均有限"时保留该患者，丢弃数
#'   记录在 \code{diagnostics} 与 \code{coverage} 中（不会按候选各自删除缺失）。
#' @param reference 固定参照的候选取名，默认 \code{"Lasso"}。
#' @param rules 要输出的规则名；\code{NULL}（默认）表示全部**可实施**规则
#'   （见 Details），不含 \code{oracle_test}。若给出未知规则名会报错，
#'   不会静默丢弃；\code{linear_only}/\code{tree_only} 需要同时提供 \code{family}。
#' @param top_k 主 top-k 集成的 k，默认 5。须为正整数；该 k 会被真正计算，
#'   并与 top-2/3/5/10 一起构成可用的集成规则集合。若候选数小于 k，
#'   集成使用 \code{top-min(k, K)}。
#' @param family 可选命名向量：候选名 -> 家族标签（\code{"linear"}/\code{"tree"}）；
#'   仅当请求 \code{linear_only}/\code{tree_only} 时使用。
#' @param seed 并列优先级所用的随机种子（绑定候选名称，与输入列顺序无关）。
#' @param seed_resample 三条重采样规则（majority_vote / worst_half / penalized）
#'   所用的独立随机种子，默认 1（与研究脚本 `rules_one()` 的 `set.seed(1)` 口径一致；
#'   把它与 \code{seed} 分开，才能在更换划分种子时仍精确复现流水线）。
#' @param deduplicate 逻辑值：若为 \code{TRUE}，先在**验证预测完全相同**上分组去重，
#'   组内取候选名称字典序最小者为代表，再在去重池上重算规则。默认 \code{FALSE}。
#'   注意判据只是验证预测相等，不代表这些候选的测试预测或拟合模型完全相同。
#' @param B_vote,B_half \code{majority_vote} / \code{worst_half} / \code{penalized}
#'   三条重采样规则的次数。
#' @return \code{sigAudit_prediction} 对象，含：
#'   \describe{
#'     \item{candidates}{逐 任务×候选 的 \code{Cv}/\code{Ct}/\code{constant}}
#'     \item{rules}{逐 任务×规则 的 \code{Cv}/\code{Ct}/\code{Gv}/\code{Gt}/\code{Orel}、
#'       \code{winner}、\code{members}、\code{evaluable}（该任务-规则是否可评估）}
#'     \item{diagnostics}{逐任务的候选池冗余诊断（唯一候选数、常数候选数、重复组数、
#'       丢弃患者数等）}
#'     \item{summary}{按规则的宏汇总；仅在**共同可评估集合**上计算，
#'       因此 \code{Orel = Gv - Gt} 在汇总层仍成立}
#'     \item{by_cohort}{按队列的汇总（提供 cohort 时），与 \code{summary} 使用同一
#'       中间层（先队列×划分、再队列），不是任务层直接平均}
#'     \item{coverage}{任务数、规则记录数、可评估/不可评估计数、丢弃患者数}
#'   }
#' @details 规则口径与配套研究脚本一致：
#'   \itemize{
#'     \item \code{none}：固定参照本身（\code{Gv = Gt = Orel = 0}）。
#'     \item \code{argmax_val}：在全池（含参照）中按验证 C-index 选优。
#'     \item \code{argmax_nonreference}：仅在非参照候选中选优（敏感性）。
#'     \item \code{runner_up}：验证 C-index 次优。
#'     \item \code{top2/top3/top5/top10_ensemble}（以及 \code{top_k} 对应的 k）：
#'       按验证 C-index 取前 k 做集成；集成用验证侧均值/标准差标准化后取平均
#'       （不使用测试侧信息）。
#'     \item \code{all_ensemble}：全池集成。
#'     \item \code{linear_only}/\code{tree_only}：仅在 \code{family} 下生成。
#'     \item \code{majority_vote}/\code{worst_half}/\code{penalized}：验证侧重采样类补救规则。
#'     \item \code{oracle_test}：用测试结局选优的**事后描述性参照**，不可实施；
#'       默认不返回，需显式请求。
#'   }
#'   汇总口径：**任务 -> (队列×划分) 平均 -> 队列内平均 -> 队列等权**；
#'   缺少 cohort/split 标签时退化为任务级等权平均，并在 \code{notes} 中提示。
#' @examples
#' ## 一个极小示例：2 个候选 + 1 个参照，验证/测试风险分数与结局
#' set.seed(1)
#' val <- cbind(Lasso = rnorm(40), A = rnorm(40), B = rnorm(40))
#' test <- cbind(Lasso = rnorm(40), A = rnorm(40), B = rnorm(40))
#' tk <- list(list(val_pred = val, test_pred = test,
#'                 time_val = rexp(40), event_val = rbinom(40, 1, .6),
#'                 time_test = rexp(40), event_test = rbinom(40, 1, .6)))
#' res <- audit_predictions(tk, reference = "Lasso")
#' res$summary
#' @export
audit_predictions <- function(tasks, reference = "Lasso", rules = NULL, top_k = 5L,
                              family = NULL, seed = 1L, seed_resample = 1L,
                              deduplicate = FALSE, B_vote = 12L, B_half = 6L) {
  if (.is_single_task(tasks)) tasks <- list(tasks)
  if (!is.list(tasks) || !length(tasks)) stop("tasks must be a task list")

  ## ---- 规则注册与实际实现统一 ----
  if (length(top_k) != 1L || is.na(top_k) || !is.finite(top_k) ||
      top_k != as.integer(top_k) || as.integer(top_k) < 1L)
    stop("top_k must be a single positive integer")
  top_k <- as.integer(top_k)
  ens_ks    <- sort(unique(c(2L, 3L, 5L, 10L, top_k)))
  ens_rules <- paste0("top", ens_ks, "_ensemble")
  USABLE <- c("none", "argmax_val", "argmax_nonreference", "runner_up",
              ens_rules, "all_ensemble",
              "majority_vote", "worst_half", "penalized")
  OPT_IN <- "oracle_test"
  FAMILY_RULES <- c("linear_only", "tree_only")
  AVAIL <- c(USABLE, OPT_IN, FAMILY_RULES)
  if (is.null(rules)) {
    rules <- USABLE
  } else {
    rules <- unique(as.character(rules))
    bad <- setdiff(rules, AVAIL)
    if (length(bad))
      stop("unknown rule(s): ", paste(bad, collapse = ", "),
           "\navailable rules: ", paste(AVAIL, collapse = ", "))
    if (any(FAMILY_RULES %in% rules) && is.null(family))
      stop("linear_only/tree_only require the `family` argument (candidate -> family label)")
    rules <- unique(c("none", rules))
  }
  if (!is.null(family)) {
    if (is.null(names(family))) stop("family must be a named vector (candidate -> label)")
    if (anyNA(family)) stop("family has missing labels")
  }

  TK <- lapply(seq_along(tasks), function(i) .norm_task(tasks[[i]], i))

  cand_rows <- list(); rule_rows <- list(); diag_rows <- list()
  for (zk in TK) {
    Pv <- zk$val_pred; Pt <- zk$test_pred; mods <- colnames(Pv)
    if (!reference %in% mods)
      stop("reference '", reference, "' is not among candidates of task ", zk$name)
    ## 全池候选指标 + 常数诊断
    Cv <- vapply(mods, function(m) .ci_direct(zk$time_val, zk$event_val, Pv[, m]), numeric(1))
    Ct <- vapply(mods, function(m) .ci_direct(zk$time_test, zk$event_test, Pt[, m]), numeric(1))
    const <- vapply(mods, function(m) length(unique(Pv[, m])) == 1L, logical(1))
    names(Cv) <- names(Ct) <- names(const) <- mods
    cand_rows[[length(cand_rows) + 1]] <- data.frame(
      task = zk$name, cohort = zk$cohort, split = zk$split, model = mods,
      Cv = Cv, Ct = Ct, constant = const, stringsAsFactors = FALSE)
    ## 冗余诊断（基于验证预测精确相等）
    key <- apply(Pv, 2, function(col) paste(format(col, digits = 17), collapse = ","))
    grp <- split(mods, key)
    n_dup_group <- sum(vapply(grp, length, integer(1)) > 1L)
    n_unique <- length(grp)
    ## 池：原池或去重池（组内名称字典序最小者为代表）
    if (deduplicate) {
      pool <- sort(unique(vapply(grp, function(g) sort(g)[1], character(1))))
    } else {
      pool <- mods
    }
    prio <- .fixed_prio(mods, seed)          # 与主分析一致：名字绑定、列顺序无关
    pv <- Cv[pool]
    pick1 <- function(cols) { cols <- cols[!is.na(cols)]
      if (!length(cols)) return(NA_character_)
      m <- suppressWarnings(max(pv[cols], na.rm = TRUE))
      if (!is.finite(m)) return(NA_character_)
      cc <- cols[is.finite(pv[cols]) & pv[cols] == m]; cc[which.max(prio[cc])] }
    member_k <- function(cols, k) {
      cols <- cols[!is.na(cols)]
      o <- order(-pv[cols], -prio[cols], na.last = NA)
      if (!length(o)) return(character(0))
      cols[o][seq_len(min(k, length(o)))] }
    ## 重采样规则（验证侧）。仅当确实请求时才计算，避免无谓开销。
    ## 重采样 RNG 独立于并列优先级 seed，默认 1，
    ## 与 scripts/_lib_analysis.R::rules_one 的 set.seed(1) 口径一致。
    need_resamp <- any(c("majority_vote", "worst_half", "penalized") %in% rules)
    if (need_resamp) {
      set.seed(seed_resample); wins <- stats::setNames(numeric(length(pool)), pool)
      for (b in seq_len(B_vote)) {
        ib <- sample(length(zk$time_val), replace = TRUE)
        if (sum(zk$event_val[ib]) < 5) next
        vb <- vapply(pool, function(m) .ci_direct(zk$time_val[ib], zk$event_val[ib], Pv[ib, m]), numeric(1))
        if (all(!is.finite(vb))) next
        jj <- pool[pick_max_local(vb, prio[pool])]; wins[jj] <- wins[jj] + 1
      }
      set.seed(seed_resample); acc <- stats::setNames(vector("list", length(pool)), pool)
      for (b in seq_len(B_half)) {
        pm <- sample(length(zk$time_val)); h1 <- pm[seq_len(floor(length(zk$time_val) / 2))]
        h2 <- setdiff(pm, h1)
        for (m in pool) {
          x1 <- .ci_direct(zk$time_val[h1], zk$event_val[h1], Pv[h1, m])
          x2 <- .ci_direct(zk$time_val[h2], zk$event_val[h2], Pv[h2, m])
          if (is.finite(x1) && is.finite(x2)) acc[[m]] <- c(acc[[m]], x1, x2)
        }
      }
      hmin <- vapply(acc, function(z2) if (length(z2) >= 4) {
        k2 <- length(z2) %/% 2 * 2; mean(pmin(z2[seq(1, k2, 2)], z2[seq(2, k2, 2)])) } else NA_real_, numeric(1))
      hpen <- vapply(acc, function(z2) if (length(z2) >= 4) mean(z2) - stats::sd(z2) else NA_real_, numeric(1))
    }

    ## ---- 只计算被请求的规则 ----
    sel <- list(none = reference)
    if ("argmax_val" %in% rules)          sel$argmax_val <- pick1(pool)
    if ("argmax_nonreference" %in% rules) sel$argmax_nonreference <- pick1(setdiff(pool, reference))
    if ("runner_up" %in% rules) {
      o <- order(-pv[pool], -prio[pool], na.last = NA)
      sel$runner_up <- if (length(o) >= 2) pool[o[2]] else NA_character_
    }
    for (nm in intersect(ens_rules, rules)) {
      k <- as.integer(sub("^top(\\d+)_ensemble$", "\\1", nm))
      mem <- member_k(pool, k)
      pr <- .ens_pred_local(Pv, Pt, mem)
      sel[[nm]] <- if (is.null(pr)) NA_character_ else list(ens = pr, mem = mem)
    }
    if ("all_ensemble" %in% rules) {
      pr <- .ens_pred_local(Pv, Pt, pool)
      sel$all_ensemble <- if (is.null(pr)) NA_character_ else list(ens = pr, mem = pool)
    }
    if (any(FAMILY_RULES %in% rules)) {
      tre <- pool[pool %in% names(family)[family == "tree"]]
      lin <- setdiff(pool, tre)
      if ("linear_only" %in% rules) sel$linear_only <- if (length(lin) >= 2) pick1(lin) else NA_character_
      if ("tree_only"   %in% rules) sel$tree_only   <- if (length(tre) >= 2) pick1(tre) else NA_character_
    }
    if (need_resamp) {
      if ("majority_vote" %in% rules && sum(wins) > 0)
        sel$majority_vote <- pool[pick_max_local(wins[pool], prio[pool])]
      if ("worst_half" %in% rules && any(is.finite(hmin)))
        sel$worst_half <- pool[pick_max_local(hmin, prio[pool])]
      if ("penalized" %in% rules && any(is.finite(hpen)))
        sel$penalized <- pool[pick_max_local(hpen, prio[pool])]
    }
    if ("oracle_test" %in% rules) {
      ot <- Ct[pool]; m <- suppressWarnings(max(ot, na.rm = TRUE))
      sel$oracle_test <- if (!is.finite(m)) NA_character_ else {
        cc <- pool[is.finite(ot) & ot == m]; cc[which.max(prio[cc])] }
    }

    ## ---- 逐规则记录（共同可评估集合：Cv/Ct/refCv/refCt 均有限）----
    refCv <- Cv[[reference]]; refCt <- Ct[[reference]]
    for (rn in rules) {
      v <- sel[[rn]]
      if (is.list(v)) {
        cvv <- .ci_direct(zk$time_val, zk$event_val, v$ens$v)
        ctt <- .ci_direct(zk$time_test, zk$event_test, v$ens$t)
        win <- "ensemble"; mem <- paste(v$mem, collapse = "|")
      } else if (is.null(v) || is.na(v) || !v %in% mods) {
        cvv <- ctt <- NA_real_; win <- mem <- NA_character_
      } else {
        cvv <- Cv[[v]]; ctt <- Ct[[v]]; win <- v; mem <- v
      }
      ev <- is.finite(cvv) && is.finite(ctt) && is.finite(refCv) && is.finite(refCt)
      rule_rows[[length(rule_rows) + 1]] <- data.frame(
        task = zk$name, cohort = zk$cohort, split = zk$split, rule = rn,
        Cv = cvv, Ct = ctt,
        Gv = if (ev) cvv - refCv else NA_real_,
        Gt = if (ev) ctt - refCt else NA_real_,
        Orel = if (ev) (cvv - refCv) - (ctt - refCt) else NA_real_,
        winner = win, members = mem, evaluable = ev, stringsAsFactors = FALSE)
    }
    ## 诊断
    tk5 <- member_k(pool, as.integer(top_k))
    aw <- if (is.null(sel$argmax_val)) NA_character_ else sel$argmax_val
    diag_rows[[length(diag_rows) + 1]] <- data.frame(
      task = zk$name, cohort = zk$cohort, split = zk$split,
      n_candidates = length(mods), n_unique = n_unique, n_constant = sum(const),
      n_dup_groups = n_dup_group,
      n_drop_val = zk$n_drop_val, n_drop_test = zk$n_drop_test,
      argmax_winner = aw,
      argmax_is_constant = !is.na(aw) && aw %in% mods && isTRUE(const[[aw]]),
      topk = as.integer(top_k), topk_members = paste(tk5, collapse = "|"),
      topk_has_duplicate = any(duplicated(key[match(tk5, mods)])),
      stringsAsFactors = FALSE)
  }
  CAND <- do.call(rbind, cand_rows); RUL <- do.call(rbind, rule_rows); DIA <- do.call(rbind, diag_rows)
  rownames(CAND) <- NULL; rownames(RUL) <- NULL; rownames(DIA) <- NULL

  ## ---- 汇总：只在共同可评估集合上计算（保证 Orel = Gv − Gt 在汇总层成立）----
  MC <- c("Cv", "Ct", "Gv", "Gt", "Orel")
  RE <- RUL[RUL$evaluable, , drop = FALSE]
  has_grp <- !anyNA(RUL$cohort) && !anyNA(RUL$split) &&
    any(nzchar(as.character(RUL$cohort))) && any(nzchar(as.character(RUL$split)))
  if (nrow(RE)) {
    if (has_grp) {
      s1 <- stats::aggregate(RE[MC], by = RE[c("cohort", "split", "rule")], FUN = mean)
      s2 <- stats::aggregate(s1[MC], by = s1[c("cohort", "rule")], FUN = mean)
      SUM <- stats::aggregate(s2[MC], by = s2["rule"], FUN = mean)
      BYC <- s2                     # 队列层：先队列×划分、再队列（不是任务层平均）
    } else {
      SUM <- stats::aggregate(RE[MC], by = RE["rule"], FUN = mean)
      BYC <- NULL
    }
  } else {
    SUM <- data.frame(rule = character(0), Cv = numeric(0), Ct = numeric(0),
                      Gv = numeric(0), Gt = numeric(0), Orel = numeric(0))
    BYC <- NULL
  }
  ntask <- stats::aggregate(list(n_tasks = RE$Cv), by = RE["rule"], FUN = length)
  SUM$n_tasks <- ntask$n_tasks[match(SUM$rule, ntask$rule)]
  SUM$n_tasks[is.na(SUM$n_tasks)] <- 0L

  notes <- c(
    "Orel is reference-adjusted selection optimism (Gv - Gt), not 'pure selection bias'",
    "Ct is a held-out test estimate, not a known truth",
    "oracle_test is a post-hoc test-set maximum and is NOT implementable (opt-in only)",
    "the reference is kept as its own output even if it duplicates another candidate",
    "summary is computed on the common evaluable set (Cv, Ct, reference Cv and Ct all finite)")
  if (!has_grp) notes <- c(notes,
    "no cohort/split labels supplied: tasks were aggregated with equal weight; if tasks share patients this is not a valid unit of replication")
  if (deduplicate) notes <- c(notes,
    "candidates were deduplicated on EXACT equality of validation predictions; representatives chosen by lexicographic candidate name (this does NOT imply identical test predictions or identical fitted models)")
  if (any(DIA$n_drop_val > 0) || any(DIA$n_drop_test > 0))
    notes <- c(notes, "some patients were dropped by the task-level common finite mask; see diagnostics/coverage")

  structure(list(module = "predictions", reference = reference, rules_run = rules,
                 deduplicate = deduplicate,
                 candidates = CAND, rules = RUL, diagnostics = DIA,
                 summary = SUM, by_cohort = BYC,
                 coverage = c(n_tasks = length(TK),
                              n_rule_records = nrow(RUL),
                              n_rule_evaluable = sum(RUL$evaluable),
                              n_rule_failed = sum(!RUL$evaluable),
                              n_drop_val = sum(DIA$n_drop_val),
                              n_drop_test = sum(DIA$n_drop_test)),
                 notes = notes),
            class = c("sigAudit_prediction", "sigAudit_module"))
}

## 内部：与 .fixed_prio 一致的取最大（返回下标），供重采样规则使用
#' @keywords internal
#' @noRd
pick_max_local <- function(v, prio) {
  m <- suppressWarnings(max(v, na.rm = TRUE))
  if (!is.finite(m)) return(NA_integer_)
  cand <- which(is.finite(v) & v == m)
  cand[which.max(prio[cand])]
}

#' 打印预测级审计结果
#' @param x 一个 sigAudit_prediction 对象
#' @param ... 其它参数（未使用）
#' @return 不可见地返回 x
#' @export
print.sigAudit_prediction <- function(x, ...) {
  cat("===============================================\n")
  cat(" sigAudit prediction-level audit\n")
  cat("===============================================\n")
  cat(sprintf("reference: %s   tasks: %d   deduplicated: %s\n\n",
              x$reference, x$coverage[["n_tasks"]], x$deduplicate))
  cat("--- candidate-pool diagnostics (per task means) ---\n")
  for (nm in c("n_candidates", "n_unique", "n_constant", "n_dup_groups",
               "n_drop_val", "n_drop_test")) {
    cat(sprintf("  %-18s %.3f\n", nm, mean(x$diagnostics[[nm]], na.rm = TRUE)))
  }
  cat(sprintf("  argmax winner constant in %.1f%% of tasks\n",
              100 * mean(x$diagnostics$argmax_is_constant, na.rm = TRUE)))
  cat(sprintf("  evaluable rule records: %d / %d\n",
              x$coverage[["n_rule_evaluable"]], x$coverage[["n_rule_records"]]))
  cat("\n--- rule summary (Gv/Gt/Orel vs the fixed reference) ---\n")
  if (nrow(x$summary)) {
    s <- x$summary[order(-x$summary$Ct), ]
    for (i in seq_len(nrow(s)))
      cat(sprintf("  %-22s Ct=%.4f  Gv=%+.4f  Gt=%+.4f  Orel=%+.4f  n=%d\n",
                  s$rule[i], s$Ct[i], s$Gv[i], s$Gt[i], s$Orel[i], s$n_tasks[i]))
  } else {
    cat("  (no evaluable rule records)\n")
  }
  cat("\n--- notes ---\n"); for (n in x$notes) cat("  - ", n, "\n", sep = "")
  invisible(x)
}
