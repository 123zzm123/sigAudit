## ============================================================
## audit_selection.R — 模块 4/5：选择乐观度核算与规则比较
##
## 本文件于 2026-10-16 从 audit_modules.R 拆出，并于 2026-10-17 按**第三轮**
## 外部评审意见再度修改：
##   · 并列处理：优先级改由 `.fixed_prio()` 在**候选名排序后的规范集合**上生成，
##     再由名字映射。此前 `setNames(seq_along(mods), mods)` 把优先级绑定在**输入列顺序**上，
##     等价于偏向排在后面的候选；交换候选列会改变获胜者与报告的测试性能（已给反例）。
##   · `linear_only`：定义统一为**全部非树候选、含参照**（与 scripts/20 一致）。
##     原文"两侧均不含参照，条件对称"与实现不符，已删除该表述；
##     家族规则改为**必须显式提供 `family`** 才生成，避免在无家族信息时
##     把 `linear_only` 退化成 `argmax_val` 的重复。
##   · 行对应：新增 `row_ids` 与 `assume_positional_pairing`，不再只比尺寸。
##   · 有限性/范围：NaN/Inf 记 NA 并计数；[0,1] 之外的值直接拒绝。
##
## 术语：不称"纯选择偏倚"，改称【参照校正的选择乐观度 Orel】；
##       Oraw = Cv − Ct 是另一指标，严格区分。
## ============================================================

#' 模块 4：选择乐观度与平均边际进入效应
#'
#' 给定"每个签名 × 每个候选"的验证集/独立测试集 C-index 矩阵，计算：
#' \itemize{
#'   \item \strong{池级选择乐观度}：在全池中按验证集选最优，其
#'         `Orel = Gv − Gt`（\code{Gv} 为验证增益、\code{Gt} 为独立测试增益，
#'         两者均相对参照模型）。
#'   \item \strong{平均边际进入效应}：把候选 \code{t} 加入一个大小为 \code{k−1}
#'         的**随机背景子池**，偏倚的变化量；对 \code{n_background} 个随机背景
#'         求平均，并给出背景抽样的 Monte Carlo 标准误，
#'         以及背景组成的分散程度（\code{sd_bg}，**不是**患者抽样的置信区间）。
#' }
#'
#' @section 使用限制（重要）:
#' 单个候选的进入效应排名在**跨数据划分时只有中等可复现性**，
#' **不足以在每个队列里可靠地点名"头号元凶"**。因此本模块的主输出是
#' 池级乐观度；进入效应应作为辅助分析，并同时报告其跨划分的变化。
#'
#' @section 行对应:
#' 只核对矩阵尺寸与列名**不足以**确认两份矩阵的行是同一批任务 —— 行顺序被调换
#' 但尺寸相同会静默产生错误配对。因此：提供 `row_ids`（唯一任务 ID）时按 ID 对齐；
#' 未提供时，仅在两份矩阵**行名完全一致且唯一**的前提下接受位置配对，
#' 否则必须显式设 `assume_positional_pairing = TRUE` 由调用者承担。
#'
#' @param val_ci,test_ci 数值矩阵，行 = 签名，列 = 候选（含参照）
#' @param reference 参照候选列名，默认第一列
#' @param n_background 每个 (目标, 规模) 的随机背景子池数，默认 50
#' @param k_set 评估的池规模（加入目标后的总数）；默认覆盖 2…K
#' @param seed 随机种子（决定背景抽样与并列优先级）
#' @param row_ids 可选，唯一的任务行 ID；用于显式校验/对齐 val_ci 与 test_ci 的行
#' @param assume_positional_pairing 未提供 row_ids 且行名不可核验时，是否声明接受位置配对
#' @return 含 pool_bias / dose / entry / entry_summary 的列表
#' @examples
#' f <- system.file("extdata", "example_ci.csv", package = "sigAudit")
#' d <- read.csv(f)
#' V  <- as.matrix(xtabs(val_ci  ~ geneset + model, data = d))
#' Tt <- as.matrix(xtabs(test_ci ~ geneset + model, data = d))
#' r  <- audit_selection(V, Tt, reference = "Lasso", n_background = 10)
#' r$pool_bias
#' head(r$entry)
#' @export
audit_selection <- function(val_ci, test_ci, reference = NULL,
                            n_background = 50, k_set = NULL, seed = 1,
                            row_ids = NULL, assume_positional_pairing = FALSE) {
  V <- as.matrix(val_ci); Tt <- as.matrix(test_ci)
  if (!identical(dim(V), dim(Tt))) stop("val_ci and test_ci have different dimensions")
  if (is.null(colnames(V))) stop("val_ci must have column names (candidate names)")
  if (!identical(colnames(V), colnames(Tt))) stop("val_ci and test_ci have different column names")

  ## 行对应 + 有限性/范围校验
  al <- .align_rows(V, Tt, row_ids, assume_positional_pairing)
  V <- al$V; Tt <- al$Tt
  cv <- .validate_ci(V, "val_ci"); V <- cv$m
  ct <- .validate_ci(Tt, "test_ci"); Tt <- ct$m
  n_nonfinite <- cv$n_nonfinite + ct$n_nonfinite

  ok <- stats::complete.cases(V) & stats::complete.cases(Tt)
  n_drop_incomplete <- sum(!ok)
  V <- V[ok, , drop = FALSE]; Tt <- Tt[ok, , drop = FALSE]
  if (nrow(V) < 10) stop("fewer than 10 complete signatures; cannot estimate")
  K <- ncol(V)
  if (K < 3) stop("fewer than 3 candidates")
  if (is.null(reference)) reference <- colnames(V)[1]
  if (!reference %in% colnames(V)) stop("reference model not among column names: ", reference)

  mods <- colnames(V)
  refv <- V[, reference]; reft <- Tt[, reference]
  # 与输入列顺序无关的固定优先级（共用于 audit_rules）
  prio <- .fixed_prio(mods, seed)
  pick <- function(v, cols) {
    m <- suppressWarnings(max(v[cols], na.rm = TRUE))
    if (!is.finite(m)) return(NA_character_)
    cand <- cols[is.finite(v[cols]) & v[cols] == m]
    cand[which.max(prio[cand])]
  }
  # 子池 S 的乐观度（逐签名 Orel 的均值）
  orel_S <- function(S) {
    sel <- vapply(seq_len(nrow(V)), function(i) pick(V[i, ], S), character(1))
    if (anyNA(sel)) return(NA_real_)
    gv <- V[cbind(seq_len(nrow(V)), match(sel, colnames(V)))] - refv
    gt <- Tt[cbind(seq_len(nrow(V)), match(sel, colnames(V)))] - reft
    mean(gv - gt, na.rm = TRUE)
  }

  pool_bias <- orel_S(mods)
  dose <- data.frame(k = 2:K,
                     candidate_added = mods[2:K],
                     orel = vapply(2:K, function(k) orel_S(mods[seq_len(k)]), numeric(1)))
  if (is.null(k_set)) k_set <- sort(unique(pmin(c(2, 3, 5, 8, 12, K), K)))
  k_set <- sort(unique(k_set[k_set >= 2 & k_set <= K]))
  set.seed(seed + 7)
  entry <- do.call(rbind, lapply(setdiff(mods, reference), function(tg) {
    bg <- sort(setdiff(mods, tg))   # 排序后抽样：使结果与输入列顺序无关
    do.call(rbind, lapply(k_set, function(k) {
      if (length(bg) < k - 1) return(NULL)
      if (k - 1 == length(bg)) {           # 全池：只有一个背景集合，不假装抽样
        d <- orel_S(c(bg, tg)) - orel_S(bg)
        return(data.frame(target = tg, k = k, entry = d, sd_bg = NA_real_,
                          mc_se = NA_real_, n_bg = 1L))
      }
      b <- vapply(seq_len(n_background), function(i) {
        S <- sample(bg, k - 1); orel_S(c(S, tg)) - orel_S(S) }, numeric(1))
      data.frame(target = tg, k = k, entry = mean(b), sd_bg = stats::sd(b),
                 mc_se = stats::sd(b) / sqrt(n_background), n_bg = n_background)
    }))
  }))
  agg <- stats::aggregate(entry ~ target + k, data = entry, FUN = function(x) mean(x, na.rm = TRUE))
  structure(list(module = "selection", reference = reference, n = nrow(V),
                 n_candidate = K, pool_bias = pool_bias, dose = dose,
                 entry = entry, entry_summary = agg,
                 coverage = c(n_signatures = nrow(V), n_dropped_incomplete = n_drop_incomplete,
                              n_nonfinite = n_nonfinite),
                 notes = c("pool_bias is the reproducible quantity",
                           "the dose/prefix curve follows the input column order; its composition depends on that order",
                           "per-candidate entry effects are marginal changes under random background pools",
                           "entry-effect rankings are only moderately reproducible across data splits")),
            class = "sigAudit_module")
}

#' 模块 5：决策规则核算
#'
#' 偏倚小不等于好——"专挑最差模型"的乐观度也是 0。评价一条规则必须同时看
#' **独立测试集上的性能**与**乐观度**。本模块把若干规则放在同一把尺子上：
#' \describe{
#'   \item{\code{none}}{不做选择，直接使用参照模型（核算基准）}
#'   \item{\code{argmax_val}}{在**全池（含参照）**中挑验证集最优}
#'   \item{\code{runner_up}}{挑验证集次优}
#'   \item{\code{argmax_nonreference}}{仅在非参照候选中挑（敏感性；回答的是
#'         "强制选一个非参照"的收益，**不是**通用的全池选优）}
#'   \item{\code{linear_only} / \code{tree_only}}{在指定子池中挑。
#'         \strong{两者不是对称的家族对照}：\code{linear_only} 是
#'         **全部非树候选、包含参照**（与 `scripts/20_analysis_corrected.R` 一致），
#'         \code{tree_only} 是全部树候选。只有当两侧候选规模与任务集一致时才可作家族比较；
#'         本包不主张家族效应。需要这两条规则时必须提供 \code{family}。}
#'   \item{\code{oracle_test}}{用测试集挑——**该候选集合内、该测试样本上的事后最大值**，
#'         不是总体上界，仅作描述性参照}
#' }
#'
#' @param val_ci,test_ci 数值矩阵，行 = 签名，列 = 候选
#' @param reference 参照候选列名，默认第一列
#' @param family 可选命名向量，给出每个候选的家族标签；提供时才会生成
#'   \code{linear_only}/\code{tree_only}
#' @param tree_label \code{family} 中代表树集成的取值，默认 \code{"tree"}
#' @param seed 并列优先级所用种子（与 \code{audit_selection} 同一函数）
#' @param row_ids,assume_positional_pairing 见 \code{\link{audit_selection}}
#' @return 含 table（每行一条规则）的 sigAudit_module
#' @examples
#' f <- system.file("extdata", "example_ci.csv", package = "sigAudit")
#' d <- read.csv(f)
#' V  <- as.matrix(xtabs(val_ci  ~ geneset + model, data = d))
#' Tt <- as.matrix(xtabs(test_ci ~ geneset + model, data = d))
#' fam <- ifelse(grepl("^(RSF|GBM)", colnames(V)), "tree", "linear")
#' names(fam) <- colnames(V)
#' audit_rules(V, Tt, reference = "Lasso", family = fam)$table
#' @export
audit_rules <- function(val_ci, test_ci, reference = NULL, family = NULL,
                        tree_label = "tree", seed = 1,
                        row_ids = NULL, assume_positional_pairing = FALSE) {
  V <- as.matrix(val_ci); Tt <- as.matrix(test_ci)
  if (!identical(dim(V), dim(Tt))) stop("val_ci and test_ci have different dimensions")
  if (is.null(colnames(V))) stop("val_ci must have column names (candidate names)")
  if (!identical(colnames(V), colnames(Tt))) stop("val_ci and test_ci have different column names")

  al <- .align_rows(V, Tt, row_ids, assume_positional_pairing)
  V <- al$V; Tt <- al$Tt
  cv <- .validate_ci(V, "val_ci"); V <- cv$m
  ct <- .validate_ci(Tt, "test_ci"); Tt <- ct$m
  n_nonfinite <- cv$n_nonfinite + ct$n_nonfinite

  ok <- stats::complete.cases(V) & stats::complete.cases(Tt)
  n_drop_incomplete <- sum(!ok)
  V <- V[ok, , drop = FALSE]; Tt <- Tt[ok, , drop = FALSE]
  n <- nrow(V)
  if (n < 10) stop("fewer than 10 complete signatures")
  if (is.null(reference)) reference <- colnames(V)[1]
  if (!reference %in% colnames(V)) stop("reference model not among column names: ", reference)
  mods <- colnames(V)
  if (length(mods) < 3) stop("fewer than 3 candidates")

  # 与输入列顺序无关的固定优先级（与 audit_selection 共用同一函数）
  prio <- .fixed_prio(mods, seed)
  pick1 <- function(cols, i) {
    v <- V[i, cols]; m <- max(v, na.rm = TRUE)
    cc <- cols[is.finite(v) & v == m]; cc[which.max(prio[cc])] }
  pickT <- function(cols, i) {
    v <- Tt[i, cols]; m <- max(v, na.rm = TRUE)
    cc <- cols[is.finite(v) & v == m]; cc[which.max(prio[cc])] }
  rank2 <- function(cols, i) {
    o <- order(-V[i, cols], -prio[cols], na.last = NA); cols[o][min(2, length(o))] }

  # 家族规则必须显式提供 family；否则 linear_only 会退化成 argmax_val 的重复
  if (is.null(family)) {
    lin <- NULL; tre <- NULL
  } else {
    if (is.null(names(family))) stop("family must be a named vector (names = candidate names)")
    tre <- mods[mods %in% names(family)[family == tree_label]]
    lin <- setdiff(mods, tre)   # 含参照，与 scripts/20 的 linear_only 定义一致
  }
  RULES <- list(
    none = list(sel = rep(reference, n)),
    argmax_val = list(sel = vapply(seq_len(n), function(i) pick1(mods, i), character(1))),
    runner_up  = list(sel = vapply(seq_len(n), function(i) rank2(mods, i), character(1))),
    argmax_nonreference = list(sel = vapply(seq_len(n), function(i)
      pick1(setdiff(mods, reference), i), character(1))),
    oracle_test = list(sel = vapply(seq_len(n), function(i) pickT(mods, i), character(1))))
  if (length(lin) >= 2) RULES$linear_only <- list(sel = vapply(seq_len(n), function(i) pick1(lin, i), character(1)))
  if (length(tre) >= 2) RULES$tree_only <- list(sel = vapply(seq_len(n), function(i) pick1(tre, i), character(1)))

  refv <- V[, reference]; reft <- Tt[, reference]
  out <- do.call(rbind, lapply(names(RULES), function(rn) {
    s <- RULES[[rn]]$sel; j <- match(s, mods)
    cvv <- V[cbind(seq_len(n), j)]; ctt <- Tt[cbind(seq_len(n), j)]
    gv <- cvv - refv; gt <- ctt - reft
    data.frame(rule = rn, n = n,
               reported = mean(cvv, na.rm = TRUE),      # Cv：验证集上的性能
               test     = mean(ctt, na.rm = TRUE),      # Ct：独立保留测试集上的性能估计
               Gv = mean(gv, na.rm = TRUE), Gt = mean(gt, na.rm = TRUE),
               Orel = mean(gv - gt, na.rm = TRUE),     # 参照校正的选择乐观度
               Oraw = mean(cvv - ctt, na.rm = TRUE))   # 原始落差（另一指标）
  }))
  structure(list(module = "rules", reference = reference, n = n, table = out,
                 coverage = c(n_signatures = n, n_dropped_incomplete = n_drop_incomplete,
                              n_nonfinite = n_nonfinite),
                 notes = c("Orel is reference-adjusted selection optimism, not 'pure selection bias'",
                           "Ct is a held-out test estimate, not a known truth",
                           "linear_only includes the reference and is NOT a symmetric family contrast",
                           "oracle_test is a post-hoc maximum within this candidate set, not a general upper bound")),
            class = "sigAudit_module")
}
