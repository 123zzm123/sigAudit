## ============================================================
## utils.R — 内部工具函数
## ============================================================

#' 固定方向的 C-index
#'
#' 方向恒定为「分数越高 = 风险越高」，使用
#' `survival::concordance(..., reverse = TRUE, timewt = "n")` **直接评价**。
#'
#' 为什么不再"自动识别方向"：旧写法
#' `summary(coxph(Surv ~ s))$concordance[1]` 会在**当前评价集**上重新估一个系数，
#' 方向由该评价集的结局决定 —— 在独立测试集上调用即构成泄漏。已实证：对 `s` 与
#' `-s` 该写法给出**完全相同**的值（如本项目的 MESO 样例 0.667077）。
#' 因此本函数**不做任何重拟合**，方向固定。
#'
#' @param min_n,min_ev 计算 C-index 所需的最少样本数/事件数。默认 15/5。
#'   **不要**调回 30/10：那会让小队列（如 MESO 的测试集只有 22 例、19 事件）
#'   静默返回 NA。
#' @keywords internal
#' @noRd
.c_index <- function(time, event, score, min_n = 15, min_ev = 5) {
  ok <- !is.na(score) & is.finite(score) & !is.na(time) & !is.na(event) & time > 0
  if (sum(ok) < min_n || sum(event[ok]) < min_ev) return(NA_real_)
  # 直接评价：不在评价集上重新拟合任何模型，方向固定为"高值 = 高风险"。
  # 旧写法 summary(coxph(Surv ~ s))$concordance[1] 会在【当前评价集】上重新估一个
  # 系数，方向由该评价集的结局决定 —— 在测试集上调用即构成泄漏（已实证：对 s 与
  # -s 给出完全相同的值）。reverse=TRUE 才对应"高值=高风险"的 Harrell C。
  cd <- tryCatch(survival::concordance(survival::Surv(time[ok], event[ok]) ~ score[ok],
                                       reverse = TRUE, timewt = "n"),
                 error = function(e) NULL)
  if (is.null(cd)) return(NA_real_)
  as.numeric(cd$concordance)
}

#' 可比对生存对数（取 concordance 的 count 字段，不用自写判据）
#'
#' **已修正（第三轮审查）**：Harrell C 的有效可比对数分母是
#' `concordant + discordant + tied.x`。
#' 旧实现把所有以 `tied` 开头的 count 都加了进去，其中
#' `tied.y`（事件时间并列）与 `tied.xy`（分数与时间都并列）**不属于可比对**，
#' 会把分母算大（本项目 MESO 测试集因此报 224，正确值为 223）。
#' 注意：这只是**计数诊断**的错误，`c_index()` 本身不受影响。
#' @keywords internal
#' @noRd
.n_pairs <- function(time, event) {
  ok <- !is.na(time) & !is.na(event) & is.finite(time) & is.finite(event)
  if (sum(ok) < 2) return(0)
  cd <- tryCatch(survival::concordance(survival::Surv(time[ok], event[ok]) ~ seq_along(time[ok]),
                                       timewt = "n"), error = function(e) NULL)
  if (is.null(cd) || is.null(cd$count)) return(NA_real_)
  cts <- cd$count
  if (is.null(names(cts))) return(as.numeric(sum(cts, na.rm = TRUE)))
  tx <- if ("tied.x" %in% names(cts)) unname(cts["tied.x"]) else 0
  if (all(c("concordant", "discordant") %in% names(cts)))
    return(as.numeric(unname(cts["concordant"]) + unname(cts["discordant"]) + tx))
  as.numeric(sum(cts, na.rm = TRUE))
}

#' 与输入列顺序【无关】的固定并列优先级
#'
#' 在**候选名字排序后的规范集合**上生成随机优先级，再按名字映射。
#' 关键是候选身份由**名字**决定，不由列位置决定，因此交换输入列顺序不改变结果。
#' `audit_selection()` 与 `audit_rules()` 共用本函数（审查要求"所有接口共用同一函数
#' 并显式传入 seed"）。旧实现 `setNames(seq_along(mods), mods)` 把优先级绑定在
#' **输入列顺序**上，等价于"偏向排在后面的候选"。
#' @keywords internal
#' @noRd
.fixed_prio <- function(mods, seed = 1) {
  mods <- sort(unique(as.character(mods)))
  set.seed(seed)
  stats::setNames(stats::runif(length(mods)), mods)
}

#' 校验 C-index 矩阵：非有限值记为 NA 并计数，越界值直接拒绝
#' @keywords internal
#' @noRd
.validate_ci <- function(m, what = "matrix") {
  m <- as.matrix(m)
  if (!is.numeric(m)) stop(what, " must be numeric")
  nf <- !is.finite(m)
  n_nf <- sum(nf)
  if (n_nf) m[nf] <- NA_real_
  bad <- is.finite(m) & (m < 0 | m > 1)
  if (any(bad)) stop(sprintf(
    "%s contains %d C-index value(s) outside [0,1]; a C-index must lie in [0,1]",
    what, sum(bad)))
  list(m = m, n_nonfinite = n_nf)
}

#' 校验并（必要时）按任务 ID 对齐两份 C-index 矩阵的行
#'
#' 审查要求：只核对矩阵尺寸与列名不够 —— 行顺序被调换但尺寸相同会**静默错配**。
#' 规则：
#' \itemize{
#'   \item 提供 \code{row_ids}：必须唯一；若两份矩阵行名是同一集合，按 \code{row_ids} 重排对齐。
#'   \item 未提供：仅当两份矩阵**行名完全一致且唯一**时才接受（此时位置对应可核验）；
#'         否则必须显式声明 \code{assume_positional_pairing = TRUE}，由调用者承担位置配对。
#' }
#' @keywords internal
#' @noRd
.align_rows <- function(V, Tt, row_ids = NULL, assume_positional_pairing = FALSE, what = "input") {
  rnV <- rownames(V); rnT <- rownames(Tt)
  if (is.null(row_ids)) {
    if (!is.null(rnV) && !is.null(rnT) && identical(rnV, rnT)) {
      if (anyDuplicated(rnV)) stop("row names must be unique to verify row correspondence")
      return(list(V = V, Tt = Tt, row_ids = rnV))
    }
    if (!assume_positional_pairing) stop(
      "cannot verify that rows of val_ci and test_ci correspond: ",
      "supply unique `row_ids`, or set `assume_positional_pairing = TRUE` to ",
      "accept positional pairing explicitly.")
    return(list(V = V, Tt = Tt, row_ids = as.character(seq_len(nrow(V)))))
  }
  row_ids <- as.character(row_ids)
  if (length(row_ids) != nrow(V)) stop("length(row_ids) != nrow(val_ci)")
  if (anyDuplicated(row_ids)) stop("row_ids must be unique")
  if (!is.null(rnV) && !is.null(rnT)) {
    if (!setequal(row_ids, rnV) || !setequal(row_ids, rnT))
      stop("row_ids is not the same set as the row names of val_ci/test_ci; cannot align tasks")
    if (!identical(rnV, rnT)) {
      V <- V[row_ids, , drop = FALSE]; Tt <- Tt[row_ids, , drop = FALSE]
    }
  }
  list(V = V, Tt = Tt, row_ids = row_ids)
}

#' NA 安全的比较：NA / NULL / 长度不为 1 一律视为"不成立"，不抛错
#'
#' 为什么需要：`if (NA > 0.05)` 在 R 里会直接报
#' "missing value where TRUE/FALSE needed"。而 NA 在这里是很现实的输入
#' ——小队列的 C-index 就是 NA。陈旧写法会让整个分级函数崩掉。
#' @keywords internal
#' @noRd
.gt  <- function(a, b) isTRUE(is.numeric(a) && length(a) == 1L && !is.na(a) && a > b)
#' @keywords internal
#' @noRd
.lt  <- function(a, b) isTRUE(is.numeric(a) && length(a) == 1L && !is.na(a) && a < b)

#' 判断矩阵方向：返回 TRUE 表示已是【基因 × 样本】
#'
#' 不用"行数<列数"这种脆弱判据 —— 只有少量基因(如 30)而样本很多(如 200)时
#' 会被误判。改用【标签匹配】：哪一侧的标签能命中签名基因，哪一侧就是基因轴。
#' 无基因信息时退回"行数远多于列数"的保守判据。
#' @keywords internal
#' @noRd
.is_genes_by_samples <- function(mat, genes = NULL) {
  if (!is.null(genes) && length(genes)) {
    r <- length(intersect(genes, rownames(mat)))
    c <- length(intersect(genes, colnames(mat)))
    if (r != c) return(r > c)
  }
  nrow(mat) > ncol(mat)          # 保守兜底：表达矩阵基因数通常远多于样本数
}

#' 在【训练侧】估计基因层 z-score 所需的全部参数
#'
#' 为什么拆出来：`audit_transport()` 需要把**训练队列**估计的标准化原样应用到
#' 外部队列，而不是在外部队列上重新估计（后者违反"禁止在任何一步使用外部数据拟合"
#' 的归纳式预测承诺 —— 这是第三轮审查指出的实际缺陷）。
#' @keywords internal
#' @noRd
.fit_zscore <- function(train, genes = NULL, max_miss = 0.20) {
  tr <- as.matrix(train)
  if (!.is_genes_by_samples(tr, genes)) tr <- t(tr)
  n0 <- nrow(tr)
  miss <- rowSums(is.na(tr)) / ncol(tr)
  keep_miss <- miss <= max_miss
  tr <- tr[keep_miss, , drop = FALSE]
  med <- apply(tr, 1, function(v) stats::median(v, na.rm = TRUE))
  sdt <- apply(tr, 1, function(v) stats::sd(v, na.rm = TRUE))
  keep_sd <- is.finite(sdt) & sdt > 0
  tr <- tr[keep_sd, , drop = FALSE]
  list(genes = rownames(tr), med = med[keep_sd],
       mu = rowMeans(tr), sd = apply(tr, 1, function(v) stats::sd(v)),
       n_gene_in = n0, n_drop_missing = sum(!keep_miss), n_drop_zerovar = sum(!keep_sd))
}

#' 用已拟合的（训练侧）参数把任意矩阵标准化；基因取交集
#' @keywords internal
#' @noRd
.apply_zscore <- function(mat, par, genes = NULL) {
  m <- as.matrix(mat)
  if (!.is_genes_by_samples(m, genes)) m <- t(m)
  g <- intersect(par$genes, rownames(m))
  if (!length(g)) stop("no shared genes between the data and the fitted z-score parameters")
  m <- m[g, , drop = FALSE]
  n_imp <- sum(is.na(m))
  if (n_imp > 0) {
    for (j in seq_len(nrow(m))) {
      v <- m[j, ]; k <- is.na(v)
      if (any(k)) { v[k] <- par$med[g[j]]; m[j, ] <- v }
    }
  }
  z <- (m - par$mu[g]) / par$sd[g]
  z[!is.finite(z)] <- 0        # 已剔除零方差基因，此处仅为数值残渣兜底
  attr(z, "n_imputed") <- n_imp
  z
}

#' 把表达矩阵统一成【基因 × 样本】并做基因层 z-score（训练侧拟合 + 应用）
#' @keywords internal
#' @noRd
.zscore_genes <- function(mat, genes = NULL, train = NULL, max_miss = 0.20) {
  # 标准化参数（均值/标准差/中位数/零方差规则）**只在训练集上估计**，
  # 再原样应用到验证与测试。旧实现用整队列估计，构成信息泄漏。
  mat <- as.matrix(mat)
  if (!.is_genes_by_samples(mat, genes)) mat <- t(mat)
  if (is.null(train)) train <- seq_len(ncol(mat))   # 未指定时退回全体（不推荐，仅为兼容）
  par <- .fit_zscore(mat[, train, drop = FALSE], genes, max_miss)
  z <- .apply_zscore(mat, par, genes)
  attr(z, "n_gene_in")     <- par$n_gene_in
  attr(z, "n_drop_missing") <- par$n_drop_missing
  attr(z, "n_drop_zerovar") <- par$n_drop_zerovar
  attr(z, "train_idx")      <- train
  z
}

#' 统一的样本 ID 规范：截断到 patient 级（12 字符）
#'
#' 表达矩阵常按 sampleId(TCGA-XX-XXXX-01) 索引，生存表按
#' patientId(TCGA-XX-XXXX)。不统一会导致交集为 0 —— 这是本项目
#' 实际踩过三次的坑（BRCA / METABRIC / LGG），故固化为标准步骤。
#' @keywords internal
#' @noRd
.harmonise_ids <- function(x) substr(as.character(x), 1, 12)

#' 按生存数据结构对齐表达矩阵与生存表
#'
#' @param event_levels 可选的**显式终点字典**：命名向量，名字是原始事件列的取值字符串、
#'   值是 0/1（1 = 事件）。**强烈建议提供**。未提供时退回启发式正则
#'   （`^1|DECEASED|Dead` -> 1，`^0|LIVING|Alive` -> 0），该启发式对
#'   非标准编码（如 `1:PROGRESSION`/`0:CENSORED` 之外的组合）可能静默误判，
#'   因此不要把启发式结果当作终点语义的证据。
#' @keywords internal
#' @noRd
.align_cohort <- function(expr, surv, time_col, event_col, genes = NULL, event_levels = NULL) {
  expr <- as.matrix(expr)
  if (!.is_genes_by_samples(expr, genes)) expr <- t(expr)   # -> 基因 × 样本
  colnames(expr) <- .harmonise_ids(colnames(expr))
  expr <- expr[, !duplicated(colnames(expr)), drop = FALSE]

  sv <- as.data.frame(surv)
  t <- suppressWarnings(as.numeric(sv[[time_col]]))
  raw <- as.character(sv[[event_col]])
  if (!is.null(event_levels)) {
    if (is.null(names(event_levels))) stop("event_levels must be a named vector (names = raw values)")
    e <- unname(event_levels[raw])
    e <- ifelse(is.na(e), NA_integer_, ifelse(e == 1, 1L, 0L))
  } else {
    e <- ifelse(grepl("^1|DECEASED|Dead", raw), 1L,
                ifelse(grepl("^0|LIVING|Alive", raw), 0L, NA_integer_))
  }
  keep <- !is.na(t) & !is.na(e) & t > 0
  sv <- data.frame(id = .harmonise_ids(rownames(sv)), t = t, e = e)[keep, ]
  sv <- sv[!duplicated(sv$id), ]

  com <- intersect(colnames(expr), sv$id)
  if (!length(com)) stop("no shared samples between expression and survival data; check sample ID formats")
  list(expr = expr[, com, drop = FALSE],
       surv = sv[match(com, sv$id), ])
}

#' 生成训练/测试划分（按事件状态分层）
#' @keywords internal
#' @noRd
.stratified_split <- function(event, prop = 0.7, seed = 20261014) {
  set.seed(seed)
  ev <- which(event == 1); ce <- which(event == 0)
  tr <- sort(c(sample(ev, floor(prop * length(ev))),
               sample(ce, floor(prop * length(ce)))))
  list(train = tr, test = setdiff(seq_along(event), tr))
}
