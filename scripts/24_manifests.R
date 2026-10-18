#!/usr/bin/env Rscript
# ============================================================
# 24 — 生成评审要求的两份清单
#   results/manifest/sample_manifest.tsv  样本清单（含每个种子的 train/val/test 归属）
#   results/manifest/geneset_manifest.tsv 基因集清单（抽样得到的 100 个基因集及各队列交集）
#   results/manifest/candidate_diversity.tsv 候选多样性（每个任务上实际有多少个互不相同的预测）
#
# 说明：
#  · 划分归属由【与 19 完全相同的种子逻辑】重新推导，因此对尚未跑完的种子也能给出；
#    跑完后可用 19 产物的 train/val/test ID 复核（_verify_19.R 已做此校验）。
#  · 样本类型由 TCGA 条形码第 14–15 位判定（01 = 原发实体瘤）。
#  · 候选多样性用于回答评审问题："18 个名字对应多少个真正不同的预测"。
# ============================================================
suppressMessages(library(jsonlite))
options(stringsAsFactors = FALSE)
`%||%` <- function(a, b) if (is.null(a)) b else a
PROJ <- Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
MAN  <- file.path(PROJ, "results/manifest"); dir.create(MAN, showWarnings = FALSE, recursive = TRUE)

ENDPOINTS <- list(
  LGG  = list(e=file.path(PROJ, "data/cohorts/LGG_expr.tsv"),
              f=file.path(PROJ, "data/cohorts/LGG_surv.tsv"),
              t="PFS_MONTHS", s="PFS_STATUS", endpoint="PFS",
              event=c("1:PROGRESSION"), censor=c("0:CENSORED")),
  PAAD = list(t="PFS_MONTHS", s="PFS_STATUS", endpoint="PFS",
              event=c("1:PROGRESSION"), censor=c("0:CENSORED")),
  BRCA = list(t="OS_MONTHS", s="OS_STATUS", endpoint="OS",
              event=c("1:DECEASED"), censor=c("0:LIVING")),
  COAD = list(t="OS_MONTHS", s="OS_STATUS", endpoint="OS",
              event=c("1:DECEASED"), censor=c("0:LIVING")),
  KIRC = list(t="OS_MONTHS", s="OS_STATUS", endpoint="OS",
              event=c("1:DECEASED"), censor=c("0:LIVING")),
  LUAD = list(t="OS_MONTHS", s="OS_STATUS", endpoint="OS",
              event=c("1:DECEASED"), censor=c("0:LIVING")),
  MESO = list(t="OS_MONTHS", s="OS_STATUS", endpoint="OS",
              event=c("1:DECEASED"), censor=c("0:LIVING")),
  STAD = list(t="OS_MONTHS", s="OS_STATUS", endpoint="OS",
              event=c("1:DECEASED"), censor=c("0:LIVING"))
)
STUDY <- function(c) if (c == "LGG") "local:LGG_broad" else
  sprintf("%s_tcga_pan_can_atlas_2018", tolower(if (c == "COAD") "coadread" else c))
SEEDS <- c(20261014, 20261016, 20261017, 20261018, 20261019)
STYPE <- c("01"="primary_solid_tumor","02"="recurrent_tumor","06"="metastatic",
           "11"="solid_tissue_normal","03"="primary_blood_cancer","07"="additional_metastatic",
           "12"="blood_derived_normal","14"="bone_marrow_normal")

## ---------- 1. 样本清单 ----------
rows <- list()
for (cn in names(ENDPOINTS)) {
  E <- ENDPOINTS[[cn]]
  cf <- E[["e"]] %||% sprintf("%s/data/cohorts/%s_expr.tsv", PROJ, cn)
  sf <- E[["f"]] %||% sprintf("%s/data/cohorts/%s_surv.tsv", PROJ, cn)
  expr <- as.matrix(read.table(cf, header=TRUE, sep="\t", row.names=1, check.names=FALSE))
  if (nrow(expr) < ncol(expr)) expr <- t(expr)
  sv <- read.table(sf, header=TRUE, sep="\t", row.names=1, check.names=FALSE)
  raw <- as.character(sv[[E[["s"]]]])
  sv$t <- suppressWarnings(as.numeric(sv[[E[["t"]]]]))
  sv$e <- ifelse(raw %in% E[["event"]], 1L, ifelse(raw %in% E[["censor"]], 0L, NA_integer_))
  sv <- sv[!is.na(sv$t) & !is.na(sv$e) & sv$t > 0, ]
  ids_full <- colnames(expr)
  ids12 <- substr(ids_full, 1, 12)
  keep <- !duplicated(ids12)
  inclusion <- ifelse(keep, "included", "excluded_duplicate_patient")
  ids_full <- ids_full[keep]; ids12 <- ids12[keep]
  in_sv <- ids12 %in% rownames(sv)
  inclusion[!in_sv] <- "excluded_no_survival"
  base <- data.frame(cohort=cn, study_id=STUDY(cn), endpoint=E[["endpoint"]],
                     sample_id=ids_full, patient_id=ids12,
                     sample_type_code=substr(ids_full, 14, 15),
                     stringsAsFactors = FALSE)
  base$sample_type <- STYPE[base$sample_type_code]
  base$inclusion <- inclusion
  # 与 19 相同的划分逻辑（用进入分析的样本）
  use <- base$inclusion == "included" & in_sv
  ord <- which(use)                       # 保持与 19 相同的样本顺序（expr 列序 ∩ sv 行序）
  svsub <- sv[base$patient_id[ord], ]
  for (s in SEEDS) {
    set.seed(s)
    ev <- which(svsub$e == 1); ce <- which(svsub$e == 0)
    tr <- sort(c(sample(ev, floor(.5*length(ev))), sample(ce, floor(.5*length(ce)))))
    rest <- setdiff(seq_along(svsub$e), tr)
    re <- rest[svsub$e[rest] == 1]; rc <- rest[svsub$e[rest] == 0]
    va <- sort(c(sample(re, floor(.5*length(re))), sample(rc, floor(.5*length(rc)))))
    te <- setdiff(rest, va)
    assign <- rep(NA_character_, nrow(base)); assign[ord[tr]] <- "train"
    assign[ord[va]] <- "val"; assign[ord[te]] <- "test"
    base[[paste0("split_", s)]] <- assign
  }
  base$event <- NA_integer_; base$event[ord] <- svsub$e
  base$time_months <- NA_real_; base$time_months[ord] <- svsub$t
  rows[[cn]] <- base
}
SM <- do.call(rbind, rows)
write.table(SM, file.path(MAN, "sample_manifest.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## ---------- 2. 基因集清单 ----------
gs <- fromJSON(file.path(PROJ,"data/signatures/msigdb_genesets.json"))
gs <- gs[sapply(gs,length)>=10 & sapply(gs,length)<=500]
set.seed(20261015); pick <- sample(names(gs), 100)
gm <- data.frame(geneset = pick, n_gene_source = sapply(gs[pick], length), stringsAsFactors = FALSE)
for (cn in names(ENDPOINTS)) {
  E <- ENDPOINTS[[cn]]
  cf <- E[["e"]] %||% sprintf("%s/data/cohorts/%s_expr.tsv", PROJ, cn)
  expr <- as.matrix(read.table(cf, header=TRUE, sep="\t", row.names=1, check.names=FALSE))
  if (nrow(expr) < ncol(expr)) expr <- t(expr)
  gm[[paste0("n_gene_in_", cn)]] <- sapply(pick, function(g) length(intersect(gs[[g]], rownames(expr))))
}
write.table(gm, file.path(MAN, "geneset_manifest.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## ---------- 3. 候选多样性（**只读最终输入清单列出的产物**）----------
## 第三轮审查：旧实现扫描 results/main/ 目录，结果只覆盖 198 个任务（当时只有
## MESO/PAAD 的 20261014），却被当成全实验统计使用。改为只读
## results/manifest/main_input_files.txt，保证覆盖恰好是最终 40 组。
## 同时区分"精确相等 / 6 位小数相等 / 仅排序相同"三档，正文需说明标准。
mi <- file.path(MAN, "main_input_files.txt")
if (!file.exists(mi)) stop("缺少 results/manifest/main_input_files.txt（请先跑 scripts/23_manifest.R）")
fs <- readLines(mi, warn = FALSE); fs <- fs[nzchar(fs)]
div <- list()
for (f in fs) {
  R <- tryCatch(readRDS(f), error=function(e) NULL)
  if (is.null(R) || !length(R$genesets)) next
  for (z in R$genesets) {
    P <- as.matrix(z$val_pred)
    cn <- colnames(P)
    k_exact <- apply(P, 2, function(v) paste(format(v, digits = 17), collapse = ","))
    k_r6    <- apply(P, 2, function(v) paste(sprintf("%.6f", v), collapse = ","))
    k_rank  <- apply(P, 2, function(v) paste(rank(v, ties.method = "min"), collapse = ","))
    # 精确重复的候选分组（用于判断 top-k 集成是否把同一预测重复加权）
    grp <- split(cn, k_exact)
    grp <- grp[vapply(grp, length, integer(1)) >= 2]
    dup_groups <- paste(vapply(grp, function(g) paste(g, collapse = ","), character(1)),
                        collapse = ";")
    div[[length(div)+1]] <- data.frame(cohort=R$cohort, seed=R$seed, geneset=z$name,
      n_models        = ncol(P),
      n_distinct_exact = length(unique(k_exact)),
      n_distinct_r6    = length(unique(k_r6)),
      n_distinct_rank  = length(unique(k_rank)),
      n_constant       = sum(apply(P, 2, function(v) length(unique(v)) == 1)),
      dup_groups       = dup_groups, stringsAsFactors = FALSE)
  }
}
if (length(div)) {
  D <- do.call(rbind, div)

  ## ---------- 4. 重复候选对集成的影响（只用验证侧信息判定，不看测试性能）----------
  rl <- file.path(PROJ, "results/tables/20_rule_members_long.tsv")
  if (file.exists(rl)) {
    RL <- read.table(rl, sep="\t", header=TRUE, quote="", comment.char="")
    t5 <- RL[RL$rule == "top5_ensemble", c("cohort","seed","geneset","members")]
    kk <- D[, c("cohort","seed","geneset","n_models","n_distinct_exact","n_distinct_r6","dup_groups")]
    t5 <- merge(t5, kk, by = c("cohort","seed","geneset"), all.x = TRUE)
    t5$top5_n        <- vapply(strsplit(t5$members, "|", fixed = TRUE), length, integer(1))
    t5$pool_has_dup  <- !is.na(t5$dup_groups) & nzchar(t5$dup_groups)
    # top5 的成员里**是否真的同时含有一对精确重复的候选**（这才叫被重复加权）
    dup_in_top5 <- function(members, dg) {
      if (is.na(dg) || !nzchar(dg)) return(FALSE)
      ms <- strsplit(members, "|", fixed = TRUE)[[1]]
      any(vapply(strsplit(dg, ";", fixed = TRUE)[[1]], function(g)
        sum(strsplit(g, ",", fixed = TRUE)[[1]] %in% ms) >= 2, logical(1)))
    }
    t5$top5_has_dup <- mapply(dup_in_top5, t5$members, t5$dup_groups)
    D <- merge(D, t5[, c("cohort","seed","geneset","pool_has_dup","top5_has_dup")],
               by = c("cohort","seed","geneset"), all.x = TRUE)
    write.table(t5, file.path(MAN, "candidate_top5_dup.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
    cat(sprintf("\n候选多样性：%d 个已拟合任务（分析集为其中 3919 个）\n", nrow(D)))
    cat(sprintf("  名义候选 %.1f 个；精确互不相同 %.1f；6 位小数互不相同 %.1f；排序互不相同 %.1f；常数 %.1f\n",
                mean(D$n_models), mean(D$n_distinct_exact), mean(D$n_distinct_r6),
                mean(D$n_distinct_rank), mean(D$n_constant)))
    cat(sprintf("  存在精确重复候选的任务占比: %.1f%%（6 位小数口径 %.1f%%）\n",
                100*mean(D$n_distinct_exact < D$n_models),
                100*mean(D$n_distinct_r6    < D$n_models)))
    cat(sprintf("  其中 top5 集成成员里【确实含有一对精确重复候选】的任务: %d / %d（%.1f%%）\n",
                sum(t5$top5_has_dup, na.rm = TRUE), nrow(t5),
                100*mean(t5$top5_has_dup, na.rm = TRUE)))
  } else {
    cat("\n[注意] 未找到 20_rule_members_long.tsv，跳过重复候选对集成的影响统计\n")
  }
  write.table(D, file.path(MAN, "candidate_diversity.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
}

cat("\n清单生成完毕：\n")
cat("  results/manifest/sample_manifest.tsv   ", nrow(SM), "行（含 5 个种子的划分归属）\n")
cat("  results/manifest/geneset_manifest.tsv ", nrow(gm), "个基因集\n")
if (length(div)) cat("  results/manifest/candidate_diversity.tsv ", nrow(D), "个任务\n")
cat("\n样本类型分布（全部队列）:\n"); print(table(SM$sample_type, useNA="ifany"))
