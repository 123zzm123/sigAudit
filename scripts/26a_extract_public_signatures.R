#!/usr/bin/env Rscript
# ============================================================
# 26a — 从 genefu 2.44.0 提取「已发表乳腺癌预后签名」的基因列表
#
# 第二轮评审建议的现实任务模块：用**可追溯的已发表签名**替代/补充 MSigDB 基因集合，
# 以缩短"基因集合"与"真实发表签名"之间的距离。
#
# 溯源原则（评审要求）：
#   · 只纳入**能恢复出干净基因符号**的签名；
#   · 明确区分"探针数"与"唯一基因数"（例如 70 基因签名实为 70 个探针对应 51 个唯一符号）；
#   · 无法恢复符号的（行名是 Affy/Agilent 探针 ID）一律排除并记录原因；
#   · 记录来源对象与原始论文，不凭记忆写基因名。
#
# 依赖：需要 genefu 源码包的数据对象（脚本会说明获取方式；若本地无则报错退出）。
# 输出：results/manifest/public_signatures_genes.tsv（长表）
#       results/manifest/public_signatures_meta.tsv （签名层面 + 溯源信息）
# ============================================================
options(stringsAsFactors = FALSE)
DATA_DIR <- Sys.getenv("SIGAUDIT_GENEFU_DATA",
                       file.path(Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit"), "third_party/genefu/data"))
OUT <- file.path(Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit"), "results/manifest"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(DATA_DIR))
  stop("找不到 genefu 数据目录: ", DATA_DIR, "\n",
       "获取方式：curl -sL https://bioconductor.org/packages/release/bioc/src/contrib/",
       "genefu_2.44.0.tar.gz | tar xz  （只需解包取 data/*.rda，不必安装整个包）")

load_obj <- function(f) { e <- new.env(); load(file.path(DATA_DIR, paste0(f, ".rda")), envir = e)
  get(ls(e)[1], envir = e) }

long <- list(); meta <- list()
add <- function(sig, genes, n_probe, hascoef, src, ref, note) {
  genes <- as.character(genes)
  raw_n <- length(genes)                       # 清洗前的条目数（探针/条目）
  genes <- sort(unique(trimws(genes[!is.na(genes) & nzchar(trimws(genes))])))
  if (length(genes) < 5) return(invisible(NULL))
  # note 一律由【清洗后的向量】生成，避免与 n_gene_unique 自相矛盾
  # （旧版把唯一符号数写进 note 时用的是未清洗的列，导致出现 "128 -> 100" 而列值是 99）
  auto <- sprintf("%d 个探针/条目 -> %d 个唯一基因符号", raw_n, length(genes))
  note <- if (missing(note) || is.null(note) || !nzchar(note)) auto else paste0(auto, "；", note)
  long[[sig]] <<- data.frame(signature = sig, gene = genes, stringsAsFactors = FALSE)
  meta[[sig]] <<- data.frame(signature = sig, n_gene_unique = length(genes),
    n_probe_or_entry = n_probe, has_coefficients = hascoef,
    source_object = src, reference = ref, note = note, stringsAsFactors = FALSE)
}

## ---------- 纳入：能恢复干净基因符号的签名 ----------
o <- load_obj("pam50.robust")
add("PAM50", rownames(o$centroids), nrow(o$centroids), FALSE, "pam50.robust",
    "Parker JS et al. J Clin Oncol 2009;27:1160-1167")

o <- load_obj("sig.ggi")
add("GGI", o$HUGO.gene.symbol, nrow(o), FALSE, "sig.ggi",
    "Sotiriou C et al. J Natl Cancer Inst 2006;98:262-272",
    NULL)

o <- load_obj("sig.gene70")
add("GENE70_MammaPrint", o$HUGO.gene.symbol, nrow(o), TRUE, "sig.gene70",
    "van 't Veer LJ et al. Nature 2002;415:530-536",
    NULL)

o <- load_obj("sig.gene76")
add("GENE76_Rotterdam", o$NCBI.gene.symbol, nrow(o), TRUE, "sig.gene76",
    "Wang Y et al. Lancet 2005;365:671-679",
    NULL)

o <- load_obj("sig.oncotypedx")
add("OncotypeDX_21", o$symbol, nrow(o), TRUE, "sig.oncotypedx",
    "Paik S et al. N Engl J Med 2004;351:2817-2826", "21 个基因")

o <- load_obj("sig.endoPredict")
add("EndoPredict_11", o$symbol, nrow(o), TRUE, "sig.endoPredict",
    "Filipits M et al. Clin Cancer Res 2011;17:6012-6020", "11 个基因")

## ---------- 排除：记录原因，不纳入 ----------
excluded <- data.frame(
  object = c("ssp2003.robust","ssp2006.robust","mod1","mod2","sig.tamr13","claudinLowData","sig.pik3cags"),
  reason = c(
    "行名是 Agilent 探针 ID（AGI_HUM1_OLIGO_...），无法可靠恢复基因符号",
    "同上（306 个探针 ID）",
    "是亚型分类模块的基因集合（并集 2733），不是单个预后签名",
    "同上（1666）",
    "对象为 11 个 pclust 矩阵，首矩阵仅 7 行，定义不清",
    "xfnames 为 Entrez ID（807 个），属分类器特征而非预后签名",
    "278 个 Entrez ID，属通路基因列表而非预后签名"),
  stringsAsFactors = FALSE)

L <- do.call(rbind, long); M <- do.call(rbind, meta)
M <- M[order(M$n_gene_unique), ]
write.table(L, file.path(OUT, "public_signatures_genes.tsv"), sep="\t", quote=FALSE, row.names=FALSE, eol="\n")
write.table(M, file.path(OUT, "public_signatures_meta.tsv"),  sep="\t", quote=FALSE, row.names=FALSE, eol="\n")
write.table(excluded, file.path(OUT, "public_signatures_excluded.tsv"), sep="\t", quote=FALSE, row.names=FALSE, eol="\n")

cat("纳入的已发表签名：\n")
for (i in seq_len(nrow(M))) cat(sprintf("  %-22s %3d 唯一基因 (探针/条目 %3d)  系数:%s\n",
  M$signature[i], M$n_gene_unique[i], M$n_probe_or_entry[i],
  ifelse(M$has_coefficients[i], "有", "无")))
cat(sprintf("\n共 %d 个签名；去重后基因 %d 个\n", nrow(M), length(unique(L$gene))))
cat("排除 ", nrow(excluded), " 个对象（原因见 public_signatures_excluded.tsv）\n")
