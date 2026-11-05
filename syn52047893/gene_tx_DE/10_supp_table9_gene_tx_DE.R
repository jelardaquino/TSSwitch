###############################
# Supplementary Table 9: every gene in Supplementary Figure 5 (long-read, 6 AD vs 6 CT), one row per gene.
#   Rows: genes that are a gene DEG, have >= 1 DE transcript, or are one of the 86 DIU genes (DESeq2 ~ condition).
# DE = FDR < 0.05 and |apeglm-shrunken log2FC| > 1, as in re-analysis_PMI_RIN/04b_venn_JA_DEG.R.
# Writes a TSV; 10b_supp_table9_xlsx.py turns it into the workbook.
###############################

S <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing"
PY <- file.path(S, "syn52065646_LongRead/pychopperfq")
RE <- file.path(PY, "re-analysis_PMI_RIN")
OUT <- file.path(S, "ERP161086/emtab14167_replication/tables")

DIU_JA <- file.path(PY, "minimap2/bam_mapq10/all_bambu_Rscript/all_bambu_python_analysis",
                    "IsoformSwitchAnalyzeR/ISO_results_v3/86_isoformswitch_JA.csv")
DIU_ANNOT <- file.path(OUT, "gene_level_86_DIU_genes.tsv")  # promoter-switching flag per DIU gene
MODELS <- "noCov"

strip_ver <- function(x) sub("\\.[0-9]+$", "", x)
rd <- function(f) read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
yn <- function(x) ifelse(is.na(x), "not tested", ifelse(x, "yes", "no"))
is_de <- function(d) !is.na(d$padj) & d$padj < 0.05 & abs(d$log2FC_shrunk) > 1

###############################
# DIU genes and switching isoforms
###############################

iso <- read.csv(DIU_JA, stringsAsFactors = FALSE)
iso$gene_id <- strip_ver(iso$gene_id)
iso$isoform_id <- strip_ver(iso$isoform_id)
diu <- unique(iso$gene_id)
stopifnot(length(diu) == 86)
annot <- rd(DIU_ANNOT)
annot$gene_id <- strip_ver(annot$gene_id)
sw_ids <- tapply(sprintf("%s (%+.2f)", iso$isoform_id, iso$dIF), iso$gene_id, paste, collapse = "; ")

###############################
# DE results per design
###############################

res <- lapply(setNames(MODELS, MODELS), function(m) {
  g <- rd(file.path(RE, "gene_DEG", sprintf("gene_DESeq2_%s_results.tsv", m)))
  t <- rd(file.path(RE, "tx_DEG", sprintf("tx_DESeq2_%s_results.tsv", m)))
  g$gene_id <- strip_ver(g$gene_id); t$gene_id <- strip_ver(t$gene_id); t$id <- strip_ver(t$id)
  g <- g[!is.na(g$pvalue), ]; t <- t[!is.na(t$pvalue) & !is.na(t$gene_id) & t$gene_id != "", ]
  g$de <- is_de(g); t$de <- is_de(t)
  list(g = g, t = t, deg = g$gene_id[g$de], txdeg = unique(t$gene_id[t$de]))
})

ids <- unique(c(diu, unlist(lapply(res, function(r) c(r$deg, r$txdeg)))))

# Gene symbol: DIU file first, then the gene and transcript tables
sym <- rbind(unique(iso[, c("gene_id", "gene_name")]),
             do.call(rbind, lapply(res, function(r) rbind(
               data.frame(gene_id = r$g$gene_id, gene_name = r$g$gene_symbol),
               data.frame(gene_id = r$t$gene_id, gene_name = r$t$gene_symbol)))))
sym <- sym[!is.na(sym$gene_name) & sym$gene_name != "" & !duplicated(sym$gene_id), ]

A <- data.frame(gene_id = ids, gene_name = sym$gene_name[match(ids, sym$gene_id)],
                DIU = yn(ids %in% diu), switch_ids = unname(sw_ids[ids]),
                promoter_switching_gene = yn(as.logical(annot$promoter_switching_gene[match(ids, annot$gene_id)])))
A$promoter_switching_gene[A$DIU == "no"] <- ""

for (m in MODELS) {
  g <- res[[m]]$g; t <- res[[m]]$t
  gi <- match(A$gene_id, g$gene_id)
  A[[paste0(m, "_log2FC")]] <- g$log2FC_shrunk[gi]
  A[[paste0(m, "_p")]] <- g$pvalue[gi]
  A[[paste0(m, "_padj")]] <- g$padj[gi]
  A[[paste0(m, "_geneDE")]] <- yn(g$de[gi])
  td <- t[t$de & t$gene_id %in% A$gene_id, ]
  A[[paste0(m, "_n_txDE")]] <- as.integer(table(factor(td$gene_id, A$gene_id)))
  A[[paste0(m, "_txDE_ids")]] <- unname(tapply(sprintf("%s (%s)", td$id, ifelse(td$log2FC_shrunk > 0, "up", "down")),
                                               factor(td$gene_id, A$gene_id), paste, collapse = "; "))
  gde <- A[[paste0(m, "_geneDE")]] == "yes"; tde <- A[[paste0(m, "_n_txDE")]] > 0
  A[[paste0(m, "_region")]] <- apply(cbind(gde, tde, A$DIU == "yes"), 1, function(z)
    if (!any(z)) "none" else paste(c("gene DEG", "transcript DEG", "DIU")[z], collapse = " + "))
}

# DIU genes first, then by gene-level FDR without covariates
A <- A[order(A$DIU != "yes", A$noCov_padj, A$noCov_p, na.last = TRUE), ]

# Region counts must reproduce Supplementary Figure 5
for (m in MODELS) { cat(m, "\n"); print(table(A[[paste0(m, "_region")]])) }
cat("rows:", nrow(A), "\n")

write.table(A, file.path(OUT, "SuppTable9_gene_tx_DE_DIU.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
