#!/usr/bin/env Rscript

# DIU with isoformSwitchTestSatuRn, using the same inputs, import settings and
# preFilter() as isoformswitchanalyzer_DEXSeq.r (the manuscript DIU run), so the
# two tests can be compared directly. All tested genes are kept
# (reduceToSwitchingGenes = FALSE) so every gene has a satuRn q-value.

suppressPackageStartupMessages({
  library(IsoformSwitchAnalyzeR)
})

set.seed(1)

BASE_DIR <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead"
INPUT_DIR <- file.path(BASE_DIR, "pychopperfq/minimap2/bam_mapq10/all_bambu_Rscript")

counts_file <- file.path(INPUT_DIR, "syn52065646LR_counts_transcript_reordered.txt")
tpm_file <- file.path(INPUT_DIR, "bambu_tx_tpm_GENEID_genesymbol.txt")
gtf_file <- file.path(INPUT_DIR, "syn52065646LR_extended_annotations.gtf")
metadata_file <- file.path(BASE_DIR, "jaquino_analysis/long_read_metadata.txt")

dexseq_top_file <- file.path(BASE_DIR, "jaquino_analysis/diu/IsoformSwitchAnalyzeR_top_switches_gene_name.tsv")
previous_saturn_top_file <- file.path(
  BASE_DIR,
  "pychopperfq/minimap2/bam_mapq10/all_bambu_Rscript/all_bambu_python_analysis/IsoformSwitchAnalyzeR/ISO_resultSatuRn/IsoformSwitchAnalyzeR_top_switches_gene_name.tsv"
)

outdir <- file.path(BASE_DIR, "jaquino_analysis/diu_satuRn")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)

ALPHA <- 0.05
DIF_CUTOFF <- 0.1

###############################
# Inputs
###############################

counts_raw <- read.delim(counts_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
tpm_raw <- read.delim(tpm_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
meta <- read.delim(metadata_file, header = TRUE, stringsAsFactors = FALSE)

sample_ids <- setdiff(names(counts_raw), c("TXNAME", "GENEID", "gene_symbol", "TX"))

if (!setequal(sample_ids, meta$sample)) {
  stop("Sample IDs in the count matrix do not match long_read_metadata.txt")
}
if (!all(sample_ids %in% names(tpm_raw))) {
  stop("TPM file is missing samples present in the count matrix")
}

isoform_counts <- counts_raw[, c("TXNAME", sample_ids), drop = FALSE]
names(isoform_counts)[1] <- "isoform_id"

isoform_tpm <- tpm_raw[, c("TXNAME", sample_ids), drop = FALSE]
names(isoform_tpm)[1] <- "isoform_id"

for (s in sample_ids) {
  isoform_counts[[s]] <- as.numeric(isoform_counts[[s]])
  isoform_tpm[[s]] <- as.numeric(isoform_tpm[[s]])
}

design_df <- data.frame(
  sampleID = sample_ids,
  condition = meta$condition[match(sample_ids, meta$sample)],
  stringsAsFactors = FALSE
)

comparisons_df <- data.frame(condition_1 = "CT", condition_2 = "AD", stringsAsFactors = FALSE)

gene_map <- unique(tpm_raw[, c("GENEID", "gene_symbol")])
names(gene_map) <- c("gene_id", "gene_name")
gene_map <- gene_map[!is.na(gene_map$gene_id) & gene_map$gene_id != "" &
                       !is.na(gene_map$gene_name) & gene_map$gene_name != "", ]
gene_map <- gene_map[!duplicated(gene_map$gene_id), ]

###############################
# Import, filter, test (same import/filter settings as the DEXSeq run)
###############################

switchList <- importRdata(
  isoformCountMatrix = isoform_counts,
  isoformRepExpression = isoform_tpm,
  designMatrix = design_df,
  isoformExonAnnoation = gtf_file,
  comparisonsToMake = comparisons_df,
  showProgress = FALSE,
  quiet = FALSE,
  removeNonConvensionalChr = FALSE,
  fixStringTieAnnotationProblem = FALSE
)

switchList <- preFilter(switchList)

pdf(file.path(outdir, "satuRn_diagnostic_plots.pdf"))
switchList <- isoformSwitchTestSatuRn(
  switchAnalyzeRlist = switchList,
  alpha = ALPHA,
  dIFcutoff = DIF_CUTOFF,
  reduceToSwitchingGenes = FALSE,
  reduceFurtherToGenesWithConsequencePotential = FALSE,
  onlySigIsoforms = FALSE,
  showProgress = FALSE,
  quiet = FALSE
)
dev.off()

switchList$isoformFeatures$gene_name <- gene_map$gene_name[match(switchList$isoformFeatures$gene_id, gene_map$gene_id)]

saveRDS(switchList, file.path(outdir, "switchList_satuRn_all_genes.rds"))
write.table(switchList$designMatrix, file.path(outdir, "design_matrix_with_surrogate_variables.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

###############################
# Isoform- and gene-level results
###############################

iso <- switchList$isoformFeatures
iso_cols <- c("isoform_id", "gene_id", "gene_name", "condition_1", "condition_2",
              "iso_value_1", "iso_value_2", "IF1", "IF2", "dIF",
              "isoform_switch_q_value", "gene_switch_q_value")
iso <- iso[, intersect(iso_cols, names(iso))]
iso$significant_isoform <- !is.na(iso$isoform_switch_q_value) &
  iso$isoform_switch_q_value < ALPHA & abs(iso$dIF) > DIF_CUTOFF
iso <- iso[order(iso$isoform_switch_q_value, -abs(iso$dIF)), ]

write.table(iso, file.path(outdir, "isoform_results_satuRn_all_tested.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(iso[iso$significant_isoform, ], file.path(outdir, "significant_isoforms_satuRn.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

min_or_na <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)

genes <- do.call(rbind, lapply(split(iso, iso$gene_id), function(g) {
  big_dif <- abs(g$dIF) > DIF_CUTOFF
  data.frame(
    gene_id = g$gene_id[1],
    gene_name = g$gene_name[1],
    n_isoforms_tested = nrow(g),
    gene_switch_q_value = g$gene_switch_q_value[1],
    min_isoform_q_with_dIF_gt_cutoff = min_or_na(g$isoform_switch_q_value[big_dif]),
    max_abs_dIF = max(abs(g$dIF), na.rm = TRUE),
    n_significant_isoforms = sum(g$significant_isoform),
    significant_gene = any(g$significant_isoform),
    stringsAsFactors = FALSE
  )
}))
genes <- genes[order(genes$min_isoform_q_with_dIF_gt_cutoff), ]

write.table(genes, file.path(outdir, "gene_results_satuRn_all_tested.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

switch_summary <- extractSwitchSummary(switchList, alpha = ALPHA, dIFcutoff = DIF_CUTOFF)
write.table(switch_summary, file.path(outdir, "IsoformSwitchAnalyzeR_switch_summary_satuRn.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

###############################
# Comparison with the manuscript DEXSeq run and the earlier satuRn run
###############################

dexseq_top <- read.delim(dexseq_top_file, stringsAsFactors = FALSE)
previous_saturn_top <- read.delim(previous_saturn_top_file, stringsAsFactors = FALSE)

saturn_sig_ids <- genes$gene_id[genes$significant_gene]
all_ids <- union(dexseq_top$gene_id, saturn_sig_ids)

q_bin <- function(q) {
  ifelse(is.na(q), "not tested",
  ifelse(q < 0.05, "< 0.05",
  ifelse(q < 0.10, "0.05-0.10",
  ifelse(q < 0.25, "0.10-0.25", ">= 0.25"))))
}

comparison <- data.frame(gene_id = all_ids, stringsAsFactors = FALSE)
comparison$gene_name <- gene_map$gene_name[match(comparison$gene_id, gene_map$gene_id)]
comparison$significant_DEXSeq_manuscript <- comparison$gene_id %in% dexseq_top$gene_id
comparison$significant_satuRn_rerun <- comparison$gene_id %in% saturn_sig_ids
comparison$significant_satuRn_previous_run <- comparison$gene_id %in% previous_saturn_top$gene_id
comparison$DEXSeq_gene_q_value <- dexseq_top$gene_switch_q_value[match(comparison$gene_id, dexseq_top$gene_id)]
m <- match(comparison$gene_id, genes$gene_id)
comparison$satuRn_gene_q_value <- genes$gene_switch_q_value[m]
comparison$satuRn_min_isoform_q_with_dIF_gt_cutoff <- genes$min_isoform_q_with_dIF_gt_cutoff[m]
comparison$satuRn_q_bin <- q_bin(comparison$satuRn_min_isoform_q_with_dIF_gt_cutoff)
comparison$satuRn_max_abs_dIF <- genes$max_abs_dIF[m]
comparison$category <- ifelse(
  comparison$significant_DEXSeq_manuscript & comparison$significant_satuRn_rerun, "Both",
  ifelse(comparison$significant_DEXSeq_manuscript, "DEXSeq only", "satuRn only")
)
comparison <- comparison[order(factor(comparison$category, levels = c("Both", "DEXSeq only", "satuRn only")),
                               comparison$satuRn_min_isoform_q_with_dIF_gt_cutoff), ]

write.table(comparison, file.path(outdir, "comparison_DEXSeq_vs_satuRn_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

###############################
# Run record
###############################

dex_ids <- dexseq_top$gene_id
prev_ids <- previous_saturn_top$gene_id
dex_only <- comparison[comparison$category == "DEXSeq only", ]

report <- c(
  paste("Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("IsoformSwitchAnalyzeR", as.character(packageVersion("IsoformSwitchAnalyzeR")),
        "| satuRn", as.character(packageVersion("satuRn"))),
  paste("alpha =", ALPHA, "| dIF cutoff =", DIF_CUTOFF),
  paste("Genes tested:", nrow(genes), "| isoforms tested:", nrow(iso)),
  "",
  "extractSwitchSummary():",
  capture.output(print(switch_summary, row.names = FALSE)),
  "",
  paste("satuRn rerun significant genes:", length(saturn_sig_ids)),
  paste("Manuscript DEXSeq significant genes:", length(dex_ids)),
  paste("Shared with DEXSeq:", length(intersect(dex_ids, saturn_sig_ids))),
  paste("Previous satuRn run significant genes:", length(prev_ids),
        "| shared with rerun:", length(intersect(prev_ids, saturn_sig_ids))),
  "",
  "satuRn q-value bins for the DEXSeq-only genes (min isoform q among isoforms with |dIF| > cutoff):",
  capture.output(print(table(dex_only$satuRn_q_bin))),
  "",
  capture.output(sessionInfo())
)
writeLines(report, file.path(outdir, "run_summary.txt"))
cat(report[1:20], sep = "\n")
