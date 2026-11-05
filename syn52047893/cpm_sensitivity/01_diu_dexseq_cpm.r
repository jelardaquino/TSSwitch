#!/usr/bin/env Rscript

# Step 1: DIU with IsoformSwitchAnalyzeR + DEXSeq using CPM as the isoform abundance.
# Copy of scripts/isoformswitchanalyzer_DEXSeq.r (manuscript DIU). The only analytical
# change is isoformRepExpression: CPM computed from the same Bambu counts, instead of TPM.
# The three surrogate variables are fixed to those of the manuscript run: importRdata()
# otherwise re-estimates them from the abundance input, and CPM-derived estimates changed
# the result far more than CPM itself (see cpm_pipeline/diagnostics).
# Also writes the downstream inputs the manuscript pipeline used (switchList.txt,
# sig_transcripts.txt, diu_sig_genes.txt) and the ORF/PTC table used for NMD calls.

suppressPackageStartupMessages({
  library(IsoformSwitchAnalyzeR)
})

set.seed(1)

BASE_DIR <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead"
INPUT_DIR <- file.path(BASE_DIR, "pychopperfq/minimap2/bam_mapq10/all_bambu_Rscript")

counts_file <- file.path(INPUT_DIR, "syn52065646LR_counts_transcript_reordered.txt")
gtf_file <- file.path(INPUT_DIR, "syn52065646LR_extended_annotations.gtf")
map_file <- file.path(INPUT_DIR, "bambu_tx_tpm_GENEID_genesymbol.txt")
metadata_file <- file.path(BASE_DIR, "jaquino_analysis/long_read_metadata.txt")
manuscript_diu_rds <- file.path(BASE_DIR, "jaquino_analysis/diu/switchList_with_gene_name.rds")

outdir <- file.path(BASE_DIR, "jaquino_analysis/cpm_pipeline/diu")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)

counts_raw <- read.delim(counts_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
meta <- read.delim(metadata_file, stringsAsFactors = FALSE)

sample_ids <- setdiff(names(counts_raw), c("TXNAME", "GENEID", "gene_symbol", "TX"))
if (!setequal(sample_ids, meta$sample)) stop("Count matrix samples do not match long_read_metadata.txt")

isoform_counts <- counts_raw[, c("TXNAME", sample_ids), drop = FALSE]
names(isoform_counts)[1] <- "isoform_id"
for (s in sample_ids) isoform_counts[[s]] <- as.numeric(isoform_counts[[s]])

count_mat <- as.matrix(isoform_counts[, sample_ids])
cpm_mat <- sweep(count_mat, 2, colSums(count_mat), FUN = "/") * 1e6
isoform_cpm <- data.frame(isoform_id = isoform_counts$isoform_id, cpm_mat, check.names = FALSE)

design_df <- data.frame(
  sampleID = sample_ids,
  condition = meta$condition[match(sample_ids, meta$sample)],
  stringsAsFactors = FALSE
)
manuscript_sv <- readRDS(manuscript_diu_rds)$designMatrix
design_df <- cbind(design_df, manuscript_sv[match(sample_ids, manuscript_sv$sampleID), c("sv1", "sv2", "sv3")])

comparisons_df <- data.frame(condition_1 = "CT", condition_2 = "AD", stringsAsFactors = FALSE)

switchList <- importRdata(
  isoformCountMatrix = isoform_counts,
  isoformRepExpression = isoform_cpm,
  designMatrix = design_df,
  isoformExonAnnoation = gtf_file,
  comparisonsToMake = comparisons_df,
  detectUnwantedEffects = FALSE,
  showProgress = FALSE,
  quiet = FALSE,
  removeNonConvensionalChr = FALSE,
  fixStringTieAnnotationProblem = FALSE
)

switchList <- preFilter(switchList)

switchList <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = switchList,
  alpha = 0.05,
  dIFcutoff = 0.1,
  reduceToSwitchingGenes = TRUE,
  reduceFurtherToGenesWithConsequencePotential = FALSE,
  onlySigIsoforms = FALSE,
  quiet = FALSE
)

gene_map <- read.delim(map_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
names(gene_map)[names(gene_map) == "GENEID"] <- "gene_id"
names(gene_map)[names(gene_map) == "gene_symbol"] <- "gene_name"
gene_map <- unique(gene_map[, c("gene_id", "gene_name")])
gene_map <- gene_map[!is.na(gene_map$gene_id) & gene_map$gene_id != "" &
                       !is.na(gene_map$gene_name) & gene_map$gene_name != "", ]

add_gene_name <- function(df, map_df) {
  if (!is.data.frame(df)) return(df)
  id_col <- intersect(c("gene_id", "geneID", "geneId"), names(df))
  if (length(id_col) == 0) return(df)
  df$gene_name <- map_df$gene_name[match(df[[id_col[1]]], map_df$gene_id)]
  df
}

switchList$isoformFeatures <- add_gene_name(switchList$isoformFeatures, gene_map)

saveRDS(switchList, file = file.path(outdir, "switchList_with_gene_name.rds"))
write.table(switchList$designMatrix, file.path(outdir, "design_matrix_with_surrogate_variables.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

switch_summary <- extractSwitchSummary(switchList, dIFcutoff = 0.1)
top_switches <- add_gene_name(extractTopSwitches(switchList, filterForConsequences = FALSE, n = Inf), gene_map)

write.table(switch_summary, file.path(outdir, "IsoformSwitchAnalyzeR_switch_summary_gene_name.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(top_switches, file.path(outdir, "IsoformSwitchAnalyzeR_top_switches_gene_name.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Downstream inputs, built the same way as the manuscript files
write.table(switchList$isoformFeatures, file.path(outdir, "switchList.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c("isoform_id", unique(switchList$isoformFeatures$isoform_id)),
           file.path(outdir, "sig_transcripts.txt"))
writeLines(unique(top_switches$gene_id), file.path(outdir, "diu_sig_genes.txt"))

# ORF prediction and PTC (NMD-sensitivity) calls, as in the manuscript
library(BSgenome.Hsapiens.UCSC.hg38)
switchList <- analyzeORF(switchList, genomeObject = BSgenome.Hsapiens.UCSC.hg38)
saveRDS(switchList, file = file.path(outdir, "switchList_with_ORF.rds"))
write.table(switchList$isoformFeatures, file.path(outdir, "isoform_PTC_status.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

iso <- switchList$isoformFeatures
sig_iso <- iso[!is.na(iso$isoform_switch_q_value) & iso$isoform_switch_q_value < 0.05 & abs(iso$dIF) > 0.1, ]
sig_iso <- sig_iso[!duplicated(sig_iso$isoform_id), ]
nmd <- c(
  paste("Switching isoforms:", nrow(sig_iso)),
  paste("PTC status assigned:", sum(!is.na(sig_iso$PTC))),
  paste("Predicted NMD-sensitive (PTC = TRUE):", sum(sig_iso$PTC %in% TRUE)),
  paste("  increased usage in AD (dIF > 0):", sum(sig_iso$PTC %in% TRUE & sig_iso$dIF > 0)),
  paste("  increased usage in CT (dIF < 0):", sum(sig_iso$PTC %in% TRUE & sig_iso$dIF < 0))
)
writeLines(nmd, file.path(outdir, "nmd_summary.txt"))

print(switch_summary)
cat(nmd, sep = "\n")
