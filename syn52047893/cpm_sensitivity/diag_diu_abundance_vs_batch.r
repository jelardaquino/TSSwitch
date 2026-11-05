#!/usr/bin/env Rscript

# Diagnostic: separates the effect of the abundance measure (TPM vs CPM) from the effect of
# the surrogate variables IsoformSwitchAnalyzeR estimates from it during import.
# Usage: Rscript diag_diu_abundance_vs_batch.r <tpm_noSV | cpm_noSV | cpm_tpmSV>
# Otherwise identical to 01_diu_dexseq_cpm.r / isoformswitchanalyzer_DEXSeq.r.

suppressPackageStartupMessages({
  library(IsoformSwitchAnalyzeR)
})

set.seed(1)
mode <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(mode %in% c("tpm_noSV", "cpm_noSV", "cpm_tpmSV"))

BASE_DIR <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead"
INPUT_DIR <- file.path(BASE_DIR, "pychopperfq/minimap2/bam_mapq10/all_bambu_Rscript")
counts_file <- file.path(INPUT_DIR, "syn52065646LR_counts_transcript_reordered.txt")
tpm_file <- file.path(INPUT_DIR, "bambu_tx_tpm_GENEID_genesymbol.txt")
gtf_file <- file.path(INPUT_DIR, "syn52065646LR_extended_annotations.gtf")
metadata_file <- file.path(BASE_DIR, "jaquino_analysis/long_read_metadata.txt")
tpm_sv_file <- file.path(BASE_DIR, "jaquino_analysis/diu_satuRn/design_matrix_with_surrogate_variables.tsv")

outdir <- file.path(BASE_DIR, "jaquino_analysis/cpm_pipeline/diagnostics", mode)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)

counts_raw <- read.delim(counts_file, stringsAsFactors = FALSE, check.names = FALSE)
meta <- read.delim(metadata_file, stringsAsFactors = FALSE)
sample_ids <- setdiff(names(counts_raw), c("TXNAME", "GENEID", "gene_symbol", "TX"))

isoform_counts <- counts_raw[, c("TXNAME", sample_ids)]
names(isoform_counts)[1] <- "isoform_id"
for (s in sample_ids) isoform_counts[[s]] <- as.numeric(isoform_counts[[s]])

if (mode == "tpm_noSV") {
  tpm_raw <- read.delim(tpm_file, stringsAsFactors = FALSE, check.names = FALSE)
  abundance <- tpm_raw[, c("TXNAME", sample_ids)]
  names(abundance)[1] <- "isoform_id"
  for (s in sample_ids) abundance[[s]] <- as.numeric(abundance[[s]])
} else {
  m <- as.matrix(isoform_counts[, sample_ids])
  abundance <- data.frame(isoform_id = isoform_counts$isoform_id,
                          sweep(m, 2, colSums(m), "/") * 1e6, check.names = FALSE)
}

design_df <- data.frame(sampleID = sample_ids,
                        condition = meta$condition[match(sample_ids, meta$sample)],
                        stringsAsFactors = FALSE)
if (mode == "cpm_tpmSV") {
  sv <- read.delim(tpm_sv_file, stringsAsFactors = FALSE)
  design_df <- cbind(design_df, sv[match(sample_ids, sv$sampleID), c("sv1", "sv2", "sv3")])
}

switchList <- importRdata(
  isoformCountMatrix = isoform_counts,
  isoformRepExpression = abundance,
  designMatrix = design_df,
  isoformExonAnnoation = gtf_file,
  comparisonsToMake = data.frame(condition_1 = "CT", condition_2 = "AD"),
  detectUnwantedEffects = FALSE,
  showProgress = FALSE,
  quiet = FALSE,
  removeNonConvensionalChr = FALSE,
  fixStringTieAnnotationProblem = FALSE
)
switchList <- preFilter(switchList)
switchList <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = switchList, alpha = 0.05, dIFcutoff = 0.1,
  reduceToSwitchingGenes = TRUE, reduceFurtherToGenesWithConsequencePotential = FALSE,
  onlySigIsoforms = FALSE, quiet = FALSE
)

summary_tbl <- extractSwitchSummary(switchList, dIFcutoff = 0.1)
write.table(summary_tbl, "switch_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(unique(extractTopSwitches(switchList, filterForConsequences = FALSE, n = Inf)$gene_id), "sig_genes.txt")
write.table(switchList$designMatrix, "design_matrix.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(summary_tbl)
