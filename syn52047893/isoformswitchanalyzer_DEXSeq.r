library(IsoformSwitchAnalyzeR)

outdir <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/diu"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)

counts_file <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/pychopperfq/minimap2/bam_mapq10/all_bambu_Rscript/syn52065646LR_counts_transcript_reordered.txt"
tpm_file    <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/pychopperfq/minimap2/bam_mapq10/all_bambu_Rscript/bambu_tx_tpm_GENEID_genesymbol.txt"
gtf_file    <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/pychopperfq/minimap2/bam_mapq10/all_bambu_Rscript/syn52065646LR_extended_annotations.gtf"

# Map file with gene_id -> gene_name
map_file <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/pychopperfq/minimap2/bam_mapq10/all_bambu_Rscript/bambu_tx_tpm_GENEID_genesymbol.txt"

counts_raw <- read.delim(counts_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
tpm_raw    <- read.delim(tpm_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

if (!"TXNAME" %in% names(counts_raw)) stop("counts file must contain TXNAME")
if (!"TXNAME" %in% names(tpm_raw)) stop("TPM file must contain TXNAME")

sample_ids <- setdiff(names(counts_raw), c("TXNAME", "GENEID", "gene_symbol", "TX"))
group <- factor(c(rep("CT", 6), rep("AD", 6)), levels = c("CT", "AD"))

if (length(sample_ids) != length(group)) {
  stop("Sample column count does not match the group labels.")
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
  condition = as.character(group),
  stringsAsFactors = FALSE
)

comparisons_df <- data.frame(
  condition_1 = "CT",
  condition_2 = "AD",
  stringsAsFactors = FALSE
)

switchList <- importRdata(
  isoformCountMatrix = isoform_counts,
  isoformRepExpression = isoform_tpm,
  designMatrix = design_df,
  isoformExonAnnoation = gtf_file,
  comparisonsToMake = comparisons_df,
  showProgress = TRUE,
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

# Read gene mapping table
gene_map <- read.delim(map_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# Standardize names
if ("GENEID" %in% names(gene_map)) names(gene_map)[names(gene_map) == "GENEID"] <- "gene_id"
if ("gene_symbol" %in% names(gene_map)) names(gene_map)[names(gene_map) == "gene_symbol"] <- "gene_name"
if ("symbol" %in% names(gene_map)) names(gene_map)[names(gene_map) == "symbol"] <- "gene_name"

if (!all(c("gene_id", "gene_name") %in% names(gene_map))) {
  stop("gene_map must contain gene_id and gene_name columns.")
}

gene_map <- unique(gene_map[, c("gene_id", "gene_name")])
gene_map <- gene_map[!is.na(gene_map$gene_id) & gene_map$gene_id != "", ]
gene_map <- gene_map[!is.na(gene_map$gene_name) & gene_map$gene_name != "", ]

add_gene_name <- function(df, map_df) {
  if (!is.data.frame(df)) return(df)

  id_col <- intersect(c("gene_id", "geneID", "geneId"), names(df))
  if (length(id_col) == 0) return(df)

  id_col <- id_col[1]
  df$gene_name <- map_df$gene_name[match(df[[id_col]], map_df$gene_id)]
  df
}

# Add gene_name to internal tables if present
if (!is.null(switchList$isoformFeatures)) {
  switchList$isoformFeatures <- add_gene_name(switchList$isoformFeatures, gene_map)
}
if (!is.null(switchList$geneFeatures)) {
  switchList$geneFeatures <- add_gene_name(switchList$geneFeatures, gene_map)
}

# Save updated object
saveRDS(switchList, file = file.path(outdir, "switchList_with_gene_name.rds"))

# Extract results
switch_summary <- extractSwitchSummary(
  switchList,
  dIFcutoff = 0.1
)
switch_summary <- add_gene_name(switch_summary, gene_map)

top_switches <- extractTopSwitches(
  switchList,
  filterForConsequences = FALSE,
  n = Inf
)
top_switches <- add_gene_name(top_switches, gene_map)

# Save TSVs
write.table(
  switch_summary,
  file = file.path(outdir, "IsoformSwitchAnalyzeR_switch_summary_gene_name.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  top_switches,
  file = file.path(outdir, "IsoformSwitchAnalyzeR_top_switches_gene_name.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Switch plot
library(BSgenome.Hsapiens.UCSC.hg38)
switchList <- analyzeORF(
    switchList,
    genomeObject = BSgenome.Hsapiens.UCSC.hg38
)

png(file.path(outdir, "switch_plot_ENSG00000183648.png"), width = 8, height = 6, units = "in", res = 300)
switchPlot(
  switchList,
  gene="ENSG00000183648"
)
dev.off()
