library(IsoformSwitchAnalyzeR)
library(rtracklayer)
library(ggrepel)
library(ggplot2)
library(dplyr)

outdir <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/PRJNA1206164/diu"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)

### Import quantifications
salmonQuant <- importIsoformExpression(
    parentDir = "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/PRJNA1206164/salmon/"
)

normalize_tx_id <- function(x) {
  x <- as.character(x)
  x <- sub("\\|.*$", "", x)   # remove text after bar
  x <- sub(" .*$", "", x)       # remove text after first space
  x <- sub("\\..*$", "", x)    # remove version suffix after period
  x
}

isoform_counts <- salmonQuant$counts
colnames(isoform_counts) <- gsub("_quant", "", colnames(isoform_counts))

isoform_tpm <- salmonQuant$abundance
colnames(isoform_tpm) <- gsub("_quant", "", colnames(isoform_tpm))

# Move isoform_id column into rownames if present
if ("isoform_id" %in% colnames(isoform_counts)) {
  rownames(isoform_counts) <- isoform_counts$isoform_id
  isoform_counts$isoform_id <- NULL
}
if ("isoform_id" %in% colnames(isoform_tpm)) {
  rownames(isoform_tpm) <- isoform_tpm$isoform_id
  isoform_tpm$isoform_id <- NULL
}

rownames(isoform_counts) <- normalize_tx_id(rownames(isoform_counts))
rownames(isoform_tpm) <- normalize_tx_id(rownames(isoform_tpm))

# Remove empty/duplicated transcript IDs consistently
keep_counts <- !is.na(rownames(isoform_counts)) & nzchar(rownames(isoform_counts)) & !duplicated(rownames(isoform_counts))
keep_tpm <- !is.na(rownames(isoform_tpm)) & nzchar(rownames(isoform_tpm)) & !duplicated(rownames(isoform_tpm))
isoform_counts <- isoform_counts[keep_counts, , drop = FALSE]
isoform_tpm <- isoform_tpm[keep_tpm, , drop = FALSE]

gtf_file <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/reproducing_nextflow_pipeline/Homo_sapiens.GRCh38.107_ERCC.gtf"
# ---- GTF sanity check + cleanup for IsoformSwitchAnalyzeR ----
gr <- rtracklayer::import(gtf_file)

feat_type <- as.character(S4Vectors::mcols(gr)$type)
tx_id <- as.character(S4Vectors::mcols(gr)$transcript_id)

keep_feat <- feat_type %in% c("transcript", "exon")
keep_tx <- !is.na(tx_id) & nzchar(tx_id)

gr2 <- gr[keep_feat & keep_tx]
feat_type2 <- as.character(S4Vectors::mcols(gr2)$type)
tx_id2 <- as.character(S4Vectors::mcols(gr2)$transcript_id)
tx_id2 <- normalize_tx_id(tx_id2)
S4Vectors::mcols(gr2)$transcript_id <- tx_id2

tx_transcript <- unique(tx_id2[feat_type2 == "transcript"])
tx_exon <- unique(tx_id2[feat_type2 == "exon"])
tx_common <- intersect(tx_transcript, tx_exon)

cat("transcript IDs in transcript rows:", length(tx_transcript), "\n")
cat("transcript IDs in exon rows:", length(tx_exon), "\n")
cat("common transcript IDs:", length(tx_common), "\n")

gr_clean <- gr2[tx_id2 %in% tx_common]
clean_gtf <- file.path(outdir, "Homo_sapiens.GRCh38.107_ERCC.clean_for_ISAR.gtf")
rtracklayer::export(gr_clean, clean_gtf, format = "gtf")

# Preflight overlap check (prevents opaque importRdata failures)
annot_tx <- unique(as.character(S4Vectors::mcols(gr_clean)$transcript_id))
quant_tx <- rownames(isoform_counts)
overlap_n <- length(intersect(quant_tx, annot_tx))
cat("quantified isoforms:", length(quant_tx), "\n")
cat("annotated isoforms:", length(annot_tx), "\n")
cat("overlap isoforms:", overlap_n, "\n")
if (overlap_n == 0) {
  stop("Zero overlap between quantification and annotation transcript IDs after normalization. Check salmon index and GTF source pairing.")
}

# Restrict to shared isoforms to avoid importRdata abort on small ID mismatches
common_tx <- intersect(rownames(isoform_counts), annot_tx)
isoform_counts <- isoform_counts[common_tx, , drop = FALSE]
isoform_tpm <- isoform_tpm[common_tx, , drop = FALSE]

cat("quantified isoforms after intersect filter:", nrow(isoform_counts), "\n")



metadata <- read.table("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/PRJNA1206164/metadata_AD_CTR.txt", header=TRUE, sep="\t")
# Ensure diagnosis is character and replace "CTR" with "CT"
metadata$diagnosis[metadata$diagnosis == "CTR"] <- "CT"


design_df <- data.frame(
  sampleID = metadata$Run,
  condition = metadata$diagnosis,
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
  isoformExonAnnoation = clean_gtf,
  comparisonsToMake = comparisons_df,
  showProgress = TRUE,
  quiet = FALSE,
  removeNonConvensionalChr = FALSE,
  fixStringTieAnnotationProblem = FALSE,
  ignoreAfterBar = TRUE,
  ignoreAfterSpace = TRUE,
  ignoreAfterPeriod = TRUE
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

switchList <- isoformSwitchTestSatuRn(
  switchAnalyzeRlist = switchList,
  alpha = 0.05,
  dIFcutoff = 0.1,
  reduceToSwitchingGenes = TRUE,
  reduceFurtherToGenesWithConsequencePotential = FALSE,
  onlySigIsoforms = FALSE,
  quiet = FALSE
)

# Read gene mapping table
# gene_map <- read.delim(map_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# # Standardize names
# if ("GENEID" %in% names(gene_map)) names(gene_map)[names(gene_map) == "GENEID"] <- "gene_id"
# if ("gene_symbol" %in% names(gene_map)) names(gene_map)[names(gene_map) == "gene_symbol"] <- "gene_name"
# if ("symbol" %in% names(gene_map)) names(gene_map)[names(gene_map) == "symbol"] <- "gene_name"

# if (!all(c("gene_id", "gene_name") %in% names(gene_map))) {
#   stop("gene_map must contain gene_id and gene_name columns.")
# }

# gene_map <- unique(gene_map[, c("gene_id", "gene_name")])
# gene_map <- gene_map[!is.na(gene_map$gene_id) & gene_map$gene_id != "", ]
# gene_map <- gene_map[!is.na(gene_map$gene_name) & gene_map$gene_name != "", ]

# add_gene_name <- function(df, map_df) {
#   if (!is.data.frame(df)) return(df)

#   id_col <- intersect(c("gene_id", "geneID", "geneId"), names(df))
#   if (length(id_col) == 0) return(df)

#   id_col <- id_col[1]
#   df$gene_name <- map_df$gene_name[match(df[[id_col]], map_df$gene_id)]
#   df
# }

# # Add gene_name to internal tables if present
# if (!is.null(switchList$isoformFeatures)) {
#   switchList$isoformFeatures <- add_gene_name(switchList$isoformFeatures, gene_map)
# }
# if (!is.null(switchList$geneFeatures)) {
#   switchList$geneFeatures <- add_gene_name(switchList$geneFeatures, gene_map)
# }

# Save updated object
saveRDS(switchList, file = file.path(outdir, "switchList.rds"))

# Extract results
switch_summary <- extractSwitchSummary(
  switchList,
  dIFcutoff = 0.1
)
# switch_summary <- add_gene_name(switch_summary, gene_map)

top_switches <- extractTopSwitches(
  switchList,
  filterForConsequences = FALSE,
  n = Inf
)


# Save TSVs
write.table(
  switch_summary,
  file = file.path(outdir, "IsoformSwitchAnalyzeR_switch_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  top_switches,
  file = file.path(outdir, "IsoformSwitchAnalyzeR_top_switches.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Switch plot
library(BSgenome.Hsapiens.UCSC.hg38)
switchList <- analyzeORF(
    switchList,
    genomeObject = BSgenome.Hsapiens.UCSC.hg38
)

png(file.path(outdir, "switch_plot_ENSG00000133318.14.png"), width = 8, height = 6, units = "in", res = 300)
switchPlot(
  switchList,
  gene="ENSG00000133318.14"
)
dev.off()

longread_isoformswitch <- read.table("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/diu/IsoformSwitchAnalyzeR_top_switches_gene_name.tsv", header=TRUE, sep="\t")

top_switches$gene_id <- sub("\\..*$", "", top_switches$gene_id)

# matching long-read and short-read top switches by gene_id
match_genes <- merge(longread_isoformswitch, top_switches, by = "gene_id")
write.table(match_genes, file = file.path(outdir, "longread_vs_shortread_top_switches_comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)


# barplot of switching features
df <- data.frame(
    category = factor(c("Isoforms", "Switches", "Genes"),
                      levels = c("Isoforms", "Switches", "Genes")),
    value = c(1708, 1667, 1598)
)

ggplot(df, aes(x = category, y = value, fill = category)) +
    geom_col() +
    geom_text(aes(label = value), vjust = -0.4, size = 3.5) +
    scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) + # color-friendly palette
    labs(x = NULL, y = "number of switching features") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12)))

ggsave("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/PRJNA1206164/diu/switching_features.png", width = 4, height = 3, dpi = 300)

# volcano plot of switches
df <- switchList$isoformFeatures

# Pick one row per gene per facet first (best q-value), then top 15 per facet
top15_labels <- df %>%
  filter(!is.na(dIF), !is.na(isoform_switch_q_value)) %>%
  group_by(condition_2, gene_name) %>%  # use gene_id if gene_name is not available
  slice_min(order_by = isoform_switch_q_value, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(condition_2) %>%
  slice_min(order_by = isoform_switch_q_value, n = 15, with_ties = FALSE) %>%
  ungroup()

ggplot(data = df, aes(x = dIF, y = -log10(isoform_switch_q_value))) +
  geom_point(
    aes(color = abs(dIF) > 0.1 & isoform_switch_q_value < 0.05),
    size = 1
  ) +
  geom_text_repel(
    data = top15_labels,
    aes(label = gene_name),  # or gene_id
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    show.legend = FALSE
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed") +
  facet_wrap(~ condition_2) +
  scale_color_manual("Signficant\nIsoform Switch", values = c("black", "red")) +
  labs(x = "dIF", y = "-Log10 ( Isoform Switch Q Value )") +
  theme_classic(base_size = 12)

ggsave("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/PRJNA1206164/diu/dIF_vs_qvalue_plot.png", width=6, height=4, dpi=300) 

#Spp1
png(file.path(outdir, "switch_plot_ENSG00000118785.15.png"), width = 8, height = 6, units = "in", res = 300)
switchPlot(
  switchList,
  gene="ENSG00000118785.15"
)
dev.off()


png(file.path(outdir, "switch_plot_ENSG00000134871.19.png"), width = 8, height = 6, units = "in", res = 300)
switchPlot(
  switchList,
  gene="ENSG00000134871.19"
)
dev.off()
