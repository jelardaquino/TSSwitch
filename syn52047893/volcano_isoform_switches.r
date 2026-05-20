library(IsoformSwitchAnalyzeR)

mySwitchList <- readRDS("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/diu/switchList_with_gene_name.rds")

extractConsequenceEnrichment(
    mySwitchList,
    analysisOppositeConsequence = TRUE,
    localTheme = theme_bw(base_size = 14), # Increase font size in vignette
    returnResult = FALSE # if TRUE returns a data.frame with the summary statistics
)

mySwitchList <- analyzeAlternativeSplicing( SwitchList )

# Switch plot
library(BSgenome.Hsapiens.UCSC.hg38)
mySwitchList <- analyzeORF(
    mySwitchList,
    genomeObject = BSgenome.Hsapiens.UCSC.hg38
)

mySwitchList <- analyzeNovelIsoformORF( mySwitchList )
mySwitchList <- analyzeCPC2(mySwitchList)
mySwitchList <- analyzeCPAT( mySwitchList, pathToCPATresultFile = "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/diu/cpat_result.txt")
mySwitchList <- analyzeSwitchConsequences( mySwitchList )

extractSplicingSummary(
    mySwitchList,
    asFractionTotal = FALSE,
    plotGenes=FALSE,
    localTheme = theme_classic(base_size = 14)
)

splicingEnrichment <- extractSplicingEnrichment(
    mySwitchList,
    splicingToAnalyze='all',
    returnResult=TRUE,
    returnSummary=TRUE
)

extractSplicingGenomeWide(
    mySwitchList,
    featureToExtract = 'all',                 # all isoforms stored in the switchAnalyzeRlist
    splicingToAnalyze = c('A3','MES','ATSS'), # Splice types significantly enriched in COAD
    plot=TRUE,
    returnResult=FALSE  # Preventing the summary statistics to be returned as a data.frame
)

library(dplyr)
library(ggplot2)
library(ggrepel)

df <- mySwitchList$isoformFeatures

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

ggsave("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/diu/dIF_vs_qvalue_plot.png", width=6, height=4, dpi=300)  
