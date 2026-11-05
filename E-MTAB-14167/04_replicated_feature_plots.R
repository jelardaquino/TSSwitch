###############################
# Per-sample plots for the long-read findings that were significant again in the
# 41-donor E-MTAB-14167 short-read cohort:
#   - isoform switches significant in the short-read IsoformSwitchAnalyzeR run
#     (CLSTN3-01, CLSTN3-02, COL4A2), shown as raw TPM fractions and as the
#     covariate-adjusted fractions IsoformSwitchAnalyzeR tests and reports (limma::removeBatchEffect on age, sex and SVs)
#   - the two promoters with the strongest concordant short-read effect (GSDMB, NTHL1)
#
# Run with the isoformswitch_env R (no data.table there, so base R only).
###############################

suppressPackageStartupMessages({
  library(IsoformSwitchAnalyzeR)
  library(ggplot2)
})

S <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing"
LR <- file.path(S, "syn52065646_LongRead/jaquino_analysis")
OUT <- file.path(S, "ERP161086/emtab14167_replication")
FIGS <- file.path(OUT, "figures")

BAMBU_GTF <- file.path(LR, "NanoporeRNASeq_extended_annotations.gtf")
LR_TPM <- file.path(LR, "bambu_tx_tpm.txt")
LR_META <- file.path(LR, "long_read_metadata.txt")
LR_SWITCHLIST_RDS <- file.path(LR, "diu/switchList_with_gene_name.rds")
LR_PROMOTER <- file.path(S, "syn52065646_LongRead/jaquino_analysis_v2/tables/promoter_usage_86_genes_12samples/promoter_usage_betareg_results_9samples.tsv")
SR_SWITCHLIST_RDS <- file.path(OUT, "diu_bambu/switchList_SR41_bambu_annotation_all_genes.rds")
SR_QUANT_DIR <- file.path(OUT, "salmon_bambu/quant")
TARGETED <- file.path(OUT, "tables/targeted_replication_all_features.tsv")

COLORS <- c(CT = "#00BFC4", AD = "#F8766D")

ISOFORMS <- data.frame(
  isoform_id = c("ENST00000266546", "ENST00000331148", "ENST00000648222"),
  gene_id = c("ENSG00000139182", "ENSG00000139182", "ENSG00000134871"),
  label = c("CLSTN3-01", "CLSTN3-02", "COL4A2 (ENST00000648222)")
)
PROMOTERS <- data.frame(
  promoter_id = c("17:39916316-39917416(-)", "16:2043079-2044179(-)"),
  gene_id = c("ENSG00000073605", "ENSG00000065057"),
  label = c("GSDMB promoter\nchr17:39,916,316-39,917,416 (-)", "NTHL1 promoter\nchr16:2,043,079-2,044,179 (-)")
)
GENES <- unique(c(ISOFORMS$gene_id, PROMOTERS$gene_id))

###############################
# Annotation for the four genes
###############################

gtf_cmd <- sprintf("awk -F'\\t' '$3 == \"transcript\" && $9 ~ /%s/' %s",
                   paste(GENES, collapse = "|"), BAMBU_GTF)
gtf <- read.delim(pipe(gtf_cmd), header = FALSE, quote = "")
tx <- data.frame(
  transcript_id = sub('.*transcript_id "([^"]+)".*', "\\1", gtf$V9),
  gene_id = sub('.*gene_id "([^"]+)".*', "\\1", gtf$V9),
  chr = gtf$V1, start = gtf$V4, end = gtf$V5, strand = gtf$V7
)
tx <- tx[tx$gene_id %in% GENES, ]
tss <- ifelse(tx$strand == "-", tx$end, tx$start)
tx$promoter_id <- paste0(tx$chr, ":", pmax(1L, tss - 1000L), "-", tss + 100L, "(", tx$strand, ")")

###############################
# Expression matrices (transcripts of the four genes)
###############################

lr_meta <- read.delim(LR_META)
lr_meta$condition <- factor(lr_meta$condition, levels = c("CT", "AD"))
lr_tpm <- as.matrix(read.delim(LR_TPM, check.names = FALSE))
colnames(lr_tpm) <- paste0("X", colnames(lr_tpm))
lr_tpm <- lr_tpm[intersect(tx$transcript_id, rownames(lr_tpm)), lr_meta$sample]

sr <- readRDS(SR_SWITCHLIST_RDS)
sr_meta <- sr$designMatrix[, c("sampleID", "condition")]
names(sr_meta)[1] <- "sample"
sr_meta$condition <- factor(sr_meta$condition, levels = c("CT", "AD"))
sr_tpm <- sapply(sr_meta$sample, function(s) {
  q <- read.delim(file.path(SR_QUANT_DIR, paste0(s, "_quant"), "quant.sf"))
  setNames(q$TPM, q$Name)[tx$transcript_id]
})
rownames(sr_tpm) <- tx$transcript_id
sr_tpm[is.na(sr_tpm)] <- 0

# Raw isoform fraction: transcript TPM / gene TPM over all Bambu transcripts, NA where gene TPM < 1
raw_fraction <- function(tpm, iso, gene) {
  gt <- colSums(tpm[intersect(tx$transcript_id[tx$gene_id == gene], rownames(tpm)), , drop = FALSE])
  f <- tpm[iso, ] / gt
  f[gt < 1] <- NA
  f
}

# Promoter usage as in 02_targeted_replication.R: promoter TPM / gene TPM, NA where gene TPM < 1
promoter_usage <- function(tpm, prom, gene) {
  g_tx <- intersect(tx$transcript_id[tx$gene_id == gene], rownames(tpm))
  p_tx <- intersect(tx$transcript_id[tx$promoter_id == prom], g_tx)
  gt <- colSums(tpm[g_tx, , drop = FALSE])
  u <- colSums(tpm[p_tx, , drop = FALSE]) / gt
  u[gt < 1] <- NA
  u
}

rep_if <- function(sl, iso) {
  d <- sl$designMatrix
  unlist(sl$isoformRepIF[sl$isoformRepIF$isoform_id == iso, d$sampleID])
}

long_df <- function(values, meta, feature, panel) {
  data.frame(sample = meta$sample, condition = meta$condition,
             value = as.numeric(values[meta$sample]), feature = feature, panel = panel)
}

###############################
# Statistics for the panel labels
###############################

lr <- readRDS(LR_SWITCHLIST_RDS)
targeted <- read.delim(TARGETED)
targeted <- targeted[targeted$set == "all41", ]
lr_prom <- read.delim(LR_PROMOTER)

fmt_p <- function(x) ifelse(x < 0.001, formatC(x, format = "e", digits = 1), formatC(x, format = "f", digits = 3))
fmt_d <- function(x) sprintf("%+.2f", x)
mean_diff <- function(df) {
  m <- tapply(df$value, df$condition, mean, na.rm = TRUE)
  m[["AD"]] - m[["CT"]]
}

###############################
# Figure: isoform switches
###############################

PANELS <- c("Long-read\nraw TPM fraction", "Long-read\nadjusted: SVs (tested)",
            "Short-read\nraw TPM fraction", "Short-read\nadjusted: age, sex, SVs (tested)")
iso_df <- list()
iso_lab <- list()
for (i in seq_len(nrow(ISOFORMS))) {
  iso <- ISOFORMS$isoform_id[i]
  gene <- ISOFORMS$gene_id[i]
  lab <- ISOFORMS$label[i]
  lr_f <- lr$isoformFeatures[lr$isoformFeatures$isoform_id == iso, ][1, ]
  sr_f <- sr$isoformFeatures[sr$isoformFeatures$isoform_id == iso, ][1, ]
  tg <- targeted[targeted$family == "isoform" & targeted$feature_id == iso, ]
  d <- list(
    long_df(raw_fraction(lr_tpm, iso, gene), lr_meta, lab, PANELS[1]),
    long_df(rep_if(lr, iso), data.frame(sample = lr$designMatrix$sampleID, condition = factor(lr$designMatrix$condition, levels = c("CT", "AD"))), lab, PANELS[2]),
    long_df(raw_fraction(sr_tpm, iso, gene), sr_meta, lab, PANELS[3]),
    long_df(rep_if(sr, iso), sr_meta, lab, PANELS[4])
  )
  iso_df <- c(iso_df, d)
  iso_lab <- c(iso_lab, list(
    data.frame(feature = lab, panel = PANELS[1], text = paste0("dIF ", fmt_d(mean_diff(d[[1]])), "\nnot tested")),
    data.frame(feature = lab, panel = PANELS[2], text = paste0("dIF ", fmt_d(lr_f$dIF), "\nFDR ", fmt_p(lr_f$isoform_switch_q_value))),
    data.frame(feature = lab, panel = PANELS[3], text = paste0("dIF ", fmt_d(mean_diff(d[[3]])), "\nbeta reg. p ", fmt_p(tg$p))),
    data.frame(feature = lab, panel = PANELS[4], text = paste0("dIF ", fmt_d(sr_f$dIF), "\nFDR ", fmt_p(sr_f$isoform_switch_q_value)))
  ))
}
iso_df <- do.call(rbind, iso_df)
iso_lab <- do.call(rbind, iso_lab)
for (x in c("iso_df", "iso_lab")) {
  y <- get(x)
  y$panel <- factor(y$panel, levels = PANELS)
  y$feature <- factor(y$feature, levels = ISOFORMS$label)
  assign(x, y)
}

strip_plot <- function(df, lab, ylab) {
  ggplot(df, aes(condition, value)) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.55, linewidth = 0.3, colour = "grey25") +
    geom_jitter(aes(fill = condition), shape = 21, colour = "white", stroke = 0.4,
                size = 2.2, width = 0.15, height = 0, alpha = 0.9, na.rm = TRUE) +
    geom_text(data = lab, aes(x = 0.5, y = 1.25, label = text), hjust = 0, vjust = 1,
              size = 2.7, lineheight = 0.9, colour = "grey15", inherit.aes = FALSE) +
    scale_fill_manual(values = COLORS, name = NULL) +
    scale_y_continuous(limits = c(0, 1.26), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0.02, 0))) +
    labs(x = NULL, y = ylab) +
    theme_bw(base_size = 9) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
          panel.border = element_rect(colour = "grey70"),
          strip.background = element_rect(fill = "grey95", colour = "grey70"),
          strip.text = element_text(size = 8), legend.position = "bottom")
}

set.seed(1)
p_iso <- strip_plot(iso_df, iso_lab, "Isoform fraction") +
  facet_grid(feature ~ panel) +
  labs(title = "Isoform switches significant in both long-read and short-read data",
       subtitle = "Long-read: 6 AD vs 6 CT. Short-read: E-MTAB-14167, 20 AD vs 21 CT. Bars are group means.") +
  theme(plot.title = element_text(size = 10, face = "bold"), plot.subtitle = element_text(size = 8))
ggsave(file.path(FIGS, "replicated_isoform_switches.png"), p_iso, width = 8.5, height = 7, dpi = 300)
ggsave(file.path(FIGS, "replicated_isoform_switches.pdf"), p_iso, width = 8.5, height = 7)

# Same figure with only the adjusted (tested) fractions
ADJ <- PANELS[c(2, 4)]
set.seed(1)
p_iso_adj <- strip_plot(droplevels(iso_df[iso_df$panel %in% ADJ, ]), droplevels(iso_lab[iso_lab$panel %in% ADJ, ]),
                        "Adjusted isoform fraction") +
  facet_grid(feature ~ panel) +
  labs(title = "Isoform switches significant in both\nlong-read and short-read data",
       subtitle = "Long-read: 6 AD vs 6 CT. Short-read: E-MTAB-14167, 20 AD vs 21 CT.\nFractions after IsoformSwitchAnalyzeR covariate adjustment. Bars are group means.") +
  theme(plot.title = element_text(size = 10, face = "bold"), plot.subtitle = element_text(size = 8))
ggsave(file.path(FIGS, "replicated_isoform_switches_adjusted.png"), p_iso_adj, width = 4.8, height = 7, dpi = 300)
ggsave(file.path(FIGS, "replicated_isoform_switches_adjusted.pdf"), p_iso_adj, width = 4.8, height = 7)

###############################
# Figure: promoters
###############################

PPANELS <- c("Long-read", "Short-read")
pr_df <- list()
pr_lab <- list()
for (i in seq_len(nrow(PROMOTERS))) {
  prom <- PROMOTERS$promoter_id[i]
  gene <- PROMOTERS$gene_id[i]
  lab <- PROMOTERS$label[i]
  lp <- lr_prom[lr_prom$promoter_id == prom, ]
  tg <- targeted[targeted$family == "promoter" & targeted$feature_id == prom, ]
  d_lr <- long_df(promoter_usage(lr_tpm, prom, gene), lr_meta, lab, PPANELS[1])
  d_sr <- long_df(promoter_usage(sr_tpm, prom, gene), sr_meta, lab, PPANELS[2])
  pr_df <- c(pr_df, list(d_lr, d_sr))
  pr_lab <- c(pr_lab, list(
    data.frame(feature = lab, panel = PPANELS[1], text = paste0("diff ", fmt_d(mean_diff(d_lr)), "\nFDR ", fmt_p(lp$padj))),
    data.frame(feature = lab, panel = PPANELS[2], text = paste0("diff ", fmt_d(mean_diff(d_sr)), "\np ", fmt_p(tg$p), ", FDR ", fmt_p(tg$fdr)))
  ))
}
pr_df <- do.call(rbind, pr_df)
pr_lab <- do.call(rbind, pr_lab)
for (x in c("pr_df", "pr_lab")) {
  y <- get(x)
  y$panel <- factor(y$panel, levels = PPANELS)
  y$feature <- factor(y$feature, levels = PROMOTERS$label)
  assign(x, y)
}

set.seed(1)
p_pr <- strip_plot(pr_df, pr_lab, "Promoter usage") +
  facet_grid(feature ~ panel) +
  labs(title = "Promoters with concordant short-read effects",
       subtitle = "Short-read effects did not reach FDR < 0.05 (beta regression ~ condition + age + sex).") +
  theme(plot.title = element_text(size = 10, face = "bold"), plot.subtitle = element_text(size = 8))
ggsave(file.path(FIGS, "replicated_promoters.png"), p_pr, width = 5.5, height = 5.5, dpi = 300)
ggsave(file.path(FIGS, "replicated_promoters.pdf"), p_pr, width = 5.5, height = 5.5)

###############################
# Values behind the figures
###############################

write.table(rbind(cbind(type = "isoform", iso_df), cbind(type = "promoter", pr_df)),
            file.path(OUT, "tables/replicated_features_per_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(iso_lab)
print(pr_lab)
