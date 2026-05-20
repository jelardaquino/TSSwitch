library(tximport)
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(rtracklayer)
library(ballgown)
library(tidyverse)
library(IsoformSwitchAnalyzeR)
library(data.table)
library(biomaRt)
output_dir_tables <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/tables/suppa_beta_regression/"
dir.create(output_dir_tables, showWarnings = FALSE)

output_dir_figures <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/figures/suppa_beta_regression/"
dir.create(output_dir_figures, showWarnings = FALSE)

metadata <- read.table("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn25671149/jaquino_analysis/tmp/long_read_metadata.txt", header=TRUE, sep="\t")
metadata$sample <- gsub("^sample_", "", metadata$sample)
metadata$sample <- gsub("_PAM\\d+$|_PAG\\d+$", "", metadata$sample)
rownames(metadata) <- metadata$sample
conditions <- metadata$condition
unique(conditions)
grpA <- unique(conditions)[1]   # "Control"
grpB <- setdiff(unique(conditions), grpA)[1]  # "Case"

psi_files <- list.files(path = "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/suppa_bambu_gtf", pattern = "*.psi", full.names = TRUE)

psi_list <- lapply(psi_files, function(f) {
  df <- read.table(f, header = TRUE, sep = "\t", check.names = FALSE)
  return(df)
})

psi_all <- do.call(rbind, psi_list)

# 2. Filter: average PSI >= 0.05 in >= 75% samples
min_avg_psi <- 0.05
min_frac_samples <- 0.75
keep_event <- apply(psi_all, 1, function(v) {
  v_num <- as.numeric(v)
  non_na <- !is.na(v_num)
  if(sum(non_na)==0) return(FALSE)
  mean_psi_all <- mean(v_num[non_na], na.rm=TRUE)
  frac_ok <- mean(v_num[non_na] >= min_avg_psi, na.rm=TRUE)
  return(frac_ok >= min_frac_samples)
})
sum(keep_event)
psi_filt <- psi_all[keep_event, , drop=FALSE]

psi_filt <- psi_filt[, sample_info$sampleID]
A3_psi <- psi_filt[grepl("A3", rownames(psi_filt)), ]
A5_psi <- psi_filt[grepl("A5", rownames(psi_filt)), ]
SE_psi <- psi_filt[grepl("SE", rownames(psi_filt)), ]
RI_psi <- psi_filt[grepl("RI", rownames(psi_filt)), ]
MX_psi <- psi_filt[grepl("MX", rownames(psi_filt)), ]
AF_psi <- psi_filt[grepl("AF", rownames(psi_filt)), ]
AL_psi <- psi_filt[grepl("AL", rownames(psi_filt)), ]
# 3. compute per-event test and deltaPSI
library(glmmTMB)
library(broom.mixed)

sex      <- metadata$sex          # vector of sex for each sample
group    <- metadata$condition        # your existing group vector (grpA / grpB)

res_list <- lapply(rownames(psi_filt), function(e) {

  vals <- as.numeric(psi_filt[e, ])

  # build data frame for the model
  df <- data.frame(
    psi   = vals,
    group = relevel(factor(group), ref = "CT"),
    sex   = factor(sex)
  )

  # remove NAs
  df <- df[!is.na(df$psi), ]

  # need at least 2/group or test fails
  if (length(unique(df$group)) < 2 || any(table(df$group) < 2)) {
    return(list(event=e, pval=NA, delta=NA, n=nrow(df)))
  }

  # PSI must be strictly between 0 and 1 for beta regression
  df$psi[df$psi == 0] <- 1e-6
  df$psi[df$psi == 1] <- 1 - 1e-6

  # fit beta regression: psi ~ group + sex
  fit <- try(
    glmmTMB(psi ~ group + sex, data=df, family=beta_family()),
    silent=TRUE
  )

  if (inherits(fit, "try-error")) {
    return(list(event=e, pval=NA, delta=NA, n=nrow(df)))
  }

  # extract p-value for group
  coef_table <- broom.mixed::tidy(fit)
  pval <- coef_table$p.value[coef_table$term == paste0("group", grpB)]  # adjust name if needed

  # compute effect size (median delta PSI like before)
  delta <- median(df$psi[df$group == grpB]) - median(df$psi[df$group == grpA]) # positive delta means higher in case

  list(event=e, pval=pval, delta=delta, n=nrow(df))
})

clean_list <- Filter(
  function(x) length(x$pval) > 0,
  res_list
)
res_df <- do.call(rbind, lapply(clean_list, as.data.frame))
res_list_fixed <- lapply(res_list, function(x) {
  if (length(x$pval) == 0) x$pval <- NA
  x
})

res_df <- do.call(rbind, lapply(res_list_fixed, as.data.frame))
res_df$p_adj <- p.adjust(res_df$pval, method="BH")
res_df <- res_df %>% filter(!is.na(pval))
res_sig <- res_df %>%
  filter(!is.na(p_adj)) %>%
  filter(p_adj <= 0.05 & abs(delta) >= 0.2) %>%
  arrange(p_adj)

write.table(paste0(output_dir_tables, "/suppa_beta_regression_significant_results.tsv"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
save.image(paste0(output_dir_tables, "/suppa_beta_regression_analysis.RData"))

res_sig <- res_sig %>%
  mutate(
    geneID = str_extract(event, "ENSG\\d+"),
    event_type = str_extract(event, "(?<=;)\\w{2}")
  )

res_df <- res_df %>%
  mutate(
    geneID = str_extract(event, "ENSG\\d+"),
    event_type = str_extract(event, "(?<=;)\\w{2}")
  )
write.table(res_df, paste0(output_dir_tables, "/suppa_beta_regression_all_results.tsv"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

# prepare counts
event_counts <- as.data.frame(table(res_sig$event_type), stringsAsFactors = FALSE)
colnames(event_counts) <- c("event_type", "count")
event_counts <- event_counts[!is.na(event_counts$event_type), ]

# order factor by descending count
event_counts$event_type <- factor(event_counts$event_type,
    levels = event_counts$event_type[order(-event_counts$count)]
)

# plot
p <- ggplot(event_counts, aes(x = event_type, y = count, fill = event_type)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(label = count), vjust = -0.5, size = 3) +
    scale_fill_brewer(palette = "Set2") +
    labs(
        x = "Event type", y = "Number of significant events",
        title = "Significant SUPPA events by type"
    ) +
    theme_minimal()

print(p)

ggsave(paste0(output_dir_figures, "/res_sig_event_type_counts.png"), plot = p, width = 7, height = 4, dpi = 300)

merged_gtf <- import("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/NanoporeRNASeq_extended_annotations.gtf")
merged_gtf <- as.data.frame(merged_gtf)
AD_genes <- read.table("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/reproducing_RNA_figures/ad_genes.tsv", header=TRUE, stringsAsFactors=FALSE, sep="\t")

# Select Ensembl database
ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror = "useast"   # or "uswest", "asia"
)

# Get mapping of Ensembl gene IDs to gene names
gene_map <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                  mart = ensembl)
merged_gtf$gene_name <- gene_map$external_gene_name[match(merged_gtf$gene_id, gene_map$ensembl_gene_id)]
res_sig$gene_name <- gene_map$external_gene_name[match(res_sig$geneID, gene_map$ensembl_gene_id)]
res_df$gene_name <- gene_map$external_gene_name[match(res_df$geneID, gene_map$ensembl_gene_id)]
res_sig_AD <- res_sig %>% filter(gene_name %in% AD_genes$gene_name) 
# genes with differential AS: TPCN1, CTSB, MS4A4A, CLU, RBCK1, TNIP1, ADAM17, SCIMP, CD2AP, MME, MYO15A
write.table(res_sig, file=paste0(output_dir_tables, "/suppa_beta_regression_significant_results.tsv"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(res_sig_AD, file=paste0(output_dir_tables, "/suppa_beta_regression_significant_results_ADgenes.tsv"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

# did not work with the gene sets
# library(clusterProfiler)
# library(org.Hs.eg.db)

# mapped <- bitr(res_sig$gene_name,
#                fromType = "SYMBOL",
#                toType = "ENTREZID",
#                OrgDb = org.Hs.eg.db)

# background <- bitr(res_df$gene_name,
#                    fromType = "SYMBOL",
#                    toType = "ENTREZID",
#                    OrgDb = org.Hs.eg.db)

# ego <- enrichGO(
#   gene          = unique(mapped$ENTREZID),
#   universe      = unique(background$ENTREZID),
#   OrgDb         = org.Hs.eg.db,
#   keyType       = "ENTREZID",
#   ont           = "MF",
#   pAdjustMethod = "BH",
#   qvalueCutoff  = 0.05
# )

# png(paste0(output_dir_figures, "/GO_enrichment_SUPPA_significant_events_barplot.png"), width = 8, height = 6, units = "in", res = 300)
# dotplot(ego, showCategory = 20)
# dev.off()

# gene_df <- bitr(res_sig$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# ekegg <- enrichKEGG(
#   gene         = mapped$ENTREZID,
#   organism     = "hsa",          # human
#   pvalueCutoff = 0.05
# )
# # View results
# head(ekegg@result)
# # Visualize
# dotplot(ekegg, showCategory = 20, title = "KEGG Pathway Enrichment")

# originial_gene_list <- res_sig$delta 
# names(originial_gene_list) <- res_sig$gene_name
# # remove NA values
# gene_list <- na.omit(originial_gene_list)
# # sort in decreasing order
# gene_list <- sort(gene_list, decreasing = TRUE)
# png(paste0(output_dir_figures, "/GO_enrichment_SUPPA_significant_events_cnetplot.png"), width = 8, height = 6, units = "in", res = 300)
# p <- cnetplot(ego,
#               categorySize = "pvalue",
#               foldChange = gene_list,  # named numeric vector: names = genes
#               showCategory = 4)
# # Change legend title to "ΔPSI"
# p + guides(color = guide_colorbar(title = expression(Delta*"PSI")))
# dev.off()





