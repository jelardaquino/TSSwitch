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

# Gather salmon quant.sf files
samples <- read.table("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/PRJNA1206164/SRR_Acc_List.txt", header = FALSE)$V1
files <- file.path(paste0("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/PRJNA1206164/salmon/", samples, "_quant"), "quant.sf")
names(files) <- samples

# Obtain transcript counts
txi <- tximport(files, type = "salmon", txOut = TRUE)
#counts <- txi$counts
counts <- txi$abundance # TPM values
#write.table(counts, "transcript_counts.txt", sep="\t", quote = FALSE, row.names = TRUE, col.names = TRUE)


metadata <- read.table("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/PRJNA1206164/metadata_AD_CTR.txt", header=TRUE, sep="\t")
rownames(metadata) <- metadata$Run
conditions <- metadata$diagnosis
unique(conditions)
grpA <- unique(conditions)[1]   # e.g., "Case"
grpB <- setdiff(unique(conditions), grpA)[1]  # e.g., "Control"

sample_info <- data.frame(sampleID = metadata$Run,
                          group = metadata$diagnosis,
                          sex = factor(metadata$sex))

counts <- counts[, sample_info$sampleID]
write.table(counts, "transcript_tpm.txt", sep="\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

library(DESeq2)

# dds: DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = sample_info,
                              design = ~1)  # no design needed for PCA

# # variance stabilizing transformation
# vsd <- vst(dds, blind = TRUE)

# # PCA plot
# png("PRJNA1206164_FC_PCA_plot.png", width = 1600, height = 1200, res = 300)
# plotPCA(vsd, intgroup = "group")
# dev.off()

# plotPCA(vsd, intgroup = "sex")

# 1. Log-transform the TPM values (recommended for RNA-seq data)
tpm_log <- log2(counts + 1)  # +1 to avoid log(0)

# Remove transcripts with very low variance across samples
row_vars <- rowVars(as.matrix(tpm_log))
keep <- row_vars > quantile(row_vars, 0.1)  # keep top 90%
tpm_filtered <- tpm_log[keep, ]

# Then proceed with t() and prcomp() as above
# 2. Transpose so samples are rows (required by prcomp)
tpm_t <- t(tpm_filtered)

# 3. Run PCA (center and scale the data)
pca_result <- prcomp(tpm_t, center = TRUE, scale. = TRUE)

# 4. View summary (variance explained by each PC)
summary(pca_result)

# Extract PCA coordinates
pca_df <- as.data.frame(pca_result$x)

# Add sample names
pca_df$sample <- rownames(pca_df)

# Calculate % variance explained
var_explained <- round(
  100 * (pca_result$sdev^2 / sum(pca_result$sdev^2)), 1
)

# Add a grouping column (e.g., treatment vs control)
pca_df$group <- metadata$diagnosis[match(pca_df$sample, metadata$Run)]
# Replace label "CTR" with "CT" in the PCA grouping column
pca_df$group <- ifelse(pca_df$group == "CTR", "CT", pca_df$group)
ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4) +
  xlab(paste0("PC1 (", var_explained[1], "% variance)")) +
  ylab(paste0("PC2 (", var_explained[2], "% variance)")) +
  theme_bw(base_size=14) +
  ggtitle("PCA of TPM Transcript Counts")
ggsave("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/PRJNA1206164/PCA_TPM_Transcript_Counts.png", width = 6, height = 4, dpi = 300)

suppa.py generateEvents -i /home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/reproducing_nextflow_pipeline/Homo_sapiens.GRCh38.107_ERCC.gtf -o FC_events -f ioe -e SE MX SS RI FL
for event in SE RI MX AL AF A5 A3
do
    suppa.py psiPerEvent \
        --ioe-file FC_events_${event}_strict.ioe \
        --expression-file ../transcript_tpm.txt \
        -o psi_${event}
done


psi_files <- list.files(path = ".", pattern = "*.psi", full.names = TRUE)

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

sex      <- sample_info$sex          # vector of sex for each sample
group    <- sample_info$group        # your existing group vector (grpA / grpB)

res_list <- lapply(rownames(psi_filt), function(e) {

  vals <- as.numeric(psi_filt[e, ])

  # build data frame for the model
  df <- data.frame(
    psi   = vals,
    group = relevel(factor(group), ref = "CTR"),
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
  pval <- coef_table$p.value[coef_table$term == paste0("group", grpA)]  # adjust name if needed

  # compute effect size (median delta PSI like before)
  delta <- median(df$psi[df$group == grpA]) - median(df$psi[df$group == grpB]) # positive delta means higher in case

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

write.table(res_sig, file="suppa_beta_regression_significant_results.tsv", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
save.image("suppa_beta_regression_analysis.RData")

res_sig <- res_sig %>%
  mutate(
    geneID = str_extract(event, "ENSG\\d+"),
    event_type = str_extract(event, "(?<=;)\\w{2}")
  )

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

ggsave("res_sig_event_type_counts.png", plot = p, width = 7, height = 4, dpi = 300)

gtf <- import("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/reproducing_nextflow_pipeline/Homo_sapiens.GRCh38.107_ERCC.gtf")
gtf_df <- as.data.frame(gtf)
geneID_genename <- gtf_df[,c("gene_id", "gene_name")]
geneID_genename <- geneID_genename %>% distinct(gene_id, .keep_all = TRUE)
res_sig <- res_sig %>%
  left_join(geneID_genename, by = c("geneID" = "gene_id"))

long_read_sig <- read.table("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn25671149/jaquino_analysis/tmp/suppa/suppa_beta_regression_significant_results.tsv", header=TRUE, sep="\t")

merge_genes <- merge(res_sig, long_read_sig, by="geneID")

write.table(res_sig, "PRJNA1206164_sig_events.tsv", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(merge_genes, "matching_long_read_sig_events.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


merged_gtf <- import("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn25671149/jaquino_analysis/tmp/NanoporeRNASeq_extended_annotations.gtf")
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
res_sig_AD <- res_sig %>% filter(gene_name %in% AD_genes$gene_name) 
# genes with differential AS: TPCN1, CTSB, MS4A4A, CLU, RBCK1, TNIP1, ADAM17, SCIMP, CD2AP, MME, MYO15A
write.table(res_sig, file="suppa_beta_regression_significant_results.tsv", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
write.table(res_sig_AD, file="suppa_beta_regression_significant_results_ADgenes.tsv", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

library(clusterProfiler)
library(org.Hs.eg.db)

ego <- enrichGO(
  gene          = res_sig$gene_name,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)
png("GO_enrichment_SUPPA_significant_events_dotplot.png", width = 8, height = 12, units = "in", res = 300)
dotplot(ego, showCategory = 20)
dev.off()

gene_df <- bitr(res_sig$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ekegg <- enrichKEGG(
  gene         = gene_df$ENTREZID,
  organism     = "hsa",          # human
  pvalueCutoff = 0.05
)
# View results
head(ekegg@result)
# Visualize
png("KEGG_SUPPA_significant_events_dotplot.png", width = 6, height = 4, units = "in", res = 300)
dotplot(ekegg, showCategory = 20, title = "KEGG Pathway Enrichment")
dev.off()

originial_gene_list <- res_sig$delta 
names(originial_gene_list) <- res_sig$gene_name
# remove NA values
gene_list <- na.omit(originial_gene_list)
# sort in decreasing order
gene_list <- sort(gene_list, decreasing = TRUE)
png("GO_enrichment_SUPPA_significant_events_cnetplot.png", width = 8, height = 6, units = "in", res = 300)
p <- cnetplot(ego,
              categorySize = "pvalue",
              foldChange = gene_list,  # named numeric vector: names = genes
              showCategory = 4)
# Change legend title to "ΔPSI"
p + guides(color = guide_colorbar(title = expression(Delta*"PSI")))
dev.off()