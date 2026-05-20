library(ggplot2)
library(matrixStats)
# PCA of TPM counts

out_dir_figures <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/figures/"
dir.create(out_dir_figures, showWarnings = FALSE)
out_dir_tables <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/tables/"
dir.create(out_dir_tables, showWarnings = FALSE)

# Load metadata and TPM counts
metadata <- read.table("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn25671149/jaquino_analysis/tmp/long_read_metadata.txt", header=TRUE, sep="\t")
metadata$sample <- gsub("^sample_", "X", metadata$sample)
metadata$sample <- gsub("_PAM\\d+$|_PAG\\d+$", "", metadata$sample)
rownames(metadata) <- metadata$sample
conditions <- metadata$condition
unique(conditions)
grpA <- unique(conditions)[1]   # "Control"
grpB <- setdiff(unique(conditions), grpA)[1]  # "Case"

# tpm counts where rows = transcripts/genes and columns = samples
counts <- read.table("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/bambu_tx_tpm.txt", header=TRUE, sep="\t", row.names=1)

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

# Plot PC1 vs PC2
ggplot(pca_df, aes(x = PC1, y = PC2, label = sample)) +
  geom_point(size = 10) +
  geom_text(vjust = -0.5, size = 10) +
  xlab(paste0("PC1 (", var_explained[1], "% variance)")) +
  ylab(paste0("PC2 (", var_explained[2], "% variance)")) +
  theme_bw() +
  ggtitle("PCA of TPM Transcript Counts")

# Add a grouping column (e.g., treatment vs control)
pca_df$group <- metadata$condition[match(pca_df$sample, metadata$sample)]

ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4) +
  xlab(paste0("PC1 (", var_explained[1], "% variance)")) +
  ylab(paste0("PC2 (", var_explained[2], "% variance)")) +
  theme_bw(base_size=14) +
  ggtitle("PCA of TPM Transcript Counts")
ggsave(paste0(out_dir_figures, "PCA_TPM_Transcript_Counts.png"), width = 6, height = 4, dpi = 300)
write.table(pca_df, paste0(out_dir_tables, "PCA_TPM_Transcript_Counts.txt"), sep="\t", quote=FALSE, row.names=FALSE)