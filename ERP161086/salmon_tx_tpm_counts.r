###############################################################################
# salmon_tx_tpm_counts.r
#
# Purpose : Import Salmon transcript quantifications, build a sample-level
#           TPM matrix, and perform PCA to assess sample quality and grouping.
#
# Inputs  : - SRR_Acc_List.txt          : SRA run accession list
#           - salmon/<sample>_quant/    : Salmon quant.sf output directories
#           - metadata_AD_CTR.txt       : Sample metadata (Run, condition)
#
# Outputs : - transcript_tpm.txt        : Transcript-level TPM matrix
#           - PCA_TPM_Transcript_Counts.png : PCA plot colored by condition
#
# Next step: Run suppa_generate_events_sr.sh to compute PSI values,
#            then beta_regression_psi.r for differential splicing analysis.
###############################################################################

# ── Libraries ─────────────────────────────────────────────────────────────────
library(tximport)
library(readr)
library(dplyr)
library(ggplot2)
library(matrixStats)

# ── Path configuration ────────────────────────────────────────────────────────
# ⚠️ Update these paths to match your local environment before running.
BASE_DIR    <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086"
SALMON_DIR  <- file.path(BASE_DIR, "salmon")
OUTPUT_DIR  <- file.path(BASE_DIR, "results")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Load sample accession list and build file paths ───────────────────────────
samples <- read.table(
  file.path(BASE_DIR, "SRR_Acc_List.txt"),
  header = FALSE,
  stringsAsFactors = FALSE
)$V1

# Construct paths to each sample's quant.sf file
quant_files <- file.path(SALMON_DIR, paste0(samples, "_quant"), "quant.sf")
names(quant_files) <- samples

# Verify all quant files exist before importing
missing_files <- quant_files[!file.exists(quant_files)]
if (length(missing_files) > 0) {
  stop("Missing Salmon quant.sf files:\n", paste(missing_files, collapse = "\n"))
}

# ── Import transcript-level TPM values via tximport ───────────────────────────
# txOut = TRUE keeps results at transcript level (not summarized to gene)
txi <- tximport(quant_files, type = "salmon", txOut = TRUE)
tpm_matrix <- txi$abundance  # rows = transcripts, cols = samples

# ── Load metadata and align sample order ─────────────────────────────────────
metadata <- read.table(
  file.path(BASE_DIR, "metadata.txt"),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)
rownames(metadata) <- metadata$Run

sample_info <- data.frame(
  sampleID  = metadata$Run,
  group     = metadata$condition,
  stringsAsFactors = FALSE
)

# Reorder TPM columns to match sample_info row order
tpm_matrix <- tpm_matrix[, sample_info$sampleID]

# ── Write TPM matrix to disk ──────────────────────────────────────────────────
write.table(
  tpm_matrix,
  file      = file.path(OUTPUT_DIR, "transcript_tpm.txt"),
  sep       = "\t",
  quote     = FALSE,
  row.names = TRUE,
  col.names = TRUE
)
message("TPM matrix written to: ", file.path(OUTPUT_DIR, "transcript_tpm.txt"))

# ── PCA of transcript-level TPM ───────────────────────────────────────────────
# Step 1: log2(TPM + 1) transformation to stabilize variance
tpm_log <- log2(tpm_matrix + 1)

# Step 2: Remove the bottom 10% of transcripts by variance (low-information features)
row_vars    <- rowVars(as.matrix(tpm_log))
tpm_filtered <- tpm_log[row_vars > quantile(row_vars, 0.1), ]

# Step 3: Transpose so rows = samples (required by prcomp)
tpm_transposed <- t(tpm_filtered)

# Step 4: Run PCA with centering and scaling
pca_result <- prcomp(tpm_transposed, center = TRUE, scale. = TRUE)

# Step 5: Build data frame of PC coordinates with group labels
pca_df <- as.data.frame(pca_result$x)
pca_df$sample <- rownames(pca_df)
pca_df$group  <- metadata$condition[match(pca_df$sample, metadata$Run)]


# Calculate % variance explained per PC
var_explained <- round(100 * (pca_result$sdev^2 / sum(pca_result$sdev^2)), 1)

# Step 6: Plot PC1 vs PC2
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4) +
  xlab(paste0("PC1 (", var_explained[1], "% variance)")) +
  ylab(paste0("PC2 (", var_explained[2], "% variance)")) +
  theme_bw(base_size = 14) +
  ggtitle("PCA of TPM Transcript Counts — ERP161086")

ggsave(
  filename = file.path(OUTPUT_DIR, "PCA_TPM_Transcript_Counts.png"),
  plot     = pca_plot,
  width    = 6,
  height   = 4,
  dpi      = 300
)
message("PCA plot saved to: ", file.path(OUTPUT_DIR, "PCA_TPM_Transcript_Counts.png"))

# ── Session info for reproducibility ─────────────────────────────────────────
message("\n── Session Info ──")
sessionInfo()