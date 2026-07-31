###############################################################################
# beta_regression_psi.r
#
# Purpose : Test for differential alternative splicing between AD and CT
#           using beta regression (glmmTMB) on SUPPA2 PSI values.
#           Results are filtered by BH-adjusted
#           p-value (≤ 0.05) and effect size (|ΔPSI| ≥ 0.2).
#
# Inputs  : - suppa/psi_*.psi files                    : SUPPA2 PSI matrices
#           - metadata_AD_CTR.txt                       : Sample metadata
#           - Homo_sapiens.GRCh38.107_ERCC.gtf          : For gene name mapping
#           - ad_genes.tsv                              : Known AD risk gene list
#           - suppa_beta_regression_significant_results.tsv : Long-read results
#                                                         for cross-dataset comparison
#
# Outputs : - suppa_beta_regression_significant_results.tsv : All significant events
#           - suppa_beta_regression_significant_results_ADgenes.tsv : AD-gene subset
#           - PRJNA1206164_sig_events.tsv               : Annotated significant events
#           - matching_long_read_sig_events.txt          : Events shared with long-read
#           - suppa_beta_regression_analysis.RData       : Full workspace checkpoint
#           - res_sig_event_type_counts.png              : Bar plot of event type counts
#
# Next step: Run go_kegg_enrichment.r for pathway enrichment analysis.
###############################################################################

# ── Libraries ─────────────────────────────────────────────────────────────────
library(dplyr)
library(stringr)
library(ggplot2)
library(rtracklayer)
library(data.table)
library(biomaRt)
library(glmmTMB)
library(broom.mixed)

# ── Path configuration ────────────────────────────────────────────────────────
# ⚠️ Update these paths to match your local environment before running.
BASE_DIR       <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086"
SUPPA_DIR      <- file.path(BASE_DIR, "suppa")
GTF_FILE       <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/reproducing_nextflow_pipeline/Homo_sapiens.GRCh38.107_ERCC.gtf"
AD_GENES_FILE  <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/reproducing_RNA_figures/ad_genes.tsv"
LONG_READ_SIG  <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/tables/suppa_beta_regression/suppa_beta_regression_significant_results.tsv"
OUTPUT_DIR     <- file.path(BASE_DIR, "results")

# dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Load sample metadata ──────────────────────────────────────────────────────
metadata <- read.table(
  file.path(BASE_DIR, "metadata.txt"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)
rownames(metadata) <- metadata$Run

sample_info <- data.frame(
  sampleID = metadata$Run,
  group    = metadata$condition,
  sex      = factor(metadata$sex),
  stringsAsFactors = FALSE
)

# Define reference group for beta regression contrast (CT = reference)
grpA <- "AD"   # case group (positive delta = higher PSI in AD)
grpB <- "CT"   # reference group

# ── Load and combine all SUPPA2 PSI files ────────────────────────────────────
psi_files <- list.files(path = SUPPA_DIR, pattern = "\\.psi$", full.names = TRUE)

if (length(psi_files) == 0) {
  stop("No .psi files found in: ", SUPPA_DIR,
       "\nRun suppa_generate_events_sr.sh first.")
}

psi_list <- lapply(psi_files, function(f) {
  read.table(f, header = TRUE, sep = "\t", check.names = FALSE,
             stringsAsFactors = FALSE)
})
psi_all <- do.call(rbind, psi_list)

# ── Filter low-coverage events ────────────────────────────────────────────────
# Keep events where PSI >= 0.05 in at least 75% of samples
min_avg_psi      <- 0.05
min_frac_samples <- 0.75

keep_event <- apply(psi_all, 1, function(v) {
  v_num   <- suppressWarnings(as.numeric(v))
  non_na  <- !is.na(v_num)
  if (sum(non_na) == 0) return(FALSE)
  frac_ok <- mean(v_num[non_na] >= min_avg_psi, na.rm = TRUE)
  frac_ok >= min_frac_samples
})

psi_filt <- psi_all[keep_event, , drop = FALSE]
message("Events passing PSI filter: ", sum(keep_event), " of ", nrow(psi_all))

# Reorder columns to match sample_info order
psi_filt <- psi_filt[, sample_info$sampleID, drop = FALSE]

# ── Split PSI matrix by event type (for inspection/downstream use) ────────────
A3_psi <- psi_filt[grepl("A3", rownames(psi_filt)), ]
A5_psi <- psi_filt[grepl("A5", rownames(psi_filt)), ]
SE_psi <- psi_filt[grepl("SE", rownames(psi_filt)), ]
RI_psi <- psi_filt[grepl("RI", rownames(psi_filt)), ]
MX_psi <- psi_filt[grepl("MX", rownames(psi_filt)), ]
AF_psi <- psi_filt[grepl("AF", rownames(psi_filt)), ]
AL_psi <- psi_filt[grepl("AL", rownames(psi_filt)), ]

# ── Beta regression: PSI ~ group, per splicing event ───────────────────
# Uses glmmTMB with beta_family() which is appropriate for bounded [0,1] PSI values
group_vec <- sample_info$group
sex_vec   <- sample_info$sex

res_list <- lapply(rownames(psi_filt), function(e) {

  vals <- as.numeric(psi_filt[e, ])

  df <- data.frame(
    psi   = vals,
    group = relevel(factor(group_vec), ref = grpB), # CT as reference
    sex   = factor(sex_vec)
  )

  # Remove missing PSI values
  df <- df[!is.na(df$psi), ]

  # Require at least 2 samples per group for a valid test
  if (length(unique(df$group)) < 2 || any(table(df$group) < 2)) {
    return(list(event = e, pval = NA, delta = NA, n = nrow(df)))
  }

  # Beta regression requires PSI strictly between 0 and 1
  df$psi[df$psi == 0] <- 1e-6
  df$psi[df$psi == 1] <- 1 - 1e-6

  fit <- tryCatch(
    glmmTMB(psi ~ group + sex, data = df, family = beta_family()),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(list(event = e, pval = NA, delta = NA, n = nrow(df)))
  }

  # Extract p-value for the AD group coefficient
  coef_table <- broom.mixed::tidy(fit)
  group_term <- paste0("group", grpA)
  pval <- coef_table$p.value[coef_table$term == group_term]
  if (length(pval) == 0) pval <- NA

  # Effect size: median ΔPSI (AD minus CT); positive = higher PSI in AD
  delta <- median(df$psi[df$group == grpA], na.rm = TRUE) -
           median(df$psi[df$group == grpB], na.rm = TRUE)

  list(event = e, pval = pval, delta = delta, n = nrow(df))
})

# ── Collect results and apply BH multiple testing correction ─────────────────
# Safely handle any entries where pval length == 0
res_list_clean <- lapply(res_list, function(x) {
  if (length(x$pval) == 0) x$pval <- NA
  x
})

res_df       <- do.call(rbind, lapply(res_list_clean, as.data.frame))
res_df$p_adj <- p.adjust(res_df$pval, method = "BH")
res_df       <- res_df %>% filter(!is.na(pval))

# Significant events: BH-adjusted p ≤ 0.05 and |ΔPSI| ≥ 0.2
res_sig <- res_df %>%
  filter(!is.na(p_adj), p_adj <= 0.05, abs(delta) >= 0.2) %>%
  arrange(p_adj)

message("Significant events (p_adj ≤ 0.05, |ΔPSI| ≥ 0.2): ", nrow(res_sig))

# ── Annotate significant events with gene IDs and event types ─────────────────
res_sig <- res_sig %>%
  mutate(
    geneID     = str_extract(event, "ENSG\\d+"),
    event_type = str_extract(event, "(?<=;)\\w{2}")
  )

# ── Map Ensembl gene IDs to gene symbols via biomaRt ─────────────────────────
ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror  = "useast"  # alternatives: "uswest", "asia"
)

gene_map <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  mart       = ensembl
)

res_sig$gene_name <- gene_map$external_gene_name[
  match(res_sig$geneID, gene_map$ensembl_gene_id)
]

# ── Subset to known AD risk genes ─────────────────────────────────────────────
AD_genes   <- read.table(AD_GENES_FILE, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)
res_sig_AD <- res_sig %>% filter(gene_name %in% AD_genes$gene_name)
message("Significant events in known AD genes: ", nrow(res_sig_AD))

# ── Cross-dataset comparison: match with long-read significant events ─────────
long_read_sig <- read.table(LONG_READ_SIG, header = TRUE, sep = "\t",
                             stringsAsFactors = FALSE)
merge_genes <- merge(res_sig, long_read_sig, by = "geneID")
message("Events significant in both short-read and long-read: ", nrow(merge_genes))

# ── Write result tables ───────────────────────────────────────────────────────
write.table(res_sig,
  file = file.path(OUTPUT_DIR, "suppa_beta_regression_significant_results.tsv"),
  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(res_sig_AD,
  file = file.path(OUTPUT_DIR, "suppa_beta_regression_significant_results_ADgenes.tsv"),
  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(res_sig,
  file = file.path(OUTPUT_DIR, "ERP161086_sig_events.tsv"),
  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(merge_genes,
  file = file.path(OUTPUT_DIR, "matching_long_read_sig_events.txt"),
  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Save full workspace as checkpoint for downstream enrichment analysis
save.image(file.path(OUTPUT_DIR, "suppa_beta_regression_analysis.RData"))
message("Workspace checkpoint saved.")

# ── Bar plot: significant event counts by type ────────────────────────────────
event_counts <- as.data.frame(table(res_sig$event_type), stringsAsFactors = FALSE)
colnames(event_counts) <- c("event_type", "count")
event_counts <- event_counts[!is.na(event_counts$event_type), ]

# Order bars by descending count
event_counts$event_type <- factor(
  event_counts$event_type,
  levels = event_counts$event_type[order(-event_counts$count)]
)

event_type_plot <- ggplot(event_counts,
    aes(x = event_type, y = count, fill = event_type)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = count), vjust = -0.5, size = 3) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x     = "Event type",
    y     = "Number of significant events",
    title = "Significant SUPPA2 Events by Type — ERP161086"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(OUTPUT_DIR, "res_sig_event_type_counts.png"),
  plot     = event_type_plot,
  width    = 7, height = 4, dpi = 300
)
message("Event type bar plot saved.")

# ── Session info for reproducibility ─────────────────────────────────────────
message("\n── Session Info ──")
sessionInfo()
