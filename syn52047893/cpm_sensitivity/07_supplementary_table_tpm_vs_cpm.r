#!/usr/bin/env Rscript

# Step 7: one supplementary table comparing the TPM (manuscript) and CPM analyses.
# One row per feature significant in either analysis, with the significance statistic and
# effect size from both where the analysis produced one.

suppressPackageStartupMessages({
  library(dplyr)
})

J <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis"
V <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis_v2"
C <- file.path(J, "cpm_pipeline")
OUT_FILE <- file.path(C, "comparison", "Supplementary_Table_X_TPM_vs_CPM_sensitivity.csv")

# PSI tables contain plot labels with embedded line breaks; rejoin physical lines.
read_tsv_multiline <- function(path) {
  lines <- readLines(path, warn = FALSE)
  header <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  n <- length(header)
  rows <- list()
  buf <- NULL
  for (ln in lines[-1]) {
    buf <- if (is.null(buf)) ln else paste(buf, ln)
    if (lengths(regmatches(buf, gregexpr("\t", buf))) + 1 >= n) {
      rows[[length(rows) + 1]] <- strsplit(buf, "\t", fixed = TRUE)[[1]][seq_len(n)]
      buf <- NULL
    }
  }
  setNames(as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE), header)
}

num <- function(x) suppressWarnings(as.numeric(x))

concordance <- function(tpm, cpm) {
  case_when(tpm & cpm ~ "Both", tpm ~ "TPM only", cpm ~ "CPM only")
}

empty_cols <- function(df) {
  for (col in c("feature_detail", "statistic_TPM", "statistic_CPM", "effect_TPM", "effect_CPM")) {
    if (!col %in% names(df)) df[[col]] <- NA
  }
  df
}

###############################
# DIU genes and isoforms
###############################

diu_tpm <- read.delim(file.path(J, "diu/IsoformSwitchAnalyzeR_top_switches_gene_name.tsv"), stringsAsFactors = FALSE)
diu_cpm <- read.delim(file.path(C, "diu/IsoformSwitchAnalyzeR_top_switches_gene_name.tsv"), stringsAsFactors = FALSE)

diu_genes <- full_join(
  diu_tpm %>% transmute(gene_id, gene_TPM = gene_name, statistic_TPM = gene_switch_q_value),
  diu_cpm %>% transmute(gene_id, gene_CPM = gene_name, statistic_CPM = gene_switch_q_value),
  by = "gene_id"
) %>%
  transmute(
    analysis = "DIU gene",
    feature_id = gene_id,
    gene_symbol = coalesce(gene_TPM, gene_CPM),
    significant_TPM = !is.na(statistic_TPM),
    significant_CPM = !is.na(statistic_CPM),
    statistic_name = "gene switch q-value",
    statistic_TPM, statistic_CPM,
    effect_name = NA_character_
  )

iso_features <- function(rds) {
  readRDS(rds)$isoformFeatures %>%
    distinct(isoform_id, .keep_all = TRUE) %>%
    transmute(isoform_id, gene_id, gene_name, q = isoform_switch_q_value, dIF,
              sig = !is.na(q) & q < 0.05 & abs(dIF) > 0.1)
}
iso_tpm <- iso_features(file.path(J, "diu/switchList_with_gene_name.rds"))
iso_cpm <- iso_features(file.path(C, "diu/switchList_with_gene_name.rds"))

diu_iso <- full_join(iso_tpm, iso_cpm, by = "isoform_id", suffix = c("_TPM", "_CPM")) %>%
  filter(sig_TPM %in% TRUE | sig_CPM %in% TRUE) %>%
  transmute(
    analysis = "DIU switching isoform",
    feature_id = isoform_id,
    gene_symbol = coalesce(gene_name_TPM, gene_name_CPM),
    feature_detail = coalesce(gene_id_TPM, gene_id_CPM),
    significant_TPM = sig_TPM %in% TRUE,
    significant_CPM = sig_CPM %in% TRUE,
    statistic_name = "isoform switch q-value",
    statistic_TPM = q_TPM, statistic_CPM = q_CPM,
    effect_name = "dIF (AD - CT)",
    effect_TPM = dIF_TPM, effect_CPM = dIF_CPM
  )

###############################
# PSI events
###############################

psi_rows <- function(label, tpm_dir, tpm_file, cpm_dir, cpm_file) {
  load <- function(path) {
    read_tsv_multiline(path) %>%
      transmute(event_id, gene_symbol, event_type,
                adj_p = num(adj_p_value), dPSI = num(dPSI),
                sig = !is.na(adj_p) & adj_p < 0.05 & !is.na(dPSI) & abs(dPSI) >= 0.10)
  }
  tpm <- load(file.path(tpm_dir, tpm_file))
  cpm <- load(file.path(cpm_dir, cpm_file))
  full_join(tpm, cpm, by = "event_id", suffix = c("_TPM", "_CPM")) %>%
    filter(sig_TPM %in% TRUE | sig_CPM %in% TRUE) %>%
    transmute(
      analysis = label,
      feature_id = event_id,
      gene_symbol = coalesce(gene_symbol_TPM, gene_symbol_CPM),
      feature_detail = coalesce(event_type_TPM, event_type_CPM),
      significant_TPM = sig_TPM %in% TRUE,
      significant_CPM = sig_CPM %in% TRUE,
      statistic_name = "BH-adjusted p (LRT)",
      statistic_TPM = adj_p_TPM, statistic_CPM = adj_p_CPM,
      effect_name = "dPSI (AD - CT)",
      effect_TPM = dPSI_TPM, effect_CPM = dPSI_CPM
    )
}

psi_additive <- psi_rows("PSI event, additive model",
                         file.path(J, "tables/suppa_beta_regression"), "all_beta_results.tsv",
                         file.path(C, "psi_additive/tables"), "all_beta_results.tsv")
psi_sex <- psi_rows("PSI event, condition x sex model",
                    file.path(J, "tables/suppa_beta_regression_sex"), "all_beta_results_sex.tsv",
                    file.path(C, "psi_condition_x_sex/tables"), "all_beta_results_sex.tsv")

###############################
# Promoters and promoter-switching genes
###############################

prom_load <- function(path) {
  read.delim(path, stringsAsFactors = FALSE) %>%
    transmute(promoter_id, gene_name, padj, logitFC, sig = !is.na(padj) & padj <= 0.05)
}
prom_tpm <- prom_load(file.path(V, "tables/promoter_usage_86_genes_12samples/promoter_usage_betareg_results_9samples.tsv"))
prom_cpm <- prom_load(file.path(C, "promoter_usage/tables/promoter_usage_betareg_results_9samples.tsv"))

promoters <- full_join(prom_tpm, prom_cpm, by = "promoter_id", suffix = c("_TPM", "_CPM")) %>%
  filter(sig_TPM %in% TRUE | sig_CPM %in% TRUE) %>%
  transmute(
    analysis = "Differentially used promoter",
    feature_id = promoter_id,
    gene_symbol = coalesce(gene_name_TPM, gene_name_CPM),
    significant_TPM = sig_TPM %in% TRUE,
    significant_CPM = sig_CPM %in% TRUE,
    statistic_name = "BH-adjusted p (beta regression)",
    statistic_TPM = padj_TPM, statistic_CPM = padj_CPM,
    effect_name = "logit fold change (AD vs CT)",
    effect_TPM = logitFC_TPM, effect_CPM = logitFC_CPM
  )

sw_tpm <- unique(read.delim(file.path(V, "tables/promoter_usage_86_genes_12samples/switching_promoter_usage_9samples.tsv"))$gene_name)
sw_cpm <- unique(read.delim(file.path(C, "promoter_usage/tables/switching_promoter_usage_9samples.tsv"))$gene_name)
switching <- data.frame(gene_symbol = union(sw_tpm, sw_cpm), stringsAsFactors = FALSE) %>%
  transmute(
    analysis = "Promoter-switching gene",
    feature_id = gene_symbol,
    gene_symbol,
    significant_TPM = gene_symbol %in% sw_tpm,
    significant_CPM = gene_symbol %in% sw_cpm,
    statistic_name = NA_character_,
    effect_name = NA_character_
  )

###############################
# Motif enrichment
###############################

motif_load <- function(dir) {
  best <- function(file, direction) {
    read.csv(file.path(dir, file), stringsAsFactors = FALSE) %>%
      filter(!is.na(target)) %>%
      arrange(p.value) %>%
      distinct(target, .keep_all = TRUE) %>%
      mutate(rank = row_number(), direction = direction)
  }
  bind_rows(
    best("motif_enrichment_promoters_up_AD_2kbps_v4.csv", "Increased usage in AD"),
    best("motif_enrichment_promoters_down_AD_2kbps_v4.csv", "Decreased usage in AD")
  ) %>% transmute(target, direction, p = p.value, rank, sig = p <= 0.05)
}
mot_tpm <- motif_load(file.path(V, "tables/promoter_usage_86_genes_12samples"))
mot_cpm <- motif_load(file.path(C, "promoter_usage/tables"))

motifs <- full_join(mot_tpm, mot_cpm, by = c("target", "direction"), suffix = c("_TPM", "_CPM")) %>%
  filter(sig_TPM %in% TRUE | sig_CPM %in% TRUE) %>%
  transmute(
    analysis = "Enriched TF motif",
    feature_id = target,
    gene_symbol = NA_character_,
    feature_detail = direction,
    significant_TPM = sig_TPM %in% TRUE,
    significant_CPM = sig_CPM %in% TRUE,
    statistic_name = "nominal p (PWMEnrich)",
    statistic_TPM = p_TPM, statistic_CPM = p_CPM,
    effect_name = "rank among TF targets",
    effect_TPM = rank_TPM, effect_CPM = rank_CPM
  )

###############################
# Combine and write
###############################

analysis_order <- c("DIU gene", "DIU switching isoform", "PSI event, additive model",
                    "PSI event, condition x sex model", "Differentially used promoter",
                    "Promoter-switching gene", "Enriched TF motif")

combined <- list(diu_genes, diu_iso, psi_additive, psi_sex, promoters, switching, motifs) %>%
  lapply(empty_cols) %>%
  bind_rows() %>%
  mutate(
    concordance = concordance(significant_TPM, significant_CPM),
    across(c(statistic_TPM, statistic_CPM, effect_TPM, effect_CPM), ~ signif(as.numeric(.x), 4))
  ) %>%
  select(analysis, feature_id, gene_symbol, feature_detail, concordance,
         significant_TPM, significant_CPM, statistic_name, statistic_TPM, statistic_CPM,
         effect_name, effect_TPM, effect_CPM) %>%
  arrange(factor(analysis, levels = analysis_order),
          factor(concordance, levels = c("Both", "TPM only", "CPM only")),
          coalesce(statistic_TPM, statistic_CPM))

title <- paste(
  "Supplementary Table X - Sensitivity of results to transcript abundance measure (TPM, primary analysis; CPM, sensitivity analysis).",
  "Rows are features significant in at least one analysis. The CPM analysis used the same surrogate variables as the TPM analysis.",
  "Gene switch q-values are reported only for genes significant in that analysis."
)
writeLines(title, OUT_FILE)
suppressWarnings(write.table(combined, OUT_FILE, sep = ",", row.names = FALSE, na = "NA",
                             append = TRUE, col.names = TRUE, quote = TRUE))

summary_counts <- combined %>%
  count(analysis = factor(analysis, levels = analysis_order), concordance) %>%
  tidyr::pivot_wider(names_from = concordance, values_from = n, values_fill = 0)
print(as.data.frame(summary_counts))
message("Wrote ", OUT_FILE, " (", nrow(combined), " rows)")
