#!/usr/bin/env Rscript

# Step 6: compare the CPM pipeline results with the manuscript (TPM) results for
# DIU, PSI (additive and condition x sex), promoter usage and motif enrichment.

suppressPackageStartupMessages({
  library(dplyr)
})

J <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis"
V <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis_v2"
C <- file.path(J, "cpm_pipeline")
OUT <- file.path(C, "comparison")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

NAMED_GENES <- c("CLSTN3", "SNCA", "BECN1", "USP10", "EGFR", "DMXL1", "NSD2", "NOTCH4",
                 "SPP1", "PSMC6", "PVT1", "COL4A2", "PLSCR3", "KDSR", "CXXC1", "C11orf80")
NAMED_MOTIFS <- list(up = c("IRF9", "IRF8", "IRF7", "STAT2::STAT1"), down = c("NFKB2", "NFKB1"))

# Some PSI tables contain plot labels with embedded line breaks; rejoin physical lines
# until a record has the full number of fields.
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
  if (length(rows) == 0) return(setNames(data.frame(matrix(character(), 0, n)), header))
  setNames(as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE), header)
}

overlap_table <- function(tpm, cpm, key) {
  ids <- union(tpm[[key]], cpm[[key]])
  data.frame(
    id = ids,
    in_TPM = ids %in% tpm[[key]],
    in_CPM = ids %in% cpm[[key]],
    stringsAsFactors = FALSE
  ) %>%
    mutate(status = case_when(in_TPM & in_CPM ~ "Both", in_TPM ~ "TPM only", TRUE ~ "CPM only"))
}

line_counts <- function(label, ov) {
  sprintf("%-38s TPM %4d | CPM %4d | shared %4d", label,
          sum(ov$in_TPM), sum(ov$in_CPM), sum(ov$in_TPM & ov$in_CPM))
}

report <- c(paste("TPM vs CPM comparison,", format(Sys.time(), "%Y-%m-%d %H:%M")), "")

###############################
# DIU
###############################

diu_tpm <- read.delim(file.path(J, "diu/IsoformSwitchAnalyzeR_top_switches_gene_name.tsv"), stringsAsFactors = FALSE)
diu_cpm <- read.delim(file.path(C, "diu/IsoformSwitchAnalyzeR_top_switches_gene_name.tsv"), stringsAsFactors = FALSE)

sig_isoforms <- function(rds) {
  f <- readRDS(rds)$isoformFeatures
  unique(f$isoform_id[!is.na(f$isoform_switch_q_value) & f$isoform_switch_q_value < 0.05 & abs(f$dIF) > 0.1])
}
iso_tpm <- data.frame(isoform_id = sig_isoforms(file.path(J, "diu/switchList_with_gene_name.rds")))
iso_cpm <- data.frame(isoform_id = sig_isoforms(file.path(C, "diu/switchList_with_gene_name.rds")))

ov_diu <- overlap_table(diu_tpm, diu_cpm, "gene_id")
gene_names <- bind_rows(diu_tpm[, c("gene_id", "gene_name")], diu_cpm[, c("gene_id", "gene_name")]) %>% distinct(gene_id, .keep_all = TRUE)
ov_diu$gene_name <- gene_names$gene_name[match(ov_diu$id, gene_names$gene_id)]
write.table(ov_diu, file.path(OUT, "diu_genes_tpm_vs_cpm.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
ov_iso <- overlap_table(iso_tpm, iso_cpm, "isoform_id")
write.table(ov_iso, file.path(OUT, "diu_isoforms_tpm_vs_cpm.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

diag_sets <- c(
  "CPM, covariates re-estimated from CPM" = file.path(C, "archive_covariates_estimated_from_cpm/diu/diu_sig_genes.txt"),
  "TPM, no covariates" = file.path(C, "diagnostics/tpm_noSV/sig_genes.txt"),
  "CPM, no covariates" = file.path(C, "diagnostics/cpm_noSV/sig_genes.txt")
)
diag_lines <- vapply(names(diag_sets), function(nm) {
  g <- if (file.exists(diag_sets[[nm]])) unique(readLines(diag_sets[[nm]])) else character()
  sprintf("  %-38s %5d genes | shared with TPM manuscript: %d", nm, length(g), sum(g %in% diu_tpm$gene_id))
}, character(1))

report <- c(report, "DIU sensitivity to the surrogate-variable (batch covariate) step",
            unname(diag_lines), "")

report <- c(report, "DIU (IsoformSwitchAnalyzeR, DEXSeq; CPM with the manuscript's surrogate variables)",
            line_counts("  genes", ov_diu),
            line_counts("  switching isoforms", ov_iso),
            paste("  named genes significant with CPM:",
                  paste(sprintf("%s=%s", NAMED_GENES, ifelse(NAMED_GENES %in% diu_cpm$gene_name, "yes", "no")), collapse = ", ")),
            paste(" ", readLines(file.path(C, "diu/nmd_summary.txt"))), "")

###############################
# PSI
###############################

psi_sets <- list(
  additive = c(file.path(J, "tables/suppa_beta_regression/all_sig_beta_events.tsv"),
               file.path(C, "psi_additive/tables/all_sig_beta_events.tsv")),
  condition_x_sex = c(file.path(J, "tables/suppa_beta_regression_sex/all_sig_beta_events_sex.tsv"),
                      file.path(C, "psi_condition_x_sex/tables/all_sig_beta_events_sex.tsv"))
)

for (m in names(psi_sets)) {
  tpm <- read_tsv_multiline(psi_sets[[m]][1])
  cpm <- read_tsv_multiline(psi_sets[[m]][2])
  ov <- overlap_table(tpm, cpm, "event_id")
  sym <- bind_rows(tpm[, c("event_id", "gene_symbol", "event_type")], cpm[, c("event_id", "gene_symbol", "event_type")]) %>%
    distinct(event_id, .keep_all = TRUE)
  ov$gene_symbol <- sym$gene_symbol[match(ov$id, sym$event_id)]
  ov$event_type <- sym$event_type[match(ov$id, sym$event_id)]
  write.table(ov, file.path(OUT, paste0("psi_events_", m, "_tpm_vs_cpm.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  report <- c(report, paste0("PSI, ", m, " model"),
              line_counts("  significant events", ov),
              paste("  genes (CPM):", paste(sort(unique(cpm$gene_symbol)), collapse = ", ")), "")
}

###############################
# Promoter usage
###############################

prom_tpm <- read.delim(file.path(V, "tables/promoter_usage_86_genes_12samples/beta_reg_promoter_sig.tsv"), stringsAsFactors = FALSE)
prom_cpm <- read.delim(file.path(C, "promoter_usage/tables/beta_reg_promoter_sig.tsv"), stringsAsFactors = FALSE)
sw_tpm <- read.delim(file.path(V, "tables/promoter_usage_86_genes_12samples/switching_promoter_usage_9samples.tsv"), stringsAsFactors = FALSE)
sw_cpm <- read.delim(file.path(C, "promoter_usage/tables/switching_promoter_usage_9samples.tsv"), stringsAsFactors = FALSE)

ov_prom <- overlap_table(prom_tpm, prom_cpm, "promoter_id")
write.table(ov_prom, file.path(OUT, "promoters_tpm_vs_cpm.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
ov_sw <- overlap_table(distinct(sw_tpm, gene_name), distinct(sw_cpm, gene_name), "gene_name")
write.table(ov_sw, file.path(OUT, "promoter_switching_genes_tpm_vs_cpm.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

report <- c(report, "Promoter usage (beta regression, padj <= 0.05)",
            line_counts("  differentially used promoters", ov_prom),
            sprintf("  CPM: %d increased / %d decreased in AD", sum(prom_cpm$logitFC > 0), sum(prom_cpm$logitFC < 0)),
            line_counts("  promoter-switching genes", ov_sw),
            paste("  switching genes (CPM):", paste(sort(unique(sw_cpm$gene_name)), collapse = ", ")), "")

###############################
# Motif enrichment
###############################

motif_files <- list(
  TPM = file.path(V, "tables/promoter_usage_86_genes_12samples"),
  CPM = file.path(C, "promoter_usage/tables")
)

read_motifs <- function(dir, label) {
  up <- read.csv(file.path(dir, "motif_enrichment_promoters_up_AD_2kbps_v4.csv"), stringsAsFactors = FALSE)
  dn <- read.csv(file.path(dir, "motif_enrichment_promoters_down_AD_2kbps_v4.csv"), stringsAsFactors = FALSE)
  rank_best <- function(d) {
    d %>% arrange(p.value) %>% distinct(target, .keep_all = TRUE) %>%
      mutate(target_rank = row_number()) %>% select(target, p.value, target_rank)
  }
  bind_rows(
    rank_best(up) %>% mutate(direction = "Up in AD"),
    rank_best(dn) %>% mutate(direction = "Down in AD")
  ) %>% mutate(source = label)
}

mot <- bind_rows(read_motifs(motif_files$TPM, "TPM"), read_motifs(motif_files$CPM, "CPM"))
mot_wide <- mot %>%
  tidyr::pivot_wider(id_cols = c(target, direction), names_from = source,
                     values_from = c(p.value, target_rank))
write.table(mot_wide, file.path(OUT, "motifs_tpm_vs_cpm.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

sig_n <- function(src, dir) sum(mot$source == src & mot$direction == dir & mot$p.value <= 0.05)
top10 <- function(src, dir) paste(head(mot$target[mot$source == src & mot$direction == dir], 10), collapse = ", ")
named_rank <- function(dir, targets) {
  paste(vapply(targets, function(t) {
    r <- mot_wide[mot_wide$target %in% t & mot_wide$direction %in% dir, ]
    if (nrow(r) == 0) return(sprintf("%s: not found", t))
    sprintf("%s rank TPM %s / CPM %s (p %.2g / %.2g)", t, r$target_rank_TPM, r$target_rank_CPM, r$p.value_TPM, r$p.value_CPM)
  }, character(1)), collapse = "; ")
}

sig_table_counts <- function(dir) {
  s <- read.delim(file.path(dir, "motif_enrichment_promoters_sig_AD_2kbps.tsv"), stringsAsFactors = FALSE)
  sprintf("%d motifs, %d unique TF targets", nrow(s), length(unique(s$target)))
}

report <- c(report, "Motif enrichment (PWMEnrich, nominal p <= 0.05 per TF target)",
            paste("  Significant motif table, TPM:", sig_table_counts(motif_files$TPM)),
            paste("  Significant motif table, CPM:", sig_table_counts(motif_files$CPM)),
            sprintf("  Up in AD:   TPM %d targets | CPM %d targets", sig_n("TPM", "Up in AD"), sig_n("CPM", "Up in AD")),
            sprintf("  Down in AD: TPM %d targets | CPM %d targets", sig_n("TPM", "Down in AD"), sig_n("CPM", "Down in AD")),
            paste("  Top 10 Up in AD (TPM):", top10("TPM", "Up in AD")),
            paste("  Top 10 Up in AD (CPM):", top10("CPM", "Up in AD")),
            paste("  Top 10 Down in AD (TPM):", top10("TPM", "Down in AD")),
            paste("  Top 10 Down in AD (CPM):", top10("CPM", "Down in AD")),
            paste("  Named, Up in AD:", named_rank("Up in AD", NAMED_MOTIFS$up)),
            paste("  Named, Down in AD:", named_rank("Down in AD", NAMED_MOTIFS$down)))

writeLines(report, file.path(OUT, "summary_tpm_vs_cpm.txt"))
cat(report, sep = "\n")
