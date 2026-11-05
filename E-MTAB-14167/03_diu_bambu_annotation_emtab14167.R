#!/usr/bin/env Rscript

# DIU in the full E-MTAB-14167 prefrontal cortex cohort (20 AD vs 21 CT) using the
# long-read Bambu extended annotation, then comparison with the long-read DIU results.
# Replaces the 7 vs 7 ERP161086 run (ERP161086/scripts/isoformswitchanalyzer_bambu_annotation.r),
# whose 14 donors are a subset of this cohort.
#
# Import, preFilter() and DEXSeq settings match the long-read DIU
# (syn52065646_LongRead/jaquino_analysis/scripts/isoformswitchanalyzer_DEXSeq.r); age and sex
# are added to the design because the cohort is unbalanced for age. Counts and TPM are read
# from the Salmon quant.sf files made by 01_salmon_bambu_quant_emtab14167.sh.

suppressPackageStartupMessages({
  library(IsoformSwitchAnalyzeR)
  library(data.table)
})

set.seed(1)

S <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing"
QUANT_DIR <- file.path(S, "ERP161086/emtab14167_replication/salmon_bambu/quant")
GTF_FILE <- file.path(S, "syn52065646_LongRead/jaquino_analysis/NanoporeRNASeq_extended_annotations.gtf")
SR_META <- file.path(S, "E-MTAB-14167/Step01_metadata/sample_metadata_matchedv2.tsv")
GENE_MAP_FILE <- file.path(S, "syn52065646_LongRead/pychopperfq/minimap2/bam_mapq10/all_bambu_Rscript/bambu_tx_tpm_GENEID_genesymbol.txt")

LR_DIU_RDS <- file.path(S, "syn52065646_LongRead/jaquino_analysis/diu/switchList_with_gene_name.rds")
LR_DIU_TOP <- file.path(S, "syn52065646_LongRead/jaquino_analysis/diu/IsoformSwitchAnalyzeR_top_switches_gene_name.tsv")
# The long-read DEXSeq object was reduced to switching genes, so the set of genes tested in
# the long-read data is taken from the satuRn rerun (same import and preFilter()).
LR_TESTED_RDS <- file.path(S, "syn52065646_LongRead/jaquino_analysis/diu_satuRn/switchList_satuRn_all_genes.rds")

OUT <- file.path(S, "ERP161086/emtab14167_replication/diu_bambu")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(OUT)

ALPHA <- 0.05
DIF_CUTOFF <- 0.1

write_tsv <- function(x, name) fwrite(as.data.table(x), file.path(OUT, name), sep = "\t", na = "NA")

###############################
# Samples and quantification
###############################

meta <- fread(SR_META, colClasses = list(character = "Individual"))[RegionShort == "PFC"]
meta <- meta[, .(sample = ENA_RUN, individual = Individual,
                 condition = ifelse(DiseaseGroup == "AD", "AD", "CT"),
                 sex = Sex, age = as.numeric(Age))]
setorder(meta, condition, sample)
stopifnot(nrow(meta) == 41, all(file.exists(file.path(QUANT_DIR, paste0(meta$sample, "_quant"), "quant.sf"))))
print(meta[, .N, by = .(condition, sex)])

read_quant <- function(run, column) {
  q <- fread(file.path(QUANT_DIR, paste0(run, "_quant"), "quant.sf"))
  setNames(q[[column]], q$Name)
}
counts_list <- lapply(meta$sample, read_quant, column = "NumReads")
tpm_list <- lapply(meta$sample, read_quant, column = "TPM")
tx_ids <- names(counts_list[[1]])
stopifnot(all(vapply(c(counts_list, tpm_list), function(x) identical(names(x), tx_ids), logical(1))))

isoform_counts <- data.frame(isoform_id = tx_ids, do.call(cbind, counts_list), check.names = FALSE)
isoform_tpm <- data.frame(isoform_id = tx_ids, do.call(cbind, tpm_list), check.names = FALSE)
colnames(isoform_counts)[-1] <- meta$sample
colnames(isoform_tpm)[-1] <- meta$sample

gene_map <- fread(GENE_MAP_FILE)[, .(gene_id = GENEID, gene_name = gene_symbol)]
gene_map <- unique(gene_map[!is.na(gene_name) & gene_name != ""], by = "gene_id")

###############################
# IsoformSwitchAnalyzeR
###############################

design_df <- data.frame(sampleID = meta$sample, condition = meta$condition,
                        age = meta$age, sex = meta$sex, stringsAsFactors = FALSE)

switchList <- importRdata(
  isoformCountMatrix = isoform_counts,
  isoformRepExpression = isoform_tpm,
  designMatrix = design_df,
  isoformExonAnnoation = GTF_FILE,
  comparisonsToMake = data.frame(condition_1 = "CT", condition_2 = "AD"),
  showProgress = FALSE,
  quiet = FALSE,
  removeNonConvensionalChr = FALSE,
  fixStringTieAnnotationProblem = FALSE
)
switchList <- preFilter(switchList)
switchList <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = switchList,
  alpha = ALPHA,
  dIFcutoff = DIF_CUTOFF,
  reduceToSwitchingGenes = FALSE,
  reduceFurtherToGenesWithConsequencePotential = FALSE,
  onlySigIsoforms = FALSE,
  quiet = FALSE
)
switchList$isoformFeatures$gene_name <- gene_map$gene_name[match(switchList$isoformFeatures$gene_id, gene_map$gene_id)]
saveRDS(switchList, file.path(OUT, "switchList_SR41_bambu_annotation_all_genes.rds"))
write_tsv(switchList$designMatrix, "design_matrix.tsv")

sr_summary <- extractSwitchSummary(switchList, alpha = ALPHA, dIFcutoff = DIF_CUTOFF)
write_tsv(sr_summary, "SR41_switch_summary.tsv")

###############################
# Short-read results
###############################

sr_iso <- as.data.table(switchList$isoformFeatures)[, .(isoform_id, gene_id, gene_name, IF1, IF2, dIF,
                                                        isoform_switch_q_value, gene_switch_q_value)]
sr_iso[, novel_bambu_isoform := !grepl("^ENST", isoform_id)]
sr_iso[, significant := !is.na(isoform_switch_q_value) & isoform_switch_q_value < ALPHA & abs(dIF) > DIF_CUTOFF]
setorder(sr_iso, isoform_switch_q_value, na.last = TRUE)
write_tsv(sr_iso, "SR41_isoform_results_all_tested.tsv")
write_tsv(sr_iso[significant == TRUE], "SR41_significant_isoforms.tsv")

sr_sig_genes <- unique(sr_iso[significant == TRUE, gene_id])
sr_tested_genes <- unique(sr_iso$gene_id)

###############################
# Comparison with the long-read DIU
###############################

lr <- as.data.table(readRDS(LR_DIU_RDS)$isoformFeatures)
lr_top <- fread(LR_DIU_TOP)
lr_sig_genes <- unique(lr_top$gene_id)
lr_sig_iso <- unique(lr[!is.na(isoform_switch_q_value) & isoform_switch_q_value < ALPHA & abs(dIF) > DIF_CUTOFF, isoform_id])
lr_tested_genes <- unique(readRDS(LR_TESTED_RDS)$isoformFeatures$gene_id)
stopifnot(length(lr_sig_genes) == 86, length(lr_sig_iso) == 96)

sr_gene_q <- unique(sr_iso[, .(gene_id, SR_gene_q_value = gene_switch_q_value)], by = "gene_id")
gene_cmp <- data.table(gene_id = union(lr_sig_genes, sr_sig_genes))
gene_cmp[, gene_name := gene_map$gene_name[match(gene_id, gene_map$gene_id)]]
gene_cmp[, `:=`(significant_LR = gene_id %in% lr_sig_genes,
                tested_SR = gene_id %in% sr_tested_genes,
                significant_SR = gene_id %in% sr_sig_genes)]
gene_cmp[, LR_gene_q_value := lr_top$gene_switch_q_value[match(gene_id, lr_top$gene_id)]]
gene_cmp <- merge(gene_cmp, sr_gene_q, by = "gene_id", all.x = TRUE)
gene_cmp[, category := fifelse(significant_LR & significant_SR, "Both",
                       fifelse(significant_LR & !tested_SR, "LR only (not tested in SR)",
                       fifelse(significant_LR, "LR only", "SR only")))]
setorder(gene_cmp, category, SR_gene_q_value, na.last = TRUE)
write_tsv(gene_cmp, "comparison_genes_LR_vs_SR41.tsv")

universe <- intersect(lr_tested_genes, sr_tested_genes)
fisher <- fisher.test(table(factor(universe %in% lr_sig_genes, c(TRUE, FALSE)),
                            factor(universe %in% sr_sig_genes, c(TRUE, FALSE))), alternative = "greater")

lr_iso <- lr[, .(isoform_id, gene_id, dIF_LR = dIF, q_LR = isoform_switch_q_value)]
sr_small <- sr_iso[, .(isoform_id, dIF_SR = dIF, q_SR = isoform_switch_q_value, significant_SR = significant)]
iso_cmp <- merge(lr_iso[isoform_id %in% lr_sig_iso], sr_small, by = "isoform_id", all.x = TRUE)
iso_cmp[, gene_name := gene_map$gene_name[match(gene_id, gene_map$gene_id)]]
iso_cmp[, tested_SR := !is.na(q_SR)]
iso_cmp[, significant_SR := significant_SR %in% TRUE]
iso_cmp[, same_direction := fifelse(tested_SR, sign(dIF_LR) == sign(dIF_SR), NA)]
setorder(iso_cmp, q_SR, na.last = TRUE)
write_tsv(iso_cmp, "comparison_LR_significant_isoforms_in_SR41.tsv")

all_iso <- merge(lr_iso, sr_small, by = "isoform_id")
all_iso[, gene_name := gene_map$gene_name[match(gene_id, gene_map$gene_id)]]
write_tsv(all_iso, "comparison_dIF_LR_vs_SR41_isoforms_in_LR_DIU_genes.tsv")
dif_cor <- suppressWarnings(cor.test(all_iso$dIF_LR, all_iso$dIF_SR, method = "spearman"))

n_tested <- sum(iso_cmp$tested_SR)
n_same <- sum(iso_cmp$same_direction, na.rm = TRUE)
binom_p <- if (n_tested > 0) binom.test(n_same, n_tested, 0.5, alternative = "greater")$p.value else NA

report <- c(
  paste("Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("IsoformSwitchAnalyzeR", as.character(packageVersion("IsoformSwitchAnalyzeR")),
        "| DEXSeq", as.character(packageVersion("DEXSeq"))),
  paste("Cohort: E-MTAB-14167 PFC,", sum(meta$condition == "AD"), "AD vs", sum(meta$condition == "CT"), "CT;",
        "design ~ condition + age + sex (+ surrogate variables added by importRdata)"),
  paste("Annotation:", GTF_FILE),
  paste("alpha =", ALPHA, "| dIF cutoff =", DIF_CUTOFF),
  "",
  "Short-read DIU (extractSwitchSummary):",
  capture.output(print(sr_summary, row.names = FALSE)),
  paste("Short-read genes tested:", length(sr_tested_genes),
        "| significant isoforms that are novel Bambu transcripts:",
        sum(sr_iso$significant & sr_iso$novel_bambu_isoform), "of", sum(sr_iso$significant)),
  "",
  "Gene level:",
  paste("  Long-read DIU genes: 86 | tested in short-read:", sum(lr_sig_genes %in% sr_tested_genes),
        "| significant in short-read:", sum(lr_sig_genes %in% sr_sig_genes)),
  paste("  Short-read DIU genes:", length(sr_sig_genes)),
  paste("  Genes tested in both:", length(universe), "| overlap:", sum(universe %in% lr_sig_genes & universe %in% sr_sig_genes),
        "| Fisher exact (one-sided) OR =", signif(fisher$estimate, 3), "p =", signif(fisher$p.value, 3)),
  "",
  "Isoform level (96 long-read switching isoforms):",
  paste("  Tested in short-read:", n_tested, "of 96 | significant in short-read:", sum(iso_cmp$significant_SR),
        "| same dIF direction:", n_same, "of", n_tested, "(one-sided binomial p =", signif(binom_p, 3), ")"),
  paste("  dIF correlation, all isoforms of long-read DIU genes tested in both (n =", nrow(all_iso), "): Spearman rho =",
        signif(dif_cor$estimate, 3), "p =", signif(dif_cor$p.value, 3)),
  "",
  "Long-read DIU genes significant in short-read:",
  capture.output(print(gene_cmp[category == "Both"], row.names = FALSE)),
  "",
  capture.output(sessionInfo())
)
writeLines(report, file.path(OUT, "run_summary.txt"))
cat(report[1:24], sep = "\n")
