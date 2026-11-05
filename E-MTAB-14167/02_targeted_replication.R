#!/usr/bin/env Rscript
# Targeted short-read replication of the long-read findings in the full
# E-MTAB-14167 prefrontal cortex cohort (20 AD vs 21 CT), quantified against the
# long-read Bambu extended annotation (01_salmon_bambu_quant_emtab14167.sh).
#
# Pre-specified targets (from the 6 vs 6 long-read analysis):
#   - 96 switching isoforms (IsoformSwitchAnalyzeR, q < 0.05, |dIF| > 0.1)
#   - 46 differentially used promoters (beta regression, padj <= 0.05, |log2FC| >= 1)
#   - 14 differential PSI events (beta regression ~ condition + sex, FDR < 0.05, |dPSI| >= 0.1)
#
# Each target is retested in the short-read data with covariates; BH is applied
# within each target family. Directional concordance with the long-read effect
# is tested with a one-sided binomial test, at the feature level and at the gene
# level (one pre-specified feature per gene, the one with the largest long-read
# effect, because features within a gene are not independent).
#
# Sample sets:
#   all41        all PFC samples                  ~ condition + age + sex
#   all41_apoe   all PFC samples (sensitivity)    ~ condition + age + sex + APOE4
#   new27        donors NOT in the manuscript 7v7 ~ condition + age + sex  (independent)
#   apoe33       APOE3/3 donors only              ~ condition + age + sex
#   manuscript14 the 7v7 donors in the manuscript ~ condition + sex        (consistency check)
#   all41_neuronal / new27_neuronal: as above + deconvolved neuronal fraction
#                (NNLS EX + INH, syn52065646_LongRead/deconvolution/scripts/01_deconvolution.R)
#
# Long-read composition sensitivity: the promoter and PSI tests are repeated on
# the long-read data with the long-read neuronal fraction as an extra covariate.
#
# Also: gene-level DESeq2 for the 86 DIU genes in long-read and short-read data.

suppressPackageStartupMessages({
  library(data.table)
  library(betareg)
  library(rtracklayer)
  library(tximport)
  library(DESeq2)
  library(ggplot2)
  library(patchwork)
})

S <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing"
LR <- file.path(S, "syn52065646_LongRead/jaquino_analysis")
OUT <- file.path(S, "ERP161086/emtab14167_replication")
QUANT_DIR <- file.path(OUT, "salmon_bambu/quant")
TABLES <- file.path(OUT, "tables")
FIGS <- file.path(OUT, "figures")
dir.create(TABLES, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGS, showWarnings = FALSE, recursive = TRUE)

BAMBU_GTF <- file.path(LR, "NanoporeRNASeq_extended_annotations.gtf")
SR_META <- file.path(S, "E-MTAB-14167/Step01_metadata/sample_metadata_matchedv2.tsv")
ERP_RUNTABLE <- file.path(S, "ERP161086/SraRunTable.csv")
ERP_QUANT_DIR <- file.path(S, "ERP161086/salmon_bambu/quant")
LR_SWITCHLIST <- file.path(LR, "diu/switchList.txt")
LR_PTC <- file.path(LR, "diu/isoform_PTC_status.txt")
LR_DIU_GENES <- file.path(LR, "diu/IsoformSwitchAnalyzeR_top_switches_gene_name.tsv")
LR_PROMOTER <- file.path(S, "syn52065646_LongRead/jaquino_analysis_v2/tables/promoter_usage_86_genes_12samples/promoter_usage_betareg_results_9samples.tsv")
LR_PSI <- file.path(LR, "tables/suppa_beta_regression/all_sig_beta_events.tsv")
LR_SUPPA_DIR <- file.path(LR, "suppa_bambu_gtf")
LR_TPM <- file.path(LR, "bambu_tx_tpm.txt")
LR_META <- file.path(LR, "long_read_metadata.txt")
LR_GENE_COUNTS <- file.path(LR, "NanoporeRNASeq_counts_gene.txt")
NEURO_FILE <- file.path(S, "syn52065646_LongRead/deconvolution/tables/neuronal_fraction.tsv")

EVENT_TYPES <- c("SE", "RI", "MX", "AL", "AF", "A5", "A3")
FDR_CUTOFF <- 0.05
MIN_PER_GROUP <- 3

log_file <- file.path(OUT, "results", "02_targeted_replication_summary.txt")
dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
sink(log_file, split = TRUE)
cat("Targeted replication in E-MTAB-14167 PFC,", format(Sys.time()), "\n\n")

###############################
# Metadata
###############################

meta <- fread(SR_META, colClasses = list(character = "Individual"))[RegionShort == "PFC"]
meta[, `:=`(
  sample = ENA_RUN,
  condition = factor(ifelse(DiseaseGroup == "AD", "AD", "CT"), levels = c("CT", "AD")),
  sex = tolower(Sex),
  age = as.numeric(Age),
  e4 = factor(E4, levels = c("E4neg", "E4pos")),
  individual = as.character(Individual)
)]
erp <- fread(ERP_RUNTABLE, colClasses = list(character = "Individual"))
meta[, in_manuscript := individual %in% erp$Individual]
neuro <- fread(NEURO_FILE, colClasses = list(character = "sample"))[method == "NNLS"]
meta[neuro[platform == "short-read"], neuronal := i.neuronal, on = "sample"]
stopifnot(!anyNA(meta$neuronal))
setorder(meta, condition, sample)

cat("PFC samples:", nrow(meta), "\n")
print(meta[, .N, by = .(condition, sex)][order(condition, sex)])
print(meta[, .(n = .N, age_mean = round(mean(age), 1), age_sd = round(sd(age), 1),
               e4pos = sum(e4 == "E4pos"), in_manuscript = sum(in_manuscript)), by = condition])

###############################
# Annotation
###############################

gtf <- import(BAMBU_GTF)
tx_gr <- gtf[gtf$type == "transcript"]
tx2gene <- data.table(transcript_id = tx_gr$transcript_id, gene_id = tx_gr$gene_id)

###############################
# Short-read quantification
###############################

files <- file.path(QUANT_DIR, paste0(meta$sample, "_quant"), "quant.sf")
names(files) <- meta$sample
stopifnot(all(file.exists(files)))

map_rate <- sapply(meta$sample, function(s) {
  jsonlite::fromJSON(file.path(QUANT_DIR, paste0(s, "_quant"), "aux_info", "meta_info.json"))$percent_mapped
})
meta[, mapping_rate := map_rate[sample]]
cat("\nSalmon mapping rate (Bambu index): min", round(min(map_rate), 1),
    "median", round(median(map_rate), 1), "max", round(max(map_rate), 1), "\n")

txi_tx <- tximport(files, type = "salmon", txOut = TRUE)
tpm <- txi_tx$abundance
tx_gene <- tx2gene$gene_id[match(rownames(tpm), tx2gene$transcript_id)]
stopifnot(!anyNA(tx_gene))
gene_tpm <- rowsum(tpm, tx_gene)

###############################
# Sex check (XIST vs chrY genes)
###############################

y_genes <- c(RPS4Y1 = "ENSG00000129824", DDX3Y = "ENSG00000067048", KDM5D = "ENSG00000012817",
             UTY = "ENSG00000183878", EIF1AY = "ENSG00000198692")
sex_dt <- data.table(
  sample = colnames(gene_tpm),
  XIST = gene_tpm["ENSG00000229807", ],
  chrY = colSums(gene_tpm[intersect(y_genes, rownames(gene_tpm)), , drop = FALSE])
)
sex_dt[, inferred_sex := ifelse(chrY > XIST, "male", "female")]
sex_dt <- merge(sex_dt, meta[, .(sample, individual, metadata_sex = sex)], by = "sample")
fwrite(sex_dt, file.path(TABLES, "qc_sex_check.tsv"), sep = "\t")
mismatch <- sex_dt[inferred_sex != metadata_sex]
cat("\nSex check: ", nrow(mismatch), " mismatch(es) between metadata sex and XIST/chrY expression\n", sep = "")
if (nrow(mismatch) > 0) {
  print(mismatch)
  cat("Using expression-inferred sex for these samples.\n")
  meta[sex_dt, sex := i.inferred_sex, on = "sample"]
}
meta[, sex := factor(sex)]

###############################
# Consistency with the manuscript 7v7 quantification
###############################

erp_files <- file.path(ERP_QUANT_DIR, paste0(erp$Run, "_quant"), "quant.sf")
if (all(file.exists(erp_files))) {
  names(erp_files) <- erp$Run
  erp_tpm <- tximport(erp_files, type = "salmon", txOut = TRUE)$abundance
  pair <- merge(erp[, .(erp_run = Run, individual = Individual)], meta[, .(sample, individual)], by = "individual")
  pair[, spearman_log_tpm := mapply(function(a, b) {
    keep <- erp_tpm[, a] > 0 | tpm[, b] > 0
    cor(log1p(erp_tpm[keep, a]), log1p(tpm[rownames(erp_tpm)[keep], b]), method = "spearman")
  }, erp_run, sample)]
  fwrite(pair, file.path(TABLES, "qc_manuscript14_vs_emtab14167_tpm_correlation.tsv"), sep = "\t")
  cat("\nTranscript TPM correlation, manuscript ERP161086 run vs matching E-MTAB-14167 run (same donor):\n")
  print(summary(pair$spearman_log_tpm))
  # Donors 802 and 827 (both female CT) have their data swapped between the two
  # depositions: ERR13257761 (802) is identical to ERR13232570 (827) and
  # ERR13257759 (827) to ERR13232568 (802). Condition and sex are unaffected.
  cat("Low-correlation pairs (see donor 802/827 label swap note in script):\n")
  print(pair[spearman_log_tpm < 0.99])
}

###############################
# Beta regression helper
###############################

# transform = "smithson": (y * (n - 1) + 0.5) / n, as in the manuscript PSI analysis
# transform = "clip": clip to [1e-4, 1 - 1e-4], as in the manuscript promoter analysis
# test = "lrt": full vs reduced (covariates only) likelihood ratio test
# test = "wald": Wald test on the condition coefficient
fit_beta <- function(y, cd, covars, transform, test, require_group_var = FALSE) {
  d <- data.table(y = as.numeric(y), cd)[is.finite(y)]
  n_ct <- d[condition == "CT", .N]
  n_ad <- d[condition == "AD", .N]
  out <- list(n_CT = n_ct, n_AD = n_ad,
              mean_CT = if (n_ct) mean(d[condition == "CT", y]) else NA_real_,
              mean_AD = if (n_ad) mean(d[condition == "AD", y]) else NA_real_,
              estimate = NA_real_, se = NA_real_, p = NA_real_, note = NA_character_)
  if (n_ct < MIN_PER_GROUP || n_ad < MIN_PER_GROUP) {
    out$note <- "fewer than 3 non-missing values in a group"
    return(out)
  }
  if (transform == "smithson") {
    n <- nrow(d)
    d[, yb := (y * (n - 1) + 0.5) / n]
  } else {
    d[, yb := pmin(pmax(y, 1e-4), 1 - 1e-4)]
  }
  if (uniqueN(d$yb) < 3) {
    out$note <- "fewer than 3 distinct values"
    return(out)
  }
  if (require_group_var && !all(d[, var(yb), by = condition]$V1 > 0)) {
    out$note <- "no variance in a group"
    return(out)
  }
  covars <- covars[vapply(covars, function(v) uniqueN(d[[v]]) > 1, logical(1))]
  f_full <- reformulate(c("condition", covars), "yb")
  f_red <- reformulate(if (length(covars)) covars else "1", "yb")
  fit <- tryCatch(betareg(f_full, data = d), error = function(e) e)
  if (inherits(fit, "error")) {
    out$note <- paste("full model failed:", conditionMessage(fit))
    return(out)
  }
  cf <- summary(fit)$coefficients$mean
  if (!"conditionAD" %in% rownames(cf)) {
    out$note <- "condition coefficient not estimable"
    return(out)
  }
  out$estimate <- cf["conditionAD", "Estimate"]
  out$se <- cf["conditionAD", "Std. Error"]
  if (test == "wald") {
    out$p <- cf["conditionAD", "Pr(>|z|)"]
  } else {
    red <- tryCatch(betareg(f_red, data = d), error = function(e) e)
    if (inherits(red, "error")) {
      out$note <- paste("reduced model failed:", conditionMessage(red))
      return(out)
    }
    lr <- 2 * (as.numeric(logLik(fit)) - as.numeric(logLik(red)))
    out$p <- pchisq(max(lr, 0), df = 1, lower.tail = FALSE)
  }
  out
}

test_features <- function(mat, cd, covars, transform, test, require_group_var = FALSE) {
  res <- lapply(rownames(mat), function(f) {
    r <- fit_beta(mat[f, cd$sample], cd, covars, transform, test, require_group_var)
    as.data.table(c(list(feature_id = f), r))
  })
  rbindlist(res)
}

###############################
# Long-read targets
###############################

sl <- fread(LR_SWITCHLIST)
ptc <- fread(LR_PTC)[, .(isoform_id, PTC)]
iso_targets <- sl[isoform_switch_q_value < 0.05 & abs(dIF) > 0.1,
                  .(feature_id = isoform_id, gene_id, gene_name,
                    lr_mean_CT = IF1, lr_mean_AD = IF2, lr_effect = dIF, lr_q = isoform_switch_q_value)]
iso_targets <- merge(iso_targets, ptc, by.x = "feature_id", by.y = "isoform_id", all.x = TRUE)
diu_genes <- fread(LR_DIU_GENES)[, .(gene_id, gene_name)]
stopifnot(nrow(iso_targets) == 96, uniqueN(iso_targets$gene_id) == 86, nrow(diu_genes) == 86)

pr_lr <- fread(LR_PROMOTER)
prom_targets <- pr_lr[padj <= 0.05 & abs(logitFC / log(2)) >= 1,
                      .(feature_id = promoter_id, gene_id, gene_name, lr_effect = logitFC, lr_q = padj)]
prom_targets[, switching_gene := any(lr_effect > 0) & any(lr_effect < 0), by = gene_id]
stopifnot(nrow(prom_targets) == 46)

psi_lr <- fread(LR_PSI)
psi_targets <- psi_lr[, .(feature_id = event_id, event_type, gene_id, gene_name = gene_symbol,
                          lr_mean_CT = mean_CT, lr_mean_AD = mean_AD, lr_effect = dPSI, lr_q = adj_p_value)]
stopifnot(nrow(psi_targets) == 14)

cat("\nTargets: ", nrow(iso_targets), " isoforms in ", uniqueN(iso_targets$gene_id), " genes; ",
    nrow(prom_targets), " promoters in ", uniqueN(prom_targets$gene_id), " genes (",
    uniqueN(prom_targets[switching_gene == TRUE]$gene_id), " switching genes); ",
    nrow(psi_targets), " PSI events in ", uniqueN(psi_targets$gene_id), " genes\n", sep = "")

###############################
# Feature matrices
###############################

# Isoform fraction: transcript TPM / gene TPM; NA where gene TPM < 1
isoform_fraction <- function(tpm_mat, ids) {
  g <- tx2gene$gene_id[match(ids, tx2gene$transcript_id)]
  gene_tx <- tx2gene[gene_id %in% g & transcript_id %in% rownames(tpm_mat)]
  gt <- rowsum(tpm_mat[gene_tx$transcript_id, , drop = FALSE], gene_tx$gene_id)
  frac <- matrix(NA_real_, length(ids), ncol(tpm_mat), dimnames = list(ids, colnames(tpm_mat)))
  ok <- ids %in% rownames(tpm_mat)
  den <- gt[g[ok], , drop = FALSE]
  num <- tpm_mat[ids[ok], , drop = FALSE] / den
  num[den < 1] <- NA
  frac[ok, ] <- num
  frac
}

# Promoter usage, reproducing promoter_usage_analysis_86_genes_12samples.r:
# one -1000/+100 bp window per transcript TSS, TPM summed over transcripts sharing
# the window, usage = promoter TPM / gene TPM, NA where gene TPM < 1,
# genes with >= 2 promoters.
strand_char <- as.character(strand(tx_gr))
tss <- ifelse(strand_char == "-", end(tx_gr), start(tx_gr))
prom_map <- data.table(
  transcript_id = tx_gr$transcript_id,
  gene_id = tx_gr$gene_id,
  promoter_id = paste0(seqnames(tx_gr), ":", pmax(1L, tss - 1000L), "-", tss + 100L, "(", strand_char, ")")
)[gene_id %in% diu_genes$gene_id]
prom_gene <- prom_map[, .(gene_id = sort(gene_id)[1]), by = promoter_id]
multi_prom_genes <- prom_gene[, .N, by = gene_id][N >= 2, gene_id]

promoter_usage <- function(tpm_mat) {
  pm <- prom_map[transcript_id %in% rownames(tpm_mat)]
  ptpm <- rowsum(tpm_mat[pm$transcript_id, , drop = FALSE], pm$promoter_id)
  pg <- prom_gene$gene_id[match(rownames(ptpm), prom_gene$promoter_id)]
  gt <- rowsum(ptpm, pg)
  den <- gt[pg, , drop = FALSE]
  u <- ptpm / den
  u[den < 1] <- NA
  u[pg %in% multi_prom_genes, , drop = FALSE]
}

# PSI from SUPPA IOE definitions: sum(alternative TPM) / sum(total TPM)
ioe <- rbindlist(lapply(EVENT_TYPES, function(e) {
  fread(file.path(LR_SUPPA_DIR, paste0("FC_events_", e, "_strict.ioe")))
}))
ioe_t <- ioe[event_id %in% psi_targets$feature_id]
stopifnot(nrow(ioe_t) == 14)

event_psi <- function(tpm_mat) {
  m <- t(vapply(seq_len(nrow(ioe_t)), function(i) {
    alt <- strsplit(ioe_t$alternative_transcripts[i], ",")[[1]]
    tot <- strsplit(ioe_t$total_transcripts[i], ",")[[1]]
    num <- colSums(tpm_mat[intersect(alt, rownames(tpm_mat)), , drop = FALSE])
    den <- colSums(tpm_mat[intersect(tot, rownames(tpm_mat)), , drop = FALSE])
    ifelse(den > 0, num / den, NA_real_)
  }, numeric(ncol(tpm_mat))))
  dimnames(m) <- list(ioe_t$event_id, colnames(tpm_mat))
  m
}

###############################
# Verify the re-implementation on the long-read data
###############################

lr_meta <- fread(LR_META)
lr_meta[, `:=`(sample = sub("^X", "", sample),
               condition = factor(condition, levels = c("CT", "AD")),
               sex = factor(tolower(sex)))]
lr_meta[neuro[platform == "long-read"], neuronal := i.neuronal, on = "sample"]
stopifnot(!anyNA(lr_meta$neuronal))
lr_tpm <- as.matrix(read.table(LR_TPM, header = TRUE, sep = "\t", check.names = FALSE))
lr_tpm <- lr_tpm[, lr_meta$sample]

lr_prom <- promoter_usage(lr_tpm)
lr_prom_res <- test_features(lr_prom[prom_targets$feature_id, ], lr_meta, character(), "clip", "wald", TRUE)
chk <- merge(lr_prom_res[, .(feature_id, estimate)], prom_targets[, .(feature_id, lr_effect)], by = "feature_id")
cat("\nCheck - promoter logitFC recomputed from long-read TPM vs manuscript: max abs difference =",
    signif(max(abs(chk$estimate - chk$lr_effect)), 3), "\n")

lr_psi <- event_psi(lr_tpm)
suppa_mat <- do.call(rbind, lapply(EVENT_TYPES, function(e) {
  x <- as.matrix(read.table(file.path(LR_SUPPA_DIR, paste0("psi_", e, ".psi")),
                            header = TRUE, sep = "\t", check.names = FALSE))
  x[rownames(x) %in% psi_targets$feature_id, lr_meta$sample, drop = FALSE]
}))
cat("Check - PSI recomputed from long-read TPM vs SUPPA psi files: max abs difference =",
    signif(max(abs(lr_psi[rownames(suppa_mat), ] - suppa_mat), na.rm = TRUE), 3), "\n")

###############################
# Long-read composition sensitivity
###############################

# Promoters: all promoters of the 86 genes, as in the manuscript (Wald, BH over
# all tested promoters), without and with the neuronal fraction.
lr_prom_all <- rbindlist(lapply(c(original = "", neuronal = "neuronal"), function(cv) {
  r <- test_features(lr_prom, lr_meta, setdiff(cv, ""), "clip", "wald", TRUE)[!is.na(p)]
  r[, `:=`(model = if (cv == "") "~ condition" else "~ condition + neuronal",
           padj = p.adjust(p, "BH"))]
  r[, sig := padj <= 0.05 & abs(estimate / log(2)) >= 1]
}))
lr_prom_cmp <- dcast(lr_prom_all[feature_id %in% prom_targets$feature_id], feature_id ~ model,
                     value.var = c("estimate", "padj", "sig"))
cat("\n==== Long-read promoters with neuronal fraction as covariate ====\n")
cat("Promoters tested:", lr_prom_all[model == "~ condition", .N],
    "| significant original:", lr_prom_all[model == "~ condition", sum(sig)],
    "| significant with neuronal covariate:", lr_prom_all[model == "~ condition + neuronal", sum(sig)],
    "| of the 46 targets still significant:", sum(lr_prom_cmp[["sig_~ condition + neuronal"]], na.rm = TRUE), "\n")
sw <- lr_prom_all[model == "~ condition + neuronal" & sig == TRUE]
sw <- merge(sw, prom_gene, by.x = "feature_id", by.y = "promoter_id")
cat("Promoter-switching genes with neuronal covariate:",
    sw[, .(up = any(estimate > 0), down = any(estimate < 0)), by = gene_id][up & down, .N], "(original 11)\n")

# PSI events: ~ condition + sex (+ neuronal), LRT against the model without condition
lr_psi_cmp <- rbindlist(lapply(list(original = c("sex"), neuronal = c("sex", "neuronal")), function(cv) {
  r <- test_features(lr_psi, lr_meta, cv, "smithson", "lrt")
  r[, `:=`(model = paste("~ condition +", paste(cv, collapse = " + ")), fdr_within_14 = p.adjust(p, "BH"))]
}))
cat("\n==== Long-read PSI events (14 targets) with neuronal fraction as covariate ====\n")
print(dcast(merge(lr_psi_cmp, psi_targets[, .(feature_id, gene_name, event_type)], by = "feature_id"),
            gene_name + event_type + feature_id ~ model, value.var = "p")[, feature_id := NULL][])
cat("FDR < 0.05 (BH within the 14):", lr_psi_cmp[, .(n_sig = sum(fdr_within_14 < 0.05, na.rm = TRUE)), by = model][, paste(model, n_sig, collapse = "; ")], "\n")
fwrite(lr_prom_all, file.path(TABLES, "longread_promoter_neuronal_sensitivity.tsv"), sep = "\t")
fwrite(lr_psi_cmp, file.path(TABLES, "longread_psi_neuronal_sensitivity.tsv"), sep = "\t")

###############################
# Short-read feature matrices and tests
###############################

missing_iso <- setdiff(iso_targets$feature_id, rownames(tpm))
if (length(missing_iso)) cat("Isoforms absent from Salmon output:", missing_iso, "\n")

sr_if <- isoform_fraction(tpm, iso_targets$feature_id)
sr_prom <- promoter_usage(tpm)
missing_prom <- setdiff(prom_targets$feature_id, rownames(sr_prom))
if (length(missing_prom)) cat("Promoters not reconstructed in short-read data:", missing_prom, "\n")
sr_prom <- sr_prom[intersect(prom_targets$feature_id, rownames(sr_prom)), , drop = FALSE]
sr_psi <- event_psi(tpm)
iso_mean_tpm <- setNames(rep(0, nrow(iso_targets)), iso_targets$feature_id)
present <- intersect(iso_targets$feature_id, rownames(tpm))
iso_mean_tpm[present] <- rowMeans(tpm[present, ])

sets <- list(
  all41        = list(samples = meta$sample, covars = c("age", "sex")),
  all41_apoe   = list(samples = meta$sample, covars = c("age", "sex", "e4")),
  new27        = list(samples = meta[in_manuscript == FALSE, sample], covars = c("age", "sex")),
  all41_neuronal = list(samples = meta$sample, covars = c("age", "sex", "neuronal")),
  new27_neuronal = list(samples = meta[in_manuscript == FALSE, sample], covars = c("age", "sex", "neuronal")),
  apoe33       = list(samples = meta[Genotype == "ApoE3,3", sample], covars = c("age", "sex")),
  manuscript14 = list(samples = meta[in_manuscript == TRUE, sample], covars = c("sex"))
)

families <- list(
  isoform  = list(mat = sr_if,   targets = iso_targets,  transform = "smithson", test = "lrt",  gv = FALSE),
  promoter = list(mat = sr_prom, targets = prom_targets, transform = "clip",     test = "wald", gv = TRUE),
  psi      = list(mat = sr_psi,  targets = psi_targets,  transform = "smithson", test = "lrt",  gv = FALSE)
)

rep_res <- rbindlist(lapply(names(sets), function(sn) {
  st <- sets[[sn]]
  cd <- meta[sample %in% st$samples, .(sample, condition, sex, age, e4, neuronal)]
  rbindlist(lapply(names(families), function(fn) {
    fam <- families[[fn]]
    r <- test_features(fam$mat, cd, st$covars, fam$transform, fam$test, fam$gv)
    r <- merge(fam$targets[, .(feature_id, gene_id, gene_name, lr_effect, lr_q)], r,
               by = "feature_id", all.x = TRUE)
    r[, `:=`(set = sn, family = fn,
             n_CT_set = cd[condition == "CT", .N], n_AD_set = cd[condition == "AD", .N],
             model = paste("~ condition +", paste(st$covars, collapse = " + ")))]
    r[, sr_delta := mean_AD - mean_CT]
    r[, fdr := p.adjust(p, method = "BH")]
    r[, concordant := ifelse(is.na(estimate), NA, sign(estimate) == sign(lr_effect))]
    r
  }))
}))
rep_res[family == "isoform", sr_mean_isoform_tpm := iso_mean_tpm[feature_id]]
rep_res[, replicated := !is.na(fdr) & fdr < FDR_CUTOFF & concordant %in% TRUE]
setcolorder(rep_res, c("set", "family", "model", "n_CT_set", "n_AD_set", "feature_id", "gene_id", "gene_name",
                       "lr_effect", "lr_q", "n_CT", "n_AD", "mean_CT", "mean_AD", "sr_delta",
                       "estimate", "se", "p", "fdr", "concordant", "replicated"))
fwrite(rep_res, file.path(TABLES, "targeted_replication_all_features.tsv"), sep = "\t")

###############################
# Summaries
###############################

binom_greater <- function(k, n) if (n > 0) binom.test(k, n, 0.5, alternative = "greater")$p.value else NA_real_

summarise_set <- function(r) {
  tested <- r[!is.na(estimate)]
  # one pre-specified feature per gene: largest |long-read effect|
  per_gene <- tested[order(-abs(lr_effect))][, .SD[1], by = gene_id]
  rho <- if (nrow(tested) >= 5) cor.test(tested$lr_effect, tested$estimate, method = "spearman", exact = FALSE) else NULL
  data.table(
    n_targets = nrow(r),
    n_tested = nrow(tested),
    n_concordant = sum(tested$concordant),
    pct_concordant = round(100 * mean(tested$concordant), 1),
    binom_p_feature = signif(binom_greater(sum(tested$concordant), nrow(tested)), 3),
    n_genes = nrow(per_gene),
    n_genes_concordant = sum(per_gene$concordant),
    binom_p_gene = signif(binom_greater(sum(per_gene$concordant), nrow(per_gene)), 3),
    spearman_rho = if (is.null(rho)) NA_real_ else round(unname(rho$estimate), 3),
    spearman_p = if (is.null(rho)) NA_real_ else signif(rho$p.value, 3),
    n_nominal_p05 = sum(tested$p < 0.05),
    n_fdr05 = sum(tested$fdr < FDR_CUTOFF),
    n_replicated = sum(tested$replicated)
  )
}

summary_tab <- rep_res[, summarise_set(.SD), by = .(set, family, model, n_CT_set, n_AD_set)]
iso_expr <- rep_res[family == "isoform" & sr_mean_isoform_tpm >= 1,
                    summarise_set(.SD), by = .(set, model, n_CT_set, n_AD_set)]
iso_expr[, family := "isoform (mean SR TPM >= 1)"]
summary_tab <- rbind(summary_tab, iso_expr, use.names = TRUE)
setorder(summary_tab, set, family)
fwrite(summary_tab, file.path(TABLES, "targeted_replication_summary.tsv"), sep = "\t")

cat("\n==== Targeted replication summary ====\n")
print(summary_tab[, .(set, family, n = paste0(n_AD_set, "v", n_CT_set), n_tested, n_concordant, pct_concordant,
                      binom_p_feature, n_genes, n_genes_concordant, binom_p_gene,
                      spearman_rho, spearman_p, n_nominal_p05, n_fdr05, n_replicated)], nrows = 100)

cat("\nReplicated features (FDR < 0.05 within family and same direction as long-read):\n")
print(rep_res[replicated == TRUE, .(set, family, gene_name, feature_id, lr_effect = signif(lr_effect, 3),
                                    sr_delta = signif(sr_delta, 3), estimate = signif(estimate, 3),
                                    fdr = signif(fdr, 3))], nrows = 200)

cat("\nNominal p < 0.05 in the independent new27 set:\n")
print(rep_res[set == "new27" & p < 0.05, .(family, gene_name, feature_id, lr_effect = signif(lr_effect, 3),
                                          estimate = signif(estimate, 3), p = signif(p, 3),
                                          fdr = signif(fdr, 3), concordant)], nrows = 200)

###############################
# Figures: long-read vs short-read effects
###############################

family_labels <- c(isoform = "Isoform fraction (96 isoforms)",
                   promoter = "Promoter usage (46 promoters)",
                   psi = "PSI events (14 events)")
x_labels <- c(isoform = "Long-read dIF (AD - CT)",
              promoter = "Long-read logit fold change",
              psi = "Long-read dPSI (AD - CT)")
col_sig <- "#2a78d6"
col_ns <- "#a3a29c"

plot_set <- function(sn) {
  st <- summary_tab[set == sn]
  plots <- lapply(names(family_labels), function(fn) {
    d <- rep_res[set == sn & family == fn & !is.na(estimate)]
    d[, status := factor(ifelse(fdr < FDR_CUTOFF, "Short-read FDR < 0.05", "Not significant"),
                         levels = c("Not significant", "Short-read FDR < 0.05"))]
    s <- st[family == fn]
    ggplot(d, aes(lr_effect, estimate)) +
      geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.4) +
      geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.4) +
      geom_point(aes(colour = status), size = 2.4, alpha = 0.85) +
      scale_colour_manual(values = c("Not significant" = col_ns, "Short-read FDR < 0.05" = col_sig),
                          drop = FALSE, name = NULL) +
      labs(title = family_labels[[fn]],
           subtitle = sprintf("Same direction: %d/%d features (binomial p = %s)\n%d/%d genes (p = %s); Spearman rho = %s",
                              s$n_concordant, s$n_tested, format(s$binom_p_feature, digits = 2),
                              s$n_genes_concordant, s$n_genes, format(s$binom_p_gene, digits = 2),
                              format(s$spearman_rho, digits = 2)),
           x = x_labels[[fn]], y = "Short-read AD coefficient (logit scale)") +
      theme_minimal(base_size = 11) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
            plot.title = element_text(face = "bold", size = 11),
            plot.subtitle = element_text(size = 9, colour = "grey30"),
            legend.position = "bottom")
  })
  n_lab <- sprintf("%s: %d AD vs %d CT, %s", sn, st$n_AD_set[1], st$n_CT_set[1], st$model[1])
  wrap_plots(plots, nrow = 1, guides = "collect") +
    plot_annotation(title = paste("Long-read targets retested in E-MTAB-14167 PFC -", n_lab)) &
    theme(legend.position = "bottom")
}

for (sn in c("all41", "new27", "apoe33")) {
  ggsave(file.path(FIGS, paste0("replication_scatter_", sn, ".png")), plot_set(sn),
         width = 15, height = 5.4, dpi = 300, bg = "white")
}

###############################
# Step 4: gene-level expression of the 86 DIU genes
###############################

# Long-read gene-level DESeq2 (6 AD vs 6 CT, ~ sex + condition)
lr_counts <- as.matrix(read.table(LR_GENE_COUNTS, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE))
lr_counts <- round(lr_counts[, lr_meta$sample])
lr_dds <- DESeqDataSetFromMatrix(lr_counts, as.data.frame(lr_meta[, .(condition, sex)]), ~ sex + condition)
lr_dds <- lr_dds[rowSums(counts(lr_dds) >= 10) >= 6, ]
lr_dds <- DESeq(lr_dds, quiet = TRUE)
lr_de <- as.data.table(as.data.frame(results(lr_dds, contrast = c("condition", "AD", "CT"))), keep.rownames = "gene_id")

# Short-read gene-level DESeq2 from Salmon (Bambu annotation)
txi_gene <- summarizeToGene(txi_tx, as.data.frame(tx2gene))
sr_deseq <- function(samples, design) {
  cd <- as.data.frame(meta[match(samples, sample), .(condition, sex, age = as.numeric(scale(age)), e4)])
  rownames(cd) <- samples
  txi <- lapply(txi_gene[c("abundance", "counts", "length")], function(m) m[, samples])
  txi$countsFromAbundance <- txi_gene$countsFromAbundance
  dds <- DESeqDataSetFromTximport(txi, cd, design)
  dds <- dds[rowSums(counts(dds) >= 10) >= min(table(cd$condition)), ]
  dds <- DESeq(dds, quiet = TRUE)
  as.data.table(as.data.frame(results(dds, contrast = c("condition", "AD", "CT"))), keep.rownames = "gene_id")
}
sr_de_all <- sr_deseq(meta$sample, ~ age + sex + condition)
sr_de_new <- sr_deseq(meta[in_manuscript == FALSE, sample], ~ age + sex + condition)

pick <- function(de, prefix) {
  x <- de[gene_id %in% diu_genes$gene_id, .(gene_id, log2FoldChange, pvalue, padj)]
  setnames(x, c("gene_id", paste0(prefix, c("_log2FC", "_p", "_padj"))))
  x
}

iso_best <- rep_res[set %in% c("all41", "new27") & family == "isoform" & !is.na(estimate)][
  order(-abs(lr_effect))][, .SD[1], by = .(set, gene_id)]
iso_best_w <- dcast(iso_best, gene_id ~ set, value.var = c("p", "fdr", "concordant"))

gene_tab <- Reduce(function(a, b) merge(a, b, by = "gene_id", all.x = TRUE), list(
  diu_genes,
  iso_targets[, .(n_switching_isoforms = .N,
                  n_ptc_isoforms = sum(PTC %in% TRUE),
                  max_abs_lr_dIF = round(max(abs(lr_effect)), 3)), by = gene_id],
  data.table(gene_id = diu_genes$gene_id,
             promoter_switching_gene = diu_genes$gene_id %in% prom_targets[switching_gene == TRUE, gene_id]),
  pick(lr_de, "LR_6v6"),
  pick(sr_de_all, "SR_all41"),
  pick(sr_de_new, "SR_new27"),
  iso_best_w
))
setorder(gene_tab, LR_6v6_padj, na.last = TRUE)
fwrite(gene_tab, file.path(TABLES, "gene_level_86_DIU_genes.tsv"), sep = "\t")

gl_summary <- function(de, label) {
  x <- de[gene_id %in% diu_genes$gene_id & !is.na(pvalue)]
  all_tested <- de[!is.na(pvalue)]
  data.table(dataset = label, genes_tested = nrow(x),
             nominal_p05 = sum(x$pvalue < 0.05),
             padj05 = sum(x$padj < 0.05, na.rm = TRUE),
             genome_wide_padj05 = sum(all_tested$padj < 0.05, na.rm = TRUE),
             genome_wide_tested = nrow(all_tested))
}
gl <- rbind(gl_summary(lr_de, "long-read 6v6 (~ sex + condition)"),
            gl_summary(sr_de_all, "short-read all41 (~ age + sex + condition)"),
            gl_summary(sr_de_new, "short-read new27 (~ age + sex + condition)"))
fwrite(gl, file.path(TABLES, "gene_level_86_DIU_genes_summary.tsv"), sep = "\t")
cat("\n==== Gene-level differential expression of the 86 DIU genes ====\n")
print(gl)
cat("\nDIU genes with gene-level padj < 0.05 in any dataset:\n")
print(gene_tab[LR_6v6_padj < 0.05 | SR_all41_padj < 0.05 | SR_new27_padj < 0.05,
               .(gene_name, LR_6v6_log2FC, LR_6v6_padj, SR_all41_log2FC, SR_all41_padj, SR_new27_padj)])

# Positive control: canonical bulk-cortex AD genes (neuronal down, glial up)
# to show whether the cohort carries a detectable AD signal at all.
pos_ctrl <- data.table(
  gene_name = c("SST", "VGF", "NPTX2", "CRH", "BDNF", "SCG2", "GFAP", "SERPINA3", "CD44", "C4B", "HSPB1", "MT1M"),
  gene_id = c("ENSG00000157005", "ENSG00000128564", "ENSG00000106236", "ENSG00000147571", "ENSG00000176697",
              "ENSG00000171951", "ENSG00000131095", "ENSG00000196136", "ENSG00000026508", "ENSG00000224389",
              "ENSG00000106211", "ENSG00000205364"),
  expected = rep(c("down", "up"), each = 6)
)
pc_tab <- rbindlist(lapply(list(all41 = sr_de_all, new27 = sr_de_new), function(de) {
  merge(pos_ctrl, de[, .(gene_id, log2FoldChange, pvalue, padj)], by = "gene_id")
}), idcol = "set")
pc_tab[, as_expected := ifelse(expected == "down", log2FoldChange < 0, log2FoldChange > 0)]
fwrite(pc_tab, file.path(TABLES, "positive_control_AD_genes.tsv"), sep = "\t")
cat("\n==== Positive control: canonical AD genes in the expected direction ====\n")
print(pc_tab[, .(n = .N, as_expected = sum(as_expected),
                 sign_test_p = signif(binom.test(sum(as_expected), .N, 0.5, alternative = "greater")$p.value, 3),
                 nominal_p05 = sum(pvalue < 0.05), padj05 = sum(padj < 0.05, na.rm = TRUE)), by = set])

saveRDS(list(meta = meta, results = rep_res, summary = summary_tab, gene_tab = gene_tab,
             lr_de = lr_de, sr_de_all = sr_de_all, sr_de_new = sr_de_new),
        file.path(TABLES, "targeted_replication.rds"))
cat("\nDone:", format(Sys.time()), "\n")
sink()
