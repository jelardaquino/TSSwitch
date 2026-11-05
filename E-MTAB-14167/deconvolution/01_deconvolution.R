#!/usr/bin/env Rscript
# Cell-type deconvolution of the long-read (6 AD vs 6 CT, Heberle et al.) and
# short-read (E-MTAB-14167 PFC, 20 AD vs 21 CT) frontal cortex samples, using
# the GSE174367 (Morabito et al. 2021) prefrontal cortex snRNA-seq as reference.
#
# Reference: per-donor pseudobulk profiles for 7 cell types (EX, INH, ASC, MG,
# ODC, OPC, PER.END), donor x cell type groups with >= 20 nuclei.
#
# Two methods, run separately per bulk platform:
#   NNLS     (primary) non-negative least squares on linear CPM at the 50 most
#            cell-type-specific reference genes per cell type (own-type mean CPM >= 20,
#            ranked by log2 ratio to the highest other type); each gene is scaled by
#            its mean reference expression so highly expressed genes do not dominate;
#            proportions normalised to sum to 1
#   dtangle  (check) marker-based, log scale
# Each estimate is validated against single canonical marker genes. In exploratory
# comparisons NNLS tracked canonical markers best in the long-read data; microglia
# were not estimated reliably by any method on either platform.
#
# Bulk input scale: long-read ONT cDNA reads count molecules, as do 10x UMIs, so
# long-read uses gene CPM; Illumina full-length short reads use gene TPM.
# Mitochondrial genes are excluded (absent from nuclear RNA).

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(rhdf5)
  library(dtangle)
  library(nnls)
  library(tximport)
  library(rtracklayer)
  library(ggplot2)
  library(patchwork)
})

S <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing"
OUT <- file.path(S, "syn52065646_LongRead/deconvolution")
TABLES <- file.path(OUT, "tables")
FIGS <- file.path(OUT, "figures")
dir.create(TABLES, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGS, showWarnings = FALSE, recursive = TRUE)

SN_DIR <- file.path(S, "GSE174367_AD_bulk_snRNA_snATAC")
SN_H5 <- file.path(SN_DIR, "snRNA/preprocessed_snRNA/GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5")
SN_META <- file.path(SN_DIR, "GSE174367_snRNA-seq_cell_meta.csv.gz")

LR <- file.path(S, "syn52065646_LongRead/jaquino_analysis")
LR_GENE_COUNTS <- file.path(LR, "NanoporeRNASeq_counts_gene.txt")
LR_META <- file.path(LR, "long_read_metadata.txt")
BAMBU_GTF <- file.path(LR, "NanoporeRNASeq_extended_annotations.gtf")

SR_RDS <- file.path(S, "ERP161086/emtab14167_replication/tables/targeted_replication.rds")
SR_QUANT_DIR <- file.path(S, "ERP161086/emtab14167_replication/salmon_bambu/quant")

MIN_NUCLEI <- 20
N_NNLS_MARKERS <- 50
MIN_MARKER_CPM <- 20
CELL_TYPES <- c("EX", "INH", "ASC", "MG", "ODC", "OPC", "PER.END")
CANONICAL <- data.table(
  gene = c("SNAP25", "SYT1", "GAD1", "GAD2", "GFAP", "AQP4", "CX3CR1", "P2RY12", "MBP", "PLP1", "PDGFRA", "CLDN5", "FLT1"),
  gene_id = c("ENSG00000132639", "ENSG00000067715", "ENSG00000128683", "ENSG00000136750", "ENSG00000131095",
              "ENSG00000171885", "ENSG00000168329", "ENSG00000169313", "ENSG00000197971", "ENSG00000123560",
              "ENSG00000134853", "ENSG00000184113", "ENSG00000102755"),
  cell_type = c("EX", "EX", "INH", "INH", "ASC", "ASC", "MG", "MG", "ODC", "ODC", "OPC", "PER.END", "PER.END")
)

log_file <- file.path(OUT, "01_deconvolution_summary.txt")
sink(log_file, split = TRUE)
cat("Deconvolution,", format(Sys.time()), "\n\n")

###############################
# snRNA-seq reference
###############################

h5 <- function(x) as.vector(h5read(SN_H5, x))
shape <- h5("/matrix/shape")
sn <- sparseMatrix(i = h5("/matrix/indices") + 1, p = h5("/matrix/indptr"), x = as.numeric(h5("/matrix/data")),
                   dims = shape, dimnames = list(sub("\\.\\d+$", "", h5("/matrix/features/id")), h5("/matrix/barcodes")))
sn_names <- setNames(h5("/matrix/features/name"), rownames(sn))
H5close()

cells <- fread(SN_META)[Barcode %in% colnames(sn) & Cell.Type %in% CELL_TYPES]
cells[, group := paste(SampleID, Cell.Type, sep = "|")]
keep_groups <- cells[, .N, by = group][N >= MIN_NUCLEI, group]
cells <- cells[group %in% keep_groups]
cat("Reference nuclei:", nrow(cells), "from", uniqueN(cells$SampleID), "donors\n")
print(cells[, .(nuclei = .N, donors = uniqueN(SampleID)), by = Cell.Type][order(Cell.Type)])

ind <- sparseMatrix(i = match(cells$Barcode, colnames(sn)), j = match(cells$group, keep_groups), x = 1,
                    dims = c(ncol(sn), length(keep_groups)))
pb <- as.matrix(sn %*% ind)
colnames(pb) <- keep_groups
pb <- rowsum(pb, rownames(pb))   # collapse duplicated Ensembl IDs
pb_type <- sub(".*\\|", "", colnames(pb))

mt_genes <- names(sn_names)[grepl("^MT-", sn_names)]

###############################
# Bulk matrices
###############################

lr_meta <- fread(LR_META)
lr_meta[, sample := sub("^X", "", sample)]
lr_counts <- as.matrix(read.table(LR_GENE_COUNTS, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE))
lr_bulk <- lr_counts[, lr_meta$sample]

sr <- readRDS(SR_RDS)
sr_meta <- sr$meta
gtf <- import(BAMBU_GTF)
tx2gene <- unique(data.frame(tx = gtf$transcript_id, gene = gtf$gene_id)[gtf$type == "transcript", ])
sr_files <- setNames(file.path(SR_QUANT_DIR, paste0(sr_meta$sample, "_quant"), "quant.sf"), sr_meta$sample)
sr_bulk <- tximport(sr_files, type = "salmon", tx2gene = tx2gene)$abundance

###############################
# Deconvolution
###############################

cpm <- function(m) t(t(m) / colSums(m)) * 1e6

deconvolve <- function(bulk, label) {
  genes <- setdiff(intersect(rownames(pb), rownames(bulk)), mt_genes)
  ref <- cpm(pb[genes, ])
  blk <- cpm(bulk[genes, , drop = FALSE])
  # genes expressed in the reference (mean CPM >= 1 in at least one cell type) and in bulk
  type_mean <- sapply(CELL_TYPES, function(ct) rowMeans(ref[, pb_type == ct, drop = FALSE]))
  expressed <- apply(type_mean, 1, max) >= 1 & rowMeans(blk) >= 1
  ref <- ref[expressed, ]
  blk <- blk[expressed, ]
  type_mean <- type_mean[expressed, ]
  cat("\n[", label, "] genes used:", sum(expressed), "\n")

  Y <- rbind(t(log2(ref + 1)), t(log2(blk + 1)))
  pure <- lapply(CELL_TYPES, function(ct) which(pb_type == ct))
  names(pure) <- CELL_TYPES
  mk <- find_markers(Y, pure_samples = pure, data_type = "rna-seq", marker_method = "ratio")
  dt <- dtangle(Y, pure_samples = pure, markers = mk$L, data_type = "rna-seq")
  p_dt <- dt$estimates[(ncol(ref) + 1):nrow(Y), , drop = FALSE]
  colnames(p_dt) <- CELL_TYPES

  spec_markers <- lapply(CELL_TYPES, function(ct) {
    other <- apply(type_mean[, setdiff(CELL_TYPES, ct), drop = FALSE], 1, max)
    sc <- log2((type_mean[, ct] + 1) / (other + 1))
    sc[type_mean[, ct] < MIN_MARKER_CPM] <- -Inf
    names(sort(sc, decreasing = TRUE))[seq_len(N_NNLS_MARKERS)]
  })
  names(spec_markers) <- CELL_TYPES
  mk_genes <- unique(unlist(spec_markers))
  A <- type_mean[mk_genes, CELL_TYPES]
  scale_f <- rowMeans(A)
  p_nn <- t(apply(blk[mk_genes, , drop = FALSE], 2, function(b) {
    x <- nnls(A / scale_f, b / scale_f)$x
    x / sum(x)
  }))
  colnames(p_nn) <- CELL_TYPES

  long <- rbind(
    data.table(method = "dtangle", sample = rownames(p_dt), as.data.table(p_dt)),
    data.table(method = "NNLS", sample = rownames(p_nn), as.data.table(p_nn))
  )
  long <- melt(long, id.vars = c("method", "sample"), variable.name = "cell_type", value.name = "proportion")
  long[, platform := label]

  agree <- sapply(CELL_TYPES, function(ct) cor(p_dt[, ct], p_nn[rownames(p_dt), ct], method = "spearman"))
  cat("[", label, "] dtangle vs NNLS Spearman rho per cell type:\n")
  print(round(agree, 2))

  # validation against single canonical marker genes (log2 CPM)
  can <- CANONICAL[gene_id %in% rownames(blk)]
  valid <- can[, .(
    NNLS = cor(p_nn[colnames(blk), cell_type], log2(blk[gene_id, ] + 1), method = "spearman"),
    dtangle = cor(p_dt[colnames(blk), cell_type], log2(blk[gene_id, ] + 1), method = "spearman")
  ), by = .(gene, cell_type)]
  valid[, platform := label]
  cat("[", label, "] Spearman rho of estimated proportion with canonical marker expression:\n")
  print(valid[, .(gene, cell_type, NNLS = round(NNLS, 2), dtangle = round(dtangle, 2))])
  list(long = long, markers = spec_markers, validation = valid)
}

lr_dec <- deconvolve(lr_bulk, "long-read")
sr_dec <- deconvolve(sr_bulk, "short-read")

lr_long <- merge(lr_dec$long, lr_meta[, .(sample, condition, sex = tolower(sex))], by = "sample")
sr_long <- merge(sr_dec$long, sr_meta[, .(sample, condition = as.character(condition), sex = as.character(sex),
                                          age, in_manuscript)], by = "sample")
props <- rbind(lr_long, sr_long, fill = TRUE)
props[, condition := factor(condition, levels = c("CT", "AD"))]
props[, cell_type := factor(cell_type, levels = CELL_TYPES)]
fwrite(props, file.path(TABLES, "cell_type_proportions.tsv"), sep = "\t")
fwrite(rbind(lr_dec$validation, sr_dec$validation), file.path(TABLES, "marker_gene_validation.tsv"), sep = "\t")

# neuronal fraction (EX + INH) per sample, used as covariate downstream
neuro <- props[cell_type %in% c("EX", "INH"), .(neuronal = sum(proportion)), by = .(platform, method, sample, condition)]
fwrite(neuro, file.path(TABLES, "neuronal_fraction.tsv"), sep = "\t")

###############################
# AD vs CT composition
###############################

test_ct <- function(d, platform) {
  if (platform == "long-read") {
    w <- wilcox.test(proportion ~ condition, data = d, exact = TRUE)
    f <- lm(proportion ~ condition + sex, data = d)
  } else {
    w <- wilcox.test(proportion ~ condition, data = d, exact = FALSE)
    f <- lm(proportion ~ condition + age + sex, data = d)
  }
  cf <- summary(f)$coefficients
  data.table(mean_CT = mean(d[condition == "CT", proportion]), mean_AD = mean(d[condition == "AD", proportion]),
             wilcox_p = w$p.value, lm_diff = cf["conditionAD", "Estimate"], lm_p = cf["conditionAD", "Pr(>|t|)"])
}
comp <- props[, test_ct(.SD, platform), by = .(platform, method, cell_type)]
comp[, lm_fdr := p.adjust(lm_p, "BH"), by = .(platform, method)]
nd <- merge(neuro, unique(props[, .(sample, sex, age)]), by = "sample")
setnames(nd, "neuronal", "proportion")
nf <- nd[, test_ct(.SD, platform), by = .(platform, method)]
nf[, cell_type := "Neuronal (EX+INH)"]
comp <- rbind(comp, nf, fill = TRUE)
fwrite(comp, file.path(TABLES, "composition_AD_vs_CT.tsv"), sep = "\t")
cat("\n==== Cell-type proportions, AD vs CT (lm: long-read ~ condition + sex; short-read ~ condition + age + sex) ====\n")
print(comp[, .(platform, method, cell_type, mean_CT = round(mean_CT, 3), mean_AD = round(mean_AD, 3),
               diff = round(lm_diff, 3), lm_p = signif(lm_p, 2), lm_fdr = signif(lm_fdr, 2),
               wilcox_p = signif(wilcox_p, 2))], nrows = 100)

###############################
# Figure
###############################

cond_cols <- c(CT = "turquoise3", AD = "salmon")
plot_platform <- function(pl, title) {
  d <- props[platform == pl & method == "NNLS"]
  ggplot(d, aes(condition, proportion, colour = condition)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.4, colour = "grey45") +
    geom_jitter(width = 0.12, height = 0, size = 1.9, alpha = 0.9) +
    facet_wrap(~ cell_type, nrow = 1, scales = "free_y") +
    scale_colour_manual(values = cond_cols, guide = "none") +
    labs(title = title, x = NULL, y = "Estimated proportion (NNLS)") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
          strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold", size = 11))
}
p <- plot_platform("long-read", "Long-read ONT (6 AD vs 6 CT)") /
  plot_platform("short-read", "Short-read E-MTAB-14167 (20 AD vs 21 CT)") +
  plot_annotation(title = "Cell-type composition estimated against GSE174367 PFC snRNA-seq")
ggsave(file.path(FIGS, "cell_type_proportions_NNLS.png"), p, width = 13, height = 7.5, dpi = 300, bg = "white")

agree_dt <- dcast(props, platform + sample + cell_type ~ method, value.var = "proportion")
p2 <- ggplot(agree_dt, aes(dtangle, NNLS)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey75", linewidth = 0.4) +
  geom_point(colour = "#2a78d6", size = 1.6, alpha = 0.8) +
  facet_grid(platform ~ cell_type, scales = "free") +
  labs(title = "Agreement between deconvolution methods", x = "dtangle", y = "NNLS") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))
ggsave(file.path(FIGS, "deconvolution_method_agreement.png"), p2, width = 14, height = 5, dpi = 300, bg = "white")

saveRDS(list(props = props, neuro = neuro, comp = comp, lr_markers = lr_dec$markers, sr_markers = sr_dec$markers),
        file.path(TABLES, "deconvolution.rds"))
cat("\nDone:", format(Sys.time()), "\n")
sink()
