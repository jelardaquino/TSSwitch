# promoter usage differential analysis using betaregression on 61 DTU significant genes

library(GenomicRanges)
library(data.table)
library(dplyr)
library(rtracklayer)
library(tidyr)
library(tibble)
library(betareg)
library(ggplot2)

gtf_file <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/reproducing_nextflow_pipeline/Homo_sapiens.GRCh38.107_ERCC.gtf"   # your GTF
tx_counts_file <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/PRJNA1206164/transcript_tpm.txt"
DTU_sig_file <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/PRJNA1206164/diu/IsoformSwitchAnalyzeR_top_switches.tsv"
DTU_61_genes <- read.table(DTU_sig_file, header = TRUE, sep = "\t")

output_dir_tables <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/PRJNA1206164/tables/"
output_dir_figures <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/PRJNA1206164/figures/"
dir.create(output_dir_tables, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir_figures, showWarnings = FALSE, recursive = TRUE)

# import GTF and extract promoter regions
g <- rtracklayer::import(gtf_file)
# parameters for promoter window (adjust as needed)
up <- 1000L
down <- 100L

# extract transcript features from the imported GTF
tx_gr <- g[g$type == "transcript"]

# locate a transcript identifier column in the GTF metadata
meta_names <- names(mcols(tx_gr))
tid_col <- meta_names[grepl("transcript_id|transcriptname|transcript_name|transcript|trans_id|tx_id", meta_names, ignore.case = TRUE)][1]
if (is.na(tid_col)) tid_col <- meta_names[1]
tx_ids <- as.character(mcols(tx_gr)[[tid_col]])

# compute TSS per transcript and create promoter ranges
strand_char <- as.character(strand(tx_gr))
tss <- ifelse(strand_char == "-", end(tx_gr), start(tx_gr))
prom_start <- pmax(1L, tss - up)
prom_end <- tss + down
prom_gr <- GRanges(seqnames = seqnames(tx_gr),
                   ranges = IRanges(prom_start, prom_end),
                   strand = strand(tx_gr))
mcols(prom_gr)$transcript_id <- tx_ids
mcols(prom_gr)$promoter_id <- paste0(seqnames(prom_gr), ":", prom_start, "-", prom_end, "(", strand_char, ")")

# read transcript-level counts
tx_counts <- fread(tx_counts_file)

# filter to only transcripts in the DTU significant genes
# map gene_id -> transcript_id from the GTF and filter counts accordingly
gene_col <- meta_names[grepl("^gene_id$|geneid|gene", meta_names, ignore.case = TRUE)][1]
if (is.na(gene_col)) stop("Could not find a gene_id column in the GTF transcript metadata.")

tx_gene_map <- data.table(
  transcript_id = tx_ids,
  gene_id = as.character(mcols(tx_gr)[[gene_col]])
)


# normalize gene_id by removing any version suffix after a dot (e.g. ENSG000001.1 -> ENSG000001)
DTU_61_genes$gene_id <- sub("\\..*$", "", as.character(DTU_61_genes$gene_id))
# keep tx->gene mappings only for the DTU genes (normalize gene ids first)

tx_gene_map <- tx_gene_map[gene_id %in% as.character(DTU_61_genes$gene_id)]
# get transcripts to retain
sig_tx <- unique(tx_gene_map$transcript_id)
# detect id column in tx_counts robustly and subset counts to those transcripts
tx_counts <- tx_counts %>% filter(V1 %in% sig_tx)


colnames(tx_counts)[1] <- "transcript_id"  # assuming first column is transcript id 
# detect transcript id column in the counts file
count_id_candidates <- c("transcript_id", "transcript", "tx_id", "transcript_name", "Name", "ID")
id_col <- intersect(names(tx_counts), count_id_candidates)[1]
if (is.na(id_col)) id_col <- names(tx_counts)[1]
tx_counts[[id_col]] <- as.character(tx_counts[[id_col]])

# remove specific samples from tx_counts (by column name, if present)
# samples_to_drop <- c("1163", "1271", "1092")
# tx_counts <- tx_counts %>%
#     select(-any_of(samples_to_drop))
# prepare promoter metadata table and merge with counts by transcript id
prom_dt <- data.table(
  promoter_id = mcols(prom_gr)$promoter_id,
  transcript_id = mcols(prom_gr)$transcript_id,
  seqnames = as.character(seqnames(prom_gr)),
  start = start(prom_gr),
  end = end(prom_gr),
  strand = as.character(strand(prom_gr))
)
setDT(prom_dt); setDT(tx_counts)
combined <- merge(prom_dt, tx_counts, by.x = "transcript_id", by.y = id_col, all.x = FALSE, all.y = FALSE)

# add gene_id column to prom_dt by joining transcript -> gene_id from the GTF
gene_id_col <- meta_names[grepl("gene_id|geneid", meta_names, ignore.case = TRUE)][1]
tx_gene_id <- if (!is.na(gene_id_col)) as.character(mcols(tx_gr)[[gene_id_col]]) else NA_character_

tx_gene_dt <- data.table(
    transcript_id = as.character(mcols(tx_gr)[[tid_col]]),
    gene_id = tx_gene_id
)

setDT(prom_dt)
prom_dt <- merge(prom_dt, tx_gene_dt, by = "transcript_id", all.x = TRUE)
# add gene_name to prom_dt (fallback to gene_id if gene_name not present/empty)
# If the currently imported GTF doesn't include gene_name, re-import an alternative GTF that does.
alt_gtf_file <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/reproducing_nextflow_pipeline/Homo_sapiens.GRCh38.107_ERCC.gtf"

g_alt <- rtracklayer::import(alt_gtf_file)
tx_alt <- g_alt[g_alt$type == "transcript"]

alt_meta <- names(mcols(tx_alt))
alt_tid_col <- alt_meta[grepl("transcript_id|transcriptname|transcript_name|transcript|trans_id|tx_id", alt_meta, ignore.case = TRUE)][1]
if (is.na(alt_tid_col)) stop("Could not find transcript_id in the alternative GTF.")

alt_gene_col <- alt_meta[grepl("^gene_id$|geneid|gene", alt_meta, ignore.case = TRUE)][1]
if (is.na(alt_gene_col)) stop("Could not find gene_id in the alternative GTF transcript metadata.")

alt_gene_name_col <- alt_meta[grepl("^gene_name$|genename|gene_name", alt_meta, ignore.case = TRUE)][1]
# gene_name might still be absent in some GTFs; we handle that below.

tx_gene_annot_dt <- data.table(
  transcript_id = as.character(mcols(tx_alt)[[alt_tid_col]]),
  gene_id_alt   = as.character(mcols(tx_alt)[[alt_gene_col]]),
  gene_name     = if (!is.na(alt_gene_name_col)) as.character(mcols(tx_alt)[[alt_gene_name_col]]) else NA_character_
)

# Keep one row per transcript_id (avoid many-to-many merges)
setkey(tx_gene_annot_dt, transcript_id)
tx_gene_annot_dt <- unique(tx_gene_annot_dt, by = "transcript_id")

# Join onto prom_dt by transcript_id
prom_dt <- merge(prom_dt, tx_gene_annot_dt, by = "transcript_id", all.x = TRUE)

# Prefer existing gene_id from your original GTF; otherwise fill from alt
prom_dt[, gene_id := fifelse(is.na(gene_id) | gene_id == "", gene_id_alt, gene_id)]
prom_dt[, gene_id_alt := NULL]

# Fill gene_name; fallback to gene_id
prom_dt[, gene_name := fifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name)]

count_cols <- names(tx_counts)[sapply(tx_counts, is.numeric)]
if (length(count_cols) == 0) stop("No numeric count columns detected in the transcript counts file.")
promoter_counts <- combined[, lapply(.SD, sum), by = .(promoter_id, seqnames, start, end, strand), .SDcols = count_cols]

# create a GRanges of promoters with aggregated counts (optional)
promoter_granges <- GRanges(
  seqnames = promoter_counts$seqnames,
  ranges = IRanges(promoter_counts$start, promoter_counts$end),
  strand = promoter_counts$strand
)
mcols(promoter_granges) <- as.data.frame(promoter_counts[, ..count_cols])
mcols(promoter_granges)$promoter_id <- promoter_counts$promoter_id

# results:
# - promoter_counts : data.table with summed counts per promoter
# - promoter_granges : GRanges with count columns in metadata
promoter_counts_dt <- promoter_counts
write.table(promoter_counts_dt, paste0(output_dir_tables, "promoter_counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


# prepare count matrix (rows = promoters, cols = samples)
counts_mat <- as.matrix(promoter_counts_dt[, ..count_cols])
rownames(counts_mat) <- promoter_counts_dt$promoter_id

sample_info_file <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/PRJNA1206164/metadata_AD_CT.txt"

    coldata <- fread(sample_info_file)
    setDF(coldata)
    rownames(coldata) <- coldata[[1]]
    coldata[[1]] <- NULL
    # drop the same samples from coldata (by rowname / sample id)
    # samples_to_drop <- c("1163", "1271", "1092")
    # coldata <- coldata[!rownames(coldata) %in% samples_to_drop, , drop = FALSE]
counts_long <- counts_mat %>%
  as.data.frame() %>%
  rownames_to_column("promoter_id") %>%
  pivot_longer(
    cols = -promoter_id,
    names_to = "sample",
    values_to = "abundance"
  )

prom_anno <- prom_dt %>%
  distinct(promoter_id, gene_id, gene_name)
# ---- clean up prom_anno: de-dup promoter_id and standardize gene_name per gene_id ----
prom_anno <- prom_anno %>%
    mutate(
        gene_id   = as.character(gene_id),
        gene_name = na_if(as.character(gene_name), "")
    )

# pick ONE "best" gene_name per gene_id: prefer a non-ENSG symbol; else any non-empty; else gene_id
gene_name_map <- prom_anno %>%
    group_by(gene_id) %>%
    summarise(
        gene_name_best = {
            vals <- unique(na.omit(gene_name))
            sym  <- vals[!grepl("^ENSG\\d+", vals)]
            if (length(sym) > 0) sym[[1]] else if (length(vals) > 0) vals[[1]] else gene_id[[1]]
        },
        .groups = "drop"
    )

# apply mapping and enforce one row per promoter_id
prom_anno <- prom_anno %>%
    left_join(gene_name_map, by = "gene_id") %>%
    mutate(gene_name = gene_name_best) %>%
    select(promoter_id, gene_id, gene_name) %>%
    distinct() %>%
    group_by(promoter_id) %>%
    arrange(gene_id, gene_name, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()

write.table(prom_anno, paste0(output_dir_tables, "promoter_annotations.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
counts_long <- counts_long %>%
  left_join(prom_anno, by = "promoter_id")

counts_long <- counts_long %>%
  left_join(
    coldata %>% 
      rownames_to_column("sample"),
    by = "sample"
  )

#compute promoter usage
# usage = promoter abundance / total gene abundance per sample
usage_long <- counts_long %>%
  group_by(gene_id, sample) %>%
  mutate(
    gene_total = sum(abundance),
    promoter_usage = abundance / gene_total
  ) %>%
  ungroup()

# filter low-information genes
usage_long <- usage_long %>%
  filter(gene_total >= 1) %>%        # tune if needed
  group_by(gene_id) %>%
  filter(n_distinct(promoter_id) >= 2) %>%
  ungroup()

#Beta regression analysis for differential promoter usage
usage_long <- usage_long %>%
  mutate(
    usage_adj = pmin(pmax(promoter_usage, 1e-4), 1 - 1e-4),
    condition = factor(condition, levels = c("CT", "AD"))
  )

write.table(usage_long, paste0(output_dir_tables, "promoter_usage_long.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
# filter promoters with low variability or missing groups
usage_filt <- usage_long %>%
  group_by(promoter_id) %>%
  filter(
    n_distinct(condition) == 2,                     # both groups present
    n_distinct(usage_adj) >= 3,                     # variability
    all(tapply(usage_adj, condition, var) > 0)      # variance in both groups
  ) %>%
  ungroup()

safe_betareg <- function(df) {
  tryCatch({
    fit <- betareg(usage_adj ~ condition, data = df)
    coef <- summary(fit)$coefficients$mean
    data.frame(
      logitFC = coef["conditionAD", "Estimate"],
      SE = coef["conditionAD", "Std. Error"],
      pvalue = coef["conditionAD", "Pr(>|z|)"]
    )
  }, error = function(e) {
    data.frame(
      logitFC = NA,
      SE = NA,
      pvalue = NA
    )
  })
}

res_usage <- usage_filt %>%
  group_by(promoter_id, gene_id, gene_name) %>%
  do(safe_betareg(.)) %>%
  ungroup() %>%
  filter(!is.na(pvalue)) %>%
  mutate(padj = p.adjust(pvalue, method = "BH"))
write.table(res_usage, paste0(output_dir_tables, "promoter_usage_betareg_results.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# true switching genes: at least one promoter upregulated and at least one downregulated
switching_genes_usage <- res_usage %>%
  filter(padj <= 0.05) %>%
  group_by(gene_id) %>%
  filter(
    any(logitFC > 0) &
    any(logitFC < 0)
  ) %>%
  ungroup()

# keep only unique genes (one row per gene) to avoid inflated gene counts downstream
switching_genes_usage_genes <- switching_genes_usage %>%
  distinct(gene_id, gene_name)

write.table(switching_genes_usage, paste0(output_dir_tables, "switching_promoter_usage.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# asymmetric usage shifts
asymmetric_usage <- res_usage %>%
  filter(padj <= 0.05) %>%
  group_by(gene_id) %>%
  filter(n_distinct(sign(logitFC)) == 1) %>%
  ungroup()
write.table(asymmetric_usage, paste0(output_dir_tables, "asymmetric_promoter_usage.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# keep only genes with >= 2 promoters AND restrict to unique genes in switching_genes_usage
switch_gene_ids <- switching_genes_usage %>%
  pull(gene_name) %>%
  unique() %>%
  as.character()

plot_df <- usage_long %>%
  filter(gene_name %in% switch_gene_ids) %>%
  group_by(gene_name) %>%
  filter(n_distinct(promoter_id) >= 2) %>%
  ungroup()

  # filter out promoters with very low adjusted usage
  plot_df <- plot_df %>%
    group_by(promoter_id) %>%
    filter(any(usage_adj >= 0.10, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(gene_name) %>%
    filter(n_distinct(promoter_id) >= 2) %>%
    ungroup()

# plot_df <- plot_df %>%
#   filter(gene_name=="NOTCH4" | gene_name=="SNCA" | gene_name=="NSD2" )


mean_df <- plot_df %>%
  group_by(gene_id, gene_name, promoter_id, condition) %>%
  summarise(
    mean_usage = mean(promoter_usage),
    .groups = "drop"
  )

p <- ggplot(
      plot_df,
      aes(
        x = promoter_id,
        y = promoter_usage,
        fill = condition
      )
      ) +
      geom_violin(
        scale = "width",
        trim = FALSE,
        alpha = 0.7
      ) +
      geom_point(
        data = mean_df,
        aes(
          x = promoter_id,
          y = mean_usage,
          group = condition
        ),
        position = position_dodge(width = 0.9),
        size = 2,
        shape = 21,
        fill = "white",
        color = "black",
        inherit.aes = FALSE
      ) +
      facet_wrap(
        ~ gene_name,
        scales = "free_x",
        # ncol = 1,
        # nrow = 3
      ) +
      labs(
        x = "Promoter",
        y = "Promoter usage",
        fill = "Condition"
      ) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"),
        panel.grid.major.x = element_blank()
      ) +
      scale_fill_manual(
      values = c("CT" = "turquoise", "AD" = "salmon"),
      breaks = c("AD", "CT")
      ) +
      ylim(0, 1)

p +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 50)
  )
ggsave(paste0(output_dir_figures, "promoter_usage_switching_genes_all.png"), width = 25, height = 22, dpi=300)
# ggsave(paste0(output_dir_figures, "promoter_usage_switching_genes_3genes.png"), width = 10, height = 5, dpi=300)
# ---- Enhanced volcano plot (res_all = res_usage; top_labels = unique gene_name from switching genes) ----
library(EnhancedVolcano)

# res_all <- as.data.frame(res_usage)

# # map promoter-level results to gene_name labels
# res_all$label <- ifelse(!is.na(res_all$gene_name) & res_all$gene_name != "",
#                         as.character(res_all$gene_name),
#                         as.character(res_all$gene_id))

# # label only switching genes (unique gene_name)
# top_labels <- unique(as.character(switching_genes_usage$gene_name))
# top_labels <- top_labels[!is.na(top_labels) & top_labels != ""]

# thresholds (adjust if desired)
p_cut <- 0.05
fc_cut <- 1

# p <- EnhancedVolcano(
#   toptable = res_all,
#   lab = res_all$label,
#   x = "logitFC",
#   y = "padj",
#   selectLab = top_labels,
#   pCutoff = p_cut,
#   FCcutoff = fc_cut,
#   xlab = bquote(~logit~" fold change"),
#   ylab = bquote(~-log[10]~" adjusted p-value"),
#   title = "Differential Promoter Usage Volcano Plot",
#   subtitle = NULL,
#   legendPosition = "right",
#   col = c("grey70", "#1f78b4", "#e31a1c", "purple"),
#   pointSize = 2.0,
#   labSize = 3.0,
#   drawConnectors = TRUE,
#   widthConnectors = 0.5
# )

# ggsave(
#   filename = "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn25671149/jaquino_analysis/tmp/promoter_usage_61_genes/promoter_usage_volcano.png",
#   plot = p,
#   width = 20,
#   height = 15
# )


res_sig <- as.data.frame(res_usage) %>%
  mutate(
    log2FC = logitFC / log(2),
    sig_dir = case_when(
      padj > p_cut | is.na(padj) ~ "NS",
      log2FC >= fc_cut           ~ "Up",
      log2FC <= -fc_cut          ~ "Down",
      TRUE                       ~ "NS"
    )
  )

keyvals <- c(
  NS   = "grey70",
  Up   = "#e31a1c",
  Down = "#1f78b4"
)

keyvals.colour <- keyvals[res_sig$sig_dir]

# thresholds (adjust if desired)
p_cut <- 0.05
fc_cut <- 1

genes_to_label <- c("PLD3", "MFN2", "CISD2", "DKK3")

label_rows <- res_sig %>%
  filter(
    gene_name %in% genes_to_label &
    sig_dir %in% c("Up", "Down")
  ) %>%
  group_by(gene_name, sig_dir) %>%
  slice_min(padj, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(promoter_id)

res_sig <- res_sig %>%
  mutate(
    label_sig = ifelse(
      promoter_id %in% label_rows$promoter_id,
      gene_name,
      NA
    )
  )

p <- EnhancedVolcano(
  toptable = res_sig,
  lab = res_sig$label_sig,
  x = "log2FC",
  y = "padj",
  pCutoff = p_cut,
  FCcutoff = fc_cut,
  colCustom = keyvals.colour,
  drawConnectors = TRUE,
  widthConnectors = 0.4,
  max.overlaps = 10,
  ylim = c(0, 15)
)

ggsave(
  filename = paste0(output_dir_figures, "promoter_usage_volcano.png"),
  plot = p,
  width = 10,
  height = 8,
  dpi = 300
)

# exploring whether there is promoter preference in AD
# asking: is promoter usage non-uniform/heterogeneous within a gene in AD samples?
# is there evidence that multiple promoters contribute significantly (and unevenly) to transcription in AD?
usage_AD <- usage_long %>%
  filter(condition == "AD") %>%
  group_by(gene_id, promoter_id) %>%
  summarise(
    total_abundance = sum(abundance),
    .groups = "drop"
  )


# gene-level Dirichlet-multinomial test
library(DirichletMultinomial)
library(dirmult)
uniform_test <- usage_AD %>%
  group_by(gene_id) %>%
  filter(n_distinct(promoter_id) >= 2) %>%
  do({
    counts <- .$total_abundance
    counts <- counts / sum(counts) * 1e6  # scale to pseudo-counts
    
    expected <- rep(mean(counts), length(counts))
    
    stat <- sum((counts - expected)^2 / expected)
    pval <- pchisq(stat, df = length(counts) - 1, lower.tail = FALSE)
    
    data.frame(
      chisq = stat,
      pvalue = pval
    )
  }) %>%
  ungroup() %>%
  mutate(padj = p.adjust(pvalue, "BH"))

heterogeneous_AD <- uniform_test %>%
  filter(padj <= 0.05)

sig_genes <- heterogeneous_AD$gene_id

plot_df <- usage_long %>%
  filter(
    condition == "AD",
    gene_id %in% sig_genes
  )

ggplot(plot_df, aes(x = promoter_id, y = usage_adj)) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    fill = "steelblue",
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.1,
    size = 0.8,
    alpha = 0.6
  ) +
  facet_wrap(
    ~ gene_id,
    scales = "free_x",
    ncol = 5
  ) +
  labs(
    x = "Promoter",
    y = "Promoter usage (PSI)",
    title = "Within-AD promoter usage heterogeneity for significant genes"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(paste0(output_dir_figures, "promoter_usage_heterogeneous_AD.png"), width = 50, height = 50, limitsize = FALSE)
# the violin plots shows: Within AD samples, promoter usage is non-uniform for a subset of genes. Some genes exhibit strong dominance of a single promoter, while others retain usage of multiple promoters, consistent with heterogeneous isoform expression in disease

