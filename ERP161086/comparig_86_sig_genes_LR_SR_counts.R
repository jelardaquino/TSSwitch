

library(dplyr)
library(data.table)
library(rtracklayer)
library(ggpattern)
sr_tpm_counts <- read.table("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/results/transcript_tpm.txt", header = TRUE, row.names = 1, sep="\t")
sr_metadata <- read.table("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/metadata.txt", header = TRUE, sep="\t")
lr_tpm_counts <- read.table("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis_v3/bambu_tx_tpm.txt", header = TRUE, row.names = 1, sep="\t")
lr_metadata <- read.table("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis_v3/metadata_pipeline.txt", header = TRUE, sep="\t")
lr_metadata <- lr_metadata %>% mutate(Run = paste0("X", Run))

sig_genes <- read.table("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis_v2/tables/diu_sig_genes.txt", header=FALSE)

gtf_file <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis_v3/NanoporeRNASeq_extended_annotations.gtf"
gtf <- rtracklayer::import(gtf_file)


library(dplyr)
library(tibble)
library(ggplot2)

# ------------------------------------------------------------
# 1. Clean inputs
# ------------------------------------------------------------

sig <- unique(as.character(sig_genes[[1]]))

lr_counts <- as.matrix(lr_tpm_counts)
sr_counts <- as.matrix(sr_tpm_counts)

storage.mode(lr_counts) <- "numeric"
storage.mode(sr_counts) <- "numeric"

# Standardize metadata column names
lr_meta <- lr_metadata
sr_meta <- sr_metadata

if ("AGE" %in% names(sr_meta) && !"age" %in% names(sr_meta)) {
  names(sr_meta)[names(sr_meta) == "AGE"] <- "age"
}

lr_meta <- lr_meta %>%
  mutate(
    Run = as.character(Run),
    condition = as.character(condition),
    sex = as.character(sex),
    age = as.numeric(age)
  ) %>%
  filter(Run %in% colnames(lr_counts)) %>%
  arrange(match(Run, colnames(lr_counts)))

sr_meta <- sr_meta %>%
  mutate(
    Run = as.character(Run),
    condition = as.character(condition),
    sex = as.character(sex),
    age = as.numeric(age)
  ) %>%
  filter(Run %in% colnames(sr_counts)) %>%
  arrange(match(Run, colnames(sr_counts)))

lr_counts <- lr_counts[, lr_meta$Run, drop = FALSE]
sr_counts <- sr_counts[, sr_meta$Run, drop = FALSE]

# ------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------

strip_version <- function(x) {
  sub("\\.\\d+$", "", as.character(x))
}

make_tx2gene <- function(gtf, count_ids, preferred = "transcript_id") {
  count_ids <- as.character(count_ids)
  count_key <- strip_version(count_ids)

  # If the count matrix is already gene-level
  if (mean(grepl("^ENSG", count_key)) > 0.8) {
    return(tibble(tx_id = count_ids, gene_id = count_key))
  }

  gtf_tx <- gtf[mcols(gtf)$type == "transcript"]
  mc <- as.data.frame(mcols(gtf_tx))

  if (!"gene_id" %in% names(mc)) {
    stop("No gene_id column found in the GTF metadata.")
  }

  candidate_cols <- setdiff(names(mc), "gene_id")

  overlaps <- sapply(candidate_cols, function(nm) {
    ids <- strip_version(mc[[nm]])
    ids <- ids[!is.na(ids) & ids != ""]
    length(intersect(unique(ids), count_key))
  })

  overlaps <- sort(overlaps, decreasing = TRUE)

  message("Top GTF columns overlapping count matrix rownames:")
  print(head(data.frame(
    gtf_column = names(overlaps),
    n_overlap = as.integer(overlaps)
  ), 10))

  positive_cols <- names(overlaps)[overlaps > 0]

  if (length(positive_cols) == 0) {
    stop(
      "No GTF metadata column overlaps the count matrix rownames. ",
      "Your count rownames may be Bambu internal IDs that need a separate rowData/transcript mapping."
    )
  }

  tx_col <- if (preferred %in% positive_cols) preferred else positive_cols[1]

  message("Using GTF column for transcript IDs: ", tx_col)

  tx_map <- tibble(
    tx_key = strip_version(mc[[tx_col]]),
    gene_id = strip_version(mc$gene_id)
  ) %>%
    filter(
      !is.na(tx_key), tx_key != "",
      !is.na(gene_id), gene_id != ""
    ) %>%
    distinct()

  # Drop ambiguous transcript IDs that map to multiple genes
  tx_map <- tx_map %>%
    add_count(tx_key, name = "n_genes") %>%
    filter(n_genes == 1) %>%
    select(-n_genes)

  tibble(
    tx_id = count_ids,
    tx_key = count_key
  ) %>%
    inner_join(tx_map, by = "tx_key") %>%
    select(tx_id, gene_id) %>%
    distinct()
}

aggregate_tpm_to_gene <- function(tpm, tx2gene) {
  tx2gene <- tx2gene %>%
    filter(tx_id %in% rownames(tpm)) %>%
    distinct(tx_id, gene_id)

  m <- tpm[tx2gene$tx_id, , drop = FALSE]

  gene_tpm <- rowsum(
    m,
    group = tx2gene$gene_id,
    reorder = FALSE
  )

  as.matrix(gene_tpm)
}

# ------------------------------------------------------------
# 3. Transcript TPM -> gene TPM
# ------------------------------------------------------------

lr_tx2gene <- make_tx2gene(gtf, rownames(lr_counts), preferred = "transcript_id")
sr_tx2gene <- make_tx2gene(gtf, rownames(sr_counts), preferred = "transcript_id")

lr_gene_tpm <- aggregate_tpm_to_gene(lr_counts, lr_tx2gene)
sr_gene_tpm <- aggregate_tpm_to_gene(sr_counts, sr_tx2gene)

cat("LR genes:", nrow(lr_gene_tpm), "\n")
cat("SR genes:", nrow(sr_gene_tpm), "\n")
cat("Significant genes:", length(sig), "\n")

common_sig <- Reduce(intersect, list(
  sig,
  rownames(lr_gene_tpm),
  rownames(sr_gene_tpm)
))

cat("Significant genes found in both LR and SR:", length(common_sig), "\n")

lr_sig <- lr_gene_tpm[common_sig, , drop = FALSE]
sr_sig <- sr_gene_tpm[common_sig, , drop = FALSE]


# ------------------------------------------------------------
# 4. Compare condition-level expression profiles
# ------------------------------------------------------------

mean_log_by_condition <- function(mat, meta) {
  stopifnot(all(colnames(mat) == meta$Run))

  conds <- sort(unique(as.character(meta$condition)))

  out <- sapply(conds, function(cond) {
    rowMeans(log2(mat[, meta$condition == cond, drop = FALSE] + 1), na.rm = TRUE)
  })

  if (is.null(dim(out))) {
    out <- matrix(out, ncol = length(conds))
    colnames(out) <- conds
    rownames(out) <- rownames(mat)
  }

  as.matrix(out)
}

safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)

  if (sum(ok) < 3) return(NA_real_)
  if (sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)

  cor(x[ok], y[ok], method = method)
}

lr_cond_mean <- mean_log_by_condition(lr_sig, lr_meta)
sr_cond_mean <- mean_log_by_condition(sr_sig, sr_meta)

common_conditions <- intersect(colnames(lr_cond_mean), colnames(sr_cond_mean))

condition_similarity <- bind_rows(lapply(common_conditions, function(cond) {
  tibble(
    condition = cond,
    n_genes = length(common_sig),
    spearman = safe_cor(lr_cond_mean[, cond], sr_cond_mean[, cond], "spearman"),
    pearson = safe_cor(lr_cond_mean[, cond], sr_cond_mean[, cond], "pearson")
  )
}))

condition_similarity

# ------------------------------------------------------------
# 5. Compare LR vs SR log2 fold-change direction
# ------------------------------------------------------------

case_label <- "AD"
control_label <- "CT"

log2fc_by_condition <- function(mat, meta, case_label, control_label) {
  conds <- unique(as.character(meta$condition))

  if (!all(c(case_label, control_label) %in% conds)) {
    stop("Could not find both case and control labels in metadata.")
  }

  case_samples <- meta$Run[meta$condition == case_label]
  control_samples <- meta$Run[meta$condition == control_label]

  rowMeans(log2(mat[, case_samples, drop = FALSE] + 1), na.rm = TRUE) -
    rowMeans(log2(mat[, control_samples, drop = FALSE] + 1), na.rm = TRUE)
}

lr_log2fc <- log2fc_by_condition(lr_sig, lr_meta, case_label, control_label)
sr_log2fc <- log2fc_by_condition(sr_sig, sr_meta, case_label, control_label)

fc_similarity <- tibble(
  gene_id = common_sig,
  lr_log2fc = lr_log2fc[common_sig],
  sr_log2fc = sr_log2fc[common_sig],
  same_direction = sign(lr_log2fc[common_sig]) == sign(sr_log2fc[common_sig]),
  abs_delta = abs(lr_log2fc[common_sig] - sr_log2fc[common_sig])
) %>%
  arrange(desc(same_direction), abs_delta)

summary_stats <- tibble(
  n_common_sig_genes = nrow(fc_similarity),
  spearman_log2fc = safe_cor(fc_similarity$lr_log2fc, fc_similarity$sr_log2fc, "spearman"),
  pearson_log2fc = safe_cor(fc_similarity$lr_log2fc, fc_similarity$sr_log2fc, "pearson"),
  percent_same_direction = mean(fc_similarity$same_direction, na.rm = TRUE) * 100
)

summary_stats

top_concordant <- fc_similarity %>%
  filter(same_direction) %>%
  arrange(abs_delta) %>%
  head(25)

top_concordant

top_discordant <- fc_similarity %>%
  filter(!same_direction) %>%
  arrange(desc(abs_delta)) %>%
  head(25)

top_discordant

ggplot(fc_similarity, aes(x = sr_log2fc, y = lr_log2fc)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = same_direction), size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  labs(
    title = "Short-read vs long-read similarity among significant genes",
    subtitle = paste0(case_label, " vs ", control_label),
    x = "Short-read log2FC",
    y = "Long-read log2FC",
    shape = "Same direction"
  )

expr_similarity <- tibble(
  gene_id = common_sig,
  lr_mean_log_tpm = rowMeans(log2(lr_sig + 1), na.rm = TRUE),
  sr_mean_log_tpm = rowMeans(log2(sr_sig + 1), na.rm = TRUE)
)

expr_cor <- tibble(
  n_common_sig_genes = nrow(expr_similarity),
  spearman_mean_expression = safe_cor(expr_similarity$lr_mean_log_tpm, expr_similarity$sr_mean_log_tpm, "spearman"),
  pearson_mean_expression = safe_cor(expr_similarity$lr_mean_log_tpm, expr_similarity$sr_mean_log_tpm, "pearson")
)

expr_cor

ggplot(expr_similarity, aes(x = sr_mean_log_tpm, y = lr_mean_log_tpm)) +
  geom_point(alpha = 0.8, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  labs(
    title = "Mean expression similarity among significant genes",
    x = "Short-read mean log2(TPM + 1)",
    y = "Long-read mean log2(TPM + 1)"
  )

write.csv(condition_similarity, "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/results/compare_longread_shortread/lr_sr_condition_similarity_sig_genes.csv", row.names = FALSE)
write.csv(summary_stats, "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/results/compare_longread_shortread/lr_sr_log2fc_similarity_summary.csv", row.names = FALSE)
write.csv(fc_similarity, "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/results/compare_longread_shortread/lr_sr_gene_level_log2fc_similarity.csv", row.names = FALSE)
write.csv(expr_similarity, "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/results/compare_longread_shortread/lr_sr_mean_expression_similarity.csv", row.names = FALSE)



###################################################################################################################################################
# Matching by isoform usage

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

# ------------------------------------------------------------
# 1. Basic cleanup
# ------------------------------------------------------------

strip_version <- function(x) {
  sub("\\.\\d+$", "", as.character(x))
}

safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)

  if (sum(ok) < 3) return(NA_real_)
  if (sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)

  cor(x[ok], y[ok], method = method)
}

sig <- unique(strip_version(sig_genes[[1]]))

lr_counts <- as.matrix(lr_tpm_counts)
sr_counts <- as.matrix(sr_tpm_counts)

storage.mode(lr_counts) <- "numeric"
storage.mode(sr_counts) <- "numeric"

lr_meta <- lr_metadata
sr_meta <- sr_metadata

if ("AGE" %in% names(sr_meta) && !"age" %in% names(sr_meta)) {
  names(sr_meta)[names(sr_meta) == "AGE"] <- "age"
}

lr_meta <- lr_meta %>%
  mutate(
    Run = as.character(Run),
    condition = as.character(condition),
    sex = as.character(sex),
    age = as.numeric(age)
  ) %>%
  filter(Run %in% colnames(lr_counts)) %>%
  arrange(match(Run, colnames(lr_counts)))

sr_meta <- sr_meta %>%
  mutate(
    Run = as.character(Run),
    condition = as.character(condition),
    sex = as.character(sex),
    age = as.numeric(age)
  ) %>%
  filter(Run %in% colnames(sr_counts)) %>%
  arrange(match(Run, colnames(sr_counts)))

lr_counts <- lr_counts[, lr_meta$Run, drop = FALSE]
sr_counts <- sr_counts[, sr_meta$Run, drop = FALSE]

# ------------------------------------------------------------
# 2. Build transcript-to-gene table from GTF
# ------------------------------------------------------------

# ------------------------------------------------------------
# Build transcript table from GTF
# ------------------------------------------------------------

gtf_tx <- gtf[as.character(mcols(gtf)$type) == "transcript"]

gtf_tx_tbl <- as.data.frame(mcols(gtf_tx)) %>%
  mutate(
    tx_index = row_number(),
    gene_id = strip_version(gene_id),
    transcript_id = strip_version(transcript_id)
  ) %>%
  dplyr::select(tx_index, gene_id, transcript_id) %>%
  filter(
    !is.na(gene_id), gene_id != "",
    !is.na(transcript_id), transcript_id != ""
  )

head(gtf_tx_tbl, 20)

# ------------------------------------------------------------
# 3. Map TPM row names to gene_id and transcript_id
# ------------------------------------------------------------

map_tpm_rows_to_gtf <- function(tpm, gtf_tx_tbl, label, allow_bambu_index_fallback = TRUE) {
  row_ids <- rownames(tpm)
  row_key <- strip_version(row_ids)

  # ------------------------------------------------------------
  # 1. Direct transcript_id match
  #    Important: after joining row_key = transcript_id,
  #    dplyr keeps row_key, not transcript_id, so we restore it.
  # ------------------------------------------------------------

  direct_tx <- tibble(
    row_id = row_ids,
    row_key = row_key
  ) %>%
    inner_join(
      gtf_tx_tbl %>% dplyr::select(gene_id, transcript_id),
      by = c("row_key" = "transcript_id")
    ) %>%
    mutate(
      transcript_id = row_key,
      match_type = "transcript_id"
    ) %>%
    dplyr::select(row_id, gene_id, transcript_id, match_type)

  out <- direct_tx %>%
    distinct(row_id, gene_id, transcript_id, .keep_all = TRUE)

  # ------------------------------------------------------------
  # 2. BambuTx1/BambuTx2 fallback for LR
  # ------------------------------------------------------------

  if (nrow(out) == 0 && allow_bambu_index_fallback) {
    if (all(grepl("^BambuTx[0-9]+$", row_ids))) {
      idx <- as.integer(sub("^BambuTx", "", row_ids))

      if (max(idx, na.rm = TRUE) <= nrow(gtf_tx_tbl)) {
        message(label, ": using BambuTx numeric index fallback.")

        out <- tibble(
          row_id = row_ids,
          tx_index = idx
        ) %>%
          left_join(gtf_tx_tbl, by = "tx_index") %>%
          filter(!is.na(gene_id), !is.na(transcript_id)) %>%
          dplyr::select(row_id, gene_id, transcript_id) %>%
          mutate(match_type = "BambuTx_index")
      }
    }
  }

  message(label, " rows in TPM matrix: ", nrow(tpm))
  message(label, " rows mapped to genes: ", nrow(out))
  message(label, " unique genes mapped: ", length(unique(out$gene_id)))
  message(label, " unique transcript IDs mapped: ", length(unique(out$transcript_id)))

  if (nrow(out) == 0) {
    stop(
      label,
      " could not be mapped. For LR BambuTx IDs, you may need the Bambu transcript metadata ",
      "that explicitly maps BambuTx IDs to transcript_id and gene_id."
    )
  }

  out
}

lr_tx_map <- map_tpm_rows_to_gtf(lr_tpm_counts, gtf_tx_tbl, "LR")
sr_tx_map <- map_tpm_rows_to_gtf(sr_tpm_counts, gtf_tx_tbl, "SR")

head(lr_tx_map, 20)
head(sr_tx_map, 20)

table(lr_tx_map$match_type)
table(sr_tx_map$match_type)

sig <- unique(strip_version(sig_genes[[1]]))

lr_sig_tx_map <- lr_tx_map %>%
  filter(gene_id %in% sig)

sr_sig_tx_map <- sr_tx_map %>%
  filter(gene_id %in% sig)

cat("LR significant-gene transcripts:", nrow(lr_sig_tx_map), "\n")
cat("SR significant-gene transcripts:", nrow(sr_sig_tx_map), "\n")

cat("LR significant genes represented:",
    length(unique(lr_sig_tx_map$gene_id)), "\n")

cat("SR significant genes represented:",
    length(unique(sr_sig_tx_map$gene_id)), "\n")

cat("Significant genes represented in both:",
    length(intersect(lr_sig_tx_map$gene_id, sr_sig_tx_map$gene_id)), "\n")

# ------------------------------------------------------------
# 4. Keep only transcripts from significant genes
# ------------------------------------------------------------

lr_sig_tx_map <- lr_tx_map %>%
  filter(gene_id %in% sig)

sr_sig_tx_map <- sr_tx_map %>%
  filter(gene_id %in% sig)

cat("LR significant-gene transcripts:", nrow(lr_sig_tx_map), "\n")
cat("SR significant-gene transcripts:", nrow(sr_sig_tx_map), "\n")

common_sig_genes_with_tx <- intersect(lr_sig_tx_map$gene_id, sr_sig_tx_map$gene_id)

cat("Significant genes with transcripts in both LR and SR:",
    length(common_sig_genes_with_tx), "\n")

# ------------------------------------------------------------
# 5. Collapse duplicate transcript IDs if needed
# ------------------------------------------------------------

collapse_to_transcript_id <- function(tpm, tx_map) {
  tx_map <- tx_map %>%
    filter(row_id %in% rownames(tpm)) %>%
    filter(!is.na(transcript_id), transcript_id != "") %>%
    distinct(row_id, transcript_id, gene_id)

  m <- tpm[tx_map$row_id, , drop = FALSE]

  tx_tpm <- rowsum(
    m,
    group = tx_map$transcript_id,
    reorder = FALSE
  )

  tx_gene_map <- tx_map %>%
    distinct(transcript_id, gene_id)

  list(
    tpm = as.matrix(tx_tpm),
    map = tx_gene_map
  )
}

lr_tx_obj <- collapse_to_transcript_id(lr_counts, lr_sig_tx_map)
sr_tx_obj <- collapse_to_transcript_id(sr_counts, sr_sig_tx_map)

lr_sig_tx_tpm <- lr_tx_obj$tpm
sr_sig_tx_tpm <- sr_tx_obj$tpm

tx_gene_map <- bind_rows(lr_tx_obj$map, sr_tx_obj$map) %>%
  distinct(transcript_id, gene_id)

common_tx <- intersect(rownames(lr_sig_tx_tpm), rownames(sr_sig_tx_tpm))

cat("Shared annotated transcripts from significant genes:",
    length(common_tx), "\n")


# ------------------------------------------------------------
# 6. Transcript-level expression similarity by condition
# ------------------------------------------------------------

condition_mean_tpm <- function(mat, meta) {
  stopifnot(all(colnames(mat) == meta$Run))

  conds <- sort(unique(as.character(meta$condition)))

  out <- sapply(conds, function(cond) {
    rowMeans(mat[, meta$condition == cond, drop = FALSE], na.rm = TRUE)
  })

  if (is.null(dim(out))) {
    out <- matrix(out, ncol = length(conds))
    colnames(out) <- conds
    rownames(out) <- rownames(mat)
  }

  as.matrix(out)
}

lr_tx_cond_tpm <- condition_mean_tpm(lr_sig_tx_tpm, lr_meta)
sr_tx_cond_tpm <- condition_mean_tpm(sr_sig_tx_tpm, sr_meta)

common_conditions <- intersect(colnames(lr_tx_cond_tpm), colnames(sr_tx_cond_tpm))

tx_expression_similarity <- bind_rows(lapply(common_conditions, function(cond) {
  lr_vec <- log2(lr_tx_cond_tpm[common_tx, cond] + 1)
  sr_vec <- log2(sr_tx_cond_tpm[common_tx, cond] + 1)

  tibble(
    condition = cond,
    n_shared_transcripts = length(common_tx),
    spearman_log_tpm = safe_cor(lr_vec, sr_vec, "spearman"),
    pearson_log_tpm = safe_cor(lr_vec, sr_vec, "pearson")
  )
}))

tx_expression_similarity

# ------------------------------------------------------------
# 7. Isoform usage similarity per significant gene
# ------------------------------------------------------------

js_similarity <- function(p, q) {
  p <- as.numeric(p)
  q <- as.numeric(q)

  if (sum(p, na.rm = TRUE) == 0 || sum(q, na.rm = TRUE) == 0) {
    return(NA_real_)
  }

  p <- p / sum(p)
  q <- q / sum(q)

  m <- 0.5 * (p + q)

  kl <- function(a, b) {
    ok <- a > 0 & b > 0 & is.finite(a) & is.finite(b)
    sum(a[ok] * log2(a[ok] / b[ok]))
  }

  jsd <- 0.5 * kl(p, m) + 0.5 * kl(q, m)

  # 1 = identical, 0 = maximally different
  1 - sqrt(jsd)
}

cosine_similarity <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)

  denom <- sqrt(sum(x^2)) * sqrt(sum(y^2))

  if (!is.finite(denom) || denom == 0) return(NA_real_)

  sum(x * y) / denom
}

usage_similarity_one_condition <- function(cond) {
  genes <- sort(unique(tx_gene_map$gene_id))
  genes <- intersect(genes, common_sig_genes_with_tx)

  bind_rows(lapply(genes, function(g) {
    gene_txs <- tx_gene_map %>%
      filter(gene_id == g) %>%
      pull(transcript_id) %>%
      unique()

    lr_vec <- rep(0, length(gene_txs))
    sr_vec <- rep(0, length(gene_txs))

    names(lr_vec) <- gene_txs
    names(sr_vec) <- gene_txs

    lr_present <- intersect(gene_txs, rownames(lr_tx_cond_tpm))
    sr_present <- intersect(gene_txs, rownames(sr_tx_cond_tpm))

    lr_vec[lr_present] <- lr_tx_cond_tpm[lr_present, cond]
    sr_vec[sr_present] <- sr_tx_cond_tpm[sr_present, cond]

    lr_total <- sum(lr_vec, na.rm = TRUE)
    sr_total <- sum(sr_vec, na.rm = TRUE)

    if (lr_total > 0) {
      lr_usage <- lr_vec / lr_total
    } else {
      lr_usage <- rep(NA_real_, length(lr_vec))
    }

    if (sr_total > 0) {
      sr_usage <- sr_vec / sr_total
    } else {
      sr_usage <- rep(NA_real_, length(sr_vec))
    }

    lr_top_tx <- if (lr_total > 0) names(which.max(lr_usage)) else NA_character_
    sr_top_tx <- if (sr_total > 0) names(which.max(sr_usage)) else NA_character_

    tibble(
      condition = cond,
      gene_id = g,
      n_union_transcripts = length(gene_txs),
      n_lr_detected = sum(lr_vec > 0, na.rm = TRUE),
      n_sr_detected = sum(sr_vec > 0, na.rm = TRUE),
      n_shared_detected = sum(lr_vec > 0 & sr_vec > 0, na.rm = TRUE),
      lr_gene_tpm = lr_total,
      sr_gene_tpm = sr_total,
      lr_top_tx = lr_top_tx,
      sr_top_tx = sr_top_tx,
      same_top_tx = lr_top_tx == sr_top_tx,
      usage_cosine = cosine_similarity(lr_usage, sr_usage),
      usage_js_similarity = js_similarity(lr_usage, sr_usage),
      usage_spearman = safe_cor(lr_usage, sr_usage, "spearman"),
      usage_pearson = safe_cor(lr_usage, sr_usage, "pearson")
    )
  }))
}

usage_similarity <- bind_rows(lapply(common_conditions, usage_similarity_one_condition))

usage_similarity <- usage_similarity %>%
  arrange(condition, desc(usage_cosine), desc(usage_js_similarity))

usage_similarity <- merge(usage_similarity, gene_annot, by= "gene_id")

head(usage_similarity, 25)

usage_summary <- usage_similarity %>%
  group_by(condition) %>%
  summarise(
    n_genes = n(),
    median_usage_cosine = median(usage_cosine, na.rm = TRUE),
    mean_usage_cosine = mean(usage_cosine, na.rm = TRUE),
    median_js_similarity = median(usage_js_similarity, na.rm = TRUE),
    mean_js_similarity = mean(usage_js_similarity, na.rm = TRUE),
    percent_same_top_tx = mean(same_top_tx, na.rm = TRUE) * 100,
    median_shared_detected = median(n_shared_detected, na.rm = TRUE),
    .groups = "drop"
  )

usage_summary

top_similar_genes <- usage_similarity %>%
  filter(!is.na(usage_cosine)) %>%
  arrange(desc(usage_cosine), desc(usage_js_similarity)) %>%
  dplyr::select(
    condition,
    gene_id,
    n_union_transcripts,
    n_lr_detected,
    n_sr_detected,
    n_shared_detected,
    same_top_tx,
    lr_top_tx,
    sr_top_tx,
    usage_cosine,
    usage_js_similarity
  ) %>%
  head(50)

top_similar_genes

top_different_genes <- usage_similarity %>%
  filter(!is.na(usage_cosine)) %>%
  arrange(usage_cosine, usage_js_similarity) %>%
  dplyr::select(
    condition,
    gene_id,
    n_union_transcripts,
    n_lr_detected,
    n_sr_detected,
    n_shared_detected,
    same_top_tx,
    lr_top_tx,
    sr_top_tx,
    usage_cosine,
    usage_js_similarity
  ) %>%
  head(50)

top_different_genes

plot_tx_expr <- bind_rows(lapply(common_conditions, function(cond) {
  tibble(
    condition = cond,
    transcript_id = common_tx,
    gene_id = tx_gene_map$gene_id[match(common_tx, tx_gene_map$transcript_id)],
    lr_log_tpm = log2(lr_tx_cond_tpm[common_tx, cond] + 1),
    sr_log_tpm = log2(sr_tx_cond_tpm[common_tx, cond] + 1)
  )
}))

ggplot(plot_tx_expr, aes(x = sr_log_tpm, y = lr_log_tpm)) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ condition) +
  theme_bw() +
  labs(
    title = "Transcript TPM similarity for significant genes",
    x = "Short-read log2(TPM + 1)",
    y = "Long-read log2(TPM + 1)"
  )
ggsave("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/results/compare_longread_shortread/sr_lr_sig_gene_transcript_expression_similarity.png", dpi = 300, width = 6, height = 5)

ggplot(usage_similarity, aes(x = usage_cosine)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ condition) +
  theme_bw() +
  labs(
    title = "Isoform usage similarity between SR and LR",
    x = "Usage cosine similarity",
    y = "Number of significant genes"
  )
ggsave("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/results/compare_longread_shortread/sr_lr_sig_gene_isoform_usage_similarity_histogram.png", dpi = 300, width = 6, height = 5)
write.csv(tx_expression_similarity,
          "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/results/compare_longread_shortread/sr_lr_sig_gene_transcript_expression_similarity.csv",
          row.names = FALSE)

write.csv(usage_summary,
          "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/results/compare_longread_shortread/sr_lr_sig_gene_isoform_usage_summary.csv",
          row.names = FALSE)

write.csv(usage_similarity,
          "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/results/compare_longread_shortread/sr_lr_sig_gene_isoform_usage_similarity_by_gene.csv",
          row.names = FALSE)

write.csv(top_similar_genes,
          "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/results/compare_longread_shortread/sr_lr_sig_gene_top_similar_isoform_usage.csv",
          row.names = FALSE)

write.csv(top_different_genes,
          "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/results/compare_longread_shortread/sr_lr_sig_gene_top_different_isoform_usage.csv",
          row.names = FALSE)


#########################################################
#plotting
#########################################################

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

case_label <- "AD"
control_label <- "CT"
top_n <- 10

# Avoid selecting transcripts/genes that look "similar" only because both are zero
min_mean_tpm <- 0.5

mean_by_condition <- function(mat, meta, platform) {
  mat <- as.matrix(mat)

  if (is.null(rownames(mat))) {
    stop("Input matrix has no rownames. Need gene_id/transcript_id rownames.")
  }

  meta <- meta %>%
    mutate(Run = as.character(Run)) %>%
    filter(Run %in% colnames(mat)) %>%
    arrange(match(Run, colnames(mat)))

  mat <- mat[, meta$Run, drop = FALSE]

  conds <- sort(unique(as.character(meta$condition)))

  bind_rows(lapply(conds, function(cond) {
    samples <- meta$Run[as.character(meta$condition) == cond]

    tibble(
      feature_id = rownames(mat),
      condition = cond,
      mean_tpm = rowMeans(mat[, samples, drop = FALSE], na.rm = TRUE),
      mean_log2_tpm = rowMeans(log2(mat[, samples, drop = FALSE] + 1), na.rm = TRUE),
      platform = platform
    )
  }))
}

tx_to_gene_tpm <- function(tx_tpm, tx_gene_map) {
  map <- tx_gene_map %>%
    filter(transcript_id %in% rownames(tx_tpm)) %>%
    distinct(transcript_id, gene_id)

  m <- tx_tpm[map$transcript_id, , drop = FALSE]

  gene_tpm <- rowsum(
    m,
    group = map$gene_id,
    reorder = FALSE
  )

  as.matrix(gene_tpm)
}

# ------------------------------------------------------------
# Pick top similar genes
# ------------------------------------------------------------

top_similar_gene_ids <- usage_similarity %>%
  filter(!is.na(usage_cosine)) %>%
  group_by(gene_id) %>%
  summarise(
    mean_usage_cosine = mean(usage_cosine, na.rm = TRUE),
    mean_js_similarity = mean(usage_js_similarity, na.rm = TRUE),
    mean_lr_gene_tpm = mean(lr_gene_tpm, na.rm = TRUE),
    mean_sr_gene_tpm = mean(sr_gene_tpm, na.rm = TRUE),
    same_top_tx_rate = mean(same_top_tx, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(mean_lr_gene_tpm >= min_mean_tpm | mean_sr_gene_tpm >= min_mean_tpm) %>%
  arrange(desc(mean_usage_cosine), desc(mean_js_similarity)) %>%
  slice_head(n = top_n) %>%
  pull(gene_id)

top_similar_gene_ids

# ------------------------------------------------------------
# Gene-level TPM from transcript TPM
# ------------------------------------------------------------
library(org.Hs.eg.db)
library(AnnotationDbi)

all_gene_ids <- unique(c(
  rownames(lr_gene_tpm_from_tx),
  rownames(sr_gene_tpm_from_tx),
  tx_gene_map$gene_id
))

all_gene_ids <- strip_version(all_gene_ids)

gene_annot <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = all_gene_ids,
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "GENENAME")
) %>%
  as_tibble() %>%
  dplyr::rename(
    gene_id = ENSEMBL,
    gene_symbol = SYMBOL,
    gene_name = GENENAME
  ) %>%
  mutate(
    gene_id = strip_version(gene_id),
    gene_symbol = ifelse(is.na(gene_symbol) | gene_symbol == "", gene_id, gene_symbol),
    gene_label = paste0(gene_symbol, " (", gene_id, ")")
  ) %>%
  distinct(gene_id, .keep_all = TRUE)

head(gene_annot)

lr_gene_tpm_from_tx <- tx_to_gene_tpm(lr_sig_tx_tpm, tx_gene_map)
sr_gene_tpm_from_tx <- tx_to_gene_tpm(sr_sig_tx_tpm, tx_gene_map)

common_top_gene_ids <- intersect(
  top_similar_gene_ids,
  intersect(rownames(lr_gene_tpm_from_tx), rownames(sr_gene_tpm_from_tx))
)

lr_gene_plot_tpm <- lr_gene_tpm_from_tx[common_top_gene_ids, , drop = FALSE]
sr_gene_plot_tpm <- sr_gene_tpm_from_tx[common_top_gene_ids, , drop = FALSE]

gene_label_levels <- tibble(
  gene_id = common_top_gene_ids
) %>%
  left_join(gene_annot, by = "gene_id") %>%
  mutate(
    gene_label = ifelse(is.na(gene_label), gene_id, gene_label)
  )

gene_expr_plot_df <- bind_rows(
  mean_by_condition(lr_gene_plot_tpm, lr_meta, "Long-read"),
  mean_by_condition(sr_gene_plot_tpm, sr_meta, "Short-read")
) %>%
  dplyr::rename(gene_id = feature_id) %>%
  left_join(gene_annot, by = "gene_id") %>%
  mutate(
    gene_label = ifelse(is.na(gene_label), gene_id, gene_label)
  ) %>%
  filter(condition %in% c(control_label, case_label)) %>%
  mutate(
    condition = factor(condition, levels = c(control_label, case_label)),
    platform_condition = paste(platform, condition, sep = " - "),
    gene_label = factor(
      gene_label,
      levels = rev(gene_label_levels$gene_label)
    )
  )

p_top_similar_genes <- ggplot(
  gene_expr_plot_df,
  aes(x = gene_label, y = mean_log2_tpm, fill = platform_condition)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Top similar significant genes: AD vs CT expression",
    subtitle = "Genes ranked by long-read vs short-read transcript usage similarity",
    x = "Gene",
    y = "Mean log2(TPM + 1)",
    fill = "Platform / condition"
  )

p_top_similar_genes
ggsave(
  "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/results/compare_longread_shortread/top_similar_genes_AD_CT_expression_barplot.png",
  p_top_similar_genes,
  width = 11,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------
# Condition means for all shared significant-gene transcripts
# ------------------------------------------------------------

common_tx <- intersect(rownames(lr_sig_tx_tpm), rownames(sr_sig_tx_tpm))

lr_tx_means <- mean_by_condition(
  lr_sig_tx_tpm[common_tx, , drop = FALSE],
  lr_meta,
  "Long-read"
)

sr_tx_means <- mean_by_condition(
  sr_sig_tx_tpm[common_tx, , drop = FALSE],
  sr_meta,
  "Short-read"
)

tx_means_wide <- bind_rows(lr_tx_means, sr_tx_means) %>%
  filter(condition %in% c(control_label, case_label)) %>%
  dplyr::select(feature_id, condition, platform, mean_log2_tpm, mean_tpm) %>%
  pivot_wider(
    names_from = c(platform, condition),
    values_from = c(mean_log2_tpm, mean_tpm),
    names_sep = "__"
  )

top_similar_tx_ids <- tx_means_wide %>%
  rowwise() %>%
  mutate(
    lr_ct = get(paste0("mean_log2_tpm__Long-read__", control_label)),
    lr_ad = get(paste0("mean_log2_tpm__Long-read__", case_label)),
    sr_ct = get(paste0("mean_log2_tpm__Short-read__", control_label)),
    sr_ad = get(paste0("mean_log2_tpm__Short-read__", case_label)),

    lr_mean_tpm = mean(c(
      get(paste0("mean_tpm__Long-read__", control_label)),
      get(paste0("mean_tpm__Long-read__", case_label))
    ), na.rm = TRUE),

    sr_mean_tpm = mean(c(
      get(paste0("mean_tpm__Short-read__", control_label)),
      get(paste0("mean_tpm__Short-read__", case_label))
    ), na.rm = TRUE),

    mean_abs_delta = mean(abs(c(lr_ct - sr_ct, lr_ad - sr_ad)), na.rm = TRUE),
    ad_ct_direction_same = sign(lr_ad - lr_ct) == sign(sr_ad - sr_ct)
  ) %>%
  ungroup() %>%
  filter(lr_mean_tpm >= min_mean_tpm | sr_mean_tpm >= min_mean_tpm) %>%
  arrange(mean_abs_delta, desc(ad_ct_direction_same)) %>%
  slice_head(n = top_n) %>%
  pull(feature_id)

top_similar_tx_ids

tx_gene_lookup <- tx_gene_map %>%
  distinct(transcript_id, gene_id) %>%
  left_join(gene_annot, by = "gene_id") %>%
  mutate(
    gene_symbol = ifelse(is.na(gene_symbol), gene_id, gene_symbol),
    gene_label = ifelse(is.na(gene_label), gene_id, gene_label)
  )

tx_expr_plot_df <- bind_rows(lr_tx_means, sr_tx_means) %>%
  filter(
    feature_id %in% top_similar_tx_ids,
    condition %in% c(control_label, case_label)
  ) %>%
  dplyr::rename(transcript_id = feature_id) %>%
  left_join(tx_gene_lookup, by = "transcript_id") %>%
  mutate(
    transcript_label = paste0(gene_symbol, " | ", transcript_id),
    condition = factor(condition, levels = c(control_label, case_label)),
    platform_condition = paste(platform, condition, sep = " - ")
  )

head(tx_expr_plot_df)

tx_label_levels <- tibble(
  transcript_id = top_similar_tx_ids
) %>%
  left_join(tx_gene_lookup, by = "transcript_id") %>%
  mutate(
    gene_symbol = ifelse(is.na(gene_symbol), gene_id, gene_symbol),
    transcript_label = paste0(gene_symbol, " | ", transcript_id)
  )

tx_expr_plot_df <- tx_expr_plot_df %>%
  mutate(
    transcript_label = factor(
      transcript_label,
      levels = rev(tx_label_levels$transcript_label)
    )
  )

p_top_similar_transcripts <- ggplot(
  tx_expr_plot_df,
  aes(x = transcript_label, y = mean_log2_tpm, fill = platform_condition)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Top similar significant-gene transcripts: AD vs CT expression",
    subtitle = "Transcripts ranked by long-read vs short-read expression-profile similarity",
    x = "Transcript",
    y = "Mean log2(TPM + 1)",
    fill = "Platform / condition"
  )

p_top_similar_transcripts
ggsave(
  "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/results/compare_longread_shortread/top_similar_transcripts_AD_CT_expression_barplot.png",
  p_top_similar_transcripts,
  width = 12,
  height = 8,
  dpi = 300
)



#################################
# Pick transcripts from sig genes
#################################

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

strip_version <- function(x) {
  sub("\\.\\d+$", "", as.character(x))
}

case_label <- "AD"
control_label <- "CT"

sig <- unique(strip_version(sig_genes[[1]]))

tx_gene_lookup <- tx_gene_map %>%
  mutate(
    transcript_id = strip_version(transcript_id),
    gene_id = strip_version(gene_id)
  ) %>%
  distinct(transcript_id, gene_id)

target_tx <- tx_gene_lookup %>%
  filter(gene_id %in% sig) %>%
  pull(transcript_id) %>%
  unique()

# Make sure rownames are version-stripped
lr_tpm <- as.matrix(lr_sig_tx_tpm)
sr_tpm <- as.matrix(sr_sig_tx_tpm)

rownames(lr_tpm) <- strip_version(rownames(lr_tpm))
rownames(sr_tpm) <- strip_version(rownames(sr_tpm))

common_target_tx <- Reduce(intersect, list(
  target_tx,
  rownames(lr_tpm),
  rownames(sr_tpm)
))

length(common_target_tx)

##############################################
# Function to run Welch t-tests per transcript
##############################################
run_tx_ttests <- function(tx_tpm, meta, platform,
                          transcripts = NULL,
                          case_label = "AD",
                          control_label = "CT") {
  tx_tpm <- as.matrix(tx_tpm)

  if (!is.null(transcripts)) {
    transcripts <- intersect(transcripts, rownames(tx_tpm))
    tx_tpm <- tx_tpm[transcripts, , drop = FALSE]
  }

  meta <- meta %>%
    mutate(
      Run = as.character(Run),
      condition = as.character(condition)
    ) %>%
    filter(
      Run %in% colnames(tx_tpm),
      condition %in% c(control_label, case_label)
    ) %>%
    arrange(match(Run, colnames(tx_tpm)))

  tx_tpm <- tx_tpm[, meta$Run, drop = FALSE]

  log_tpm <- log2(tx_tpm + 1)

  res <- lapply(rownames(log_tpm), function(tx) {
    ad_samples <- meta$Run[meta$condition == case_label]
    ct_samples <- meta$Run[meta$condition == control_label]

    ad_vals <- as.numeric(log_tpm[tx, ad_samples])
    ct_vals <- as.numeric(log_tpm[tx, ct_samples])

    ad_raw <- as.numeric(tx_tpm[tx, ad_samples])
    ct_raw <- as.numeric(tx_tpm[tx, ct_samples])

    # Need at least 2 samples per group
    if (length(ad_vals) < 2 || length(ct_vals) < 2) {
      pval <- NA_real_
      tstat <- NA_real_
    } else if (sd(ad_vals, na.rm = TRUE) == 0 && sd(ct_vals, na.rm = TRUE) == 0) {
      pval <- NA_real_
      tstat <- NA_real_
    } else {
      tt <- t.test(ad_vals, ct_vals, var.equal = FALSE)
      pval <- tt$p.value
      tstat <- unname(tt$statistic)
    }

    tibble(
      transcript_id = tx,
      platform = platform,
      mean_CT_tpm = mean(ct_raw, na.rm = TRUE),
      mean_AD_tpm = mean(ad_raw, na.rm = TRUE),
      mean_CT_log2_tpm = mean(ct_vals, na.rm = TRUE),
      mean_AD_log2_tpm = mean(ad_vals, na.rm = TRUE),
      log2FC_AD_vs_CT = mean(ad_vals, na.rm = TRUE) - mean(ct_vals, na.rm = TRUE),
      t_statistic = tstat,
      p_value = pval,
      n_CT = length(ct_vals),
      n_AD = length(ad_vals)
    )
  })

  bind_rows(res) %>%
    mutate(
      padj = p.adjust(p_value, method = "BH")
    )
}

######################################
# Run t-tests separately for LR and SR
######################################

lr_ttest <- run_tx_ttests(
  lr_tpm,
  lr_meta,
  platform = "Long-read",
  transcripts = common_target_tx,
  case_label = case_label,
  control_label = control_label
)

sr_ttest <- run_tx_ttests(
  sr_tpm,
  sr_meta,
  platform = "Short-read",
  transcripts = common_target_tx,
  case_label = case_label,
  control_label = control_label
)

###############################
# Join LR and SR t-test results
###############################
alpha <- 0.05
lfc_cutoff <- 0       # change to 0.25 or 0.5 if you want a minimum effect size

lr_wide <- lr_ttest %>%
  dplyr::select(
    transcript_id,
    lr_log2FC = log2FC_AD_vs_CT,
    lr_p = p_value,
    lr_padj = padj,
    lr_mean_CT_tpm = mean_CT_tpm,
    lr_mean_AD_tpm = mean_AD_tpm,
    lr_mean_CT_log2_tpm = mean_CT_log2_tpm,
    lr_mean_AD_log2_tpm = mean_AD_log2_tpm
  )

sr_wide <- sr_ttest %>%
  dplyr::select(
    transcript_id,
    sr_log2FC = log2FC_AD_vs_CT,
    sr_p = p_value,
    sr_padj = padj,
    sr_mean_CT_tpm = mean_CT_tpm,
    sr_mean_AD_tpm = mean_AD_tpm,
    sr_mean_CT_log2_tpm = mean_CT_log2_tpm,
    sr_mean_AD_log2_tpm = mean_AD_log2_tpm
  )

tx_ttest_comparison <- inner_join(
  lr_wide,
  sr_wide,
  by = "transcript_id"
) %>%
  left_join(tx_gene_lookup, by = "transcript_id") %>%
  mutate(
    lr_significant = lr_padj < alpha & abs(lr_log2FC) >= lfc_cutoff,
    sr_significant = sr_padj < alpha & abs(sr_log2FC) >= lfc_cutoff,
    significant_in_both = lr_significant & sr_significant,

    same_direction = sign(lr_log2FC) == sign(sr_log2FC),

    direction = case_when(
      lr_log2FC > 0 & sr_log2FC > 0 ~ "Up in AD in both",
      lr_log2FC < 0 & sr_log2FC < 0 ~ "Down in AD in both",
      lr_log2FC > 0 & sr_log2FC < 0 ~ "LR up, SR down",
      lr_log2FC < 0 & sr_log2FC > 0 ~ "LR down, SR up",
      TRUE ~ "No clear direction"
    ),

    abs_log2FC_delta = abs(lr_log2FC - sr_log2FC),

    effect_size_ratio = pmin(abs(lr_log2FC), abs(sr_log2FC)) /
      pmax(abs(lr_log2FC), abs(sr_log2FC)),

    effect_size_ratio = ifelse(is.nan(effect_size_ratio), NA_real_, effect_size_ratio),

    similar_difference = significant_in_both &
      same_direction &
      effect_size_ratio >= 0.5
  ) %>%
  arrange(
    desc(similar_difference),
    desc(significant_in_both),
    abs_log2FC_delta
  )

head(tx_ttest_comparison)

sig_in_both_similar <- tx_ttest_comparison %>%
  filter(similar_difference) %>%
  arrange(abs_log2FC_delta, desc(effect_size_ratio))

sig_in_both_similar

tx_ttest_comparison %>%
  select(
    gene_id,
    transcript_id,
    lr_log2FC,
    sr_log2FC,
    lr_padj,
    sr_padj,
    significant_in_both,
    same_direction,
    effect_size_ratio,
    abs_log2FC_delta,
    similar_difference,
    direction
  )

tx_ttest_summary <- tx_ttest_comparison %>%
  summarise(
    n_tested = n(),
    n_lr_significant = sum(lr_significant, na.rm = TRUE),
    n_sr_significant = sum(sr_significant, na.rm = TRUE),
    n_significant_in_both = sum(significant_in_both, na.rm = TRUE),
    n_sig_both_same_direction = sum(significant_in_both & same_direction, na.rm = TRUE),
    n_sig_both_similar = sum(similar_difference, na.rm = TRUE),
    spearman_log2FC = cor(lr_log2FC, sr_log2FC, method = "spearman", use = "complete.obs"),
    pearson_log2FC = cor(lr_log2FC, sr_log2FC, method = "pearson", use = "complete.obs")
  )

tx_ttest_summary

#############################
# Plot LR vs SR AD-CT effects
#############################

p_tx_ttest_scatter <- ggplot(
  tx_ttest_comparison,
  aes(x = sr_log2FC, y = lr_log2FC)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(aes(shape = significant_in_both), size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  labs(
    title = "Transcript AD vs CT effects in short-read and long-read",
    subtitle = "Welch t-tests on log2(TPM + 1)",
    x = "Short-read log2FC: AD - CT",
    y = "Long-read log2FC: AD - CT",
    shape = "Significant in both"
  )

p_tx_ttest_scatter

write.csv(
  tx_ttest_comparison,
  "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/results/compare_longread_shortread/transcript_AD_vs_CT_ttest_LR_SR_comparison.csv",
  row.names = FALSE
)

write.csv(
  sig_in_both_similar,
  "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/results/compare_longread_shortread/transcripts_significant_similar_AD_vs_CT_ttest_LR_SR.csv",
  row.names = FALSE
)

write.csv(
  tx_ttest_summary,
  "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/results/compare_longread_shortread/transcript_AD_vs_CT_ttest_LR_SR_summary.csv",
  row.names = FALSE
)

###################################################################################
# Large difference
###################################################################################
strip_version <- function(x) {
  sub("\\.\\d+$", "", as.character(x))
}

case_label <- "AD"
control_label <- "CT"

# Set these based on how strict you want to be
min_abs_log2FC <- 0.5       # 0.5 = moderate; 1 = stronger, about 2-fold on log2 scale
min_mean_tpm <- 0.5         # avoid transcripts that are basically absent
min_effect_ratio <- 0.5     # 0.5 means smaller effect must be at least half of larger effect

calc_tx_effects <- function(tx_tpm, meta, platform,
                            case_label = "AD",
                            control_label = "CT") {
  tx_tpm <- as.matrix(tx_tpm)

  rownames(tx_tpm) <- strip_version(rownames(tx_tpm))

  meta <- meta %>%
    mutate(
      Run = as.character(Run),
      condition = as.character(condition)
    ) %>%
    filter(
      Run %in% colnames(tx_tpm),
      condition %in% c(control_label, case_label)
    ) %>%
    arrange(match(Run, colnames(tx_tpm)))

  tx_tpm <- tx_tpm[, meta$Run, drop = FALSE]

  ct_samples <- meta$Run[meta$condition == control_label]
  ad_samples <- meta$Run[meta$condition == case_label]

  log_tpm <- log2(tx_tpm + 1)

  tibble(
    transcript_id = rownames(tx_tpm),
    platform = platform,

    mean_CT_tpm = rowMeans(tx_tpm[, ct_samples, drop = FALSE], na.rm = TRUE),
    mean_AD_tpm = rowMeans(tx_tpm[, ad_samples, drop = FALSE], na.rm = TRUE),

    mean_CT_log2_tpm = rowMeans(log_tpm[, ct_samples, drop = FALSE], na.rm = TRUE),
    mean_AD_log2_tpm = rowMeans(log_tpm[, ad_samples, drop = FALSE], na.rm = TRUE),

    log2FC_AD_vs_CT =
      rowMeans(log_tpm[, ad_samples, drop = FALSE], na.rm = TRUE) -
      rowMeans(log_tpm[, ct_samples, drop = FALSE], na.rm = TRUE),

    mean_tpm = rowMeans(tx_tpm, na.rm = TRUE),
    max_tpm = apply(tx_tpm, 1, max, na.rm = TRUE)
  )
}

lr_tpm <- as.matrix(lr_sig_tx_tpm)
sr_tpm <- as.matrix(sr_sig_tx_tpm)

rownames(lr_tpm) <- strip_version(rownames(lr_tpm))
rownames(sr_tpm) <- strip_version(rownames(sr_tpm))

common_tx <- intersect(rownames(lr_tpm), rownames(sr_tpm))

lr_effects <- calc_tx_effects(
  lr_tpm[common_tx, , drop = FALSE],
  lr_meta,
  platform = "Long-read",
  case_label = case_label,
  control_label = control_label
)

sr_effects <- calc_tx_effects(
  sr_tpm[common_tx, , drop = FALSE],
  sr_meta,
  platform = "Short-read",
  case_label = case_label,
  control_label = control_label
)

tx_gene_lookup <- tx_gene_map %>%
  mutate(
    transcript_id = strip_version(transcript_id),
    gene_id = strip_version(gene_id)
  ) %>%
  distinct(transcript_id, gene_id)

lr_wide <- lr_effects %>%
  dplyr::select(
    transcript_id,
    lr_log2FC = log2FC_AD_vs_CT,
    lr_mean_CT_tpm = mean_CT_tpm,
    lr_mean_AD_tpm = mean_AD_tpm,
    lr_mean_CT_log2_tpm = mean_CT_log2_tpm,
    lr_mean_AD_log2_tpm = mean_AD_log2_tpm,
    lr_mean_tpm = mean_tpm,
    lr_max_tpm = max_tpm
  )

sr_wide <- sr_effects %>%
  dplyr::select(
    transcript_id,
    sr_log2FC = log2FC_AD_vs_CT,
    sr_mean_CT_tpm = mean_CT_tpm,
    sr_mean_AD_tpm = mean_AD_tpm,
    sr_mean_CT_log2_tpm = mean_CT_log2_tpm,
    sr_mean_AD_log2_tpm = mean_AD_log2_tpm,
    sr_mean_tpm = mean_tpm,
    sr_max_tpm = max_tpm
  )

tx_similarity_effects <- inner_join(
  lr_wide,
  sr_wide,
  by = "transcript_id"
) %>%
  left_join(tx_gene_lookup, by = "transcript_id") %>%
  mutate(
    same_direction = sign(lr_log2FC) == sign(sr_log2FC),

    direction = case_when(
      lr_log2FC > 0 & sr_log2FC > 0 ~ "Higher in AD in both",
      lr_log2FC < 0 & sr_log2FC < 0 ~ "Lower in AD in both",
      lr_log2FC > 0 & sr_log2FC < 0 ~ "LR higher in AD, SR lower in AD",
      lr_log2FC < 0 & sr_log2FC > 0 ~ "LR lower in AD, SR higher in AD",
      TRUE ~ "No clear direction"
    ),

    lr_abs_log2FC = abs(lr_log2FC),
    sr_abs_log2FC = abs(sr_log2FC),

    min_abs_log2FC_both = pmin(lr_abs_log2FC, sr_abs_log2FC),
    max_abs_log2FC_both = pmax(lr_abs_log2FC, sr_abs_log2FC),

    abs_log2FC_delta = abs(lr_log2FC - sr_log2FC),

    effect_size_ratio = min_abs_log2FC_both / max_abs_log2FC_both,
    effect_size_ratio = ifelse(is.nan(effect_size_ratio), NA_real_, effect_size_ratio),

    mean_abs_log2FC = mean(c(lr_abs_log2FC, sr_abs_log2FC), na.rm = TRUE),

    expressed_enough =
      lr_mean_tpm >= min_mean_tpm | sr_mean_tpm >= min_mean_tpm,

    large_difference_in_both =
      lr_abs_log2FC >= min_abs_log2FC &
      sr_abs_log2FC >= min_abs_log2FC,

    similar_large_difference =
      expressed_enough &
      large_difference_in_both &
      same_direction &
      effect_size_ratio >= min_effect_ratio,

    similarity_rank_score =
      ifelse(
        same_direction,
        min_abs_log2FC_both * effect_size_ratio,
        -min_abs_log2FC_both
      )
  ) %>%
  arrange(
    desc(similar_large_difference),
    desc(similarity_rank_score),
    abs_log2FC_delta
  )

tx_similarity_effects <- merge(tx_similarity_effects, gene_annot, by="gene_id")

similar_large_tx <- tx_similarity_effects %>%
  filter(similar_large_difference) %>%
  arrange(desc(similarity_rank_score), abs_log2FC_delta)

similar_large_tx %>%
  dplyr::select(
    gene_id,
    gene_label,
    transcript_id,
    direction,
    lr_log2FC,
    sr_log2FC,
    lr_mean_CT_tpm,
    lr_mean_AD_tpm,
    sr_mean_CT_tpm,
    sr_mean_AD_tpm,
    effect_size_ratio,
    abs_log2FC_delta,
    similarity_rank_score
  )

tx_similarity_summary <- tx_similarity_effects %>%
  summarise(
    n_transcripts_compared = n(),

    n_same_direction = sum(same_direction, na.rm = TRUE),

    n_large_difference_lr = sum(lr_abs_log2FC >= min_abs_log2FC, na.rm = TRUE),
    n_large_difference_sr = sum(sr_abs_log2FC >= min_abs_log2FC, na.rm = TRUE),

    n_large_difference_both = sum(large_difference_in_both, na.rm = TRUE),

    n_large_same_direction = sum(
      large_difference_in_both & same_direction,
      na.rm = TRUE
    ),

    n_large_similar_difference = sum(similar_large_difference, na.rm = TRUE),

    percent_large_both_same_direction =
      100 * mean(same_direction[large_difference_in_both], na.rm = TRUE),

    spearman_log2FC_all = cor(
      lr_log2FC,
      sr_log2FC,
      method = "spearman",
      use = "complete.obs"
    ),

    pearson_log2FC_all = cor(
      lr_log2FC,
      sr_log2FC,
      method = "pearson",
      use = "complete.obs"
    )
  )

tx_similarity_summary

p_tx_similarity_scatter <- ggplot(
  tx_similarity_effects,
  aes(x = sr_log2FC, y = lr_log2FC)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(
    aes(shape = similar_large_difference),
    alpha = 0.8,
    size = 2
  ) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  labs(
    title = "Transcript AD vs CT effects in short-read and long-read",
    subtitle = paste0(
      "Highlighted transcripts have |log2FC| >= ",
      min_abs_log2FC,
      " in both platforms and similar direction"
    ),
    x = "Short-read log2FC: AD - CT",
    y = "Long-read log2FC: AD - CT",
    shape = "Large and similar"
  )

p_tx_similarity_scatter

top_large_similar_tx <- similar_large_tx %>%
  slice_head(n = 20) %>%
  pull(transcript_id)

plot_bar_df <- tx_similarity_effects %>%
  filter(transcript_id %in% top_large_similar_tx) %>%
  dplyr::select(
    transcript_id,
    gene_id,
    gene_label,
    lr_mean_CT_log2_tpm,
    lr_mean_AD_log2_tpm,
    sr_mean_CT_log2_tpm,
    sr_mean_AD_log2_tpm
  ) %>%
  pivot_longer(
    cols = c(
      lr_mean_CT_log2_tpm,
      lr_mean_AD_log2_tpm,
      sr_mean_CT_log2_tpm,
      sr_mean_AD_log2_tpm
    ),
    names_to = "measurement",
    values_to = "mean_log2_tpm"
  ) %>%
  mutate(
    platform = case_when(
      grepl("^lr_", measurement) ~ "Long-read",
      grepl("^sr_", measurement) ~ "Short-read"
    ),
    condition = case_when(
      grepl("_CT_", measurement) ~ "CT",
      grepl("_AD_", measurement) ~ "AD"
    ),
    platform_condition = paste(platform, condition, sep = " - "),
    transcript_label = paste0(transcript_id, "\n", gene_label),
    transcript_label = factor(
      transcript_label,
      levels = rev(unique(transcript_label[match(top_large_similar_tx, transcript_id)]))
    )
  )

p_top_large_similar_tx <- ggplot(
  plot_bar_df,
  aes(
    x = transcript_label,
    y = mean_log2_tpm,
    fill = condition,
    pattern = platform,
    group = interaction(platform, condition)
  )
) +
  geom_col_pattern(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black",
    pattern_angle = 45,
    pattern_density = 0.35,
    pattern_spacing = 0.03,
    pattern_fill = "black",
    pattern_colour = "black"
  ) +
  scale_fill_manual(
    values = c(
      "AD" = "salmon",
      "CT" = "turquoise"
    )
  ) +
  scale_pattern_manual(
    values = c(
      "Long-read" = "none",
      "Short-read" = "stripe"
    )
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        pattern = "none"
      )
    ),
    pattern = guide_legend(
      override.aes = list(
        fill = "white"
      )
    )
  ) +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Top transcripts with large, similar AD-vs-CT differences",
    subtitle = "Ranked by concordant LR/SR log2FC magnitude",
    x = "Transcript / gene",
    y = "Mean log2(TPM + 1)",
    fill = "Condition",
    pattern = "Platform"
  )

p_top_large_similar_tx

plot_fc_df <- tx_similarity_effects %>%
  filter(transcript_id %in% top_large_similar_tx) %>%
  dplyr::select(
    transcript_id,
    gene_id,
    gene_label,
    lr_log2FC,
    sr_log2FC
  ) %>%
  pivot_longer(
    cols = c(lr_log2FC, sr_log2FC),
    names_to = "platform",
    values_to = "log2FC"
  ) %>%
  mutate(
    platform = case_when(
      platform == "lr_log2FC" ~ "Long-read",
      platform == "sr_log2FC" ~ "Short-read"
    ),
    transcript_label = paste0(transcript_id, "\n", gene_label),
    transcript_label = factor(
      transcript_label,
      levels = rev(unique(transcript_label[match(top_large_similar_tx, transcript_id)]))
    )
  )

p_top_large_similar_tx_fc <- ggplot(
  plot_fc_df,
  aes(x = transcript_label, y = log2FC, fill = platform)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  coord_flip() +
  theme_bw() +
  labs(
    title = "AD-vs-CT effect sizes for top similar transcripts",
    subtitle = "Positive = higher in AD; negative = lower in AD",
    x = "Transcript / gene",
    y = "log2FC: AD - CT",
    fill = "Platform"
  )

p_top_large_similar_tx_fc

write.csv(
  tx_similarity_effects,
  "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/results/compare_longread_shortread/transcript_LR_SR_AD_CT_similarity_effects_all.csv",
  row.names = FALSE
)

write.csv(
  similar_large_tx,
  "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/results/compare_longread_shortread/transcripts_large_similar_AD_CT_effects_LR_SR.csv",
  row.names = FALSE
)

write.csv(
  tx_similarity_summary,
  "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/results/compare_longread_shortread/transcript_LR_SR_AD_CT_similarity_summary.csv",
  row.names = FALSE
)

ggsave(
  "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/results/compare_longread_shortread/transcript_LR_SR_AD_CT_effect_similarity_scatter.png",
  p_tx_similarity_scatter,
  width = 8,
  height = 7,
  dpi = 300
)

ggsave(
  "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/results/compare_longread_shortread/top_large_similar_transcripts_AD_CT_expression_barplot.png",
  p_top_large_similar_tx,
  width = 12,
  height = 9,
  dpi = 300
)

ggsave(
  "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/results/compare_longread_shortread/top_large_similar_transcripts_AD_CT_log2FC_barplot.png",
  p_top_large_similar_tx_fc,
  width = 12,
  height = 9,
  dpi = 300
)