# long distance Motif Enrichment (2kb upstream) with PWMEnrich for AD promoter switching analysis
library(data.table)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(PWMEnrich)
library(PWMEnrich.Hsapiens.background)
library(ggplot2)

output_dir_tables <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/PRJNA1206164/tables/"
output_dir_figures <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/PRJNA1206164/figures/"
# promoter switching table 
prom_dt <- read.table(paste0(output_dir_tables, "promoter_usage_betareg_results.tsv"), header = TRUE, sep = "\t")
# split promoter_id into seqnames/start/end/strand (e.g., 19:18587685-18588785(+))
# handle promoter_id entries with strand "*" as well as "+" / "-"
    sig <- prom_dt %>%
        tidyr::extract(
            col    = promoter_id,
            into   = c("seqnames", "start", "end", "strand"),
            regex  = "^([^:]+):([0-9]+)-([0-9]+)\\(([+\\-*])\\)$",
            remove = FALSE
        ) %>%
        dplyr::mutate(
            seqnames = as.character(seqnames),
            start    = as.integer(start),
            end      = as.integer(end),
            strand   = as.character(strand)
        )

prom_gr <- GRanges(
  seqnames = sig$seqnames,
  ranges   = IRanges(start = sig$start, end = sig$end),
  strand   = sig$strand
)
mcols(prom_gr) <- sig

bad_names <- c("seqnames","ranges","strand","seqlevels","seqlengths",
               "isCircular","start","end","width","element")

mc <- mcols(prom_gr)
mcols(prom_gr) <- mc[, setdiff(colnames(mc), bad_names), drop = FALSE]
sig_prom  <- prom_gr[mcols(prom_gr)$padj <= 0.05]

sig_prom <- as.data.frame(sig_prom) %>%
    group_by(gene_name) %>%
    filter(n_distinct(promoter_id) >= 2) %>%
    ungroup() %>%
    group_by(gene_name) %>%
    filter(n_distinct(promoter_id) >= 2) %>%
    ungroup()

sig_prom <- sig_prom %>%
  group_by(gene_name) %>%
  filter(
    any(logitFC > 0, na.rm = TRUE) &
    any(logitFC < 0, na.rm = TRUE)
  ) %>%
  ungroup()
sig_prom <- makeGRangesFromDataFrame(sig_prom, keep.extra.columns = TRUE)
unique(mcols(sig_prom)$gene_name)
prom_up   <- sig_prom[mcols(sig_prom)$logitFC > 0]
prom_down <- sig_prom[mcols(sig_prom)$logitFC < 0]

# prom_up   <- sig_prom[mcols(sig_prom)$logitFC > 1]
# prom_down <- sig_prom[mcols(sig_prom)$logitFC < -1]

length(prom_up)
length(prom_down)

# function to get TSS per GRanges
get_promoter_TSS <- function(gr) {
  is_plus  <- strand(gr) == "+"
  is_minus <- strand(gr) == "-"

  tss <- ifelse(is_plus, start(gr), end(gr))
  GRanges(seqnames = seqnames(gr),
          ranges   = IRanges(start = tss, width = 1),
          strand   = strand(gr),
          mcols(gr))
}



prom_up_tss   <- get_promoter_TSS(prom_up)
prom_down_tss <- get_promoter_TSS(prom_down)

long_upstream  <- 2000L
long_downstream <- 2000L

# Significant promoters: TSS-centered ±2kb
prom_up_win   <- promoters(prom_up_tss,   upstream = long_upstream, downstream = long_downstream)
prom_down_win <- promoters(prom_down_tss, upstream = long_upstream, downstream = long_downstream)

# background from nonsig promoters
nonsig_prom <- prom_gr[mcols(prom_gr)$padj > 0.2 | is.na(mcols(prom_gr)$padj)]
# Background: same window ±2kb
nonsig_tss  <- get_promoter_TSS(nonsig_prom)
bg_win      <- promoters(nonsig_tss, upstream = long_upstream, downstream = long_downstream)

seqlevelsStyle(prom_up_win) <- "UCSC"
seqlevelsStyle(prom_down_win) <- "UCSC"
seqlevelsStyle(bg_win) <- "UCSC"
# sequences
hg38      <- BSgenome.Hsapiens.UCSC.hg38
seqs_up   <- getSeq(hg38, prom_up_win)
seqs_down <- getSeq(hg38, prom_down_win)

strand(bg_win) <- "*"
seqs_bg <- getSeq(hg38, bg_win)

library(PWMEnrich)
library(PWMEnrich.Hsapiens.background)

data(PWMLogn.hg19.MotifDb.Hsap)


res_up     <- motifEnrichment(seqs_up, PWMLogn.hg19.MotifDb.Hsap)
res_up_sum <- groupReport(res_up)

res_down     <- motifEnrichment(seqs_down, PWMLogn.hg19.MotifDb.Hsap)
res_down_sum <- groupReport(res_down)

write.csv(as.data.frame(res_up_sum),   paste0(output_dir_tables, "motif_enrichment_promoters_up_AD_2kbps.csv"),   row.names = FALSE)
write.csv(as.data.frame(res_down_sum), paste0(output_dir_tables, "motif_enrichment_promoters_down_AD_2kbps.csv"), row.names = FALSE)

df_up   <- as.data.frame(res_up_sum)
df_down <- as.data.frame(res_down_sum)

# peek at columns
colnames(df_up)

library(dplyr)
library(ggplot2)

topN <- 20

df_up <- as.data.frame(res_up_sum) %>%
  arrange(p.value) %>%
  slice_head(n = topN) %>%
  mutate(direction = "Down in AD")

df_down <- as.data.frame(res_down_sum) %>%
  arrange(p.value) %>%
  slice_head(n = topN) %>%
  mutate(direction = "Up in AD")

df_motifs <- bind_rows(df_up, df_down)

# collapse duplicated motifs (keep only best)
df_motifs_unique <- df_motifs %>%
  group_by(target, direction) %>%
  slice_min(p.value, n = 1) %>% 
  ungroup()

target_levels <- df_motifs_unique %>%
  arrange(p.value) %>%
  pull(target) %>%
  as.character() %>%
  unique()

df_motifs_unique$target <- factor(
  df_motifs_unique$target,
  levels = target_levels
)

ggplot(df_motifs_unique,
       aes(x = direction, y = target)) +
  geom_point(aes(size = -log10(p.value), color = raw.score)) +
  scale_color_viridis_c() +
  scale_size_continuous(name = "-log10(p.value)") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12)) +
  labs(title = "Motif enrichment in AD promoter switching",
       x = "", y = "Motif / TF") +
  theme(
    plot.margin = margin(t = 10, r = 100, b = 10, l = 10)
  )
  
ggsave(paste0(output_dir_figures, "motif_enrichment_promoter_switching_AD_2kbps.png"), width = 6, height = 8, dpi = 300)



library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(MotifDb)
library(TFBSTools)
library(motifmatchr)
library(circlize)
library(dplyr)

# get motifs from MotifDb
motifs_all <- MotifDb
motifs_hs_all <- motifs_all[grepl("^Hsapiens-", names(motifs_all))]

ad_tfs <- readLines("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn25671149/jaquino_analysis/tmp/promoter_usage_61_genes/AD_TFs.txt")
# keep motifs whose name contains any of your TFs
ad_pat <- paste0("\\b(", paste(ad_tfs, collapse="|"), ")\\b")
motifs_ad <- motifs_all[grepl(ad_pat, names(motifs_all), ignore.case = TRUE)]

length(motifs_ad)
head(names(motifs_ad), 3)

library(MotifDb)
library(TFBSTools)
library(S4Vectors)

# keep human motifs
motifs_hs <- motifs_ad[grepl("^Hsapiens-", names(motifs_ad))]

# helper: turn a single motif into a PWMatrix
# Robust path: coerce to integer PFMatrix first, then to PWM.
motif_to_pwm <- function(m, nm, pseudocount = 0.8, count_scale = 1000L) {
  mat <- as.matrix(m)
  storage.mode(mat) <- "double"

  rn <- toupper(rownames(mat))
  if (!all(c("A", "C", "G", "T") %in% rn)) {
    stop("Motif does not have A/C/G/T rownames: ", nm)
  }

  rownames(mat) <- rn
  mat <- mat[c("A", "C", "G", "T"), , drop = FALSE]
  mat[!is.finite(mat)] <- 0
  mat[mat < 0] <- 0

  cs <- colSums(mat)
  valid_cols <- cs > 0
  if (!all(valid_cols)) {
    mat <- mat[, valid_cols, drop = FALSE]
    cs <- cs[valid_cols]
  }
  if (ncol(mat) == 0) {
    stop("Motif has no valid columns after cleaning: ", nm)
  }

  # If matrix looks like probabilities, scale to pseudo-counts.
  if (max(mat, na.rm = TRUE) <= 1.0) {
    mat <- mat * as.numeric(count_scale)
  }

  # Build strictly integer PFM to satisfy TFBSTools validation.
  pfm_counts <- round(mat)
  pfm_counts[!is.finite(pfm_counts)] <- 0
  pfm_counts[pfm_counts < 0] <- 0
  storage.mode(pfm_counts) <- "integer"

  if (sum(pfm_counts) == 0) {
    stop("Motif has zero total counts after integer conversion: ", nm)
  }

  pfm <- TFBSTools::PFMatrix(
    ID = nm,
    name = nm,
    profileMatrix = pfm_counts,
    bg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25)
  )

  TFBSTools::toPWM(
    pfm,
    type = "log2probratio",
    pseudocounts = pseudocount,
    bg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25)
  )
}

# build PWMatrixList
pw_list <- mapply(
  FUN = function(m, nm) {
    tryCatch(
      motif_to_pwm(m, nm),
      error = function(e) {
        message("Skipping motif '", nm, "' during PWMatrix conversion: ", conditionMessage(e))
        NULL
      }
    )
  },
  m  = as.list(motifs_hs),
  nm = names(motifs_hs),
  SIMPLIFY = FALSE
)

pw_list <- Filter(Negate(is.null), pw_list)

if (length(pw_list) == 0) {
  stop("No motifs could be converted to PWMatrix objects; check motif set/mapping.")
}

pfms_ad <- do.call(TFBSTools::PWMatrixList, pw_list)

pfms_ad

library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

hg38 <- BSgenome.Hsapiens.UCSC.hg38

hits_up   <- matchMotifs(pfms_ad, prom_up_win,   genome = hg38, out = "scores")
hits_down <- matchMotifs(pfms_ad, prom_down_win, genome = hg38, out = "scores")

score_mat_up   <- motifmatchr::motifScores(hits_up)
score_mat_down <- motifmatchr::motifScores(hits_down)


score_mat_up   <- motifmatchr::motifScores(hits_up)   # promoters x motifs
score_mat_down <- motifmatchr::motifScores(hits_down)

score_mat_up_dense <- as.matrix(score_mat_up)
score_mat_down_dense <- as.matrix(score_mat_down)

# ==========================================
# Motif/TF proximity to significant promoters
# ==========================================
library(Biostrings)

# Extract a TF-like symbol from MotifDb names (e.g.,
# Hsapiens-HOCOMOCOv13-BCL6.H13CORE.0.PSM.A -> BCL6)
extract_tf_symbol <- function(x) {
  x <- as.character(x)
  x <- toupper(x)

  # Keep the last '-' token, then trim trailing motif metadata
  x <- sub("^.*-", "", x)
  x <- sub("[\\._ ].*$", "", x)
  x <- sub("-MA[0-9.]+.*$", "", x)
  x <- sub("-[0-9]+$", "", x)

  x
}

# Define significant motifs from enrichment summaries
sig_motif_tbl <- bind_rows(
  as.data.frame(res_up_sum) %>% mutate(direction = "Down in AD"),
  as.data.frame(res_down_sum) %>% mutate(direction = "Up in AD")
) %>%
  filter(!is.na(p.value), p.value <= 0.05)

# Map significant PWMEnrich TF targets -> MotifDb motif IDs
motif_tf_map <- data.frame(
  motif_name = names(motifs_hs),
  TF = vapply(names(motifs_hs), extract_tf_symbol, character(1)),
  stringsAsFactors = FALSE
) %>%
  distinct(motif_name, TF)

# Global human motif map (not restricted to AD TF list), used for top-10 heatmap
motif_tf_map_all <- data.frame(
  motif_name = names(motifs_hs_all),
  TF = vapply(names(motifs_hs_all), extract_tf_symbol, character(1)),
  stringsAsFactors = FALSE
) %>%
  distinct(motif_name, TF)

sig_tf_targets <- unique(toupper(as.character(sig_motif_tbl$target)))
sig_motif_names <- motif_tf_map %>%
  filter(TF %in% sig_tf_targets) %>%
  pull(motif_name) %>%
  unique()

top10_hs_tf_targets <- sig_motif_tbl %>%
  mutate(TF = toupper(as.character(target))) %>%
  filter(TF %in% motif_tf_map_all$TF) %>%
  group_by(TF) %>%
  summarise(best_p = min(p.value, na.rm = TRUE), .groups = "drop") %>%
  arrange(best_p) %>%
  slice_head(n = 10) %>%
  pull(TF)

top10_hs_motif_names <- motif_tf_map_all %>%
  filter(TF %in% top10_hs_tf_targets) %>%
  pull(motif_name) %>%
  unique()

message(
  "Significant PWMEnrich targets (p<=0.05): ", length(sig_tf_targets),
  "; matched MotifDb motifs: ", length(sig_motif_names),
  "; matched TFs: ", dplyr::n_distinct(motif_tf_map$TF[motif_tf_map$motif_name %in% sig_motif_names])
)

message(
  "Top-10 significant TF targets with human motifs: ", length(top10_hs_tf_targets),
  "; mapped human motifs: ", length(top10_hs_motif_names)
)

if (length(sig_motif_names) == 0) {
  message("No significant motifs (p.value <= 0.05) were found in MotifDb motif set; skipping motif-distance plots.")
} else {
  # Reuse already-validated PWMs from pw_list to avoid reconversion failures
  pwm_lookup <- stats::setNames(pw_list, names(motifs_hs))

  # Build PWM lookup for top-10 significant TF targets with human motifs
  pw_list_top10 <- mapply(
    FUN = function(m, nm) {
      tryCatch(
        motif_to_pwm(m, nm),
        error = function(e) {
          message("Skipping motif '", nm, "' during top10 PWMatrix conversion: ", conditionMessage(e))
          NULL
        }
      )
    },
    m = as.list(motifs_hs_all[top10_hs_motif_names]),
    nm = top10_hs_motif_names,
    SIMPLIFY = FALSE
  )
  pw_list_top10 <- Filter(Negate(is.null), pw_list_top10)
  pwm_lookup_top10 <- stats::setNames(pw_list_top10, names(pw_list_top10))

  # Scan one promoter sequence for motif hits on both strands
  scan_one_sequence <- function(seq_obj, pwm, min_score = "85%") {
    if (is.null(pwm) || length(seq_obj) == 0) return(data.frame())

    seq_len <- nchar(as.character(seq_obj))
    pwm_mat <- if (methods::is(pwm, "PWMatrix")) TFBSTools::Matrix(pwm) else as.matrix(pwm)

    hits_fwd <- Biostrings::matchPWM(pwm_mat, seq_obj, min.score = min_score)
    hits_rev <- Biostrings::matchPWM(pwm_mat, Biostrings::reverseComplement(seq_obj), min.score = min_score)

    df_fwd <- if (length(hits_fwd) > 0) {
      data.frame(start = start(hits_fwd), end = end(hits_fwd), strand_hit = "+")
    } else {
      data.frame()
    }

    df_rev <- if (length(hits_rev) > 0) {
      data.frame(
        start = seq_len - end(hits_rev) + 1,
        end = seq_len - start(hits_rev) + 1,
        strand_hit = "-"
      )
    } else {
      data.frame()
    }

    bind_rows(df_fwd, df_rev)
  }

  tss_idx <- long_upstream + 1L

  build_motif_distance_df <- function(seqs, prom_win, direction_label, motif_names, min_score = "85%") {
    if (length(seqs) == 0 || length(motif_names) == 0) return(data.frame())

    out <- vector("list", length = length(motif_names))

    for (mi in seq_along(motif_names)) {
      motif_name <- motif_names[[mi]]
      pwm <- pwm_lookup[[motif_name]]
      if (is.null(pwm)) next

      one_motif <- vector("list", length = length(seqs))

      for (i in seq_along(seqs)) {
        hits_df <- scan_one_sequence(seqs[[i]], pwm, min_score = min_score)
        if (nrow(hits_df) == 0) next

        hit_mid <- floor((hits_df$start + hits_df$end) / 2)
        dist_to_tss <- hit_mid - tss_idx

        one_motif[[i]] <- data.frame(
          promoter_id = as.character(mcols(prom_win)$promoter_id[[i]]),
          gene_name = as.character(mcols(prom_win)$gene_name[[i]]),
          direction = direction_label,
          motif_name = motif_name,
          TF = extract_tf_symbol(motif_name),
          hit_start = hits_df$start,
          hit_end = hits_df$end,
          hit_strand = hits_df$strand_hit,
          dist_to_tss_bp = dist_to_tss,
          abs_dist_to_tss_bp = abs(dist_to_tss)
        )
      }

      out[[mi]] <- bind_rows(one_motif)
    }

    bind_rows(out)
  }

  build_motif_distance_df_with_lookup <- function(seqs, prom_win, direction_label, motif_names, pwm_lookup_local, min_score = "85%") {
    if (length(seqs) == 0 || length(motif_names) == 0) return(data.frame())

    out <- vector("list", length = length(motif_names))

    for (mi in seq_along(motif_names)) {
      motif_name <- motif_names[[mi]]
      pwm <- pwm_lookup_local[[motif_name]]
      if (is.null(pwm)) next

      one_motif <- vector("list", length = length(seqs))

      for (i in seq_along(seqs)) {
        hits_df <- scan_one_sequence(seqs[[i]], pwm, min_score = min_score)
        if (nrow(hits_df) == 0) next

        hit_mid <- floor((hits_df$start + hits_df$end) / 2)
        dist_to_tss <- hit_mid - tss_idx

        one_motif[[i]] <- data.frame(
          promoter_id = as.character(mcols(prom_win)$promoter_id[[i]]),
          gene_name = as.character(mcols(prom_win)$gene_name[[i]]),
          direction = direction_label,
          motif_name = motif_name,
          TF = extract_tf_symbol(motif_name),
          hit_start = hits_df$start,
          hit_end = hits_df$end,
          hit_strand = hits_df$strand_hit,
          dist_to_tss_bp = dist_to_tss,
          abs_dist_to_tss_bp = abs(dist_to_tss)
        )
      }

      out[[mi]] <- bind_rows(one_motif)
    }

    bind_rows(out)
  }

  plot_gene_tf_heatmap <- function(motif_nearest_df, title_txt, out_file) {
    if (nrow(motif_nearest_df) == 0) {
      message("No motif-nearest rows for ", out_file, "; skipping.")
      return(invisible(NULL))
    }

    gene_tf_dist <- motif_nearest_df %>%
      group_by(gene_name, TF, direction) %>%
      summarise(min_dist_bp = min(nearest_dist_bp), .groups = "drop")

    tf_order <- gene_tf_dist %>%
      group_by(TF) %>%
      summarise(med = median(min_dist_bp), .groups = "drop") %>%
      arrange(med) %>%
      pull(TF)

    gene_order_dist <- gene_tf_dist %>%
      group_by(gene_name) %>%
      summarise(med = median(min_dist_bp), .groups = "drop") %>%
      arrange(med) %>%
      pull(gene_name)

    p_gene_tf <- ggplot(
      gene_tf_dist %>%
        mutate(
          TF = factor(TF, levels = tf_order),
          gene_name = factor(gene_name, levels = gene_order_dist)
        ),
      aes(x = TF, y = gene_name, fill = min_dist_bp)
    ) +
      geom_tile(color = "grey90") +
      facet_wrap(~direction) +
      scale_fill_viridis_c(option = "C", trans = "sqrt") +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12)
      ) +
      labs(
        title = title_txt,
        x = "TF",
        y = "Gene",
        fill = "Min distance\n(bp)"
      )

    ggsave(
      out_file,
      plot = p_gene_tf,
      width = 10,
      height = 8,
      dpi = 300
    )
  }

  motif_dist_up <- build_motif_distance_df(
    seqs = seqs_up,
    prom_win = prom_up_win,
    direction_label = "Down in AD",
    motif_names = sig_motif_names,
    min_score = "85%"
  )

  motif_dist_down <- build_motif_distance_df(
    seqs = seqs_down,
    prom_win = prom_down_win,
    direction_label = "Up in AD",
    motif_names = sig_motif_names,
    min_score = "85%"
  )

  motif_dist_all <- bind_rows(motif_dist_up, motif_dist_down)

  if (nrow(motif_dist_all) == 0) {
    message("No motif hits found at min.score = 85%; skipping motif-distance plots.")
  } else {
    write.csv(
      motif_dist_all,
      paste0(output_dir_tables, "significant_motif_hits_distance_to_promoter_tss.csv"),
      row.names = FALSE
    )

    # Keep nearest hit per promoter x TF for cleaner summary plots
    motif_nearest <- motif_dist_all %>%
      group_by(gene_name, promoter_id, TF, motif_name, direction) %>%
      summarise(nearest_dist_bp = min(abs_dist_to_tss_bp), .groups = "drop")

    top_tfs_dist <- motif_nearest %>%
      count(TF, sort = TRUE) %>%
      slice_head(n = 20) %>%
      pull(TF)

    motif_nearest_top <- motif_nearest %>%
      filter(TF %in% top_tfs_dist)

    p_dist_tf <- ggplot(motif_nearest_top,
                        aes(x = reorder(TF, nearest_dist_bp, FUN = median),
                            y = nearest_dist_bp,
                            fill = direction)) +
      geom_boxplot(outlier.alpha = 0.2, alpha = 0.85) +
      coord_flip() +
      theme_bw(base_size = 12) +
      labs(
        title = "Significant motif/TF proximity to significant promoter TSS",
        x = "TF",
        y = "Nearest motif hit distance to promoter TSS (bp)",
        fill = "Direction"
      )

    ggsave(
      paste0(output_dir_figures, "significant_TF_distance_to_promoter_TSS_boxplot.png"),
      plot = p_dist_tf,
      width = 10,
      height = 7,
      dpi = 300
    )

    plot_gene_tf_heatmap(
      motif_nearest_df = motif_nearest,
      title_txt = "Closest motif/TF hit to significant promoters (AD TF list)",
      out_file = paste0(output_dir_figures, "significant_TF_gene_promoter_distance_heatmap_in_AD_list.png")
    )

    # Keep legacy filename for compatibility
    plot_gene_tf_heatmap(
      motif_nearest_df = motif_nearest,
      title_txt = "Closest motif/TF hit to significant promoters (AD TF list)",
      out_file = paste0(output_dir_figures, "significant_TF_gene_promoter_distance_heatmap.png")
    )

    # Top-10 significant TF targets with any human motif
    motif_dist_up_top10 <- build_motif_distance_df_with_lookup(
      seqs = seqs_up,
      prom_win = prom_up_win,
      direction_label = "Down in AD",
      motif_names = top10_hs_motif_names,
      pwm_lookup_local = pwm_lookup_top10,
      min_score = "85%"
    )

    motif_dist_down_top10 <- build_motif_distance_df_with_lookup(
      seqs = seqs_down,
      prom_win = prom_down_win,
      direction_label = "Up in AD",
      motif_names = top10_hs_motif_names,
      pwm_lookup_local = pwm_lookup_top10,
      min_score = "85%"
    )

    motif_dist_all_top10 <- bind_rows(motif_dist_up_top10, motif_dist_down_top10)

    if (nrow(motif_dist_all_top10) == 0) {
      message("No motif hits found for top-10 human motif TF set at min.score = 85%; skipping top-10 heatmap.")
    } else {
      motif_nearest_top10 <- motif_dist_all_top10 %>%
        group_by(gene_name, promoter_id, TF, motif_name, direction) %>%
        summarise(nearest_dist_bp = min(abs_dist_to_tss_bp), .groups = "drop")

      plot_gene_tf_heatmap(
        motif_nearest_df = motif_nearest_top10,
        title_txt = "Closest motif/TF hit to significant promoters (Top-10 TFs with human motifs)",
        out_file = paste0(output_dir_figures, "significant_TF_gene_promoter_distance_heatmap_top10_human_motif.png")
      )
    }

  }
}

mat_to_edges <- function(score_mat, prom_gr, direction_label, score_thresh = 6) {
  df <- as.data.frame(score_mat)
  df$promoter_id <- mcols(prom_gr)$promoter_id

  long <- tidyr::pivot_longer(df,
                             cols = -promoter_id,
                             names_to = "motif_name",
                             values_to = "score")

  long %>%
    filter(score >= score_thresh) %>%
    mutate(direction = direction_label) %>%
    # crude TF label extraction: first token up to _ or . or space
    mutate(TF = toupper(gsub("[\\._ ].*$", "", motif_name))) %>%
    dplyr::select(TF, promoter_id, score, direction, motif_name)
}
edges_up   <- mat_to_edges(score_mat_up_dense,   prom_up_win,   "Down in AD", score_thresh = 6)
edges_down <- mat_to_edges(score_mat_down_dense, prom_down_win, "Up in AD",   score_thresh = 6)
edges <- bind_rows(edges_up, edges_down)
head(edges)


top_tfs <- edges %>% count(TF, sort=TRUE) %>% slice_head(n=12) %>% pull(TF)

edges_f <- edges %>%
  #filter(TF %in% top_tfs) %>%
  group_by(TF, promoter_id) %>%
  summarise(score = max(score), direction = first(direction), .groups="drop")


# =========================
# A) Build TF -> gene edges
# =========================
library(dplyr)
library(tidyr)
library(circlize)

# 1) Start from promoter-level edges
chord_df <- edges_f %>%
  transmute(from = TF, to = promoter_id, value = score, direction = direction)

# 2) Clean TF labels (keep heterodimers; drop MA IDs; drop one junk label)
extract_TF <- function(x) {
  x <- toupper(x)
  x <- sub("^HSAPIENS-[^-]+-", "", x)     # remove source prefix
  x <- sub("-MA[0-9]+.*$", "", x)         # remove JASPAR MA IDs
  x <- sub("-[0-9]+$", "", x)             # collapse SOX9-2 -> SOX9, etc.
  x
}

normalize_ad_tf <- function(x) {
  x <- toupper(as.character(x))
  x <- sub("::.*$", "", x)                # keep first TF in heterodimers
  x <- sub("-[0-9]+$", "", x)
  x
}

ad_tfs_norm <- unique(normalize_ad_tf(ad_tfs))

chord_df2 <- chord_df %>%
  mutate(TF = vapply(from, extract_TF, character(1))) %>%
  filter(!is.na(TF), TF != "HSAPIENS-JASPAR") %>%
  filter(TF %in% ad_tfs_norm)

if (nrow(chord_df2) == 0) {
  stop("No TF-promoter edges remain for circos after filtering to ad_tfs list.")
}

# 3) Map promoter_id -> gene_name
chord_df2 <- chord_df2 %>%
  left_join(prom_dt %>% dplyr::select(promoter_id, gene_name), by = c("to" = "promoter_id")) %>%
  mutate(gene_name = ifelse(is.na(gene_name) | gene_name == "", to, gene_name))

# 4) Collapse promoter-level links to gene-level links
#    Choose max(value) for strongest motif evidence per TF–gene–direction
edges_gene <- chord_df2 %>%
  group_by(TF, gene_name, direction) %>%
  summarise(value = max(value), .groups = "drop")

# Link colors (direction overlay)
link_cols <- ifelse(edges_gene$direction == "Up in AD", "red",
                    ifelse(edges_gene$direction == "Down in AD", "blue", "grey70"))


# =============================================
# B) Hierarchical clustering of genes (k = 4)
# =============================================
# Build gene x TF matrix from edges_gene
mat <- edges_gene %>%
  group_by(gene_name, TF) %>%
  summarise(w = sum(value), .groups = "drop") %>%
  pivot_wider(names_from = TF, values_from = w, values_fill = 0)

gene_names <- mat$gene_name
X <- as.matrix(mat[, -1, drop = FALSE])
rownames(X) <- gene_names

hc <- hclust(dist(X), method = "ward.D2")
k <- 4
cl <- cutree(hc, k = k)                         # named vector: gene -> cluster
gene_order <- hc$labels[hc$order]               # dendrogram order


# ===========================================================
# C) Ensure sector order + gaps match EXACT plotted sectors
# ===========================================================
# Keep only genes/TFs actually present in edges_gene
tf_levels   <- sort(unique(edges_gene$TF))
gene_levels <- gene_order[gene_order %in% unique(edges_gene$gene_name)]

# Refactor edges to enforce sector order
edges_plot <- edges_gene %>%
  mutate(
    TF = factor(TF, levels = tf_levels),
    gene_name = factor(gene_name, levels = gene_levels)
  ) %>%
  filter(!is.na(TF), !is.na(gene_name))

# Drop unused levels so sector definitions exactly match plotted sectors
edges_plot$TF <- droplevels(edges_plot$TF)
edges_plot$gene_name <- droplevels(edges_plot$gene_name)

tf_levels <- levels(edges_plot$TF)
gene_levels <- levels(edges_plot$gene_name)

# Cluster-based gaps for genes (big gap when cluster changes)
cl_used <- cl[gene_levels]
gap_gene <- ifelse(c(TRUE, diff(cl_used) != 0), 10, 1)

# Gaps: within TFs + big separator + within genes with cluster breaks
gap_tf <- c(rep(2, length(tf_levels) - 1), 12)
gap_after <- c(gap_tf, gap_gene)

sector_ids <- unique(c(as.character(edges_plot$TF), as.character(edges_plot$gene_name)))
if (length(gap_after) != length(sector_ids)) {
  message(
    "Gap vector length (", length(gap_after),
    ") != sector count (", length(sector_ids),
    "); using uniform gap.after = 1 (possible TF/gene name overlap)."
  )
  gap_after <- 1
}


# =====================================
# D) Plot: blocks + links + big labels
# =====================================
# Sector block colors: TF grey, genes black
grid_col <- c(
  setNames(rep("grey70", length(tf_levels)), tf_levels),
  setNames(rep("black",  length(gene_levels)), gene_levels)
)


png(paste0(output_dir_figures, "TF_gene_circos_AD.png"), width = 3000, height = 3500, res = 300)
circos.clear()
circos.par(
  start.degree = 90,
  gap.after = gap_after,
  track.margin = c(0.02, 0.02)
)

chordDiagram(
  x = edges_plot %>% transmute(from = TF, to = gene_name, value = value),
  grid.col = grid_col,
  col = link_cols,
  transparency = 0.30,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.30)   # more space for bigger text
)

# Big, rotated labels on the same outer ring (aligned)
circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector <- get.cell.meta.data("sector.index")
    xlim   <- get.cell.meta.data("xlim")
    ylim   <- get.cell.meta.data("ylim")

    is_tf <- sector %in% tf_levels

    # TF labels (inside, bigger)
    if (is_tf) {
      circos.text(
        x = mean(xlim),
        y = 0.25,
        labels = sector,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = 1.0,
        col = "goldenrod4",
        font = 2
      )
    }

    # Gene labels (pushed outward)
    if (!is_tf) {
      circos.text(
        x = mean(xlim),
        y = 0.25,   # <-- push outside the track
        labels = sector,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = 1.0,
        col = "black",
        font = 2
      )
    }
  },
  bg.border = NA
)

dev.off()

# ==========================================
# E) Heatmap: genes with mixed Up/Down TF hits
# ==========================================
library(dplyr)
library(tidyr)
library(pheatmap)

# 1) Identify genes that have BOTH directions across TF links
mixed_genes <- edges_gene %>%
  group_by(gene_name) %>%
  summarise(has_up = any(direction == "Up in AD"),
            has_down = any(direction == "Down in AD"),
            .groups = "drop") %>%
  filter(has_up & has_down) %>%
  pull(gene_name)

heat_df <- chord_df2 %>%
  mutate(dir_num = case_when(
    direction == "Up in AD"   ~  1,
    direction == "Down in AD" ~ -1
  ))

heat_df <- heat_df %>%
  filter(gene_name %in% mixed_genes)

heat_df_collapsed <- heat_df %>%
  group_by(to, TF) %>%
  summarise(
    dir_num = ifelse(
      any(dir_num == 1) & any(dir_num == -1),
      0,
      sum(dir_num)
    ),
    gene_name = first(gene_name),
    .groups = "drop"
  )

heat_mat_df <- heat_df_collapsed %>%
  pivot_wider(
    names_from = TF,
    values_from = dir_num,
    values_fill = 0
  )

heat_mat <- as.matrix(heat_mat_df[, 3:ncol(heat_mat_df)])
rownames(heat_mat) <- heat_mat_df$to

# Binarize/clip values in heat_mat: > 1 -> 1, < -1 -> -1 (leave -1/0/1 as-is)
heat_mat[!is.na(heat_mat) & heat_mat > 1]  <- 1
heat_mat[!is.na(heat_mat) & heat_mat < -1] <- -1

# Remove TF columns that are empty (all 0 / NA) across all rows
keep_cols <- colSums(abs(replace(heat_mat, is.na(heat_mat), 0))) > 0
heat_mat <- heat_mat[, keep_cols, drop = FALSE]
row_anno <- data.frame(
  gene = heat_mat_df$gene_name
)

rownames(row_anno) <- rownames(heat_mat)

ph <- pheatmap(
  heat_mat,
  color = c("blue", "white", "red"), 
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_row = row_anno,
  show_rownames = TRUE,
  fontsize_col = 12,
  #border_color = NA,
  main = "TF-associated promoter usage changes in AD",
  silent = TRUE
)

# Add explicit blank space on the right so legend/labels are not clipped
ph$gtable <- gtable::gtable_add_cols(ph$gtable, widths = grid::unit(1.4, "in"), pos = 0)
ph$gtable <- gtable::gtable_add_cols(ph$gtable, widths = grid::unit(1.4, "in"))

ggsave(
  filename = paste0(output_dir_figures, "mixed_direction_genes_TF_heatmap.png"),
  plot = ph$gtable,
  width = 10,
  height = max(6, 0.2 * nrow(heat_mat) + 3),
  units = "in",
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)



# ==========================================
# ATAC-seq overlap with Promoters
# ==========================================
library(GenomicRanges)
library(stringr)

DAC_peaks <- read.csv("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn25671149/jaquino_analysis/tmp/promoter_usage_61_genes/Differential_ATAC.csv", header=TRUE)

dac_gr <- GRanges(
  seqnames = sub(":.*", "", DAC_peaks$Peak),
  ranges   = IRanges(
    start = as.integer(sub(".*:(\\d+)-.*", "\\1", DAC_peaks$Peak)),
    end   = as.integer(sub(".*-(\\d+)$", "\\1", DAC_peaks$Peak))
  )
)

ov_up <- findOverlaps(dac_gr, prom_up_win)
ov_down <- findOverlaps(dac_gr, prom_down_win)

overlap_up_table <- data.frame(
  dac_peak = DAC_peaks[queryHits(ov_up),],
  prom_up = prom_up_win[subjectHits(ov_up)]
)

overlap_down_table <- data.frame(
  dac_peak = DAC_peaks[queryHits(ov_down),],
  prom_down = prom_down_win[subjectHits(ov_down)]
)

tf_up_window_gr <- promoters(
  prom_up_win,
  upstream = 2000,
  downstream = 2000
)

tf_down_window_gr <- promoters(
  prom_down_win,
  upstream = 2000,
  downstream = 2000
)

ov_up <- findOverlaps(dac_gr, tf_up_window_gr)
ov_down <- findOverlaps(dac_gr, tf_down_window_gr)

overlap_up_table <- data.frame(
  dac_peak = DAC_peaks[queryHits(ov_up),],
  prom_up = tf_up_window_gr[subjectHits(ov_up)]
)

overlap_down_table <- data.frame(
  dac_peak = DAC_peaks[queryHits(ov_down),],
  prom_down = tf_down_window_gr[subjectHits(ov_down)]
)

tf_up_window_gr <- resize(
  prom_up_win,
  width = width(prom_up_win) + 4000,
  fix = "center"
)

tf_down_window_gr <- resize(
  prom_down_win,
  width = width(prom_down_win) + 4000,
  fix = "center"
)

ov_up <- findOverlaps(dac_gr, tf_up_window_gr)
ov_down <- findOverlaps(dac_gr, tf_down_window_gr)

overlap_up_table <- data.frame(
  dac_peak = DAC_peaks[queryHits(ov_up),],
  prom_up = tf_up_window_gr[subjectHits(ov_up)]
)

overlap_down_table <- data.frame(
  dac_peak = DAC_peaks[queryHits(ov_down),],
  prom_down = tf_down_window_gr[subjectHits(ov_down)]
)

overlap_up_table$direction <- "Promoter Up"
overlap_down_table$direction <- "Promoter Down"

df <- bind_rows(
  overlap_up_table %>%
    dplyr::select(gene = prom_up.gene_name,
           promoter_fc = prom_up.logitFC,
           atac_fc = dac_peak.avg_logFC,
           cell_type = dac_peak.cell_type,
           direction),

  overlap_down_table %>%
    dplyr::select(gene = prom_down.gene_name,
           promoter_fc = prom_down.logitFC,
           atac_fc = dac_peak.avg_logFC,
           cell_type = dac_peak.cell_type,
           direction)
) %>%
  distinct(gene, cell_type, .keep_all = TRUE)

ggplot(df, aes(x = promoter_fc,
               y = atac_fc,
               color = cell_type,
               shape = gene)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_shape_manual(values = c(
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
    10, 11, 12, 13, 14,
    15, 16, 17, 18, 19
  )) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  theme(
    plot.margin = margin(t = 30, r = 10, b = 30, l = 10)
  ) +
  guides(
    color = guide_legend(ncol = 2),
    shape = guide_legend(ncol = 2)
  )

ggsave(paste0(output_dir_figures, "promoter_usage_ATAC_overlap_scatter.png"), width = 8, height = 6, dpi = 300)

df$gene_group <- ifelse(df$gene %in% c("PPM1B","PAXX","REX1BD"),
                        df$gene,
                        "Other")

ggplot(df, aes(x = promoter_fc,
               y = atac_fc,
               color = cell_type,
               shape = gene_group)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()
ggsave(paste0(output_dir_figures, "promoter_usage_ATAC_overlap_scatter_subset.png"), width = 8, height = 6, dpi = 300)