#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(betareg)
  library(broom)
  library(ggplot2)
})

# Differential PSI beta regression for SUPPA event types:
# SE, RI, MX, AL, AF, A5, A3.
#
# The script:
#   1. Matches significant ENST IDs to SUPPA IOE events.
#   2. Runs event-level beta regression of PSI ~ condition + sex.
#   3. Adds compact plot labels based on SUPPA event coordinates.
#   4. Writes per-event-type and combined tables/figures.

###############################
# User configuration
###############################

BASE_DIR <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis"
SUPPA_DIR <- file.path(BASE_DIR, "suppa_bambu_gtf")

SIG_TX_FILE <- file.path(BASE_DIR, "diu", "sig_transcripts.txt")
METADATA_FILE <- file.path(BASE_DIR, "long_read_metadata.txt")

IOE_TEMPLATE <- file.path(SUPPA_DIR, "FC_events_{EVENT}_strict.ioe")
PSI_TEMPLATE <- file.path(SUPPA_DIR, "psi_{EVENT}.psi")

REFERENCE_GTF_FILE <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/reproducing_nextflow_pipeline/Homo_sapiens.GRCh38.107_ERCC.gtf"
BAMBU_GTF_FILE <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/NanoporeRNASeq_extended_annotations.gtf"

FIGURE_DIR <- file.path(BASE_DIR, "figures", "suppa_beta_regression")
TABLE_DIR <- file.path(BASE_DIR, "tables", "suppa_beta_regression")

EVENT_TYPES <- c("SE", "RI", "MX", "AL", "AF", "A5", "A3")

CONDITION_REFERENCE <- "CT"
CONDITION_TEST <- "AD"
MIN_NON_MISSING_PER_GROUP <- 3
FDR_CUTOFF <- 0.05
DPSI_CUTOFF <- 0.10

# These defaults are set to reproduce the behavior of your original SE script:
#   sig_tx <- read_tsv(...)
#   inner_join(sig_tx, by = "transcript_id")
SIG_TX_HAS_HEADER <- TRUE
MATCH_TRANSCRIPT_ID_BY_BASE <- FALSE

# TRUE reproduces the coordinate style you asked for, e.g.
# TRIT1 | SE | exon 39844532-39844639
# These are the internal affected bases between SUPPA splice-site coordinates.
USE_INTERNAL_EVENT_INTERVALS <- TRUE

# Set TRUE if you want labels like chr1:39844532-39844639.
INCLUDE_CHR_IN_PLOT_LABEL <- FALSE

# Keep this TRUE if the beta regression should only test events linked to your
# significant IsoformSwitchAnalyzeR transcripts. Set FALSE to test every SUPPA
# event with PSI values.
ANALYZE_ONLY_EVENTS_LINKED_TO_SIG_TX <- TRUE

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

###############################
# Helpers
###############################

path_for_event <- function(template, event_type) {
  stringr::str_replace(template, "\\{EVENT\\}", event_type)
}

strip_tx_version <- function(x) {
  stringr::str_remove(x, "\\.\\d+$")
}

clean_sample_name <- function(x) {
  make.names(x)
}

split_transcripts <- function(x) {
  if (is.na(x) || x == "") {
    return(character())
  }

  stringr::str_split(x, ",", simplify = FALSE)[[1]] %>%
    stringr::str_trim() %>%
    discard(~ .x == "")
}

read_sig_transcripts <- function(path) {
  sig_tx <- readr::read_tsv(path, col_names = SIG_TX_HAS_HEADER, show_col_types = FALSE)

  sig_tx %>%
    dplyr::select(transcript_id = 1) %>%
    mutate(transcript_id = as.character(transcript_id)) %>%
    filter(!is.na(transcript_id), transcript_id != "", transcript_id != "transcript_id") %>%
    mutate(
      tx_key = if (MATCH_TRANSCRIPT_ID_BY_BASE) {
        strip_tx_version(transcript_id)
      } else {
        transcript_id
      }
    ) %>%
    distinct(tx_key, .keep_all = TRUE)
}

read_metadata <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE) %>%
    mutate(
      sample = clean_sample_name(sample),
      condition = factor(condition, levels = c(CONDITION_REFERENCE, CONDITION_TEST)),
      sex = factor(sex)
    )
}

read_psi <- function(path) {
  psi <- read.table(
    path,
    row.names = 1,
    header = TRUE,
    sep = "\t",
    check.names = TRUE,
    stringsAsFactors = FALSE
  )

  psi <- as_tibble(psi, rownames = "event_id")
  colnames(psi) <- clean_sample_name(colnames(psi))
  colnames(psi)[colnames(psi) == "event_id"] <- "event_id"

  psi %>%
    mutate(event_id = as.character(event_id))
}

read_gene_symbols <- function(reference_gtf_file, bambu_gtf_file) {
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    warning("Package 'rtracklayer' is not installed; plot labels will use gene_id.")
    return(tibble(gene_id = character(), gene_symbol = character()))
  }

  symbol_tbl <- tibble(gene_id = character(), gene_symbol = character())

  if (file.exists(reference_gtf_file)) {
    message("Reading reference GTF gene symbols: ", reference_gtf_file)
    ref_gtf <- rtracklayer::import(reference_gtf_file) %>%
      as.data.frame() %>%
      as_tibble()

    if ("gene_name" %in% colnames(ref_gtf)) {
      symbol_tbl <- ref_gtf %>%
        filter(!is.na(gene_id), !is.na(gene_name)) %>%
        distinct(gene_id, gene_symbol = gene_name)
    } else if ("gene_symbol" %in% colnames(ref_gtf)) {
      symbol_tbl <- ref_gtf %>%
        filter(!is.na(gene_id), !is.na(gene_symbol)) %>%
        distinct(gene_id, gene_symbol)
    }
  }

  if (file.exists(bambu_gtf_file)) {
    message("Reading Bambu GTF for gene IDs: ", bambu_gtf_file)
    bambu_gtf <- rtracklayer::import(bambu_gtf_file) %>%
      as.data.frame() %>%
      as_tibble()

    if ("gene_name" %in% colnames(bambu_gtf)) {
      bambu_symbols <- bambu_gtf %>%
        filter(!is.na(gene_id), !is.na(gene_name)) %>%
        distinct(gene_id, gene_symbol = gene_name)
    } else if ("gene_symbol" %in% colnames(bambu_gtf)) {
      bambu_symbols <- bambu_gtf %>%
        filter(!is.na(gene_id), !is.na(gene_symbol)) %>%
        distinct(gene_id, gene_symbol)
    } else {
      bambu_symbols <- bambu_gtf %>%
        filter(!is.na(gene_id)) %>%
        distinct(gene_id) %>%
        left_join(symbol_tbl, by = "gene_id")
    }

    symbol_tbl <- bind_rows(symbol_tbl, bambu_symbols) %>%
      filter(!is.na(gene_id), !is.na(gene_symbol)) %>%
      distinct(gene_id, .keep_all = TRUE)
  }

  symbol_tbl
}

format_chr <- function(chr) {
  ifelse(stringr::str_detect(chr, "^chr"), chr, paste0("chr", chr))
}

format_region <- function(chr, start, end, include_chr = INCLUDE_CHR_IN_PLOT_LABEL) {
  if (is.na(start) || is.na(end)) {
    return("coord?")
  }

  start <- as.integer(start)
  end <- as.integer(end)

  if (start > end) {
    tmp <- start
    start <- end
    end <- tmp
  }

  region <- paste0(start, "-", end)

  if (isTRUE(include_chr)) {
    region <- paste0(format_chr(chr), ":", region)
  }

  region
}

parse_coord_token <- function(token) {
  pieces <- strsplit(token, "-", fixed = TRUE)[[1]]
  nums <- suppressWarnings(as.integer(pieces))

  if (length(nums) == 1) {
    tibble(raw = token, start = nums[1], end = nums[1], is_range = FALSE)
  } else if (length(nums) == 2) {
    tibble(raw = token, start = nums[1], end = nums[2], is_range = TRUE)
  } else {
    tibble(raw = token, start = NA_integer_, end = NA_integer_, is_range = NA)
  }
}

empty_event_detail <- function(reason) {
  tibble(
    region_type = NA_character_,
    region_1_start = NA_integer_,
    region_1_end = NA_integer_,
    region_2_start = NA_integer_,
    region_2_end = NA_integer_,
    plot_detail = "coords?",
    parse_note = reason
  )
}

make_se_detail <- function(chr, coords) {
  if (length(coords) != 2) {
    return(empty_event_detail("SE event did not have 2 coordinate tokens."))
  }

  c1 <- parse_coord_token(coords[1])
  c2 <- parse_coord_token(coords[2])

  if (USE_INTERNAL_EVENT_INTERVALS) {
    region_start <- c1$end + 1
    region_end <- c2$start - 1
  } else {
    region_start <- c1$end
    region_end <- c2$start
  }

  tibble(
    region_type = "skipped_exon",
    region_1_start = region_start,
    region_1_end = region_end,
    region_2_start = NA_integer_,
    region_2_end = NA_integer_,
    plot_detail = paste0("exon ", format_region(chr, region_start, region_end)),
    parse_note = NA_character_
  )
}

make_ri_detail <- function(chr, coords) {
  if (length(coords) != 3) {
    return(empty_event_detail("RI event did not have 3 coordinate tokens."))
  }

  intron <- parse_coord_token(coords[2])

  if (USE_INTERNAL_EVENT_INTERVALS) {
    region_start <- intron$start + 1
    region_end <- intron$end - 1
  } else {
    region_start <- intron$start
    region_end <- intron$end
  }

  tibble(
    region_type = "retained_intron",
    region_1_start = region_start,
    region_1_end = region_end,
    region_2_start = NA_integer_,
    region_2_end = NA_integer_,
    plot_detail = paste0("intron ", format_region(chr, region_start, region_end)),
    parse_note = NA_character_
  )
}

make_mx_detail <- function(chr, coords) {
  if (length(coords) != 4) {
    return(empty_event_detail("MX event did not have 4 coordinate tokens."))
  }

  c1 <- parse_coord_token(coords[1])
  c2 <- parse_coord_token(coords[2])
  c3 <- parse_coord_token(coords[3])
  c4 <- parse_coord_token(coords[4])

  if (USE_INTERNAL_EVENT_INTERVALS) {
    exon_1_start <- c1$end + 1
    exon_1_end <- c2$start - 1
    exon_2_start <- c3$end + 1
    exon_2_end <- c4$start - 1
  } else {
    exon_1_start <- c1$end
    exon_1_end <- c2$start
    exon_2_start <- c3$end
    exon_2_end <- c4$start
  }

  tibble(
    region_type = "mutually_exclusive_exons",
    region_1_start = exon_1_start,
    region_1_end = exon_1_end,
    region_2_start = exon_2_start,
    region_2_end = exon_2_end,
    plot_detail = paste0(
      "exons ",
      format_region(chr, exon_1_start, exon_1_end),
      " vs ",
      format_region(chr, exon_2_start, exon_2_end)
    ),
    parse_note = NA_character_
  )
}

make_a35_detail <- function(chr, coords, event_type) {
  if (length(coords) != 2) {
    return(empty_event_detail(paste0(event_type, " event did not have 2 coordinate tokens.")))
  }

  c1 <- parse_coord_token(coords[1])
  c2 <- parse_coord_token(coords[2])

  if (isTRUE(c1$start != c2$start) && isTRUE(c1$end == c2$end)) {
    site_1 <- c1$start
    site_2 <- c2$start
    region_start <- min(site_1, site_2) + 1
    region_end <- max(site_1, site_2)
    note <- NA_character_
  } else if (isTRUE(c1$start == c2$start) && isTRUE(c1$end != c2$end)) {
    site_1 <- c1$end
    site_2 <- c2$end
    region_start <- min(site_1, site_2)
    region_end <- max(site_1, site_2) - 1
    note <- NA_character_
  } else {
    site_1 <- NA_integer_
    site_2 <- NA_integer_
    region_start <- min(c(c1$start, c1$end, c2$start, c2$end), na.rm = TRUE)
    region_end <- max(c(c1$start, c1$end, c2$start, c2$end), na.rm = TRUE)
    note <- "Both splice-site coordinates differed; label uses full coordinate span."
  }

  if (!is.na(region_start) && !is.na(region_end) && region_start <= region_end) {
    detail <- paste0("splice region ", format_region(chr, region_start, region_end))
  } else {
    sites <- paste(sort(unique(c(site_1, site_2))), collapse = "/")
    detail <- paste0("splice sites ", sites)
  }

  tibble(
    region_type = "alternative_splice_site",
    region_1_start = region_start,
    region_1_end = region_end,
    region_2_start = NA_integer_,
    region_2_end = NA_integer_,
    plot_detail = detail,
    parse_note = note
  )
}

make_fl_detail <- function(chr, coords, event_type) {
  if (length(coords) != 4) {
    return(empty_event_detail(paste0(event_type, " event did not have 4 coordinate tokens.")))
  }

  c1 <- parse_coord_token(coords[1])
  c2 <- parse_coord_token(coords[2])
  c3 <- parse_coord_token(coords[3])
  c4 <- parse_coord_token(coords[4])

  pattern <- paste(
    ifelse(c(c1$is_range, c2$is_range, c3$is_range, c4$is_range), "R", "S"),
    collapse = ""
  )

  if (pattern == "SRSR") {
    terminal_1_start <- c1$start
    terminal_1_end <- c2$start
    terminal_2_start <- c3$start
    terminal_2_end <- c4$start
    note <- NA_character_
  } else if (pattern == "RSRS") {
    terminal_1_start <- c1$end
    terminal_1_end <- c2$start
    terminal_2_start <- c3$end
    terminal_2_end <- c4$start
    note <- NA_character_
  } else {
    terminal_1_start <- NA_integer_
    terminal_1_end <- NA_integer_
    terminal_2_start <- NA_integer_
    terminal_2_end <- NA_integer_
    note <- paste0("Unexpected ", event_type, " coordinate pattern: ", pattern)
  }

  if (is.na(terminal_1_start) || is.na(terminal_2_start)) {
    detail <- paste0("terminal coords ", paste(coords, collapse = " vs "))
  } else {
    detail <- paste0(
      "terminal exons ",
      format_region(chr, terminal_1_start, terminal_1_end),
      " vs ",
      format_region(chr, terminal_2_start, terminal_2_end)
    )
  }

  tibble(
    region_type = "alternative_terminal_exons",
    region_1_start = terminal_1_start,
    region_1_end = terminal_1_end,
    region_2_start = terminal_2_start,
    region_2_end = terminal_2_end,
    plot_detail = detail,
    parse_note = note
  )
}

make_event_detail <- function(event_type, chr, coords) {
  if (event_type == "SE") {
    make_se_detail(chr, coords)
  } else if (event_type == "RI") {
    make_ri_detail(chr, coords)
  } else if (event_type == "MX") {
    make_mx_detail(chr, coords)
  } else if (event_type %in% c("A3", "A5")) {
    make_a35_detail(chr, coords, event_type)
  } else if (event_type %in% c("AF", "AL")) {
    make_fl_detail(chr, coords, event_type)
  } else {
    empty_event_detail(paste0("Unsupported event type: ", event_type))
  }
}

parse_suppa_event_one <- function(event_id) {
  event_split <- strsplit(event_id, ";", fixed = TRUE)[[1]]

  if (length(event_split) < 2) {
    return(tibble(
      event_id = event_id,
      gene_id_from_event = NA_character_,
      event_type = NA_character_,
      chr = NA_character_,
      strand = NA_character_,
      coord_string = NA_character_
    ) %>%
      bind_cols(empty_event_detail("Could not split event_id on ';'.")))
  }

  gene_id_from_event <- event_split[1]
  event_body <- paste(event_split[-1], collapse = ";")
  tokens <- strsplit(event_body, ":", fixed = TRUE)[[1]]

  if (length(tokens) < 4) {
    return(tibble(
      event_id = event_id,
      gene_id_from_event = gene_id_from_event,
      event_type = NA_character_,
      chr = NA_character_,
      strand = NA_character_,
      coord_string = event_body
    ) %>%
      bind_cols(empty_event_detail("Could not parse event body.")))
  }

  event_type <- tokens[1]
  chr <- tokens[2]
  strand <- tokens[length(tokens)]
  coords <- tokens[3:(length(tokens) - 1)]
  detail <- make_event_detail(event_type, chr, coords)

  tibble(
    event_id = event_id,
    gene_id_from_event = gene_id_from_event,
    event_type = event_type,
    chr = chr,
    strand = strand,
    coord_string = paste(coords, collapse = ":")
  ) %>%
    bind_cols(detail)
}

parse_suppa_events <- function(event_ids) {
  purrr::map_dfr(unique(event_ids), parse_suppa_event_one)
}

make_event_transcript_matches <- function(ioe, event_type, sig_tx) {
  event_tx <- ioe %>%
    mutate(
      event_type = event_type,
      alternative_list = purrr::map(alternative_transcripts, split_transcripts),
      total_list = purrr::map(total_transcripts, split_transcripts),
      opposite_list = purrr::map2(total_list, alternative_list, setdiff)
    )

  alt_long <- event_tx %>%
    dplyr::select(event_type, gene_id, event_id, transcript_id = alternative_list) %>%
    tidyr::unnest(transcript_id) %>%
    mutate(role = "alternative")

  opposite_long <- event_tx %>%
    dplyr::select(event_type, gene_id, event_id, transcript_id = opposite_list) %>%
    tidyr::unnest(transcript_id) %>%
    mutate(role = "opposite")

  bind_rows(alt_long, opposite_long) %>%
    mutate(
      tx_key = if (MATCH_TRANSCRIPT_ID_BY_BASE) {
        strip_tx_version(transcript_id)
      } else {
        transcript_id
      }
    ) %>%
    inner_join(
      sig_tx %>% dplyr::select(sig_transcript_id = transcript_id, tx_key),
      by = "tx_key"
    ) %>%
    distinct(event_type, gene_id, event_id, transcript_id, sig_transcript_id, tx_key, role, .keep_all = TRUE)
}

run_beta_one_event <- function(df_event, meta, sample_cols) {
  long <- df_event %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "sample",
      values_to = "psi"
    ) %>%
    left_join(meta, by = "sample") %>%
    mutate(psi = as.numeric(psi)) %>%
    filter(is.finite(psi), !is.na(condition))

  n_reference <- sum(long$condition == CONDITION_REFERENCE)
  n_test <- sum(long$condition == CONDITION_TEST)

  if (n_reference < MIN_NON_MISSING_PER_GROUP || n_test < MIN_NON_MISSING_PER_GROUP) {
    return(tibble(
      term = paste0("condition", CONDITION_TEST),
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      n_reference = n_reference,
      n_test = n_test,
      model_used = NA_character_,
      model_note = "Too few non-missing PSI values in one or both groups."
    ))
  }

  n <- nrow(long)
  long <- long %>%
    mutate(psi_beta = (psi * (n - 1) + 0.5) / n)

  has_sex <- n_distinct(long$sex) >= 2

  full_formula <- if (has_sex) {
    psi_beta ~ condition + sex
  } else {
    psi_beta ~ condition
  }

  reduced_formula <- if (has_sex) {
    psi_beta ~ sex
  } else {
    psi_beta ~ 1
  }

  model_used <- paste(
    "full:", paste(deparse(full_formula), collapse = " "),
    "| reduced:", paste(deparse(reduced_formula), collapse = " ")
  )

  fit_full <- tryCatch(
    betareg::betareg(full_formula, data = long),
    error = function(e) e
  )

  if (inherits(fit_full, "error")) {
    return(tibble(
      term = paste0("condition", CONDITION_TEST),
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      n_reference = n_reference,
      n_test = n_test,
      model_used = model_used,
      model_note = paste0("Full model failed: ", fit_full$message)
    ))
  }

  fit_reduced <- tryCatch(
    betareg::betareg(reduced_formula, data = long),
    error = function(e) e
  )

  if (inherits(fit_reduced, "error")) {
    return(tibble(
      term = paste0("condition", CONDITION_TEST),
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      n_reference = n_reference,
      n_test = n_test,
      model_used = model_used,
      model_note = paste0("Reduced model failed: ", fit_reduced$message)
    ))
  }

  fit_tidy <- broom::tidy(fit_full)

  if ("component" %in% colnames(fit_tidy)) {
    fit_tidy <- fit_tidy %>%
      filter(component == "mean")
  }

  condition_row <- fit_tidy %>%
    filter(term == paste0("condition", CONDITION_TEST)) %>%
    dplyr::select(any_of(c("term", "estimate", "std.error"))) %>%
    slice(1)

  if (nrow(condition_row) == 0) {
    return(tibble(
      term = paste0("condition", CONDITION_TEST),
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      n_reference = n_reference,
      n_test = n_test,
      model_used = model_used,
      model_note = "condition term not found in beta regression output."
    ))
  }

  # Likelihood ratio test: full model (condition + sex) vs reduced model (sex only)
  df_diff <- length(coef(fit_full)) - length(coef(fit_reduced))
  lr_statistic <- 2 * (as.numeric(logLik(fit_full)) - as.numeric(logLik(fit_reduced)))
  lr_p_value <- pchisq(lr_statistic, df = df_diff, lower.tail = FALSE)

  condition_row %>%
    mutate(
      statistic = lr_statistic,
      p.value = lr_p_value,
      n_reference = n_reference,
      n_test = n_test,
      model_used = model_used,
      model_note = NA_character_
    )
}

run_beta_on_events <- function(psi_event, meta, sample_cols) {
  if (nrow(psi_event) == 0) {
    return(tibble())
  }

  beta_results <- psi_event %>%
    group_by(event_type, gene_id, event_id) %>%
    group_modify(~ run_beta_one_event(.x, meta, sample_cols)) %>%
    ungroup() %>%
    group_by(event_type) %>%
    mutate(adj_p_value = p.adjust(p.value, method = "BH")) %>%
    ungroup()

  ct_samples <- meta %>%
    filter(condition == CONDITION_REFERENCE) %>%
    pull(sample) %>%
    intersect(sample_cols)

  test_samples <- meta %>%
    filter(condition == CONDITION_TEST) %>%
    pull(sample) %>%
    intersect(sample_cols)

  effect_sizes <- psi_event %>%
    rowwise() %>%
    mutate(
      mean_CT = mean(c_across(all_of(ct_samples)), na.rm = TRUE),
      mean_AD = mean(c_across(all_of(test_samples)), na.rm = TRUE),
      dPSI = mean_AD - mean_CT
    ) %>%
    ungroup() %>%
    dplyr::select(event_type, gene_id, event_id, mean_CT, mean_AD, dPSI)

  beta_results %>%
    left_join(effect_sizes, by = c("event_type", "gene_id", "event_id")) %>%
    relocate(mean_CT, mean_AD, dPSI, .after = p.value)
}

plot_sig_events <- function(plot_df, event_type, figure_path) {
  if (nrow(plot_df) == 0) {
    return(invisible(NULL))
  }

  n_facets <- n_distinct(plot_df$plot_label)
  n_cols <- min(3, n_facets)
  n_rows <- ceiling(n_facets / n_cols)
  fig_width <- max(8, 5 * n_cols)
  fig_height <- max(4.5, 3.2 * n_rows)

  p <- ggplot(plot_df, aes(x = condition, y = psi, shape = sex)) +
    geom_point(
      size = 4,
      alpha = 0.85,
      position = position_jitterdodge(
        jitter.width = 0.12,
        jitter.height = 0,
        dodge.width = 0.45
      )
    ) +
    scale_shape_manual(values = c("Male" = "\u2642", "Female" = "\u2641")) +
    stat_summary(
      aes(group = condition),
      fun = mean,
      geom = "crossbar",
      width = 0.35,
      color = "black",
      linewidth = 0.25
    ) +
    facet_wrap(~ plot_label, scales = "fixed", ncol = n_cols) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
    theme_bw(base_size = 12) +
    theme(
      strip.text = element_text(size = 9),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = paste0(event_type, " significant differential PSI events"),
      x = NULL,
      y = "PSI",
      shape = "Sex"
    )

  ggsave(
    filename = figure_path,
    plot = p,
    width = fig_width,
    height = fig_height,
    dpi = 300,
    limitsize = FALSE
  )

  invisible(p)
}

process_event_type <- function(event_type, sig_tx, meta, gene_symbols) {
  message("Processing ", event_type, " events")

  ioe_path <- path_for_event(IOE_TEMPLATE, event_type)
  psi_path <- path_for_event(PSI_TEMPLATE, event_type)

  if (!file.exists(ioe_path)) {
    warning("Missing IOE file for ", event_type, ": ", ioe_path)
    return(NULL)
  }

  if (!file.exists(psi_path)) {
    warning("Missing PSI file for ", event_type, ": ", psi_path)
    return(NULL)
  }

  ioe <- readr::read_tsv(ioe_path, show_col_types = FALSE) %>%
    mutate(
      event_type = event_type,
      gene_id = as.character(gene_id),
      event_id = as.character(event_id)
    )

  psi <- read_psi(psi_path)

  sample_cols <- intersect(meta$sample, colnames(psi))
  missing_samples <- setdiff(meta$sample, colnames(psi))

  if (length(sample_cols) == 0) {
    stop("No metadata samples were found in PSI file: ", psi_path)
  }

  if (length(missing_samples) > 0) {
    warning(
      "These metadata samples were not found in ",
      event_type,
      " PSI file and will be ignored: ",
      paste(missing_samples, collapse = ", ")
    )
  }

  meta_event <- meta %>%
    filter(sample %in% sample_cols) %>%
    arrange(match(sample, sample_cols))

  sample_cols <- meta_event$sample

  psi_event_all <- ioe %>%
    dplyr::select(event_type, gene_id, event_id) %>%
    distinct() %>%
    inner_join(psi, by = "event_id") %>%
    mutate(across(all_of(sample_cols), as.numeric))

  tx_matches <- make_event_transcript_matches(ioe, event_type, sig_tx)

  readr::write_tsv(
    tx_matches,
    file.path(TABLE_DIR, paste0(event_type, "_transcript_matches.tsv"))
  )

  if (ANALYZE_ONLY_EVENTS_LINKED_TO_SIG_TX) {
    psi_event <- psi_event_all %>%
      semi_join(tx_matches, by = c("event_type", "gene_id", "event_id"))
  } else {
    psi_event <- psi_event_all
  }

  beta_results <- run_beta_on_events(psi_event, meta_event, sample_cols)

  if (nrow(beta_results) == 0) {
    readr::write_tsv(
      beta_results,
      file.path(TABLE_DIR, paste0(event_type, "_beta_results.tsv"))
    )

    return(list(
      event_type = event_type,
      matches = tx_matches,
      beta_results = beta_results,
      sig_beta_events = beta_results,
      sig_event_tx = tibble(),
      labels = tibble(),
      plot_df = tibble()
    ))
  }

  event_regions <- parse_suppa_events(beta_results$event_id) %>%
    dplyr::select(
      event_id,
      parsed_event_type = event_type,
      chr,
      strand,
      coord_string,
      region_type,
      region_1_start,
      region_1_end,
      region_2_start,
      region_2_end,
      plot_detail,
      parse_note
    )

  event_labels <- beta_results %>%
    distinct(event_type, gene_id, event_id) %>%
    left_join(event_regions, by = "event_id") %>%
    left_join(gene_symbols, by = "gene_id") %>%
    mutate(
      gene_symbol = coalesce(gene_symbol, gene_id),
      plot_label = paste0(gene_symbol, " | ", event_type, "\n", plot_detail)
    )

  beta_results <- beta_results %>%
    left_join(
      event_labels %>%
        dplyr::select(
          event_type,
          gene_id,
          event_id,
          gene_symbol,
          plot_label,
          chr,
          strand,
          coord_string,
          region_type,
          region_1_start,
          region_1_end,
          region_2_start,
          region_2_end,
          parse_note
        ),
      by = c("event_type", "gene_id", "event_id")
    ) %>%
    relocate(gene_symbol, plot_label, .after = event_id)

  sig_beta_events <- beta_results %>%
    filter(
      !is.na(adj_p_value),
      adj_p_value < FDR_CUTOFF,
      abs(dPSI) >= DPSI_CUTOFF
    ) %>%
    arrange(adj_p_value, desc(abs(dPSI)))

  sig_event_tx <- tx_matches %>%
    inner_join(
      sig_beta_events %>%
        dplyr::select(event_type, gene_id, event_id, gene_symbol, plot_label, mean_CT, mean_AD, dPSI, p.value, adj_p_value),
      by = c("event_type", "gene_id", "event_id")
    ) %>%
    arrange(adj_p_value, event_id, role, transcript_id)

  plot_df <- psi_event %>%
    filter(event_id %in% sig_beta_events$event_id) %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "sample",
      values_to = "psi"
    ) %>%
    left_join(meta_event, by = "sample") %>%
    left_join(
      sig_beta_events %>%
        dplyr::select(event_type, gene_id, event_id, gene_symbol, plot_label, mean_CT, mean_AD, dPSI, adj_p_value),
      by = c("event_type", "gene_id", "event_id")
    ) %>%
    mutate(
      psi = as.numeric(psi),
      plot_label = factor(plot_label, levels = unique(sig_beta_events$plot_label))
    )

  readr::write_tsv(
    beta_results,
    file.path(TABLE_DIR, paste0(event_type, "_beta_results.tsv"))
  )

  readr::write_tsv(
    sig_beta_events,
    file.path(TABLE_DIR, paste0(event_type, "_sig_beta_events.tsv"))
  )

  readr::write_tsv(
    sig_event_tx,
    file.path(TABLE_DIR, paste0(event_type, "_sig_beta_event_transcript_matches.tsv"))
  )

  readr::write_tsv(
    event_labels,
    file.path(TABLE_DIR, paste0(event_type, "_event_label_lookup.tsv"))
  )

  readr::write_tsv(
    plot_df,
    file.path(TABLE_DIR, paste0(event_type, "_sig_beta_events_plot_data.tsv"))
  )

  if (nrow(sig_beta_events) > 0) {
    plot_sig_events(
      plot_df = plot_df,
      event_type = event_type,
      figure_path = file.path(FIGURE_DIR, paste0(event_type, "_sig_beta_events.png"))
    )
  }

  list(
    event_type = event_type,
    matches = tx_matches,
    beta_results = beta_results,
    sig_beta_events = sig_beta_events,
    sig_event_tx = sig_event_tx,
    labels = event_labels,
    plot_df = plot_df
  )
}

###############################
# Main analysis
###############################

message("Reading significant transcripts")
sig_tx <- read_sig_transcripts(SIG_TX_FILE)

message("Reading metadata")
meta <- read_metadata(METADATA_FILE)

message("Reading gene symbols")
gene_symbols <- read_gene_symbols(REFERENCE_GTF_FILE, BAMBU_GTF_FILE)

results <- EVENT_TYPES %>%
  set_names() %>%
  purrr::map(~ process_event_type(.x, sig_tx, meta, gene_symbols)) %>%
  compact()

all_matches <- purrr::map_dfr(results, "matches")
all_beta_results <- purrr::map_dfr(results, "beta_results")
all_sig_beta_events <- purrr::map_dfr(results, "sig_beta_events")
all_sig_event_tx <- purrr::map_dfr(results, "sig_event_tx")
all_labels <- purrr::map_dfr(results, "labels")
all_plot_df <- purrr::map_dfr(results, "plot_df")

readr::write_tsv(
  all_matches,
  file.path(TABLE_DIR, "all_event_transcript_matches.tsv")
)

readr::write_tsv(
  all_beta_results,
  file.path(TABLE_DIR, "all_beta_results.tsv")
)

readr::write_tsv(
  all_sig_beta_events,
  file.path(TABLE_DIR, "all_sig_beta_events.tsv")
)

readr::write_tsv(
  all_sig_event_tx,
  file.path(TABLE_DIR, "all_sig_beta_event_transcript_matches.tsv")
)

readr::write_tsv(
  all_labels,
  file.path(TABLE_DIR, "all_event_label_lookup.tsv")
)

readr::write_tsv(
  all_plot_df,
  file.path(TABLE_DIR, "all_sig_beta_events_plot_data.tsv")
)

if (nrow(all_plot_df) > 0) {
  n_facets <- n_distinct(all_plot_df$plot_label)
  n_cols <- min(3, n_facets)
  n_rows <- ceiling(n_facets / n_cols)

  combined_plot <- ggplot(all_plot_df, aes(x = condition, y = psi, shape = sex)) +
    geom_point(
      size = 4,
      alpha = 0.85,
      position = position_jitterdodge(
        jitter.width = 0.12,
        jitter.height = 0,
        dodge.width = 0.45
      )
    ) +
    scale_shape_manual(values = c("Male" = "\u2642", "Female" = "\u2641")) +
    stat_summary(
      aes(group = condition),
      fun = mean,
      geom = "crossbar",
      width = 0.35,
      color = "black",
      linewidth = 0.25
    ) +
    facet_wrap(~ plot_label, scales = "fixed", ncol = n_cols) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
    theme_bw(base_size = 12) +
    theme(
      strip.text = element_text(size = 8),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "Significant differential PSI events",
      x = NULL,
      y = "PSI",
      shape = "Sex"
    )

  ggsave(
    filename = file.path(FIGURE_DIR, "all_sig_beta_events.png"),
    plot = combined_plot,
    width = max(8, 5 * n_cols),
    height = max(4.5, 3.2 * n_rows),
    dpi = 300,
    limitsize = FALSE
  )
}

summary_table <- tibble(event_type = EVENT_TYPES) %>%
  left_join(
    all_beta_results %>%
      group_by(event_type) %>%
      summarise(
        tested_events = n(),
        significant_events = sum(!is.na(adj_p_value) & adj_p_value < FDR_CUTOFF & abs(dPSI) >= DPSI_CUTOFF),
        .groups = "drop"
      ),
    by = "event_type"
  ) %>%
  left_join(
    all_matches %>%
      group_by(event_type) %>%
      summarise(
        transcript_matched_events = n_distinct(event_id),
        matched_transcript_rows = n(),
        .groups = "drop"
      ),
    by = "event_type"
  ) %>%
  mutate(
    across(
      c(tested_events, significant_events, transcript_matched_events, matched_transcript_rows),
      ~ replace_na(.x, 0L)
    )
  )

readr::write_tsv(
  summary_table,
  file.path(TABLE_DIR, "summary_by_event_type.tsv")
)

message("Done.")
message("Tables written to: ", TABLE_DIR)
message("Figures written to: ", FIGURE_DIR)
