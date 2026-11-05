suppressPackageStartupMessages({
# CPM pipeline copy of jaquino_analysis_v2/scripts/promoterswitch_isoformswitch_merge.r (CPM inputs).
  library(GenomicRanges)
  library(rtracklayer)
})

# -----------------------------
# Input files
# -----------------------------

isoform_switch_file <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/cpm_pipeline/diu/switchList.txt"

gtf_file <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/NanoporeRNASeq_extended_annotations.gtf"

switching_promoters_file <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/cpm_pipeline/promoter_usage/tables/switching_promoter_usage_9samples.tsv"

out_dir <- dirname(switching_promoters_file)

annotated_promoters_out <- file.path(out_dir, "switching_promoters_with_tss_transcripts.tsv")
combined_matches_out <- file.path(out_dir, "switching_promoters_isoform_switch_matches.tsv")
combined_all_out <- file.path(out_dir, "switching_promoters_isoform_switch_left_join.tsv")


# -----------------------------
# Parameters
# -----------------------------

# 0 = require exact transcript TSS to fall inside promoter interval.
# Use 50 or 100 if you want a small tolerance around the TSS.
TSS_FLANK_BP <- 0

REQUIRE_SAME_STRAND <- TRUE
REQUIRE_SAME_GENE <- TRUE


# -----------------------------
# Helper functions
# -----------------------------

read_tsv <- function(path) {
  read.table(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

first_present <- function(available, candidates) {
  hit <- candidates[candidates %in% available]
  if (length(hit) == 0) NA_character_ else hit[1]
}

collapse_ids <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (length(x) == 0) NA_character_ else paste(x, collapse = ";")
}

split_semicolon_ids <- function(x) {
  lapply(strsplit(as.character(x), ";", fixed = TRUE), function(v) {
    v <- trimws(v)
    v[!is.na(v) & nzchar(v)]
  })
}

parse_promoter_ids <- function(promoter_ids) {
  # Example promoter_id:
  # "10:14950417-14951517(-)"

  pat <- "^([^:]+):([0-9]+)-([0-9]+)\\(([+-])\\)$"
  regs <- regmatches(promoter_ids, regexec(pat, promoter_ids))

  parsed <- t(vapply(regs, function(x) {
    if (length(x) == 5) {
      x[2:5]
    } else {
      c(NA_character_, NA_character_, NA_character_, NA_character_)
    }
  }, character(4)))

  parsed <- as.data.frame(parsed, stringsAsFactors = FALSE)
  colnames(parsed) <- c("seqnames", "start", "end", "strand")

  parsed$start <- as.integer(parsed$start)
  parsed$end <- as.integer(parsed$end)

  bad_rows <- which(is.na(parsed$seqnames) | is.na(parsed$start) | is.na(parsed$end) | is.na(parsed$strand))

  if (length(bad_rows) > 0) {
    warning(
      "Some promoter_id values could not be parsed. Rows: ",
      paste(bad_rows, collapse = ", ")
    )
  }

  GRanges(
    seqnames = parsed$seqnames,
    ranges = IRanges(start = parsed$start, end = parsed$end),
    strand = parsed$strand
  )
}

harmonize_seqnames <- function(query_gr, reference_gr) {
  ref_seqs <- seqlevels(reference_gr)
  query_seqs <- unique(as.character(seqnames(query_gr)))
  query_seqs <- query_seqs[!is.na(query_seqs)]

  if (all(query_seqs %in% ref_seqs)) {
    return(query_gr)
  }

  # Case 1: promoter IDs are "10", GTF uses "chr10"
  if (all(paste0("chr", query_seqs) %in% ref_seqs)) {
    seqnames(query_gr) <- Rle(paste0("chr", as.character(seqnames(query_gr))))
    return(query_gr)
  }

  # Case 2: promoter IDs are "chr10", GTF uses "10"
  stripped <- sub("^chr", "", query_seqs)

  if (all(stripped %in% ref_seqs)) {
    seqnames(query_gr) <- Rle(sub("^chr", "", as.character(seqnames(query_gr))))
    return(query_gr)
  }

  warning(
    "Could not fully harmonize chromosome naming between promoter table and GTF. ",
    "Some overlaps may be missed."
  )

  query_gr
}

make_transcript_tss_granges <- function(transcripts, transcript_ids, gene_ids = NULL, flank_bp = 0) {
  tx_strand <- as.character(strand(transcripts))

  tss_pos <- ifelse(
    tx_strand == "-",
    end(transcripts),
    start(transcripts)
  )

  if (flank_bp > 0) {
    tss_start <- pmax(1L, tss_pos - flank_bp)
    tss_end <- tss_pos + flank_bp
  } else {
    tss_start <- tss_pos
    tss_end <- tss_pos
  }

  tx_tss <- GRanges(
    seqnames = seqnames(transcripts),
    ranges = IRanges(start = tss_start, end = tss_end),
    strand = strand(transcripts)
  )

  mcols(tx_tss)$transcript_id <- as.character(transcript_ids)

  if (!is.null(gene_ids)) {
    mcols(tx_tss)$gene_id_for_matching <- as.character(gene_ids)
  }

  tx_tss
}

annotate_promoters_with_transcripts <- function(
  switching_promoters,
  prom_gr,
  transcripts,
  transcript_ids,
  require_same_strand = TRUE,
  require_same_gene = TRUE,
  tss_flank_bp = 0
) {
  transcript_meta_cols <- names(mcols(transcripts))
  promoter_cols <- colnames(switching_promoters)

  common_gene_col <- first_present(
    intersect(promoter_cols, transcript_meta_cols),
    c("gene_id", "gene_name", "gene", "GeneID")
  )

  use_gene_filter <- require_same_gene && !is.na(common_gene_col)

  if (require_same_gene && !use_gene_filter) {
    warning(
      "Same-gene filtering requested, but no shared gene column was found between ",
      "switching_promoters and the GTF transcript metadata. Proceeding without gene filter."
    )
  }

  gene_ids <- NULL

  if (use_gene_filter) {
    gene_ids <- mcols(transcripts)[[common_gene_col]]
    message("Using same-gene filter with column: ", common_gene_col)
  }

  tx_tss <- make_transcript_tss_granges(
    transcripts = transcripts,
    transcript_ids = transcript_ids,
    gene_ids = gene_ids,
    flank_bp = tss_flank_bp
  )

  hits <- findOverlaps(
    prom_gr,
    tx_tss,
    ignore.strand = !require_same_strand
  )

  if (use_gene_filter && length(hits) > 0) {
    q <- queryHits(hits)
    s <- subjectHits(hits)

    promoter_gene <- as.character(switching_promoters[[common_gene_col]][q])
    transcript_gene <- as.character(mcols(tx_tss)$gene_id_for_matching[s])

    keep <- !is.na(promoter_gene) &
      !is.na(transcript_gene) &
      promoter_gene == transcript_gene

    hits <- hits[keep]
  }

  hits_list <- split(
    mcols(tx_tss)$transcript_id[subjectHits(hits)],
    queryHits(hits)
  )

  transcript_id <- vapply(seq_along(prom_gr), function(i) {
    ids <- hits_list[[as.character(i)]]
    if (is.null(ids)) return(NA_character_)
    collapse_ids(as.character(ids))
  }, character(1))

  n_linked_transcripts <- vapply(split_semicolon_ids(transcript_id), length, integer(1))

  switching_promoters$transcript_id <- transcript_id
  switching_promoters$n_linked_transcripts <- n_linked_transcripts

  switching_promoters
}

add_isoform_switch_matches <- function(switching_promoters, isoform_switch_list) {
  if (!"isoform_id" %in% colnames(isoform_switch_list)) {
    stop("isoform_switch_list must contain an isoform_id column.")
  }

  isoform_switch_list$isoform_id <- as.character(isoform_switch_list$isoform_id)
  iso_ids <- unique(isoform_switch_list$isoform_id)

  promoter_transcript_ids <- split_semicolon_ids(switching_promoters$transcript_id)

  matches_list <- lapply(promoter_transcript_ids, function(ids) {
    intersect(ids, iso_ids)
  })

  switching_promoters$isoform_match <- vapply(matches_list, function(x) {
    length(x) > 0
  }, logical(1))

  switching_promoters$isoform_match_ids <- vapply(matches_list, collapse_ids, character(1))

  list(
    switching_promoters = switching_promoters,
    matches_list = matches_list
  )
}

build_combined_table <- function(
  switching_promoters,
  isoform_switch_list,
  matches_list,
  keep_unmatched = FALSE
) {
  isoform_prefixed <- isoform_switch_list
  colnames(isoform_prefixed) <- paste0("switch_", colnames(isoform_prefixed))

  res_list <- lapply(seq_len(nrow(switching_promoters)), function(i) {
    ids <- matches_list[[i]]
    promoter_row <- switching_promoters[i, , drop = FALSE]

    if (length(ids) == 0) {
      if (!keep_unmatched) return(NULL)

      na_iso <- as.data.frame(
        matrix(NA, nrow = 1, ncol = ncol(isoform_prefixed)),
        stringsAsFactors = FALSE
      )
      colnames(na_iso) <- colnames(isoform_prefixed)

      return(cbind(promoter_row, na_iso, stringsAsFactors = FALSE))
    }

    matched <- isoform_prefixed[
      isoform_switch_list$isoform_id %in% ids,
      ,
      drop = FALSE
    ]

    if (nrow(matched) == 0) {
      return(NULL)
    }

    promoter_rows <- promoter_row[rep(1, nrow(matched)), , drop = FALSE]

    cbind(promoter_rows, matched, stringsAsFactors = FALSE)
  })

  res_list <- res_list[!vapply(res_list, is.null, logical(1))]

  if (length(res_list) == 0) {
    return(data.frame())
  }

  combined <- do.call(rbind, res_list)
  rownames(combined) <- NULL

  combined
}


# -----------------------------
# Load data
# -----------------------------

message("Loading isoform switch list...")
isoform_switch_list <- read_tsv(isoform_switch_file)

message("Loading switching promoter table...")
switching_promoters <- read_tsv(switching_promoters_file)

if (!"promoter_id" %in% colnames(switching_promoters)) {
  stop("switching_promoters table must contain a promoter_id column.")
}

message("Importing GTF...")
gtf <- import(gtf_file)

transcripts <- gtf[gtf$type %in% c("transcript", "mRNA")]

if (length(transcripts) == 0) {
  stop("No transcript or mRNA features found in the GTF.")
}


# -----------------------------
# Identify transcript ID column
# -----------------------------

transcript_id_col <- first_present(
  names(mcols(transcripts)),
  c("transcript_id", "transcriptId", "ID", "Name")
)

if (is.na(transcript_id_col)) {
  stop("Could not find a transcript ID column in the GTF metadata.")
}

message("Using transcript ID column from GTF: ", transcript_id_col)

transcript_ids <- as.character(mcols(transcripts)[[transcript_id_col]])


# -----------------------------
# Convert promoter IDs to GRanges
# -----------------------------

message("Parsing promoter_id values...")
prom_gr <- parse_promoter_ids(switching_promoters$promoter_id)

message("Harmonizing chromosome names...")
prom_gr <- harmonize_seqnames(prom_gr, transcripts)


# -----------------------------
# Annotate promoters by transcript TSS
# -----------------------------

message("Annotating promoters with transcript TSS overlaps...")

switching_promoters <- annotate_promoters_with_transcripts(
  switching_promoters = switching_promoters,
  prom_gr = prom_gr,
  transcripts = transcripts,
  transcript_ids = transcript_ids,
  require_same_strand = REQUIRE_SAME_STRAND,
  require_same_gene = REQUIRE_SAME_GENE,
  tss_flank_bp = TSS_FLANK_BP
)

message(
  "Promoter rows with >=1 linked transcript: ",
  sum(switching_promoters$n_linked_transcripts > 0, na.rm = TRUE),
  " / ",
  nrow(switching_promoters)
)

message(
  "Maximum linked transcripts for a promoter: ",
  max(switching_promoters$n_linked_transcripts, na.rm = TRUE)
)


# -----------------------------
# Match promoter-linked transcripts to switchList isoforms
# -----------------------------

message("Matching promoter-linked transcript IDs to isoform switch list...")

match_result <- add_isoform_switch_matches(
  switching_promoters = switching_promoters,
  isoform_switch_list = isoform_switch_list
)

switching_promoters <- match_result$switching_promoters
matches_list <- match_result$matches_list

message(
  "Promoter rows with >=1 isoform switch match: ",
  sum(switching_promoters$isoform_match, na.rm = TRUE)
)

message(
  "Total unique matching isoform IDs: ",
  length(unique(unlist(matches_list)))
)


# -----------------------------
# Build combined result tables
# -----------------------------

combined_matches <- build_combined_table(
  switching_promoters = switching_promoters,
  isoform_switch_list = isoform_switch_list,
  matches_list = matches_list,
  keep_unmatched = FALSE
)

combined_all <- build_combined_table(
  switching_promoters = switching_promoters,
  isoform_switch_list = isoform_switch_list,
  matches_list = matches_list,
  keep_unmatched = TRUE
)


# -----------------------------
# Write outputs
# -----------------------------

message("Writing outputs...")

write.table(
  switching_promoters,
  annotated_promoters_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  combined_matches,
  combined_matches_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  combined_all,
  combined_all_out,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Done.")
message("Annotated promoters: ", annotated_promoters_out)
message("Promoter/isoform matches only: ", combined_matches_out)
message("Promoter/isoform left join: ", combined_all_out)