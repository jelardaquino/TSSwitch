library(GenomicFeatures)
library(GenomicRanges)
library(SplicingGraphs)
library(dplyr)
library(tidyr)

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("No geneID provided")
}

geneID <- args[1]

cat("Running analysis for gene:", geneID, "\n")

output_dir_figures <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis_v2/figures/splicing_graph/"
dir.create(output_dir_figures, showWarnings = FALSE, recursive = TRUE)

load("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn25671149/jaquino_analysis/tmp/bambu.RData")
txdb <- loadDb("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn25671149/jaquino_analysis/tmp/splicing_graph/Homo_sapiens.GRCh38.107_ERC_txdb.sqlite")


# GRangesList: one element per gene
tx_by_gene <- transcriptsBy(txdb, "gene")

# Transcripts just for your gene
gene_tx_gr <- tx_by_gene[[geneID]]

gene_tx <- data.frame(
  tx_id   = mcols(gene_tx_gr)$tx_id,
  tx_name = mcols(gene_tx_gr)$tx_name,
  seqnames = as.character(seqnames(gene_tx_gr)),
  start    = start(gene_tx_gr),
  end      = end(gene_tx_gr),
  strand   = as.character(strand(gene_tx_gr))
)

head(gene_tx)

tx_map <- gene_tx[, c("tx_id", "tx_name", "strand")]  # numeric → ENST mapping (+ strand)
tx_map$tx_id <- as.integer(tx_map$tx_id)

head(tx_map)

sample_info_file <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/tables/promoter_usage_61_genes/sample_info.txt"
sample_info <- read.table(sample_info_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Prefix an "X" to each sampleID (but avoid double-prefixing if already present)
sample_info$sampleID <- as.character(sample_info$sampleID)
sample_info$sampleID <- ifelse(grepl("^X", sample_info$sampleID), sample_info$sampleID, paste0("X", sample_info$sampleID))
## tx_counts: rows = transcript IDs, cols = samples
## Make sure rownames(tx_counts) are transcript IDs that match tx_list
# e.g. rownames(tx_counts) <- c("ENST0000...", ...) 
#tx_counts <- read.table("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn25671149/jaquino_analysis/tmp/NanoporeRNASeq_counts_transcript.txt", header = TRUE, row.names = 1)
tx_counts <- read.table("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/pychopperfq/minimap2/bam_mapq10/all_bambu_Rscript/bambu_tx_tpm_GENEID_genesymbol.txt", header = TRUE, row.names = 1)
# txID_geneID <- data.frame(tx_id = rownames(tx_counts), gene_id = tx_counts$GENEID)

# tx_counts$GENEID <- NULL
# colnames(tx_counts) <- gsub("_pychopper.fq.gz.filtered", "", colnames(tx_counts))
tx_counts <- tx_counts[, sample_info$sampleID]
# keep only rows whose rownames start with ENST
if (!is.null(rownames(tx_counts))) {
    tx_counts <- tx_counts[grepl("^ENST", rownames(tx_counts)), , drop = FALSE]
}

# Keep only transcripts that are both in tx_counts and this gene
common_tx <- intersect(rownames(tx_counts), tx_map$tx_name)
tx_counts_sub <- tx_counts[common_tx, , drop = FALSE]

# Remove low-expression transcripts by total abundance across samples
tx_counts_sub <- tx_counts_sub[rowMeans(tx_counts_sub, na.rm = TRUE) > 1, , drop = FALSE]
if (nrow(tx_counts_sub) == 0) {
  stop("No transcripts with rowMeans > 1 for gene: ", geneID)
}

# Long format with conditions
expr_long <- tx_counts_sub %>%
  as.data.frame() %>%
  tibble::rownames_to_column("tx_name") %>%
  pivot_longer(
    cols = -tx_name,
    names_to = "sample",
    values_to = "expr"
  ) %>%
  mutate(condition = sample_info$condition[match(sample, sample_info$sampleID)])

# Attach numeric tx_id from tx_map
expr_long <- expr_long %>%
  left_join(tx_map, by = "tx_name")

# Summarise mean expression per transcript & condition
expr_summary <- expr_long %>%
  group_by(tx_id, tx_name, condition) %>%
  summarise(mean_expr = mean(expr), .groups = "drop") %>%
  pivot_wider(
    names_from  = condition,
    values_from = mean_expr,
    values_fill = 0
  )

expr_summary

# Filter transcripts by expression
expr_cutoff <- 0       # tweak as you like

cond_cols <- setdiff(colnames(expr_summary), c("tx_id", "tx_name"))
expr_summary$max_expr <- apply(expr_summary[, cond_cols, drop = FALSE], 1, max)

expr_summary_filt <- expr_summary %>%
  filter(max_expr >= expr_cutoff)

expr_summary_filt

sg <- SplicingGraphs(txdb)

# All edges for the gene (exons + introns + R/L)
gene_edges <- as.data.frame(sgedges(sg[geneID]))

# Inspect columns once:
# names(gene_edges)

# Keep only exons
exon_edges <- gene_edges %>%
  filter(ex_or_in == "ex")

exon_edges_long <- exon_edges %>%
  tidyr::unnest_longer(tx_id)   # or `unnest(tx_id)` depending on your tidyr version

head(exon_edges_long)

exon_edges_long <- exon_edges_long %>%
  mutate(
    from = as.numeric(as.character(from)),
    to   = as.numeric(as.character(to))
  )

node_labels <- sgnodes(sg[geneID])  # c("R","1","2",...,"11","L")
n_nodes <- length(node_labels)
x_nodes <- seq_len(n_nodes)

gene_edge_gr <- sgedgesByGene(sg[geneID])[[1]] 

node_pos_df <- data.frame(
  node = c(mcols(gene_edge_gr)$from, mcols(gene_edge_gr)$to),
  pos  = c(start(gene_edge_gr),      end(gene_edge_gr)),
  stringsAsFactors = FALSE
) %>%
  group_by(node) %>%
  summarise(pos = median(pos), .groups = "drop")

# Keep only nodes that exist in the plot ordering
node_pos_df <- node_pos_df %>%
  filter(node %in% node_labels)

# Transcript order & labels
tx_plot <- expr_summary_filt
tx_plot_ids  <- tx_plot$tx_name
n_tx         <- length(tx_plot_ids)

# Node x positions where plotted transcripts start (used to align promoters)
gene_strand <- unique(gene_tx$strand)[1]
tx_start_nodes <- exon_edges_long %>%
  filter(tx_id %in% tx_plot_ids) %>%
  group_by(tx_id) %>%
  summarise(
    start_node = if (gene_strand == "-") max(pmax(from, to), na.rm = TRUE) else min(pmin(from, to), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(start_x = start_node + 1)

tx_start_x_candidates <- sort(unique(tx_start_nodes$start_x))

# One base y per transcript
tx_base_y <- seq_len(n_tx)

# Conditions and colors
cond_levels <- cond_cols
n_cond      <- length(cond_levels)
cond_colors <- setNames(
  c("salmon", "turquoise")[seq_len(n_cond)],
  cond_levels
)

# Slight y-offset for each condition so CT/AD don't sit exactly on top of each other
cond_offsets <- seq(from = -0.15, to = 0.15, length.out = n_cond)

# Helper: y position for transcript index t_idx and condition index c_idx
tx_y_for_cond <- function(t_idx, c_idx) {
  tx_base_y[t_idx] + cond_offsets[c_idx]
}

# For scaling line widths per condition
max_expr_per_cond <- sapply(
  cond_levels,
  function(cc) max(tx_plot[[cc]], na.rm = TRUE)
)

# ---- Get gene edges as GRanges with from/to and genomic ranges ----
gene_edge_gr <- sgedgesByGene(sg[geneID])[[1]]  # GRanges for this gene

# Build node -> genomic coordinate lookup from edge boundaries
node_pos_df <- data.frame(
  node = c(mcols(gene_edge_gr)$from, mcols(gene_edge_gr)$to),
  pos  = c(start(gene_edge_gr), end(gene_edge_gr)),
  stringsAsFactors = FALSE
) %>%
  mutate(node_num = suppressWarnings(as.integer(node))) %>%
  filter(!is.na(node_num)) %>%
  group_by(node_num) %>%
  summarise(pos = median(pos), .groups = "drop") %>%
  arrange(node_num)

# Node coordinate vector indexed by node number
node_nums   <- node_pos_df$node_num
node_coords <- node_pos_df$pos



# ---- Junction parsing helpers (regtools BED12-ish) ----
read_regtools_junctions <- function(bed_file, sample_id = NA, condition = NA) {
  j <- read.table(bed_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "")
  if (ncol(j) < 12) stop("Expected 12 columns in regtools junction BED: ", bed_file)

  colnames(j)[1:12] <- c(
    "chrom","chromStart","chromEnd","name","count","strand",
    "thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts"
  )

  # Parse blockSizes "225,70" -> c(225, 70)
  parse_sizes <- function(x) as.integer(strsplit(gsub(",+$","", x), ",")[[1]])

  sizes <- lapply(j$blockSizes, parse_sizes)
  left_size  <- vapply(sizes, `[[`, integer(1), 1)
  right_size <- vapply(sizes, `[[`, integer(1), 2)

  # Exact splice junction (intron) boundaries:
  # donor   = chromStart + left_block_size
  # acceptor= chromEnd   - right_block_size
  j$donor    <- j$chromStart + left_size
  j$acceptor <- j$chromEnd   - right_size

  j$junction_id <- paste0(j$chrom, ":", j$donor, ":", j$acceptor)
  j$sample <- sample_id
  j$condition <- condition

  j
}
# ---- Load junctions from regtools BEDs ----
junction_dir <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/regtools"

junction_list <- lapply(sample_info$sampleID, function(s) {
  bed <- file.path(junction_dir, paste0(s, "_junctions.bed"))
  if (!file.exists(bed)) return(NULL)
  read_regtools_junctions(
    bed_file = bed,
    sample_id = s,
    condition = sample_info$condition[match(s, sample_info$sampleID)]
  )
})
junction_df <- bind_rows(junction_list)

# If nothing loaded, keep going without junctions
if (nrow(junction_df) == 0) {
  warning("No junction BEDs loaded; skipping junction overlay.")
}

# Keep junctions on same chrom and overlapping gene bounds
gene_chr   <- as.character(unique(seqnames(gene_edge_gr)))[1]
gene_start <- min(start(gene_edge_gr))
gene_end   <- max(end(gene_edge_gr))

junction_df_gene <- junction_df %>%
  filter(
    chrom == gene_chr,
    donor >= gene_start,
    acceptor <= gene_end
  )

# Map genomic coordinate -> nearest node number
map_coord_to_node_label <- function(coord, node_pos_df) {
  as.character(node_pos_df$node_num[which.min(abs(node_pos_df$pos - coord))])
}

if (exists("junction_df_gene") && nrow(junction_df_gene) > 0) {
  junction_df_gene$from_node_label <- vapply(
    junction_df_gene$donor,
    map_coord_to_node_label,
    character(1),
    node_pos_df = node_pos_df
  )

  junction_df_gene$to_node_label <- vapply(
    junction_df_gene$acceptor,
    map_coord_to_node_label,
    character(1),
    node_pos_df = node_pos_df
  )

  # Convert node labels -> plot x positions
  junction_df_gene$from_x <- match(junction_df_gene$from_node_label, node_labels)
  junction_df_gene$to_x   <- match(junction_df_gene$to_node_label,   node_labels)

  # Keep only junctions that match intron edges from transcripts that survived filtering
  intron_edges_long <- gene_edges %>%
    filter(ex_or_in == "in") %>%
    tidyr::unnest_longer(tx_id) %>%
    mutate(
      from_x = match(from, node_labels),
      to_x   = match(to,   node_labels),
      jkey   = paste(pmin(from_x, to_x), pmax(from_x, to_x), sep = "__")
    ) %>%
    filter(
      tx_id %in% tx_plot_ids,
      !is.na(from_x), !is.na(to_x),
      from_x != to_x
    ) %>%
    distinct(jkey)

  allowed_jkeys <- intron_edges_long$jkey

  junction_sum <- junction_df_gene %>%
    filter(!is.na(from_x), !is.na(to_x), from_x != to_x) %>%
    mutate(jkey = paste(pmin(from_x, to_x), pmax(from_x, to_x), sep = "__")) %>%
    filter(jkey %in% allowed_jkeys) %>%
    group_by(junction_id, chrom, donor, acceptor, from_x, to_x, condition, jkey) %>%
    summarise(junc_count = sum(count), .groups = "drop")
} else {
  junction_sum <- NULL
}

# (optional) drop tiny counts
junc_cutoff <- 2

if (!is.null(junction_sum) && nrow(junction_sum) > 0) {
  junction_sum <- junction_sum %>%
    mutate(condition = as.character(condition)) %>%
    filter(junc_count >= junc_cutoff)

  missing_conds <- setdiff(unique(junction_sum$condition), names(cond_colors))
  if (length(missing_conds) > 0) {
    warning("junction_sum has conditions not in cond_colors: ", paste(missing_conds, collapse=","))
  }

  if (nrow(junction_sum) == 0) {
    junction_sum <- NULL
  }
}

# ==== PROMOTERS: filter to gene, summarize by condition, map to plot x ====
prom_anno <- read.table("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis_v2/tables/promoter_usage_86_genes/promoter_annotations.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
promoter_counts_dt <- read.table("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis_v2/tables/promoter_usage_86_genes/promoter_counts_9samples.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

prom_gene_ids <- prom_anno %>%
  filter(gene_id == geneID) %>%
  distinct(promoter_id)

# If no promoters, we'll just skip plotting them
has_promoters <- nrow(prom_gene_ids) > 0

prom_sum <- NULL

if (has_promoters) {
  # 2) subset counts table to promoters for this gene
  prom_dt <- promoter_counts_dt %>%
    as.data.frame() %>%
    dplyr::inner_join(prom_gene_ids, by = "promoter_id")

  # 3) gather sample columns -> long
  # (assumes sample_info$sampleID matches your promoter count column names)
  sample_cols <- intersect(sample_info$sampleID, colnames(prom_dt))

  prom_long <- prom_dt %>%
    tidyr::pivot_longer(
      cols = all_of(sample_cols),
      names_to = "sample",
      values_to = "count"
    ) %>%
    dplyr::mutate(condition = sample_info$condition[match(sample, sample_info$sampleID)])

  # 4) summarize mean promoter usage per condition
  prom_sum <- prom_long %>%
    dplyr::group_by(promoter_id, seqnames, start, end, strand, condition) %>%
    dplyr::summarise(mean_count = mean(count, na.rm = TRUE), .groups = "drop")

  # 5) keep promoters on same chrom and near the gene region (allow some upstream padding)
  upstream_pad <- 5000
  prom_sum <- prom_sum %>%
    dplyr::filter(
      as.character(seqnames) == gene_chr,
      start <= (gene_end + 1000),
      end   >= (gene_start - upstream_pad)
    )

  # 6) map promoter midpoints to plot x via interpolation using node genomic positions
  #    node_pos_df has: node_num, pos ; x_nodes are 1..n_nodes
  coord_to_x_interp <- function(coord, node_pos_df, x_nodes) {
    ok <- !is.na(node_pos_df$pos) & !is.na(node_pos_df$node_num)
    xs <- x_nodes[node_pos_df$node_num[ok]]  # node_num is 1..(n_nodes-2) typically
    ps <- node_pos_df$pos[ok]
    # ensure sorted by genomic pos
    o <- order(ps)
    ps <- ps[o]; xs <- xs[o]
    # interpolate (rule=2 clamps outside range)
    as.numeric(approx(x = ps, y = xs, xout = coord, rule = 2)$y)
  }

  prom_sum <- prom_sum %>%
    dplyr::mutate(
      mid   = (start + end) / 2,
      x_raw = coord_to_x_interp(mid, node_pos_df, x_nodes),
      # snap promoter x to nearest transcript start node
      x = if (length(tx_start_x_candidates) > 0) {
        vapply(
          x_raw,
          function(xx) tx_start_x_candidates[which.min(abs(tx_start_x_candidates - xx))],
          numeric(1)
        )
      } else {
        x_raw
      }
    )
}


png(paste0(output_dir_figures, "/", geneID, "_splicing_graph.png"),
    width = 2400, height = 1600, res = 300)
par(mar = c(6, 8, 4, 12), bty = "n")

arc_room <- 1.0
# plot(
#   0, 0,
#   type = "n", xaxt = "n", yaxt = "n",
#   xlab = geneID, ylab = "",
#   xlim = c(1, n_nodes),
#   ylim = c(0, n_tx + 1)
# )

plot(
  0, 0,
  type = "n", xaxt = "n", yaxt = "n",
  xlab = geneID, ylab = "",
  xlim = c(1, n_nodes),
  ylim = c(-arc_room, n_tx + 1.6)
)
# Draw splice graph nodes at the bottom
points(
  x = x_nodes,
  y = rep(0.5, n_nodes),
  pch = 16,
  col = "orange",
  cex = 1.2
)

# Reduce node label density for large graphs to avoid overlap
node_plot_labels <- rep("", n_nodes)
node_plot_labels[1] <- "R"
node_plot_labels[n_nodes] <- "L"

numeric_node_vals <- seq_len(n_nodes - 2)
if (length(numeric_node_vals) > 0) {
  if (n_nodes > 10) {
    show_idx <- which(numeric_node_vals %% 5 == 0)
  } else {
    show_idx <- seq_along(numeric_node_vals)
  }
  node_plot_labels[show_idx + 1] <- as.character(numeric_node_vals[show_idx])
}

text(
  x = x_nodes,
  y = 0.6,
  labels = node_plot_labels,
  cex = 0.6
)

# --- Exon segments on the bottom node track (connect nodes) ---
edge_df <- as.data.frame(gene_edge_gr)

# Keep only exons
exon_edge_df <- edge_df %>% filter(ex_or_in == "ex")

# Map node labels -> x positions in your plot
exon_edge_df$from_x <- match(exon_edge_df$from, node_labels)
exon_edge_df$to_x   <- match(exon_edge_df$to,   node_labels)

# Draw as horizontal segments on the same y as the nodes
y_nodes <- 0.5

# Optional: small vertical offset so the segment doesn't sit exactly on top of the points
y_exon <- y_nodes

# Draw
ok <- !is.na(exon_edge_df$from_x) & !is.na(exon_edge_df$to_x) & exon_edge_df$from_x != exon_edge_df$to_x
apply(exon_edge_df[ok, ], 1, function(r) {
  segments(
    x0 = as.numeric(r["from_x"]),
    y0 = y_exon,
    x1 = as.numeric(r["to_x"]),
    y1 = y_exon,
    col = "purple",   # pick what you like
    lwd = 5,
    lend = "butt"
  )
})

# ==== DRAW PROMOTERS (lane above transcripts) ====
if (!is.null(prom_sum) && nrow(prom_sum) > 0) {

  # lane location
  y_prom_base <- n_tx + 0.8

  # reuse condition offsets but smaller so points don't overlap
  prom_offsets <- seq(from = -0.12, to = 0.12, length.out = length(cond_levels))
  names(prom_offsets) <- cond_levels

  # size scaling
  prom_sum$w <- log10(prom_sum$mean_count + 1)
  wmin <- min(prom_sum$w, na.rm = TRUE)
  wmax <- max(prom_sum$w, na.rm = TRUE)

  # draw a thin baseline to show it's a separate track (optional)
  segments(1, y_prom_base, n_nodes, y_prom_base, col = "gray85", lwd = 1)

  for (i in seq_len(nrow(prom_sum))) {
    r <- prom_sum[i, ]
    cond_key <- as.character(r$condition)
    if (!cond_key %in% names(cond_colors)) next

    # size: 0.6..1.8-ish
    cex <- 0.6 + 1.2 * ((r$w - wmin) / (wmax - wmin + 1e-6))

    # y with condition offset
    y <- y_prom_base + prom_offsets[[cond_key]]

    # promoter marker (triangle)
    col <- grDevices::adjustcolor(cond_colors[[cond_key]], alpha.f = 0.85)

    points(r$x, y, pch = 17, col = col, cex = cex)

    # optional: small stem down toward node track (visual association)
    segments(r$x, y - 0.05, r$x, 0.65, col = grDevices::adjustcolor(col, alpha.f = 0.35), lwd = 1, lty = 2)
  }

  # label the lane
  text(x = 1, y = y_prom_base + 0.75, labels = "Promoters (mean by condition)", adj = 0, cex = 0.7)
}

# Transcript labels on the left (ENST IDs)
axis(
  side = 2,
  at   = tx_base_y,
  labels = tx_plot_ids,
  las  = 2,
  tick = FALSE,
  cex.axis = 0.6
)

# Precompute a log-scaled global range
expr_log <- log10(tx_plot[, cond_levels, drop = FALSE] + 1)
global_min_log <- min(expr_log, na.rm = TRUE)
global_max_log <- max(expr_log, na.rm = TRUE)

# Draw exon segments per transcript & condition
for (t_idx in seq_len(n_tx)) {
  this_tx <- tx_plot_ids[t_idx]

  # All edges for this transcript
  tx_edges <- exon_edges_long %>% filter(tx_id == this_tx)
  if (nrow(tx_edges) == 0) next

  for (c_idx in seq_along(cond_levels)) {
    cond_name <- cond_levels[c_idx]
    expr_val  <- tx_plot[[cond_name]][t_idx]

    if (expr_val <= 0) next

    # ---- line width & y ----
    expr_val_log <- log10(expr_val + 1)
    lwd <- 0.5 + 6 * ( (expr_val_log - global_min_log) /
                   (global_max_log - global_min_log + 1e-6) )
    y   <- tx_y_for_cond(t_idx, c_idx)

    # ---- draw exon segments ----
    segments(
      x0 = tx_edges$from + 1,
      y0 = y,
      x1 = tx_edges$to + 1,
      y1 = y,
      col = cond_colors[cond_name],
      lwd = lwd
    )

    # ==== NEW: expression label for this transcript & condition ====
    # choose x position: center of the transcript's exonic span
    # x_label <- (min(tx_edges$from + 1) + max(tx_edges$to + 1)) / 2
    # or left of first exon:
    x_label <- min(tx_edges$from + 1) - 0.5

    text(
      x = x_label,
      y = y + 0.12,                            # small offset above the line
      labels = round(expr_val, 2),            # format as you like
      col = cond_colors[cond_name],
      cex = 0.6
    )
    # ===============================================================
  }
}

# ---- Draw junction arcs at the bottom (by condition) ----
if (!is.null(junction_sum) && nrow(junction_sum) > 0) {

  # scaling for arc line width
  junction_sum$w <- log10(junction_sum$junc_count + 1)
  w_min <- min(junction_sum$w, na.rm = TRUE)
  w_max <- max(junction_sum$w, na.rm = TRUE)

  # arc drawing helper
  draw_arc <- function(x1, x2, y = 0.5, height = 0.35, mid_offset = 0,
                     col = "gray40", lwd = 1, shape = 0.6) {

  xm <- (x1 + x2) / 2
  xs <- c(x1, xm, x2)
  ys <- c(y, y - height + mid_offset, y)   # <-- flipped

  xspline(xs, ys, shape = shape, open = TRUE, border = col, col = NA, lwd = lwd)
}
#   draw_arc <- function(x1, x2, y = 0.5, height = 0.35, mid_offset = 0,
#                      col = "gray40", lwd = 1, shape = 0.6) {

#   xm <- (x1 + x2) / 2
#   xs <- c(x1, xm, x2)
#   ys <- c(y, y + height + mid_offset, y)

#   xspline(
#     xs, ys,
#     shape = shape,
#     open  = TRUE,
#     border = col,   # <-- THIS is the stroke color
#     col    = NA,    # no fill
#     lwd    = lwd
#   )
# }


  # small vertical offsets so condition arcs don’t fully overlap
  # arc_offsets <- setNames(seq(from = -0.06, to = 0.06, length.out = length(cond_levels)), cond_levels)
arc_offsets <- setNames(seq(from = -0.10, to = -0.02, length.out = length(cond_levels)), cond_levels)

  for (i in seq_len(nrow(junction_sum))) {
    row <- junction_sum[i, ]

    # skip weird mappings
    if (row$from_x < 1 || row$to_x < 1 || row$from_x > n_nodes || row$to_x > n_nodes) next
    if (row$from_x == row$to_x) next

    # scale lwd
    lwd <- 0.5 + 6 * ((row$w - w_min) / (w_max - w_min + 1e-6))

    # arc height: proportional to distance (keeps short introns from being huge)
    dist <- abs(row$to_x - row$from_x)
    height <- min(1.0, 0.15 + 0.05 * dist)

    # condition color (with a bit of transparency)
    cond_key <- as.character(row$condition)
    base_col <- cond_colors[[cond_key]]
    if (is.null(base_col) || is.na(base_col)) base_col <- "gray40"
    col <- grDevices::adjustcolor(base_col, alpha.f = 0.6)

    mid_offset <- arc_offsets[[cond_key]]
    if (is.null(mid_offset) || is.na(mid_offset)) mid_offset <- 0

    draw_arc(row$from_x, row$to_x, y = 0.5, height = height, mid_offset = mid_offset,
         col = col, lwd = lwd)
  }

  # optional legend note
  mtext(paste0("Junction arcs: regtools counts (cutoff ≥ ", junc_cutoff, ")"), side = 1, line = 4, cex = 0.7)
}

# Legend for conditions
legend(
  "top",
  inset  = 0.02,
  legend = cond_levels,
  col    = cond_colors[cond_levels],
  lwd    = 3,
  horiz  = TRUE,
  bty    = "n",
  cex    = 0.8
)

# Right-side feature legend (shifted farther right)
feature_labels <- c(
  "Promoters",
  "Splicing graph exons",
  "Splice junctions",
  "Splice junction counts",
  "AD Transcript exons",
  "CT Transcript exons"
)

feature_pch <- c(17, NA, 16, NA, NA, NA)
feature_lty <- c(NA, 1, NA, NA, 1, 1)
feature_lwd <- c(NA, 5, NA, NA, 3, 3)
feature_col <- c("black", "purple", "orange", "gray40", "salmon", "turquoise")
feature_pt_cex <- c(1.1, NA, 1.1, NA, NA, NA)

usr <- par("usr")
legend_x <- usr[2] + 0.35
legend_y <- usr[4]

legend_info <- legend(
  x = legend_x,
  y = legend_y,
  xjust = 0,
  yjust = 1,
  legend = feature_labels,
  pch = feature_pch,
  lty = feature_lty,
  lwd = feature_lwd,
  col = feature_col,
  pt.cex = feature_pt_cex,
  bty = "n",
  cex = 0.72,
  y.intersp = 1.1,
  seg.len = 2.2,
  xpd = NA,
  plot = FALSE
)

legend(
  x = legend_x,
  y = legend_y,
  xjust = 0,
  yjust = 1,
  legend = feature_labels,
  pch = feature_pch,
  lty = feature_lty,
  lwd = feature_lwd,
  col = feature_col,
  pt.cex = feature_pt_cex,
  bty = "n",
  cex = 0.72,
  y.intersp = 1.1,
  seg.len = 2.2,
  xpd = NA
)

# Draw a concave-up arc glyph for "Splice junction counts" in the legend
idx_junc_counts <- match("Splice junction counts", feature_labels)
if (!is.na(idx_junc_counts)) {
  y_arc <- legend_info$text$y[idx_junc_counts]
  x1 <- legend_info$rect$left + 0.03 * legend_info$rect$w
  x2 <- x1 + 0.11 * legend_info$rect$w
  xm <- (x1 + x2) / 2
  h  <- 0.02 * diff(par("usr")[3:4])

  xspline(
    x = c(x1, xm, x2),
    y = c(y_arc, y_arc - h, y_arc),
    shape = 0.6,
    open = TRUE,
    border = "gray40",
    col = NA,
    lwd = 2,
    xpd = NA
  )
}

dev.off()
