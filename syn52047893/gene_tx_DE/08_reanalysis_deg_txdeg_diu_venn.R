###############################
# Area-proportional 3-set Venn of re-analysis (PMI_RIN) gene DEGs, genes with >= 1 transcript DEG,
# and the original 86 DIU genes (6 AD vs 6 CT). Style follows 06_deg_diu_venn.R; sets follow
# re-analysis_PMI_RIN/04b_venn_JA_DEG.R.
#   Gene DEG: gene_DEG/gene_DESeq2_<model>_sig_padj0.05_absLFC1.tsv (padj < 0.05, |shrunk log2FC| > 1)
#   Tx DEG:   tx_DEG/tx_DESeq2_<model>_sig_padj0.05_absLFC1.tsv, collapsed to parent gene_id
#   DIU:      86_isoformswitch_JA.csv (IsoformSwitchAnalyzeR, FDR < 0.05, |dIF| > 0.1)
# Circle areas match set sizes and each pairwise lens matches the pairwise overlap;
# the triple region is not fitted. Gene IDs compared with the Ensembl version stripped.
# Usage: Rscript 08_reanalysis_deg_txdeg_diu_venn.R [model]     (default noCov)
###############################

suppressPackageStartupMessages({
  library(ggplot2)
})

# When sourced (09_supp_figure_deg_txdeg_diu.R), `model` and `save_plots` are set by the caller
args <- commandArgs(trailingOnly = TRUE)
if (!exists("model")) model <- if (length(args) >= 1) args[1] else "noCov"
if (!exists("save_plots")) save_plots <- TRUE
stopifnot(model %in% c("noCov", "withCov"))

S <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing"
PY <- file.path(S, "syn52065646_LongRead/pychopperfq")
RE <- file.path(PY, "re-analysis_PMI_RIN")
OUT <- file.path(S, "ERP161086/emtab14167_replication")
FIGS <- file.path(OUT, "figures")

GENE_SIG <- file.path(RE, "gene_DEG", sprintf("gene_DESeq2_%s_sig_padj0.05_absLFC1.tsv", model))
TX_SIG <- file.path(RE, "tx_DEG", sprintf("tx_DESeq2_%s_sig_padj0.05_absLFC1.tsv", model))
DIU_JA <- file.path(PY, "minimap2/bam_mapq10/all_bambu_Rscript/all_bambu_python_analysis",
                    "IsoformSwitchAnalyzeR/ISO_results_v3/86_isoformswitch_JA.csv")

COLS <- c(DEG = "#3A6FB0", TX = "#1BAF7A", DIU = "#D9822B")

strip_ver <- function(x) sub("\\.[0-9]+$", "", x)

###############################
# Gene sets
###############################

g <- read.delim(GENE_SIG, stringsAsFactors = FALSE)
t <- read.delim(TX_SIG, stringsAsFactors = FALSE)
diu <- read.csv(DIU_JA, stringsAsFactors = FALSE)
diu$gene_id <- strip_ver(diu$gene_id)
diu <- unique(diu[, c("gene_id", "gene_name")])
stopifnot(nrow(diu) == 86)

sets <- list(DEG = unique(strip_ver(g$gene_id)),
             TX = unique(strip_ver(t$gene_id[!is.na(t$gene_id) & t$gene_id != ""])),
             DIU = diu$gene_id)
sym <- rbind(data.frame(gene_id = diu$gene_id, gene_name = diu$gene_name),
             data.frame(gene_id = strip_ver(g$gene_id), gene_name = g$gene_symbol),
             data.frame(gene_id = strip_ver(t$gene_id), gene_name = t$gene_symbol))
sym <- sym[!is.na(sym$gene_name) & sym$gene_name != "" & !duplicated(sym$gene_id), ]

# Exclusive regions, keyed by membership in DEG / TX / DIU
all_ids <- unique(unlist(sets))
mem <- sapply(sets, function(s) all_ids %in% s)
key <- apply(mem, 1, function(r) paste(as.integer(r), collapse = ""))
region <- split(all_ids, key)
n_reg <- function(k) length(region[[k]])
n <- lengths(sets)
cat("model:", model, "| gene DEG:", n[["DEG"]], "| tx-DEG genes:", n[["TX"]], "| DIU:", n[["DIU"]], "\n")
for (k in c("100", "010", "001", "110", "101", "011", "111")) cat("  region", k, ":", n_reg(k), "\n")

###############################
# Circle geometry: areas proportional to set sizes, each pairwise lens area = pairwise overlap
###############################

lens_area <- function(d, r1, r2) {
  if (d >= r1 + r2) return(0)
  if (d <= abs(r1 - r2)) return(pi * min(r1, r2)^2)
  a1 <- r1^2 * acos((d^2 + r1^2 - r2^2) / (2 * d * r1))
  a2 <- r2^2 * acos((d^2 + r2^2 - r1^2) / (2 * d * r2))
  a3 <- 0.5 * sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))
  a1 + a2 - a3
}
r <- sqrt(n / pi)
k <- r[["DEG"]] / sqrt(1568 / pi)  # scale label offsets to the 06_deg_diu_venn.R layout
pair_d <- function(a, b) {
  ov <- length(intersect(sets[[a]], sets[[b]]))
  if (ov == 0) return(r[[a]] + r[[b]] + 1 * k)
  uniroot(function(x) lens_area(x, r[[a]], r[[b]]) - ov,
          c(abs(r[[a]] - r[[b]]) + 1e-6, r[[a]] + r[[b]] - 1e-6))$root
}
d_dd <- pair_d("DEG", "DIU"); d_dt <- pair_d("DEG", "TX"); d_td <- pair_d("TX", "DIU")

# DEG at the origin, DIU on the x axis, TX triangulated below
ctr <- data.frame(set = c("DEG", "DIU", "TX"), x = c(0, d_dd, NA), y = c(0, 0, NA))
tx_x <- (d_dt^2 - d_td^2 + d_dd^2) / (2 * d_dd)
ctr[3, c("x", "y")] <- c(tx_x, -sqrt(max(d_dt^2 - tx_x^2, 0)))
rownames(ctr) <- ctr$set

circle <- function(set) {
  th <- seq(0, 2 * pi, length.out = 400)
  data.frame(x = ctr[set, "x"] + r[[set]] * cos(th), y = ctr[set, "y"] + r[[set]] * sin(th), set = set)
}
circles <- do.call(rbind, lapply(c("DEG", "TX", "DIU"), circle))

# Region anchor = mean of grid points falling in that exclusive region
xr <- range(circles$x); yr <- range(circles$y)
grid <- expand.grid(x = seq(xr[1], xr[2], length.out = 600), y = seq(yr[1], yr[2], length.out = 600))
inside <- sapply(c("DEG", "TX", "DIU"), function(s) (grid$x - ctr[s, "x"])^2 + (grid$y - ctr[s, "y"])^2 <= r[[s]]^2)
grid$key <- apply(inside, 1, function(z) paste(as.integer(z), collapse = ""))
anchor <- function(kk, target = NULL) {
  p <- grid[grid$key == kk, ]
  if (nrow(p) == 0) return(NULL)
  if (is.null(target)) target <- c(mean(p$x), mean(p$y))
  i <- which.min((p$x - target[1])^2 + (p$y - target[2])^2)  # nearest grid point to the target stays inside the region
  c(p$x[i], p$y[i])
}
# DIU-only callout points at the upper right of the DIU circle so its leader line clears the overlaps
anchor_call <- function(kk) {
  if (kk == "001") anchor(kk, c(ctr["DIU", "x"] + 0.35 * r[["DIU"]], ctr["DIU", "y"] + 0.55 * r[["DIU"]])) else anchor(kk)
}

###############################
# Labels: large exclusive regions inside their circle, overlaps as callouts on the right
###############################

names_of <- function(kk) {
  nm <- sym$gene_name[match(region[[kk]], sym$gene_id)]
  nm[is.na(nm)] <- region[[kk]][is.na(nm)]
  sort(nm)
}
callouts <- list(
  list(key = "001", txt = sprintf("%d of %d DIU genes (%.0f%%)\nwith no gene- or\ntranscript-level change",
                                  n_reg("001"), n[["DIU"]], 100 * n_reg("001") / n[["DIU"]])),
  list(key = "011", txt = "DIU and transcript DEG,\nno gene-level change"),
  list(key = "101", txt = "DIU and gene DEG"),
  list(key = "111", txt = "in all three"))
label_x <- max(circles$x) + 4 * k
call_y <- seq(r[["DEG"]] * 0.75, -r[["DEG"]] * 0.75, length.out = length(callouts))
# Slot order for callouts with a leader line: fewest crossings, then shortest total length; empty regions last
crosses <- function(p1, p2, q1, q2) {
  o <- function(a, b, c) sign((b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]))
  o(p1, p2, q1) != o(p1, p2, q2) && o(q1, q2, p1) != o(q1, q2, p2)
}
perms <- function(v) if (length(v) <= 1) list(v) else
  do.call(c, lapply(seq_along(v), function(i) lapply(perms(v[-i]), function(p) c(v[i], p))))
has_line <- sapply(callouts, function(cc) n_reg(cc$key) > 0 && !is.null(anchor_call(cc$key)))
lined <- which(has_line)
best <- NULL; best_score <- Inf
for (pp in perms(lined)) {
  ends <- lapply(seq_along(pp), function(s) list(a = anchor_call(callouts[[pp[s]]]$key),
                                                 b = c(label_x - 0.5 * k, call_y[s])))
  nx <- 0
  for (i in seq_along(ends)) for (j in seq_along(ends)) if (i < j)
    nx <- nx + crosses(ends[[i]]$a, ends[[i]]$b, ends[[j]]$a, ends[[j]]$b)
  len <- sum(sapply(ends, function(e) sqrt(sum((e$a - e$b)^2))))
  score <- nx * 1e6 + len
  if (score < best_score) { best_score <- score; best <- pp }
}
callouts <- callouts[c(best, which(!has_line))]
seg <- txt <- NULL
for (i in seq_along(callouts)) {
  cc <- callouts[[i]]; m <- n_reg(cc$key)
  lab <- if (cc$key == "001") cc$txt else sprintf("%d %s", m, cc$txt)
  if (m > 0 && m <= 5) lab <- sprintf("%s\n(%s)", lab, paste(names_of(cc$key), collapse = ", "))
  txt <- rbind(txt, data.frame(x = label_x, y = call_y[i], label = lab))
  a <- anchor_call(cc$key)
  if (m > 0 && !is.null(a)) seg <- rbind(seg, data.frame(x = a[1], y = a[2], xend = label_x - 0.5 * k, yend = call_y[i]))
}

big <- data.frame(set = c("DEG", "TX"), key = c("100", "010"),
                  label = c("gene DEGs only", "genes with a transcript\nDEG only"))
big$n <- sapply(big$key, n_reg)
big_xy <- t(sapply(big$key, anchor))
big$x <- big_xy[, 1]; big$y <- big_xy[, 2]
# gene DEG and transcript DEG lens is labelled in place
lens_xy <- anchor("110")
lens <- data.frame(x = lens_xy[1], y = lens_xy[2], n = n_reg("110"))

model_txt <- c(noCov = "~ condition", withCov = "~ Sex + Age + PMI + RIN + condition")[[model]]
leg <- c(DEG = sprintf("Differentially expressed genes (n = %s)", format(n[["DEG"]], big.mark = ",")),
         TX = sprintf("Genes with a differentially expressed transcript (n = %s)", format(n[["TX"]], big.mark = ",")),
         DIU = sprintf("Genes with differential isoform usage (n = %d)", n[["DIU"]]))

###############################
# Plot
###############################

p <- ggplot() +
  geom_polygon(data = circles, aes(x, y, fill = set, colour = set), alpha = 0.35, linewidth = 0.6) +
  geom_text(data = big, aes(x, y + 1.5 * k, label = format(n, big.mark = ",")),
            size = 7, fontface = "bold", colour = "grey15") +
  geom_text(data = big, aes(x, y - 2.5 * k, label = label), size = 3.2, colour = "grey25", lineheight = 0.95) +
  geom_text(data = lens, aes(x, y + 1.6 * k, label = n), size = 5, fontface = "bold", colour = "grey15") +
  geom_text(data = lens, aes(x, y - 1.4 * k), label = "gene and\ntranscript DEG",
            size = 2.6, colour = "grey25", lineheight = 0.9) +
  geom_segment(data = seg, aes(x, y, xend = xend, yend = yend), colour = "grey45", linewidth = 0.3) +
  geom_text(data = txt, aes(x, y, label = label), hjust = 0, vjust = 0.5, lineheight = 0.95,
            colour = "grey15", size = 3.2) +
  scale_fill_manual(values = COLS, breaks = c("DEG", "TX", "DIU"), labels = leg, name = NULL) +
  scale_colour_manual(values = COLS, guide = "none") +
  coord_equal(xlim = c(min(circles$x) - 1 * k, label_x + 26 * k),
              ylim = c(min(circles$y) - 1 * k, max(circles$y) + 1 * k), expand = FALSE) +
  guides(fill = guide_legend(ncol = 1, override.aes = list(alpha = 0.6))) +
  labs(title = "Most genes with isoform switching show no change in total expression",
       subtitle = sprintf(paste("Long-read RNA-seq, 6 AD vs 6 CT. Gene and transcript DEG: DESeq2 %s, FDR < 0.05, |log2FC| > 1.",
                                "\nDIU: IsoformSwitchAnalyzeR, FDR < 0.05, |dIF| > 0.1. Circle areas and pairwise overlaps",
                                "are proportional to gene counts."),
                          model_txt)) +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(size = 10.5, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey30", lineheight = 1.05),
        legend.position = "bottom", legend.justification = "left",
        plot.margin = margin(8, 8, 8, 8), plot.background = element_rect(fill = "white", colour = NA))

if (save_plots) {
  stem <- file.path(FIGS, paste0("reanalysis_deg_txdeg_diu_venn_", model))
  ggsave(paste0(stem, ".png"), p, width = 7.5, height = 6.5, dpi = 300)
  ggsave(paste0(stem, ".pdf"), p, width = 7.5, height = 6.5)
}
