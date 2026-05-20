###############################################################################
# Example Usage of Custom Switch Plot Helpers
###############################################################################

library(IsoformSwitchAnalyzeR)
library(ggplot2)

setwd("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/diu")
switchListAnalyzed <- readRDS("switchList_with_gene_name.rds")
outdir <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/diu/custom_plots"
# Source the custom helper functions
source("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/scripts/custom_switch_plot_helpers.R")

# ─── Define your custom colors ──────────────────────────────────────────────
# Named vector: condition name -> color
my_colors <- c("CT" = "turquoise", "AD" = "salmon")   # turquoise vs salmon
# my_colors <- c("CT" = "#1B9E77", "AD" = "#D95F02")  # teal vs orange
# my_colors <- c("CT" = "#7570B3", "AD" = "#E7298A")  # purple vs pink

# ─── Rename condition labels on plots ────────────────────────────────────────
# Named vector: original condition name -> display label
# my_labels <- c("CT" = "CT", "AD" = "AD")
my_labels <- NULL

# ─── Transcript colors ───────────────────────────────────────────────────────
# Set to NULL to auto-assign a publishable color-blind friendly palette.
# Or supply a named vector manually (names must match isoform IDs).
my_transcript_colors <- NULL


# ═══════════════════════════════════════════════════════════════════════════════
# OPTION 1: All-in-one combined plot (like switchPlot but with your colors)
# ═══════════════════════════════════════════════════════════════════════════════
result <- customSwitchPlot(
    switchAnalyzeRlist = switchListAnalyzed,
    gene = "NDUFB1",
    condition1 = "CT",
    condition2 = "AD",
    condition_colors = my_colors,
    condition_labels = my_labels,       # <-- rename conditions on plots
    transcript_colors = my_transcript_colors,  # <-- uncomment to color transcripts
    localTheme = theme_bw(base_size = 5),
    reverseMinus = FALSE,               # keep native strand direction for arrows
    plotTopology = FALSE
)

# Optional: inspect generated transcript names and colors used in the top panel
print(result$transcript_name_map)
print(result$transcript_color_map)

# Save the combined figure
ggsave(file.path(outdir, "NDUFB1_custom_switch_combined.pdf"), plot = result$combined,
       width = 14, height = 10, dpi = 300)
ggsave(file.path(outdir, "NDUFB1_custom_switch_combined.png"), plot = result$combined,
       width = 10, height = 10, dpi = 300)

# ═══════════════════════════════════════════════════════════════════════════════
# OPTION 2: Build each panel separately for full customization
# ═══════════════════════════════════════════════════════════════════════════════

# Step 1: Extract the data
dat <- extractSwitchData(
    switchListAnalyzed,
    gene = "NDUFB1",
    condition1 = "CT",
    condition2 = "AD"
)

# Step 2: Transcript structure plot (top panel)
p_transcripts <- plotTranscriptStructure(
    switchListAnalyzed,
    gene = "NDUFB1",
    condition1 = "CT",
    condition2 = "AD",
    transcript_colors = my_transcript_colors,  # <-- uncomment for custom colors
    localTheme = theme_bw(base_size = 13)
)
# You can add layers to this ggplot:
# p_transcripts <- p_transcripts + ggtitle("My Custom Transcript Title")

# Step 3: Gene expression barplot
p_gene <- plotGeneExpression(
    dat,
    condition_colors = my_colors,
    condition_labels = my_labels,
    localTheme = theme_bw(base_size = 13)
)
# Customize further:
# p_gene <- p_gene + theme(plot.title = element_text(face = "italic"))

# Step 4: Isoform expression barplot
p_isoexp <- plotIsoformExpression(
    dat,
    condition_colors = my_colors,
    condition_labels = my_labels,
    localTheme = theme_bw(base_size = 13)
)

# Step 5: Isoform usage barplot
p_usage <- plotIsoformUsage(
    dat,
    condition_colors = my_colors,
    condition_labels = my_labels,
    localTheme = theme_bw(base_size = 13)
)

# ─── Save individual panels ─────────────────────────────────────────────────
ggsave(file.path(outdir, "NDUFB1_transcripts.png"), p_transcripts, width = 12, height = 5, dpi = 300)
ggsave(file.path(outdir, "NDUFB1_gene_exp.png"), p_gene, width = 4, height = 4, dpi = 300)
ggsave(file.path(outdir, "NDUFB1_iso_exp.png"), p_isoexp, width = 6, height = 4, dpi = 300)
ggsave(file.path(outdir, "NDUFB1_iso_usage.png"), p_usage, width = 6, height = 4, dpi = 300)

# ─── Or manually arrange them however you like ──────────────────────────────
library(gridExtra)

bottom <- arrangeGrob(p_gene, p_isoexp, p_usage, ncol = 3, widths = c(1, 2, 2))
full <- arrangeGrob(p_transcripts, bottom, nrow = 2, heights = c(3, 2))
ggsave(file.path(outdir, "NDUFB1_custom_full.pdf"), full, width = 14, height = 10)


# ═══════════════════════════════════════════════════════════════════════════════
# OPTION 3: Heavy customization example
# ═══════════════════════════════════════════════════════════════════════════════

# Use a completely custom theme and color palette
fancy_colors <- c("CT" = "#2166AC", "AD" = "#B2182B")

fancy_theme <- theme_minimal(base_size = 14) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        panel.grid.minor = element_blank(),
        legend.position = "bottom"
    )

p_gene_fancy <- plotGeneExpression(
    dat,
    condition_colors = fancy_colors,
    localTheme = fancy_theme,
    title = expression(bold("Gene Expression: ") * italic("NDUFB1"))
)

p_isoexp_fancy <- plotIsoformExpression(
    dat,
    condition_colors = fancy_colors,
    localTheme = fancy_theme,
    title = expression(bold("Isoform Expression: ") * italic("NDUFB1"))
)

p_usage_fancy <- plotIsoformUsage(
    dat,
    condition_colors = fancy_colors,
    localTheme = fancy_theme,
    title = expression(bold("Isoform Usage: ") * italic("NDUFB1"))
)

ggsave(file.path(outdir, "NDUFB1_gene_exp_fancy.pdf"), p_gene_fancy, width = 5, height = 5)
ggsave(file.path(outdir, "NDUFB1_iso_exp_fancy.pdf"), p_isoexp_fancy, width = 7, height = 5)
ggsave(file.path(outdir, "NDUFB1_iso_usage_fancy.pdf"), p_usage_fancy, width = 7, height = 5)
