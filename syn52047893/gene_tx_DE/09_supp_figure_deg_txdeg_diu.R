###############################
# Supplementary Figure 5: gene DEG vs genes with a transcript DEG vs the 86 DIU genes, long-read data, DESeq2 ~ condition.
# Drawn by 08_reanalysis_deg_txdeg_diu_venn.R; title and subtitle go in the figure legend.
###############################

suppressPackageStartupMessages({
  library(ggplot2)
})

SCRIPTS <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/emtab14167_replication/scripts"
FIGS <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/ERP161086/emtab14167_replication/figures"

env <- new.env()
assign("model", "noCov", envir = env)
assign("save_plots", FALSE, envir = env)
sys.source(file.path(SCRIPTS, "08_reanalysis_deg_txdeg_diu_venn.R"), envir = env)
fig <- env$p + labs(title = NULL, subtitle = NULL)

stem <- file.path(FIGS, "SuppFig_gene_tx_DEG_DIU_venn")
ggsave(paste0(stem, ".pdf"), fig, width = 7, height = 5.8)
ggsave(paste0(stem, ".tiff"), fig, width = 7, height = 5.8, dpi = 600, compression = "lzw")
ggsave(paste0(stem, ".png"), fig, width = 7, height = 5.8, dpi = 300)
