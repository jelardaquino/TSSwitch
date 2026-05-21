###############################################################################
# go_kegg_enrichment.r
#
# Purpose : Run GO (Biological Process) and KEGG pathway enrichment analysis
#           on genes with significant differential alternative splicing
#           (identified by beta_regression_psi.r). Produces dot plots and a
#           cnetplot showing the relationship between top GO terms and genes.
#
# Inputs  : - suppa_beta_regression_analysis.RData : Workspace checkpoint from
#                                                    beta_regression_psi.r
#                                                    (contains res_sig data frame)
#
# Outputs : - GO_enrichment_SUPPA_significant_events_dotplot.png : GO BP dot plot
#           - KEGG_SUPPA_significant_events_dotplot.png           : KEGG dot plot
#           - GO_enrichment_SUPPA_significant_events_cnetplot.png : GO cnetplot
#                                                                   (top 4 terms)
#
# Prerequisite: Run beta_regression_psi.r first to generate the RData checkpoint.
###############################################################################

# ── Libraries ─────────────────────────────────────────────────────────────────
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)

# ── Path configuration ────────────────────────────────────────────────────────
# ⚠️ Update these paths to match your local environment before running.
BASE_DIR   <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/PRJNA1206164"
OUTPUT_DIR <- BASE_DIR

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Load results from beta regression checkpoint ──────────────────────────────
rdata_file <- file.path(BASE_DIR, "suppa_beta_regression_analysis.RData")
if (!file.exists(rdata_file)) {
  stop("RData checkpoint not found: ", rdata_file,
       "\nRun beta_regression_psi.r first.")
}
load(rdata_file)

# Verify res_sig is available and has gene_name column
stopifnot(
  exists("res_sig"),
  "gene_name" %in% colnames(res_sig)
)

# Use gene symbols for enrichment; drop missing gene names
query_genes <- unique(res_sig$gene_name[!is.na(res_sig$gene_name) & nzchar(res_sig$gene_name)])
message("Genes submitted to enrichment: ", length(query_genes))

# ── GO Biological Process enrichment ─────────────────────────────────────────
ego <- enrichGO(
  gene          = query_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",            # Biological Process
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

if (nrow(as.data.frame(ego)) == 0) {
  message("No significant GO terms found at q ≤ 0.05.")
} else {
  message("Significant GO BP terms: ", nrow(as.data.frame(ego)))

  # Dot plot: top 20 GO BP terms
  png(
    file.path(OUTPUT_DIR, "GO_enrichment_SUPPA_significant_events_dotplot.png"),
    width = 8, height = 12, units = "in", res = 300
  )
  print(dotplot(ego, showCategory = 20,
                title = "GO Biological Process — PRJNA1206164 Significant SUPPA2 Events"))
  dev.off()
  message("GO dot plot saved.")

  # ── cnetplot: gene-to-term network for top 4 GO terms ──────────────────────
  # Uses ΔPSI as fold change values to color genes in the network
  delta_vec        <- res_sig$delta
  names(delta_vec) <- res_sig$gene_name
  delta_vec        <- sort(na.omit(delta_vec), decreasing = TRUE)

  png(
    file.path(OUTPUT_DIR, "GO_enrichment_SUPPA_significant_events_cnetplot.png"),
    width = 8, height = 6, units = "in", res = 300
  )
  cnet_plot <- cnetplot(
    ego,
    categorySize = "pvalue",
    foldChange   = delta_vec,  # named numeric vector (gene → ΔPSI)
    showCategory = 4
  )
  # Relabel color legend as ΔPSI
  print(cnet_plot + guides(color = guide_colorbar(title = expression(Delta * "PSI"))))
  dev.off()
  message("GO cnetplot saved.")
}

# ── KEGG pathway enrichment ───────────────────────────────────────────────────
# Convert gene symbols to Entrez IDs (required by enrichKEGG)
gene_entrez <- bitr(
  query_genes,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)
message("Genes successfully mapped to Entrez IDs: ", nrow(gene_entrez))

ekegg <- enrichKEGG(
  gene         = gene_entrez$ENTREZID,
  organism     = "hsa",   # Homo sapiens
  pvalueCutoff = 0.05
)

if (nrow(as.data.frame(ekegg)) == 0) {
  message("No significant KEGG pathways found at p ≤ 0.05.")
} else {
  message("Significant KEGG pathways: ", nrow(as.data.frame(ekegg)))

  png(
    file.path(OUTPUT_DIR, "KEGG_SUPPA_significant_events_dotplot.png"),
    width = 6, height = 4, units = "in", res = 300
  )
  print(dotplot(ekegg, showCategory = 20,
                title = "KEGG Pathway Enrichment — PRJNA1206164 Significant SUPPA2 Events"))
  dev.off()
  message("KEGG dot plot saved.")
}

# ── Session info for reproducibility ─────────────────────────────────────────
message("\n── Session Info ──")
sessionInfo()
