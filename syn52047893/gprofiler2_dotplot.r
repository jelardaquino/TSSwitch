library(clusterProfiler)
library(org.Hs.eg.db)
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyverse)
library(gprofiler2)
library(biomaRt)
library(forcats)

res_df <- read.table("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/tables/suppa_beta_regression/suppa_beta_regression_all_results.tsv", header=TRUE, stringsAsFactors=FALSE, sep="\t")
# res_sig <- res_df %>% filter(p_adj < 0.05)
res_sig <- res_df %>% filter(p_adj < 0.05, abs(delta) > 0.2)

# ── g:Profiler GO enrichment ──────────────────────────────────────────────────
query_genes <- unique(res_sig$gene_name[!is.na(res_sig$gene_name)])
gost_res <- gost(
  query             = query_genes,
  organism          = "hsapiens")
# View results
gost_res$result %>%
  filter(p_value < 0.05) %>%
  arrange(p_value) %>%
  dplyr::select(source, term_name, p_value, term_size, intersection_size)
write.table(gost_res$result, file = "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/tables/suppa_beta_regression/gprofiler_results.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
png("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/figures/suppa_beta_regression/gprofiler_dotplot.png", width = 10, height = 8, units = "in", res = 300)   
gostplot(gost_res, capped = FALSE, interactive = FALSE)
dev.off()

dot_df <- gost_res$result %>%
  filter(p_value < 0.05) %>%
  filter(source != "TF") %>%
  group_by(source) %>%
  arrange(p_value) %>%
  slice_head(n = 15) %>%         # top 15 per source
  ungroup() %>%
  mutate(
    gene_ratio = intersection_size / term_size,
    term = fct_reorder(term_name, gene_ratio)
  )



library(ggplot2)
library(patchwork)
library(stringr)

dot_df$term_wrapped <- str_wrap(dot_df$term, width = 30)
size_range <- range(dot_df$intersection_size, na.rm = TRUE)

p_left <- ggplot(subset(dot_df, source == "GO:BP"),
                 aes(x = gene_ratio, y = reorder(term_wrapped, gene_ratio))) +
  geom_point(aes(size = intersection_size, color = p_value)) +
  scale_color_viridis_c(direction = -1, option = "magma") +
  facet_wrap(~ source, scales = "free_y") +  # gives real strip
  labs(y = "Terms", x = "Gene Ratio") +
  geom_point(aes(size = intersection_size, color = p_value)) +
  scale_size(range = c(2, 8), limits = size_range) +
  theme_bw(base_size = 12)

p_right <- ggplot(subset(dot_df, source != "GO:BP"),
                  aes(x = gene_ratio, y = reorder(term_wrapped, gene_ratio))) +
  geom_point(aes(size = intersection_size, color = p_value)) +
  facet_wrap(~ source, scales = "free_y", ncol = 1) +  # stack vertically
  scale_color_viridis_c(direction = -1, option = "magma") +
  geom_point(aes(size = intersection_size, color = p_value)) +
  scale_size(range = c(2, 8), limits = size_range) +
  labs(y = "Terms", x = "Gene Ratio") +
  theme_bw(base_size = 12)

(p_left | p_right) + plot_layout(guides = "collect")

ggsave("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis/figures/suppa_beta_regression/gprofiler_dotplot_facetwrap.png", width = 12, height = 8, dpi = 300)