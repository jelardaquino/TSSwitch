###############################
# Summary figure for the short-read (E-MTAB-14167, 41 donors) DIU replication:
#   A  overlap of long-read and short-read DIU genes vs chance, at three short-read thresholds
#   B  long-read vs short-read dIF for isoforms of long-read DIU genes tested in both
#   C  long-read and short-read dIF for switching isoforms significant in both
# Short-read dIF values are IsoformSwitchAnalyzeR's covariate-adjusted (age, sex, SV) values,
# as reported in the text.
#
# Run with the isoformswitch_env R (base R + ggplot2 + gridExtra).
###############################

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
})

S <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing"
OUT <- file.path(S, "ERP161086/emtab14167_replication")
DIU <- file.path(OUT, "diu_bambu")
FIGS <- file.path(OUT, "figures")

LR_DIU_TOP <- file.path(S, "syn52065646_LongRead/jaquino_analysis/diu/IsoformSwitchAnalyzeR_top_switches_gene_name.tsv")
LR_TESTED_RDS <- file.path(S, "syn52065646_LongRead/jaquino_analysis/diu_satuRn/switchList_satuRn_all_genes.rds")

COL_LR <- "#3A6FB0"
COL_SR <- "#D9822B"
COL_BG <- "#C9CDD3"

ISO_LABELS <- c(ENST00000266546 = "CLSTN3-01", ENST00000331148 = "CLSTN3-02",
                ENST00000648222 = "COL4A2", ENST00000307765 = "RXFP1")

base_theme <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey93", linewidth = 0.3),
        panel.border = element_rect(colour = "grey70"),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        legend.position = "bottom", legend.key.size = unit(0.35, "cm"))

fmt_p <- function(x) formatC(x, format = "f", digits = 2)

###############################
# A: gene overlap vs chance
###############################

lr_sig <- unique(read.delim(LR_DIU_TOP)$gene_id)
lr_tested <- unique(readRDS(LR_TESTED_RDS)$isoformFeatures$gene_id)
sr <- read.delim(file.path(DIU, "SR41_isoform_results_all_tested.tsv"))
universe <- intersect(lr_tested, unique(sr$gene_id))

thresholds <- data.frame(fdr = c(0.05, 0.1, 0.1), dif = c(0.1, 0.1, 0.05))
ov <- do.call(rbind, lapply(seq_len(nrow(thresholds)), function(i) {
  q <- sr$isoform_switch_q_value
  sr_sig <- unique(sr$gene_id[!is.na(q) & q < thresholds$fdr[i] & abs(sr$dIF) > thresholds$dif[i]])
  a <- universe %in% lr_sig
  b <- universe %in% sr_sig
  f <- fisher.test(table(factor(a, c(TRUE, FALSE)), factor(b, c(TRUE, FALSE))), alternative = "greater")
  data.frame(threshold = sprintf("FDR < %s\n|dIF| > %s", thresholds$fdr[i], thresholds$dif[i]),
             n_sr = sum(b), observed = sum(a & b), expected = sum(a) * sum(b) / length(universe),
             p = f$p.value)
}))
ov$threshold <- factor(ov$threshold, levels = ov$threshold)
ov_long <- rbind(data.frame(threshold = ov$threshold, what = "Observed", n = ov$observed),
                 data.frame(threshold = ov$threshold, what = "Expected by chance", n = ov$expected))
ov_long$what <- factor(ov_long$what, levels = c("Observed", "Expected by chance"))

pA <- ggplot(ov_long, aes(threshold, n, fill = what)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, colour = "white", linewidth = 0.6) +
  geom_text(aes(label = ifelse(what == "Observed", n, sprintf("%.1f", n))),
            position = position_dodge(width = 0.75), vjust = -0.4, size = 2.8, colour = "grey15") +
  geom_text(data = ov, aes(threshold, pmax(observed, expected) + 1.6, label = paste0("Fisher p = ", fmt_p(p))),
            inherit.aes = FALSE, size = 2.7, colour = "grey30") +
  scale_fill_manual(values = c("Observed" = COL_LR, "Expected by chance" = COL_BG), name = NULL) +
  scale_y_continuous(limits = c(0, max(ov$observed, ov$expected) + 2.5), expand = expansion(mult = c(0, 0))) +
  labs(title = "A  Long-read DIU genes also significant in short-read data",
       subtitle = sprintf("%d of 86 long-read DIU genes tested; %s genes tested in both",
                          sum(universe %in% lr_sig), format(length(universe), big.mark = ",")),
       x = "Short-read significance threshold", y = "Overlapping genes") +
  base_theme + theme(panel.grid.major.x = element_blank())

###############################
# B: dIF scatter
###############################

cmp <- read.delim(file.path(DIU, "comparison_dIF_LR_vs_SR41_isoforms_in_LR_DIU_genes.tsv"))
cmp$lr_switch <- !is.na(cmp$q_LR) & cmp$q_LR < 0.05 & abs(cmp$dIF_LR) > 0.1
cmp$both <- cmp$lr_switch & cmp$significant_SR
cmp$cls <- ifelse(cmp$both, "Significant in both",
                  ifelse(cmp$lr_switch, "Long-read switching isoform", "Other isoform of a DIU gene"))
cmp$cls <- factor(cmp$cls, levels = c("Other isoform of a DIU gene", "Long-read switching isoform", "Significant in both"))
cmp <- cmp[order(cmp$cls), ]

rho <- cor(cmp$dIF_LR, cmp$dIF_SR, method = "spearman")
sw <- cmp[cmp$lr_switch, ]
same <- sum(sign(sw$dIF_LR) == sign(sw$dIF_SR))
binom_p <- binom.test(same, nrow(sw), 0.5, alternative = "greater")$p.value

lab_b <- cmp[cmp$both, ]
lab_b$label <- ISO_LABELS[lab_b$isoform_id]
lab_b$hjust <- ifelse(lab_b$dIF_LR > 0, 1.15, -0.15)

lim <- c(-0.8, 0.8)
pB <- ggplot(cmp, aes(dIF_LR, dIF_SR)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, colour = "grey75", linewidth = 0.3, linetype = "dashed") +
  geom_point(aes(fill = cls, size = cls), shape = 21, colour = "white", stroke = 0.4) +
  geom_text(data = lab_b, aes(label = label, hjust = hjust), size = 2.8, colour = "grey15") +
  annotate("text", x = lim[1], y = lim[2], hjust = 0, vjust = 1, size = 2.8, colour = "grey15", lineheight = 0.95,
           label = sprintf("Spearman rho = %.3f (n = %d isoforms)\nSwitching isoforms, same direction: %d of %d\n(one-sided binomial p = %s)",
                           rho, nrow(cmp), same, nrow(sw), fmt_p(binom_p))) +
  scale_fill_manual(values = c("Other isoform of a DIU gene" = COL_BG, "Long-read switching isoform" = COL_LR,
                               "Significant in both" = COL_SR), name = NULL) +
  scale_size_manual(values = c(1.6, 2.2, 2.8), guide = "none") +
  coord_equal(xlim = lim, ylim = lim) +
  guides(fill = guide_legend(nrow = 2, override.aes = list(size = 2.5))) +
  labs(title = "B  Change in isoform fraction, long-read vs short-read",
       subtitle = "Isoforms of the 86 long-read DIU genes tested in both datasets",
       x = "Long-read dIF (AD - CT)", y = "Short-read dIF (AD - CT)") +
  base_theme

###############################
# C: isoforms significant in both
###############################

sig <- read.delim(file.path(DIU, "comparison_LR_significant_isoforms_in_SR41.tsv"))
sig <- sig[sig$significant_SR %in% TRUE, ]
sig$label <- ISO_LABELS[sig$isoform_id]
sig$label <- factor(sig$label, levels = rev(ISO_LABELS[ISO_LABELS %in% sig$label]))
sig_long <- rbind(
  data.frame(label = sig$label, platform = "Long-read (6 vs 6)", dIF = sig$dIF_LR, q = sig$q_LR),
  data.frame(label = sig$label, platform = "Short-read (20 vs 21)", dIF = sig$dIF_SR, q = sig$q_SR)
)
sig_long$platform <- factor(sig_long$platform, levels = c("Long-read (6 vs 6)", "Short-read (20 vs 21)"))
sig_long$q_lab <- ifelse(sig_long$q < 0.001, sprintf("FDR %.0e", sig_long$q), sprintf("FDR %.3f", sig_long$q))
dir_lab <- data.frame(label = sig$label, text = ifelse(sig$same_direction, "same direction", "opposite direction"))

pC <- ggplot(sig_long, aes(dIF, label)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_line(aes(group = label), colour = "grey70", linewidth = 0.6) +
  geom_point(aes(fill = platform), shape = 21, colour = "white", stroke = 0.4, size = 3) +
  geom_text(aes(label = q_lab, vjust = ifelse(platform == "Long-read (6 vs 6)", -1.1, 2)),
            size = 2.4, colour = "grey30") +
  geom_text(data = dir_lab, aes(x = 0.88, y = label, label = text), hjust = 0, size = 2.6, colour = "grey15") +
  scale_fill_manual(values = c("Long-read (6 vs 6)" = COL_LR, "Short-read (20 vs 21)" = COL_SR), name = NULL) +
  scale_x_continuous(limits = c(-0.8, 1.45), breaks = seq(-0.75, 0.75, 0.25)) +
  scale_y_discrete(expand = expansion(add = 0.7)) +
  labs(title = "C  Switching isoforms significant in both datasets",
       subtitle = "Short-read FDR < 0.05 and |dIF| > 0.1",
       x = "dIF (AD - CT)", y = NULL) +
  base_theme + theme(panel.grid.major.y = element_blank())

###############################
# Assemble
###############################

caption <- grid::textGrob(
  paste("Short-read: E-MTAB-14167 prefrontal cortex, quantified against the long-read Bambu annotation and analyzed with the same",
        "IsoformSwitchAnalyzeR/DEXSeq workflow (~ condition + age + sex + surrogate variables).",
        "dIF values are IsoformSwitchAnalyzeR's covariate-adjusted estimates (expression adjusted for age, sex and surrogate variables before computing isoform fractions).", sep = "\n"),
  x = 0.01, hjust = 0, gp = grid::gpar(fontsize = 7, col = "grey30"))

fig <- arrangeGrob(
  arrangeGrob(pA, pC, ncol = 1, heights = c(1, 1)),
  pB, ncol = 2, widths = c(1, 1.05),
  bottom = caption
)
ggsave(file.path(FIGS, "diu_replication_summary.png"), fig, width = 10, height = 7.2, dpi = 300)
ggsave(file.path(FIGS, "diu_replication_summary.pdf"), fig, width = 10, height = 7.2)

print(ov)
cat("Spearman rho", rho, "n", nrow(cmp), "| same direction", same, "of", nrow(sw), "binomial p", binom_p, "\n")
print(sig[, c("label", "dIF_LR", "q_LR", "dIF_SR", "q_SR", "same_direction")])
