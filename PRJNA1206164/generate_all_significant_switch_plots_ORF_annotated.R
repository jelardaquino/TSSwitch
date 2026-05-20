###############################################################################
# Generate custom switch plots for all significant genes in switchListAnalyzed
###############################################################################

suppressPackageStartupMessages({
    library(IsoformSwitchAnalyzeR)
    library(ggplot2)
    library(dplyr)
    library(BSgenome.Hsapiens.UCSC.hg38)
})

base_dir <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/PRJNA1206164"
diu_dir <- file.path(base_dir, "diu")

switch_rds <- file.path(diu_dir, "switchList.rds")
helper_script <- file.path(diu_dir, "custom_switch_plot_helpers.R")
outdir <- diu_dir

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

source(helper_script)
switchListAnalyzed <- readRDS(switch_rds)

# Annotate ORFs so transcript structures can indicate coding/NMD status
switchListAnalyzed <- analyzeORF(
    switchListAnalyzed,
    genomeObject = BSgenome.Hsapiens.UCSC.hg38
)

# Plot settings
condition1 <- "CT"
condition2 <- "AD"
q_cutoff <- 0.05
# Optional: limit number of genes for quick test runs (set to NA for all)
max_genes <- NA_integer_

condition_colors <- c("CT" = "turquoise", "AD" = "salmon")
condition_labels <- NULL
transcript_colors <- NULL  # auto color-blind palette from helper

features <- switchListAnalyzed$isoformFeatures

sig_tbl <- features %>%
    filter(
        condition_1 == condition1,
        condition_2 == condition2,
        !is.na(isoform_switch_q_value),
        isoform_switch_q_value < q_cutoff
    ) %>%
    mutate(
        gene_label = ifelse(!is.na(gene_name) & nzchar(gene_name), gene_name, gene_id)
    ) %>%
    group_by(gene_id, gene_label) %>%
    summarize(min_q = min(isoform_switch_q_value, na.rm = TRUE), .groups = "drop") %>%
    arrange(min_q, gene_label)

if (nrow(sig_tbl) == 0) {
    message("No significant genes found for ", condition1, " vs ", condition2, " at q < ", q_cutoff)
    quit(save = "no", status = 0)
}

if (!is.na(max_genes)) {
    sig_tbl <- head(sig_tbl, max_genes)
}

safe_name <- function(x) {
    gsub("[^A-Za-z0-9._-]", "_", x)
}

status <- vector("list", nrow(sig_tbl))
combined_pages <- list()
combined_pdf_file <- file.path(outdir, "all_significant_switchPlots.pdf")

sig_tbl <- sig_tbl %>%
    filter(gene_label == "COL4A2" | gene_label == "SPP1")

for (i in seq_len(nrow(sig_tbl))) {
    gene_id <- sig_tbl$gene_id[i]
    gene_label <- sig_tbl$gene_label[i]
    min_q <- sig_tbl$min_q[i]

    file_prefix <- sprintf("%04d_%s", i, safe_name(gene_label))

    res <- tryCatch(
        {
            # customSwitchPlot() draws by design; route that draw to a null PDF
            # device so we can safely control final combined-PDF writing later.
            grDevices::pdf(file = NULL)

            customSwitchPlot(
                switchAnalyzeRlist = switchListAnalyzed,
                gene = gene_id,
                condition1 = condition1,
                condition2 = condition2,
                condition_colors = condition_colors,
                condition_labels = condition_labels,
                transcript_colors = transcript_colors,
                localTheme = theme_bw(base_size = 12),
                reverseMinus = FALSE,
                plotTopology = FALSE
            )
        },
        error = function(e) e,
        finally = {
            if (grDevices::dev.cur() > 1) {
                try(grDevices::dev.off(), silent = TRUE)
            }
        }
    )

    if (inherits(res, "error")) {
        status[[i]] <- data.frame(
            gene_id = gene_id,
            gene_label = gene_label,
            min_q = min_q,
            ok = FALSE,
            message = conditionMessage(res),
            stringsAsFactors = FALSE
        )
        next
    }

    png_file <- file.path(outdir, paste0(file_prefix, "_switchPlot.png"))

    ggsave(png_file, plot = res$combined, width = 14, height = 10, dpi = 300)
    combined_pages[[length(combined_pages) + 1]] <- res$combined

    status[[i]] <- data.frame(
        gene_id = gene_id,
        gene_label = gene_label,
        min_q = min_q,
        ok = TRUE,
        message = "ok",
        stringsAsFactors = FALSE
    )

    message(sprintf("[%d/%d] %s (%s) done", i, nrow(sig_tbl), gene_label, gene_id))
}

status_df <- bind_rows(status)
write.csv(status_df, file.path(outdir, "plot_generation_status.csv"), row.names = FALSE)

if (length(combined_pages) > 0) {
    tmp_pdf <- paste0(combined_pdf_file, ".tmp")
    grDevices::pdf(
        file = tmp_pdf,
        width = 14,
        height = 10,
        onefile = TRUE,
        paper = "special"
    )
    for (pg in combined_pages) {
        grid::grid.newpage()
        grid::grid.draw(pg)
    }
    grDevices::dev.off()

    if (file.exists(combined_pdf_file)) {
        file.remove(combined_pdf_file)
    }
    ok_rename <- file.rename(tmp_pdf, combined_pdf_file)
    if (!ok_rename) {
        file.copy(tmp_pdf, combined_pdf_file, overwrite = TRUE)
        file.remove(tmp_pdf)
    }
}

message("Done. Success: ", sum(status_df$ok), " / ", nrow(status_df))
message("Combined PDF: ", combined_pdf_file)
