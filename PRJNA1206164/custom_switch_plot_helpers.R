###############################################################################
# Custom Helper Functions for IsoformSwitchAnalyzeR Plotting
# 
# These functions extract data from a switchAnalyzeRlist and produce individual
# ggplot2 objects that you can customize with your own colors, themes, etc.
#
# The main functions:
#   1. extractSwitchData()       - Extracts all data needed for plotting
#   2. plotTranscriptStructure() - Top panel: transcript/exon structure
#   3. plotGeneExpression()      - Bottom-left: gene expression barplot
#   4. plotIsoformExpression()   - Bottom-middle: isoform expression barplot
#   5. plotIsoformUsage()        - Bottom-right: isoform usage (IF) barplot
#   6. customSwitchPlot()        - Combines all panels into one figure
###############################################################################

library(IsoformSwitchAnalyzeR)
library(ggplot2)
library(dplyr)
library(reshape2)
library(grid)
library(gridExtra)

###############################################################################
# 1. extractSwitchData: Extract expression and usage data from switchAnalyzeRlist
###############################################################################
extractSwitchData <- function(
    switchAnalyzeRlist,
    gene,
    condition1,
    condition2,
    IFcutoff = 0.05
) {
    features <- switchAnalyzeRlist$isoformFeatures
    
    # Resolve gene to gene_id
    if (tolower(gene) %in% tolower(features$gene_id)) {
        gene_id <- gene
    } else if (tolower(gene) %in% tolower(features$gene_name)) {
        gene_id <- unique(features$gene_id[which(tolower(features$gene_name) %in% tolower(gene))])
        if (length(gene_id) > 1) stop("Gene name maps to multiple gene_ids. Use gene_id directly.")
    } else {
        stop(paste("Gene", gene, "not found in switchAnalyzeRlist."))
    }
    
    # Get isoform IDs for this gene
    isoform_ids <- unique(features$isoform_id[features$gene_id == gene_id])
    
    # Get gene name
    gene_name <- unique(features$gene_name[features$gene_id == gene_id])
    gene_name <- gene_name[!is.na(gene_name)][1]
    if (is.na(gene_name)) gene_name <- gene_id
    
    # Get rows for this comparison
    rows <- which(
        features$isoform_id %in% isoform_ids &
        features$condition_1 == condition1 &
        features$condition_2 == condition2
    )
    
    if (length(rows) == 0) {
        # Try reversed conditions
        rows <- which(
            features$isoform_id %in% isoform_ids &
            features$condition_1 == condition2 &
            features$condition_2 == condition1
        )
        if (length(rows) > 0) {
            temp <- condition1
            condition1 <- condition2
            condition2 <- temp
        } else {
            stop("No data found for this gene/condition combination.")
        }
    }
    
    sub <- features[rows, ]
    
    # Filter by IFcutoff
    sub <- sub[which(pmax(sub$IF1, sub$IF2, na.rm = TRUE) >= IFcutoff), ]
    if (nrow(sub) == 0) stop("No isoforms left after IFcutoff filter.")
    
    isoform_ids <- unique(sub$isoform_id)
    short_prefix <- if (!is.na(gene_name) && nzchar(gene_name)) gene_name else gene_id
    isoform_short_labels <- setNames(
        paste0(short_prefix, "-", sprintf("%02d", seq_along(isoform_ids))),
        isoform_ids
    )
    
    # --- Gene Expression Data ---
    gene_exp <- data.frame(
        Condition = c(condition1, condition2),
        Expression = c(sub$gene_value_1[1], sub$gene_value_2[1]),
        StdErr = c(sub$gene_stderr_1[1], sub$gene_stderr_2[1]),
        stringsAsFactors = FALSE
    )
    gene_q <- sub$gene_q_value[1]
    
    # --- Isoform Expression Data ---
    iso_exp <- data.frame(
        isoform_id = rep(sub$isoform_id, 2),
        Condition = c(rep(condition1, nrow(sub)), rep(condition2, nrow(sub))),
        Expression = c(sub$iso_value_1, sub$iso_value_2),
        StdErr = c(sub$iso_stderr_1, sub$iso_stderr_2),
        stringsAsFactors = FALSE
    )
    iso_q <- data.frame(
        isoform_id = sub$isoform_id,
        q_value = sub$iso_q_value,
        stringsAsFactors = FALSE
    )
    
    # --- Isoform Usage (IF) Data ---
    iso_usage <- data.frame(
        isoform_id = rep(sub$isoform_id, 2),
        Condition = c(rep(condition1, nrow(sub)), rep(condition2, nrow(sub))),
        IF = c(sub$IF1, sub$IF2),
        stringsAsFactors = FALSE
    )
    switch_q <- data.frame(
        isoform_id = sub$isoform_id,
        q_value = sub$isoform_switch_q_value,
        dIF = sub$dIF,
        stringsAsFactors = FALSE
    )
    
    # Number of replicates per condition (for CI calculation)
    nrReplicates <- switchAnalyzeRlist$conditions
    
    return(list(
        gene_id = gene_id,
        gene_name = gene_name,
        isoform_ids = isoform_ids,
    isoform_short_labels = isoform_short_labels,
        condition1 = condition1,
        condition2 = condition2,
        gene_expression = gene_exp,
        gene_q_value = gene_q,
        isoform_expression = iso_exp,
        isoform_q_values = iso_q,
        isoform_usage = iso_usage,
        switch_q_values = switch_q,
        nrReplicates = nrReplicates
    ))
}


###############################################################################
# Helper: significance label from q-value
###############################################################################
.sigLabel <- function(pValue, alphas = c(0.05, 0.001)) {
    sapply(pValue, function(x) {
        if (is.na(x)) return("NA")
        if (x < min(alphas)) return("***")
        if (x < max(alphas)) return("*")
        return("ns")
    })
}


###############################################################################
# Helpers: transcript label/color handling
###############################################################################
.okabeItoPalette <- function(n) {
    base_cols <- c(
        "#000000", "#E69F00", "#56B4E9", "#009E73",
        "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
    )

    if (n <= length(base_cols)) {
        return(base_cols[seq_len(n)])
    }

    extra_n <- n - length(base_cols)
    extra_cols <- grDevices::hcl(
        h = seq(15, 375, length.out = extra_n + 1)[1:extra_n],
        c = 100,
        l = 60
    )
    c(base_cols, extra_cols)
}

.resolveGeneDisplayName <- function(switchAnalyzeRlist, gene) {
    features <- switchAnalyzeRlist$isoformFeatures

    if (tolower(gene) %in% tolower(features$gene_id)) {
        gene_id <- unique(features$gene_id[tolower(features$gene_id) == tolower(gene)])[1]
    } else if (tolower(gene) %in% tolower(features$gene_name)) {
        gene_id <- unique(features$gene_id[tolower(features$gene_name) == tolower(gene)])[1]
    } else {
        return(gene)
    }

    gene_name <- unique(features$gene_name[features$gene_id == gene_id])
    gene_name <- gene_name[!is.na(gene_name)][1]
    if (is.na(gene_name) || !nzchar(gene_name)) gene_name <- gene_id
    gene_name
}

.resolveGeneId <- function(switchAnalyzeRlist, gene) {
    features <- switchAnalyzeRlist$isoformFeatures

    if (tolower(gene) %in% tolower(features$gene_id)) {
        local_gene_id <- unique(features$gene_id[tolower(features$gene_id) == tolower(gene)])[1]
        local_gene_name <- unique(features$gene_name[features$gene_id == local_gene_id])
        local_gene_name <- local_gene_name[!is.na(local_gene_name)][1]
        if (!is.na(local_gene_name) && nzchar(local_gene_name)) return(local_gene_name)
        return(local_gene_id)
    } else if (tolower(gene) %in% tolower(features$gene_name)) {
        local_gene_name <- unique(features$gene_name[tolower(features$gene_name) == tolower(gene)])[1]
        if (!is.na(local_gene_name) && nzchar(local_gene_name)) return(local_gene_name)
        return(unique(features$gene_id[tolower(features$gene_name) == tolower(gene)])[1])
    }

    gene
}


###############################################################################
# 2. plotTranscriptStructure: Transcript map (top panel)
#    This is a thin wrapper around switchPlotTranscript() that returns a
#    ggplot you can further modify (add colors, themes, etc.)
###############################################################################
plotTranscriptStructure <- function(
    switchAnalyzeRlist,
    gene,
    condition1,
    condition2,
    IFcutoff = 0.05,
    dIFcutoff = 0.1,
    rescaleTranscripts = TRUE,
    reverseMinus = FALSE,
    plotORF = NULL,
    plotTopology = TRUE,
    transcript_colors = NULL,    # named vector: c("ENST00000340900" = "red", ...)
    localTheme = theme_bw(base_size = 13),
    ...
) {
    # Work around IsoformSwitchAnalyzeR edge-cases in ORF-less objects where
    # switchPlotTranscript() may fail when ORF rows are missing for isoforms.
    required_orf_cols <- c(
        "isoform_id", "orfStartGenomic", "orfEndGenomic",
        "wasTrimmed", "trimmedStartGenomic"
    )

    if (is.null(switchAnalyzeRlist$orfAnalysis)) {
        switchAnalyzeRlist$orfAnalysis <- data.frame(stringsAsFactors = FALSE)
    }

    n_orf_rows <- nrow(switchAnalyzeRlist$orfAnalysis)

    # Ensure required columns exist with sensible defaults.
    if (!"isoform_id" %in% colnames(switchAnalyzeRlist$orfAnalysis)) {
        switchAnalyzeRlist$orfAnalysis$isoform_id <- rep(NA_character_, n_orf_rows)
    }
    if (!"orfStartGenomic" %in% colnames(switchAnalyzeRlist$orfAnalysis)) {
        switchAnalyzeRlist$orfAnalysis$orfStartGenomic <- rep(NA_real_, n_orf_rows)
    }
    if (!"orfEndGenomic" %in% colnames(switchAnalyzeRlist$orfAnalysis)) {
        switchAnalyzeRlist$orfAnalysis$orfEndGenomic <- rep(NA_real_, n_orf_rows)
    }
    if (!"wasTrimmed" %in% colnames(switchAnalyzeRlist$orfAnalysis)) {
        switchAnalyzeRlist$orfAnalysis$wasTrimmed <- rep(FALSE, n_orf_rows)
    }
    if (!"trimmedStartGenomic" %in% colnames(switchAnalyzeRlist$orfAnalysis)) {
        switchAnalyzeRlist$orfAnalysis$trimmedStartGenomic <- rep(NA_real_, n_orf_rows)
    }

    # Identify isoforms likely to be plotted for this gene/condition pair.
    features <- switchAnalyzeRlist$isoformFeatures
    if (tolower(gene) %in% tolower(features$gene_id)) {
        gene_id <- gene
    } else if (tolower(gene) %in% tolower(features$gene_name)) {
        gene_id <- unique(features$gene_id[which(tolower(features$gene_name) %in% tolower(gene))])[1]
    } else {
        gene_id <- gene
    }

    local_rows <- which(
        features$gene_id == gene_id &
        features$condition_1 == condition1 &
        features$condition_2 == condition2
    )
    if (!length(local_rows)) {
        local_rows <- which(features$gene_id == gene_id)
    }
    local_isoforms <- unique(features$isoform_id[local_rows])

    # Add dummy non-coding ORF rows for any missing isoforms so downstream code
    # gets nrow(orfInfo) > 0 and handles them as non-coding.
    missing_isoforms <- setdiff(local_isoforms, switchAnalyzeRlist$orfAnalysis$isoform_id)
    if (length(missing_isoforms)) {
        dummy_orf <- data.frame(
            isoform_id = missing_isoforms,
            orfStartGenomic = rep(NA_real_, length(missing_isoforms)),
            orfEndGenomic = rep(NA_real_, length(missing_isoforms)),
            wasTrimmed = rep(FALSE, length(missing_isoforms)),
            trimmedStartGenomic = rep(NA_real_, length(missing_isoforms)),
            stringsAsFactors = FALSE
        )

        switchAnalyzeRlist$orfAnalysis <- rbind(
            switchAnalyzeRlist$orfAnalysis[, required_orf_cols, drop = FALSE],
            dummy_orf
        )
    }

    # NOTE: 'plotORF' is kept only for backwards-compatible helper calls.
    # switchPlotTranscript() in current IsoformSwitchAnalyzeR versions does
    # not expose a plotORF argument.
    if (!is.null(plotORF)) {
        warning(
            "'plotORF' is not a switchPlotTranscript() argument in this IsoformSwitchAnalyzeR version and will be ignored."
        )
    }

    # switchPlotTranscript already returns a ggplot object
    p <- switchPlotTranscript(
        switchAnalyzeRlist = switchAnalyzeRlist,
        gene = gene,
        condition1 = condition1,
        condition2 = condition2,
        IFcutoff = IFcutoff,
        dIFcutoff = dIFcutoff,
        rescaleTranscripts = rescaleTranscripts,
        reverseMinus = reverseMinus,
        plotTopology = plotTopology,
        localTheme = localTheme,
        ...
    )
    
    # Append short transcript labels (GeneID-01, GeneID-02, ...) as an
    # additional row and apply color-blind friendly transcript colors
    # (or user-supplied colors).
    transcript_name_map <- NULL
    transcript_color_map <- NULL
    transcript_order_ids <- NULL

    # Find the rect layer that draws exons (has a "transcript" column)
    for (i in seq_along(p$layers)) {
        ld <- p$layers[[i]]$data
        if (!is.null(ld) && is.data.frame(ld) && "transcript" %in% names(ld)) {
            full_labels <- unique(as.character(ld$transcript))
            short_ids <- sub("\\n.*$", "", full_labels)
            transcript_order_ids <- short_ids

            gene_id_display <- .resolveGeneId(switchAnalyzeRlist, gene)
            short_labels <- paste0(gene_id_display, "-", sprintf("%02d", seq_along(full_labels)))
            names(short_labels) <- full_labels

            # Keep original label and add short label on a new line.
            display_labels <- paste0(full_labels, "\n", short_labels)
            names(display_labels) <- full_labels
            transcript_name_map <- setNames(short_labels, short_ids)

            # Update y-axis labels used by switchPlotTranscript.
            for (s_idx in seq_along(p$scales$scales)) {
                sc <- p$scales$scales[[s_idx]]
                if ("y" %in% sc$aesthetics && !is.null(sc$labels)) {
                    sc_labels <- as.character(sc$labels)
                    mapped <- display_labels[sc_labels]
                    mapped[is.na(mapped)] <- sc_labels[is.na(mapped)]
                    p$scales$scales[[s_idx]]$labels <- unname(mapped)
                }
            }

            # Build mapping from full transcript label -> color
            color_map <- character(length(full_labels))
            names(color_map) <- full_labels

            if (is.null(transcript_colors)) {
                auto_cols <- .okabeItoPalette(length(full_labels))
                color_map[full_labels] <- auto_cols
            } else {
                for (fl in full_labels) {
                    matched <- FALSE
                    for (nm in names(transcript_colors)) {
                        if (startsWith(as.character(fl), nm)) {
                            color_map[fl] <- transcript_colors[nm]
                            matched <- TRUE
                            break
                        }
                    }
                    if (!matched) color_map[fl] <- "#161616"
                }
            }

            # Apply renamed transcript labels and fill mapping for exon layers
            mapped_labels <- unname(display_labels[as.character(ld$transcript)])
            p$layers[[i]]$data$transcript <- factor(mapped_labels, levels = unname(display_labels[full_labels]))
            p$layers[[i]]$mapping$fill <- ggplot2::aes(fill = transcript)$fill

            # scale_fill_manual must use the renamed levels as names
            renamed_color_map <- color_map
            names(renamed_color_map) <- unname(display_labels[names(color_map)])
            transcript_color_map <- setNames(unname(renamed_color_map), unname(display_labels[full_labels]))

            suppressMessages(
                p <- p + scale_fill_manual(values = renamed_color_map, guide = "none")
            )
            break
        }
    }

    attr(p, "transcript_name_map") <- transcript_name_map
    attr(p, "transcript_color_map") <- transcript_color_map
    attr(p, "transcript_order_ids") <- transcript_order_ids
    
    return(p)
}


###############################################################################
# 3. plotGeneExpression: Gene-level expression barplot
###############################################################################
plotGeneExpression <- function(
    data,                        # output from extractSwitchData()
    condition_colors = NULL,     # named vector: c("CT" = "#E41A1C", "AD" = "#377EB8")
    condition_labels = NULL,     # named vector: c("CTRL" = "Control", "AD" = "Alzheimer's")
    alphas = c(0.05, 0.001),
    addErrorbars = TRUE,
    confidenceInterval = 0.95,
    localTheme = theme_bw(base_size = 13),
    title = NULL
) {
    df <- data$gene_expression
    df$Condition <- factor(df$Condition, levels = c(data$condition1, data$condition2))
    
    # Default colors
    if (is.null(condition_colors)) {
        condition_colors <- setNames(c("darkgrey", "#333333"), c(data$condition1, data$condition2))
    }
    
    # Rename condition labels if provided
    if (!is.null(condition_labels)) {
        old_levels <- levels(df$Condition)
        new_levels <- ifelse(old_levels %in% names(condition_labels),
                             condition_labels[old_levels], old_levels)
        levels(df$Condition) <- new_levels
        old_cnames <- names(condition_colors)
        new_cnames <- ifelse(old_cnames %in% names(condition_labels),
                             condition_labels[old_cnames], old_cnames)
        names(condition_colors) <- new_cnames
    }
    
    # CI calculation
    if (addErrorbars) {
        repNrs <- data$nrReplicates
        n1 <- repNrs$nrReplicates[repNrs$condition == data$condition1]
        n2 <- repNrs$nrReplicates[repNrs$condition == data$condition2]
        df$nrRep <- c(n1, n2)
        df$CI <- df$StdErr * qnorm(confidenceInterval / 2 + 0.5)
        df$CI_up <- df$Expression + df$CI
        df$CI_down <- pmax(0, df$Expression - df$CI)
    }
    
    # Significance
    sig <- .sigLabel(data$gene_q_value, alphas)
    if (!is.na(data$gene_q_value)) {
        if (addErrorbars) {
            sigY <- max(c(df$Expression, df$CI_up), na.rm = TRUE) * 1.05
        } else {
            sigY <- max(df$Expression, na.rm = TRUE) * 1.05
        }
    }
    
    # Plot
    g <- ggplot(df, aes(x = Condition, y = Expression, fill = Condition)) +
        geom_bar(stat = "identity", position = "dodge", width = 0.7) +
        scale_fill_manual(values = condition_colors) +
        localTheme +
        theme(legend.position = "none")
    
    if (addErrorbars) {
        g <- g + geom_errorbar(aes(ymin = CI_down, ymax = CI_up), width = 0.2)
    }
    
    # Add significance bracket
    if (!is.na(data$gene_q_value) && sig != "NA") {
        g <- g +
            geom_segment(aes(x = 1, xend = 2, y = sigY, yend = sigY), inherit.aes = FALSE) +
            annotate("text", x = 1.5, y = sigY, label = sig, vjust = -0.3,
                     size = localTheme$text$size * 0.3)
    }
    
    if (is.null(title)) title <- paste0("Gene Expression: ", data$gene_name)
    g <- g + labs(x = "Condition", y = "Gene Expression", title = title)
    
    return(g)
}


###############################################################################
# 4. plotIsoformExpression: Isoform-level expression barplot
###############################################################################
plotIsoformExpression <- function(
    data,                        # output from extractSwitchData()
    condition_colors = NULL,     # named vector
    condition_labels = NULL,     # named vector: c("CTRL" = "Control", "AD" = "Alzheimer's")
    alphas = c(0.05, 0.001),
    addErrorbars = TRUE,
    confidenceInterval = 0.95,
    localTheme = theme_bw(base_size = 13),
    title = NULL
) {
    df <- data$isoform_expression
    if (!is.null(data$isoform_order)) {
        df$isoform_id <- factor(df$isoform_id, levels = data$isoform_order)
    } else {
        df$isoform_id <- factor(df$isoform_id, levels = unique(df$isoform_id))
    }
    df$Condition <- factor(df$Condition, levels = c(data$condition1, data$condition2))
    
    if (is.null(condition_colors)) {
        condition_colors <- setNames(c("darkgrey", "#333333"), c(data$condition1, data$condition2))
    }
    
    # CI (must be done BEFORE renaming labels, since nrReplicates uses original names)
    if (addErrorbars) {
        repNrs <- data$nrReplicates
        df <- merge(df, repNrs, by.x = "Condition", by.y = "condition")
        df$Condition <- factor(df$Condition, levels = c(data$condition1, data$condition2))
        df$CI <- df$StdErr * qnorm(confidenceInterval / 2 + 0.5)
        df$CI_up <- df$Expression + df$CI
        df$CI_down <- pmax(0, df$Expression - df$CI)
    }
    
    # Rename condition labels if provided (AFTER merge)
    if (!is.null(condition_labels)) {
        old_levels <- levels(df$Condition)
        new_levels <- ifelse(old_levels %in% names(condition_labels),
                             condition_labels[old_levels], old_levels)
        levels(df$Condition) <- new_levels
        old_cnames <- names(condition_colors)
        new_cnames <- ifelse(old_cnames %in% names(condition_labels),
                             condition_labels[old_cnames], old_cnames)
        names(condition_colors) <- new_cnames
    }
    
    # Significance per isoform
    sigDF <- data$isoform_q_values
    sigDF$sig <- .sigLabel(sigDF$q_value, alphas)
    
    # Calculate y-position for significance labels
    sigDF <- merge(sigDF, 
        aggregate(
            if (addErrorbars) df$CI_up else df$Expression,
            by = list(isoform_id = df$isoform_id),
            FUN = max, na.rm = TRUE
        ),
        by = "isoform_id"
    )
    colnames(sigDF)[ncol(sigDF)] <- "ymax"
    sigDF$ymax <- sigDF$ymax * 1.05
    
    g <- ggplot(df, aes(x = isoform_id, y = Expression, fill = Condition)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
        scale_fill_manual(values = condition_colors) +
        scale_x_discrete(labels = data$isoform_short_labels) +
        localTheme +
        theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))
    
    if (addErrorbars) {
        g <- g + geom_errorbar(
            aes(ymin = CI_down, ymax = CI_up, group = Condition),
            position = position_dodge(width = 0.9), width = 0.2
        )
    }
    
    # Add sig labels
    sigDF2 <- sigDF[!is.na(sigDF$q_value) & sigDF$sig != "NA", ]
    if (nrow(sigDF2) > 0) {
    sigDF2$isoform_id <- factor(sigDF2$isoform_id, levels = levels(df$isoform_id))
        sigDF2$x <- as.numeric(sigDF2$isoform_id)
        sigDF2$y_line <- sigDF2$ymax * 0.985
        sigDF2$y_label <- sigDF2$ymax
        g <- g +
            geom_text(data = sigDF2, aes(x = isoform_id, y = y_label, label = sig),
                      inherit.aes = FALSE, vjust = -0.2, size = localTheme$text$size * 0.3) +
            geom_segment(data = sigDF2,
                         aes(x = x - 0.25, xend = x + 0.25, y = y_line, yend = y_line),
                         inherit.aes = FALSE)
    }
    
    if (is.null(title)) title <- paste0("Isoform Expression: ", data$gene_name)
    g <- g + labs(x = "Isoform", y = "Isoform Expression", title = title)
    
    return(g)
}


###############################################################################
# 5. plotIsoformUsage: Isoform usage (IF) barplot
###############################################################################
plotIsoformUsage <- function(
    data,                        # output from extractSwitchData()
    condition_colors = NULL,     # named vector
    condition_labels = NULL,     # named vector: c("CTRL" = "Control", "AD" = "Alzheimer's")
    alphas = c(0.05, 0.001),
    localTheme = theme_bw(base_size = 13),
    title = NULL
) {
    df <- data$isoform_usage
    if (!is.null(data$isoform_order)) {
        df$isoform_id <- factor(df$isoform_id, levels = data$isoform_order)
    } else {
        df$isoform_id <- factor(df$isoform_id, levels = unique(df$isoform_id))
    }
    df$Condition <- factor(df$Condition, levels = c(data$condition1, data$condition2))
    
    if (is.null(condition_colors)) {
        condition_colors <- setNames(c("darkgrey", "#333333"), c(data$condition1, data$condition2))
    }
    
    # Rename condition labels if provided
    if (!is.null(condition_labels)) {
        old_levels <- levels(df$Condition)
        new_levels <- ifelse(old_levels %in% names(condition_labels),
                             condition_labels[old_levels], old_levels)
        levels(df$Condition) <- new_levels
        old_cnames <- names(condition_colors)
        new_cnames <- ifelse(old_cnames %in% names(condition_labels),
                             condition_labels[old_cnames], old_cnames)
        names(condition_colors) <- new_cnames
    }
    
    # Significance per isoform
    sigDF <- data$switch_q_values
    sigDF$sig <- .sigLabel(sigDF$q_value, alphas)
    
    sigDF <- merge(sigDF,
        aggregate(df$IF, by = list(isoform_id = df$isoform_id), FUN = max, na.rm = TRUE),
        by = "isoform_id"
    )
    colnames(sigDF)[ncol(sigDF)] <- "ymax"
    sigDF$ymax <- sigDF$ymax * 1.05
    
    g <- ggplot(df, aes(x = isoform_id, y = IF, fill = Condition)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
        scale_fill_manual(values = condition_colors) +
        scale_x_discrete(labels = data$isoform_short_labels) +
        localTheme +
        theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1))
    
    # Add sig labels
    sigDF2 <- sigDF[!is.na(sigDF$q_value) & sigDF$sig != "NA", ]
    if (nrow(sigDF2) > 0) {
    sigDF2$isoform_id <- factor(sigDF2$isoform_id, levels = levels(df$isoform_id))
        sigDF2$x <- as.numeric(sigDF2$isoform_id)
        sigDF2$y_line <- sigDF2$ymax * 0.985
        sigDF2$y_label <- sigDF2$ymax
        g <- g +
            geom_text(data = sigDF2, aes(x = isoform_id, y = y_label, label = sig),
                      inherit.aes = FALSE, vjust = -0.15, size = localTheme$text$size * 0.3) +
            geom_segment(data = sigDF2,
                         aes(x = x - 0.25, xend = x + 0.25, y = y_line, yend = y_line),
                         inherit.aes = FALSE)
    }
    
    if (is.null(title)) title <- paste0("Isoform Usage: ", data$gene_name)
    g <- g + labs(x = "Isoform", y = "Isoform Fraction (IF)", title = title)
    
    return(g)
}


###############################################################################
# 6. customSwitchPlot: Combine all panels into one figure
#    Top:    transcript structure
#    Bottom: gene exp | isoform exp | isoform usage
###############################################################################
customSwitchPlot <- function(
    switchAnalyzeRlist,
    gene,
    condition1,
    condition2,
    condition_colors = NULL,       # e.g. c("CT" = "#E41A1C", "AD" = "#377EB8")
    condition_labels = NULL,       # e.g. c("CTRL" = "Control", "AD" = "Alzheimer's")
    transcript_colors = NULL,      # named vector: c("ENST000001" = "#E41A1C", "ENST000002" = "#377EB8")
    IFcutoff = 0.05,
    dIFcutoff = 0.1,
    alphas = c(0.05, 0.001),
    addErrorbars = TRUE,
    rescaleTranscripts = TRUE,
    reverseMinus = FALSE,
    plotORF = NULL,
    plotTopology = TRUE,
    localTheme = theme_bw(base_size = 11),
    transcript_height = 3,         # relative height of transcript panel
    expression_height = 2          # relative height of expression panels
) {
    # Extract data
    dat <- extractSwitchData(switchAnalyzeRlist, gene, condition1, condition2, IFcutoff)
    
    # Default colors if not supplied
    if (is.null(condition_colors)) {
        condition_colors <- setNames(c("#4393C3", "#D6604D"), c(condition1, condition2))
    }
    
    # Build individual plots
    p_transcript <- plotTranscriptStructure(
        switchAnalyzeRlist, gene, condition1, condition2,
        IFcutoff = IFcutoff, dIFcutoff = dIFcutoff,
        rescaleTranscripts = rescaleTranscripts,
        reverseMinus = reverseMinus,
        plotORF = plotORF,
        plotTopology = plotTopology,
        transcript_colors = transcript_colors,
        localTheme = localTheme
    )

    # Keep isoform x-axis labels in expression/usage panels synchronized with
    # transcript labels from the top panel (important when ORF annotation adds
    # extra transcript metadata in labels/order).
    tx_name_map <- attr(p_transcript, "transcript_name_map")
    tx_order_ids <- attr(p_transcript, "transcript_order_ids")
    if (!is.null(tx_name_map) && !is.null(dat$isoform_short_labels)) {
        common_ids <- intersect(names(dat$isoform_short_labels), names(tx_name_map))
        if (length(common_ids)) {
            dat$isoform_short_labels[common_ids] <- tx_name_map[common_ids]
        }
    }
    if (!is.null(tx_order_ids)) {
        dat$isoform_order <- tx_order_ids[tx_order_ids %in% dat$isoform_ids]
    }
    
    p_gene_exp <- plotGeneExpression(
        dat, condition_colors = condition_colors,
        condition_labels = condition_labels,
        alphas = alphas, addErrorbars = addErrorbars,
        localTheme = localTheme, title = "Gene Exp"
    )
    
    p_iso_exp <- plotIsoformExpression(
        dat, condition_colors = condition_colors,
        condition_labels = condition_labels,
        alphas = alphas, addErrorbars = addErrorbars,
        localTheme = localTheme, title = "Isoform Exp"
    ) + theme(legend.position = "none")
    
    p_iso_usage <- plotIsoformUsage(
        dat, condition_colors = condition_colors,
        condition_labels = condition_labels,
        alphas = alphas,
        localTheme = localTheme, title = "Isoform Usage"
    )
    
    # Title (use renamed labels if provided)
    lbl1 <- if (!is.null(condition_labels)) condition_labels[condition1] else condition1
    lbl2 <- if (!is.null(condition_labels)) condition_labels[condition2] else condition2
    title_grob <- textGrob(
        paste0("Isoform Switch: ", dat$gene_name, " (", lbl1, " vs ", lbl2, ")"),
        gp = gpar(fontsize = 14, fontface = "bold")
    )
    
    # Assemble using gridExtra (align grob heights to equalize panel heights)
    g_gene <- ggplotGrob(p_gene_exp)
    g_iso <- ggplotGrob(p_iso_exp)
    g_usage <- ggplotGrob(p_iso_usage)
    max_heights <- unit.pmax(g_gene$heights, g_iso$heights, g_usage$heights)
    g_gene$heights <- max_heights
    g_iso$heights <- max_heights
    g_usage$heights <- max_heights

    bottom_row <- arrangeGrob(
        g_gene, g_iso, g_usage,
        ncol = 3, widths = c(1, 2, 2)
    )
    
    combined <- arrangeGrob(
        title_grob,
        p_transcript,
        bottom_row,
        nrow = 3,
        heights = c(0.3, transcript_height, expression_height)
    )
    
    grid.newpage()
    grid.draw(combined)
    
    # Return invisibly so user can save with ggsave or further modify
    invisible(list(
        transcript = p_transcript,
        transcript_name_map = attr(p_transcript, "transcript_name_map"),
        transcript_color_map = attr(p_transcript, "transcript_color_map"),
        gene_expression = p_gene_exp,
        isoform_expression = p_iso_exp,
        isoform_usage = p_iso_usage,
        combined = combined,
        data = dat
    ))
}

