#!/usr/bin/env Rscript

# Converts Bambu transcript counts to counts per million (CPM). The output has the
# same layout as bambu_tx_tpm.txt so it can be passed to SUPPA2 psiPerEvent.

BASE_DIR <- "/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/jaquino_analysis"
COUNTS_FILE <- file.path(BASE_DIR, "NanoporeRNASeq_counts_transcript.txt")
CPM_FILE <- file.path(BASE_DIR, "bambu_tx_cpm.txt")

counts <- read.delim(COUNTS_FILE, row.names = 1, check.names = FALSE)
counts$GENEID <- NULL
counts <- as.matrix(counts)

cpm <- sweep(counts, 2, colSums(counts), FUN = "/") * 1e6

write.table(cpm, file = CPM_FILE, sep = "\t", quote = FALSE)

message("Library sizes (reads): ", paste(colnames(counts), round(colSums(counts)), sep = "=", collapse = ", "))
message("Wrote ", CPM_FILE)
