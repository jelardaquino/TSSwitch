library(bambu)
library(GenomicRanges)
library(SummarizedExperiment)
load("/home/AD.UNLV.EDU/Shared_Data/AlternativeSplicing/syn52065646_LongRead/pychopperfq/minimap2/bam_mapq10/all_bambu_Rscript/bambu.RData")

se <- bambu_results
head(assays(se)$counts)
rowRanges(se)
tx_length <- width(rowRanges(se))
tx_length_vec <- sum(tx_length)

names(tx_length_vec) <- names(tx_length)

tx_counts <- assays(se)$counts

counts_to_tpm <- function(counts, lengths) {
  lengths_kb <- lengths / 1000
  rate <- counts / lengths_kb
  sweep(rate, 2, colSums(rate), FUN = "/") * 1e6
}

tx_tpm <- counts_to_tpm(tx_counts, tx_length_vec)
write.table(tx_tpm, file = "bambu_tx_tpm.txt", sep = "\t", quote = FALSE)
