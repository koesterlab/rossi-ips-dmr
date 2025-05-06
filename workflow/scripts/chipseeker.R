log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(GenomicRanges)
library(ChIPseeker)
library(dplyr)
library(AnnotationDbi)

metilene <- read.table(snakemake@input[['metilene']], header=FALSE, sep="\t")
colnames(metilene) <- c("chr", "start", "end", "q_value", "mean_methylation_difference", "num_CpGs", "p_MWU", "p_2D_KS", "mean_g1", "mean_g2")
metilene <- metilene[order(metilene$chr, metilene$start), ]

# Create GRanges Object from DMRs
gr <- GRanges(seqnames = Rle(metilene$chr),
              ranges = IRanges(start = metilene$start, end = metilene$end),
              strand = Rle("*"))

# Annotate with chipseeker
txdb <- loadDb(snakemake@input[["txdb"]])
txnames <- readRDS(snakemake@input[["txnames"]])

chipseeker_output <- annotatePeak(gr, TxDb = txdb)
chipseeker_output_df <- as.data.frame(chipseeker_output)

# Postprocess output
chipseeker_output_df <- chipseeker_output_df %>%
  rename(chr = seqnames)
output <- merge(chipseeker_output_df, metilene, by = c("chr", "start", "end")) %>%
  filter(q_value <= 0.25) %>%
  mutate(
    absolute_signed_pi_val = ifelse(q_value == 0, NaN, abs(-log10(q_value) * mean_methylation_difference))
  ) %>%
  arrange(desc(absolute_signed_pi_val)) %>%
  rename(start_dmr = start) %>%
  rename(end_dmr = end) 
  # %>%
  # dplyr::select(chr, start_dmr, end_dmr, mean_methylation_difference, annotation, transcriptId, q_value, p_MWU, absolute_signed_pi_val)

write.table(output, file = snakemake@output[['chipseeker']], sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


