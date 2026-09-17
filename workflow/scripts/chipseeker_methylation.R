log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")
library(GenomicRanges)
library(ChIPseeker)
library(dplyr)
library(AnnotationDbi)
library(arrow)

# Read methylation data from parquet
methylation <- read_parquet(snakemake@input[['methylation']])
methylation <- methylation[order(methylation$chromosome, methylation$position), ]

# Create GRanges Object from single CpG positions (start == end)
gr <- GRanges(seqnames = Rle(methylation$chromosome),
              ranges = IRanges(start = methylation$position, end = methylation$position),
              strand = Rle("*"))

# Annotate with chipseeker
txdb <- loadDb(snakemake@input[["txdb"]])
txnames <- readRDS(snakemake@input[["txnames"]])
chipseeker_output <- annotatePeak(gr, TxDb = txdb)
chipseeker_output_df <- as.data.frame(chipseeker_output)
# Postprocess: rename + merge back original CpG-level data
chipseeker_output_df <- chipseeker_output_df %>%
  rename(chromosome = seqnames, position = start)

cpg_annotated <- merge(chipseeker_output_df, methylation, by = c("chromosome", "position"))
# Aggregate: average methylation per gene, separately per germ layer

write.table(cpg_annotated, file = snakemake@output[['chipseeker']], sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
