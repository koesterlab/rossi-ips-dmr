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
print(head(chipseeker_output_df), width = Inf)
# Postprocess: rename + merge back original CpG-level data
chipseeker_output_df <- chipseeker_output_df %>%
  rename(chromosome = seqnames, position = start)

cpg_annotated <- merge(chipseeker_output_df, methylation, by = c("chromosome", "position"))
print(head(cpg_annotated), width = Inf)
# Aggregate: average methylation per gene, separately per germ layer
output <- cpg_annotated %>%
  group_by(geneId) %>%
  summarise(
    chromosome = first(chromosome),
    annotation = first(annotation),
    transcriptId = first(transcriptId),
    num_CpGs = n(),
    mean_psc_methylation = mean(psc_methylation, na.rm = TRUE),
    mean_endoderm_methylation = mean(endoderm_methylation, na.rm = TRUE),
    mean_mesoderm_methylation = mean(mesoderm_methylation, na.rm = TRUE),
    mean_ectoderm_methylation = mean(ectoderm_methylation, na.rm = TRUE),
    .groups = "drop"
  )
print(head(output), width = Inf)
write.table(output, file = snakemake@output[['chipseeker']], sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
