log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(GenomicFeatures)
library(AnnotationDbi)
library(rtracklayer)

annotation_gtf <- snakemake@input[[1]]
txdb <- makeTxDbFromGFF(annotation_gtf)

saveDb(txdb, snakemake@output[["txdb"]])
gtf <- import(annotation_gtf, "gtf", colnames=c("transcript_id", "transcript_name"))
txnames <- gtf$transcript_name
names(txnames) <- gtf$transcript_id

saveRDS(txnames, snakemake@output[["txnames"]])