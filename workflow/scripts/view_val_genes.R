#!/usr/bin/env Rscript

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript view_val_genes.R <input.rds> <output.tsv>")
}

input_file <- args[1]
output_file <- args[2]

# Read RDS file
data <- readRDS(input_file)

# Write to TSV
write.table(data, file = output_file, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

cat("Successfully wrote validation genes to", output_file, "\n")
