#!/usr/bin/env Rscript
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]

# Read RDS file
data <- readRDS(input_file)

print(data)
print(typeof(data))

# Convert character vector to data frame with column name 'gene_id'
data_df <- data.frame(ext_gene = data, stringsAsFactors = FALSE)

# Write to TSV
write.table(data_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

cat("Successfully wrote validation genes to", output_file, "\n")
