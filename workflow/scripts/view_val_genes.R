#!/usr/bin/env Rscript
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

input_file  <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]

# Read RDS file
data <- readRDS(input_file)

# Convert character vector to data frame with column name 'ext_gene'
data_df <- data.frame(val_gene = data, stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

syn_list <- lapply(data_df$val_gene, function(sym) {
  # 1) Zuerst versuchen: sym als offizielles Symbol
  res_sym <- tryCatch(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys    = sym,
      keytype = "SYMBOL",
      columns = c("SYMBOL", "ALIAS")
    ),
    error = function(e) NULL
  )

  # Falls kein Treffer: sym als ALIAS interpretieren
  if (is.null(res_sym) || nrow(res_sym) == 0) {
    res_alias <- tryCatch(
      AnnotationDbi::select(
        org.Hs.eg.db,
        keys    = sym,
        keytype = "ALIAS",
        columns = c("SYMBOL")
      ),
      error = function(e) NULL
    )
    if (is.null(res_alias) || nrow(res_alias) == 0) {
      return(NA_character_)
    }

    # 2) Das/die zugehörigen offiziellen Symbole holen
    symbols <- unique(na.omit(res_alias$SYMBOL))

    # 3) Für diese Symbole nochmal als SYMBOL alle Aliase holen
    res_sym <- tryCatch(
      AnnotationDbi::select(
        org.Hs.eg.db,
        keys    = symbols,
        keytype = "SYMBOL",
        columns = c("SYMBOL", "ALIAS")
      ),
      error = function(e) NULL
    )
    if (is.null(res_sym) || nrow(res_sym) == 0) {
      return(NA_character_)
    }
  }

  # Jetzt ist res_sym eine Tabelle mit SYMBOL + allen ALIAS für das/die Symbol(e)
  aliases <- unique(c(res_sym$SYMBOL, res_sym$ALIAS))
  aliases <- aliases[!is.na(aliases)]
  if (length(aliases) == 0) return(NA_character_)

  paste(sort(unique(aliases)), collapse = ",")
})

data_df$synonyms <- unlist(syn_list)

# Dubletten nach identischer Synonymliste entfernen:
# behalte nur die erste Zeile pro 'synonyms'
data_df_dedup <- data_df[!duplicated(data_df$synonyms), ]

# Optional: Auch Zeilen mit komplett NA-Synonyms entfernen
# data_df_dedup <- subset(data_df_dedup, !is.na(synonyms))

write.table(
  data_df_dedup,
  file      = output_file,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# # Write to TSV
# write.table(data_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# cat("Successfully wrote validation genes to", output_file, "\n")
