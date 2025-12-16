log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(biomaRt)
options(biomaRt.cache = FALSE)
tmp_cache_dir <- tempfile("biomaRt_cache_")
dir.create(tmp_cache_dir)
Sys.setenv(BIOMART_CACHE = tmp_cache_dir)

library(tidyverse)
library(cli)
library(R.utils)


tryCatch({
  biomaRt::biomartCacheClear()
}, error = function(e) {
  warning("Cache konnte nicht gelöscht werden: ", e$message)
})

data <- tryCatch({
  readr::read_tsv(snakemake@input[[1]], show_col_types = FALSE)
}, error = function(e) {
  stop(paste("Fehler beim Einlesen der Datei:", e$message))
})

mart <- "useast"
rounds <- 0
while (class(mart)[[1]] != "Mart") {
  mart <- tryCatch(
    {
      if (mart == "www") rounds <- rounds + 1 
        withTimeout({
          biomaRt::useEnsembl(
            biomart = "ENSEMBL_MART_ENSEMBL",
            dataset = paste0(snakemake@params[["species"]], "_gene_ensembl"),
            mirror = mart
          )
        }, timeout = 60, onTimeout = "error")
      
    },
    error = function(e) {
      if (rounds >= 3) {
        stop(paste("Alle Ensembl-Mirrors durchlaufen,", rounds, "Runden. Kein Erfolg. Letzter Fehler:", e$message))
      }
      mart <- switch(mart,
        useast = "uswest",
        uswest = "asia",
        asia = "www",
        www = {
          Sys.sleep(30)
          "useast"
        }
      )
    }
  )
}

print("Verbindung zu Ensembl erfolgreich.")

# Geninformationen abfragen
print("Frage Geninformationen ab...")
show(mart)
gene_info <- tryCatch({
  withTimeout({
    getBM(
      attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"),
      filters = "ensembl_transcript_id",
      values = data$transcriptId,
      mart = mart
    )
  }, timeout = 3600, onTimeout = "error")

}, error = function(e) {
  stop(paste("Fehler bei getBM():", e$message))
})

print(paste("getBM() erfolgreich, Zeilen:", nrow(gene_info)))

# Schreibe Ergebnis
print("Schreibe Ausgabe-Datei...")
write.table(gene_info, file = snakemake@output[[1]], sep = "\t", row.names = FALSE, quote = FALSE)

print("Script abgeschlossen.")

