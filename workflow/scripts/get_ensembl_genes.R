#  log <- file(snakemake@log[[1]], open="wt")
#  sink(log)
#  sink(log, type="message")

library(biomaRt)
library("tidyverse")
library("cli")

data <- read.table(snakemake@input[[1]], sep="\t", header=TRUE)

mart <- "useast"
rounds <- 0
while (class(mart)[[1]] != "Mart") {
  mart <- tryCatch(
    {
      # done here, because error function does not
      # modify outer scope variables, I tried
      if (mart == "www") rounds <- rounds + 1
      # equivalent to useMart, but you can choose
      # the mirror instead of specifying a host
      biomaRt::useEnsembl(
        biomart = "ENSEMBL_MART_ENSEMBL",
        dataset = str_c(snakemake@params[["species"]], "_gene_ensembl"),
        version = snakemake@params[["version"]],
        mirror = mart
      )
    },
    error = function(e) {
      # change or make configurable if you want more or
      # less rounds of tries of all the mirrors
      if (rounds >= 3) {
        cli_abort(
          str_c(
            "Have tried all 4 available Ensembl biomaRt mirrors ",
            rounds,
            " times. You might have a connection problem, or no mirror is responsive.\n",
            "The last error message was:\n",
            message(e)
          )
        )
      }
      # hop to next mirror
      mart <- switch(mart,
        useast = "uswest",
        uswest = "asia",
        asia = "www",
        www = {
          # wait before starting another round through the mirrors,
          # hoping that intermittent problems disappear
          Sys.sleep(30)
          "useast"
        }
      )
    }
  )
}

gene_info <- getBM(attributes=c('ensembl_transcript_id', 'ensembl_gene_id', 'external_gene_name'),
                   filters='ensembl_transcript_id', 
                   values=data$transcriptId, 
                   mart=mart)


# Join does not work as intended
# data_with_genes <- 
#                     inner_join(
#                     data,
#                     gene_info,
#                     by = join_by(transcriptId == ensembl_transcript_id))


write.table(gene_info, file=snakemake@output[[1]], sep="\t", row.names=FALSE, quote=FALSE)