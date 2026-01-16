library(decoupleR)
library(tidyverse)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)


# Transform diffexp data into matrices for decoupleR
prepare_diffexp <- function() {
  diffexp <- read_tsv(snakemake@input[[1]])

  diffexp_count <- diffexp %>%
    dplyr::filter(!is.na(ext_gene)) %>%
    dplyr::mutate(across(where(is.numeric), ~ if_else(is.na(.x), 0, .x))) %>%
    dplyr::group_by(ext_gene) %>%
    # TODO: duplicate 'row.names' are not allowed. Right now we compute means
    dplyr::summarise(across(all_of(samples), mean), .groups = "drop") %>%
    tibble::column_to_rownames("ext_gene") %>%
    as.matrix()

    deg <- diffexp %>%
    dplyr::select(
      ext_gene,
      logFC = all_of(paste0("b_condition", layer)),
      t = all_of(paste0("signed_pi_value_condition", layer)),
      P.Value = all_of("pval")
    ) %>%
    dplyr::filter(!is.na(t) & !is.na(ext_gene)) %>%
    dplyr::group_by(ext_gene) %>%
    # TODO: duplicate 'row.names' are not allowed. Right now we compute means
    dplyr::summarise(
      logFC = mean(logFC, na.rm = TRUE),
      t = mean(t, na.rm = TRUE),
      P.Value = mean(P.Value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tibble::column_to_rownames("ext_gene") %>%
    as.matrix()
    diffexp <- list(deg = deg, diffexp_count = diffexp_count)
  return(diffexp)
}


target_genes_by_effect <- function(max_p_value, min_contribution) {
  tf_gene_table <- net %>%
    inner_join(
      as.data.frame(deg) %>%
        rownames_to_column("target"),
      by = "target"
    ) %>%
    mutate(
      contribution = mor * t,
      direction = case_when(
        contribution > 0 ~ "supports_TF_activity",
        contribution < 0 ~ "opposes_TF_activity",
        TRUE ~ "neutral"
      )
    ) %>%
    arrange(source, desc(abs(contribution))) %>%
    filter(
      P.Value < max_p_value,
      abs(contribution) > min_contribution
    )

  gene_tf_summary <- tf_gene_table %>%
    group_by(target) %>%
    summarise(
      TFs = paste(sort(unique(source)), collapse = ", "),
      n_TFs = n_distinct(source),
      .groups = "drop"
    )
  return(gene_tf_summary)
}

target_genes_by_significant_tfs <- function(minsize, max_p_value, min_contribution) {
    # Run ULM to find active TFs (contrast analysis)
  contrast_acts <- decoupleR::run_ulm(
    mat = deg[, 't', drop = FALSE],
    net = net,
    .source = 'source',
    .target = 'target',
    .mor = 'mor',
    minsize = minsize
  )

  # Filter significant TFs
  sig_TFs <- contrast_acts %>%
    filter(p_value < max_p_value) %>%
    pull(source)

  # Get target genes of significant TFs
  tf_targets <- net %>%
    filter(source %in% sig_TFs) %>%
    select(source, target, mor)

  # Combine with DEG statistics and calculate contribution / direction
  tf_gene_table <- tf_targets %>%
    inner_join(
      as.data.frame(deg) %>% rownames_to_column("target"),
      by = "target"
    ) %>%
    mutate(
      contribution = mor * t,
      direction = case_when(
        contribution > 0 ~ "supports_TF_activity",
        contribution < 0 ~ "opposes_TF_activity",
        TRUE ~ "neutral"
      )
    ) %>%
    arrange(source, desc(abs(contribution)))

  # Summarise TFs per target gene
  gene_tf_summary <- tf_gene_table %>%
    group_by(target) %>%
    summarise(
      TFs = paste(sort(unique(source)), collapse = ", "),
      n_TFs = n_distinct(source),
      .groups = "drop"
    )
  return(gene_tf_summary)
}


# Maximum p-value to consider a TF significant
max_p_value = snakemake@params[["max_p_value"]]

# Compute samples for the given layer
layer = snakemake@params[["layer"]]
samples <- snakemake@input[["samples"]] %>%
  read_tsv() %>%
  filter(condition %in% c(layer, "psc")) %>%
  pull(sample) %>%
  paste0("-1")

# Prepare differential expression data
diffexp <- prepare_diffexp()
deg <- diffexp$deg
diffexp_count <- diffexp$diffexp_count

# Download possible TFs to target genes from CollecTRI via decoupleR
net <- decoupleR::get_collectri(organism = 'human', 
                                split_complexes = FALSE)


# Compute TF target genes contributions by effect size alone without using the ulm of decoupleR

min_contribution = snakemake@params[["min_contribution"]] # Minimum signed_pi_value to consider a TF-target gene interaction
gene_tf_summary <- target_genes_by_effect(max_p_value = max_p_value, min_contribution = min_contribution)
write_csv(gene_tf_summary, snakemake@output[[1]])


# Compute TF target genes contributions using only significant TFs from the ulm of decoupleR

minsize = snakemake@params[["minsize"]] # Minimal size of TF affected target genes to be considered
gene_tf_summary <- target_genes_by_significant_tfs(minsize = minsize, max_p_value = max_p_value, min_contribution = min_contribution)
write_csv(gene_tf_summary, snakemake@output[[2]])

