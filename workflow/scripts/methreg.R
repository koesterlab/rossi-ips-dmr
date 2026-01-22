
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# JASPAR2022 is needed for methreg but does not work with conda
if (!requireNamespace("JASPAR2022", quietly = TRUE)) {
  BiocManager::install(
    "JASPAR2022",
    ask = FALSE,
    update = FALSE
  )
}

suppressPackageStartupMessages({
    library(MethReg)
    library(arrow)
    library(dplyr)
    library(GenomicRanges)
    library(SummarizedExperiment)
})
# calls <- snakemake@input[["calls"]]

# 4.1.1 Input DNA methylation dataset
meth_df <- read_parquet(snakemake@input[["meth"]])


meth_df <- meth_df %>%
  mutate(
    chromosome = paste0("chr", chromosome)
  )

meth_mat <- meth_df %>%
  select(
    psc_methylation,
    endoderm_methylation,
    mesoderm_methylation,
    ectoderm_methylation
  ) %>%
  as.matrix()

colnames(meth_mat) <- c("psc", "endoderm", "mesoderm", "ectoderm")
meth_mat <- meth_mat / 100
# We have to build the Grange object ourselves since make_dnam_se expects array data
gr <- GRanges(
  seqnames = meth_df$chromosome,
  ranges = IRanges(
    start = meth_df$position,
    end = meth_df$position+1
  )
)

names(gr) <- paste0(
  meth_df$chromosome, ":",
  meth_df$position, "-",
  meth_df$position + 1
)
meth_se <- SummarizedExperiment(
  assays = list(meth = meth_mat),
  rowRanges = gr
)

cat("\n\n#############################################\n")
print("METHYLATION DATA IN GRANGES OBJECT:")
meth_se
assay(meth_se)[1:4, 1:4]
SummarizedExperiment::rowRanges(meth_se)
cat("#############################################\n\n")

# 4.1.2 Input gene expression dataset
gene_exp_df <- read.table(
  snakemake@input[["gene_exp"]],
  sep = "\t",
  header = TRUE,
)

# We have more samples in our gene expression than in our meth calling. Since we need to have the same samples in both datasets, we compute the average expression per germ layer. Each gene becomes the name of the row.
gene_exp_df <- gene_exp_df %>%
  mutate(
    psc      = rowMeans(select(., starts_with("psc_")),      na.rm = TRUE),
    endoderm = rowMeans(select(., starts_with("endoderm_")), na.rm = TRUE),
    mesoderm = rowMeans(select(., starts_with("mesoderm_")), na.rm = TRUE),
    ectoderm = rowMeans(select(., starts_with("ectoderm_")), na.rm = TRUE)
  ) %>%
  select(ens_gene, psc, endoderm, mesoderm, ectoderm)
gene_exp_mat <- as.matrix(gene_exp_df[ , -1])
rownames(gene_exp_mat) <- gene_exp_df$ens_gene
gene_exp_df_se <- make_exp_se(
  exp = gene_exp_mat,
  genome = "hg38",
  verbose = FALSE
)
cat("\n\n#############################################\n")
print("GENE EXPRESSION DATA:")
gene_exp_df_se
SummarizedExperiment::rowRanges(gene_exp_df_se)[1:5,]
cat("#############################################\n\n")


# 4.1.3 Creating triplet dataset
# MethReg offer two methods for enhancer linking 1) a window-based method similar to the ones in the GWAS studies, 2) a fixed number of genes upstream and downstream of the DNA methylation loci similar to the one suggested by Yao et al. (2015), and one method for promoter linking, which maps to the gene of the promoter region.
# For the analysis of probes in gene promoter region, we recommend setting method = "genes.promoter.overlap", or method = "closest.gene". For the analysis of probes in distal regions, we recommend setting either method = "window" or method = "nearby.genes". Note that the distal analysis will be more time and resource consuming.
# Why doesn't it need the gene expression data here?
triplet.promoter <- create_triplet_distance_based(
  region = meth_se,
  target.method = "genes.promoter.overlap",
  genome = "hg38",
  target.promoter.upstream.dist.tss = 2000,
  target.promoter.downstream.dist.tss = 2000,
  motif.search.window.size = 400,
  motif.search.p.cutoff  = 1e-08,
  cores = 1  
)


triplet.promoter

print("TEST TRIPLET DATASET")
nrow(triplet.promoter)
sum(triplet.promoter$target %in% rownames(gene_exp_df_se))
sum(triplet.promoter$regionID %in% rownames(meth_se))
print("TEST TRIPLET DATASET")




library(dplyr)

# 1. Identifiziere Regionen, die bei deinem Threshold zu wenig Samples haben
# Bei 4 Samples und 0.25 Threshold bleibt nur 1 Sample pro Gruppe -> Absturz.
# Wir filtern Regionen raus, die nicht mindestens 2 unterschiedliche Werte in den Extremen haben.

valid_regions <- apply(assay(meth_se), 1, function(x) {
  x <- x[!is.na(x)]
  if(length(unique(x)) < 2) return(FALSE) # Keine Variation
  
  # Simuliere die Gruppenbildung
  low_cut <- quantile(x, 0.25, na.rm = TRUE)
  high_cut <- quantile(x, 0.75, na.rm = TRUE)
  
  n_low <- sum(x <= low_cut)
  n_high <- sum(x >= high_cut)
  
  # Wilcoxon braucht n > 1 pro Gruppe
  return(n_low > 1 && n_high > 1)
})

# 2. Filter die Namen der validen Regionen
keep_region_ids <- names(valid_regions)[valid_regions]

# 3. Triplets filtern
triplet.filtered <- triplet.promoter %>% 
  filter(regionID %in% keep_region_ids)
tail(triplet.filtered)
tail(triplet.promoter)
cat(paste0("Regionen vorher: ", nrow(meth_se), "\n"))
cat(paste0("Regionen nach Filter: ", length(keep_region_ids), "\n"))
cat(paste0("Triplets nach Filter: ", nrow(triplet.filtered), "\n"))




















# 4.2.1 Analysis using model with methylation by TF interaction
results.interaction.model <- interaction_model(
    triplet = triplet.filtered, 
    dnam = meth_se,
    exp = gene_exp_df_se,
    dnam.group.threshold = 0.5,
    sig.threshold = 0.5,
    fdr = T,
    stage.wise.analysis = FALSE,
    filter.correlated.tf.exp.dnam = F,
    filter.triplet.by.sig.term = T
)

results.interaction.model





















# ############################################# Example
# load(snakemake@input[["gene_exp_official"]])
# load(snakemake@input[["meth_official"]])
# meth_vals <- dna.met.chr21
# meth_vals[1:5,1:5]


# meth_vals_se <- make_dnam_se(
#   dnam = meth_vals,
#   genome = "hg38",
#   arrayType = "450k",
#   betaToM = FALSE, # transform beta to m-values 
#   verbose = FALSE # hide informative messages
# )
# meth_vals_se
# SummarizedExperiment::rowRanges(meth_vals_se)[1:4,1:4]


# gene_vals <- gene.exp.chr21.log2
# gene_vals[1:5,1:5]

# gene_vals_se <- make_exp_se(
#   exp = gene_vals,
#   genome = "hg38",
#   verbose = FALSE
# )
# gene_vals_se
# SummarizedExperiment::rowRanges(gene_vals_se)[1:5,]


# # Map probes to 5 genes upstream and 5 downstream
# triplet.distal.nearby.genes <- create_triplet_distance_based(
#   region = meth_vals_se,
#     genome = "hg38", 
#     target.method = "nearby.genes",
#     target.num.flanking.genes = 5,
#     target.window.size = 500 * 10^3,
#     target.rm.promoter.regions.from.distal.linking = TRUE,
#     motif.search.window.size = 400,
#     motif.search.p.cutoff  = 1e-08,
#     cores = 1  
# )

# str(triplet.distal.nearby.genes)