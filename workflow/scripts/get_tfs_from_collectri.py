import decoupler as dc

net = dc.op.collectri(organism="human")
# net <- decoupleR::get_collectri(
#   organism = 9606,  # Homo sapiens
#   split_complexes = FALSE
# )

# net
net.to_csv(snakemake.output[0])

# write.csv(net, file = snakemake@output[[1]], row.names = FALSE)

# # Compute TF target genes contributions by effect size alone without using the ulm of decoupleR
# min_contribution = snakemake@params[["min_contribution"]] # Minimum igned_pi_value to consider a TF-target gene interaction
# gene_tf_summary <- target_genes_by_effect(max_p_value = max_p_value, min_contribution = min_contribution)
# write_csv(gene_tf_summary, snakemake@output[[1]])


# # Compute TF target genes contributions using only significant TFs from the ulm of decoupleR
# minsize = snakemake@params[["minsize"]] # Minimal size of TF affected target genes to be considered
# gene_tf_summary <- target_genes_by_significant_tfs(minsize = minsize, max_p_value = max_p_value, min_contribution = min_contribution)
# write_csv(gene_tf_summary, snakemake@output[[2]])
