import pandas as pd
import sys

# Define file paths from Snakemake
chipseeker_file = snakemake.input["chipseeker"]
ensembl_genes_file = snakemake.input["genes"]
output_file = snakemake.output[0]

# Load the input files
chipseeker_df = pd.read_csv(chipseeker_file, sep="\t")
ensembl_genes_df = pd.read_csv(ensembl_genes_file, sep="\t")

# Rename columns to match for merging
ensembl_genes_df.rename(columns={"ensembl_transcript_id": "transcriptId"}, inplace=True)

# Perform the inner join
merged_df = pd.merge(ensembl_genes_df, chipseeker_df, on="transcriptId", how="inner")

# Drop rows with missing values in the specified columns
filtered_df = merged_df.dropna(
    subset=["external_gene_name", "transcriptId", "ensembl_gene_id"]
)

filtered_df["annotation_type"] = filtered_df["annotation"].str.replace(
    r"\s*\([^)]*\)", "", regex=True
)

filtered_df = filtered_df.rename(
    columns={
        "q_value": "qval",
        "ensembl_gene_id": "ens_gene",
        "external_gene_name": "ext_gene",
        "p_MWU": "pval",
    }
)

column_order = [
    "ext_gene",
    "transcriptId",
    "ens_gene",
    "annotation_type",
    "mean_methylation_difference",
    "absolute_signed_pi_val",
    "qval",
    "pval",
    "chr",
    "start_dmr",
    "end_dmr",
    "width",
    "strand",
    "annotation",
    "geneChr",
    "geneStart",
    "geneEnd",
    "geneLength",
    "geneStrand",
    "geneId",
    "distanceToTSS",
    "num_CpGs",
    "p_2D_KS",
    "mean_g1",
    "mean_g2",
]
filtered_df = filtered_df[column_order]

filtered_df = filtered_df.sort_values(by="absolute_signed_pi_val", ascending=False)

# Save the processed data to the output file
filtered_df.to_csv(output_file, sep="\t", index=False)
