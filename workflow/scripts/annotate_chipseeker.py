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
merged_df = pd.merge(chipseeker_df, ensembl_genes_df, on="transcriptId", how="inner")

# Drop rows with missing values in the specified columns
filtered_df = merged_df.dropna(
    subset=["transcriptId", "ensembl_gene_id", "external_gene_name"]
)

filtered_df = filtered_df.rename(
    columns={
        "q_value": "qval",
        "ensembl_gene_id": "ens_gene",
        "external_gene_name": "ext_gene",
        "p_MWU": "pval",
    }
)

# Save the processed data to the output file
filtered_df.to_csv(output_file, sep="\t", index=False)
