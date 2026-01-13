import pandas as pd
import sys
import numpy as np
import openpyxl

# sys.stderr = open(snakemake.log[0], "w", buffering=1)

germ_layer = snakemake.params["germ_layer"]
top_n = snakemake.params.get("top_n_genes", 10)
germ_layer = germ_layer.replace("derm", "")

genes_df = pd.read_csv(snakemake.input["genes_transcripts"], sep="\t")
genes_df = genes_df[genes_df["annotation_type"] == "Promoter"]
log2FC = f"log2fc_{germ_layer}_avg_combined"

##### Jochens RNA data
rna_df = pd.read_excel(snakemake.input["rna_seq"])
rna_df = rna_df[rna_df["test"] == germ_layer].reset_index(drop=True)
# Compute signed pi-value for rna_df
rna_df["absolute_signed_pi_value"] = np.abs(-np.log10(rna_df["FDR"] * rna_df["logFC"]))
rna_df = rna_df.sort_values(by="absolute_signed_pi_value", ascending=False)
# Filter to top N genes
# rna_df = rna_df.head(top_n)
rna_df = rna_df[["genes", "absolute_signed_pi_value", "logFC", "FDR"]].rename(
    columns={
        "genes": "ext_gene",
        "absolute_signed_pi_value": "rna_signed_pi",
        "logFC": "rna_log2FC",
        "FDR": "rna_FDR",
    }
)


##### DMR data
dmr_df = pd.read_csv(snakemake.input["genes_transcripts"], sep="\t")
dmr_df = dmr_df[dmr_df["annotation_type"] == "Promoter"]
dmr_df = dmr_df[
    ["ext_gene", "absolute_signed_pi_val", "mean_methylation_difference", "qval"]
].rename(
    columns={
        "absolute_signed_pi_val": "dmr_signed_pi",
        "mean_methylation_difference": "dmr_mean_methylation_difference",
        "qval": "dmr_qval",
    }
)


###### Diffexp data
diffexp_df = pd.read_csv(snakemake.input["diffexp"], sep="\t")
diffexp_df = diffexp_df.sort_values(
    by=f"signed_pi_value_condition{germ_layer}derm", ascending=False
)
diffexp_df = diffexp_df[
    [
        "ext_gene",
        f"signed_pi_value_condition{germ_layer}derm",
        "qval",
        f"b_condition{germ_layer}derm",
    ]
].rename(
    columns={
        f"signed_pi_value_condition{germ_layer}derm": "diffexp_signed_pi",
        "qval": "diffexp_qval",
        f"b_condition{germ_layer}derm": "diffexp_log2FC",
    }
)

# Rank dfs
# Sort + index-based ranking
rna_df = rna_df.sort_values("rna_signed_pi", ascending=False).reset_index(drop=True)
rna_df["rna_rank"] = rna_df.index + 1

dmr_df = dmr_df.sort_values("dmr_signed_pi", ascending=False).reset_index(drop=True)
dmr_df["dmr_rank"] = dmr_df.index + 1

diffexp_df = diffexp_df.sort_values("diffexp_signed_pi", ascending=False).reset_index(
    drop=True
)
diffexp_df["diffexp_rank"] = diffexp_df.index + 1


merged_all = rna_df.merge(dmr_df, on="ext_gene", how="outer").merge(
    diffexp_df, on="ext_gene", how="outer"
)


final_genes = merged_all[
    (merged_all["rna_rank"] <= top_n)
    | (merged_all["dmr_rank"] <= top_n)
    | (merged_all["diffexp_rank"] <= top_n)
].copy()


final_genes = final_genes.sort_values(
    by=["rna_rank", "dmr_rank", "diffexp_rank"]
).reset_index(drop=True)
print(final_genes)
final_genes.to_excel(snakemake.output["comparison"], index=False)
final_genes["ext_gene"].to_csv(snakemake.output["top_genes"], index=False)
