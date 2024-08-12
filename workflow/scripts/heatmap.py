import pandas as pd
import altair as alt
import os
from sklearn import cluster
import seaborn as sns

filename_to_name = {
    "distal_intergenic": "Distal Intergenic",
    "promoter": "Promoter",
    "intron": "Intron",
    "exon": "Exon",
    "3_utr": "3' UTR",
    "5_utr": "5' UTR",
    "downstream": "Downstream",
}


# Group df by gene regions
def aggregate_by_gene_region(df):
    df["transcriptId"] = df["transcriptId"].astype(str)
    df["annotation"] = df["annotation"].astype(str)
    # For every gene region take the maximum of methylation differences
    df_grouped = (
        df.groupby(["transcriptId", "annotation"])["mean_methylation_difference"]
        .max()
        .reset_index()
    )
    df_grouped["annotation_type"] = df_grouped["annotation"].str.replace(
        "\s*\([^)]*\)", "", regex=True
    )
    df_grouped["region"] = df_grouped.apply(
        lambda row: f"{row['transcriptId']}:{row['annotation']}", axis=1
    )
    return df_grouped[["region", "annotation_type", "mean_methylation_difference"]]

input_files = snakemake.input
output = snakemake.output[0]

# Merge single dfs to one df
aggregated_data = [
    aggregate_by_gene_region(pd.read_csv(file, sep="\t")) for file in input_files
]
all_regions = set(region for df in aggregated_data for region in df["region"])
for idx, df in enumerate(aggregated_data):
    aggregated_data[idx] = df.rename(
        columns={"mean_methylation_difference": f"BC0{idx+2}"}
    )
heatmap_data = (
    aggregated_data[0]
    .merge(aggregated_data[1], on=["region", "annotation_type"], how="outer")
    .merge(aggregated_data[2], on=["region", "annotation_type"], how="outer")
)
heatmap_data = heatmap_data.fillna(0)


# Filter on annotation_type
name = os.path.basename(output).replace(".png", "")
annotation_type = filename_to_name[name]
df_filtered = heatmap_data[heatmap_data["annotation_type"] == annotation_type]
df_filtered = df_filtered[["BC02", "BC03", "BC04"]]

# Create heatmap with seaborn
heatmap = sns.clustermap(
    df_filtered,
    z_score=0,
    cmap="vlag_r",
    center=0,
    row_cluster=False,
    col_cluster=False,
)
heatmap.ax_heatmap.set_title(
    f"DMRs between Samples and BC01 for annotation type {annotation_type}"
)
heatmap.ax_heatmap.set_xlabel("Samples")
heatmap.ax_heatmap.set_ylabel("Generegions")
heatmap.ax_heatmap.yaxis.set_label_position("left")
heatmap.ax_heatmap.yaxis.set_ticks([])


heatmap.savefig(output, format="png")
