import pandas as pd
import os
import sys
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

# We need this to cluster huge data
sys.stderr = open(snakemake.log[0], "w", buffering=1)
pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)
sys.setrecursionlimit(100000)

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
    # For every gene region take the mean of methylation differences
    df_grouped = (
        df.groupby(["transcriptId", "annotation"])["mean_methylation_difference"]
        .mean()
        .reset_index()
    )
    df_grouped["annotation_type"] = df_grouped["annotation"].str.replace(
        r"\s*\([^)]*\)", "", regex=True
    )
    df_grouped["region"] = df_grouped.apply(
        lambda row: f"{row['transcriptId']}:{row['annotation_type']}", axis=1
    )
    return df_grouped[["region", "annotation_type", "mean_methylation_difference"]]


input_files = snakemake.input
output = snakemake.output[0]

sample_names = [
    os.path.basename(os.path.dirname(os.path.dirname(file))) for file in input_files
]

aggregated_data = []
for file, sample_name in zip(input_files, sample_names):
    df = pd.read_csv(
        file, sep="\t", dtype={"chr": str, "transcriptId": str, "annotation": str}
    )
    agg_df = aggregate_by_gene_region(df)
    agg_df = agg_df.rename(columns={"mean_methylation_difference": sample_name})
    aggregated_data.append(agg_df)

heatmap_data = aggregated_data[0]
for df in aggregated_data[1:]:
    heatmap_data = heatmap_data.merge(df, on=["region", "annotation_type"], how="outer")
# Drop rows with all NaN values

heatmap_data = heatmap_data.dropna()


name = os.path.basename(output).replace(".png", "")
annotation_type = filename_to_name[name]
df_filtered = heatmap_data[heatmap_data["annotation_type"] == annotation_type]
df_filtered = df_filtered[sample_names]
df_filtered = df_filtered.replace([np.inf, -np.inf], np.nan).dropna()
# If df_filtered has less than 2 rows, we cannot cluster it, so we will just plot a heatmap without clustering
if df_filtered.shape[0] < 2:
    fig, ax = plt.subplots(figsize=(6, 4))

    ax.set_facecolor("#f5f5f5")
    ax.text(
        0.5,
        0.5,
        f"No common DMR of type {annotation_type} between samples",
        ha="center",
        va="center",
        fontsize=14,
        fontweight="bold",
    )

    ax.set_xticks([])
    ax.set_yticks([])

    ax.set_title(f"Annotation type: {annotation_type}")

    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()
else:

    heatmap = sns.clustermap(
        df_filtered,
        cmap="vlag_r",
        center=0,
        col_cluster=False,
    )

    heatmap.ax_heatmap.set_title(
        f"DMRs between Samples and psc for annotation type {annotation_type}"
    )
    heatmap.ax_row_dendrogram.set_visible(False)
    heatmap.ax_col_dendrogram.set_visible(False)
    heatmap.ax_heatmap.set_xlabel("Samples")
    heatmap.ax_heatmap.set_ylabel("Generegions")
    heatmap.ax_heatmap.yaxis.set_label_position("left")
    heatmap.ax_heatmap.yaxis.set_ticks([])

    heatmap.savefig(output, format="png")
