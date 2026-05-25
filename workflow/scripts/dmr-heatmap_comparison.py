import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage

# We need this to cluster huge data
sys.stderr = open(snakemake.log[0], "w", buffering=1)
pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)
sys.setrecursionlimit(100000)
base = snakemake.params["base"]

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


pacbio_input_files = snakemake.input.pacbio
nanopore_input_files = snakemake.input.nanopore
output = snakemake.output[0]
dfs = []
for input_files in [pacbio_input_files, nanopore_input_files]:
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
        heatmap_data = heatmap_data.merge(
            df, on=["region", "annotation_type"], how="outer"
        )
    # Drop rows with all NaN values
    # heatmap_data = heatmap_data.dropna()

    name = os.path.basename(output).replace(".png", "")
    annotation_type = filename_to_name[name]
    # df_filtered = heatmap_data[heatmap_data["annotation_type"] == annotation_type]
    df_filtered = heatmap_data[sample_names + ["region"]]
    df_filtered = df_filtered.replace([np.inf, -np.inf], np.nan)
    dfs.append(df_filtered)

dfs[0]["region"] = dfs[0]["region"].astype(str)
dfs[1]["region"] = dfs[1]["region"].astype(str)
number_pacbio_genes = dfs[0].shape[0]
number_nanopore_genes = dfs[1].shape[0]
df_complete = pd.merge(
    dfs[0], dfs[1], on="region", how="outer", suffixes=("_pacbio", "_nanopore")
)
# Count non-NaN values per row to sort by data completeness
pacbio_cols = [col for col in df_complete.columns if col.endswith("_pacbio")]
nanopore_cols = [col for col in df_complete.columns if col.endswith("_nanopore")]


# df_complete = df_complete.drop("non_nan_count", axis=1)
df_complete = df_complete.set_index("region")

vmin = df_complete.min().min()
vmax = df_complete.max().max()

fig, axes = plt.subplots(1, 3, figsize=(15, 6))

layers = ["endoderm", "mesoderm", "ectoderm"]

for idx, layer in enumerate(layers):
    number_nanopore_genes = df_complete[f"{layer}_nanopore"].notna().sum()
    number_pacbio_genes = df_complete[f"{layer}_pacbio"].notna().sum()
    number_common_genes = (
        df_complete[[f"{layer}_nanopore", f"{layer}_pacbio"]].dropna().shape[0]
    )
    print(layer, number_nanopore_genes, number_pacbio_genes, number_common_genes)
    df_temp = df_complete[[f"{layer}_nanopore", f"{layer}_pacbio"]].dropna()
    df_temp = df_temp.rename(
        columns={f"{layer}_nanopore": "nanopore", f"{layer}_pacbio": "pacbio"}
    )
    if len(df_temp) > 1:
        # Clustering durchführen
        row_linkage = linkage(df_temp, method="ward")
        col_linkage = linkage(df_temp.T, method="ward")

        # Daten nach Clustering sortieren
        row_order = dendrogram(row_linkage, no_plot=True)["leaves"]
        df_sorted = df_temp.iloc[row_order]
    else:
        df_sorted = df_temp
    corr = df_sorted["nanopore"].corr(df_sorted["pacbio"])
    corr_pear = df_sorted.corr(method="pearson")
    covariance = df_sorted.cov()["nanopore"]["pacbio"]
    variance_nanopore = df_sorted["nanopore"].var()
    variance_pacbio = df_sorted["pacbio"].var()
    own_corr = covariance / ((variance_nanopore**0.5) * (variance_pacbio**0.5))
    print(layer)
    print("COrr", corr)
    print("Pearson", corr_pear)
    print("Own corr", own_corr)
    # Heatmap zeichnen mit gemeinsamer Skala
    sns.heatmap(
        df_sorted,
        cmap="vlag_r",
        center=0,
        vmin=vmin,
        vmax=vmax,
        ax=axes[idx],
        cbar=(idx == 2),  # Nur rechts eine Colorbar
        cbar_kws={"label": "Value"},
    )

    axes[idx].set_title(f"{layer} ({number_common_genes})")
    # X-ticks mit Nummern
    x_labels = [
        f"nanopore\n{number_nanopore_genes}",
        f"pacbio\n{number_pacbio_genes}",
    ]
    axes[idx].set_xticklabels(x_labels, fontsize=9)

    # Y-ticks verstecken
    axes[idx].set_yticks([])

    # Y-Achsen-Label mit common gene count
    # axes[idx].set_ylabel(f"Gene regions ({number_common_genes})")
    # X-ticks drehen
    # axes[idx].set_xticklabels(axes[idx].get_xticklabels(), rotation=45, ha="right")

    # # Y-ticks verstecken
    # axes[idx].set_yticks([])

    if idx == 0:
        axes[idx].set_ylabel("Gene regions")
    else:
        axes[idx].set_ylabel("")

plt.tight_layout()
plt.savefig(output, format="png", dpi=300, bbox_inches="tight")
plt.close()
