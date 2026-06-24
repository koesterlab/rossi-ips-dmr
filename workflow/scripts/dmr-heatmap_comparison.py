import os
import sys

import altair as alt

# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def bin_methylation(series: pd.Series, bin_size: int) -> pd.Series:
    """Round methylation values to nearest bin_size and cast to int."""
    return (np.round(series / bin_size) * bin_size).astype(int)


# We need this to cluster huge data
sys.stderr = open(snakemake.log[0], "w", buffering=1)
pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)
alt.data_transformers.enable("vegafusion")

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
    df_grouped = (
        df.groupby(["transcriptId", "annotation_type"])["mean_methylation_difference"]
        .mean()
        .reset_index()
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

    name = os.path.basename(output).replace(".png", "")
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

# fig, axes = plt.subplots(1, 3, figsize=(15, 6))

layers = ["endoderm", "mesoderm", "ectoderm"]
charts = []
for idx, layer in enumerate(layers):
    number_nanopore_genes = df_complete[f"{layer}_nanopore"].notna().sum()
    number_pacbio_genes = df_complete[f"{layer}_pacbio"].notna().sum()
    number_common_genes = (
        df_complete[[f"{layer}_nanopore", f"{layer}_pacbio"]].dropna().shape[0]
    )
    df_temp = df_complete[[f"{layer}_nanopore", f"{layer}_pacbio"]].dropna()
    df_temp = df_temp.rename(
        columns={f"{layer}_nanopore": "nanopore", f"{layer}_pacbio": "pacbio"}
    )

    df_sorted = df_temp
    corr = df_sorted["nanopore"].corr(df_sorted["pacbio"])
    print(layer)
    print("Pearson", corr)

    df_sorted["pacbio"] = df_sorted["pacbio"] * 100
    df_sorted["nanopore"] = df_sorted["nanopore"] * 100
    df_sorted = df_sorted.assign(
        pacbio_bin=bin_methylation(df_sorted["pacbio"], 10),
        nanopore_bin=bin_methylation(df_sorted["nanopore"], 10),
    )

    counts = (
        pd.crosstab(df_sorted["pacbio_bin"], df_sorted["nanopore_bin"])
        .stack()
        .reset_index(name="count")
    )

    all_bins = range(-100, 101, 10)

    # create complete grid
    full_index = pd.MultiIndex.from_product(
        [all_bins, all_bins],
        names=["pacbio_bin", "nanopore_bin"],
    )

    # fill missing combinations with 0
    counts = (
        counts.set_index(["pacbio_bin", "nanopore_bin"])
        .reindex(full_index, fill_value=0)
        .reset_index()
    )

    counts["pacbio_bin"] = counts["pacbio_bin"] / 100
    counts["nanopore_bin"] = counts["nanopore_bin"] / 100

    plot = (
        alt.Chart(
            counts,
            title=alt.Title(text=f"{layer}: {number_common_genes}", fontSize=20),
        )
        .mark_rect()
        .encode(
            x=alt.X(
                "pacbio_bin:Q",
                bin=alt.Bin(step=0.1),
                sort=alt.SortOrder("ascending"),
                title=alt.Title(text=f"PacBio:\n{number_pacbio_genes}", fontSize=12),
            ),
            y=alt.Y(
                "nanopore_bin:Q",
                bin=alt.Bin(step=0.1),
                sort=alt.SortOrder("ascending"),
                title=alt.Title(
                    text=f"Nanopore:\n{number_nanopore_genes}", fontSize=12
                ),
            ),
            color=alt.Color(
                "count:Q",
                scale=alt.Scale(
                    type="log", scheme="viridis", domain=[1, counts["count"].max() + 1]
                ),
            ),
        )
    ).properties(width=200, height=200)

    v_line = (
        alt.Chart(pd.DataFrame({"x": [0]}))
        .mark_rule(color="red", size=1)
        .encode(x=alt.X("x:O", title=None, axis=None))
    )
    h_line = (
        alt.Chart(pd.DataFrame({"y": [0]}))
        .mark_rule(color="red", size=1)
        .encode(y=alt.Y("y:O", title=None, axis=None))
    )
    plot = plot + v_line + h_line

    charts.append(plot)

chart = alt.hconcat(*charts).resolve_scale("shared")

chart.save(snakemake.output[0])
# plt.tight_layout()
# plt.savefig(output, format="png", dpi=300, bbox_inches="tight")
# plt.close()
